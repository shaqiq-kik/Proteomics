#!/usr/bin/env python3
"""One-command runner for the `proteomics_de/` SILAC differential-expression pipeline.

Why this file exists
--------------------
The pipeline is 13 scripts that must run in a specific order. Until now that
order lived in comments, in `BUILD_LOG.md` prose, and in whoever last ran it.
This module makes the order **executable data**: `STAGES` below is the single
source of truth for what runs, in what sequence, what each stage produces, and
what "produced correctly" means for each output.

    .venv/bin/python run_pipeline.py --all

The design constraint that shapes everything: **emptiness is not failure.**
Three of this pipeline's outputs are header-only *on purpose* and a naive
"the file is non-empty, ship it" check would false-positive on all three:

  * `results/ipa_input_significant.csv` — 0/1938 proteins pass FDR<0.05
    (min adj.p = 0.305). The expected n=2 technical-replicate ceiling. See
    DECISIONS_LOG D2.
  * `results/enrichment/ora_up.csv`, `ora_down.csv` — 0 GO/KEGG/Reactome terms
    survive the honest detected-proteome background. See DECISIONS_LOG D6. This
    is the report's headline finding, not a bug.

So every output declares its own expected shape, and a stage resolves to one of
`PASS` / `PASS_EMPTY_EXPECTED` / `FAIL`. A table with `expected_rows == 0` that
returns 0 rows is `PASS_EMPTY_EXPECTED` — printed distinctly, with the reason.
Emptiness is a failure only where `expected_rows > 0`. Row counts are read from
`proteomics_de/tests/expected/frozen_counts.json`, never hardcoded here, so the
expectations live in exactly one place across the whole repo.

Deliberate non-goals
--------------------
* **`--jobs` / parallel execution.** Stages 3-6 (the four viz scripts) are
  mutually independent and would parallelise cleanly, saving perhaps 20 seconds
  on a run measured in minutes. Execution is nevertheless strictly **serial**:
  a deterministic, linear, greppable log — where a traceback belongs to
  obviously one stage and stage N's output was demonstrably written before
  stage N+1 read it — is worth far more than those 20 seconds on a pipeline a
  human runs a handful of times. Not an oversight; a choice.
* **Resuming from a crash.** Use `--from <stage>` explicitly. Automatic
  "figure out what's stale and redo it" needs a real dependency-on-content
  model; `--from` is honest about being manual.
* **Rewriting the stage scripts.** The runner invokes them exactly as they are,
  from `cwd=<repo root>` with `sys.executable`, because stage 1 is cwd-locked to
  the repo root. Nothing here monkeypatches or imports pipeline internals.

Network
-------
Stages 7 (STRING), 9 (g:Profiler) and 11 (Enrichr, via gseapy) make **live
outbound API calls** with the mouse gene list. `--skip-network` skips them and
reports them as SKIPPED, not FAILED; their committed artifacts stay on disk, so
downstream stages still run against them.

Paths
-----
Everything resolves from `Path(__file__)`, never from the current working
directory. `cd /tmp && /path/to/.venv/bin/python /path/to/run_pipeline.py --list`
behaves identically to running it from the repo root.

CLI
---
    python run_pipeline.py --list                 # the stage table, run nothing
    python run_pipeline.py --all                  # every stage, in order
    python run_pipeline.py --all --skip-network   # offline stages only
    python run_pipeline.py --from viz_volcano     # that stage and everything after
    python run_pipeline.py --only viz_volcano,viz_ma_plot
    python run_pipeline.py --all --dry-run        # print the plan, execute nothing
    python run_pipeline.py --verify-frozen        # sha256 the byte-frozen outputs

Exit code is 0 only when nothing FAILed and nothing was BLOCKED by a failure.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Paths. Resolved from __file__ so cwd is irrelevant (see module docstring).
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve()
ROOT = HERE.parent
PDE = ROOT / "proteomics_de"
FROZEN_COUNTS_PATH = PDE / "tests" / "expected" / "frozen_counts.json"
STATUS_TOOL_PATH = ROOT / "tools" / "status.py"

# --------------------------------------------------------------------------
# Outcomes. Ordered worst-first for `max`-style aggregation.
# --------------------------------------------------------------------------
PASS = "PASS"
PASS_EMPTY_EXPECTED = "PASS_EMPTY_EXPECTED"
FAIL = "FAIL"
SKIPPED = "SKIPPED"
BLOCKED = "BLOCKED"
PLANNED = "PLANNED"  # --dry-run only

# Higher = worse. Used to fold per-output verdicts into a stage verdict.
_SEVERITY = {PASS: 0, PASS_EMPTY_EXPECTED: 1, FAIL: 2}


# --------------------------------------------------------------------------
# Declarative stage table
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Output:
    """One artifact a stage is contracted to produce.

    ``expected_rows`` is the whole point of this class:

    * ``None``      — existence + non-empty file is the contract (figures, JSON,
                      logs). No row count is asserted.
    * ``int``       — a literal expected data-row count (header excluded).
    * ``str``       — a KEY into ``tests/expected/frozen_counts.json``, resolved
                      at import time. Prefer this: it keeps the numbers in one
                      place repo-wide.

    ``expected_rows == 0`` means *header-only by design*. ``empty_reason`` is
    then mandatory and is printed whenever the file verifies, so nobody ever
    reads a 0 and files a bug. See module docstring.
    """

    path: str                                   # repo-root-relative
    kind: str = "table"                         # table | figure | json | text
    expected_rows: int | str | None = None
    empty_reason: str | None = None


@dataclass(frozen=True)
class Stage:
    id: str
    description: str
    script: str                                 # repo-root-relative
    outputs: tuple[Output, ...] = ()
    depends_on: tuple[str, ...] = ()
    network: bool = False                       # makes live outbound HTTP calls
    slow: bool = False                          # > ~1 min
    requires_r: bool = False                    # shells out to Rscript


# Reasons for the three by-design header-only outputs. Written out here rather
# than inline so the wording is identical everywhere it is printed.
_WHY_NO_SIG = (
    "0/1938 proteins pass FDR<0.05 at n=2 technical replicates "
    "(min adj.p = 0.305) — DECISIONS_LOG D2"
)
_WHY_NO_ORA = (
    "0 GO/KEGG/Reactome terms survive the honest detected-proteome background "
    "— DECISIONS_LOG D6; this is the report's headline finding"
)

STAGES: tuple[Stage, ...] = (
    Stage(
        id="foldchange",
        description=(
            "SILAC log2 fold change from both Protein Report sheets; chains "
            "centering_check, replicate_check and limma_test (R) in-process"
        ),
        script="proteomics_de/foldchange.py",
        requires_r=True,
        slow=True,
        outputs=(
            Output("proteomics_de/results/foldchange_all.csv",
                   expected_rows="foldchange_all_rows"),
            Output("proteomics_de/results/single_condition_proteins.csv",
                   expected_rows="single_condition_rows"),
            Output("proteomics_de/results/onoff_proteins.csv",
                   expected_rows="onoff_rows"),
            Output("proteomics_de/results/ipa_input.csv",
                   expected_rows="ipa_input_rows"),
            Output("proteomics_de/results/qc_centering.csv"),
            Output("proteomics_de/results/qc_replicate_correlation.csv"),
            Output("proteomics_de/results/replicate_correlation.png", kind="figure"),
            Output("proteomics_de/results/qc_limma.csv",
                   expected_rows="qc_limma_rows"),
            Output("proteomics_de/results/ipa_input_significant.csv",
                   expected_rows="n_significant_fdr05",
                   empty_reason=_WHY_NO_SIG),
        ),
    ),
    Stage(
        id="ipa_export",
        description=(
            "IPA deliverables with limma p-values/FDR + provenance sidecars "
            "(research1.md section 3)"
        ),
        # Its absence from this table is exactly why ipa_input_full.csv went
        # stale: --all refreshed qc_limma.csv under D9 and never rebuilt the
        # file that quotes it, so the QIAGEN deliverable carried vanilla
        # p-values while the report quoted trend/robust. Caught by
        # test_golden_outputs.py, not by the byte-freeze -- the manifest is
        # happy to freeze a self-consistent but stale file.
        script="proteomics_de/export/ipa_export.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/ipa_input_full.csv",
                   kind="table", expected_rows="ipa_input_rows"),
            Output("proteomics_de/results/ipa_input_full.txt", kind="text"),
        ),
    ),
    Stage(
        id="regulated_lists",
        description=(
            "Split ipa_input_full.csv into UP/DOWN CSVs with a linear "
            "fold_change column, sorted by descending magnitude of change "
            "-- for handing to a person rather than QIAGEN"
        ),
        script="proteomics_de/export/regulated_lists.py",
        depends_on=("ipa_export",),
        outputs=(
            Output("proteomics_de/results/regulated_up.csv",
                   kind="table", expected_rows="n_up"),
            Output("proteomics_de/results/regulated_down.csv",
                   kind="table", expected_rows="n_down"),
        ),
    ),
    Stage(
        id="regulated_lists_supplementary",
        description=(
            "Two more files closing the gap regulated_up/down.csv structurally "
            "cannot: tier-3 partial-missingness proteins reclassified from "
            "limma's imputed log2FC, and tier-1/2 fully-undetected proteins "
            "unioned with a qualitative direction call (DECISIONS_LOG D17)"
        ),
        script="proteomics_de/export/supplementary_lists.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/regulated_up_partial.csv",
                   kind="table", expected_rows="n_up_partial"),
            Output("proteomics_de/results/regulated_down_partial.csv",
                   kind="table", expected_rows="n_down_partial"),
            Output("proteomics_de/results/qualitative_changes.csv",
                   kind="table", expected_rows="n_qualitative"),
        ),
    ),
    Stage(
        id="ipa_export_extended",
        description=(
            "QIAGEN uploads that also carry the 862 proteins D17 recovered: a "
            "963-row quantitative file (715 core + 248 tier-3 partial, tagged "
            "by tier) plus two ID-only lists for the 614 qualitative proteins "
            "that have no valid fold change (DECISIONS_LOG D18), plus the "
            "2552-accession measured-proteome Reference Set (D19)"
        ),
        # Registered here for the reason D13 gives: an export that is not in
        # this table goes stale the moment an upstream stage is re-run and
        # nothing notices. ipa_input_full.csv reached QIAGEN carrying vanilla
        # p-values exactly that way. Depends on BOTH producers because it reads
        # ipa_export's ipa_input_full.csv and supplementary_lists'
        # regulated_{up,down}_partial.csv + qualitative_changes.csv. `foldchange`
        # is named explicitly as well, though it is already a transitive
        # dependency of both: D19's Reference Set reads foldchange_all.csv and
        # single_condition_proteins.csv DIRECTLY, and a dependency that is only
        # true by accident of another stage's edges is one refactor from wrong.
        script="proteomics_de/export/ipa_extended.py",
        depends_on=("foldchange", "ipa_export", "regulated_lists_supplementary"),
        outputs=(
            Output("proteomics_de/results/ipa_input_extended.csv",
                   kind="table", expected_rows="n_ipa_extended"),
            # kind="text", following ipa_input_full.txt: the runner's
            # data_row_count() reads .txt with a comma separator, so a
            # tab-delimited multi-column file only counts right by accident.
            # The row count of this twin is asserted properly by
            # ipa_extended.assert_twins_agree, which knows the delimiter.
            Output("proteomics_de/results/ipa_input_extended.txt", kind="text"),
            # These two ARE row-counted despite the .txt extension: they are
            # single-column, so there is no delimiter to get wrong.
            Output("proteomics_de/results/ipa_qualitative_up.txt",
                   kind="table", expected_rows="n_ipa_qualitative_up"),
            Output("proteomics_de/results/ipa_qualitative_down.txt",
                   kind="table", expected_rows="n_ipa_qualitative_down"),
            # The Reference Set (D19). Also single-column, so row-countable.
            # Belongs to this stage rather than a new one: it is written by the
            # same script, and its contract is a cross-file one -- it must be a
            # strict superset of the three files above, which only holds if all
            # four are regenerated together.
            Output("proteomics_de/results/ipa_background_measured.txt",
                   kind="table", expected_rows="n_ipa_background_measured"),
        ),
    ),
    Stage(
        id="qc_validate",
        description="pandera schema validation of the DE outputs; writes results/qc/qc_report.{json,md}",
        script="proteomics_de/qc/validate.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/qc/qc_report.json", kind="json"),
            Output("proteomics_de/results/qc/qc_report.md", kind="text"),
        ),
    ),
    Stage(
        id="viz_qc_plots",
        description="QC figures: sample correlation, intensity distributions, missingness, rank-abundance",
        script="proteomics_de/viz/qc_plots.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/figures/sample_correlation.png", kind="figure"),
            Output("proteomics_de/results/figures/intensity_distributions.png", kind="figure"),
            Output("proteomics_de/results/figures/missing_values.png", kind="figure"),
            Output("proteomics_de/results/figures/rank_abundance.png", kind="figure"),
            Output("proteomics_de/results/figures/figures_manifest.json", kind="json"),
        ),
    ),
    Stage(
        id="viz_volcano",
        description="Volcano plot — limma log2FC vs -log10(raw p), UP/DOWN/NO CHANGE at |log2FC|>=0.585",
        script="proteomics_de/viz/volcano.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/figures/volcano.png", kind="figure"),
            Output("proteomics_de/results/figures/volcano.svg", kind="figure"),
        ),
    ),
    Stage(
        id="viz_ma_plot",
        description="MA plot — mean log2 intensity vs limma log2 fold change",
        script="proteomics_de/viz/ma_plot.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/figures/ma_plot.png", kind="figure"),
            Output("proteomics_de/results/figures/ma_plot.svg", kind="figure"),
        ),
    ),
    Stage(
        id="viz_heatmap",
        description="Clustered heatmap of the top 40 proteins by raw limma p-value",
        script="proteomics_de/viz/heatmap.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/figures/heatmap_top_de.png", kind="figure"),
            Output("proteomics_de/results/figures/heatmap_top_de.svg", kind="figure"),
        ),
    ),
    Stage(
        id="enrich_string_ppi",
        description="STRING v12 REST (species 10090) PPI network over the 715 regulated seeds + networkx metrics",
        script="proteomics_de/enrich/string_ppi.py",
        depends_on=("foldchange",),
        network=True,
        outputs=(
            Output("proteomics_de/results/enrichment/string_node_metrics.csv",
                   expected_rows="string_nodes"),
            Output("proteomics_de/results/enrichment/string_edges.tsv"),
            Output("proteomics_de/results/enrichment/string_meta.json", kind="json"),
            Output("proteomics_de/results/enrichment/raw/string_get_ids.json", kind="json"),
            Output("proteomics_de/results/enrichment/raw/string_network.tsv"),
        ),
    ),
    Stage(
        id="enrich_network_figure",
        description="Publication PPI figure from the STRING tables (log2FC node colouring, hub labels)",
        script="proteomics_de/enrich/network_figure.py",
        depends_on=("enrich_string_ppi",),
        outputs=(
            Output("proteomics_de/results/figures/ppi_network.png", kind="figure"),
            Output("proteomics_de/results/figures/ppi_network.svg", kind="figure"),
            Output("proteomics_de/results/figures/figures_manifest_network.json", kind="json"),
        ),
    ),
    Stage(
        id="enrich_ora",
        description="g:Profiler g:GOSt over-representation (mmusculus, CUSTOM detected-proteome background) + dotplot",
        script="proteomics_de/enrich/ora.py",
        depends_on=("foldchange",),
        network=True,
        outputs=(
            # Literal 0, not a frozen_counts key: frozen_counts.json has no entry
            # for the ORA term counts. The 0 is the D6 finding, and if it ever
            # stops being 0 that is a scientific event a human must sign off on.
            Output("proteomics_de/results/enrichment/ora_up.csv",
                   expected_rows=0, empty_reason=_WHY_NO_ORA),
            Output("proteomics_de/results/enrichment/ora_down.csv",
                   expected_rows=0, empty_reason=_WHY_NO_ORA),
            Output("proteomics_de/results/enrichment/ora_meta.json", kind="json"),
            Output("proteomics_de/results/enrichment/ora_top_terms_detail.json", kind="json"),
            Output("proteomics_de/results/enrichment/raw/gprofiler_up.json", kind="json"),
            Output("proteomics_de/results/enrichment/raw/gprofiler_down.json", kind="json"),
            Output("proteomics_de/results/figures/ora_dotplot.png", kind="figure"),
        ),
    ),
    Stage(
        id="enrich_upset",
        description="UpSet plot of gene overlap across the top ORA leads (UP + DOWN)",
        script="proteomics_de/enrich/upset.py",
        depends_on=("enrich_ora",),
        outputs=(
            Output("proteomics_de/results/figures/upset.png", kind="figure"),
            Output("proteomics_de/results/figures/upset.svg", kind="figure"),
        ),
    ),
    Stage(
        id="enrich_gsea",
        description="gseapy prerank against mouse Enrichr libraries (1000 permutations, seed 42)",
        script="proteomics_de/enrich/gsea.py",
        depends_on=("foldchange",),
        network=True,
        slow=True,
        outputs=(
            Output("proteomics_de/results/enrichment/gsea_results.csv",
                   expected_rows="gsea_terms"),
            Output("proteomics_de/results/enrichment/gsea_meta.json", kind="json"),
            Output("proteomics_de/results/figures/gsea_top.png", kind="figure"),
            Output("proteomics_de/results/figures/figures_manifest_enrich.json", kind="json"),
        ),
    ),
    Stage(
        id="gated_pca_cluster",
        description="PCA + hierarchical clustering behind the n>=6 gate (fires at n=4: QC-only, writes skip_log)",
        script="proteomics_de/gated/pca_cluster.py",
        depends_on=("foldchange",),
        outputs=(
            Output("proteomics_de/results/gated/skip_log.csv"),
            Output("proteomics_de/results/gated/pca_coords.csv"),
            Output("proteomics_de/results/gated/pca_variance.csv"),
            Output("proteomics_de/results/figures/pca_qc.png", kind="figure"),
            Output("proteomics_de/results/figures/sample_dendrogram.png", kind="figure"),
            Output("proteomics_de/results/figures/figures_manifest_gated.json", kind="json"),
        ),
    ),
    Stage(
        id="report",
        description="Assemble the single self-contained interactive HTML report (must run last)",
        script="proteomics_de/report/build_report.py",
        depends_on=(
            "qc_validate",
            "viz_qc_plots", "viz_volcano", "viz_ma_plot", "viz_heatmap",
            "enrich_network_figure", "enrich_upset", "enrich_gsea",
            "gated_pca_cluster",
        ),
        outputs=(
            Output("proteomics_de/report/report.html", kind="text"),
        ),
    ),
)

STAGE_BY_ID = {s.id: s for s in STAGES}
STAGE_ORDER = {s.id: i for i, s in enumerate(STAGES)}


# --------------------------------------------------------------------------
# Stage-table integrity. Enforced at import so the table can never silently
# drift into a cycle or a forward reference; also exercised directly by the
# test suite.
# --------------------------------------------------------------------------
def validate_stage_table(stages: tuple[Stage, ...] = STAGES) -> None:
    """Raise if the declared dependency graph is not a well-ordered DAG.

    Because execution is serial and in list order, "no forward references" is
    the property that matters: every dependency must appear EARLIER in the
    list. That is strictly stronger than acyclicity and needs no toposort.
    """
    seen: set[str] = set()
    for i, stage in enumerate(stages):
        if stage.id in seen:
            raise ValueError(f"duplicate stage id: {stage.id!r}")
        for dep in stage.depends_on:
            if dep not in {s.id for s in stages}:
                raise ValueError(f"stage {stage.id!r} depends on unknown stage {dep!r}")
            if dep not in seen:
                raise ValueError(
                    f"stage {stage.id!r} (position {i + 1}) declares a FORWARD "
                    f"dependency on {dep!r}, which runs later. The stage list "
                    f"order is the execution order; dependencies must precede."
                )
        if not (ROOT / stage.script).name.endswith(".py"):
            raise ValueError(f"stage {stage.id!r}: script must be a .py file")
        for out in stage.outputs:
            if out.expected_rows == 0 and not out.empty_reason:
                raise ValueError(
                    f"stage {stage.id!r} output {out.path!r} expects 0 rows but "
                    f"declares no empty_reason. An intentionally header-only "
                    f"output must say WHY, or the next reader will 'fix' it."
                )
            if out.kind not in ("table", "figure", "json", "text"):
                raise ValueError(f"stage {stage.id!r}: unknown output kind {out.kind!r}")
        seen.add(stage.id)


validate_stage_table()


# --------------------------------------------------------------------------
# frozen_counts.json — the single source of expected row counts
# --------------------------------------------------------------------------
def load_frozen_counts(path: Path = FROZEN_COUNTS_PATH) -> dict:
    """Read tests/expected/frozen_counts.json. Keys starting with '_' are prose."""
    with path.open(encoding="utf-8") as fh:
        return {k: v for k, v in json.load(fh).items() if not k.startswith("_")}


def resolve_expected_rows(out: Output, frozen: dict) -> int | None:
    """Turn ``Output.expected_rows`` into an int (or None) using frozen_counts."""
    if out.expected_rows is None:
        return None
    if isinstance(out.expected_rows, int):
        return out.expected_rows
    if out.expected_rows not in frozen:
        raise KeyError(
            f"{out.path}: expected_rows key {out.expected_rows!r} is not in "
            f"{FROZEN_COUNTS_PATH.name}. Add it there rather than hardcoding."
        )
    return int(frozen[out.expected_rows])


# --------------------------------------------------------------------------
# Output verification
# --------------------------------------------------------------------------
def data_row_count(path: Path) -> int:
    """Rows EXCLUDING the header.

    pandas, not ``wc -l``: several of these CSVs carry quoted embedded newlines
    (some accession fields run to 32 KB), which would inflate a line count.
    """
    import pandas as pd  # local: keeps --list/--help usable without the sci stack

    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    return len(pd.read_csv(path, sep=sep, low_memory=False))


def check_output(out: Output, frozen: dict, since: float | None = None) -> tuple[str, str]:
    """Verify one artifact. Returns ``(verdict, human-readable detail)``.

    ``since`` (a ``time.time()`` stamp) enables a staleness NOTE — an output
    older than the stage that supposedly wrote it is suspicious but is NOT a
    failure here, because a stage may legitimately decide an artifact is
    already correct.
    """
    path = ROOT / out.path
    if not path.exists():
        return FAIL, "missing"

    size = path.stat().st_size
    if size == 0:
        # A truly zero-BYTE file has no header either -- that is never by design.
        return FAIL, "zero bytes on disk (not even a header)"

    stale = since is not None and path.stat().st_mtime < since

    if out.kind == "json":
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except Exception as exc:
            return FAIL, f"unparsable JSON ({type(exc).__name__})"
        return PASS, _with_stale(f"{size:,} B", stale)

    if out.kind in ("figure", "text"):
        return PASS, _with_stale(f"{size:,} B", stale)

    # kind == "table"
    try:
        n = data_row_count(path)
    except Exception as exc:
        return FAIL, f"unparsable table ({type(exc).__name__}: {exc})"

    expected = resolve_expected_rows(out, frozen)

    if expected is None:
        return PASS, _with_stale(f"{n:,} rows", stale)

    if expected == 0:
        if n == 0:
            return PASS_EMPTY_EXPECTED, _with_stale(
                f"0 rows, EXPECTED — {out.empty_reason}", stale
            )
        return FAIL, (
            f"declared header-only by design ({out.empty_reason}) but found "
            f"{n:,} rows. That is a SCIENTIFIC change, not a crash — a human "
            f"must update the stage table deliberately."
        )

    if n == expected:
        return PASS, _with_stale(f"{n:,} rows (matches expected)", stale)
    return FAIL, f"expected {expected:,} rows, found {n:,}"


def _with_stale(detail: str, stale: bool) -> str:
    return f"{detail} [NOTE: not rewritten by this run]" if stale else detail


# --------------------------------------------------------------------------
# --verify-frozen. Reuses tools/status.py rather than reimplementing sha256
# bookkeeping, so the runner and the status report can never disagree.
# --------------------------------------------------------------------------
def _load_status_tool():
    """Import tools/status.py by path (``tools/`` is not a package)."""
    spec = importlib.util.spec_from_file_location("_pipeline_status_tool", STATUS_TOOL_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {STATUS_TOOL_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def verify_frozen(allow_drift: bool = False, stream=None) -> int:
    """sha256 every frozen SCIENTIFIC OUTPUT (tests/expected/outputs.sha256).

    Delegates to ``tools.status.freeze_check()``, which owns the manifest
    parsing and hashing; this function only presents the result and decides the
    exit code. Same verdict as ``tools/status.py --check`` by construction:
    CHANGED *or* MISSING is a failure (a vanished frozen output is strictly
    worse than a modified one), downgraded to a warning by ``--allow-drift``.
    """
    out = stream or sys.stdout
    print(_rule(), file=out)
    print("FROZEN-OUTPUT VERIFICATION", file=out)
    print(f"  manifest: {STATUS_TOOL_PATH.parent.parent / 'proteomics_de/tests/expected/outputs.sha256'}",
          file=out)
    print(f"  hashing via: {STATUS_TOOL_PATH} (freeze_check)", file=out)
    print(_rule(), file=out)

    status = _load_status_tool()
    rows, counts = status.freeze_check()

    if rows is None:
        print("  manifest not found — cannot verify.", file=out)
        return 0 if allow_drift else 1

    width = max((len(rel) for rel, _ in rows), default=10)
    for rel, verdict in rows:
        mark = {"OK": "ok    ", "CHANGED": "CHANGED", "MISSING": "MISSING"}[verdict]
        print(f"  {mark}  {rel:<{width}}", file=out)

    total = len(rows)
    print(_rule("-"), file=out)
    print(
        f"  {counts['OK']} OK · {counts['CHANGED']} CHANGED · "
        f"{counts['MISSING']} MISSING   (of {total} frozen files)",
        file=out,
    )

    bad = counts["CHANGED"] + counts["MISSING"]
    if bad == 0:
        print("  No drift: every byte-frozen output is identical to its baseline.", file=out)
        return 0
    if allow_drift:
        print(f"  {bad} file(s) drifted — tolerated because --allow-drift was passed.", file=out)
        return 0
    print(
        f"  DRIFT: {bad} file(s) changed or vanished. Re-run with --allow-drift "
        f"only if you intended it, and update the manifest in the same commit.",
        file=out,
    )
    return 1


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------
def select_stages(all_: bool = False, from_: str | None = None,
                  only: str | None = None) -> list[Stage]:
    """Resolve the CLI selection to stages, always in canonical execution order."""
    if all_:
        return list(STAGES)
    if from_ is not None:
        if from_ not in STAGE_BY_ID:
            raise SystemExit(f"unknown stage id: {from_!r}\n{_id_hint()}")
        return list(STAGES[STAGE_ORDER[from_]:])
    if only is not None:
        wanted = [tok.strip() for tok in only.split(",") if tok.strip()]
        unknown = [w for w in wanted if w not in STAGE_BY_ID]
        if unknown:
            raise SystemExit(f"unknown stage id(s): {', '.join(unknown)}\n{_id_hint()}")
        # Canonical order, de-duplicated: `--only viz_volcano,foldchange` runs
        # foldchange first. The stage list order IS the correct order; letting
        # the CLI reorder it would be a footgun, not a feature.
        return [STAGE_BY_ID[i] for i in sorted(set(wanted), key=lambda i: STAGE_ORDER[i])]
    return []


def _id_hint() -> str:
    return "valid ids: " + ", ".join(s.id for s in STAGES)


# --------------------------------------------------------------------------
# Presentation
# --------------------------------------------------------------------------
def _rule(ch: str = "=") -> str:
    return ch * 78


def _flags(stage: Stage) -> str:
    bits = []
    if stage.network:
        bits.append("network")
    if stage.slow:
        bits.append("slow")
    if stage.requires_r:
        bits.append("R")
    return ",".join(bits) or "-"


def _describe_output(out: Output, frozen: dict) -> str:
    rel = out.path.replace("proteomics_de/", "")
    try:
        expected = resolve_expected_rows(out, frozen)
    except KeyError:
        expected = None
    if out.kind != "table":
        return f"{rel} [{out.kind}]"
    if expected is None:
        return f"{rel} [rows: not asserted]"
    if expected == 0:
        return f"{rel} [0 rows EXPECTED — {out.empty_reason}]"
    return f"{rel} [{expected:,} rows]"


def print_list(frozen: dict, stream=None) -> None:
    out = stream or sys.stdout
    print(_rule(), file=out)
    print(f"proteomics_de pipeline — {len(STAGES)} stages, serial execution order", file=out)
    print(f"repo root: {ROOT}", file=out)
    print(_rule(), file=out)
    for i, stage in enumerate(STAGES, start=1):
        deps = ", ".join(stage.depends_on) or "—"
        print(f"\n{i:>2}. {stage.id}   [{_flags(stage)}]", file=out)
        print(f"    {stage.description}", file=out)
        print(f"    script:  {stage.script}", file=out)
        print(f"    after:   {deps}", file=out)
        for o in stage.outputs:
            print(f"      -> {_describe_output(o, frozen)}", file=out)
    print(f"\n{_rule('-')}", file=out)
    print(
        "flags: network = live outbound API call (skippable with --skip-network) · "
        "slow = > ~1 min · R = shells out to Rscript",
        file=out,
    )
    print(
        "Three outputs are header-only BY DESIGN (marked '0 rows EXPECTED'). They "
        "verify as PASS_EMPTY_EXPECTED, never as failures.",
        file=out,
    )
    print(_rule(), file=out)


def print_plan(stages: list[Stage], skip_network: bool, frozen: dict, stream=None) -> None:
    out = stream or sys.stdout
    print(_rule(), file=out)
    print(f"DRY RUN — plan only. {len(stages)} stage(s) selected. NOTHING WILL EXECUTE.", file=out)
    print(_rule(), file=out)
    for i, stage in enumerate(stages, start=1):
        if skip_network and stage.network:
            print(f"\n{i:>2}. {stage.id}  ->  SKIPPED (--skip-network; live API stage)", file=out)
            continue
        print(f"\n{i:>2}. {stage.id}  [{_flags(stage)}]", file=out)
        print(f"    would run: {sys.executable} {stage.script}", file=out)
        print(f"    with cwd:  {ROOT}", file=out)
        print(f"    then verify {len(stage.outputs)} output(s):", file=out)
        for o in stage.outputs:
            print(f"      -> {_describe_output(o, frozen)}", file=out)
    print(f"\n{_rule('-')}", file=out)
    print("Dry run complete. No subprocess was spawned and no file was written.", file=out)
    print(_rule(), file=out)


def _fmt_seconds(sec: float) -> str:
    if sec < 60:
        return f"{sec:5.1f}s"
    return f"{int(sec // 60):d}m{sec % 60:04.1f}s"


# --------------------------------------------------------------------------
# Execution
# --------------------------------------------------------------------------
@dataclass
class StageResult:
    stage: Stage
    outcome: str
    seconds: float = 0.0
    detail: str = ""
    output_lines: list[str] = field(default_factory=list)


def run_stage(stage: Stage, frozen: dict, stream=None) -> StageResult:
    """Run one stage as a subprocess, then verify its declared outputs."""
    out = stream or sys.stdout
    script = ROOT / stage.script
    if not script.exists():
        return StageResult(stage, FAIL, 0.0, f"script not found: {stage.script}")

    print(f"\n{_rule()}", file=out)
    print(f"RUN  {stage.id}  [{_flags(stage)}]", file=out)
    print(f"     {stage.description}", file=out)
    print(f"     {sys.executable} {stage.script}   (cwd={ROOT})", file=out)
    print(_rule(), file=out)
    out.flush()

    started_wall = time.time()
    t0 = time.perf_counter()
    # No capture: the stage's own stdout/stderr streams straight through, which
    # is the whole point of serial execution (see the module docstring).
    completed = subprocess.run([sys.executable, str(script)], cwd=str(ROOT))
    elapsed = time.perf_counter() - t0

    if completed.returncode != 0:
        return StageResult(
            stage, FAIL, elapsed,
            f"exit code {completed.returncode}; outputs NOT verified",
        )

    verdicts: list[str] = []
    lines: list[str] = []
    for o in stage.outputs:
        verdict, detail = check_output(o, frozen, since=started_wall)
        verdicts.append(verdict)
        rel = o.path.replace("proteomics_de/", "")
        lines.append(f"  [{verdict:<20}] {rel}: {detail}")

    outcome = PASS
    for v in verdicts:
        if _SEVERITY[v] > _SEVERITY[outcome]:
            outcome = v

    n_empty = sum(v == PASS_EMPTY_EXPECTED for v in verdicts)
    n_fail = sum(v == FAIL for v in verdicts)
    if n_fail:
        detail = f"{n_fail}/{len(verdicts)} output(s) failed verification"
    elif n_empty:
        detail = f"{len(verdicts)} output(s) OK, {n_empty} header-only as expected"
    else:
        detail = f"{len(verdicts)} output(s) verified"

    print(f"\n  -- output verification for {stage.id} --", file=out)
    for line in lines:
        print(line, file=out)

    return StageResult(stage, outcome, elapsed, detail, lines)


def print_summary(results: list[StageResult], total_seconds: float, stream=None) -> None:
    out = stream or sys.stdout
    print(f"\n{_rule()}", file=out)
    print("SUMMARY", file=out)
    print(_rule(), file=out)
    width = max((len(r.stage.id) for r in results), default=10)
    for i, r in enumerate(results, start=1):
        t = "     —" if r.outcome in (SKIPPED, BLOCKED) else _fmt_seconds(r.seconds)
        print(f"{i:>3}. {r.stage.id:<{width}}  {r.outcome:<20}  {t}  {r.detail}", file=out)

    counts = {k: sum(r.outcome == k for r in results)
              for k in (PASS, PASS_EMPTY_EXPECTED, FAIL, SKIPPED, BLOCKED)}
    print(_rule("-"), file=out)
    print(
        f"{counts[PASS]} PASS · {counts[PASS_EMPTY_EXPECTED]} PASS_EMPTY_EXPECTED · "
        f"{counts[FAIL]} FAIL · {counts[SKIPPED]} SKIPPED · {counts[BLOCKED]} BLOCKED"
        f"      total {_fmt_seconds(total_seconds)}",
        file=out,
    )
    if counts[PASS_EMPTY_EXPECTED]:
        print(
            "PASS_EMPTY_EXPECTED is a SUCCESS: the stage produced a header-only "
            "file it was contracted to produce. See DECISIONS_LOG D2/D6.",
            file=out,
        )
    print(_rule(), file=out)


def execute(stages: list[Stage], skip_network: bool, frozen: dict, stream=None) -> int:
    out = stream or sys.stdout
    results: list[StageResult] = []
    failed: set[str] = set()
    t0 = time.perf_counter()

    for stage in stages:
        blocking = [d for d in stage.depends_on if d in failed]
        if blocking:
            # Note the asymmetry with SKIPPED below: a FAILED dependency means
            # this stage's inputs may be wrong, so we refuse to run it. A
            # deliberately SKIPPED dependency leaves the committed artifact in
            # place, so downstream stages are still meaningful.
            results.append(StageResult(
                stage, BLOCKED, 0.0,
                f"upstream failure: {', '.join(blocking)}",
            ))
            failed.add(stage.id)
            print(f"\nBLOCKED  {stage.id}  (upstream failure: {', '.join(blocking)})", file=out)
            continue

        if skip_network and stage.network:
            results.append(StageResult(
                stage, SKIPPED, 0.0,
                "--skip-network: live outbound API stage; committed artifacts left in place",
            ))
            print(f"\nSKIPPED  {stage.id}  (--skip-network)", file=out)
            continue

        result = run_stage(stage, frozen, stream=out)
        results.append(result)
        if result.outcome == FAIL:
            failed.add(stage.id)
            print(f"\nFAIL  {stage.id}: {result.detail}", file=out)

    print_summary(results, time.perf_counter() - t0, stream=out)
    return 1 if any(r.outcome in (FAIL, BLOCKED) for r in results) else 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="run_pipeline.py",
        description=(
            "Run the proteomics_de SILAC DE pipeline end to end, in the one "
            "order that works, verifying each stage's declared outputs."
        ),
        epilog=(
            "Three outputs are header-only BY DESIGN and verify as "
            "PASS_EMPTY_EXPECTED, not FAIL: results/ipa_input_significant.csv "
            "(D2) and results/enrichment/ora_{up,down}.csv (D6). "
            "Execution is serial on purpose; there is no --jobs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    mode = p.add_mutually_exclusive_group()
    mode.add_argument("--all", action="store_true", help="run every stage, in order")
    mode.add_argument("--from", dest="from_", metavar="STAGE",
                      help="run STAGE and every stage after it")
    mode.add_argument("--only", metavar="a,b,c",
                      help="run exactly these stages (comma-separated; canonical order)")
    mode.add_argument("--list", action="store_true",
                      help="print the stage table and exit")

    p.add_argument("--dry-run", action="store_true",
                   help="print the execution plan; spawn nothing, write nothing")
    p.add_argument("--skip-network", action="store_true",
                   help="skip the stages that make live outbound API calls "
                        "(STRING, g:Profiler, Enrichr); report them SKIPPED")
    p.add_argument("--verify-frozen", action="store_true",
                   help="sha256 every frozen output (tests/expected/outputs.sha256) "
                        "and report OK/CHANGED/MISSING per file")
    p.add_argument("--allow-drift", action="store_true",
                   help="with --verify-frozen: report drift but exit 0")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.allow_drift and not args.verify_frozen:
        print("--allow-drift has no effect without --verify-frozen", file=sys.stderr)

    frozen = load_frozen_counts()

    if args.list:
        print_list(frozen)
        return 0

    rc = 0
    if args.verify_frozen:
        rc = verify_frozen(allow_drift=args.allow_drift)

    selected = select_stages(all_=args.all, from_=args.from_, only=args.only)
    if not selected:
        if args.verify_frozen:
            return rc
        build_parser().print_help()
        print("\nPick one of --all / --from / --only / --list / --verify-frozen.",
              file=sys.stderr)
        return 2

    if args.dry_run:
        print_plan(selected, args.skip_network, frozen)
        return rc

    return max(rc, execute(selected, args.skip_network, frozen))


if __name__ == "__main__":
    raise SystemExit(main())
