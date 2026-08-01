#!/usr/bin/env python3
"""Generate ``proteomics_de/STATUS.md`` -- an inventory of what actually exists.

Why this script exists
----------------------
Hand-written status documents go stale the moment someone adds a file. This one
is regenerated *from the filesystem*, so "what exists" is always observed, never
remembered. Nothing here interprets the science; it reports presence, size, row
counts, hashes, and versions.

Companion documents (all three are meant to be read together):

* ``proteomics_de/DECISIONS_LOG.md`` -- human decisions (D1..D6). Hand-written.
* ``proteomics_de/BUILD_LOG.md``     -- per-work-package build history. Hand-written,
                                        append-only.
* ``proteomics_de/STATUS.md``        -- THIS script's output. Generated. Disposable.

Design constraints
------------------
* All paths resolve from ``Path(__file__)``, never from the current working
  directory, so the script behaves identically from any cwd (and from CI).
* Standard library + pandas only (pandas is already a hard pipeline dependency).
* Degrades gracefully: a missing ``Rscript``, an unreadable CSV, or a missing
  freeze manifest produces a note in the report, never a traceback.

CLI
---
    python tools/status.py            # write proteomics_de/STATUS.md
    python tools/status.py --stdout   # print the report instead of writing it
    python tools/status.py --check    # exit 1 if any byte-frozen file drifted
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# --------------------------------------------------------------------------
# Path resolution -- cwd-independent by construction.
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve()
ROOT = HERE.parent.parent                 # repo root
PDE = ROOT / "proteomics_de"
RESULTS = PDE / "results"
TESTS = PDE / "tests"
# The gate watches SCIENTIFIC OUTPUTS, not source code. The original manifest
# froze all 93 tracked files including 21 .py/.R sources, which made the gate
# fail by construction the moment anyone refactored a script -- the very
# activity the gate exists to make safe. See tools/freeze.py for the rationale
# and for why SVGs are hashed canonically rather than raw.
PROTECTED_MANIFEST = TESTS / "expected" / "outputs.sha256"
#: Retained for reference only; NOT a gate. Source files are versioned by git.
SOURCES_MANIFEST = TESTS / "expected" / "protected.sha256"
OUTPUT = PDE / "STATUS.md"

STATUS_OK = "✅ implemented"
STATUS_ELSEWHERE = "↔ implemented elsewhere"
STATUS_MISSING = "❌ missing"


# --------------------------------------------------------------------------
# (a) research1.md Section 6 master build list, items 1-20.
#
# ``doc_path`` / ``lang`` / ``desc`` are transcribed verbatim from the table at
# research1.md lines 257-279. ``actual`` lists the real file(s) on disk that
# DELIVER the item (relative to proteomics_de/) and drives the status; ``extra``
# lists supporting files that are reported but do not affect status. ``note``
# explains any divergence from the doc's proposal. Status is computed against the
# filesystem, not asserted here.
# --------------------------------------------------------------------------
SECTION6 = [
    {
        "n": 1,
        "doc_path": "config/sample_sheet.tsv + config.yaml",
        "lang": "YAML/TSV",
        "desc": "Declarative samples, groups, contrasts, thresholds, paths",
        "actual": ["config/config.yaml", "config/sample_sheet.tsv"],
        "note": "Descriptive only -- the frozen scripts do not read it (DECISIONS_LOG D1).",
    },
    {
        "n": 2,
        "doc_path": "io/load.py",
        "lang": "Python/pandas",
        "desc": "Read both SILAC sheets, normalize column names",
        "actual": ["foldchange.py"],
        "note": "Inlined: foldchange.py reads 'Protein Report L'/'H' directly via pd.read_excel.",
    },
    {
        "n": 3,
        "doc_path": "qc/schema.py (pandera)",
        "lang": "Python/pandera",
        "desc": "Validate schema, types, UniProt regex, uniqueness",
        "actual": ["qc/schema.py"],
        "extra": ["qc/validate.py", "qc/__init__.py"],
        "note": "Built as specified; validate.py is the runner that emits results/qc/qc_report.*.",
    },
    {
        "n": 4,
        "doc_path": "etl/merge.py",
        "lang": "Python/pandas",
        "desc": "Safe join w/ cardinality guard + dropped-exclusives log",
        "actual": ["foldchange.py"],
        "note": "Inlined: outer merge with indicator; exclusives -> single_condition_proteins.csv.",
    },
    {
        "n": 5,
        "doc_path": "etl/foldchange.py",
        "lang": "Python/numpy",
        "desc": "log2-space ratios, symmetric thresholds, inf/NaN masking (fixes bugs 1-3)",
        "actual": ["foldchange.py"],
        "note": "Same module, flat layout instead of etl/ (DECISIONS_LOG D1).",
    },
    {
        "n": 6,
        "doc_path": "qc/replicate_qc.py",
        "lang": "Python",
        "desc": "Replicate correlation, distribution, missingness checks",
        "actual": ["replicate_check.py", "centering_check.py"],
        "note": "Split into two flat modules: replicate correlation (bug 6) + centering check.",
    },
    {
        "n": 7,
        "doc_path": "etl/build_matrix.py",
        "lang": "Python",
        "desc": "Emit intensity_matrix.tsv + design.tsv (limma contract)",
        "actual": ["limma_test.py"],
        "note": "Inlined: limma_test.py writes the handoff CSV (_limma_input.csv) directly.",
    },
    {
        "n": 8,
        "doc_path": "de/run_limma.R",
        "lang": "R/limma+DEP",
        "desc": "MinProb impute + limma; emits limma_results.tsv",
        "actual": ["limma_test.R"],
        "note": "Same role, flat layout; emits _limma_output.csv + _limma_versions.txt.",
    },
    {
        "n": 9,
        "doc_path": "de/run_limma.py",
        "lang": "Python/subprocess",
        "desc": "Orchestrate Rscript, check exit code, read results",
        "actual": ["limma_test.py"],
        "note": "Same role, flat layout.",
    },
    {
        "n": 10,
        "doc_path": "export/ipa_export.py",
        "lang": "Python",
        "desc": "Write IPA CSV (UniProt + log2FC + FDR), export validation",
        "actual": ["foldchange.py", "limma_test.py"],
        "note": "Inlined: foldchange.py writes ipa_input.csv; limma_test.py writes ipa_input_significant.csv.",
    },
    {
        "n": 11,
        "doc_path": "viz/volcano.py / .R",
        "lang": "Python/R",
        "desc": "Volcano plot",
        "actual": ["viz/volcano.py"],
        "note": "Python variant built (the doc allowed either).",
    },
    {
        "n": 12,
        "doc_path": "viz/ma_plot.py",
        "lang": "Python",
        "desc": "MA plot",
        "actual": ["viz/ma_plot.py"],
        "note": "",
    },
    {
        "n": 13,
        "doc_path": "viz/heatmap.R (pheatmap)",
        "lang": "R",
        "desc": "Top-variable/DE protein heatmap",
        "actual": ["viz/heatmap.py"],
        "note": "Built in Python (matplotlib/scipy) instead of R/pheatmap -- no extra R dependency.",
    },
    {
        "n": 14,
        "doc_path": "viz/qc_plots.py",
        "lang": "Python/seaborn",
        "desc": "Correlation matrix, boxplots, density, missing-value heatmap",
        "actual": ["viz/qc_plots.py"],
        "note": "",
    },
    {
        "n": 15,
        "doc_path": "enrich/string_ppi.py",
        "lang": "Python/requests+igraph",
        "desc": "STRING API -> network -> igraph metrics",
        "actual": ["enrich/string_ppi.py"],
        "extra": ["enrich/enrich_common.py"],
        "note": "networkx used for graph metrics instead of igraph; STRING species=10090 (D5).",
    },
    {
        "n": 16,
        "doc_path": "enrich/ora_gsea.R",
        "lang": "R/clusterProfiler+fgsea",
        "desc": "GO/KEGG/Reactome ORA + ranked GSEA",
        "actual": ["enrich/ora.py", "enrich/gsea.py"],
        "note": "Python + web APIs (g:Profiler, gseapy/Enrichr) instead of R Bioconductor (D3).",
    },
    {
        "n": 17,
        "doc_path": "viz/network_cytoscape.py",
        "lang": "py4cytoscape",
        "desc": "Publication network with log2FC node coloring",
        "actual": ["enrich/network_figure.py"],
        "note": "Static networkx figure -- py4cytoscape needs a running Cytoscape desktop (D3).",
    },
    {
        "n": 18,
        "doc_path": "viz/upset.py",
        "lang": "Python/upsetplot",
        "desc": "Set-intersection plots",
        "actual": ["enrich/upset.py"],
        "note": "Lives under enrich/ rather than viz/ (it consumes enrichment sets).",
    },
    {
        "n": 19,
        "doc_path": "gated/pca_cluster.py",
        "lang": "Python",
        "desc": "PCA/clustering -- runs only if n>=6/group (QC-only otherwise)",
        "actual": ["gated/pca_cluster.py"],
        "note": "Gate fires: n=4 < 6, so output is QC-only + results/gated/skip_log.csv.",
    },
    {
        "n": 20,
        "doc_path": "report/build_report.py",
        "lang": "Python",
        "desc": "Assemble HTML report + provenance",
        "actual": ["report/build_report.py"],
        "extra": ["report/report_template.html", "report/report.html", "report/report_facts.json"],
        "note": "Self-contained single-file HTML report.",
    },
]

# Files that are header-only ON PURPOSE. A naive reader sees "0 rows" and assumes
# a broken pipeline; these are the honest scientific result. Keys are relative to
# proteomics_de/.
EXPECTED_EMPTY = {
    "results/ipa_input_significant.csv":
        "expected -- 0/1938 proteins pass FDR<0.05 at n=2 technical replicates "
        "(DECISIONS_LOG D2)",
    "results/enrichment/ora_up.csv":
        "expected -- 0 GO/KEGG/Reactome terms survive the honest detected-proteome "
        "background (DECISIONS_LOG D6)",
    "results/enrichment/ora_down.csv":
        "expected -- 0 GO/KEGG/Reactome terms survive the honest detected-proteome "
        "background (DECISIONS_LOG D6)",
}

# Known-good row counts, from config/config.yaml's file_contracts block and the
# verified build. A mismatch here means an artifact moved and is flagged loudly.
EXPECTED_ROWS = {
    "results/foldchange_all.csv": 1948,
    # 604 since DECISIONS_LOG D11: the 2 junk MaxQuant row-index-list
    # accessions are quarantined to results/qc/quarantine_accessions.csv.
    "results/single_condition_proteins.csv": 604,
    "results/onoff_proteins.csv": 10,
    "results/qc_limma.csv": 1938,
    "results/ipa_input.csv": 715,
    "results/enrichment/gsea_results.csv": 568,
    "results/enrichment/string_node_metrics.csv": 694,
}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def human_size(nbytes: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if nbytes < 1024 or unit == "GB":
            return f"{nbytes:,.0f} {unit}" if unit == "B" else f"{nbytes:,.1f} {unit}"
        nbytes /= 1024.0
    return f"{nbytes} B"


def data_row_count(path: Path):
    """Row count EXCLUDING the header, or None if the file is not tabular/parsable."""
    sep = "\t" if path.suffix.lower() == ".tsv" else ","
    try:
        return len(pd.read_csv(path, sep=sep, low_memory=False))
    except Exception as exc:  # unreadable/odd file -> report, don't crash
        return f"unparsable ({type(exc).__name__})"


def git_head() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        if out.returncode == 0:
            return out.stdout.strip()
        return "unknown"
    except Exception:
        return "unknown"


def rel_pde(path: Path) -> str:
    return path.relative_to(PDE).as_posix()


# --------------------------------------------------------------------------
# (a) Section 6 build-list status
# --------------------------------------------------------------------------
def section6_status():
    rows = []
    for item in SECTION6:
        present, absent = [], []
        for rel in item["actual"]:
            (present if (PDE / rel).exists() else absent).append(rel)

        if not present:
            status = STATUS_MISSING
        else:
            # "Implemented as proposed" means every actual path matches the doc's
            # proposed path (compared on basename+dir, ignoring the doc's " + " lists).
            doc_targets = {p.strip() for p in item["doc_path"].replace(" / ", " + ").split(" + ")}
            doc_targets = {t.split(" ")[0] for t in doc_targets}  # drop "(pandera)" etc.
            # config item 1 lists "config.yaml" bare; normalise to its real home.
            doc_targets = {"config/config.yaml" if t == "config.yaml" else t for t in doc_targets}
            status = STATUS_OK if set(present) <= doc_targets else STATUS_ELSEWHERE
        extra_present = [r for r in item.get("extra", []) if (PDE / r).exists()]
        rows.append({**item, "status": status, "present": present,
                     "absent": absent, "extra_present": extra_present})
    return rows


def render_section6(rows) -> list[str]:
    n_ok = sum(r["status"] == STATUS_OK for r in rows)
    n_el = sum(r["status"] == STATUS_ELSEWHERE for r in rows)
    n_missing = sum(r["status"] == STATUS_MISSING for r in rows)

    out = [
        "## a) research1.md Section 6 build list (items 1-20)",
        "",
        f"**{n_ok} implemented as proposed · {n_el} implemented elsewhere · "
        f"{n_missing} missing** (of 20).",
        "",
        "`↔ implemented elsewhere` is not a gap: the build kept the flat, "
        "individually-verified script layout instead of the doc's `io/ etl/ de/` "
        "tree (DECISIONS_LOG D1), and substituted Python for R in three places.",
        "",
        "| # | Doc proposes | Actual file(s) on disk | Status | Note |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        actual = ", ".join(f"`{p}`" for p in r["present"]) or "—"
        if r["absent"]:
            actual += " (missing: " + ", ".join(f"`{p}`" for p in r["absent"]) + ")"
        if r["extra_present"]:
            actual += "<br>+ " + ", ".join(f"`{p}`" for p in r["extra_present"])
        out.append(
            f"| {r['n']} | `{r['doc_path']}` | {actual} | {r['status']} | {r['note']} |"
        )
    out.append("")
    return out


# --------------------------------------------------------------------------
# (b) Artifact inventory
# --------------------------------------------------------------------------
def render_artifacts() -> list[str]:
    out = ["## b) Artifact inventory (`proteomics_de/results/`)", ""]
    if not RESULTS.is_dir():
        out += ["_`results/` directory not found._", ""]
        return out

    files = sorted(p for p in RESULTS.rglob("*") if p.is_file())
    total_bytes = sum(p.stat().st_size for p in files)
    tabular = [p for p in files if p.suffix.lower() in (".csv", ".tsv")]

    out += [
        f"**{len(files)} files**, {human_size(total_bytes)} total, "
        f"{len(tabular)} tabular (`.csv`/`.tsv`).",
        "",
        "Row counts EXCLUDE the header. Three files are header-only **by design** — "
        "they are the honest scientific result, not a failure. See "
        "`DECISIONS_LOG.md` D2 and D6.",
        "",
        "| File | Size | Data rows | Note |",
        "|---|---|---|---|",
    ]

    mismatches = []
    for p in files:
        rel = rel_pde(p)
        size = human_size(p.stat().st_size)
        note = ""
        if p.suffix.lower() in (".csv", ".tsv"):
            n = data_row_count(p)
            if isinstance(n, int):
                rel_key = rel
                if rel_key in EXPECTED_EMPTY and n == 0:
                    rows_cell = f"**0 rows ({EXPECTED_EMPTY[rel_key]})**"
                else:
                    rows_cell = f"{n:,}"
                    if n == 0:
                        rows_cell += " ⚠️ unexpectedly empty"
                exp = EXPECTED_ROWS.get(rel_key)
                if exp is not None:
                    if n == exp:
                        note = f"matches expected {exp:,}"
                    else:
                        note = f"⚠️ **MISMATCH — expected {exp:,}, got {n:,}**"
                        mismatches.append((rel_key, exp, n))
            else:
                rows_cell = str(n)
        else:
            rows_cell = "—"
        out.append(f"| `{rel}` | {size} | {rows_cell} | {note} |")

    out.append("")
    if mismatches:
        out += ["> ⚠️ **Row-count mismatches detected:**", ""]
        out += [f"> * `{k}`: expected {e:,}, found {g}" for k, e, g in mismatches]
        out.append("")
    else:
        out += [
            "All 7 headline row counts match the contract in `config/config.yaml`.",
            "",
        ]
    return out


# --------------------------------------------------------------------------
# (c) Byte-freeze drift
# --------------------------------------------------------------------------
def freeze_check():
    """Return (rows, counts) where rows are (relpath, status)."""
    if not PROTECTED_MANIFEST.exists():
        return None, None
    # Delegate hashing to tools/freeze.py so the manifest's per-file mode
    # (raw vs. svg-canon) is honoured. A raw byte compare on SVGs would report
    # spurious drift: matplotlib stamps a wall-clock <dc:date> and random
    # element-id salts into every SVG, so two runs of identical code differ.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import freeze as _freeze

    ok, changed, missing, _extra = _freeze.verify(PROTECTED_MANIFEST, ROOT)
    rows = (
        [(rel, "OK") for rel in ok]
        + [(rel, "CHANGED") for rel in changed]
        + [(rel, "MISSING") for rel in missing]
    )
    rows.sort(key=lambda r: r[0])
    counts = {
        "OK": sum(s == "OK" for _, s in rows),
        "CHANGED": sum(s == "CHANGED" for _, s in rows),
        "MISSING": sum(s == "MISSING" for _, s in rows),
    }
    return rows, counts


def render_freeze(rows, counts) -> list[str]:
    out = ["## c) Byte-freeze drift", ""]
    if rows is None:
        out += [
            f"_Freeze manifest not found at `{PROTECTED_MANIFEST.relative_to(ROOT)}` — "
            "drift cannot be checked._",
            "",
        ]
        return out
    total = len(rows)
    verdict = (
        "✅ **No drift.** Every frozen file is byte-identical to its baseline."
        if counts["CHANGED"] == 0 and counts["MISSING"] == 0
        else "🚨 **DRIFT DETECTED.**"
    )
    out += [
        f"Manifest: `{PROTECTED_MANIFEST.relative_to(ROOT)}` ({total} files).",
        "",
        f"**{counts['OK']} OK · {counts['CHANGED']} CHANGED · {counts['MISSING']} MISSING**",
        "",
        verdict,
        "",
    ]
    bad = [(rel, st) for rel, st in rows if st != "OK"]
    if bad:
        out += ["| File | Status |", "|---|---|"]
        out += [f"| `{rel}` | **{st}** |" for rel, st in bad]
        out.append("")
    return out


# --------------------------------------------------------------------------
# (d) Test coverage
# --------------------------------------------------------------------------
def pipeline_modules() -> list[Path]:
    skip_dirs = {"tests", "__pycache__", "results", "config", "assets"}
    mods = []
    for p in sorted(PDE.rglob("*.py")):
        if any(part in skip_dirs for part in p.relative_to(PDE).parts[:-1]):
            continue
        if p.name == "__init__.py":
            continue
        mods.append(p)
    return mods


def test_functions(path: Path):
    """Names of test functions/methods in a file, via ast (never executed)."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    except Exception:
        return None
    names = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test"):
            names.append(node.name)
    return names


def render_tests() -> list[str]:
    out = ["## d) Test coverage", ""]
    if not TESTS.is_dir():
        out += ["_No `proteomics_de/tests/` directory._", ""]
        return out

    test_files = sorted(p for p in TESTS.rglob("test_*.py"))
    if not test_files:
        out += [
            "_No `test_*.py` files found under `proteomics_de/tests/` "
            "(the directory exists but holds no tests yet)._",
            "",
        ]
    else:
        total = 0
        out += ["| Test file | Test functions |", "|---|---|"]
        for p in test_files:
            fns = test_functions(p)
            if fns is None:
                out.append(f"| `{rel_pde(p)}` | (parse error) |")
                continue
            total += len(fns)
            out.append(f"| `{rel_pde(p)}` | {len(fns)} |")
        out += ["", f"**{len(test_files)} test files · {total} test functions.**", ""]

    # module -> has a test file naming it?
    stems = {p.stem for p in test_files}
    blob = " ".join(sorted(stems))
    out += [
        "Module coverage — a module counts as covered if a test file is named "
        "after it (`test_<module>.py`) or mentions it in its filename:",
        "",
        "| Pipeline module | Test file? |",
        "|---|---|",
    ]
    covered = 0
    for m in pipeline_modules():
        hit = f"test_{m.stem}" in stems or m.stem in blob
        covered += bool(hit)
        out.append(f"| `{rel_pde(m)}` | {'✅ yes' if hit else '❌ none'} |")
    mods = pipeline_modules()
    out += [
        "",
        f"**{covered}/{len(mods)} pipeline modules have a matching test file.**",
        "",
    ]
    return out


# --------------------------------------------------------------------------
# (e) Environment
# --------------------------------------------------------------------------
def r_info() -> list[str]:
    """Rscript version + limma/imputeLCMD presence. Never raises."""
    lines = []
    try:
        ver = subprocess.run(["Rscript", "--version"], capture_output=True, text=True, timeout=30)
        text = (ver.stdout + ver.stderr).strip().splitlines()
        lines.append(f"* Rscript: `{text[0] if text else 'present (version unknown)'}`")
    except FileNotFoundError:
        return ["* Rscript: ❌ not on PATH — the limma step (item 8) cannot run here."]
    except Exception as exc:
        return [f"* Rscript: ⚠️ probe failed ({type(exc).__name__})"]

    probe = (
        'cat(paste(sapply(c("limma","imputeLCMD"), function(p) '
        'paste0(p, "=", if (requireNamespace(p, quietly=TRUE)) '
        'as.character(utils::packageVersion(p)) else "MISSING")), collapse=" "))'
    )
    try:
        out = subprocess.run(["Rscript", "-e", probe], capture_output=True, text=True, timeout=90)
        payload = out.stdout.strip()
        if payload:
            for tok in payload.split():
                name, _, val = tok.partition("=")
                mark = "❌" if val == "MISSING" else "✅"
                lines.append(f"* R package `{name}`: {mark} {val}")
        else:
            lines.append("* R packages: ⚠️ probe returned nothing")
    except Exception as exc:
        lines.append(f"* R packages: ⚠️ probe failed ({type(exc).__name__})")
    return lines


def render_env() -> list[str]:
    venv = ROOT / ".venv"
    out = [
        "## e) Environment",
        "",
        f"* Python: `{sys.version.split()[0]}` (`{sys.executable}`)",
        f"* pandas: `{pd.__version__}`",
        f"* `.venv/`: {'✅ present' if venv.is_dir() else '❌ absent'}"
        + (f" (`{venv}`)" if venv.is_dir() else ""),
    ]
    for name in ("pyproject.toml", "requirements-dev.txt", "requirements-lock.txt"):
        p = ROOT / name
        out.append(f"* `{name}`: {'✅ present' if p.exists() else '❌ absent'}")
    out += r_info()
    out.append("")
    return out


# --------------------------------------------------------------------------
# Report assembly
# --------------------------------------------------------------------------
def build_report() -> tuple[str, dict]:
    s6 = section6_status()
    freeze_rows, freeze_counts = freeze_check()

    now = datetime.now(timezone.utc).astimezone()
    lines = [
        "# STATUS — what actually exists in `proteomics_de/`",
        "",
        "> **This file is GENERATED. Do not hand-edit.**",
        "> Regenerate with `.venv/bin/python tools/status.py` "
        "(source: `tools/status.py`). Every number below is read off the "
        "filesystem at generation time, so it cannot silently go stale.",
        "",
        f"* Generated: **{now.strftime('%Y-%m-%d %H:%M:%S %Z')}**",
        f"* Git HEAD: **`{git_head()}`**",
        f"* Repo root: `{ROOT}`",
        "",
        "Companions: `DECISIONS_LOG.md` (human decisions D1–D6) · "
        "`BUILD_LOG.md` (per-work-package history) · "
        "`../research1.md` §Build Log (per-bug narrative).",
        "",
        "---",
        "",
    ]
    lines += render_section6(s6)
    lines += ["---", ""]
    lines += render_artifacts()
    lines += ["---", ""]
    lines += render_freeze(freeze_rows, freeze_counts)
    lines += ["---", ""]
    lines += render_tests()
    lines += ["---", ""]
    lines += render_env()

    summary = {
        "ok": sum(r["status"] == STATUS_OK for r in s6),
        "elsewhere": sum(r["status"] == STATUS_ELSEWHERE for r in s6),
        "missing": sum(r["status"] == STATUS_MISSING for r in s6),
        "freeze": freeze_counts,
    }
    return "\n".join(lines) + "\n", summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stdout", action="store_true",
                    help="print the report instead of writing STATUS.md")
    ap.add_argument("--check", action="store_true",
                    help="exit nonzero if any byte-frozen file CHANGED or is MISSING")
    args = ap.parse_args(argv)

    report, summary = build_report()

    if args.stdout:
        sys.stdout.write(report)
    else:
        OUTPUT.write_text(report, encoding="utf-8")
        print(f"wrote {OUTPUT}")

    print(
        f"Section 6: {summary['ok']} implemented, {summary['elsewhere']} elsewhere, "
        f"{summary['missing']} missing",
        file=sys.stderr,
    )

    fc = summary["freeze"]
    if fc is None:
        print("freeze: manifest not found", file=sys.stderr)
        return 1 if args.check else 0
    print(
        f"freeze: {fc['OK']} OK, {fc['CHANGED']} CHANGED, {fc['MISSING']} MISSING",
        file=sys.stderr,
    )
    if args.check and (fc["CHANGED"] or fc["MISSING"]):
        print("FAIL: byte-frozen outputs have drifted", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
