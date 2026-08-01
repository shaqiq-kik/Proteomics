"""Bug 7 — per-protein statistical significance via limma (R) + MinProb imputation.

The earlier bugs gave every both-condition protein a fold-change and a +/-0.585
"regulated" call, but a fold-change alone says nothing about whether the change is
*reproducible* — a protein measured twice can look "UP" by chance. Bug 7 adds the
missing half: a per-protein p-value and a Benjamini-Hochberg adjusted p-value
(FDR), computed with **limma** (a moderated t-test that borrows variance across
proteins) after **MinProb imputation** of missing intensities, following Peng et
al. 2024 (Nat Commun, DOI 10.1038/s41467-024-47899-w).

The statistics live in R (``limma_test.R``); this module owns the subprocess
boundary and the CSV file handoff. The split is deliberate:

  * Fork 1 (handoff)  — Python <-> R talk over CSV files, not stdin pipes, so the
                        exact matrix R saw is auditable on disk.
  * Fork 2 (eligibility) — only the ~1938 genuine both-condition proteins are
                        tested; the 10 ON_OFF proteins are excluded (testing them
                        would invent an entire absent condition) and the 606
                        single-condition proteins are not in foldchange_all.csv.
  * Fork 3 (transform) — log2 + imputation happen in R, close to limma, per Peng.

It only ADDS files (qc_limma.csv, qc_limma_vanilla.csv, de/limma_results.tsv,
ipa_input_significant.csv) and the auditable intermediates; it overwrites NONE of
the existing pipeline outputs. It reuses the existing ``regulated`` column for the
+/-0.585 part (it does not recompute the threshold) and imports nothing from
``foldchange`` (the module graph stays acyclic; it reads ``foldchange_all.csv``
from disk).

Two user decisions shape the current output (see ``DECISIONS_LOG.md``):

  * **D9** — ``eBayes(trend=TRUE, robust=TRUE)`` is the default flavour, per
    research1.md line 124 and proteomics practice. Vanilla is still run on every
    invocation and preserved as ``results/qc_limma_vanilla.csv``, and the two are
    cross-checked: the switch may move p-values, never logFC.
  * **D10** — the limma columns ``n_imputed, AveExpr, t, B`` are restored
    (research1.md line 169), and ``results/de/limma_results.tsv`` carries the
    documented contract schema including ``fold_change`` and ``direction``.
    ``n_imputed`` is the important one: MinProb is stochastic and this study is
    n=2, so without it nothing distinguishes a measured value from an invented
    one.

If the R worker fails for any reason, this module raises loudly and writes NO
result files — it never proceeds on empty or partial p-values.
"""

import os
import subprocess
import sys

import numpy as np
import pandas as pd

# --- Constants (single source of truth) -------------------------------------
ADJ_P_THRESHOLD = 0.05          # FDR cutoff for "significant"
R_SEED          = 42            # passed to R so MinProb imputation is reproducible
LIMMA_INPUT     = "_limma_input.csv"
LIMMA_OUTPUT    = "_limma_output.csv"
LIMMA_VERSIONS  = "_limma_versions.txt"
LIMMA_DESIGN    = "_limma_design.tsv"   # design handed to R via --design

#: eBayes flavour used by default. DECISIONS_LOG **D9** flipped this from
#: "vanilla" to "trend_robust": research1.md line 124 specifies
#: ``eBayes(trend=TRUE, robust=TRUE)`` and it is the proteomics field standard
#: (intensity-dependent prior variance + outlier-robust moderation). The build
#: originally shipped vanilla purely as a byte-reproducible baseline.
#:
#: The switch is confined to the variance model: logFC is BIT-IDENTICAL between
#: the two flavours (eBayes moderates variances, it does not refit
#: coefficients), only the moderated t and the p-values move. That invariant is
#: not assumed -- :func:`run_limma_test` re-runs vanilla on every default
#: invocation and *asserts* it (see ``vanilla_companion`` below).
DEFAULT_EBAYES_MODE = "trend_robust"

#: The vanilla companion, preserved so both flavours stay comparable (D9).
VANILLA_EBAYES_MODE = "vanilla"
VANILLA_QC_FILENAME = "qc_limma_vanilla.csv"
VANILLA_OUTPUT      = "_limma_output_vanilla.csv"

#: research1.md Section 1 file contract, alongside intensity_matrix.tsv/design.tsv.
DE_SUBDIR           = "de"
LIMMA_RESULTS_NAME  = "limma_results.tsv"

#: The four columns DECISIONS_LOG **D10** restores, appended (never interleaved)
#: to the R worker's original ``id, limma_log2FC, p_value, adj_p_value`` contract.
D10_COLUMNS = ["n_imputed", "AveExpr", "t", "B"]

#: ``results/de/limma_results.tsv`` — research1.md line 169, verbatim.
LIMMA_RESULTS_COLUMNS = [
    "accession", "gene", "logFC", "fold_change", "AveExpr", "t",
    "P.Value", "adj.P.Val", "B", "n_imputed", "direction",
]

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_R_SCRIPT = "limma_test.R"  # resolved relative to cwd=_HERE in the subprocess call

if _ROOT not in sys.path:  # importable however the caller was launched
    sys.path.insert(0, _ROOT)

from proteomics_de.config import constants       # noqa: E402
from proteomics_de.config import design          # noqa: E402
from proteomics_de.etl import build_matrix       # noqa: E402

#: |log2FC| boundary for the ``direction`` call in ``results/de/limma_results.tsv``
#: (research1.md line 140-141). Taken from ``config.constants`` rather than
#: re-typed, so this module does not become yet another copy of 0.585.
LOG2_THRESHOLD = constants.LOG2_THRESHOLD

# The FDR cutoff above predates config/, and test_limma_contract imports it from
# here. Rather than move it (and break that import) or let the two silently
# drift, assert they agree at import time.
assert ADJ_P_THRESHOLD == constants.ADJ_P_THRESHOLD, (
    f"ADJ_P_THRESHOLD disagrees with config.constants: "
    f"{ADJ_P_THRESHOLD} != {constants.ADJ_P_THRESHOLD}"
)

DEFAULT_RESULTS_DIR = os.path.join(_HERE, "results")
DEFAULT_FOLDCHANGE_CSV = os.path.join(DEFAULT_RESULTS_DIR, "foldchange_all.csv")

# The design, read from config/sample_sheet.tsv rather than hardcoded here.
#
# COLUMN ORDER IS DELIBERATELY *PHYSICAL* (acquisition order), NOT grouped.
# ------------------------------------------------------------------------
# It is tempting to hand R the columns in canonical control-then-treated order
# and let position imply the design. Doing so makes the matrix layout depend on
# the condition labels, and MinProb imputation is stochastic *per column*: with
# set.seed(42), reordering the columns hands different random draws to different
# samples. Measured during the D7 relabelling: the 1578 fully-observed proteins
# inverted exactly (max |new + old| = 0.0) but the 360 proteins with a missing
# value did not, because they had been imputed from a different draw -- p-values
# moved by up to 0.80 for a change that is supposed to be a pure relabelling.
#
# Keeping the matrix in acquisition order and letting design.tsv carry the group
# labels makes the imputation invariant, so relabelling the conditions negates
# logFC and leaves every p-value untouched -- which is the mathematically correct
# behaviour for swapping the levels of a two-group contrast.
_CTRL_COLS = design.control_columns()
_TRT_COLS = design.treated_columns()

#: Acquisition order: the order the channels appear in the workbook.
_PHYSICAL_COLS = ["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581"]
_INTENSITY_COLS = _PHYSICAL_COLS
_COL_TO_GROUP = {c: "control" for c in _CTRL_COLS}
_COL_TO_GROUP.update({c: "treated" for c in _TRT_COLS})
_GROUP_VECTOR = [_COL_TO_GROUP[c] for c in _INTENSITY_COLS]

# Handoff names follow the physical order and encode each column's group, so the
# CSV header stays self-describing after a relabelling.
_HANDOFF_NAMES = [
    ("ctrl_" if _COL_TO_GROUP[c] == "control" else "trt_") + c.split()[-1]
    for c in _INTENSITY_COLS
]
_HANDOFF_COLS = ["id", "gene"] + _HANDOFF_NAMES

# The sheet must describe exactly this workbook's four channels.
if sorted(_CTRL_COLS + _TRT_COLS) != sorted(_PHYSICAL_COLS):
    raise ValueError(
        "BUG7: sample_sheet.tsv does not describe this workbook.\n"
        f"  sheet resolves to : {_CTRL_COLS + _TRT_COLS}\n"
        f"  workbook provides : {_PHYSICAL_COLS}\n"
        "See proteomics_de/config/design.py."
    )
# R sets the reference level from the first group it sees, so control must be
# present and must be the first label in the design file (see design.py's
# canonical sort and tests/test_design.py).
if "control" not in _GROUP_VECTOR or "treated" not in _GROUP_VECTOR:
    raise ValueError(f"BUG7: design must contain both groups, got {_GROUP_VECTOR}")


def _missing_to_blank(series):
    """Raw intensity -> float, with blank/non-numeric/<=0 turned into NaN.

    NaN is written by ``to_csv`` as an empty cell (its default ``na_rep``), which
    R's ``read.csv`` reads back as ``NA`` so MinProb fills it. This is the single
    missing-value rule for the handoff -- shared with ``intensity_matrix.tsv``
    via ``etl.build_matrix`` so the two cannot drift apart.
    """
    return build_matrix.intensity_series(series)


def _write_design_handoff(path):
    """Write the ``--design`` TSV the R worker reads (columns: sample, group).

    Distinct from the ``results/de/design.tsv`` contract file: R needs `sample`
    to *name the intensity columns of the input CSV*, so this one carries the
    handoff names (``ctrl_31578``), not the bare sample ids (``31578``). See the
    "Known contract wart" note in ``etl/build_matrix``.
    """
    pd.DataFrame({"sample": _HANDOFF_NAMES, "group": _GROUP_VECTOR}).to_csv(
        path, sep="\t", index=False, encoding="utf-8", lineterminator="\n"
    )
    return path


def run_limma_test(foldchange_csv=DEFAULT_FOLDCHANGE_CSV, outdir=DEFAULT_RESULTS_DIR,
                   ebayes_mode=DEFAULT_EBAYES_MODE, qc_filename="qc_limma.csv",
                   output_name=LIMMA_OUTPUT, reuse_input=False, write_ipa=True,
                   write_contract=True, vanilla_companion=True):
    """Add limma p-values / FDR to the both-condition proteins.

    Reads ``foldchange_csv`` (the both-condition group), tests the eligible rows
    (``onoff`` empty) in R via MinProb + limma, and writes NEW files into
    ``outdir``: ``qc_limma.csv`` (full audit), ``de/limma_results.tsv``
    (research1.md's contract schema) and ``ipa_input_significant.csv`` (the
    honest "regulated AND significant" IPA list). Returns a summary-counts dict.
    Raises ``RuntimeError`` if the R worker fails (no partial outputs).

      * ``ebayes_mode``  — "trend_robust" (default, DECISIONS_LOG D9) or
                           "vanilla"; forwarded to R as ``--mode``.
      * ``qc_filename``  — name of the qc CSV under ``outdir``.
      * ``output_name``  — name of the R worker's output intermediate in
                           ``proteomics_de/``.
      * ``reuse_input``  — if True, reuse the existing ``_limma_input.csv``
                           instead of regenerating it (proves identical data
                           going in — this is how the vanilla companion below
                           guarantees it is comparing like with like).
      * ``write_ipa``    — if False, skip the IPA file entirely.
      * ``write_contract`` — if False, skip ``results/de/`` (the companion run
                           must not overwrite the primary run's contract files).
      * ``vanilla_companion`` — if True and the mode is "trend_robust", ALSO run
                           vanilla over the identical input and write
                           ``results/qc_limma_vanilla.csv``.

    **Why the companion runs every time.** D9 made trend/robust the default on
    the strength of one claim: the change is confined to the variance model, so
    logFC cannot move. Preserving vanilla as a file someone regenerated once
    would let that claim rot. Running both and asserting
    ``logFC_trend == logFC_vanilla`` on every invocation turns it into a
    standing check, and costs one extra R call over an input that is reused
    byte-for-byte rather than rebuilt.
    """
    if ebayes_mode not in (VANILLA_EBAYES_MODE, "trend_robust"):
        raise ValueError(
            f"ebayes_mode must be 'vanilla' or 'trend_robust', got: {ebayes_mode!r}"
        )

    # 1) Eligibility (Fork 2). Read everything as strings + keep empty cells empty
    #    so the raw intensities pass through untouched and no NaN literal sneaks in.
    fc = pd.read_csv(foldchange_csv, dtype=str, keep_default_na=False)
    onoff_nonempty = fc["onoff"].str.strip() != ""
    n_onoff = int(onoff_nonempty.sum())
    eligible = fc[~onoff_nonempty].reset_index(drop=True)
    eligible_count = len(eligible)

    assert eligible_count == len(fc) - n_onoff, (
        f"eligible_count {eligible_count} != both-group {len(fc)} - ON_OFF {n_onoff}"
    )
    print("[Bug 7] limma statistical testing (R + MinProb imputation)")
    print(f"  eligible_count : {eligible_count}  "
          f"(both-group {len(fc)} - ON_OFF {n_onoff})")

    input_path = os.path.join(_HERE, LIMMA_INPUT)
    design_path = os.path.join(_HERE, LIMMA_DESIGN)
    output_path = os.path.join(_HERE, output_name)
    # Versions filename is derived from the output name (mirrors the R worker), so a
    # non-default output writes its own versions file rather than the committed one.
    versions_name = os.path.splitext(output_name)[0].replace(
        "_limma_output", "_limma_versions") + ".txt"
    versions_path = os.path.join(_HERE, versions_name)

    # 2) Write the handoff CSV (Fork 1) in the EXACT column order — unless reusing
    #    the existing intermediate (so the data going in is provably the same).
    if reuse_input:
        if not os.path.exists(input_path):
            raise RuntimeError(
                f"BUG7: reuse_input=True but {LIMMA_INPUT} is missing — run the "
                "baseline first to produce it."
            )
        print(f"  reusing existing {LIMMA_INPUT} (not regenerated)")
    else:
        inp = pd.DataFrame(
            {
                "id": eligible["UniProt Accession Number"].values,
                "gene": eligible["Gene names"].values,
            }
        )
        for handoff_name, raw_col in zip(_HANDOFF_NAMES, _INTENSITY_COLS):
            inp[handoff_name] = _missing_to_blank(eligible[raw_col])
        inp[_HANDOFF_COLS].to_csv(input_path, index=False, encoding="utf-8")

    # 2b) The design the R worker reads (--design). Derived from the sample sheet,
    #     not from the data, so it is written even when the input is reused.
    _write_design_handoff(design_path)

    # 2c) research1.md's Section 1 file contract: results/de/intensity_matrix.tsv
    #     + results/de/design.tsv. Additive and replicate-count-agnostic; nothing
    #     downstream reads them yet, so they cannot perturb the frozen numbers.
    #     Skipped for the companion run, which must not restate the primary
    #     run's contract.
    if write_contract:
        build_matrix.build(eligible, outdir)

    # Remove any stale worker artifacts so a leftover from a previous run can never
    # be mistaken for a fresh success (critical for the fail-loud guarantee).
    for stale in (output_path, versions_path):
        if os.path.exists(stale):
            os.remove(stale)

    # 3) Call the R worker. cwd=_HERE so the script + relative file names resolve
    #    and the intermediates land in proteomics_de/.
    try:
        result = subprocess.run(
            ["Rscript", _R_SCRIPT,
             "--in", LIMMA_INPUT, "--out", output_name,
             "--seed", str(R_SEED), "--mode", ebayes_mode,
             "--design", LIMMA_DESIGN],
            capture_output=True, text=True, cwd=_HERE,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "BUG7: 'Rscript' not found on PATH — cannot run the limma worker. "
            "Install R (and the limma + imputeLCMD packages). No outputs written."
        ) from exc

    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0 or not os.path.exists(output_path):
        raise RuntimeError(
            "BUG7: limma R worker failed — refusing to continue with no/partial "
            f"p-values. returncode={result.returncode}\n"
            f"--- R stderr ---\n{result.stderr}\n--- R stdout ---\n{result.stdout}"
        )

    # 4) Read results and validate before trusting them.
    res = pd.read_csv(output_path)
    for col in ("p_value", "adj_p_value", "limma_log2FC", "AveExpr", "t", "B"):
        vals = pd.to_numeric(res[col], errors="coerce")
        assert vals.notna().all() and np.isfinite(vals).all(), (
            f"BUG7: non-finite values found in '{col}' of {output_name}"
        )
    assert len(res) == eligible_count, (
        f"BUG7: R returned {len(res)} rows, expected eligible_count {eligible_count}"
    )

    # n_imputed (D10) is the one column whose truth lives on THIS side of the
    # boundary too: R counted NAs on its log2 matrix, and the same matrix was
    # written here by _missing_to_blank. Recompute and compare, so a silent
    # divergence in the missing-value rule surfaces as a hard failure rather
    # than as a column nobody can trust.
    n_imputed = pd.to_numeric(res["n_imputed"], errors="coerce")
    assert n_imputed.notna().all(), "BUG7: non-numeric n_imputed from R"
    assert n_imputed.between(0, len(_INTENSITY_COLS)).all(), (
        f"BUG7: n_imputed outside [0, {len(_INTENSITY_COLS)}]"
    )
    expected_missing = pd.concat(
        [_missing_to_blank(eligible[c]).isna() for c in _INTENSITY_COLS], axis=1
    ).sum(axis=1)
    if not reuse_input:
        # (only meaningful when this call wrote the handoff itself; with
        #  reuse_input the CSV on disk is the authority, not `eligible`)
        assert res["id"].tolist() == eligible["UniProt Accession Number"].tolist(), (
            "BUG7: R returned the proteins in a different order than they were "
            "sent — the per-row joins below would be silently misaligned."
        )
        assert n_imputed.tolist() == expected_missing.tolist(), (
            "BUG7: R's n_imputed disagrees with the missing-value count Python "
            "wrote into the handoff — the two sides no longer share one "
            "missing-value rule."
        )

    # 5) Merge back onto the eligible rows by id (ids are unique here, so 1:1).
    merged = eligible.merge(
        res[["id", "limma_log2FC", "p_value", "adj_p_value"] + D10_COLUMNS],
        left_on="UniProt Accession Number", right_on="id", how="left",
    )
    assert merged["adj_p_value"].notna().all(), (
        "BUG7: an eligible protein has no limma result after the merge."
    )
    merged["significant"] = merged["adj_p_value"] < ADJ_P_THRESHOLD
    is_regulated = merged["regulated"].isin(["UP", "DOWN"])
    merged["regulated_significant"] = is_regulated & merged["significant"]

    os.makedirs(outdir, exist_ok=True)

    # 6a) qc_limma.csv — full audit of every eligible protein.
    #      The raw intensity columns are echoed in canonical sample-sheet order
    #      (dicts preserve insertion order, so the CSV layout is unchanged).
    qc_cols = {
        "id": merged["UniProt Accession Number"],
        "gene": merged["Gene names"],
    }
    for raw_col in _INTENSITY_COLS:
        qc_cols[raw_col] = merged[raw_col]
    qc_cols.update(
        {
            "limma_log2FC": merged["limma_log2FC"].round(6),
            "p_value": merged["p_value"],
            "adj_p_value": merged["adj_p_value"],
            "significant": merged["significant"],
            "regulated": merged["regulated"],
        }
    )
    #      D10: the restored limma columns are APPENDED after the original 11.
    #      Never inserted, never reordered — viz/*, gated/ and enrich/gsea.py
    #      select by name, but the 11-column prefix is also what every committed
    #      review of this file describes.
    qc_cols["n_imputed"] = merged["n_imputed"].astype(int)
    for col in ("AveExpr", "t", "B"):
        qc_cols[col] = merged[col]
    qc = pd.DataFrame(qc_cols)
    qc_path = os.path.join(outdir, qc_filename)
    qc.to_csv(qc_path, index=False, encoding="utf-8")

    # 6a-bis) results/de/limma_results.tsv — research1.md line 169's schema,
    #      verbatim and in its order: accession, gene, logFC, fold_change,
    #      AveExpr, t, P.Value, adj.P.Val, B, n_imputed, direction. This is a
    #      CONTRACT file (limma's own column names), distinct from qc_limma.csv
    #      (the pipeline's audit view, with its own historical names). Nothing
    #      reads it yet; it exists so the contract is on disk and testable.
    de_results_path = None
    if write_contract:
        de_dir = os.path.join(outdir, DE_SUBDIR)
        os.makedirs(de_dir, exist_ok=True)
        logfc = pd.to_numeric(merged["limma_log2FC"])
        adj_p = pd.to_numeric(merged["adj_p_value"])
        de_results = pd.DataFrame(
            {
                "accession": merged["UniProt Accession Number"],
                "gene": merged["Gene names"],
                "logFC": logfc,
                # research1.md line 132: fold_change = 2^logFC. Derived from the
                # logFC written in the same row, so the two can never disagree.
                "fold_change": np.power(2.0, logfc),
                "AveExpr": merged["AveExpr"],
                "t": merged["t"],
                "P.Value": merged["p_value"],
                "adj.P.Val": adj_p,
                "B": merged["B"],
                "n_imputed": merged["n_imputed"].astype(int),
                # research1.md lines 140-141, verbatim: the direction call
                # requires BOTH significance and the effect-size boundary, so on
                # this dataset every row is "NS" (0/1938 pass FDR<0.05). That is
                # the honest reading and it must not be softened into an
                # effect-size-only call — see DECISIONS_LOG D2.
                "direction": np.where(
                    (adj_p < ADJ_P_THRESHOLD) & (logfc >= LOG2_THRESHOLD), "UP",
                    np.where(
                        (adj_p < ADJ_P_THRESHOLD) & (logfc <= -LOG2_THRESHOLD),
                        "DOWN", "NS",
                    ),
                ),
            }
        )
        assert list(de_results.columns) == LIMMA_RESULTS_COLUMNS, (
            f"BUG7: limma_results.tsv schema drifted: {list(de_results.columns)}"
        )
        de_results_path = os.path.join(de_dir, LIMMA_RESULTS_NAME)
        de_results.to_csv(
            de_results_path, sep="\t", index=False, encoding="utf-8",
            lineterminator="\n",
        )

    # 6b) ipa_input_significant.csv — only "regulated AND significant" rows, in the
    #     existing ipa_input.csv column layout (UniProt / Gene / log2FC / regulated)
    #     plus adj_p_value, for drop-in use. log2FC here is the pipeline's existing
    #     fold-change (UP/DOWN rows are always complete, so it is finite). Skipped
    #     entirely when write_ipa is False (the side experiment never touches IPA).
    ipa_sig_path = None
    if write_ipa:
        ipa_sig_rows = merged[merged["regulated_significant"]]
        ipa_sig = pd.DataFrame(
            {
                "UniProt Accession Number": ipa_sig_rows["UniProt Accession Number"],
                "Gene names": ipa_sig_rows["Gene names"],
                "log2FC": ipa_sig_rows["log2FC"],
                "regulated": ipa_sig_rows["regulated"],
                "adj_p_value": ipa_sig_rows["adj_p_value"],
            }
        )
        ipa_sig_path = os.path.join(outdir, "ipa_input_significant.csv")
        ipa_sig.to_csv(ipa_sig_path, index=False, encoding="utf-8")

    # 7) Summary: of the existing regulated proteins, how many are also significant
    #    (trustworthy) vs not (probably noise)?
    n_significant = int(merged["significant"].sum())
    n_regulated = int(is_regulated.sum())
    n_reg_sig = int(merged["regulated_significant"].sum())
    n_reg_not_sig = n_regulated - n_reg_sig

    min_adj_p = float(pd.to_numeric(merged["adj_p_value"]).min())
    n_raw_p_lt_05 = int((pd.to_numeric(merged["p_value"]) < 0.05).sum())

    print(f"  eBayes mode                : {ebayes_mode}")
    print(f"  proteins tested            : {eligible_count}")
    print(f"  significant (adj p < {ADJ_P_THRESHOLD}) : {n_significant}")
    print(f"  min adj p / raw p < 0.05   : {min_adj_p:.6f} / {n_raw_p_lt_05}")
    print(f"  imputed values (n_imputed) : "
          f"{dict(sorted(n_imputed.astype(int).value_counts().items()))}")
    print(f"  of {n_regulated} regulated (+/-0.585) proteins:")
    print(f"    also significant (trustworthy)   : {n_reg_sig}")
    print(f"    not significant (probably noise) : {n_reg_not_sig}")
    written = os.path.relpath(qc_path, _ROOT)
    if de_results_path is not None:
        written += f", {os.path.relpath(de_results_path, _ROOT)}"
    if ipa_sig_path is not None:
        written += f", {os.path.relpath(ipa_sig_path, _ROOT)}"
    print(f"  wrote: {written}")
    if os.path.exists(versions_path):
        print(f"  versions: {os.path.relpath(versions_path, _ROOT)}")

    summary = {
        "eligible_count": eligible_count,
        "n_tested": eligible_count,
        "ebayes_mode": ebayes_mode,
        "n_significant": n_significant,
        "n_regulated": n_regulated,
        "n_regulated_significant": n_reg_sig,
        "n_regulated_not_significant": n_reg_not_sig,
        "min_adj_p_value": min_adj_p,
        "n_raw_p_lt_0.05": n_raw_p_lt_05,
        "n_imputed_total": int(n_imputed.sum()),
    }

    # 8) The vanilla companion (D9). Same input file, same seed, same design —
    #    only the eBayes call differs — so the two runs are directly comparable
    #    and `qc_limma_vanilla.csv` stays in lockstep with `qc_limma.csv`.
    if vanilla_companion and ebayes_mode == "trend_robust":
        print(f"[Bug 7] vanilla companion (D9) -> {VANILLA_QC_FILENAME}")
        vanilla = run_limma_test(
            foldchange_csv=foldchange_csv,
            outdir=outdir,
            ebayes_mode=VANILLA_EBAYES_MODE,
            qc_filename=VANILLA_QC_FILENAME,
            output_name=VANILLA_OUTPUT,
            reuse_input=True,        # byte-identical data going in, by construction
            write_ipa=False,         # one IPA list only, from the default flavour
            write_contract=False,    # must not restate results/de/
            vanilla_companion=False,
        )

        # THE invariant D9 rests on: eBayes moderates the variance, it does not
        # refit the model, so switching flavours cannot move a single logFC. If
        # this ever trips, the change is NOT confined to the variance model and
        # every downstream fold-change conclusion is in play.
        trend_fc = pd.read_csv(qc_path)["limma_log2FC"].to_numpy(dtype=float)
        vanilla_fc = pd.read_csv(
            os.path.join(outdir, VANILLA_QC_FILENAME)
        )["limma_log2FC"].to_numpy(dtype=float)
        if not np.array_equal(trend_fc, vanilla_fc):
            worst = float(np.nanmax(np.abs(trend_fc - vanilla_fc)))
            raise RuntimeError(
                "BUG7: eBayes(trend, robust) moved logFC (max |diff| = "
                f"{worst:.3e}). It must only change the moderated t and the "
                "p-values — DECISIONS_LOG D9 does not hold and the results are "
                "not comparable."
            )

        summary["vanilla"] = vanilla
        print(f"  logFC identical across flavours : True ({len(trend_fc)} proteins)")
        print(f"  min adj p  trend_robust / vanilla : {min_adj_p:.6f} / "
              f"{vanilla['min_adj_p_value']:.6f}")
        print(f"  raw p<0.05 trend_robust / vanilla : {n_raw_p_lt_05} / "
              f"{vanilla['n_raw_p_lt_0.05']}")

    return summary


if __name__ == "__main__":
    run_limma_test()
