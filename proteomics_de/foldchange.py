"""Corrected SILAC fold-change analysis (proteomics_de).

Rebuilds the fold-change logic from the legacy ``proteomics_analysis.py`` while
fixing three correctness bugs:

  Bug 1 — log2 each replicate ratio FIRST, then average the logs (the legacy
          script averaged the raw ratios and logged once).
  Bug 2 — symmetric up/down cutoffs in log2 space (the legacy script used
          ratio >= 1.5 for UP but ratio <= 0.5 for DOWN, which is asymmetric).
  Bug 3 — never emit inf/-inf: ratios are only computed on complete rows, so
          division by zero can't occur, and we assert no inf/NaN survives.

The legacy script is the frozen baseline and is NOT modified.

This file is **wiring**, not logic. Every transformation lives in
:mod:`proteomics_de.etl.foldchange_core` as an ordinary function over
DataFrames, and the merge cardinality guards live in
:mod:`proteomics_de.etl.merge_guard`; ``main`` reads the workbook, calls them in
order, and writes the CSVs. Previously all of it sat inside ``if __name__ ==
"__main__"``, which is why none of it could be tested.

Running it
----------
Every path is resolved from ``__file__``, so the script works from any working
directory::

    python proteomics_de/foldchange.py            # from the repo root
    python -m proteomics_de.foldchange            # as a module
    python /abs/path/proteomics_de/foldchange.py  # from /tmp, or anywhere

``--input`` / ``--results-dir`` / ``--sample-sheet`` override the defaults; the
defaults reproduce the committed run exactly.

Assertions
----------
Two kinds, deliberately kept apart:

* **Structural guards** (``etl/merge_guard.py``) are derived from the data in
  hand and hold for any dataset.
* **Frozen expectations** (UP=206, DOWN=509, ...) are specific to the committed
  workbook and are read from ``tests/expected/frozen_counts.json``, not typed
  into this file. They are ON by default because they have caught real
  regressions; a genuinely new dataset sets ``PDE_EXPECT_BASELINE=0`` and
  re-derives the JSON.
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np

# --- Path resolution --------------------------------------------------------
# Resolved from __file__, never from the cwd. The pipeline is invoked from the
# repo root, from proteomics_de/, and (in tests) from a temp dir; the previous
# cwd-relative INPUT_FILE / RESULTS_DIR silently produced a FileNotFoundError
# anywhere but the repo root.
_HERE = Path(__file__).resolve().parent          # proteomics_de/
_ROOT = _HERE.parent                             # repo root

# `_ROOT` makes `proteomics_de.*` importable; `_HERE` makes the flat sibling
# imports at the bottom of main (`from centering_check import ...`) resolve under
# `-m` and from any cwd -- a bare `import centering_check` only ever worked by
# accident of sys.path[0] being the script's directory. Inserted in this order so
# the result is [_ROOT, _HERE, ...], which is exactly the sys.path the committed
# run ends up with today (script dir is _HERE; centering_check.py prepends _ROOT).
for _entry in (str(_HERE), str(_ROOT)):
    if _entry not in sys.path:
        sys.path.insert(0, _entry)

# LOG2_THRESHOLD (= log2(1.5), with the symmetric -log2(1.5) on the down side) is
# re-exported from this module: centering_check.py imports the cutoff from here
# rather than re-typing 0.585. The value itself now comes from config/constants.py,
# so there is exactly one literal for it in the tree.
from proteomics_de.config import design  # noqa: E402
from proteomics_de.config.constants import LOG2_THRESHOLD  # noqa: E402,F401
from proteomics_de.etl import foldchange_core as core  # noqa: E402
from proteomics_de.etl import merge_guard  # noqa: E402
from proteomics_de.export.ipa_export import write_ipa  # noqa: E402
from proteomics_de.qc import boundaries  # noqa: E402

INPUT_FILE = str(_ROOT / "Copy of General Sheet.xlsx")
RESULTS_DIR = str(_HERE / "results")

# Two different things used to be conflated in one pair of constants, and
# separating them is what makes DECISIONS_LOG D7 a one-line change:
#
#   * WHICH SHEET a sample lives in is a fact about the workbook. Samples 31578
#     and 31580 are in "Protein Report L", 31579 and 31581 are in "Protein
#     Report H". No experimental relabelling can move them, so COLS_L/COLS_H
#     stay literal.
#   * WHICH CONDITION a sample is is an experimental fact, and it belongs to
#     config/sample_sheet.tsv. D7 records that this assignment shipped
#     inverted, so it must be derivable, not typed in here.
#
# PHYSICAL_COLS is the acquisition order and is used only for file I/O, so
# output column order is stable across a relabelling -- a D7-style flip changes
# the sign of log2FC and the UP/DOWN labels, not the shape of any CSV.
PHYSICAL_COLS = ["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581"]
INTENSITY_COLS = PHYSICAL_COLS  # back-compat alias for importers

SHEET_L = "Protein Report L"
SHEET_H = "Protein Report H"

SHEET_L_COLS = ["Intensity 31578", "Intensity 31580"]  # channels living in sheet L
SHEET_H_COLS = ["Intensity 31579", "Intensity 31581"]  # channels living in sheet H

COLS_L = ["UniProt Accession Number", "Gene names"] + SHEET_L_COLS
COLS_H = ["UniProt Accession Number", "Gene names"] + SHEET_H_COLS

OUT_COLS = [
    "UniProt Accession Number", "Gene names",
    "Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581",
    "ratio_rep1", "ratio_rep2", "log2_rep1", "log2_rep2", "log2FC",
    "complete", "regulated", "onoff",
]
IPA_COLS = ["UniProt Accession Number", "Gene names", "log2FC", "regulated"]
SINGLE_COLS = [
    "accession", "gene", "detected_in",
    "Intensity 31578", "Intensity 31580",  # treated (testosterone) replicates
    "Intensity 31579", "Intensity 31581",  # control (vehicle) replicates
]
ONOFF_COLS = [
    "accession", "gene", "onoff",
    "Intensity 31578", "Intensity 31580",  # treated (testosterone) replicates
    "Intensity 31579", "Intensity 31581",  # control (vehicle) replicates
]


def _display(path) -> str:
    """Repo-relative rendering of `path` for the console, absolute if outside.

    Keeps the "Saved ..." lines reading ``proteomics_de/results/...`` now that
    the paths themselves are absolute, and matches how ``centering_check.py`` and
    ``replicate_check.py`` already report their outputs.
    """
    try:
        rel = os.path.relpath(str(path), str(_ROOT))
    except ValueError:  # different drive (Windows); no sensible relative form
        return str(path)
    return str(path) if rel.startswith(os.pardir) else rel


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--input", default=INPUT_FILE,
                    help="input workbook (default: %(default)s)")
    ap.add_argument("--results-dir", default=RESULTS_DIR,
                    help="directory for the output CSVs (default: %(default)s)")
    ap.add_argument("--sample-sheet", default=None,
                    help="sample sheet to validate the intensity columns against "
                         "(default: proteomics_de/config/sample_sheet.tsv)")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    input_file = args.input
    results_dir = args.results_dir

    # The boundary hooks below are called as `check(stage, df)` -- the Wave-0
    # two-argument contract, unchanged. This is how they learn about
    # --results-dir, so a run into a temp directory does not write its QC
    # records into the committed tree.
    boundaries.set_results_dir(results_dir)

    # 0) Resolve the condition assignment from config/sample_sheet.tsv. The
    #    SHAPE of this stage is still not design-driven -- the L/H sheet split
    #    and the left-order/Heavy-dtype restore below are inherently 2-channel
    #    SILAC, so a genuinely new design means regenerating this stage. But
    #    WHICH samples are control and which are treated is read from the sheet,
    #    which is what makes the D7 correction a one-line edit to that TSV.
    sheet = design.read_sample_sheet(args.sample_sheet)
    control_cols = design.control_columns(sheet)
    treated_cols = design.treated_columns(sheet)

    # The sheet must still describe THIS workbook: same four channels, and two
    # per group so the ratio pairing below is well defined. assert_matches turns
    # "silently analysing the wrong four columns" into a loud, actionable failure.
    if sorted(control_cols + treated_cols) != sorted(PHYSICAL_COLS):
        raise ValueError(
            f"design drift in foldchange.py: sample sheet resolves to "
            f"{control_cols + treated_cols}, workbook provides {PHYSICAL_COLS}"
        )
    if len(control_cols) != 2 or len(treated_cols) != 2:
        raise ValueError(
            "foldchange.py is 2-channel SILAC-specific: it pairs replicate i of "
            f"control with replicate i of treated, so it needs exactly 2 samples "
            f"per group, got {len(control_cols)} control / {len(treated_cols)} treated."
        )
    print(f"Control replicates: {control_cols}")
    print(f"Treated replicates: {treated_cols}")

    # 1) Read both protein report sheets, keeping only the relevant columns
    df_L, df_H = core.read_sheets(input_file, COLS_L, COLS_H,
                                  sheet_L=SHEET_L, sheet_H=SHEET_H)
    boundaries.check("after_load", df_L)
    boundaries.check("after_load", df_H)

    print(f"Proteins in Protein Report L: {len(df_L)}")
    print(f"Proteins in Protein Report H: {len(df_H)}")

    # 2) Outer-merge on accession with an indicator so proteins detected in only ONE
    #    sheet are NOT silently dropped (Bug 4). Coalesce Gene names (L first, then H).
    merged = core.merge_with_indicator(df_L, df_H)
    boundaries.check("after_merge", merged)

    # Split: "both" feeds the EXISTING fold-change pipeline unchanged; single-condition
    # rows (left_only / right_only) are rescued to their own file (no fold change).
    both, single_cond = core.split_both_single(merged)

    # Cardinality guard (research1.md lines 52, 183): duplicate accessions on
    # either side turn the merge into a cross product, and the extra rows are
    # indistinguishable from real proteins downstream. Derived entirely from the
    # frames in hand, so it stays valid for any dataset.
    merge_guard.assert_merge_safe(df_L, df_H, merged, both)

    n_both = (merged["_merge"] == "both").sum()
    n_control_only = (merged["_merge"] == "left_only").sum()
    n_treated_only = (merged["_merge"] == "right_only").sum()
    print(f"Proteins found in both sheets: {n_both}")
    print(f"Single-condition proteins (control_only): {n_control_only}")
    print(f"Single-condition proteins (treated_only): {n_treated_only}")
    print(f"Single-condition proteins (total): {len(single_cond)}")

    # The fold-change pipeline below operates ONLY on the `both` group, so every
    # current result stays identical to the previous inner-join behaviour. The outer
    # join has two side effects we must undo so the existing outputs stay byte-for-byte
    # unchanged: (a) it reorders the matched rows, and (b) its NaNs (from the rescued
    # single-condition rows) coerce the Heavy intensity columns to float. Restore the
    # original Light-sheet row order and the original Heavy dtypes.
    df = core.restore_left_order(
        both, df_L, df_H,
        out_cols=["UniProt Accession Number", "Gene names"] + INTENSITY_COLS,
        dtype_cols=SHEET_H_COLS,  # sheet-H channels: dtype, not condition
    )

    # 3) Completeness: a row is INCOMPLETE if any of the 4 intensities is 0 or NaN
    df = core.mark_complete(df, INTENSITY_COLS)
    print(f"Complete proteins (all 4 intensities present & non-zero): {df['complete'].sum()}")

    # 4) Per-replicate SILAC ratios (Heavy / Light), only on complete rows.
    #    Restricting to complete rows means denominators are never 0, so no inf can
    #    be produced (Bug 3).
    mask = df["complete"]
    df, ratio_cols = core.compute_ratios(df, control_cols, treated_cols, mask)

    # 5) Bug 1 fix: log2 each ratio first, THEN average the logs.
    df, _log_cols = core.compute_log2fc(df, ratio_cols)

    # 6) Bug 2 fix: symmetric cutoffs in log2 space (only on complete rows)
    df = core.classify_regulated(df, mask, LOG2_THRESHOLD)

    # 7) Bug 4b: "cousin" on/off proteins. A protein present in BOTH sheets but whose
    #    entire control OR treated side is absent (both intensities 0/NaN) is a real
    #    on/off signal, yet the Bug 3 completeness logic parks it as incomplete and it
    #    falls into NO CHANGE — a wrong label. Detect these and relabel them ON_OFF.
    #    A protein absent on BOTH sides is fully empty and stays NO CHANGE.
    df, _control_off, _treated_off = core.detect_onoff(df, control_cols, treated_cols)
    boundaries.check("after_foldchange", df)

    # Sanity checks
    counts = core.summarize(df)
    n_up = counts["n_up"]
    n_down = counts["n_down"]
    n_nc = counts["n_nochange"]
    n_onoff = counts["n_onoff"]
    n_on = counts["n_on"]
    n_off = counts["n_off"]
    print(f"on_with_treatment: {n_on}")
    print(f"off_with_treatment: {n_off}")
    print(f"total on/off: {n_onoff}")
    print(f"UP: {n_up}")
    print(f"DOWN: {n_down}")
    print(f"NO CHANGE: {n_nc}")
    print(f"ON_OFF: {n_onoff}")

    # Dataset-specific expectations, read from tests/expected/frozen_counts.json
    # rather than typed here. `expect` is None when PDE_EXPECT_BASELINE=0.
    expect = core.load_frozen_counts() if core.baseline_checks_enabled() else None

    # Bug 4 + Bug 4b asserts: UP/DOWN are untouched, and on/off proteins ONLY move from
    # NO CHANGE to ON_OFF — nothing else should shift between buckets.
    if expect is not None:
        assert n_up == expect["n_up"], f"UP changed to {n_up}, expected {expect['n_up']}"
        assert n_down == expect["n_down"], f"DOWN changed to {n_down}, expected {expect['n_down']}"
    assert n_onoff == n_on + n_off, "ON_OFF total != on_with_treatment + off_with_treatment"
    if expect is not None:
        assert n_nc + n_onoff == expect["n_nochange_plus_onoff"], (
            f"NO CHANGE + ON_OFF = {n_nc + n_onoff}, expected "
            f"{expect['n_nochange_plus_onoff']}; the NO CHANGE drop "
            "does not equal the on/off total, so something other than on/off proteins moved."
        )
    # On/off proteins must NOT carry a numeric log2FC.
    assert df.loc[df["regulated"] == "ON_OFF", "log2FC"].isnull().all(), (
        "An ON_OFF protein has a numeric log2FC."
    )

    # Bug 3 assert: no inf or NaN in log2FC for complete rows
    complete_log2fc = df.loc[mask, "log2FC"]
    assert not np.isinf(complete_log2fc).any(), "Found inf in log2FC for complete rows"
    assert not complete_log2fc.isnull().any(), "Found NaN in log2FC for complete rows"

    # Every row carries exactly one of UP / DOWN / NO CHANGE / ON_OFF. Catches a
    # stray or missing label, which would otherwise surface only as a downstream
    # file being mysteriously short (research1.md line 183).
    merge_guard.assert_classification_partition(df)

    # Output
    os.makedirs(results_dir, exist_ok=True)

    foldchange_path = os.path.join(results_dir, "foldchange_all.csv")
    df[OUT_COLS].to_csv(foldchange_path, index=False)
    print(f"\nSaved {_display(foldchange_path)} ({len(df)} rows)")

    # IPA input: complete AND regulated; accession leftmost. UTF-8, no BOM.
    # The LIVE frame is handed to the writer — never a re-read of the CSV just
    # written. Round-tripping through text perturbs ~30 rows in the last float
    # ULP, which is invisible in a diff and fatal to the byte-freeze.
    df_ipa = core.build_ipa_frame(df)
    ipa_path = os.path.join(results_dir, "ipa_input.csv")
    write_ipa(df_ipa, ipa_path, IPA_COLS)
    print(f"Saved {_display(ipa_path)} ({len(df_ipa)} rows)")

    if expect is not None:
        assert len(df_ipa) == expect["ipa_input_rows"], (
            f"ipa_input.csv has {len(df_ipa)} rows, expected {expect['ipa_input_rows']}"
        )

    # Bug 4 rescue file: single-condition proteins (no fold change computed). The
    # absent condition's intensity columns stay blank/NaN, which is expected and
    # visually shows the protein was off in that condition.
    # Sheet L holds SHEET_L_COLS and sheet H holds SHEET_H_COLS; which CONDITION
    # each of those is comes from the sample sheet, so the detected_in labels are
    # resolved rather than assumed (see D7).
    left_condition = "control" if SHEET_L_COLS[0] in control_cols else "treated"
    right_condition = "treated" if left_condition == "control" else "control"
    single_cond = core.build_single_condition_frame(
        single_cond, left_condition=left_condition, right_condition=right_condition
    )

    # DECISIONS_LOG D11 — quarantine, don't ship. Two rows of the raw Light
    # sheet carry a ';'-joined list of bare MaxQuant row indices where an
    # accession belongs (32,759 and 681 characters). They are not accessions,
    # and qc/schema.py used to carve an exception so they would PASS
    # validation; they are now written, in full and with a reason, to
    # results/qc/quarantine_accessions.csv and dropped from this file
    # (606 -> 604 rows).
    #
    # No numeric logic above changes: these rows are single-condition, so they
    # never had a fold change, never entered limma, and never reached the IPA
    # export. Only the enrichment BACKGROUND shrinks with them (2554 -> 2552),
    # which is read from tests/expected/frozen_counts.json, not from here.
    single_cond, quarantined = boundaries.quarantine_junk_accessions(
        single_cond,
        "accession",
        source="single_condition_proteins.csv",
        results_dir=results_dir,
    )
    if len(quarantined):
        print(
            f"Quarantined {len(quarantined)} junk accession(s) -> "
            f"{_display(boundaries.quarantine_path(results_dir))}"
        )

    single_path = os.path.join(results_dir, "single_condition_proteins.csv")
    single_cond[SINGLE_COLS].to_csv(single_path, index=False, encoding="utf-8")
    print(f"Saved {_display(single_path)} ({len(single_cond)} rows)")

    # Bug 4b file: "cousin" on/off proteins (present in both sheets, one side absent).
    onoff_df = core.build_onoff_frame(df)
    onoff_path = os.path.join(results_dir, "onoff_proteins.csv")
    onoff_df[ONOFF_COLS].to_csv(onoff_path, index=False, encoding="utf-8")
    print(f"Saved {_display(onoff_path)} ({len(onoff_df)} rows)")

    # Bug 5 — SILAC centering QC. Runs AFTER foldchange_all.csv is written and
    # only ever writes NEW files (qc_centering.csv, and on WARN
    # foldchange_all_centered.csv); the four outputs above are untouched.
    from centering_check import run_centering_check
    run_centering_check(foldchange_path, results_dir)

    # Bug 6 — replicate correlation QC. Pure read of foldchange_all.csv; writes
    # only NEW files (qc_replicate_correlation.csv, replicate_correlation.png).
    # Local import keeps the module-level import graph acyclic (same as Bug 5).
    from replicate_check import run_replicate_correlation
    run_replicate_correlation(foldchange_path, results_dir)

    # Bug 7 — per-protein statistical testing via limma (R) + MinProb imputation.
    # Reads foldchange_all.csv from disk, shells out to limma_test.R for the stats,
    # and writes only NEW files (qc_limma.csv, ipa_input_significant.csv). The local
    # import keeps the module graph acyclic (same as Bug 5/6); limma_test.py does
    # not import from foldchange.
    from limma_test import run_limma_test
    run_limma_test(foldchange_csv=foldchange_path, outdir=results_dir)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
