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
"""

import os

import numpy as np
import pandas as pd

INPUT_FILE = "Copy of General Sheet.xlsx"
RESULTS_DIR = os.path.join("proteomics_de", "results")

LOG2_THRESHOLD = 0.585  # = log2(1.5); the down side is the symmetric -log2(1.5)

INTENSITY_COLS = ["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581"]


if __name__ == "__main__":
    # 1) Read both protein report sheets, keeping only the relevant columns
    cols_L = ["UniProt Accession Number", "Gene names", "Intensity 31578", "Intensity 31580"]
    cols_H = ["UniProt Accession Number", "Gene names", "Intensity 31579", "Intensity 31581"]

    df_L = pd.read_excel(INPUT_FILE, sheet_name="Protein Report L")[cols_L]
    df_H = pd.read_excel(INPUT_FILE, sheet_name="Protein Report H")[cols_H]

    print(f"Proteins in Protein Report L: {len(df_L)}")
    print(f"Proteins in Protein Report H: {len(df_H)}")

    # 2) Outer-merge on accession with an indicator so proteins detected in only ONE
    #    sheet are NOT silently dropped (Bug 4). Coalesce Gene names (L first, then H).
    merged = pd.merge(
        df_L, df_H, on="UniProt Accession Number", how="outer",
        indicator=True, suffixes=("_L", "_H"),
    )
    merged["Gene names"] = merged["Gene names_L"].combine_first(merged["Gene names_H"])

    # Split: "both" feeds the EXISTING fold-change pipeline unchanged; single-condition
    # rows (left_only / right_only) are rescued to their own file (no fold change).
    both = merged[merged["_merge"] == "both"].copy()
    single_cond = merged[merged["_merge"] != "both"].copy()

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
    left_order = {acc: i for i, acc in enumerate(df_L["UniProt Accession Number"])}
    both = both.sort_values(
        "UniProt Accession Number", key=lambda s: s.map(left_order), kind="stable"
    )
    df = both[["UniProt Accession Number", "Gene names"] + INTENSITY_COLS].reset_index(drop=True)
    for col in ["Intensity 31579", "Intensity 31581"]:  # Heavy/treated replicates
        df[col] = df[col].astype(df_H[col].dtype)

    # 3) Completeness: a row is INCOMPLETE if any of the 4 intensities is 0 or NaN
    df["complete"] = ~(
        df[INTENSITY_COLS].isnull().any(axis=1) | (df[INTENSITY_COLS] == 0).any(axis=1)
    )
    print(f"Complete proteins (all 4 intensities present & non-zero): {df['complete'].sum()}")

    # 4) Per-replicate SILAC ratios (Heavy / Light), only on complete rows.
    #    Restricting to complete rows means denominators are never 0, so no inf can
    #    be produced (Bug 3).
    df["ratio_rep1"] = np.nan
    df["ratio_rep2"] = np.nan
    mask = df["complete"]
    df.loc[mask, "ratio_rep1"] = df.loc[mask, "Intensity 31579"] / df.loc[mask, "Intensity 31578"]
    df.loc[mask, "ratio_rep2"] = df.loc[mask, "Intensity 31581"] / df.loc[mask, "Intensity 31580"]

    # 5) Bug 1 fix: log2 each ratio first, THEN average the logs.
    df["log2_rep1"] = np.log2(df["ratio_rep1"])
    df["log2_rep2"] = np.log2(df["ratio_rep2"])
    df["log2FC"] = df[["log2_rep1", "log2_rep2"]].mean(axis=1)

    # 6) Bug 2 fix: symmetric cutoffs in log2 space (only on complete rows)
    df["regulated"] = "NO CHANGE"
    df.loc[mask & (df["log2FC"] >= LOG2_THRESHOLD), "regulated"] = "UP"
    df.loc[mask & (df["log2FC"] <= -LOG2_THRESHOLD), "regulated"] = "DOWN"

    # 7) Bug 4b: "cousin" on/off proteins. A protein present in BOTH sheets but whose
    #    entire control OR treated side is absent (both intensities 0/NaN) is a real
    #    on/off signal, yet the Bug 3 completeness logic parks it as incomplete and it
    #    falls into NO CHANGE — a wrong label. Detect these and relabel them ON_OFF.
    #    A protein absent on BOTH sides is fully empty and stays NO CHANGE.
    CONTROL_COLS = ["Intensity 31578", "Intensity 31580"]  # light / control replicates
    TREATED_COLS = ["Intensity 31579", "Intensity 31581"]  # heavy / treated replicates
    control_absent = (df[CONTROL_COLS].isnull() | (df[CONTROL_COLS] == 0)).all(axis=1)
    treated_absent = (df[TREATED_COLS].isnull() | (df[TREATED_COLS] == 0)).all(axis=1)
    control_off = control_absent & ~treated_absent  # present in treated only -> "on"
    treated_off = treated_absent & ~control_absent  # present in control only -> "off"

    df["onoff"] = ""
    df.loc[control_off, "onoff"] = "on_with_treatment"
    df.loc[treated_off, "onoff"] = "off_with_treatment"

    # Relabel out of NO CHANGE. No log2FC is computed for on/off proteins (they are
    # incomplete, so log2FC is already NaN).
    df.loc[control_off | treated_off, "regulated"] = "ON_OFF"

    # Sanity checks
    n_up = (df["regulated"] == "UP").sum()
    n_down = (df["regulated"] == "DOWN").sum()
    n_nc = (df["regulated"] == "NO CHANGE").sum()
    n_onoff = (df["regulated"] == "ON_OFF").sum()
    n_on = (df["onoff"] == "on_with_treatment").sum()
    n_off = (df["onoff"] == "off_with_treatment").sum()
    print(f"on_with_treatment: {n_on}")
    print(f"off_with_treatment: {n_off}")
    print(f"total on/off: {n_onoff}")
    print(f"UP: {n_up}")
    print(f"DOWN: {n_down}")
    print(f"NO CHANGE: {n_nc}")
    print(f"ON_OFF: {n_onoff}")

    # Bug 4 + Bug 4b asserts: UP/DOWN are untouched, and on/off proteins ONLY move from
    # NO CHANGE to ON_OFF — nothing else should shift between buckets.
    assert n_up == 206, f"UP changed to {n_up}, expected 206"
    assert n_down == 509, f"DOWN changed to {n_down}, expected 509"
    assert n_onoff == n_on + n_off, "ON_OFF total != on_with_treatment + off_with_treatment"
    assert n_nc + n_onoff == 1233, (
        f"NO CHANGE + ON_OFF = {n_nc + n_onoff}, expected 1233; the NO CHANGE drop "
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

    # Output
    os.makedirs(RESULTS_DIR, exist_ok=True)

    out_cols = [
        "UniProt Accession Number", "Gene names",
        "Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581",
        "ratio_rep1", "ratio_rep2", "log2_rep1", "log2_rep2", "log2FC",
        "complete", "regulated", "onoff",
    ]
    foldchange_path = os.path.join(RESULTS_DIR, "foldchange_all.csv")
    df[out_cols].to_csv(foldchange_path, index=False)
    print(f"\nSaved {foldchange_path} ({len(df)} rows)")

    # IPA input: complete AND regulated; accession leftmost. UTF-8, no BOM.
    df_ipa = df[df["complete"] & (df["regulated"] != "NO CHANGE")].copy()
    ipa_cols = ["UniProt Accession Number", "Gene names", "log2FC", "regulated"]
    ipa_path = os.path.join(RESULTS_DIR, "ipa_input.csv")
    df_ipa[ipa_cols].to_csv(ipa_path, index=False, encoding="utf-8")
    print(f"Saved {ipa_path} ({len(df_ipa)} rows)")

    assert len(df_ipa) == 715, f"ipa_input.csv has {len(df_ipa)} rows, expected 715"

    # Bug 4 rescue file: single-condition proteins (no fold change computed). The
    # absent condition's intensity columns stay blank/NaN, which is expected and
    # visually shows the protein was off in that condition.
    single_cond["accession"] = single_cond["UniProt Accession Number"]
    single_cond["gene"] = single_cond["Gene names"]
    single_cond["detected_in"] = np.where(
        single_cond["_merge"] == "left_only", "control_only", "treated_only"
    )
    single_cols = [
        "accession", "gene", "detected_in",
        "Intensity 31578", "Intensity 31580",  # control replicates
        "Intensity 31579", "Intensity 31581",  # treated replicates
    ]
    single_path = os.path.join(RESULTS_DIR, "single_condition_proteins.csv")
    single_cond[single_cols].to_csv(single_path, index=False, encoding="utf-8")
    print(f"Saved {single_path} ({len(single_cond)} rows)")

    # Bug 4b file: "cousin" on/off proteins (present in both sheets, one side absent).
    onoff_df = df[df["regulated"] == "ON_OFF"].copy()
    onoff_df["accession"] = onoff_df["UniProt Accession Number"]
    onoff_df["gene"] = onoff_df["Gene names"]
    onoff_cols = [
        "accession", "gene", "onoff",
        "Intensity 31578", "Intensity 31580",  # control replicates
        "Intensity 31579", "Intensity 31581",  # treated replicates
    ]
    onoff_path = os.path.join(RESULTS_DIR, "onoff_proteins.csv")
    onoff_df[onoff_cols].to_csv(onoff_path, index=False, encoding="utf-8")
    print(f"Saved {onoff_path} ({len(onoff_df)} rows)")

    # Bug 5 — SILAC centering QC. Runs AFTER foldchange_all.csv is written and
    # only ever writes NEW files (qc_centering.csv, and on WARN
    # foldchange_all_centered.csv); the four outputs above are untouched.
    from centering_check import run_centering_check
    run_centering_check(foldchange_path, RESULTS_DIR)
