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

# Sanity checks
n_up = (df["regulated"] == "UP").sum()
n_down = (df["regulated"] == "DOWN").sum()
n_nc = (df["regulated"] == "NO CHANGE").sum()
print(f"UP: {n_up}")
print(f"DOWN: {n_down}")
print(f"NO CHANGE: {n_nc}")

# Bug 4 assert: rescuing single-condition proteins must NOT leak any of them into
# the fold-change logic, so the `both` group's classification must be unchanged.
assert (n_up, n_down, n_nc) == (206, 509, 1233), (
    f"Classification on `both` changed to {n_up}/{n_down}/{n_nc}; "
    "single-condition rows likely leaked into the fold-change split."
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
    "complete", "regulated",
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
