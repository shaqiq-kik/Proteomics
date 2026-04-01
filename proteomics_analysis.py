import pandas as pd
import numpy as np


INPUT_FILE = "Copy of General Sheet.xlsx"
OUTPUT_FILE = "alligned_proteins.xlsx"

# 1) Read both protein report sheets
df_L = pd.read_excel(INPUT_FILE, sheet_name="Protein Report L")
df_H = pd.read_excel(INPUT_FILE, sheet_name="Protein Report H")

# 2) Keep only the relev    ant columns from each sheet
cols_L = ["UniProt Accession Number", "Gene names", "Intensity 31578", "Intensity 31580"]
cols_H = ["UniProt Accession Number", "Gene names", "Intensity 31579", "Intensity 31581"]

df_L = df_L[cols_L]
df_H = df_H[cols_H]

print(f"Proteins in Protein Report L: {len(df_L)}")
print(f"Proteins in Protein Report H: {len(df_H)}")

# 3) Merge on UniProt Accession Number (inner join = found in both)
df_merged = pd.merge(
    df_L, df_H,
    on="UniProt Accession Number",
    how="inner",
    suffixes=("_L", "_H")
)

# Use Gene names from L, fall back to H if missing
df_merged["Gene names"] = df_merged["Gene names_L"].combine_first(df_merged["Gene names_H"])
df_merged = df_merged.drop(columns=["Gene names_L", "Gene names_H"])

# Reorder columns
df_merged = df_merged[
    ["UniProt Accession Number", "Gene names",
     "Intensity 31578", "Intensity 31580",
     "Intensity 31579", "Intensity 31581"]
]

print(f"Proteins found in both sheets: {len(df_merged)}")

# 4) Flag rows where any intensity value is 0 or NaN
intensity_cols = ["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581"]
df_merged["has_missing_data"] = (
    df_merged[intensity_cols].isnull().any(axis=1) |
    (df_merged[intensity_cols] == 0).any(axis=1)
)

n_missing = df_merged["has_missing_data"].sum()
print(f"Proteins with missing/zero intensity data: {n_missing}")

# 5) Add fold change columns
df_merged["fold_change_exp1"] = df_merged["Intensity 31579"] / df_merged["Intensity 31578"]
df_merged["fold_change_exp2"] = df_merged["Intensity 31581"] / df_merged["Intensity 31580"]

# 6) Calculate mean, log2, and regulation status
df_merged["mean_fold_change"] = (df_merged["fold_change_exp1"] + df_merged["fold_change_exp2"]) / 2
df_merged["log2_fold_change"] = np.log2(df_merged["mean_fold_change"])

# Set 'regulated' (UP, DOWN, or Stable)
df_merged["regulated"] = "NO CHANGE"
complete_data = ~df_merged["has_missing_data"]
df_merged.loc[complete_data & (df_merged["mean_fold_change"] >= 1.5), "regulated"] = "UP"
df_merged.loc[complete_data & (df_merged["mean_fold_change"] <= 0.5), "regulated"] = "DOWN"

# 7) Print protein regulation counts
print("UP proteins:", (df_merged["regulated"] == "UP").sum())
print("DOWN proteins:", (df_merged["regulated"] == "DOWN").sum())
print("NO CHANGE proteins:", (df_merged["regulated"] == "NO CHANGE").sum())

# 8) Sort: False (complete) first, True (missing) last, then save
df_merged = df_merged.sort_values("has_missing_data", ascending=True)
df_merged.to_excel(OUTPUT_FILE, index=False)
print(f"\nResults saved to {OUTPUT_FILE}")

# filter out missing data and update IPA input file
IPA_FILE = "IPA_input.xlsx"

df_ipa = df_merged[
    (~df_merged["has_missing_data"]) & (df_merged["regulated"] != "NO CHANGE")
].copy()

df_ipa.to_excel(IPA_FILE, index=False)

print(f"\nFiltered data saved to {IPA_FILE}")
print(f"Rown in new IPA input file: {len(df_ipa)}")

