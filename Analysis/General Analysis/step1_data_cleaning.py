"""
STEP 1.1: Data Cleaning & Formatting for Proteomics Analysis
Processes SILAC proteomics data comparing Vehicle Control vs Testosterone treatment
"""

import pandas as pd
import numpy as np
import re

print("="*80)
print("STEP 1.1: PROTEOMICS DATA CLEANING & FORMATTING")
print("="*80)

# ============================================================================
# STEP 1: DATA LOADING AND CLEANING
# ============================================================================
print("\n[STEP 1] Loading and cleaning data...")

# Load the Excel file with proper header row
df = pd.read_excel('../CLEANED Silac Proteomics Soluble Factors.xlsx', header=0)

print(f"Initial shape: {df.shape}")
print(f"Initial columns: {list(df.columns)}")

# Function to clean numeric values
def clean_numeric(val):
    """Remove commas, spaces, handle special values"""
    if pd.isna(val):
        return np.nan

    # Convert to string
    val_str = str(val).strip()

    # Handle special cases
    if val_str in ['', '-', ' - ', '#DIV/0!', 'nan', 'NaN']:
        return np.nan

    # Remove commas
    val_str = val_str.replace(',', '')

    # Try to convert to float
    try:
        return float(val_str)
    except:
        return np.nan

# Function to clean text values
def clean_text(val):
    """Strip whitespace, handle empty/missing values"""
    if pd.isna(val):
        return np.nan

    val_str = str(val).strip()

    if val_str in ['', '-', ' - ', 'nan', 'NaN']:
        return np.nan

    return val_str

# Define numeric columns (using actual column names from file)
numeric_cols = [
    'log_10_fold_change', 'log_2_fold_change',
    'Vehicle_Rep1_31579', 'Vehicle_Rep2_31581',
    'Vehicle_Mean', 'Vehicle_SD', 'Vehicle_SD_Percent',
    'Testosterone_Rep1_31578', 'Testosterone_Rep2_31580',
    'Testosterone_Mean', 'Testosterone_SD', 'Testosterone_SD_Percent',
    'Fold_Change'
]

# Define text columns
text_cols = [
    'UniProt_Accession', 'Gene', 'Protein_Description', 'UniProt_ID',
    'Alternate_Names', 'Mouse_Gene_ID', 'Cellular_Location', 'Functional_Class'
]

# Clean numeric columns
print("  - Cleaning numeric columns (removing commas, handling missing values)...")
for col in numeric_cols:
    if col in df.columns:
        df[col] = df[col].apply(clean_numeric)

# Clean text columns
print("  - Cleaning text columns (stripping whitespace)...")
for col in text_cols:
    if col in df.columns:
        df[col] = df[col].apply(clean_text)

# Standardize gene symbols to uppercase
print("  - Standardizing gene symbols to uppercase...")
if 'Gene' in df.columns:
    df['Gene'] = df['Gene'].apply(lambda x: x.upper() if pd.notna(x) else np.nan)

# Validate UniProt IDs
print("  - Validating UniProt ID format...")
if 'UniProt_Accession' in df.columns:
    df['UniProt_Valid'] = df['UniProt_Accession'].apply(
        lambda x: bool(re.match(r'^[PQ]\d{5}$', str(x))) if pd.notna(x) else False
    )
    valid_count = df['UniProt_Valid'].sum()
    print(f"    Valid UniProt IDs: {valid_count}/{len(df)}")

print(f"✓ Data cleaned. Final shape: {df.shape}")

# ============================================================================
# STEP 2: SUMMARY STATISTICS
# ============================================================================
print("\n[STEP 2] Calculating summary statistics...")

total_proteins = len(df)
upregulated = len(df[df['Fold_Change'] > 1.0])
downregulated = len(df[df['Fold_Change'] < 1.0])
unchanged = len(df[(df['Fold_Change'] >= 0.95) & (df['Fold_Change'] <= 1.05)])

# Fold change statistics
fc_stats = df['Fold_Change'].describe()

# Create summary dictionary
summary_stats = {
    'Metric': [
        'Total Proteins',
        'Upregulated (FC > 1)',
        'Downregulated (FC < 1)',
        'Unchanged (FC ≈ 1)',
        'Fold Change - Mean',
        'Fold Change - Median',
        'Fold Change - Std Dev',
        'Fold Change - Min',
        'Fold Change - Max',
    ],
    'Value': [
        total_proteins,
        upregulated,
        downregulated,
        unchanged,
        fc_stats['mean'],
        fc_stats['50%'],
        fc_stats['std'],
        fc_stats['min'],
        fc_stats['max'],
    ]
}

summary_df = pd.DataFrame(summary_stats)

print("\nSUMMARY STATISTICS:")
print(summary_df.to_string(index=False))

# Functional class distribution
print("\nFUNCTIONAL CLASS DISTRIBUTION:")
func_dist = df['Functional_Class'].value_counts()
func_dist_pct = df['Functional_Class'].value_counts(normalize=True) * 100

for func_class, count in func_dist.items():
    pct = func_dist_pct[func_class]
    print(f"  {func_class}: {count} ({pct:.1f}%)")

# ============================================================================
# STEP 3: CREATE CLASSIFICATION FILES
# ============================================================================
print("\n[STEP 3] Creating classification dataframes...")

# Upregulated proteins
df_upregulated = df[df['Fold_Change'] > 1.0].copy()
df_upregulated = df_upregulated.sort_values('Fold_Change', ascending=False)
print(f"  - Upregulated proteins: {len(df_upregulated)}")

# Downregulated proteins
df_downregulated = df[df['Fold_Change'] < 1.0].copy()
df_downregulated = df_downregulated.sort_values('Fold_Change', ascending=True)
print(f"  - Downregulated proteins: {len(df_downregulated)}")

# Top 10 upregulated
df_top10_up = df_upregulated.nlargest(10, 'Fold_Change') if len(df_upregulated) >= 10 else df_upregulated
print(f"  - Top upregulated: {len(df_top10_up)}")

# Top 10 downregulated
df_top10_down = df_downregulated.nsmallest(10, 'Fold_Change') if len(df_downregulated) >= 10 else df_downregulated
print(f"  - Top downregulated: {len(df_top10_down)}")

# ============================================================================
# STEP 4: EXPORT CLEAN FILES
# ============================================================================
print("\n[STEP 4] Exporting CSV files...")

# 1. Main cleaned dataset
df.to_csv('cleaned_proteomics_data.csv', index=False)
print("  ✓ cleaned_proteomics_data.csv")

# 2. Summary statistics
summary_df.to_csv('summary_statistics.csv', index=False)
print("  ✓ summary_statistics.csv")

# 3. Upregulated proteins
df_upregulated.to_csv('upregulated_proteins.csv', index=False)
print("  ✓ upregulated_proteins.csv")

# 4. Downregulated proteins
df_downregulated.to_csv('downregulated_proteins.csv', index=False)
print("  ✓ downregulated_proteins.csv")

# 5. Top changers (combined top 10 up and down)
df_top10_up_subset = df_top10_up[['Gene', 'Protein_Description', 'Fold_Change',
                                    'log_2_fold_change', 'Functional_Class']].copy()
df_top10_up_subset['Regulation'] = 'Upregulated'

df_top10_down_subset = df_top10_down[['Gene', 'Protein_Description', 'Fold_Change',
                                        'log_2_fold_change', 'Functional_Class']].copy()
df_top10_down_subset['Regulation'] = 'Downregulated'

df_top_changers = pd.concat([df_top10_up_subset, df_top10_down_subset])
df_top_changers.to_csv('top_changers.csv', index=False)
print("  ✓ top_changers.csv")

# 6. Functional class distribution with regulation counts
func_class_stats = []
for func_class in df['Functional_Class'].dropna().unique():
    subset = df[df['Functional_Class'] == func_class]
    total_count = len(subset)
    up_count = len(subset[subset['Fold_Change'] > 1.0])
    down_count = len(subset[subset['Fold_Change'] < 1.0])
    percentage = (total_count / total_proteins) * 100

    func_class_stats.append({
        'Functional_Class': func_class,
        'Count': total_count,
        'Percentage': percentage,
        'Upregulated_Count': up_count,
        'Downregulated_Count': down_count
    })

df_func_dist = pd.DataFrame(func_class_stats)
df_func_dist = df_func_dist.sort_values('Count', ascending=False)
df_func_dist.to_csv('functional_class_distribution.csv', index=False)
print("  ✓ functional_class_distribution.csv")

# ============================================================================
# STEP 5: VALIDATION & REPORT
# ============================================================================
print("\n[STEP 5] Validation & Quality Report")
print("="*80)

print("\n1. FIRST 5 ROWS OF CLEANED DATA:")
print(df[['Gene', 'Protein_Description', 'Fold_Change', 'log_2_fold_change',
          'Functional_Class']].head().to_string())

print("\n2. DATA QUALITY REPORT:")
print(f"\n   Total proteins: {len(df)}")
print(f"   Total columns: {len(df.columns)}")

print("\n   Missing values per column:")
missing_counts = df.isna().sum()
missing_counts = missing_counts[missing_counts > 0].sort_values(ascending=False)
if len(missing_counts) > 0:
    for col, count in missing_counts.items():
        pct = (count / len(df)) * 100
        print(f"     {col}: {count} ({pct:.1f}%)")
else:
    print("     No missing values!")

# Complete data check
complete_data = df.dropna()
print(f"\n   Proteins with complete data (no missing values): {len(complete_data)}")
print(f"   Proteins with some missing data: {len(df) - len(complete_data)}")

print("\n3. FUNCTIONAL CLASS SUMMARY:")
print(df_func_dist.to_string(index=False))

print("\n4. MEAN FOLD CHANGE PER FUNCTIONAL CLASS:")
fc_by_class = df.groupby('Functional_Class')['Fold_Change'].agg(['mean', 'median', 'count'])
fc_by_class = fc_by_class.sort_values('mean', ascending=False)
print(fc_by_class.to_string())

print("\n" + "="*80)
print("✓ DATA CLEANING COMPLETE!")
print("="*80)
import os
print(f"\nAll output files saved to: {os.getcwd()}")
print("\nFiles created:")
print("  1. cleaned_proteomics_data.csv")
print("  2. summary_statistics.csv")
print("  3. upregulated_proteins.csv")
print("  4. downregulated_proteins.csv")
print("  5. top_changers.csv")
print("  6. functional_class_distribution.csv")
