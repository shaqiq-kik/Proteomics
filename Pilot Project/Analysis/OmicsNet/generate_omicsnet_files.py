#!/usr/bin/env python3
"""
Generate OmicsNet 2.0 formatted input files from proteomics data
"""

import pandas as pd
import numpy as np
import os

# Read the data
input_file = "../General Analysis/cleaned_proteomics_data_with_QC_flags.csv"
output_dir = "."

print("Reading proteomics data...")
df = pd.read_csv(input_file)

print(f"Total proteins in dataset: {len(df)}")
print(f"\nFC_Type distribution:")
print(df['FC_Type'].value_counts())
print(f"\nConfidence_Level distribution:")
print(df['Confidence_Level'].value_counts())

# ============================================================================
# FILTERING FUNCTIONS
# ============================================================================

def get_filtered_data(df, confidence_levels=['High', 'Medium']):
    """
    Filter data according to OmicsNet requirements:
    - EXCLUDE FC_Type == 'Cannot_Calculate'
    - EXCLUDE NaN Fold_Change
    - INCLUDE 'Complete_Suppression' (they have valid log2FC)
    - Filter by confidence level
    """
    filtered = df[
        (df['FC_Type'] != 'Cannot_Calculate') &  # Exclude Cannot_Calculate
        (~df['Fold_Change'].isna()) &  # Exclude NaN fold changes
        (df['Confidence_Level'].isin(confidence_levels))  # Confidence filter
    ].copy()

    return filtered

def sort_by_abs_fc(df):
    """Sort by absolute log2 fold change (most changed first)"""
    df['abs_log2fc'] = df['log_2_fold_change'].abs()
    df_sorted = df.sort_values('abs_log2fc', ascending=False)
    df_sorted = df_sorted.drop('abs_log2fc', axis=1)
    return df_sorted

# ============================================================================
# PART 1: MAIN OMICSNET INPUT FILE (WITH EXPRESSION)
# ============================================================================

print("\n" + "="*80)
print("PART 1: Creating main OmicsNet input file with expression data")
print("="*80)

main_data = get_filtered_data(df, ['High', 'Medium'])
main_data = sort_by_abs_fc(main_data)

print(f"Proteins after filtering: {len(main_data)}")
print(f"  - High confidence: {len(main_data[main_data['Confidence_Level'] == 'High'])}")
print(f"  - Medium confidence: {len(main_data[main_data['Confidence_Level'] == 'Medium'])}")

# Create 3-column format: Gene, log2FC, "Protein"
omicsnet_main = main_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_main['Type'] = 'Protein'

# Save WITHOUT header
output_file = os.path.join(output_dir, "omicsnet_input_with_expression.txt")
omicsnet_main.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Format: Tab-delimited, 3 columns (Gene, log2FC, 'Protein'), NO header")
print(f"  Proteins: {len(omicsnet_main)}")
print(f"\nTop 5 most changed proteins:")
for idx, row in omicsnet_main.head().iterrows():
    print(f"  {row['Gene']:<15} {row['log_2_fold_change']:>6.2f}")

# ============================================================================
# PART 2: SIMPLE GENE LIST (NO EXPRESSION)
# ============================================================================

print("\n" + "="*80)
print("PART 2: Creating simple gene list (no expression)")
print("="*80)

gene_list = main_data['Gene'].tolist()

output_file = os.path.join(output_dir, "omicsnet_genelist_only.txt")
with open(output_file, 'w') as f:
    for gene in gene_list:
        f.write(f"{gene}\n")

print(f"\n✓ Created: {output_file}")
print(f"  Format: One gene per line, NO header")
print(f"  Genes: {len(gene_list)}")

# ============================================================================
# PART 3: HIGH-CONFIDENCE ONLY
# ============================================================================

print("\n" + "="*80)
print("PART 3: Creating high-confidence only file")
print("="*80)

high_conf_data = get_filtered_data(df, ['High'])
high_conf_data = sort_by_abs_fc(high_conf_data)

omicsnet_high = high_conf_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_high['Type'] = 'Protein'

output_file = os.path.join(output_dir, "omicsnet_high_confidence_only.txt")
omicsnet_high.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Format: Tab-delimited, 3 columns (Gene, log2FC, 'Protein'), NO header")
print(f"  Proteins: {len(omicsnet_high)} (High confidence only)")

# ============================================================================
# PART 4: FUNCTIONAL CLASS-SPECIFIC NETWORKS
# ============================================================================

print("\n" + "="*80)
print("PART 4: Creating functional class-specific files")
print("="*80)

# Cytokines only
cytokines_data = main_data[main_data['Functional_Class'] == 'cytokine'].copy()
cytokines_data = sort_by_abs_fc(cytokines_data)
omicsnet_cytokines = cytokines_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_cytokines['Type'] = 'Protein'

output_file = os.path.join(output_dir, "omicsnet_cytokines.txt")
omicsnet_cytokines.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Cytokines: {len(omicsnet_cytokines)}")
print(f"  Genes: {', '.join(omicsnet_cytokines['Gene'].tolist())}")

# Growth factors only
growth_factors_data = main_data[main_data['Functional_Class'] == 'growth factor'].copy()
growth_factors_data = sort_by_abs_fc(growth_factors_data)
omicsnet_gf = growth_factors_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_gf['Type'] = 'Protein'

output_file = os.path.join(output_dir, "omicsnet_growth_factors.txt")
omicsnet_gf.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Growth factors: {len(omicsnet_gf)}")
print(f"  Genes: {', '.join(omicsnet_gf['Gene'].tolist())}")

# Signaling proteins (cytokines + growth factors)
signaling_data = main_data[
    main_data['Functional_Class'].isin(['cytokine', 'growth factor'])
].copy()
signaling_data = sort_by_abs_fc(signaling_data)
omicsnet_signaling = signaling_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_signaling['Type'] = 'Protein'

output_file = os.path.join(output_dir, "omicsnet_signaling_proteins.txt")
omicsnet_signaling.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Signaling proteins (cytokines + growth factors): {len(omicsnet_signaling)}")

# ============================================================================
# PART 5: UPREGULATED AND DOWNREGULATED FILES
# ============================================================================

print("\n" + "="*80)
print("PART 5: Creating upregulated and downregulated files")
print("="*80)

# Upregulated (log2FC > 0)
upregulated_data = main_data[main_data['log_2_fold_change'] > 0].copy()
upregulated_data = sort_by_abs_fc(upregulated_data)
omicsnet_up = upregulated_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_up['Type'] = 'Protein'

output_file = os.path.join(output_dir, "omicsnet_upregulated.txt")
omicsnet_up.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Upregulated proteins: {len(omicsnet_up)}")
print(f"  Log2FC range: {omicsnet_up['log_2_fold_change'].min():.2f} to {omicsnet_up['log_2_fold_change'].max():.2f}")

# Downregulated (log2FC < 0)
downregulated_data = main_data[main_data['log_2_fold_change'] < 0].copy()
downregulated_data = sort_by_abs_fc(downregulated_data)
omicsnet_down = downregulated_data[['Gene', 'log_2_fold_change']].copy()
omicsnet_down['Type'] = 'Protein'

output_file = os.path.join(output_dir, "omicsnet_downregulated.txt")
omicsnet_down.to_csv(output_file, sep='\t', index=False, header=False)
print(f"\n✓ Created: {output_file}")
print(f"  Downregulated proteins: {len(omicsnet_down)}")
print(f"  Log2FC range: {omicsnet_down['log_2_fold_change'].min():.2f} to {omicsnet_down['log_2_fold_change'].max():.2f}")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*80)
print("SUMMARY OF GENERATED FILES")
print("="*80)

summary_data = {
    'File': [
        'omicsnet_input_with_expression.txt',
        'omicsnet_high_confidence_only.txt',
        'omicsnet_genelist_only.txt',
        'omicsnet_cytokines.txt',
        'omicsnet_growth_factors.txt',
        'omicsnet_signaling_proteins.txt',
        'omicsnet_upregulated.txt',
        'omicsnet_downregulated.txt'
    ],
    'Proteins': [
        len(omicsnet_main),
        len(omicsnet_high),
        len(gene_list),
        len(omicsnet_cytokines),
        len(omicsnet_gf),
        len(omicsnet_signaling),
        len(omicsnet_up),
        len(omicsnet_down)
    ],
    'Description': [
        'Main file (High + Medium confidence)',
        'High confidence only',
        'Simple gene list (no expression)',
        'Cytokines only',
        'Growth factors only',
        'Cytokines + Growth factors',
        'Upregulated proteins',
        'Downregulated proteins'
    ]
}

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

print("\n" + "="*80)
print("FILE GENERATION COMPLETE!")
print("="*80)
print("\nNext step: Review OMICSNET_UPLOAD_INSTRUCTIONS.md for upload guidance")
