#!/usr/bin/env python3
"""
Generate PaintOmics 4 formatted input files from proteomics data
PaintOmics paints expression data onto KEGG pathway diagrams
"""

import pandas as pd
import numpy as np
import os

# Read the data
input_file = "../General Analysis/cleaned_proteomics_data_with_QC_flags.csv"
output_dir = "."

print("="*80)
print("PAINTOMICS 4 FILE GENERATOR")
print("="*80)
print("\nReading proteomics data...")
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
    Filter data according to PaintOmics requirements:
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

def add_regulation_column(df):
    """Add Regulation column: Up, Down, or Unchanged"""
    def classify_regulation(log2fc):
        if log2fc > 0:
            return "Up"
        elif log2fc < 0:
            return "Down"
        else:
            return "Unchanged"

    df['Regulation'] = df['log_2_fold_change'].apply(classify_regulation)
    return df

# ============================================================================
# PART 1: MAIN PAINTOMICS INPUT FILE (FULL FORMAT)
# ============================================================================

print("\n" + "="*80)
print("PART 1: Creating main PaintOmics input file (full format)")
print("="*80)

main_data = get_filtered_data(df, ['High', 'Medium'])
main_data = sort_by_abs_fc(main_data)
main_data = add_regulation_column(main_data)

print(f"Proteins after filtering: {len(main_data)}")
print(f"  - High confidence: {len(main_data[main_data['Confidence_Level'] == 'High'])}")
print(f"  - Medium confidence: {len(main_data[main_data['Confidence_Level'] == 'Medium'])}")
print(f"  - Upregulated: {len(main_data[main_data['Regulation'] == 'Up'])}")
print(f"  - Downregulated: {len(main_data[main_data['Regulation'] == 'Down'])}")

# Create PaintOmics format with all columns
paintomics_main = pd.DataFrame({
    'Name': main_data['Gene'],
    'Testosterone_vs_Vehicle': main_data['log_2_fold_change'],
    'Fold_Change': main_data['Fold_Change'],
    'Regulation': main_data['Regulation'],
    'Confidence': main_data['Confidence_Level'],
    'Function': main_data['Functional_Class'],
    'Description': main_data['Protein_Description']
})

# Save WITH header (required for PaintOmics)
output_file = os.path.join(output_dir, "paintomics_input_expression_data.csv")
paintomics_main.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  Format: CSV with header")
print(f"  Columns: Name, Testosterone_vs_Vehicle, Fold_Change, Regulation, Confidence, Function, Description")
print(f"  Proteins: {len(paintomics_main)}")
print(f"\nTop 5 most changed proteins:")
for idx, row in paintomics_main.head().iterrows():
    print(f"  {row['Name']:<15} {row['Testosterone_vs_Vehicle']:>6.2f}  {row['Regulation']:<5}  {row['Function']}")

# ============================================================================
# PART 2: MINIMAL FORMAT (TWO COLUMNS ONLY)
# ============================================================================

print("\n" + "="*80)
print("PART 2: Creating minimal input file (2 columns only)")
print("="*80)

paintomics_minimal = pd.DataFrame({
    'Name': main_data['Gene'],
    'Testosterone_vs_Vehicle': main_data['log_2_fold_change']
})

output_file = os.path.join(output_dir, "paintomics_minimal_input.csv")
paintomics_minimal.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  Format: CSV with header (minimal)")
print(f"  Columns: Name, Testosterone_vs_Vehicle")
print(f"  Proteins: {len(paintomics_minimal)}")
print(f"  Use this for: Quick analysis with basic coloring")

# ============================================================================
# PART 3: HIGH-CONFIDENCE ONLY
# ============================================================================

print("\n" + "="*80)
print("PART 3: Creating high-confidence only file")
print("="*80)

high_conf_data = get_filtered_data(df, ['High'])
high_conf_data = sort_by_abs_fc(high_conf_data)
high_conf_data = add_regulation_column(high_conf_data)

paintomics_high = pd.DataFrame({
    'Name': high_conf_data['Gene'],
    'Testosterone_vs_Vehicle': high_conf_data['log_2_fold_change'],
    'Fold_Change': high_conf_data['Fold_Change'],
    'Regulation': high_conf_data['Regulation'],
    'Confidence': high_conf_data['Confidence_Level'],
    'Function': high_conf_data['Functional_Class'],
    'Description': high_conf_data['Protein_Description']
})

output_file = os.path.join(output_dir, "paintomics_high_confidence.csv")
paintomics_high.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  Format: CSV with header (full format)")
print(f"  Proteins: {len(paintomics_high)} (High confidence only)")

# ============================================================================
# PART 4: FUNCTIONAL CLASS-SPECIFIC FILES
# ============================================================================

print("\n" + "="*80)
print("PART 4: Creating functional class-specific files")
print("="*80)

# 4.1 Cytokines
cytokines_data = main_data[main_data['Functional_Class'] == 'cytokine'].copy()
cytokines_data = sort_by_abs_fc(cytokines_data)

paintomics_cytokines = pd.DataFrame({
    'Name': cytokines_data['Gene'],
    'Testosterone_vs_Vehicle': cytokines_data['log_2_fold_change'],
    'Fold_Change': cytokines_data['Fold_Change'],
    'Regulation': cytokines_data['Regulation'],
    'Confidence': cytokines_data['Confidence_Level'],
    'Function': cytokines_data['Functional_Class'],
    'Description': cytokines_data['Protein_Description']
})

output_file = os.path.join(output_dir, "paintomics_cytokines.csv")
paintomics_cytokines.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  Cytokines: {len(paintomics_cytokines)}")
print(f"  Genes: {', '.join(paintomics_cytokines['Name'].tolist())}")
print(f"  Use for: Inflammatory pathway visualization (cytokine-cytokine receptor interaction)")

# 4.2 Growth factors
growth_factors_data = main_data[main_data['Functional_Class'] == 'growth factor'].copy()
growth_factors_data = sort_by_abs_fc(growth_factors_data)

paintomics_gf = pd.DataFrame({
    'Name': growth_factors_data['Gene'],
    'Testosterone_vs_Vehicle': growth_factors_data['log_2_fold_change'],
    'Fold_Change': growth_factors_data['Fold_Change'],
    'Regulation': growth_factors_data['Regulation'],
    'Confidence': growth_factors_data['Confidence_Level'],
    'Function': growth_factors_data['Functional_Class'],
    'Description': growth_factors_data['Protein_Description']
})

output_file = os.path.join(output_dir, "paintomics_growth_factors.csv")
paintomics_gf.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  Growth factors: {len(paintomics_gf)}")
print(f"  Genes: {', '.join(paintomics_gf['Name'].tolist())}")
print(f"  Use for: TGF-beta, PI3K-Akt, MAPK pathway visualization")

# 4.3 ECM proteins
# Include specific ECM-related proteins from "other" category
ecm_protein_names = ['BGN', 'VCAN', 'LUM', 'TGFBI', 'SERPINH1', 'SERPINE2', 'SERPINF1', 'PTX3', 'FRZB', 'PLXDC2', 'SLIT3']
ecm_data = main_data[main_data['Gene'].isin(ecm_protein_names)].copy()
ecm_data = sort_by_abs_fc(ecm_data)

paintomics_ecm = pd.DataFrame({
    'Name': ecm_data['Gene'],
    'Testosterone_vs_Vehicle': ecm_data['log_2_fold_change'],
    'Fold_Change': ecm_data['Fold_Change'],
    'Regulation': ecm_data['Regulation'],
    'Confidence': ecm_data['Confidence_Level'],
    'Function': ecm_data['Functional_Class'],
    'Description': ecm_data['Protein_Description']
})

output_file = os.path.join(output_dir, "paintomics_ecm_proteins.csv")
paintomics_ecm.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  ECM-related proteins: {len(paintomics_ecm)}")
print(f"  Genes: {', '.join(paintomics_ecm['Name'].tolist())}")
print(f"  Use for: ECM-receptor interaction, focal adhesion pathway visualization")

# ============================================================================
# PART 5: TIME-COURSE TEMPLATE (FOR FUTURE USE)
# ============================================================================

print("\n" + "="*80)
print("PART 5: Creating time-course template")
print("="*80)

# For now, use Vehicle as baseline (0) and Testosterone as treatment (log2FC)
paintomics_timecourse = pd.DataFrame({
    'Name': main_data['Gene'],
    'Vehicle': 0.0,  # Baseline
    'Testosterone': main_data['log_2_fold_change']
})

output_file = os.path.join(output_dir, "paintomics_timecourse_template.csv")
paintomics_timecourse.to_csv(output_file, index=False)
print(f"\n✓ Created: {output_file}")
print(f"  Format: CSV with header")
print(f"  Columns: Name, Vehicle, Testosterone")
print(f"  Note: Vehicle = 0.0 (baseline), Testosterone = log2FC")
print(f"  Use this as: Template for future time-course experiments")
print(f"  Future format example: Name, T0h, T6h, T24h, T48h")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

print("\n" + "="*80)
print("SUMMARY OF GENERATED FILES")
print("="*80)

summary_data = {
    'File': [
        'paintomics_input_expression_data.csv',
        'paintomics_minimal_input.csv',
        'paintomics_high_confidence.csv',
        'paintomics_cytokines.csv',
        'paintomics_growth_factors.csv',
        'paintomics_ecm_proteins.csv',
        'paintomics_timecourse_template.csv'
    ],
    'Proteins': [
        len(paintomics_main),
        len(paintomics_minimal),
        len(paintomics_high),
        len(paintomics_cytokines),
        len(paintomics_gf),
        len(paintomics_ecm),
        len(paintomics_timecourse)
    ],
    'Format': [
        'Full (7 cols)',
        'Minimal (2 cols)',
        'Full (7 cols)',
        'Full (7 cols)',
        'Full (7 cols)',
        'Full (7 cols)',
        'Timecourse (3 cols)'
    ],
    'Description': [
        'Main file - RECOMMENDED',
        'Quick analysis',
        'High confidence only',
        'Inflammatory pathways',
        'Growth signaling pathways',
        'ECM/structural pathways',
        'Template for future use'
    ]
}

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

# ============================================================================
# DETAILED PATHWAY PREDICTIONS
# ============================================================================

print("\n" + "="*80)
print("EXPECTED KEGG PATHWAYS TO EXPLORE")
print("="*80)

print("\n1. CYTOKINE-CYTOKINE RECEPTOR INTERACTION (mmu04060)")
print("   Proteins likely mapped: CCL7, CCL2, CXCL5")
print("   Regulation: All downregulated")
print("   Implication: Reduced inflammatory signaling")

print("\n2. TGF-BETA SIGNALING PATHWAY (mmu04350)")
print("   Proteins likely mapped: TGFB1, TGFB3, GDF15, AMH")
print("   Regulation: Mixed (TGFB1/3/GDF15 down, AMH up)")
print("   Implication: Complex TGF-beta pathway modulation")

print("\n3. PI3K-AKT SIGNALING PATHWAY (mmu04151)")
print("   Proteins likely mapped: AGT, GAS6, multiple growth factors")
print("   Regulation: Mixed")
print("   Implication: Cell survival and proliferation changes")

print("\n4. ECM-RECEPTOR INTERACTION (mmu04512)")
print("   Proteins likely mapped: VCAN, LUM, BGN, TGFBI")
print("   Regulation: Mostly upregulated")
print("   Implication: ECM remodeling and tissue structure changes")

print("\n5. COMPLEMENT AND COAGULATION CASCADES (mmu04610)")
print("   Proteins likely mapped: C4A/C4B, SERPIN family")
print("   Regulation: C4A/C4B upregulated, SERPINs mixed")
print("   Implication: Immune system modulation")

print("\n6. FOCAL ADHESION (mmu04510)")
print("   Proteins likely mapped: VCAN, BGN, multiple ECM proteins")
print("   Regulation: Mostly upregulated")
print("   Implication: Cell-matrix interaction changes")

print("\n" + "="*80)
print("FILE GENERATION COMPLETE!")
print("="*80)
print("\nAll files are CSV format with HEADERS (required by PaintOmics)")
print("Next step: Review PAINTOMICS_UPLOAD_INSTRUCTIONS.md for detailed guidance")
print("\nPaintOmics 4 URL: http://www.paintomics.org/")
