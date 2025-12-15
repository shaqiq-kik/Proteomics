"""
STEP 1.2: Data Quality Investigation & Resolution
Comprehensive analysis of missing data, edge cases, and validation
"""

import pandas as pd
import numpy as np
import re

print("="*80)
print("STEP 1.2: DATA QUALITY INVESTIGATION")
print("="*80)

# Load cleaned data
df = pd.read_csv('cleaned_proteomics_data.csv')
print(f"\nLoaded {len(df)} proteins from cleaned_proteomics_data.csv")

# ============================================================================
# TASK 1: MISSING REPLICATE DATA ANALYSIS
# ============================================================================
print("\n" + "="*80)
print("TASK 1: MISSING REPLICATE DATA ANALYSIS")
print("="*80)

replicate_cols = {
    'Vehicle_Rep1': 'Vehicle_Rep1_31579',
    'Vehicle_Rep2': 'Vehicle_Rep2_31581',
    'Testosterone_Rep1': 'Testosterone_Rep1_31578',
    'Testosterone_Rep2': 'Testosterone_Rep2_31580'
}

# Create missing data report
missing_data_report = []

for idx, row in df.iterrows():
    gene = row['Gene']

    rep_status = {}
    missing_count = 0

    for short_name, col_name in replicate_cols.items():
        is_present = pd.notna(row[col_name]) and row[col_name] != 0
        rep_status[short_name] = 'Present' if is_present else 'Missing'
        if not is_present:
            missing_count += 1

    # Determine completeness status
    if missing_count == 0:
        status = 'Complete (4/4)'
    elif missing_count == 1:
        status = 'Missing 1 (3/4)'
    elif missing_count == 2:
        status = 'Missing 2 (2/4)'
    elif missing_count == 3:
        status = 'Missing 3 (1/4)'
    else:
        status = 'Missing All (0/4)'

    missing_data_report.append({
        'Gene': gene,
        'Vehicle_Rep1': rep_status['Vehicle_Rep1'],
        'Vehicle_Rep2': rep_status['Vehicle_Rep2'],
        'Testosterone_Rep1': rep_status['Testosterone_Rep1'],
        'Testosterone_Rep2': rep_status['Testosterone_Rep2'],
        'Data_Completeness_Status': status,
        'Missing_Count': missing_count
    })

df_missing = pd.DataFrame(missing_data_report)

# Count completeness categories
complete_data = len(df_missing[df_missing['Missing_Count'] == 0])
missing_1 = len(df_missing[df_missing['Missing_Count'] == 1])
missing_2_plus = len(df_missing[df_missing['Missing_Count'] >= 2])

print(f"\nREPLICATE DATA COMPLETENESS SUMMARY:")
print(f"  Complete data (all 4 replicates): {complete_data} proteins ({complete_data/len(df)*100:.1f}%)")
print(f"  Missing 1 replicate: {missing_1} proteins ({missing_1/len(df)*100:.1f}%)")
print(f"  Missing 2+ replicates: {missing_2_plus} proteins ({missing_2_plus/len(df)*100:.1f}%)")

print(f"\nPROTEINS WITH INCOMPLETE DATA:")
incomplete = df_missing[df_missing['Missing_Count'] > 0].sort_values('Missing_Count', ascending=False)
if len(incomplete) > 0:
    print(incomplete[['Gene', 'Vehicle_Rep1', 'Vehicle_Rep2', 'Testosterone_Rep1',
                      'Testosterone_Rep2', 'Data_Completeness_Status']].to_string(index=False))
else:
    print("  All proteins have complete replicate data!")

# ============================================================================
# TASK 2: ZERO FOLD CHANGE INVESTIGATION
# ============================================================================
print("\n" + "="*80)
print("TASK 2: ZERO FOLD CHANGE INVESTIGATION")
print("="*80)

zero_fc_proteins = df[df['Fold_Change'] == 0.0].copy()
print(f"\nFound {len(zero_fc_proteins)} proteins with Fold_Change = 0")

if len(zero_fc_proteins) > 0:
    print("\nDETAILED ANALYSIS:")
    for idx, row in zero_fc_proteins.iterrows():
        print(f"\n  Protein: {row['Gene']} ({row['Protein_Description']})")
        print(f"    Vehicle:")
        print(f"      Rep1: {row['Vehicle_Rep1_31579']}")
        print(f"      Rep2: {row['Vehicle_Rep2_31581']}")
        print(f"      Mean: {row['Vehicle_Mean']}")
        print(f"    Testosterone:")
        print(f"      Rep1: {row['Testosterone_Rep1_31578']}")
        print(f"      Rep2: {row['Testosterone_Rep2_31580']}")
        print(f"      Mean: {row['Testosterone_Mean']}")
        print(f"    Fold Change: {row['Fold_Change']}")
        print(f"    Log2 FC: {row['log_2_fold_change']}")

        # Determine cause
        v_mean = row['Vehicle_Mean']
        t_mean = row['Testosterone_Mean']

        if pd.isna(t_mean) or t_mean == 0:
            if pd.notna(v_mean) and v_mean > 0:
                print(f"    DIAGNOSIS: Complete suppression (Testosterone = 0, Vehicle = {v_mean:.2e})")
                print(f"    STATUS: Biologically valid - protein completely suppressed by testosterone")
            elif pd.isna(v_mean) or v_mean == 0:
                print(f"    DIAGNOSIS: No detection in either condition")
                print(f"    STATUS: Protein not detected - may need to exclude from analysis")

        print(f"    ---")

# ============================================================================
# TASK 3: UNIPROT ID VALIDATION
# ============================================================================
print("\n" + "="*80)
print("TASK 3: UNIPROT ID VALIDATION")
print("="*80)

# Standard UniProt pattern: [PQ][0-9]{5}
# Extended pattern: [OPQ][0-9][A-Z0-9]{3}[0-9]
standard_pattern = r'^[PQ]\d{5}$'
extended_pattern = r'^[OPQ]\d[A-Z0-9]{3}\d$'

df['UniProt_Standard_Format'] = df['UniProt_Accession'].apply(
    lambda x: bool(re.match(standard_pattern, str(x))) if pd.notna(x) else False
)

df['UniProt_Extended_Format'] = df['UniProt_Accession'].apply(
    lambda x: bool(re.match(extended_pattern, str(x))) if pd.notna(x) else False
)

df['UniProt_Valid'] = df['UniProt_Standard_Format'] | df['UniProt_Extended_Format']

standard_count = df['UniProt_Standard_Format'].sum()
extended_count = df['UniProt_Extended_Format'].sum() - standard_count
valid_count = df['UniProt_Valid'].sum()
invalid_count = len(df) - valid_count

print(f"\nUNIPROT ID VALIDATION RESULTS:")
print(f"  Standard format ([PQ]#####): {standard_count} proteins")
print(f"  Extended format ([OPQ]#XXX#): {extended_count} proteins")
print(f"  Total valid: {valid_count} proteins ({valid_count/len(df)*100:.1f}%)")
print(f"  Invalid/Missing: {invalid_count} proteins ({invalid_count/len(df)*100:.1f}%)")

# Check non-standard IDs
non_standard = df[~df['UniProt_Valid']][['Gene', 'UniProt_Accession']]
if len(non_standard) > 0:
    print(f"\nNON-STANDARD UNIPROT IDs:")
    print(non_standard.to_string(index=False))

# Special check for AOC1 (Q8JZQ5)
aoc1 = df[df['Gene'] == 'AOC1']
if len(aoc1) > 0:
    aoc1_id = aoc1.iloc[0]['UniProt_Accession']
    print(f"\nSPECIAL CASE - AOC1:")
    print(f"  UniProt ID: {aoc1_id}")
    print(f"  Format: Extended/Non-standard")
    print(f"  Analysis: Q8JZQ5 is a valid extended UniProt format")
    print(f"  Recommendation: Keep as-is (valid ID)")

# ============================================================================
# TASK 4: OVERALL DATA QUALITY ASSESSMENT
# ============================================================================
print("\n" + "="*80)
print("TASK 4: OVERALL DATA QUALITY ASSESSMENT")
print("="*80)

# Missing data statistics per column
missing_stats = []
for col in df.columns:
    missing_count = df[col].isna().sum()
    missing_pct = (missing_count / len(df)) * 100
    if missing_count > 0:
        missing_stats.append({
            'Column': col,
            'Missing_Count': missing_count,
            'Missing_Percent': missing_pct,
            'Present_Count': len(df) - missing_count
        })

df_missing_stats = pd.DataFrame(missing_stats).sort_values('Missing_Percent', ascending=False)

print(f"\nMISSING DATA BY COLUMN:")
if len(df_missing_stats) > 0:
    print(df_missing_stats.to_string(index=False))
else:
    print("  No missing data!")

# Total data completeness
total_cells = df.shape[0] * df.shape[1]
missing_cells = df.isna().sum().sum()
completeness = ((total_cells - missing_cells) / total_cells) * 100

print(f"\nOVERALL DATA COMPLETENESS: {completeness:.2f}%")
print(f"  Total cells: {total_cells}")
print(f"  Missing cells: {missing_cells}")
print(f"  Complete cells: {total_cells - missing_cells}")

# ============================================================================
# TASK 5: CREATE QUALITY FLAGS
# ============================================================================
print("\n" + "="*80)
print("TASK 5: CREATING QUALITY FLAGS")
print("="*80)

# Flag 1: Replicate_Completeness
def determine_replicate_completeness(row):
    reps = [row['Vehicle_Rep1_31579'], row['Vehicle_Rep2_31581'],
            row['Testosterone_Rep1_31578'], row['Testosterone_Rep2_31580']]
    missing = sum([1 for r in reps if pd.isna(r) or r == 0])

    if missing == 0:
        return "Complete"
    elif missing == 1:
        return "Missing_1"
    else:
        return "Missing_2+"

df['Replicate_Completeness'] = df.apply(determine_replicate_completeness, axis=1)

# Flag 2: FC_Type
def determine_fc_type(row):
    v_mean = row['Vehicle_Mean']
    t_mean = row['Testosterone_Mean']
    fc = row['Fold_Change']

    # Check for missing data
    if pd.isna(fc) or pd.isna(v_mean) or pd.isna(t_mean):
        return "Cannot_Calculate"

    # Check for both zero
    if v_mean == 0 and t_mean == 0:
        return "No_Detection"

    # Complete suppression (testosterone = 0, vehicle > 0)
    if t_mean == 0 and v_mean > 0:
        return "Complete_Suppression"

    # Complete induction (vehicle = 0, testosterone > 0)
    if v_mean == 0 and t_mean > 0:
        return "Complete_Induction"

    # Normal calculation
    return "Normal"

df['FC_Type'] = df.apply(determine_fc_type, axis=1)

# Flag 3: Confidence_Level
def determine_confidence(row):
    completeness = row['Replicate_Completeness']
    fc_type = row['FC_Type']

    # Low confidence conditions
    if completeness == "Missing_2+":
        return "Low"
    if fc_type in ["Cannot_Calculate", "No_Detection"]:
        return "Low"

    # High confidence
    if completeness == "Complete" and fc_type == "Normal":
        return "High"

    # Medium confidence (everything else)
    return "Medium"

df['Confidence_Level'] = df.apply(determine_confidence, axis=1)

# Display flag distribution
print("\nQUALITY FLAG DISTRIBUTION:")
print("\n1. Replicate_Completeness:")
print(df['Replicate_Completeness'].value_counts().to_string())

print("\n2. FC_Type:")
print(df['FC_Type'].value_counts().to_string())

print("\n3. Confidence_Level:")
print(df['Confidence_Level'].value_counts().to_string())

# ============================================================================
# EXPORT RESULTS
# ============================================================================
print("\n" + "="*80)
print("EXPORTING RESULTS")
print("="*80)

# 1. Main dataset with QC flags
df.to_csv('cleaned_proteomics_data_with_QC_flags.csv', index=False)
print("✓ cleaned_proteomics_data_with_QC_flags.csv")

# 2. Missing data summary
df_missing.to_csv('missing_data_summary.csv', index=False)
print("✓ missing_data_summary.csv")

# 3. Proteins flagged for review
flagged_proteins = []

for idx, row in df.iterrows():
    issues = []
    severity = None
    recommendation = None

    # Check for issues
    if row['Confidence_Level'] == 'Low':
        if row['Replicate_Completeness'] == 'Missing_2+':
            issues.append('Missing 2+ replicates')
            severity = 'High'
            recommendation = 'Exclude'
        if row['FC_Type'] == 'No_Detection':
            issues.append('Not detected in either condition')
            severity = 'High'
            recommendation = 'Exclude'
        if row['FC_Type'] == 'Cannot_Calculate':
            issues.append('Cannot calculate fold change')
            severity = 'High'
            recommendation = 'Exclude'

    if row['Confidence_Level'] == 'Medium':
        if row['Replicate_Completeness'] == 'Missing_1':
            issues.append('Missing 1 replicate')
            severity = 'Medium'
            recommendation = 'Flag'
        if row['FC_Type'] in ['Complete_Suppression', 'Complete_Induction']:
            issues.append(f'{row["FC_Type"]} - edge case')
            severity = 'Medium'
            recommendation = 'Flag'

    if not row['UniProt_Valid']:
        issues.append('Non-standard UniProt ID')
        if severity is None:
            severity = 'Low'
        if recommendation is None:
            recommendation = 'Flag'

    if issues:
        flagged_proteins.append({
            'Gene': row['Gene'],
            'Protein_Description': row['Protein_Description'],
            'Issues': '; '.join(issues),
            'Severity': severity,
            'Recommendation': recommendation,
            'Confidence_Level': row['Confidence_Level'],
            'FC_Type': row['FC_Type'],
            'Replicate_Completeness': row['Replicate_Completeness']
        })

df_flagged = pd.DataFrame(flagged_proteins)
df_flagged.to_csv('proteins_flagged_for_review.csv', index=False)
print("✓ proteins_flagged_for_review.csv")

print(f"\n  Total proteins flagged: {len(df_flagged)}")
if len(df_flagged) > 0:
    print(f"  Severity breakdown:")
    print(f"    High: {len(df_flagged[df_flagged['Severity'] == 'High'])}")
    print(f"    Medium: {len(df_flagged[df_flagged['Severity'] == 'Medium'])}")
    print(f"    Low: {len(df_flagged[df_flagged['Severity'] == 'Low'])}")

# 4. High confidence proteins
high_conf = df[df['Confidence_Level'] == 'High'].copy()
high_conf.to_csv('high_confidence_proteins.csv', index=False)
print(f"✓ high_confidence_proteins.csv ({len(high_conf)} proteins)")

# 5. Data quality report (text file)
with open('data_quality_report.txt', 'w') as f:
    f.write("="*80 + "\n")
    f.write("DATA QUALITY REPORT - PROTEOMICS ANALYSIS\n")
    f.write("="*80 + "\n\n")

    f.write("SUMMARY\n")
    f.write("-" * 40 + "\n")
    f.write(f"Total proteins analyzed: {len(df)}\n")
    f.write(f"Overall data completeness: {completeness:.2f}%\n")
    f.write(f"High confidence proteins: {len(high_conf)} ({len(high_conf)/len(df)*100:.1f}%)\n")
    f.write(f"Flagged for review: {len(df_flagged)} ({len(df_flagged)/len(df)*100:.1f}%)\n\n")

    f.write("TASK 1: REPLICATE DATA COMPLETENESS\n")
    f.write("-" * 40 + "\n")
    f.write(f"Complete data (4/4): {complete_data} proteins ({complete_data/len(df)*100:.1f}%)\n")
    f.write(f"Missing 1 replicate: {missing_1} proteins ({missing_1/len(df)*100:.1f}%)\n")
    f.write(f"Missing 2+ replicates: {missing_2_plus} proteins ({missing_2_plus/len(df)*100:.1f}%)\n\n")

    if len(incomplete) > 0:
        f.write("Proteins with incomplete data:\n")
        f.write(incomplete[['Gene', 'Data_Completeness_Status']].to_string(index=False))
        f.write("\n\n")

    f.write("TASK 2: ZERO FOLD CHANGE PROTEINS\n")
    f.write("-" * 40 + "\n")
    f.write(f"Proteins with FC = 0: {len(zero_fc_proteins)}\n")
    if len(zero_fc_proteins) > 0:
        for idx, row in zero_fc_proteins.iterrows():
            # Get FC_Type from main dataframe
            fc_type = df.loc[idx, 'FC_Type']
            f.write(f"\n  {row['Gene']}:\n")
            f.write(f"    Vehicle Mean: {row['Vehicle_Mean']}\n")
            f.write(f"    Testosterone Mean: {row['Testosterone_Mean']}\n")
            f.write(f"    FC Type: {fc_type}\n")
    f.write("\n")

    f.write("TASK 3: UNIPROT ID VALIDATION\n")
    f.write("-" * 40 + "\n")
    f.write(f"Valid UniProt IDs: {valid_count}/{len(df)} ({valid_count/len(df)*100:.1f}%)\n")
    f.write(f"  Standard format: {standard_count}\n")
    f.write(f"  Extended format: {extended_count}\n")
    f.write(f"Invalid/Missing: {invalid_count}\n\n")

    if len(non_standard) > 0:
        f.write("Non-standard IDs:\n")
        f.write(non_standard.to_string(index=False))
        f.write("\n\n")

    f.write("TASK 4: MISSING DATA STATISTICS\n")
    f.write("-" * 40 + "\n")
    if len(df_missing_stats) > 0:
        f.write(df_missing_stats.to_string(index=False))
    else:
        f.write("No missing data!\n")
    f.write("\n\n")

    f.write("TASK 5: QUALITY FLAGS\n")
    f.write("-" * 40 + "\n")
    f.write("Replicate Completeness:\n")
    f.write(df['Replicate_Completeness'].value_counts().to_string())
    f.write("\n\nFC Type:\n")
    f.write(df['FC_Type'].value_counts().to_string())
    f.write("\n\nConfidence Level:\n")
    f.write(df['Confidence_Level'].value_counts().to_string())
    f.write("\n\n")

    f.write("RECOMMENDATIONS\n")
    f.write("-" * 40 + "\n")
    f.write(f"1. USE FOR ALL ANALYSES: {len(high_conf)} high-confidence proteins\n")
    f.write(f"   - Complete replicate data\n")
    f.write(f"   - Normal fold change calculations\n")
    f.write(f"   - No data quality concerns\n\n")

    medium_conf = df[df['Confidence_Level'] == 'Medium']
    f.write(f"2. USE WITH CAUTION: {len(medium_conf)} medium-confidence proteins\n")
    f.write(f"   - May have 1 missing replicate OR edge case FC\n")
    f.write(f"   - Include but flag in results\n\n")

    low_conf = df[df['Confidence_Level'] == 'Low']
    f.write(f"3. EXCLUDE FROM ANALYSIS: {len(low_conf)} low-confidence proteins\n")
    f.write(f"   - Missing 2+ replicates OR cannot calculate FC\n")
    f.write(f"   - Not reliable for statistical analysis\n\n")

    f.write("="*80 + "\n")
    f.write("END OF REPORT\n")
    f.write("="*80 + "\n")

print("✓ data_quality_report.txt")

print("\n" + "="*80)
print("DATA QUALITY INVESTIGATION COMPLETE!")
print("="*80)
print(f"\nFiles created:")
print("  1. cleaned_proteomics_data_with_QC_flags.csv")
print("  2. data_quality_report.txt")
print("  3. missing_data_summary.csv")
print("  4. proteins_flagged_for_review.csv")
print("  5. high_confidence_proteins.csv")

print(f"\nKEY FINDINGS:")
print(f"  High confidence proteins: {len(high_conf)}/{len(df)} ({len(high_conf)/len(df)*100:.1f}%)")
print(f"  Medium confidence: {len(medium_conf)}/{len(df)} ({len(medium_conf)/len(df)*100:.1f}%)")
print(f"  Low confidence: {len(low_conf)}/{len(df)} ({len(low_conf)/len(df)*100:.1f}%)")
print(f"  Overall completeness: {completeness:.2f}%")
