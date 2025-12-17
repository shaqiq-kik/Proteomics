#!/usr/bin/env python3
"""
Create TSV (tab-separated) versions for PaintOmics
Some bioinformatics tools prefer TSV even when labeled as CSV
"""

import pandas as pd
import os
import glob

print("="*80)
print("CREATING TSV (TAB-SEPARATED) VERSIONS FOR PAINTOMICS")
print("="*80)
print("\nPaintOmics may expect TAB-separated files instead of comma-separated")
print("Creating both CSV and TSV versions...\n")

# Find all CSV files
csv_files = glob.glob("paintomics_*.csv")

for input_file in csv_files:
    print(f"Processing: {input_file}")

    # Read the CSV file
    df = pd.read_csv(input_file)

    # Verify format
    if list(df.columns) != ['ID', 'Vehicle', 'Testosterone']:
        print(f"  WARNING: Unexpected columns in {input_file}")
        print(f"  Columns: {list(df.columns)}")
        continue

    print(f"  Rows: {len(df)}")
    print(f"  Columns: {list(df.columns)}")

    # Create TSV version
    tsv_file = input_file.replace('.csv', '.txt')

    # Save as TSV (tab-separated)
    df.to_csv(tsv_file, sep='\t', index=False, encoding='utf-8', lineterminator='\n')
    print(f"  ✓ Created TSV: {tsv_file}")

    # Also recreate CSV with explicit settings to ensure no BOM
    df.to_csv(input_file, sep=',', index=False, encoding='utf-8', lineterminator='\n')
    print(f"  ✓ Recreated CSV: {input_file}")

    print()

print("="*80)
print("CONVERSION COMPLETE")
print("="*80)

print("\nFiles created:")
print("\nCSV files (comma-separated):")
for f in sorted(glob.glob("paintomics_*.csv")):
    print(f"  - {f}")

print("\nTXT files (tab-separated):")
for f in sorted(glob.glob("paintomics_*.txt")):
    print(f"  - {f}")

print("\n" + "="*80)
print("UPLOAD INSTRUCTIONS")
print("="*80)
print("\nTry uploading in this order:")
print("1. First try: paintomics_input_expression_data.txt (TAB-separated)")
print("2. If that fails: paintomics_input_expression_data.csv (COMMA-separated)")
print("\nSettings:")
print("  - Organism: Mouse (Mus musculus)")
print("  - Has header: YES")
print("  - ID column: ID")
print("  - Conditions: Select BOTH 'Vehicle' and 'Testosterone'")
print("\n" + "="*80)

# Show example of both formats
print("\nEXAMPLE FORMATS:")
print("\nCSV format (comma-separated):")
csv_example = pd.read_csv("paintomics_input_expression_data.csv")
print(csv_example.head(3).to_csv(index=False))

print("TXT/TSV format (tab-separated):")
tsv_example = pd.read_csv("paintomics_input_expression_data.txt", sep='\t')
print(tsv_example.head(3).to_csv(sep='\t', index=False))
