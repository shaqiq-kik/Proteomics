# PaintOmics 4 - Quick Start Guide

## ✅ FILES ARE READY

All files have been converted to PaintOmics-compatible format.

**Status**: ✅ READY FOR UPLOAD
**Error fixed**: "Expected at least 2 columns" - RESOLVED ✅

---

## Upload in 5 Steps

### 1. Go to PaintOmics
**URL**: http://www.paintomics.org/
Click: **"Start a new job"**

### 2. Job Settings
- **Job name**: `Testosterone_Proteomics`
- **Organism**: **Mouse (Mus musculus)** ⚠️ CRITICAL
- **Description**: "SILAC proteomics: testosterone vs. vehicle"

### 3. Upload File
**RECOMMENDED FILE**: `paintomics_input_expression_data_FIXED.csv`

- Click "Upload data"
- Select file: `paintomics_input_expression_data_FIXED.csv`
- Data type: Proteomics
- Format: CSV

### 4. Column Settings
⚠️ **CRITICAL SETTINGS**:
- **Has header**: ✅ YES (check this!)
- **ID column**: Select **"ID"**
- **Condition columns**: Select **BOTH "Vehicle" AND "Testosterone"**
  - Must select BOTH columns for PaintOmics to work
  - Vehicle = baseline (0)
  - Testosterone = log2FC

### 5. Gene ID Settings
- **ID type**: Gene Symbol
- Click **"Submit"**
- Wait 2-10 minutes

---

## File Selection Guide

| File | Proteins | Use When |
|------|----------|----------|
| **paintomics_input_expression_data_FIXED.csv** ⭐ | 33 | MAIN FILE - Start here |
| paintomics_high_confidence_FIXED.csv | 25 | Want cleanest data only |
| paintomics_cytokines_FIXED.csv | 3 | Focus on inflammatory pathways |
| paintomics_growth_factors_FIXED.csv | 7 | Focus on growth signaling |
| paintomics_ecm_proteins_FIXED.csv | 10 | Focus on ECM remodeling |

---

## What the Format Looks Like

### ✅ CORRECT (FIXED files):
```csv
ID,Vehicle,Testosterone
FRZB,0,4.42
GAS6,0,3.04
CCL7,0,-2.86
```
- 3 columns
- 2 numeric columns (Vehicle, Testosterone)
- Vehicle = 0 (baseline)
- Testosterone = log2FC

### ❌ INCORRECT (Original files):
```csv
Name,Testosterone_vs_Vehicle
FRZB,4.42
GAS6,3.04
```
- Only 1 numeric column
- PaintOmics error: "Expected at least 2 columns"

---

## Expected Results

### Pathways You Should See:

1. **Cytokine-cytokine receptor interaction** (mmu04060)
   - CCL7, CCL2, CXCL5 (all downregulated)

2. **TGF-beta signaling** (mmu04350)
   - TGFB1, TGFB3, GDF15, AMH (mixed)

3. **PI3K-Akt signaling** (mmu04151)
   - GAS6, AGT, growth factors

4. **ECM-receptor interaction** (mmu04512)
   - VCAN, LUM, BGN, TGFBI (mostly upregulated)

5. **Focal adhesion** (mmu04510)
   - ECM proteins

6. **Complement cascades** (mmu04610)
   - C4A/C4B, SERPINs

---

## Troubleshooting

### "Expected at least 2 columns" error
✅ **FIXED** - Should not happen with _FIXED files
- If it does: Check you selected BOTH "Vehicle" AND "Testosterone" columns

### "Gene IDs not mapped"
- Ensure organism = Mouse (Mus musculus)
- Try high-confidence file instead
- Some genes (C4A/C4B) may not map perfectly

### "No pathways enriched"
- Normal for small gene lists
- Lower p-value threshold to 0.1
- Try main file (33 proteins) instead of specialized files

---

## Top 5 Proteins to Watch

### Upregulated:
1. **FRZB** (+4.42) - Wnt signaling inhibitor
2. **GAS6** (+3.04) - TAM receptor ligand
3. **C4A/C4B** (+2.73) - Complement cascade

### Downregulated:
1. **CCL7** (-2.86) - Chemokine, inflammation
2. **CCL2** (-2.35) - Chemokine, monocyte recruitment

---

## Final Confirmation

✅ All _FIXED.csv files are valid PaintOmics inputs
✅ Will NOT trigger "expected at least 2 columns" error
✅ Format: ID, Vehicle, Testosterone
✅ CSV, UTF-8, comma-separated
✅ Ready for upload NOW

---

## Need Help?

See detailed guides:
- `FIXED_FILES_SUMMARY.md` - Full documentation
- `PAINTOMICS_UPLOAD_INSTRUCTIONS.md` - Comprehensive guide
- `FILE_SUMMARY.md` - File descriptions

---

**Generated**: 2025-12-16
**Status**: ✅ READY FOR PAINTOMICS 4
**Upload**: http://www.paintomics.org/
