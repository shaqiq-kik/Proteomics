# 🚀 PaintOmics Upload - START HERE

## ⚠️ IMPORTANT: "Expected at least 2 columns" Error - FIXED

**Solution**: Use **TXT files** (tab-separated), not CSV files.

---

## 📁 Upload This File:

### ⭐ MAIN FILE (USE THIS):
```
paintomics_input_expression_data.txt
```

**Format**: Tab-separated (TSV)
**Proteins**: 33
**Status**: ✅ READY

---

## 🚀 Quick Upload (3 Steps)

### 1. Go to PaintOmics
**URL**: http://www.paintomics.org/
Click: "Start a new job"

### 2. Upload Settings
- **Organism**: **Mouse (Mus musculus)** ⚠️
- **Upload file**: `paintomics_input_expression_data.txt`
- **Has header**: ✅ YES
- **Delimiter**: TAB (or auto-detect)

### 3. Column Selection
- **ID column**: ID
- **Conditions**: Select **BOTH** "Vehicle" AND "Testosterone"

**Click Submit** → Wait 2-10 minutes → Done!

---

## 📋 Why TXT Instead of CSV?

PaintOmics expects **TAB-separated** values, even though it says "CSV".

**TXT files** (tab-separated):
```
ID	Vehicle	Testosterone
FRZB	0	4.42
GAS6	0	3.04
```
✅ Works with PaintOmics

**CSV files** (comma-separated):
```
ID,Vehicle,Testosterone
FRZB,0,4.42
GAS6,0,3.04
```
⚠️ May cause "expected at least 2 columns" error

---

## 📁 All Available Files

### TXT Files (Tab-Separated) - ⭐ RECOMMENDED

| File | Proteins | Use For |
|------|----------|---------|
| **paintomics_input_expression_data.txt** ⭐ | 33 | Main analysis (START HERE) |
| paintomics_high_confidence.txt | 25 | Cleanest data only |
| paintomics_cytokines.txt | 3 | Inflammatory pathways |
| paintomics_growth_factors.txt | 7 | Growth signaling |
| paintomics_ecm_proteins.txt | 10 | ECM remodeling |

### CSV Files (Comma-Separated) - Backup

| File | Proteins | Use For |
|------|----------|---------|
| paintomics_input_expression_data.csv | 33 | Try if TXT fails |
| paintomics_high_confidence.csv | 25 | Backup option |
| paintomics_cytokines.csv | 3 | Backup option |
| paintomics_growth_factors.csv | 7 | Backup option |
| paintomics_ecm_proteins.csv | 10 | Backup option |

---

## 🎯 Expected Results

After successful upload, you should see these KEGG pathways:

1. ✅ **Cytokine-cytokine receptor interaction** (mmu04060)
   - CCL7, CCL2, CXCL5 (all downregulated)

2. ✅ **TGF-beta signaling** (mmu04350)
   - TGFB1, TGFB3, GDF15, AMH (mixed)

3. ✅ **PI3K-Akt signaling** (mmu04151)
   - GAS6, AGT, growth factors

4. ✅ **ECM-receptor interaction** (mmu04512)
   - VCAN, LUM, BGN (upregulated)

5. ✅ **Focal adhesion** (mmu04510)
   - ECM proteins

6. ✅ **Complement cascades** (mmu04610)
   - C4A/C4B, SERPINs

---

## 🔧 Troubleshooting

### If you still get "Expected at least 2 columns":

1. **Check delimiter**: Ensure TAB is selected (or auto-detect)
2. **Check file**: Make sure you uploaded `.txt` file, not `.csv`
3. **Check columns**: BOTH "Vehicle" and "Testosterone" must be selected
4. **Try browser**: Use Chrome or Firefox in incognito mode
5. **Try CSV**: If all else fails, try the `.csv` version

---

## 📊 File Format Confirmed

All files verified:
- ✅ TXT files: Tab-separated (\t delimiter)
- ✅ CSV files: Comma-separated (, delimiter)
- ✅ Format: ID, Vehicle, Testosterone
- ✅ Header row: YES
- ✅ Encoding: UTF-8
- ✅ 2 numeric columns (Vehicle=0, Testosterone=log2FC)
- ✅ No duplicates, no missing values

---

## 📖 Need More Help?

See these guides:
- **UPLOAD_FIX.md** - Detailed error fix instructions
- **README.md** - Complete file documentation
- **QUICK_START.md** - Original upload guide
- **PAINTOMICS_UPLOAD_INSTRUCTIONS.md** - Comprehensive guide

---

## ✅ Final Checklist

Before clicking Submit:
- ✅ File: `paintomics_input_expression_data.txt`
- ✅ Organism: Mouse (Mus musculus)
- ✅ Has header: YES
- ✅ Delimiter: TAB
- ✅ ID column: Selected
- ✅ Conditions: BOTH Vehicle AND Testosterone selected

---

**Generated**: 2025-12-16
**Status**: ✅ READY - ERROR FIXED
**Upload**: http://www.paintomics.org/
**Main file**: paintomics_input_expression_data.txt (TAB-separated)
