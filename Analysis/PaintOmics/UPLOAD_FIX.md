# PaintOmics Upload Fix - "Expected at least 2 columns" Error

## ✅ SOLUTION: Use TXT Files (Tab-Separated)

The error occurs because PaintOmics may be expecting **TAB-separated** values, not comma-separated, despite calling them "CSV" files.

---

## 📁 Files Available

### TXT Files (Tab-Separated) - ⭐ TRY THESE FIRST

```
paintomics_input_expression_data.txt  ⭐ MAIN FILE - USE THIS
paintomics_high_confidence.txt
paintomics_minimal_input.txt
paintomics_cytokines.txt
paintomics_growth_factors.txt
paintomics_ecm_proteins.txt
```

### CSV Files (Comma-Separated) - Try if TXT fails

```
paintomics_input_expression_data.csv
paintomics_high_confidence.csv
paintomics_minimal_input.csv
paintomics_cytokines.csv
paintomics_growth_factors.csv
paintomics_ecm_proteins.csv
```

---

## 🚀 Upload Instructions (UPDATED)

### Step 1: Go to PaintOmics
**URL**: http://www.paintomics.org/
Click: "Start a new job"

### Step 2: Job Settings
- **Job name**: "Testosterone_Proteomics"
- **Organism**: **Mouse (Mus musculus)** ⚠️ CRITICAL
- **Description**: "SILAC proteomics"

### Step 3: Upload File ⚠️ IMPORTANT CHANGE
**Upload this file**: `paintomics_input_expression_data.txt` (NOT .csv)

- Click "Upload data"
- Select file: **paintomics_input_expression_data.txt** ⭐
- Data type: Proteomics
- Format: CSV/TSV (or Text file)

### Step 4: Column Settings
- **Has header**: ✅ YES (check this!)
- **Delimiter**: If asked, select **TAB** or **Auto-detect**
- **ID column**: Select "ID"
- **Condition columns**: Select **BOTH "Vehicle" AND "Testosterone"**

### Step 5: Submit
- Click "Submit"
- Wait 2-10 minutes

---

## 📊 Format Comparison

### TXT Format (Tab-Separated) - ⭐ USE THIS:
```
ID	Vehicle	Testosterone
FRZB	0	4.42
GAS6	0	3.04
CCL7	0	-2.86
```
**Delimiter**: TAB character (not visible, but columns are separated by tabs)

### CSV Format (Comma-Separated) - Backup:
```
ID,Vehicle,Testosterone
FRZB,0,4.42
GAS6,0,3.04
CCL7,0,-2.86
```
**Delimiter**: Comma (,)

---

## 🔧 If TXT Files Work

If the `.txt` files upload successfully:
✅ You're done! Proceed with analysis
✅ The "expected at least 2 columns" error is fixed
✅ PaintOmics will generate pathway diagrams

---

## 🔧 If TXT Files Still Fail

Try these alternatives in order:

### Option 1: Try CSV Files
Upload: `paintomics_input_expression_data.csv`
- Sometimes PaintOmics auto-detects the delimiter

### Option 2: Check Delimiter Setting
- In upload dialog, look for "Delimiter" or "Separator" option
- Try selecting: TAB, Comma, or Auto-detect

### Option 3: Try Different Browser
- Some browsers handle file uploads differently
- Try: Chrome, Firefox, or Safari
- Clear cache and cookies first

### Option 4: Manual Column Specification
- After upload, if it asks about columns:
  - Manually specify that columns are separated by TAB or COMMA
  - Explicitly select ID, Vehicle, Testosterone columns

---

## 📋 File Format Verification

### TXT Files (Tab-Separated):
- ✅ 3 columns: ID, Vehicle, Testosterone
- ✅ Delimiter: TAB (\t)
- ✅ Header: YES
- ✅ Encoding: UTF-8
- ✅ Line endings: Unix (\n)
- ✅ No BOM
- ✅ 2 numeric columns (Vehicle, Testosterone)

### CSV Files (Comma-Separated):
- ✅ 3 columns: ID, Vehicle, Testosterone
- ✅ Delimiter: COMMA (,)
- ✅ Header: YES
- ✅ Encoding: UTF-8
- ✅ Line endings: Unix (\n)
- ✅ No BOM
- ✅ 2 numeric columns (Vehicle, Testosterone)

---

## 🎯 Why TXT Files Might Work Better

Many bioinformatics tools (including PaintOmics) prefer **TAB-separated** files because:

1. **Gene names may contain commas** (e.g., "protein, isoform 1")
   - Tabs avoid this conflict
   - Commas can cause parsing errors

2. **TAB is unambiguous**
   - No confusion with decimal separators
   - No need for quoting

3. **Standard in bioinformatics**
   - Many tools default to TAB-separated
   - Even when labeled as "CSV"

---

## 📊 File Statistics (Main File)

**File**: paintomics_input_expression_data.txt

- **Proteins**: 33
- **Format**: Tab-separated (TSV)
- **Columns**: ID, Vehicle, Testosterone
- **Testosterone range**: -2.859 to +4.420
- **Upregulated**: 21 proteins
- **Downregulated**: 12 proteins

**Top proteins**:
- FRZB: +4.42 (most upregulated)
- GAS6: +3.04
- C4A/C4B: +2.73
- CCL7: -2.86 (most downregulated)
- CCL2: -2.35

---

## 🎯 Expected Pathways (After Upload)

If upload succeeds, you should see these KEGG pathways:

1. **Cytokine-cytokine receptor interaction** (mmu04060)
2. **TGF-beta signaling** (mmu04350)
3. **PI3K-Akt signaling** (mmu04151)
4. **ECM-receptor interaction** (mmu04512)
5. **Focal adhesion** (mmu04510)
6. **Complement and coagulation cascades** (mmu04610)

---

## 📖 Quick Reference

| Try This | File | Format | When |
|----------|------|--------|------|
| 1st ⭐ | paintomics_input_expression_data.txt | TAB | Default choice |
| 2nd | paintomics_input_expression_data.csv | COMMA | If TXT fails |
| 3rd | paintomics_high_confidence.txt | TAB | For cleanest data |
| 4th | paintomics_high_confidence.csv | COMMA | Backup option |

---

## ✅ Final Checklist

Before uploading, verify:

- ✅ Using `.txt` file (tab-separated)
- ✅ Organism set to **Mouse (Mus musculus)**
- ✅ "Has header" is checked
- ✅ ID column selected
- ✅ **BOTH** Vehicle and Testosterone selected as conditions
- ✅ Delimiter set to TAB (if asked)

---

## 🔍 Still Getting Errors?

If you still get "Expected at least 2 columns" error:

1. **Check organism**: Must be Mouse (Mus musculus)
2. **Check delimiter**: Try manually selecting TAB
3. **Check browser**: Try different browser or incognito mode
4. **Check file**: Ensure you uploaded .txt file, not .csv
5. **Check columns**: Ensure BOTH conditions are selected

---

## 📞 Contact Info

If problems persist:
- Check PaintOmics documentation: http://www.paintomics.org/tutorial/
- Contact PaintOmics support through their website
- Verify file format using a text editor (should show tabs between columns)

---

**Generated**: 2025-12-16
**Status**: ✅ BOTH TXT AND CSV FORMATS AVAILABLE
**Recommendation**: ⭐ TRY TXT FILES FIRST (Tab-separated)
**Main file**: paintomics_input_expression_data.txt
