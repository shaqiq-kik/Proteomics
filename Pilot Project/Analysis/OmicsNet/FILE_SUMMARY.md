# OmicsNet 2.0 Input Files - Summary

## Quick Start

**RECOMMENDED FILE**: `omicsnet_input_with_expression.txt`
- Upload this file first to OmicsNet 2.0
- 33 proteins with expression data
- Best for comprehensive network analysis

**INSTRUCTIONS**: See `OMICSNET_UPLOAD_INSTRUCTIONS.md` for detailed step-by-step guide

---

## All Files Generated

### Primary Analysis Files

| File | Proteins | Description | Format |
|------|----------|-------------|--------|
| `omicsnet_input_with_expression.txt` | 33 | Main file (High + Medium confidence) | 3-column with expression |
| `omicsnet_high_confidence_only.txt` | 25 | High confidence only | 3-column with expression |
| `omicsnet_genelist_only.txt` | 33 | Gene list only (no expression) | 1-column gene list |

### Specialized Analysis Files

| File | Proteins | Description | Format |
|------|----------|-------------|--------|
| `omicsnet_cytokines.txt` | 3 | Cytokines only (CCL7, CCL2, CXCL5) | 3-column with expression |
| `omicsnet_growth_factors.txt` | 7 | Growth factors only | 3-column with expression |
| `omicsnet_signaling_proteins.txt` | 10 | Cytokines + Growth factors | 3-column with expression |
| `omicsnet_upregulated.txt` | 21 | Upregulated proteins only | 3-column with expression |
| `omicsnet_downregulated.txt` | 12 | Downregulated proteins only | 3-column with expression |

---

## File Format Details

### 3-Column Format (with expression)
- **Column 1**: Gene symbol
- **Column 2**: log2 fold change (numeric)
- **Column 3**: "Protein" (literal text)
- **Delimiter**: Tab
- **Header**: NO HEADER

Example:
```
FRZB	4.42	Protein
GAS6	3.04	Protein
CCL7	-2.86	Protein
```

### 1-Column Format (gene list only)
- One gene symbol per line
- No header
- Plain text

Example:
```
FRZB
GAS6
CCL7
```

---

## Data Filtering Applied

✅ **INCLUDED**:
- High confidence proteins (25)
- Medium confidence proteins (8)
- Complete_Suppression proteins (valid log2FC values)

❌ **EXCLUDED**:
- Low confidence proteins
- Proteins with FC_Type == 'Cannot_Calculate'
- Proteins with NaN fold changes

---

## Top Proteins by Absolute Fold Change

### Most Changed (Top 10):
1. **FRZB** - log2FC: +4.42 (21.4x upregulated)
2. **CRISP1/CRISP3** - log2FC: +3.69 (12.9x upregulated)
3. **GAS6** - log2FC: +3.04 (8.2x upregulated)
4. **CCL7** - log2FC: -2.86 (7.3x downregulated)
5. **C4A/C4B** - log2FC: +2.73 (6.7x upregulated)
6. **CCL2** - log2FC: -2.35 (5.1x downregulated)
7. **CST12** - log2FC: +2.35 (5.1x upregulated)
8. **SLIT3** - log2FC: +1.95 (3.9x upregulated)
9. **VCAN** - log2FC: +1.81 (3.5x upregulated)
10. **LUM** - log2FC: +1.80 (3.5x upregulated)

---

## Expected Network Insights

### Key Pathways Expected:
1. **TGF-beta signaling** (TGFB1, TGFB3, GDF15, TGFBI)
2. **Cytokine signaling** (CCL7, CCL2, CXCL5)
3. **Growth factor signaling** (GAS6, AGT, AMH, NENF)
4. **ECM remodeling** (VCAN, LUM, BGN, SERPINE2)
5. **Complement cascade** (C4A/C4B)

### Regulation Pattern:
- **Upregulated**: 21 proteins (mostly growth factors, ECM proteins)
- **Downregulated**: 12 proteins (mostly cytokines, chemokines)

---

## Upload to OmicsNet 2.0

1. Go to: https://www.omicsnet.ca/
2. Click "Get Started"
3. Select organism: **Mouse (Mus musculus)**
4. Select ID type: **Gene Symbol**
5. Upload your chosen file
6. Data type: **Protein**
7. See `OMICSNET_UPLOAD_INSTRUCTIONS.md` for detailed settings

---

## File Verification

All files have been verified to:
- ✅ Have correct format (tab-delimited for 3-column, plain text for gene list)
- ✅ Have NO header row
- ✅ Contain valid gene symbols
- ✅ Be sorted by absolute fold change (for expression files)
- ✅ Include only filtered, high-quality data

---

## Regenerating Files

If you need to regenerate these files (e.g., with different filtering):

```bash
cd Analysis/OmicsNet
python3 generate_omicsnet_files.py
```

The script reads from:
```
../General Analysis/cleaned_proteomics_data_with_QC_flags.csv
```

---

**Generated**: 2025-12-16
**Source**: Testosterone vs. Vehicle SILAC proteomics
**Total proteins analyzed**: 33 (from 42 total, after QC filtering)
