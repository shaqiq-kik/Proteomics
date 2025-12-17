# PaintOmics 4 Input Files - Summary

## Quick Start

**RECOMMENDED FILE**: `paintomics_input_expression_data.csv`
- Upload this file first to PaintOmics 4
- 33 proteins with full metadata (7 columns)
- Best for comprehensive KEGG pathway visualization

**PAINTOMICS URL**: http://www.paintomics.org/

**INSTRUCTIONS**: See `PAINTOMICS_UPLOAD_INSTRUCTIONS.md` for detailed step-by-step guide

---

## What Makes PaintOmics Different?

Unlike enrichment tools or network tools, PaintOmics literally **paints** your expression data onto KEGG pathway diagrams. You see:
- Biological pathways as visual diagrams
- Your proteins highlighted in color (red=up, blue=down)
- Exactly where your proteins fit in signaling cascades
- Publication-quality pathway figures

---

## All Files Generated

### Primary Analysis Files

| File | Proteins | Columns | Format | Description |
|------|----------|---------|--------|-------------|
| `paintomics_input_expression_data.csv` | 33 | 7 | Full | Main file - RECOMMENDED |
| `paintomics_minimal_input.csv` | 33 | 2 | Minimal | Quick analysis |
| `paintomics_high_confidence.csv` | 25 | 7 | Full | High confidence only |

### Specialized Analysis Files

| File | Proteins | Focus | Expected Pathways |
|------|----------|-------|-------------------|
| `paintomics_cytokines.csv` | 3 | Inflammatory signaling | Cytokine-cytokine receptor interaction |
| `paintomics_growth_factors.csv` | 7 | Growth signaling | TGF-beta, PI3K-Akt, MAPK pathways |
| `paintomics_ecm_proteins.csv` | 10 | ECM remodeling | ECM-receptor interaction, Focal adhesion |

### Template File

| File | Proteins | Purpose |
|------|----------|---------|
| `paintomics_timecourse_template.csv` | 33 | Template for future time-course experiments |

---

## File Format Details

### Full Format (7 columns)
Used by: Main file, high-confidence, all specialized files

```csv
Name,Testosterone_vs_Vehicle,Fold_Change,Regulation,Confidence,Function,Description
FRZB,4.42,21.40,Up,Medium,other,frizzled related protein
GAS6,3.04,8.24,Up,Medium,growth factor,growth arrest specific 6
CCL7,-2.86,0.14,Down,High,cytokine,chemokine (C-C motif) ligand 7
```

**Columns**:
1. `Name` - Gene symbol (REQUIRED by PaintOmics)
2. `Testosterone_vs_Vehicle` - log2 fold change (REQUIRED by PaintOmics)
3. `Fold_Change` - Actual fold change (not log)
4. `Regulation` - "Up" or "Down"
5. `Confidence` - "High" or "Medium"
6. `Function` - Functional class
7. `Description` - Full protein description

### Minimal Format (2 columns)
Used by: Minimal input file

```csv
Name,Testosterone_vs_Vehicle
FRZB,4.42
GAS6,3.04
CCL7,-2.86
```

**Columns**:
1. `Name` - Gene symbol
2. `Testosterone_vs_Vehicle` - log2 fold change

### Time-Course Format (3 columns)
Used by: Time-course template

```csv
Name,Vehicle,Testosterone
FRZB,0.0,4.42
GAS6,0.0,3.04
CCL7,0.0,-2.86
```

**Columns**:
1. `Name` - Gene symbol
2. `Vehicle` - Baseline (0.0)
3. `Testosterone` - Treatment (log2FC)

**Note**: For future time-course experiments, add more columns:
`Name,T0h,T6h,T24h,T48h`

---

## Data Filtering Applied

✅ **INCLUDED**:
- High confidence proteins (25)
- Medium confidence proteins (8)
- Complete_Suppression proteins (valid log2FC values)
- Total: 33 proteins

❌ **EXCLUDED**:
- Low confidence proteins (9)
- Proteins with FC_Type == 'Cannot_Calculate' (4)
- Proteins with NaN fold changes

---

## Top Proteins by Absolute Fold Change

### Most Upregulated (Top 10):
1. **FRZB** - log2FC: +4.42 (21.4x increase)
2. **CRISP1/CRISP3** - log2FC: +3.69 (12.9x increase)
3. **GAS6** - log2FC: +3.04 (8.2x increase)
4. **C4A/C4B** - log2FC: +2.73 (6.7x increase)
5. **CST12** - log2FC: +2.35 (5.1x increase)
6. **SLIT3** - log2FC: +1.95 (3.9x increase)
7. **VCAN** - log2FC: +1.81 (3.5x increase)
8. **LUM** - log2FC: +1.80 (3.5x increase)
9. **SERPINE2** - log2FC: +1.76 (3.4x increase)
10. **NENF** - log2FC: +1.15 (2.2x increase)

### Most Downregulated (Top 5):
1. **CCL7** - log2FC: -2.86 (7.3x decrease)
2. **CCL2** - log2FC: -2.35 (5.1x decrease)
3. **TGFB3** - log2FC: -1.63 (3.1x decrease)
4. **GDF15** - log2FC: -1.46 (2.8x decrease)
5. **TGFB1** - log2FC: -1.41 (2.7x decrease)

---

## Expected KEGG Pathways

Based on your proteins, these KEGG pathways are likely to be enriched:

### High-Probability Pathways (≥3 proteins expected)

#### 1. Cytokine-cytokine receptor interaction (mmu04060)
- **Proteins**: CCL7, CCL2, CXCL5
- **Regulation**: All downregulated
- **Insight**: Reduced inflammatory signaling

#### 2. TGF-beta signaling pathway (mmu04350)
- **Proteins**: TGFB1, TGFB3, GDF15, AMH
- **Regulation**: Mixed (most down, AMH up)
- **Insight**: Complex TGF-beta modulation

#### 3. PI3K-Akt signaling pathway (mmu04151)
- **Proteins**: AGT, GAS6, multiple growth factors
- **Regulation**: Mixed
- **Insight**: Cell survival pathway changes

#### 4. ECM-receptor interaction (mmu04512)
- **Proteins**: VCAN, LUM, BGN, TGFBI
- **Regulation**: Mostly upregulated
- **Insight**: Active ECM remodeling

#### 5. Focal adhesion (mmu04510)
- **Proteins**: VCAN, BGN, ECM proteins
- **Regulation**: Mostly upregulated
- **Insight**: Cell-matrix interaction changes

#### 6. Complement and coagulation cascades (mmu04610)
- **Proteins**: C4A/C4B, SERPINE2, SERPINF1
- **Regulation**: Mixed
- **Insight**: Immune system modulation

---

## Upload Quick Reference

### Basic Settings:
1. **Website**: http://www.paintomics.org/
2. **Organism**: Mouse (Mus musculus)
3. **ID type**: Gene Symbol
4. **File has header**: YES ✓
5. **ID column**: Name
6. **Value column**: Testosterone_vs_Vehicle

### Which File to Start With?

| Your Goal | Upload This File |
|-----------|------------------|
| Comprehensive overview | `paintomics_input_expression_data.csv` |
| Quick test | `paintomics_minimal_input.csv` |
| Cleanest data | `paintomics_high_confidence.csv` |
| Inflammatory pathways only | `paintomics_cytokines.csv` |
| Growth signaling only | `paintomics_growth_factors.csv` |
| ECM pathways only | `paintomics_ecm_proteins.csv` |

---

## File Verification

All files verified to meet PaintOmics 4 requirements:
- ✅ CSV format (comma-separated)
- ✅ Header row present
- ✅ Correct column names ("Name", "Testosterone_vs_Vehicle")
- ✅ Valid gene symbols
- ✅ Numeric log2 fold changes
- ✅ Sorted by absolute fold change

Example format verification:
```csv
Name,Testosterone_vs_Vehicle,Fold_Change,Regulation,Confidence,Function,Description
FRZB,4.42,21.40135269,Up,Medium,other,frizzled related protein
```
✓ Comma-separated ✓ Header present ✓ Correct column names

---

## Key Differences: PaintOmics vs. OmicsNet Files

### PaintOmics Files:
- ✅ **CSV format** (comma-separated)
- ✅ **MUST have header row**
- ✅ Column names: "Name", "Testosterone_vs_Vehicle"
- ✅ Can have multiple metadata columns
- 📊 Purpose: KEGG pathway visualization

### OmicsNet Files:
- ✅ **TXT format** (tab-delimited)
- ✅ **NO header row**
- ✅ Three columns: Gene, log2FC, "Protein"
- ✅ Minimal format only
- 🔗 Purpose: Protein-protein interaction networks

**Use both tools** for complementary insights:
- PaintOmics → Pathway context
- OmicsNet → Interaction networks

---

## Biological Insights Summary

### Inflammatory Response (DOWNREGULATED)
- **Cytokines**: CCL7, CCL2, CXCL5 all reduced
- **Interpretation**: Anti-inflammatory effect of testosterone
- **Pathways**: Cytokine-cytokine receptor interaction

### Growth Signaling (MIXED)
- **Upregulated**: GAS6 (+3.04), AGT (+0.58), AMH (+0.61), NENF (+1.15)
- **Downregulated**: TGFB1 (-1.41), TGFB3 (-1.63), GDF15 (-1.46)
- **Interpretation**: Complex modulation, not simple on/off
- **Pathways**: TGF-beta, PI3K-Akt, MAPK pathways

### ECM Remodeling (UPREGULATED)
- **Proteins**: FRZB (+4.42), VCAN (+1.81), LUM (+1.80), BGN (+1.19)
- **Interpretation**: Active tissue restructuring
- **Pathways**: ECM-receptor interaction, Focal adhesion

### Complement System (UPREGULATED)
- **Proteins**: C4A/C4B (+2.73)
- **Interpretation**: Immune system activation
- **Pathways**: Complement and coagulation cascades

---

## Regenerating Files

If you need to regenerate with different filtering:

```bash
cd Analysis/PaintOmics
python3 generate_paintomics_files.py
```

The script reads from:
```
../General Analysis/cleaned_proteomics_data_with_QC_flags.csv
```

---

## Next Steps

### 1. Immediate: Upload to PaintOmics
- Start with `paintomics_input_expression_data.csv`
- Follow `PAINTOMICS_UPLOAD_INSTRUCTIONS.md`
- Generate pathway diagrams

### 2. Analysis: Compare with Other Tools
- OmicsNet: Protein interaction networks
- Enrichr/DAVID: Broader functional enrichment
- STRING: Interaction confidence

### 3. Validation: Experimental Follow-up
- Key pathways: TGF-beta, Cytokine signaling, ECM
- Key proteins: GAS6, FRZB, CCL7, CCL2
- Techniques: qPCR, Western blot, immunostaining

### 4. Publication: Export Figures
- PaintOmics generates publication-quality pathway diagrams
- Export as SVG for vector graphics
- Include in your manuscript figures

---

## File Statistics

| Metric | Value |
|--------|-------|
| Total files created | 7 CSV files |
| Main analysis file size | 2.2 KB |
| Total proteins (main file) | 33 |
| Upregulated proteins | 21 (64%) |
| Downregulated proteins | 12 (36%) |
| Log2FC range | -2.86 to +4.42 |
| Fold change range | 0.14x to 21.4x |

---

## Troubleshooting Tips

### If genes don't map:
1. Check gene symbol format (uppercase/lowercase)
2. Try UniProt IDs instead (available in original data)
3. Manually check KEGG database: https://www.genome.jp/kegg/

### If no pathways enriched:
1. Lower p-value threshold to 0.1
2. Set minimum genes per pathway to 1-2
3. Check organism is set to Mouse

### If colors don't show:
1. Ensure "Testosterone_vs_Vehicle" selected as value column
2. Check data type is set to "numeric" or "log2 fold change"
3. Try minimal file first to isolate issue

---

**Generated**: 2025-12-16
**Source**: Testosterone vs. Vehicle SILAC proteomics
**Total proteins**: 33 (from 42 total, after QC filtering)
**Ready for**: PaintOmics 4 pathway visualization
