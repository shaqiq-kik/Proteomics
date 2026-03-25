# PaintOmics 4 Pathway Visualization Guide

## What is PaintOmics?

**PaintOmics 4** is a web-based tool that takes your expression data and literally **PAINTS** it onto KEGG pathway diagrams. Instead of seeing abstract lists of proteins, you see:

- **Visual pathway diagrams** with your proteins highlighted
- **Color-coded nodes** showing upregulation (red) and downregulation (blue/green)
- **Pathway components** in biological context (receptors, enzymes, signaling cascades)
- **Integrated multi-omics** (can combine proteomics, metabolomics, transcriptomics)

**PaintOmics 4 URL**: http://www.paintomics.org/

---

## Files Created

### Primary Analysis Files

#### 1. **paintomics_input_expression_data.csv** ⭐ RECOMMENDED - MAIN FILE
- **Proteins**: 33 (High + Medium confidence)
- **Format**: CSV with 7 columns (full metadata)
- **Columns**:
  - `Name` - Gene symbol
  - `Testosterone_vs_Vehicle` - log2 fold change
  - `Fold_Change` - Actual fold change (not log)
  - `Regulation` - "Up" or "Down"
  - `Confidence` - Confidence level
  - `Function` - Functional class
  - `Description` - Protein description
- **Best for**: Comprehensive pathway visualization with full annotations

#### 2. **paintomics_minimal_input.csv**
- **Proteins**: 33
- **Format**: CSV with 2 columns (simplest format)
- **Columns**: `Name`, `Testosterone_vs_Vehicle`
- **Best for**: Quick analysis, testing, or when you just want basic coloring

#### 3. **paintomics_high_confidence.csv**
- **Proteins**: 25 (High confidence only)
- **Format**: CSV with 7 columns (full metadata)
- **Best for**: Ultra-clean analysis if main file has too many unmapped genes

---

### Specialized Analysis Files

#### 4. **paintomics_cytokines.csv**
- **Proteins**: 3 cytokines (CCL7, CCL2, CXCL5)
- **All downregulated** (log2FC: -2.86 to -1.00)
- **Best for**: Cytokine-cytokine receptor interaction pathway (mmu04060)
- **Expected insight**: Reduced inflammatory signaling

#### 5. **paintomics_growth_factors.csv**
- **Proteins**: 7 growth factors (GAS6, TGFB3, GDF15, TGFB1, NENF, AMH, AGT)
- **Mixed regulation**
- **Best for**: TGF-beta signaling (mmu04350), PI3K-Akt (mmu04151), MAPK pathways
- **Expected insight**: Complex growth factor signaling modulation

#### 6. **paintomics_ecm_proteins.csv**
- **Proteins**: 10 ECM-related proteins (FRZB, VCAN, LUM, BGN, SLIT3, etc.)
- **Mostly upregulated**
- **Best for**: ECM-receptor interaction (mmu04512), Focal adhesion (mmu04510)
- **Expected insight**: Tissue remodeling and ECM reorganization

---

### Template for Future Use

#### 7. **paintomics_timecourse_template.csv**
- **Format**: CSV with 3 columns (Name, Vehicle, Testosterone)
- **Current use**: Shows baseline (Vehicle=0) vs. treatment (Testosterone=log2FC)
- **Future use**: Template for time-course experiments
  - Example future format: `Name, T0h, T6h, T24h, T48h`
  - Can track how expression changes over time

---

## Step-by-Step Upload Instructions

### Step 1: Navigate to PaintOmics 4

1. Go to: **http://www.paintomics.org/**
2. Click **"Start a new job"** or **"Upload your data"**

### Step 2: Create a New Job

1. **Job Name**: Enter a descriptive name
   - Example: `Testosterone_Proteomics_KEGG_Pathways`

2. **Job Description** (optional but recommended):
   - Example: "SILAC proteomics comparing testosterone vs. vehicle treatment in mouse prostate tissue. 33 proteins analyzed."

3. **Select organism**: **Mouse (Mus musculus)**
   - CRITICAL: Make sure you select the correct organism!
   - Look for "mmu" prefix for mouse KEGG pathways

### Step 3: Upload Expression Data

#### Option A: Main File (RECOMMENDED for first analysis)
1. **Upload file**: `paintomics_input_expression_data.csv`
2. **Data type**: Proteomics / Protein Expression
3. **File format**: CSV (comma-separated)
4. **Has header**: YES (check this box!)
5. **ID column**: Select "Name"
6. **Value columns**: Select "Testosterone_vs_Vehicle"
   - You can also select additional columns for metadata viewing

#### Option B: Minimal File (For quick testing)
1. **Upload file**: `paintomics_minimal_input.csv`
2. Settings same as above
3. Only has 2 columns, simpler

#### Option C: Specialized Files
- Upload `paintomics_cytokines.csv` for inflammatory pathway focus
- Upload `paintomics_growth_factors.csv` for growth signaling focus
- Upload `paintomics_ecm_proteins.csv` for ECM pathway focus

### Step 4: Data Processing Settings

**Gene/Protein ID Mapping:**
- **ID type**: Gene Symbol (or Gene Name)
- **Alternative**: If mapping fails, try UniProt ID (you have this in your original data)
- PaintOmics will attempt to map your gene symbols to KEGG IDs

**Expression value interpretation:**
- **Data type**: Log2 fold change
- **Scale**: Already log-transformed (check this if asked)
- **Reference**: Vehicle (control)
- **Treatment**: Testosterone

### Step 5: Pathway Analysis Settings

**Recommended Settings:**

1. **Pathway database**: KEGG Pathways (default)

2. **Pathway enrichment**:
   - Enable pathway enrichment analysis
   - Method: Fisher's exact test or Hypergeometric test
   - P-value threshold: 0.05
   - Correction: Benjamini-Hochberg (FDR)

3. **Visualization threshold**:
   - Display pathways with at least **2-3 mapped genes**
   - For small gene lists (e.g., cytokines file with 3 genes), set minimum to 1

4. **Color scheme**:
   - Upregulated: Red/Orange
   - Downregulated: Blue/Green
   - Not significant: Gray
   - Gradient scale: Automatic (based on log2FC range)

### Step 6: Run Analysis

1. Click **"Submit"** or **"Start Analysis"**
2. Wait for processing (usually 2-10 minutes)
3. You'll receive a Job ID - save this to access results later

---

## Understanding Your Results

### Result Sections in PaintOmics

#### 1. **Pathway Overview**
- **Enriched Pathways Table**: Lists pathways with significant protein mapping
- **Columns**:
  - Pathway name
  - Number of matched proteins
  - P-value (enrichment significance)
  - Adjusted P-value (FDR)
  - List of matched proteins

#### 2. **Interactive Pathway Diagrams**
- Click on any pathway to view the colored diagram
- **Node colors**:
  - **Red**: Upregulated proteins
  - **Blue/Green**: Downregulated proteins
  - **Gray**: Proteins in pathway but not in your dataset
  - **White**: Other pathway components (not proteins)

#### 3. **Comparative Views**
- Can compare different conditions (if you uploaded time-course data)
- Side-by-side pathway visualization

---

## Expected KEGG Pathways for Your Data

### High-Probability Pathways (Should appear with multiple mapped proteins)

#### 1. **Cytokine-cytokine receptor interaction** (mmu04060)
**Matched proteins expected**: CCL7, CCL2, CXCL5
**Regulation**: All downregulated (-2.86 to -1.00 log2FC)
**Biological interpretation**:
- Reduced chemokine signaling
- Decreased monocyte and neutrophil recruitment
- Lower inflammatory signaling
**Why important**: Suggests anti-inflammatory effects of testosterone

#### 2. **TGF-beta signaling pathway** (mmu04350)
**Matched proteins expected**: TGFB1, TGFB3, GDF15, AMH
**Regulation**: Mixed (TGFB1/3/GDF15 downregulated, AMH upregulated)
**Biological interpretation**:
- Complex modulation of TGF-beta pathway
- TGF-beta ligands reduced but AMH increased
- May affect cell differentiation and fibrosis
**Why important**: TGF-beta is critical for prostate development and homeostasis

#### 3. **PI3K-Akt signaling pathway** (mmu04151)
**Matched proteins expected**: AGT, GAS6, multiple growth factors
**Regulation**: Mixed
**Biological interpretation**:
- Cell survival and proliferation pathway activation
- GAS6 strongly upregulated (+3.04 log2FC) - TAM receptor ligand
- AGT upregulated (+0.58 log2FC)
**Why important**: Key pathway for cell growth and survival

#### 4. **ECM-receptor interaction** (mmu04512)
**Matched proteins expected**: VCAN, LUM, BGN, TGFBI
**Regulation**: Mostly upregulated
**Biological interpretation**:
- Extracellular matrix remodeling
- Proteoglycans increased (VCAN, LUM, BGN)
- Tissue structure reorganization
**Why important**: ECM changes affect tissue architecture

#### 5. **Focal adhesion** (mmu04510)
**Matched proteins expected**: VCAN, BGN, multiple ECM proteins
**Regulation**: Mostly upregulated
**Biological interpretation**:
- Cell-matrix adhesion changes
- Integrin signaling modulation
**Why important**: Affects cell migration and tissue integrity

#### 6. **Complement and coagulation cascades** (mmu04610)
**Matched proteins expected**: C4A/C4B, SERPINE2, SERPINF1
**Regulation**: C4A/C4B upregulated (+2.73), SERPINs mixed
**Biological interpretation**:
- Complement activation
- Protease inhibitor regulation
**Why important**: Immune system modulation

### Moderate-Probability Pathways (May appear with 1-2 mapped proteins)

7. **MAPK signaling pathway** (mmu04010)
8. **Ras signaling pathway** (mmu04014)
9. **Rap1 signaling pathway** (mmu04015)
10. **HIF-1 signaling pathway** (mmu04066)
11. **Protein digestion and absorption** (mmu04974)
12. **Proteoglycans in cancer** (mmu05205)

---

## Tips for Best Results

### Data Quality Tips

1. **Gene symbol formatting**:
   - PaintOmics is sensitive to gene symbol formatting
   - If mapping fails, try:
     - Converting to all uppercase
     - Removing special characters
     - Using UniProt IDs instead

2. **Problematic gene names**:
   - `C4A/C4B` - May need to split into separate entries
   - `CRISP1/CRISP3` - May need to choose one or split
   - If unmapped, check KEGG database for correct mouse gene symbols

3. **Minimum dataset size**:
   - PaintOmics works best with at least 10-20 genes
   - For smaller lists (like cytokines file with 3 genes), lower enrichment thresholds

### Visualization Tips

1. **Pathway selection**:
   - Start with pathways showing highest enrichment (lowest p-value)
   - Focus on pathways with ≥3 mapped proteins for robust conclusions

2. **Color scale adjustment**:
   - If all proteins are similar fold change, adjust color scale manually
   - Helps highlight subtle differences

3. **Multiple pathway views**:
   - Look at related pathways together
   - Example: Compare TGF-beta pathway with PI3K-Akt pathway to see crosstalk

4. **Export options**:
   - Export pathway diagrams as high-resolution images (PNG, SVG)
   - Export pathway tables as CSV for further analysis
   - Save session for later review (use your Job ID)

---

## Comparison: PaintOmics vs. Other Tools

### PaintOmics Strengths:
✅ Visual pathway diagrams (KEGG integration)
✅ See exactly where proteins fit in biological pathways
✅ Multi-omics integration capability
✅ Publication-quality pathway figures
✅ Easy to understand for non-bioinformaticians

### When to use other tools instead:

**Use STRING/OmicsNet when you want**:
- Protein-protein interaction networks
- Hub protein identification
- Network topology analysis
- Interaction confidence scores

**Use Enrichr/DAVID when you want**:
- Broader database coverage (not just KEGG)
- Gene Ontology (GO) enrichment
- Transcription factor analysis
- Disease association analysis

**Use both PaintOmics AND other tools**:
- Best practice: Use multiple complementary tools
- PaintOmics shows pathways visually
- Enrichr/DAVID provide broader functional context
- OmicsNet shows protein interactions

---

## Troubleshooting

### Issue: Many genes not mapped

**Possible causes**:
1. Gene symbol format mismatch
2. Mouse vs. human gene symbols
3. Outdated gene symbols

**Solutions**:
1. Check gene symbols match KEGG database
   - Visit: https://www.genome.jp/kegg/genes.html
   - Search for your genes manually
2. Try uploading with UniProt IDs instead
3. Use `paintomics_high_confidence.csv` (cleaner dataset)
4. Remove problematic genes (C4A/C4B, CRISP1/CRISP3) and retry

### Issue: No enriched pathways found

**Possible causes**:
1. Gene list too small
2. Genes don't share common pathways
3. Enrichment threshold too stringent

**Solutions**:
1. Lower p-value threshold to 0.1
2. Use main file (33 proteins) instead of specialized files
3. Try "Show all pathways with ≥1 mapped gene" option
4. Check if genes mapped successfully in first place

### Issue: Pathways look empty or sparse

**Explanation**:
- This is normal for proteomics data (limited coverage)
- Most pathways have 50-200 components
- Your 33 proteins won't cover entire pathways

**What to look for**:
- Focus on pathways with ≥3 of your proteins
- Even 2-3 proteins can give biological insights
- Look for clustering (multiple proteins in same pathway region)

### Issue: Colors don't show up

**Possible causes**:
1. Data uploaded as categorical instead of numeric
2. Column not selected as "value" column
3. Data not recognized as log2 fold change

**Solutions**:
1. Re-upload and ensure "Testosterone_vs_Vehicle" is selected as numeric value
2. Check "log-transformed" option if available
3. Try `paintomics_minimal_input.csv` (simpler format)

---

## Advanced Features

### Multi-Omics Integration

PaintOmics allows you to integrate multiple omics layers:

**Example workflow for future experiments**:
1. Upload proteomics (this file)
2. Upload transcriptomics (RNA-seq, if available)
3. Upload metabolomics (if available)
4. PaintOmics will show all three on same pathway diagram

**File format for multi-omics**:
- Create separate CSV for each omics type
- Use same gene/metabolite naming convention
- Upload sequentially in PaintOmics

### Time-Course Analysis

If you have multiple timepoints:

1. Use the template: `paintomics_timecourse_template.csv`
2. Modify to include your timepoints:
   ```
   Name, T0h, T6h, T24h, T48h
   FRZB, 0, 1.2, 3.5, 4.4
   ```
3. PaintOmics will show temporal changes in pathway activity

### Custom Pathways

PaintOmics allows you to:
- Create custom pathway diagrams
- Combine multiple KEGG pathways
- Annotate pathways with your own data

---

## Publication-Quality Figures

### Exporting Pathway Diagrams

1. **Image format**:
   - **PNG**: Good for presentations, PowerPoint
   - **SVG**: Best for publications (vector format, editable)
   - **PDF**: Good for supplementary materials

2. **Resolution**:
   - Use at least 300 DPI for publications
   - SVG is resolution-independent (recommended)

3. **Color scheme**:
   - Ensure red/green colorblind-friendly options if needed
   - Use blue-orange diverging scale instead of red-green

### Figure Legends (Example)

> **Figure X. KEGG pathway analysis of testosterone-regulated proteins.**
> Expression data from 33 quantified proteins (log2 fold change, testosterone vs. vehicle) were visualized using PaintOmics 4 on KEGG pathway diagrams. Red indicates upregulated proteins, blue indicates downregulated proteins. (A) Cytokine-cytokine receptor interaction pathway showing downregulation of chemokines CCL7, CCL2, and CXCL5. (B) TGF-beta signaling pathway showing complex modulation with both up- and downregulated components. (C) ECM-receptor interaction pathway showing upregulation of multiple proteoglycans.

---

## Data Interpretation Guide

### What do the pathways tell you?

#### Downregulated Inflammatory Pathways
**Observation**: CCL7, CCL2, CXCL5 all downregulated
**Pathways affected**: Cytokine-cytokine receptor interaction
**Interpretation**:
- Testosterone reduces pro-inflammatory chemokine expression
- May decrease immune cell recruitment
- Anti-inflammatory effect

**Next steps**:
- Validate with immunostaining for immune cells
- Check for reduced macrophage/neutrophil infiltration
- Measure inflammatory cytokines (IL-6, TNF-alpha)

#### Upregulated ECM Remodeling
**Observation**: VCAN, LUM, BGN upregulated
**Pathways affected**: ECM-receptor interaction, Focal adhesion
**Interpretation**:
- Active tissue remodeling
- Extracellular matrix reorganization
- Potential fibrotic response

**Next steps**:
- Histological analysis of ECM structure
- Fibrosis markers (collagen staining)
- Functional consequences on tissue stiffness

#### Mixed TGF-beta Signaling
**Observation**: TGFB1/3 down, AMH up
**Pathways affected**: TGF-beta signaling
**Interpretation**:
- Complex regulation of TGF-beta pathway
- Not simple activation or suppression
- Context-dependent signaling

**Next steps**:
- Measure TGF-beta target genes (SMAD signaling)
- Check downstream effects (EMT markers, proliferation)
- Time-course to understand dynamics

---

## Recommended Analysis Workflow

### Step-by-Step Analysis Plan

**Week 1: Initial Exploration**
1. Upload `paintomics_input_expression_data.csv` (main file)
2. Identify top 5 enriched pathways
3. Export pathway diagrams for each
4. Take notes on observed patterns

**Week 2: Targeted Analysis**
1. Upload specialized files:
   - `paintomics_cytokines.csv` for inflammatory pathways
   - `paintomics_growth_factors.csv` for growth signaling
   - `paintomics_ecm_proteins.csv` for ECM pathways
2. Compare results across different subsets
3. Look for consistent patterns

**Week 3: Integration and Validation**
1. Compare PaintOmics results with:
   - OmicsNet interaction networks
   - Enrichr/DAVID functional enrichment
   - Literature on testosterone effects
2. Identify key pathways for experimental validation
3. Plan follow-up experiments

---

## Citation

If you use PaintOmics 4 in your research, please cite:

> Hernández-de-Diego R, Tarazona S, Martínez-Mira C, Balzano-Nogueira L, Furió-Tarí P, Pappas GJ Jr, Conesa A. PaintOmics 3: a web resource for the pathway analysis and visualization of multi-omics data. *Nucleic Acids Research*, 2018, 46(W1):W503-W509.

**Note**: Check the PaintOmics website for the most current citation.

---

## Additional Resources

### PaintOmics Help
- **User Manual**: http://www.paintomics.org/tutorial/
- **Video Tutorials**: http://www.paintomics.org/videos/
- **FAQs**: http://www.paintomics.org/faq/

### KEGG Database
- **KEGG Pathways**: https://www.genome.jp/kegg/pathway.html
- **Mouse pathways**: https://www.genome.jp/kegg-bin/show_organism?org=mmu
- **Gene search**: https://www.genome.jp/kegg/genes.html

### Related Tools
- **STRING** (protein interactions): https://string-db.org/
- **OmicsNet** (network analysis): https://www.omicsnet.ca/
- **Enrichr** (functional enrichment): https://maayanlab.cloud/Enrichr/
- **DAVID** (functional annotation): https://david.ncifcrf.gov/

---

## File Locations

All PaintOmics input files are located in:
```
Analysis/PaintOmics/
```

Generated from:
```
Analysis/General Analysis/cleaned_proteomics_data_with_QC_flags.csv
```

Using script:
```
Analysis/PaintOmics/generate_paintomics_files.py
```

---

## Quick Reference Card

### File Selection Guide

| Your Goal | Use This File | Why |
|-----------|--------------|-----|
| Comprehensive pathway view | `paintomics_input_expression_data.csv` | Most complete, 33 proteins |
| Quick test/validation | `paintomics_minimal_input.csv` | Simplest format |
| Ultra-clean analysis | `paintomics_high_confidence.csv` | Only highest quality data |
| Inflammatory pathways | `paintomics_cytokines.csv` | Focused on chemokines |
| Growth signaling | `paintomics_growth_factors.csv` | TGF-beta, PI3K-Akt pathways |
| ECM/tissue structure | `paintomics_ecm_proteins.csv` | ECM remodeling focus |

### Key Settings Checklist

- ✅ Organism: Mouse (Mus musculus)
- ✅ ID type: Gene Symbol
- ✅ Has header: YES
- ✅ Value column: Testosterone_vs_Vehicle
- ✅ Data type: Log2 fold change
- ✅ P-value threshold: 0.05 (or 0.1 for small datasets)
- ✅ Minimum genes per pathway: 2-3

---

**Last updated**: 2025-12-16
**Generated for**: Testosterone vs. Vehicle SILAC proteomics
**PaintOmics version**: 4 (check website for latest version)
