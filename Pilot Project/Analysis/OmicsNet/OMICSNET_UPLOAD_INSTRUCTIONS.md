# OmicsNet 2.0 Network Analysis Guide

## Overview

This guide provides step-by-step instructions for uploading your proteomics data to OmicsNet 2.0 for protein-protein interaction (PPI) network analysis.

**OmicsNet 2.0 URL**: https://www.omicsnet.ca/

---

## Files Created

### Primary Analysis Files

#### 1. **omicsnet_input_with_expression.txt** ⭐ RECOMMENDED
- **Proteins**: 33 (High + Medium confidence)
- **Format**: Tab-delimited, 3 columns (Gene, log2FC, "Protein")
- **Contains**: Expression data (log2 fold changes)
- **Best for**: Colored networks showing upregulation (red) and downregulation (blue)
- **Use this for**: Your main comprehensive network analysis

#### 2. **omicsnet_high_confidence_only.txt**
- **Proteins**: 25 (High confidence only)
- **Format**: Tab-delimited, 3 columns (Gene, log2FC, "Protein")
- **Contains**: Expression data (log2 fold changes)
- **Use this for**: Cleaner networks if the main file produces too complex/crowded networks

#### 3. **omicsnet_genelist_only.txt**
- **Proteins**: 33 genes
- **Format**: Plain text, one gene per line
- **Contains**: Gene symbols only (no expression data)
- **Use this for**: Quick network generation without expression coloring

---

### Specialized Analysis Files

#### 4. **omicsnet_cytokines.txt**
- **Proteins**: 3 cytokines (CCL7, CCL2, CXCL5)
- **Use this for**: Focused inflammatory signaling network analysis
- **Expected network**: Small, focused on chemokine signaling pathways

#### 5. **omicsnet_growth_factors.txt**
- **Proteins**: 7 growth factors (GAS6, TGFB3, GDF15, TGFB1, NENF, AMH, AGT)
- **Use this for**: Growth factor signaling network analysis
- **Expected network**: TGF-beta pathway, angiogenesis, growth signaling

#### 6. **omicsnet_signaling_proteins.txt**
- **Proteins**: 10 proteins (cytokines + growth factors combined)
- **Use this for**: Comprehensive signaling network analysis
- **Expected network**: Combined inflammatory and growth factor signaling

#### 7. **omicsnet_upregulated.txt**
- **Proteins**: 21 upregulated proteins
- **Log2FC range**: +0.58 to +4.42
- **Use this for**: Understanding what pathways are activated/increased
- **Top upregulated**: FRZB (4.42), CRISP1/CRISP3 (3.69), GAS6 (3.04)

#### 8. **omicsnet_downregulated.txt**
- **Proteins**: 12 downregulated proteins
- **Log2FC range**: -2.86 to -1.00
- **Use this for**: Understanding what pathways are suppressed/decreased
- **Top downregulated**: CCL7 (-2.86), CCL2 (-2.35), CXCL5 (-1.00)

---

## Step-by-Step Upload Instructions

### Step 1: Navigate to OmicsNet 2.0
1. Go to: **https://www.omicsnet.ca/**
2. Click **"Get Started"** or **"Upload Data"**

### Step 2: Select Input Type
Choose one of these options:
- **Option A**: "Upload data with expression values" (RECOMMENDED)
  - Use for files with expression data (.txt files with log2FC values)
  - Allows colored networks

- **Option B**: "Upload gene list only"
  - Use for `omicsnet_genelist_only.txt`
  - Simple network without expression coloring

### Step 3: Configure Upload Settings

#### For Expression Data Upload (Option A):
1. **Select organism**: Mouse (*Mus musculus*)
2. **ID type**: Gene Symbol
3. **Upload file**: Choose your .txt file (e.g., `omicsnet_input_with_expression.txt`)
4. **Data type**: Protein
5. **Expression format**: Log2 fold change (if asked)

#### For Gene List Upload (Option B):
1. **Select organism**: Mouse (*Mus musculus*)
2. **ID type**: Gene Symbol
3. **Upload file**: `omicsnet_genelist_only.txt`
4. **Data type**: Protein

### Step 4: Data Processing
- OmicsNet will map your gene symbols to protein IDs
- Review the mapping results
- Check for any unmapped genes
  - Common issue: `C4A/C4B` may need to be split into separate entries if mapping fails
  - Alternative names may be suggested

### Step 5: Network Construction Settings

**Recommended Settings:**
- **Network type**: Protein-Protein Interaction (PPI)
- **Confidence level**:
  - Start with **"High confidence" (700+)** for cleaner networks
  - If network is too sparse, lower to **"Medium confidence" (400+)**
- **Network degree**:
  - Start with **1st order** (direct interactions only)
  - Expand to **2nd order** if you want indirect interactions
- **Minimum network size**: Default (usually 5-10 nodes)
- **Maximum network size**:
  - For main file (33 proteins): Set to 500-1000 nodes
  - For focused files (3-10 proteins): Set to 200-500 nodes

### Step 6: Visualization Options

**Node coloring:**
- If you uploaded expression data, nodes will be colored:
  - **Red**: Upregulated proteins
  - **Blue**: Downregulated proteins
  - **Gray**: Interacting partners (not in your input list)

**Layout options:**
- **Force-directed layout**: Good for most networks
- **Circular layout**: Good for small, focused networks
- **Hierarchical layout**: Good for pathway visualization

### Step 7: Analysis Features to Explore

1. **Pathway Enrichment**
   - OmicsNet will automatically perform pathway enrichment
   - Look for KEGG pathways, GO terms, Reactome pathways
   - Expected enriched pathways:
     - Cytokine-cytokine receptor interaction
     - TGF-beta signaling pathway
     - ECM-receptor interaction
     - PI3K-Akt signaling pathway

2. **Hub Analysis**
   - Identify key hub proteins in your network
   - These are proteins with many interactions
   - May represent key regulatory nodes

3. **Module Detection**
   - OmicsNet can identify network modules/clusters
   - Useful for finding functional groups of proteins

4. **Export Options**
   - Export network as image (PNG, SVG)
   - Export network data (GraphML, SIF, JSON)
   - Export enrichment results (CSV, TSV)

---

## Recommended Analysis Strategy

### Strategy 1: Comprehensive Analysis
1. **Start with**: `omicsnet_input_with_expression.txt` (33 proteins)
2. **Goal**: Get overall view of protein interactions
3. **Settings**: 1st order network, high confidence
4. **Expected outcome**: Medium-sized network showing major interaction hubs

### Strategy 2: High-Quality Analysis
1. **Start with**: `omicsnet_high_confidence_only.txt` (25 proteins)
2. **Goal**: Most reliable interactions only
3. **Settings**: 1st order network, high confidence
4. **Expected outcome**: Cleaner, more focused network

### Strategy 3: Functional Focus
1. **Start with**: `omicsnet_signaling_proteins.txt` (10 proteins)
2. **Goal**: Understand cytokine and growth factor crosstalk
3. **Settings**: 1st or 2nd order network, medium-high confidence
4. **Expected outcome**: Signaling pathway network

### Strategy 4: Direction-Specific Analysis
1. **Upload separately**:
   - `omicsnet_upregulated.txt` (21 proteins)
   - `omicsnet_downregulated.txt` (12 proteins)
2. **Goal**: Understand activation vs. suppression networks
3. **Compare**: Are different pathways upregulated vs. downregulated?

---

## Troubleshooting

### Issue: Too many unmapped genes
**Solution**:
- Check if gene symbols are current
- Try alternative gene symbols from the `Alternate_Names` column
- Use UniProt IDs instead (from `UniProt_ID` column)

### Issue: Network too large/complex
**Solutions**:
1. Increase confidence threshold (700 → 900)
2. Use only 1st order interactions
3. Switch to `omicsnet_high_confidence_only.txt`
4. Use one of the focused files (cytokines, growth factors)

### Issue: Network too small/sparse
**Solutions**:
1. Lower confidence threshold (700 → 400)
2. Include 2nd order interactions
3. Use `omicsnet_input_with_expression.txt` instead of high-confidence only
4. Check if gene mapping was successful

### Issue: No pathway enrichment
**Possible reasons**:
- Network too small (need at least 5-10 mapped proteins)
- Proteins don't share common pathways
- Try expanding to 2nd order interactions

---

## Key Proteins to Watch For

Based on your data, these proteins are likely to be network hubs:

### Top Upregulated (Network Hubs):
1. **FRZB** (log2FC: +4.42) - Wnt signaling inhibitor
2. **CRISP1/CRISP3** (log2FC: +3.69) - Signaling proteins
3. **GAS6** (log2FC: +3.04) - Growth arrest-specific protein, TAM receptor ligand
4. **C4A/C4B** (log2FC: +2.73) - Complement cascade
5. **TGFB1, TGFB3** - TGF-beta signaling pathway

### Top Downregulated (Network Hubs):
1. **CCL7** (log2FC: -2.86) - Chemokine, inflammatory signaling
2. **CCL2** (log2FC: -2.35) - Chemokine, monocyte recruitment
3. **CXCL5** (log2FC: -1.00) - Chemokine, neutrophil recruitment

---

## Expected Biological Insights

Based on your dataset, you should expect to find:

1. **TGF-beta signaling network**
   - TGFB1, TGFB3, GDF15, TGFBI
   - Implications: Tissue remodeling, fibrosis, cell growth

2. **Chemokine/cytokine signaling**
   - CCL7, CCL2, CXCL5 (all downregulated)
   - Implications: Reduced inflammatory signaling

3. **Growth factor signaling**
   - GAS6, AGT, AMH, NENF
   - Implications: Cell survival, proliferation pathways

4. **Extracellular matrix (ECM) remodeling**
   - VCAN, LUM, BGN, SERPINE2
   - Implications: Tissue structure changes

5. **Complement cascade activation**
   - C4A/C4B upregulated
   - Implications: Immune system modulation

---

## Data Quality Notes

**Your data filtering:**
- ✅ Excluded proteins with `FC_Type = 'Cannot_Calculate'` (4 proteins removed)
- ✅ Excluded proteins with NaN fold changes
- ✅ Included `Complete_Suppression` proteins (0 values in treatment → very strong downregulation)
- ✅ Only High and Medium confidence proteins included

**Confidence levels:**
- **High confidence (25 proteins)**: Complete data in all replicates
- **Medium confidence (8 proteins)**: Missing 1 replicate but valid fold change
- **Total in main file**: 33 proteins

---

## Next Steps After OmicsNet Analysis

1. **Export your networks**
   - Save as high-resolution images for publications
   - Export network data for further analysis

2. **Validate key interactions**
   - Literature search for top hub proteins
   - Consider experimental validation of key interactions

3. **Integrate with pathway analysis**
   - Compare OmicsNet results with your previous enrichment analyses
   - Look for concordance between methods

4. **Generate publication figures**
   - OmicsNet networks make excellent figures for papers
   - Combine with your volcano plots and heatmaps

---

## Citation

If you use OmicsNet 2.0 in your research, please cite:

> Zhou G, Xia J. OmicsNet 2.0: a web-based platform for network-based multi-omics integration and interpretation. *Nucleic Acids Research*, 2022, 50(W1):W317-W324.

---

## Additional Resources

- **OmicsNet Tutorial**: https://www.omicsnet.ca/docs/
- **FAQ**: https://www.omicsnet.ca/docs/FAQs.xhtml
- **Video Tutorials**: https://www.omicsnet.ca/docs/Tutorial.xhtml

---

## File Location

All OmicsNet input files are located in:
```
Analysis/OmicsNet/
```

Generated from:
```
Analysis/General Analysis/cleaned_proteomics_data_with_QC_flags.csv
```

Using script:
```
Analysis/OmicsNet/generate_omicsnet_files.py
```

---

**Last updated**: 2025-12-16
**Generated for**: Testosterone vs. Vehicle SILAC proteomics comparison
