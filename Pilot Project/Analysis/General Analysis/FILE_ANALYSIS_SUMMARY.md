# General Analysis Folder - Comprehensive File Analysis

## Overview
This folder contains the complete data processing pipeline for SILAC proteomics analysis comparing **Testosterone vs Vehicle Control** conditions. The analysis processes 42 proteins with 2 biological replicates per condition.

---

## 📁 File Structure & Contents

### 🔧 **Python Scripts (Processing Pipeline)**

#### 1. `step1_data_cleaning.py` (285 lines)
**Purpose:** Initial data cleaning and formatting pipeline

**Key Functions:**
- Loads raw Excel data from `../CLEANED Silac Proteomics Soluble Factors.xlsx`
- Cleans numeric columns (removes commas, handles special values like `#DIV/0!`)
- Cleans text columns (strips whitespace, standardizes missing values)
- Standardizes gene symbols to uppercase
- Validates UniProt ID format (standard: `[PQ]#####`, extended: `[OPQ]#XXX#`)
- Calculates summary statistics (total proteins, upregulated/downregulated counts)
- Creates functional class distribution analysis

**Output Files Generated:**
- `cleaned_proteomics_data.csv` (main cleaned dataset)
- `summary_statistics.csv`
- `upregulated_proteins.csv`
- `downregulated_proteins.csv`
- `top_changers.csv`
- `functional_class_distribution.csv`

**Key Statistics:**
- Total proteins: 42
- Upregulated (FC > 1): 21
- Downregulated (FC < 1): 17
- Unchanged (FC ≈ 1): 1

---

#### 2. `step1_2_data_quality_investigation.py` (483 lines)
**Purpose:** Comprehensive data quality assessment and QC flagging

**Key Functions:**
- **Task 1:** Missing replicate data analysis
  - Tracks data completeness for 4 replicates (2 Vehicle, 2 Testosterone)
  - Categorizes: Complete (4/4), Missing 1 (3/4), Missing 2+ (2/4 or less)
  
- **Task 2:** Zero fold change investigation
  - Identifies proteins with FC = 0
  - Diagnoses causes: Complete suppression, No detection, etc.
  
- **Task 3:** UniProt ID validation
  - Validates standard and extended UniProt formats
  - Flags non-standard IDs
  
- **Task 4:** Overall data quality assessment
  - Missing data statistics per column
  - Overall completeness percentage (94.35%)
  
- **Task 5:** Quality flag creation
  - `Replicate_Completeness`: Complete, Missing_1, Missing_2+
  - `FC_Type`: Normal, Complete_Suppression, Complete_Induction, Cannot_Calculate, No_Detection
  - `Confidence_Level`: High, Medium, Low

**Output Files Generated:**
- `cleaned_proteomics_data_with_QC_flags.csv` (main dataset with QC flags)
- `missing_data_summary.csv`
- `proteins_flagged_for_review.csv`
- `high_confidence_proteins.csv`
- `data_quality_report.txt`

**Key Findings:**
- High confidence proteins: 25 (59.5%)
- Medium confidence: 8 (19.0%)
- Low confidence: 9 (21.4%)
- Complete data (4/4 replicates): 25 proteins (59.5%)

---

### 📊 **Data Files (CSV)**

#### 3. `cleaned_proteomics_data_with_QC_flags.csv` (43 rows, 27 columns)
**Purpose:** Main processed dataset with all QC flags

**Columns:**
- **Identifiers:** UniProt_Accession, Gene, Protein_Description, UniProt_ID, Alternate_Names, Mouse_Gene_ID
- **Location/Function:** Cellular_Location, Functional_Class
- **Replicate Data:** Vehicle_Rep1_31579, Vehicle_Rep2_31581, Testosterone_Rep1_31578, Testosterone_Rep2_31580
- **Statistics:** Vehicle_Mean, Vehicle_SD, Vehicle_SD_Percent, Testosterone_Mean, Testosterone_SD, Testosterone_SD_Percent
- **Fold Changes:** log_10_fold_change, log_2_fold_change, Fold_Change
- **QC Flags:** UniProt_Valid, UniProt_Standard_Format, UniProt_Extended_Format, Replicate_Completeness, FC_Type, Confidence_Level

**Key Features:**
- All 42 proteins with complete metadata
- QC flags for data quality assessment
- Ready for downstream analysis

---

#### 4. `high_confidence_proteins.csv` (27 rows)
**Purpose:** Filtered dataset containing only high-confidence proteins

**Criteria for High Confidence:**
- Complete replicate data (all 4 replicates present)
- Normal fold change calculations (not edge cases)
- No data quality concerns

**Contents:**
- 25 high-confidence proteins (59.5% of total)
- All columns from main dataset
- Recommended for primary statistical analyses

**Top Proteins:**
- **Upregulated:** FRZB (21.4×), GAS6 (8.2×), C4A/C4B (6.7×)
- **Downregulated:** CCL7 (0.14×), CCL2 (0.40×), TGFB3 (0.32×)

---

#### 5. `summary_statistics.csv` (11 rows)
**Purpose:** Overall dataset statistics

**Metrics:**
- Total Proteins: 42
- Upregulated (FC > 1): 21
- Downregulated (FC < 1): 17
- Unchanged (FC ≈ 1): 1
- Fold Change - Mean: 2.23
- Fold Change - Median: 1.34
- Fold Change - Std Dev: 3.67
- Fold Change - Min: 0.0
- Fold Change - Max: 21.40

---

#### 6. `functional_class_distribution.csv` (7 rows)
**Purpose:** Distribution of proteins by functional class

**Functional Classes:**
1. **other** (22 proteins, 52.4%)
   - Upregulated: 16, Downregulated: 4
2. **growth factor** (9 proteins, 21.4%)
   - Upregulated: 4, Downregulated: 5
3. **enzyme** (5 proteins, 11.9%)
   - Upregulated: 1, Downregulated: 2
4. **cytokine** (4 proteins, 9.5%)
   - Upregulated: 0, Downregulated: 4
5. **peptidase** (2 proteins, 4.8%)
   - Upregulated: 0, Downregulated: 2

**Key Insight:** Cytokines are exclusively downregulated by testosterone treatment.

---

#### 7. `upregulated_proteins.csv` (23 rows)
**Purpose:** List of all upregulated proteins (FC > 1)

**Contents:**
- 21 upregulated proteins
- Sorted by Fold_Change (descending)
- Includes all metadata and QC flags

**Top 5 Upregulated:**
1. FRZB: 21.4× (log₂FC = 4.42)
2. GAS6: 8.2× (log₂FC = 3.04)
3. C4A/C4B: 6.7× (log₂FC = 2.73)
4. CST12: 5.1× (log₂FC = 2.35)
5. VCAN: 3.5× (log₂FC = 1.81)

---

#### 8. `downregulated_proteins.csv` (19 rows)
**Purpose:** List of all downregulated proteins (FC < 1)

**Contents:**
- 17 downregulated proteins
- Sorted by Fold_Change (ascending)
- Includes complete suppression cases (FC = 0)

**Top 5 Downregulated:**
1. EGF: 0.001× (Complete suppression)
2. EREG: 0.001× (Complete suppression)
3. PENK: 0.001× (Complete suppression)
4. CXCL2: 0.001× (Complete suppression)
5. CCL7: 0.14× (log₂FC = -2.86)

---

#### 9. `top_changers.csv` (22 rows)
**Purpose:** Combined list of top 10 upregulated + top 10 downregulated

**Contents:**
- Top 10 upregulated proteins
- Top 10 downregulated proteins
- Columns: Gene, Protein_Description, Fold_Change, log_2_fold_change, Functional_Class, Confidence_Level, Replicate_Completeness, FC_Type, Regulation

**Use Case:** Quick reference for most significant changes

---

#### 10. `missing_data_summary.csv` (44 rows)
**Purpose:** Detailed missing data analysis per protein

**Columns:**
- Gene
- Vehicle_Rep1, Vehicle_Rep2, Testosterone_Rep1, Testosterone_Rep2 (Present/Missing status)
- Data_Completeness_Status (Complete, Missing 1, Missing 2, Missing 3, Missing All)
- Missing_Count

**Key Statistics:**
- Complete (4/4): 25 proteins
- Missing 1 (3/4): 8 proteins
- Missing 2+ (2/4 or less): 9 proteins
- Missing All (0/4): 4 proteins (LYZ, SERPINH1, NAXE, CYRIB)

---

#### 11. `proteins_flagged_for_review.csv` (19 rows)
**Purpose:** Proteins with data quality issues requiring review

**Columns:**
- Gene, Protein_Description
- Issues (comma-separated list)
- Severity (High, Medium, Low)
- Recommendation (Exclude, Flag)
- Confidence_Level, FC_Type, Replicate_Completeness

**Severity Breakdown:**
- **High (5 proteins):** Missing 2+ replicates or cannot calculate FC → **Exclude**
- **Medium (8 proteins):** Missing 1 replicate → **Flag**
- **Low (0 proteins):** Minor issues → **Flag**

**Proteins to Exclude:**
- EGF, EREG, PENK, CXCL2 (Complete suppression, missing 2+ replicates)
- AOC1 (Missing 2+ replicates)
- LYZ, SERPINH1, NAXE, CYRIB (Cannot calculate FC, missing all data)

---

### 📄 **Report Files**

#### 12. `data_quality_report.txt` (123 lines)
**Purpose:** Comprehensive text report of data quality assessment

**Sections:**
1. **Summary:** Overall statistics and completeness
2. **Task 1:** Replicate data completeness analysis
3. **Task 2:** Zero fold change proteins investigation
4. **Task 3:** UniProt ID validation results
5. **Task 4:** Missing data statistics by column
6. **Task 5:** Quality flag distributions
7. **Recommendations:** Usage guidelines for each confidence level

**Key Recommendations:**
- **USE FOR ALL ANALYSES:** 25 high-confidence proteins
- **USE WITH CAUTION:** 8 medium-confidence proteins
- **EXCLUDE FROM ANALYSIS:** 9 low-confidence proteins

---

## 🔄 **Data Processing Workflow**

```
1. step1_data_cleaning.py
   ↓
   Raw Excel → Cleaned CSV + Summary Files
   
2. step1_2_data_quality_investigation.py
   ↓
   Cleaned CSV → QC-Flagged Dataset + Quality Reports
   
3. Output Files Ready for:
   - Statistical analysis
   - Visualization
   - Publication
```

---

## 📈 **Key Dataset Characteristics**

### **Sample Information:**
- **Total Proteins:** 42
- **Conditions:** Vehicle Control vs Testosterone
- **Replicates:** 2 biological replicates per condition (4 total)
- **Data Completeness:** 94.35% overall

### **Quality Distribution:**
- **High Confidence:** 25 proteins (59.5%) ✅
- **Medium Confidence:** 8 proteins (19.0%) ⚠️
- **Low Confidence:** 9 proteins (21.4%) ❌

### **Regulation Summary:**
- **Upregulated:** 21 proteins (50%)
- **Downregulated:** 17 proteins (40.5%)
- **Unchanged:** 1 protein (2.4%)
- **Complete Suppression:** 4 proteins (9.5%)

### **Functional Class Distribution:**
- **Other:** 22 proteins (52.4%)
- **Growth Factor:** 9 proteins (21.4%)
- **Enzyme:** 5 proteins (11.9%)
- **Cytokine:** 4 proteins (9.5%)
- **Peptidase:** 2 proteins (4.8%)

---

## 🎯 **Recommended Usage**

### **For Primary Analysis:**
- Use `high_confidence_proteins.csv` (25 proteins)
- Complete replicate data
- Reliable fold change calculations

### **For Exploratory Analysis:**
- Use `cleaned_proteomics_data_with_QC_flags.csv` (all 42 proteins)
- Filter by `Confidence_Level` as needed
- Review `proteins_flagged_for_review.csv` for issues

### **For Reporting:**
- Reference `data_quality_report.txt` for methodology
- Use `summary_statistics.csv` for overall metrics
- Use `functional_class_distribution.csv` for class breakdowns

---

## ⚠️ **Important Notes**

1. **Complete Suppression Cases:** 4 proteins (EGF, EREG, PENK, CXCL2) show complete suppression (Testosterone = 0). These are biologically valid but flagged as low confidence due to missing replicates.

2. **Missing All Data:** 4 proteins (LYZ, SERPINH1, NAXE, CYRIB) have no replicate data and should be excluded.

3. **UniProt IDs:** All 42 proteins have valid UniProt IDs (100% validation rate).

4. **Replicate Completeness:** 59.5% of proteins have complete data across all 4 replicates.

---

## 📝 **File Dependencies**

```
step1_data_cleaning.py
  → Requires: ../CLEANED Silac Proteomics Soluble Factors.xlsx
  → Generates: cleaned_proteomics_data.csv + 5 summary CSVs

step1_2_data_quality_investigation.py
  → Requires: cleaned_proteomics_data.csv
  → Generates: cleaned_proteomics_data_with_QC_flags.csv + 4 QC files
```

---

**Last Updated:** Analysis complete
**Dataset:** SILAC Proteomics - Testosterone vs Vehicle Control
**Total Files:** 12 (2 Python scripts + 10 data/report files)

