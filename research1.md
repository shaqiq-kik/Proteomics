# Technical Design Doc: Proteomics Differential Expression Pipeline (v2)

## TL;DR
- The existing `proteomics_analysis.py` is a working 88-line ETL script (SILAC Light vs. Heavy intensities, inner-join on UniProt Accession, ratio-based fold-change classification) but has three correctness bugs that will produce wrong "UP/DOWN" calls: (1) it averages raw ratios arithmetically then log-transforms, (2) its UP (≥1.5) and DOWN (≤0.5) thresholds are asymmetric on the log2 scale, and (3) it has no statistical testing. Build all new work in a separate `proteomics_de/` module and leave the legacy script untouched.
- The statistically-defensible core, per Peng, H., Wang, H., Kong, W., Li, J. & Goh, W. W. B. "Optimizing differential expression analysis for proteomics data via high-performing rules and ensemble inference." *Nat. Commun.* **15**, 3922 (2024) (DOI 10.1038/s41467-024-47899-w), is: build a protein intensity matrix → **no extra normalization** → **MinProb** imputation (`imputeLCMD`/`DEP`) → **limma** moderated t-test, invoked from Python as a subprocess over a clean TSV file contract.
- With n=2 and zero biological replicates, p-values are technical-noise p-values only; limma will run and is the best available option, but PCA, clustering, WGCNA, UMAP, and most ML are statistically degenerate at 4 samples and must be excluded. The architecture should be config-driven (sample sheet + parameterized design matrix) so the future addition of biological replicates is a config change, not a rewrite.

## Key Findings

1. **The repo is real and SILAC-based.** The repository at `shaqiq-kik/Proteomics` contains `proteomics_analysis.py` (88 lines, 3.17 KB), a `Pilot Project` folder, `Copy of General Sheet.xlsx`, `General Sheet.xlsx`, `IPA_input.csv`, `IPA_input.xlsx`, `Limma.pdf`, and `alligned_proteins.xlsx`. GitHub reports the repo as 60.7% R / 39.3% Python, confirming R work exists (consistent with `Limma.pdf` and the `Pilot Project` folder). The two input sheets are named **"Protein Report L"** and **"Protein Report H"** — Light and Heavy SILAC channels, i.e. vehicle control (L) vs. testosterone-treated (H).

2. **Three concrete correctness bugs** in the fold-change logic (averaging ratios before log; asymmetric 1.5/0.5 thresholds; inf/NaN propagation from division), plus a silent data-loss issue from the inner join dropping proteins detected in only one condition.

3. **Peng et al. 2024 gives a directly applicable recipe** for label-free proteomics. The paper's verbatim finding: "high-performing workflows prefer no normalization and incline MinProb (probabilistic minimum) for imputation while eschewing simple statistical tools (e.g., ANOVA, SAM and t-test) for DEA." Its general-case recommendation, verbatim: "For expression matrices without acquisition platform or quantification information, we recommend a workflow combining no normalization, loess, Rlr, center.median for normalization options, MinProb, SeqKNN, Impseq or MinDet for MVI, and limma or ROTS for DEA."

4. **QIAGEN IPA needs only an identifier column** to run; fold change and p-value are optional but unlock cutoffs and statistical power. IPA accepts UniProt accession directly. A flat tab/CSV file avoids the Excel header-recognition problem.

5. **n=2 is a hard architectural constraint.** It is below the documented minimums for WGCNA. The official WGCNA FAQ (Langfelder & Horvath) states verbatim: "We do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples. In a typical high-throughput setting, correlations on fewer than 15 samples will simply be too noisy for the network to be biologically meaningful. If at all possible, one should have at least 20 samples." n=2 also makes PCA on 4 samples rank-deficient and clustering unstable. These analyses must be gated behind a replicate-count check.

## Details

### Audit of the existing code (`proteomics_analysis.py`)

The script does the following, in order (referencing the actual variable names):

```python
INPUT_FILE = "Copy of General Sheet.xlsx"
df_L = pd.read_excel(INPUT_FILE, sheet_name="Protein Report L")   # control (Light)
df_H = pd.read_excel(INPUT_FILE, sheet_name="Protein Report H")   # testosterone (Heavy)
cols_L = ["UniProt Accession Number","Gene names","Intensity 31578","Intensity 31580"]
cols_H = ["UniProt Accession Number","Gene names","Intensity 31579","Intensity 31581"]
df_merged = pd.merge(df_L, df_H, on="UniProt Accession Number", how="inner", suffixes=("_L","_H"))
df_merged["Gene names"] = df_merged["Gene names_L"].combine_first(df_merged["Gene names_H"])
df_merged["has_missing_data"] = (df_merged[intensity_cols].isnull().any(axis=1) | (df_merged[intensity_cols]==0).any(axis=1))
df_merged["fold_change_exp1"] = df_merged["Intensity 31579"] / df_merged["Intensity 31578"]
df_merged["fold_change_exp2"] = df_merged["Intensity 31581"] / df_merged["Intensity 31580"]
df_merged["mean_fold_change"] = (df_merged["fold_change_exp1"] + df_merged["fold_change_exp2"]) / 2
df_merged["log2_fold_change"] = np.log2(df_merged["mean_fold_change"])
df_merged.loc[complete_data & (df_merged["mean_fold_change"]>=1.5), "regulated"] = "UP"
df_merged.loc[complete_data & (df_merged["mean_fold_change"]<=0.5), "regulated"] = "DOWN"
df_ipa = df_merged[(~df_merged["has_missing_data"]) & (df_merged["regulated"]!="NO CHANGE")].copy()
```

So run 1 = (Intensity 31578 Light + Intensity 31579 Heavy); run 2 = (31580 Light + 31581 Heavy). `fold_change_exp1` and `fold_change_exp2` are the two H/L SILAC ratios.

**Bug 1 — Arithmetic mean of ratios then log (statistical bias).** `mean_fold_change = (fc1 + fc2)/2` averages raw ratios. Ratios are multiplicative; the correct central tendency is the geometric mean (equivalently, the mean of the log2 ratios). Worked example: a protein up 2× in run 1 (ratio 2.0) and down 2× in run 2 (ratio 0.5) yields arithmetic mean 1.25 → log2 = +0.32 (looks "up"), but the geometric mean is 1.0 → log2 = 0 (no change). The current code systematically inflates apparent up-regulation. **Fix:** compute `log2(fc1)` and `log2(fc2)` first, then average.

**Bug 2 — Asymmetric thresholds on the log2 scale.** UP is `mean_fold_change >= 1.5` (log2 = +0.585) while DOWN is `<= 0.5` (log2 = −1.0). A protein must rise only 1.5× to be UP but must fall a full 2× to be DOWN. The symmetric down-threshold for a 1.5× cut is `1/1.5 = 0.667` (log2 = −0.585). As written, the classifier is biased against detecting down-regulation. **Fix:** use symmetric log2 cutoffs (e.g. |log2FC| ≥ 0.585).

**Bug 3 — inf/NaN propagation from division and log.** Fold changes are computed for every row before the missing/zero filter is applied, so rows with a zero Light denominator produce `inf`, and `np.log2(0)` produces `−inf`. These survive in the saved `alligned_proteins.xlsx` columns. They are correctly excluded from UP/DOWN (gated by `complete_data`), but the persisted artifact contains non-finite values that will break any downstream consumer that does not re-filter. **Fix:** mask non-finite values and compute fold changes only on complete rows.

**Bug 4 — Inner join silently drops the most interesting proteins.** `how="inner"` keeps only the ~1,948 accessions present in both sheets. Proteins detected only in control or only in treated ("on/off" proteins — often the strongest biological signal) are discarded with no log. Additional join hazards on `UniProt Accession Number`: isoform suffixes (`P12345-2` vs `P12345`), multi-protein groups (`P12345;Q67890` semicolon lists), and duplicate accessions causing many-to-many row explosion. **Fix:** validate accession cardinality before the merge; record dropped exclusives to a side file; decide isoform collapse policy explicitly.

**Bug 5 — No normalization sanity check.** SILAC log-ratios should be centered on zero; the code never checks this. **Bug 6 — No replicate correlation check** between `fold_change_exp1` and `fold_change_exp2`. **Bug 7 — No statistical testing**, so no p-values. Minor: misspelled artifacts `alligned_proteins.xlsx` and the `"Rown in new IPA"` print string.

**Validated-correct:** the IPA export filter `(~has_missing_data) & (regulated != "NO CHANGE")` is now correct, addressing the earlier row-count bug where the IPA file had contained all proteins. The presence of both `IPA_input.csv` and `IPA_input.xlsx` is consistent with the documented fix of switching to CSV after IPA failed to read Excel headers.

---

### SECTION 1 — Statistical Testing Module (limma, per Peng et al. 2024)

**Why limma, not a t-test.** With n=2 a per-protein t-test has 2 degrees of freedom split across groups and essentially no power; the variance estimate is wildly unstable. limma's empirical-Bayes step shrinks each protein's variance toward a global trend estimated across all ~1,948 proteins, which is exactly what rescues inference at small n. Peng et al. 2024 found limma and ROTS "consistently rank amongst the top 4 options across all 6 settings," while plain t-test/ANOVA/SAM were enriched in low-performing workflows. Per the paper's label-free recommendation: **directLFQ intensity matrix, no normalization, MinProb imputation, limma DEA.**

**How MinProb works (plain English).** Missing values in MS proteomics are mostly missing-not-at-random (MNAR): a protein is absent because it was below the detection limit, not by random chance. MinProb (`imputeLCMD::impute.MinProb`, default `q = 0.01, tune.sigma = 1`) fills each missing value with a random draw from a Gaussian centered near the low end of each sample's observed-intensity distribution (the q-th quantile), with width set from the median per-protein standard deviation. In effect it says "this protein was probably present at a low, near-detection-limit level," which is the biologically correct assumption for left-censored data. Implement via `DEP::impute(..., fun="MinProb", q=0.01)` (a thin wrapper over `imputeLCMD`) on log2-transformed data.

**Design matrix for a two-condition comparison.** Treat the four intensity columns as two groups: control = {31578, 31580}, treated = {31579, 31581}. `group <- factor(c("control","control","treated","treated"), levels=c("control","treated"))`; `design <- model.matrix(~group)`; the coefficient `grouptreated` is the log2 fold change (treated − control). `lmFit → eBayes → topTable(coef="grouptreated", number=Inf, adjust.method="BH")`.

**The runnable R script** (`run_limma.R`):

```r
#!/usr/bin/env Rscript
# run_limma.R — two-condition limma DE for SILAC intensity data.
# Per Peng et al. 2024 (Nat Commun 15:3922; DOI 10.1038/s41467-024-47899-w):
#   no extra normalization; MinProb imputation; limma moderated t-test.
# File contract (CLI args):
#   --in  <intensity_matrix.tsv>  rows=proteins, cols=accession,gene + 4 intensity cols
#   --design <design.tsv>         columns: sample, group   (group in {control,treated})
#   --out <results.tsv>           limma topTable + provenance
suppressPackageStartupMessages({
  library(optparse); library(limma); library(DEP)
  library(SummarizedExperiment); library(matrixStats)
})

opt <- parse_args(OptionParser(option_list=list(
  make_option("--in",     dest="infile"),
  make_option("--design", dest="design"),
  make_option("--out",    dest="outfile"),
  make_option("--minprob_q", default=0.01),
  make_option("--seed",   default=42L)
)))
set.seed(opt$seed)   # MinProb is stochastic — seed for reproducibility

# 1) READ -------------------------------------------------------------
mat_df <- read.delim(opt$infile, check.names=FALSE, stringsAsFactors=FALSE)
design_df <- read.delim(opt$design, stringsAsFactors=FALSE)
stopifnot(c("accession","gene") %in% colnames(mat_df))
intensity_cols <- design_df$sample
stopifnot(all(intensity_cols %in% colnames(mat_df)))

# 2) MATRIX + LOG2 (zeros -> NA so MinProb treats them as MNAR) --------
M <- as.matrix(mat_df[, intensity_cols]); M[M == 0] <- NA
logM <- log2(M)

# 3) NO NORMALIZATION (Peng et al. 2024 high-performing rule) ----------
#    Intentionally skipped. Diagnostic only:
boxplot(logM, main="log2 intensity (no normalization)", las=2)

# 4) MinProb IMPUTATION via DEP (wraps imputeLCMD::impute.MinProb) -----
rowData_df <- data.frame(name=mat_df$accession, ID=mat_df$accession,
                         gene=mat_df$gene, stringsAsFactors=FALSE)
se <- SummarizedExperiment(assays=list(logM),
        rowData=rowData_df,
        colData=DataFrame(label=intensity_cols,
                          condition=design_df$group,
                          replicate=ave(design_df$group, design_df$group,
                                        FUN=seq_along)))
se_imp <- DEP::impute(se, fun="MinProb", q=opt$minprob_q)
logM_imp <- assay(se_imp)

# 5) DESIGN MATRIX + limma --------------------------------------------
grp <- factor(design_df$group, levels=c("control","treated"))
design <- model.matrix(~grp)            # coef 2 = treated - control (log2FC)
fit <- lmFit(logM_imp, design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
tt  <- topTable(fit, coef=2, number=Inf, adjust.method="BH", sort.by="none")

# 6) ASSEMBLE OUTPUT (stable schema) ----------------------------------
out <- data.frame(
  accession   = mat_df$accession,
  gene        = mat_df$gene,
  logFC       = tt$logFC,
  fold_change = 2^tt$logFC,
  AveExpr     = tt$AveExpr,
  t           = tt$t,
  P.Value     = tt$P.Value,
  adj.P.Val   = tt$adj.P.Val,
  B           = tt$B,
  n_imputed   = rowSums(is.na(logM)),
  stringsAsFactors=FALSE)
out$direction <- ifelse(out$adj.P.Val < 0.05 & out$logFC >=  0.585, "UP",
                 ifelse(out$adj.P.Val < 0.05 & out$logFC <= -0.585, "DOWN", "NS"))

write.table(out, opt$outfile, sep="\t", quote=FALSE, row.names=FALSE)
writeLines(c(
  paste0("# limma_version: ", as.character(packageVersion("limma"))),
  paste0("# minprob_q: ", opt$minprob_q),
  paste0("# n_control: ", sum(grp=="control"), " n_treated: ", sum(grp=="treated")),
  paste0("# generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
), con=paste0(opt$outfile, ".meta"))
```

**Integration / subprocess / file contract.** Python orchestrates; R computes. Python writes a deterministic TSV matrix and a design TSV, calls `Rscript run_limma.R ...` via `subprocess.run([...], check=True, capture_output=True)`, checks the return code, then reads the results TSV back into pandas. The contract is three flat files with fixed schemas (below). This keeps the languages decoupled: neither `rpy2` nor a shared process, just files and exit codes — easy to test, cache, and reproduce.

```python
# de/run_limma.py (sketch)
import subprocess, pandas as pd, pathlib
def run_limma(matrix_tsv, design_tsv, out_tsv, rscript="run_limma.R"):
    r = subprocess.run(["Rscript", rscript, "--in", matrix_tsv,
                        "--design", design_tsv, "--out", out_tsv],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"limma failed:\n{r.stderr}")
    return pd.read_csv(out_tsv, sep="\t")
```

**Output schemas / naming conventions.**
- `intensity_matrix.tsv`: `accession, gene, <sample_1..n>` (one row per protein).
- `design.tsv`: `sample, group` (group ∈ {control, treated}).
- `limma_results.tsv`: `accession, gene, logFC, fold_change, AveExpr, t, P.Value, adj.P.Val, B, n_imputed, direction`.
- Naming: `de_<contrast>_<YYYYMMDD>.tsv`, e.g. `de_testosterone_vs_vehicle_20260611.tsv`; sidecar `.meta` carries package versions and parameters for reproducibility.

**Alternatives.** ROTS (`ROTS` Bioconductor) is the paper's co-recommended method and a good cross-check; it optimizes a modified t-statistic for reproducibility and estimates FDR by permutation, but permutation-based FDR is itself unreliable at n=2. DEqMS adds PSM-count-aware variance but requires a peptide/PSM count column not present here. For a pure SILAC-ratio formulation, `Proteus::limmaRatioDE` runs a one-sample limma test against the null log-ratio = 0; this is a valid secondary analysis but the two-group intensity design above matches Peng et al. more directly.

---

### SECTION 2 — Data Quality & Validation System

Build a dedicated `qc/` module. Use **pandera** for schema/range validation (lightweight, ~12 dependencies, DataFrameModel API, integrates as decorators on pipeline functions) and **pytest** for unit tests; reserve **great_expectations** for later if a data-warehouse-style profiling layer is wanted (it pulls in 100+ dependencies and is overkill here). Validate at every boundary: after load, after merge, after fold-change, before IPA export.

Concrete checks:
- **Schema validation (pandera):** required columns present and typed; intensity columns numeric; `accession` matches UniProt regex `^[A-NR-Z0-9]([A-Z0-9]{5}|[A-Z0-9]{9})(-\d+)?$`; `accession` unique (`Field(unique=True)`); `strict=True` to catch upstream schema drift.
- **Range checks:** intensities ≥ 0; flag negatives; fold changes finite; log2FC within a sane bound.
- **Duplicate / cardinality detection:** count duplicate accessions pre-merge; assert merge does not increase row count beyond `min(len(L),len(H))` (guards against many-to-many explosion); detect and strip isoform suffixes per policy.
- **Missingness audit:** per-sample and per-protein missing/zero counts; a missing-value heatmap; assert imputation only touches values flagged missing.
- **Replicate correlation QC:** Pearson/Spearman of `log2(fc_exp1)` vs `log2(fc_exp2)` (and raw log2 intensities run1 vs run2); fail/warn below a threshold (e.g. r < 0.7) since at n=2 the two replicates are the entire basis of the variance estimate.
- **Distribution checks:** SILAC log-ratios centered near zero (assert |median(log2FC)| < 0.2, else normalization may be needed despite the "no normalization" rule); intensity distributions comparable across samples (boxplot).
- **Export validation:** assert IPA file row count == count of (complete & UP|DOWN); assert no NA/inf in exported fold-change/p-value columns; assert identifier column is leftmost and free of stray characters.

Each check emits a structured pass/fail record to a `qc_report.json` + human-readable `qc_report.html`, and the pipeline halts (or quarantines invalid rows with logging) on hard failures.

---

### SECTION 3 — QIAGEN IPA Integration (ETL contract)

Treat IPA as an external sink with a strict file contract. From QIAGEN's official documentation ("Data Upload definitions"), verbatim: "the dataset should start with at least one column of molecular identifiers (for genes, proteins, or chemicals) and then ideally have additional columns for 'measurements' such as fold changes and p-values. However, you can upload a list of molecules without [any measurement values]." So the **identifier column is the only strictly required element**; fold change, p-value, and FDR are optional. IPA accepts a wide range of IDs **including UniProt accession and gene symbol** directly (it cannot accept metabolite names). From QIAGEN's "Formatting Your Dataset" KB, verbatim: "Your file may contain up to 20 observations (experimental conditions) with up to eight measurement value types per observation. Measurement value types are ratio, fold change, log ratio, p-value, FDR, Intensity, etc." Accepted formats: `.xlsx, .xls, .txt (tab-delimited), .csv, .diff`; from QIAGEN's "File Format Specification" KB, verbatim: "We strongly recommend that the dataset files be in tab-delimited text format (.txt, rather than Excel-based formats) for faster upload." For Excel uploads the data must be on the **first worksheet**. This is the documented basis for the Excel→CSV fix: switching to a flat CSV resolved IPA's failure to recognize the column headers.

**What the pipeline must produce for IPA:**
- File `ipa_<contrast>_<date>.csv`: columns `UniProtAccession, GeneSymbol, ExprLogRatio (or FoldChange), pValue, FDR`. Identifier column leftmost. UTF-8, no BOM, single header row, no merged cells.
- Use **log2 fold change** mapped to IPA's "Expr Log Ratio" column type (cleanest, sign-symmetric), plus the limma `adj.P.Val` as FDR.

**Which IPA modules need what (run order):**
1. **Core Analysis** — the parent analysis; requires identifiers + a measurement (fold change/expression value). Run first; everything else is computed from it. p-value optional but enables expression cutoffs (Standard/Bonferroni/FDR).
2. **Canonical Pathways** — enrichment via right-tailed Fisher's Exact Test; IPA computes its own overlap p-value. Needs the molecule set (fold change to color/direction). z-scores for pathway activation require directional (red/green) fold-change overlay. (Already populated in the existing analysis.)
3. **Upstream Regulator Analysis** — uses the **direction** of fold change of target molecules to compute activation z-scores; overlap p-value computed internally. Fold change sufficient; p-value improves the molecule set quality.
4. **Causal Network Analysis** — extends upstream regulators into multi-tier networks scored by activation z-score; same directional fold-change requirement.
5. **Downstream Effects / Diseases & Functions** — predicts increased/decreased biological functions from fold-change direction; z-score based.
6. **Regulator Effects** — links upstream regulators to downstream functions via a Consistency Score; needs the Core Analysis with directional fold changes.
7. **Comparison Analysis** — only meaningful across ≥2 Core Analyses (e.g. multiple contrasts or dose levels); not applicable to a single contrast now.

**Does p-value unlock features?** Not for *running* the analytics (direction of fold change drives all z-score predictions; enrichment p-values are computed by IPA). It unlocks **cutoff-based filtering** (only push significant molecules) and improves the reliability/defensibility of the analyzed set. So: ship fold change + adj.P.Val once limma is in place.

---

### SECTION 4 — Full Analysis & Visualization Engineering Plan

For ~1,948 proteins with fold changes (now) and p-values (after limma). "P?" = requires p-values.

**Buildable now / appropriate:**
- **Volcano plot** — *P? Yes.* x=log2FC, y=−log10(adj.P.Val). Python `matplotlib`/`seaborn` or R `EnhancedVolcano`/`ggplot2`. Input: `limma_results.tsv`.
- **MA plot** — *P? No (color by sig if available).* x=AveExpr, y=log2FC. `matplotlib` / limma `plotMD`.
- **Heatmap / clustermap** — *P? No.* log2 intensities of top-N variable or DE proteins; R `pheatmap` or Python `seaborn.clustermap`. Note: only 4 columns; cluster proteins (rows), not samples.
- **Correlation matrices** — *P? No.* replicate–replicate and sample–sample Pearson/Spearman; `seaborn.heatmap` / `corrplot`. Critical QC at n=2.
- **Intensity distribution plots** — *P? No.* boxplots / KDE density per sample; `seaborn`. Pre/post imputation overlay (`DEP::plot_imputation`).
- **Missing-value heatmap** — *P? No.* binary present/absent matrix; `DEP::plot_missval` / `seaborn`.
- **Rank-abundance plot** — *P? No.* sorted log2 intensity vs rank; `matplotlib`.
- **STRING PPI network** — *P? No (use DE list to seed).* Map DE accessions via STRING REST API (`/api/tsv/get_string_ids` then `/api/tsv/network`), species = 9606 (human), `required_score` ≥ 400/700; or R `STRINGdb` (caches network as igraph) / `rbioapi`. Respect rate limit (≥1 s between calls).
- **Network construction & visualization** — *P? No.* `igraph`/`networkx` for graph metrics (degree, betweenness, clusters); `py4cytoscape` (CyREST, 250+ functions, RdBu node coloring by log2FC) or Cytoscape for publication figures; map log2FC onto node color.
- **GO / pathway over-representation (ORA)** — *P? Yes, to define the DE gene set.* `clusterProfiler` (`enrichGO`, `enrichKEGG`), `ReactomePA`, `gprofiler2` (g:Profiler REST), `enrichR` (Enrichr API), `g:Profiler`/`PANTHER`/`DAVID`/`Reactome` web APIs. Input: DE accession/symbol list + background of all 1,948 detected.
- **GSEA (rank-based)** — *P? Yes (needs ranked metric).* rank all proteins by `sign(logFC)*-log10(P.Value)`; `fgsea`, `clusterProfiler::gseGO/gseKEGG`, `ReactomePA::gsePathway`. Uses all proteins, not a thresholded list.
- **UpSet plots** — *P? No.* intersections of DE sets / pathway memberships; Python `upsetplot`, R `UpSetR`, or `clusterProfiler::upsetplot`.
- **PCA** — *P? No, but see Section 5.* technically runs on 4 samples but is degenerate; build it gated and clearly captioned as QC-only.

**Inappropriate for ~1,948 proteins × n=2 — and why:**
- **UMAP / t-SNE** — designed for many samples/cells; 4 points yield meaningless embeddings. Exclude.
- **WGCNA** — official WGCNA FAQ (Langfelder & Horvath), verbatim: "We do not recommend attempting WGCNA on a data set consisting of fewer than 15 samples. In a typical high-throughput setting, correlations on fewer than 15 samples will simply be too noisy for the network to be biologically meaningful. If at all possible, one should have at least 20 samples." Hard exclude.
- **mixOmics (sPLS-DA etc.)** — multivariate methods need samples ≫ classes for cross-validation; 4 samples cannot support it.
- **Variational Autoencoders / Graph Neural Networks** — need thousands of training samples; 4 samples → instant overfit. Exclude.
- **AlphaFold / structure prediction** — orthogonal to differential expression; per-protein structure does not answer the abundance question and is enormous compute for no DE insight. Exclude (consider only for a single follow-up protein of interest).

---

### SECTION 5 — Limitations as Engineering Constraints

- **Variance estimation at n=2:** each group has 2 values; the per-protein sample variance has 1 degree of freedom and is extremely unstable. limma's empirical-Bayes shrinkage is what makes inference possible, but it cannot manufacture replication. Failure mode: unmoderated t-tests give near-random p-values.
- **Technical vs biological replicates:** the two runs are technical replicates (same biological material, re-measured). They capture instrument/prep noise, **not** biological variability. Therefore p-values describe "is this ratio reproducible on this machine," not "does testosterone change this protein in the population." This is the single most important interpretive caveat and must be stated on every output.
- **PCA on 4 samples:** the sample covariance matrix is rank-deficient (rank ≤ n−1 = 3); only 3 PCs exist; eigenvalues are biased (first ones overestimated, per the n<p literature). Failure mode: PC1 can be driven entirely by one sample. Gate behind `n_samples ≥ 6` per group or label as QC-only.
- **Clustering stability:** hierarchical clustering / k-means on 4 samples has no stability (no bootstrap support possible). Dendrograms are not interpretable.
- **Multiple-testing power:** BH-FDR across ~1,948 tests with n=2 has very low power; few or zero proteins may survive `adj.P.Val < 0.05`. Plan to also report ranked nominal p-values and effect sizes, and treat results as hypothesis-generating.
- **WGCNA / UMAP / ML minimums:** WGCNA ≥15 (ideally ≥20) samples; UMAP/t-SNE dozens+; VAEs/GNNs thousands. All are hard-gated off.

Architecturally, every analysis registers a `min_samples` and `min_replicates_per_group` requirement; a dispatcher runs only those whose requirements are met and emits a skip-log for the rest.

---

### SECTION 6 — Master Build List (dependency-ordered roadmap)

| # | Module / script | Lang/tool | Depends on | P? | One-line description |
|---|---|---|---|---|---|
| 1 | `config/sample_sheet.tsv` + `config.yaml` | YAML/TSV | — | No | Declarative samples, groups, contrasts, thresholds, paths |
| 2 | `io/load.py` | Python/pandas | 1 | No | Read both SILAC sheets, normalize column names |
| 3 | `qc/schema.py` (pandera) | Python/pandera | 2 | No | Validate schema, types, UniProt regex, uniqueness |
| 4 | `etl/merge.py` | Python/pandas | 2,3 | No | Safe join w/ cardinality guard + dropped-exclusives log |
| 5 | `etl/foldchange.py` | Python/numpy | 4 | No | log2-space ratios, symmetric thresholds, inf/NaN masking (fixes bugs 1–3) |
| 6 | `qc/replicate_qc.py` | Python | 5 | No | Replicate correlation, distribution, missingness checks |
| 7 | `etl/build_matrix.py` | Python | 4 | No | Emit `intensity_matrix.tsv` + `design.tsv` (limma contract) |
| 8 | `de/run_limma.R` | R/limma+DEP | 7 | →makes | MinProb impute + limma; emits `limma_results.tsv` |
| 9 | `de/run_limma.py` | Python/subprocess | 7,8 | No | Orchestrate Rscript, check exit code, read results |
| 10 | `export/ipa_export.py` | Python | 8,9 | Yes | Write IPA CSV (UniProt + log2FC + FDR), export validation |
| 11 | `viz/volcano.py` / `.R` | Python/R | 8 | Yes | Volcano plot |
| 12 | `viz/ma_plot.py` | Python | 5/8 | No | MA plot |
| 13 | `viz/heatmap.R` (pheatmap) | R | 7 | No | Top-variable/DE protein heatmap |
| 14 | `viz/qc_plots.py` | Python/seaborn | 6 | No | Correlation matrix, boxplots, density, missing-value heatmap |
| 15 | `enrich/string_ppi.py` | Python/requests+igraph | 8/10 | No | STRING API → network → igraph metrics |
| 16 | `enrich/ora_gsea.R` | R/clusterProfiler+fgsea | 8 | Yes | GO/KEGG/Reactome ORA + ranked GSEA |
| 17 | `viz/network_cytoscape.py` | py4cytoscape | 15 | No | Publication network with log2FC node coloring |
| 18 | `viz/upset.py` | Python/upsetplot | 16 | No | Set-intersection plots |
| 19 | `gated/pca_cluster.py` | Python | 7 | No | PCA/clustering — runs only if `n≥6/group` (QC-only otherwise) |
| 20 | `report/build_report.py` | Python | all | mixed | Assemble HTML report + provenance |

---

### FORWARD-PATH SECTION — If biological replicates arrive

**Design to make the transition a config change, not a rewrite.** Everything keys off `config/sample_sheet.tsv` and a parameterized `model.matrix(~group)` (plus optional batch/covariate terms). Adding biological replicates means appending rows to the sample sheet and re-running; the design matrix regenerates automatically. The file contract (`intensity_matrix.tsv` + `design.tsv` + `limma_results.tsv`) is replicate-count-agnostic by construction.

What changes statistically and which analyses unlock:
- **Design matrix:** extend to `~ batch + group` (or `~ group` with biological replicates), enabling blocking on batch and the `duplicateCorrelation` approach if technical replicates are nested within biological ones.
- **Imputation reconsideration:** with more samples, missingness patterns are estimable; consider MAR-aware methods (kNN/SeqKNN, Impseq — also high-performing in Peng et al.) alongside MinProb, and re-evaluate per the paper's setting-specific guidance.
- **Now-reliable analyses:** PCA (≥6 samples), hierarchical clustering with bootstrap support, real biological-variance p-values with meaningful FDR power.
- **Newly appropriate tools:** at ≥15–20 samples, **WGCNA** co-expression modules; **mixOmics** sPLS-DA; **UMAP** for sample structure; potentially ML classifiers with proper cross-validation. Each is registered with its `min_samples` gate so the dispatcher turns them on automatically once the sample sheet grows.
- **Interpretation shift:** p-values graduate from "technical reproducibility" to "population-level biological significance" — the caveat banner is removed only when biological replicates exist.

Recommended folder structure for the new module:

```
proteomics_de/
├── config/
│   ├── config.yaml
│   └── sample_sheet.tsv
├── io/            load.py
├── qc/            schema.py  replicate_qc.py
├── etl/           merge.py  foldchange.py  build_matrix.py
├── de/            run_limma.R  run_limma.py
├── enrich/        string_ppi.py  ora_gsea.R
├── export/        ipa_export.py
├── viz/           volcano.py  ma_plot.py  heatmap.R  qc_plots.py
│                  network_cytoscape.py  upset.py
├── gated/         pca_cluster.py
├── report/        build_report.py
├── tests/         test_*.py  (pytest)
├── data/          raw/  interim/  processed/
├── results/       de/  figures/  ipa/  qc/
└── README.md
```

## Recommendations

1. **Freeze the legacy script; build `proteomics_de/` fresh.** Do not edit `proteomics_analysis.py`; reproduce its corrected logic in `etl/foldchange.py`. Benchmark: new pipeline reproduces the same protein list minus the three bug-induced miscalls.
2. **Fix the three fold-change bugs first** (log-space averaging, symmetric thresholds, inf/NaN masking) — these change results today, before any new analysis.
3. **Stand up the limma module (items 7–9) next.** This is the highest-value addition: it produces the p-values that unlock the volcano plot, ORA/GSEA, and IPA's full statistical features. Threshold to watch: if `0` proteins survive `adj.P.Val<0.05`, report ranked nominal p-values + effect sizes and explicitly frame results as hypothesis-generating.
4. **Switch IPA input to the validated CSV contract** (UniProt + Expr Log Ratio + FDR), single header row, UTF-8 no BOM, identifier leftmost.
5. **Gate every multivariate/ML analysis behind a replicate-count check.** Build PCA/clustering as QC-only now; leave WGCNA/UMAP/mixOmics unbuilt until the sample sheet reaches ≥15 biological samples.
6. **Make reproducibility first-class:** pin R/Python package versions, emit `.meta` sidecars, set seeds (MinProb is stochastic), and run `pytest` + pandera in CI.

## Caveats
- I read `proteomics_analysis.py` in full; I could not fetch the contents of the `Pilot Project` folder or `Limma.pdf` directly (GitHub tree pages blocked automated access), so the existing R work is inferred from the repo's 60.7% R language share and the `Limma.pdf` artifact rather than line-audited. Re-audit those files before finalizing the R module to avoid duplicating existing code.
- The n=2 technical-replicate limitation is fundamental: no statistical method, including limma, converts technical replicates into population-level biological inference. All p-values must be labeled accordingly.
- QIAGEN IPA's Salesforce KB pages block automated fetching; the format facts (20-observation limit, optional p-value, accepted IDs, tab-delimited recommendation) come from indexed captures of official QIAGEN pages and corroborating manuals, not a live fetch. The verbatim quotes above were cross-checked against multiple captures.
- Peng et al. 2024's recommendations were derived on spike-in DDA/DIA label-free benchmarks; this dataset is SILAC, so "no normalization" should be verified against the centered-log-ratio QC check rather than assumed.

## Build Log

### Bug 6 — Replicate correlation check ✅ (2026-06-27)
- New module replicate_check.py. Reports fold-change agreement between runs
  (Pearson r of log2_rep1 vs log2_rep2 + sign agreement) and raw reproducibility
  per condition. Verdict on the fold-change r vs REPLICATE_FC_R_MIN=0.50.
- Result: fold-change r = 0.2703 -> WARN; sign agreement 889/1578 (56.3%).
  Raw reproducibility good: control r = 0.862 (n=1723), treated r = 0.841 (n=1656).
- Lesson: raw readings reproduce well, but the fold-change (a difference) is noisy
  because subtracting two large similar numbers amplifies small wobbles. So
  per-protein up/down calls at n=2 are not individually trustworthy.
- This is the on-paper justification for Bug 7: limma borrows information across
  proteins to stabilize estimates precisely when single-protein replicates wobble.
- Outputs: qc_replicate_correlation.csv + replicate_correlation.png. Check-only;
  six protected outputs byte-identical (sha256). Bug numbering canonical.

### Bug 7 — Statistical testing via limma + MinProb ✅ (2026-06-27)

What: Per-protein significance. limma (R) called from Python over a CSV/subprocess
file handoff; MinProb imputation for missing cells. Adds raw p-value, BH-adjusted
p-value, and a combined "regulated AND significant" call.

Design (locked before spec):
- Fork 1 — file/CSV handoff; auditable intermediates kept; R fails loud; versions pinned.
- Fork 3 — input = 4 raw intensities (31578/31580 control, 31579/31581 treated),
  log2, un-centered (per Peng).
- Fork 2 — limma on the both-group only. Eligible = 1938 (1948 − 10 ON_OFF). 606
  single-condition + 10 ON_OFF stay in their own CSVs, untested (an entire condition
  would be imputed = manufactured data).

Files: limma_test.R (worker), limma_test.py (orchestrator), foldchange.py (acyclic
wiring after Bug 5/6). Outputs: results/qc_limma.csv (1938 rows),
results/ipa_input_significant.csv; intermediates _limma_input.csv / _limma_output.csv
+ _limma_versions.txt in proteomics_de/.

Result: 0/1938 significant at FDR<0.05. Smallest raw p = 1.6e-4 (63 at raw p<0.05);
smallest BH-adjusted p = 0.305. corr(limma_log2FC, pipeline log2FC) = 1.0000 → contrast
correct, not a wiring error. All 715 fold-change "regulated" calls land in probably-noise
→ ipa_input_significant.csv is header-only by design.

Lesson: This is the finding Bug 7 was built to surface, not a failure. At n=2 with
technical-only replicates (no biological replicates), large fold-changes do not survive
multiple-testing correction; the ±0.585 list overstated confidence. Consistent with Bug 6
(replicate FC r=0.27, WARN).

Env: R 4.6.1, limma 3.68.4, imputeLCMD 2.1 (Homebrew site-library, found by plain Rscript).

Open follow-up: Re-run with eBayes(trend=TRUE, robust=TRUE) — the field-standard proteomics
refinement deliberately deferred from the first verified run — before treating "zero
significant" as final.

### Bug 7 follow-up — trend/robust eBayes experiment ✅ (2026-06-28)

Re-ran limma statistical testing with eBayes(trend=TRUE, robust=TRUE), the
field-standard proteomics refinement deferred from the first run, to confirm the
"zero significant" baseline before treating it as final. Parameterized
limma_test.R/limma_test.py with an optional ebayes_mode arg (default vanilla,
verified byte-identical to the committed baseline by sha256); experiment reused the
existing _limma_input.csv so the input is provably identical. Result: still 0/1938
significant at FDR < 0.05. The refinement sharpened the extremes — smallest
BH-adjusted p fell 0.305 → 0.116 (top candidate's gap to the 0.05 line roughly
halved, but still non-significant), raw-p<0.05 count trimmed 63 → 55 as robust
moderation down-weighted outlier variances. Conclusion: the n=2 technical-replicate
ceiling dominates regardless of eBayes flavor; trend/robust corroborates the baseline
rather than overturning it. Output in results/qc_limma_trend.csv (canonical
qc_limma.csv untouched).