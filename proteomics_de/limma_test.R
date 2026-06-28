#!/usr/bin/env Rscript
# Bug 7 — per-protein statistical testing via limma, with MinProb imputation.
#
# This is the R "worker" half of Bug 7. It owns the statistics; the Python
# orchestrator (limma_test.py) owns the subprocess boundary and the file handoff.
#
# Invocation:
#   Rscript limma_test.R <input_csv> <output_csv> <seed>
#
# Input  (written by Python): columns id, gene, ctrl_31578, ctrl_31580,
#         trt_31579, trt_31581 — raw intensities; blank cell == missing.
# Output (read by Python):    columns id, limma_log2FC, p_value, adj_p_value.
#
# Pipeline (Peng et al. 2024, Nat Commun, DOI 10.1038/s41467-024-47899-w):
#   raw intensity -> (<=0 -> NA) -> log2 (un-centered, no normalization)
#   -> MinProb imputation -> limma moderated t-test -> BH adjusted p-values.
#
# Reproducibility: set.seed() before MinProb (it draws random values); R + package
# versions are written to _limma_versions.txt next to the output.
#
# Errors: any failure stops with a "BUG7 R ERROR: ..." message on stderr and a
# nonzero exit code, so the Python side can fail loudly and never proceed on
# empty/partial p-values.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("BUG7 R ERROR: usage: Rscript limma_test.R <input_csv> <output_csv> <seed>")
}
input_csv  <- args[[1L]]
output_csv <- args[[2L]]
seed       <- suppressWarnings(as.integer(args[[3L]]))
if (is.na(seed)) {
  stop("BUG7 R ERROR: seed must be an integer, got: ", args[[3L]])
}

# Required packages: limma (Bioconductor), imputeLCMD (CRAN).
if (!requireNamespace("limma", quietly = TRUE) ||
    !requireNamespace("imputeLCMD", quietly = TRUE)) {
  stop(paste0(
    "BUG7 R ERROR: required R packages missing (need 'limma' and 'imputeLCMD'). ",
    "Install with:\n",
    "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install('limma'); install.packages('imputeLCMD')"
  ))
}

run <- function() {
  df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

  required_cols <- c("id", "gene", "ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0L) {
    stop("BUG7 R ERROR: input missing columns: ", paste(missing, collapse = ", "))
  }
  if (nrow(df) == 0L) {
    stop("BUG7 R ERROR: input has zero rows.")
  }

  # Build the numeric matrix in the FIXED column order: control, control,
  # treated, treated (= 31578, 31580, 31579, 31581).
  intensity_cols <- c("ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581")
  mat <- sapply(
    df[, intensity_cols, drop = FALSE],
    function(x) suppressWarnings(as.numeric(as.character(x)))
  )
  mat <- as.matrix(mat)

  # Any value <= 0 is missing; then log2-transform (un-centered, no normalization).
  mat[!is.na(mat) & mat <= 0] <- NA
  mat <- log2(mat)

  # MinProb imputation on the log2 matrix. set.seed() immediately before it,
  # because MinProb draws random values (reproducibility hinges on this).
  # impute.MinProb prints an internal value to stdout; capture.output swallows it
  # so the worker's stdout stays clean (the assignment still happens).
  set.seed(seed)
  invisible(utils::capture.output(
    mat <- imputeLCMD::impute.MinProb(mat, q = 0.01, tune.sigma = 1)
  ))

  if (any(!is.finite(mat))) {
    stop("BUG7 R ERROR: non-finite values remain after imputation.")
  }

  # Two-group design; coef 2 is treated - control.
  group <- c("control", "control", "treated", "treated")
  design <- model.matrix(~ factor(group, levels = c("control", "treated")))

  fit <- limma::lmFit(mat, design)
  # Plain eBayes for this first verified run, so results are reproducible. A common
  # proteomics refinement is eBayes(fit, trend = TRUE, robust = TRUE) (intensity-
  # dependent prior variance + outlier-robust moderation); leave it as a future
  # toggle rather than changing the baseline. (renv for full R-environment pinning
  # is likewise a future refinement; _limma_versions.txt is enough for now.)
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(
    fit, coef = 2L, number = Inf, sort.by = "none", adjust.method = "BH"
  )

  if (nrow(tt) != nrow(df)) {
    stop("BUG7 R ERROR: topTable returned ", nrow(tt),
         " rows for ", nrow(df), " input proteins.")
  }
  if (any(!is.finite(tt$P.Value)) || any(!is.finite(tt$adj.P.Val)) ||
      any(!is.finite(tt$logFC))) {
    stop("BUG7 R ERROR: non-finite logFC/p-value/adj-p produced.")
  }

  out <- data.frame(
    id           = df$id,
    limma_log2FC = tt$logFC,
    p_value      = tt$P.Value,
    adj_p_value  = tt$adj.P.Val,
    stringsAsFactors = FALSE
  )
  write.csv(out, output_csv, row.names = FALSE)

  # Reproducibility record, next to the output file.
  versions_path <- file.path(dirname(output_csv), "_limma_versions.txt")
  writeLines(
    c(
      R.version.string,
      paste0("limma ", as.character(utils::packageVersion("limma"))),
      paste0("imputeLCMD ", as.character(utils::packageVersion("imputeLCMD"))),
      paste0("seed ", seed)
    ),
    versions_path
  )

  invisible(TRUE)
}

tryCatch(
  run(),
  error = function(e) {
    msg <- conditionMessage(e)
    if (!grepl("^BUG7 R ERROR", msg)) {
      msg <- paste0("BUG7 R ERROR: ", msg)
    }
    stop(msg, call. = FALSE)
  }
)
