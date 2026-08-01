#!/usr/bin/env Rscript
# Bug 7 — per-protein statistical testing via limma, with MinProb imputation.
#
# This is the R "worker" half of Bug 7. It owns the statistics; the Python
# orchestrator (limma_test.py) owns the subprocess boundary and the file handoff.
#
# Invocation (two equivalent forms):
#   Rscript limma_test.R <input_csv> <output_csv> <seed> [ebayes_mode]
#   Rscript limma_test.R --in <csv> --out <csv> --seed <int> [--mode <m>] [--design <tsv>]
#
#   The positional form is the original one and is kept working byte-for-byte.
#   The flag form adds --design, which is the point of the exercise: with it the
#   group vector and the expected intensity columns are READ FROM THE DESIGN FILE
#   instead of being hardcoded here, so adding replicates is an edit to
#   config/sample_sheet.tsv and nothing else.
#
#   ebayes_mode / --mode (optional): "vanilla" (default) or "trend_robust". The
#   mode only switches the eBayes call; everything upstream/downstream is
#   identical, so the default invocation reproduces the committed baseline exactly.
#
#   Flags are parsed by hand (below) rather than with optparse: optparse is not
#   installed here, and five flags do not justify a new dependency in a script
#   whose whole value is that it reproduces a frozen numeric result.
#
# Input  (written by Python): columns id, gene, then one column per sample —
#         ctrl_31578, ctrl_31580, trt_31579, trt_31581 for today's design.
#         Raw intensities; blank cell == missing.
# Design (optional, --design): TAB-separated, columns sample, group. `sample`
#         names the INPUT CSV's intensity columns (i.e. the handoff names), in
#         the order the design matrix expects: control group first.
# Output (read by Python):    columns id, limma_log2FC, p_value, adj_p_value.
#
# Pipeline (Peng et al. 2024, Nat Commun, DOI 10.1038/s41467-024-47899-w):
#   raw intensity -> (<=0 -> NA) -> log2 (un-centered, no normalization)
#   -> MinProb imputation -> limma moderated t-test -> BH adjusted p-values.
#
# Reproducibility: set.seed() before MinProb (it draws random values); R + package
# versions are written to _limma_versions.txt next to the output.
#
# Errors: any failure prints a message whose FIRST characters on stderr are
# "BUG7 R ERROR" and exits nonzero, so the Python side can fail loudly and never
# proceed on empty/partial p-values. bug7_abort() exists because a bare stop()
# would let R prepend its own "Error:" / "Error in f(x) :" banner, which breaks
# that prefix contract for anything trying to detect our errors specifically.

bug7_abort <- function(...) {
  msg <- paste0(...)
  if (!grepl("^BUG7 R ERROR", msg)) {
    msg <- paste0("BUG7 R ERROR: ", msg)
  }
  cat(msg, "\n", sep = "", file = stderr())
  quit(save = "no", status = 1L)
}

args <- commandArgs(trailingOnly = TRUE)

# --- Argument parsing: flag form if any arg starts with "--", else positional --
# The two forms are never mixed; that keeps the positional path exactly as it was.
parse_flags <- function(args) {
  keys <- c("--in" = "input_csv", "--out" = "output_csv", "--design" = "design_tsv",
            "--seed" = "seed", "--mode" = "ebayes_mode")
  opt <- list(input_csv = NA_character_, output_csv = NA_character_,
              design_tsv = NA_character_, seed = NA_character_,
              ebayes_mode = "vanilla")
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% names(keys)) {
      bug7_abort("BUG7 R ERROR: unknown flag: ", key, " (expected one of ",
                 paste(names(keys), collapse = ", "), ")")
    }
    if (i + 1L > length(args)) {
      bug7_abort("BUG7 R ERROR: flag ", key, " requires a value.")
    }
    opt[[keys[[key]]]] <- args[[i + 1L]]
    i <- i + 2L
  }
  opt
}

if (any(startsWith(args, "--"))) {
  opt <- parse_flags(args)
  if (is.na(opt$input_csv) || is.na(opt$output_csv) || is.na(opt$seed)) {
    bug7_abort("BUG7 R ERROR: usage: Rscript limma_test.R --in <csv> --out <csv> ",
               "--seed <int> [--mode <mode>] [--design <tsv>]")
  }
  input_csv  <- opt$input_csv
  output_csv <- opt$output_csv
  seed_arg   <- opt$seed
  ebayes_mode <- opt$ebayes_mode
  design_tsv <- opt$design_tsv
} else {
  if (length(args) < 3L) {
    bug7_abort("BUG7 R ERROR: usage: Rscript limma_test.R <input_csv> <output_csv> <seed>")
  }
  input_csv  <- args[[1L]]
  output_csv <- args[[2L]]
  seed_arg   <- args[[3L]]
  # Optional 4th arg: eBayes mode. Default keeps the committed baseline behavior.
  ebayes_mode <- if (length(args) >= 4L) args[[4L]] else "vanilla"
  design_tsv <- NA_character_
}

seed <- suppressWarnings(as.integer(seed_arg))
if (is.na(seed)) {
  bug7_abort("BUG7 R ERROR: seed must be an integer, got: ", seed_arg)
}
if (!ebayes_mode %in% c("vanilla", "trend_robust")) {
  bug7_abort("BUG7 R ERROR: ebayes_mode must be 'vanilla' or 'trend_robust', got: ",
             ebayes_mode)
}

# Required packages: limma (Bioconductor), imputeLCMD (CRAN).
if (!requireNamespace("limma", quietly = TRUE) ||
    !requireNamespace("imputeLCMD", quietly = TRUE)) {
  bug7_abort(paste0(
    "BUG7 R ERROR: required R packages missing (need 'limma' and 'imputeLCMD'). ",
    "Install with:\n",
    "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
    "  BiocManager::install('limma'); install.packages('imputeLCMD')"
  ))
}

# Resolve the experimental design: from --design when given, else the historical
# hardcoded 2x2 SILAC layout. Both paths produce the same three values, so
# everything downstream (including model.matrix) is a single code path.
resolve_design <- function(design_tsv) {
  if (is.na(design_tsv)) {
    return(list(
      intensity_cols = c("ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581"),
      group          = c("control", "control", "treated", "treated"),
      group_levels   = c("control", "treated")
    ))
  }
  if (!file.exists(design_tsv)) {
    stop("BUG7 R ERROR: design file not found: ", design_tsv)
  }
  dz <- read.delim(design_tsv, stringsAsFactors = FALSE, check.names = FALSE,
                   colClasses = "character")
  missing_design <- setdiff(c("sample", "group"), colnames(dz))
  if (length(missing_design) > 0L) {
    stop("BUG7 R ERROR: design file missing columns: ",
         paste(missing_design, collapse = ", "))
  }
  if (nrow(dz) == 0L) {
    stop("BUG7 R ERROR: design file has zero rows.")
  }
  group <- trimws(as.character(dz$group))
  if (length(unique(group)) < 2L) {
    stop("BUG7 R ERROR: design file has fewer than two groups: ",
         paste(unique(group), collapse = ", "))
  }
  list(
    intensity_cols = trimws(as.character(dz$sample)),
    group          = group,
    # Levels in order of first appearance. The design file is written in
    # canonical order (control group first), so level 1 is the reference and
    # coef 2 stays "treated - control" -- the sign of every logFC depends on it.
    group_levels   = unique(group)
  )
}

run <- function() {
  df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

  dsn <- resolve_design(design_tsv)
  intensity_cols <- dsn$intensity_cols
  group          <- dsn$group
  group_levels   <- dsn$group_levels

  required_cols <- c("id", "gene", intensity_cols)
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0L) {
    stop("BUG7 R ERROR: input missing columns: ", paste(missing, collapse = ", "))
  }
  if (nrow(df) == 0L) {
    stop("BUG7 R ERROR: input has zero rows.")
  }
  if (length(group) != length(intensity_cols)) {
    stop("BUG7 R ERROR: design has ", length(group), " rows but ",
         length(intensity_cols), " sample columns.")
  }

  # Build the numeric matrix in the design's column order: control group first,
  # treated second (= 31578, 31580, 31579, 31581 for today's sheet).
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

  # Two-group design; coef 2 is treated - control. `group` / `group_levels` come
  # from resolve_design() above, so this line is identical for both invocation
  # forms -- which is exactly why --design cannot perturb the design matrix.
  design <- model.matrix(~ factor(group, levels = group_levels))

  fit <- limma::lmFit(mat, design)
  # eBayes mode toggle. "vanilla" (the committed baseline) keeps the plain call so
  # results stay reproducible. "trend_robust" applies the common proteomics
  # refinement: intensity-dependent prior variance (trend) + outlier-robust
  # moderation (robust). Only this call differs between the two modes.
  # (renv for full R-environment pinning remains a future refinement;
  # _limma_versions.txt is enough for now.)
  if (ebayes_mode == "trend_robust") {
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  } else {
    fit <- limma::eBayes(fit)
  }
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

  # Reproducibility record, next to the output file. The versions filename is
  # derived from the output filename (s/_limma_output/_limma_versions/) so a
  # non-default output (e.g. the trend/robust experiment) writes its own versions
  # file instead of overwriting the committed _limma_versions.txt. For the default
  # output (_limma_output.csv) this resolves to _limma_versions.txt unchanged.
  out_base <- tools::file_path_sans_ext(basename(output_csv))
  ver_name <- sub("_limma_output", "_limma_versions", out_base)
  versions_path <- file.path(dirname(output_csv), paste0(ver_name, ".txt"))
  # Content stays exactly as the committed baseline (4 lines) so a default re-run
  # reproduces _limma_versions.txt byte-for-byte; the mode is captured by the
  # filename instead.
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
  error = function(e) bug7_abort(conditionMessage(e))
)
