suppressMessages({library(limma); library(imputeLCMD)})
R <- "/Users/abrarshakik/Documents/Proteomics/proteomics_de/"
df <- read.csv(paste0(R,"_limma_input.csv"), check.names=FALSE, stringsAsFactors=FALSE)
dz <- read.delim(paste0(R,"_limma_design.tsv"), colClasses="character")
cols <- dz$sample; group <- dz$group
mat <- as.matrix(sapply(df[,cols], function(x) suppressWarnings(as.numeric(as.character(x)))))
mat[!is.na(mat) & mat<=0] <- NA; mat <- log2(mat)
nimp <- rowSums(is.na(mat))
set.seed(42); invisible(capture.output(mat <- impute.MinProb(mat, q=0.01, tune.sigma=1)))
g <- factor(group, levels=c("control","treated"))
# pairing: exp1 = 31579+31578, exp2 = 31581+31580
expt <- factor(ifelse(grepl("31578|31579", cols), "E1", "E2"))
print(data.frame(cols, group, expt))
fitit <- function(design, mode, coef) {
  f <- lmFit(mat, design)
  f <- if (mode=="trend_robust") eBayes(f, trend=TRUE, robust=TRUE) else eBayes(f)
  tt <- topTable(f, coef=coef, number=Inf, sort.by="none", adjust.method="BH")
  cat(mode, colnames(design)[coef], "resid df:", unique(f$df.residual), " prior df (median):", median(f$df.prior), "\n")
  data.frame(id=df$id, gene=df$gene, logFC=tt$logFC, p=tt$P.Value, adj=tt$adj.P.Val, t=tt$t, n_imputed=nimp)
}
dU <- model.matrix(~g); dP <- model.matrix(~expt+g)
for (m in c("trend_robust","vanilla")) {
  u <- fitit(dU, m, 2); p <- fitit(dP, m, "gtreated")
  write.csv(u, sprintf("unpaired_%s.csv", m), row.names=FALSE)
  write.csv(p, sprintf("paired_%s.csv", m), row.names=FALSE)
}
