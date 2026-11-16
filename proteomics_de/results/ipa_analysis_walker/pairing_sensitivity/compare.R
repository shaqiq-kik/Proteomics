R <- "/Users/abrarshakik/Documents/Proteomics/proteomics_de/"
for (m in c("trend_robust","vanilla")) {
  ref <- read.csv(paste0(R, if (m=="vanilla") "_limma_output_vanilla.csv" else "_limma_output.csv"))
  u <- read.csv(sprintf("unpaired_%s.csv",m)); p <- read.csv(sprintf("paired_%s.csv",m))
  stopifnot(all(ref$id==u$id))
  cat("\n==", m, "== reproduce committed: max|dlogFC|", max(abs(ref$limma_log2FC-u$logFC)),
      " max|dp|", max(abs(ref$p_value-u$p)), " max|dadj|", max(abs(ref$adj_p_value-u$adj)),
      " max|dt|", max(abs(ref$t-u$t)), "\n")
  cat("paired vs unpaired max|dlogFC|", max(abs(p$logFC-u$logFC)), "\n")
  for (x in list(list("unpaired",u),list("paired",p))) { d <- x[[2]]
    cat(sprintf("%-9s FDR<0.05: %d | FDR<0.10: %d | raw p<0.05: %d | min adj p: %.4g | min p: %.3g\n", x[[1]],
      sum(d$adj<0.05), sum(d$adj<0.1), sum(d$p<0.05), min(d$adj), min(d$p))) }
  f <- which(u$gene=="Frzb"); cat("Frzb unpaired: logFC", u$logFC[f], "p", u$p[f], "adj", u$adj[f], "t", u$t[f], " n_imp", u$n_imputed[f], "\n")
  cat("Frzb paired:   logFC", p$logFC[f], "p", p$p[f], "adj", p$adj[f], "t", p$t[f], "\n")
  cat("Spearman p unpaired vs paired:", cor(u$p, p$p, method="spearman"), "\n")
  for (x in list(list("unpaired",u),list("paired",p))) { d <- x[[2]]; o <- head(order(d$p),10)
    cat("top10", x[[1]], ":\n"); print(data.frame(gene=d$gene[o], logFC=round(d$logFC[o],3), p=signif(d$p[o],3), adj=signif(d$adj[o],3), n_imp=d$n_imputed[o])) }
  cat("top10 overlap:", length(intersect(head(u$id[order(u$p)],10), head(p$id[order(p$p)],10))), "\n")
}
