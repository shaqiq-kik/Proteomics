# ============================================================================
# Publication-grade Clustered Heatmap for SILAC Proteomics Data
# High Confidence Secreted Proteins
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
  library(RColorBrewer)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "high_confidence_proteins.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Select only the 25 high-confidence proteins (all in the file are high confidence)
# Select intensity columns
intensity_cols <- c("Vehicle_Rep1_31579", "Vehicle_Rep2_31581",
                    "Testosterone_Rep1_31578", "Testosterone_Rep2_31580")

# Create matrix for heatmap
intensity_matrix <- as.matrix(df[, intensity_cols])
rownames(intensity_matrix) <- df$Gene

# Log10 transform (add small value to avoid log(0))
log_matrix <- log10(intensity_matrix + 1)

# Z-score normalize across rows (per protein)
zscore_matrix <- t(scale(t(log_matrix)))

# Clip z-scores to -2.5 to +2.5 for visualization
zscore_clipped <- zscore_matrix
zscore_clipped[zscore_clipped > 2.5] <- 2.5
zscore_clipped[zscore_clipped < -2.5] <- -2.5

# Rename columns for display
colnames(zscore_clipped) <- c("Vehicle Rep1", "Vehicle Rep2",
                               "Testosterone Rep1", "Testosterone Rep2")

# Create row annotation for Functional Class
functional_class_display <- case_when(
  df$Functional_Class == "cytokine" ~ "Cytokine",
  df$Functional_Class == "growth factor" ~ "Growth Factor",
  df$Functional_Class %in% c("enzyme", "peptidase") ~ "Enzyme",
  TRUE ~ "Other"
)

row_annotation <- data.frame(
  `Functional Class` = functional_class_display,
  row.names = df$Gene,
  check.names = FALSE
)

# Define annotation colors
annotation_colors <- list(
  `Functional Class` = c(
    "Cytokine" = "#E41A1C",
    "Growth Factor" = "#4DAF4A",
    "Enzyme" = "#377EB8",
    "Other" = "#808080"
  )
)

# Define color palette (RdBu reversed: red=high, blue=low)
color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Create breaks for color scale (-2.5 to +2.5)
breaks <- seq(-2.5, 2.5, length.out = 101)

# Generate heatmap and save to PNG
png(
  filename = "clustered_heatmap.png",
  width = 12,
  height = 16,
  units = "in",
  res = 300
)

pheatmap(
  zscore_clipped,
  # Clustering parameters
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  # Dendrogram
  treeheight_row = 50,
  treeheight_col = 30,
  # Colors
  color = color_palette,
  breaks = breaks,
  # Annotations
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  # Labels
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  # Cell styling
  border_color = "white",
  cellwidth = NA,
  cellheight = NA,
  # Title
  main = "Hierarchical Clustering - High Confidence Secreted Proteins",
  # Legend
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_row = TRUE,
  # Gaps
  gaps_col = 2  # Gap between Vehicle and Testosterone groups
)

dev.off()

# Generate PDF version
pdf(
  file = "clustered_heatmap.pdf",
  width = 12,
  height = 16
)

pheatmap(
  zscore_clipped,
  # Clustering parameters
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  # Dendrogram
  treeheight_row = 50,
  treeheight_col = 30,
  # Colors
  color = color_palette,
  breaks = breaks,
  # Annotations
  annotation_row = row_annotation,
  annotation_colors = annotation_colors,
  # Labels
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  fontsize = 10,
  # Cell styling
  border_color = "white",
  cellwidth = NA,
  cellheight = NA,
  # Title
  main = "Hierarchical Clustering - High Confidence Secreted Proteins",
  # Legend
  legend = TRUE,
  annotation_legend = TRUE,
  annotation_names_row = TRUE,
  # Gaps
  gaps_col = 2  # Gap between Vehicle and Testosterone groups
)

dev.off()

cat("Clustered heatmap saved to:\n")
cat("  - clustered_heatmap.png\n")
cat("  - clustered_heatmap.pdf\n")
