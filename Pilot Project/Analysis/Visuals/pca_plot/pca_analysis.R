# ============================================================================
# Publication-grade PCA Plot for SILAC Proteomics Data
# Vehicle Control vs Testosterone Treatment
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(ggrepel)
  library(grid)
})

# ============================================================================
# Load and prepare data
# ============================================================================

script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "high_confidence_proteins.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Select intensity columns
intensity_cols <- c("Vehicle_Rep1_31579", "Vehicle_Rep2_31581",
                    "Testosterone_Rep1_31578", "Testosterone_Rep2_31580")

# Extract intensity data (proteins as rows, samples as columns)
intensity_matrix <- as.matrix(df[, intensity_cols])
rownames(intensity_matrix) <- df$Gene

# Remove rows with missing or zero values
complete_mask <- complete.cases(intensity_matrix) &
                 rowSums(intensity_matrix == 0, na.rm = TRUE) == 0
intensity_complete <- intensity_matrix[complete_mask, ]

n_proteins <- nrow(intensity_complete)

# Log10 transform
log_intensity <- log10(intensity_complete + 1)

# Z-score normalize across samples for each protein (row-wise)
zscore_matrix <- t(scale(t(log_intensity)))

# Transpose: samples as rows, proteins as columns
data_for_pca <- t(zscore_matrix)

# ============================================================================
# Run PCA
# ============================================================================

pca_result <- prcomp(data_for_pca, center = TRUE, scale. = FALSE)

# Extract scores (sample coordinates)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Sample <- rownames(pca_scores)

# Create sample metadata
pca_scores <- pca_scores %>%
  mutate(
    Condition = case_when(
      grepl("Vehicle", Sample) ~ "Vehicle",
      grepl("Testosterone", Sample) ~ "Testosterone"
    ),
    Replicate = case_when(
      grepl("Rep1", Sample) ~ "Rep 1",
      grepl("Rep2", Sample) ~ "Rep 2"
    ),
    Display_Name = case_when(
      Sample == "Vehicle_Rep1_31579" ~ "Vehicle Rep1",
      Sample == "Vehicle_Rep2_31581" ~ "Vehicle Rep2",
      Sample == "Testosterone_Rep1_31578" ~ "Testosterone Rep1",
      Sample == "Testosterone_Rep2_31580" ~ "Testosterone Rep2"
    )
  )

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2, ] * 100
cumvar_explained <- summary(pca_result)$importance[3, ] * 100

pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)
total_var <- round(pc1_var + pc2_var, 1)

# Extract loadings (protein contributions)
loadings <- as.data.frame(pca_result$rotation)
loadings$Gene <- rownames(loadings)

# Get top proteins contributing to PC1 and PC2
top_pc1 <- loadings %>%
  arrange(desc(abs(PC1))) %>%
  head(5)

top_pc2 <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  head(5)

# Combine unique proteins for biplot
top_proteins <- unique(c(top_pc1$Gene, top_pc2$Gene))
loadings_subset <- loadings %>% filter(Gene %in% top_proteins)

# Scale loadings for visualization
loading_scale <- min(max(abs(pca_scores$PC1)), max(abs(pca_scores$PC2))) * 0.8
loadings_subset <- loadings_subset %>%
  mutate(
    PC1_scaled = PC1 * loading_scale / max(abs(PC1)),
    PC2_scaled = PC2 * loading_scale / max(abs(PC2))
  )

# ============================================================================
# Calculate ellipses for confidence regions
# ============================================================================

# Function to calculate ellipse points
calc_ellipse <- function(data, x_col, y_col, level = 0.95) {
  x <- data[[x_col]]
  y <- data[[y_col]]

  if (length(x) < 2) return(NULL)

  center_x <- mean(x)
  center_y <- mean(y)

  # For only 2 points, create a simple ellipse around the center
  if (length(x) == 2) {
    diff_x <- abs(x[1] - x[2]) / 2 + 0.5
    diff_y <- abs(y[1] - y[2]) / 2 + 0.5

    theta <- seq(0, 2*pi, length.out = 100)
    ellipse_x <- center_x + diff_x * cos(theta)
    ellipse_y <- center_y + diff_y * sin(theta)

    return(data.frame(x = ellipse_x, y = ellipse_y))
  }

  # For more points, use covariance
  cov_mat <- cov(cbind(x, y))
  eig <- eigen(cov_mat)

  # Chi-square value for confidence level
  chi_sq <- qchisq(level, df = 2)

  # Ellipse parameters
  a <- sqrt(chi_sq * eig$values[1])
  b <- sqrt(chi_sq * eig$values[2])
  angle <- atan2(eig$vectors[2, 1], eig$vectors[1, 1])

  theta <- seq(0, 2*pi, length.out = 100)
  ellipse_x <- center_x + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
  ellipse_y <- center_y + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)

  return(data.frame(x = ellipse_x, y = ellipse_y))
}

# Calculate ellipses for each condition
vehicle_ellipse <- calc_ellipse(
  pca_scores %>% filter(Condition == "Vehicle"),
  "PC1", "PC2", level = 0.95
)
vehicle_ellipse$Condition <- "Vehicle"

testosterone_ellipse <- calc_ellipse(
  pca_scores %>% filter(Condition == "Testosterone"),
  "PC1", "PC2", level = 0.95
)
testosterone_ellipse$Condition <- "Testosterone"

ellipse_data <- rbind(vehicle_ellipse, testosterone_ellipse)

# ============================================================================
# Assess separation quality
# ============================================================================

# Calculate centroid distance vs within-group variance
vehicle_centroid <- c(
  mean(pca_scores$PC1[pca_scores$Condition == "Vehicle"]),
  mean(pca_scores$PC2[pca_scores$Condition == "Vehicle"])
)
testosterone_centroid <- c(
  mean(pca_scores$PC1[pca_scores$Condition == "Testosterone"]),
  mean(pca_scores$PC2[pca_scores$Condition == "Testosterone"])
)

centroid_dist <- sqrt(sum((vehicle_centroid - testosterone_centroid)^2))

# Simple separation assessment
if (pc1_var > 50 && centroid_dist > 1) {
  separation_quality <- "Good"
  separation_symbol <- "\u2713"  # Checkmark
  separation_color <- "#2ECC71"
} else {
  separation_quality <- "Moderate"
  separation_symbol <- "\u25CB"  # Circle
  separation_color <- "#F39C12"
}

# ============================================================================
# Create main PCA plot
# ============================================================================

# Define colors and shapes
condition_colors <- c("Vehicle" = "#3498DB", "Testosterone" = "#E74C3C")
replicate_shapes <- c("Rep 1" = 21, "Rep 2" = 22)  # Circle, Square

main_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2)) +

  # Reference lines at origin
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +

  # Confidence ellipses
  geom_polygon(data = ellipse_data,
               aes(x = x, y = y, fill = Condition),
               alpha = 0.1, color = NA) +
  geom_path(data = ellipse_data,
            aes(x = x, y = y, color = Condition),
            linetype = "dashed", linewidth = 1, alpha = 0.6) +

  # Loading arrows (biplot)
  geom_segment(data = loadings_subset,
               aes(x = 0, y = 0, xend = PC1_scaled, yend = PC2_scaled),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
               color = "gray40", linewidth = 0.8, alpha = 0.7) +

  # Loading labels
  geom_text_repel(data = loadings_subset,
                  aes(x = PC1_scaled, y = PC2_scaled, label = Gene),
                  size = 2.8, color = "gray30", fontface = "italic",
                  box.padding = 0.3, point.padding = 0.2,
                  max.overlaps = 15) +

  # Sample points
  geom_point(aes(fill = Condition, shape = Replicate),
             size = 8, color = "black", stroke = 1.5, alpha = 0.9) +

  # Sample labels
  geom_text_repel(aes(label = Display_Name),
                  size = 3.5, fontface = "bold",
                  box.padding = 0.8, point.padding = 0.5,
                  max.overlaps = 10) +

  # Scales
  scale_fill_manual(values = condition_colors, name = "Condition") +
  scale_color_manual(values = condition_colors, guide = "none") +
  scale_shape_manual(values = replicate_shapes, name = "Replicate") +

  # Labels
  labs(
    title = "Principal Component Analysis - Secreted Proteome",
    subtitle = "Vehicle Control vs Testosterone Treatment",
    x = paste0("PC1 (", pc1_var, "% variance explained)"),
    y = paste0("PC2 (", pc2_var, "% variance explained)")
  ) +

  # Annotation boxes
  annotate("label",
           x = min(pca_scores$PC1) - 0.5, y = max(pca_scores$PC2) + 0.3,
           label = paste0("Separation: ", separation_quality, " ", separation_symbol, "\n",
                         "PC1+PC2 variance: ", total_var, "%"),
           hjust = 0, vjust = 1, size = 3.2,
           fill = "white", color = separation_color,
           label.padding = unit(0.4, "lines"),
           fontface = "bold") +

  annotate("label",
           x = max(pca_scores$PC1) + 0.5, y = min(pca_scores$PC2) - 0.3,
           label = paste0("n = ", n_proteins, " proteins\n",
                         "4 samples (2 conditions \u00D7 2 replicates)"),
           hjust = 1, vjust = 0, size = 3,
           fill = "white", color = "gray40",
           label.padding = unit(0.4, "lines")) +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 15)),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.margin = margin(5, 8, 5, 8),
    plot.margin = margin(15, 15, 15, 15)
  ) +

  # Coordinate expansion
  coord_cartesian(clip = "off")

# ============================================================================
# Create scree plot inset
# ============================================================================

scree_data <- data.frame(
  PC = paste0("PC", 1:4),
  Variance = var_explained[1:4],
  Highlighted = c(TRUE, TRUE, FALSE, FALSE)
)
scree_data$PC <- factor(scree_data$PC, levels = scree_data$PC)

scree_plot <- ggplot(scree_data, aes(x = PC, y = Variance, fill = Highlighted)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  geom_text(aes(label = paste0(round(Variance, 0), "%")),
            vjust = -0.3, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#3498DB", "FALSE" = "gray70"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Variance Explained", x = NULL, y = "%") +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = "gray50", linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# Combine main plot with inset
# ============================================================================

final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(scree_plot, x = 0.68, y = 0.65, width = 0.28, height = 0.28)

# ============================================================================
# Save outputs
# ============================================================================

ggsave(
  filename = "pca_analysis.png",
  plot = final_plot,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "pca_analysis.pdf",
  plot = final_plot,
  width = 12,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

cat("PCA analysis figure saved to:\n")
cat("  - pca_analysis.png\n")
cat("  - pca_analysis.pdf\n")
cat("\nPCA Summary:\n")
cat(paste0("  Proteins analyzed: ", n_proteins, "\n"))
cat(paste0("  PC1 variance: ", pc1_var, "%\n"))
cat(paste0("  PC2 variance: ", pc2_var, "%\n"))
cat(paste0("  Total (PC1+PC2): ", total_var, "%\n"))
cat(paste0("  Separation quality: ", separation_quality, "\n"))
cat("\nTop proteins contributing to PC1:\n")
print(top_pc1[, c("Gene", "PC1")])
cat("\nTop proteins contributing to PC2:\n")
print(top_pc2[, c("Gene", "PC2")])
