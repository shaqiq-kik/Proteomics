# ============================================================================
# Publication-grade Correlation Heatmap for Technical Replicates
# SILAC Proteomics Data Quality Assessment
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "high_confidence_proteins.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Select intensity columns
intensity_cols <- c("Vehicle_Rep1_31579", "Vehicle_Rep2_31581",
                    "Testosterone_Rep1_31578", "Testosterone_Rep2_31580")

# Extract intensity data
intensity_data <- df[, intensity_cols]

# Remove rows with any missing or zero values
complete_mask <- complete.cases(intensity_data) &
                 rowSums(intensity_data == 0, na.rm = TRUE) == 0
intensity_complete <- intensity_data[complete_mask, ]

n_proteins <- nrow(intensity_complete)

# Log10 transform
log_intensity <- log10(intensity_complete)

# Calculate Pearson correlation matrix
cor_matrix <- cor(log_intensity, method = "pearson")

# Square the correlations to get r²
cor_matrix_r2 <- cor_matrix^2

# Create nice labels
labels <- c("Vehicle Rep 1", "Vehicle Rep 2", "Testosterone Rep 1", "Testosterone Rep 2")
rownames(cor_matrix_r2) <- labels
colnames(cor_matrix_r2) <- labels

# Calculate quality metrics
# Within-condition correlations
vehicle_within <- cor_matrix_r2[1, 2]
testosterone_within <- cor_matrix_r2[3, 4]
mean_within <- (vehicle_within + testosterone_within) / 2

# Between-condition correlations (all pairwise between Vehicle and Testosterone)
between_cors <- c(cor_matrix_r2[1, 3], cor_matrix_r2[1, 4],
                  cor_matrix_r2[2, 3], cor_matrix_r2[2, 4])
mean_between <- mean(between_cors)

# Overall mean (excluding diagonal)
all_cors <- cor_matrix_r2[upper.tri(cor_matrix_r2)]
mean_overall <- mean(all_cors)

# Convert to long format for ggplot
cor_long <- as.data.frame(cor_matrix_r2)
cor_long$Row <- rownames(cor_matrix_r2)
cor_long <- cor_long %>%
  pivot_longer(cols = -Row, names_to = "Column", values_to = "Correlation")

# Set factor levels for correct ordering
cor_long$Row <- factor(cor_long$Row, levels = rev(labels))
cor_long$Column <- factor(cor_long$Column, levels = labels)

# Quality metrics text
quality_text <- paste0(
  "Quality Metrics:\n",
  "Within Vehicle r\u00B2: ", sprintf("%.3f", vehicle_within), "\n",
  "Within Testosterone r\u00B2: ", sprintf("%.3f", testosterone_within), "\n",
  "Mean within-condition: ", sprintf("%.3f", mean_within), "\n",
  "Mean between-condition: ", sprintf("%.3f", mean_between), "\n",
  "Overall mean r\u00B2: ", sprintf("%.3f", mean_overall)
)

# Create heatmap
p <- ggplot(cor_long, aes(x = Column, y = Row, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 1) +

  # Add correlation values as text
  geom_text(aes(label = sprintf("%.2f", Correlation)),
            color = "white", fontface = "bold", size = 5) +

  # Color scale (YlOrRd - yellow to red)
  scale_fill_gradientn(
    colors = c("#FFFFB2", "#FECC5C", "#FD8D3C", "#F03B20", "#BD0026"),
    limits = c(0.7, 1.0),
    oob = scales::squish,
    name = expression("Pearson r"^2)
  ) +

  # Square aspect ratio
  coord_fixed() +

  # Labels
  labs(
    title = "Technical Replicate Correlation Analysis",
    subtitle = paste0("n = ", n_proteins, " proteins with complete data"),
    x = NULL,
    y = NULL
  ) +

  # Quality metrics annotation
  annotate("label",
           x = 4.7, y = 1.5,
           label = quality_text,
           hjust = 0, vjust = 0.5,
           size = 3,
           fill = "white",
           color = "black",
           label.padding = unit(0.4, "lines"),
           label.r = unit(0.1, "lines"),
           fontface = "plain") +

  # Theme
  theme_minimal(base_size = 11) +
  theme(
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),

    # Axes
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1,
                               face = "bold"),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),

    # Title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray30",
                                 margin = margin(b = 15)),

    # Legend
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.5, "cm"),

    # Margins
    plot.margin = margin(15, 100, 15, 15)
  )

# Save PNG
ggsave(
  filename = "correlation_heatmap.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "correlation_heatmap.pdf",
  plot = p,
  width = 10,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)

cat("Correlation heatmap saved to:\n")
cat("  - correlation_heatmap.png\n")
cat("  - correlation_heatmap.pdf\n")
cat("\nSummary:\n")
cat(paste0("  Proteins analyzed: ", n_proteins, "\n"))
cat(paste0("  Mean within-condition r²: ", sprintf("%.3f", mean_within), "\n"))
cat(paste0("  Mean between-condition r²: ", sprintf("%.3f", mean_between), "\n"))
cat(paste0("  Overall mean r²: ", sprintf("%.3f", mean_overall), "\n"))
