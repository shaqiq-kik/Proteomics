# ============================================================================
# Publication-grade Replicate Correlation Analysis Figure
# High-confidence proteins only
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(grid)
  library(gridExtra)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "high_confidence_proteins.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Filter for high-confidence proteins only
df_high <- df %>%
  filter(Confidence_Level == "High")

# Select intensity columns
intensity_cols <- c("Vehicle_Rep1_31579", "Vehicle_Rep2_31581",
                    "Testosterone_Rep1_31578", "Testosterone_Rep2_31580")

# Extract and filter complete data
intensity_data <- df_high[, intensity_cols]
complete_mask <- complete.cases(intensity_data) &
                 rowSums(intensity_data == 0, na.rm = TRUE) == 0
intensity_complete <- intensity_data[complete_mask, ]

n_proteins <- nrow(intensity_complete)

# Log10 transform
log_intensity <- log10(intensity_complete + 1)

# Rename columns for display
display_names <- c("Vehicle\nRep 1", "Vehicle\nRep 2",
                   "Testosterone\nRep 1", "Testosterone\nRep 2")
colnames(log_intensity) <- display_names

# ============================================================================
# Calculate correlations with p-values
# ============================================================================

n_cols <- ncol(log_intensity)
cor_matrix <- matrix(0, nrow = n_cols, ncol = n_cols)
r2_matrix <- matrix(0, nrow = n_cols, ncol = n_cols)
p_matrix <- matrix(0, nrow = n_cols, ncol = n_cols)

for (i in 1:n_cols) {
  for (j in 1:n_cols) {
    if (i == j) {
      cor_matrix[i, j] <- 1
      r2_matrix[i, j] <- 1
      p_matrix[i, j] <- 0
    } else {
      test <- cor.test(log_intensity[[i]], log_intensity[[j]], method = "pearson")
      cor_matrix[i, j] <- test$estimate
      r2_matrix[i, j] <- test$estimate^2
      p_matrix[i, j] <- test$p.value
    }
  }
}

rownames(r2_matrix) <- display_names
colnames(r2_matrix) <- display_names
rownames(p_matrix) <- display_names
colnames(p_matrix) <- display_names

# ============================================================================
# Create significance labels
# ============================================================================

sig_labels <- matrix("", nrow = n_cols, ncol = n_cols)
for (i in 1:n_cols) {
  for (j in 1:n_cols) {
    if (i != j) {
      p <- p_matrix[i, j]
      if (p < 0.001) {
        sig_labels[i, j] <- "***"
      } else if (p < 0.01) {
        sig_labels[i, j] <- "**"
      } else if (p < 0.05) {
        sig_labels[i, j] <- "*"
      }
    }
  }
}

# ============================================================================
# Calculate summary statistics
# ============================================================================

# Within-condition correlations
vehicle_r2 <- r2_matrix[1, 2]
testosterone_r2 <- r2_matrix[3, 4]
mean_within <- (vehicle_r2 + testosterone_r2) / 2

# Between-condition correlations
between_cors <- c(r2_matrix[1, 3], r2_matrix[1, 4],
                  r2_matrix[2, 3], r2_matrix[2, 4])
mean_between <- mean(between_cors)

# Overall mean (excluding diagonal)
all_cors <- r2_matrix[upper.tri(r2_matrix)]
mean_overall <- mean(all_cors)

# Quality assessment
quality_score <- mean_within
if (quality_score > 0.95) {
  quality_level <- "Excellent"
  quality_color <- "#2ECC71"
} else if (quality_score > 0.90) {
  quality_level <- "Good"
  quality_color <- "#F1C40F"
} else if (quality_score > 0.85) {
  quality_level <- "Acceptable"
  quality_color <- "#E67E22"
} else {
  quality_level <- "Poor"
  quality_color <- "#E74C3C"
}

# ============================================================================
# Convert to long format for ggplot heatmap
# ============================================================================

r2_long <- as.data.frame(r2_matrix)
r2_long$Row <- rownames(r2_matrix)
r2_long <- r2_long %>%
  pivot_longer(cols = -Row, names_to = "Column", values_to = "r2")

# Add significance labels
sig_long <- as.data.frame(sig_labels)
colnames(sig_long) <- display_names
sig_long$Row <- display_names
sig_long <- sig_long %>%
  pivot_longer(cols = -Row, names_to = "Column", values_to = "sig")

# Merge
r2_long <- r2_long %>%
  left_join(sig_long, by = c("Row", "Column"))

# Create combined label
r2_long <- r2_long %>%
  mutate(
    label = ifelse(Row == Column, "1.00", sprintf("%.2f", r2)),
    full_label = ifelse(sig == "", label, paste0(label, "\n", sig))
  )

# Set factor levels for correct ordering
r2_long$Row <- factor(r2_long$Row, levels = rev(display_names))
r2_long$Column <- factor(r2_long$Column, levels = display_names)

# ============================================================================
# Create main heatmap
# ============================================================================

heatmap_plot <- ggplot(r2_long, aes(x = Column, y = Row, fill = r2)) +
  geom_tile(color = "white", linewidth = 2) +

  # Add correlation values
  geom_text(aes(label = full_label),
            color = "black", fontface = "bold", size = 4.5, lineheight = 0.85) +

  # Color scale (RdYlGn)
  scale_fill_gradientn(
    colors = c("#D73027", "#FC8D59", "#FEE08B", "#D9EF8B", "#91CF60", "#1A9850"),
    limits = c(0.5, 1.0),
    oob = scales::squish,
    name = expression("Pearson r"^2)
  ) +

  # Square aspect
  coord_fixed() +

  # Labels
  labs(
    title = "Technical Replicate Correlation Analysis",
    subtitle = paste0("High-confidence proteins only (n=", n_proteins, ")"),
    x = NULL,
    y = NULL
  ) +

  # Theme
  theme_minimal(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, color = "black", face = "bold",
                               lineheight = 0.9),
    axis.text.y = element_text(size = 10, color = "black", face = "bold",
                               lineheight = 0.9),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.5, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================================
# Create summary statistics panel
# ============================================================================

summary_text <- paste0(
  "CORRELATION SUMMARY\n",
  "\n",
  "Within-Condition:\n",
  "  Vehicle r\u00B2: ", sprintf("%.3f", vehicle_r2), "\n",
  "  Testosterone r\u00B2: ", sprintf("%.3f", testosterone_r2), "\n",
  "  Mean within: ", sprintf("%.3f", mean_within), "\n",
  "\n",
  "Between-Condition:\n",
  "  Mean r\u00B2: ", sprintf("%.3f", mean_between), "\n",
  "\n",
  "Overall:\n",
  "  Mean r\u00B2: ", sprintf("%.3f", mean_overall), "\n",
  "\n",
  "Proteins: n = ", n_proteins
)

summary_panel <- ggplot() +
  annotate("label",
           x = 0.5, y = 0.65,
           label = summary_text,
           hjust = 0.5, vjust = 0.5,
           size = 3.5,
           fill = "gray95",
           color = "black",
           label.padding = unit(0.6, "lines"),
           label.r = unit(0.15, "lines"),
           fontface = "plain",
           lineheight = 1.1) +

  # Quality indicator box
  annotate("rect",
           xmin = 0.15, xmax = 0.85,
           ymin = 0.15, ymax = 0.35,
           fill = quality_color, color = "black", linewidth = 1) +
  annotate("text",
           x = 0.5, y = 0.25,
           label = paste0("Quality: ", quality_level),
           size = 4.5, fontface = "bold", color = "white") +

  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================================
# Create footer with statistical note
# ============================================================================

footer <- ggdraw() +
  draw_label(
    "All pairwise correlations significant at p < 0.001; Log\u2081\u2080-transformed intensities | *** p<0.001, ** p<0.01, * p<0.05",
    size = 9, x = 0.5, hjust = 0.5, fontface = "italic", color = "gray50"
  )

# ============================================================================
# Combine all elements
# ============================================================================

# Main content: heatmap + summary panel
main_content <- plot_grid(
  heatmap_plot,
  summary_panel,
  ncol = 2,
  rel_widths = c(3, 1.2),
  align = "h"
)

# Final assembly with footer
final_plot <- plot_grid(
  main_content,
  footer,
  ncol = 1,
  rel_heights = c(1, 0.05)
)

# ============================================================================
# Save outputs
# ============================================================================

ggsave(
  filename = "replicate_correlation.png",
  plot = final_plot,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "replicate_correlation.pdf",
  plot = final_plot,
  width = 12,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

cat("Replicate correlation analysis figure saved to:\n")
cat("  - replicate_correlation.png\n")
cat("  - replicate_correlation.pdf\n")
cat("\nSummary Statistics:\n")
cat(paste0("  Proteins analyzed: ", n_proteins, "\n"))
cat(paste0("  Vehicle within r\u00B2: ", sprintf("%.3f", vehicle_r2), "\n"))
cat(paste0("  Testosterone within r\u00B2: ", sprintf("%.3f", testosterone_r2), "\n"))
cat(paste0("  Mean within-condition r\u00B2: ", sprintf("%.3f", mean_within), "\n"))
cat(paste0("  Mean between-condition r\u00B2: ", sprintf("%.3f", mean_between), "\n"))
cat(paste0("  Overall quality: ", quality_level, "\n"))
