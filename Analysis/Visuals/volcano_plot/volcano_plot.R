# ============================================================================
# Publication-grade Volcano Plot for SILAC Proteomics Data
# Testosterone vs Vehicle Control - Secreted Proteome
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

# Load data - use absolute path derived from script location
script_dir <- getwd()
data_path <- file.path(dirname(dirname(script_dir)), "General Analysis", "high_confidence_proteins.csv")

# If running from script directory
if (!file.exists(data_path)) {
  data_path <- file.path(script_dir, "..", "..", "General Analysis", "high_confidence_proteins.csv")
}
# Normalize path
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Calculate composite significance score
# Map confidence levels
confidence_map <- c("High" = 3, "Medium" = 2, "Low" = 1)
df$confidence_score <- confidence_map[df$Confidence_Level]

# Replicate completeness weight
df$completeness_weight <- ifelse(df$Replicate_Completeness == "Complete", 1.0, 0.5)

# Calculate average SD percent and inverse variability score
df$avg_sd_percent <- (df$Vehicle_SD_Percent + df$Testosterone_SD_Percent) / 2
max_sd <- max(df$avg_sd_percent, na.rm = TRUE)
df$variability_score <- 1 - (df$avg_sd_percent / max_sd)

# Composite significance score
df$significance_score <- df$confidence_score * df$completeness_weight * (1 + df$variability_score)

# Calculate -log10 of inverse significance (higher score = more significant)
df$neg_log_sig <- -log10(1 / (df$significance_score + 0.1))

# Calculate average intensity for point sizing
df$avg_intensity <- (df$Vehicle_Mean + df$Testosterone_Mean) / 2

# Normalize intensity for point sizes
min_size <- 1
max_size <- 8
df$point_size <- min_size + (max_size - min_size) * (
  (log10(df$avg_intensity + 1) - log10(min(df$avg_intensity) + 1)) /
    (log10(max(df$avg_intensity) + 1) - log10(min(df$avg_intensity) + 1))
)

# Color mapping by functional class
df$Functional_Class_Display <- case_when(
  df$Functional_Class == "cytokine" ~ "Cytokine",
  df$Functional_Class == "growth factor" ~ "Growth Factor",
  df$Functional_Class %in% c("enzyme", "peptidase") ~ "Enzyme",
  TRUE ~ "Other"
)

# Define color palette
color_palette <- c(
  "Cytokine" = "#E41A1C",
  "Growth Factor" = "#4DAF4A",
  "Enzyme" = "#377EB8",
  "Other" = "#808080"
)

# Proteins to annotate
upregulated <- c("FRZB", "GAS6", "C4A/C4B", "SLIT3", "VCAN")
downregulated <- c("CCL7", "CCL2", "TGFB3", "GDF15", "TGFB1")
proteins_to_label <- c(upregulated, downregulated)

# Add label column
df$label <- ifelse(df$Gene %in% proteins_to_label, df$Gene, "")

# Calculate significance threshold (median)
sig_threshold <- median(df$neg_log_sig, na.rm = TRUE)

# Create the plot
p <- ggplot(df, aes(x = log_2_fold_change, y = neg_log_sig)) +
  # Add threshold lines first (behind points)
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.5) +

  # Add points
  geom_point(
    aes(fill = Functional_Class_Display, size = point_size),
    shape = 21,
    color = "black",
    stroke = 0.5,
    alpha = 0.7
  ) +

  # Add labels with ggrepel
  geom_text_repel(
    aes(label = label),
    size = 3.2,
    fontface = "bold",
    family = "Arial",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    force = 2
  ) +

  # Set scales
  scale_fill_manual(
    values = color_palette,
    name = "Functional Class"
  ) +
  scale_size_identity() +  # Use actual point_size values

  # Set axis limits
  scale_x_continuous(
    limits = c(-3, 5),
    breaks = seq(-3, 5, 1),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    limits = c(0, max(df$neg_log_sig) * 1.15),
    expand = c(0.02, 0)
  ) +

  # Labels
  labs(
    title = "Testosterone vs Vehicle Control - Secreted Proteome",
    x = expression("Log"[2]*" Fold Change (Testosterone/Vehicle)"),
    y = "Composite Significance Score"
  ) +

  # Theme (Nature/Cell style)
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    # Plot background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Grid
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),

    # Axes
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12, face = "bold"),

    # Title
    plot.title = element_text(
      hjust = 0.5,
      size = 14,
      face = "bold",
      margin = margin(b = 15)
    ),

    # Legend
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines"),
    legend.margin = margin(5, 5, 5, 5),

    # Margins
    plot.margin = margin(15, 15, 15, 15)
  ) +

  # Guide for legend
  guides(
    fill = guide_legend(
      override.aes = list(size = 4)
    )
  )

# Save outputs
ggsave(
  filename = "volcano_plot.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "volcano_plot.pdf",
  plot = p,
  width = 10,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)

cat("Volcano plot saved to:\n")
cat("  - volcano_plot.png\n")
cat("  - volcano_plot.pdf\n")
