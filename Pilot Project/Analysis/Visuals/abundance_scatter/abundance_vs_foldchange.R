# ============================================================================
# Publication-grade Scatter Plot: Protein Abundance vs Fold Change
# With LOWESS smoothing and correlation analysis
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "cleaned_proteomics_data_with_QC_flags.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Filter out proteins where fold change cannot be calculated
df_filtered <- df %>%
  filter(FC_Type != "Cannot_Calculate") %>%
  filter(!is.na(Vehicle_Mean) & !is.na(Testosterone_Mean))

# Calculate average intensity
df_filtered <- df_filtered %>%
  mutate(
    avg_intensity = (Vehicle_Mean + Testosterone_Mean) / 2,
    log10_avg_intensity = log10(avg_intensity + 1)
  )

# Define point sizes based on confidence level
size_map <- c("High" = 5, "Medium" = 3.5, "Low" = 2)
df_filtered$point_size <- size_map[df_filtered$Confidence_Level]

# Define colors for functional class
color_palette <- c(
  "cytokine" = "#E74C3C",
  "growth factor" = "#2ECC71",
  "enzyme" = "#3498DB",
  "peptidase" = "#3498DB",
  "other" = "#95A5A6"
)

# Create display labels for functional class
df_filtered <- df_filtered %>%
  mutate(
    Functional_Class_Display = case_when(
      Functional_Class == "cytokine" ~ "Cytokine",
      Functional_Class == "growth factor" ~ "Growth Factor",
      Functional_Class %in% c("enzyme", "peptidase") ~ "Enzyme",
      TRUE ~ "Other"
    )
  )

# Proteins to label
proteins_to_label <- c("FRZB", "GAS6", "C4A/C4B", "CCL7", "CCL2", "TGFB3", "VCAN", "SLIT3")
df_filtered$label <- ifelse(df_filtered$Gene %in% proteins_to_label, df_filtered$Gene, "")

# Calculate Spearman correlation
cor_test <- cor.test(df_filtered$log10_avg_intensity, df_filtered$log_2_fold_change,
                      method = "spearman", exact = FALSE)
spearman_r <- round(cor_test$estimate, 3)
p_value <- cor_test$p.value
n_proteins <- nrow(df_filtered)

# Format p-value for display
if (p_value < 0.001) {
  p_display <- "p < 0.001"
} else {
  p_display <- paste0("p = ", format(round(p_value, 3), nsmall = 3))
}

# Correlation stats text
stats_text <- paste0(
  "Spearman \u03C1 = ", spearman_r, "\n",
  p_display, "\n",
  "n = ", n_proteins, " proteins"
)

# Create the plot
p <- ggplot(df_filtered, aes(x = log10_avg_intensity, y = log_2_fold_change)) +

  # Reference lines (behind everything)
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "gray50", linewidth = 0.6) +

  # LOWESS smoothing curve with confidence interval
  geom_smooth(method = "loess", se = TRUE,
              color = "#8B0000", fill = "#FFB6C1",
              linewidth = 2, alpha = 0.3, span = 0.75) +

  # Scatter points
  geom_point(aes(fill = Functional_Class_Display, size = Confidence_Level),
             shape = 21, color = "black", stroke = 0.8, alpha = 0.7) +

  # Labels for key proteins
  geom_text_repel(
    aes(label = label),
    size = 3.2,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    max.overlaps = 20,
    min.segment.length = 0,
    force = 3,
    arrow = arrow(length = unit(0.015, "npc"), type = "closed")
  ) +

  # Reference line labels
  annotate("text", x = 7.2, y = 0.25, label = "No change",
           size = 2.8, hjust = 0, fontface = "italic", color = "black") +
  annotate("text", x = 7.2, y = 1.25, label = "2-fold up",
           size = 2.5, hjust = 0, fontface = "italic", color = "gray40") +
  annotate("text", x = 7.2, y = -0.75, label = "2-fold down",
           size = 2.5, hjust = 0, fontface = "italic", color = "gray40") +
  annotate("text", x = 7.2, y = 2.25, label = "4-fold up",
           size = 2.5, hjust = 0, fontface = "italic", color = "gray40") +
  annotate("text", x = 7.2, y = -1.75, label = "4-fold down",
           size = 2.5, hjust = 0, fontface = "italic", color = "gray40") +

  # Correlation stats box
  annotate("label",
           x = max(df_filtered$log10_avg_intensity) - 0.3,
           y = max(df_filtered$log_2_fold_change) - 0.3,
           label = stats_text,
           hjust = 1, vjust = 1,
           size = 3.2,
           fill = "white",
           color = "black",
           label.padding = unit(0.4, "lines"),
           label.r = unit(0.15, "lines"),
           fontface = "plain") +

  # Scales
  scale_fill_manual(
    values = c("Cytokine" = "#E74C3C", "Growth Factor" = "#2ECC71",
               "Enzyme" = "#3498DB", "Other" = "#95A5A6"),
    name = "Functional Class"
  ) +
  scale_size_manual(
    values = c("High" = 5, "Medium" = 3.5, "Low" = 2),
    name = "Confidence Level"
  ) +
  scale_x_continuous(
    limits = c(6.5, 11),
    breaks = seq(7, 11, 1)
  ) +
  scale_y_continuous(
    limits = c(-4, 5),
    breaks = seq(-4, 5, 1)
  ) +

  # Labels
  labs(
    title = "Protein Abundance vs Differential Expression",
    subtitle = "Point size indicates data confidence level",
    x = expression("Mean Protein Intensity (Log"[10]*")"),
    y = expression("Log"[2]*" Fold Change (Testosterone/Vehicle)")
  ) +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Grid
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),

    # Axes
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),

    # Title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 15)),

    # Legend
    legend.position = "right",
    legend.box = "vertical",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines"),
    legend.margin = margin(5, 8, 5, 8),
    legend.spacing.y = unit(0.3, "cm"),

    # Margins
    plot.margin = margin(15, 15, 15, 15)
  ) +

  # Guide adjustments
  guides(
    fill = guide_legend(order = 1, override.aes = list(size = 4)),
    size = guide_legend(order = 2)
  )

# Save PNG
ggsave(
  filename = "abundance_vs_foldchange.png",
  plot = p,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "abundance_vs_foldchange.pdf",
  plot = p,
  width = 12,
  height = 9,
  device = cairo_pdf,
  bg = "white"
)

cat("Abundance vs fold change scatter plot saved to:\n")
cat("  - abundance_vs_foldchange.png\n")
cat("  - abundance_vs_foldchange.pdf\n")
cat("\nCorrelation Analysis:\n")
cat(paste0("  Spearman rho: ", spearman_r, "\n"))
cat(paste0("  P-value: ", format(p_value, scientific = TRUE, digits = 3), "\n"))
cat(paste0("  N proteins: ", n_proteins, "\n"))
