# ============================================================================
# Publication-grade Data Quality Summary Figure
# Multi-panel figure (2x2 grid)
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(gridExtra)
  library(grid)
})

# ============================================================================
# Data for plots (from data_quality_report.txt)
# ============================================================================

# Panel A: Confidence Level Distribution
confidence_data <- data.frame(
  Category = c("High", "Medium", "Low"),
  Count = c(25, 8, 9),
  stringsAsFactors = FALSE
)
confidence_data$Percentage <- round(100 * confidence_data$Count / sum(confidence_data$Count), 1)
confidence_data$Label <- paste0(confidence_data$Category, "\n", confidence_data$Count,
                                 " (", confidence_data$Percentage, "%)")
confidence_data$Category <- factor(confidence_data$Category, levels = c("High", "Medium", "Low"))

# Panel B: Replicate Completeness
completeness_data <- data.frame(
  Category = c("Complete (4/4)", "Missing 1 (3/4)", "Missing 2+ (\u22642/4)"),
  Count = c(25, 8, 9),
  stringsAsFactors = FALSE
)
completeness_data$Category <- factor(completeness_data$Category,
                                      levels = rev(completeness_data$Category))

# Panel C: Fold Change Type Distribution
fc_type_data <- data.frame(
  Category = c("Normal", "Complete Suppression", "Cannot Calculate"),
  Count = c(34, 4, 4),
  stringsAsFactors = FALSE
)
fc_type_data$Percentage <- round(100 * fc_type_data$Count / sum(fc_type_data$Count), 1)
fc_type_data$Category <- factor(fc_type_data$Category,
                                 levels = c("Normal", "Complete Suppression", "Cannot Calculate"))

# ============================================================================
# Panel A: Pie Chart - Confidence Level Distribution
# ============================================================================

# Calculate positions for pie chart labels
confidence_data <- confidence_data %>%
  arrange(desc(Category)) %>%
  mutate(
    ypos = cumsum(Count) - 0.5 * Count
  )

# Colors
confidence_colors <- c("High" = "#2ECC71", "Medium" = "#F39C12", "Low" = "#E74C3C")

panel_a <- ggplot(confidence_data, aes(x = "", y = Count, fill = Category)) +
  geom_col(width = 1, color = "white", linewidth = 1) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = confidence_colors, name = "Confidence") +
  geom_text(aes(y = ypos, label = paste0(Count, "\n(", Percentage, "%)")),
            color = "white", fontface = "bold", size = 3.5) +
  labs(title = "Confidence Level Distribution") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# Panel B: Horizontal Bar Chart - Replicate Completeness
# ============================================================================

completeness_colors <- c("Complete (4/4)" = "#2ECC71",
                          "Missing 1 (3/4)" = "#F1C40F",
                          "Missing 2+ (\u22642/4)" = "#E74C3C")

panel_b <- ggplot(completeness_data, aes(x = Count, y = Category, fill = Category)) +
  geom_col(color = "black", linewidth = 0.5, width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.3, fontface = "bold", size = 4) +
  scale_fill_manual(values = completeness_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 30)) +
  labs(
    title = "Replicate Completeness",
    x = "Number of Proteins",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
    axis.text.y = element_text(face = "bold", size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 5)),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    plot.margin = margin(5, 15, 5, 5)
  )

# ============================================================================
# Panel C: Donut Chart - Fold Change Type Distribution
# ============================================================================

fc_colors <- c("Normal" = "#3498DB",
               "Complete Suppression" = "#9B59B6",
               "Cannot Calculate" = "#95A5A6")

# Calculate positions
fc_type_data <- fc_type_data %>%
  mutate(
    fraction = Count / sum(Count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2,
    label = paste0(Count, " (", Percentage, "%)")
  )

panel_c <- ggplot(fc_type_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = Category)) +
  geom_rect(color = "white", linewidth = 1) +
  geom_text(aes(x = 3.25, y = labelPosition, label = label),
            color = "white", fontface = "bold", size = 3) +
  coord_polar(theta = "y") +
  xlim(c(1, 4)) +
  scale_fill_manual(values = fc_colors, name = "FC Type") +
  labs(title = "Fold Change Type Distribution") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# Panel D: Summary Statistics Table
# ============================================================================

summary_table <- data.frame(
  Metric = c("Total proteins analyzed",
             "High confidence (usable)",
             "Upregulated proteins",
             "Downregulated proteins",
             "Complete suppression",
             "Overall data completeness"),
  Value = c("42",
            "25 (59.5%)",
            "21",
            "17",
            "4",
            "94.4%"),
  stringsAsFactors = FALSE
)

# Create table grob
table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 10, fontface = "plain", hjust = 0, x = 0.05),
    bg_params = list(fill = c("gray95", "white"))
  ),
  colhead = list(
    fg_params = list(fontsize = 11, fontface = "bold", hjust = 0, x = 0.05),
    bg_params = list(fill = "gray80")
  )
)

table_grob <- tableGrob(summary_table, rows = NULL, theme = table_theme)

# Wrap in ggplot for consistent sizing
panel_d <- ggplot() +
  annotation_custom(table_grob) +
  labs(title = "Summary Statistics") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5, margin = margin(b = 10)),
    plot.margin = margin(5, 5, 5, 5)
  )

# ============================================================================
# Combine all panels
# ============================================================================

# Add panel labels
panel_a_labeled <- panel_a + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 14))
panel_b_labeled <- panel_b + labs(tag = "B") +
  theme(plot.tag = element_text(face = "bold", size = 14))
panel_c_labeled <- panel_c + labs(tag = "C") +
  theme(plot.tag = element_text(face = "bold", size = 14))
panel_d_labeled <- panel_d + labs(tag = "D") +
  theme(plot.tag = element_text(face = "bold", size = 14))

# Create main title
main_title <- ggdraw() +
  draw_label("Data Quality Assessment - SILAC Proteomics Analysis",
             fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

# Combine panels in 2x2 grid
panel_grid <- plot_grid(
  panel_a_labeled, panel_b_labeled,
  panel_c_labeled, panel_d_labeled,
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  rel_heights = c(1, 1),
  align = "hv"
)

# Add main title
final_plot <- plot_grid(
  main_title,
  panel_grid,
  ncol = 1,
  rel_heights = c(0.06, 1)
)

# Save PNG
ggsave(
  filename = "data_quality_summary.png",
  plot = final_plot,
  width = 14,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "data_quality_summary.pdf",
  plot = final_plot,
  width = 14,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

cat("Data quality summary figure saved to:\n")
cat("  - data_quality_summary.png\n")
cat("  - data_quality_summary.pdf\n")
