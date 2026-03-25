# ============================================================================
# Publication-grade Multi-Panel Data Quality Figure
# 2x2 grid: Pie, Bar, Donut, Table
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
# Load and process data
# ============================================================================

script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "cleaned_proteomics_data_with_QC_flags.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# ============================================================================
# PANEL A: Confidence Level Distribution - PIE CHART
# ============================================================================

confidence_data <- df %>%
  count(Confidence_Level) %>%
  mutate(
    Percentage = round(100 * n / sum(n), 1),
    Label = paste0(n, " (", Percentage, "%)")
  )

# Reorder levels
confidence_data$Confidence_Level <- factor(
  confidence_data$Confidence_Level,
  levels = c("High", "Medium", "Low")
)
confidence_data <- confidence_data %>% arrange(Confidence_Level)

# Colors
conf_colors <- c("High" = "#2ECC71", "Medium" = "#F39C12", "Low" = "#E74C3C")

# Calculate positions for labels
confidence_data <- confidence_data %>%
  arrange(desc(Confidence_Level)) %>%
  mutate(
    ypos = cumsum(n) - 0.5 * n,
    pct_label = paste0(Percentage, "%")
  )

panel_a <- ggplot(confidence_data, aes(x = "", y = n, fill = Confidence_Level)) +

geom_col(width = 1, color = "white", linewidth = 1.5) +
  coord_polar(theta = "y", start = pi/2) +
  scale_fill_manual(values = conf_colors, name = "Confidence") +
  geom_text(aes(y = ypos, label = pct_label),
            color = "white", fontface = "bold", size = 4.5) +
  labs(title = "A. Data Confidence Distribution") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

# ============================================================================
# PANEL B: Replicate Completeness - HORIZONTAL BAR CHART
# ============================================================================

completeness_data <- df %>%
  count(Replicate_Completeness) %>%
  mutate(
    Percentage = round(100 * n / sum(n), 1),
    Pct_Label = paste0("(", Percentage, "%)")
  )

# Rename for display
completeness_data <- completeness_data %>%
  mutate(
    Display_Name = case_when(
      Replicate_Completeness == "Complete" ~ "Complete (4/4)",
      Replicate_Completeness == "Missing_1" ~ "Missing 1 (3/4)",
      Replicate_Completeness == "Missing_2+" ~ "Missing 2+ (\u22642/4)",
      TRUE ~ Replicate_Completeness
    )
  )

# Set order
completeness_data$Display_Name <- factor(
  completeness_data$Display_Name,
  levels = rev(c("Complete (4/4)", "Missing 1 (3/4)", "Missing 2+ (\u22642/4)"))
)

# Colors
comp_colors <- c(
  "Complete (4/4)" = "#2ECC71",
  "Missing 1 (3/4)" = "#F1C40F",
  "Missing 2+ (\u22642/4)" = "#E74C3C"
)

panel_b <- ggplot(completeness_data, aes(x = n, y = Display_Name, fill = Display_Name)) +
  geom_col(color = "black", linewidth = 0.6, width = 0.6) +
  geom_text(aes(label = n, x = n - 1.5),
            color = "white", fontface = "bold", size = 4, hjust = 1) +
  geom_text(aes(label = Pct_Label, x = n + 0.8),
            color = "black", size = 3.5, hjust = 0) +
  scale_fill_manual(values = comp_colors, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15)), limits = c(0, 32)) +
  labs(
    title = "B. Replicate Data Completeness",
    x = "Number of Proteins",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 8)),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )

# ============================================================================
# PANEL C: FC Type Distribution - DONUT CHART
# ============================================================================

fc_type_data <- df %>%
  count(FC_Type) %>%
  mutate(
    Percentage = round(100 * n / sum(n), 1),
    Label = paste0(n, "\n(", Percentage, "%)")
  )

# Rename for display
fc_type_data <- fc_type_data %>%
  mutate(
    Display_Name = case_when(
      FC_Type == "Normal" ~ "Normal",
      FC_Type == "Complete_Suppression" ~ "Complete\nSuppression",
      FC_Type == "Cannot_Calculate" ~ "Cannot\nCalculate",
      TRUE ~ FC_Type
    )
  )

# Colors
fc_colors <- c(
  "Normal" = "#3498DB",
  "Complete\nSuppression" = "#9B59B6",
  "Cannot\nCalculate" = "#E67E22"
)

# Calculate positions for donut
fc_type_data <- fc_type_data %>%
  mutate(
    fraction = n / sum(n),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2
  )

panel_c <- ggplot(fc_type_data, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.4, fill = Display_Name)) +
  geom_rect(color = "white", linewidth = 1) +
  geom_text(aes(x = 3.2, y = labelPosition, label = Label),
            color = "white", fontface = "bold", size = 3.2, lineheight = 0.9) +
  # Center text
  annotate("text", x = 0, y = 0, label = "42\nTotal\nProteins",
           size = 5, fontface = "bold", lineheight = 0.9) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  scale_fill_manual(values = fc_colors, name = "FC Type") +
  labs(title = "C. Fold Change Calculation Status") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

# ============================================================================
# PANEL D: Summary Statistics TABLE
# ============================================================================

# Calculate summary stats
total_proteins <- nrow(df)
high_conf <- sum(df$Confidence_Level == "High")
high_conf_pct <- round(100 * high_conf / total_proteins, 1)

# Calculate upregulated/downregulated from the data
df_with_fc <- df %>% filter(FC_Type != "Cannot_Calculate" & !is.na(Fold_Change))
upregulated <- sum(df_with_fc$Fold_Change > 1, na.rm = TRUE)
downregulated <- sum(df_with_fc$Fold_Change < 1, na.rm = TRUE)
complete_supp <- sum(df$FC_Type == "Complete_Suppression")

mean_fc <- round(mean(df_with_fc$Fold_Change[df_with_fc$Fold_Change > 0], na.rm = TRUE), 2)
median_fc <- round(median(df_with_fc$Fold_Change[df_with_fc$Fold_Change > 0], na.rm = TRUE), 2)

# Create table data
summary_table <- data.frame(
  Metric = c(
    "Total Proteins",
    "High Confidence",
    "Usable for Analysis",
    "Upregulated",
    "Downregulated",
    "Complete Suppression",
    "Mean Fold Change",
    "Median Fold Change"
  ),
  Value = c(
    as.character(total_proteins),
    paste0(high_conf, " (", high_conf_pct, "%)"),
    as.character(high_conf),
    paste0(upregulated, " (", round(100*upregulated/total_proteins, 1), "%)"),
    paste0(downregulated, " (", round(100*downregulated/total_proteins, 1), "%)"),
    paste0(complete_supp, " (", round(100*complete_supp/total_proteins, 1), "%)"),
    paste0(mean_fc, "\u00D7"),
    paste0(median_fc, "\u00D7")
  ),
  stringsAsFactors = FALSE
)

# Create table grob with professional styling
table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 10, fontface = "plain", hjust = c(0, 1), x = c(0.05, 0.95)),
    bg_params = list(fill = c("gray95", "white"), col = "gray70", lwd = 0.5)
  ),
  colhead = list(
    fg_params = list(fontsize = 11, fontface = "bold", hjust = c(0, 1), x = c(0.05, 0.95)),
    bg_params = list(fill = "gray80", col = "gray50", lwd = 1)
  )
)

table_grob <- tableGrob(summary_table, rows = NULL, theme = table_theme)

# Add border
table_grob <- gtable::gtable_add_grob(
  table_grob,
  grobs = rectGrob(gp = gpar(fill = NA, lwd = 2, col = "gray50")),
  t = 1, b = nrow(table_grob), l = 1, r = ncol(table_grob)
)

# Wrap in ggplot
panel_d <- ggplot() +
  annotation_custom(table_grob) +
  labs(title = "D. Dataset Summary Statistics") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
    plot.margin = margin(10, 10, 10, 10)
  )

# ============================================================================
# Combine all panels
# ============================================================================

# Create main title
main_title <- ggdraw() +
  draw_label(
    "Data Quality Assessment - SILAC Proteomics Analysis",
    fontface = "bold", size = 18, x = 0.5, hjust = 0.5
  )

# Create subtitle
subtitle <- ggdraw() +
  draw_label(
    "Testosterone vs Vehicle Control (n=42 secreted proteins)",
    size = 12, x = 0.5, hjust = 0.5, color = "gray40"
  )

# Create footer
footer <- ggdraw() +
  draw_label(
    "Quality control metrics for downstream analysis selection",
    size = 10, x = 0.5, hjust = 0.5, fontface = "italic", color = "gray50"
  )

# Combine panels in 2x2 grid
panel_grid <- plot_grid(
  panel_a, panel_b,
  panel_c, panel_d,
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1.2),
  rel_heights = c(1, 1),
  align = "hv"
)

# Final assembly
final_plot <- plot_grid(
  main_title,
  subtitle,
  panel_grid,
  footer,
  ncol = 1,
  rel_heights = c(0.05, 0.03, 1, 0.03)
)

# Save PNG
ggsave(
  filename = "quality_multipanel.png",
  plot = final_plot,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "quality_multipanel.pdf",
  plot = final_plot,
  width = 16,
  height = 12,
  device = cairo_pdf,
  bg = "white"
)

cat("Multi-panel quality figure saved to:\n")
cat("  - quality_multipanel.png\n")
cat("  - quality_multipanel.pdf\n")
