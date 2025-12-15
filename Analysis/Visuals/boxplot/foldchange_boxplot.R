# ============================================================================
# Publication-grade Box Plot with Strip Plot Overlay
# Fold Change Distribution by Protein Function
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "high_confidence_proteins.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Filter for high-confidence proteins only
df_high <- df %>%
  filter(Confidence_Level == "High")

# Standardize functional class names and combine enzyme/peptidase
df_high <- df_high %>%
  mutate(
    Functional_Class_Display = case_when(
      Functional_Class == "cytokine" ~ "Cytokine",
      Functional_Class == "growth factor" ~ "Growth Factor",
      Functional_Class %in% c("enzyme", "peptidase") ~ "Enzyme",
      TRUE ~ "Other"
    )
  )

# Set factor order
df_high$Functional_Class_Display <- factor(
  df_high$Functional_Class_Display,
  levels = c("Cytokine", "Growth Factor", "Enzyme", "Other")
)

# Define color palette
color_palette <- c(
  "Cytokine" = "#E74C3C",
  "Growth Factor" = "#2ECC71",
  "Enzyme" = "#3498DB",
  "Other" = "#95A5A6"
)

# Calculate summary statistics per class
summary_stats <- df_high %>%
  group_by(Functional_Class_Display) %>%
  summarise(
    n = n(),
    mean_fc = round(mean(log_2_fold_change, na.rm = TRUE), 2),
    median_fc = round(median(log_2_fold_change, na.rm = TRUE), 2),
    min_fc = round(min(log_2_fold_change, na.rm = TRUE), 2),
    max_fc = round(max(log_2_fold_change, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    range_fc = paste0("[", min_fc, ", ", max_fc, "]")
  )

# Create labels for n= and median
label_data <- summary_stats %>%
  mutate(
    n_label = paste0("n=", n),
    med_label = paste0("med=", median_fc)
  )

# Y position for labels (above the max value or box)
y_max <- 5
label_data$y_n <- y_max - 0.3
label_data$y_med <- y_max - 0.8

# Create the main plot
p <- ggplot(df_high, aes(x = Functional_Class_Display, y = log_2_fold_change)) +

  # Reference lines (behind everything)
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.7) +

  # Box plot (gray boxes)
  geom_boxplot(
    width = 0.6,
    fill = "gray85",
    alpha = 0.6,
    outlier.shape = 18,  # Diamond
    outlier.size = 3,
    outlier.color = "black",
    coef = 1.5  # 1.5 IQR for whiskers
  ) +

  # Overlay with colored strip plot
  geom_jitter(
    aes(fill = Functional_Class_Display),
    shape = 21,
    size = 4,
    width = 0.2,
    alpha = 0.8,
    color = "black",
    stroke = 0.5
  ) +

  # Add median line enhancement (red thick line)
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    color = "#C0392B",
    linewidth = 0.8,
    fatten = 0
  ) +

  # n= labels
  geom_text(
    data = label_data,
    aes(x = Functional_Class_Display, y = y_n, label = n_label),
    size = 3.5,
    fontface = "bold"
  ) +

  # Median labels
  geom_text(
    data = label_data,
    aes(x = Functional_Class_Display, y = y_med, label = med_label),
    size = 3,
    fontface = "italic",
    color = "#C0392B"
  ) +

  # Reference line labels
  annotate("text", x = 4.4, y = 0.15, label = "No change",
           size = 2.8, hjust = 0, fontface = "italic", color = "black") +
  annotate("text", x = 4.4, y = 1.15, label = "2-fold up",
           size = 2.8, hjust = 0, fontface = "italic", color = "gray50") +
  annotate("text", x = 4.4, y = -0.85, label = "2-fold down",
           size = 2.8, hjust = 0, fontface = "italic", color = "gray50") +

  # Scales
  scale_fill_manual(values = color_palette, guide = "none") +
  scale_y_continuous(
    limits = c(-3, 5),
    breaks = seq(-3, 5, 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +

  # Labels
  labs(
    title = "Fold Change Distribution by Protein Function",
    x = "Functional Class",
    y = expression("Log"[2]*" Fold Change (Testosterone/Vehicle)")
  ) +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Grid
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),

    # Axes
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),

    # Title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 15)),

    # Margins
    plot.margin = margin(15, 50, 80, 15)
  )

# Create summary table
table_data <- summary_stats %>%
  select(Functional_Class_Display, n, mean_fc, median_fc, range_fc) %>%
  rename(
    `Class` = Functional_Class_Display,
    `N` = n,
    `Mean` = mean_fc,
    `Median` = median_fc,
    `Range` = range_fc
  )

# Create table grob
table_theme <- ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 9, fontface = "plain"),
    bg_params = list(fill = c("gray95", "white"))
  ),
  colhead = list(
    fg_params = list(fontsize = 9, fontface = "bold"),
    bg_params = list(fill = "gray80")
  )
)

table_grob <- tableGrob(table_data, rows = NULL, theme = table_theme)

# Add title to table
table_title <- textGrob("Summary Statistics",
                        gp = gpar(fontsize = 10, fontface = "bold"))

# Combine plot and table
combined <- arrangeGrob(
  p,
  table_title,
  table_grob,
  ncol = 1,
  heights = c(6, 0.3, 1.2)
)

# Save PNG
ggsave(
  filename = "foldchange_boxplot.png",
  plot = combined,
  width = 10,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "foldchange_boxplot.pdf",
  plot = combined,
  width = 10,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

cat("Fold change box plot saved to:\n")
cat("  - foldchange_boxplot.png\n")
cat("  - foldchange_boxplot.pdf\n")
