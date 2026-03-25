# ============================================================================
# Publication-grade Horizontal Bar Chart - Top 20 Protein Changes
# Top 10 Upregulated and Top 10 Downregulated Proteins
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "top_changers.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Separate upregulated and downregulated
upregulated <- df %>%
  filter(Regulation == "Upregulated") %>%
  arrange(desc(Fold_Change)) %>%
  head(10)

downregulated <- df %>%
  filter(Regulation == "Downregulated") %>%
  arrange(Fold_Change) %>%
  head(10)

# Define color palette for functional classes
color_palette <- c(
  "cytokine" = "#E74C3C",
  "growth factor" = "#2ECC71",
  "enzyme" = "#3498DB",
  "other" = "#95A5A6"
)

# Create fold change labels (format: "2.3x" or "0.3x")
upregulated$fc_label <- sprintf("%.1f\u00D7", upregulated$Fold_Change)
downregulated$fc_label <- sprintf("%.1f\u00D7", downregulated$Fold_Change)

# For very small values, use scientific notation
downregulated$fc_label[downregulated$Fold_Change == 0] <- "0\u00D7"
downregulated$fc_label[downregulated$Fold_Change > 0 & downregulated$Fold_Change < 0.1] <-
  sprintf("%.2f\u00D7", downregulated$Fold_Change[downregulated$Fold_Change > 0 & downregulated$Fold_Change < 0.1])

# Order genes by fold change for plotting
upregulated$Gene <- factor(upregulated$Gene, levels = rev(upregulated$Gene))
downregulated$Gene <- factor(downregulated$Gene, levels = rev(downregulated$Gene))

# Calculate label positions
upregulated$label_x <- upregulated$log_2_fold_change + 0.15
downregulated$label_x <- downregulated$log_2_fold_change - 0.15

# Create upregulated plot
p_up <- ggplot(upregulated, aes(x = log_2_fold_change, y = Gene, fill = Functional_Class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(aes(x = label_x, label = fc_label),
            hjust = 0, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = color_palette, name = "Functional Class",
                    labels = c("cytokine" = "Cytokine",
                              "growth factor" = "Growth Factor",
                              "enzyme" = "Enzyme",
                              "other" = "Other")) +
  scale_x_continuous(limits = c(0, max(upregulated$log_2_fold_change) * 1.25),
                     expand = c(0, 0)) +
  labs(title = "Top 10 Upregulated Proteins",
       x = NULL,
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 20, 5, 10)
  )

# Create downregulated plot
p_down <- ggplot(downregulated, aes(x = log_2_fold_change, y = Gene, fill = Functional_Class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
  geom_text(aes(x = label_x, label = fc_label),
            hjust = 1, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = color_palette, name = "Functional Class",
                    labels = c("cytokine" = "Cytokine",
                              "growth factor" = "Growth Factor",
                              "enzyme" = "Enzyme",
                              "other" = "Other")) +
  scale_x_continuous(limits = c(min(downregulated$log_2_fold_change) * 1.15, 0),
                     expand = c(0, 0)) +
  labs(title = "Top 10 Downregulated Proteins",
       x = expression("Log"[2]*" Fold Change (Testosterone/Vehicle)"),
       y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "bold", size = 10, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 20, 10, 10)
  )

# Create a separate legend
legend_data <- data.frame(
  Functional_Class = c("cytokine", "growth factor", "enzyme", "other"),
  x = 1:4,
  y = 1:4
)

p_legend <- ggplot(legend_data, aes(x = x, y = y, fill = Functional_Class)) +
  geom_point(shape = 22, size = 5, color = "black", stroke = 0.5) +
  scale_fill_manual(values = color_palette, name = "Functional Class",
                    labels = c("cytokine" = "Cytokine",
                              "growth factor" = "Growth Factor",
                              "enzyme" = "Enzyme",
                              "other" = "Other")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines")
  ) +
  guides(fill = guide_legend(nrow = 1))

# Extract legend
legend <- get_legend(p_legend)

# Combine plots
combined_plot <- plot_grid(
  p_up,
  p_down,
  legend,
  ncol = 1,
  rel_heights = c(1, 1, 0.15),
  align = "v",
  axis = "lr"
)

# Save PNG
ggsave(
  filename = "top_changers_barplot.png",
  plot = combined_plot,
  width = 10,
  height = 12,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "top_changers_barplot.pdf",
  plot = combined_plot,
  width = 10,
  height = 12,
  device = cairo_pdf,
  bg = "white"
)

cat("Top changers bar plot saved to:\n")
cat("  - top_changers_barplot.png\n")
cat("  - top_changers_barplot.pdf\n")
