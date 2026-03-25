# ============================================================================
# Publication-grade Waterfall Plot - Complete Proteome Response
# All 42 proteins ranked by fold change
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(scales)
})

# ============================================================================
# Load and prepare data
# ============================================================================

script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "cleaned_proteomics_data_with_QC_flags.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Sort by fold change (ascending)
df <- df %>%
  arrange(log_2_fold_change) %>%
  mutate(Rank = row_number())

n_proteins <- nrow(df)

# ============================================================================
# Create color mapping based on fold change magnitude and direction
# ============================================================================

df <- df %>%
  mutate(
    # Determine regulation category
    Regulation = case_when(
      FC_Type == "Complete_Suppression" ~ "Complete Suppression",
      FC_Type == "Cannot_Calculate" ~ "Cannot Calculate",
      log_2_fold_change > 0.5 ~ "Upregulated",
      log_2_fold_change < -0.5 ~ "Downregulated",
      TRUE ~ "Unchanged"
    ),

    # Color intensity based on magnitude
    Color_Intensity = abs(log_2_fold_change) / max(abs(log_2_fold_change), na.rm = TRUE),

    # Alpha based on confidence level
    Alpha_Value = case_when(
      Confidence_Level == "High" ~ 1.0,
      Confidence_Level == "Medium" ~ 0.8,
      Confidence_Level == "Low" ~ 0.5,
      TRUE ~ 0.5
    ),

    # Bar color
    Bar_Color = case_when(
      FC_Type == "Complete_Suppression" ~ "#7D3C98",
      FC_Type == "Cannot_Calculate" ~ "#BDC3C7",
      log_2_fold_change > 0 ~ colorRampPalette(c("#FADBD8", "#E74C3C", "#922B21"))(100)[
        pmin(100, pmax(1, round(abs(log_2_fold_change) / max(abs(log_2_fold_change[log_2_fold_change > 0]), na.rm = TRUE) * 100)))
      ],
      log_2_fold_change < 0 ~ colorRampPalette(c("#D4E6F1", "#3498DB", "#1A5276"))(100)[
        pmin(100, pmax(1, round(abs(log_2_fold_change) / max(abs(log_2_fold_change[log_2_fold_change < 0]), na.rm = TRUE) * 100)))
      ],
      TRUE ~ "#95A5A6"
    ),

    # Functional class display
    Functional_Class_Display = case_when(
      Functional_Class == "cytokine" ~ "Cytokine",
      Functional_Class == "growth factor" ~ "Growth Factor",
      Functional_Class %in% c("enzyme", "peptidase") ~ "Enzyme",
      TRUE ~ "Other"
    )
  )

# Manual color assignment for each bar
up_colors <- colorRampPalette(c("#FADBD8", "#E74C3C", "#922B21"))(100)
down_colors <- colorRampPalette(c("#D4E6F1", "#3498DB", "#1A5276"))(100)

# Calculate color for each protein
df <- df %>%
  mutate(
    Bar_Fill = case_when(
      FC_Type == "Complete_Suppression" ~ "#7D3C98",
      FC_Type == "Cannot_Calculate" ~ "#BDC3C7",
      is.na(log_2_fold_change) ~ "#BDC3C7",
      log_2_fold_change > 0 ~ {
        max_up <- max(log_2_fold_change[log_2_fold_change > 0], na.rm = TRUE)
        idx <- pmin(100, pmax(1, round(log_2_fold_change / max_up * 100)))
        up_colors[idx]
      },
      log_2_fold_change < 0 ~ {
        max_down <- max(abs(log_2_fold_change[log_2_fold_change < 0]), na.rm = TRUE)
        idx <- pmin(100, pmax(1, round(abs(log_2_fold_change) / max_down * 100)))
        down_colors[idx]
      },
      TRUE ~ "#95A5A6"
    )
  )

# Recalculate fills properly
max_up <- max(df$log_2_fold_change[df$log_2_fold_change > 0 & !is.na(df$log_2_fold_change)], na.rm = TRUE)
max_down <- max(abs(df$log_2_fold_change[df$log_2_fold_change < 0 & !is.na(df$log_2_fold_change)]), na.rm = TRUE)

df$Bar_Fill <- sapply(1:nrow(df), function(i) {
  row <- df[i, ]
  if (row$FC_Type == "Complete_Suppression") return("#7D3C98")
  if (row$FC_Type == "Cannot_Calculate") return("#BDC3C7")
  if (is.na(row$log_2_fold_change)) return("#BDC3C7")
  if (row$log_2_fold_change > 0) {
    idx <- pmin(100, pmax(1, round(row$log_2_fold_change / max_up * 100)))
    return(up_colors[idx])
  }
  if (row$log_2_fold_change < 0) {
    idx <- pmin(100, pmax(1, round(abs(row$log_2_fold_change) / max_down * 100)))
    return(down_colors[idx])
  }
  return("#95A5A6")
})

# ============================================================================
# Calculate summary statistics
# ============================================================================

# Handle NA values for calculations
df_valid <- df %>% filter(!is.na(Fold_Change) & FC_Type != "Cannot_Calculate")

total_proteins <- nrow(df)
upregulated <- sum(df_valid$Fold_Change > 1, na.rm = TRUE)
downregulated <- sum(df_valid$Fold_Change < 1 & df_valid$Fold_Change > 0, na.rm = TRUE)
complete_supp <- sum(df$FC_Type == "Complete_Suppression")
cannot_calc <- sum(df$FC_Type == "Cannot_Calculate")
unchanged <- sum(df_valid$Fold_Change >= 0.8 & df_valid$Fold_Change <= 1.2, na.rm = TRUE)

fc_range_max <- max(df_valid$Fold_Change, na.rm = TRUE)
median_fc <- median(df_valid$Fold_Change[df_valid$Fold_Change > 0], na.rm = TRUE)

# ============================================================================
# Proteins to label
# ============================================================================

top_upregulated <- c("FRZB", "GAS6", "C4A/C4B", "SLIT3", "VCAN")
top_downregulated <- c("CCL7", "CCL2", "TGFB3", "GDF15", "TGFB1")
complete_suppression_genes <- c("EGF", "EREG", "PENK", "CXCL2")

proteins_to_label <- c(top_upregulated, top_downregulated, complete_suppression_genes)

df <- df %>%
  mutate(
    Label = ifelse(Gene %in% proteins_to_label, Gene, ""),
    Label_Y = ifelse(log_2_fold_change >= 0 | FC_Type == "Cannot_Calculate",
                     log_2_fold_change + 0.5,
                     log_2_fold_change - 0.5)
  )

# Fix label positions for NA/special cases
df$Label_Y[is.na(df$Label_Y)] <- 0.5

# ============================================================================
# Functional class colors
# ============================================================================

func_colors <- c(
  "Cytokine" = "#E74C3C",
  "Growth Factor" = "#2ECC71",
  "Enzyme" = "#3498DB",
  "Other" = "#95A5A6"
)

# ============================================================================
# Create the waterfall plot
# ============================================================================

# Y-axis limits
y_min <- min(df$log_2_fold_change, na.rm = TRUE) - 1.5
y_max <- max(df$log_2_fold_change, na.rm = TRUE) + 1.5

# Summary stats text
stats_text <- paste0(
  "SUMMARY STATISTICS\n",
  "Total proteins: ", total_proteins, "\n",
  "Upregulated: ", upregulated, " (", round(100*upregulated/total_proteins, 1), "%)\n",
  "Downregulated: ", downregulated, " (", round(100*downregulated/total_proteins, 1), "%)\n",
  "Complete suppression: ", complete_supp, "\n",
  "Cannot calculate: ", cannot_calc, "\n",
  "FC range: 0 to ", round(fc_range_max, 1), "\u00D7\n",
  "Median FC: ", round(median_fc, 2), "\u00D7"
)

main_plot <- ggplot(df, aes(x = Rank, y = log_2_fold_change)) +

  # Background shading for zones
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -1,
           fill = "#D4E6F1", alpha = 0.3) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -0.5, ymax = 0.5,
           fill = "gray90", alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "#FADBD8", alpha = 0.3) +

  # Reference lines
  geom_hline(yintercept = 0, color = "black", linewidth = 1.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray40", linewidth = 1) +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", color = "gray50", linewidth = 0.8) +

  # Bars
  geom_col(aes(alpha = Alpha_Value), fill = df$Bar_Fill,
           color = "black", linewidth = 0.3, width = 0.8) +

  # Functional class dots on top of bars
  geom_point(aes(y = ifelse(log_2_fold_change >= 0, log_2_fold_change + 0.3, log_2_fold_change - 0.3),
                 color = Functional_Class_Display),
             size = 2.5, shape = 19) +

  # Labels for key proteins
  geom_text(data = df %>% filter(Label != ""),
            aes(y = Label_Y, label = Label),
            angle = 90, hjust = ifelse(df$log_2_fold_change[df$Label != ""] >= 0, 0, 1),
            size = 2.5, fontface = "bold") +

  # Zone labels
  annotate("text", x = 3, y = y_min + 0.8, label = "Strong\nDownregulation",
           size = 3, fontface = "italic", color = "#1A5276", hjust = 0) +
  annotate("text", x = n_proteins - 2, y = y_max - 0.8, label = "Strong\nUpregulation",
           size = 3, fontface = "italic", color = "#922B21", hjust = 1) +

  # Reference line labels
  annotate("text", x = n_proteins + 1, y = 1, label = "2-fold", size = 2.5, hjust = 0, color = "gray40") +
  annotate("text", x = n_proteins + 1, y = -1, label = "2-fold", size = 2.5, hjust = 0, color = "gray40") +
  annotate("text", x = n_proteins + 1, y = 2, label = "4-fold", size = 2.5, hjust = 0, color = "gray50") +
  annotate("text", x = n_proteins + 1, y = -2, label = "4-fold", size = 2.5, hjust = 0, color = "gray50") +

  # Summary stats box
  annotate("label", x = n_proteins - 5, y = y_max - 0.3,
           label = stats_text,
           hjust = 1, vjust = 1, size = 2.8,
           fill = "white", color = "black",
           label.padding = unit(0.4, "lines"),
           fontface = "plain", lineheight = 1.1) +

  # Scales
  scale_alpha_identity() +
  scale_color_manual(values = func_colors, name = "Functional Class") +
  scale_x_continuous(breaks = seq(5, n_proteins, 5), expand = c(0.02, 0)) +
  scale_y_continuous(breaks = seq(-10, 5, 1),
                     limits = c(y_min, y_max)) +
  coord_cartesian(clip = "off") +

  # Labels
  labs(
    title = "Waterfall Plot - Complete Proteome Response to Testosterone",
    subtitle = "All 42 secreted proteins ranked by differential expression",
    x = "Protein Rank (Sorted by Fold Change)",
    y = expression("Log"[2]*" Fold Change (Testosterone/Vehicle)")
  ) +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40",
                                 margin = margin(b = 15)),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "gray80", linewidth = 0.5),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.margin = margin(5, 8, 5, 8),
    plot.margin = margin(15, 60, 15, 15)
  )

# ============================================================================
# Create legend for color scheme
# ============================================================================

# Custom legend explanation
legend_data <- data.frame(
  x = c(1, 2, 3, 4),
  y = c(1, 1, 1, 1),
  label = c("Upregulated", "Downregulated", "Complete\nSuppression", "Cannot\nCalculate"),
  color = c("#E74C3C", "#3498DB", "#7D3C98", "#BDC3C7")
)

color_legend <- ggplot(legend_data) +
  geom_tile(aes(x = x, y = y, fill = color), color = "black", width = 0.8, height = 0.5) +
  geom_text(aes(x = x, y = y - 0.5, label = label), size = 2.5, lineheight = 0.9) +
  scale_fill_identity() +
  labs(title = "Regulation") +
  theme_void() +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  coord_fixed(ratio = 1)

# ============================================================================
# Save outputs
# ============================================================================

ggsave(
  filename = "waterfall_plot.png",
  plot = main_plot,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "waterfall_plot.pdf",
  plot = main_plot,
  width = 16,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)

cat("Waterfall plot saved to:\n")
cat("  - waterfall_plot.png\n")
cat("  - waterfall_plot.pdf\n")
cat("\nSummary:\n")
cat(paste0("  Total proteins: ", total_proteins, "\n"))
cat(paste0("  Upregulated (FC > 1): ", upregulated, "\n"))
cat(paste0("  Downregulated (FC < 1): ", downregulated, "\n"))
cat(paste0("  Complete suppression: ", complete_supp, "\n"))
cat(paste0("  Cannot calculate: ", cannot_calc, "\n"))
cat(paste0("  Max fold change: ", round(fc_range_max, 2), "\n"))
cat(paste0("  Median fold change: ", round(median_fc, 2), "\n"))
