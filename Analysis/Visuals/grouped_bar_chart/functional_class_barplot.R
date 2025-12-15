# ============================================================================
# Publication-grade Grouped Bar Chart - Regulation Patterns by Functional Class
# Style: Nature/Cell publication standards
# ============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Load data
script_dir <- getwd()
data_path <- file.path(script_dir, "..", "..", "General Analysis", "functional_class_distribution.csv")
data_path <- normalizePath(data_path, mustWork = FALSE)

df <- read.csv(data_path, stringsAsFactors = FALSE)

# Order by total count descending
df <- df %>%
  arrange(desc(Count))

# Capitalize functional class names
df$Functional_Class_Display <- tools::toTitleCase(df$Functional_Class)

# Create factor with correct order
df$Functional_Class_Display <- factor(df$Functional_Class_Display,
                                       levels = df$Functional_Class_Display)

# Calculate totals for stats box
total_proteins <- sum(df$Count)
total_upregulated <- sum(df$Upregulated_Count)
total_downregulated <- sum(df$Downregulated_Count)

# Reshape data for grouped bar chart
df_long <- df %>%
  select(Functional_Class_Display, Upregulated_Count, Downregulated_Count, Percentage) %>%
  pivot_longer(
    cols = c(Upregulated_Count, Downregulated_Count),
    names_to = "Regulation",
    values_to = "Count"
  ) %>%
  mutate(
    Regulation = case_when(
      Regulation == "Upregulated_Count" ~ "Upregulated",
      Regulation == "Downregulated_Count" ~ "Downregulated"
    )
  )

# Set regulation factor order
df_long$Regulation <- factor(df_long$Regulation, levels = c("Upregulated", "Downregulated"))

# Define colors
colors <- c("Upregulated" = "#E74C3C", "Downregulated" = "#3498DB")

# Create percentage labels for x-axis
x_labels <- paste0(df$Functional_Class_Display, "\n(", sprintf("%.1f%%", df$Percentage), ")")
names(x_labels) <- df$Functional_Class_Display

# Stats box text
stats_text <- paste0(
  "Total proteins: ", total_proteins, "\n",
  "Total upregulated: ", total_upregulated, "\n",
  "Total downregulated: ", total_downregulated, "\n",
  "Pattern: Cytokines exclusively\ndownregulated"
)

# Create the plot
p <- ggplot(df_long, aes(x = Functional_Class_Display, y = Count, fill = Regulation)) +
  # Bars
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "black", linewidth = 0.5) +

  # Count labels on bars
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5, fontface = "bold") +

  # Add asterisk for cytokines (0 upregulated)
  annotate("text", x = 4, y = 0.3, label = "*",
           size = 6, fontface = "bold", color = "#E74C3C") +

  # Stats box
  annotate("label",
           x = 5.3, y = max(df_long$Count) * 0.95,
           label = stats_text,
           hjust = 1, vjust = 1,
           size = 3.2,
           fill = "white",
           color = "black",
           label.padding = unit(0.5, "lines"),
           label.r = unit(0.15, "lines"),
           fontface = "plain") +

  # Note about asterisk
  annotate("text",
           x = 5.3, y = 1.5,
           label = "* No upregulated proteins",
           hjust = 1, size = 3, fontface = "italic", color = "#666666") +

  # Scales
  scale_fill_manual(values = colors, name = NULL) +
  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     breaks = seq(0, 20, 4)) +

  # Labels
  labs(
    title = "Testosterone Regulation Patterns by Protein Functional Class",
    x = NULL,
    y = "Number of Proteins"
  ) +

  # Theme
  theme_classic(base_size = 11) +
  theme(
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),

    # Grid
    panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
    panel.grid.major.x = element_blank(),

    # Axes
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(size = 10, color = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),

    # Title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5,
                              margin = margin(b = 15)),

    # Legend
    legend.position = c(0.15, 0.9),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "lines"),
    legend.margin = margin(5, 8, 5, 8),

    # Margins
    plot.margin = margin(15, 20, 15, 15)
  )

# Save PNG
ggsave(
  filename = "functional_class_barplot.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Save PDF
ggsave(
  filename = "functional_class_barplot.pdf",
  plot = p,
  width = 12,
  height = 8,
  device = cairo_pdf,
  bg = "white"
)

cat("Functional class bar plot saved to:\n")
cat("  - functional_class_barplot.png\n")
cat("  - functional_class_barplot.pdf\n")
