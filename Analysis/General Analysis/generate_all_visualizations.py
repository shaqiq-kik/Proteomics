"""
COMPREHENSIVE PROTEOMICS VISUALIZATION SUITE
Generates 12 publication-quality figures for proteomics analysis
Testosterone vs Vehicle Control SILAC experiment
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import pearsonr, zscore, gaussian_kde
from adjustText import adjust_text
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import warnings
import os
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

print("="*80)
print("COMPREHENSIVE PROTEOMICS VISUALIZATION SUITE")
print("="*80)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Color scheme
COLOR_UP = '#d62728'      # Red for upregulated
COLOR_DOWN = '#1f77b4'    # Blue for downregulated
COLOR_NS = '#7f7f7f'      # Gray for unchanged

# Thresholds
FC_THRESHOLD = 1.5        # 1.5-fold change cutoff
LOG2FC_THRESHOLD = np.log2(FC_THRESHOLD)

# Figure settings
FIGURE_DPI = 300
FIGURE_SIZE_STANDARD = (10, 8)
FIGURE_SIZE_WIDE = (14, 8)
FIGURE_SIZE_TALL = (10, 12)

# Font sizes
TITLE_SIZE = 16
LABEL_SIZE = 14
TICK_SIZE = 12
LEGEND_SIZE = 11

print("\nConfiguration:")
print(f"  Fold change threshold: {FC_THRESHOLD}")
print(f"  Log2 FC threshold: ±{LOG2FC_THRESHOLD:.3f}")
print(f"  Figure DPI: {FIGURE_DPI}")
print(f"  Color scheme: Up={COLOR_UP}, Down={COLOR_DOWN}, NS={COLOR_NS}")

# ============================================================================
# CREATE FOLDER STRUCTURE
# ============================================================================

base_dir = "../Visuals"
folders = {
    'volcano': os.path.join(base_dir, 'Volcano Plot'),
    'heatmap': os.path.join(base_dir, 'Heatmaps'),
    'bars': os.path.join(base_dir, 'Bar Charts'),
    'distributions': os.path.join(base_dir, 'Distributions'),
    'correlations': os.path.join(base_dir, 'Correlations'),
    'qc': os.path.join(base_dir, 'QC Plots'),
    'publication': os.path.join(base_dir, 'Publication Figures')
}

for folder in folders.values():
    os.makedirs(folder, exist_ok=True)

print("\nFolder structure created in Visuals/")

# ============================================================================
# LOAD DATA
# ============================================================================

print("\nLoading data files...")
df = pd.read_csv('cleaned_proteomics_data_with_QC_flags.csv')
df_high = pd.read_csv('high_confidence_proteins.csv')
func_dist = pd.read_csv('functional_class_distribution.csv')
missing_data = pd.read_csv('missing_data_summary.csv')

print(f"  Main dataset: {len(df)} proteins")
print(f"  High-confidence: {len(df_high)} proteins")

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Calculate magnitude score for volcano plot
df['magnitude_score'] = -np.log10(np.abs(df['log_2_fold_change']) + 0.001)
df_high['magnitude_score'] = -np.log10(np.abs(df_high['log_2_fold_change']) + 0.001)

# Assign colors based on fold change
def assign_color(fc):
    if fc > FC_THRESHOLD:
        return COLOR_UP
    elif fc < (1/FC_THRESHOLD):
        return COLOR_DOWN
    else:
        return COLOR_NS

df['color'] = df['Fold_Change'].apply(assign_color)
df_high['color'] = df_high['Fold_Change'].apply(assign_color)

# Regulation status
def assign_regulation(fc):
    if fc > FC_THRESHOLD:
        return 'Upregulated'
    elif fc < (1/FC_THRESHOLD):
        return 'Downregulated'
    else:
        return 'Unchanged'

df['Regulation_Status'] = df['Fold_Change'].apply(assign_regulation)
df_high['Regulation_Status'] = df_high['Fold_Change'].apply(assign_regulation)

print("\nData preparation complete!")

# ============================================================================
# PLOT 1: VOLCANO PLOT (Enhanced version)
# ============================================================================

print("\n[1/12] Creating Volcano Plot...")

fig, ax = plt.subplots(figsize=FIGURE_SIZE_STANDARD, dpi=FIGURE_DPI)

# Scatter points
ax.scatter(df_high['log_2_fold_change'], df_high['magnitude_score'],
           c=df_high['color'], alpha=0.6, s=80, edgecolors='black', linewidth=0.5, zorder=2)

# Add threshold lines
ax.axvline(LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, zorder=1)
ax.axvline(-LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, zorder=1)
ax.axhline(0, color='black', linestyle='-', alpha=0.2, linewidth=0.8, zorder=1)
ax.axvline(0, color='black', linestyle='-', alpha=0.2, linewidth=0.8, zorder=1)

# Label top proteins
df_high_valid = df_high.dropna(subset=['Fold_Change'])
top_up = df_high_valid.nlargest(10, 'Fold_Change')
top_down = df_high_valid.nsmallest(10, 'Fold_Change')
to_label = pd.concat([top_up, top_down])

texts = []
for idx, row in to_label.iterrows():
    texts.append(ax.text(row['log_2_fold_change'], row['magnitude_score'],
                        row['Gene'], fontsize=8, fontweight='bold', zorder=3))

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.7))

# Labels and styling
ax.set_xlabel('Log₂ Fold Change (Testosterone vs Vehicle)', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('-Log₁₀(|Log₂ FC|) [Magnitude Score]', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('Volcano Plot: Testosterone vs Vehicle Control\n(High-Confidence Proteins)',
             fontsize=TITLE_SIZE, fontweight='bold', pad=20)
ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)

# Legend
up_count = len(df_high[df_high['Regulation_Status'] == 'Upregulated'])
down_count = len(df_high[df_high['Regulation_Status'] == 'Downregulated'])
ns_count = len(df_high[df_high['Regulation_Status'] == 'Unchanged'])

legend_elements = [
    Patch(facecolor=COLOR_UP, edgecolor='black', label=f'Upregulated (n={up_count})'),
    Patch(facecolor=COLOR_DOWN, edgecolor='black', label=f'Downregulated (n={down_count})'),
    Patch(facecolor=COLOR_NS, edgecolor='black', label=f'Unchanged (n={ns_count})')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=LEGEND_SIZE, framealpha=0.9)

plt.tight_layout()
plt.savefig(f"{folders['volcano']}/01_volcano_plot.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['volcano']}/01_volcano_plot.pdf", bbox_inches='tight')
plt.close()

# Interactive version
fig_int = go.Figure()

for status in ['Upregulated', 'Downregulated', 'Unchanged']:
    df_subset = df_high[df_high['Regulation_Status'] == status]
    color_map = {'Upregulated': COLOR_UP, 'Downregulated': COLOR_DOWN, 'Unchanged': COLOR_NS}

    fig_int.add_trace(go.Scatter(
        x=df_subset['log_2_fold_change'],
        y=df_subset['magnitude_score'],
        mode='markers',
        name=f'{status} (n={len(df_subset)})',
        marker=dict(size=10, color=color_map[status], line=dict(width=1, color='black')),
        text=[f"<b>{row['Gene']}</b><br>{row['Protein_Description']}<br>FC: {row['Fold_Change']:.2f}"
              for _, row in df_subset.iterrows()],
        hovertemplate='%{text}<extra></extra>'
    ))

fig_int.update_layout(
    title='Interactive Volcano Plot: Testosterone vs Vehicle',
    xaxis_title='Log₂ Fold Change',
    yaxis_title='-Log₁₀(Magnitude)',
    hovermode='closest',
    template='plotly_white',
    width=1000,
    height=700
)

fig_int.write_html(f"{folders['volcano']}/01_volcano_plot_interactive.html")

print("  ✓ Volcano plot saved (static + interactive)")

# ============================================================================
# PLOT 2: CLUSTERED HEATMAP
# ============================================================================

print("[2/12] Creating Clustered Heatmap...")

# Prepare heatmap data
heatmap_data = df_high[['Gene', 'Vehicle_Mean', 'Testosterone_Mean', 'Functional_Class']].copy()
heatmap_data = heatmap_data.set_index('Gene')

# Log transform
expression_matrix = heatmap_data[['Vehicle_Mean', 'Testosterone_Mean']].apply(lambda x: np.log2(x + 1))
expression_matrix = expression_matrix.fillna(0)

# Z-score normalization with careful handling of edge cases
expression_matrix_zscore = expression_matrix.apply(lambda x: (x - x.mean()) / (x.std() if x.std() > 0 else 1), axis=1)
expression_matrix_zscore = expression_matrix_zscore.fillna(0).replace([np.inf, -np.inf], 0)

# Functional class colors
func_classes = df_high['Functional_Class'].unique()
class_colors = dict(zip(func_classes, sns.color_palette("Set2", len(func_classes))))
row_colors = df_high.set_index('Gene')['Functional_Class'].map(class_colors)

# Create clustered heatmap
g = sns.clustermap(expression_matrix_zscore,
                   cmap='RdBu_r',
                   center=0,
                   row_colors=row_colors,
                   figsize=(10, 14),
                   dendrogram_ratio=0.15,
                   cbar_pos=(0.02, 0.8, 0.03, 0.15),
                   linewidths=0.5,
                   yticklabels=True,
                   xticklabels=['Vehicle', 'Testosterone'])

g.ax_heatmap.set_xlabel('Condition', fontsize=LABEL_SIZE, fontweight='bold')
g.ax_heatmap.set_ylabel('Protein', fontsize=LABEL_SIZE, fontweight='bold')
plt.suptitle('Hierarchically Clustered Heatmap\n(Z-score normalized, High-Confidence Proteins)',
             fontsize=TITLE_SIZE, fontweight='bold', y=0.98)

plt.savefig(f"{folders['heatmap']}/02_clustered_heatmap.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['heatmap']}/02_clustered_heatmap.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Clustered heatmap saved")

# ============================================================================
# PLOT 3: TOP CHANGERS BAR CHART
# ============================================================================

print("[3/12] Creating Top Changers Bar Chart...")

top_up = df_high.dropna(subset=['Fold_Change']).nlargest(10, 'Fold_Change')[['Gene', 'Fold_Change', 'Functional_Class']]
top_down = df_high.dropna(subset=['Fold_Change']).nsmallest(10, 'Fold_Change')[['Gene', 'Fold_Change', 'Functional_Class']]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8), dpi=FIGURE_DPI)

# Upregulated
ax1.barh(range(len(top_up)), top_up['Fold_Change'],
         color=COLOR_UP, alpha=0.7, edgecolor='black', linewidth=1)
ax1.set_yticks(range(len(top_up)))
ax1.set_yticklabels(top_up['Gene'], fontsize=12)
ax1.set_xlabel('Fold Change', fontsize=LABEL_SIZE, fontweight='bold')
ax1.set_title('Top 10 Upregulated Proteins', fontsize=TITLE_SIZE, fontweight='bold')
ax1.invert_yaxis()
ax1.grid(axis='x', alpha=0.3)

for i, (idx, row) in enumerate(top_up.iterrows()):
    ax1.text(row['Fold_Change'] + 0.2, i, f"{row['Fold_Change']:.2f}×",
             va='center', fontsize=10, fontweight='bold')

# Downregulated
ax2.barh(range(len(top_down)), top_down['Fold_Change'],
         color=COLOR_DOWN, alpha=0.7, edgecolor='black', linewidth=1)
ax2.set_yticks(range(len(top_down)))
ax2.set_yticklabels(top_down['Gene'], fontsize=12)
ax2.set_xlabel('Fold Change', fontsize=LABEL_SIZE, fontweight='bold')
ax2.set_title('Top 10 Downregulated Proteins', fontsize=TITLE_SIZE, fontweight='bold')
ax2.invert_yaxis()
ax2.grid(axis='x', alpha=0.3)

for i, (idx, row) in enumerate(top_down.iterrows()):
    ax2.text(row['Fold_Change'] + 0.01, i, f"{row['Fold_Change']:.2f}×",
             va='center', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.savefig(f"{folders['bars']}/03_top_changers_bars.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['bars']}/03_top_changers_bars.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Top changers bar chart saved")

# ============================================================================
# PLOT 4: FUNCTIONAL CLASS DISTRIBUTION
# ============================================================================

print("[4/12] Creating Functional Class Distribution...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), dpi=FIGURE_DPI)

# Left: Pie chart
colors_pie = sns.color_palette("Set2", len(func_dist))
wedges, texts, autotexts = ax1.pie(func_dist['Count'],
                                     labels=func_dist['Functional_Class'],
                                     autopct='%1.1f%%',
                                     colors=colors_pie,
                                     startangle=90,
                                     wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))

for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
    autotext.set_fontsize(11)

ax1.set_title('Protein Distribution by Functional Class',
              fontsize=TITLE_SIZE, fontweight='bold')

# Right: Up/Down stacked bars
x = np.arange(len(func_dist))
width = 0.35

ax2.bar(x - width/2, func_dist['Upregulated_Count'], width,
        label='Upregulated', color=COLOR_UP, alpha=0.7, edgecolor='black')
ax2.bar(x + width/2, func_dist['Downregulated_Count'], width,
        label='Downregulated', color=COLOR_DOWN, alpha=0.7, edgecolor='black')

ax2.set_xlabel('Functional Class', fontsize=LABEL_SIZE, fontweight='bold')
ax2.set_ylabel('Number of Proteins', fontsize=LABEL_SIZE, fontweight='bold')
ax2.set_title('Regulation Direction by Functional Class',
              fontsize=TITLE_SIZE, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(func_dist['Functional_Class'], rotation=45, ha='right', fontsize=11)
ax2.legend(fontsize=LEGEND_SIZE)
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{folders['distributions']}/04_functional_class_distribution.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['distributions']}/04_functional_class_distribution.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Functional class distribution saved")

# ============================================================================
# PLOT 5: FOLD CHANGE DISTRIBUTION
# ============================================================================

print("[5/12] Creating Fold Change Distribution...")

fc_nonzero = df_high[df_high['Fold_Change'] > 0]['Fold_Change']

fig, ax = plt.subplots(figsize=FIGURE_SIZE_STANDARD, dpi=FIGURE_DPI)

# Histogram
ax.hist(fc_nonzero, bins=20, alpha=0.6, color='steelblue',
        edgecolor='black', density=True, label='Histogram', linewidth=1.5)

# Density overlay
density = gaussian_kde(fc_nonzero)
xs = np.linspace(fc_nonzero.min(), fc_nonzero.max(), 200)
ax.plot(xs, density(xs), 'r-', linewidth=3, label='Density', alpha=0.8)

# Statistics lines
mean_fc = fc_nonzero.mean()
median_fc = fc_nonzero.median()
ax.axvline(mean_fc, color='green', linestyle='--', linewidth=2.5, label=f'Mean: {mean_fc:.2f}', alpha=0.8)
ax.axvline(median_fc, color='orange', linestyle='--', linewidth=2.5, label=f'Median: {median_fc:.2f}', alpha=0.8)

ax.set_xlabel('Fold Change', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('Density', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('Fold Change Distribution\n(Testosterone vs Vehicle, High-Confidence Proteins)',
             fontsize=TITLE_SIZE, fontweight='bold')
ax.legend(fontsize=LEGEND_SIZE, loc='upper right')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{folders['distributions']}/05_fold_change_distribution.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['distributions']}/05_fold_change_distribution.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Fold change distribution saved")

# ============================================================================
# PLOT 6: MA PLOT
# ============================================================================

print("[6/12] Creating MA Plot...")

df_high['A_value'] = (np.log2(df_high['Vehicle_Mean'] + 1) + np.log2(df_high['Testosterone_Mean'] + 1)) / 2
df_high['M_value'] = df_high['log_2_fold_change']

fig, ax = plt.subplots(figsize=FIGURE_SIZE_STANDARD, dpi=FIGURE_DPI)

ax.scatter(df_high['A_value'], df_high['M_value'], c=df_high['color'],
           alpha=0.6, s=60, edgecolors='black', linewidth=0.5)

# Reference lines
ax.axhline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)
ax.axhline(LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax.axhline(-LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)

ax.set_xlabel('A: Average Log₂ Intensity', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('M: Log₂ Fold Change', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('MA Plot: Intensity vs Fold Change\n(High-Confidence Proteins)',
             fontsize=TITLE_SIZE, fontweight='bold')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{folders['correlations']}/06_ma_plot.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['correlations']}/06_ma_plot.pdf", bbox_inches='tight')
plt.close()

print("  ✓ MA plot saved")

# ============================================================================
# PLOT 7: VEHICLE VS TESTOSTERONE SCATTER
# ============================================================================

print("[7/12] Creating Vehicle vs Testosterone Scatter...")

vehicle_log = np.log2(df_high['Vehicle_Mean'] + 1)
testo_log = np.log2(df_high['Testosterone_Mean'] + 1)

# Calculate correlation
valid_data = pd.DataFrame({'v': vehicle_log, 't': testo_log}).dropna()
r, p = pearsonr(valid_data['v'], valid_data['t'])

fig, ax = plt.subplots(figsize=FIGURE_SIZE_STANDARD, dpi=FIGURE_DPI)

ax.scatter(vehicle_log, testo_log, c=df_high['color'],
           alpha=0.6, s=70, edgecolors='black', linewidth=0.5)

# Diagonal line
lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0, linewidth=2, label='No Change (y=x)')

# Correlation text
ax.text(0.05, 0.95, f'Pearson R = {r:.3f}\np-value = {p:.2e}\nn = {len(valid_data)}',
        transform=ax.transAxes, fontsize=12, verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8, edgecolor='black', linewidth=1.5))

ax.set_xlabel('Vehicle: Log₂ Mean Intensity', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('Testosterone: Log₂ Mean Intensity', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('Condition Correlation: Vehicle vs Testosterone\n(High-Confidence Proteins)',
             fontsize=TITLE_SIZE, fontweight='bold')
ax.legend(fontsize=LEGEND_SIZE, loc='lower right')
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(f"{folders['correlations']}/07_scatter_vehicle_vs_testosterone.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['correlations']}/07_scatter_vehicle_vs_testosterone.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Scatter correlation plot saved")

# ============================================================================
# PLOT 8: BOX PLOTS BY FUNCTIONAL CLASS
# ============================================================================

print("[8/12] Creating Box Plots by Functional Class...")

fig, ax = plt.subplots(figsize=(14, 8), dpi=FIGURE_DPI)

# Prepare data
df_high_melted = df_high[['Gene', 'Functional_Class', 'Vehicle_Mean', 'Testosterone_Mean']].melt(
    id_vars=['Gene', 'Functional_Class'],
    value_vars=['Vehicle_Mean', 'Testosterone_Mean'],
    var_name='Condition',
    value_name='Intensity'
)
df_high_melted['Log2_Intensity'] = np.log2(df_high_melted['Intensity'] + 1)

# Violin + box + strip plot
sns.violinplot(data=df_high_melted, x='Functional_Class', y='Log2_Intensity',
               hue='Condition', split=False, inner=None, alpha=0.4, ax=ax)
sns.boxplot(data=df_high_melted, x='Functional_Class', y='Log2_Intensity',
            hue='Condition', width=0.3, ax=ax, showcaps=True, boxprops={'alpha': 0.7},
            whiskerprops={'linewidth': 1.5}, medianprops={'color': 'red', 'linewidth': 2})
sns.stripplot(data=df_high_melted, x='Functional_Class', y='Log2_Intensity',
              hue='Condition', dodge=True, alpha=0.5, size=3, ax=ax, legend=False)

# Remove duplicate legends
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], fontsize=LEGEND_SIZE, title='Condition', loc='upper right')

ax.set_xlabel('Functional Class', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('Log₂ Intensity', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('Protein Expression Distribution by Functional Class\n(Violin + Box + Strip Plot)',
             fontsize=TITLE_SIZE, fontweight='bold')
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(f"{folders['distributions']}/08_boxplots_by_functional_class.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['distributions']}/08_boxplots_by_functional_class.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Box plots by functional class saved")

# ============================================================================
# PLOT 9: REPLICATE CORRELATION GRID (1x2)
# ============================================================================

print("[9/12] Creating Replicate Correlation Grid...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=FIGURE_DPI)

# Actual replicate column names
vehicle_cols = ['Vehicle_Rep1_31579', 'Vehicle_Rep2_31581']
testo_cols = ['Testosterone_Rep1_31578', 'Testosterone_Rep2_31580']

# Plot 1: Vehicle rep1 vs rep2
ax = axes[0]
df_v12 = df_high[[vehicle_cols[0], vehicle_cols[1]]].dropna()
if len(df_v12) > 0:
    r12, p12 = pearsonr(df_v12[vehicle_cols[0]], df_v12[vehicle_cols[1]])
    ax.scatter(df_v12[vehicle_cols[0]], df_v12[vehicle_cols[1]],
               alpha=0.6, s=60, edgecolors='black', linewidth=0.5, color='steelblue')
    lims = [0, max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=2, label='Perfect Correlation')
    ax.text(0.05, 0.95, f'Pearson R = {r12:.3f}\np-value = {p12:.2e}\nn = {len(df_v12)}',
            transform=ax.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    ax.set_xlabel('Vehicle Rep1', fontsize=LABEL_SIZE, fontweight='bold')
    ax.set_ylabel('Vehicle Rep2', fontsize=LABEL_SIZE, fontweight='bold')
    ax.set_title('Vehicle Control: Rep1 vs Rep2', fontsize=TITLE_SIZE, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(alpha=0.3)

# Plot 2: Testosterone rep1 vs rep2
ax = axes[1]
df_t12 = df_high[[testo_cols[0], testo_cols[1]]].dropna()
if len(df_t12) > 0:
    rt12, pt12 = pearsonr(df_t12[testo_cols[0]], df_t12[testo_cols[1]])
    ax.scatter(df_t12[testo_cols[0]], df_t12[testo_cols[1]],
               alpha=0.6, s=60, edgecolors='black', linewidth=0.5, color='coral')
    lims = [0, max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=2, label='Perfect Correlation')
    ax.text(0.05, 0.95, f'Pearson R = {rt12:.3f}\np-value = {pt12:.2e}\nn = {len(df_t12)}',
            transform=ax.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))
    ax.set_xlabel('Testosterone Rep1', fontsize=LABEL_SIZE, fontweight='bold')
    ax.set_ylabel('Testosterone Rep2', fontsize=LABEL_SIZE, fontweight='bold')
    ax.set_title('Testosterone: Rep1 vs Rep2', fontsize=TITLE_SIZE, fontweight='bold')
    ax.legend(fontsize=10, loc='lower right')
    ax.grid(alpha=0.3)

plt.suptitle('Replicate Correlation Analysis\nHigh-Confidence Proteins (2 biological replicates per condition)',
             fontsize=TITLE_SIZE+2, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(f"{folders['correlations']}/09_replicate_correlation_grid.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['correlations']}/09_replicate_correlation_grid.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Replicate correlation grid saved")

# ============================================================================
# PLOT 10: MISSING DATA HEATMAP
# ============================================================================

print("[10/12] Creating Missing Data Heatmap...")

# Create binary matrix (1 = data present, 0 = missing)
replicate_cols = vehicle_cols + testo_cols
missing_matrix = df[['Gene'] + replicate_cols].set_index('Gene')
missing_binary = (~missing_matrix.isna()).astype(int)

fig, ax = plt.subplots(figsize=(8, 12), dpi=FIGURE_DPI)

# Create heatmap
sns.heatmap(missing_binary, cmap=['#d62728', '#2ca02c'], cbar_kws={'label': 'Data Availability'},
            linewidths=0.5, linecolor='gray', ax=ax, yticklabels=True,
            cbar=True)

# Customize colorbar
cbar = ax.collections[0].colorbar
cbar.set_ticks([0.25, 0.75])
cbar.set_ticklabels(['Missing', 'Present'])

ax.set_xlabel('Replicate', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('Protein', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('Missing Data Pattern Across Replicates\n(All 42 Proteins)',
             fontsize=TITLE_SIZE, fontweight='bold', pad=20)
ax.set_xticklabels(['V-Rep1', 'V-Rep2', 'T-Rep1', 'T-Rep2'],
                   rotation=45, ha='right', fontsize=11)

plt.tight_layout()
plt.savefig(f"{folders['qc']}/10_missing_data_heatmap.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['qc']}/10_missing_data_heatmap.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Missing data heatmap saved")

# ============================================================================
# PLOT 11: QC CONFIDENCE DISTRIBUTION
# ============================================================================

print("[11/12] Creating QC Confidence Distribution...")

confidence_counts = df['Confidence_Level'].value_counts()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7), dpi=FIGURE_DPI)

# Pie chart
colors_conf = ['#2ca02c', '#ff7f0e', '#d62728']  # Green, Orange, Red
wedges, texts, autotexts = ax1.pie(confidence_counts,
                                    labels=confidence_counts.index,
                                    autopct=lambda pct: f'{pct:.1f}%\n({int(pct/100.*confidence_counts.sum())})',
                                    colors=colors_conf,
                                    startangle=90,
                                    explode=(0.05, 0.05, 0.05),
                                    wedgeprops=dict(width=0.6, edgecolor='white', linewidth=3))

for text in texts:
    text.set_fontsize(13)
    text.set_fontweight('bold')

for autotext in autotexts:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
    autotext.set_fontsize(11)

ax1.set_title('Quality Control Confidence Distribution\n(All 42 Proteins)',
              fontsize=TITLE_SIZE, fontweight='bold', pad=20)

# Bar chart with confidence breakdown
ax2.bar(confidence_counts.index, confidence_counts.values,
        color=colors_conf, alpha=0.7, edgecolor='black', linewidth=2)

for i, (idx, val) in enumerate(confidence_counts.items()):
    ax2.text(i, val + 0.5, str(val), ha='center', fontsize=13, fontweight='bold')

ax2.set_xlabel('Confidence Level', fontsize=LABEL_SIZE, fontweight='bold')
ax2.set_ylabel('Number of Proteins', fontsize=LABEL_SIZE, fontweight='bold')
ax2.set_title('QC Confidence Breakdown', fontsize=TITLE_SIZE, fontweight='bold')
ax2.grid(axis='y', alpha=0.3)
ax2.set_ylim(0, max(confidence_counts.values) * 1.15)

plt.tight_layout()
plt.savefig(f"{folders['qc']}/11_qc_confidence_distribution.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['qc']}/11_qc_confidence_distribution.pdf", bbox_inches='tight')
plt.close()

print("  ✓ QC confidence distribution saved")

# ============================================================================
# PLOT 12: MULTI-PANEL PUBLICATION FIGURE (6 panels)
# ============================================================================

print("[12/12] Creating Multi-Panel Publication Figure...")

fig = plt.figure(figsize=(18, 24), dpi=FIGURE_DPI)
gs = fig.add_gridspec(6, 2, hspace=0.35, wspace=0.3)

# Panel A: Volcano plot
ax_a = fig.add_subplot(gs[0, :])
ax_a.scatter(df_high['log_2_fold_change'], df_high['magnitude_score'],
            c=df_high['color'], alpha=0.6, s=60, edgecolors='black', linewidth=0.5)
ax_a.axvline(LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax_a.axvline(-LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax_a.set_xlabel('Log₂ Fold Change', fontsize=12, fontweight='bold')
ax_a.set_ylabel('-Log₁₀(Magnitude)', fontsize=12, fontweight='bold')
ax_a.set_title('A. Volcano Plot', fontsize=14, fontweight='bold', loc='left')
ax_a.grid(alpha=0.3)

# Panel B: Top upregulated
ax_b = fig.add_subplot(gs[1, 0])
top_up_pub = df_high.nlargest(8, 'Fold_Change')[['Gene', 'Fold_Change']]
ax_b.barh(range(len(top_up_pub)), top_up_pub['Fold_Change'],
         color=COLOR_UP, alpha=0.7, edgecolor='black')
ax_b.set_yticks(range(len(top_up_pub)))
ax_b.set_yticklabels(top_up_pub['Gene'], fontsize=10)
ax_b.invert_yaxis()
ax_b.set_xlabel('Fold Change', fontsize=12, fontweight='bold')
ax_b.set_title('B. Top 8 Upregulated', fontsize=14, fontweight='bold', loc='left')
ax_b.grid(axis='x', alpha=0.3)

# Panel C: Top downregulated
ax_c = fig.add_subplot(gs[1, 1])
top_down_pub = df_high.nsmallest(8, 'Fold_Change')[['Gene', 'Fold_Change']]
ax_c.barh(range(len(top_down_pub)), top_down_pub['Fold_Change'],
         color=COLOR_DOWN, alpha=0.7, edgecolor='black')
ax_c.set_yticks(range(len(top_down_pub)))
ax_c.set_yticklabels(top_down_pub['Gene'], fontsize=10)
ax_c.invert_yaxis()
ax_c.set_xlabel('Fold Change', fontsize=12, fontweight='bold')
ax_c.set_title('C. Top 8 Downregulated', fontsize=14, fontweight='bold', loc='left')
ax_c.grid(axis='x', alpha=0.3)

# Panel D: Functional class distribution
ax_d = fig.add_subplot(gs[2, :])
x_func = np.arange(len(func_dist))
width_func = 0.35
ax_d.bar(x_func - width_func/2, func_dist['Upregulated_Count'], width_func,
        label='Upregulated', color=COLOR_UP, alpha=0.7, edgecolor='black')
ax_d.bar(x_func + width_func/2, func_dist['Downregulated_Count'], width_func,
        label='Downregulated', color=COLOR_DOWN, alpha=0.7, edgecolor='black')
ax_d.set_xlabel('Functional Class', fontsize=12, fontweight='bold')
ax_d.set_ylabel('Count', fontsize=12, fontweight='bold')
ax_d.set_title('D. Regulation by Functional Class', fontsize=14, fontweight='bold', loc='left')
ax_d.set_xticks(x_func)
ax_d.set_xticklabels(func_dist['Functional_Class'], rotation=45, ha='right', fontsize=10)
ax_d.legend(fontsize=11)
ax_d.grid(axis='y', alpha=0.3)

# Panel E: MA plot
ax_e = fig.add_subplot(gs[3, 0])
ax_e.scatter(df_high['A_value'], df_high['M_value'], c=df_high['color'],
            alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
ax_e.axhline(0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)
ax_e.axhline(LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5)
ax_e.axhline(-LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5)
ax_e.set_xlabel('Average Log₂ Intensity (A)', fontsize=12, fontweight='bold')
ax_e.set_ylabel('Log₂ Fold Change (M)', fontsize=12, fontweight='bold')
ax_e.set_title('E. MA Plot', fontsize=14, fontweight='bold', loc='left')
ax_e.grid(alpha=0.3)

# Panel F: Fold change distribution
ax_f = fig.add_subplot(gs[3, 1])
fc_pub = df_high[df_high['Fold_Change'] > 0]['Fold_Change']
ax_f.hist(fc_pub, bins=15, alpha=0.6, color='steelblue',
         edgecolor='black', density=True)
density_pub = gaussian_kde(fc_pub)
xs_pub = np.linspace(fc_pub.min(), fc_pub.max(), 200)
ax_f.plot(xs_pub, density_pub(xs_pub), 'r-', linewidth=2.5)
ax_f.axvline(fc_pub.mean(), color='green', linestyle='--', linewidth=2)
ax_f.set_xlabel('Fold Change', fontsize=12, fontweight='bold')
ax_f.set_ylabel('Density', fontsize=12, fontweight='bold')
ax_f.set_title('F. FC Distribution', fontsize=14, fontweight='bold', loc='left')
ax_f.grid(alpha=0.3)

# Panel G: QC Confidence (pie chart)
ax_g = fig.add_subplot(gs[4, 0])
colors_g = ['#2ca02c', '#ff7f0e', '#d62728']
wedges_g, texts_g, autotexts_g = ax_g.pie(confidence_counts,
                                            labels=confidence_counts.index,
                                            autopct='%1.1f%%',
                                            colors=colors_g,
                                            startangle=90,
                                            wedgeprops=dict(width=0.5, edgecolor='white', linewidth=2))
for autotext in autotexts_g:
    autotext.set_color('white')
    autotext.set_fontweight('bold')
ax_g.set_title('G. QC Confidence', fontsize=14, fontweight='bold', loc='left')

# Panel H: Vehicle vs Testosterone scatter
ax_h = fig.add_subplot(gs[4, 1])
ax_h.scatter(vehicle_log, testo_log, c=df_high['color'],
            alpha=0.6, s=50, edgecolors='black', linewidth=0.5)
lims_h = [np.min([ax_h.get_xlim(), ax_h.get_ylim()]), np.max([ax_h.get_xlim(), ax_h.get_ylim()])]
ax_h.plot(lims_h, lims_h, 'k--', alpha=0.5, linewidth=2)
ax_h.text(0.05, 0.95, f'R = {r:.3f}', transform=ax_h.transAxes, fontsize=11,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
ax_h.set_xlabel('Vehicle (Log₂)', fontsize=12, fontweight='bold')
ax_h.set_ylabel('Testosterone (Log₂)', fontsize=12, fontweight='bold')
ax_h.set_title('H. Condition Correlation', fontsize=14, fontweight='bold', loc='left')
ax_h.grid(alpha=0.3)

# Panel I: Regulation status summary (large text box)
ax_i = fig.add_subplot(gs[5, :])
ax_i.axis('off')

summary_text = f"""
PROTEOMICS ANALYSIS SUMMARY
Testosterone vs Vehicle Control (SILAC)

Dataset Overview:
• Total proteins analyzed: {len(df)}
• High-confidence proteins: {len(df_high)} ({len(df_high)/len(df)*100:.1f}%)
• Medium confidence: {len(df[df['Confidence_Level']=='Medium'])} ({len(df[df['Confidence_Level']=='Medium'])/len(df)*100:.1f}%)
• Low confidence: {len(df[df['Confidence_Level']=='Low'])} ({len(df[df['Confidence_Level']=='Low'])/len(df)*100:.1f}%)

Differential Expression (FC threshold = {FC_THRESHOLD}):
• Upregulated proteins: {up_count} ({up_count/len(df_high)*100:.1f}%)
• Downregulated proteins: {down_count} ({down_count/len(df_high)*100:.1f}%)
• Unchanged proteins: {ns_count} ({ns_count/len(df_high)*100:.1f}%)

Top Upregulated: {', '.join(df_high.nlargest(5, 'Fold_Change')['Gene'].tolist())}
Top Downregulated: {', '.join(df_high.nsmallest(5, 'Fold_Change')['Gene'].tolist())}

Functional Classes: {len(func_dist)} categories identified
Mean Fold Change: {df_high['Fold_Change'].mean():.2f}
Median Fold Change: {df_high['Fold_Change'].median():.2f}
"""

ax_i.text(0.5, 0.5, summary_text, transform=ax_i.transAxes,
         fontsize=11, verticalalignment='center', horizontalalignment='center',
         fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='#f0f0f0', alpha=0.9,
                  edgecolor='black', linewidth=2, pad=20))

plt.suptitle('Comprehensive Proteomics Analysis: Testosterone vs Vehicle Control\nSILAC Experiment - Publication Figure',
             fontsize=18, fontweight='bold', y=0.995)

plt.savefig(f"{folders['publication']}/12_multipanel_publication_figure.png", dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig(f"{folders['publication']}/12_multipanel_publication_figure.pdf", bbox_inches='tight')
plt.close()

print("  ✓ Multi-panel publication figure saved")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print("\n" + "="*80)
print("VISUALIZATION SUITE GENERATION COMPLETE!")
print("="*80)

print("\n📊 FILES CREATED:")
print("\n1. VOLCANO PLOT FOLDER:")
print(f"   • {folders['volcano']}/01_volcano_plot.png")
print(f"   • {folders['volcano']}/01_volcano_plot.pdf")
print(f"   • {folders['volcano']}/01_volcano_plot_interactive.html")

print("\n2. HEATMAPS FOLDER:")
print(f"   • {folders['heatmap']}/02_clustered_heatmap.png")
print(f"   • {folders['heatmap']}/02_clustered_heatmap.pdf")

print("\n3. BAR CHARTS FOLDER:")
print(f"   • {folders['bars']}/03_top_changers_bars.png")
print(f"   • {folders['bars']}/03_top_changers_bars.pdf")

print("\n4. DISTRIBUTIONS FOLDER:")
print(f"   • {folders['distributions']}/04_functional_class_distribution.png")
print(f"   • {folders['distributions']}/04_functional_class_distribution.pdf")
print(f"   • {folders['distributions']}/05_fold_change_distribution.png")
print(f"   • {folders['distributions']}/05_fold_change_distribution.pdf")
print(f"   • {folders['distributions']}/08_boxplots_by_functional_class.png")
print(f"   • {folders['distributions']}/08_boxplots_by_functional_class.pdf")

print("\n5. CORRELATIONS FOLDER:")
print(f"   • {folders['correlations']}/06_ma_plot.png")
print(f"   • {folders['correlations']}/06_ma_plot.pdf")
print(f"   • {folders['correlations']}/07_scatter_vehicle_vs_testosterone.png")
print(f"   • {folders['correlations']}/07_scatter_vehicle_vs_testosterone.pdf")
print(f"   • {folders['correlations']}/09_replicate_correlation_grid.png")
print(f"   • {folders['correlations']}/09_replicate_correlation_grid.pdf")

print("\n6. QC PLOTS FOLDER:")
print(f"   • {folders['qc']}/10_missing_data_heatmap.png")
print(f"   • {folders['qc']}/10_missing_data_heatmap.pdf")
print(f"   • {folders['qc']}/11_qc_confidence_distribution.png")
print(f"   • {folders['qc']}/11_qc_confidence_distribution.pdf")

print("\n7. PUBLICATION FIGURES FOLDER:")
print(f"   • {folders['publication']}/12_multipanel_publication_figure.png")
print(f"   • {folders['publication']}/12_multipanel_publication_figure.pdf")

print("\n" + "="*80)
print("📈 ANALYSIS STATISTICS:")
print("="*80)
print(f"Total proteins analyzed: {len(df)}")
print(f"High-confidence proteins: {len(df_high)}")
print(f"Upregulated (FC > {FC_THRESHOLD}): {up_count}")
print(f"Downregulated (FC < {1/FC_THRESHOLD:.2f}): {down_count}")
print(f"Unchanged: {ns_count}")
print(f"Functional classes: {len(func_dist)}")

print("\n" + "="*80)
print("✅ All 12 publication-quality figures generated successfully!")
print("="*80)
