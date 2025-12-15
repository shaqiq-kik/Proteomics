"""
PLOT 1: VOLCANO PLOT
Professional visualization of fold change vs magnitude for proteomics data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from adjustText import adjust_text
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("PLOT 1: VOLCANO PLOT GENERATION")
print("="*80)

# ============================================================================
# CONFIGURATION
# ============================================================================

# Figure settings
FIGURE_SIZE_STANDARD = (10, 8)
FIGURE_DPI = 300
LABEL_SIZE = 12
TITLE_SIZE = 14
LEGEND_SIZE = 10

# Color scheme
COLOR_UP = '#E63946'      # Red for upregulated
COLOR_DOWN = '#457B9D'    # Blue for downregulated
COLOR_NS = '#CCCCCC'      # Gray for unchanged

# Thresholds
FC_THRESHOLD = 1.5        # Fold change threshold
LOG2FC_THRESHOLD = np.log2(FC_THRESHOLD)  # ±0.585

print(f"\nSettings:")
print(f"  Fold change threshold: {FC_THRESHOLD}")
print(f"  Log2 FC threshold: ±{LOG2FC_THRESHOLD:.3f}")
print(f"  Colors: Up={COLOR_UP}, Down={COLOR_DOWN}, NS={COLOR_NS}")

# ============================================================================
# LOAD DATA
# ============================================================================

# Use high confidence proteins for cleaner visualization
df = pd.read_csv('high_confidence_proteins.csv')
print(f"\nLoaded {len(df)} high-confidence proteins")

# Also prepare full dataset version
df_full = pd.read_csv('cleaned_proteomics_data_with_QC_flags.csv')
print(f"Loaded {len(df_full)} total proteins (for comparison)")

# ============================================================================
# PREPARE DATA FOR VOLCANO PLOT
# ============================================================================

# Calculate magnitude score (since no p-values available)
# Using -log10 of absolute log2 fold change as a proxy for significance
# Add small constant to avoid log(0)
df['magnitude_score'] = -np.log10(np.abs(df['log_2_fold_change']) + 0.001)
df_full['magnitude_score'] = -np.log10(np.abs(df_full['log_2_fold_change']) + 0.001)

# Assign colors based on fold change thresholds
def assign_color(row):
    if row['Fold_Change'] > FC_THRESHOLD:
        return COLOR_UP
    elif row['Fold_Change'] < (1/FC_THRESHOLD):
        return COLOR_DOWN
    else:
        return COLOR_NS

def assign_regulation(row):
    if row['Fold_Change'] > FC_THRESHOLD:
        return 'Upregulated'
    elif row['Fold_Change'] < (1/FC_THRESHOLD):
        return 'Downregulated'
    else:
        return 'Unchanged'

df['color'] = df.apply(assign_color, axis=1)
df['Regulation_Status'] = df.apply(assign_regulation, axis=1)

df_full['color'] = df_full.apply(assign_color, axis=1)
df_full['Regulation_Status'] = df_full.apply(assign_regulation, axis=1)

# Count regulation categories
up_count = len(df[df['Regulation_Status'] == 'Upregulated'])
down_count = len(df[df['Regulation_Status'] == 'Downregulated'])
ns_count = len(df[df['Regulation_Status'] == 'Unchanged'])

print(f"\nHigh-confidence protein distribution:")
print(f"  Upregulated (FC > {FC_THRESHOLD}): {up_count}")
print(f"  Downregulated (FC < {1/FC_THRESHOLD:.2f}): {down_count}")
print(f"  Unchanged: {ns_count}")

# ============================================================================
# STATIC VOLCANO PLOT (MATPLOTLIB)
# ============================================================================

print("\n[1] Creating static volcano plot...")

fig, ax = plt.subplots(figsize=FIGURE_SIZE_STANDARD, dpi=FIGURE_DPI)

# Scatter plot
scatter = ax.scatter(df['log_2_fold_change'],
                     df['magnitude_score'],
                     c=df['color'],
                     alpha=0.6,
                     s=80,
                     edgecolors='black',
                     linewidth=0.5,
                     zorder=2)

# Add threshold lines
ax.axvline(LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, zorder=1)
ax.axvline(-LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, zorder=1)

# Add horizontal line at y=0 for reference
ax.axhline(0, color='black', linestyle='-', alpha=0.2, linewidth=0.8, zorder=1)
ax.axvline(0, color='black', linestyle='-', alpha=0.2, linewidth=0.8, zorder=1)

# Label top proteins (top 10 up and top 10 down)
# Remove rows with NaN fold change
df_valid = df.dropna(subset=['Fold_Change'])

top_up = df_valid.nlargest(10, 'Fold_Change')
top_down = df_valid.nsmallest(10, 'Fold_Change')
to_label = pd.concat([top_up, top_down])

# Create text labels for adjustText
texts = []
for idx, row in to_label.iterrows():
    texts.append(ax.text(row['log_2_fold_change'],
                        row['magnitude_score'],
                        row['Gene'],
                        fontsize=8,
                        fontweight='bold',
                        zorder=3))

# Adjust text to avoid overlaps
adjust_text(texts,
           arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.7),
           expand_points=(1.5, 1.5),
           force_points=(0.5, 0.5))

# Labels and title
ax.set_xlabel('Log₂ Fold Change (Testosterone vs Vehicle)', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_ylabel('-Log₁₀(|Log₂ FC|) [Magnitude Score]', fontsize=LABEL_SIZE, fontweight='bold')
ax.set_title('Volcano Plot: Testosterone vs Vehicle Control\n(High-Confidence Proteins, n=25)',
             fontsize=TITLE_SIZE, fontweight='bold', pad=20)

# Add grid
ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLOR_UP, edgecolor='black', label=f'Upregulated (n={up_count})'),
    Patch(facecolor=COLOR_DOWN, edgecolor='black', label=f'Downregulated (n={down_count})'),
    Patch(facecolor=COLOR_NS, edgecolor='black', label=f'Unchanged (n={ns_count})')
]
ax.legend(handles=legend_elements,
         loc='upper left',
         fontsize=LEGEND_SIZE,
         framealpha=0.9,
         edgecolor='black')

# Add annotation about thresholds
threshold_text = f'Thresholds: FC > {FC_THRESHOLD} or < {1/FC_THRESHOLD:.2f}\n(Log₂ FC > ±{LOG2FC_THRESHOLD:.2f})'
ax.text(0.98, 0.02, threshold_text,
        transform=ax.transAxes,
        fontsize=8,
        verticalalignment='bottom',
        horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

plt.tight_layout()

# Save static versions
plt.savefig('../Visuals/Volcano Plot/01_volcano_plot.png', dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig('../Visuals/Volcano Plot/01_volcano_plot.pdf', bbox_inches='tight')
print("  ✓ Saved: ../Visuals/Volcano Plot/01_volcano_plot.png")
print("  ✓ Saved: ../Visuals/Volcano Plot/01_volcano_plot.pdf")

plt.show()
plt.close()

# ============================================================================
# INTERACTIVE VOLCANO PLOT (PLOTLY)
# ============================================================================

print("\n[2] Creating interactive volcano plot...")

# Prepare hover text
df['hover_text'] = (
    '<b>' + df['Gene'] + '</b><br>' +
    df['Protein_Description'] + '<br>' +
    'Fold Change: ' + df['Fold_Change'].round(2).astype(str) + '<br>' +
    'Log₂ FC: ' + df['log_2_fold_change'].round(3).astype(str) + '<br>' +
    'Confidence: ' + df['Confidence_Level']
)

# Create figure
fig = go.Figure()

# Add points for each regulation status
for status in ['Upregulated', 'Downregulated', 'Unchanged']:
    df_subset = df[df['Regulation_Status'] == status]

    color_map = {
        'Upregulated': COLOR_UP,
        'Downregulated': COLOR_DOWN,
        'Unchanged': COLOR_NS
    }

    fig.add_trace(go.Scatter(
        x=df_subset['log_2_fold_change'],
        y=df_subset['magnitude_score'],
        mode='markers',
        name=f'{status} (n={len(df_subset)})',
        marker=dict(
            size=10,
            color=color_map[status],
            line=dict(width=1, color='black'),
            opacity=0.7
        ),
        text=df_subset['hover_text'],
        hovertemplate='%{text}<extra></extra>',
        customdata=df_subset['Gene']
    ))

# Add threshold lines
fig.add_vline(x=LOG2FC_THRESHOLD, line_dash="dash", line_color="gray", opacity=0.5)
fig.add_vline(x=-LOG2FC_THRESHOLD, line_dash="dash", line_color="gray", opacity=0.5)
fig.add_hline(y=0, line_dash="dot", line_color="black", opacity=0.3)
fig.add_vline(x=0, line_dash="dot", line_color="black", opacity=0.3)

# Update layout
fig.update_layout(
    title={
        'text': 'Interactive Volcano Plot: Testosterone vs Vehicle Control<br><sub>High-Confidence Proteins (n=25)</sub>',
        'x': 0.5,
        'xanchor': 'center',
        'font': {'size': 16, 'family': 'Arial, sans-serif'}
    },
    xaxis_title='Log₂ Fold Change (Testosterone vs Vehicle)',
    yaxis_title='-Log₁₀(|Log₂ FC|) [Magnitude Score]',
    hovermode='closest',
    width=1000,
    height=700,
    template='plotly_white',
    font=dict(size=12),
    legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01,
        bgcolor='rgba(255,255,255,0.8)',
        bordercolor='gray',
        borderwidth=1
    ),
    xaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='black'
    ),
    yaxis=dict(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        zeroline=True,
        zerolinewidth=2,
        zerolinecolor='black'
    )
)

# Save interactive version
fig.write_html('../Visuals/Volcano Plot/01_volcano_plot_interactive.html')
print("  ✓ Saved: ../Visuals/Volcano Plot/01_volcano_plot_interactive.html")

# ============================================================================
# BONUS: ALL PROTEINS VERSION (with QC flags indicated)
# ============================================================================

print("\n[3] Creating full dataset version with QC indicators...")

fig_full, ax_full = plt.subplots(figsize=FIGURE_SIZE_STANDARD, dpi=FIGURE_DPI)

# Plot all proteins with different markers for confidence levels
for conf_level, marker, alpha in [('Low', 'x', 0.3), ('Medium', 's', 0.5), ('High', 'o', 0.7)]:
    df_conf = df_full[df_full['Confidence_Level'] == conf_level]
    if len(df_conf) > 0:
        ax_full.scatter(df_conf['log_2_fold_change'],
                       df_conf['magnitude_score'],
                       c=df_conf['color'],
                       marker=marker,
                       s=80 if marker == 'o' else 60,
                       alpha=alpha,
                       edgecolors='black',
                       linewidth=0.5,
                       label=f'{conf_level} Confidence (n={len(df_conf)})',
                       zorder=2 if conf_level == 'High' else 1)

# Add threshold lines
ax_full.axvline(LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax_full.axvline(-LOG2FC_THRESHOLD, color='gray', linestyle='--', alpha=0.5, linewidth=1.5)
ax_full.axhline(0, color='black', linestyle='-', alpha=0.2, linewidth=0.8)
ax_full.axvline(0, color='black', linestyle='-', alpha=0.2, linewidth=0.8)

# Labels and title
ax_full.set_xlabel('Log₂ Fold Change (Testosterone vs Vehicle)', fontsize=LABEL_SIZE, fontweight='bold')
ax_full.set_ylabel('-Log₁₀(|Log₂ FC|) [Magnitude Score]', fontsize=LABEL_SIZE, fontweight='bold')
ax_full.set_title('Volcano Plot: All Proteins with QC Indicators\n(n=42, marker shape indicates confidence)',
                  fontsize=TITLE_SIZE, fontweight='bold', pad=20)

ax_full.grid(True, alpha=0.2, linestyle=':', linewidth=0.5)
ax_full.legend(loc='upper left', fontsize=LEGEND_SIZE, framealpha=0.9, edgecolor='black')

# Add color legend manually
from matplotlib.lines import Line2D
color_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_UP, markersize=8, label='Upregulated', markeredgecolor='black'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_DOWN, markersize=8, label='Downregulated', markeredgecolor='black'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor=COLOR_NS, markersize=8, label='Unchanged', markeredgecolor='black')
]
ax_full.legend(handles=color_elements, loc='upper right', fontsize=LEGEND_SIZE,
               title='Regulation Status', framealpha=0.9, edgecolor='black')

plt.tight_layout()

plt.savefig('../Visuals/Volcano Plot/01_volcano_plot_all_proteins.png', dpi=FIGURE_DPI, bbox_inches='tight')
plt.savefig('../Visuals/Volcano Plot/01_volcano_plot_all_proteins.pdf', bbox_inches='tight')
print("  ✓ Saved: ../Visuals/Volcano Plot/01_volcano_plot_all_proteins.png")
print("  ✓ Saved: ../Visuals/Volcano Plot/01_volcano_plot_all_proteins.pdf")

plt.show()
plt.close()

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "="*80)
print("VOLCANO PLOT GENERATION COMPLETE!")
print("="*80)

print("\nFiles created in '../Visuals/Volcano Plot/' folder:")
print("  1. 01_volcano_plot.png - High-confidence proteins (static)")
print("  2. 01_volcano_plot.pdf - High-confidence proteins (static, vector)")
print("  3. 01_volcano_plot_interactive.html - High-confidence proteins (interactive)")
print("  4. 01_volcano_plot_all_proteins.png - All proteins with QC markers")
print("  5. 01_volcano_plot_all_proteins.pdf - All proteins with QC markers (vector)")

print("\nKey Statistics:")
print(f"  High-confidence proteins plotted: {len(df)}")
print(f"  Upregulated (FC > {FC_THRESHOLD}): {up_count} proteins")
print(f"  Downregulated (FC < {1/FC_THRESHOLD:.2f}): {down_count} proteins")
print(f"  Unchanged: {ns_count} proteins")
print(f"  Log₂ FC threshold: ±{LOG2FC_THRESHOLD:.3f}")

print("\nTop 5 Upregulated:")
for i, (idx, row) in enumerate(df.nlargest(5, 'Fold_Change').iterrows(), 1):
    print(f"  {i}. {row['Gene']:10s} - FC: {row['Fold_Change']:6.2f}, Log₂ FC: {row['log_2_fold_change']:6.2f}")

print("\nTop 5 Downregulated:")
for i, (idx, row) in enumerate(df.nsmallest(5, 'Fold_Change').iterrows(), 1):
    print(f"  {i}. {row['Gene']:10s} - FC: {row['Fold_Change']:6.2f}, Log₂ FC: {row['log_2_fold_change']:6.2f}")

print("\n" + "="*80)
