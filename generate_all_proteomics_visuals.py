"""
Comprehensive Proteomics Visualization Script
Generates all required visualizations for SILAC proteomics analysis
Author: Computational Proteomics Analyst
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set publication-ready style
plt.style.use('default')
sns.set_palette("husl")

# Color standards
COLOR_UPREGULATED = '#d62728'  # red
COLOR_DOWNREGULATED = '#1f77b4'  # blue
CONFIDENCE_COLORS = {
    'High': '#2ca02c',  # dark green
    'Medium': '#ff7f0e',  # orange
    'Low': '#7f7f7f'  # gray
}
FUNCTIONAL_CLASS_PALETTE = {
    'cytokine': '#e377c2',
    'growth factor': '#8c564b',
    'enzyme': '#17becf',
    'peptidase': '#bcbd22',
    'other': '#9467bd'
}

print("="*70)
print("PROTEOMICS VISUALIZATION PIPELINE")
print("="*70)

# ============================================================================
# 1. LOAD AND PREPARE DATA
# ============================================================================
print("\n[1/17] Loading dataset...")
data_path = "Analysis/General Analysis/cleaned_proteomics_data_with_QC_flags.csv"
df = pd.read_csv(data_path)

# Replace empty values with NaN
df = df.replace('', np.nan)

# Define replicate columns
replicate_cols = [
    'Vehicle_Rep1_31579',
    'Vehicle_Rep2_31581',
    'Testosterone_Rep1_31578',
    'Testosterone_Rep2_31580'
]

print(f"   Loaded {len(df)} proteins")
print(f"   Columns: {len(df.columns)}")

# ============================================================================
# 2. CREATE FOLDER STRUCTURE
# ============================================================================
print("\n[2/17] Creating folder structure...")
base_dir = Path("visuals")
folders = {
    'core_results': base_dir / 'core_results',
    'heatmaps': base_dir / 'heatmaps',
    'qc': base_dir / 'qc',
    'statistical_diagnostics': base_dir / 'statistical_diagnostics',
    'summaries': base_dir / 'summaries'
}

for folder_name, folder_path in folders.items():
    folder_path.mkdir(parents=True, exist_ok=True)
    print(f"   Created: {folder_path}")

# ============================================================================
# 3. VOLCANO-STYLE SCATTER PLOT
# ============================================================================
print("\n[3/17] Generating volcano-style scatter plot...")

fig, ax = plt.subplots(figsize=(12, 8))

# Prepare data
plot_df = df.dropna(subset=['log_2_fold_change']).copy()
plot_df['abs_log2FC'] = np.abs(plot_df['log_2_fold_change'])

# Plot by confidence level
for conf_level in ['High', 'Medium', 'Low']:
    subset = plot_df[plot_df['Confidence_Level'] == conf_level]

    # Different markers for FC_Type
    for fc_type in subset['FC_Type'].unique():
        data = subset[subset['FC_Type'] == fc_type]
        marker = 'o' if fc_type == 'Normal' else 's'
        edgecolor = 'black' if fc_type == 'Complete_Suppression' else 'none'

        ax.scatter(
            data['log_2_fold_change'],
            data['abs_log2FC'],
            c=CONFIDENCE_COLORS[conf_level],
            marker=marker,
            s=100,
            alpha=0.7,
            edgecolors=edgecolor,
            linewidths=1.5,
            label=f"{conf_level} ({fc_type})" if fc_type != 'Normal' else conf_level
        )

# Add vertical line at log2FC = 0
ax.axvline(x=0, color='black', linestyle='--', linewidth=1.5, alpha=0.5)

ax.set_xlabel('log₂ Fold Change (Testosterone/Vehicle)', fontsize=14, fontweight='bold')
ax.set_ylabel('|log₂ Fold Change|', fontsize=14, fontweight='bold')
ax.set_title('Differential Protein Regulation by Testosterone', fontsize=16, fontweight='bold', pad=20)
ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig(folders['core_results'] / 'volcano_log2FC_confidence.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['core_results'] / 'volcano_log2FC_confidence.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: volcano_log2FC_confidence.[png/pdf]")

# ============================================================================
# 4. RANKED WATERFALL PLOT
# ============================================================================
print("\n[4/17] Generating ranked waterfall plot...")

fig, ax = plt.subplots(figsize=(14, 7))

# Sort by log2FC and prepare data
waterfall_df = df.dropna(subset=['log_2_fold_change']).copy()
waterfall_df = waterfall_df.sort_values('log_2_fold_change', ascending=True).reset_index(drop=True)

# Create color array
colors = [COLOR_DOWNREGULATED if x < 0 else COLOR_UPREGULATED for x in waterfall_df['log_2_fold_change']]

# Create bar plot
bars = ax.bar(range(len(waterfall_df)), waterfall_df['log_2_fold_change'], color=colors, alpha=0.8, edgecolor='black', linewidth=0.5)

ax.axhline(y=0, color='black', linestyle='-', linewidth=2)
ax.set_xlabel('Proteins (ranked by log₂ FC)', fontsize=14, fontweight='bold')
ax.set_ylabel('log₂ Fold Change', fontsize=14, fontweight='bold')
ax.set_title('Ranked Protein Regulation Profile', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, axis='y', linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLOR_UPREGULATED, label='Upregulated (log₂FC > 0)'),
    Patch(facecolor=COLOR_DOWNREGULATED, label='Downregulated (log₂FC < 0)')
]
ax.legend(handles=legend_elements, loc='upper left', frameon=True, fancybox=True, shadow=True)

plt.tight_layout()
plt.savefig(folders['core_results'] / 'ranked_log2FC_waterfall.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['core_results'] / 'ranked_log2FC_waterfall.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: ranked_log2FC_waterfall.[png/pdf]")

# ============================================================================
# 5. TOP 20 CHANGERS BAR CHART
# ============================================================================
print("\n[5/17] Generating top 20 changers bar chart...")

fig, ax = plt.subplots(figsize=(12, 10))

# Get top changers
top_changers_df = df.dropna(subset=['log_2_fold_change']).copy()
top_up = top_changers_df.nlargest(10, 'log_2_fold_change')
top_down = top_changers_df.nsmallest(10, 'log_2_fold_change')
top_20 = pd.concat([top_down, top_up]).reset_index(drop=True)

# Create color mapping
colors_top20 = [CONFIDENCE_COLORS[conf] for conf in top_20['Confidence_Level']]

# Create horizontal bar plot
y_pos = np.arange(len(top_20))
bars = ax.barh(y_pos, top_20['log_2_fold_change'], color=colors_top20, alpha=0.8, edgecolor='black', linewidth=1)

ax.set_yticks(y_pos)
ax.set_yticklabels(top_20['Gene'], fontsize=10)
ax.axvline(x=0, color='black', linestyle='-', linewidth=2)
ax.set_xlabel('log₂ Fold Change', fontsize=14, fontweight='bold')
ax.set_ylabel('Gene', fontsize=14, fontweight='bold')
ax.set_title('Top 20 Differentially Regulated Proteins', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, axis='x', linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Add legend
legend_elements = [Patch(facecolor=CONFIDENCE_COLORS[conf], label=conf) for conf in ['High', 'Medium', 'Low']]
ax.legend(handles=legend_elements, loc='lower right', title='Confidence Level', frameon=True, fancybox=True, shadow=True)

plt.tight_layout()
plt.savefig(folders['core_results'] / 'top_20_changers.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['core_results'] / 'top_20_changers.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: top_20_changers.[png/pdf]")

# ============================================================================
# 6. MEAN EXPRESSION HEATMAP
# ============================================================================
print("\n[6/17] Generating mean expression heatmap...")

fig, ax = plt.subplots(figsize=(8, 14))

# Prepare data for heatmap
heatmap_df = df[['Gene', 'Vehicle_Mean', 'Testosterone_Mean', 'Functional_Class']].copy()
heatmap_df = heatmap_df.dropna(subset=['Vehicle_Mean', 'Testosterone_Mean'])

# Apply log10 transform
heatmap_data = heatmap_df[['Vehicle_Mean', 'Testosterone_Mean']].apply(lambda x: np.log10(x + 1))
heatmap_data.index = heatmap_df['Gene'].values

# Create heatmap with clustering
sns.clustermap(
    heatmap_data,
    cmap='RdYlBu_r',
    center=None,
    col_cluster=False,
    row_cluster=True,
    figsize=(8, 14),
    cbar_kws={'label': 'log₁₀(Intensity + 1)'},
    xticklabels=['Vehicle', 'Testosterone'],
    yticklabels=True,
    linewidths=0.5,
    linecolor='gray',
    dendrogram_ratio=0.15,
    cbar_pos=(0.02, 0.8, 0.03, 0.15)
)

plt.suptitle('Condition Mean Expression Profiles', fontsize=14, fontweight='bold', y=0.98)
plt.savefig(folders['heatmaps'] / 'heatmap_condition_means.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['heatmaps'] / 'heatmap_condition_means.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: heatmap_condition_means.[png/pdf]")

# ============================================================================
# 7. REPLICATE EXPRESSION HEATMAP
# ============================================================================
print("\n[7/17] Generating replicate expression heatmap...")

# Prepare replicate data
rep_heatmap_df = df[['Gene'] + replicate_cols].copy()
rep_heatmap_df = rep_heatmap_df.dropna(subset=replicate_cols, how='all')

# Apply log10 transform
rep_data = rep_heatmap_df[replicate_cols].apply(lambda x: np.log10(x + 1))
rep_data.index = rep_heatmap_df['Gene'].values

# Rename columns for clarity
rep_data.columns = ['Vehicle Rep1', 'Vehicle Rep2', 'Testosterone Rep1', 'Testosterone Rep2']

# Create clustered heatmap
sns.clustermap(
    rep_data,
    cmap='viridis',
    col_cluster=False,
    row_cluster=True,
    figsize=(10, 14),
    cbar_kws={'label': 'log₁₀(Intensity + 1)'},
    yticklabels=True,
    linewidths=0.5,
    linecolor='gray',
    dendrogram_ratio=0.15,
    cbar_pos=(0.02, 0.8, 0.03, 0.15)
)

plt.suptitle('Replicate Expression Profiles', fontsize=14, fontweight='bold', y=0.98)
plt.savefig(folders['heatmaps'] / 'heatmap_replicates.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['heatmaps'] / 'heatmap_replicates.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: heatmap_replicates.[png/pdf]")

# ============================================================================
# 8. CONFIDENCE LEVEL DISTRIBUTION
# ============================================================================
print("\n[8/17] Generating confidence level distribution...")

fig, ax = plt.subplots(figsize=(10, 6))

conf_counts = df['Confidence_Level'].value_counts()
colors_conf = [CONFIDENCE_COLORS[level] for level in conf_counts.index]

bars = ax.bar(conf_counts.index, conf_counts.values, color=colors_conf, alpha=0.8, edgecolor='black', linewidth=2)

# Add count labels on bars
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(height)}',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_xlabel('Confidence Level', fontsize=14, fontweight='bold')
ax.set_ylabel('Number of Proteins', fontsize=14, fontweight='bold')
ax.set_title('Protein Confidence Level Distribution', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, axis='y', linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig(folders['qc'] / 'confidence_level_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['qc'] / 'confidence_level_distribution.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: confidence_level_distribution.[png/pdf]")

# ============================================================================
# 9. REPLICATE COMPLETENESS DISTRIBUTION
# ============================================================================
print("\n[9/17] Generating replicate completeness distribution...")

fig, ax = plt.subplots(figsize=(10, 6))

rep_comp_counts = df['Replicate_Completeness'].value_counts()
colors_rep = ['#2ca02c', '#ff7f0e', '#d62728'][:len(rep_comp_counts)]

bars = ax.bar(rep_comp_counts.index, rep_comp_counts.values, color=colors_rep, alpha=0.8, edgecolor='black', linewidth=2)

# Add count labels on bars
for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height,
            f'{int(height)}',
            ha='center', va='bottom', fontsize=12, fontweight='bold')

ax.set_xlabel('Replicate Completeness', fontsize=14, fontweight='bold')
ax.set_ylabel('Number of Proteins', fontsize=14, fontweight='bold')
ax.set_title('Replicate Data Completeness', fontsize=16, fontweight='bold', pad=20)
ax.grid(True, alpha=0.3, axis='y', linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig(folders['qc'] / 'replicate_completeness_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['qc'] / 'replicate_completeness_distribution.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: replicate_completeness_distribution.[png/pdf]")

# ============================================================================
# 10. MISSING DATA MATRIX
# ============================================================================
print("\n[10/17] Generating missing data matrix...")

fig, ax = plt.subplots(figsize=(8, 14))

# Create binary matrix (1 = present, 0 = missing)
missing_df = df[['Gene'] + replicate_cols].copy()
missing_matrix = (~missing_df[replicate_cols].isna()).astype(int)
missing_matrix.index = missing_df['Gene'].values

# Rename columns
missing_matrix.columns = ['Vehicle\nRep1', 'Vehicle\nRep2', 'Testosterone\nRep1', 'Testosterone\nRep2']

sns.heatmap(
    missing_matrix,
    cmap=['#d62728', '#2ca02c'],  # red for missing, green for present
    cbar_kws={'label': 'Data Presence', 'ticks': [0, 1]},
    linewidths=1,
    linecolor='white',
    yticklabels=True,
    ax=ax,
    vmin=0,
    vmax=1
)

ax.set_xlabel('Replicates', fontsize=14, fontweight='bold')
ax.set_ylabel('Gene', fontsize=10)
ax.set_title('Missing Data Matrix', fontsize=16, fontweight='bold', pad=20)

# Customize colorbar
cbar = ax.collections[0].colorbar
cbar.set_ticks([0.25, 0.75])
cbar.set_ticklabels(['Missing (0)', 'Present (1)'])

plt.tight_layout()
plt.savefig(folders['qc'] / 'missing_data_matrix.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['qc'] / 'missing_data_matrix.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: missing_data_matrix.[png/pdf]")

# ============================================================================
# 11. COMPLETE SUPPRESSION PANEL
# ============================================================================
print("\n[11/17] Generating complete suppression panel...")

# Get proteins with complete suppression
suppressed = df[df['FC_Type'] == 'Complete_Suppression'].copy()

if len(suppressed) > 0:
    fig, ax = plt.subplots(figsize=(12, 6))

    x_pos = np.arange(len(suppressed))
    bars = ax.bar(x_pos, suppressed['Vehicle_Mean'], color='#1f77b4', alpha=0.8, edgecolor='black', linewidth=1.5, label='Vehicle Mean')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(suppressed['Gene'], rotation=45, ha='right', fontsize=10)
    ax.set_xlabel('Gene', fontsize=14, fontweight='bold')
    ax.set_ylabel('Vehicle Mean Intensity', fontsize=14, fontweight='bold')
    ax.set_title('Completely Suppressed Proteins (Testosterone = 0)', fontsize=16, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, axis='y', linestyle=':')
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')

    # Add annotation
    ax.text(0.98, 0.97, 'Testosterone: 0 for all proteins shown',
            transform=ax.transAxes, fontsize=12,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(folders['qc'] / 'complete_suppression_proteins.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(folders['qc'] / 'complete_suppression_proteins.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   ✓ Saved: complete_suppression_proteins.[png/pdf] ({len(suppressed)} proteins)")
else:
    print("   ⚠ No completely suppressed proteins found")

# ============================================================================
# 12. MA PLOT
# ============================================================================
print("\n[12/17] Generating MA plot...")

fig, ax = plt.subplots(figsize=(12, 8))

# Calculate A (average) and M (log2FC)
ma_df = df.dropna(subset=['Vehicle_Mean', 'Testosterone_Mean', 'log_2_fold_change']).copy()
ma_df['A'] = np.log10((ma_df['Vehicle_Mean'] + ma_df['Testosterone_Mean']) / 2 + 1)
ma_df['M'] = ma_df['log_2_fold_change']

# Plot by confidence level
for conf_level in ['High', 'Medium', 'Low']:
    subset = ma_df[ma_df['Confidence_Level'] == conf_level]
    ax.scatter(
        subset['A'],
        subset['M'],
        c=CONFIDENCE_COLORS[conf_level],
        s=100,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5,
        label=conf_level
    )

# Add horizontal line at M = 0
ax.axhline(y=0, color='black', linestyle='--', linewidth=2, alpha=0.7)

ax.set_xlabel('A: log₁₀(Mean Intensity)', fontsize=14, fontweight='bold')
ax.set_ylabel('M: log₂ Fold Change', fontsize=14, fontweight='bold')
ax.set_title('MA Plot: Intensity vs Fold Change', fontsize=16, fontweight='bold', pad=20)
ax.legend(title='Confidence Level', loc='upper right', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

plt.tight_layout()
plt.savefig(folders['statistical_diagnostics'] / 'MA_plot.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['statistical_diagnostics'] / 'MA_plot.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: MA_plot.[png/pdf]")

# ============================================================================
# 13. REPLICATE CORRELATION - VEHICLE
# ============================================================================
print("\n[13/17] Generating replicate correlation (Vehicle)...")

fig, ax = plt.subplots(figsize=(10, 10))

# Get data
veh_corr_df = df[['Vehicle_Rep1_31579', 'Vehicle_Rep2_31581']].dropna()

# Calculate correlation
from scipy.stats import pearsonr
if len(veh_corr_df) > 1:
    r, p = pearsonr(veh_corr_df['Vehicle_Rep1_31579'], veh_corr_df['Vehicle_Rep2_31581'])

    # Apply log10 for visualization
    x_log = np.log10(veh_corr_df['Vehicle_Rep1_31579'] + 1)
    y_log = np.log10(veh_corr_df['Vehicle_Rep2_31581'] + 1)

    ax.scatter(x_log, y_log, alpha=0.6, s=100, c='#1f77b4', edgecolors='black', linewidths=0.5)

    # Add diagonal line
    min_val = min(x_log.min(), y_log.min())
    max_val = max(x_log.max(), y_log.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, alpha=0.7)

    ax.set_xlabel('Vehicle Rep1 [log₁₀(Intensity + 1)]', fontsize=14, fontweight='bold')
    ax.set_ylabel('Vehicle Rep2 [log₁₀(Intensity + 1)]', fontsize=14, fontweight='bold')
    ax.set_title(f'Vehicle Replicate Correlation (Pearson r = {r:.3f}, p = {p:.2e})',
                 fontsize=16, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(folders['statistical_diagnostics'] / 'replicate_correlation_vehicle.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(folders['statistical_diagnostics'] / 'replicate_correlation_vehicle.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   ✓ Saved: replicate_correlation_vehicle.[png/pdf] (r={r:.3f})")

# ============================================================================
# 14. REPLICATE CORRELATION - TESTOSTERONE
# ============================================================================
print("\n[14/17] Generating replicate correlation (Testosterone)...")

fig, ax = plt.subplots(figsize=(10, 10))

# Get data
test_corr_df = df[['Testosterone_Rep1_31578', 'Testosterone_Rep2_31580']].dropna()

# Calculate correlation
if len(test_corr_df) > 1:
    r, p = pearsonr(test_corr_df['Testosterone_Rep1_31578'], test_corr_df['Testosterone_Rep2_31580'])

    # Apply log10 for visualization
    x_log = np.log10(test_corr_df['Testosterone_Rep1_31578'] + 1)
    y_log = np.log10(test_corr_df['Testosterone_Rep2_31580'] + 1)

    ax.scatter(x_log, y_log, alpha=0.6, s=100, c='#d62728', edgecolors='black', linewidths=0.5)

    # Add diagonal line
    min_val = min(x_log.min(), y_log.min())
    max_val = max(x_log.max(), y_log.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, alpha=0.7)

    ax.set_xlabel('Testosterone Rep1 [log₁₀(Intensity + 1)]', fontsize=14, fontweight='bold')
    ax.set_ylabel('Testosterone Rep2 [log₁₀(Intensity + 1)]', fontsize=14, fontweight='bold')
    ax.set_title(f'Testosterone Replicate Correlation (Pearson r = {r:.3f}, p = {p:.2e})',
                 fontsize=16, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    plt.savefig(folders['statistical_diagnostics'] / 'replicate_correlation_testosterone.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(folders['statistical_diagnostics'] / 'replicate_correlation_testosterone.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"   ✓ Saved: replicate_correlation_testosterone.[png/pdf] (r={r:.3f})")

# ============================================================================
# 15. SD PERCENT COMPARISON
# ============================================================================
print("\n[15/17] Generating SD percent comparison...")

fig, ax = plt.subplots(figsize=(10, 10))

# Get data
sd_df = df[['Vehicle_SD_Percent', 'Testosterone_SD_Percent', 'Confidence_Level']].dropna()

# Plot by confidence level
for conf_level in ['High', 'Medium', 'Low']:
    subset = sd_df[sd_df['Confidence_Level'] == conf_level]
    ax.scatter(
        subset['Vehicle_SD_Percent'],
        subset['Testosterone_SD_Percent'],
        c=CONFIDENCE_COLORS[conf_level],
        s=100,
        alpha=0.7,
        edgecolors='black',
        linewidths=0.5,
        label=conf_level
    )

# Add diagonal line
max_val = max(sd_df['Vehicle_SD_Percent'].max(), sd_df['Testosterone_SD_Percent'].max())
ax.plot([0, max_val], [0, max_val], 'k--', linewidth=2, alpha=0.5)

ax.set_xlabel('Vehicle SD (%)', fontsize=14, fontweight='bold')
ax.set_ylabel('Testosterone SD (%)', fontsize=14, fontweight='bold')
ax.set_title('Standard Deviation Comparison Between Conditions', fontsize=16, fontweight='bold', pad=20)
ax.legend(title='Confidence Level', loc='upper left', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')
ax.set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.savefig(folders['statistical_diagnostics'] / 'sd_percent_comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['statistical_diagnostics'] / 'sd_percent_comparison.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: sd_percent_comparison.[png/pdf]")

# ============================================================================
# 16. FUNCTIONAL CLASS STACKED BAR
# ============================================================================
print("\n[16/17] Generating functional class stacked bar...")

fig, ax = plt.subplots(figsize=(12, 7))

# Prepare data
fc_class_df = df[df['log_2_fold_change'].notna()].copy()
fc_class_df['Regulation'] = fc_class_df['log_2_fold_change'].apply(lambda x: 'Upregulated' if x > 0 else 'Downregulated')

# Count by functional class and regulation
fc_counts = fc_class_df.groupby(['Functional_Class', 'Regulation']).size().unstack(fill_value=0)

# Create stacked bar plot
fc_counts.plot(
    kind='bar',
    stacked=True,
    color=[COLOR_DOWNREGULATED, COLOR_UPREGULATED],
    alpha=0.8,
    edgecolor='black',
    linewidth=1.5,
    ax=ax
)

ax.set_xlabel('Functional Class', fontsize=14, fontweight='bold')
ax.set_ylabel('Number of Proteins', fontsize=14, fontweight='bold')
ax.set_title('Protein Regulation by Functional Class', fontsize=16, fontweight='bold', pad=20)
ax.legend(title='Regulation', frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, axis='y', linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.savefig(folders['summaries'] / 'functional_class_up_down.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['summaries'] / 'functional_class_up_down.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: functional_class_up_down.[png/pdf]")

# ============================================================================
# 17. FOLD-CHANGE DISTRIBUTION
# ============================================================================
print("\n[17/17] Generating fold-change distribution...")

fig, ax = plt.subplots(figsize=(12, 7))

# Get data
fc_dist_df = df['log_2_fold_change'].dropna()

# Separate up and down
up_fc = fc_dist_df[fc_dist_df > 0]
down_fc = fc_dist_df[fc_dist_df < 0]

# Create histogram
bins = np.linspace(fc_dist_df.min(), fc_dist_df.max(), 30)
ax.hist(down_fc, bins=bins, color=COLOR_DOWNREGULATED, alpha=0.7, edgecolor='black', linewidth=1, label='Downregulated')
ax.hist(up_fc, bins=bins, color=COLOR_UPREGULATED, alpha=0.7, edgecolor='black', linewidth=1, label='Upregulated')

# Add vertical line at 0
ax.axvline(x=0, color='black', linestyle='--', linewidth=2, alpha=0.7)

ax.set_xlabel('log₂ Fold Change', fontsize=14, fontweight='bold')
ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
ax.set_title('Distribution of log₂ Fold Changes', fontsize=16, fontweight='bold', pad=20)
ax.legend(frameon=True, fancybox=True, shadow=True)
ax.grid(True, alpha=0.3, axis='y', linestyle=':')
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Add statistics
stats_text = f'Total proteins: {len(fc_dist_df)}\n'
stats_text += f'Upregulated: {len(up_fc)} ({len(up_fc)/len(fc_dist_df)*100:.1f}%)\n'
stats_text += f'Downregulated: {len(down_fc)} ({len(down_fc)/len(fc_dist_df)*100:.1f}%)'

ax.text(0.98, 0.97, stats_text,
        transform=ax.transAxes, fontsize=11,
        verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

plt.tight_layout()
plt.savefig(folders['summaries'] / 'log2FC_distribution.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(folders['summaries'] / 'log2FC_distribution.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("   ✓ Saved: log2FC_distribution.[png/pdf]")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print("\n" + "="*70)
print("VISUALIZATION PIPELINE COMPLETE")
print("="*70)
print(f"\nAll visualizations saved to: {base_dir}/")
print("\nFolder structure:")
print(f"  {base_dir}/")
for folder_name, folder_path in folders.items():
    file_count = len(list(folder_path.glob('*.png')))
    print(f"  ├── {folder_name}/ ({file_count} plots)")

print("\n" + "="*70)
print("DONE!")
print("="*70)
