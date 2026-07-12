"""
Figure 13: Clustered heatmap of the top differentially-expressed proteins.

Reads proteomics_de/results/qc_limma.csv (read-only) and writes
proteomics_de/results/figures/heatmap_top_de.png + .svg.

Method (documented, as required by the task spec): the top 40 proteins are
selected by SMALLEST RAW p_value, not by largest |log2FC|. Rationale: limma's
moderated t-statistic already balances effect size against the precision
available from the two technical replicates, so ranking on p_value is a less
noise-prone selection than ranking on raw fold-change magnitude alone (a huge
log2FC estimated from very low/near-zero intensities is often the least
reliable, not the most interesting, protein - see e.g. the -10 log2FC outliers
visible near the origin-adjacent cloud in the MA plot).

Each protein's 4 raw sample intensities are log2-transformed (0/absent -> NaN)
and z-scored across the 4 samples (row-wise), then clustered by protein
(rows) only - columns stay in the fixed control_1/control_2/treated_1/treated_2
order so the two conditions remain directly comparable left-to-right.
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist

import style

N_TOP = 40


def short_gene(g):
    if not isinstance(g, str):
        return str(g)
    parts = g.split(";")
    return parts[0] + (f" (+{len(parts)-1})" if len(parts) > 1 else "")


def load_top(n=N_TOP):
    df = pd.read_csv(style.RESULTS_DIR / "qc_limma.csv")
    top = df.nsmallest(n, "p_value").copy()

    log2 = pd.DataFrame(index=top.index)
    for c in style.SAMPLE_COLS:
        log2[style.SAMPLE_LABELS[c]] = style.safe_log2(top[c])

    # Row-wise z-score across the 4 samples (ddof=0), NaN-safe.
    row_mean = log2.mean(axis=1, skipna=True)
    row_std = log2.std(axis=1, ddof=0, skipna=True).replace(0, np.nan)
    z = log2.sub(row_mean, axis=0).div(row_std, axis=0)

    n_missing_cells = int(z.isna().sum().sum())
    # Columns are already in control_1/control_2/treated_1/treated_2 order
    # because `style.SAMPLE_COLS` is iterated in that fixed order above.
    z.index = [short_gene(g) for g in top["gene"]]
    # De-duplicate any repeated short labels (protein groups sharing a lead gene).
    z.index = pd.Index(_dedupe(z.index))

    return top, z, n_missing_cells


def _dedupe(labels):
    seen = {}
    out = []
    for lab in labels:
        seen[lab] = seen.get(lab, 0) + 1
        out.append(lab if seen[lab] == 1 else f"{lab}.{seen[lab]}")
    return out


def plot(top, z):
    max_abs = np.nanmax(np.abs(z.values))
    row_colors = top["regulated"].map(style.REGULATED_COLORS).values
    row_colors = pd.Series(row_colors, index=z.index, name="regulated")

    # A handful of rows have exactly 1 of 4 intensities missing (0 in the raw
    # data), which leaves a NaN z-score cell that scipy's linkage() can't take
    # a distance on. Cluster on a 0-imputed copy (0 = the row mean after
    # z-scoring, i.e. "no information") but DISPLAY the true z matrix, so
    # missing cells render as blank rather than a fabricated color.
    row_link = linkage(pdist(z.fillna(0.0).values, metric="euclidean"), method="average")

    g = sns.clustermap(
        z,
        row_cluster=True,
        col_cluster=False,
        row_linkage=row_link,
        cmap=style.DIVERGING_CMAP,
        vmin=-max_abs,
        vmax=max_abs,
        center=0,
        row_colors=row_colors,
        figsize=(7.8, 11.5),
        dendrogram_ratio=(0.16, 0.02),
        cbar_pos=(0.02, 0.86, 0.03, 0.12),
        linewidths=0.4,
        linecolor=style.CHROME["surface"],
        xticklabels=True,
        yticklabels=True,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), rotation=35, ha="right")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), rotation=0, fontsize=7.6)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.cax.set_ylabel("row z-score\n(log2 intensity)", fontsize=8)
    g.cax.tick_params(labelsize=7.5)
    g.fig.suptitle(
        f"Top {len(z)} proteins by smallest raw p-value — row-standardized "
        "log2 intensity, clustered by protein (columns fixed)",
        fontsize=11.5, y=1.005, color=style.CHROME["ink_primary"],
    )
    g.fig.text(
        0.5, -0.01,
        style.CAVEAT_TEXT, ha="center", va="top", fontsize=7.8,
        style="italic", color=style.CHROME["ink_muted"], wrap=True,
    )

    # Legend for row_colors (UP/DOWN)
    handles = [
        plt.Line2D([0], [0], marker="s", linestyle="", markersize=9,
                   markerfacecolor=style.REGULATED_COLORS[c], markeredgewidth=0)
        for c in ["UP", "DOWN"] if c in top["regulated"].unique()
    ]
    labels = [c.title() for c in ["UP", "DOWN"] if c in top["regulated"].unique()]
    g.ax_heatmap.legend(
        handles, labels, title="Regulated call", loc="upper left",
        bbox_to_anchor=(1.06, 1.0), frameon=False, fontsize=8.5, title_fontsize=9,
    )
    return g, max_abs


def main():
    top, z, n_missing_cells = load_top()
    g, max_abs = plot(top, z)

    png_path = style.FIGURES_DIR / "heatmap_top_de.png"
    svg_path = style.FIGURES_DIR / "heatmap_top_de.svg"
    g.fig.savefig(png_path, dpi=200, bbox_inches="tight", facecolor=style.CHROME["surface"])
    g.fig.savefig(svg_path, bbox_inches="tight", facecolor=style.CHROME["surface"])
    plt.close(g.fig)
    png = str(png_path.relative_to(style.RESULTS_DIR))
    svg = str(svg_path.relative_to(style.RESULTS_DIR))

    counts = top["regulated"].value_counts()
    key_numbers = {
        "n_proteins": int(len(top)),
        "selection_method": "smallest raw p_value (limma), NOT largest |log2FC|",
        "n_up": int(counts.get("UP", 0)),
        "n_down": int(counts.get("DOWN", 0)),
        "n_missing_intensity_cells": n_missing_cells,
        "z_score_abs_max": round(float(max_abs), 3),
        "p_value_range": [
            round(float(top["p_value"].min()), 6),
            round(float(top["p_value"].max()), 6),
        ],
        "genes": [short_gene(x) for x in top["gene"].tolist()],
    }
    style.record_manifest([{
        "file": png,
        "title": "Clustered heatmap — top 40 DE proteins by raw p-value",
        "caption": (
            "Row-standardized (z-scored) log2 intensities for the 40 proteins "
            "with the smallest raw limma p-values (not the largest |log2FC| — "
            "see script docstring for rationale), clustered by protein only; "
            "sample columns stay in a fixed control/treated order. The row-color "
            "sidebar marks each protein's UP/DOWN call. As with all statistical "
            "figures here, n=2 technical replicates means none of these proteins "
            "are FDR-significant; this view is for pattern inspection only."
        ),
        "key_numbers": key_numbers,
    }])
    print("heatmap.py:", key_numbers)
    print("wrote:", png, svg)


if __name__ == "__main__":
    main()
