"""
Figure 12: MA plot (mean log2 intensity vs. limma log2FC).

Reads proteomics_de/results/qc_limma.csv (read-only) and writes
proteomics_de/results/figures/ma_plot.png + .svg.

x = mean log2 intensity across the 4 raw samples (log2 of each intensity,
    treating values <= 0 as missing, then row-wise nanmean).
y = limma_log2FC.
Color follows the same diverging regulated scheme as the volcano plot, so the
two statistical figures share one visual language.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import style


def load_data():
    df = pd.read_csv(style.RESULTS_DIR / "qc_limma.csv")
    log2_cols = []
    for c in style.SAMPLE_COLS:
        lc = f"log2_{c}"
        df[lc] = style.safe_log2(df[c])
        log2_cols.append(lc)
    df["n_valid_samples"] = df[log2_cols].notna().sum(axis=1)
    df["mean_log2_intensity"] = df[log2_cols].mean(axis=1, skipna=True)
    return df, log2_cols


def plot(df):
    plotted = df[df["n_valid_samples"] > 0].copy()
    excluded = df[df["n_valid_samples"] == 0]

    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    counts = plotted["regulated"].value_counts()

    for cls in style.REGULATED_ORDER:
        if cls not in plotted["regulated"].unique():
            continue
        sub = plotted[plotted["regulated"] == cls]
        ax.scatter(
            sub["mean_log2_intensity"], sub["limma_log2FC"],
            s=14 if cls == "NO CHANGE" else 18,
            c=style.REGULATED_COLORS[cls],
            alpha=style.REGULATED_ALPHA[cls],
            linewidths=0,
            label=f"{cls.title() if cls != 'NO CHANGE' else 'No change'} (n={len(sub)})",
            zorder=3 if cls == "NO CHANGE" else 4,
        )

    ax.axhline(0, color=style.CHROME["baseline"], lw=1, ls="-", zorder=2)
    ax.axhline(style.FC_THRESHOLD, color=style.CHROME["ink_muted"], lw=1, ls="--", zorder=2)
    ax.axhline(-style.FC_THRESHOLD, color=style.CHROME["ink_muted"], lw=1, ls="--", zorder=2)

    style.recede_spines(ax)
    ax.set_xlabel("Mean log2 intensity (4 samples; 0 treated as missing)")
    ax.set_ylabel("limma log2 fold change (treated / control)")
    ax.set_title("MA plot — treated vs. control (SILAC)")
    ax.legend(loc="upper right", markerscale=1.6)
    fig.suptitle(
        "n = 2 technical replicates; color = regulated call at |log2FC| >= 0.585",
        fontsize=9.5, color=style.CHROME["ink_secondary"], y=0.955,
    )
    style.add_caveat(fig)
    fig.tight_layout(rect=(0, 0.035, 1, 0.94))
    return fig, counts, plotted, excluded


def main():
    df, log2_cols = load_data()
    fig, counts, plotted, excluded = plot(df)
    png, svg = style.save_fig(fig, "ma_plot")

    key_numbers = {
        "n_points_total": int(len(df)),
        "n_points_plotted": int(len(plotted)),
        "n_excluded_all_samples_missing": int(len(excluded)),
        "n_up": int(counts.get("UP", 0)),
        "n_down": int(counts.get("DOWN", 0)),
        "n_no_change": int(counts.get("NO CHANGE", 0)),
        "fc_threshold": style.FC_THRESHOLD,
        "mean_log2_intensity_range": [
            round(float(plotted["mean_log2_intensity"].min()), 3),
            round(float(plotted["mean_log2_intensity"].max()), 3),
        ],
    }
    style.record_manifest([{
        "file": png,
        "title": "MA plot — treated vs. control (SILAC)",
        "caption": (
            "Mean log2 intensity across the 4 raw samples (x) vs. limma log2 "
            "fold change (y) for the 1938 both-condition proteins; per-sample "
            "zero/absent intensities are treated as missing before averaging. "
            "Colored UP/DOWN/NO CHANGE at |log2FC| >= 0.585, dashed reference "
            "lines at 0 and +/-0.585. Same n=2 technical-replicate caveat as "
            "the volcano plot applies — no FDR-significant calls exist in this "
            "dataset."
        ),
        "key_numbers": key_numbers,
    }])
    print("ma_plot.py:", key_numbers)
    print("wrote:", png, svg)


if __name__ == "__main__":
    main()
