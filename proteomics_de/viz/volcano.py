"""
Figure 11: Volcano plot (limma log2FC vs. -log10 raw p-value).

Reads proteomics_de/results/qc_limma.csv (read-only) and writes
proteomics_de/results/figures/volcano.png + .svg.

Design notes (see style.py for the shared system):
  - Colored by `regulated` (UP/DOWN/NO CHANGE) as a diverging encoding:
    red = UP, blue = DOWN, recessive gray = NO CHANGE.
  - Dashed vertical lines at the |log2FC| = 0.585 up/down call boundary.
  - Dashed horizontal line at raw p = 0.05 (uncorrected). No FDR<0.05 line is
    drawn because nothing crosses it (min adj. p = 0.305 here) - that fact is
    stated explicitly in an annotation instead of a misleading absent line.
  - Top 10 proteins by smallest p_value are labeled with gene names using a
    small deterministic manual-offset scheme (no adjustText dependency).
  - Carries the n=2 technical-replicate caveat (style.add_caveat).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import style


def load_data():
    df = pd.read_csv(style.RESULTS_DIR / "qc_limma.csv")
    df["neg_log10_p"] = -np.log10(df["p_value"])
    return df


def label_offsets(n):
    """Deterministic (dx, dy) offsets in points, cycling to reduce overlap."""
    base = [
        (24, 10), (24, -12), (-24, 10), (-24, -12),
        (30, 22), (30, -24), (-30, 22), (-30, -24),
        (0, 26), (0, -26),
    ]
    return [base[i % len(base)] for i in range(n)]


def annotate_top_hits(ax, df, n=10):
    top = df.nsmallest(n, "p_value").copy()
    # Order by y (most significant first) so offsets fan out predictably.
    top = top.sort_values("neg_log10_p", ascending=False).reset_index(drop=True)
    offsets = label_offsets(len(top))
    for (dx, dy), (_, row) in zip(offsets, top.iterrows()):
        gene = row["gene"] if isinstance(row["gene"], str) else row["id"]
        gene_short = gene.split(";")[0] + ("…" if ";" in gene else "")
        ha = "left" if dx >= 0 else "right"
        ax.annotate(
            gene_short,
            xy=(row["limma_log2FC"], row["neg_log10_p"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=7.6,
            color=style.CHROME["ink_primary"],
            ha=ha,
            va="center",
            arrowprops=dict(
                arrowstyle="-", color=style.CHROME["ink_muted"],
                lw=0.6, alpha=0.7, shrinkA=2, shrinkB=2,
            ),
            zorder=6,
        )
    return top


def plot(df):
    fig, ax = plt.subplots(figsize=(9.6, 6.6))

    counts = df["regulated"].value_counts()
    for cls in style.REGULATED_ORDER:
        if cls not in df["regulated"].unique():
            continue
        sub = df[df["regulated"] == cls]
        n = len(sub)
        label = f"{cls.title() if cls != 'NO CHANGE' else 'No change'} (n={n})"
        ax.scatter(
            sub["limma_log2FC"], sub["neg_log10_p"],
            s=14 if cls == "NO CHANGE" else 18,
            c=style.REGULATED_COLORS[cls],
            alpha=style.REGULATED_ALPHA[cls],
            linewidths=0,
            label=label,
            zorder=3 if cls == "NO CHANGE" else 4,
        )

    # Reference lines
    ax.axvline(style.FC_THRESHOLD, color=style.CHROME["ink_muted"], lw=1, ls="--", zorder=2)
    ax.axvline(-style.FC_THRESHOLD, color=style.CHROME["ink_muted"], lw=1, ls="--", zorder=2)
    ax.axhline(-np.log10(style.RAW_P_THRESHOLD), color=style.CHROME["ink_muted"], lw=1, ls="--", zorder=2)
    ax.text(
        ax.get_xlim()[0], -np.log10(style.RAW_P_THRESHOLD),
        "raw p = 0.05 (uncorrected) ", fontsize=7.6, color=style.CHROME["ink_secondary"],
        va="bottom", ha="left",
    )

    min_adj_p = df["adj_p_value"].min()
    n_fdr_sig = int((df["adj_p_value"] < 0.05).sum())

    top = annotate_top_hits(ax, df, n=10)

    style.recede_spines(ax)
    ax.set_ylim(top=ax.get_ylim()[1] * 1.08)  # headroom for top-left labels
    ax.set_xlabel("limma log2 fold change (treated / control)")
    ax.set_ylabel("-log10(raw p-value)")
    ax.set_title("Volcano plot — treated vs. control (SILAC)", pad=12)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), markerscale=1.6,
              borderaxespad=0.0, title="Regulated call")
    fig.suptitle(
        f"n = 2 technical replicates — 0/{len(df)} proteins significant at "
        f"FDR < 0.05 (min adj. p = {min_adj_p:.3f}); no FDR line is drawn "
        f"because nothing crosses it",
        fontsize=9.3, color=style.CHROME["ink_secondary"], y=0.965,
    )
    style.add_caveat(fig)
    fig.tight_layout(rect=(0, 0.035, 0.83, 0.93))
    return fig, counts, top, min_adj_p, n_fdr_sig


def main():
    df = load_data()
    fig, counts, top, min_adj_p, n_fdr_sig = plot(df)
    png, svg = style.save_fig(fig, "volcano")

    key_numbers = {
        "n_points": int(len(df)),
        "n_up": int(counts.get("UP", 0)),
        "n_down": int(counts.get("DOWN", 0)),
        "n_no_change": int(counts.get("NO CHANGE", 0)),
        "n_raw_p_lt_0.05": int((df["p_value"] < 0.05).sum()),
        "n_significant_fdr_lt_0.05": n_fdr_sig,
        "min_adj_p_value": round(float(min_adj_p), 6),
        "fc_threshold": style.FC_THRESHOLD,
        "top_10_by_p_value_genes": top["gene"].fillna(top["id"]).tolist(),
        "top_10_by_p_value_p": [round(float(p), 6) for p in top["p_value"].tolist()],
    }
    style.record_manifest([{
        "file": png,
        "title": "Volcano plot — treated vs. control (SILAC)",
        "caption": (
            "limma log2 fold change vs. -log10(raw p-value) for the 1938 "
            "proteins tested in both conditions, colored UP/DOWN/NO CHANGE at "
            "|log2FC| >= 0.585. With n=2 technical replicates per condition "
            "(no biological replication), 0 of 1938 proteins pass FDR < 0.05 "
            "(minimum adjusted p = 0.305) — no dashed FDR line is drawn because "
            "nothing crosses it; the raw p = 0.05 line is shown instead for "
            "orientation only. Top 10 proteins by smallest raw p-value are "
            "labeled; treat all calls here as hypothesis-generating, not "
            "confirmed significant."
        ),
        "key_numbers": key_numbers,
    }])
    print("volcano.py:", key_numbers)
    print("wrote:", png, svg)


if __name__ == "__main__":
    main()
