"""
Figure 14 (a-d): QC figures.

Reads proteomics_de/results/{foldchange_all.csv, single_condition_proteins.csv}
(read-only) and writes 4 figures under proteomics_de/results/figures/:
  a. sample_correlation.png/.svg   - 4x4 Pearson correlation of log2 intensity
  b. intensity_distributions.png/.svg - per-sample log2 intensity distributions
  c. missing_values.png/.svg       - missingness bar chart + present/absent matrix
  d. rank_abundance.png/.svg       - sorted mean log2 intensity vs. rank

All four draw from the combined "whole proteome" universe: the 1948 rows in
foldchange_all.csv (both-condition-tested, including the 10 ON_OFF proteins)
unioned with the 606 rows in single_condition_proteins.csv (seen in only one
condition) = 2554 unique proteins total. This is a superset of the 1938-row
qc_limma.csv used in the volcano/MA/heatmap figures (those three are
restricted to proteins limma could actually test in both conditions).
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import style


def load_universe():
    fc = pd.read_csv(style.RESULTS_DIR / "foldchange_all.csv")
    sc = pd.read_csv(style.RESULTS_DIR / "single_condition_proteins.csv")

    fc = fc.rename(columns={"UniProt Accession Number": "accession", "Gene names": "gene"})
    fc_u = fc[["accession", "gene"] + style.SAMPLE_COLS].copy()
    fc_u["group"] = np.where(
        fc["regulated"] == "ON_OFF", "on_off",
        np.where(fc["complete"], "complete", "partial"),
    )

    sc_u = sc[["accession", "gene"] + style.SAMPLE_COLS].copy()
    sc_u["group"] = "single_condition"

    universe = pd.concat([fc_u, sc_u], ignore_index=True)
    return universe


GROUP_ORDER = ["complete", "partial", "on_off", "single_condition"]
GROUP_LABELS = {
    "complete": "Both-condition, complete\n(4/4 intensities)",
    "partial": "Both-condition, partial\n(>=1 of 4 missing)",
    "on_off": "ON/OFF\n(present in 1 condition only)",
    "single_condition": "Single-condition only\n(not in both-condition table)",
}


# ---------------------------------------------------------------------------
# a. Sample correlation matrix
# ---------------------------------------------------------------------------
def fig_sample_correlation(universe):
    log2_df = pd.DataFrame({
        style.SAMPLE_LABELS[c]: style.safe_log2(universe[c]) for c in style.SAMPLE_COLS
    })
    corr = log2_df.corr(method="pearson")  # pairwise, NaN-safe

    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    sns.heatmap(
        corr, annot=True, fmt=".3f", cmap=style.SEQUENTIAL_CMAP,
        vmin=corr.values[~np.eye(4, dtype=bool)].min() - 0.01, vmax=1.0,
        square=True, linewidths=1.5, linecolor=style.CHROME["surface"],
        cbar_kws={"label": "Pearson r (log2 intensity)"}, ax=ax,
        annot_kws={"fontsize": 9.5, "color": style.CHROME["ink_primary"]},
    )
    ax.set_title("Sample correlation — log2 intensity", pad=12)
    fig.suptitle(
        "Pairwise-complete Pearson r, 4 raw samples", fontsize=9,
        color=style.CHROME["ink_secondary"], y=0.965,
    )
    ax.tick_params(left=False, bottom=False)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    plt.setp(ax.get_yticklabels(), rotation=0)
    style.add_caveat(fig, y=-0.02)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    png, svg = style.save_fig(fig, "sample_correlation")

    key_numbers = {
        "n_proteins_used": int(log2_df.notna().all(axis=1).sum()),
        "n_proteins_universe": int(len(universe)),
        "control1_vs_control2_r": round(float(corr.loc["control_1", "control_2"]), 4),
        "treated1_vs_treated2_r": round(float(corr.loc["treated_1", "treated_2"]), 4),
        "control_vs_treated_r_range": [
            round(float(corr.loc[["control_1", "control_2"], ["treated_1", "treated_2"]].values.min()), 4),
            round(float(corr.loc[["control_1", "control_2"], ["treated_1", "treated_2"]].values.max()), 4),
        ],
        "full_matrix": corr.round(4).to_dict(),
    }
    style.record_manifest([{
        "file": png,
        "title": "Sample correlation matrix (log2 intensity)",
        "caption": (
            "Pairwise Pearson correlation of log2 intensity across the 4 raw "
            "samples, computed over the full 2554-protein identified universe "
            "(pairwise-complete, zero/absent treated as missing). Within-"
            "condition replicate pairs (control_1 vs control_2; treated_1 vs "
            "treated_2) are the QC check for technical reproducibility; "
            "cf. qc_replicate_correlation.csv, which reports a similar check "
            "on the derived fold-change ratios rather than raw intensities."
        ),
        "key_numbers": key_numbers,
    }])
    print("sample_correlation:", {k: v for k, v in key_numbers.items() if k != "full_matrix"})
    return png


# ---------------------------------------------------------------------------
# b. Intensity distributions
# ---------------------------------------------------------------------------
def fig_intensity_distributions(universe):
    long_rows = []
    for c in style.SAMPLE_COLS:
        vals = style.safe_log2(universe[c]).dropna()
        label = style.SAMPLE_LABELS[c]
        cond = style.SAMPLE_CONDITION[c]
        for v in vals:
            long_rows.append((label, cond, v))
    long_df = pd.DataFrame(long_rows, columns=["sample", "condition", "log2_intensity"])

    fig, ax = plt.subplots(figsize=(7.6, 5.6))
    order = [style.SAMPLE_LABELS[c] for c in style.SAMPLE_COLS]
    sns.violinplot(
        data=long_df, x="sample", y="log2_intensity", order=order,
        hue="condition", palette=style.CONDITION_COLORS, dodge=False,
        cut=0, inner=None, linewidth=0, alpha=0.35, ax=ax, legend=False,
    )
    sns.boxplot(
        data=long_df, x="sample", y="log2_intensity", order=order,
        width=0.18, showfliers=False, ax=ax,
        boxprops=dict(facecolor=style.CHROME["surface"], edgecolor=style.CHROME["ink_secondary"]),
        medianprops=dict(color=style.CHROME["ink_primary"], linewidth=1.6),
        whiskerprops=dict(color=style.CHROME["ink_secondary"]),
        capprops=dict(color=style.CHROME["ink_secondary"]),
    )
    medians = long_df.groupby("sample")["log2_intensity"].median().reindex(order)
    means = long_df.groupby("sample")["log2_intensity"].mean().reindex(order)

    style.recede_spines(ax)
    ax.set_xlabel("")
    ax.set_ylabel("log2 intensity (0 treated as missing)")
    ax.set_title("Per-sample intensity distributions")
    handles = [
        plt.Line2D([0], [0], color=style.CONDITION_COLORS["control"], lw=6, alpha=0.5),
        plt.Line2D([0], [0], color=style.CONDITION_COLORS["treated"], lw=6, alpha=0.5),
    ]
    ax.legend(handles, ["Control", "Treated"], loc="upper right", title="Condition")
    fig.suptitle(
        "Comparable medians/spread across samples support skipping an extra "
        "normalization step", fontsize=9.3, color=style.CHROME["ink_secondary"], y=0.965,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    png, svg = style.save_fig(fig, "intensity_distributions")

    key_numbers = {
        "n_values_per_sample": {s: int((long_df["sample"] == s).sum()) for s in order},
        "median_log2_intensity": {s: round(float(medians[s]), 3) for s in order},
        "mean_log2_intensity": {s: round(float(means[s]), 3) for s in order},
    }
    style.record_manifest([{
        "file": png,
        "title": "Per-sample log2 intensity distributions",
        "caption": (
            "Violin + box plots of log2 intensity for each of the 4 raw "
            "samples (zero/absent values excluded), colored by condition "
            "(control = blue, treated = red). Medians and spread are similar "
            "across all 4 samples, which is the QC basis for not applying an "
            "additional cross-sample normalization step upstream."
        ),
        "key_numbers": key_numbers,
    }])
    print("intensity_distributions:", key_numbers)
    return png


# ---------------------------------------------------------------------------
# c. Missing values
# ---------------------------------------------------------------------------
def fig_missing_values(universe):
    present = pd.DataFrame({
        style.SAMPLE_LABELS[c]: (universe[c].fillna(0) > 0) for c in style.SAMPLE_COLS
    })
    n_total = len(universe)
    missing_counts = (~present).sum()

    group_sizes = universe["group"].value_counts().reindex(GROUP_ORDER).fillna(0).astype(int)
    pct_present = pd.DataFrame({
        style.SAMPLE_LABELS[c]: universe.assign(_p=(universe[c].fillna(0) > 0)).groupby("group")["_p"].mean() * 100
        for c in style.SAMPLE_COLS
    }).reindex(GROUP_ORDER)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.2), gridspec_kw={"width_ratios": [1, 1.3]})

    # Left: per-sample missing/zero counts
    ax = axes[0]
    order = [style.SAMPLE_LABELS[c] for c in style.SAMPLE_COLS]
    colors = [style.CONDITION_COLORS[style.SAMPLE_CONDITION[c]] for c in style.SAMPLE_COLS]
    bars = ax.bar(order, [missing_counts[s] for s in order], color=colors, alpha=0.85, width=0.6)
    for b, s in zip(bars, order):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + n_total * 0.01,
                f"{missing_counts[s]:,}\n({missing_counts[s]/n_total:.0%})",
                ha="center", va="bottom", fontsize=8.5, color=style.CHROME["ink_secondary"])
    style.recede_spines(ax)
    ax.set_ylabel(f"Missing / zero-intensity proteins (of {n_total:,})")
    ax.set_title("Missing values per sample")
    ax.set_ylim(0, max(missing_counts) * 1.28)

    # Right: present/absent matrix by category
    ax2 = axes[1]
    mat = pct_present.values
    im = ax2.imshow(mat, cmap=style.SEQUENTIAL_CMAP, vmin=0, vmax=100, aspect="auto")
    ax2.set_xticks(range(len(order)))
    ax2.set_xticklabels(order, rotation=30, ha="right")
    ax2.set_yticks(range(len(GROUP_ORDER)))
    ax2.set_yticklabels([f"{GROUP_LABELS[g]}\n(n={group_sizes[g]:,})" for g in GROUP_ORDER], fontsize=8.5)
    for i in range(len(GROUP_ORDER)):
        for j in range(len(order)):
            v = mat[i, j]
            txt_color = "white" if v > 55 else style.CHROME["ink_primary"]
            ax2.text(j, i, f"{v:.0f}%", ha="center", va="center", fontsize=9.5, color=txt_color)
    ax2.set_title("% present, by protein category x sample")
    cbar = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    cbar.set_label("% present (non-zero)")
    for spine in ax2.spines.values():
        spine.set_visible(False)

    fig.suptitle(
        f"Whole identified proteome, n={n_total:,} (both-condition-tested + single-condition)",
        fontsize=9.3, color=style.CHROME["ink_secondary"], y=1.0,
    )
    style.add_caveat(fig, y=-0.02)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    png, svg = style.save_fig(fig, "missing_values")

    key_numbers = {
        "n_total_universe": int(n_total),
        "missing_or_zero_per_sample": {s: int(missing_counts[s]) for s in order},
        "group_sizes": {g: int(group_sizes[g]) for g in GROUP_ORDER},
        "pct_present_matrix": {g: {s: round(float(pct_present.loc[g, s]), 1) for s in order} for g in GROUP_ORDER},
    }
    style.record_manifest([{
        "file": png,
        "title": "Missingness summary",
        "caption": (
            "Left: count of missing/zero-intensity proteins per sample across "
            "the full 2554-protein identified universe. Right: % of proteins "
            "present (non-zero) per sample, broken out by protein category — "
            "complete both-condition (n={complete}), partial both-condition "
            "(n={partial}), ON/OFF (n={on_off}), and single-condition-only "
            "(n={single_condition}) — showing that most 'missingness' is "
            "structural (proteins genuinely absent from one condition) rather "
            "than random dropout."
        ).format(**{g: group_sizes[g] for g in GROUP_ORDER}),
        "key_numbers": key_numbers,
    }])
    print("missing_values:", {k: v for k, v in key_numbers.items() if k != "pct_present_matrix"})
    return png


# ---------------------------------------------------------------------------
# d. Rank abundance
# ---------------------------------------------------------------------------
def fig_rank_abundance(universe):
    log2_df = pd.DataFrame({c: style.safe_log2(universe[c]) for c in style.SAMPLE_COLS})
    mean_log2 = log2_df.mean(axis=1, skipna=True)
    valid = mean_log2.dropna().sort_values(ascending=False).reset_index(drop=True)
    rank = np.arange(1, len(valid) + 1)

    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    ax.plot(rank, valid.values, color=style.CATEGORICAL["blue"], lw=2)
    ax.fill_between(rank, valid.values, valid.min(), color=style.CATEGORICAL["blue"], alpha=0.08)

    top_idx, bot_idx = 0, len(valid) - 1
    top_gene = universe.loc[mean_log2.sort_values(ascending=False).index[0], "gene"]
    bot_gene = universe.loc[mean_log2.dropna().sort_values(ascending=True).index[0], "gene"]
    for idx, gene, xy_off in [
        (top_idx, top_gene, (12, -4)),
        (bot_idx, bot_gene, (-70, 8)),
    ]:
        ax.annotate(
            str(gene).split(";")[0], xy=(rank[idx], valid.values[idx]),
            xytext=xy_off, textcoords="offset points", fontsize=8.5,
            color=style.CHROME["ink_primary"],
            arrowprops=dict(arrowstyle="-", color=style.CHROME["ink_muted"], lw=0.7, alpha=0.7),
        )

    style.recede_spines(ax)
    ax.set_xlabel("Abundance rank")
    ax.set_ylabel("Mean log2 intensity (4 samples)")
    ax.set_title("Rank-abundance — dynamic range of the identified proteome")
    fig.suptitle(
        f"n = {len(valid):,} proteins with >=1 valid sample intensity; "
        f"dynamic range = {valid.max() - valid.min():.1f} log2 units",
        fontsize=9.3, color=style.CHROME["ink_secondary"], y=0.965,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    png, svg = style.save_fig(fig, "rank_abundance")

    key_numbers = {
        "n_points": int(len(valid)),
        "max_mean_log2_intensity": round(float(valid.max()), 3),
        "min_mean_log2_intensity": round(float(valid.min()), 3),
        "dynamic_range_log2_units": round(float(valid.max() - valid.min()), 3),
        "most_abundant_gene": str(top_gene),
        "least_abundant_gene": str(bot_gene),
    }
    style.record_manifest([{
        "file": png,
        "title": "Rank-abundance plot",
        "caption": (
            "Mean log2 intensity (across up to 4 samples, zero/absent excluded) "
            "sorted from most to least abundant, showing the proteome's dynamic "
            f"range (~{key_numbers['dynamic_range_log2_units']:.0f} log2 units, "
            "i.e. several orders of magnitude) — context for interpreting "
            "fold-change estimates, which are noisier near the low-abundance tail."
        ),
        "key_numbers": key_numbers,
    }])
    print("rank_abundance:", key_numbers)
    return png


def main():
    universe = load_universe()
    print(f"universe size: {len(universe)}")
    fig_sample_correlation(universe)
    fig_intensity_distributions(universe)
    fig_missing_values(universe)
    fig_rank_abundance(universe)


if __name__ == "__main__":
    main()
