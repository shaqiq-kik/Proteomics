"""
Section 6, item 19 -- gated PCA / clustering QC layer.

Reads proteomics_de/results/qc_limma.csv (read-only, 1938 both-condition
proteins) and writes NEW files only, under proteomics_de/results/gated/ and
proteomics_de/results/figures/ (a separate figures_manifest_gated.json; the
existing report manifests are never touched).

-----------------------------------------------------------------------------
SCOPE -- read this before trusting anything this script prints
-----------------------------------------------------------------------------
This is a SILAC design with n=2 technical replicates per condition, 4 samples
total. Per research1.md SS5 ("Limitations as Engineering Constraints"):

  - PCA on 4 samples: the sample covariance matrix is rank-deficient
    (rank <= n-1 = 3); only 3 PCs exist; eigenvalues are biased (the first
    ones overestimated, per the n<p literature); PC1 can be driven entirely
    by a single sample.
  - Clustering stability: hierarchical clustering on 4 samples has no
    bootstrap support possible -- dendrogram topology is not statistically
    stable.
  - The "trustworthy" gate for both is n_samples >= 6 (research1.md).

At n=4 we are BELOW that gate. This script still RUNS PCA and hierarchical
clustering -- they are useful as a QC sanity check (e.g. "do the two
technical replicates of each condition sit near each other?") -- but every
output is tagged QC_ONLY and captioned accordingly. Nothing here should be
read as evidence of biological structure or reproducible clusters.

-----------------------------------------------------------------------------
Forward-path dispatcher
-----------------------------------------------------------------------------
`REGISTRY` below declares, for every candidate analysis, the sample-count
gates from research1.md SS5 / SS "FORWARD-PATH SECTION". `dispatch()` compares
those gates against the CURRENT sample sheet (4 samples, 2 per group) and
assigns a status (QC_ONLY or SKIP) per analysis; the two QC_ONLY analyses
(pca, hierarchical_clustering) are then actually executed below, tagged
accordingly. When the sample sheet grows (config/sample_sheet.tsv in the
forward-path design), only N_SAMPLES / N_REPLICATES_PER_GROUP need to change
for the gates to flip automatically -- no other code changes required.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Import the shared report style module (proteomics_de/viz/style.py)
# ----------------------------------------------------------------------------
GATED_DIR = Path(__file__).resolve().parent
PROTEOMICS_DE_DIR = GATED_DIR.parent
VIZ_DIR = PROTEOMICS_DE_DIR / "viz"
sys.path.insert(0, str(VIZ_DIR))
import style  # noqa: E402  (palette, apply_style(), save_fig, CHROME, etc.)

RESULTS_DIR = style.RESULTS_DIR
FIGURES_DIR = style.FIGURES_DIR
GATED_RESULTS_DIR = RESULTS_DIR / "gated"
GATED_RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# NEW manifest -- never touches the existing figures_manifest*.json files
GATED_MANIFEST_PATH = FIGURES_DIR / "figures_manifest_gated.json"

SEED = style.SEED  # 42, set on import of style


def record_gated_manifest(entries):
    """Merge manifest entries into figures_manifest_gated.json (NEW file,
    keyed by `file`). Mirrors style.record_manifest / enrich_common's
    record_enrich_manifest but targets the gated-layer manifest path."""
    import json

    existing = []
    if GATED_MANIFEST_PATH.exists():
        try:
            existing = json.loads(GATED_MANIFEST_PATH.read_text())
        except json.JSONDecodeError:
            existing = []
    by_file = {e["file"]: e for e in existing}
    for e in entries:
        by_file[e["file"]] = e
    merged = list(by_file.values())
    GATED_MANIFEST_PATH.write_text(
        json.dumps(merged, indent=2, default=style._json_default)
    )
    return merged


# -----------------------------------------------------------------------
# Current sample sheet (hardcoded today; forward-path = read from
# config/sample_sheet.tsv once it exists -- see research1.md SS"FORWARD-PATH").
# -----------------------------------------------------------------------
N_SAMPLES = 4
N_GROUPS = 2
N_REPLICATES_PER_GROUP = N_SAMPLES // N_GROUPS  # 2

# Gate for "trustworthy" (not just runnable) PCA/clustering, per research1.md SS5.
TRUSTWORTHY_MIN_SAMPLES = 6
TRUSTWORTHY_MIN_REPLICATES_PER_GROUP = 3

# -----------------------------------------------------------------------
# Registry: every candidate analysis declares min_samples and
# min_replicates_per_group. This IS the forward-path mechanism: grow
# N_SAMPLES / N_REPLICATES_PER_GROUP above and status flips automatically,
# no code changes required elsewhere.
# -----------------------------------------------------------------------
REGISTRY = [
    dict(
        analysis="pca",
        min_samples_to_run=4,
        min_replicates_per_group_to_run=2,
        trustworthy_min_samples=TRUSTWORTHY_MIN_SAMPLES,
        trustworthy_min_replicates_per_group=TRUSTWORTHY_MIN_REPLICATES_PER_GROUP,
        reason_qc_only=(
            "n=4 < trustworthy gate n>=6 (research1.md SS5): the sample "
            "covariance matrix is rank-deficient (rank<=n-1=3, exactly 3 PCs "
            "exist by construction), eigenvalues are biased, and PC1 can be "
            "driven by a single sample. Runs for QC purposes only (e.g. do "
            "technical replicate pairs land near each other); NOT evidence "
            "of biological structure."
        ),
        reason_skip=None,
    ),
    dict(
        analysis="hierarchical_clustering",
        min_samples_to_run=4,
        min_replicates_per_group_to_run=2,
        trustworthy_min_samples=TRUSTWORTHY_MIN_SAMPLES,
        trustworthy_min_replicates_per_group=TRUSTWORTHY_MIN_REPLICATES_PER_GROUP,
        reason_qc_only=(
            "n=4 < trustworthy gate n>=6 (research1.md SS5): hierarchical "
            "clustering of 4 samples has no bootstrap support possible -- "
            "dendrogram topology is not statistically stable. Runs for QC "
            "purposes only; NOT a claim of reproducible sample clusters."
        ),
        reason_skip=None,
    ),
    dict(
        analysis="wgcna",
        min_samples_to_run=15,
        min_replicates_per_group_to_run=8,
        trustworthy_min_samples=20,
        trustworthy_min_replicates_per_group=10,
        reason_qc_only=None,
        reason_skip=(
            "WGCNA co-expression modules require >=15 samples (ideally "
            ">=20) for stable network topology; have 4."
        ),
    ),
    dict(
        analysis="umap_tsne",
        min_samples_to_run=24,
        min_replicates_per_group_to_run=12,
        trustworthy_min_samples=24,
        trustworthy_min_replicates_per_group=12,
        reason_qc_only=None,
        reason_skip=(
            "UMAP/t-SNE require 'dozens+' samples for a stable "
            "low-dimensional embedding (approximated here as >=24); have 4."
        ),
    ),
    dict(
        analysis="mixomics_splsda",
        min_samples_to_run=10 * N_GROUPS,
        min_replicates_per_group_to_run=10,
        trustworthy_min_samples=10 * N_GROUPS,
        trustworthy_min_replicates_per_group=10,
        reason_qc_only=None,
        reason_skip=(
            "mixOmics sPLS-DA requires samples >> classes (approximated "
            f"here as >=10x the {N_GROUPS} classes = {10 * N_GROUPS}); have 4."
        ),
    ),
    dict(
        analysis="vae_gnn",
        min_samples_to_run=1000,
        min_replicates_per_group_to_run=500,
        trustworthy_min_samples=1000,
        trustworthy_min_replicates_per_group=500,
        reason_qc_only=None,
        reason_skip=(
            "VAE/GNN architectures require thousands of samples to train "
            "without overfitting (approximated here as >=1000); have 4."
        ),
    ),
]


def dispatch():
    """Compare REGISTRY gates against the current sample sheet (N_SAMPLES,
    N_REPLICATES_PER_GROUP) and assign each analysis a status. Returns a
    list of row dicts with the skip_log.csv schema: analysis, min_samples,
    have_samples, status, reason.

    `min_samples` in the emitted row is the gate that is actually binding
    for that row's status: the trustworthy-gate for QC_ONLY rows (the
    number that, once met, would upgrade the analysis out of QC-only), and
    the run-gate for SKIP rows (the number that must be met before the
    analysis executes at all).
    """
    rows = []
    for entry in REGISTRY:
        can_run = (
            N_SAMPLES >= entry["min_samples_to_run"]
            and N_REPLICATES_PER_GROUP >= entry["min_replicates_per_group_to_run"]
        )
        if not can_run:
            status = "SKIP"
            min_samples = entry["min_samples_to_run"]
            reason = entry["reason_skip"]
        else:
            trustworthy = (
                N_SAMPLES >= entry["trustworthy_min_samples"]
                and N_REPLICATES_PER_GROUP >= entry["trustworthy_min_replicates_per_group"]
            )
            status = "TRUSTWORTHY" if trustworthy else "QC_ONLY"
            min_samples = entry["trustworthy_min_samples"]
            reason = entry["reason_qc_only"] if not trustworthy else (
                "n_samples/replicates meet the trustworthy gate."
            )
        rows.append(dict(
            analysis=entry["analysis"],
            min_samples=min_samples,
            have_samples=N_SAMPLES,
            status=status,
            reason=reason,
        ))
    return rows


def write_skip_log(rows):
    df = pd.DataFrame(rows, columns=["analysis", "min_samples", "have_samples", "status", "reason"])
    out_path = GATED_RESULTS_DIR / "skip_log.csv"
    df.to_csv(out_path, index=False)
    print(f"\n[dispatch] wrote {out_path} ({len(df)} rows)")
    print(df.to_string(index=False))
    return df, out_path


# -----------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------
def load_complete_matrix():
    """Build the log2 sample x protein matrix from qc_limma.csv.

    Zero/absent intensities are treated as missing (style.safe_log2). Only
    proteins with all 4 samples present ("complete") are used for PCA /
    clustering, per the task spec -- drop rows with any missing among the
    4 and report how many were used.
    """
    df = pd.read_csv(RESULTS_DIR / "qc_limma.csv")
    n_total = len(df)

    log2 = pd.DataFrame(
        {style.SAMPLE_LABELS[c]: style.safe_log2(df[c]) for c in style.SAMPLE_COLS}
    )
    log2.index = df["id"]

    complete_mask = log2.notna().all(axis=1)
    n_complete = int(complete_mask.sum())
    log2_complete = log2[complete_mask]

    print(f"[load] qc_limma.csv: {n_total} both-condition proteins")
    print(f"[load] complete (4/4 intensities present): {n_complete} / {n_total} used for PCA/clustering")

    return df, log2_complete, n_total, n_complete


SAMPLE_ORDER = [style.SAMPLE_LABELS[c] for c in style.SAMPLE_COLS]  # control_1, control_2, treated_1, treated_2
LABEL_CONDITION = {style.SAMPLE_LABELS[c]: style.SAMPLE_CONDITION[c] for c in style.SAMPLE_COLS}


# -----------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------
def run_pca(log2_complete, n_complete, n_total):
    """PCA on the 4 samples.

    Design choice (documented): samples are the 4 observations, proteins
    are the features. Each protein is standardized (z-scored across the 4
    samples, zero mean / unit variance) before PCA -- this is the
    correlation-based variant of PCA rather than raw-covariance PCA, and is
    the standard choice for expression data where proteins span several
    orders of magnitude in absolute intensity; without it, PCA would be
    dominated by the handful of highest-abundance proteins rather than by
    which proteins actually vary in a coordinated way across samples. This
    choice does not fix the fundamental n=4 degeneracy from research1.md
    SS5 (rank<=3, biased eigenvalues) -- it only controls for scale.
    """
    X = log2_complete.T  # 4 samples (rows) x n_complete proteins (columns)
    X = X.loc[SAMPLE_ORDER]  # fixed sample order

    scaler = StandardScaler()  # per-protein (per-column) z-score across the 4 samples
    X_scaled = scaler.fit_transform(X.values)

    n_components = min(3, X_scaled.shape[0] - 1)  # rank <= n_samples - 1 = 3
    pca = PCA(n_components=n_components, random_state=SEED)
    coords = pca.fit_transform(X_scaled)
    variance_ratio = pca.explained_variance_ratio_

    pc_cols = [f"PC{i+1}" for i in range(n_components)]
    coords_df = pd.DataFrame(coords, columns=pc_cols, index=SAMPLE_ORDER)
    coords_df.insert(0, "sample", SAMPLE_ORDER)
    coords_df = coords_df.reset_index(drop=True)

    variance_df = pd.DataFrame({
        "PC": pc_cols,
        "variance_explained": variance_ratio,
    })

    coords_path = GATED_RESULTS_DIR / "pca_coords.csv"
    variance_path = GATED_RESULTS_DIR / "pca_variance.csv"
    coords_df.to_csv(coords_path, index=False)
    variance_df.to_csv(variance_path, index=False)

    print(f"\n[pca] n_proteins_used={n_complete} (of {n_total} both-condition proteins)")
    print(f"[pca] wrote {coords_path}")
    print(coords_df.to_string(index=False))
    print(f"[pca] wrote {variance_path}")
    print(variance_df.to_string(index=False))
    print(f"[pca] cumulative variance explained by {n_components} PCs: {variance_ratio.sum():.4f} "
          f"(rank<=n-1={X_scaled.shape[0]-1} -> {'100%' if abs(variance_ratio.sum()-1.0)<1e-6 else f'{variance_ratio.sum():.1%}'} "
          f"by construction, direct evidence of the n=4 degeneracy)")

    return coords_df, variance_df, pc_cols


def fig_pca(coords_df, variance_df, pc_cols, n_complete, n_total):
    pct = {row["PC"]: row["variance_explained"] * 100 for _, row in variance_df.iterrows()}

    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    for _, row in coords_df.iterrows():
        cond = LABEL_CONDITION[row["sample"]]
        ax.scatter(
            row["PC1"], row["PC2"], s=200,
            color=style.CONDITION_COLORS[cond],
            edgecolor=style.CHROME["ink_primary"], linewidth=1.3,
            zorder=3, alpha=0.9,
        )
        ax.annotate(
            row["sample"], (row["PC1"], row["PC2"]),
            xytext=(9, 7), textcoords="offset points",
            fontsize=9.5, color=style.CHROME["ink_primary"],
        )

    style.recede_spines(ax)
    ax.axhline(0, color=style.CHROME["gridline"], lw=0.8, zorder=1)
    ax.axvline(0, color=style.CHROME["gridline"], lw=0.8, zorder=1)
    ax.set_xlabel(f"PC1 ({pct['PC1']:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pct['PC2']:.1f}% variance)")
    ax.set_title("PCA of 4 samples")

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=10,
                   markerfacecolor=style.CONDITION_COLORS["control"],
                   markeredgecolor=style.CHROME["ink_primary"]),
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=10,
                   markerfacecolor=style.CONDITION_COLORS["treated"],
                   markeredgecolor=style.CHROME["ink_primary"]),
    ]
    ax.legend(handles, ["Control", "Treated"], loc="best", title="Condition")

    fig.text(
        0.5, 0.975,
        f"n_proteins_used = {n_complete:,} of {n_total:,} (complete across all 4 samples); "
        f"trustworthy gate = n≥6 samples (research1.md §5)",
        ha="center", va="top", fontsize=8.8, color=style.CHROME["ink_secondary"],
    )
    fig.text(
        0.5, 0.925,
        "QC-ONLY — n=4 is statistically degenerate for PCA: rank≤n−1=3, "
        "only 3 PCs exist, eigenvalues biased, PC1 can be driven by 1 sample. "
        "NOT evidence of biological structure.",
        ha="center", va="top", fontsize=8.3, color="white", wrap=True,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=style.STATUS["critical"], edgecolor="none"),
    )
    style.add_caveat(fig, y=-0.02)
    fig.tight_layout(rect=(0, 0.02, 1, 0.83))
    png, svg = style.save_fig(fig, "pca_qc")
    print(f"[pca] wrote figures/pca_qc.png (+ .svg)")
    return png, svg, pct


# -----------------------------------------------------------------------
# Hierarchical clustering
# -----------------------------------------------------------------------
def run_clustering(log2_complete, n_complete):
    """Hierarchical clustering of the 4 samples.

    Distance: 1 - Pearson correlation between samples' log2 intensity
    vectors (scipy pdist metric='correlation'), computed on the same
    complete-protein matrix used for PCA.
    Linkage: average (UPGMA) -- the standard choice for correlation
    distances; Ward's method assumes (squared) Euclidean distance and is
    not appropriate here, and single linkage is prone to chaining on only
    4 observations.
    """
    X = log2_complete.T.loc[SAMPLE_ORDER]  # 4 samples x n_complete proteins
    dist_condensed = pdist(X.values, metric="correlation")  # 1 - Pearson r
    dist_square = squareform(dist_condensed)
    dist_df = pd.DataFrame(dist_square, index=SAMPLE_ORDER, columns=SAMPLE_ORDER)

    Z = linkage(dist_condensed, method="average")

    print(f"\n[cluster] n_proteins_used={n_complete}")
    print("[cluster] pairwise 1-Pearson-correlation distance matrix:")
    print(dist_df.round(4).to_string())

    return Z, dist_df


def fig_dendrogram(Z, dist_df, n_complete):
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ddata = dendrogram(
        Z, labels=SAMPLE_ORDER, ax=ax,
        color_threshold=0, above_threshold_color=style.CHROME["ink_secondary"],
    )
    for lbl in ax.get_xticklabels():
        cond = LABEL_CONDITION[lbl.get_text()]
        lbl.set_color(style.CONDITION_COLORS[cond])
        lbl.set_fontweight("semibold")
        lbl.set_fontsize(10)

    style.recede_spines(ax)
    ax.set_ylabel("1 − Pearson correlation distance")
    ax.set_xlabel("")
    ax.set_title("Hierarchical clustering of 4 samples (average linkage)")

    fig.text(
        0.5, 0.975,
        f"n_proteins_used = {n_complete:,}; distance = 1−Pearson r on log2 intensity; "
        f"trustworthy gate = n≥6 samples (research1.md §5)",
        ha="center", va="top", fontsize=8.8, color=style.CHROME["ink_secondary"],
    )
    fig.text(
        0.5, 0.925,
        "QC-ONLY — no bootstrap support is possible at n=4; dendrogram "
        "topology is not statistically stable. NOT a claim of reproducible "
        "sample clusters.",
        ha="center", va="top", fontsize=8.3, color="white", wrap=True,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=style.STATUS["critical"], edgecolor="none"),
    )
    style.add_caveat(fig, y=-0.02)
    fig.tight_layout(rect=(0, 0.02, 1, 0.80))
    png, svg = style.save_fig(fig, "sample_dendrogram")
    print(f"[cluster] wrote figures/sample_dendrogram.png (+ .svg)")
    return png, svg


# -----------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------
def record_manifest_entries(n_complete, n_total, pct, dist_df, skip_rows):
    skip_list = [r["analysis"] for r in skip_rows if r["status"] == "SKIP"]
    qc_only_list = [r["analysis"] for r in skip_rows if r["status"] == "QC_ONLY"]

    pca_entry = {
        "file": "figures/pca_qc.png",
        "title": "PCA of 4 samples (QC only)",
        "caption": (
            "PC1 vs PC2 from PCA on the 4 SILAC samples (2 control, 2 "
            "treated), standardized per-protein across samples before "
            "decomposition. QC-ONLY / NOT INTERPRETABLE AS BIOLOGICAL "
            "STRUCTURE: at n=4 the sample covariance matrix is "
            "rank-deficient (rank<=n-1=3, exactly 3 PCs exist by "
            "construction), eigenvalues are biased, and PC1 can be driven "
            "by a single sample (research1.md SS5). The trustworthy gate is "
            "n>=6 samples; this dataset is below it, so PCA is shown only "
            "as a QC sanity check (do technical replicate pairs sit near "
            "each other), not as evidence of biological grouping."
        ),
        "key_numbers": {
            "n_proteins_used": n_complete,
            "n_proteins_total_qc_limma": n_total,
            "pc_variance_explained_pct": {k: round(v, 2) for k, v in pct.items()},
            "gate_status": "QC_ONLY",
            "trustworthy_min_samples": TRUSTWORTHY_MIN_SAMPLES,
            "have_samples": N_SAMPLES,
            "skip_list": skip_list,
            "qc_only_list": qc_only_list,
        },
    }

    dendro_entry = {
        "file": "figures/sample_dendrogram.png",
        "title": "Hierarchical clustering of 4 samples (QC only)",
        "caption": (
            "Dendrogram of the 4 SILAC samples, average-linkage clustering "
            "on 1-minus-Pearson-correlation distance of log2 intensity "
            "(complete proteins only). QC-ONLY / NOT INTERPRETABLE: no "
            "bootstrap support is possible with 4 observations, so "
            "dendrogram topology is not statistically stable (research1.md "
            "SS5). The trustworthy gate is n>=6 samples; shown only as a QC "
            "sanity check that technical replicate pairs are each other's "
            "nearest neighbor, not as a claim of reproducible clusters."
        ),
        "key_numbers": {
            "n_proteins_used": n_complete,
            "n_proteins_total_qc_limma": n_total,
            "distance_metric": "1 - Pearson correlation (log2 intensity)",
            "linkage_method": "average",
            "pairwise_correlation_distance": dist_df.round(4).to_dict(),
            "gate_status": "QC_ONLY",
            "trustworthy_min_samples": TRUSTWORTHY_MIN_SAMPLES,
            "have_samples": N_SAMPLES,
            "skip_list": skip_list,
            "qc_only_list": qc_only_list,
        },
    }

    record_gated_manifest([pca_entry, dendro_entry])
    print(f"\n[manifest] wrote {GATED_MANIFEST_PATH}")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
def main():
    print("=" * 78)
    print("Gated PCA / clustering QC layer (Section 6, item 19)")
    print(f"Current sample sheet: N_SAMPLES={N_SAMPLES}, N_GROUPS={N_GROUPS}, "
          f"N_REPLICATES_PER_GROUP={N_REPLICATES_PER_GROUP}")
    print("=" * 78)

    skip_rows = dispatch()
    skip_df, skip_log_path = write_skip_log(skip_rows)

    df, log2_complete, n_total, n_complete = load_complete_matrix()

    coords_df, variance_df, pc_cols = run_pca(log2_complete, n_complete, n_total)
    _, _, pct = fig_pca(coords_df, variance_df, pc_cols, n_complete, n_total)

    Z, dist_df = run_clustering(log2_complete, n_complete)
    fig_dendrogram(Z, dist_df, n_complete)

    record_manifest_entries(n_complete, n_total, pct, dist_df, skip_rows)

    print("\n" + "=" * 78)
    print("DONE. QC-ONLY outputs -- PCA/clustering below the n>=6 trustworthy "
          "gate are NOT evidence of biological structure. See skip_log.csv "
          "for the forward-path dispatcher's full status table.")
    print("=" * 78)


if __name__ == "__main__":
    main()
