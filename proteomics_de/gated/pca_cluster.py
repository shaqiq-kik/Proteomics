"""
Section 6, item 19 -- gated PCA / clustering QC layer.

Reads proteomics_de/results/qc_limma.csv (read-only, 1938 both-condition
proteins) and writes NEW files only, under proteomics_de/results/gated/ and
proteomics_de/results/figures/ (a separate figures_manifest_gated.json; the
existing report manifests are never touched).

-----------------------------------------------------------------------------
SCOPE -- read this before trusting anything this script prints
-----------------------------------------------------------------------------
The design is read from config/sample_sheet.tsv via config.design -- nothing
about the sample count is hardcoded here. For the committed sheet that is a
SILAC design with n=2 technical replicates per condition, 4 samples total, and
the paragraphs below describe that case. Per research1.md SS5 ("Limitations as
Engineering Constraints"):

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
`REGISTRY` declares, for every candidate analysis, the sample-count gates from
research1.md SS5 / SS "FORWARD-PATH SECTION". `evaluate_registry()` is a pure
function of (registry, n_samples, n_replicates_per_group) and assigns each
analysis one of three statuses:

  SKIP     -- below the run gate; the analysis does not execute at all.
  QC_ONLY  -- runs, but below the trustworthy gate; output is captioned as a
              sanity check and is NOT interpretable as biology.
  RUN      -- meets the trustworthy gate; output is interpretable on its own
              terms.

`dispatch()` is a thin wrapper that evaluates `REGISTRY` against the CURRENT
sample sheet. N_SAMPLES / N_GROUPS / N_REPLICATES_PER_GROUP are read from
config/sample_sheet.tsv through config.design, so appending replicate rows to
that TSV is sufficient to flip a gate -- no edit to this file is required.
`main()` then executes exactly the analyses whose status is RUN or QC_ONLY.

Scope of "switches itself on", stated precisely so this docstring does not
overclaim: only `pca` and `hierarchical_clustering` are IMPLEMENTED here, and
those two do switch themselves on and off from the sheet. The other four
registry entries (wgcna, umap_tsne, mixomics_splsda, vae_gnn) have no
implementation yet -- growing the sheet moves them out of SKIP in skip_log.csv,
which is the signal to go build them, but it does not conjure the analysis.
`tests/test_gating.py` exercises the dispatcher at n=4, n=6 and n=20 so the
flip is verified rather than assumed.
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
REPO_ROOT = PROTEOMICS_DE_DIR.parent
VIZ_DIR = PROTEOMICS_DE_DIR / "viz"
sys.path.insert(0, str(VIZ_DIR))
import style  # noqa: E402  (palette, apply_style(), save_fig, CHROME, etc.)

# ----------------------------------------------------------------------------
# config package bootstrap (same pattern as viz/style.py)
# ----------------------------------------------------------------------------
# This script is run as a file path (`python proteomics_de/gated/pca_cluster.py`),
# so sys.path[0] is *this* directory and the repo root is not importable. Try the
# normal import first, and only fall back to extending sys.path. The path comes
# from `__file__`, never from the cwd.
try:  # pragma: no cover - exercised by both branches across entry points
    from proteomics_de.config import design as design_cfg  # noqa: E402
except ImportError:  # pragma: no cover
    if str(REPO_ROOT) not in sys.path:
        sys.path.append(str(REPO_ROOT))
    from proteomics_de.config import design as design_cfg  # noqa: E402

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
# The current design, read from config/sample_sheet.tsv (config.design resolves
# the sheet from its own __file__, never from the cwd). For the committed sheet
# these are 4 / 2 / 2. Adding biological replicates means appending rows to that
# TSV; nothing in this file needs editing for the gates below to re-evaluate.
# -----------------------------------------------------------------------
N_SAMPLES = design_cfg.n_samples()
N_GROUPS = design_cfg.n_groups()
N_REPLICATES_PER_GROUP = design_cfg.replicates_per_group()

# Gate for "trustworthy" (not just runnable) PCA/clustering, per research1.md SS5.
TRUSTWORTHY_MIN_SAMPLES = 6
TRUSTWORTHY_MIN_REPLICATES_PER_GROUP = 3

# Statuses emitted by evaluate_registry(), in increasing order of confidence.
SKIP = "SKIP"          # below the run gate -- the analysis does not execute
QC_ONLY = "QC_ONLY"    # runs, but below the trustworthy gate -- not interpretable
RUN = "RUN"            # meets the trustworthy gate -- interpretable
STATUSES = (SKIP, QC_ONLY, RUN)

#: Statuses for which main() actually executes the analysis.
EXECUTING_STATUSES = (QC_ONLY, RUN)

# -----------------------------------------------------------------------
# Registry: every candidate analysis declares min_samples and
# min_replicates_per_group. This IS the forward-path mechanism: grow the
# sample sheet and the statuses re-evaluate, with no code change here.
#
# The `reason_*` fields are FORMAT TEMPLATES, not finished sentences. They are
# rendered by evaluate_registry() against the n actually being evaluated, so a
# reason can never go stale by claiming "have 4" at n=20. Placeholders resolve
# against the entry's own keys plus `n_samples`, `n_replicates_per_group` and
# `rank` (= n_samples - 1).
# -----------------------------------------------------------------------
#: sPLS-DA wants samples >> classes; approximated as this many per class.
MIXOMICS_MIN_SAMPLES_PER_CLASS = 10

REASON_MEETS_GATE = "n_samples/replicates meet the trustworthy gate."


def build_registry(n_groups=None):
    """The candidate-analysis registry, as a list of gate declarations.

    Only `n_groups` is baked in (sPLS-DA's gate scales with the number of
    classes); the sample-count gates are pure research1.md thresholds. Kept a
    function rather than a literal so tests can build a registry for a design
    with a different number of groups.
    """
    if n_groups is None:
        n_groups = N_GROUPS
    mixomics_min = MIXOMICS_MIN_SAMPLES_PER_CLASS * n_groups
    return [
        dict(
            analysis="pca",
            min_samples_to_run=4,
            min_replicates_per_group_to_run=2,
            trustworthy_min_samples=TRUSTWORTHY_MIN_SAMPLES,
            trustworthy_min_replicates_per_group=TRUSTWORTHY_MIN_REPLICATES_PER_GROUP,
            reason_qc_only=(
                "n={n_samples} < trustworthy gate n>={trustworthy_min_samples} "
                "(research1.md SS5): the sample "
                "covariance matrix is rank-deficient (rank<=n-1={rank}, exactly {rank} PCs "
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
                "n={n_samples} < trustworthy gate n>={trustworthy_min_samples} "
                "(research1.md SS5): hierarchical "
                "clustering of {n_samples} samples has no bootstrap support possible -- "
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
                "WGCNA co-expression modules require >={min_samples_to_run} samples (ideally "
                ">={trustworthy_min_samples}) for stable network topology; have {n_samples}."
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
                "low-dimensional embedding (approximated here as "
                ">={min_samples_to_run}); have {n_samples}."
            ),
        ),
        dict(
            analysis="mixomics_splsda",
            min_samples_to_run=mixomics_min,
            min_replicates_per_group_to_run=MIXOMICS_MIN_SAMPLES_PER_CLASS,
            trustworthy_min_samples=mixomics_min,
            trustworthy_min_replicates_per_group=MIXOMICS_MIN_SAMPLES_PER_CLASS,
            n_groups=n_groups,
            min_samples_per_class=MIXOMICS_MIN_SAMPLES_PER_CLASS,
            reason_qc_only=None,
            reason_skip=(
                "mixOmics sPLS-DA requires samples >> classes (approximated "
                "here as >={min_samples_per_class}x the {n_groups} classes = "
                "{min_samples_to_run}); have {n_samples}."
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
                "without overfitting (approximated here as "
                ">={min_samples_to_run}); have {n_samples}."
            ),
        ),
    ]


#: The registry for the committed design. `dispatch()` evaluates this one.
REGISTRY = build_registry()


def _render_reason(template, entry, n_samples, n_replicates_per_group):
    """Fill a registry reason template. Never returns an empty reason.

    A registry entry may legitimately have no template for a status it was
    never expected to reach (e.g. `pca` has no `reason_skip`, because at any
    n>=4 it runs). Rather than emit a blank cell in skip_log.csv, fall back to
    a reason generated from the entry's own gates -- so every row of the log
    always says why, at any sample count.
    """
    context = dict(
        entry,
        n_samples=n_samples,
        n_replicates_per_group=n_replicates_per_group,
        rank=max(n_samples - 1, 0),
    )
    if template:
        return template.format(**context)
    return (
        "{analysis} requires >={min_samples_to_run} samples and "
        ">={min_replicates_per_group_to_run} replicates per group; have "
        "{n_samples} samples ({n_replicates_per_group} per group)."
    ).format(**context)


def evaluate_registry(registry, n_samples, n_replicates_per_group):
    """Assign every registry entry a status at the given sample counts.

    Pure function -- no I/O, no module state. This is what makes the
    forward-path claim testable: `tests/test_gating.py` calls it at n=6 and
    n=20 to prove the gates really do flip, which cannot be observed from the
    committed n=4 design alone.

    Returns a list of row dicts with the skip_log.csv schema: analysis,
    min_samples, have_samples, status, reason.

    `min_samples` in the emitted row is the gate that is actually binding
    for that row's status: the trustworthy-gate for QC_ONLY rows (the
    number that, once met, would upgrade the analysis out of QC-only), and
    the run-gate for SKIP rows (the number that must be met before the
    analysis executes at all). For RUN rows both gates are already met, so
    the trustworthy gate is reported as the one that was cleared last.
    """
    rows = []
    for entry in registry:
        can_run = (
            n_samples >= entry["min_samples_to_run"]
            and n_replicates_per_group >= entry["min_replicates_per_group_to_run"]
        )
        if not can_run:
            status = SKIP
            min_samples = entry["min_samples_to_run"]
            template = entry.get("reason_skip")
        else:
            trustworthy = (
                n_samples >= entry["trustworthy_min_samples"]
                and n_replicates_per_group >= entry["trustworthy_min_replicates_per_group"]
            )
            status = RUN if trustworthy else QC_ONLY
            min_samples = entry["trustworthy_min_samples"]
            template = REASON_MEETS_GATE if trustworthy else entry.get("reason_qc_only")
        rows.append(dict(
            analysis=entry["analysis"],
            min_samples=min_samples,
            have_samples=n_samples,
            status=status,
            reason=_render_reason(template, entry, n_samples, n_replicates_per_group),
        ))
    return rows


def dispatch():
    """Evaluate REGISTRY against the current sample sheet.

    Thin wrapper -- all the logic lives in :func:`evaluate_registry`; the only
    thing this adds is the design read from config/sample_sheet.tsv.
    """
    return evaluate_registry(REGISTRY, N_SAMPLES, N_REPLICATES_PER_GROUP)


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
    proteins present in ALL samples ("complete") are used for PCA /
    clustering, per the task spec -- drop rows with any missing among the
    N_SAMPLES columns and report how many were used.
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
    print(f"[load] complete ({N_SAMPLES}/{N_SAMPLES} intensities present): "
          f"{n_complete} / {n_total} used for PCA/clustering")

    return df, log2_complete, n_total, n_complete


SAMPLE_ORDER = [style.SAMPLE_LABELS[c] for c in style.SAMPLE_COLS]  # control_1, control_2, treated_1, treated_2
LABEL_CONDITION = {style.SAMPLE_LABELS[c]: style.SAMPLE_CONDITION[c] for c in style.SAMPLE_COLS}


# -----------------------------------------------------------------------
# PCA
# -----------------------------------------------------------------------
def run_pca(log2_complete, n_complete, n_total):
    """PCA on the samples named by the sample sheet.

    Design choice (documented): samples are the observations, proteins are
    the features. Each protein is standardized (z-scored across the
    samples, zero mean / unit variance) before PCA -- this is the
    correlation-based variant of PCA rather than raw-covariance PCA, and is
    the standard choice for expression data where proteins span several
    orders of magnitude in absolute intensity; without it, PCA would be
    dominated by the handful of highest-abundance proteins rather than by
    which proteins actually vary in a coordinated way across samples. This
    choice does not fix the small-n degeneracy from research1.md SS5
    (rank<=n-1, biased eigenvalues) -- it only controls for scale.
    """
    X = log2_complete.T  # n_samples (rows) x n_complete proteins (columns)
    X = X.loc[SAMPLE_ORDER]  # fixed sample order

    scaler = StandardScaler()  # per-protein (per-column) z-score across samples
    X_scaled = scaler.fit_transform(X.values)

    # Centering costs one degree of freedom, so the sample covariance matrix has
    # rank <= n_samples - 1: that many PCs exist, and no more (research1.md:246).
    # For the committed 4-sample sheet this is exactly 3.
    n_components = min(X_scaled.shape[0] - 1, X_scaled.shape[1])
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
          f"by construction, direct evidence of the n={N_SAMPLES} degeneracy)")

    return coords_df, variance_df, pc_cols


def gate_banner(status, qc_only_text):
    """(text, facecolor) for a figure's gate banner.

    QC_ONLY figures carry the analysis-specific degeneracy warning; RUN figures
    say so plainly instead of shouting a caveat that no longer applies. Keeping
    this in one place means the banner can never contradict skip_log.csv.
    """
    if status == RUN:
        return (
            f"n={N_SAMPLES} meets the trustworthy gate "
            f"(n≥{TRUSTWORTHY_MIN_SAMPLES} samples, "
            f"≥{TRUSTWORTHY_MIN_REPLICATES_PER_GROUP} per group; research1.md §5).",
            style.STATUS["good"],
        )
    return qc_only_text, style.STATUS["critical"]


def fig_pca(coords_df, variance_df, pc_cols, n_complete, n_total, status=QC_ONLY):
    pct = {row["PC"]: row["variance_explained"] * 100 for _, row in variance_df.iterrows()}
    rank = len(pc_cols)

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
    ax.set_title(f"PCA of {N_SAMPLES} samples")

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=10,
                   markerfacecolor=style.CONDITION_COLORS["control"],
                   markeredgecolor=style.CHROME["ink_primary"]),
        plt.Line2D([0], [0], marker="o", linestyle="", markersize=10,
                   markerfacecolor=style.CONDITION_COLORS["treated"],
                   markeredgecolor=style.CHROME["ink_primary"]),
    ]
    ax.legend(handles, ["Control", "Treated"], loc="best", title="Condition")

    banner_text, banner_color = gate_banner(
        status,
        f"QC-ONLY — n={N_SAMPLES} is statistically degenerate for PCA: rank≤n−1={rank}, "
        f"only {rank} PCs exist, eigenvalues biased, PC1 can be driven by 1 sample. "
        f"NOT evidence of biological structure.",
    )
    fig.text(
        0.5, 0.975,
        f"n_proteins_used = {n_complete:,} of {n_total:,} (complete across all {N_SAMPLES} samples); "
        f"trustworthy gate = n≥{TRUSTWORTHY_MIN_SAMPLES} samples (research1.md §5)",
        ha="center", va="top", fontsize=8.8, color=style.CHROME["ink_secondary"],
    )
    fig.text(
        0.5, 0.925,
        banner_text,
        ha="center", va="top", fontsize=8.3, color="white", wrap=True,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=banner_color, edgecolor="none"),
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
    """Hierarchical clustering of the samples named by the sample sheet.

    Distance: 1 - Pearson correlation between samples' log2 intensity
    vectors (scipy pdist metric='correlation'), computed on the same
    complete-protein matrix used for PCA.
    Linkage: average (UPGMA) -- the standard choice for correlation
    distances; Ward's method assumes (squared) Euclidean distance and is
    not appropriate here, and single linkage is prone to chaining at the
    handful of observations this design provides.
    """
    X = log2_complete.T.loc[SAMPLE_ORDER]  # n_samples x n_complete proteins
    dist_condensed = pdist(X.values, metric="correlation")  # 1 - Pearson r
    dist_square = squareform(dist_condensed)
    dist_df = pd.DataFrame(dist_square, index=SAMPLE_ORDER, columns=SAMPLE_ORDER)

    Z = linkage(dist_condensed, method="average")

    print(f"\n[cluster] n_proteins_used={n_complete}")
    print("[cluster] pairwise 1-Pearson-correlation distance matrix:")
    print(dist_df.round(4).to_string())

    return Z, dist_df


def fig_dendrogram(Z, dist_df, n_complete, status=QC_ONLY):
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
    ax.set_title(f"Hierarchical clustering of {N_SAMPLES} samples (average linkage)")

    banner_text, banner_color = gate_banner(
        status,
        f"QC-ONLY — no bootstrap support is possible at n={N_SAMPLES}; dendrogram "
        f"topology is not statistically stable. NOT a claim of reproducible "
        f"sample clusters.",
    )
    fig.text(
        0.5, 0.975,
        f"n_proteins_used = {n_complete:,}; distance = 1−Pearson r on log2 intensity; "
        f"trustworthy gate = n≥{TRUSTWORTHY_MIN_SAMPLES} samples (research1.md §5)",
        ha="center", va="top", fontsize=8.8, color=style.CHROME["ink_secondary"],
    )
    fig.text(
        0.5, 0.925,
        banner_text,
        ha="center", va="top", fontsize=8.3, color="white", wrap=True,
        bbox=dict(boxstyle="round,pad=0.5", facecolor=banner_color, edgecolor="none"),
    )
    style.add_caveat(fig, y=-0.02)
    fig.tight_layout(rect=(0, 0.02, 1, 0.80))
    png, svg = style.save_fig(fig, "sample_dendrogram")
    print(f"[cluster] wrote figures/sample_dendrogram.png (+ .svg)")
    return png, svg


# -----------------------------------------------------------------------
# Manifest
# -----------------------------------------------------------------------
#: Human-readable suffix for a figure title, per gate status.
TITLE_SUFFIX = {QC_ONLY: "QC only", RUN: "trustworthy"}


def _group_composition():
    """e.g. "2 control, 2 treated" -- read from the sample sheet, not assumed."""
    return ", ".join(
        f"{N_REPLICATES_PER_GROUP} {group}" for group in design_cfg.group_names()
    )


def record_manifest_entries(n_complete, n_total, pct, dist_df, skip_rows):
    status_by_analysis = {r["analysis"]: r["status"] for r in skip_rows}
    skip_list = [r["analysis"] for r in skip_rows if r["status"] == SKIP]
    qc_only_list = [r["analysis"] for r in skip_rows if r["status"] == QC_ONLY]

    pca_status = status_by_analysis["pca"]
    cluster_status = status_by_analysis["hierarchical_clustering"]
    rank = len(pct)
    composition = _group_composition()

    # Captions are assembled per status so a figure can never carry a caveat
    # that its own gate status contradicts.
    if pca_status == RUN:
        pca_tail = (
            f"Meets the trustworthy gate (n>={TRUSTWORTHY_MIN_SAMPLES} samples, "
            f">={TRUSTWORTHY_MIN_REPLICATES_PER_GROUP} per group; research1.md "
            f"SS5), so the decomposition is interpretable on its own terms. "
            f"Note that rank<=n-1 still holds: exactly {rank} PCs exist."
        )
    else:
        pca_tail = (
            f"QC-ONLY / NOT INTERPRETABLE AS BIOLOGICAL "
            f"STRUCTURE: at n={N_SAMPLES} the sample covariance matrix is "
            f"rank-deficient (rank<=n-1={rank}, exactly {rank} PCs exist by "
            f"construction), eigenvalues are biased, and PC1 can be driven "
            f"by a single sample (research1.md SS5). The trustworthy gate is "
            f"n>={TRUSTWORTHY_MIN_SAMPLES} samples; this dataset is below it, so PCA is shown only "
            f"as a QC sanity check (do technical replicate pairs sit near "
            f"each other), not as evidence of biological grouping."
        )

    pca_entry = {
        "file": "figures/pca_qc.png",
        "title": f"PCA of {N_SAMPLES} samples ({TITLE_SUFFIX[pca_status]})",
        "caption": (
            f"PC1 vs PC2 from PCA on the {N_SAMPLES} SILAC samples ({composition}"
            f"), standardized per-protein across samples before "
            f"decomposition. " + pca_tail
        ),
        "key_numbers": {
            "n_proteins_used": n_complete,
            "n_proteins_total_qc_limma": n_total,
            "pc_variance_explained_pct": {k: round(v, 2) for k, v in pct.items()},
            "gate_status": pca_status,
            "trustworthy_min_samples": TRUSTWORTHY_MIN_SAMPLES,
            "have_samples": N_SAMPLES,
            "skip_list": skip_list,
            "qc_only_list": qc_only_list,
        },
    }

    if cluster_status == RUN:
        cluster_tail = (
            f"Meets the trustworthy gate (n>={TRUSTWORTHY_MIN_SAMPLES} samples, "
            f">={TRUSTWORTHY_MIN_REPLICATES_PER_GROUP} per group; research1.md "
            f"SS5), so the topology may be read as sample structure."
        )
    else:
        cluster_tail = (
            f"QC-ONLY / NOT INTERPRETABLE: no "
            f"bootstrap support is possible with {N_SAMPLES} observations, so "
            f"dendrogram topology is not statistically stable (research1.md "
            f"SS5). The trustworthy gate is n>={TRUSTWORTHY_MIN_SAMPLES} samples; shown only as a QC "
            f"sanity check that technical replicate pairs are each other's "
            f"nearest neighbor, not as a claim of reproducible clusters."
        )

    dendro_entry = {
        "file": "figures/sample_dendrogram.png",
        "title": (
            f"Hierarchical clustering of {N_SAMPLES} samples "
            f"({TITLE_SUFFIX[cluster_status]})"
        ),
        "caption": (
            f"Dendrogram of the {N_SAMPLES} SILAC samples, average-linkage clustering "
            f"on 1-minus-Pearson-correlation distance of log2 intensity "
            f"(complete proteins only). " + cluster_tail
        ),
        "key_numbers": {
            "n_proteins_used": n_complete,
            "n_proteins_total_qc_limma": n_total,
            "distance_metric": "1 - Pearson correlation (log2 intensity)",
            "linkage_method": "average",
            "pairwise_correlation_distance": dist_df.round(4).to_dict(),
            "gate_status": cluster_status,
            "trustworthy_min_samples": TRUSTWORTHY_MIN_SAMPLES,
            "have_samples": N_SAMPLES,
            "skip_list": skip_list,
            "qc_only_list": qc_only_list,
        },
    }

    # Only describe figures this run actually produced. A SKIPped analysis
    # leaves whatever the manifest already said about it untouched.
    entries = []
    if pca_status in EXECUTING_STATUSES:
        entries.append(pca_entry)
    if cluster_status in EXECUTING_STATUSES:
        entries.append(dendro_entry)

    record_gated_manifest(entries)
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

    # The dispatcher decides what runs. Nothing below is unconditional: grow
    # config/sample_sheet.tsv past a gate and the corresponding analysis starts
    # executing (and loses its QC-only caveats) without an edit to this file.
    status_by_analysis = {r["analysis"]: r["status"] for r in skip_rows}
    pca_status = status_by_analysis["pca"]
    cluster_status = status_by_analysis["hierarchical_clustering"]

    df, log2_complete, n_total, n_complete = load_complete_matrix()

    pct = {}
    dist_df = pd.DataFrame()

    if pca_status in EXECUTING_STATUSES:
        coords_df, variance_df, pc_cols = run_pca(log2_complete, n_complete, n_total)
        _, _, pct = fig_pca(
            coords_df, variance_df, pc_cols, n_complete, n_total, status=pca_status
        )
    else:
        print(f"\n[pca] SKIPPED by the dispatcher ({status_by_analysis['pca']}).")

    if cluster_status in EXECUTING_STATUSES:
        Z, dist_df = run_clustering(log2_complete, n_complete)
        fig_dendrogram(Z, dist_df, n_complete, status=cluster_status)
    else:
        print(f"\n[cluster] SKIPPED by the dispatcher ({cluster_status}).")

    record_manifest_entries(n_complete, n_total, pct, dist_df, skip_rows)

    print("\n" + "=" * 78)
    if QC_ONLY in (pca_status, cluster_status):
        print(f"DONE. QC-ONLY outputs -- PCA/clustering below the "
              f"n>={TRUSTWORTHY_MIN_SAMPLES} trustworthy "
              f"gate are NOT evidence of biological structure. See skip_log.csv "
              f"for the forward-path dispatcher's full status table.")
    else:
        print(f"DONE. PCA/clustering meet the n>={TRUSTWORTHY_MIN_SAMPLES} "
              f"trustworthy gate. See skip_log.csv for the forward-path "
              f"dispatcher's full status table.")
    print("=" * 78)


if __name__ == "__main__":
    main()
