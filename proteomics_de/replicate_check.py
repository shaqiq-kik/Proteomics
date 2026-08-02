"""Bug 6 — Replicate correlation check (do the two technical runs agree?).

Everything was measured twice (two technical replicates). This module asks, in
two complementary ways, whether the two runs agree:

  A. Fold-change agreement (the Bug 6 headline) — do the two runs' treated-vs-
     control answers (``log2_rep1`` vs ``log2_rep2``) line up? This is what drives
     the up/down calls, so the PASS/WARN verdict is taken from it.
  B. Raw replicate reproducibility (reassuring context) — do the raw intensity
     readings reproduce between the paired runs (log10 control rep1 vs rep2, and
     treated rep1 vs rep2)?

It only *measures and flags*. It changes no values and writes no corrected copy,
consistent with the Bug 5 "keep the pipeline un-normalized" stance (Peng et al.
2024). Fold-change correlations run lower than raw-intensity correlations because
taking a difference amplifies noise — so the fold-change bar below is deliberately
lower than a raw-reproducibility bar would be.

Guardrails honored here:
  * The frozen ``proteomics_analysis.py`` and the existing pipeline outputs
    (foldchange_all.csv, foldchange_all_centered.csv, qc_centering.csv,
    ipa_input.csv, single_condition_proteins.csv, onoff_proteins.csv) are never
    touched. This is a pure read of foldchange_all.csv.
  * No inf/NaN literals are written; absent values are empty cells (Bug 3).
  * Check only — no normalization, no sliding, no re-labeling.
"""

import os
import sys

import numpy as np
import pandas as pd

# Tunable judgment knob: the "good enough" bar for the FOLD-CHANGE correlation
# between the two runs. Fold-change correlations run lower than raw-intensity
# correlations (subtracting control from treated amplifies noise), so this bar is
# deliberately lower than a raw-reproducibility bar. No config.py exists yet, so
# the constant lives at the module top; move it to config.py if one is added.
REPLICATE_FC_R_MIN = 0.50

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:  # importable however the caller was launched
    sys.path.insert(0, _ROOT)

from proteomics_de.config import design  # noqa: E402

DEFAULT_RESULTS_DIR = os.path.join(_HERE, "results")
DEFAULT_FOLDCHANGE_CSV = os.path.join(DEFAULT_RESULTS_DIR, "foldchange_all.csv")

# Raw-intensity column mapping, DERIVED from config/sample_sheet.tsv.
#
# These were hardcoded as CONTROL = 31578/31580, TREATED = 31579/31581 -- the
# assignment DECISIONS_LOG D7 established is INVERTED. Because this module only
# ever *labels* two correlations (it does not compute a contrast), the error was
# silent: the numbers were right and the names were swapped, so
# qc_replicate_correlation.csv reported the treated pair's r = 0.8624 (n=1723)
# under `control_raw_r`, and the report repeated it.
#
# Nothing here needs to know which group is which -- it correlates replicate 1
# against replicate 2 within each group -- so deriving the mapping removes the
# only thing that could be wrong.
CONTROL_RAW = tuple(design.control_columns())  # control rep1, rep2
TREATED_RAW = tuple(design.treated_columns())  # treated rep1, rep2


def _pearson(a, b):
    """Pearson r of two finite, equal-length arrays; None if not computable."""
    if len(a) < 2:
        return None
    r = np.corrcoef(a, b)[0, 1]
    return float(r) if np.isfinite(r) else None


def _raw_log10_r(df, col1, col2):
    """log10 Pearson r over rows where BOTH paired raw intensities are > 0.

    Returns ``(r, n)`` where r is None if it cannot be computed. The > 0 filter
    is per-pair, so n differs between the control and treated pairs.
    """
    x = pd.to_numeric(df[col1], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(df[col2], errors="coerce").to_numpy(dtype=float)
    mask = (x > 0) & (y > 0)
    n = int(mask.sum())
    r = _raw = None
    if n >= 2:
        r = _pearson(np.log10(x[mask]), np.log10(y[mask]))
    return r, n


def _fmt(value, ndigits):
    """Round to ndigits, or "" if value is None/non-finite (Bug 3: no inf/NaN)."""
    if value is None or not np.isfinite(value):
        return ""
    return round(float(value), ndigits)


def _make_figure(rep1, rep2, fc_pearson_r, results_dir):
    """Scatter log2_rep1 (x) vs log2_rep2 (y) with a dashed y=x reference line."""
    import matplotlib

    matplotlib.use("Agg")  # headless / reproducible
    import matplotlib.pyplot as plt

    lo = float(min(rep1.min(), rep2.min()))
    hi = float(max(rep1.max(), rep2.max()))

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(rep1, rep2, s=8, alpha=0.4, edgecolors="none")
    ax.plot([lo, hi], [lo, hi], "k--", linewidth=1, label="y = x (perfect agreement)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel("log2 fold-change, run 1 (log2_rep1)")
    ax.set_ylabel("log2 fold-change, run 2 (log2_rep2)")
    ax.set_title("Replicate fold-change agreement")
    r_text = "n/a" if fc_pearson_r is None else f"{fc_pearson_r:.3f}"
    ax.annotate(
        f"Pearson r = {r_text}\nn = {len(rep1)}",
        xy=(0.04, 0.96), xycoords="axes fraction",
        va="top", ha="left",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()

    fig_path = os.path.join(results_dir, "replicate_correlation.png")
    fig.savefig(fig_path, dpi=150)
    plt.close(fig)
    return fig_path


def run_replicate_correlation(foldchange_csv_path, results_dir):
    """Measure whether the two technical runs agree; record a QC note + figure.

    Reads ``foldchange_csv_path`` (pure read). Always writes
    ``qc_replicate_correlation.csv``; also writes ``replicate_correlation.png``.
    Mutates none of the existing pipeline outputs. Returns the metrics dict.
    """
    df = pd.read_csv(foldchange_csv_path)

    # A. Fold-change agreement — over rows where BOTH per-run log2 FCs are finite
    #    (the complete rows). This is the headline the verdict is taken from.
    rep1_all = pd.to_numeric(df["log2_rep1"], errors="coerce")
    rep2_all = pd.to_numeric(df["log2_rep2"], errors="coerce")
    fc_mask = np.isfinite(rep1_all) & np.isfinite(rep2_all)
    rep1 = rep1_all[fc_mask].to_numpy(dtype=float)
    rep2 = rep2_all[fc_mask].to_numpy(dtype=float)

    n_fc = int(fc_mask.sum())
    fc_pearson_r = _pearson(rep1, rep2)

    same_direction = np.sign(rep1) == np.sign(rep2)
    fc_sign_agree_n = int(same_direction.sum())
    fc_sign_agree_pct = (100.0 * fc_sign_agree_n / n_fc) if n_fc else None

    # B. Raw reproducibility — context only.
    control_raw_r, control_raw_n = _raw_log10_r(df, *CONTROL_RAW)
    treated_raw_r, treated_raw_n = _raw_log10_r(df, *TREATED_RAW)

    # Verdict is on the fold-change correlation (section A), since that drives the
    # up/down calls. A non-computable r is treated as a WARN (cannot confirm PASS).
    verdict = "PASS" if (fc_pearson_r is not None and fc_pearson_r >= REPLICATE_FC_R_MIN) else "WARN"

    os.makedirs(results_dir, exist_ok=True)

    # Output 1 (always): the persisted QC record. Floats to 4 dp, percent to 1 dp;
    # empty cells (never inf/NaN literals) where a value is not computable.
    qc_path = os.path.join(results_dir, "qc_replicate_correlation.csv")
    pd.DataFrame(
        [{
            "n_fc": n_fc,
            "fc_pearson_r": _fmt(fc_pearson_r, 4),
            "fc_sign_agree_n": fc_sign_agree_n,
            "fc_sign_agree_pct": _fmt(fc_sign_agree_pct, 1),
            "control_raw_r": _fmt(control_raw_r, 4),
            "control_raw_n": control_raw_n,
            "treated_raw_r": _fmt(treated_raw_r, 4),
            "treated_raw_n": treated_raw_n,
            "fc_r_min": REPLICATE_FC_R_MIN,
            "verdict": verdict,
        }]
    ).to_csv(qc_path, index=False, encoding="utf-8")

    # Output 2 (recommended): the QC scatter figure.
    fig_path = _make_figure(rep1, rep2, fc_pearson_r, results_dir)

    # Console summary (shape mirrors the design doc).
    r_str = "n/a" if fc_pearson_r is None else f"{fc_pearson_r:.3f}"
    pct_str = "n/a" if fc_sign_agree_pct is None else f"{fc_sign_agree_pct:.1f}%"
    if verdict == "PASS":
        verdict_note = "runs agree"
    else:
        verdict_note = "runs agree only weakly"

    def _ctx(r, n):
        return f"r = {r:.3f}  (n={n})" if r is not None else f"r = n/a  (n={n})"

    print("[Bug 6] Replicate correlation check")
    print("  fold-change agreement (run1 vs run2)")
    print(f"    proteins         : {n_fc}")
    print(f"    Pearson r        : {r_str}")
    print(f"    same direction   : {fc_sign_agree_n} / {n_fc}  ({pct_str})")
    print(f"    threshold        : {REPLICATE_FC_R_MIN:.2f}")
    print(f"    verdict          : {verdict}  ({verdict_note})")
    print("  raw reproducibility (context)")
    print(f"    control rep1 vs rep2 : {_ctx(control_raw_r, control_raw_n)}")
    print(f"    treated rep1 vs rep2 : {_ctx(treated_raw_r, treated_raw_n)}")
    print(
        f"  wrote: {os.path.relpath(qc_path, _ROOT)}, "
        f"{os.path.relpath(fig_path, _ROOT)}"
    )

    return {
        "n_fc": n_fc,
        "fc_pearson_r": fc_pearson_r,
        "fc_sign_agree_n": fc_sign_agree_n,
        "fc_sign_agree_pct": fc_sign_agree_pct,
        "control_raw_r": control_raw_r,
        "control_raw_n": control_raw_n,
        "treated_raw_r": treated_raw_r,
        "treated_raw_n": treated_raw_n,
        "fc_r_min": REPLICATE_FC_R_MIN,
        "verdict": verdict,
    }


if __name__ == "__main__":
    run_replicate_correlation(DEFAULT_FOLDCHANGE_CSV, DEFAULT_RESULTS_DIR)
