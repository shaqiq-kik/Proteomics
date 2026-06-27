"""Bug 5 — SILAC centering check (QC note + optional 'slid onto zero' copy).

Most proteins should not change between conditions, so the pile of complete-row
log2 fold-changes ought to sit on zero ("no change"). This module *measures*
where that pile actually sits (its median), prints a PASS/WARN note, and — only
when the pile is off-centre — also writes a SEPARATE "centered" copy for visual
inspection. It never edits the real pipeline outputs.

Per Peng et al. 2024 the primary pipeline stays un-normalized: the centered copy
is a "what-if" artifact for judging whether the shift is a technical loading
offset, NOT an input to IPA or to limma (Bug 7).

Guardrails honored here:
  * The frozen ``proteomics_analysis.py`` and the existing outputs
    (foldchange_all.csv, ipa_input.csv, single_condition_proteins.csv,
    onoff_proteins.csv) are never touched.
  * The +/-0.585 regulation threshold is NOT redefined — it is imported from the
    Bug 2 module (single source of truth).
  * No inf/NaN literals are written; absent values are empty cells (Bug 3).
  * Only complete rows (finite log2FC) are centered; ON_OFF and single-condition
    proteins are left alone.
"""

import os
import sys

import numpy as np
import pandas as pd

# Single source of truth for the regulation cutoff: import it from the Bug 2
# module rather than re-typing 0.585. foldchange.py guards its body behind a
# __main__ check, so importing it is side-effect free (no pipeline re-run).
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from proteomics_de.foldchange import LOG2_THRESHOLD as REG_THRESHOLD  # noqa: E402

# Tunable judgment knob: how close to zero counts as "already centered".
# 0.10 in log2 space ~= a 1.07x (~7%) imbalance between the channels.
CENTERING_TOLERANCE_LOG2 = 0.10

DEFAULT_RESULTS_DIR = os.path.join(_HERE, "results")
DEFAULT_FOLDCHANGE_CSV = os.path.join(DEFAULT_RESULTS_DIR, "foldchange_all.csv")


def _classify(log2fc_value):
    """UP / DOWN / NO CHANGE for a finite log2FC, using the shared threshold."""
    if log2fc_value >= REG_THRESHOLD:
        return "UP"
    if log2fc_value <= -REG_THRESHOLD:
        return "DOWN"
    return "NO CHANGE"


def run_centering_check(foldchange_csv_path, results_dir):
    """Measure where the complete-row log2FC pile sits and record a QC note.

    Reads ``foldchange_csv_path``; always writes ``qc_centering.csv`` and, only
    when the pile is off-centre (verdict == WARN), also writes
    ``foldchange_all_centered.csv`` (inspection only). Pure with respect to the
    existing outputs: ``foldchange_all.csv`` is never mutated.

    Returns the metrics dict.
    """
    # Read every cell as a string so the original columns survive byte-for-byte
    # when we re-emit them in the centered copy (keep_default_na=False prevents
    # empty cells from becoming NaN literals — the Bug 3 invariant).
    raw = pd.read_csv(foldchange_csv_path, dtype=str, keep_default_na=False)

    complete_mask = raw["complete"] == "True"
    log2fc = pd.to_numeric(raw.loc[complete_mask, "log2FC"], errors="coerce")

    n_complete = int(complete_mask.sum())
    center_mean = float(log2fc.mean())
    center_median = float(log2fc.median())

    if not np.isfinite(center_median):
        raise ValueError(
            "center_median is not finite; refusing to write any QC output."
        )

    verdict = "PASS" if abs(center_median) <= CENTERING_TOLERANCE_LOG2 else "WARN"
    centered_copy_written = verdict == "WARN"

    os.makedirs(results_dir, exist_ok=True)

    # Output 1 (always): the persisted QC note.
    qc_path = os.path.join(results_dir, "qc_centering.csv")
    pd.DataFrame(
        [{
            "n_complete": n_complete,
            "mean_log2FC": round(center_mean, 4),
            "median_log2FC": round(center_median, 4),
            "tolerance": CENTERING_TOLERANCE_LOG2,
            "verdict": verdict,
            "centered_copy_written": centered_copy_written,
        }]
    ).to_csv(qc_path, index=False, encoding="utf-8")

    # Console summary (shape mirrors the design doc).
    print("[Bug 5] SILAC centering check")
    print(f"  complete rows : {n_complete}")
    print(f"  mean  log2FC  : {center_mean:.4f}")
    print(f"  median log2FC : {center_median:.4f}   (centering value)")
    print(f"  tolerance     : {CENTERING_TOLERANCE_LOG2:.2f}")

    before_counts = {"UP": 0, "DOWN": 0, "NO CHANGE": 0}
    after_counts = {"UP": 0, "DOWN": 0, "NO CHANGE": 0}

    if verdict == "PASS":
        print("  verdict       : PASS  (pile already centered)")
        print("  already centered; no copy written")
        metrics = {
            "n_complete": n_complete,
            "mean_log2FC": round(center_mean, 4),
            "median_log2FC": round(center_median, 4),
            "tolerance": CENTERING_TOLERANCE_LOG2,
            "verdict": verdict,
            "centered_copy_written": centered_copy_written,
        }
        return metrics

    # verdict == WARN -> build and write the centered copy (Output 2).
    log2fc_centered = log2fc - center_median  # finite, complete rows only

    # New column 1: log2FC_centered. Empty string for non-complete rows so no
    # inf/NaN literal is ever written.
    centered_col = pd.Series([""] * len(raw), index=raw.index, dtype=object)
    centered_col.loc[complete_mask] = log2fc_centered.map(lambda v: repr(float(v)))

    # New column 2: regulated_centered. Complete rows are re-classified against
    # the shared threshold; non-complete rows keep their original `regulated`
    # label (ON_OFF stays ON_OFF, NO CHANGE stays NO CHANGE).
    regulated_centered = raw["regulated"].astype(object).copy()
    regulated_centered.loc[complete_mask] = log2fc_centered.map(_classify)

    out = raw.copy()
    out["log2FC_centered"] = centered_col
    out["regulated_centered"] = regulated_centered

    centered_path = os.path.join(results_dir, "foldchange_all_centered.csv")
    out.to_csv(centered_path, index=False, encoding="utf-8")

    # Before/after regulation calls on the complete rows, for the console.
    before = raw.loc[complete_mask, "regulated"]
    after = regulated_centered.loc[complete_mask]
    for key in before_counts:
        before_counts[key] = int((before == key).sum())
        after_counts[key] = int((after == key).sum())

    print("  verdict       : WARN  (pile is off-center)")
    rel = os.path.relpath(centered_path, _ROOT)
    print(f"  wrote: {rel}  (inspection only)")
    print("")
    print("  regulated calls on complete rows   BEFORE   AFTER")
    print(f"    UP                                {before_counts['UP']:>5}   {after_counts['UP']:>5}")
    print(f"    DOWN                              {before_counts['DOWN']:>5}   {after_counts['DOWN']:>5}")
    print(f"    NO CHANGE                         {before_counts['NO CHANGE']:>5}   {after_counts['NO CHANGE']:>5}")

    metrics = {
        "n_complete": n_complete,
        "mean_log2FC": round(center_mean, 4),
        "median_log2FC": round(center_median, 4),
        "tolerance": CENTERING_TOLERANCE_LOG2,
        "verdict": verdict,
        "centered_copy_written": centered_copy_written,
        "before_counts": before_counts,
        "after_counts": after_counts,
    }
    return metrics


if __name__ == "__main__":
    run_centering_check(DEFAULT_FOLDCHANGE_CSV, DEFAULT_RESULTS_DIR)
