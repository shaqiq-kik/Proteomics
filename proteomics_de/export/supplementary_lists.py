"""Supplementary regulated-adjacent CSVs for the two missing-data tiers that
``regulated_up.csv``/``regulated_down.csv`` structurally -- and correctly --
exclude (see DECISIONS_LOG D17).

``regulated_up.csv``/``regulated_down.csv`` are gated on ``complete=True``
(``etl.foldchange_core.mark_complete``), decided BEFORE limma ever runs. Two
tiers of real signal never cross that gate:

* Tier 3 ("partial missingness"): ``complete=False`` because 1-2 of 4 raw
  replicates are zero/missing, so the raw log2FC/regulated call never runs --
  but limma DOES test these rows after standard MinProb imputation, and
  ``qc_limma.csv`` already carries a legitimate ``limma_log2FC``/``p_value``/
  ``adj_p_value``/``n_imputed`` for every one of them.
* Tier 1/2 ("fully undetected"): the protein was never identified at all on
  one side -- either its whole row is absent from one raw MaxQuant sheet
  (``single_condition_proteins.csv``) or every replicate on one side reads
  null/zero while the row exists structurally (``onoff_proteins.csv``). These
  are correctly excluded from limma too (testing them would invent an entire
  absent condition) and so never carry a fold change or p-value at all.

This module does NOT touch ``ipa_input.csv``, ``ipa_input_full.csv``,
``regulated_up.csv`` or ``regulated_down.csv`` -- it is additive only, reading
their upstream sources (``qc_limma.csv``, ``single_condition_proteins.csv``,
``onoff_proteins.csv``) and writing three new files:

``regulated_up_partial.csv`` / ``regulated_down_partial.csv``
    The tier-3 rows, reclassified from ``limma_log2FC`` at the existing
    ``config.constants.LOG2_THRESHOLD``. Every row carries ``n_imputed`` (always
    >= 1 by construction of the filter) so a reader can never mistake an
    imputed call for a complete-case measurement -- this is a DIFFERENT
    quantity than the core files' ``log2FC`` (raw mean-of-log2-ratios on
    unimputed data), not just a re-run of the same one.

``qualitative_changes.csv``
    The union of ``single_condition_proteins.csv`` + ``onoff_proteins.csv``,
    with a computed ``direction`` column and NO fold-change or p-value
    column -- neither is statistically valid when one whole condition has
    zero data points (tier 1) or zero variance (tier 2).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# proteomics_de/export/supplementary_lists.py -> export -> proteomics_de -> repo root
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_ROOT = _PKG_DIR.parent

if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from proteomics_de.export.regulated_lists import add_fold_change, sort_by_magnitude  # noqa: E402

#: Output column order for the two tier-3 files. Mirrors
#: ``regulated_lists.REGULATED_OUTPUT_COLS`` with ``n_imputed`` appended --
#: mandatory on every row here, absent there, so the two are never confusable
#: by column shape alone.
PARTIAL_OUTPUT_COLS = [
    "Gene names",
    "UniProt Accession Number",
    "log2FC",
    "fold_change",
    "p_value",
    "adj_p_value",
    "n_imputed",
]

#: Output column order for the qualitative file. Deliberately carries no
#: log2FC/fold_change/p_value/adj_p_value -- none is statistically valid for
#: these rows (see module docstring) -- and carries ``basis``/``note`` instead
#: so the omission reads as a decision, not a missing column.
QUALITATIVE_OUTPUT_COLS = [
    "Gene names",
    "UniProt Accession Number",
    "direction",
    "basis",
    "Intensity 31578",
    "Intensity 31580",
    "Intensity 31579",
    "Intensity 31581",
    "note",
]

REGULATED_UP_PARTIAL_FILENAME = "regulated_up_partial.csv"
REGULATED_DOWN_PARTIAL_FILENAME = "regulated_down_partial.csv"
QUALITATIVE_CHANGES_FILENAME = "qualitative_changes.csv"

#: Caveat text stamped into every qualitative_changes.csv row's ``note``
#: column, split by which of the two source files a row came from.
_SINGLE_CONDITION_NOTE = (
    "qualitative: identified in only one condition -- the protein's entire "
    "row is absent from the other MaxQuant sheet, not just zero. No fold "
    "change or p-value is statistically valid when one condition has zero "
    "data points."
)
_ONOFF_NOTE = (
    "qualitative: identified in both MaxQuant sheets, but every replicate on "
    "one side reads null/zero. No fold change or p-value is statistically "
    "valid when one condition has zero variance."
)

#: Raw intensity columns carried through into qualitative_changes.csv,
#: physical acquisition order (matches limma_test.py's _PHYSICAL_COLS).
_INTENSITY_COLS = ["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581"]


class SupplementaryListsError(AssertionError):
    """A supplementary output violates its contract. Mirrors RegulatedListsError."""


def _frozen_counts(path=None):
    path = Path(path) if path else _PKG_DIR / "tests" / "expected" / "frozen_counts.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# Tier 3: regulated_up_partial.csv / regulated_down_partial.csv
# ---------------------------------------------------------------------------
def select_partial_regulated(qc_limma, threshold):
    """Tier-3 rows (``n_imputed`` > 0) reclassified at `threshold` from ``limma_log2FC``.

    Returns ``(up, down)``, each with columns renamed to the core-file
    convention (``id``->``UniProt Accession Number``, ``gene``->``Gene names``,
    ``limma_log2FC``->``log2FC``) but not yet fold-changed or sorted.
    """
    required = {"id", "gene", "limma_log2FC", "p_value", "adj_p_value", "n_imputed"}
    missing = required - set(qc_limma.columns)
    if missing:
        raise SupplementaryListsError(f"qc_limma frame is missing columns: {sorted(missing)}")

    tier3 = qc_limma[qc_limma["n_imputed"] > 0].copy()
    tier3 = tier3.rename(columns={
        "id": "UniProt Accession Number",
        "gene": "Gene names",
        "limma_log2FC": "log2FC",
    })
    up = tier3[tier3["log2FC"] >= threshold].copy()
    down = tier3[tier3["log2FC"] <= -threshold].copy()
    return up, down


def build_partial_regulated_lists(results_dir=None, *, threshold=None,
                                   validate=True, counts_path=None):
    """Read ``qc_limma.csv``, write ``regulated_up_partial.csv`` / ``regulated_down_partial.csv``.

    Returns
    -------
    list[pathlib.Path]
        ``[regulated_up_partial.csv, regulated_down_partial.csv]``.
    """
    results_dir = Path(results_dir) if results_dir else _PKG_DIR / "results"
    limma_path = results_dir / "qc_limma.csv"
    if not limma_path.exists():
        raise FileNotFoundError(
            f"{limma_path} not found -- run foldchange.py (which runs limma_test.py) first."
        )

    if threshold is None:
        from proteomics_de.config.constants import LOG2_THRESHOLD
        threshold = LOG2_THRESHOLD

    qc = pd.read_csv(limma_path, float_precision="round_trip")
    up, down = select_partial_regulated(qc, threshold)
    up = sort_by_magnitude(add_fold_change(up))
    down = sort_by_magnitude(add_fold_change(down))

    up_path = results_dir / REGULATED_UP_PARTIAL_FILENAME
    down_path = results_dir / REGULATED_DOWN_PARTIAL_FILENAME
    up[PARTIAL_OUTPUT_COLS].to_csv(up_path, index=False, encoding="utf-8")
    down[PARTIAL_OUTPUT_COLS].to_csv(down_path, index=False, encoding="utf-8")

    if validate:
        counts = _frozen_counts(counts_path)
        if len(up) != counts["n_up_partial"]:
            raise SupplementaryListsError(
                f"{up_path}: {len(up)} rows, expected {counts['n_up_partial']} (frozen_counts n_up_partial)"
            )
        if len(down) != counts["n_down_partial"]:
            raise SupplementaryListsError(
                f"{down_path}: {len(down)} rows, expected {counts['n_down_partial']} (frozen_counts n_down_partial)"
            )
        _assert_disjoint_from_core(
            results_dir,
            set(up["UniProt Accession Number"]) | set(down["UniProt Accession Number"]),
            "tier-3 partial", up_path,
        )

    return [up_path, down_path]


# ---------------------------------------------------------------------------
# Tier 1/2: qualitative_changes.csv
# ---------------------------------------------------------------------------
def build_qualitative_changes(results_dir=None, *, validate=True, counts_path=None):
    """Read ``single_condition_proteins.csv`` + ``onoff_proteins.csv``, write
    ``qualitative_changes.csv``.

    Returns
    -------
    list[pathlib.Path]
        ``[qualitative_changes.csv]``.
    """
    results_dir = Path(results_dir) if results_dir else _PKG_DIR / "results"
    sc_path = results_dir / "single_condition_proteins.csv"
    oo_path = results_dir / "onoff_proteins.csv"
    for p in (sc_path, oo_path):
        if not p.exists():
            raise FileNotFoundError(f"{p} not found -- run foldchange.py first.")

    sc = pd.read_csv(sc_path, float_precision="round_trip")
    oo = pd.read_csv(oo_path, float_precision="round_trip")

    sc = sc.rename(columns={"accession": "UniProt Accession Number", "gene": "Gene names"})
    sc["direction"] = np.where(sc["detected_in"] == "treated_only", "UP", "DOWN")
    sc["basis"] = sc["detected_in"]
    sc["note"] = _SINGLE_CONDITION_NOTE

    oo = oo.rename(columns={"accession": "UniProt Accession Number", "gene": "Gene names"})
    oo["direction"] = np.where(oo["onoff"] == "on_with_treatment", "UP", "DOWN")
    oo["basis"] = oo["onoff"]
    oo["note"] = _ONOFF_NOTE

    combined = pd.concat(
        [sc[QUALITATIVE_OUTPUT_COLS], oo[QUALITATIVE_OUTPUT_COLS]], ignore_index=True
    )

    # Deterministic sort: UP before DOWN, then gene name case-insensitive
    # (mergesort for byte-freeze reproducibility -- DECISIONS_LOG D14). Some
    # `Gene names` values are blank/NaN -- filled to "" so they sort first
    # rather than raising.
    combined["_dir_rank"] = combined["direction"].map({"UP": 0, "DOWN": 1})
    combined["_gene_key"] = combined["Gene names"].fillna("").astype(str).str.upper()
    combined = (
        combined.sort_values(["_dir_rank", "_gene_key"], kind="mergesort")
        .drop(columns=["_dir_rank", "_gene_key"])
        .reset_index(drop=True)
    )

    out_path = results_dir / QUALITATIVE_CHANGES_FILENAME
    combined[QUALITATIVE_OUTPUT_COLS].to_csv(out_path, index=False, encoding="utf-8")

    if validate:
        counts = _frozen_counts(counts_path)
        if len(combined) != counts["n_qualitative"]:
            raise SupplementaryListsError(
                f"{out_path}: {len(combined)} rows, expected {counts['n_qualitative']} (frozen_counts n_qualitative)"
            )
        n_up = int((combined["direction"] == "UP").sum())
        n_down = int((combined["direction"] == "DOWN").sum())
        if n_up != counts["n_qualitative_up"] or n_down != counts["n_qualitative_down"]:
            raise SupplementaryListsError(
                f"{out_path}: UP={n_up}/DOWN={n_down}, expected "
                f"UP={counts['n_qualitative_up']}/DOWN={counts['n_qualitative_down']}"
            )
        _assert_disjoint_from_core(
            results_dir, set(combined["UniProt Accession Number"]), "qualitative", out_path
        )

    return [out_path]


def _assert_disjoint_from_core(results_dir, ids, label, path):
    """Best-effort cross-check against ``regulated_up.csv``/``regulated_down.csv``.

    By construction this can never overlap (``complete=False`` rows and
    single/onoff rows never reach ``ipa_input.csv``), but "can never" is
    exactly the kind of claim that deserves a runtime assertion, not just a
    comment. Skipped silently if the core files aren't present yet (e.g. a
    partial results directory in a test fixture).
    """
    up_path = results_dir / "regulated_up.csv"
    down_path = results_dir / "regulated_down.csv"
    if not (up_path.exists() and down_path.exists()):
        return
    core = (
        set(pd.read_csv(up_path, float_precision="round_trip")["UniProt Accession Number"])
        | set(pd.read_csv(down_path, float_precision="round_trip")["UniProt Accession Number"])
    )
    overlap = ids & core
    if overlap:
        raise SupplementaryListsError(
            f"{path}: {label} rows overlap the core regulated set (should be "
            f"impossible by construction): {sorted(overlap)[:10]}"
        )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", default=None,
                     help="directory holding qc_limma.csv / single_condition_proteins.csv "
                          "/ onoff_proteins.csv (default: proteomics_de/results)")
    ap.add_argument("--no-sidecars", action="store_true",
                     help="skip writing the .provenance.json sidecars")
    args = ap.parse_args(argv)

    results_dir = Path(args.results_dir) if args.results_dir else _PKG_DIR / "results"

    written = build_partial_regulated_lists(results_dir) + build_qualitative_changes(results_dir)
    for p in written:
        n_rows = sum(1 for _ in open(p, encoding="utf-8")) - 1
        print(f"Saved {p} ({n_rows} rows)")

    if not args.no_sidecars:
        from proteomics_de.provenance import sidecar

        for p in written:
            print(f"Saved {sidecar(p)}")

    return 0


__all__ = [
    "PARTIAL_OUTPUT_COLS",
    "QUALITATIVE_OUTPUT_COLS",
    "REGULATED_UP_PARTIAL_FILENAME",
    "REGULATED_DOWN_PARTIAL_FILENAME",
    "QUALITATIVE_CHANGES_FILENAME",
    "SupplementaryListsError",
    "select_partial_regulated",
    "build_partial_regulated_lists",
    "build_qualitative_changes",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
