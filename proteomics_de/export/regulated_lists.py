"""UP/DOWN regulated-protein CSVs, for handing to a person rather than QIAGEN.

``ipa_input_full.csv`` (:mod:`proteomics_de.export.ipa_export`) already carries
everything a collaborator asking "what went up and what went down, and by how
much" needs -- but it is one 715-row file sorted by accession, with only
``log2FC`` (not a plain multiplicative fold change) to describe magnitude.
This module reads that file and writes two smaller, human-oriented ones:

``regulated_up.csv`` / ``regulated_down.csv``
    The UP-only and DOWN-only rows, each with a ``fold_change = 2**log2FC``
    column appended and sorted by descending magnitude of change (biggest
    movers first). Gene name is the leftmost column here, not the accession --
    unlike the IPA files, nothing consumes these as a QIAGEN upload, so
    ordering for a human reader wins.

No new logic runs on the underlying numbers: this is a split-and-reformat of
data ``ipa_export.py`` and ``limma_test.R`` already computed and validated.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

# proteomics_de/export/regulated_lists.py -> export -> proteomics_de -> repo root
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_ROOT = _PKG_DIR.parent

#: Output column order: gene name leftmost, for a human reader (not an IPA
#: upload -- contrast with ipa_export.IPA_ID_COL being required leftmost there).
REGULATED_OUTPUT_COLS = [
    "Gene names",
    "UniProt Accession Number",
    "log2FC",
    "fold_change",
    "p_value",
    "adj_p_value",
]

REGULATED_UP_FILENAME = "regulated_up.csv"
REGULATED_DOWN_FILENAME = "regulated_down.csv"


class RegulatedListsError(AssertionError):
    """A regulated-list output violates its contract.

    Subclasses ``AssertionError`` on purpose, matching
    :class:`proteomics_de.export.ipa_export.IpaValidationError`.
    """


def _frozen_counts(path=None):
    path = Path(path) if path else _PKG_DIR / "tests" / "expected" / "frozen_counts.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def add_fold_change(df):
    """Append ``fold_change = 2**log2FC`` -- the same formula limma_test.py
    already uses for its own ``fold_change`` column.
    """
    out = df.copy()
    out["fold_change"] = 2.0 ** out["log2FC"]
    return out


def split_regulated(df_full):
    """Partition an ``ipa_input_full``-shaped frame into (up, down).

    Raises if any row's ``regulated`` is not exactly ``UP``/``DOWN``, or if the
    two outputs are not a total partition of the input (no row dropped, none
    duplicated).
    """
    labels = set(df_full["regulated"].unique())
    if not labels <= {"UP", "DOWN"}:
        raise RegulatedListsError(
            f"unexpected 'regulated' values, expected only UP/DOWN: {labels - {'UP', 'DOWN'}}"
        )

    up = df_full[df_full["regulated"] == "UP"]
    down = df_full[df_full["regulated"] == "DOWN"]
    if len(up) + len(down) != len(df_full):
        raise RegulatedListsError(
            f"UP ({len(up)}) + DOWN ({len(down)}) != input rows ({len(df_full)})"
        )
    return up, down


def sort_by_magnitude(df):
    """Descending |log2FC|, stable sort.

    Stability (``kind="mergesort"``) matters here, not just for tidiness: the
    byte-freeze gate (DECISIONS_LOG D14) requires re-running the pipeline to
    reproduce identical bytes, and the default quicksort is not guaranteed
    stable across ties.
    """
    return (
        df.assign(_abs=df["log2FC"].abs())
        .sort_values("_abs", ascending=False, kind="mergesort")
        .drop(columns="_abs")
        .reset_index(drop=True)
    )


def write_regulated(df, path, columns=REGULATED_OUTPUT_COLS):
    df[columns].to_csv(path, index=False, encoding="utf-8")


def build_regulated_lists(results_dir=None, *, validate=True, counts_path=None):
    """Read ``ipa_input_full.csv``, split, sort, write both CSVs.

    Returns
    -------
    list[pathlib.Path]
        ``[regulated_up.csv, regulated_down.csv]``.
    """
    results_dir = Path(results_dir) if results_dir else _PKG_DIR / "results"
    full_path = results_dir / "ipa_input_full.csv"
    if not full_path.exists():
        raise FileNotFoundError(
            f"{full_path} not found -- run ipa_export.py (which runs foldchange.py "
            "and limma_test.py) first."
        )

    df_full = pd.read_csv(full_path, float_precision="round_trip")
    df_full = add_fold_change(df_full)
    up, down = split_regulated(df_full)
    up = sort_by_magnitude(up)
    down = sort_by_magnitude(down)

    up_path = results_dir / REGULATED_UP_FILENAME
    down_path = results_dir / REGULATED_DOWN_FILENAME
    write_regulated(up, up_path)
    write_regulated(down, down_path)

    if validate:
        counts = _frozen_counts(counts_path)
        if len(up) != counts["n_up"]:
            raise RegulatedListsError(
                f"{up_path}: {len(up)} rows, expected {counts['n_up']} (frozen_counts n_up)"
            )
        if len(down) != counts["n_down"]:
            raise RegulatedListsError(
                f"{down_path}: {len(down)} rows, expected {counts['n_down']} (frozen_counts n_down)"
            )

    return [up_path, down_path]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", default=None,
                     help="directory holding ipa_input_full.csv "
                          "(default: proteomics_de/results)")
    ap.add_argument("--no-sidecars", action="store_true",
                     help="skip writing the .provenance.json sidecars")
    args = ap.parse_args(argv)

    results_dir = Path(args.results_dir) if args.results_dir else _PKG_DIR / "results"

    written = build_regulated_lists(results_dir)
    for p in written:
        n_rows = sum(1 for _ in open(p, encoding="utf-8")) - 1
        print(f"Saved {p} ({n_rows} rows)")

    if not args.no_sidecars:
        if str(_ROOT) not in sys.path:
            sys.path.insert(0, str(_ROOT))
        from proteomics_de.provenance import sidecar

        for p in written:
            print(f"Saved {sidecar(p)}")

    return 0


__all__ = [
    "REGULATED_OUTPUT_COLS",
    "REGULATED_UP_FILENAME",
    "REGULATED_DOWN_FILENAME",
    "RegulatedListsError",
    "add_fold_change",
    "split_regulated",
    "sort_by_magnitude",
    "write_regulated",
    "build_regulated_lists",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
