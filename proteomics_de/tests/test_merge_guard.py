"""The merge cardinality guard research1.md asks for twice (lines 52 and 183).

Two halves:

* **Negative tests** inject a duplicated accession into a clean fixture and
  assert the guard trips. Without them the guard is untested code that happens
  never to fire, which is indistinguishable from no guard at all.
* **Positive tests** run the guard against the real ``Copy of General
  Sheet.xlsx`` and assert it passes, so a future edit to the guard cannot start
  rejecting the committed dataset.

The real-sheet read is ~0.4 s (two ``read_excel`` calls limited to the accession
and gene columns), so it stays in the default, un-marked test run.
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.etl import foldchange_core as core  # noqa: E402
from proteomics_de.etl import merge_guard  # noqa: E402

ACC = merge_guard.ACCESSION_COL
# Which channels live in which SHEET of the workbook -- a fact about the file,
# independent of which condition each one is. (DECISIONS_LOG D7 corrected the
# condition assignment to 31579/31581 = control, 31578/31580 = treated; it did
# not, and could not, move a channel between sheets. The merge only cares about
# sheet position, so nothing here changes with it.)
SHEET_L_COLS = ["Intensity 31578", "Intensity 31580"]
SHEET_H_COLS = ["Intensity 31579", "Intensity 31581"]

_WORKBOOK = _REPO_ROOT / "Copy of General Sheet.xlsx"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def sheet(accessions, cols, start=100.0):
    """A minimal sheet: one row per accession, with `cols` as its intensities."""
    rows = []
    for i, acc in enumerate(accessions):
        value = start + i
        rows.append((acc, f"Gene{acc}", value, value))
    return pd.DataFrame(rows, columns=[ACC, "Gene names"] + list(cols))


@pytest.fixture
def clean_sheets():
    """Three shared accessions plus one exclusive to each side. No duplicates."""
    df_L = sheet(["A", "B", "C", "L_ONLY"], SHEET_L_COLS)
    df_H = sheet(["A", "B", "C", "H_ONLY"], SHEET_H_COLS)
    return df_L, df_H


def merge_and_split(df_L, df_H):
    merged = core.merge_with_indicator(df_L, df_H)
    both, single = core.split_both_single(merged)
    return merged, both, single


# ---------------------------------------------------------------------------
# duplicate detection
# ---------------------------------------------------------------------------
def test_clean_sheets_report_zero_duplicates(clean_sheets):
    df_L, df_H = clean_sheets

    assert merge_guard.duplicate_accession_count(df_L) == 0
    assert merge_guard.duplicate_accession_count(df_H) == 0
    assert merge_guard.duplicate_accessions(df_L).empty


def test_duplicate_accessions_reports_the_accession_and_its_count():
    df = sheet(["A", "B", "A", "C", "A"], SHEET_L_COLS)
    dupes = merge_guard.duplicate_accessions(df)

    assert dupes.to_dict() == {"A": 3}
    assert merge_guard.duplicate_accession_count(df) == 1


def test_missing_accessions_are_not_counted_as_duplicates():
    """Two NaN keys are a schema problem, not a cardinality problem.

    ``pd.merge`` does not join NaN to NaN, so missing keys cannot cause the row
    explosion this guard exists to catch. Counting them here would fire the wrong
    alarm and point at the wrong fix.
    """
    df = pd.DataFrame(
        [(np.nan, "g1", 1.0, 1.0), (np.nan, "g2", 2.0, 2.0), ("A", "g3", 3.0, 3.0)],
        columns=[ACC, "Gene names"] + SHEET_L_COLS,
    )
    assert merge_guard.duplicate_accession_count(df) == 0


# ---------------------------------------------------------------------------
# assert_merge_safe -- the guard trips
# ---------------------------------------------------------------------------
def test_guard_trips_on_a_duplicate_in_the_left_sheet(clean_sheets):
    df_L, df_H = clean_sheets
    df_L = pd.concat([df_L, df_L.iloc[[0]]], ignore_index=True)  # inject "A" twice
    merged, both, _single = merge_and_split(df_L, df_H)

    with pytest.raises(AssertionError, match=r"duplicate .* in the LEFT sheet"):
        merge_guard.assert_merge_safe(df_L, df_H, merged, both)


def test_guard_trips_on_a_duplicate_in_the_right_sheet(clean_sheets):
    df_L, df_H = clean_sheets
    df_H = pd.concat([df_H, df_H.iloc[[1]]], ignore_index=True)  # inject "B" twice
    merged, both, _single = merge_and_split(df_L, df_H)

    with pytest.raises(AssertionError, match=r"duplicate .* in the RIGHT sheet"):
        merge_guard.assert_merge_safe(df_L, df_H, merged, both)


def test_the_error_message_names_the_offending_accession(clean_sheets):
    df_L, df_H = clean_sheets
    df_L = pd.concat([df_L, df_L.iloc[[2]]], ignore_index=True)  # "C" twice
    merged, both, _single = merge_and_split(df_L, df_H)

    with pytest.raises(AssertionError) as excinfo:
        merge_guard.assert_merge_safe(df_L, df_H, merged, both)

    message = str(excinfo.value)
    assert "'C' x2" in message
    assert "allow_duplicates=True" in message  # says what NOT to do, and why


def test_a_duplicate_on_both_sides_really_does_explode_the_merge(clean_sheets):
    """The failure the guard exists for, demonstrated rather than asserted.

    "A" twice on each side yields 2x2 = 4 matched rows for one protein. The
    matched-row count then exceeds ``min(len(L), len(H))``, which is the
    research1.md line-183 bound. Both checks are exercised: the duplicate check
    fires first, and with it waived the cardinality check still catches it.
    """
    df_L, df_H = clean_sheets
    df_L = pd.concat([df_L, df_L.iloc[[0]]], ignore_index=True)
    df_H = pd.concat([df_H] + [df_H.iloc[[0]]] * 3, ignore_index=True)
    merged, both, _single = merge_and_split(df_L, df_H)

    # One protein now occupies 2 x 4 = 8 rows.
    assert int((both[ACC] == "A").sum()) == 8
    assert len(both) > min(len(df_L), len(df_H))

    with pytest.raises(AssertionError, match="duplicate"):
        merge_guard.assert_merge_safe(df_L, df_H, merged, both)

    with pytest.raises(AssertionError, match="merge cardinality violated"):
        merge_guard.assert_merge_safe(df_L, df_H, merged, both, allow_duplicates=True)


def test_allow_duplicates_still_enforces_the_cardinality_bound(clean_sheets):
    """The escape hatch downgrades check 1 only; checks 2 and 3 always run."""
    df_L, df_H = clean_sheets
    df_L = pd.concat([df_L, df_L.iloc[[0]]], ignore_index=True)
    merged, both, _single = merge_and_split(df_L, df_H)

    # A duplicate on ONE side fans out but does not exceed min(len(L), len(H)),
    # so with the duplicate check waived this particular case is allowed through.
    stats = merge_guard.assert_merge_safe(df_L, df_H, merged, both, allow_duplicates=True)
    assert stats["n_dup_L"] == 1
    assert stats["n_dup_H"] == 0


def test_guard_returns_stats_for_logging(clean_sheets):
    df_L, df_H = clean_sheets
    merged, both, _single = merge_and_split(df_L, df_H)
    stats = merge_guard.assert_merge_safe(df_L, df_H, merged, both)

    assert stats == {
        "n_L": 4, "n_H": 4, "n_merged": 5, "n_both": 3,
        "n_dup_L": 0, "n_dup_H": 0, "max_both": 4,
    }


def test_guard_hardcodes_no_dataset_numbers(frozen_counts):
    """The bounds are derived from the inputs, so they survive a data change.

    A guard that knows 2315, 2187 or 1948 is a frozen expectation wearing a
    guard's clothes; those numbers belong in ``tests/expected/frozen_counts.json``.
    Numeric *literals* are what is checked -- the module's prose is free to quote
    today's counts, and does.
    """
    tree = ast.parse(Path(merge_guard.__file__).read_text(encoding="utf-8"))
    literals = {
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float))
        and not isinstance(node.value, bool)
    }
    # Small numbers (0, 1, 5) are ordinary programming, not dataset facts.
    dataset_numbers = {v for v in frozen_counts.values() if isinstance(v, int) and v > 100}

    assert dataset_numbers  # the fixture really did load
    assert literals & dataset_numbers == set(), (
        f"merge_guard.py hardcodes dataset-specific number(s): "
        f"{sorted(literals & dataset_numbers)}"
    )


# ---------------------------------------------------------------------------
# assert_classification_partition
# ---------------------------------------------------------------------------
def test_partition_holds_for_a_well_labelled_frame():
    df = pd.DataFrame({"regulated": ["UP", "DOWN", "NO CHANGE", "ON_OFF", "UP"]})
    counts = merge_guard.assert_classification_partition(df)

    assert counts == {"UP": 2, "DOWN": 1, "NO CHANGE": 1, "ON_OFF": 1}
    assert sum(counts.values()) == len(df)


def test_partition_trips_on_an_unknown_label():
    df = pd.DataFrame({"regulated": ["UP", "DOWN", "MAYBE"]})

    with pytest.raises(AssertionError) as excinfo:
        merge_guard.assert_classification_partition(df)

    message = str(excinfo.value)
    assert "does not partition" in message
    assert "'MAYBE'" in message


def test_partition_trips_on_a_null_label():
    df = pd.DataFrame({"regulated": ["UP", None, "DOWN"]})

    with pytest.raises(AssertionError, match=r"null labels\s+: 1"):
        merge_guard.assert_classification_partition(df)


def test_partition_is_checked_against_the_real_class_set():
    assert set(merge_guard.REGULATED_CLASSES) == set(core.REGULATED_CLASSES)


# ---------------------------------------------------------------------------
# The real sheets pass
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def real_sheets():
    if not _WORKBOOK.exists():  # pragma: no cover - the workbook is committed
        pytest.skip(f"input workbook not found: {_WORKBOOK}")
    cols_L = [ACC, "Gene names", "Intensity 31578", "Intensity 31580"]
    cols_H = [ACC, "Gene names", "Intensity 31579", "Intensity 31581"]
    return core.read_sheets(_WORKBOOK, cols_L, cols_H)


def test_the_committed_sheets_have_no_duplicate_accessions(real_sheets, frozen_counts):
    """L: 2315/2315 unique, H: 2187/2187. The guard is a no-op today, on purpose."""
    df_L, df_H = real_sheets

    assert len(df_L) == frozen_counts["raw_sheet_L_rows"]
    assert len(df_H) == frozen_counts["raw_sheet_H_rows"]
    assert merge_guard.duplicate_accession_count(df_L) == 0
    assert merge_guard.duplicate_accession_count(df_H) == 0
    assert df_L[ACC].nunique() == len(df_L)
    assert df_H[ACC].nunique() == len(df_H)


def test_the_real_merge_passes_the_guard(real_sheets, frozen_counts):
    df_L, df_H = real_sheets
    merged, both, single = merge_and_split(df_L, df_H)
    stats = merge_guard.assert_merge_safe(df_L, df_H, merged, both)

    assert stats["n_both"] == frozen_counts["foldchange_all_rows"]
    assert len(single) == frozen_counts["single_condition_rows"]
    assert stats["n_both"] <= stats["max_both"]
    assert stats["n_merged"] == stats["n_both"] + len(single)
