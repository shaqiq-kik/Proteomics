"""Merge cardinality and classification-partition guards.

research1.md asks for these twice and the pipeline had neither:

* line 52 (Bug 4) -- *"Additional join hazards on ``UniProt Accession Number``:
  ... duplicate accessions causing many-to-many row explosion. **Fix:** validate
  accession cardinality before the merge"*.
* line 183 (Section 2, Duplicate / cardinality detection) -- *"count duplicate
  accessions pre-merge; assert merge does not increase row count beyond
  ``min(len(L), len(H))``"*.

Why it matters. ``pd.merge`` on a key that repeats in **both** frames is a
cross product: an accession appearing twice on each side yields four rows, and
the extra rows are indistinguishable from real proteins downstream. Nothing in
the pipeline would notice -- the row count would simply be wrong, every count
derived from it would be wrong, and the fold changes would be silently
duplicated into the statistics. Today's sheets are clean (L: 2315/2315 unique,
H: 2187/2187), so this guard is currently a no-op; that is exactly when it is
worth installing, because the failure mode it catches is invisible.

**Dataset-agnostic by construction.** Nothing here knows 206, 509, 715 or 1948.
Every bound is derived from the frames it is handed, so the guards keep working
when the input data changes -- which is the difference between a guard and a
frozen expectation. The dataset-specific counts live in
``tests/expected/frozen_counts.json`` and are asserted separately.

The two entry points are :func:`assert_merge_safe` (called right after the
merge) and :func:`assert_classification_partition` (called once ``regulated`` is
populated). Both raise ``AssertionError`` with a message that names the
offending values, because a guard whose failure you cannot diagnose is only
marginally better than no guard.
"""

from __future__ import annotations

import pandas as pd

#: Default merge key. Matches ``foldchange_core.ACCESSION_COL``.
ACCESSION_COL = "UniProt Accession Number"

#: The exhaustive, mutually exclusive set of ``regulated`` labels.
REGULATED_CLASSES = ("UP", "DOWN", "NO CHANGE", "ON_OFF")

#: How many offending values to name in an error message.
_MAX_EXAMPLES = 5


def duplicate_accessions(df, key=ACCESSION_COL) -> pd.Series:
    """Accessions occurring more than once in `df`, as ``{accession: count}``.

    Missing keys are excluded: a NaN accession is a separate data-quality
    problem (the schema layer reports it) and ``pd.merge`` does not join on it,
    so it cannot cause the row explosion this module guards against.
    """
    counts = df[key].dropna().value_counts()
    return counts[counts > 1]


def duplicate_accession_count(df, key=ACCESSION_COL) -> int:
    """Number of *distinct* accessions that occur more than once in `df`."""
    return int(len(duplicate_accessions(df, key)))


def _describe_duplicates(dupes: pd.Series) -> str:
    head = dupes.head(_MAX_EXAMPLES)
    shown = ", ".join(f"{acc!r} x{int(n)}" for acc, n in head.items())
    if len(dupes) > _MAX_EXAMPLES:
        shown += f", ... ({len(dupes) - _MAX_EXAMPLES} more)"
    return shown


def assert_merge_safe(df_L, df_H, merged, both, key=ACCESSION_COL,
                      allow_duplicates=False) -> dict:
    """Guard the L/H merge against a many-to-many row explosion.

    Three checks, in the order a failure is most likely to be diagnosable:

    1. **Duplicate accessions per input sheet.** Counted for both sheets and
       returned in the stats dict. A duplicate on *one* side only fans rows out;
       a duplicate on *both* sides multiplies them. Either way the merge is no
       longer a protein-to-protein correspondence, so both are rejected unless
       `allow_duplicates` is set (which the caller must then justify in writing).
    2. **Matched-row bound.** ``len(both) <= min(len(df_L), len(df_H))`` -- the
       matched group cannot be larger than the smaller sheet, because each
       matched row consumes at least one row from each side. This is the
       research1.md line-183 assertion, and it is the one that trips on an
       explosion even if the duplicates check were somehow bypassed.
    3. **Total-row bound.** ``len(merged) <= len(df_L) + len(df_H)`` -- an outer
       join over unique keys emits at most one row per input row.

    Parameters
    ----------
    df_L, df_H :
        The two input sheets, *before* the merge.
    merged :
        The full outer-merge result.
    both :
        The ``_merge == "both"`` subset of `merged`.
    key :
        Merge key column name.
    allow_duplicates :
        Downgrade check 1 to a reported statistic. Checks 2 and 3 still run.

    Returns
    -------
    dict
        ``n_L``, ``n_H``, ``n_merged``, ``n_both``, ``n_dup_L``, ``n_dup_H``,
        ``max_both``. Useful for logging and for tests that want the numbers
        rather than just the absence of an exception.

    Raises
    ------
    AssertionError
        On any failed check, naming the offending accessions or counts.
    """
    dup_L = duplicate_accessions(df_L, key)
    dup_H = duplicate_accessions(df_H, key)

    stats = {
        "n_L": int(len(df_L)),
        "n_H": int(len(df_H)),
        "n_merged": int(len(merged)),
        "n_both": int(len(both)),
        "n_dup_L": int(len(dup_L)),
        "n_dup_H": int(len(dup_H)),
        "max_both": int(min(len(df_L), len(df_H))),
    }

    if not allow_duplicates:
        if len(dup_L):
            raise AssertionError(
                f"duplicate {key!r} in the LEFT sheet: {len(dup_L)} accession(s) "
                f"occur more than once ({_describe_duplicates(dup_L)}). "
                "pd.merge would fan these out into extra rows that are "
                "indistinguishable from real proteins downstream. De-duplicate "
                "the input, or decide explicitly how repeated accessions should "
                "be collapsed -- do not pass allow_duplicates=True to silence it."
            )
        if len(dup_H):
            raise AssertionError(
                f"duplicate {key!r} in the RIGHT sheet: {len(dup_H)} accession(s) "
                f"occur more than once ({_describe_duplicates(dup_H)}). "
                "pd.merge would fan these out into extra rows that are "
                "indistinguishable from real proteins downstream. De-duplicate "
                "the input, or decide explicitly how repeated accessions should "
                "be collapsed -- do not pass allow_duplicates=True to silence it."
            )

    if len(both) > stats["max_both"]:
        raise AssertionError(
            f"merge cardinality violated: {len(both)} matched rows exceed "
            f"min(len(L), len(H)) = {stats['max_both']} "
            f"(L={stats['n_L']}, H={stats['n_H']}). The join produced more "
            "matches than either sheet can supply, which means a many-to-many "
            f"explosion on {key!r} (research1.md lines 52, 183)."
        )

    if len(merged) > stats["n_L"] + stats["n_H"]:
        raise AssertionError(
            f"merge cardinality violated: the outer join emitted {len(merged)} "
            f"rows from {stats['n_L']} + {stats['n_H']} = "
            f"{stats['n_L'] + stats['n_H']} input rows. An outer join over "
            f"unique keys emits at most one row per input row, so {key!r} "
            "repeats on both sides."
        )

    return stats


def assert_classification_partition(df, column="regulated",
                                    classes=REGULATED_CLASSES) -> dict:
    """Assert that `column` partitions `df`: the class counts sum to ``len(df)``.

    ``n_up + n_down + n_nochange + n_onoff == len(df)``. This catches two things
    the individual count assertions cannot: a row that acquired a label outside
    the known set (a typo, or a new class added upstream without updating the
    consumers), and a row left unlabelled. Both would otherwise show up only as a
    downstream file being mysteriously short.

    Returns the per-class counts. Raises ``AssertionError`` naming the unexpected
    labels when the partition does not hold.
    """
    counts = {cls: int((df[column] == cls).sum()) for cls in classes}
    total = sum(counts.values())
    if total != len(df):
        seen = set(df[column].dropna().unique())
        unexpected = sorted(seen - set(classes))
        n_missing_label = int(df[column].isnull().sum())
        raise AssertionError(
            f"{column!r} does not partition the frame: class counts sum to "
            f"{total} but the frame has {len(df)} rows.\n"
            f"  expected classes : {list(classes)}\n"
            f"  counts           : {counts}\n"
            f"  unexpected labels: {unexpected}\n"
            f"  null labels      : {n_missing_label}\n"
            "Every row must carry exactly one known class; a row that does not "
            "will be dropped or double-counted by the downstream exports."
        )
    return counts


__all__ = [
    "ACCESSION_COL",
    "REGULATED_CLASSES",
    "duplicate_accessions",
    "duplicate_accession_count",
    "assert_merge_safe",
    "assert_classification_partition",
]
