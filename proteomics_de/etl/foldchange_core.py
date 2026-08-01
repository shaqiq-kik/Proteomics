"""Pure, testable stages of the SILAC fold-change computation.

Everything in ``foldchange.py`` used to live inside its ``if __name__ ==
"__main__"`` block, which meant not one line of the fold-change logic could be
imported, exercised on a hand-built frame, or reasoned about in isolation. This
module holds those stages as ordinary functions over DataFrames; ``foldchange.py``
is now just the wiring that reads the workbook, calls these in order, and writes
the CSVs.

**Byte-identity contract.** Each function reproduces the corresponding block of
the pre-refactor script *statement for statement* -- same pandas calls, same
order, same dtypes. The refactor's acceptance test is that re-running the
pipeline leaves all 13 committed outputs sha256-identical, so "equivalent but
tidier" is not good enough here; the operations themselves are the specification.
Two places are subtle enough to call out:

* :func:`restore_left_order` undoes the two side effects of switching the merge
  from ``how="inner"`` to ``how="outer"`` (Bug 4): the outer join reorders the
  matched rows, and its NaNs coerce the Heavy intensity columns to float. Both
  the row order and the Heavy dtypes are restored from the input sheets.
* :func:`compute_ratios` creates *all* ratio columns as NaN before filling any of
  them, matching the original's two-phase write.

**Generality.** The functions take their column lists as arguments rather than
reading module-level literals, so they are exercisable on three-row fixtures and
readable without the sample sheet in hand. They are *not* a claim that the stage
is design-agnostic: the L/H sheet split and the left-order/dtype restore are
inherently 2-channel SILAC, which is why ``foldchange.py`` keeps its explicit
column lists and guards them with ``config.design.assert_matches``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

#: Merge key. Whole accession string, never a split token -- see
#: :mod:`proteomics_de.etl.accessions` for why a protein group is an identity.
ACCESSION_COL = "UniProt Accession Number"

#: Gene-symbol column, coalesced L-first across the merge.
GENE_COL = "Gene names"

#: pandas' ``indicator=True`` column name.
INDICATOR_COL = "_merge"

#: Class labels written to ``regulated``. The four are mutually exclusive and
#: exhaustive; :func:`proteomics_de.etl.merge_guard.assert_classification_partition`
#: enforces that.
REGULATED_CLASSES = ("UP", "DOWN", "NO CHANGE", "ON_OFF")

#: Dataset-specific expectations, lifted out of inline literals in the source.
FROZEN_COUNTS_PATH = (
    Path(__file__).resolve().parent.parent / "tests" / "expected" / "frozen_counts.json"
)

#: Environment variable that turns the frozen-count assertions off. They are ON
#: by default: they have caught real regressions, and the only thing wrong with
#: them was that the numbers were typed into the source.
BASELINE_ENV_VAR = "PDE_EXPECT_BASELINE"


# ---------------------------------------------------------------------------
# Frozen-count expectations
# ---------------------------------------------------------------------------
def baseline_checks_enabled(env=None) -> bool:
    """True unless ``PDE_EXPECT_BASELINE=0``.

    The dataset-specific row/class counts are asserted by default. A future
    dataset legitimately produces different numbers; that run sets
    ``PDE_EXPECT_BASELINE=0`` (and then re-derives ``frozen_counts.json``)
    instead of deleting the assertions.
    """
    env = os.environ if env is None else env
    return str(env.get(BASELINE_ENV_VAR, "1")).strip() != "0"


def load_frozen_counts(path=None) -> dict:
    """Read ``tests/expected/frozen_counts.json``.

    Resolved from ``__file__``, never from the cwd. Keys beginning with ``_``
    are commentary and are left in place; callers ask for the keys they need.
    """
    p = Path(path) if path is not None else FROZEN_COUNTS_PATH
    with p.open(encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# 1) Load
# ---------------------------------------------------------------------------
def read_sheets(input_file, cols_L, cols_H,
                sheet_L="Protein Report L", sheet_H="Protein Report H"):
    """Read the Light and Heavy protein report sheets, keeping only `cols_*`.

    Returns ``(df_L, df_H)``. Column selection happens here (not downstream) so
    the frames that cross the ``after_load`` boundary are already narrowed to the
    accession, the gene symbol, and that channel's intensity columns.
    """
    df_L = pd.read_excel(input_file, sheet_name=sheet_L)[list(cols_L)]
    df_H = pd.read_excel(input_file, sheet_name=sheet_H)[list(cols_H)]
    return df_L, df_H


# ---------------------------------------------------------------------------
# 2) Merge
# ---------------------------------------------------------------------------
def merge_with_indicator(df_L, df_H, key=ACCESSION_COL, gene_col=GENE_COL):
    """Outer-merge L and H on `key`, keeping pandas' ``_merge`` indicator.

    Bug 4: an inner join silently discards every protein seen in only one
    channel -- often the strongest biological signal. The outer join keeps them
    and the indicator says which side they came from. `gene_col` is coalesced
    L-first so a protein that has a symbol in either sheet keeps it.
    """
    merged = pd.merge(
        df_L, df_H, on=key, how="outer",
        indicator=True, suffixes=("_L", "_H"),
    )
    merged[gene_col] = merged[f"{gene_col}_L"].combine_first(merged[f"{gene_col}_H"])
    return merged


def split_both_single(merged, indicator=INDICATOR_COL):
    """Split the merged frame into ``(both, single_cond)``.

    ``both`` (present in each sheet) feeds the fold-change pipeline unchanged;
    ``single_cond`` (``left_only`` / ``right_only``) is rescued to its own file
    with no fold change. The two are disjoint and cover `merged` exactly.
    """
    both = merged[merged[indicator] == "both"].copy()
    single_cond = merged[merged[indicator] != "both"].copy()
    return both, single_cond


def restore_left_order(both, df_L, df_H, out_cols, dtype_cols, key=ACCESSION_COL):
    """Undo the outer join's two cosmetic side effects, then project `out_cols`.

    The fold-change pipeline operates only on `both`, so every number is
    unaffected by the inner→outer switch -- but the *file* would still change,
    because the outer join (a) reorders the matched rows and (b) lets the NaNs of
    the rescued single-condition rows coerce the Heavy intensity columns to
    float. This restores the original Light-sheet row order (a stable sort on the
    L-sheet position of each accession) and re-casts `dtype_cols` back to the
    dtype they carry in `df_H`.

    Do not "simplify" this into a reindex: the stable sort is what keeps
    duplicate-free ties in their original relative order, and the dtype restore
    must read the dtype from `df_H`, not assume ``int64``.
    """
    left_order = {acc: i for i, acc in enumerate(df_L[key])}
    both = both.sort_values(key, key=lambda s: s.map(left_order), kind="stable")
    df = both[list(out_cols)].reset_index(drop=True)
    for col in dtype_cols:
        df[col] = df[col].astype(df_H[col].dtype)
    return df


# ---------------------------------------------------------------------------
# 3) Completeness
# ---------------------------------------------------------------------------
def mark_complete(df, intensity_cols, column="complete"):
    """Flag rows whose every intensity is present and non-zero.

    Bug 3: ratios are computed only on complete rows, so a zero denominator can
    never reach the division and no ``inf`` can be produced. A row is INCOMPLETE
    if any intensity is NaN **or** 0. Mutates and returns `df`.
    """
    cols = list(intensity_cols)
    df[column] = ~(
        df[cols].isnull().any(axis=1) | (df[cols] == 0).any(axis=1)
    )
    return df


# ---------------------------------------------------------------------------
# 4-5) Ratios and log2 fold change
# ---------------------------------------------------------------------------
def compute_ratios(df, control_cols, treated_cols, mask, prefix="ratio_rep"):
    """Per-replicate SILAC ratios (treated / control), on `mask` rows only.

    Replicate *i* pairs ``treated_cols[i]`` over ``control_cols[i]``, producing
    ``ratio_rep1``, ``ratio_rep2``, ... Every column is created as NaN first and
    filled second, so incomplete rows carry NaN rather than a stale value and the
    guarded division never sees a zero denominator (Bug 3).

    Returns ``(df, ratio_cols)``.
    """
    control_cols = list(control_cols)
    treated_cols = list(treated_cols)
    if len(control_cols) != len(treated_cols):
        raise ValueError(
            f"control/treated replicate counts differ: {len(control_cols)} vs "
            f"{len(treated_cols)}. This stage assumes a balanced, paired design."
        )

    ratio_cols = [f"{prefix}{i}" for i in range(1, len(control_cols) + 1)]
    for col in ratio_cols:
        df[col] = np.nan
    for col, num, den in zip(ratio_cols, treated_cols, control_cols):
        df.loc[mask, col] = df.loc[mask, num] / df.loc[mask, den]
    return df, ratio_cols


def compute_log2fc(df, ratio_cols, ratio_prefix="ratio_rep", prefix="log2_rep",
                   out="log2FC"):
    """Bug 1 fix: log2 each ratio FIRST, then average the logs.

    Ratios are multiplicative, so their central tendency is the geometric mean --
    equivalently, the arithmetic mean of the log2 ratios. The legacy script
    averaged the raw ratios and logged once, which systematically inflates
    apparent up-regulation: a protein up 2x in one run and down 2x in the other
    gives ``mean(2.0, 0.5) = 1.25 -> log2 = +0.32`` ("up") where the honest answer
    is ``mean(log2 2.0, log2 0.5) = 0.0`` (no change).

    Returns ``(df, log_cols)``.
    """
    log_cols = []
    for rc in ratio_cols:
        lc = rc.replace(ratio_prefix, prefix, 1)
        df[lc] = np.log2(df[rc])
        log_cols.append(lc)
    df[out] = df[log_cols].mean(axis=1)
    return df, log_cols


def classify_regulated(df, mask, threshold, column="regulated", log2fc_col="log2FC"):
    """Bug 2 fix: symmetric UP/DOWN cutoffs in log2 space, on `mask` rows only.

    The legacy rule was ``ratio >= 1.5`` for UP but ``ratio <= 0.5`` for DOWN: a
    protein had to rise 1.5x to be called UP but fall a full 2x to be called
    DOWN, which biases the classifier against down-regulation. The symmetric
    partner of a 1.5x rise is a fall to ``1/1.5 = 0.667``, i.e. ``|log2FC| >=
    log2(1.5) = 0.585``. Mutates and returns `df`.
    """
    df[column] = "NO CHANGE"
    df.loc[mask & (df[log2fc_col] >= threshold), column] = "UP"
    df.loc[mask & (df[log2fc_col] <= -threshold), column] = "DOWN"
    return df


# ---------------------------------------------------------------------------
# 6) On/off ("cousin") proteins
# ---------------------------------------------------------------------------
def detect_onoff(df, control_cols, treated_cols,
                 column="onoff", regulated_col="regulated"):
    """Bug 4b: relabel one-sided proteins out of NO CHANGE into ON_OFF.

    A protein present in *both* sheets but whose entire control side (or entire
    treated side) is absent is a real on/off signal -- yet the Bug 3
    completeness rule parks it as incomplete, where it lands in NO CHANGE. That
    label is simply wrong.

    * control absent, treated present -> ``on_with_treatment``
    * treated absent, control present -> ``off_with_treatment``
    * absent on **both** sides -> fully empty, stays NO CHANGE (an empty row is
      not evidence of anything, and calling it ON_OFF would invent a signal)

    No log2FC is computed for on/off proteins: they are incomplete, so their
    log2FC is already NaN and stays that way. Mutates `df`; returns
    ``(df, control_off, treated_off)``.
    """
    control_cols = list(control_cols)
    treated_cols = list(treated_cols)
    control_absent = (df[control_cols].isnull() | (df[control_cols] == 0)).all(axis=1)
    treated_absent = (df[treated_cols].isnull() | (df[treated_cols] == 0)).all(axis=1)
    control_off = control_absent & ~treated_absent  # present in treated only -> "on"
    treated_off = treated_absent & ~control_absent  # present in control only -> "off"

    df[column] = ""
    df.loc[control_off, column] = "on_with_treatment"
    df.loc[treated_off, column] = "off_with_treatment"
    df.loc[control_off | treated_off, regulated_col] = "ON_OFF"
    return df, control_off, treated_off


# ---------------------------------------------------------------------------
# 7) Summary
# ---------------------------------------------------------------------------
def summarize(df, regulated_col="regulated", onoff_col="onoff",
              complete_col="complete") -> dict:
    """Class counts for the console summary and the sanity assertions.

    Returned as plain ``int`` so the values are JSON-serialisable and print
    identically to the ``numpy.int64`` they replace.
    """
    return {
        "n_rows": int(len(df)),
        "n_complete": int(df[complete_col].sum()),
        "n_up": int((df[regulated_col] == "UP").sum()),
        "n_down": int((df[regulated_col] == "DOWN").sum()),
        "n_nochange": int((df[regulated_col] == "NO CHANGE").sum()),
        "n_onoff": int((df[regulated_col] == "ON_OFF").sum()),
        "n_on": int((df[onoff_col] == "on_with_treatment").sum()),
        "n_off": int((df[onoff_col] == "off_with_treatment").sum()),
    }


# ---------------------------------------------------------------------------
# 8) Output frames
# ---------------------------------------------------------------------------
def build_ipa_frame(df, complete_col="complete", regulated_col="regulated"):
    """The IPA export selection: complete AND regulated (UP / DOWN / ON_OFF).

    Row selection only -- no reordering, no column projection. The caller owns
    the column list, which puts the identifier leftmost per QIAGEN's contract.
    """
    return df[df[complete_col] & (df[regulated_col] != "NO CHANGE")].copy()


def build_single_condition_frame(single_cond, key=ACCESSION_COL, gene_col=GENE_COL,
                                 indicator=INDICATOR_COL,
                                 left_condition="control", right_condition="treated"):
    """Bug 4 rescue frame: proteins seen in exactly one sheet.

    Adds ``accession`` / ``gene`` aliases and a ``detected_in`` label derived
    from the merge indicator. The absent condition's intensity columns stay
    blank, which is the point: the blank *is* the finding. Mutates and returns
    the frame it was given (the caller passes a copy).

    ``left_condition`` / ``right_condition`` name the conditions that the LEFT
    ("Protein Report L") and RIGHT ("Protein Report H") sheets hold. These are
    parameters rather than literals because sheet position and experimental
    condition are independent facts: which sheet a channel sits in is fixed by
    the workbook, but which condition it is comes from
    ``config/sample_sheet.tsv``. Hardcoding ``left_only -> "control_only"``
    silently mislabelled all 606 rescued proteins the moment DECISIONS_LOG D7
    corrected the assignment -- the label said ``control_only`` for proteins
    detected only in the *treated* channels.
    """
    single_cond["accession"] = single_cond[key]
    single_cond["gene"] = single_cond[gene_col]
    single_cond["detected_in"] = np.where(
        single_cond[indicator] == "left_only",
        f"{left_condition}_only",
        f"{right_condition}_only",
    )
    return single_cond


def build_onoff_frame(df, key=ACCESSION_COL, gene_col=GENE_COL,
                      regulated_col="regulated"):
    """Bug 4b file: the ON_OFF rows, with ``accession`` / ``gene`` aliases."""
    onoff_df = df[df[regulated_col] == "ON_OFF"].copy()
    onoff_df["accession"] = onoff_df[key]
    onoff_df["gene"] = onoff_df[gene_col]
    return onoff_df


__all__ = [
    "ACCESSION_COL",
    "GENE_COL",
    "INDICATOR_COL",
    "REGULATED_CLASSES",
    "FROZEN_COUNTS_PATH",
    "BASELINE_ENV_VAR",
    "baseline_checks_enabled",
    "load_frozen_counts",
    "read_sheets",
    "merge_with_indicator",
    "split_both_single",
    "restore_left_order",
    "mark_complete",
    "compute_ratios",
    "compute_log2fc",
    "classify_regulated",
    "detect_onoff",
    "summarize",
    "build_ipa_frame",
    "build_single_condition_frame",
    "build_onoff_frame",
]
