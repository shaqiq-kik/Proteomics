"""The experimental design, read from ``config/sample_sheet.tsv``.

This module is the single source of truth for the design. Today the same facts
are hardcoded as literals in nine files -- ``foldchange.py``, ``limma_test.py``,
``replicate_check.py``, ``centering_check.py``, ``qc/schema.py``,
``viz/style.py``, ``gated/pca_cluster.py`` and friends all carry their own copy
of ``["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity
31581"]`` and of ``["control", "control", "treated", "treated"]``. That is the
pipeline's central architectural flaw: the sample sheet drives nothing.

The sheet has four columns::

    sample    group    channel             replicate
    31578     control  Intensity 31578     1
    31580     control  Intensity 31580     2
    31579     treated  Intensity 31579     1
    31581     treated  Intensity 31581     2

``channel`` is the literal intensity column name as it appears in the input
workbook, so ``sample_columns()`` is just the ``channel`` column read out in the
canonical order.

**Canonical order** (relied on by limma's ``model.matrix(~group)`` and by every
downstream handoff): control group first, treated group second; within a group,
ascending ``replicate``. The order is derived from the ``group`` column, not from
file row order, so re-sorting the TSV cannot silently permute the design matrix.

**Forward path.** Adding biological replicates means appending rows to the sheet
(new sample ids, ``replicate`` 3, 4, ...). Every function here is
replicate-count- and group-count-agnostic; nothing needs editing. Stages that
still carry hardcoded literals should call :func:`assert_matches` so that they
fail loudly instead of silently analysing the wrong four columns.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

#: Path to the committed sample sheet. Resolved from ``__file__``, never from the
#: process's cwd -- the pipeline is invoked from several different directories.
DEFAULT_SAMPLE_SHEET = Path(__file__).resolve().parent / "sample_sheet.tsv"

REQUIRED_COLUMNS = ["sample", "group", "channel", "replicate"]

#: Group ordering. Groups named here sort in this order; any other group name
#: sorts after them, in order of first appearance in the sheet.
CANONICAL_GROUP_ORDER = ["control", "treated"]

#: Prefix used when a group's samples cross the limma handoff boundary
#: (``limma_test.py:56`` writes ``ctrl_31578`` / ``trt_31579``).
GROUP_HANDOFF_PREFIX = {"control": "ctrl", "treated": "trt"}


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def read_sample_sheet(path=None) -> pd.DataFrame:
    """Read the sample sheet and return it in :ref:`canonical order <order>`.

    Parameters
    ----------
    path :
        Override the sheet location. Defaults to :data:`DEFAULT_SAMPLE_SHEET`,
        which is resolved relative to this file -- never relative to the cwd.

    Returns
    -------
    pandas.DataFrame
        Columns ``sample``, ``group``, ``channel``, ``replicate``, sorted control
        group first then treated, and by ascending ``replicate`` within a group.
        ``sample`` and ``channel`` are strings; ``replicate`` is an int.
    """
    sheet_path = Path(path) if path is not None else DEFAULT_SAMPLE_SHEET
    if not sheet_path.exists():
        raise FileNotFoundError(f"sample sheet not found: {sheet_path}")

    df = pd.read_csv(sheet_path, sep="\t", dtype=str)

    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"sample sheet {sheet_path} is missing required column(s): {missing}. "
            f"Expected columns: {REQUIRED_COLUMNS}"
        )
    if df.empty:
        raise ValueError(f"sample sheet {sheet_path} has no rows.")

    for col in ("sample", "group", "channel"):
        df[col] = df[col].astype(str).str.strip()
    df["replicate"] = df["replicate"].astype(int)

    dupes = df["channel"][df["channel"].duplicated()].tolist()
    if dupes:
        raise ValueError(
            f"sample sheet {sheet_path} has duplicate channel(s): {dupes}. "
            "Each sample must map to exactly one intensity column."
        )

    return _canonical_sort(df).reset_index(drop=True)


def _canonical_sort(df: pd.DataFrame) -> pd.DataFrame:
    """Control group first, treated second, ascending ``replicate`` within group."""
    order = list(CANONICAL_GROUP_ORDER)
    for g in df["group"]:
        if g not in order:
            order.append(g)
    rank = {g: i for i, g in enumerate(order)}
    return df.assign(_grank=df["group"].map(rank)).sort_values(
        ["_grank", "replicate", "sample"], kind="mergesort"
    ).drop(columns="_grank")


def _sheet(sheet=None) -> pd.DataFrame:
    """Accept a preloaded sheet, a path, or nothing (load the default)."""
    if sheet is None:
        return read_sample_sheet()
    if isinstance(sheet, pd.DataFrame):
        return _canonical_sort(sheet).reset_index(drop=True)
    return read_sample_sheet(sheet)


# ---------------------------------------------------------------------------
# Derived views
# ---------------------------------------------------------------------------
def sample_columns(sheet=None) -> list[str]:
    """Intensity column names in canonical (control-then-treated) order.

    For today's sheet this is exactly::

        ["Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581"]
    """
    return _sheet(sheet)["channel"].tolist()


def group_names(sheet=None) -> list[str]:
    """Distinct group labels, in canonical order (``["control", "treated"]``)."""
    df = _sheet(sheet)
    seen: list[str] = []
    for g in df["group"]:
        if g not in seen:
            seen.append(g)
    return seen


def columns_for_group(group: str, sheet=None) -> list[str]:
    """Intensity column names belonging to `group`, ascending by replicate."""
    df = _sheet(sheet)
    cols = df.loc[df["group"] == group, "channel"].tolist()
    if not cols:
        raise KeyError(
            f"group {group!r} is not present in the sample sheet; "
            f"known groups: {group_names(df)}"
        )
    return cols


def control_columns(sheet=None) -> list[str]:
    """Intensity columns of the control (denominator) group."""
    return columns_for_group("control", sheet)


def treated_columns(sheet=None) -> list[str]:
    """Intensity columns of the treated (numerator) group."""
    return columns_for_group("treated", sheet)


def group_vector(sheet=None) -> list[str]:
    """Per-sample group labels aligned with :func:`sample_columns`.

    This is the vector handed to R as ``group <- c(...)`` for
    ``model.matrix(~ factor(group, levels = c("control", "treated")))``.
    For today's sheet: ``["control", "control", "treated", "treated"]``.
    """
    return _sheet(sheet)["group"].tolist()


def sample_ids(sheet=None) -> list[str]:
    """Sample identifiers aligned with :func:`sample_columns`."""
    return _sheet(sheet)["sample"].tolist()


def handoff_names(sheet=None) -> list[str]:
    """Column names used across the limma handoff, aligned with sample order.

    Pattern is ``<group prefix>_<sample id>``, matching ``limma_test.py:56``.
    For today's sheet: ``["ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581"]``.
    """
    df = _sheet(sheet)
    names = []
    for group, sample in zip(df["group"], df["sample"]):
        prefix = GROUP_HANDOFF_PREFIX.get(group, group)
        names.append(f"{prefix}_{sample}")
    return names


def n_samples(sheet=None) -> int:
    """Total number of samples (rows in the sheet)."""
    return len(_sheet(sheet))


def n_groups(sheet=None) -> int:
    """Number of distinct groups."""
    return len(group_names(sheet))


def replicates_per_group(sheet=None) -> int:
    """Replicate count per group.

    Raises
    ------
    ValueError
        If the design is unbalanced. limma tolerates unbalanced designs, but
        every stage in this pipeline assumes a balanced one, so an unbalanced
        sheet is reported rather than silently averaged away.
    """
    df = _sheet(sheet)
    counts = df.groupby("group", sort=False).size()
    unique = set(counts.tolist())
    if len(unique) != 1:
        raise ValueError(
            "unbalanced design: replicate counts per group are "
            f"{counts.to_dict()}. This pipeline assumes a balanced design; "
            "regenerate the affected stages before proceeding."
        )
    return unique.pop()


# ---------------------------------------------------------------------------
# Emitters / guards
# ---------------------------------------------------------------------------
def write_design_tsv(path, sheet=None) -> Path:
    """Write the ``design.tsv`` limma contract from research1.md Section 1.

    Format: TAB-separated, header ``sample<TAB>group``, one row per sample, in
    the same order as :func:`sample_columns`. Trailing newline, no index.
    """
    df = _sheet(sheet)
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df[["sample", "group"]].to_csv(out, sep="\t", index=False, lineterminator="\n")
    return out


def _caller_stage() -> str:
    """Best-effort name of the module that called :func:`assert_matches`."""
    try:
        frame = sys._getframe(2)
    except ValueError:  # pragma: no cover - only if called at stack bottom
        return "<unknown stage>"
    return Path(frame.f_globals.get("__file__", "<unknown stage>")).name


def assert_matches(expected_columns, sheet=None, stage=None) -> None:
    """Guard a stage's hardcoded column list against the sample sheet.

    Tier-4 stages (``foldchange.py``, ``replicate_check.py``, ``qc/schema.py``)
    still carry their own literal column lists. They call this so that a change
    to the sample sheet cannot silently leave them analysing the wrong columns:
    the run fails immediately with a message that says what to do.

    Parameters
    ----------
    expected_columns :
        The stage's hardcoded list, in the stage's own order.
    sheet :
        Optional sheet / path override, forwarded to :func:`read_sample_sheet`.
    stage :
        Human-readable stage name for the error message. Defaults to the calling
        module's filename.

    Raises
    ------
    ValueError
        If `expected_columns` differs from :func:`sample_columns` in content or
        in order.
    """
    stage_name = stage or _caller_stage()
    actual = sample_columns(sheet)
    expected = list(expected_columns)
    if expected == actual:
        return

    missing = [c for c in actual if c not in expected]
    unexpected = [c for c in expected if c not in actual]
    if missing or unexpected:
        detail = (
            f"  columns in the sample sheet but not in {stage_name}: {missing}\n"
            f"  columns in {stage_name} but not in the sample sheet: {unexpected}\n"
        )
    else:
        detail = "  same columns, different order\n"

    raise ValueError(
        f"design drift in {stage_name}: its hardcoded intensity columns no "
        f"longer match config/sample_sheet.tsv.\n"
        f"  sample sheet : {actual}\n"
        f"  {stage_name} : {expected}\n"
        f"{detail}"
        f"{stage_name} is 2-channel-SILAC-specific: it assumes exactly two "
        f"groups of two technical replicates, in control-then-treated order, "
        f"and it must be REGENERATED for a new design -- do not patch the "
        f"literal list to silence this error. See "
        f"proteomics_de/config/design.py and research1.md's FORWARD-PATH "
        f"section."
    )


__all__ = [
    "DEFAULT_SAMPLE_SHEET",
    "REQUIRED_COLUMNS",
    "CANONICAL_GROUP_ORDER",
    "GROUP_HANDOFF_PREFIX",
    "read_sample_sheet",
    "sample_columns",
    "group_names",
    "columns_for_group",
    "control_columns",
    "treated_columns",
    "group_vector",
    "sample_ids",
    "handoff_names",
    "n_samples",
    "n_groups",
    "replicates_per_group",
    "write_design_tsv",
    "assert_matches",
]
