"""Stage-boundary validation hooks (research1.md Section 6, line 178).

WHAT CHANGED
------------
This was a Wave-0 no-op shim: ``check(stage, df, schema=None)`` returned ``df``
and did nothing. The call sites in ``foldchange.py`` (``after_load`` x2,
``after_merge``, ``after_foldchange``) were installed against that contract and
are unchanged by filling the body in -- :func:`check` still returns *the same
object* it was handed, still copies nothing, and still mutates nothing.

WHY TWO POLICIES
----------------
A single "validate and raise" rule would be wrong here, and the reason is
specific rather than a matter of taste:

* ``after_load`` sees the RAW MaxQuant sheets. Two rows of the Light sheet carry
  a ';'-joined list of bare MaxQuant row indices where an accession belongs
  (32,759 and 681 characters). That is real, it is upstream of anything this
  repo controls, and it is *expected* -- so an accession-regex failure there is
  routed to QUARANTINE (:func:`quarantine_junk_accessions`), not raised.
  Raising would mean the pipeline could not read its own input.
* ``after_merge`` is the same frame plus a merge indicator; the junk rows are
  still in it (they are ``left_only``), so it is permissive for the same reason.
* ``after_foldchange`` and ``before_ipa_export`` are STRICT. By then the junk
  rows are gone (they never had a fold change to compute) and every value has
  been produced by code in this repo. A failure there is a bug in the pipeline,
  not a fact about the input, and it must stop the run.

The distinction is the whole point: permissive-everywhere is a validator that
cannot fail, strict-everywhere is a validator that cannot run.

THE RECORD
----------
Every call appends a structured record to ``results/qc/qc_boundaries.json``.
The file is rewritten at the first :func:`check` of each process, so it always
describes the most recent run rather than growing without bound.

It carries **no timestamp on purpose**. It lives under ``results/`` and is
therefore covered by the byte-freeze gate (``tools/freeze.py``); a wall-clock
field would make the file drift on every run for no informational gain, since
git already records when it changed. ``qc_report.json`` does stamp a time and
does drift -- that is the behaviour not being copied here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pandera.pandas as pa

_HERE = Path(__file__).resolve().parent          # proteomics_de/qc
_PKG_DIR = _HERE.parent                          # proteomics_de
_REPO_ROOT = _PKG_DIR.parent                     # repo root

# Same dual-import accommodation as schema.py: this module is reachable both as
# `from proteomics_de.qc import boundaries` and, in the flat-script world, with
# proteomics_de/qc on sys.path.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.qc import schema as qc_schema  # noqa: E402

#: The boundaries a stage may validate at.
#:
#: * ``after_load``        -- raw sheets read from the input workbook
#: * ``after_merge``       -- L/H outer merge, before fold changes exist
#: * ``after_foldchange``  -- log2FC / regulated / onoff columns populated
#: * ``before_ipa_export`` -- the filtered frame handed to the IPA writer
STAGES = ("after_load", "after_merge", "after_foldchange", "before_ipa_export")

#: Stages where a schema failure is a bug and must stop the run.
STRICT_STAGES = frozenset({"after_foldchange", "before_ipa_export"})

#: Stages that see raw upstream data, where a malformed accession is routed to
#: quarantine rather than raised.
PERMISSIVE_STAGES = frozenset({"after_load", "after_merge"})

DEFAULT_RESULTS_DIR = _PKG_DIR / "results"
BOUNDARY_RECORD_NAME = "qc_boundaries.json"
QUARANTINE_NAME = "quarantine_accessions.csv"

MAX_FAILURE_CASES_RECORDED = 25

#: Longest failing VALUE reproduced verbatim in a boundary record. The junk
#: accessions are 32,759 and 681 characters, and embedding them whole made this
#: committed JSON 70 KB of mostly one repeated string. The full value is not
#: lost -- it is recorded in results/qc/quarantine_accessions.csv, which is what
#: `quarantine_file` on the record points at.
MAX_FAILURE_VALUE_CHARS = 120

QUARANTINE_COLUMNS = [
    "source",
    "stage",
    "row_index",
    "column",
    "reason",
    "n_chars",
    "n_tokens",
    "value_preview",
    "value",
]

JUNK_INDEX_LIST_REASON = (
    "accession is a MaxQuant row-index list (every ';'-token is a bare "
    "integer), not a UniProt accession or protein group -- DECISIONS_LOG D11"
)


class BoundaryValidationError(AssertionError):
    """A strict boundary rejected the frame crossing it.

    Subclasses AssertionError so it reads like the structural guards in
    ``etl/merge_guard.py`` that it sits alongside, and so a bare
    ``except AssertionError`` in a caller still catches it.
    """


# ---------------------------------------------------------------------------
# Stage schemas
#
# These describe frames IN FLIGHT, so they differ from the file contracts in
# schema.py in two deliberate ways:
#   * strict=False -- an in-flight frame legitimately carries scratch columns
#     (`_merge`, `Gene names_L/_H`, the ratio intermediates) that never reach a
#     file. Missing or renamed DECLARED columns are still failures.
#   * no row-count checks -- 1948/606/715 are frozen facts about one dataset
#     and live in tests/expected/frozen_counts.json. A stage schema must hold
#     for any dataset.
# ---------------------------------------------------------------------------

# Both raw sheets share an accession and a gene column but carry DIFFERENT
# intensity columns (L: 31578/31580, H: 31579/31581), and foldchange.py calls
# after_load once for each. A regex column validates whichever pair is present
# instead of forcing two near-identical schemas.
_INTENSITY_COLUMN_PATTERN = r"^Intensity \d+$"

RAW_SHEET_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": qc_schema.accession_column(),
        "Gene names": qc_schema.gene_column(),
        _INTENSITY_COLUMN_PATTERN: pa.Column(
            float, pa.Check.ge(0), nullable=True, coerce=True, regex=True
        ),
    },
    strict=False,
    coerce=True,
    name="RAW_SHEET_SCHEMA",
)

MERGED_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": qc_schema.accession_column(),
        "Gene names": qc_schema.gene_column(),
        # `_merge` arrives as a pandas Categorical, so no dtype is declared --
        # only the vocabulary, which is what actually matters downstream
        # (split_both_single keys off exactly these three values).
        "_merge": pa.Column(
            checks=pa.Check.isin(["left_only", "right_only", "both"]),
        ),
        _INTENSITY_COLUMN_PATTERN: pa.Column(
            float, pa.Check.ge(0), nullable=True, coerce=True, regex=True
        ),
    },
    strict=False,
    coerce=True,
    name="MERGED_SCHEMA",
)

FOLDCHANGE_STAGE_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": qc_schema.accession_column(),
        "Gene names": qc_schema.gene_column(),
        _INTENSITY_COLUMN_PATTERN: pa.Column(
            float, pa.Check.ge(0), nullable=True, coerce=True, regex=True
        ),
        "ratio_rep1": pa.Column(float, nullable=True, coerce=True),
        "ratio_rep2": pa.Column(float, nullable=True, coerce=True),
        "log2_rep1": pa.Column(float, nullable=True, coerce=True),
        "log2_rep2": pa.Column(float, nullable=True, coerce=True),
        "log2FC": pa.Column(float, nullable=True, coerce=True),
        "complete": pa.Column(bool, coerce=True),
        "regulated": pa.Column(
            str, pa.Check.isin(qc_schema.REGULATED_CATEGORIES), coerce=True
        ),
        # IN FLIGHT the "not an on/off protein" sentinel is the EMPTY STRING,
        # not NaN: etl/foldchange_core.detect_onoff seeds the column with
        # `df[column] = ""` and only fills the two real labels. It becomes an
        # empty CSV field on write, which pandas reads back as NaN -- which is
        # why the file schema (schema.py) declares this column nullable and
        # this one does not. Declaring the sentinel is describing the contract,
        # not relaxing it; the dataframe-level check below then pins the
        # sentinel to exactly the non-ON_OFF rows, which the file schema cannot
        # express as tightly.
        "onoff": pa.Column(
            str, pa.Check.isin([*qc_schema.ONOFF_CATEGORIES, ""]), coerce=True
        ),
    },
    checks=[
        pa.Check(
            qc_schema._log2fc_finite_when_complete,
            error="log2FC must be finite for complete==True rows, NaN otherwise",
        ),
        pa.Check(
            lambda df: (df["onoff"].astype(str) != "").eq(df["regulated"] == "ON_OFF").all(),
            error=(
                "the onoff label must be set on exactly the ON_OFF rows: a "
                "labelled row that is not ON_OFF, or an ON_OFF row with no "
                "label, means detect_onoff and classify_regulated disagree"
            ),
        ),
    ],
    strict=False,
    coerce=True,
    name="FOLDCHANGE_STAGE_SCHEMA",
)

IPA_EXPORT_STAGE_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": qc_schema.accession_column(),
        "Gene names": qc_schema.gene_column(),
        # Every exported row is regulated with a real ratio behind it, so
        # unlike the stage above, log2FC here is non-nullable and must be
        # finite. A NaN reaching the IPA upload is silent data loss.
        "log2FC": pa.Column(float, qc_schema._isfinite_check(), coerce=True),
        "regulated": pa.Column(str, pa.Check.isin(["UP", "DOWN"]), coerce=True),
    },
    strict=False,
    coerce=True,
    name="IPA_EXPORT_STAGE_SCHEMA",
)

STAGE_SCHEMAS: dict[str, pa.DataFrameSchema] = {
    "after_load": RAW_SHEET_SCHEMA,
    "after_merge": MERGED_SCHEMA,
    "after_foldchange": FOLDCHANGE_STAGE_SCHEMA,
    "before_ipa_export": IPA_EXPORT_STAGE_SCHEMA,
}


# ---------------------------------------------------------------------------
# The run record
# ---------------------------------------------------------------------------

_records: list[dict] = []
_results_dir_override: Path | None = None


def set_results_dir(path) -> None:
    """Point the boundary record and quarantine files at `path` for this process.

    The four ``foldchange.py`` call sites pass only ``(stage, df)`` -- that is
    the Wave-0 contract and it is deliberately not being changed. But
    ``foldchange.py`` accepts ``--results-dir``, and without this the boundary
    record would be written to the PACKAGE's ``results/qc/`` even when the run's
    outputs are going somewhere else. That is not a cosmetic inconsistency: it
    means a run into a temp directory silently writes into the committed tree.
    ``main`` calls this once, so the whole stage layer follows ``--results-dir``
    while every call site stays a two-argument call.

    Pass ``None`` to restore the package default.
    """
    global _results_dir_override
    _results_dir_override = Path(path) if path is not None else None


def reset() -> None:
    """Drop the in-memory records and any results-dir override."""
    _records.clear()
    set_results_dir(None)


def records() -> list[dict]:
    """The records accumulated so far in this process (a copy)."""
    return [dict(r) for r in _records]


def _resolve_results_dir(results_dir=None) -> Path:
    if results_dir is not None:
        return Path(results_dir)
    if _results_dir_override is not None:
        return _results_dir_override
    return DEFAULT_RESULTS_DIR


def _record_path(results_dir=None) -> Path:
    return _resolve_results_dir(results_dir) / "qc" / BOUNDARY_RECORD_NAME


def quarantine_path(results_dir=None) -> Path:
    return _resolve_results_dir(results_dir) / "qc" / QUARANTINE_NAME


def _write_records(results_dir=None) -> Path:
    path = _record_path(results_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "note": (
            "Stage-boundary validation records for the most recent pipeline run "
            "(proteomics_de/qc/boundaries.py). Permissive stages (after_load, "
            "after_merge) see raw MaxQuant data and route malformed accessions to "
            "results/qc/quarantine_accessions.csv; strict stages "
            "(after_foldchange, before_ipa_export) raise. Deliberately carries no "
            "timestamp so the file is byte-reproducible under the freeze gate."
        ),
        "stages": list(STAGES),
        "strict_stages": sorted(STRICT_STAGES),
        "permissive_stages": sorted(PERMISSIVE_STAGES),
        "records": _records,
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


def _truncate(value):
    """Shorten an over-long failing value, marking that it was shortened."""
    if isinstance(value, str) and len(value) > MAX_FAILURE_VALUE_CHARS:
        return (
            f"{value[:MAX_FAILURE_VALUE_CHARS]}... "
            f"[truncated, {len(value)} chars total; full value in "
            f"qc/{QUARANTINE_NAME}]"
        )
    return value


def _failure_records(failure_cases: pd.DataFrame) -> list[dict]:
    head = failure_cases.head(MAX_FAILURE_CASES_RECORDED)
    records_ = json.loads(head.to_json(orient="records", date_format="iso"))
    return [{k: _truncate(v) for k, v in row.items()} for row in records_]


def _accession_columns(df: pd.DataFrame) -> list[str]:
    """Accession-bearing columns present in `df`, in a fixed order."""
    return [c for c in ("UniProt Accession Number", "accession", "id") if c in df.columns]


def _count_quarantinable(df: pd.DataFrame) -> int:
    total = 0
    for col in _accession_columns(df):
        total += int(df[col].map(qc_schema.is_quarantinable).sum())
    return total


# ---------------------------------------------------------------------------
# Quarantine
# ---------------------------------------------------------------------------


def find_junk_accessions(df: pd.DataFrame, column: str) -> pd.Series:
    """Boolean mask of rows whose `column` value is a MaxQuant row-index list.

    Discrimination is by token shape (all tokens are bare integers), never by
    length -- delegated to ``etl/accessions.is_junk_index_list`` so the
    69-character legitimate protein groups
    (``P08752;P20612;Q9DC51;...``) and the 681-character junk list cannot be
    confused. See DECISIONS_LOG D11.
    """
    return df[column].map(qc_schema.is_quarantinable).astype(bool)


def quarantine_junk_accessions(
    df: pd.DataFrame,
    column: str,
    *,
    source: str,
    stage: str = "after_load",
    results_dir=None,
    write: bool = True,
):
    """Split `df` into (kept, quarantined) on junk accessions in `column`.

    The quarantined rows are recorded IN FULL in
    ``results/qc/quarantine_accessions.csv`` -- quarantine means set aside with
    a reason, not deleted, so the original value stays recoverable from the
    record. Row order of the kept frame is untouched.

    Returns ``(kept, quarantined)``. `write` is provided so tests can exercise
    the split without touching the committed artifact.
    """
    mask = find_junk_accessions(df, column)
    quarantined = df.loc[mask]
    kept = df.loc[~mask]

    if write:
        write_quarantine(
            quarantined, column, source=source, stage=stage, results_dir=results_dir
        )
    return kept, quarantined


def build_quarantine_frame(
    quarantined: pd.DataFrame, column: str, *, source: str, stage: str
) -> pd.DataFrame:
    """The audit rows for `quarantined`, one per quarantined value."""
    rows = []
    for position, (_, row) in enumerate(quarantined.iterrows()):
        value = "" if pd.isna(row[column]) else str(row[column])
        tokens = value.split(";") if value else []
        rows.append(
            {
                "source": source,
                "stage": stage,
                "row_index": position,
                "column": column,
                "reason": JUNK_INDEX_LIST_REASON,
                "n_chars": len(value),
                "n_tokens": len(tokens),
                "value_preview": (value[:60] + "...") if len(value) > 63 else value,
                "value": value,
            }
        )
    return pd.DataFrame(rows, columns=QUARANTINE_COLUMNS)


def write_quarantine(
    quarantined: pd.DataFrame, column: str, *, source: str, stage: str, results_dir=None
) -> Path:
    """Write the quarantine record for `quarantined` (header-only if empty)."""
    path = quarantine_path(results_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = build_quarantine_frame(quarantined, column, source=source, stage=stage)
    frame.to_csv(path, index=False, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The hook itself
# ---------------------------------------------------------------------------


def check(stage, df, schema=None, results_dir=None):
    """Validate `df` at the named `stage` and return it.

    Returns THE SAME OBJECT it was passed, never a copy: the Wave-0 shim
    promised that and the ``foldchange.py`` call sites were written against it.
    Coercion happens inside pandera on its own copy and is discarded; nothing
    here mutates the caller's frame.

    Parameters
    ----------
    stage :
        One of :data:`STAGES`.
    df :
        The frame crossing the boundary.
    schema :
        Optional explicit schema override; defaults to :data:`STAGE_SCHEMAS`.
    results_dir :
        Where ``qc/qc_boundaries.json`` is written; defaults to the package's
        own ``results/``.

    Raises
    ------
    BoundaryValidationError
        If `stage` is strict (:data:`STRICT_STAGES`) and validation failed.
    ValueError
        If `stage` is not a known boundary -- a typo'd stage name would
        otherwise silently validate nothing at all.
    """
    if stage not in STAGES:
        raise ValueError(f"unknown boundary stage {stage!r}; expected one of {STAGES}")

    active = schema if schema is not None else STAGE_SCHEMAS[stage]
    strict = stage in STRICT_STAGES

    record = {
        "seq": len(_records),
        "stage": stage,
        "policy": "strict" if strict else "permissive",
        "schema": getattr(active, "name", None) or type(active).__name__,
        "rows": int(len(df)),
        "columns": [str(c) for c in df.columns],
        "passed": True,
        "n_failures": 0,
        "n_quarantinable_rows": _count_quarantinable(df),
        "failure_cases": [],
        "raised": False,
    }

    error = None
    try:
        active.validate(df, lazy=True)
    except pa.errors.SchemaErrors as exc:
        failure_cases = exc.failure_cases
        record["passed"] = False
        record["n_failures"] = int(len(failure_cases))
        record["failure_cases"] = _failure_records(failure_cases)
        error = exc

    if error is not None and not strict:
        # Permissive: the failures are recorded and, for accessions, handled by
        # quarantine_junk_accessions at the point the frame is written out.
        # The path is recorded RELATIVE to the results dir, not absolutely: this
        # file is committed, and an absolute path would bake one machine's
        # directory layout (and one --results-dir choice) into it.
        record["routed_to_quarantine"] = record["n_quarantinable_rows"] > 0
        record["quarantine_file"] = f"qc/{QUARANTINE_NAME}"

    if error is not None and strict:
        record["raised"] = True

    _records.append(record)
    _write_records(results_dir)

    if error is not None and strict:
        raise BoundaryValidationError(
            f"boundary {stage!r} rejected the frame: {record['n_failures']} "
            f"schema failure(s) against {record['schema']}. See "
            f"{_record_path(results_dir)} for the recorded failure cases.\n{error}"
        ) from error

    return df


__all__ = [
    "STAGES",
    "STRICT_STAGES",
    "PERMISSIVE_STAGES",
    "STAGE_SCHEMAS",
    "RAW_SHEET_SCHEMA",
    "MERGED_SCHEMA",
    "FOLDCHANGE_STAGE_SCHEMA",
    "IPA_EXPORT_STAGE_SCHEMA",
    "BoundaryValidationError",
    "QUARANTINE_COLUMNS",
    "check",
    "reset",
    "set_results_dir",
    "records",
    "find_junk_accessions",
    "quarantine_junk_accessions",
    "build_quarantine_frame",
    "write_quarantine",
    "quarantine_path",
]
