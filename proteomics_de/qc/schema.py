"""
proteomics_de/qc/schema.py

Pandera schemas encoding the REAL, observed contracts of the five committed,
frozen pipeline output CSVs (research1.md Section 6, item 3; DECISIONS_LOG.md
D1). This module is READ-ONLY documentation-as-code: it never writes to the
CSVs it describes. Run `python proteomics_de/qc/validate.py` to check the
committed files against these schemas.

Every schema here was built by first inspecting the actual committed data
(not just the task's abstract spec), because a schema that fails on
legitimate data is wrong. Three adaptations were necessary and are documented
below and at their point of use:

1. UniProt accession regex, extended for O/P/Q-prefixed accessions.
   The literal regex `^[A-NR-Z0-9]([A-Z0-9]{5}|[A-Z0-9]{9})(-\\d+)?$` is
   exactly the NON-[OPQ] branch of the two-branch official UniProt accession
   format (its first-character class `[A-NR-Z0-9]` deliberately excludes the
   letters O, P, Q). But O/P/Q-prefixed accessions are common and legitimate
   in this real dataset -- P19137 (Lama1), Q9JHU4 (Dync1h1), O08528 (Hk2),
   etc. -- and would ALL be rejected by the literal regex. UNIPROT_TOKEN_RE
   below adds the missing `[OPQ][A-Z0-9]{5}` branch (OPQ-prefixed UniProt
   accessions are always 6 characters, never 10) so the schema validates the
   real accession format rather than a subset of it.

2. Multi-protein-group entries. MaxQuant-style output frequently cannot
   assign a peptide to one protein unambiguously and instead reports several
   UniProt accessions joined by ';' (e.g. 'P05132;P68181', or
   'O70555;O70562;O70558;O70554;O70559'). foldchange_all.csv has 48 such
   rows, qc_limma.csv 48, ipa_input.csv 21, single_condition_proteins.csv 67.
   Accession validation therefore splits on ';' and validates each token.

3. NaN gene names. Some UniProt entries in this dataset (often the
   multi-protein-group ones) have no resolved gene symbol: foldchange_all.csv
   has 8 NaN 'Gene names', qc_limma.csv 8, ipa_input.csv 2,
   single_condition_proteins.csv 15, onoff_proteins.csv 0. Gene-name columns
   are therefore nullable=True everywhere.

4. single_condition_proteins.csv numeric-junk accessions are QUARANTINED,
   not excepted (DECISIONS_LOG D11).
   Exactly 2 rows of the raw SILAC Light sheet have an accession that is NOT
   a UniProt accession at all -- it is a ';'-joined list of bare MaxQuant row
   indices (32,759 and 681 characters), paired with a NaN gene and NaN in all
   four intensity columns. They reach single_condition_proteins.csv through
   foldchange.py's single-condition rescue path, which passes the raw sheets'
   accession column through unmodified.

   An EARLIER version of this module carved a narrowly-scoped exception in
   `_single_condition_accession_ok` so those 2 rows would PASS validation.
   That is precisely what D11 rejects: a validator that is taught to accept
   its own known-bad input stops being a validator. The exception is gone.
   The 2 rows are now removed from single_condition_proteins.csv (606 -> 604)
   and recorded, in full and with a reason, in
   results/qc/quarantine_accessions.csv. Nothing is deleted -- the quarantine
   file is the auditable record, and `accession` here is an ordinary strict
   accession column again.

   Discrimination is by TOKEN SHAPE, never by length, and is delegated to
   `proteomics_de.etl.accessions.is_junk_index_list` so exactly one
   implementation exists: the two 69-character
   'P08752;P20612;Q9DC51;...' rows are legitimate multi-protein groups and
   stay, while a 681-character string of nothing but digits and semicolons
   does not.

5. Contract files beyond the original five.
   `results/de/{intensity_matrix,design,limma_results}.tsv` and
   `results/ipa_input_full.csv` were added by later layers and had no schema,
   so they could drift unobserved -- which is the exact failure this module
   had already suffered once (D10 appended four columns to qc_limma.csv and
   `strict=True` turned that into a hard validation failure). They are
   schema'd here and covered by the schema-drift test in
   tests/test_schemas.py, which asserts every schema's declared column set
   equals its file's actual column set.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pandera.pandas as pa

# The repo root has to be importable for the `proteomics_de.etl.*` package form.
# This module is loaded two ways -- `import schema` after validate.py prepends
# proteomics_de/qc/ to sys.path, and `from proteomics_de.qc import schema` -- and
# only the second puts the repo root on sys.path by itself.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

# The one documented accession policy (etl/accessions.py). Junk-index-list
# detection is NOT reimplemented here: D11 turns on the difference between a
# 681-character bare-integer list (junk) and a 69-character protein group
# (legitimate), and two implementations of that test would eventually disagree.
from proteomics_de.etl.accessions import is_junk_index_list  # noqa: E402

# ---------------------------------------------------------------------------
# Shared vocabulary (mirrors config/config.yaml; kept independent on purpose
# -- this module must be able to validate the CSVs even if config.yaml is
# absent or malformed, since config.yaml is documentation, not a dependency
# of the frozen pipeline).
# ---------------------------------------------------------------------------

REGULATED_CATEGORIES = ["UP", "DOWN", "NO CHANGE", "ON_OFF"]
ONOFF_CATEGORIES = ["on_with_treatment", "off_with_treatment"]
DETECTED_IN_CATEGORIES = ["control_only", "treated_only"]
GROUP_CATEGORIES = ["control", "treated"]
DIRECTION_CATEGORIES = ["UP", "DOWN", "NS"]

INTENSITY_COLUMNS = [
    "Intensity 31578",  # treated (testosterone) / Light, replicate 1
    "Intensity 31580",  # treated (testosterone) / Light, replicate 2
    "Intensity 31579",  # control (vehicle) / Heavy, replicate 1
    "Intensity 31581",  # control (vehicle) / Heavy, replicate 2
]

# ---------------------------------------------------------------------------
# UniProt accession validation
# ---------------------------------------------------------------------------
# Branch 1 is the task spec's literal regex, verbatim (non-[OPQ] accessions,
# 6 or 10 chars, optional '-N' isoform suffix).
# Branch 2 extends it to [OPQ]-prefixed accessions (always 6 chars) -- see
# module docstring point 1 for why this is necessary on real data.
UNIPROT_TOKEN_RE = re.compile(
    r"^(?:[A-NR-Z0-9](?:[A-Z0-9]{5}|[A-Z0-9]{9})|[OPQ][A-Z0-9]{5})(?:-\d+)?$"
)


def _token_ok(token: str) -> bool:
    return bool(UNIPROT_TOKEN_RE.match(token))


def _accession_field_ok(value: object) -> bool:
    """True if `value` is one UniProt accession, or several joined by ';'
    (a multi-protein-group entry -- see module docstring point 2). NaN is
    invalid here; NaN accessions do not occur in any of the five files
    (confirmed by inspection), so accession columns are non-nullable and NaN
    is correctly rejected as a validation failure, not silently skipped.
    """
    if pd.isna(value):
        return False
    return all(_token_ok(tok) for tok in str(value).split(";"))


def is_quarantinable(value: object) -> bool:
    """True if `value` is a MaxQuant row-index list masquerading as an accession.

    Thin, deliberate re-export of the one implementation in
    ``etl/accessions.py``. A value that answers True here is not "an accession
    the schema tolerates" -- it is not an accession at all, and belongs in
    results/qc/quarantine_accessions.csv rather than in a pipeline output
    (DECISIONS_LOG D11).
    """
    return is_junk_index_list(value)


def accession_column() -> pa.Column:
    return pa.Column(
        str,
        pa.Check(
            _accession_field_ok,
            element_wise=True,
            error="not a valid UniProt accession (or ';'-joined protein group)",
        ),
        nullable=False,
    )


def gene_column() -> pa.Column:
    # NaN genes occur (see module docstring point 3); free text otherwise,
    # itself sometimes ';'-joined in lockstep with a multi-accession group
    # (e.g. 'Sprr2d;Sprr2k;Sprr2g;Sprr2b;Sprr2h').
    return pa.Column(str, nullable=True)


def intensity_column() -> pa.Column:
    # "intensities >= 0 or NaN": NaN means genuinely absent/not measured in
    # that channel (e.g. single_condition_proteins.csv), 0 means measured-
    # and-zero (e.g. onoff_proteins.csv); both are valid, negative is not.
    return pa.Column(float, pa.Check.ge(0), nullable=True, coerce=True)


def _isfinite_check() -> pa.Check:
    return pa.Check(lambda s: np.isfinite(s), element_wise=False, error="must be finite (no NaN/inf)")


def _log2fc_finite_when_complete(df: pd.DataFrame) -> pd.Series:
    """log2FC must be finite exactly on rows where complete==True; rows
    where complete==False are single/partial-condition rows and are
    EXPECTED to carry a NaN log2FC (no ratio can be computed)."""
    ok = pd.Series(True, index=df.index)
    complete_mask = df["complete"].astype(bool)
    ok.loc[complete_mask] = np.isfinite(df.loc[complete_mask, "log2FC"])
    return ok


def _row_count_check(expected: int) -> pa.Check:
    return pa.Check(lambda df: len(df) == expected, error=f"expected {expected} rows")


# ---------------------------------------------------------------------------
# foldchange_all.csv
# ---------------------------------------------------------------------------

FOLDCHANGE_ALL_EXPECTED_ROWS = 1948

FOLDCHANGE_ALL_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": accession_column(),
        "Gene names": gene_column(),
        "Intensity 31578": intensity_column(),
        "Intensity 31580": intensity_column(),
        "Intensity 31579": intensity_column(),
        "Intensity 31581": intensity_column(),
        "ratio_rep1": pa.Column(float, nullable=True, coerce=True),
        "ratio_rep2": pa.Column(float, nullable=True, coerce=True),
        "log2_rep1": pa.Column(float, nullable=True, coerce=True),
        "log2_rep2": pa.Column(float, nullable=True, coerce=True),
        "log2FC": pa.Column(float, nullable=True, coerce=True),
        "complete": pa.Column(bool, coerce=True),
        "regulated": pa.Column(str, pa.Check.isin(REGULATED_CATEGORIES), coerce=True),
        "onoff": pa.Column(str, pa.Check.isin(ONOFF_CATEGORIES), nullable=True, coerce=True),
    },
    checks=[
        pa.Check(
            _log2fc_finite_when_complete,
            error="log2FC must be finite for complete==True rows, NaN otherwise",
        ),
        _row_count_check(FOLDCHANGE_ALL_EXPECTED_ROWS),
    ],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# qc_limma.csv
# ---------------------------------------------------------------------------

QC_LIMMA_EXPECTED_ROWS = 1938

QC_LIMMA_SCHEMA = pa.DataFrameSchema(
    {
        "id": accession_column(),
        "gene": gene_column(),
        "Intensity 31578": intensity_column(),
        "Intensity 31580": intensity_column(),
        "Intensity 31579": intensity_column(),
        "Intensity 31581": intensity_column(),
        "limma_log2FC": pa.Column(float, _isfinite_check(), coerce=True),
        "p_value": pa.Column(float, pa.Check.in_range(0, 1), coerce=True),
        "adj_p_value": pa.Column(float, pa.Check.in_range(0, 1), coerce=True),
        "significant": pa.Column(bool, coerce=True),
        "regulated": pa.Column(str, pa.Check.isin(REGULATED_CATEGORIES), coerce=True),
        # --- DECISIONS_LOG D10: restored limma columns -------------------
        # These four were appended to qc_limma.csv by the D10 change and, with
        # strict=True, immediately turned every validation run into a
        # `column_in_schema` failure. That is the schema doing its job; the fix
        # is to declare them, which is what the schema-drift test in
        # tests/test_schemas.py now enforces for every file.
        #
        # n_imputed counts how many of the 4 intensities were invented by
        # MinProb rather than measured, so its ceiling is the number of
        # samples. in_range(0, 4) is that structural bound, not a description
        # of today's data (which only reaches 2: 1578/218/142 rows at 0/1/2).
        # Pinning it to 2 would make the schema fail on a legitimately
        # sparser protein.
        "n_imputed": pa.Column(int, pa.Check.in_range(0, 4), coerce=True),
        "AveExpr": pa.Column(float, _isfinite_check(), coerce=True),
        "t": pa.Column(float, _isfinite_check(), coerce=True),
        "B": pa.Column(float, _isfinite_check(), coerce=True),
    },
    checks=[_row_count_check(QC_LIMMA_EXPECTED_ROWS)],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# ipa_input.csv
# ---------------------------------------------------------------------------

IPA_INPUT_EXPECTED_ROWS = 715

IPA_INPUT_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": accession_column(),
        "Gene names": gene_column(),
        "log2FC": pa.Column(float, _isfinite_check(), coerce=True),
        "regulated": pa.Column(str, pa.Check.isin(REGULATED_CATEGORIES), coerce=True),
    },
    checks=[_row_count_check(IPA_INPUT_EXPECTED_ROWS)],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# single_condition_proteins.csv
# ---------------------------------------------------------------------------

# 606 before DECISIONS_LOG D11, 604 after: the 2 bare-integer MaxQuant
# row-index lists are quarantined to results/qc/quarantine_accessions.csv
# instead of being waved through by a schema exception. Kept in lockstep with
# tests/expected/frozen_counts.json:single_condition_rows.
SINGLE_CONDITION_EXPECTED_ROWS = 604

SINGLE_CONDITION_SCHEMA = pa.DataFrameSchema(
    {
        # D11: an ordinary strict accession column again. The former
        # dataframe-level `_single_condition_accession_ok` exception -- which
        # forgave the 2 junk rows when they came paired with a NaN gene and
        # all-NaN intensities -- is deleted, not relaxed. Junk is quarantined
        # upstream (qc/boundaries.py), so by the time a row reaches this file
        # its accession must be a real UniProt accession or protein group.
        "accession": accession_column(),
        "gene": gene_column(),
        "detected_in": pa.Column(str, pa.Check.isin(DETECTED_IN_CATEGORIES), coerce=True),
        "Intensity 31578": intensity_column(),
        "Intensity 31580": intensity_column(),
        "Intensity 31579": intensity_column(),
        "Intensity 31581": intensity_column(),
    },
    checks=[_row_count_check(SINGLE_CONDITION_EXPECTED_ROWS)],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# onoff_proteins.csv
# ---------------------------------------------------------------------------

ONOFF_EXPECTED_ROWS = 10

ONOFF_SCHEMA = pa.DataFrameSchema(
    {
        "accession": accession_column(),
        "gene": gene_column(),
        "onoff": pa.Column(str, pa.Check.isin(ONOFF_CATEGORIES), coerce=True),
        "Intensity 31578": intensity_column(),
        "Intensity 31580": intensity_column(),
        "Intensity 31579": intensity_column(),
        "Intensity 31581": intensity_column(),
    },
    checks=[_row_count_check(ONOFF_EXPECTED_ROWS)],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# results/de/intensity_matrix.tsv  (the frame handed to limma)
# ---------------------------------------------------------------------------

INTENSITY_MATRIX_EXPECTED_ROWS = 1938

# Column ORDER here is the acquisition-independent order the DE layer writes:
# treated channels first is NOT assumed -- the file lists 31579/31581 before
# 31578/31580, and after the D7 correction those are the CONTROL channels.
# Order is not enforced by pandera (a schema is a set of named columns), but
# the schema-drift test compares sets, so a reordering is deliberately not a
# failure while an added/removed/renamed column is.
INTENSITY_MATRIX_SCHEMA = pa.DataFrameSchema(
    {
        "accession": accession_column(),
        "gene": gene_column(),
        "Intensity 31579": intensity_column(),
        "Intensity 31581": intensity_column(),
        "Intensity 31578": intensity_column(),
        "Intensity 31580": intensity_column(),
    },
    checks=[_row_count_check(INTENSITY_MATRIX_EXPECTED_ROWS)],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# results/de/design.tsv  (the model.matrix(~group) source of truth)
# ---------------------------------------------------------------------------

DESIGN_EXPECTED_ROWS = 4

DESIGN_SCHEMA = pa.DataFrameSchema(
    {
        # Sample ids are the numeric MaxQuant run ids (31578..31581), written
        # bare, so pandas reads them as int64. They are identifiers, not
        # measurements: uniqueness is the contract that matters.
        "sample": pa.Column(int, unique=True, coerce=True),
        "group": pa.Column(str, pa.Check.isin(GROUP_CATEGORIES), coerce=True),
    },
    checks=[
        _row_count_check(DESIGN_EXPECTED_ROWS),
        pa.Check(
            lambda df: df["group"].nunique() == 2,
            error="design.tsv must describe exactly two groups (a two-group contrast)",
        ),
        pa.Check(
            lambda df: df["group"].value_counts().nunique() == 1,
            error="design.tsv groups must be balanced (equal samples per group)",
        ),
    ],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# results/de/limma_results.tsv  (limma's own output, before renaming)
# ---------------------------------------------------------------------------

LIMMA_RESULTS_EXPECTED_ROWS = 1938

LIMMA_RESULTS_SCHEMA = pa.DataFrameSchema(
    {
        "accession": accession_column(),
        "gene": gene_column(),
        "logFC": pa.Column(float, _isfinite_check(), coerce=True),
        # fold_change == 2**logFC, so it is strictly positive and never zero.
        "fold_change": pa.Column(float, pa.Check.gt(0), coerce=True),
        "AveExpr": pa.Column(float, _isfinite_check(), coerce=True),
        "t": pa.Column(float, _isfinite_check(), coerce=True),
        "P.Value": pa.Column(float, pa.Check.in_range(0, 1), coerce=True),
        "adj.P.Val": pa.Column(float, pa.Check.in_range(0, 1), coerce=True),
        "B": pa.Column(float, _isfinite_check(), coerce=True),
        "n_imputed": pa.Column(int, pa.Check.in_range(0, 4), coerce=True),
        # Today every row is "NS" (0/1938 survive FDR<0.05, DECISIONS_LOG D2).
        # The schema still admits UP/DOWN: pinning it to {"NS"} would encode
        # the null result as a contract and fail the moment the study gains
        # biological replicates, which is the FORWARD-PATH case config.yaml
        # is explicitly designed for.
        "direction": pa.Column(str, pa.Check.isin(DIRECTION_CATEGORIES), coerce=True),
    },
    checks=[
        _row_count_check(LIMMA_RESULTS_EXPECTED_ROWS),
        pa.Check(
            lambda df: np.allclose(df["fold_change"], np.power(2.0, df["logFC"])),
            error="fold_change must equal 2**logFC",
        ),
        pa.Check(
            lambda df: (df["adj.P.Val"] >= df["P.Value"] - 1e-12).all(),
            error="BH-adjusted p-values must be >= their raw p-values",
        ),
    ],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# results/ipa_input_full.csv  (IPA upload + the p-values IPA can weight by)
# ---------------------------------------------------------------------------

IPA_INPUT_FULL_EXPECTED_ROWS = 715

IPA_INPUT_FULL_SCHEMA = pa.DataFrameSchema(
    {
        "UniProt Accession Number": accession_column(),
        "Gene names": gene_column(),
        "log2FC": pa.Column(float, _isfinite_check(), coerce=True),
        # Only regulated proteins are exported, so NO CHANGE / ON_OFF cannot
        # appear -- a narrower vocabulary than REGULATED_CATEGORIES on purpose.
        "regulated": pa.Column(str, pa.Check.isin(["UP", "DOWN"]), coerce=True),
        "p_value": pa.Column(float, pa.Check.in_range(0, 1), coerce=True),
        "adj_p_value": pa.Column(float, pa.Check.in_range(0, 1), coerce=True),
    },
    checks=[
        _row_count_check(IPA_INPUT_FULL_EXPECTED_ROWS),
        pa.Check(
            lambda df: (df["adj_p_value"] >= df["p_value"] - 1e-12).all(),
            error="BH-adjusted p-values must be >= their raw p-values",
        ),
    ],
    strict=True,
    coerce=True,
)

# ---------------------------------------------------------------------------
# Registry consumed by validate.py
#
# Keys are paths RELATIVE TO results/, because the later layers write into
# subdirectories (results/de/). `READERS` records the separator, since the DE
# layer emits TSV and everything else CSV -- a validator that assumed ',' would
# read a TSV as one giant column and report a meaningless failure.
# ---------------------------------------------------------------------------

FILE_SCHEMAS: dict[str, pa.DataFrameSchema] = {
    "foldchange_all.csv": FOLDCHANGE_ALL_SCHEMA,
    "qc_limma.csv": QC_LIMMA_SCHEMA,
    "ipa_input.csv": IPA_INPUT_SCHEMA,
    "single_condition_proteins.csv": SINGLE_CONDITION_SCHEMA,
    "onoff_proteins.csv": ONOFF_SCHEMA,
    "de/intensity_matrix.tsv": INTENSITY_MATRIX_SCHEMA,
    "de/design.tsv": DESIGN_SCHEMA,
    "de/limma_results.tsv": LIMMA_RESULTS_SCHEMA,
    "ipa_input_full.csv": IPA_INPUT_FULL_SCHEMA,
}

EXPECTED_ROWS: dict[str, int] = {
    "foldchange_all.csv": FOLDCHANGE_ALL_EXPECTED_ROWS,
    "qc_limma.csv": QC_LIMMA_EXPECTED_ROWS,
    "ipa_input.csv": IPA_INPUT_EXPECTED_ROWS,
    "single_condition_proteins.csv": SINGLE_CONDITION_EXPECTED_ROWS,
    "onoff_proteins.csv": ONOFF_EXPECTED_ROWS,
    "de/intensity_matrix.tsv": INTENSITY_MATRIX_EXPECTED_ROWS,
    "de/design.tsv": DESIGN_EXPECTED_ROWS,
    "de/limma_results.tsv": LIMMA_RESULTS_EXPECTED_ROWS,
    "ipa_input_full.csv": IPA_INPUT_FULL_EXPECTED_ROWS,
}

#: Separator each contract file is written with.
SEPARATORS: dict[str, str] = {name: "\t" if name.endswith(".tsv") else "," for name in FILE_SCHEMAS}


def declared_columns(name: str) -> list[str]:
    """Column names a schema declares, in declaration order.

    Used by the schema-drift test: a schema whose declared column set has
    fallen behind its file's actual columns is the exact failure that blocked
    this package (D10 appended four columns to qc_limma.csv and `strict=True`
    turned every run red).
    """
    return list(FILE_SCHEMAS[name].columns)


if __name__ == "__main__":
    # Convenience alias: `python -m proteomics_de.qc.schema` (or a direct
    # `python proteomics_de/qc/schema.py` run) delegates to the full
    # validator in validate.py rather than duplicating its CLI logic.
    from validate import main
    main()
