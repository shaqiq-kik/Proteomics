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

4. single_condition_proteins.csv numeric-junk accession quirk (2 rows).
   Exactly 2 rows in the committed single_condition_proteins.csv have an
   'accession' value that is NOT a UniProt accession at all -- it is a long
   ';'-joined list of bare integers (e.g. '73;74;75;...;140;3904;...'),
   paired with a NaN 'gene' and NaN in ALL FOUR intensity columns. This
   traces to foldchange.py's single-condition rescue path, which passes the
   raw sheets' 'UniProt Accession Number' / 'Gene names' columns straight
   through unmodified; these 2 rows are themselves malformed in the ORIGINAL
   raw SILAC sheet (almost certainly an unresolved MaxQuant protein-group
   row-index list rather than a resolved accession), not something the
   frozen script introduced. Per the guardrail that existing scripts/outputs
   are frozen and must not be modified, this schema does NOT loosen the
   accession regex globally to paper over this -- it instead carves out a
   narrowly-scoped, auditable exception (see
   `_single_condition_accession_ok`) that only forgives an accession-regex
   failure when it is ALSO paired with a NaN gene AND all-NaN intensities
   AND looks like the specific digits-and-semicolons pattern observed. Any
   other malformed accession still fails validation.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import pandera.pandas as pa

# ---------------------------------------------------------------------------
# Shared vocabulary (mirrors config/config.yaml; kept independent on purpose
# -- this module must be able to validate the CSVs even if config.yaml is
# absent or malformed, since config.yaml is documentation, not a dependency
# of the frozen pipeline).
# ---------------------------------------------------------------------------

REGULATED_CATEGORIES = ["UP", "DOWN", "NO CHANGE", "ON_OFF"]
ONOFF_CATEGORIES = ["on_with_treatment", "off_with_treatment"]
DETECTED_IN_CATEGORIES = ["control_only", "treated_only"]

INTENSITY_COLUMNS = [
    "Intensity 31578",  # control / Light, replicate 1
    "Intensity 31580",  # control / Light, replicate 2
    "Intensity 31579",  # treated / Heavy, replicate 1
    "Intensity 31581",  # treated / Heavy, replicate 2
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


def _single_condition_accession_ok(df: pd.DataFrame) -> pd.Series:
    """Row-wise accession validity for single_condition_proteins.csv, with
    the narrowly-scoped documented-quirk exception from module docstring
    point 4. A row passes if EITHER:
      (a) its accession is a valid UniProt accession / protein group, OR
      (b) it matches the exact observed quirk signature: gene is NaN, all
          four intensity columns are NaN, and accession is composed only of
          digits and ';' (i.e. clearly not a UniProt accession at all, and
          co-occurring with a row that carries no other usable data).
    Any malformed accession that does NOT also satisfy (b) still fails.
    """
    acc = df["accession"].astype(str)
    is_valid_accession = df["accession"].apply(_accession_field_ok)
    is_documented_quirk = (
        df["gene"].isna()
        & df[INTENSITY_COLUMNS].isna().all(axis=1)
        & acc.str.match(r"^[0-9;]+$")
    )
    return is_valid_accession | is_documented_quirk


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

SINGLE_CONDITION_EXPECTED_ROWS = 606

SINGLE_CONDITION_SCHEMA = pa.DataFrameSchema(
    {
        # accession is validated at the dataframe level below (needs gene +
        # intensity columns in scope for the documented-quirk exception).
        "accession": pa.Column(str, nullable=False),
        "gene": gene_column(),
        "detected_in": pa.Column(str, pa.Check.isin(DETECTED_IN_CATEGORIES), coerce=True),
        "Intensity 31578": intensity_column(),
        "Intensity 31580": intensity_column(),
        "Intensity 31579": intensity_column(),
        "Intensity 31581": intensity_column(),
    },
    checks=[
        pa.Check(
            _single_condition_accession_ok,
            error=(
                "accession is not a valid UniProt accession/protein-group and "
                "does not match the documented all-NaN numeric-index quirk "
                "(see schema.py module docstring point 4)"
            ),
        ),
        _row_count_check(SINGLE_CONDITION_EXPECTED_ROWS),
    ],
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
# Registry consumed by validate.py
# ---------------------------------------------------------------------------

FILE_SCHEMAS: dict[str, pa.DataFrameSchema] = {
    "foldchange_all.csv": FOLDCHANGE_ALL_SCHEMA,
    "qc_limma.csv": QC_LIMMA_SCHEMA,
    "ipa_input.csv": IPA_INPUT_SCHEMA,
    "single_condition_proteins.csv": SINGLE_CONDITION_SCHEMA,
    "onoff_proteins.csv": ONOFF_SCHEMA,
}

EXPECTED_ROWS: dict[str, int] = {
    "foldchange_all.csv": FOLDCHANGE_ALL_EXPECTED_ROWS,
    "qc_limma.csv": QC_LIMMA_EXPECTED_ROWS,
    "ipa_input.csv": IPA_INPUT_EXPECTED_ROWS,
    "single_condition_proteins.csv": SINGLE_CONDITION_EXPECTED_ROWS,
    "onoff_proteins.csv": ONOFF_EXPECTED_ROWS,
}


if __name__ == "__main__":
    # Convenience alias: `python -m proteomics_de.qc.schema` (or a direct
    # `python proteomics_de/qc/schema.py` run) delegates to the full
    # validator in validate.py rather than duplicating its CLI logic.
    from validate import main
    main()
