"""Lock for the accession-handling policy in ``proteomics_de/etl/accessions.py``.

The real committed data is loaded where the assertion is about the real data --
the isoform count and the junk-row count are facts about today's inputs, and the
point of the test is to notice when they stop being true.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.etl import accessions  # noqa: E402

_RESULTS = _REPO_ROOT / "proteomics_de" / "results"
FOLDCHANGE_CSV = _RESULTS / "foldchange_all.csv"
SINGLE_CONDITION_CSV = _RESULTS / "single_condition_proteins.csv"
QUARANTINE_CSV = _RESULTS / "qc" / "quarantine_accessions.csv"


# ---------------------------------------------------------------------------
# is_valid_uniprot -- including the accessions the spec's regex wrongly rejects
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "token",
    [
        "P12345",       # canonical
        "Q9JHU4",       # Q-prefixed; rejected by the task spec's regex
        "O08528",       # O-prefixed; rejected by the task spec's regex
        "P19137",       # P-prefixed; rejected by the task spec's regex
        "A0A023GPI8",   # 10-character accession
        "P12345-2",     # isoform suffix is still a valid accession
    ],
)
def test_valid_uniprot_tokens(token):
    assert accessions.is_valid_uniprot(token) is True


@pytest.mark.parametrize(
    "token",
    ["", "73", "P1234", "P123456", "p12345", "P05132;P68181", "nan", "-2"],
)
def test_invalid_uniprot_tokens(token):
    assert accessions.is_valid_uniprot(token) is False


def test_group_validity_checks_every_token():
    assert accessions.is_valid_group("P05132;P68181") is True
    assert accessions.is_valid_group("P05132;NOTANACC") is False
    assert accessions.is_valid_group("73;74;75") is False
    assert accessions.is_valid_group(float("nan")) is False


# ---------------------------------------------------------------------------
# first_token / split_group -- must stay behaviour-identical to enrich_common
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "value",
    ["P12345", "P05132;P68181", "73;74;75", "", float("nan"), " P12345 ; P68181 "],
)
def test_first_token_matches_legacy_expression(value):
    """Byte-for-byte the same as ``str(acc).split(';')[0].strip()``."""
    assert accessions.first_token(value) == str(value).split(";")[0].strip()


def test_first_token_examples():
    assert accessions.first_token("P12345") == "P12345"
    assert accessions.first_token("P05132;P68181") == "P05132"
    assert accessions.first_token("73;74;75") == "73"
    assert accessions.first_token("") == ""
    assert accessions.first_token(float("nan")) == "nan"


def test_split_group():
    assert accessions.split_group("P05132;P68181") == ["P05132", "P68181"]
    assert accessions.split_group("P12345") == ["P12345"]
    assert accessions.split_group(" P05132 ; P68181 ") == ["P05132", "P68181"]
    assert accessions.split_group("73;74;75") == ["73", "74", "75"]
    assert accessions.split_group("") == [""]


# ---------------------------------------------------------------------------
# Isoforms: detected and reported, NEVER silently stripped
# ---------------------------------------------------------------------------
def test_strip_isoform():
    assert accessions.strip_isoform("P12345-2") == "P12345"
    assert accessions.strip_isoform("P12345") == "P12345"
    assert accessions.strip_isoform("A0A023GPI8-11") == "A0A023GPI8"


def test_has_isoform():
    assert accessions.has_isoform("P12345-2") is True
    assert accessions.has_isoform("P12345") is False


def test_detect_isoforms_synthetic_row_is_flagged():
    found = accessions.detect_isoforms(["P12345", "P05132;P68181-3", float("nan")])
    assert found == ["P05132;P68181-3"]
    assert len(found) > 0


def test_detect_isoforms_is_empty_on_the_real_foldchange_data():
    """Today's committed data carries ZERO isoform suffixes.

    If this ever fails, the input changed: decide explicitly what to do with
    isoforms before letting them through. Do NOT start stripping them.
    """
    df = pd.read_csv(FOLDCHANGE_CSV)
    found = accessions.detect_isoforms(df["UniProt Accession Number"])
    assert found == [], f"{len(found)} isoform-bearing accession(s) appeared: {found[:5]}"


def test_detect_isoforms_does_not_mutate_or_strip():
    values = ["P12345-2"]
    accessions.detect_isoforms(values)
    assert values == ["P12345-2"]


# ---------------------------------------------------------------------------
# Junk index lists: quarantined, and NOT confused with long protein groups
# ---------------------------------------------------------------------------
def test_is_junk_index_list_units():
    assert accessions.is_junk_index_list("73;74;75") is True
    assert accessions.is_junk_index_list("P05132;P68181") is False
    assert accessions.is_junk_index_list("P12345") is False
    assert accessions.is_junk_index_list("") is False
    assert accessions.is_junk_index_list(float("nan")) is False
    # A single bare integer is a malformed accession, not a junk index list --
    # the accession validator should report it rather than have it quarantined.
    assert accessions.is_junk_index_list("73") is False


def test_junk_index_lists_in_the_real_single_condition_file():
    """The junk rows are OUT of the shipped file; the protein groups stayed IN.

    Rewritten at DECISIONS_LOG D11. This used to assert the 2 junk rows were
    *present* in single_condition_proteins.csv -- true when qc/schema.py carved
    an exception so they passed validation. They are now quarantined to
    results/qc/quarantine_accessions.csv and dropped from the shipped file, so
    the assertion inverts: the shipped file must be clean, and the quarantine
    file must hold exactly the 2, with their full values intact.
    """
    df = pd.read_csv(SINGLE_CONDITION_CSV)
    acc = df["accession"].astype(str)

    assert not acc.map(accessions.is_junk_index_list).any(), (
        "a junk MaxQuant row-index list is still shipped in "
        "single_condition_proteins.csv"
    )
    assert acc.map(accessions.is_valid_group).all()

    # The two 69-character legitimate protein groups must have SURVIVED. This
    # is the over-quarantining guard: a length-based rule would have taken them
    # too, leaving 602 rows instead of 604.
    long_legit = acc[acc.str.len() == 69]
    assert len(long_legit) == 2, "the 2 known 69-char protein groups were dropped"
    assert not long_legit.map(accessions.is_junk_index_list).any()
    assert long_legit.map(accessions.is_valid_group).all()

    # The junk is set aside with a reason, not deleted: 32,759 and 681 chars.
    quarantined = pd.read_csv(QUARANTINE_CSV)
    assert len(quarantined) == 2
    assert sorted(quarantined["value"].str.len()) == [681, 32759]
    assert quarantined["value"].map(accessions.is_junk_index_list).all()
    assert not quarantined["value"].map(accessions.is_valid_group).any()
    assert quarantined["reason"].notna().all()


def test_regex_is_the_same_pattern_qc_schema_settled_on():
    from proteomics_de.qc import schema

    assert accessions.UNIPROT_TOKEN_RE.pattern == schema.UNIPROT_TOKEN_RE.pattern
