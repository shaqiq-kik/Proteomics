"""Tests for the pandera file contracts in ``qc/schema.py``.

THE TEST THAT MATTERS MOST IS THE FIRST ONE.

This package was blocked by a schema silently falling behind its file:
DECISIONS_LOG D10 appended ``n_imputed``, ``AveExpr``, ``t`` and ``B`` to
``results/qc_limma.csv``, ``QC_LIMMA_SCHEMA`` is ``strict=True``, and every
validation run from that moment on died with ``column_in_schema``. Nothing
tested that a schema's declared columns still matched its file's actual
columns, so the drift was only discovered when someone ran the validator by
hand.

``test_schema_columns_match_the_file`` closes that. It is parametrised over
every entry in ``FILE_SCHEMAS``, so a new contract file cannot be added without
being covered, and an added/removed/renamed column fails immediately and names
the delta.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pandera.pandas as pa
import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPO_ROOT = _PKG_DIR.parent

# The pipeline scripts are flat modules; mirror the sys.path layout they expect
# (conftest.py does this too, but this file must also work standalone -- and
# conftest's fixture runs after collection imports this module).
for _entry in (_PKG_DIR / "qc", _PKG_DIR, _REPO_ROOT):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))

from proteomics_de.etl import accessions  # noqa: E402
from proteomics_de.qc import schema as qc_schema  # noqa: E402

_RESULTS = _PKG_DIR / "results"

ALL_FILES = sorted(qc_schema.FILE_SCHEMAS)


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(_RESULTS / name, sep=qc_schema.SEPARATORS[name])


# ---------------------------------------------------------------------------
# 1. Schema drift -- the failure this package was blocked on
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", ALL_FILES)
def test_schema_columns_match_the_file(name):
    """Declared columns == actual columns, for every contract file.

    Compared as SETS with the symmetric difference reported, so the message
    names exactly which columns drifted rather than dumping two lists. Order is
    deliberately not asserted: a pandera schema is a mapping of names to
    columns, and reordering a file's columns is not a contract break.
    """
    declared = set(qc_schema.declared_columns(name))
    actual = set(_read(name).columns)

    missing_from_schema = actual - declared
    missing_from_file = declared - actual
    assert not missing_from_schema and not missing_from_file, (
        f"{name}: schema drift.\n"
        f"  in the file but NOT declared: {sorted(missing_from_schema)}\n"
        f"  declared but NOT in the file: {sorted(missing_from_file)}"
    )


@pytest.mark.parametrize("name", ALL_FILES)
def test_every_contract_file_validates(name):
    """The committed file passes its own schema, with the failures reported."""
    try:
        qc_schema.FILE_SCHEMAS[name].validate(_read(name), lazy=True)
    except pa.errors.SchemaErrors as exc:  # pragma: no cover - failure path
        cases = exc.failure_cases[["column", "check", "failure_case"]].head(10)
        pytest.fail(f"{name} failed its schema:\n{cases.to_string()}")


@pytest.mark.parametrize("name", ALL_FILES)
def test_expected_rows_matches_the_file(name):
    assert len(_read(name)) == qc_schema.EXPECTED_ROWS[name]


def test_every_schema_is_registered_with_a_row_count():
    """No schema may exist in the registry without an expected row count."""
    assert set(qc_schema.FILE_SCHEMAS) == set(qc_schema.EXPECTED_ROWS)
    assert set(qc_schema.FILE_SCHEMAS) == set(qc_schema.SEPARATORS)


def test_the_d10_columns_are_declared():
    """Regression pin for the exact drift that blocked this package.

    A generic drift test would go green again if someone deleted the columns
    from BOTH the file and the schema. D10 decided these four ship; this says so
    by name.
    """
    declared = qc_schema.declared_columns("qc_limma.csv")
    for column in ("n_imputed", "AveExpr", "t", "B"):
        assert column in declared, f"D10 column {column!r} is not declared"


def test_config_yaml_key_columns_agree_with_the_schemas():
    """config.yaml documents the same contract the schemas enforce.

    config.yaml is descriptive documentation, so it can drift silently in a way
    no validator would catch -- which is exactly how its qc_limma key_columns
    list came to be missing the four D10 columns.
    """
    yaml = pytest.importorskip("yaml")
    config = yaml.safe_load(
        (_PKG_DIR / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    contracts = config["file_contracts"]
    by_path = {
        Path(entry["path"]).name: entry
        for entry in contracts.values()
    }
    for name in ALL_FILES:
        entry = by_path.get(Path(name).name)
        assert entry is not None, f"{name} has a schema but no config.yaml file_contract"
        assert set(entry["key_columns"]) == set(qc_schema.declared_columns(name)), (
            f"{name}: config.yaml key_columns disagree with qc/schema.py"
        )
        assert entry["expected_rows"] == qc_schema.EXPECTED_ROWS[name], (
            f"{name}: config.yaml expected_rows disagrees with qc/schema.py"
        )


# ---------------------------------------------------------------------------
# 2. D11 -- the junk accessions are quarantined, not excepted
# ---------------------------------------------------------------------------
def test_the_junk_accession_exception_is_gone():
    """The exception that let the 2 junk rows PASS must not come back.

    D11's objection was not "the schema is wrong" but "the schema was taught to
    accept its own known-bad input". Checked against the module's public
    surface AND its source text, because the failure mode is someone
    reintroducing an equivalent helper under a different name.
    """
    assert not hasattr(qc_schema, "_single_condition_accession_ok")
    source = (_PKG_DIR / "qc" / "schema.py").read_text(encoding="utf-8")
    assert "is_documented_quirk" not in source


def test_a_junk_accession_now_fails_the_single_condition_schema():
    """Positive proof the exception is really gone.

    Deleting a helper proves nothing on its own; this feeds the schema a row
    shaped exactly like the ones that used to be forgiven -- junk accession,
    NaN gene, all-NaN intensities -- and requires it to FAIL.
    """
    junk = pd.DataFrame(
        [
            {
                "accession": "73;74;75;76;77",
                "gene": None,
                "detected_in": "treated_only",
                "Intensity 31578": None,
                "Intensity 31580": None,
                "Intensity 31579": None,
                "Intensity 31581": None,
            }
        ]
    )
    with pytest.raises(pa.errors.SchemaErrors):
        qc_schema.SINGLE_CONDITION_SCHEMA.validate(junk, lazy=True)


def test_quarantine_discriminates_on_token_shape_not_length():
    """The 69-char protein groups stay; the 681-char index list goes.

    This is the whole substance of D11. A length threshold anywhere between
    them would 'work' on today's data and be wrong in principle, so both
    directions are pinned.
    """
    legit = "P08752;P20612;Q9DC51;Q3V3I2;P50149;P18872;B2RSH2;Q8CGK7;P63094;Q6R0H7"
    junk = ";".join(str(i) for i in range(200))
    assert len(legit) == 69
    assert len(junk) > len(legit)

    assert not qc_schema.is_quarantinable(legit)
    assert qc_schema.is_quarantinable(junk)
    # ...and the short junk list is caught too, so length plays no part at all.
    assert qc_schema.is_quarantinable("73;74;75")


def test_is_quarantinable_is_the_shared_helper():
    """One implementation, not two that can drift apart."""
    assert qc_schema.is_quarantinable("1;2;3") == accessions.is_junk_index_list("1;2;3")
    source = (_PKG_DIR / "qc" / "schema.py").read_text(encoding="utf-8")
    assert "from proteomics_de.etl.accessions import is_junk_index_list" in source


def test_committed_single_condition_file_holds_no_junk():
    scp = _read("single_condition_proteins.csv")
    junk = [a for a in scp["accession"] if qc_schema.is_quarantinable(a)]
    assert junk == [], f"{len(junk)} junk accession(s) still shipped"


def test_the_legitimate_protein_groups_survived_the_quarantine():
    """Guard against over-quarantining: 604, not 602."""
    scp = _read("single_condition_proteins.csv")
    long_groups = [a for a in scp["accession"] if len(str(a)) == 69]
    assert len(long_groups) == 2, "the 69-character protein groups were dropped"
    assert all(accessions.is_valid_group(a) for a in long_groups)


def test_quarantine_file_records_both_rows_with_a_reason():
    quarantine = pd.read_csv(_RESULTS / "qc" / "quarantine_accessions.csv")
    assert len(quarantine) == 2
    assert quarantine["reason"].notna().all()
    # Quarantine means set aside, not deleted: the original value is recoverable.
    assert sorted(quarantine["n_chars"]) == [681, 32759]
    assert (quarantine["value"].str.len() == quarantine["n_chars"]).all()
    assert quarantine["source"].eq("single_condition_proteins.csv").all()


# ---------------------------------------------------------------------------
# 3. The counts move together
# ---------------------------------------------------------------------------
def test_frozen_counts_track_the_quarantine(frozen_counts):
    assert frozen_counts["single_condition_rows"] == 604
    assert frozen_counts["background_union"] == 2552
    assert (
        frozen_counts["foldchange_all_rows"] + frozen_counts["single_condition_rows"]
        == frozen_counts["background_union"]
    )


def test_the_quarantine_moved_nothing_else(frozen_counts):
    """D11 touches the rescue file and the background, and nothing else.

    The quarantined rows are single-condition, so they never had a fold change,
    never entered limma and never reached the IPA export. If any of these
    counts moved, something other than the quarantine did it.
    """
    assert frozen_counts["foldchange_all_rows"] == 1948
    assert frozen_counts["qc_limma_rows"] == 1938
    assert frozen_counts["ipa_input_rows"] == 715
    assert frozen_counts["onoff_rows"] == 10
    assert frozen_counts["n_up"] == 509
    assert frozen_counts["n_down"] == 206


def test_config_yaml_background_agrees_with_frozen_counts(frozen_counts):
    yaml = pytest.importorskip("yaml")
    config = yaml.safe_load(
        (_PKG_DIR / "config" / "config.yaml").read_text(encoding="utf-8")
    )
    assert (
        config["enrichment"]["gprofiler"]["background_n"]
        == frozen_counts["background_union"]
    )


def test_validate_py_does_not_hardcode_the_background():
    """The background union is read from frozen_counts.json, not typed.

    Same rule ``test_no_dataset_counts_are_hardcoded_in_enrich_common``
    applies to enrich_common.py: 2554 was a literal here and D11 moved it.
    AST-based so prose in docstrings and comments is not mistaken for code.
    """
    import ast

    tree = ast.parse((_PKG_DIR / "qc" / "validate.py").read_text(encoding="utf-8"))
    literals = [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, int)
        and not isinstance(node.value, bool)
        and node.value in {2552, 2554}
    ]
    assert literals == [], f"background union hardcoded in validate.py: {literals}"


def test_qc_report_is_green(results_dir):
    """The committed validation report reflects a passing run."""
    report = json.loads((results_dir / "qc" / "qc_report.json").read_text())
    assert report["overall_passed"] is True
    assert {f["file"] for f in report["files"]} == set(ALL_FILES)
    assert all(f["passed"] for f in report["files"])
    assert all(c["passed"] for c in report["cross_file_checks"])
