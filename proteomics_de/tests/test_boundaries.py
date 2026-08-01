"""Tests for the stage-boundary validation hooks (``qc/boundaries.py``).

``boundaries.check`` was a Wave-0 no-op shim that returned its argument. Three
call sites in ``foldchange.py`` were installed against that contract before the
body existed, so the tests here fall into two groups:

1. **The shim's contract still holds.** ``check`` returns the SAME object, takes
   the same two positional arguments, and mutates nothing. If any of that broke,
   ``foldchange.py`` would have to change -- and it must not.

2. **The two policies are real and are not each other.** ``after_load`` /
   ``after_merge`` see raw MaxQuant sheets, where a malformed accession is a
   fact about the input: they record it and route to quarantine.
   ``after_foldchange`` / ``before_ipa_export`` see frames this repo produced,
   where a schema failure is a bug: they raise. Both directions are asserted,
   because a validator that only ever passes and a validator that always raises
   are equally useless and neither is caught by testing one stage.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pandera.pandas as pa
import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPO_ROOT = _PKG_DIR.parent

for _entry in (_PKG_DIR / "qc", _PKG_DIR, _REPO_ROOT):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))

from proteomics_de.qc import boundaries  # noqa: E402

_RESULTS = _PKG_DIR / "results"

JUNK_ACCESSION = ";".join(str(i) for i in range(120))
LEGIT_GROUP = "P08752;P20612;Q9DC51;Q3V3I2;P50149;P18872;B2RSH2;Q8CGK7;P63094;Q6R0H7"


@pytest.fixture(autouse=True)
def _isolated_records():
    """Every test starts with an empty record list."""
    boundaries.reset()
    yield
    boundaries.reset()


def raw_sheet(accessions=("P19137", "Q9JHU4")) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "UniProt Accession Number": list(accessions),
            "Gene names": ["Lama1", None][: len(accessions)],
            "Intensity 31578": [1.0] * len(accessions),
            "Intensity 31580": [2.0] * len(accessions),
        }
    )


def foldchange_frame() -> pd.DataFrame:
    """A minimal frame shaped like the one crossing ``after_foldchange``."""
    return pd.DataFrame(
        {
            "UniProt Accession Number": ["P19137", "Q9JHU4"],
            "Gene names": ["Lama1", "Dync1h1"],
            "Intensity 31578": [1.0, 2.0],
            "Intensity 31580": [1.0, 2.0],
            "Intensity 31579": [2.0, 4.0],
            "Intensity 31581": [2.0, 4.0],
            "ratio_rep1": [2.0, 2.0],
            "ratio_rep2": [2.0, 2.0],
            "log2_rep1": [1.0, 1.0],
            "log2_rep2": [1.0, 1.0],
            "log2FC": [1.0, 1.0],
            "complete": [True, True],
            "regulated": ["UP", "UP"],
            # In flight the "not on/off" sentinel is "", not NaN.
            "onoff": ["", ""],
        }
    )


# ---------------------------------------------------------------------------
# 1. The Wave-0 shim contract survives the body being filled in
# ---------------------------------------------------------------------------
def test_check_returns_the_same_object(tmp_path):
    df = foldchange_frame()
    assert boundaries.check("after_foldchange", df, results_dir=tmp_path) is df


def test_check_does_not_mutate_the_frame(tmp_path):
    """Not even the dtypes.

    The stage schemas use ``coerce=True`` (the Heavy intensity columns arrive
    as int64 in flight), and a coercion that leaked back onto the caller's
    frame would silently change the dtypes ``restore_left_order`` works so hard
    to preserve.
    """
    df = foldchange_frame()
    df["Intensity 31579"] = df["Intensity 31579"].astype("int64")
    before = df.copy(deep=True)
    dtypes_before = df.dtypes.to_dict()

    boundaries.check("after_foldchange", df, results_dir=tmp_path)

    pd.testing.assert_frame_equal(df, before)
    assert df.dtypes.to_dict() == dtypes_before


def test_the_foldchange_call_sites_still_use_two_positional_args():
    """``boundaries.check("after_load", df_L)`` must keep working verbatim."""
    source = (_PKG_DIR / "foldchange.py").read_text(encoding="utf-8")
    for call in (
        'boundaries.check("after_load", df_L)',
        'boundaries.check("after_load", df_H)',
        'boundaries.check("after_merge", merged)',
        'boundaries.check("after_foldchange", df)',
    ):
        assert call in source, f"call site changed: {call}"


def test_stages_are_unchanged():
    assert boundaries.STAGES == (
        "after_load",
        "after_merge",
        "after_foldchange",
        "before_ipa_export",
    )
    assert boundaries.STRICT_STAGES | boundaries.PERMISSIVE_STAGES == set(
        boundaries.STAGES
    )
    assert not (boundaries.STRICT_STAGES & boundaries.PERMISSIVE_STAGES)


def test_every_stage_has_a_schema():
    assert set(boundaries.STAGE_SCHEMAS) == set(boundaries.STAGES)


def test_an_unknown_stage_is_an_error(tmp_path):
    """A typo'd stage name must not silently validate nothing."""
    with pytest.raises(ValueError, match="unknown boundary stage"):
        boundaries.check("after_lod", raw_sheet(), results_dir=tmp_path)


# ---------------------------------------------------------------------------
# 2. Permissive stages record; they do not raise
# ---------------------------------------------------------------------------
def test_after_load_does_not_raise_on_a_junk_accession(tmp_path):
    """The raw sheets legitimately contain the junk rows.

    Raising here would mean the pipeline cannot read its own input.
    """
    df = raw_sheet([JUNK_ACCESSION, "P19137"])
    result = boundaries.check("after_load", df, results_dir=tmp_path)

    assert result is df
    record = boundaries.records()[0]
    assert record["policy"] == "permissive"
    assert record["passed"] is False
    assert record["n_failures"] >= 1
    assert record["n_quarantinable_rows"] == 1
    assert record["routed_to_quarantine"] is True
    assert record["raised"] is False


def test_after_load_passes_on_clean_input(tmp_path):
    """Permissive must not mean 'never reports a failure'."""
    boundaries.check("after_load", raw_sheet(), results_dir=tmp_path)
    record = boundaries.records()[0]
    assert record["passed"] is True
    assert record["n_failures"] == 0
    assert record["n_quarantinable_rows"] == 0


def test_a_legitimate_protein_group_is_not_a_boundary_failure(tmp_path):
    """The 69-character group must sail through untouched."""
    boundaries.check("after_load", raw_sheet([LEGIT_GROUP, "P19137"]), results_dir=tmp_path)
    record = boundaries.records()[0]
    assert record["passed"] is True
    assert record["n_quarantinable_rows"] == 0


def test_after_merge_is_permissive_too(tmp_path):
    merged = raw_sheet([JUNK_ACCESSION, "P19137"])
    merged["_merge"] = ["left_only", "both"]
    boundaries.check("after_merge", merged, results_dir=tmp_path)
    record = boundaries.records()[0]
    assert record["policy"] == "permissive"
    assert record["raised"] is False


def test_after_merge_rejects_an_unknown_merge_indicator(tmp_path):
    """Permissive still RECORDS a failure -- it just does not raise."""
    merged = raw_sheet()
    merged["_merge"] = ["left_only", "sideways"]
    boundaries.check("after_merge", merged, results_dir=tmp_path)
    assert boundaries.records()[0]["passed"] is False


# ---------------------------------------------------------------------------
# 3. Strict stages raise
# ---------------------------------------------------------------------------
def test_after_foldchange_passes_on_a_well_formed_frame(tmp_path):
    boundaries.check("after_foldchange", foldchange_frame(), results_dir=tmp_path)
    assert boundaries.records()[0]["passed"] is True


def test_after_foldchange_raises_on_a_junk_accession(tmp_path):
    """The same value that is merely RECORDED at after_load stops the run here.

    This pair is the point of having two policies at all.
    """
    df = foldchange_frame()
    df.loc[0, "UniProt Accession Number"] = JUNK_ACCESSION
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("after_foldchange", df, results_dir=tmp_path)
    assert boundaries.records()[0]["raised"] is True


def test_after_foldchange_raises_on_an_unknown_regulated_label(tmp_path):
    df = foldchange_frame()
    df.loc[0, "regulated"] = "UPWARDS"
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("after_foldchange", df, results_dir=tmp_path)


def test_after_foldchange_raises_on_a_nan_log2fc_for_a_complete_row(tmp_path):
    """Bug 3's guarantee, enforced at the boundary."""
    df = foldchange_frame()
    df.loc[0, "log2FC"] = np.nan
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("after_foldchange", df, results_dir=tmp_path)


def test_after_foldchange_raises_when_onoff_and_regulated_disagree(tmp_path):
    """An onoff label on a row that is not ON_OFF is a real inconsistency."""
    df = foldchange_frame()
    df.loc[0, "onoff"] = "on_with_treatment"  # but regulated is still "UP"
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("after_foldchange", df, results_dir=tmp_path)


def test_before_ipa_export_raises_on_a_nan_log2fc(tmp_path):
    """A NaN reaching the IPA upload is silent data loss."""
    df = foldchange_frame()
    df.loc[0, "log2FC"] = np.nan
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("before_ipa_export", df, results_dir=tmp_path)


def test_before_ipa_export_rejects_unregulated_rows(tmp_path):
    """Only UP/DOWN are exported; NO CHANGE must never reach the writer."""
    df = foldchange_frame()
    df.loc[0, "regulated"] = "NO CHANGE"
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("before_ipa_export", df, results_dir=tmp_path)


def test_a_missing_declared_column_is_a_failure(tmp_path):
    """strict=False tolerates EXTRA columns, never missing ones."""
    df = foldchange_frame().drop(columns=["log2FC"])
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("after_foldchange", df, results_dir=tmp_path)


def test_extra_scratch_columns_are_tolerated(tmp_path):
    """In-flight frames carry scratch columns that never reach a file."""
    df = foldchange_frame()
    df["_merge"] = "both"
    df["Gene names_L"] = "Lama1"
    boundaries.check("after_foldchange", df, results_dir=tmp_path)
    assert boundaries.records()[0]["passed"] is True


def test_an_explicit_schema_overrides_the_stage_default(tmp_path):
    schema = pa.DataFrameSchema(
        {"log2FC": pa.Column(float, pa.Check.lt(0))}, strict=False, name="CUSTOM"
    )
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check(
            "after_foldchange", foldchange_frame(), schema=schema, results_dir=tmp_path
        )
    assert boundaries.records()[0]["schema"] == "CUSTOM"


# ---------------------------------------------------------------------------
# 4. The record file
# ---------------------------------------------------------------------------
def test_records_are_written_as_structured_json(tmp_path):
    boundaries.check("after_load", raw_sheet(), results_dir=tmp_path)
    boundaries.check("after_foldchange", foldchange_frame(), results_dir=tmp_path)

    payload = json.loads((tmp_path / "qc" / "qc_boundaries.json").read_text())
    assert payload["stages"] == list(boundaries.STAGES)
    assert [r["stage"] for r in payload["records"]] == [
        "after_load",
        "after_foldchange",
    ]
    assert [r["seq"] for r in payload["records"]] == [0, 1]
    for record in payload["records"]:
        assert {"stage", "policy", "schema", "rows", "columns", "passed"} <= set(record)


def test_a_raised_failure_is_recorded_before_it_raises(tmp_path):
    """The record must survive the exception, or the failure is undiagnosable."""
    df = foldchange_frame()
    df.loc[0, "regulated"] = "UPWARDS"
    with pytest.raises(boundaries.BoundaryValidationError):
        boundaries.check("after_foldchange", df, results_dir=tmp_path)

    payload = json.loads((tmp_path / "qc" / "qc_boundaries.json").read_text())
    assert len(payload["records"]) == 1
    assert payload["records"][0]["raised"] is True
    assert payload["records"][0]["failure_cases"]


def test_an_over_long_failing_value_is_truncated_in_the_record(tmp_path):
    """The 32,759-character junk accession must not be embedded whole.

    Recording it verbatim made the committed qc_boundaries.json 70 KB of one
    repeated string. The full value lives in the quarantine CSV, which the
    record points at.
    """
    boundaries.check("after_load", raw_sheet([JUNK_ACCESSION, "P19137"]), results_dir=tmp_path)
    record = boundaries.records()[0]

    cases = [c for c in record["failure_cases"] if isinstance(c.get("failure_case"), str)]
    assert cases, "expected a string failure case"
    for case in cases:
        assert len(case["failure_case"]) < len(JUNK_ACCESSION)
        assert "truncated" in case["failure_case"]
    assert record["quarantine_file"] == "qc/quarantine_accessions.csv"


def test_the_recorded_quarantine_path_is_not_absolute(tmp_path):
    """The record is committed; an absolute path bakes in one machine's layout."""
    boundaries.check("after_load", raw_sheet([JUNK_ACCESSION]), results_dir=tmp_path)
    assert not Path(boundaries.records()[0]["quarantine_file"]).is_absolute()


def test_the_record_file_is_byte_reproducible(tmp_path):
    """No wall-clock stamp: the file lives under the byte-freeze gate.

    ``qc_report.json`` carries a ``generated_at`` and therefore drifts on every
    run. That is the behaviour deliberately not copied here.
    """
    boundaries.check("after_load", raw_sheet(), results_dir=tmp_path)
    first = (tmp_path / "qc" / "qc_boundaries.json").read_bytes()

    boundaries.reset()
    boundaries.check("after_load", raw_sheet(), results_dir=tmp_path)
    second = (tmp_path / "qc" / "qc_boundaries.json").read_bytes()

    assert first == second
    assert b"generated_at" not in first


def test_check_writes_nowhere_but_the_given_results_dir(tmp_path):
    """A test run must not touch the committed results/qc artifacts."""
    committed = _RESULTS / "qc" / "qc_boundaries.json"
    before = committed.read_bytes() if committed.exists() else None

    boundaries.check("after_load", raw_sheet(), results_dir=tmp_path)

    after = committed.read_bytes() if committed.exists() else None
    assert after == before


def test_set_results_dir_redirects_the_two_argument_call(tmp_path):
    """Regression: ``--results-dir`` used to leak QC records into the repo.

    The four foldchange.py call sites pass only ``(stage, df)``, so before
    ``set_results_dir`` existed, running the pipeline with
    ``--results-dir /tmp/somewhere`` sent every output to the temp directory
    EXCEPT the boundary record, which was written into the committed
    ``proteomics_de/results/qc/``. Caught by an actual ``--results-dir`` run.
    """
    committed = _RESULTS / "qc" / "qc_boundaries.json"
    before = committed.read_bytes() if committed.exists() else None

    boundaries.set_results_dir(tmp_path)
    boundaries.check("after_load", raw_sheet())  # two positional args only

    assert (tmp_path / "qc" / "qc_boundaries.json").is_file()
    after = committed.read_bytes() if committed.exists() else None
    assert after == before, "boundary record leaked into the committed results dir"


def test_reset_restores_the_package_default_results_dir(tmp_path):
    boundaries.set_results_dir(tmp_path)
    assert boundaries.quarantine_path().parent.parent == tmp_path
    boundaries.reset()
    assert boundaries.quarantine_path() == (
        boundaries.DEFAULT_RESULTS_DIR / "qc" / boundaries.QUARANTINE_NAME
    )


def test_foldchange_configures_the_boundary_results_dir():
    source = (_PKG_DIR / "foldchange.py").read_text(encoding="utf-8")
    assert "boundaries.set_results_dir(results_dir)" in source


# ---------------------------------------------------------------------------
# 5. Quarantine (DECISIONS_LOG D11)
# ---------------------------------------------------------------------------
def test_quarantine_splits_junk_from_legitimate(tmp_path):
    df = pd.DataFrame(
        {
            "accession": [JUNK_ACCESSION, LEGIT_GROUP, "P19137", "73;74;75"],
            "gene": [None, "Sprr2d", "Lama1", None],
        }
    )
    kept, quarantined = boundaries.quarantine_junk_accessions(
        df, "accession", source="unit-test", results_dir=tmp_path
    )
    assert list(kept["accession"]) == [LEGIT_GROUP, "P19137"]
    assert list(quarantined["accession"]) == [JUNK_ACCESSION, "73;74;75"]


def test_quarantine_preserves_the_order_of_kept_rows(tmp_path):
    df = pd.DataFrame({"accession": ["P19137", JUNK_ACCESSION, "Q9JHU4", "O08528"]})
    kept, _ = boundaries.quarantine_junk_accessions(
        df, "accession", source="unit-test", results_dir=tmp_path
    )
    assert list(kept["accession"]) == ["P19137", "Q9JHU4", "O08528"]


def test_quarantine_record_keeps_the_value_recoverable(tmp_path):
    """Quarantine means set aside with a reason, not deleted."""
    df = pd.DataFrame({"accession": [JUNK_ACCESSION, "P19137"]})
    boundaries.quarantine_junk_accessions(
        df, "accession", source="unit-test", results_dir=tmp_path
    )
    written = pd.read_csv(boundaries.quarantine_path(tmp_path))

    assert list(written.columns) == boundaries.QUARANTINE_COLUMNS
    assert len(written) == 1
    assert written.loc[0, "value"] == JUNK_ACCESSION
    assert written.loc[0, "n_chars"] == len(JUNK_ACCESSION)
    assert written.loc[0, "n_tokens"] == 120
    assert "row-index list" in written.loc[0, "reason"]


def test_a_single_bare_integer_is_not_quarantined(tmp_path):
    """It is a malformed accession, and must be REPORTED as one.

    Quietly quarantining it would hide a different bug behind D11's mechanism.
    """
    df = pd.DataFrame({"accession": ["12345"]})
    kept, quarantined = boundaries.quarantine_junk_accessions(
        df, "accession", source="unit-test", results_dir=tmp_path
    )
    assert len(quarantined) == 0
    assert len(kept) == 1


def test_quarantine_is_a_no_op_on_clean_data(tmp_path):
    df = pd.DataFrame({"accession": ["P19137", LEGIT_GROUP]})
    kept, quarantined = boundaries.quarantine_junk_accessions(
        df, "accession", source="unit-test", results_dir=tmp_path
    )
    assert len(kept) == 2
    assert len(quarantined) == 0
    # Still writes a header-only record, so "no quarantine file" never has to be
    # guessed at as either "clean run" or "step did not run".
    assert list(pd.read_csv(boundaries.quarantine_path(tmp_path)).columns) == (
        boundaries.QUARANTINE_COLUMNS
    )


def test_the_committed_quarantine_matches_the_committed_outputs(frozen_counts):
    """End state: 604 shipped + 2 quarantined == the original 606."""
    scp = pd.read_csv(_RESULTS / "single_condition_proteins.csv")
    quarantined = pd.read_csv(_RESULTS / "qc" / "quarantine_accessions.csv")

    assert len(scp) == frozen_counts["single_condition_rows"] == 604
    assert len(quarantined) == 2
    assert len(scp) + len(quarantined) == 606
    assert not boundaries.find_junk_accessions(scp, "accession").any()
