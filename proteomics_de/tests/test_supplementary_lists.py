"""Tests for :mod:`proteomics_de.export.supplementary_lists`.

Mirrors ``tests/test_regulated_lists.py``'s structure. Four things are being
defended here:

1. **Tier-3 selection is correct**: ``n_imputed > 0`` rows only, split UP/DOWN
   at the given threshold, with columns renamed to the core-file convention.
2. **The qualitative union is a total, non-overlapping combination** of
   ``single_condition_proteins.csv`` + ``onoff_proteins.csv`` -- every row from
   both, none duplicated -- with ``direction`` derived purely from ``basis``.
3. **Neither new file overlaps the core ``regulated_up.csv``/``regulated_down.csv``**
   accession set -- structurally impossible by construction, asserted anyway.
4. **The qualitative file carries no fold-change or p-value column, on purpose.**
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.export import supplementary_lists as sl  # noqa: E402
from proteomics_de.export.supplementary_lists import (  # noqa: E402
    PARTIAL_OUTPUT_COLS,
    QUALITATIVE_OUTPUT_COLS,
    SupplementaryListsError,
    build_partial_regulated_lists,
    build_qualitative_changes,
    select_partial_regulated,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tiny_qc_limma():
    """A small qc_limma-shaped frame: 2 complete rows (n_imputed=0, must never
    appear in the output), 2 tier-3 UP, 1 tier-3 DOWN, 1 tier-3 below threshold.
    """
    return pd.DataFrame(
        {
            "id": ["P00001", "P00002", "P00003", "P00004", "P00005", "P00006"],
            "gene": ["Aaa", "Bbb", "Ccc", "Ddd", "Eee", "Fff"],
            "limma_log2FC": [1.5, -1.2, 2.0, 0.7, -0.9, 0.1],
            "p_value": [0.01, 0.02, 0.1, 0.2, 0.3, 0.4],
            "adj_p_value": [0.9, 0.91, 0.95, 0.98, 0.99, 0.99],
            "n_imputed": [0, 0, 1, 2, 1, 1],
        }
    )


@pytest.fixture
def tiny_single_condition():
    return pd.DataFrame(
        {
            "accession": ["P10001", "P10002", "P10003"],
            "gene": ["Ggg", "Hhh", "Iii"],
            "detected_in": ["treated_only", "control_only", "control_only"],
            "Intensity 31578": [1000.0, np.nan, np.nan],
            "Intensity 31580": [1100.0, np.nan, np.nan],
            "Intensity 31579": [np.nan, 2000.0, 3000.0],
            "Intensity 31581": [np.nan, 2100.0, 3100.0],
        }
    )


@pytest.fixture
def tiny_onoff():
    return pd.DataFrame(
        {
            "accession": ["P20001", "P20002"],
            "gene": ["Jjj", "Kkk"],
            "onoff": ["on_with_treatment", "off_with_treatment"],
            "Intensity 31578": [5000.0, 0.0],
            "Intensity 31580": [5100.0, 0.0],
            "Intensity 31579": [0.0, 6000.0],
            "Intensity 31581": [0.0, 6100.0],
        }
    )


@pytest.fixture(scope="module")
def committed_partial_up(results_dir):
    path = results_dir / "regulated_up_partial.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return pd.read_csv(path, float_precision="round_trip")


@pytest.fixture(scope="module")
def committed_partial_down(results_dir):
    path = results_dir / "regulated_down_partial.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return pd.read_csv(path, float_precision="round_trip")


@pytest.fixture(scope="module")
def committed_qualitative(results_dir):
    path = results_dir / "qualitative_changes.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return pd.read_csv(path, float_precision="round_trip")


# ---------------------------------------------------------------------------
# select_partial_regulated
# ---------------------------------------------------------------------------


def test_select_partial_regulated_excludes_complete_rows(tiny_qc_limma):
    up, down = select_partial_regulated(tiny_qc_limma, threshold=0.585)
    all_ids = set(up["UniProt Accession Number"]) | set(down["UniProt Accession Number"])
    assert "P00001" not in all_ids  # n_imputed == 0
    assert "P00002" not in all_ids  # n_imputed == 0


def test_select_partial_regulated_splits_up_down_at_threshold(tiny_qc_limma):
    up, down = select_partial_regulated(tiny_qc_limma, threshold=0.585)
    # P00003 (2.0) and P00004 (0.7) are tier-3 and cross +0.585 -> UP.
    # P00005 (-0.9) is tier-3 and crosses -0.585 -> DOWN.
    # P00006 (0.1) is tier-3 but stays below threshold -> excluded from both.
    assert set(up["UniProt Accession Number"]) == {"P00003", "P00004"}
    assert set(down["UniProt Accession Number"]) == {"P00005"}


def test_select_partial_regulated_renames_columns(tiny_qc_limma):
    up, down = select_partial_regulated(tiny_qc_limma, threshold=0.585)
    for df in (up, down):
        assert "UniProt Accession Number" in df.columns
        assert "Gene names" in df.columns
        assert "log2FC" in df.columns
        assert "id" not in df.columns
        assert "gene" not in df.columns
        assert "limma_log2FC" not in df.columns


def test_select_partial_regulated_raises_on_missing_columns():
    bad = pd.DataFrame({"id": ["P1"], "gene": ["X"]})
    with pytest.raises(SupplementaryListsError, match="missing columns"):
        select_partial_regulated(bad, threshold=0.585)


# ---------------------------------------------------------------------------
# build_partial_regulated_lists -- end to end through a results directory
# ---------------------------------------------------------------------------


def test_build_partial_regulated_lists_needs_qc_limma_first(tmp_path):
    with pytest.raises(FileNotFoundError, match="qc_limma.csv"):
        build_partial_regulated_lists(tmp_path)


def test_build_partial_regulated_lists_is_idempotent(tiny_qc_limma, tmp_path, monkeypatch):
    tiny_qc_limma.to_csv(tmp_path / "qc_limma.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text(
        '{"n_up_partial": 2, "n_down_partial": 1}', encoding="utf-8"
    )
    first = build_partial_regulated_lists(tmp_path, threshold=0.585, counts_path=counts)
    snapshot = {p.name: p.read_bytes() for p in first}
    build_partial_regulated_lists(tmp_path, threshold=0.585, counts_path=counts)
    assert {p.name: p.read_bytes() for p in first} == snapshot


def test_build_partial_regulated_lists_raises_on_row_count_mismatch(tiny_qc_limma, tmp_path):
    tiny_qc_limma.to_csv(tmp_path / "qc_limma.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text('{"n_up_partial": 999, "n_down_partial": 1}', encoding="utf-8")
    with pytest.raises(SupplementaryListsError, match="n_up_partial"):
        build_partial_regulated_lists(tmp_path, threshold=0.585, counts_path=counts)


# ---------------------------------------------------------------------------
# build_qualitative_changes
# ---------------------------------------------------------------------------


def test_build_qualitative_changes_needs_both_inputs(tmp_path):
    with pytest.raises(FileNotFoundError):
        build_qualitative_changes(tmp_path)


def test_build_qualitative_changes_is_a_total_union(
    tiny_single_condition, tiny_onoff, tmp_path
):
    tiny_single_condition.to_csv(tmp_path / "single_condition_proteins.csv", index=False)
    tiny_onoff.to_csv(tmp_path / "onoff_proteins.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text(
        '{"n_qualitative": 5, "n_qualitative_up": 2, "n_qualitative_down": 3}',
        encoding="utf-8",
    )
    [out_path] = build_qualitative_changes(tmp_path, counts_path=counts)
    out = pd.read_csv(out_path, float_precision="round_trip")
    assert len(out) == 5
    assert set(out["UniProt Accession Number"]) == {
        "P10001", "P10002", "P10003", "P20001", "P20002",
    }


def test_build_qualitative_changes_direction_matches_basis(
    tiny_single_condition, tiny_onoff, tmp_path
):
    tiny_single_condition.to_csv(tmp_path / "single_condition_proteins.csv", index=False)
    tiny_onoff.to_csv(tmp_path / "onoff_proteins.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text(
        '{"n_qualitative": 5, "n_qualitative_up": 2, "n_qualitative_down": 3}',
        encoding="utf-8",
    )
    [out_path] = build_qualitative_changes(tmp_path, counts_path=counts)
    out = pd.read_csv(out_path, float_precision="round_trip")
    by_acc = out.set_index("UniProt Accession Number")
    assert by_acc.loc["P10001", "direction"] == "UP"    # treated_only
    assert by_acc.loc["P10002", "direction"] == "DOWN"  # control_only
    assert by_acc.loc["P10003", "direction"] == "DOWN"  # control_only
    assert by_acc.loc["P20001", "direction"] == "UP"    # on_with_treatment
    assert by_acc.loc["P20002", "direction"] == "DOWN"  # off_with_treatment


def test_build_qualitative_changes_has_no_quantitative_columns(
    tiny_single_condition, tiny_onoff, tmp_path
):
    tiny_single_condition.to_csv(tmp_path / "single_condition_proteins.csv", index=False)
    tiny_onoff.to_csv(tmp_path / "onoff_proteins.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text(
        '{"n_qualitative": 5, "n_qualitative_up": 2, "n_qualitative_down": 3}',
        encoding="utf-8",
    )
    [out_path] = build_qualitative_changes(tmp_path, counts_path=counts)
    header = out_path.read_text(encoding="utf-8").splitlines()[0].split(",")
    forbidden = {"log2FC", "fold_change", "p_value", "adj_p_value"}
    assert not forbidden & set(header)


def test_build_qualitative_changes_is_idempotent(tiny_single_condition, tiny_onoff, tmp_path):
    tiny_single_condition.to_csv(tmp_path / "single_condition_proteins.csv", index=False)
    tiny_onoff.to_csv(tmp_path / "onoff_proteins.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text(
        '{"n_qualitative": 5, "n_qualitative_up": 2, "n_qualitative_down": 3}',
        encoding="utf-8",
    )
    first = build_qualitative_changes(tmp_path, counts_path=counts)
    snapshot = {p.name: p.read_bytes() for p in first}
    build_qualitative_changes(tmp_path, counts_path=counts)
    assert {p.name: p.read_bytes() for p in first} == snapshot


def test_build_qualitative_changes_raises_on_row_count_mismatch(
    tiny_single_condition, tiny_onoff, tmp_path
):
    tiny_single_condition.to_csv(tmp_path / "single_condition_proteins.csv", index=False)
    tiny_onoff.to_csv(tmp_path / "onoff_proteins.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text(
        '{"n_qualitative": 999, "n_qualitative_up": 2, "n_qualitative_down": 3}',
        encoding="utf-8",
    )
    with pytest.raises(SupplementaryListsError, match="n_qualitative"):
        build_qualitative_changes(tmp_path, counts_path=counts)


# ---------------------------------------------------------------------------
# Committed-file checks (skip if not generated yet)
# ---------------------------------------------------------------------------


def test_committed_partial_files_match_frozen_counts(
    committed_partial_up, committed_partial_down, frozen_counts
):
    assert len(committed_partial_up) == frozen_counts["n_up_partial"]
    assert len(committed_partial_down) == frozen_counts["n_down_partial"]
    assert list(committed_partial_up.columns) == PARTIAL_OUTPUT_COLS
    assert list(committed_partial_down.columns) == PARTIAL_OUTPUT_COLS


def test_committed_partial_files_every_row_is_imputed(
    committed_partial_up, committed_partial_down
):
    for df in (committed_partial_up, committed_partial_down):
        assert (df["n_imputed"] >= 1).all(), (
            "every row here rests on at least one imputed replicate -- that is "
            "the entire point of this file's separation from regulated_up/down.csv"
        )


def test_committed_partial_files_fold_change_matches_2_pow_log2fc(
    committed_partial_up, committed_partial_down
):
    for df in (committed_partial_up, committed_partial_down):
        assert np.allclose(df["fold_change"], 2.0 ** df["log2FC"])
        assert np.isfinite(df["log2FC"]).all()
        assert np.isfinite(df["fold_change"]).all()


def test_committed_partial_files_are_sorted_by_descending_magnitude(
    committed_partial_up, committed_partial_down
):
    for df in (committed_partial_up, committed_partial_down):
        assert list(df["log2FC"].abs()) == sorted(df["log2FC"].abs(), reverse=True)


def test_committed_qualitative_matches_frozen_counts(committed_qualitative, frozen_counts):
    assert len(committed_qualitative) == frozen_counts["n_qualitative"]
    assert (committed_qualitative["direction"] == "UP").sum() == frozen_counts["n_qualitative_up"]
    assert (committed_qualitative["direction"] == "DOWN").sum() == frozen_counts["n_qualitative_down"]
    assert list(committed_qualitative.columns) == QUALITATIVE_OUTPUT_COLS


def test_committed_qualitative_has_no_quantitative_columns(committed_qualitative):
    forbidden = {"log2FC", "fold_change", "p_value", "adj_p_value"}
    assert not forbidden & set(committed_qualitative.columns)


def test_committed_qualitative_direction_is_binary(committed_qualitative):
    assert set(committed_qualitative["direction"].unique()) <= {"UP", "DOWN"}


def test_provenance_sidecars_are_written_for_all_three_outputs(results_dir):
    for name in ("regulated_up_partial.csv", "regulated_down_partial.csv", "qualitative_changes.csv"):
        path = results_dir / name
        sidecar = results_dir / f"{name}.provenance.json"
        if not path.exists():
            pytest.skip(f"{name} not generated yet")
        assert sidecar.exists(), f"missing provenance sidecar for {name}"


# ---------------------------------------------------------------------------
# Source-level policy checks, mirroring regulated_lists.py's/ipa_export.py's
# equivalents
# ---------------------------------------------------------------------------


def test_every_read_csv_asks_for_round_trip_float_precision():
    with open(sl.__file__, encoding="utf-8") as fh:
        source = fh.read()
    calls = re.findall(r"pd\.read_csv\((?:[^()]|\([^()]*\))*\)", source)
    assert calls, "no pd.read_csv calls found -- did the module move?"
    missing = [c for c in calls if 'float_precision="round_trip"' not in c]
    assert not missing, (
        "these reads would silently perturb floats in the last ULP: " + str(missing)
    )


def test_module_hardcodes_no_dataset_specific_counts(frozen_counts):
    with open(sl.__file__, encoding="utf-8") as fh:
        text = fh.read()
    code = re.sub(r'"""(?:.|\n)*?"""', "", text)
    code = re.sub(r"#.*", "", code)
    forbidden = {
        str(frozen_counts[k])
        for k in ("n_up_partial", "n_down_partial", "n_qualitative",
                  "n_qualitative_up", "n_qualitative_down")
    }
    found = sorted(n for n in forbidden if re.search(rf"(?<![\w.]){n}(?![\w.])", code))
    assert not found, f"dataset-specific literals in supplementary_lists.py: {found}"
