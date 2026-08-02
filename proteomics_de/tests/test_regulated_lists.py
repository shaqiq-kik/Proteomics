"""Tests for :mod:`proteomics_de.export.regulated_lists`.

Three things are being defended here:

1. **The UP/DOWN split is a total, non-overlapping partition** of
   ``ipa_input_full.csv`` -- no protein lost, none duplicated, none appearing
   in both files.
2. **``fold_change`` matches ``2**log2FC``** exactly, the same formula
   ``limma_test.py`` already uses for its own ``fold_change`` column.
3. **The sort is stable across repeated builds.** The byte-freeze gate
   (DECISIONS_LOG D14) requires re-running the pipeline to reproduce identical
   bytes; an unstable sort on tied ``|log2FC|`` values would break that
   silently.
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

from proteomics_de.export import regulated_lists as rl  # noqa: E402
from proteomics_de.export.regulated_lists import (  # noqa: E402
    REGULATED_OUTPUT_COLS,
    RegulatedListsError,
    add_fold_change,
    build_regulated_lists,
    sort_by_magnitude,
    split_regulated,
    write_regulated,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def tiny_full():
    """A small ipa_input_full-shaped frame, including a near-tie for sorting."""
    return pd.DataFrame(
        {
            "UniProt Accession Number": ["P00001", "P00002", "P00003", "P00004", "P00005"],
            "Gene names": ["Aaa", "Bbb", "Ccc", "Ddd", "Eee"],
            "log2FC": [1.5, 0.6, -0.9, -2.0, 0.6],
            "regulated": ["UP", "UP", "DOWN", "DOWN", "UP"],
            "p_value": [0.01, 0.2, 0.3, 0.02, 0.4],
            "adj_p_value": [0.9, 0.95, 0.99, 0.91, 0.99],
        }
    )


@pytest.fixture(scope="module")
def full_csv(results_dir):
    path = results_dir / "ipa_input_full.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return path


@pytest.fixture(scope="module")
def committed_full(full_csv):
    return pd.read_csv(full_csv, float_precision="round_trip")


# ---------------------------------------------------------------------------
# add_fold_change
# ---------------------------------------------------------------------------


def test_add_fold_change_matches_2_pow_log2fc(tiny_full):
    out = add_fold_change(tiny_full)
    assert np.allclose(out["fold_change"], 2.0 ** tiny_full["log2FC"])


def test_add_fold_change_is_always_positive_and_finite(tiny_full):
    out = add_fold_change(tiny_full)
    assert (out["fold_change"] > 0).all()
    assert np.isfinite(out["fold_change"]).all()


def test_add_fold_change_does_not_mutate_the_input(tiny_full):
    before = tiny_full.copy()
    add_fold_change(tiny_full)
    pd.testing.assert_frame_equal(tiny_full, before)


# ---------------------------------------------------------------------------
# split_regulated
# ---------------------------------------------------------------------------


def test_split_regulated_is_a_total_partition(tiny_full):
    up, down = split_regulated(tiny_full)
    assert len(up) + len(down) == len(tiny_full)
    assert set(up["UniProt Accession Number"]) & set(down["UniProt Accession Number"]) == set()
    assert (up["regulated"] == "UP").all()
    assert (down["regulated"] == "DOWN").all()


def test_split_regulated_rejects_an_unexpected_label(tiny_full):
    bad = tiny_full.copy()
    bad.loc[0, "regulated"] = "NO CHANGE"
    with pytest.raises(RegulatedListsError, match="NO CHANGE"):
        split_regulated(bad)


# ---------------------------------------------------------------------------
# sort_by_magnitude
# ---------------------------------------------------------------------------


def test_sort_order_is_descending_by_magnitude_of_change(tiny_full):
    up, down = split_regulated(tiny_full)
    up_sorted = sort_by_magnitude(up)
    down_sorted = sort_by_magnitude(down)
    assert list(up_sorted["log2FC"].abs()) == sorted(up["log2FC"].abs(), reverse=True)
    assert list(down_sorted["log2FC"].abs()) == sorted(down["log2FC"].abs(), reverse=True)


def test_sort_is_stable_on_ties(tiny_full):
    """Bbb and Eee both have log2FC == 0.6; stable sort keeps input order."""
    up, _ = split_regulated(tiny_full)
    up_sorted = sort_by_magnitude(up)
    tied = up_sorted[up_sorted["log2FC"].abs() == 0.6]
    assert list(tied["Gene names"]) == ["Bbb", "Eee"]


# ---------------------------------------------------------------------------
# write_regulated
# ---------------------------------------------------------------------------


def test_write_regulated_emits_no_bom(tiny_full, tmp_path):
    out = tmp_path / "x.csv"
    df = add_fold_change(tiny_full)
    write_regulated(df, out)
    assert out.read_bytes()[:3] != b"\xef\xbb\xbf"


def test_write_regulated_writes_no_index_column(tiny_full, tmp_path):
    out = tmp_path / "x.csv"
    df = add_fold_change(tiny_full).iloc[::-1]  # non-trivial index
    write_regulated(df, out)
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert header.split(",")[0] == "Gene names"


def test_write_regulated_columns_match_the_declared_contract(tiny_full, tmp_path):
    out = tmp_path / "x.csv"
    df = add_fold_change(tiny_full)
    write_regulated(df, out)
    back = pd.read_csv(out, float_precision="round_trip")
    assert list(back.columns) == REGULATED_OUTPUT_COLS


# ---------------------------------------------------------------------------
# build_regulated_lists -- end to end through a results directory
# ---------------------------------------------------------------------------


def test_build_regulated_lists_is_idempotent(full_csv, tmp_path):
    """Running twice produces identical bytes -- no timestamps, no ordering drift."""
    (tmp_path / "ipa_input_full.csv").write_bytes(full_csv.read_bytes())

    first = build_regulated_lists(tmp_path)
    snapshot = {p.name: p.read_bytes() for p in first}
    build_regulated_lists(tmp_path)
    assert {p.name: p.read_bytes() for p in first} == snapshot


def test_build_regulated_lists_leaves_ipa_input_full_untouched(full_csv, tmp_path):
    (tmp_path / "ipa_input_full.csv").write_bytes(full_csv.read_bytes())
    before = (tmp_path / "ipa_input_full.csv").read_bytes()
    build_regulated_lists(tmp_path)
    assert (tmp_path / "ipa_input_full.csv").read_bytes() == before


def test_build_regulated_lists_needs_ipa_input_full_first(tmp_path):
    with pytest.raises(FileNotFoundError, match="ipa_input_full.csv"):
        build_regulated_lists(tmp_path)


def test_build_regulated_lists_raises_on_row_count_mismatch(tiny_full, tmp_path, monkeypatch):
    tiny_full.to_csv(tmp_path / "ipa_input_full.csv", index=False)
    counts = tmp_path / "counts.json"
    counts.write_text('{"n_up": 999, "n_down": 2}', encoding="utf-8")
    with pytest.raises(RegulatedListsError, match="n_up"):
        build_regulated_lists(tmp_path, counts_path=counts)


def test_committed_regulated_files_if_present(results_dir, frozen_counts):
    up_path = results_dir / "regulated_up.csv"
    down_path = results_dir / "regulated_down.csv"
    if not up_path.exists():
        pytest.skip("regulated_up.csv not generated yet")

    up = pd.read_csv(up_path, float_precision="round_trip")
    down = pd.read_csv(down_path, float_precision="round_trip")

    assert list(up.columns) == REGULATED_OUTPUT_COLS
    assert list(down.columns) == REGULATED_OUTPUT_COLS
    assert len(up) == frozen_counts["n_up"]
    assert len(down) == frozen_counts["n_down"]

    for df in (up, down):
        assert np.isfinite(df["log2FC"]).all()
        assert np.isfinite(df["fold_change"]).all()
        assert np.isfinite(df["p_value"]).all()
        assert np.isfinite(df["adj_p_value"]).all()
        assert np.allclose(df["fold_change"], 2.0 ** df["log2FC"])

    assert set(up["UniProt Accession Number"]) & set(down["UniProt Accession Number"]) == set()
    assert list(up["log2FC"].abs()) == sorted(up["log2FC"].abs(), reverse=True)
    assert list(down["log2FC"].abs()) == sorted(down["log2FC"].abs(), reverse=True)


def test_regulated_files_partition_ipa_input_full(committed_full, results_dir):
    up_path = results_dir / "regulated_up.csv"
    down_path = results_dir / "regulated_down.csv"
    if not up_path.exists():
        pytest.skip("regulated_up.csv not generated yet")

    up = pd.read_csv(up_path, float_precision="round_trip")
    down = pd.read_csv(down_path, float_precision="round_trip")
    combined = set(up["UniProt Accession Number"]) | set(down["UniProt Accession Number"])
    assert combined == set(committed_full["UniProt Accession Number"])


def test_provenance_sidecars_are_written_for_both_outputs(results_dir):
    for name in ("regulated_up.csv", "regulated_down.csv"):
        path = results_dir / name
        sidecar = results_dir / f"{name}.provenance.json"
        if not path.exists():
            pytest.skip(f"{name} not generated yet")
        assert sidecar.exists(), f"missing provenance sidecar for {name}"


# ---------------------------------------------------------------------------
# Source-level policy checks, mirroring ipa_export.py's equivalents
# ---------------------------------------------------------------------------


def test_every_read_csv_asks_for_round_trip_float_precision():
    with open(rl.__file__, encoding="utf-8") as fh:
        source = fh.read()
    calls = re.findall(r"pd\.read_csv\((?:[^()]|\([^()]*\))*\)", source)
    assert calls, "no pd.read_csv calls found -- did the module move?"
    missing = [c for c in calls if 'float_precision="round_trip"' not in c]
    assert not missing, (
        "these reads would silently perturb floats in the last ULP: " + str(missing)
    )


def test_module_hardcodes_no_dataset_specific_counts(frozen_counts):
    with open(rl.__file__, encoding="utf-8") as fh:
        text = fh.read()
    code = re.sub(r'"""(?:.|\n)*?"""', "", text)
    code = re.sub(r"#.*", "", code)
    forbidden = {
        str(frozen_counts[k])
        for k in ("ipa_input_rows", "n_up", "n_down")
    }
    found = sorted(n for n in forbidden if re.search(rf"(?<![\w.]){n}(?![\w.])", code))
    assert not found, f"dataset-specific literals in regulated_lists.py: {found}"
