"""Non-regression lock for the config-driven design layer.

These assertions are the contract that lets the rest of the refactor replace
hardcoded literals with calls into ``config.design``. If any of them fails, the
sample sheet and the frozen scripts have diverged, and the correct response is to
investigate the sheet -- NOT to relax the assertion.

The literals asserted below are copied from the frozen scripts:
``foldchange.py:26`` (intensity columns), ``limma_test.R`` (group vector) and
``limma_test.py:56`` (handoff names).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.config import design  # noqa: E402

# The literals as they appear in the frozen scripts.
EXPECTED_SAMPLE_COLUMNS = [
    "Intensity 31578",
    "Intensity 31580",
    "Intensity 31579",
    "Intensity 31581",
]
EXPECTED_GROUP_VECTOR = ["control", "control", "treated", "treated"]
EXPECTED_HANDOFF_NAMES = ["ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581"]


# ---------------------------------------------------------------------------
# The committed 2x2 SILAC design
# ---------------------------------------------------------------------------
def test_sample_columns_exact_order():
    assert design.sample_columns() == EXPECTED_SAMPLE_COLUMNS


def test_group_vector_exact():
    assert design.group_vector() == EXPECTED_GROUP_VECTOR


def test_handoff_names_exact():
    assert design.handoff_names() == EXPECTED_HANDOFF_NAMES


def test_control_and_treated_columns():
    assert design.control_columns() == ["Intensity 31578", "Intensity 31580"]
    assert design.treated_columns() == ["Intensity 31579", "Intensity 31581"]
    # control block first, then treated -- the order limma's design matrix assumes
    assert design.control_columns() + design.treated_columns() == EXPECTED_SAMPLE_COLUMNS


def test_counts():
    assert design.n_samples() == 4
    assert design.n_groups() == 2
    assert design.replicates_per_group() == 2


def test_sample_ids_align_with_columns():
    ids = design.sample_ids()
    assert ids == ["31578", "31580", "31579", "31581"]
    for sample_id, column in zip(ids, design.sample_columns()):
        assert column.endswith(sample_id)


def test_group_names_control_first():
    assert design.group_names() == ["control", "treated"]


def test_default_path_is_file_relative_not_cwd(tmp_path, monkeypatch):
    """read_sample_sheet() must not depend on the process's cwd."""
    monkeypatch.chdir(tmp_path)
    assert design.sample_columns() == EXPECTED_SAMPLE_COLUMNS


def test_row_order_in_the_tsv_does_not_change_the_design(tmp_path):
    """Canonical order is derived from `group`, not from file row order."""
    shuffled = tmp_path / "shuffled.tsv"
    shuffled.write_text(
        "sample\tgroup\tchannel\treplicate\n"
        "31581\ttreated\tIntensity 31581\t2\n"
        "31580\tcontrol\tIntensity 31580\t2\n"
        "31579\ttreated\tIntensity 31579\t1\n"
        "31578\tcontrol\tIntensity 31578\t1\n",
        encoding="utf-8",
    )
    assert design.sample_columns(shuffled) == EXPECTED_SAMPLE_COLUMNS
    assert design.group_vector(shuffled) == EXPECTED_GROUP_VECTOR
    assert design.handoff_names(shuffled) == EXPECTED_HANDOFF_NAMES


# ---------------------------------------------------------------------------
# design.tsv contract (research1.md Section 1)
# ---------------------------------------------------------------------------
def test_write_design_tsv_round_trips(tmp_path):
    out = tmp_path / "design.tsv"
    design.write_design_tsv(out)

    text = out.read_text(encoding="utf-8")
    lines = text.splitlines()

    assert lines[0] == "sample\tgroup"
    assert len(lines) == 1 + design.n_samples() == 5
    assert text.endswith("\n")

    rows = [line.split("\t") for line in lines[1:]]
    assert all(len(r) == 2 for r in rows), "design.tsv must be TAB-separated, 2 columns"

    # Row order matches sample_columns() order.
    assert [r[0] for r in rows] == design.sample_ids()
    assert [r[1] for r in rows] == EXPECTED_GROUP_VECTOR
    for (sample_id, _group), column in zip(rows, design.sample_columns()):
        assert column.endswith(sample_id)


# ---------------------------------------------------------------------------
# Forward path: the mechanism must generalise to biological replicates
# ---------------------------------------------------------------------------
@pytest.fixture
def six_sample_sheet(tmp_path):
    """3 control + 3 treated -- what the sheet looks like after bio replicates."""
    path = tmp_path / "sample_sheet_6.tsv"
    rows = ["sample\tgroup\tchannel\treplicate"]
    for i, sid in enumerate(["A1", "A2", "A3"], start=1):
        rows.append(f"{sid}\tcontrol\tIntensity {sid}\t{i}")
    for i, sid in enumerate(["B1", "B2", "B3"], start=1):
        rows.append(f"{sid}\ttreated\tIntensity {sid}\t{i}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_forward_path_six_samples(six_sample_sheet):
    sheet = six_sample_sheet
    assert design.n_samples(sheet) == 6
    assert design.n_groups(sheet) == 2
    assert design.replicates_per_group(sheet) == 3
    assert design.group_vector(sheet) == ["control"] * 3 + ["treated"] * 3
    assert design.sample_columns(sheet) == [
        "Intensity A1", "Intensity A2", "Intensity A3",
        "Intensity B1", "Intensity B2", "Intensity B3",
    ]
    assert len(design.sample_columns(sheet)) == 6
    assert design.handoff_names(sheet) == [
        "ctrl_A1", "ctrl_A2", "ctrl_A3", "trt_B1", "trt_B2", "trt_B3",
    ]


def test_forward_path_design_tsv(six_sample_sheet, tmp_path):
    out = tmp_path / "design6.tsv"
    design.write_design_tsv(out, six_sample_sheet)
    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "sample\tgroup"
    assert len(lines) == 7
    assert [line.split("\t")[1] for line in lines[1:]] == ["control"] * 3 + ["treated"] * 3


def test_unbalanced_design_is_reported_not_averaged(tmp_path):
    path = tmp_path / "unbalanced.tsv"
    path.write_text(
        "sample\tgroup\tchannel\treplicate\n"
        "A1\tcontrol\tIntensity A1\t1\n"
        "A2\tcontrol\tIntensity A2\t2\n"
        "B1\ttreated\tIntensity B1\t1\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="unbalanced design"):
        design.replicates_per_group(path)


# ---------------------------------------------------------------------------
# The drift guard
# ---------------------------------------------------------------------------
def test_assert_matches_passes_on_the_real_columns():
    design.assert_matches(EXPECTED_SAMPLE_COLUMNS)  # must not raise


def test_assert_matches_raises_on_a_wrong_list():
    with pytest.raises(ValueError) as exc:
        design.assert_matches(["Intensity 31578", "Intensity 99999"], stage="fake_stage.py")
    msg = str(exc.value)
    assert "fake_stage.py" in msg
    assert "2-channel-SILAC-specific" in msg
    assert "REGENERATED" in msg


def test_assert_matches_raises_on_wrong_order():
    reordered = [
        "Intensity 31579", "Intensity 31581", "Intensity 31578", "Intensity 31580",
    ]
    with pytest.raises(ValueError, match="different order"):
        design.assert_matches(reordered, stage="fake_stage.py")


def test_assert_matches_names_the_calling_module_by_default():
    with pytest.raises(ValueError) as exc:
        design.assert_matches(["nope"])
    assert "test_design.py" in str(exc.value)


# ---------------------------------------------------------------------------
# Sheet hygiene
# ---------------------------------------------------------------------------
def test_missing_sheet_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        design.read_sample_sheet(tmp_path / "nope.tsv")


def test_duplicate_channel_is_rejected(tmp_path):
    path = tmp_path / "dupe.tsv"
    path.write_text(
        "sample\tgroup\tchannel\treplicate\n"
        "A1\tcontrol\tIntensity A1\t1\n"
        "A2\tcontrol\tIntensity A1\t2\n"
        "B1\ttreated\tIntensity B1\t1\n"
        "B2\ttreated\tIntensity B2\t2\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate channel"):
        design.read_sample_sheet(path)


def test_missing_required_column_is_rejected(tmp_path):
    path = tmp_path / "nogroup.tsv"
    path.write_text(
        "sample\tchannel\treplicate\nA1\tIntensity A1\t1\n", encoding="utf-8"
    )
    with pytest.raises(ValueError, match="missing required column"):
        design.read_sample_sheet(path)
