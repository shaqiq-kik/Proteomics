"""Contract tests for ``etl/build_matrix.py`` (research1.md Section 6 item 7).

The point of these is the *forward path*: the two contract files must be shaped
by ``config/sample_sheet.tsv`` and by nothing else, so that adding biological
replicates is an edit to the sheet rather than a code change. The synthetic
6-sample test at the bottom is what actually proves that; the rest lock today's
2x2 layout.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.config import design  # noqa: E402
from proteomics_de.etl import build_matrix  # noqa: E402

# Canonical (control-first) order, D7-corrected: 31579/31581 are the controls
# (Vehicle), 31578/31580 the treated (Testosterone). See DECISIONS_LOG D7.
EXPECTED_SAMPLE_COLUMNS = [
    "Intensity 31579",
    "Intensity 31581",
    "Intensity 31578",
    "Intensity 31580",
]


def _fake_eligible(sample_cols, n=5):
    """A minimal stand-in for the eligible rows of foldchange_all.csv."""
    data = {
        build_matrix.ACCESSION_COL: [f"P{i:05d}" for i in range(n)],
        build_matrix.GENE_COL: [f"Gene{i}" for i in range(n)],
    }
    for j, col in enumerate(sample_cols):
        data[col] = [str(1000.0 * (i + 1) + j) for i in range(n)]
    return pd.DataFrame(data)


# ---------------------------------------------------------------------------
# intensity_matrix.tsv
# ---------------------------------------------------------------------------
def test_intensity_matrix_columns_and_order(tmp_path):
    eligible = _fake_eligible(EXPECTED_SAMPLE_COLUMNS)
    written = build_matrix.build(eligible, tmp_path)

    mat = pd.read_csv(written["intensity_matrix"], sep="\t")
    assert list(mat.columns) == ["accession", "gene"] + EXPECTED_SAMPLE_COLUMNS
    # The sample columns are the sample sheet's, in the sample sheet's order.
    assert list(mat.columns)[2:] == design.sample_columns()


def test_intensity_matrix_is_tab_separated_one_row_per_protein(tmp_path):
    eligible = _fake_eligible(EXPECTED_SAMPLE_COLUMNS, n=7)
    written = build_matrix.build(eligible, tmp_path)

    text = written["intensity_matrix"].read_text(encoding="utf-8")
    header = text.splitlines()[0]
    assert "\t" in header and "," not in header
    assert len(text.splitlines()) == 1 + 7  # header + one row per protein


def test_intensity_matrix_lands_under_de_subdir(tmp_path):
    written = build_matrix.build(_fake_eligible(EXPECTED_SAMPLE_COLUMNS), tmp_path)
    assert written["intensity_matrix"] == tmp_path / "de" / "intensity_matrix.tsv"
    assert written["design"] == tmp_path / "de" / "design.tsv"


def test_row_order_is_preserved_verbatim(tmp_path):
    eligible = _fake_eligible(EXPECTED_SAMPLE_COLUMNS, n=6)
    written = build_matrix.build(eligible, tmp_path)

    mat = pd.read_csv(written["intensity_matrix"], sep="\t")
    assert mat["accession"].tolist() == eligible[build_matrix.ACCESSION_COL].tolist()


def test_nonpositive_and_blank_intensities_become_empty_cells(tmp_path):
    """0 means 'below detection limit' (MNAR), not 'measured as zero'."""
    eligible = _fake_eligible(EXPECTED_SAMPLE_COLUMNS, n=3)
    eligible.loc[0, EXPECTED_SAMPLE_COLUMNS[0]] = "0"
    eligible.loc[1, EXPECTED_SAMPLE_COLUMNS[1]] = ""
    eligible.loc[2, EXPECTED_SAMPLE_COLUMNS[2]] = "-5"

    written = build_matrix.build(eligible, tmp_path)
    lines = written["intensity_matrix"].read_text(encoding="utf-8").splitlines()

    assert lines[1].split("\t")[2] == ""  # the 0
    assert lines[2].split("\t")[3] == ""  # the blank
    assert lines[3].split("\t")[4] == ""  # the negative


# ---------------------------------------------------------------------------
# design.tsv
# ---------------------------------------------------------------------------
def test_design_tsv_schema_and_groups(tmp_path):
    written = build_matrix.build(_fake_eligible(EXPECTED_SAMPLE_COLUMNS), tmp_path)

    dsn = pd.read_csv(written["design"], sep="\t", dtype=str)
    assert list(dsn.columns) == ["sample", "group"]
    assert set(dsn["group"]) <= {"control", "treated"}
    assert dsn["group"].tolist() == design.group_vector()
    assert len(dsn) == design.n_samples()


def test_design_tsv_lists_controls_before_treated(tmp_path):
    """Canonical order is load-bearing: it fixes the sign of every logFC."""
    written = build_matrix.build(_fake_eligible(EXPECTED_SAMPLE_COLUMNS), tmp_path)
    groups = pd.read_csv(written["design"], sep="\t", dtype=str)["group"].tolist()
    assert groups.index("treated") > max(
        i for i, g in enumerate(groups) if g == "control"
    )


# ---------------------------------------------------------------------------
# Fail-loud
# ---------------------------------------------------------------------------
def test_empty_eligible_df_raises(tmp_path):
    empty = _fake_eligible(EXPECTED_SAMPLE_COLUMNS).iloc[0:0]
    with pytest.raises(ValueError, match="zero rows"):
        build_matrix.build(empty, tmp_path)


def test_missing_sample_column_raises(tmp_path):
    eligible = _fake_eligible(EXPECTED_SAMPLE_COLUMNS).drop(
        columns=[EXPECTED_SAMPLE_COLUMNS[-1]]
    )
    with pytest.raises(ValueError, match="missing column"):
        build_matrix.build(eligible, tmp_path)


# ---------------------------------------------------------------------------
# Forward path: more replicates must be a sheet edit, not a code edit
# ---------------------------------------------------------------------------
def _six_sample_sheet():
    rows = [
        ("A1", "control", "Intensity A1", 1),
        ("A2", "control", "Intensity A2", 2),
        ("A3", "control", "Intensity A3", 3),
        ("B1", "treated", "Intensity B1", 1),
        ("B2", "treated", "Intensity B2", 2),
        ("B3", "treated", "Intensity B3", 3),
    ]
    return pd.DataFrame(rows, columns=["sample", "group", "channel", "replicate"])


def test_forward_path_six_samples(tmp_path):
    sheet = _six_sample_sheet()
    cols = design.sample_columns(sheet)
    assert len(cols) == 6

    eligible = _fake_eligible(cols, n=4)
    written = build_matrix.build(eligible, tmp_path, sheet=sheet)

    mat = pd.read_csv(written["intensity_matrix"], sep="\t")
    assert list(mat.columns) == ["accession", "gene"] + cols
    assert len(mat) == 4

    dsn = pd.read_csv(written["design"], sep="\t", dtype=str)
    assert len(dsn) == 6
    assert dsn["group"].tolist() == ["control"] * 3 + ["treated"] * 3
    assert set(dsn["group"]) <= {"control", "treated"}


def test_forward_path_ignores_sheet_row_order(tmp_path):
    """Re-sorting the TSV must not permute the design matrix."""
    shuffled = _six_sample_sheet().iloc[[4, 0, 3, 2, 5, 1]].reset_index(drop=True)
    cols = design.sample_columns(shuffled)

    eligible = _fake_eligible(cols, n=3)
    written = build_matrix.build(eligible, tmp_path, sheet=shuffled)

    mat = pd.read_csv(written["intensity_matrix"], sep="\t")
    assert list(mat.columns)[2:] == [
        "Intensity A1", "Intensity A2", "Intensity A3",
        "Intensity B1", "Intensity B2", "Intensity B3",
    ]
