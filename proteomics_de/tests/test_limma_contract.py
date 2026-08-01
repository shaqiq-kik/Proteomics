"""Schema contract for the committed ``results/qc_limma.csv``.

Shape and statistical sanity only -- these assertions must hold for *any* valid
limma run, not just today's. Byte-identity is the freeze manifest's job
(``tools/status.py --check``); duplicating it here would make the suite fail for
the wrong reason whenever the numbers are legitimately regenerated.

The intensity columns are asserted against ``config/sample_sheet.tsv`` rather
than against literals, so a sheet change surfaces here as a clear failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.config import design  # noqa: E402

QC_LIMMA = _REPO_ROOT / "proteomics_de" / "results" / "qc_limma.csv"

pytestmark = pytest.mark.skipif(
    not QC_LIMMA.exists(), reason="results/qc_limma.csv not present"
)


@pytest.fixture(scope="module")
def qc() -> pd.DataFrame:
    return pd.read_csv(QC_LIMMA)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------
def test_columns_exact_and_ordered(qc):
    expected = (
        ["id", "gene"]
        + design.sample_columns()
        + ["limma_log2FC", "p_value", "adj_p_value", "significant", "regulated"]
    )
    assert list(qc.columns) == expected


def test_row_count(qc, frozen_counts):
    assert len(qc) == frozen_counts["qc_limma_rows"] == 1938


def test_ids_are_unique_and_present(qc):
    assert qc["id"].notna().all()
    assert qc["id"].is_unique


def test_intensity_columns_come_from_the_sample_sheet(qc):
    for col in design.sample_columns():
        assert col in qc.columns


# ---------------------------------------------------------------------------
# Statistical sanity
# ---------------------------------------------------------------------------
def test_pvalues_finite_and_in_unit_interval(qc):
    for col in ("p_value", "adj_p_value"):
        vals = pd.to_numeric(qc[col], errors="coerce")
        assert vals.notna().all(), f"{col} has non-numeric entries"
        assert np.isfinite(vals).all(), f"{col} has non-finite entries"
        assert vals.between(0.0, 1.0).all(), f"{col} outside [0, 1]"


def test_adjusted_p_is_never_below_raw_p(qc):
    """BH only ever raises a p-value."""
    raw = pd.to_numeric(qc["p_value"])
    adj = pd.to_numeric(qc["adj_p_value"])
    assert (adj >= raw - 1e-12).all()


def test_adjusted_p_is_bh_monotone(qc):
    adj = pd.to_numeric(qc["adj_p_value"])
    ordered = adj[pd.to_numeric(qc["p_value"]).sort_values().index].to_numpy()
    assert np.all(np.diff(ordered) >= -1e-12)


def test_logfc_is_finite(qc):
    vals = pd.to_numeric(qc["limma_log2FC"], errors="coerce")
    assert vals.notna().all() and np.isfinite(vals).all()


# ---------------------------------------------------------------------------
# Derived-column consistency
# ---------------------------------------------------------------------------
def test_significant_flag_matches_the_fdr_threshold(qc):
    import limma_test

    adj = pd.to_numeric(qc["adj_p_value"])
    expected = adj < limma_test.ADJ_P_THRESHOLD
    assert qc["significant"].astype(bool).tolist() == expected.tolist()


def test_regulated_values_are_from_the_known_vocabulary(qc):
    assert set(qc["regulated"].dropna()) <= {"UP", "DOWN", "NO CHANGE", "ON_OFF"}


def test_no_onoff_proteins_were_tested(qc):
    """Fork 2: ON_OFF proteins are excluded -- testing them invents a condition."""
    assert "ON_OFF" not in set(qc["regulated"].dropna())
