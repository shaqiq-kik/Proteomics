"""Schema contract for the committed ``results/qc_limma.csv``.

Shape and statistical sanity only -- these assertions must hold for *any* valid
limma run, not just today's. Byte-identity is the freeze manifest's job
(``tools/status.py --check``); duplicating it here would make the suite fail for
the wrong reason whenever the numbers are legitimately regenerated.

The intensity columns are cross-checked against ``config/sample_sheet.tsv``, so a
sheet change surfaces here as a clear failure. Note they appear in ACQUISITION
order, not the sheet's canonical control-first order -- see
``ACQUISITION_ORDER`` below.

Also covered here, from the two decisions applied in this wave:

* **D9** -- ``eBayes(trend=TRUE, robust=TRUE)`` is the default, and the vanilla
  flavour is preserved as ``results/qc_limma_vanilla.csv``. The pair is asserted
  to differ in p-values and NOT in logFC: that is the whole justification for
  changing the default, so it is checked rather than believed.
* **D10** -- ``n_imputed, AveExpr, t, B`` are restored as APPENDED columns, and
  ``results/de/limma_results.tsv`` carries research1.md line 169's schema with
  ``fold_change`` and ``direction``. ``n_imputed`` is validated against the raw
  intensities in ``foldchange_all.csv``, not merely range-checked -- a count
  that is plausible but wrong would be worse than no column at all.
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

_RESULTS = _REPO_ROOT / "proteomics_de" / "results"
QC_LIMMA = _RESULTS / "qc_limma.csv"
QC_LIMMA_VANILLA = _RESULTS / "qc_limma_vanilla.csv"
LIMMA_RESULTS = _RESULTS / "de" / "limma_results.tsv"
FOLDCHANGE_ALL = _RESULTS / "foldchange_all.csv"

#: The order ``limma_test.py`` hands the intensity matrix to R, and therefore the
#: order these columns appear in ``qc_limma.csv``: PHYSICAL (acquisition) order,
#: i.e. the order the channels sit in the workbook -- deliberately NOT the sample
#: sheet's canonical control-first order. MinProb imputation is stochastic per
#: column, so a layout that depends on the condition labels would hand different
#: random draws to different samples the moment the labels changed; keeping the
#: matrix in acquisition order made the DECISIONS_LOG D7 relabelling a pure sign
#: inversion with every p-value invariant. Same four columns as the sheet, one
#: fixed permutation of them (asserted below).
ACQUISITION_ORDER = [
    "Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581",
]

#: The 11 columns qc_limma.csv shipped with, in order. DECISIONS_LOG D10 appends
#: to this list; it may never reorder, rename or drop any of it. Downstream
#: (``viz/*``, ``gated/``, ``enrich/gsea.py``) selects by name, but a stable
#: prefix is also what every committed review of this file describes.
ORIGINAL_COLUMNS = (
    ["id", "gene"]
    + ACQUISITION_ORDER
    + ["limma_log2FC", "p_value", "adj_p_value", "significant", "regulated"]
)

#: Appended by D10 (research1.md line 169), in this order.
D10_COLUMNS = ["n_imputed", "AveExpr", "t", "B"]

#: ``results/de/limma_results.tsv`` -- research1.md line 169, verbatim.
LIMMA_RESULTS_COLUMNS = [
    "accession", "gene", "logFC", "fold_change", "AveExpr", "t",
    "P.Value", "adj.P.Val", "B", "n_imputed", "direction",
]

pytestmark = pytest.mark.skipif(
    not QC_LIMMA.exists(), reason="results/qc_limma.csv not present"
)


@pytest.fixture(scope="module")
def qc() -> pd.DataFrame:
    return pd.read_csv(QC_LIMMA)


@pytest.fixture(scope="module")
def qc_vanilla() -> pd.DataFrame:
    if not QC_LIMMA_VANILLA.exists():
        pytest.skip("results/qc_limma_vanilla.csv not present")
    return pd.read_csv(QC_LIMMA_VANILLA)


@pytest.fixture(scope="module")
def limma_results() -> pd.DataFrame:
    if not LIMMA_RESULTS.exists():
        pytest.skip("results/de/limma_results.tsv not present")
    return pd.read_csv(LIMMA_RESULTS, sep="\t")


@pytest.fixture(scope="module")
def eligible() -> pd.DataFrame:
    """The rows of ``foldchange_all.csv`` that limma actually tested.

    Read as strings with blanks preserved, exactly as ``limma_test.py`` does, so
    the missing-value truth below is derived from the RAW intensities rather
    than from anything the limma path already computed.
    """
    if not FOLDCHANGE_ALL.exists():
        pytest.skip("results/foldchange_all.csv not present")
    fc = pd.read_csv(FOLDCHANGE_ALL, dtype=str, keep_default_na=False)
    return fc[fc["onoff"].str.strip() == ""].reset_index(drop=True)


def _true_missing_counts(eligible: pd.DataFrame) -> pd.Series:
    """Per row, how many intensity cells are NOT a positive number.

    That is the pipeline's missing-value rule (``etl.build_matrix``):
    blank, non-numeric and ``<= 0`` all mean "below the detection limit", and
    all three become NA before MinProb sees the matrix.
    """
    num = eligible[ACQUISITION_ORDER].apply(pd.to_numeric, errors="coerce")
    return (~(num > 0)).sum(axis=1)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------
def test_columns_exact_and_ordered(qc):
    assert list(qc.columns) == ORIGINAL_COLUMNS + D10_COLUMNS


def test_original_columns_are_an_unchanged_prefix(qc):
    """D10 APPENDS. It must not reorder, rename or drop what was already there.

    Stated separately from the exact-columns test so that a failure says which
    of the two rules broke: adding a column at the end is a schema extension,
    moving one is a breaking change for every downstream consumer.
    """
    assert list(qc.columns)[: len(ORIGINAL_COLUMNS)] == ORIGINAL_COLUMNS


def test_acquisition_order_is_a_permutation_of_the_sample_sheet(qc):
    """Acquisition order may differ from canonical order, but not in CONTENT.

    Without this, a channel added to / removed from the sheet would slip past
    ``test_columns_exact_and_ordered``'s literal list unnoticed.
    """
    assert sorted(ACQUISITION_ORDER) == sorted(design.sample_columns())


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


# ---------------------------------------------------------------------------
# D10 -- the restored limma columns
# ---------------------------------------------------------------------------
def test_n_imputed_is_a_count_within_the_sample_count(qc):
    n = pd.to_numeric(qc["n_imputed"], errors="coerce")
    assert n.notna().all(), "n_imputed has non-numeric entries"
    assert (n == n.round()).all(), "n_imputed is not integral"
    assert n.between(0, len(ACQUISITION_ORDER)).all(), (
        f"n_imputed outside [0, {len(ACQUISITION_ORDER)}]"
    )


def test_n_imputed_equals_the_true_missing_count(qc, eligible):
    """The column must count REAL missingness, not merely look plausible.

    This is the whole point of D10's ``n_imputed``: MinProb is stochastic and
    this study is n=2, so a consumer's only defence against treating an invented
    value as a measurement is this count being right. Derived here straight from
    ``foldchange_all.csv``'s raw intensity cells -- a row with k non-positive /
    blank intensities must report ``n_imputed == k``.
    """
    assert qc["id"].tolist() == eligible["UniProt Accession Number"].tolist(), (
        "qc_limma.csv and the eligible rows of foldchange_all.csv are no longer "
        "row-aligned; the comparison below would be meaningless"
    )
    truth = _true_missing_counts(eligible)
    got = pd.to_numeric(qc["n_imputed"]).astype(int)

    mismatched = truth.index[truth.to_numpy() != got.to_numpy()]
    assert len(mismatched) == 0, (
        f"{len(mismatched)} rows disagree, e.g. "
        f"{[(qc['id'][i], int(got[i]), int(truth[i])) for i in mismatched[:5]]} "
        "(id, n_imputed, true missing count)"
    )


def test_n_imputed_is_not_trivially_constant(qc):
    """Guard against a column that passes every check by always being 0."""
    n = pd.to_numeric(qc["n_imputed"])
    assert n.nunique() > 1, "n_imputed carries no information (all rows equal)"
    assert (n > 0).any(), "no protein has an imputed value -- suspicious for MNAR data"


def test_avexpr_t_and_b_are_finite(qc):
    for col in ("AveExpr", "t", "B"):
        vals = pd.to_numeric(qc[col], errors="coerce")
        assert vals.notna().all(), f"{col} has non-numeric entries"
        assert np.isfinite(vals).all(), f"{col} has non-finite entries"


def test_t_statistic_agrees_in_sign_with_logfc(qc):
    """t and logFC are the same contrast; disagreeing signs would be a bug."""
    logfc = pd.to_numeric(qc["limma_log2FC"])
    tstat = pd.to_numeric(qc["t"])
    disagree = (np.sign(logfc) != np.sign(tstat)) & (logfc != 0)
    assert not disagree.any(), f"{int(disagree.sum())} rows have sign(t) != sign(logFC)"


# ---------------------------------------------------------------------------
# D10 -- results/de/limma_results.tsv, research1.md's contract schema
# ---------------------------------------------------------------------------
def test_limma_results_schema_is_the_documented_one(limma_results):
    assert list(limma_results.columns) == LIMMA_RESULTS_COLUMNS


def test_limma_results_row_count(limma_results, frozen_counts):
    assert len(limma_results) == frozen_counts["qc_limma_rows"]


def test_limma_results_fold_change_is_two_to_the_logfc(limma_results):
    logfc = pd.to_numeric(limma_results["logFC"]).to_numpy()
    fold = pd.to_numeric(limma_results["fold_change"]).to_numpy()
    assert np.allclose(fold, np.power(2.0, logfc), rtol=1e-12, atol=1e-12), (
        "fold_change is not 2**logFC; max |rel err| = "
        f"{np.max(np.abs(fold - 2.0 ** logfc) / np.abs(2.0 ** logfc)):.3e}"
    )


def test_limma_results_direction_follows_the_documented_rule(limma_results):
    """research1.md lines 140-141, verbatim.

    UP/DOWN require BOTH ``adj.P.Val < 0.05`` AND ``|logFC| >= 0.585``; anything
    else is NS. The significance half is not optional -- dropping it would turn
    this into an effect-size-only call and manufacture 715 "regulated" hits out
    of a study where nothing survives correction (DECISIONS_LOG D2).
    """
    import limma_test

    logfc = pd.to_numeric(limma_results["logFC"])
    adj = pd.to_numeric(limma_results["adj.P.Val"])
    sig = adj < limma_test.ADJ_P_THRESHOLD
    expected = np.where(
        sig & (logfc >= limma_test.LOG2_THRESHOLD), "UP",
        np.where(sig & (logfc <= -limma_test.LOG2_THRESHOLD), "DOWN", "NS"),
    )
    assert limma_results["direction"].tolist() == list(expected)
    assert set(limma_results["direction"]) <= {"UP", "DOWN", "NS"}


def test_limma_results_agrees_with_qc_limma(qc, limma_results):
    """The contract file is a VIEW of the same run, not a second computation."""
    assert limma_results["accession"].tolist() == qc["id"].tolist()

    # qc_limma.csv rounds logFC to 6 dp (it has since Bug 7 shipped); the
    # contract file keeps full precision because ``fold_change = 2**logFC`` is
    # derived from it. Assert the EXACT relationship rather than picking a
    # tolerance loose enough to hide a real disagreement.
    assert (
        pd.to_numeric(limma_results["logFC"]).round(6).tolist()
        == pd.to_numeric(qc["limma_log2FC"]).tolist()
    ), "limma_results.tsv logFC is not qc_limma.csv's limma_log2FC at 6 dp"

    for contract_col, qc_col in (
        ("P.Value", "p_value"), ("adj.P.Val", "adj_p_value"),
        ("AveExpr", "AveExpr"), ("t", "t"), ("B", "B"),
        ("n_imputed", "n_imputed"),
    ):
        assert np.allclose(
            pd.to_numeric(limma_results[contract_col]),
            pd.to_numeric(qc[qc_col]),
            rtol=0, atol=1e-9,
        ), f"{contract_col} disagrees with qc_limma.csv's {qc_col}"


# ---------------------------------------------------------------------------
# D9 -- trend/robust is the default; vanilla is preserved and comparable
# ---------------------------------------------------------------------------
def test_vanilla_companion_has_the_same_schema(qc, qc_vanilla):
    assert list(qc_vanilla.columns) == list(qc.columns)
    assert len(qc_vanilla) == len(qc)
    assert qc_vanilla["id"].tolist() == qc["id"].tolist()


def test_the_two_flavours_have_bit_identical_logfc(qc, qc_vanilla):
    """THE assertion behind D9.

    ``eBayes`` moderates the residual variances; it never refits the linear
    model. Switching to ``trend=TRUE, robust=TRUE`` therefore cannot move a
    single fold change -- only the moderated t and the p-values derived from it.
    If this fails, the change is NOT confined to the variance model and every
    downstream fold-change conclusion (volcano, MA plot, IPA list, GSEA ranking)
    is back in question.
    """
    assert np.array_equal(
        pd.to_numeric(qc["limma_log2FC"]).to_numpy(),
        pd.to_numeric(qc_vanilla["limma_log2FC"]).to_numpy(),
    ), (
        "logFC moved between the vanilla and trend/robust flavours; max |diff| = "
        f"{np.abs(pd.to_numeric(qc['limma_log2FC']) - pd.to_numeric(qc_vanilla['limma_log2FC'])).max():.3e}"
    )


def test_the_two_flavours_share_every_non_pvalue_column(qc, qc_vanilla):
    """Only p_value and adj_p_value may differ -- nothing else in the file."""
    for col in ORIGINAL_COLUMNS + ["n_imputed"]:
        if col in ("p_value", "adj_p_value"):
            continue
        assert qc[col].equals(qc_vanilla[col]), (
            f"'{col}' differs between the flavours; eBayes must not touch it"
        )


def test_trend_robust_actually_changed_the_pvalues(qc, qc_vanilla):
    """The control for the test above: the switch must DO something.

    Without this, a default that silently fell back to vanilla would pass every
    other assertion here.
    """
    trend_p = pd.to_numeric(qc["p_value"]).to_numpy()
    vanilla_p = pd.to_numeric(qc_vanilla["p_value"]).to_numpy()
    assert not np.array_equal(trend_p, vanilla_p), (
        "p-values are identical across flavours -- trend/robust is not in force"
    )
    assert pd.to_numeric(qc["adj_p_value"]).min() < pd.to_numeric(
        qc_vanilla["adj_p_value"]
    ).min(), "trend/robust did not sharpen the minimum adjusted p-value"


def test_neither_flavour_finds_anything_significant(qc, qc_vanilla, frozen_counts):
    """DECISIONS_LOG D2/D9: the headline result survives the change of default.

    0/1938 at FDR<0.05 either way. D9 was accepted precisely BECAUSE it does not
    manufacture a result; asserting it here keeps that honest.
    """
    expected = frozen_counts["n_significant_fdr05"]
    for name, frame in (("trend_robust", qc), ("vanilla", qc_vanilla)):
        n_sig = int((pd.to_numeric(frame["adj_p_value"]) < 0.05).sum())
        assert n_sig == expected == 0, f"{name}: {n_sig} significant, expected 0"
        assert not frame["significant"].astype(bool).any()
