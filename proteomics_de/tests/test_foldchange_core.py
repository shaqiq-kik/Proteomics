"""Executable versions of research1.md's Bug 1-4b arguments.

Every test here is a hand-built three-to-five-row DataFrame -- no workbook, no
CSV, no results directory. That is the whole point of extracting
``etl/foldchange_core.py``: before it, the fold-change logic lived inside
``foldchange.py``'s ``if __name__ == "__main__"`` block and could only be
exercised by running the real 2,315-protein pipeline and squinting at the output.

The bug numbers match research1.md's audit of the legacy
``proteomics_analysis.py``. Where a test asserts what the *legacy* code would
have produced, it recomputes that arithmetic inline rather than importing it, so
the comparison is visible in the test body instead of being taken on faith.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.config.constants import LOG2_THRESHOLD  # noqa: E402
from proteomics_de.etl import foldchange_core as core  # noqa: E402

# The 2-channel SILAC layout. Sheet position and experimental condition are kept
# as separate facts here, because they came apart at DECISIONS_LOG D7 -- and
# conflating them is precisely the bug D7 exposed in
# ``build_single_condition_frame`` (``left_only`` was hardcoded to
# ``"control_only"``, mislabelling all 606 rescued proteins the moment the
# assignment was corrected).

#: Which channels live in which sheet of the workbook. A fact about the file;
#: D7 did not change it, and neither does relabelling the conditions.
SHEET_L_COLS = ["Intensity 31578", "Intensity 31580"]
SHEET_H_COLS = ["Intensity 31579", "Intensity 31581"]
#: Acquisition order, as ``foldchange.py``'s ``PHYSICAL_COLS``/``out_cols``.
PHYSICAL_COLS = SHEET_L_COLS + SHEET_H_COLS

#: Which of those is which condition, D7-corrected to match
#: ``config/sample_sheet.tsv``: the lab's Pilot Project labels 31579/31581 =
#: Vehicle (control) and 31578/31580 = Testosterone (treated), so sheet H holds
#: the controls and sheet L the treated.
CONTROL_COLS = SHEET_H_COLS
TREATED_COLS = SHEET_L_COLS
INTENSITY_COLS = CONTROL_COLS + TREATED_COLS
BASE_COLS = ["UniProt Accession Number", "Gene names"] + INTENSITY_COLS


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------
def make_frame(rows):
    """Build a `both`-shaped frame from ``(acc, gene, c1, c2, t1, t2)`` tuples."""
    return pd.DataFrame(list(rows), columns=BASE_COLS)


def run_core(df, threshold=LOG2_THRESHOLD):
    """Run the full in-memory chain exactly as ``foldchange.main`` wires it."""
    df = core.mark_complete(df, INTENSITY_COLS)
    mask = df["complete"]
    df, ratio_cols = core.compute_ratios(df, CONTROL_COLS, TREATED_COLS, mask)
    df, _log_cols = core.compute_log2fc(df, ratio_cols)
    df = core.classify_regulated(df, mask, threshold)
    df, _control_off, _treated_off = core.detect_onoff(df, CONTROL_COLS, TREATED_COLS)
    return df


def classify_log2fc(values, threshold=LOG2_THRESHOLD):
    """Label a list of log2FC values, treating every row as complete."""
    df = pd.DataFrame({"log2FC": list(values)})
    df["complete"] = True
    core.classify_regulated(df, df["complete"], threshold)
    return df["regulated"].tolist()


def legacy_log2fc(ratios):
    """The legacy rule: arithmetic mean of the raw ratios, logged once.

    ``proteomics_analysis.py``:
        mean_fold_change = (fc1 + fc2) / 2
        log2_fold_change = np.log2(mean_fold_change)
    """
    return float(np.log2(sum(ratios) / len(ratios)))


def legacy_regulated(ratios):
    """The legacy asymmetric call: UP at mean ratio >= 1.5, DOWN at <= 0.5."""
    mean_ratio = sum(ratios) / len(ratios)
    if mean_ratio >= 1.5:
        return "UP"
    if mean_ratio <= 0.5:
        return "DOWN"
    return "NO CHANGE"


# ---------------------------------------------------------------------------
# Bug 1 -- geometric, not arithmetic, mean of the replicate ratios
# ---------------------------------------------------------------------------
def test_bug1_up_twofold_and_down_twofold_cancel_to_exactly_zero():
    """research1.md line 46's worked example: 2.0 and 0.5 must give log2FC 0.0.

    A protein that doubles in one run and halves in the other has not changed.
    Logging each ratio first makes that exact -- +1.0 and -1.0 average to 0.0
    with no floating-point residue at all.
    """
    df = run_core(make_frame([("P1", "GeneA", 100.0, 100.0, 200.0, 50.0)]))

    assert df.loc[0, "ratio_rep1"] == 2.0
    assert df.loc[0, "ratio_rep2"] == 0.5
    assert df.loc[0, "log2_rep1"] == 1.0
    assert df.loc[0, "log2_rep2"] == -1.0
    # Exactly zero, not approximately -- this is the claim being made.
    assert df.loc[0, "log2FC"] == 0.0
    assert df.loc[0, "regulated"] == "NO CHANGE"


def test_bug1_legacy_arithmetic_mean_then_log_would_have_said_plus_0_3219():
    """The same protein, scored the legacy way, drifts upward by +0.32 log2.

    ``mean(2.0, 0.5) = 1.25 -> log2 = +0.3219``. The bias is not a rounding
    artefact: it is systematic, one-directional, and it is why the legacy
    volcano leaned up.
    """
    ratios = (2.0, 0.5)

    assert legacy_log2fc(ratios) == pytest.approx(0.3219, abs=1e-4)

    df = run_core(make_frame([("P1", "GeneA", 100.0, 100.0, 200.0, 50.0)]))
    corrected = float(df.loc[0, "log2FC"])

    assert corrected == 0.0
    assert legacy_log2fc(ratios) - corrected == pytest.approx(0.3219, abs=1e-4)


@pytest.mark.parametrize("ratios", [(2.0, 0.5), (4.0, 0.25), (3.0, 1.0 / 3.0), (10.0, 0.1)])
def test_bug1_the_legacy_bias_is_always_upward(ratios):
    """For any reciprocal pair the true answer is 0 and the legacy answer is > 0.

    The arithmetic mean of a number and its reciprocal is >= 1 (AM-GM), so the
    legacy estimator can only ever overstate up-regulation, never understate it.
    """
    c1 = c2 = 100.0
    df = run_core(make_frame([("P1", "G", c1, c2, c1 * ratios[0], c2 * ratios[1])]))

    assert float(df.loc[0, "log2FC"]) == pytest.approx(0.0, abs=1e-12)
    assert legacy_log2fc(ratios) > 0.0


# ---------------------------------------------------------------------------
# Bug 2 -- symmetric cutoffs in log2 space
# ---------------------------------------------------------------------------
def test_bug2_boundary_is_symmetric_and_inclusive():
    """+0.585 -> UP, -0.585 -> DOWN, +/-0.584 -> NO CHANGE.

    The comparison is ``>=`` / ``<=``, so the threshold value itself is inside
    the regulated set, and the two sides sit at the same distance from zero.
    """
    assert classify_log2fc([0.585, -0.585, 0.584, -0.584]) == [
        "UP", "DOWN", "NO CHANGE", "NO CHANGE",
    ]


@pytest.mark.parametrize("magnitude", [0.0, 0.1, 0.584, 0.585, 0.5851, 1.0, 3.7])
def test_bug2_up_and_down_are_mirror_images(magnitude):
    """``classify(+x) == UP`` if and only if ``classify(-x) == DOWN``.

    This is the actual content of the Bug 2 fix: not the value 0.585, but that
    the same |log2FC| produces the same verdict in both directions. The legacy
    rule failed this for every magnitude between log2(1.5) and 1.0.
    """
    up, down = classify_log2fc([magnitude, -magnitude])
    assert (up == "UP") == (down == "DOWN")


def test_bug2_a_two_thirds_fall_is_now_called_down():
    """A protein that falls to ~0.666x control is DOWN; the legacy rule missed it.

    Legacy required a full 2x fall (``ratio <= 0.5``) to call DOWN while calling
    UP at a 1.5x rise. Everything in between -- the entire 0.5 to 0.667 band --
    was invisible on the down side only.
    """
    ratio = 0.666
    df = run_core(make_frame([("P1", "G", 1000.0, 1000.0, 666.0, 666.0)]))

    assert float(df.loc[0, "ratio_rep1"]) == pytest.approx(ratio)
    assert df.loc[0, "regulated"] == "DOWN"
    assert legacy_regulated((ratio, ratio)) == "NO CHANGE"


def test_bug2_threshold_is_a_rounded_log2_of_1_point_5():
    """A ratio of *exactly* 1.5 (or 1/1.5) is NO CHANGE, by 7.2e-5 of log2.

    ``LOG2_THRESHOLD`` is 0.585, the 3-decimal rounding of
    ``log2(1.5) = 0.5849625...``, so the cutoff is a hair stricter than the
    nominal "1.5-fold" it is named for. This is recorded, not fixed: changing the
    constant would move real proteins across the UP/DOWN boundary and rewrite
    every committed output. The property that matters -- symmetry -- is
    unaffected, because the same rounded value bounds both sides.
    """
    assert LOG2_THRESHOLD == 0.585
    assert LOG2_THRESHOLD > np.log2(1.5)

    exactly_1_5 = float(np.log2(1.5))
    assert classify_log2fc([exactly_1_5, -exactly_1_5]) == ["NO CHANGE", "NO CHANGE"]


def test_bug2_incomplete_rows_are_never_regulated():
    """The UP/DOWN masks are ANDed with ``complete``; a NaN log2FC stays NO CHANGE."""
    df = run_core(make_frame([
        ("P1", "G1", 100.0, 100.0, 400.0, 400.0),   # complete, clearly UP
        ("P2", "G2", 0.0, 100.0, 400.0, 400.0),     # incomplete (a zero control)
    ]))

    assert df.loc[0, "regulated"] == "UP"
    assert df.loc[1, "regulated"] == "NO CHANGE"
    assert bool(df.loc[1, "complete"]) is False


# ---------------------------------------------------------------------------
# Bug 3 -- no inf, ever
# ---------------------------------------------------------------------------
def test_bug3_zero_denominator_gives_nan_not_inf():
    """A zero control intensity must never reach the division.

    The legacy script divided first and filtered afterwards, so ``200 / 0``
    became ``inf`` and ``np.log2(0)`` became ``-inf``, and both were persisted to
    the saved artifact where any downstream consumer that did not re-filter would
    choke on them.
    """
    df = run_core(make_frame([("P1", "G", 0.0, 100.0, 200.0, 200.0)]))

    assert bool(df.loc[0, "complete"]) is False
    for col in ("ratio_rep1", "ratio_rep2", "log2_rep1", "log2_rep2", "log2FC"):
        value = df.loc[0, col]
        assert pd.isna(value), f"{col} = {value!r}, expected NaN"
        assert not np.isinf(value), f"{col} is infinite"

    # What the legacy ordering (divide first, filter later) would have produced.
    with np.errstate(divide="ignore"):
        assert np.isinf(np.divide(200.0, 0.0))
        assert np.isneginf(np.log2(0.0))


def test_bug3_missing_and_zero_intensities_are_both_incomplete():
    """NaN and 0 are treated identically: absent is absent."""
    df = run_core(make_frame([
        ("P1", "G1", 100.0, 100.0, 100.0, 100.0),      # all present
        ("P2", "G2", np.nan, 100.0, 100.0, 100.0),     # NaN control
        ("P3", "G3", 100.0, 100.0, 0.0, 100.0),        # zero treated
    ]))

    assert df["complete"].tolist() == [True, False, False]


def test_bug3_no_non_finite_value_survives_in_any_complete_row():
    """The invariant ``foldchange.py`` asserts on the real data, on a fixture."""
    df = run_core(make_frame([
        ("P1", "G1", 100.0, 100.0, 400.0, 400.0),
        ("P2", "G2", 100.0, 100.0, 25.0, 25.0),
        ("P3", "G3", 0.0, 0.0, 50.0, 50.0),        # incomplete, excluded
    ]))
    complete_log2fc = df.loc[df["complete"], "log2FC"]

    assert len(complete_log2fc) == 2
    assert not np.isinf(complete_log2fc).any()
    assert not complete_log2fc.isnull().any()


# ---------------------------------------------------------------------------
# Bug 4 -- the outer merge, and what it must not disturb
# ---------------------------------------------------------------------------
@pytest.fixture
def two_sheets():
    """L and H sheets that overlap on three accessions and differ on two.

    L's row order (C, A, B) is deliberately *not* sorted, because the merge sorts
    by key and the pipeline has to put the order back.
    """
    df_L = pd.DataFrame(
        [
            ("C", "GeneC", 300.0, 300.0),
            ("A", "GeneA", 100.0, 100.0),
            ("B", "GeneB", 200.0, 200.0),
            ("L_ONLY", "GeneL", 900.0, 900.0),
        ],
        columns=["UniProt Accession Number", "Gene names"] + SHEET_L_COLS,
    )
    df_H = pd.DataFrame(
        [
            ("A", "GeneA", 150, 150),
            ("B", "GeneB", 400, 400),
            ("C", "GeneC", 300, 300),
            ("H_ONLY", "GeneH", 700, 700),
        ],
        columns=["UniProt Accession Number", "Gene names"] + SHEET_H_COLS,
    )
    return df_L, df_H


def test_bug4_single_condition_rows_are_rescued_with_the_right_label(two_sheets):
    """Proteins seen in exactly one sheet reach their own file, correctly labelled.

    The inner join discarded these silently. ``detected_in`` records which side
    they came from, which is the only thing that distinguishes "off in treated"
    from "off in control".

    D7-corrected: sheet L holds the Testosterone channels, so a protein seen only
    in sheet L is ``treated_only``. (On the real data that is the 367/239 split
    in ``single_condition_proteins.csv``, which D7 inverted.)
    """
    df_L, df_H = two_sheets
    merged = core.merge_with_indicator(df_L, df_H)
    _both, single = core.split_both_single(merged)
    single = core.build_single_condition_frame(
        single, left_condition="treated", right_condition="control",
    )

    labels = dict(zip(single["accession"], single["detected_in"]))
    assert labels == {"L_ONLY": "treated_only", "H_ONLY": "control_only"}
    assert set(single["gene"]) == {"GeneL", "GeneH"}
    # The absent sheet's intensities stay blank -- the blank is the finding.
    l_only = single.set_index("accession").loc["L_ONLY"]
    assert pd.isna(l_only[SHEET_H_COLS]).all()
    h_only = single.set_index("accession").loc["H_ONLY"]
    assert pd.isna(h_only[SHEET_L_COLS]).all()


def test_bug4_detected_in_follows_the_conditions_not_the_sheet_position(two_sheets):
    """Same merge, opposite condition arguments: every label must flip.

    ``left_only -> "control_only"`` used to be hardcoded, so the label tracked
    SHEET POSITION and silently contradicted the sample sheet once D7 corrected
    the assignment. This is the assertion that catches a regression to that.
    """
    df_L, df_H = two_sheets
    merged = core.merge_with_indicator(df_L, df_H)
    _both, single = core.split_both_single(merged)
    single = core.build_single_condition_frame(
        single, left_condition="control", right_condition="treated",
    )

    assert dict(zip(single["accession"], single["detected_in"])) == {
        "L_ONLY": "control_only", "H_ONLY": "treated_only",
    }


def test_bug4_both_plus_single_covers_the_merge_exactly(two_sheets):
    """``len(both) + len(single) == len(merged)`` -- the split loses nothing."""
    df_L, df_H = two_sheets
    merged = core.merge_with_indicator(df_L, df_H)
    both, single = core.split_both_single(merged)

    assert len(both) + len(single) == len(merged)
    assert len(both) == 3
    assert len(single) == 2
    # Disjoint as well as exhaustive.
    assert set(both["UniProt Accession Number"]) & set(single["UniProt Accession Number"]) == set()


def test_bug4_matched_rows_are_restored_to_left_sheet_order(two_sheets):
    """``both`` comes back in L-sheet order, not the merge's sorted-key order.

    ``pd.merge`` returns keys sorted (A, B, C); the L sheet has them as C, A, B.
    Without the restore, every committed output would be re-sorted -- same
    numbers, different file.
    """
    df_L, df_H = two_sheets
    merged = core.merge_with_indicator(df_L, df_H)
    both, _single = core.split_both_single(merged)

    assert both["UniProt Accession Number"].tolist() == ["A", "B", "C"]  # merge order

    df = core.restore_left_order(
        both, df_L, df_H,
        out_cols=["UniProt Accession Number", "Gene names"] + PHYSICAL_COLS,
        dtype_cols=SHEET_H_COLS,  # sheet-H channels: dtype, not condition
    )
    assert df["UniProt Accession Number"].tolist() == ["C", "A", "B"]
    assert df.index.tolist() == [0, 1, 2]  # reset, not merely reordered


def test_bug4_heavy_dtypes_are_restored_after_the_outer_join(two_sheets):
    """The outer join's NaNs coerce the Heavy columns to float; restore undoes it.

    Nothing numeric changes -- but ``400`` and ``400.0`` are different bytes in a
    CSV, so this is the difference between the byte-freeze passing and failing.
    """
    df_L, df_H = two_sheets
    assert all(df_H[c].dtype == np.int64 for c in SHEET_H_COLS)

    merged = core.merge_with_indicator(df_L, df_H)
    both, _single = core.split_both_single(merged)
    assert all(both[c].dtype == np.float64 for c in SHEET_H_COLS)  # coerced by the join

    df = core.restore_left_order(
        both, df_L, df_H,
        out_cols=["UniProt Accession Number", "Gene names"] + PHYSICAL_COLS,
        dtype_cols=SHEET_H_COLS,
    )
    assert all(df[c].dtype == df_H[c].dtype for c in SHEET_H_COLS)


def test_bug4_gene_names_are_coalesced_left_first():
    """A protein with a symbol in only one sheet keeps it after the merge."""
    df_L = pd.DataFrame(
        [("A", "FromL", 1.0, 1.0), ("B", np.nan, 1.0, 1.0)],
        columns=["UniProt Accession Number", "Gene names"] + SHEET_L_COLS,
    )
    df_H = pd.DataFrame(
        [("A", "FromH", 1.0, 1.0), ("B", "FromH", 1.0, 1.0)],
        columns=["UniProt Accession Number", "Gene names"] + SHEET_H_COLS,
    )
    merged = core.merge_with_indicator(df_L, df_H)
    genes = dict(zip(merged["UniProt Accession Number"], merged["Gene names"]))

    assert genes == {"A": "FromL", "B": "FromH"}


# ---------------------------------------------------------------------------
# Bug 4b -- the on/off truth table
# ---------------------------------------------------------------------------
@pytest.fixture
def onoff_truth_table():
    """One row per cell of the control-absent x treated-absent truth table."""
    return run_core(make_frame([
        # control absent, treated present -> switched ON by the treatment
        ("ON", "GeneOn", 0.0, np.nan, 500.0, 400.0),
        # treated absent, control present -> switched OFF by the treatment
        ("OFF", "GeneOff", 500.0, 400.0, 0.0, np.nan),
        # absent on BOTH sides -> an empty row, not a signal
        ("EMPTY", "GeneEmpty", 0.0, np.nan, np.nan, 0.0),
        # present on both sides -> ordinary fold change, untouched
        ("BOTH", "GeneBoth", 100.0, 100.0, 100.0, 100.0),
    ]))


def test_bug4b_control_absent_treated_present_is_on_with_treatment(onoff_truth_table):
    row = onoff_truth_table.set_index("UniProt Accession Number").loc["ON"]
    assert row["onoff"] == "on_with_treatment"
    assert row["regulated"] == "ON_OFF"


def test_bug4b_treated_absent_control_present_is_off_with_treatment(onoff_truth_table):
    row = onoff_truth_table.set_index("UniProt Accession Number").loc["OFF"]
    assert row["onoff"] == "off_with_treatment"
    assert row["regulated"] == "ON_OFF"


def test_bug4b_absent_on_both_sides_stays_no_change(onoff_truth_table):
    """The load-bearing negative case.

    A protein with nothing on either side is an empty row. Calling it ON_OFF
    would manufacture a biological claim out of missing data -- and because
    ON_OFF rows are exported to IPA, that claim would leave the pipeline.
    """
    row = onoff_truth_table.set_index("UniProt Accession Number").loc["EMPTY"]
    assert row["onoff"] == ""
    assert row["regulated"] == "NO CHANGE"
    assert row["regulated"] != "ON_OFF"


def test_bug4b_a_fully_present_protein_is_untouched(onoff_truth_table):
    row = onoff_truth_table.set_index("UniProt Accession Number").loc["BOTH"]
    assert row["onoff"] == ""
    assert row["regulated"] == "NO CHANGE"
    assert float(row["log2FC"]) == 0.0


def test_bug4b_every_onoff_row_has_nan_log2fc(onoff_truth_table):
    """ON_OFF proteins are incomplete by construction, so no ratio exists.

    Emitting a number here would be inventing one: with one side entirely absent
    there is no denominator, and a "fold change" against nothing is not a
    quantity. ``foldchange.py`` asserts the same thing on the real data.
    """
    onoff_rows = onoff_truth_table[onoff_truth_table["regulated"] == "ON_OFF"]

    assert len(onoff_rows) == 2
    assert onoff_rows["log2FC"].isnull().all()
    assert (~onoff_rows["complete"]).all()


def test_bug4b_onoff_only_ever_relabels_out_of_no_change():
    """UP and DOWN calls are never overwritten by the on/off pass.

    A row that is complete cannot be one-sided, so the two masks are disjoint --
    but the assertion is worth pinning, because the on/off pass writes into the
    same ``regulated`` column the classifier just filled.
    """
    df = run_core(make_frame([
        ("UPP", "G1", 100.0, 100.0, 400.0, 400.0),
        ("DWN", "G2", 400.0, 400.0, 100.0, 100.0),
        ("ON", "G3", 0.0, 0.0, 400.0, 400.0),
    ]))
    labels = dict(zip(df["UniProt Accession Number"], df["regulated"]))

    assert labels == {"UPP": "UP", "DWN": "DOWN", "ON": "ON_OFF"}


# ---------------------------------------------------------------------------
# Output frames and the summary
# ---------------------------------------------------------------------------
def test_build_ipa_frame_keeps_complete_and_regulated_only():
    """IPA gets complete UP/DOWN rows -- not NO CHANGE, and not incomplete ON_OFF."""
    df = run_core(make_frame([
        ("UPP", "G1", 100.0, 100.0, 400.0, 400.0),    # complete, UP    -> kept
        ("DWN", "G2", 400.0, 400.0, 100.0, 100.0),    # complete, DOWN  -> kept
        ("NCH", "G3", 100.0, 100.0, 100.0, 100.0),    # complete, NC    -> dropped
        ("ON", "G4", 0.0, 0.0, 400.0, 400.0),         # ON_OFF, incomplete -> dropped
    ]))
    ipa = core.build_ipa_frame(df)

    assert ipa["UniProt Accession Number"].tolist() == ["UPP", "DWN"]


def test_build_onoff_frame_aliases_accession_and_gene(onoff_truth_table):
    onoff = core.build_onoff_frame(onoff_truth_table)

    assert onoff["accession"].tolist() == ["ON", "OFF"]
    assert onoff["gene"].tolist() == ["GeneOn", "GeneOff"]
    assert onoff["onoff"].tolist() == ["on_with_treatment", "off_with_treatment"]


def test_summarize_counts_partition_the_frame(onoff_truth_table):
    counts = core.summarize(onoff_truth_table)

    assert counts["n_rows"] == 4
    assert counts["n_onoff"] == counts["n_on"] + counts["n_off"] == 2
    assert (
        counts["n_up"] + counts["n_down"] + counts["n_nochange"] + counts["n_onoff"]
        == counts["n_rows"]
    )
    # Plain ints, so the values are JSON-serialisable and print like the
    # numpy.int64 they replaced.
    assert all(isinstance(v, int) for v in counts.values())


# ---------------------------------------------------------------------------
# Frozen-count expectations
# ---------------------------------------------------------------------------
def test_frozen_counts_are_loaded_from_json_not_from_source():
    """The four numbers ``foldchange.py`` used to hardcode live in one file now."""
    counts = core.load_frozen_counts()

    for key in ("n_up", "n_down", "n_nochange_plus_onoff", "ipa_input_rows"):
        assert key in counts, f"frozen_counts.json is missing {key!r}"
        assert isinstance(counts[key], int)


def test_frozen_counts_path_is_resolved_from_file_not_cwd():
    assert core.FROZEN_COUNTS_PATH.is_absolute()
    assert core.FROZEN_COUNTS_PATH.exists()


@pytest.mark.parametrize(
    "env, expected",
    [({}, True), ({"PDE_EXPECT_BASELINE": "1"}, True), ({"PDE_EXPECT_BASELINE": "0"}, False)],
)
def test_baseline_checks_are_on_unless_explicitly_disabled(env, expected):
    """Default ON. The escape hatch exists for a new dataset, not for a red build."""
    assert core.baseline_checks_enabled(env) is expected
