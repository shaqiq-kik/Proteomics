"""Cross-file invariants: the committed outputs must be mutually *consistent*.

NOT A SECOND FREEZE TEST
------------------------
``tests/test_freeze.py`` already proves byte-identity against
``tests/expected/outputs.sha256``. This module proves something byte-identity
cannot: that the thirteen artifacts under ``results/`` still agree *with each
other*.

The difference matters because of how the freeze gate is maintained. Its
manifest is re-baselined by hand whenever a change is intentional. A
regeneration that was coherent but **wrong** — a re-inverted contrast, an IPA
export filtered on the wrong mask, a limma run over the wrong row set — produces
a new set of internally-inconsistent files, and a re-baselined manifest would
call them all correct. Everything below is derived from one file and checked
against another, so no amount of re-baselining can make it pass.

RULES OF THIS FILE
------------------
1. **No dataset number is ever inlined.** Counts come from
   ``tests/expected/frozen_counts.json``, the same file the pipeline sources
   itself asserts against. A literal ``1948`` here would just be a ninth copy of
   the number this repo already spent a wave consolidating.
2. **Nothing here is a tolerance to be widened.** Each assertion states a
   structural fact. If one fails, a file is wrong — that is a finding to
   investigate, never a test to relax. The header-only checks in particular are
   a CONTRACT (see :func:`test_ipa_input_significant_is_header_only`).
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

# The four labels ``regulated`` may take. Duplicated from
# ``etl/foldchange_core.REGULATED_CLASSES`` on purpose: a test that imports its
# expectation from the code under test cannot detect that code changing.
REGULATED_CLASSES = {"UP", "DOWN", "NO CHANGE", "ON_OFF"}

#: The IPA upload carries only the two directional classes. ON_OFF proteins are
#: excluded because they are incomplete by construction (one condition is
#: entirely absent), so they never have a numeric fold change to upload.
IPA_CLASSES = {"UP", "DOWN"}

#: Floor for corr(pipeline log2FC, limma logFC). Not a tuning knob -- the two
#: quantities are the same contrast computed by two independent code paths (a
#: pandas mean-of-log2-ratios and an R linear model on a MinProb-imputed log2
#: matrix), so anything below this is a wiring fault, not noise. Observed:
#: 0.9999999999999573.
LOGFC_CORR_FLOOR = 0.9999


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def _read(path) -> pd.DataFrame:
    """Read a committed CSV with the maximal-precision float parser.

    ``float_precision="round_trip"`` everywhere, without exception. Several
    assertions below compare floats for EXACT equality across two files, and
    pandas' default parser is fast rather than exact -- it can land one ULP away
    from the value ``to_csv`` wrote. Mixing the two parsers therefore
    manufactures a difference where the files agree perfectly, which is both a
    false alarm and (worse) a temptation to replace the exact comparison with a
    tolerance. ``export/ipa_export.py`` makes the same choice for the same
    reason.
    """
    return pd.read_csv(path, float_precision="round_trip")


@pytest.fixture(scope="module")
def foldchange_all(results_dir) -> pd.DataFrame:
    return _read(results_dir / "foldchange_all.csv")


@pytest.fixture(scope="module")
def qc_limma(results_dir) -> pd.DataFrame:
    return _read(results_dir / "qc_limma.csv")


@pytest.fixture(scope="module")
def ipa_input(results_dir) -> pd.DataFrame:
    return _read(results_dir / "ipa_input.csv")


@pytest.fixture(scope="module")
def single_condition(results_dir) -> pd.DataFrame:
    return _read(results_dir / "single_condition_proteins.csv")


@pytest.fixture(scope="module")
def onoff_proteins(results_dir) -> pd.DataFrame:
    return _read(results_dir / "onoff_proteins.csv")


# ---------------------------------------------------------------------------
# 1. Population invariants -- who is in which file
# ---------------------------------------------------------------------------

def test_qc_limma_ids_are_a_subset_of_foldchange_accessions(qc_limma, foldchange_all):
    """limma may only test proteins the fold-change stage produced.

    An id in ``qc_limma.csv`` that is absent from ``foldchange_all.csv`` means
    the statistics were computed over a row set the fold changes never saw --
    e.g. limma reading a stale file, or a single-condition protein leaking into
    the test.
    """
    fc_accessions = set(foldchange_all["UniProt Accession Number"])
    orphans = sorted(set(qc_limma["id"]) - fc_accessions)
    assert not orphans, (
        f"{len(orphans)} protein(s) in qc_limma.csv are not in "
        f"foldchange_all.csv: {orphans[:5]}"
    )


def test_qc_limma_row_count_is_foldchange_minus_onoff(
    qc_limma, foldchange_all, frozen_counts
):
    """Eligibility rule, stated across two files.

    ``limma_test.py`` excludes ON_OFF proteins: testing them would mean
    inventing an entire absent condition. So the tested population is exactly
    the fold-change population minus the on/off ones. The subtraction is spelled
    out rather than compared to ``qc_limma_rows`` alone, because the point is
    that the three numbers are *linked*, not that each matches a stored value.
    """
    n_onoff = frozen_counts["onoff_rows"]
    assert len(qc_limma) == len(foldchange_all) - n_onoff, (
        f"qc_limma has {len(qc_limma)} rows; foldchange_all has "
        f"{len(foldchange_all)} and {n_onoff} of them are ON_OFF, so "
        f"{len(foldchange_all) - n_onoff} were expected."
    )
    # ... and the ON_OFF count itself must be the one in the fold-change table,
    # not merely the number recorded in frozen_counts.json.
    assert int((foldchange_all["regulated"] == "ON_OFF").sum()) == n_onoff
    assert len(qc_limma) == frozen_counts["qc_limma_rows"]


def test_onoff_file_is_exactly_the_onoff_rows_of_foldchange(
    onoff_proteins, foldchange_all, frozen_counts
):
    """``onoff_proteins.csv`` is a view, so it must not drift from its source."""
    fc_onoff = foldchange_all[foldchange_all["regulated"] == "ON_OFF"]
    assert len(onoff_proteins) == frozen_counts["onoff_rows"]
    assert list(onoff_proteins["accession"]) == list(
        fc_onoff["UniProt Accession Number"]
    ), "onoff_proteins.csv and foldchange_all.csv disagree on which rows are ON_OFF"
    assert list(onoff_proteins["onoff"]) == list(fc_onoff["onoff"])


def test_regulated_column_is_a_partition(foldchange_all, frozen_counts):
    """Every row carries exactly one of the four labels, and the counts add up."""
    labels = set(foldchange_all["regulated"].unique())
    assert labels <= REGULATED_CLASSES, f"unexpected label(s): {labels - REGULATED_CLASSES}"
    counts = foldchange_all["regulated"].value_counts()
    assert counts.sum() == len(foldchange_all)
    assert int(counts.get("UP", 0)) == frozen_counts["n_up"]
    assert int(counts.get("DOWN", 0)) == frozen_counts["n_down"]
    assert (
        int(counts.get("NO CHANGE", 0)) + int(counts.get("ON_OFF", 0))
        == frozen_counts["n_nochange_plus_onoff"]
    )


# ---------------------------------------------------------------------------
# 2. The IPA export is a recomputable selection
# ---------------------------------------------------------------------------

def test_ipa_input_rows_equal_the_recomputed_selection(ipa_input, foldchange_all):
    """``ipa_input.csv`` == (complete AND regulated), recomputed from source.

    This is the invariant that a re-baselined hash manifest cannot protect: the
    export could be filtered on the wrong mask (dropping the completeness
    requirement, or including ON_OFF) and still be perfectly byte-stable. So the
    selection is recomputed here from ``foldchange_all.csv`` and compared
    row-for-row, not merely count-for-count.
    """
    mask = foldchange_all["complete"] & foldchange_all["regulated"].isin(IPA_CLASSES)
    expected = foldchange_all[mask]
    assert len(ipa_input) == int(mask.sum()), (
        f"ipa_input.csv has {len(ipa_input)} rows but recomputing "
        f"(complete & regulated in {sorted(IPA_CLASSES)}) over foldchange_all.csv "
        f"gives {int(mask.sum())}"
    )
    assert list(ipa_input["UniProt Accession Number"]) == list(
        expected["UniProt Accession Number"]
    ), "ipa_input.csv selects different proteins than the recomputed mask"
    # Values, not just membership: a correct selection carrying stale numbers
    # would otherwise pass.
    assert np.array_equal(
        ipa_input["log2FC"].to_numpy(float), expected["log2FC"].to_numpy(float)
    ), "ipa_input.csv log2FC values differ from foldchange_all.csv"
    assert list(ipa_input["regulated"]) == list(expected["regulated"])


def test_ipa_input_carries_no_incomplete_or_onoff_row(ipa_input):
    """Restated from the IPA side, so a failure names the offending file."""
    assert set(ipa_input["regulated"]) <= IPA_CLASSES, (
        "ipa_input.csv contains a class other than UP/DOWN -- an ON_OFF or "
        "NO CHANGE protein reached the QIAGEN upload."
    )
    assert ipa_input["log2FC"].notna().all(), "a NaN log2FC reached the IPA upload"
    assert np.isfinite(ipa_input["log2FC"].to_numpy(float)).all()


# ---------------------------------------------------------------------------
# 3. The enrichment background spans both detection files
# ---------------------------------------------------------------------------

def test_background_union_spans_foldchange_plus_single_condition(
    foldchange_all, single_condition, frozen_counts
):
    """The detected proteome is every protein seen at all, in either file.

    DECISIONS_LOG **D6** is the report's headline finding and it rests entirely
    on this number: against the honest detected-proteome background, 0 GO/KEGG/
    Reactome terms survive correction, while g:Profiler's default whole-genome
    background manufactures 196. If the background silently loses the
    single-condition proteins it becomes a different (smaller, more permissive)
    test and D6's conclusion is no longer the one that was checked.

    Also pins **D11**: the two quarantined junk accessions must be gone from
    ``single_condition_proteins.csv`` *and* out of the background.
    """
    assert (
        len(foldchange_all) + len(single_condition)
        == frozen_counts["background_union"]
    ), (
        f"background_union is {frozen_counts['background_union']} but "
        f"foldchange_all ({len(foldchange_all)}) + single_condition "
        f"({len(single_condition)}) = "
        f"{len(foldchange_all) + len(single_condition)}"
    )
    assert len(single_condition) == frozen_counts["single_condition_rows"]


def test_quarantine_accounts_for_the_dropped_single_condition_rows(
    results_dir, single_condition, frozen_counts
):
    """D11: quarantine means *set aside with a reason*, never *deleted*.

    The rows removed from ``single_condition_proteins.csv`` must reappear, in
    full, in the quarantine record -- otherwise "quarantined" is just a nicer
    word for data loss.
    """
    quarantine = _read(results_dir / "qc" / "quarantine_accessions.csv")
    assert len(quarantine) == frozen_counts["quarantined_accessions"]
    assert (
        len(single_condition) + len(quarantine)
        == frozen_counts["single_condition_rows_pre_quarantine"]
    )
    # The record must carry the value itself, not only a preview of it.
    assert quarantine["value"].notna().all()
    assert (quarantine["n_chars"] == quarantine["value"].str.len()).all()
    quarantined_accessions = set(quarantine["value"])
    assert not (quarantined_accessions & set(single_condition["accession"])), (
        "a quarantined accession is still present in single_condition_proteins.csv"
    )


def test_single_condition_and_foldchange_populations_are_disjoint(
    foldchange_all, single_condition
):
    """A protein is either seen in both sheets or in one, never both.

    Overlap would mean the merge indicator was mis-split, and would double-count
    those proteins in the enrichment background.
    """
    overlap = set(single_condition["accession"]) & set(
        foldchange_all["UniProt Accession Number"]
    )
    assert not overlap, f"{len(overlap)} protein(s) in both files: {sorted(overlap)[:5]}"


def test_single_condition_labels_use_both_directions(single_condition):
    """``detected_in`` must be resolved from the sample sheet, not hardcoded.

    The label was once literally ``left_only -> "control_only"``, which silently
    mislabelled every rescued protein the moment D7 corrected the assignment.
    Both directions must be present and no third value may appear.
    """
    labels = set(single_condition["detected_in"])
    assert labels == {"control_only", "treated_only"}, labels


# ---------------------------------------------------------------------------
# 4. The wiring check -- the two fold changes are the same contrast
# ---------------------------------------------------------------------------

def test_pipeline_and_limma_log2fc_agree(foldchange_all, qc_limma):
    """corr(foldchange log2FC, limma logFC) must be ~1 over the paired subset.

    THE most load-bearing number in the repo. The pipeline's log2FC is a pandas
    mean of per-replicate log2 ratios; limma's is the coefficient of a linear
    model fitted in R over a MinProb-imputed log2 matrix. Nothing structural
    forces them to agree -- they agree only if the contrast is wired the same
    way on both sides. This is the check that proved the D7 control/treated
    correction was a genuine relabelling and not a coding error, and a **stale
    NaN** for it once shipped in the report, which is why it is asserted here
    rather than merely computed somewhere.

    Correlation, not equality: limma imputes missing values, so the ~360
    proteins with a missing intensity legitimately differ a little. The subset
    is the paired non-null rows.
    """
    merged = foldchange_all[["UniProt Accession Number", "log2FC"]].merge(
        qc_limma[["id", "limma_log2FC"]],
        left_on="UniProt Accession Number",
        right_on="id",
        how="inner",
        validate="one_to_one",
    )
    paired = merged.dropna(subset=["log2FC", "limma_log2FC"])
    assert len(paired) > 0, "no paired non-null log2FC values to correlate"

    r = float(np.corrcoef(paired["log2FC"], paired["limma_log2FC"])[0, 1])
    assert np.isfinite(r), (
        "corr(log2FC, limma_log2FC) is NaN. A NaN here is not 'no data' -- it "
        "means one of the two columns is constant or non-numeric, and it is "
        "exactly the failure that once shipped into report_facts.json."
    )
    assert r > LOGFC_CORR_FLOOR, (
        f"corr(pipeline log2FC, limma logFC) = {r:.6f} over {len(paired)} "
        f"paired proteins, below {LOGFC_CORR_FLOOR}. The two sides are no "
        "longer computing the same contrast -- do NOT lower this floor."
    )


def test_limma_flavours_share_one_set_of_fold_changes(results_dir, qc_limma):
    """DECISIONS_LOG **D9**, checked across the two files it produced.

    ``eBayes(trend, robust)`` moderates the variance; it does not refit the
    model, so it cannot move a single logFC. ``limma_test.py`` asserts this
    inside the run; asserting it again *between the committed files* catches the
    case where the two were produced by different runs.
    """
    vanilla = _read(results_dir / "qc_limma_vanilla.csv")
    assert list(vanilla["id"]) == list(qc_limma["id"])
    assert np.array_equal(
        qc_limma["limma_log2FC"].to_numpy(float),
        vanilla["limma_log2FC"].to_numpy(float),
    ), (
        "qc_limma.csv and qc_limma_vanilla.csv disagree on logFC. The eBayes "
        "flavour is supposed to change only the moderated t and the p-values."
    )
    # ... and the flavours must genuinely differ somewhere, or the 'vanilla
    # companion' is silently running the same model twice and proving nothing.
    assert not np.array_equal(
        qc_limma["adj_p_value"].to_numpy(float),
        vanilla["adj_p_value"].to_numpy(float),
    ), "the two eBayes flavours produced identical p-values"


def test_limma_results_contract_file_matches_qc_limma(results_dir, qc_limma):
    """``results/de/limma_results.tsv`` restates qc_limma.csv in limma's names."""
    de = pd.read_csv(results_dir / "de" / "limma_results.tsv", sep="\t", float_precision="round_trip")
    assert list(de["accession"]) == list(qc_limma["id"])
    assert np.allclose(de["logFC"], qc_limma["limma_log2FC"], rtol=0, atol=1e-6)
    # research1.md line 132: fold_change = 2^logFC, derived in the same row.
    assert np.allclose(de["fold_change"], np.power(2.0, de["logFC"]), rtol=1e-12)
    assert np.array_equal(de["n_imputed"].to_numpy(int), qc_limma["n_imputed"].to_numpy(int))


def test_intensity_matrix_covers_exactly_the_tested_proteins(results_dir, qc_limma):
    """The Section-1 contract matrix and the limma audit must describe one run."""
    matrix = pd.read_csv(results_dir / "de" / "intensity_matrix.tsv", sep="\t")
    assert list(matrix["accession"]) == list(qc_limma["id"])


def test_ipa_input_full_selects_the_same_proteins(results_dir, ipa_input):
    """``ipa_input_full.csv`` is ``ipa_input.csv`` plus limma columns.

    Row set and fold changes only -- the p-value half is asserted separately
    below, because it is currently broken (see that test).
    """
    full = _read(results_dir / "ipa_input_full.csv")
    assert list(full["UniProt Accession Number"]) == list(
        ipa_input["UniProt Accession Number"]
    )
    assert np.array_equal(
        full["log2FC"].to_numpy(float), ipa_input["log2FC"].to_numpy(float)
    ), "ipa_input_full.csv and ipa_input.csv disagree on log2FC"
    assert list(full["regulated"]) == list(ipa_input["regulated"])


def test_ipa_input_full_carries_the_default_flavours_p_values(
    results_dir, qc_limma
):
    """The enriched IPA deliverable must quote the p-values the report quotes.

    ``ipa_input_full.csv`` is what a human uploads to QIAGEN. If its FDR column
    comes from a different variance model than ``qc_limma.csv``, the deliverable
    and the report disagree about the same proteins -- here by a factor of ~2.6
    at the minimum (0.116 trend/robust vs 0.305 vanilla). Both still say "0
    significant", so no conclusion changes, but the numbers on the two documents
    are not the same numbers.

    Byte-identity cannot catch this: the stale file is perfectly stable, and its
    hash was re-baselined along with everything else.
    """
    full = _read(results_dir / "ipa_input_full.csv")
    limma_by_id = qc_limma.set_index("id")
    expected_adj = limma_by_id.loc[full["UniProt Accession Number"], "adj_p_value"]
    expected_p = limma_by_id.loc[full["UniProt Accession Number"], "p_value"]
    assert np.allclose(full["p_value"], expected_p.to_numpy(float), rtol=1e-12), (
        "ipa_input_full.csv p_value does not come from qc_limma.csv"
    )
    assert np.allclose(full["adj_p_value"], expected_adj.to_numpy(float), rtol=1e-12), (
        "ipa_input_full.csv adj_p_value does not come from qc_limma.csv"
    )


def test_ipa_full_pvalues_are_not_the_vanilla_run(results_dir):
    """The complement of the test above: rule out the specific way this broke.

    ``ipa_input_full.csv`` shipped for one commit carrying the VANILLA eBayes
    p-values -- they matched ``qc_limma_vanilla.csv`` to the bit -- while the
    report and every figure quoted the trend/robust numbers D9 made the default.
    The code was correct; ``export/ipa_export.py`` simply was not one of
    ``run_pipeline.py``'s stages, so ``--all`` refreshed ``qc_limma.csv`` and
    never rebuilt the file that quotes it. It is a stage now.

    Asserting only "matches qc_limma.csv" would not catch a revert that ALSO
    reverted qc_limma. Naming the wrong source is what makes this specific.
    """
    full = _read(results_dir / "ipa_input_full.csv")
    vanilla = _read(results_dir / "qc_limma_vanilla.csv").set_index("id")
    merged = full.join(vanilla["adj_p_value"], on="UniProt Accession Number",
                       rsuffix="_vanilla")
    delta = (merged["adj_p_value"] - merged["adj_p_value_vanilla"]).abs().max()
    assert delta > 1e-6, (
        "ipa_input_full.csv's adj_p_value matches qc_limma_vanilla.csv, meaning "
        "the QIAGEN deliverable has reverted to the non-default eBayes flavour. "
        "Rebuild with `python -m proteomics_de.export.ipa_export`."
    )


# ---------------------------------------------------------------------------
# 5. Direction sanity (post-D7)
# ---------------------------------------------------------------------------

def test_regulated_label_agrees_with_the_sign_of_log2fc(foldchange_all):
    """Every UP row has log2FC > 0 and every DOWN row has log2FC < 0.

    DECISIONS_LOG **D7** was a *sign* error: 31578/31580 are testosterone and
    31579/31581 are vehicle, and the pipeline had them the other way round.
    A future re-inversion would flip 509 UP and 206 DOWN into each other while
    keeping every count identical -- so a count check cannot see it. This can:
    the labels are derived from the fold changes, so a relabelling that does not
    also move the numbers is incoherent.
    """
    up = foldchange_all[foldchange_all["regulated"] == "UP"]["log2FC"]
    down = foldchange_all[foldchange_all["regulated"] == "DOWN"]["log2FC"]
    assert len(up) and len(down)
    assert (up > 0).all(), (
        f"{int((up <= 0).sum())} UP row(s) have a non-positive log2FC "
        f"(min = {up.min()})"
    )
    assert (down < 0).all(), (
        f"{int((down >= 0).sum())} DOWN row(s) have a non-negative log2FC "
        f"(max = {down.max()})"
    )


def test_up_and_down_clear_the_regulation_threshold(foldchange_all):
    """The +/- cutoff is symmetric (Bug 2), so |log2FC| must clear it both ways."""
    from proteomics_de.config.constants import LOG2_THRESHOLD

    up = foldchange_all[foldchange_all["regulated"] == "UP"]["log2FC"]
    down = foldchange_all[foldchange_all["regulated"] == "DOWN"]["log2FC"]
    assert up.min() >= LOG2_THRESHOLD
    assert down.max() <= -LOG2_THRESHOLD
    nochange = foldchange_all[
        foldchange_all["complete"] & (foldchange_all["regulated"] == "NO CHANGE")
    ]["log2FC"]
    assert (nochange.abs() < LOG2_THRESHOLD).all(), (
        "a complete NO CHANGE row is outside the +/-threshold band"
    )


def test_limma_logfc_signs_agree_with_the_pipeline_labels(foldchange_all, qc_limma):
    """The direction check, restated across the Python/R boundary.

    If R were handed the design with the groups the other way round, the
    correlation test above would go to -1 and this one names the reason.
    """
    merged = foldchange_all[["UniProt Accession Number", "regulated", "log2FC"]].merge(
        qc_limma[["id", "limma_log2FC"]],
        left_on="UniProt Accession Number", right_on="id", how="inner",
    )
    directional = merged[merged["regulated"].isin(IPA_CLASSES)]
    disagreeing = directional[
        np.sign(directional["log2FC"]) != np.sign(directional["limma_log2FC"])
    ]
    assert disagreeing.empty, (
        f"{len(disagreeing)} regulated protein(s) have opposite signs in "
        f"foldchange_all.csv and qc_limma.csv, e.g. "
        f"{disagreeing['UniProt Accession Number'].head(3).tolist()}"
    )


def test_onoff_rows_have_no_numeric_fold_change(foldchange_all):
    """An ON_OFF protein is incomplete by definition, so its log2FC is NaN."""
    onoff = foldchange_all[foldchange_all["regulated"] == "ON_OFF"]
    assert onoff["log2FC"].isna().all(), "an ON_OFF protein carries a numeric log2FC"
    assert (~onoff["complete"]).all(), "an ON_OFF protein is marked complete"
    assert set(onoff["onoff"]) == {"on_with_treatment", "off_with_treatment"}


def test_complete_rows_carry_a_finite_fold_change(foldchange_all):
    """Bug 3: no inf, no NaN, ever -- restated on the committed file."""
    complete = foldchange_all[foldchange_all["complete"]]["log2FC"].to_numpy(float)
    assert np.isfinite(complete).all(), "inf/NaN log2FC on a complete row"
    incomplete = foldchange_all[~foldchange_all["complete"]]["log2FC"]
    assert incomplete.isna().all(), "an incomplete row carries a fold change"


# ---------------------------------------------------------------------------
# 6. Header-only outputs are a CONTRACT, not a tolerated failure
# ---------------------------------------------------------------------------

_HEADER_ONLY_REASON = {
    "ipa_input_significant.csv": (
        "DECISIONS_LOG D2: 0 of 1938 proteins pass FDR<0.05 (min adj.p = 0.116 "
        "under the D9 trend/robust default). That is the expected ceiling of an "
        "n=2 TECHNICAL-replicate design, not a bug and not a missing file."
    ),
    "enrichment/ora_up.csv": (
        "DECISIONS_LOG D6: 0 GO/KEGG/Reactome terms survive g:Profiler g:SCS "
        "correction against the honest detected-proteome background. The same "
        "query against the DEFAULT whole-genome background returns 196 "
        "'significant' terms -- a textbook background-inflation artifact the "
        "pipeline deliberately does not fall into."
    ),
    "enrichment/ora_down.csv": (
        "DECISIONS_LOG D6, down direction -- same reasoning as ora_up.csv."
    ),
}

_HEADER_ONLY_COLUMNS = {
    "ipa_input_significant.csv": [
        "UniProt Accession Number", "Gene names", "log2FC", "regulated", "adj_p_value",
    ],
    "enrichment/ora_up.csv": [
        "source", "term_id", "term_name", "p_value", "term_size", "query_size",
        "intersection_size", "intersecting_genes",
    ],
    "enrichment/ora_down.csv": [
        "source", "term_id", "term_name", "p_value", "term_size", "query_size",
        "intersection_size", "intersecting_genes",
    ],
}


@pytest.mark.parametrize("relpath", sorted(_HEADER_ONLY_COLUMNS))
def test_header_only_output_has_its_columns_and_zero_data_rows(results_dir, relpath):
    """These three files are *supposed* to be empty, and the header proves it ran.

    Emptiness is the scientific result here, so the file must be present, must
    carry its full column header (an empty file would be indistinguishable from
    a stage that crashed before writing), and must carry exactly zero data rows.

    IF A ROW EVER APPEARS HERE, THAT IS A FINDING TO INVESTIGATE, NOT A TEST TO
    RELAX. A non-empty ``ipa_input_significant.csv`` means something passed
    FDR<0.05 in a design that cannot support it; a non-empty ``ora_*.csv`` means
    the enrichment ran against a different (probably whole-genome) background.
    Both are claims the report would then be making, and both need a human to
    look before the number is believed -- deleting the assertion would let the
    pipeline start asserting significance nobody checked.
    """
    path = results_dir / relpath
    assert path.is_file(), f"{relpath} is missing entirely"
    df = pd.read_csv(path)
    assert list(df.columns) == _HEADER_ONLY_COLUMNS[relpath], (
        f"{relpath} column header drifted: {list(df.columns)}"
    )
    assert len(df) == 0, (
        f"{relpath} now has {len(df)} data row(s), and it must have 0.\n"
        f"REASON IT IS EMPTY: {_HEADER_ONLY_REASON[relpath]}\n"
        "A row appearing here is a finding to investigate, not a test to relax."
    )


def test_no_protein_is_both_regulated_and_significant(results_dir, qc_limma):
    """``ipa_input_significant.csv`` is empty *because the data says so*.

    Recomputed from ``qc_limma.csv`` rather than trusting the empty file: if the
    significance filter were simply broken (always False), the file would be
    empty for the wrong reason and the test above would still pass.
    """
    regulated_and_significant = qc_limma[
        qc_limma["regulated"].isin(IPA_CLASSES) & qc_limma["significant"]
    ]
    n_sig = len(pd.read_csv(results_dir / "ipa_input_significant.csv"))
    assert len(regulated_and_significant) == n_sig
    # The filter must be live, not stuck: `significant` is a real comparison
    # against adj_p_value, so it has to agree with recomputing that comparison.
    from proteomics_de.config.constants import ADJ_P_THRESHOLD

    recomputed = qc_limma["adj_p_value"] < ADJ_P_THRESHOLD
    assert np.array_equal(
        qc_limma["significant"].to_numpy(bool), recomputed.to_numpy(bool)
    ), "the `significant` column disagrees with adj_p_value < ADJ_P_THRESHOLD"


def test_gsea_results_are_present_but_none_survive_fdr(results_dir, frozen_counts):
    """D6's other half: GSEA ran over a real ranking and still found nothing.

    Unlike ORA this file is NOT header-only -- 568 terms were scored. The
    finding is that none pass FDR<0.05, which is a different (and stronger)
    statement than 'the file is empty'.
    """
    gsea = pd.read_csv(results_dir / "enrichment" / "gsea_results.csv")
    assert len(gsea) == frozen_counts["gsea_terms"]
    assert "fdr_q_value" in gsea.columns, list(gsea.columns)
    assert (gsea["fdr_q_value"] >= 0.05).all(), (
        "a GSEA term now passes FDR<0.05. Per DECISIONS_LOG D6 none did; this "
        "is a finding to investigate, not a threshold to move."
    )


# ---------------------------------------------------------------------------
# 7. The frozen-count file is itself consistent
# ---------------------------------------------------------------------------

def test_frozen_counts_are_internally_consistent(frozen_counts):
    """The single source of truth must not contradict itself.

    Cheap, and it is the file every other assertion in this module trusts.
    """
    c = frozen_counts
    assert c["qc_limma_rows"] == c["foldchange_all_rows"] - c["onoff_rows"]
    assert c["ipa_input_rows"] == c["n_up"] + c["n_down"]
    assert (
        c["n_up"] + c["n_down"] + c["n_nochange_plus_onoff"] == c["foldchange_all_rows"]
    )
    assert (
        c["single_condition_rows"] + c["quarantined_accessions"]
        == c["single_condition_rows_pre_quarantine"]
    )
    assert (
        c["background_union"] == c["foldchange_all_rows"] + c["single_condition_rows"]
    )
    assert c["n_significant_fdr05"] == 0, (
        "frozen_counts.json now claims a significant protein -- see D2 before "
        "believing it"
    )


#: The one non-count, non-comment key: the commit the counts were derived at.
_PROVENANCE_KEYS = {"baseline_git_head"}


def test_frozen_counts_json_documents_every_number(repo_root):
    """Every key is a count, the provenance key, or a leading-underscore comment.

    Guards against a bare number being added with no explanation of where it
    came from, which is how the eight scattered literals this file replaced got
    there in the first place.
    """
    path = repo_root / "proteomics_de" / "tests" / "expected" / "frozen_counts.json"
    raw = json.loads(path.read_text(encoding="utf-8"))
    for key, value in raw.items():
        if key.startswith("_"):
            assert isinstance(value, str) and value.strip(), key
        elif key in _PROVENANCE_KEYS:
            assert isinstance(value, str) and len(value) == 40, (
                f"{key} should be a full git sha, got {value!r}"
            )
            int(value, 16)  # raises if it is not hex
        else:
            assert isinstance(value, int), f"{key} is not an int count: {value!r}"
    assert raw["_comment"].strip()
