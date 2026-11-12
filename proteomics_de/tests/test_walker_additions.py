"""Tests for :mod:`proteomics_de.export.walker_additions` (DECISIONS_LOG D20).

Three things are being defended:

1. **The additions are exactly what Dr. Walker approved**: 715 core rows + FRZB
   (REAL values, ``value_source=measured``) + 5 control-only proteins with an
   ASSIGNED 100-fold decrease (``value_source=assigned``, blank p-values) = 721.
   Since D24 (his final answer) the recovered-protein rule is the original STRICT
   one again (all four readings seen, ``n_imputed == 0``), which admits nothing, so
   the 247 partial-tier proteins the user's WIDENED rule (D20) had admitted are
   gone; that widened code path is still tested, by setting the two config
   switches to false in a tmp copy of the config (``widened_cfg``). The file
   WITHOUT the 5 assigned rows (716, the IPA sensitivity dataset) is tested too.
2. **Assigned numbers can never pass as measurements**: labelled on every row,
   p-values blank, stated in the provenance sidecar and the manifest.
3. **Nothing that already exists moves**: the originals stay byte-identical and
   a re-run reproduces identical bytes.

Tests that need the committed results use the ``results_dir`` fixture; the rule
and invariant tests run on small synthetic frames so they cannot be satisfied by
accident of the data.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
sys.path.insert(0, str(_REPO_ROOT / "tools"))

import freeze  # noqa: E402
from proteomics_de.config import constants  # noqa: E402
from proteomics_de.config import design  # noqa: E402
from proteomics_de.export import walker_additions as wa  # noqa: E402
from proteomics_de.export.walker_additions import WalkerAdditionsError  # noqa: E402

ID = "UniProt Accession Number"
ASSIGNED_LOG2FC = -6.643856
MANIFEST = _REPO_ROOT / "proteomics_de" / "tests" / "expected" / "outputs.sha256"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _read(path):
    return pd.read_csv(path, float_precision="round_trip")


@pytest.fixture(scope="module")
def cfg():
    return wa.load_config()


@pytest.fixture(scope="module")
def widened_cfg(cfg):
    """The committed config with both 'no problems' switches OFF, in a copy (D24)."""
    w = copy.deepcopy(cfg)
    w["recovered_rule"].update(require_all_replicates_seen=False,
                               require_no_imputed_values=False)
    return w


@pytest.fixture(scope="module")
def without(results_dir):
    path = results_dir / wa.WALKER_NOASSIGNED_CSV_FILENAME
    if not path.exists():
        pytest.skip(f"{path} not generated yet")
    return _read(path)


@pytest.fixture(scope="module")
def walker(results_dir):
    path = results_dir / wa.WALKER_CSV_FILENAME
    if not path.exists():
        pytest.skip(f"{path} not generated yet")
    return _read(path)


@pytest.fixture(scope="module")
def manifest_df(results_dir):
    path = results_dir / wa.WALKER_MANIFEST_FILENAME
    if not path.exists():
        pytest.skip(f"{path} not generated yet")
    return _read(path)


@pytest.fixture(scope="module")
def full(results_dir):
    return _read(results_dir / "ipa_input_full.csv")


@pytest.fixture
def inputs_dir(results_dir, tmp_path):
    """A scratch results dir holding only the inputs the module reads."""
    names = set(wa._INPUTS) | set(wa.BYTE_FROZEN_INPUTS)
    for name in names:
        src = results_dir / name
        if src.exists():
            shutil.copy(src, tmp_path / name)
    return tmp_path


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Counts: 715 -> +FRZB -> +5 assigned -> +recovered
# ---------------------------------------------------------------------------
def test_row_count_is_core_plus_additions(walker, cfg, frozen_counts):
    n_core = frozen_counts["ipa_input_rows"]
    n_measured = len(cfg["measured_additions"])
    n_assigned = len(cfg["assigned_proteins"])
    n_recovered = frozen_counts["n_walker_recovered"]
    assert len(walker) == n_core + n_measured + n_assigned + n_recovered
    assert len(walker) == frozen_counts["n_ipa_walker"]
    assert (n_measured, n_assigned) == (1, 5)


def test_counts_at_each_step(walker, frozen_counts, cfg):
    """715 -> +1 (FRZB) -> +5 (assigned) -> +recovered, in file order."""
    n_core = frozen_counts["ipa_input_rows"]
    assert list(walker["tier"].iloc[:n_core].unique()) == ["core"]
    after_frzb = walker.iloc[: n_core + 1]
    assert after_frzb.iloc[-1][ID] == "P97401"
    after_assigned = walker.iloc[: n_core + 1 + len(cfg["assigned_proteins"])]
    assert list(after_assigned.iloc[n_core + 1:][ID]) == [
        e["accession"] for e in cfg["assigned_proteins"]]
    assert len(walker) - len(after_assigned) == frozen_counts["n_walker_recovered"]


def test_tier_and_value_source_breakdown(walker, frozen_counts):
    tiers = walker["tier"].value_counts().to_dict()
    n_rec = frozen_counts["n_walker_recovered"]
    assert tiers == {"core": frozen_counts["ipa_input_rows"], "partial": 1 + n_rec,
                     "qualitative": 5}
    assert walker["value_source"].value_counts().to_dict() == {
        "measured": frozen_counts["ipa_input_rows"] + 1 + n_rec, "assigned": 5}


def test_columns_are_the_declared_contract(walker):
    assert list(walker.columns) == wa.WALKER_COLS


def test_no_duplicate_accessions(walker):
    assert walker[ID].is_unique


def test_manifest_lists_exactly_the_additions(manifest_df, walker, full, frozen_counts):
    assert list(manifest_df.columns) == wa.MANIFEST_COLS
    assert len(manifest_df) == frozen_counts["n_walker_additions_listed"]
    added = set(walker[ID]) - set(full[ID])
    assert set(manifest_df[ID]) == added
    assert manifest_df[ID].is_unique
    assert manifest_df["reason"].notna().all() and manifest_df["source"].notna().all()


# ---------------------------------------------------------------------------
# Assigned rows
# ---------------------------------------------------------------------------
def test_assigned_rows_carry_the_assigned_value_and_label(walker, cfg):
    acc = [e["accession"] for e in cfg["assigned_proteins"]]
    assert acc == ["P01132", "Q61521", "O08542", "P55144", "Q9D1H9"]
    rows = walker[walker[ID].isin(acc)]
    assert len(rows) == 5
    assert (rows["log2FC"] == ASSIGNED_LOG2FC).all()
    assert (rows["value_source"] == "assigned").all()
    assert (rows["regulated"] == "DOWN").all()
    # 100-fold decrease: fold_change = 2**log2FC = 0.01 (to the 6-decimal rounding).
    assert np.allclose(2.0 ** rows["log2FC"], 0.01, rtol=1e-5)


def test_assigned_rows_have_blank_p_values_in_the_raw_text(results_dir, cfg):
    """Blank in the file itself, not merely NaN after parsing."""
    import csv

    acc = {e["accession"] for e in cfg["assigned_proteins"]}
    for name, delim in ((wa.WALKER_CSV_FILENAME, ","), (wa.WALKER_TXT_FILENAME, "\t")):
        with open(results_dir / name, newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter=delim))
        hit = [r for r in rows if r[ID] in acc]
        assert len(hit) == 5, name
        for r in hit:
            assert r["p_value"] == "" and r["adj_p_value"] == "", (name, r)
            assert r["value_source"] == "assigned"


def test_assigned_log2fc_is_derived_from_config_not_typed(cfg):
    assert wa.assigned_log2fc(cfg) == ASSIGNED_LOG2FC
    assert math.isclose(wa.assigned_log2fc(cfg), math.log2(0.01), abs_tol=1e-6)
    assert wa.assigned_fold_change(cfg) == 0.01
    other = copy.deepcopy(cfg)
    other["assigned_fold_change"]["decrease_factor"] = 10
    assert wa.assigned_log2fc(other) == round(-math.log2(10), 6)


def test_the_five_are_control_only_in_the_data(results_dir, cfg):
    """Detected in both control readings, in neither testosterone reading."""
    single = _read(results_dir / "single_condition_proteins.csv").set_index("accession")
    ctrl, trt = design.control_columns(), design.treated_columns()
    for e in cfg["assigned_proteins"]:
        row = single.loc[e["accession"]]
        assert row["gene"] == e["gene"]
        assert row["detected_in"] == "control_only"
        assert all(row[c] > 0 for c in ctrl), e
        assert all(not (row[c] > 0) for c in trt), e


def test_control_and_treated_columns_follow_the_sample_sheet():
    """D7: control = 31579/31581, testosterone = 31578/31580."""
    assert design.control_columns() == ["Intensity 31579", "Intensity 31581"]
    assert design.treated_columns() == ["Intensity 31578", "Intensity 31580"]


# ---------------------------------------------------------------------------
# FRZB: real values, not the 100-fold rule
# ---------------------------------------------------------------------------
def test_frzb_keeps_its_real_measured_values(walker, results_dir):
    row = walker[walker[ID] == "P97401"]
    assert len(row) == 1
    row = row.iloc[0]
    assert row["value_source"] == "measured"
    assert row["regulated"] == "UP"
    assert row["tier"] == "partial"
    assert row["log2FC"] == pytest.approx(4.40696, abs=1e-9)
    assert 2.0 ** row["log2FC"] == pytest.approx(21.2, abs=0.1)
    assert row["p_value"] == pytest.approx(0.0203, abs=5e-5)
    assert row["adj_p_value"] == pytest.approx(0.998, abs=5e-4)
    assert row["log2FC"] != ASSIGNED_LOG2FC
    # And it is the very same row the partial tier already has.
    src = _read(results_dir / "regulated_up_partial.csv")
    src = src[src[ID] == "P97401"].iloc[0]
    for col in ("log2FC", "p_value", "adj_p_value"):
        assert row[col] == src[col]
    assert int(src["n_imputed"]) == 1


# ---------------------------------------------------------------------------
# Originals untouched, idempotent
# ---------------------------------------------------------------------------
def test_core_rows_are_ipa_input_full_verbatim(walker, full):
    core = walker.iloc[: len(full)]
    assert list(core[ID]) == list(full[ID])
    for col in ("Gene names", "log2FC", "regulated", "p_value", "adj_p_value"):
        assert list(core[col]) == list(full[col]), col


def test_original_files_still_match_the_freeze_manifest(results_dir):
    entries = freeze.read_manifest(MANIFEST)
    for name in ("ipa_input.csv", "ipa_input_full.csv", "ipa_input_full.txt",
                 "ipa_input_extended.csv", "ipa_input_extended.txt",
                 "regulated_up.csv", "regulated_down.csv",
                 "regulated_up_partial.csv", "regulated_down_partial.csv"):
        rel = f"proteomics_de/results/{name}"
        sha, _mode = entries[rel]
        assert freeze.digest(results_dir / name)[0] == sha, f"{name} drifted"


def test_building_does_not_touch_any_input(inputs_dir):
    before = {n: _sha(inputs_dir / n) for n in wa.BYTE_FROZEN_INPUTS
              if (inputs_dir / n).exists()}
    assert "ipa_input.csv" in before and "ipa_input_extended.csv" in before
    wa.build_walker_from_results(inputs_dir)
    after = {n: _sha(inputs_dir / n) for n in before}
    assert before == after


def test_rerun_is_byte_identical_and_matches_the_committed_files(inputs_dir, results_dir):
    wa.build_walker_from_results(inputs_dir)
    first = {n: (inputs_dir / n).read_bytes() for n in
             (wa.WALKER_CSV_FILENAME, wa.WALKER_TXT_FILENAME, wa.WALKER_MANIFEST_FILENAME)}
    wa.build_walker_from_results(inputs_dir)
    for n, data in first.items():
        assert (inputs_dir / n).read_bytes() == data, f"{n} not idempotent"
        committed = results_dir / n
        if committed.exists():
            assert committed.read_bytes() == data, f"{n} differs from the committed file"


def test_csv_and_txt_twins_agree(results_dir):
    from proteomics_de.export.ipa_extended import assert_twins_agree

    assert_twins_agree(results_dir / wa.WALKER_CSV_FILENAME,
                       results_dir / wa.WALKER_TXT_FILENAME)


def test_no_bom_and_identifier_leftmost(results_dir):
    for name in (wa.WALKER_CSV_FILENAME, wa.WALKER_TXT_FILENAME, wa.WALKER_MANIFEST_FILENAME):
        raw = (results_dir / name).read_bytes()
        assert raw[:3] != b"\xef\xbb\xbf"
        assert raw.decode("utf-8").splitlines()[0].startswith(ID)


# ---------------------------------------------------------------------------
# Direction consistency
# ---------------------------------------------------------------------------
def test_regulated_column_agrees_with_the_sign_of_log2fc(walker):
    assert ((walker["log2FC"] > 0) == (walker["regulated"] == "UP")).all()
    assert (walker["log2FC"] != 0).all()
    assert set(walker["regulated"]) <= {"UP", "DOWN"}


def test_every_regulated_row_clears_the_threshold(walker):
    assert (walker["log2FC"].abs() >= constants.LOG2_THRESHOLD).all()


# ---------------------------------------------------------------------------
# The "recovered" rule -- on the committed data
# ---------------------------------------------------------------------------
def test_the_recovered_set_is_862_and_the_committed_strict_rule_admits_none(results_dir, cfg, frozen_counts):
    """What the "~860" is, and what the committed (strict, D24) rule keeps of it: nothing."""
    result = wa.compute_walker(results_dir, cfg)
    pool = result.counts["recovered_pool"]
    assert pool == (frozen_counts["n_up_partial"] + frozen_counts["n_down_partial"]
                    + frozen_counts["n_qualitative"]) == 862
    assert result.counts["recovered"] == frozen_counts["n_walker_recovered"] == 0
    assert result.counts["total"] == 715 + 1 + 5 + 0 == 721 == frozen_counts["n_ipa_walker"]
    counts = [n for _label, n in result.steps]
    assert counts == [862, 856, 247, 247, 0, 0]
    assert counts == sorted(counts, reverse=True), "each step may only remove rows"


def test_the_widened_rule_still_admits_247_when_the_switches_are_set_false(results_dir, widened_cfg):
    """The D20 code path stays covered: switches false in a tmp config -> 247 recovered, 968 rows."""
    result = wa.compute_walker(results_dir, widened_cfg)
    assert result.counts["recovered"] == 247
    assert result.counts["total"] == 715 + 1 + 5 + 247 == 968
    assert result.counts["without_assigned"] == 963
    counts = [n for _label, n in result.steps]
    assert counts == [862, 856, 247, 247]
    assert counts == sorted(counts, reverse=True), "each step may only remove rows"


def test_the_original_strict_rule_still_admits_none(results_dir, cfg):
    """Flipping the two switches (to true, as committed) gives the strict rule: 0 pass, by construction."""
    strict = copy.deepcopy(cfg)
    strict["recovered_rule"].update(require_all_replicates_seen=True,
                                    require_no_imputed_values=True)
    result = wa.compute_walker(results_dir, strict)
    assert [n for _l, n in result.steps] == [862, 856, 247, 247, 0, 0]
    assert result.counts["recovered"] == 0 and result.counts["total"] == 721


def test_the_247_recovered_proteins_are_the_partial_tier_minus_frzb(inputs_dir, widened_cfg, results_dir, frozen_counts):
    """The widened rule's 247, checked from a tmp-config run (they are NOT in the committed file)."""
    up = _read(results_dir / "regulated_up_partial.csv")
    down = _read(results_dir / "regulated_down_partial.csv")
    partial = pd.concat([up.assign(d="UP"), down.assign(d="DOWN")], ignore_index=True)
    partial = partial[partial[ID] != "P97401"]
    walker = wa.compute_walker(inputs_dir, widened_cfg).frame
    n_core = frozen_counts["ipa_input_rows"]
    rec = walker.iloc[n_core + 1 + 5:]
    assert set(rec[ID]) == set(partial[ID]) and len(rec) == len(partial) == 247
    src = partial.set_index(ID).loc[rec[ID]]
    for col in ("log2FC", "p_value", "adj_p_value"):
        assert list(rec[col]) == list(src[col]), col          # real values, verbatim
    assert list(rec["regulated"]) == list(src["d"])
    assert (rec["tier"] == "partial").all() and (rec["value_source"] == "measured").all()
    abs_lfc = list(rec["log2FC"].abs())
    assert abs_lfc == sorted(abs_lfc, reverse=True)           # descending |log2FC|
    # n_imputed breakdown: 178 with one imputed reading, 69 with two; 151 UP, 96 DOWN.
    assert src["n_imputed"].value_counts().to_dict() == {1: 178, 2: 69}
    assert rec["regulated"].value_counts().to_dict() == {"UP": 151, "DOWN": 96}


def test_the_committed_walker_file_contains_none_of_the_247(walker, results_dir):
    """Dr. Walker's answer (D24): the 247 are dropped. FRZB is the only partial-tier row left."""
    up = _read(results_dir / "regulated_up_partial.csv")
    down = _read(results_dir / "regulated_down_partial.csv")
    partial_ids = (set(up[ID]) | set(down[ID])) - {"P97401"}
    assert len(partial_ids) == 247
    assert not partial_ids & set(walker[ID])
    assert list(walker.loc[walker["tier"] == "partial", ID]) == ["P97401"]


def test_the_recovered_proteins_all_have_adj_p_near_one_like_the_core(inputs_dir, widened_cfg, walker, full):
    """The caveat the docs state, checked here: none of them is FDR-significant.

    Checked on the widened frame (the 247 plus FRZB), since the committed file now holds only FRZB.
    """
    rec = wa.compute_walker(inputs_dir, widened_cfg).frame.query("tier == 'partial'")
    assert len(rec) == 248
    assert (walker.loc[walker["tier"] == "partial", "adj_p_value"] >= 0.05).all()
    assert (rec["adj_p_value"] >= 0.05).all() and rec["adj_p_value"].min() > 0.1
    assert (rec["adj_p_value"] > 0.9).mean() > 0.99
    assert (full["adj_p_value"] >= 0.05).all()                # 0 of the core 715 pass


def test_every_tier3_row_has_imputed_values_and_no_qualitative_row_has_a_fold_change(results_dir):
    """The structural reason the STRICT rule is empty -- not a coincidence of thresholds."""
    for name in ("regulated_up_partial.csv", "regulated_down_partial.csv"):
        assert (_read(results_dir / name)["n_imputed"] >= 1).all()
    qual = _read(results_dir / "qualitative_changes.csv")
    assert "log2FC" not in qual.columns and "p_value" not in qual.columns


def test_rule_threshold_matches_the_pipeline_threshold(cfg):
    rule = cfg["recovered_rule"]
    assert rule["min_abs_log2fc"] == constants.LOG2_THRESHOLD
    assert rule["require_partial_tier"] is True
    # D24 (Dr. Walker: drop the 247): both "no problems" switches are TRUE again, i.e. the
    # strict rule. They stay config switches, so the widened rule is one edit away.
    assert rule["require_all_replicates_seen"] is True
    assert rule["require_no_imputed_values"] is True


# ---------------------------------------------------------------------------
# The recovered rule on synthetic frames (positive AND negative)
# ---------------------------------------------------------------------------
RULE = {"min_abs_log2fc": 0.585, "require_partial_tier": True,
        "require_all_replicates_seen": False, "require_no_imputed_values": False}
STRICT = dict(RULE, require_all_replicates_seen=True, require_no_imputed_values=True)


def _pool(rows):
    cols = wa.POOL_COLS
    base = {"Gene names": "G", "regulated": "UP", "p_value": 0.1, "adj_p_value": 0.9,
            "n_imputed": 1, "n_seen": 3, "source": "partial", "log2FC": 1.0}
    frame = pd.DataFrame([{**base, ID: acc, **kw} for acc, kw in rows])
    return frame[cols]


def test_widened_rule_admits_a_partial_protein_with_an_imputed_value():
    passed, steps = wa.apply_recovered_rule(_pool([("PASS1", {})]), RULE, set(), 4)
    assert list(passed[ID]) == ["PASS1"]
    assert [n for _l, n in steps] == [1, 1, 1, 1]


def test_widened_rule_boundary_is_inclusive_at_the_threshold():
    pool = _pool([("AT", {"log2FC": 0.585}), ("BELOW", {"log2FC": 0.5849}),
                  ("NEGAT", {"log2FC": -0.585, "regulated": "DOWN"})])
    passed, _ = wa.apply_recovered_rule(pool, RULE, set(), 4)
    assert sorted(passed[ID]) == ["AT", "NEGAT"]


@pytest.mark.parametrize("label,kw", [
    ("effect below 50 percent", {"log2FC": 0.4}),
    ("no fold change at all", {"log2FC": np.nan, "n_imputed": np.nan}),
    ("qualitative source tier", {"log2FC": np.nan, "n_imputed": np.nan, "source": "qualitative"}),
    ("qualitative source even with a fold change", {"source": "qualitative"}),
])
def test_widened_rule_still_excludes_what_has_no_usable_fold_change(label, kw):
    pool = _pool([("BAD", kw), ("GOOD", {})])
    passed, steps = wa.apply_recovered_rule(pool, RULE, set(), 4)
    assert list(passed[ID]) == ["GOOD"], label
    assert steps[0][1] == 2 and steps[-1][1] == 1


@pytest.mark.parametrize("label,kw", [
    ("one imputed value", {"n_imputed": 1, "n_seen": 4}),
    ("two imputed values", {"n_imputed": 2, "n_seen": 4}),
    ("a replicate missing", {"n_imputed": 0, "n_seen": 3}),
    ("a whole group missing", {"n_imputed": 0, "n_seen": 2}),
    ("effect below 50 percent", {"n_imputed": 0, "n_seen": 4, "log2FC": 0.4}),
    ("no fold change at all", {"log2FC": np.nan, "n_imputed": np.nan, "n_seen": 4}),
    ("imputation count unknown", {"n_imputed": np.nan, "n_seen": 4}),
])
def test_strict_switches_restore_the_original_rule(label, kw):
    pool = _pool([("BAD", kw), ("GOOD", {"n_imputed": 0, "n_seen": 4})])
    passed, steps = wa.apply_recovered_rule(pool, STRICT, set(), 4)
    assert list(passed[ID]) == ["GOOD"], label
    assert steps[0][1] == 2 and steps[-1][1] == 1


@pytest.mark.parametrize("switch", ["require_all_replicates_seen", "require_no_imputed_values"])
def test_each_switch_works_on_its_own(switch):
    pool = _pool([("MISS", {"n_seen": 3, "n_imputed": 0}),
                  ("IMP", {"n_seen": 4, "n_imputed": 1}),
                  ("OK", {"n_seen": 4, "n_imputed": 0})])
    one = dict(RULE, **{switch: True})
    passed, steps = wa.apply_recovered_rule(pool, one, set(), 4)
    expected = {"require_all_replicates_seen": {"IMP", "OK"},
                "require_no_imputed_values": {"MISS", "OK"}}[switch]
    assert set(passed[ID]) == expected
    assert len(steps) == 5


def test_rule_never_duplicates_a_protein_already_on_the_list():
    pool = _pool([("DUP", {}), ("NEW", {})])
    passed, steps = wa.apply_recovered_rule(pool, RULE, {"DUP"}, 4)
    assert list(passed[ID]) == ["NEW"]
    assert dict(steps)["not already on the list"] == 1


def test_rule_counts_every_step_in_order():
    pool = _pool([("A", {}), ("B", {"log2FC": 0.1}), ("C", {"n_seen": 2}),
                  ("D", {"n_imputed": 2}), ("E", {}), ("F", {})])
    _passed, widened = wa.apply_recovered_rule(pool, RULE, {"F"}, 4)
    assert [n for _l, n in widened] == [6, 5, 5, 4]
    _passed, strict = wa.apply_recovered_rule(pool, STRICT, {"F"}, 4)
    assert [n for _l, n in strict] == [6, 5, 5, 4, 0, 0]


def test_a_partial_row_with_imputed_values_reaches_the_file_end_to_end(inputs_dir, widened_cfg):
    """Widened rule (tmp config): an n_imputed=1 or 2 partial row is added, measured, labelled."""
    result = wa.compute_walker(inputs_dir, widened_cfg)
    n_fixed = result.counts["core"] + result.counts["measured"] + result.counts["assigned"]
    rec = result.frame.iloc[n_fixed:]
    assert len(rec) == 247 and result.counts["total"] == 968
    man = result.manifest.set_index(ID)
    assert set(man.loc[rec[ID], "n_imputed"].astype(int)) == {1, 2}
    assert (rec["tier"] == "partial").all() and (rec["value_source"] == "measured").all()
    assert rec[["p_value", "adj_p_value"]].notna().all().all()
    assert man.loc[rec[ID], "reason"].str.contains("WIDENED").all()


def test_strict_config_makes_the_same_data_admit_nothing_end_to_end(inputs_dir, cfg):
    strict = copy.deepcopy(cfg)
    strict["recovered_rule"].update(require_all_replicates_seen=True,
                                    require_no_imputed_values=True)
    result = wa.compute_walker(inputs_dir, strict)
    assert result.counts["recovered"] == 0 and result.counts["total"] == 721


def test_a_clean_partial_row_passes_even_the_strict_rule_end_to_end(inputs_dir, cfg):
    """Make one partial row genuinely clean: strict admits exactly it."""
    path = inputs_dir / "regulated_up_partial.csv"
    partial = _read(path)
    target = partial.loc[partial[ID] != "P97401", ID].iloc[0]
    partial.loc[partial[ID] == target, "n_imputed"] = 0
    partial.to_csv(path, index=False)
    limma_path = inputs_dir / "qc_limma.csv"
    limma = _read(limma_path)
    for col in design.control_columns() + design.treated_columns():
        limma.loc[limma["id"] == target, col] = 1e6
    limma.to_csv(limma_path, index=False)

    strict = copy.deepcopy(cfg)
    strict["recovered_rule"].update(require_all_replicates_seen=True,
                                    require_no_imputed_values=True)
    result = wa.compute_walker(inputs_dir, strict)
    assert result.counts["recovered"] == 1
    row = result.frame[result.frame[ID] == target].iloc[0]
    assert row["tier"] == "partial" and row["value_source"] == "measured"
    assert row["regulated"] == "UP" and not math.isnan(row["p_value"])
    assert result.frame[ID].is_unique
    assert result.steps[-1][1] == 1
    man = result.manifest.set_index(ID).loc[target]
    assert man["n_imputed"] == 0 and "WIDENED" not in man["reason"]


def test_a_partial_row_with_an_imputed_value_is_not_admitted_by_the_strict_rule(inputs_dir, cfg):
    """All four readings 'seen' but n_imputed left at >= 1: strict must still refuse it."""
    path = inputs_dir / "regulated_up_partial.csv"
    partial = _read(path)
    target = partial.loc[partial[ID] != "P97401", ID].iloc[0]
    limma_path = inputs_dir / "qc_limma.csv"
    limma = _read(limma_path)
    for col in design.control_columns() + design.treated_columns():
        limma.loc[limma["id"] == target, col] = 1e6
    limma.to_csv(limma_path, index=False)
    strict = copy.deepcopy(cfg)
    strict["recovered_rule"].update(require_all_replicates_seen=True,
                                    require_no_imputed_values=True)
    result = wa.compute_walker(inputs_dir, strict)
    assert result.counts["recovered"] == 0
    assert target not in set(result.frame[ID])


# ---------------------------------------------------------------------------
# Invariants and input checks (negative tests)
# ---------------------------------------------------------------------------
def _good_frame(full_like=True):
    core = pd.DataFrame({
        ID: ["C1", "C2"], "Gene names": ["a", "b"], "log2FC": [1.0, -1.0],
        "regulated": ["UP", "DOWN"], "p_value": [0.1, 0.2], "adj_p_value": [0.9, 0.9]})
    add = pd.DataFrame({
        ID: ["A1"], "Gene names": ["x"], "log2FC": [ASSIGNED_LOG2FC], "regulated": ["DOWN"],
        "p_value": [np.nan], "adj_p_value": [np.nan], "tier": ["qualitative"],
        "value_source": ["assigned"]})
    return wa.build_walker_frame(core, [add])


def test_a_clean_frame_passes():
    wa.check_walker_frame(_good_frame())


@pytest.mark.parametrize("mutate,match", [
    (lambda f: f.assign(**{"p_value": f["p_value"].where(f["value_source"] != "assigned", 0.5)}),
     "assigned row carries a p-value"),
    (lambda f: f.assign(regulated=["DOWN", "DOWN", "DOWN"]), "disagrees with the sign"),
    (lambda f: f.assign(value_source=["measured", "measured", "guess"]), "value_source"),
    (lambda f: f.assign(regulated=["UP", "DOWN", "SIDEWAYS"]), "regulated"),
    (lambda f: f.assign(log2FC=[1.0, -1.0, np.inf]), "non-finite"),
    (lambda f: f.assign(log2FC=[1.0, -1.0, np.nan]), "non-finite"),
    (lambda f: f.assign(**{ID: ["C1", "C2", "C1"]}), "duplicate"),
    (lambda f: f.assign(p_value=[0.1, np.nan, np.nan]), "measured row is missing"),
])
def test_frame_invariants_reject_bad_input(mutate, match):
    with pytest.raises(WalkerAdditionsError, match=match):
        wa.check_walker_frame(mutate(_good_frame()))


def test_building_rejects_a_core_addition_overlap(inputs_dir, cfg):
    bad = copy.deepcopy(cfg)
    core_acc = _read(inputs_dir / "ipa_input_full.csv")[ID].iloc[0]
    gene = _read(inputs_dir / "ipa_input_full.csv")["Gene names"].iloc[0]
    bad["assigned_proteins"] = [{"accession": core_acc, "gene": gene}]
    with pytest.raises(WalkerAdditionsError):
        wa.compute_walker(inputs_dir, bad)


def _single_row(**kw):
    row = {"accession": "X1", "gene": "Gx", "detected_in": "control_only",
           "Intensity 31578": np.nan, "Intensity 31580": np.nan,
           "Intensity 31579": 5.0, "Intensity 31581": 7.0}
    row.update(kw)
    return pd.DataFrame([row])


def _qual(direction="DOWN"):
    return pd.DataFrame({ID: ["X1"], "direction": [direction]})


def _verify(single, qual=None, gene="Gx"):
    return wa.verify_control_only(
        "X1", gene, single, qual if qual is not None else _qual(),
        design.control_columns(), design.treated_columns(), "control_only")


def test_verify_control_only_accepts_a_control_only_protein():
    assert _verify(_single_row())["accession"] == "X1"


@pytest.mark.parametrize("single,qual,gene,match", [
    (_single_row(**{"Intensity 31578": 3.0}), None, "Gx", "testosterone-side"),
    (_single_row(**{"Intensity 31581": np.nan}), None, "Gx", "both control readings"),
    (_single_row(**{"Intensity 31579": 0.0}), None, "Gx", "both control readings"),
    (_single_row(detected_in="treated_only"), None, "Gx", "NOT control-only"),
    (_single_row(), None, "Wrong", "config says"),
    (_single_row(), _qual("UP"), "Gx", "not a DOWN row"),
    (_single_row(accession="OTHER"), None, "Gx", "expected one row"),
])
def test_verify_control_only_rejects_everything_else(single, qual, gene, match):
    with pytest.raises(WalkerAdditionsError, match=match):
        _verify(single, qual, gene)


def test_a_protein_with_a_testosterone_reading_is_refused_by_the_builder(inputs_dir, cfg):
    path = inputs_dir / "single_condition_proteins.csv"
    single = _read(path)
    single.loc[single["accession"] == "P01132", "Intensity 31578"] = 123456.0
    single.to_csv(path, index=False)
    with pytest.raises(WalkerAdditionsError, match="testosterone-side"):
        wa.compute_walker(inputs_dir, cfg)


def test_measured_addition_must_exist_in_the_partial_tier(cfg):
    empty = pd.DataFrame(columns=[ID, "Gene names", "log2FC", "p_value",
                                  "adj_p_value", "n_imputed"])
    with pytest.raises(WalkerAdditionsError, match="exactly one partial-tier file"):
        wa.build_measured_rows(cfg, [(empty, "UP"), (empty, "DOWN")])


def test_the_hundred_fold_rule_is_not_applied_to_a_measured_row(cfg):
    partial = pd.DataFrame({ID: ["P97401"], "Gene names": ["Frzb"], "log2FC": [4.40696],
                            "p_value": [0.02], "adj_p_value": [0.99], "n_imputed": [1]})
    rows = wa.build_measured_rows(cfg, [(partial, "UP"), (partial.iloc[0:0], "DOWN")])
    assert rows.iloc[0]["log2FC"] == 4.40696
    assert rows.iloc[0]["value_source"] == "measured"


def test_load_config_rejects_a_missing_or_malformed_section(tmp_path):
    bad = tmp_path / "c.yaml"
    bad.write_text("organism: {}\n", encoding="utf-8")
    with pytest.raises(WalkerAdditionsError, match="walker_additions"):
        wa.load_config(bad)
    bad.write_text("walker_additions:\n  assigned_fold_change: {decrease_factor: 1}\n",
                   encoding="utf-8")
    with pytest.raises(WalkerAdditionsError):
        wa.load_config(bad)


# ---------------------------------------------------------------------------
# Provenance and labelling
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", [wa.WALKER_CSV_FILENAME, wa.WALKER_TXT_FILENAME,
                                  wa.WALKER_MANIFEST_FILENAME,
                                  wa.WALKER_NOASSIGNED_CSV_FILENAME,
                                  wa.WALKER_NOASSIGNED_TXT_FILENAME])
def test_sidecars_exist_say_assigned_and_match_the_file(results_dir, name):
    path = results_dir / name
    side = results_dir / f"{name}.provenance.json"
    if not side.exists():
        pytest.skip("sidecar not generated yet")
    rec = json.loads(side.read_text(encoding="utf-8"))
    assert rec["sha256"] == _sha(path)
    assert rec["assigned_values_are_not_measurements"] is True
    assert "not measured" in rec["assignment_note"]
    assert "below the detection limit" in rec["assignment_note"]
    assert rec["assigned_log2FC"] == ASSIGNED_LOG2FC
    assert [s["n"] for s in rec["recovered_rule_steps"]][0] == 862
    # D24: the committed rule is the strict one, so the last step admits nobody.
    assert rec["recovered_rule_steps"][-1]["n"] == 0
    assert rec["recovered_rule"]["require_all_replicates_seen"] is True
    assert rec["recovered_rule"]["require_no_imputed_values"] is True
    assert "STRICT" in rec["recovered_rule_note"] and "WIDENED" not in rec["recovered_rule_note"]
    assert "D24" in rec["recovered_rule_note"]
    # value_source_counts must describe THIS file, not the list the others describe.
    if name == wa.WALKER_MANIFEST_FILENAME:
        assert rec["value_source_counts"] == {"measured": 1, "assigned": 5}
    elif name.startswith("ipa_input_walker_without_assigned"):
        assert rec["value_source_counts"] == {"measured": 716}
        assert "assigned" not in rec["value_source_counts"]
        assert "without_assigned" in rec["value_source_counts_scope"]
        assert "without_assigned_note" in rec
    else:
        assert rec["value_source_counts"] == {"measured": 715 + 1, "assigned": 5}
    assert "added_measured_n_imputed_counts" not in rec


def test_manifest_labels_every_assigned_row_as_not_measured(manifest_df):
    assigned = manifest_df[manifest_df["value_source"] == "assigned"]
    assert len(assigned) == 5
    assert assigned["reason"].str.contains("NOT measured").all()
    assert assigned["reason"].str.contains("detection limit").all()
    frzb = manifest_df[manifest_df[ID] == "P97401"].iloc[0]
    assert frzb["value_source"] == "measured"
    assert "NOT applied" in frzb["reason"]


# ---------------------------------------------------------------------------
# Source-level policy checks, mirroring the neighbours'
# ---------------------------------------------------------------------------
def _code():
    text = Path(wa.__file__).read_text(encoding="utf-8")
    code = re.sub(r'"""(?:.|\n)*?"""', "", text)
    return re.sub(r"#.*", "", code)


def test_every_read_csv_asks_for_round_trip_float_precision():
    calls = re.findall(r"pd\.read_csv\((?:[^()]|\([^()]*\))*\)", _code())
    assert calls, "no pd.read_csv calls found -- did the module move?"
    missing = [c for c in calls if 'float_precision="round_trip"' not in c]
    assert not missing, f"these reads would perturb floats in the last ULP: {missing}"


def test_module_hardcodes_no_magic_numbers_or_accessions(frozen_counts):
    code = _code()
    forbidden = {str(frozen_counts[k]) for k in
                 ("ipa_input_rows", "n_up", "n_down", "n_ipa_walker")}
    forbidden |= {"0.585", "6.643856", "0.01", "100"}
    found = sorted(n for n in forbidden if re.search(rf"(?<![\w.]){re.escape(n)}(?![\w.])", code))
    assert not found, f"magic numbers in walker_additions.py: {found}"
    for acc in ("P01132", "Q61521", "O08542", "P55144", "Q9D1H9", "P97401"):
        assert acc not in code, f"accession {acc} hardcoded; it belongs in config.yaml"


def test_pipeline_stage_is_registered_and_pinned_to_frozen_counts(frozen_counts):
    import run_pipeline as rp

    stage = rp.STAGE_BY_ID["ipa_export_walker"]
    assert stage.script == "proteomics_de/export/walker_additions.py"
    assert {"foldchange", "ipa_export", "regulated_lists_supplementary"} <= set(stage.depends_on)
    paths = {o.path: o for o in stage.outputs}
    csv_out = paths["proteomics_de/results/ipa_input_walker.csv"]
    assert rp.resolve_expected_rows(csv_out, frozen_counts) == frozen_counts["n_ipa_walker"] == 721
    manifest_out = paths["proteomics_de/results/walker_additions_manifest.csv"]
    assert rp.resolve_expected_rows(manifest_out, frozen_counts) == 6
    # D24: the list without the five assigned rows is a registered, row-counted output too.
    na_out = paths["proteomics_de/results/ipa_input_walker_without_assigned.csv"]
    assert rp.resolve_expected_rows(na_out, frozen_counts) == frozen_counts[
        "n_ipa_walker_without_assigned"] == 716
    assert "proteomics_de/results/ipa_input_walker_without_assigned.txt" in paths


def test_new_outputs_are_frozen():
    entries = freeze.read_manifest(MANIFEST)
    for name in (wa.WALKER_CSV_FILENAME, wa.WALKER_TXT_FILENAME, wa.WALKER_MANIFEST_FILENAME,
                 wa.WALKER_NOASSIGNED_CSV_FILENAME, wa.WALKER_NOASSIGNED_TXT_FILENAME):
        assert f"proteomics_de/results/{name}" in entries, f"{name} is not in the manifest"


# ---------------------------------------------------------------------------
# Round 3: sidecar labelling is scoped per file, and FRZB is counted separately
# ---------------------------------------------------------------------------
def _facts(results_dir, which):
    cfg = wa.load_config()
    return wa.sidecar_facts(wa.compute_walker(results_dir, cfg), cfg, which)


def test_value_source_counts_are_scoped_to_the_file_each_sidecar_describes(results_dir, walker, manifest_df):
    lst, man = _facts(results_dir, "list"), _facts(results_dir, "manifest")
    assert lst["value_source_counts"] == walker["value_source"].value_counts().to_dict()
    assert man["value_source_counts"] == manifest_df["value_source"].value_counts().to_dict()
    assert lst["value_source_counts"] == {"measured": 716, "assigned": 5}
    assert man["value_source_counts"] == {"measured": 1, "assigned": 5}
    assert lst["value_source_counts"] != man["value_source_counts"]
    assert "walker.csv" in lst["value_source_counts_scope"]
    assert "manifest" in man["value_source_counts_scope"]


def test_the_committed_sidecars_carry_their_own_value_source_counts(results_dir, walker, manifest_df, without):
    for name, frame in ((wa.WALKER_CSV_FILENAME, walker), (wa.WALKER_TXT_FILENAME, walker),
                        (wa.WALKER_MANIFEST_FILENAME, manifest_df),
                        (wa.WALKER_NOASSIGNED_CSV_FILENAME, without),
                        (wa.WALKER_NOASSIGNED_TXT_FILENAME, without)):
        side = results_dir / f"{name}.provenance.json"
        if not side.exists():
            pytest.skip("sidecar not generated yet")
        rec = json.loads(side.read_text(encoding="utf-8"))
        assert rec["value_source_counts"] == frame["value_source"].value_counts().to_dict(), name
        assert sum(rec["value_source_counts"].values()) == rec["n_rows"], name


def test_an_unknown_sidecar_scope_is_refused(results_dir):
    cfg = wa.load_config()
    result = wa.compute_walker(results_dir, cfg)
    with pytest.raises(WalkerAdditionsError, match="scope"):
        wa.sidecar_facts(result, cfg, "everything")


def test_frzb_n_imputed_is_reported_separately_from_the_recovered_counts(results_dir, widened_cfg):
    cfg = wa.load_config()
    for which in ("list", "manifest", "without_assigned"):
        # committed (strict) config: FRZB only, no recovered protein.
        facts = wa.sidecar_facts(wa.compute_walker(results_dir, cfg), cfg, which)
        assert facts["measured_additions_n_imputed"] == {"P97401": 1}
        assert facts["recovered_n_imputed_counts"] == {}
        assert "added_measured_n_imputed_counts" not in facts
        # widened rule (tmp config): the 247 are counted on their own, FRZB is not among them.
        wide = wa.sidecar_facts(wa.compute_walker(results_dir, widened_cfg), widened_cfg, which)
        assert wide["measured_additions_n_imputed"] == {"P97401": 1}
        assert wide["recovered_n_imputed_counts"] == {"1": 178, "2": 69}
        assert sum(wide["recovered_n_imputed_counts"].values()) == 247


def test_the_sidecar_n_imputed_numbers_agree_with_the_manifest_and_the_docs(results_dir, manifest_df, frozen_counts, widened_cfg):
    rec = manifest_df[manifest_df[ID] != "P97401"]
    rec = rec[rec["value_source"] == "measured"]
    assert len(rec) == frozen_counts["n_walker_recovered"] == 0     # D24: none recovered
    for name in (wa.WALKER_CSV_FILENAME, wa.WALKER_MANIFEST_FILENAME):
        side = results_dir / f"{name}.provenance.json"
        if side.exists():
            r = json.loads(side.read_text(encoding="utf-8"))
            assert r["recovered_n_imputed_counts"] == {}
            assert r["measured_additions_n_imputed"] == {"P97401": 1}
    # The widened rule's manifest, from a tmp config: 178 with one imputed reading, 69 with two.
    wide = wa.compute_walker(results_dir, widened_cfg).manifest
    wide = wide[(wide[ID] != "P97401") & (wide["value_source"] == "measured")]
    assert wide["n_imputed"].astype(int).value_counts().to_dict() == {1: 178, 2: 69}
    assert len(wide) == 247


# ---------------------------------------------------------------------------
# D24: the 721 list, checked against the OTHER committed files, and the 716
# "without the five assigned" dataset for the IPA sensitivity run
# ---------------------------------------------------------------------------
def test_the_721_list_is_715_core_then_frzb_then_the_five_assigned(walker, results_dir, cfg, frozen_counts):
    """Independent of the builder: compare with ipa_input.csv and regulated_up_partial.csv."""
    assert len(walker) == 721
    ipa = pd.read_csv(results_dir / "ipa_input.csv", dtype=str, keep_default_na=False)
    got = pd.read_csv(results_dir / wa.WALKER_CSV_FILENAME, dtype=str, keep_default_na=False)
    shared = list(ipa.columns)
    assert got.iloc[:715][shared].reset_index(drop=True).equals(ipa[shared])    # rows 1-715
    up = pd.read_csv(results_dir / "regulated_up_partial.csv", dtype=str, keep_default_na=False)
    frzb = up[up[ID] == "P97401"].iloc[0]
    row = got.iloc[715]                                                           # row 716
    assert (row[ID], row["Gene names"], row["tier"], row["value_source"]) == (
        "P97401", "Frzb", "partial", "measured")
    for col in ("log2FC", "p_value", "adj_p_value"):
        assert row[col] == frzb[col], col                                         # real values
    tail = got.iloc[716:]                                                         # rows 717-721
    assert list(tail[ID]) == ["P01132", "Q61521", "O08542", "P55144", "Q9D1H9"]
    assert list(tail[ID]) == [e["accession"] for e in cfg["assigned_proteins"]]
    assert (tail["log2FC"] == "-6.643856").all() and (tail["regulated"] == "DOWN").all()
    assert (tail["p_value"] == "").all() and (tail["adj_p_value"] == "").all()
    assert (tail["value_source"] == "assigned").all() and (tail["tier"] == "qualitative").all()
    assert walker["regulated"].value_counts().to_dict() == {"UP": 510, "DOWN": 211}
    assert (frozen_counts["n_up"] + 1, frozen_counts["n_down"] + 5) == (510, 211)


def test_the_walker_manifest_has_six_rows_and_names_no_dropped_protein(manifest_df, cfg):
    assert list(manifest_df[ID]) == ["P97401"] + [e["accession"] for e in cfg["assigned_proteins"]]
    assert len(manifest_df) == 6
    text = " ".join(manifest_df["reason"].astype(str)).lower()
    assert "recovered" not in text and "widened" not in text
    assert "technical replicate" not in text and "one biological replicate" not in text


def test_the_without_assigned_list_is_the_721_minus_the_five_assigned(walker, without, frozen_counts):
    assert len(without) == frozen_counts["n_ipa_walker_without_assigned"] == 716
    assert len(without) == len(walker) - frozen_counts["n_walker_assigned"]
    assert list(without.columns) == wa.WALKER_COLS == list(walker.columns)
    pd.testing.assert_frame_equal(without, walker[walker["value_source"] != "assigned"].reset_index(drop=True))
    assert (without["value_source"] == "measured").all()
    assert without["regulated"].value_counts().to_dict() == {"UP": 510, "DOWN": 206}
    assert list(without[ID].iloc[:715]) == list(walker[ID].iloc[:715])
    assert without[ID].iloc[715] == "P97401"                   # the last row is FRZB
    assert not {"P01132", "Q61521", "O08542", "P55144", "Q9D1H9"} & set(without[ID])
    assert (without["adj_p_value"].notna() & without["p_value"].notna()).all()


def test_the_without_assigned_files_are_byte_for_byte_the_walker_rows_and_their_twins_agree(results_dir):
    full = (results_dir / wa.WALKER_CSV_FILENAME).read_text(encoding="utf-8").splitlines()
    no = (results_dir / wa.WALKER_NOASSIGNED_CSV_FILENAME).read_text(encoding="utf-8").splitlines()
    kept = [full[0]] + [ln for ln in full[1:] if not ln.endswith(",qualitative,assigned")]
    assert no == kept and len(no) == 717                        # header + 716
    txt = pd.read_csv(results_dir / wa.WALKER_NOASSIGNED_TXT_FILENAME, sep="\t", dtype=str,
                      keep_default_na=False)
    csv_ = pd.read_csv(results_dir / wa.WALKER_NOASSIGNED_CSV_FILENAME, dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(txt, csv_)
    for name in (wa.WALKER_NOASSIGNED_CSV_FILENAME, wa.WALKER_NOASSIGNED_TXT_FILENAME):
        raw = (results_dir / name).read_bytes()
        assert not raw.startswith(b"\xef\xbb\xbf") and b"\r" not in raw


def test_all_721_and_all_716_accessions_are_in_the_ipa_reference_set(walker, without, results_dir):
    background = set((results_dir / "ipa_background_measured.txt").read_text(encoding="utf-8").split("\n")[1:]) - {""}
    assert len(background) == 2552
    assert set(walker[ID]) <= background and len(set(walker[ID])) == 721
    assert set(without[ID]) <= background and len(set(without[ID])) == 716


def test_the_pipeline_check_for_counts_is_consistent_for_the_new_key(frozen_counts):
    assert frozen_counts["n_ipa_walker_without_assigned"] == (
        frozen_counts["n_ipa_walker"] - frozen_counts["n_walker_assigned"])
    assert frozen_counts["n_ipa_walker"] == (
        frozen_counts["ipa_input_rows"] + frozen_counts["n_walker_measured_added"]
        + frozen_counts["n_walker_assigned"] + frozen_counts["n_walker_recovered"])
