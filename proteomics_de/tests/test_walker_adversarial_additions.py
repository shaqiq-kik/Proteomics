"""Adversarial tests for :mod:`proteomics_de.export.walker_additions` (D20).

Written by an independent tester, deliberately NOT sharing code with
``test_walker_additions.py``. Every test works in ``tmp_path`` copies of the
committed inputs (or on synthetic frames); nothing here writes to the committed
``results/`` directory.

Round 2: the two former ``DEFECT-CANDIDATE`` xfails (assigned |log2FC| below the
regulation cutoff; a failed validation leaving deliverables behind) are fixed in
the source and now run as ordinary tests. The recovered-protein tests encode the
user's WIDENED rule (partial-tier proteins with |log2FC| >= 0.585 are added even
with imputed values); the strict rule is exercised by flipping the two config
switches ``require_all_replicates_seen`` / ``require_no_imputed_values``.
"""

from __future__ import annotations

import builtins
import copy
import csv
import hashlib
import math
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.config import design  # noqa: E402
from proteomics_de.export import walker_additions as wa  # noqa: E402
from proteomics_de.export.walker_additions import WalkerAdditionsError  # noqa: E402

ID = "UniProt Accession Number"
CTRL = design.control_columns()
TRT = design.treated_columns()
READINGS = CTRL + TRT
OUTPUT_NAMES = (wa.WALKER_CSV_FILENAME, wa.WALKER_TXT_FILENAME, wa.WALKER_MANIFEST_FILENAME)
REAL_PROTECTED = ("ipa_input.csv", "ipa_input_extended.csv", "ipa_input_extended.txt",
                  "ipa_input_full.csv", "ipa_input_full.txt")


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------
def _read(path):
    return pd.read_csv(path, float_precision="round_trip")


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def cfg():
    return wa.load_config()


@pytest.fixture
def inputs_dir(results_dir, tmp_path):
    """Scratch results dir holding every file the module reads (and the frozen ones)."""
    d = tmp_path / "res"
    d.mkdir()
    for name in set(wa._INPUTS) | set(wa.BYTE_FROZEN_INPUTS):
        src = results_dir / name
        if src.exists():
            shutil.copy(src, d / name)
    return d


@pytest.fixture(autouse=True)
def _real_results_never_touched(results_dir):
    """Whatever a test does, the committed originals must be byte-identical after."""
    names = sorted(set(REAL_PROTECTED) | set(OUTPUT_NAMES) | set(wa.BYTE_FROZEN_INPUTS))
    before = {n: _sha(results_dir / n) for n in names if (results_dir / n).exists()}
    yield
    after = {n: _sha(results_dir / n) for n in before}
    assert before == after, [n for n in before if before[n] != after[n]]


def _cfg_with(cfg, **sections):
    out = copy.deepcopy(cfg)
    out.update(sections)
    return out


def _strict(cfg):
    """The original strict rule, restored purely through the config switches."""
    c = copy.deepcopy(cfg)
    c["recovered_rule"].update(require_all_replicates_seen=True,
                               require_no_imputed_values=True)
    return c


def _recovered(res):
    """Recovered rows = additions beyond the 715 core, FRZB and the assigned five."""
    n_fixed = res.counts["core"] + res.counts["measured"] + res.counts["assigned"]
    return res.frame.iloc[n_fixed:]


def _assigned_only(cfg, accession, gene):
    return _cfg_with(cfg, assigned_proteins=[{"accession": accession, "gene": gene}])


def _no_outputs(directory):
    return [n for n in OUTPUT_NAMES if (Path(directory) / n).exists()]


def _control_only_row(accession, gene, ctrl=(5.0, 7.0), trt=(np.nan, np.nan)):
    row = {"accession": accession, "gene": gene, "detected_in": "control_only"}
    row.update(dict(zip(CTRL, ctrl)))
    row.update(dict(zip(TRT, trt)))
    return row


def _append(path, rows):
    df = _read(path)
    df = pd.concat([df, pd.DataFrame(rows)], ignore_index=True)
    df.to_csv(path, index=False)


def _make_clean(inputs_dir, accession, *, n_imputed=0, log2fc=None, reading=1e6,
                file="regulated_up_partial.csv"):
    """Make partial-tier protein `accession` look like a clean, fully-seen protein."""
    p = inputs_dir / file
    partial = _read(p)
    mask = partial[ID] == accession
    assert mask.sum() == 1
    partial.loc[mask, "n_imputed"] = n_imputed
    if log2fc is not None:
        partial.loc[mask, "log2FC"] = log2fc
    partial.to_csv(p, index=False)
    q = inputs_dir / "qc_limma.csv"
    limma = _read(q)
    for col in READINGS:
        limma.loc[limma["id"] == accession, col] = reading
    limma.to_csv(q, index=False)


def _partial_ids(inputs_dir, file="regulated_up_partial.csv"):
    ids = list(_read(inputs_dir / file)[ID])
    return [i for i in ids if i != "P97401"]


# ---------------------------------------------------------------------------
# 1. The assigned log2FC follows the configured factor
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("factor", [2, 10, 100, 1000, 3.5])
def test_assigned_value_follows_the_configured_factor(inputs_dir, cfg, factor):
    c = copy.deepcopy(cfg)
    c["assigned_fold_change"]["decrease_factor"] = factor
    res = wa.compute_walker(inputs_dir, c)
    acc = [e["accession"] for e in c["assigned_proteins"]]
    rows = res.frame[res.frame[ID].isin(acc)]
    assert len(rows) == len(acc)
    expected = round(-math.log2(factor), 6)
    assert (rows["log2FC"] == expected).all()
    assert (rows["regulated"] == "DOWN").all()
    assert (rows["value_source"] == "assigned").all()
    assert rows[["p_value", "adj_p_value"]].isna().all().all()
    assert np.allclose(2.0 ** rows["log2FC"], 1.0 / factor, rtol=1e-5)
    # the manifest reason must quote the factor actually used, not a stale 100
    reasons = res.manifest.loc[res.manifest[ID].isin(acc), "reason"]
    assert reasons.str.contains(f"{factor:g}-fold decrease").all()
    # FRZB is untouched by the factor
    frzb = res.frame[res.frame[ID] == "P97401"].iloc[0]
    assert frzb["value_source"] == "measured" and frzb["log2FC"] != expected


def test_factor_change_through_a_config_file_reaches_the_written_files(inputs_dir):
    full = yaml.safe_load(wa._PKG_DIR.joinpath("config", "config.yaml").read_text(encoding="utf-8"))
    full["walker_additions"]["assigned_fold_change"]["decrease_factor"] = 10
    conf = inputs_dir / "alt_config.yaml"
    conf.write_text(yaml.safe_dump(full), encoding="utf-8")
    wa.build_walker_from_results(inputs_dir, config_path=conf)   # validate=True: counts unaffected
    for name, delim in ((wa.WALKER_CSV_FILENAME, ","), (wa.WALKER_TXT_FILENAME, "\t")):
        with open(inputs_dir / name, newline="", encoding="utf-8") as fh:
            rows = [r for r in csv.DictReader(fh, delimiter=delim) if r["value_source"] == "assigned"]
        assert len(rows) == 5
        assert {r["log2FC"] for r in rows} == {"-3.321928"}
        assert all(r["p_value"] == "" and r["adj_p_value"] == "" for r in rows)


@pytest.mark.parametrize("bad", [1, 0.5, 0, -10, "100", None, True, float("nan")])
def test_load_config_rejects_a_nonsense_factor(tmp_path, bad):
    full = yaml.safe_load(wa._PKG_DIR.joinpath("config", "config.yaml").read_text(encoding="utf-8"))
    full["walker_additions"]["assigned_fold_change"]["decrease_factor"] = bad
    conf = tmp_path / "c.yaml"
    conf.write_text(yaml.safe_dump(full), encoding="utf-8")
    with pytest.raises(WalkerAdditionsError):
        wa.load_config(conf)


def test_an_infinite_factor_never_yields_a_row(inputs_dir, cfg):
    """inf passes load_config's `> 1`; something downstream must still refuse it."""
    c = copy.deepcopy(cfg)
    c["assigned_fold_change"]["decrease_factor"] = float("inf")
    with pytest.raises(WalkerAdditionsError):
        wa.compute_walker(inputs_dir, c)


def test_direction_up_is_refused_because_a_decrease_cannot_be_up(inputs_dir, cfg):
    c = copy.deepcopy(cfg)
    c["assigned_fold_change"]["direction"] = "UP"
    with pytest.raises(WalkerAdditionsError, match="inconsistent with direction"):
        wa.compute_walker(inputs_dir, c)


def test_a_factor_whose_log2fc_is_below_the_regulation_cutoff_is_refused(inputs_dir, cfg):
    """decrease_factor=1.2 (log2FC -0.263) must not become a 'regulated DOWN' row."""
    c = copy.deepcopy(cfg)
    c["assigned_fold_change"]["decrease_factor"] = 1.2
    with pytest.raises(WalkerAdditionsError, match="regulation cutoff"):
        wa.compute_walker(inputs_dir, c)
    assert _no_outputs(inputs_dir) == []


@pytest.mark.parametrize("factor,ok", [(1.4, False), (1.5, False), (1.51, True), (2, True)])
def test_the_assigned_cutoff_check_sits_at_the_pipeline_threshold(inputs_dir, cfg, factor, ok):
    """log2(1.5) = 0.58496 rounds to 0.584963 < 0.585, so 1.5 is refused; 1.51 clears it."""
    c = copy.deepcopy(cfg)
    c["assigned_fold_change"]["decrease_factor"] = factor
    if ok:
        res = wa.compute_walker(inputs_dir, c)
        assert (res.frame.loc[res.frame["value_source"] == "assigned", "log2FC"].abs()
                >= 0.585).all()
    else:
        with pytest.raises(WalkerAdditionsError, match="regulation cutoff"):
            wa.compute_walker(inputs_dir, c)


# ---------------------------------------------------------------------------
# 2. An accession that is NOT control-only must raise, never get a value
# ---------------------------------------------------------------------------
def _case_treated_only(d):
    return "A2A9R3", "Mageb4"                      # detected_in=treated_only


def _case_one_control_zero(d):
    return "A1A547", "Pglyrp3"                     # control 31579 == 0.0


def _case_core_protein_not_in_single(d):
    full = _read(d / "ipa_input_full.csv")
    return full[ID].iloc[0], full["Gene names"].iloc[0]


def _case_flipped_in_qualitative(d):
    q = _read(d / "qualitative_changes.csv")
    q.loc[q[ID] == "P01132", "direction"] = "UP"
    q.to_csv(d / "qualitative_changes.csv", index=False)
    return "P01132", "Egf"


def _case_missing_from_qualitative(d):
    q = _read(d / "qualitative_changes.csv")
    q[q[ID] != "P01132"].to_csv(d / "qualitative_changes.csv", index=False)
    return "P01132", "Egf"


def _set_single(d, accession, **cols):
    s = _read(d / "single_condition_proteins.csv")
    for k, v in cols.items():
        s.loc[s["accession"] == accession, k] = v
    s.to_csv(d / "single_condition_proteins.csv", index=False)


def _case_tiny_treated_reading(d):
    _set_single(d, "P01132", **{TRT[0]: 1e-9})
    return "P01132", "Egf"


def _case_zero_and_real_treated_reading(d):
    _set_single(d, "P01132", **{TRT[0]: 0.0, TRT[1]: 1.0})
    return "P01132", "Egf"


def _case_label_says_treated(d):
    _set_single(d, "P01132", detected_in="treated_only")
    return "P01132", "Egf"


def _case_wrong_gene(d):
    return "P01132", "NotEgf"


def _case_unknown_accession(d):
    return "Z9Z9Z9", "Nope"


@pytest.mark.parametrize("case", [
    _case_treated_only, _case_one_control_zero, _case_core_protein_not_in_single,
    _case_flipped_in_qualitative, _case_missing_from_qualitative, _case_tiny_treated_reading,
    _case_zero_and_real_treated_reading, _case_label_says_treated, _case_wrong_gene,
    _case_unknown_accession,
], ids=lambda f: f.__name__[6:])
def test_a_non_control_only_accession_is_refused_and_nothing_is_written(inputs_dir, cfg, case):
    acc, gene = case(inputs_dir)
    c = _assigned_only(cfg, acc, gene)
    with pytest.raises(WalkerAdditionsError):
        wa.compute_walker(inputs_dir, c)
    # through the public builder too: it must fail BEFORE writing a deliverable
    conf = inputs_dir / "c.yaml"
    full = yaml.safe_load(wa._PKG_DIR.joinpath("config", "config.yaml").read_text(encoding="utf-8"))
    full["walker_additions"]["assigned_proteins"] = [{"accession": acc, "gene": gene}]
    conf.write_text(yaml.safe_dump(full), encoding="utf-8")
    with pytest.raises(WalkerAdditionsError):
        wa.build_walker_from_results(inputs_dir, config_path=conf, validate=False)
    assert _no_outputs(inputs_dir) == []


def test_frzb_listed_as_assigned_is_refused_and_never_gets_the_assigned_value(inputs_dir, cfg):
    c = _assigned_only(cfg, "P97401", "Frzb")
    with pytest.raises(WalkerAdditionsError):
        wa.compute_walker(inputs_dir, c)


def test_every_assigned_row_in_the_output_is_control_only_in_the_raw_data(inputs_dir, cfg):
    """The converse guarantee: whatever carries value_source=assigned is provably control-only."""
    res = wa.compute_walker(inputs_dir, cfg)
    single = _read(inputs_dir / "single_condition_proteins.csv").set_index("accession")
    assigned = res.frame[res.frame["value_source"] == "assigned"]
    assert len(assigned) == len(cfg["assigned_proteins"])
    for acc in assigned[ID]:
        row = single.loc[acc]
        assert row["detected_in"] == "control_only"
        assert all(row[c] > 0 for c in CTRL) and not any(row[c] > 0 for c in TRT)


# ---------------------------------------------------------------------------
# 3. Duplicates: core overlap, config duplicates, measured/assigned collision
# ---------------------------------------------------------------------------
def test_an_assigned_protein_already_in_the_core_is_rejected_by_the_overlap_guard(inputs_dir, cfg):
    """Make a CORE protein look control-only everywhere so only the overlap guard can stop it."""
    full = _read(inputs_dir / "ipa_input_full.csv")
    acc, gene = full[ID].iloc[0], full["Gene names"].iloc[0]
    _append(inputs_dir / "single_condition_proteins.csv", [_control_only_row(acc, gene)])
    _append(inputs_dir / "qualitative_changes.csv", [{
        "Gene names": gene, ID: acc, "direction": "DOWN", "basis": "control_only",
        **{c: 1.0 for c in CTRL}, **{c: np.nan for c in TRT}, "note": "x"}])
    with pytest.raises(WalkerAdditionsError, match="already in the core"):
        wa.compute_walker(inputs_dir, _assigned_only(cfg, acc, gene))


def test_a_measured_addition_already_in_the_core_is_rejected_by_the_overlap_guard(inputs_dir, cfg):
    full = _read(inputs_dir / "ipa_input_full.csv")
    row = full.iloc[0]
    _append(inputs_dir / "regulated_up_partial.csv", [{
        "Gene names": row["Gene names"], ID: row[ID], "log2FC": 2.0, "fold_change": 4.0,
        "p_value": 0.01, "adj_p_value": 0.5, "n_imputed": 1}])
    c = _cfg_with(cfg, measured_additions=[{"accession": row[ID], "gene": row["Gene names"]}])
    with pytest.raises(WalkerAdditionsError, match="already in the core"):
        wa.compute_walker(inputs_dir, c)


def test_the_same_protein_both_measured_and_assigned_is_rejected(inputs_dir, cfg):
    gene = "Frzb"
    _append(inputs_dir / "single_condition_proteins.csv", [_control_only_row("P97401", gene)])
    _append(inputs_dir / "qualitative_changes.csv", [{
        "Gene names": gene, ID: "P97401", "direction": "DOWN", "basis": "control_only",
        **{c: 1.0 for c in CTRL}, **{c: np.nan for c in TRT}, "note": "x"}])
    c = _assigned_only(cfg, "P97401", gene)
    with pytest.raises(WalkerAdditionsError, match="both measured-added and assigned"):
        wa.compute_walker(inputs_dir, c)


def test_a_duplicated_config_entry_never_produces_a_duplicate_row(inputs_dir, cfg):
    for key, entry in (("assigned_proteins", cfg["assigned_proteins"][0]),
                       ("measured_additions", cfg["measured_additions"][0])):
        c = copy.deepcopy(cfg)
        c[key] = c[key] + [dict(entry)]
        with pytest.raises(WalkerAdditionsError, match="duplicate"):
            wa.compute_walker(inputs_dir, c)


def test_a_clean_frzb_is_not_added_a_second_time_as_recovered(inputs_dir, cfg):
    _make_clean(inputs_dir, "P97401")
    res = wa.compute_walker(inputs_dir, cfg)
    assert "P97401" not in set(_recovered(res)[ID])
    assert (res.frame[ID] == "P97401").sum() == 1
    row = res.frame[res.frame[ID] == "P97401"].iloc[0]
    assert row["tier"] == "partial" and row["value_source"] == "measured"
    assert res.frame[ID].is_unique and res.manifest[ID].is_unique


def test_a_protein_in_both_partial_files_is_a_loud_error(inputs_dir, cfg):
    up = _read(inputs_dir / "regulated_up_partial.csv")
    down = _read(inputs_dir / "regulated_down_partial.csv")
    dup = up[up[ID] != "P97401"].iloc[[0]].copy()
    dup["log2FC"] = -abs(dup["log2FC"])
    pd.concat([down, dup], ignore_index=True).to_csv(
        inputs_dir / "regulated_down_partial.csv", index=False)
    with pytest.raises(WalkerAdditionsError, match="duplicate"):
        wa.compute_walker(inputs_dir, cfg)


# ---------------------------------------------------------------------------
# 4. Missing / malformed inputs
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("name", wa._INPUTS)
def test_a_missing_input_file_is_a_clear_error_naming_the_file(inputs_dir, name):
    (inputs_dir / name).unlink()
    with pytest.raises(FileNotFoundError, match=name):
        wa.build_walker_from_results(inputs_dir)
    assert _no_outputs(inputs_dir) == []


def test_a_nonexistent_or_empty_results_dir_is_a_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="ipa_input_full.csv"):
        wa.build_walker_from_results(tmp_path / "does_not_exist")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        wa.build_walker_from_results(empty)
    assert _no_outputs(empty) == []


def test_a_missing_config_file_is_an_error_not_a_silent_default(tmp_path, inputs_dir):
    with pytest.raises(FileNotFoundError):
        wa.build_walker_from_results(inputs_dir, config_path=tmp_path / "nope.yaml")
    assert _no_outputs(inputs_dir) == []


@pytest.mark.parametrize("drop", [
    ("qc_limma.csv", "id"),
    ("single_condition_proteins.csv", "detected_in"),
    ("qualitative_changes.csv", "direction"),
    ("regulated_up_partial.csv", "n_imputed"),
    ("ipa_input_full.csv", "log2FC"),
])
def test_a_malformed_input_fails_before_any_deliverable_is_written(inputs_dir, drop):
    name, col = drop
    df = _read(inputs_dir / name)
    df.drop(columns=[col]).to_csv(inputs_dir / name, index=False)
    with pytest.raises(Exception):
        wa.build_walker_from_results(inputs_dir, validate=False)
    assert _no_outputs(inputs_dir) == []


def test_an_empty_assigned_list_in_config_is_rejected(tmp_path):
    full = yaml.safe_load(wa._PKG_DIR.joinpath("config", "config.yaml").read_text(encoding="utf-8"))
    full["walker_additions"]["assigned_proteins"] = []
    conf = tmp_path / "c.yaml"
    conf.write_text(yaml.safe_dump(full), encoding="utf-8")
    with pytest.raises(WalkerAdditionsError, match="non-empty"):
        wa.load_config(conf)


# ---------------------------------------------------------------------------
# 5. The recovered rule, end to end through the real file readers
# ---------------------------------------------------------------------------
def _widened_n(res):
    return res.counts["recovered"]


def test_a_clean_recovered_protein_is_included_as_measured(inputs_dir, cfg):
    """A fully clean partial row passes both the widened and the strict rule."""
    base = wa.compute_walker(inputs_dir, cfg).counts["recovered"]
    target = _partial_ids(inputs_dir)[0]
    _make_clean(inputs_dir, target)
    for conf, expect in ((cfg, base), (_strict(cfg), 1)):
        res = wa.compute_walker(inputs_dir, conf)
        assert res.counts["recovered"] == expect
        row = res.frame[res.frame[ID] == target].iloc[0]
        assert (row["tier"], row["value_source"], row["regulated"]) == ("partial", "measured", "UP")
        assert pd.notna(row["p_value"]) and pd.notna(row["adj_p_value"])
        src = _read(inputs_dir / "regulated_up_partial.csv").set_index(ID).loc[target]
        assert row["log2FC"] == src["log2FC"] and row["p_value"] == src["p_value"]
        m = res.manifest[res.manifest[ID] == target].iloc[0]
        assert m["value_source"] == "measured" and m["tier"] == "partial"
        assert m["n_imputed"] == 0
        assert (res.frame["value_source"] == "assigned").sum() == 5  # assigned set unaffected
        assert res.counts["total"] == res.counts["core"] + 6 + expect


def test_the_widened_rule_adds_every_partial_protein_that_clears_the_cutoff(inputs_dir, cfg):
    res = wa.compute_walker(inputs_dir, cfg)
    rec = _recovered(res)
    assert len(rec) == res.counts["recovered"] == 247
    assert res.counts["total"] == 715 + 1 + 5 + 247 == 968
    assert (rec["tier"] == "partial").all() and (rec["value_source"] == "measured").all()
    assert (rec["log2FC"].abs() >= 0.585).all()
    assert rec[["p_value", "adj_p_value"]].notna().all().all()
    assert "P97401" not in set(rec[ID])
    m = res.manifest.set_index(ID).loc[rec[ID]]
    assert set(m["n_imputed"].astype(int)) == {1, 2}
    assert (m["n_imputed"] >= 1).all()
    assert m["reason"].str.contains("WIDENED").all()
    assert m["reason"].str.contains("no problems").all()
    assert m["reason"].str.contains("imputed").all()
    assert [n for _l, n in res.steps] == [862, 856, 247, 247]


def test_the_strict_rule_is_restored_by_flipping_the_config_switches(inputs_dir, cfg):
    strict = wa.compute_walker(inputs_dir, _strict(cfg))
    assert strict.counts["recovered"] == 0 and strict.counts["total"] == 721
    assert [n for _l, n in strict.steps] == [862, 856, 247, 247, 0, 0]
    for switch in ("require_all_replicates_seen", "require_no_imputed_values"):
        one = copy.deepcopy(cfg)
        one["recovered_rule"][switch] = True
        assert wa.compute_walker(inputs_dir, one).counts["recovered"] == 0, switch


@pytest.mark.parametrize("label,kwargs,limma_edit", [
    ("log2FC just under 0.585", {"log2fc": 0.5849}, None),
], ids=["below_cut"])
def test_a_recovered_candidate_with_a_problem_is_excluded(inputs_dir, cfg, label, kwargs, limma_edit):
    """Under the WIDENED rule only the effect-size condition still excludes a partial row."""
    target = _partial_ids(inputs_dir)[0]
    base = wa.compute_walker(inputs_dir, cfg)
    assert target in set(base.frame[ID])
    _make_clean(inputs_dir, target, **kwargs)
    res = wa.compute_walker(inputs_dir, cfg)
    assert target not in set(res.frame[ID]), label
    assert res.counts["recovered"] == base.counts["recovered"] - 1


@pytest.mark.parametrize("label,kwargs,limma_edit", [
    ("n_imputed=1", {"n_imputed": 1}, None),
    ("log2FC just under 0.585", {"log2fc": 0.5849}, None),
    ("one reading NaN", {}, ("nan", 0)),
    ("one reading zero", {}, ("zero", 3)),
    ("whole treated group missing", {}, ("nan_group", None)),
], ids=["imputed", "below_cut", "reading_nan", "reading_zero", "group_missing"])
def test_under_the_strict_switches_a_candidate_with_a_problem_is_excluded(
        inputs_dir, cfg, label, kwargs, limma_edit):
    target = _partial_ids(inputs_dir)[0]
    _make_clean(inputs_dir, target, **kwargs)
    if limma_edit:
        kind, idx = limma_edit
        q = inputs_dir / "qc_limma.csv"
        limma = _read(q)
        m = limma["id"] == target
        if kind == "nan":
            limma.loc[m, READINGS[idx]] = np.nan
        elif kind == "zero":
            limma.loc[m, READINGS[idx]] = 0.0
        else:
            for c in TRT:
                limma.loc[m, c] = np.nan
        limma.to_csv(q, index=False)
    res = wa.compute_walker(inputs_dir, _strict(cfg))
    assert target not in set(res.frame[ID]), label
    assert res.counts["recovered"] == 0


def test_the_log2fc_boundary_is_inclusive_end_to_end(inputs_dir, cfg):
    target = _partial_ids(inputs_dir)[0]
    base = wa.compute_walker(inputs_dir, cfg).counts["recovered"]
    _make_clean(inputs_dir, target, log2fc=0.585)
    res = wa.compute_walker(inputs_dir, cfg)
    assert res.counts["recovered"] == base and target in set(res.frame[ID])
    _make_clean(inputs_dir, target, log2fc=0.5849)
    res = wa.compute_walker(inputs_dir, cfg)
    assert res.counts["recovered"] == base - 1 and target not in set(res.frame[ID])


def test_a_clean_down_protein_is_recovered_as_down_with_its_own_sign(inputs_dir, cfg):
    down_ids = _partial_ids(inputs_dir, "regulated_down_partial.csv")
    target = down_ids[0]
    _make_clean(inputs_dir, target, file="regulated_down_partial.csv", log2fc=-1.7)
    res = wa.compute_walker(inputs_dir, cfg)
    row = res.frame[res.frame[ID] == target].iloc[0]
    assert (row["regulated"], row["log2FC"], row["value_source"]) == ("DOWN", -1.7, "measured")


def test_a_label_that_contradicts_the_sign_is_a_loud_error_not_a_silent_row(inputs_dir, cfg):
    target = _partial_ids(inputs_dir)[0]                       # in the UP file ...
    _make_clean(inputs_dir, target, log2fc=-2.0)               # ... but negative
    with pytest.raises(WalkerAdditionsError, match="disagrees with the sign"):
        wa.compute_walker(inputs_dir, cfg)


def test_a_clean_candidate_without_p_values_is_a_loud_error(inputs_dir, cfg):
    target = _partial_ids(inputs_dir)[0]
    _make_clean(inputs_dir, target)
    p = inputs_dir / "regulated_up_partial.csv"
    df = _read(p)
    df.loc[df[ID] == target, "p_value"] = np.nan
    df.to_csv(p, index=False)
    with pytest.raises(WalkerAdditionsError, match="missing its p-value"):
        wa.compute_walker(inputs_dir, cfg)


def test_a_qualitative_protein_can_never_be_recovered_even_with_all_readings(inputs_dir, cfg):
    q = _read(inputs_dir / "qualitative_changes.csv")
    for c in READINGS:
        q[c] = 1e6                                             # "seen" everywhere
    q.to_csv(inputs_dir / "qualitative_changes.csv", index=False)
    for conf in (cfg, _strict(cfg)):
        res = wa.compute_walker(inputs_dir, conf)
        assigned = set(res.frame.loc[res.frame["value_source"] == "assigned", ID])
        # The five assigned proteins are the only qualitative rows on the list.
        assert (set(q[ID]) & set(res.frame[ID])) == assigned
        assert not set(q[ID]) & set(_recovered(res)[ID])


def test_recovered_rows_are_ordered_by_abs_log2fc_with_stable_ties(inputs_dir, cfg):
    ups = _partial_ids(inputs_dir)[:2]
    down = _partial_ids(inputs_dir, "regulated_down_partial.csv")[0]
    _make_clean(inputs_dir, ups[0], log2fc=20.0)
    _make_clean(inputs_dir, ups[1], log2fc=20.0)               # exact tie with ups[0]
    _make_clean(inputs_dir, down, file="regulated_down_partial.csv", log2fc=-30.0)
    res = wa.compute_walker(inputs_dir, cfg)
    rec = _recovered(res)
    assert list(rec[ID].iloc[:3]) == [down, ups[0], ups[1]]
    abs_lfc = list(rec["log2FC"].abs())
    assert abs_lfc == sorted(abs_lfc, reverse=True)
    # and the answer does not depend on how many times we ask
    assert list(_recovered(wa.compute_walker(inputs_dir, cfg))[ID]) == list(rec[ID])


def test_the_n_imputed_cap_is_not_a_config_key_any_more_only_the_two_switches(cfg):
    assert set(cfg["recovered_rule"]) == {
        "min_abs_log2fc", "require_partial_tier",
        "require_all_replicates_seen", "require_no_imputed_values"}
    assert cfg["recovered_rule"]["require_all_replicates_seen"] is False
    assert cfg["recovered_rule"]["require_no_imputed_values"] is False
    assert cfg["recovered_rule"]["require_partial_tier"] is True


def test_the_widened_rule_adds_only_measured_rows_and_never_touches_the_assigned_set(inputs_dir, cfg):
    res = wa.compute_walker(inputs_dir, cfg)
    rec = res.frame[res.frame["tier"] == "partial"].iloc[1:]    # after FRZB
    assert len(rec) == res.counts["recovered"] > 0
    assert (rec["value_source"] == "measured").all()
    assert rec[["p_value", "adj_p_value"]].notna().all().all()
    assert (rec["log2FC"].abs() >= 0.585).all()
    assert "P97401" not in set(rec[ID])
    assert res.frame[ID].is_unique
    assert (res.frame["value_source"] == "assigned").sum() == 5
    assert len(res.manifest) == res.counts["additions_listed"] == 6 + res.counts["recovered"]


def _counts_file(tmp_path, **overrides):
    import json
    counts = wa._frozen_counts()
    counts.update(overrides)
    path = tmp_path / "counts.json"
    path.write_text(json.dumps(counts), encoding="utf-8")
    return path


def test_recovered_rows_failing_the_frozen_count_check_raise(inputs_dir, tmp_path):
    counts = _counts_file(tmp_path, n_walker_recovered=0)
    with pytest.raises(WalkerAdditionsError, match="recovered"):
        wa.build_walker_from_results(inputs_dir, counts_path=counts)   # data says 247


def test_a_failed_validation_leaves_no_deliverable_behind(inputs_dir, tmp_path):
    """A rejected count check must not leave a deliverable (or scratch dir) behind."""
    counts = _counts_file(tmp_path, n_walker_recovered=0)
    with pytest.raises(WalkerAdditionsError):
        wa.build_walker_from_results(inputs_dir, counts_path=counts)
    assert _no_outputs(inputs_dir) == []
    assert not [p for p in inputs_dir.iterdir() if p.name.startswith(".walker_stage_")]


@pytest.mark.parametrize("failing", ["validate_ipa", "assert_twins_agree",
                                     "assert_core_rows_are_verbatim"])
def test_a_failure_after_staging_leaves_no_deliverable_and_no_scratch_dir(
        inputs_dir, monkeypatch, failing):
    """The files are staged and validated BEFORE any rename into the official names."""
    def boom(*a, **kw):
        raise WalkerAdditionsError(f"injected {failing} failure")

    monkeypatch.setattr(wa, failing, boom)
    with pytest.raises(WalkerAdditionsError, match="injected"):
        wa.build_walker_from_results(inputs_dir)
    assert _no_outputs(inputs_dir) == []
    assert not [p for p in inputs_dir.iterdir() if p.name.startswith(".walker_stage_")]


def test_a_failed_rebuild_does_not_clobber_the_previous_good_deliverables(inputs_dir, tmp_path):
    wa.build_walker_from_results(inputs_dir)
    before = {n: (inputs_dir / n).read_bytes() for n in OUTPUT_NAMES}
    counts = _counts_file(tmp_path, n_walker_recovered=0)
    with pytest.raises(WalkerAdditionsError):
        wa.build_walker_from_results(inputs_dir, counts_path=counts)
    assert before == {n: (inputs_dir / n).read_bytes() for n in OUTPUT_NAMES}


def test_a_successful_build_leaves_no_scratch_dir(inputs_dir):
    wa.build_walker_from_results(inputs_dir)
    assert not [p for p in inputs_dir.iterdir() if p.name.startswith(".walker_stage_")]


# ---------------------------------------------------------------------------
# 6. Assigned rows never carry p-values; FRZB keeps real ones
# ---------------------------------------------------------------------------
def test_planted_p_value_columns_in_the_single_file_never_reach_an_assigned_row(inputs_dir, cfg):
    s = _read(inputs_dir / "single_condition_proteins.csv")
    s["p_value"] = 0.01
    s["adj_p_value"] = 0.02
    s.to_csv(inputs_dir / "single_condition_proteins.csv", index=False)
    wa.build_walker_from_results(inputs_dir, validate=True)
    with open(inputs_dir / wa.WALKER_CSV_FILENAME, newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh) if r["value_source"] == "assigned"]
    assert len(rows) == 5
    assert all(r["p_value"] == "" and r["adj_p_value"] == "" for r in rows)


@pytest.mark.parametrize("col", ["p_value", "adj_p_value"])
def test_a_single_planted_p_value_on_an_assigned_row_is_rejected(col):
    core = pd.DataFrame({ID: ["C1"], "Gene names": ["a"], "log2FC": [1.0],
                         "regulated": ["UP"], "p_value": [0.1], "adj_p_value": [0.9]})
    add = pd.DataFrame({ID: ["A1"], "Gene names": ["x"], "log2FC": [-6.643856],
                        "regulated": ["DOWN"], "p_value": [np.nan], "adj_p_value": [np.nan],
                        "tier": ["qualitative"], "value_source": ["assigned"]})
    add[col] = 0.5
    with pytest.raises(WalkerAdditionsError, match="assigned row carries a p-value"):
        wa.build_walker_frame(core, [add])


def test_an_assigned_label_cannot_be_dropped_by_the_writers(inputs_dir):
    wa.build_walker_from_results(inputs_dir)
    for name in (wa.WALKER_CSV_FILENAME, wa.WALKER_TXT_FILENAME):
        df = pd.read_csv(inputs_dir / name, sep="," if name.endswith("csv") else "\t")
        a = df[df["log2FC"] == -6.643856]
        assert len(a) == 5 and (a["value_source"] == "assigned").all()
        assert (df["value_source"] == "assigned").sum() == 5


def test_frzb_keeps_real_values_whatever_the_assigned_factor(inputs_dir, cfg):
    src = _read(inputs_dir / "regulated_up_partial.csv").set_index(ID).loc["P97401"]
    for factor in (2, 100, 12345):
        c = copy.deepcopy(cfg)
        c["assigned_fold_change"]["decrease_factor"] = factor
        row = wa.compute_walker(inputs_dir, c).frame.set_index(ID).loc["P97401"]
        assert row["value_source"] == "measured" and row["tier"] == "partial"
        for col in ("log2FC", "p_value", "adj_p_value"):
            assert row[col] == src[col]


# ---------------------------------------------------------------------------
# 7. Ordering
# ---------------------------------------------------------------------------
def test_core_order_follows_the_core_file_and_additions_follow_config_order(inputs_dir, cfg):
    p = inputs_dir / "ipa_input_full.csv"
    shuffled = _read(p).sample(frac=1.0, random_state=7).reset_index(drop=True)
    shuffled.to_csv(p, index=False)
    c = copy.deepcopy(cfg)
    c["assigned_proteins"] = list(reversed(c["assigned_proteins"]))
    res = wa.compute_walker(inputs_dir, c)
    n = len(shuffled)
    assert list(res.frame[ID].iloc[:n]) == list(shuffled[ID])
    assert list(res.frame[ID].iloc[n]) == list("P97401")
    assert list(res.frame[ID].iloc[n + 1:n + 6]) == [e["accession"] for e in c["assigned_proteins"]]
    assert list(res.manifest[ID].iloc[:6]) == ["P97401"] + [e["accession"] for e in c["assigned_proteins"]]
    assert len(res.manifest) == 6 + res.counts["recovered"]
    assert list(res.manifest[ID].iloc[6:]) == list(res.frame[ID].iloc[n + 6:])


def test_additions_do_not_depend_on_the_row_order_of_the_other_inputs(inputs_dir, cfg):
    base = wa.compute_walker(inputs_dir, cfg)
    for name in ("qc_limma.csv", "single_condition_proteins.csv", "qualitative_changes.csv",
                 "regulated_up_partial.csv", "regulated_down_partial.csv"):
        df = _read(inputs_dir / name)
        df.sample(frac=1.0, random_state=3).to_csv(inputs_dir / name, index=False)
    again = wa.compute_walker(inputs_dir, cfg)
    pd.testing.assert_frame_equal(base.frame, again.frame)
    pd.testing.assert_frame_equal(base.manifest, again.manifest)
    assert base.steps == again.steps


def test_two_builds_give_identical_bytes_and_row_order(inputs_dir):
    wa.build_walker_from_results(inputs_dir)
    first = {n: (inputs_dir / n).read_bytes() for n in OUTPUT_NAMES}
    wa.build_walker_from_results(inputs_dir)
    assert first == {n: (inputs_dir / n).read_bytes() for n in OUTPUT_NAMES}


# ---------------------------------------------------------------------------
# 8. No mutation of inputs (frames, config, files) and nothing else written
# ---------------------------------------------------------------------------
def test_pure_functions_do_not_mutate_their_dataframe_arguments(inputs_dir, cfg):
    partial = [(_read(inputs_dir / n), lbl) for n, lbl in wa.PARTIAL_SOURCES]
    limma = _read(inputs_dir / "qc_limma.csv")
    qual = _read(inputs_dir / "qualitative_changes.csv")
    single = _read(inputs_dir / "single_condition_proteins.csv")
    full = _read(inputs_dir / "ipa_input_full.csv")

    def snap():
        return ([f.copy(deep=True) for f, _ in partial], limma.copy(deep=True),
                qual.copy(deep=True), single.copy(deep=True), full.copy(deep=True),
                copy.deepcopy(cfg))

    before = snap()
    measured = wa.build_measured_rows(cfg, partial)
    assigned = wa.build_assigned_rows(cfg, single, qual, CTRL, TRT)
    pool = wa.build_recovered_pool(partial, limma, qual, READINGS)
    pool_before = pool.copy(deep=True)
    listed = set(measured[ID]) | set(assigned[ID]) | set(full[ID])
    passed, _ = wa.apply_recovered_rule(pool, cfg["recovered_rule"], listed, 4)
    passed_before = passed.copy(deep=True)
    rec = wa.build_recovered_rows(passed)
    adds = [measured.copy(deep=True), assigned.copy(deep=True), rec.copy(deep=True)]
    adds_before = [a.copy(deep=True) for a in adds]
    wa.build_walker_frame(full, adds)

    after = snap()
    for b, a in zip(before[0], after[0]):
        pd.testing.assert_frame_equal(b, a)
    for b, a in zip(before[1:5], after[1:5]):
        pd.testing.assert_frame_equal(b, a)
    assert before[5] == after[5]
    pd.testing.assert_frame_equal(pool, pool_before)
    pd.testing.assert_frame_equal(passed, passed_before)
    for b, a in zip(adds_before, adds):
        pd.testing.assert_frame_equal(b, a)


def test_a_recovered_survivor_frame_is_not_mutated_by_build_recovered_rows():
    pool = pd.DataFrame([{ID: "P1", "Gene names": "g", "log2FC": 1.0, "regulated": "UP",
                          "p_value": 0.1, "adj_p_value": 0.9, "n_imputed": 0, "n_seen": 4,
                          "source": "partial"}])
    before = pool.copy(deep=True)
    out = wa.build_recovered_rows(pool)
    pd.testing.assert_frame_equal(pool, before)
    assert out.iloc[0]["tier"] == "partial" and "tier" not in pool.columns
    assert out.iloc[0]["value_source"] == "measured"


def test_compute_does_not_mutate_cfg_or_input_files_and_writes_nothing(inputs_dir, cfg):
    cfg_before = copy.deepcopy(cfg)
    files_before = {p.name: _sha(p) for p in inputs_dir.iterdir()}
    wa.compute_walker(inputs_dir, cfg)
    assert cfg == cfg_before
    assert {p.name: _sha(p) for p in inputs_dir.iterdir()} == files_before   # no new files either


def test_compute_on_the_committed_results_writes_nothing(results_dir, cfg):
    names_before = sorted(os.listdir(results_dir))
    wa.compute_walker(results_dir, cfg)
    assert sorted(os.listdir(results_dir)) == names_before


class _WriteSpy:
    """Record every path opened for writing, via open(), to_csv and Path.write_*."""

    def __init__(self, monkeypatch):
        self.paths = []
        real_open = builtins.open

        def spy_open(file, mode="r", *a, **kw):
            if isinstance(file, (str, bytes, os.PathLike)) and any(m in mode for m in "wax+"):
                self.paths.append(Path(os.fsdecode(file)).resolve())
            return real_open(file, mode, *a, **kw)

        monkeypatch.setattr(builtins, "open", spy_open)
        real_to_csv = pd.DataFrame.to_csv

        def spy_to_csv(df, path_or_buf=None, *a, **kw):
            if isinstance(path_or_buf, (str, os.PathLike)):
                self.paths.append(Path(path_or_buf).resolve())
            return real_to_csv(df, path_or_buf, *a, **kw)

        monkeypatch.setattr(pd.DataFrame, "to_csv", spy_to_csv)
        for attr in ("write_text", "write_bytes"):
            real = getattr(Path, attr)

            def make(real):
                def spy(self_, *a, **kw):
                    self.paths.append(self_.resolve())
                    return real(self_, *a, **kw)
                return spy

            monkeypatch.setattr(Path, attr, make(real))


def test_a_build_writes_exactly_the_three_new_files_and_nothing_else(inputs_dir, monkeypatch):
    spy = _WriteSpy(monkeypatch)
    wa.build_walker_from_results(inputs_dir)
    here = inputs_dir.resolve()
    inside = {p for p in spy.paths if here in p.parents}
    assert {p.name for p in inside} == set(OUTPUT_NAMES)
    outside = {p for p in spy.paths if here not in p.parents}
    repo_results = (Path(wa.__file__).resolve().parents[1] / "results")
    assert not [p for p in outside if repo_results in p.parents], outside
    for name in REAL_PROTECTED:
        assert not [p for p in spy.paths if p.name == name], f"{name} was opened for writing"


def test_originals_in_the_scratch_dir_are_byte_identical_after_a_build(inputs_dir):
    protected = [n for n in set(wa.BYTE_FROZEN_INPUTS) | set(wa._INPUTS) | set(REAL_PROTECTED)
                 if (inputs_dir / n).exists()]
    assert {"ipa_input.csv", "ipa_input_extended.csv"} <= set(protected)
    before = {n: _sha(inputs_dir / n) for n in protected}
    wa.build_walker_from_results(inputs_dir)
    assert {n: _sha(inputs_dir / n) for n in protected} == before


def test_the_freeze_guard_detects_a_modified_or_deleted_original(inputs_dir):
    before = wa._hash_frozen_inputs(inputs_dir)
    wa.assert_frozen_inputs_untouched(inputs_dir, before)       # unchanged -> fine
    with open(inputs_dir / "ipa_input.csv", "ab") as fh:
        fh.write(b"\n")
    with pytest.raises(WalkerAdditionsError, match="ipa_input.csv"):
        wa.assert_frozen_inputs_untouched(inputs_dir, before)
    (inputs_dir / "ipa_input_extended.csv").unlink()
    with pytest.raises(WalkerAdditionsError, match="ipa_input_extended.csv"):
        wa.assert_frozen_inputs_untouched(inputs_dir, before)


def test_the_original_ipa_inputs_are_declared_byte_frozen():
    assert {"ipa_input.csv", "ipa_input_extended.csv", "ipa_input_full.csv",
            "ipa_input_full.txt"} <= set(wa.BYTE_FROZEN_INPUTS)
    assert not set(OUTPUT_NAMES) & set(wa.BYTE_FROZEN_INPUTS)


def test_the_cli_writes_only_into_the_requested_dir_and_leaves_the_committed_files_alone(
        inputs_dir, results_dir):
    sidecars_before = {p.name: p.read_bytes() for p in results_dir.glob("*walker*.provenance.json")}
    assert wa.main(["--results-dir", str(inputs_dir)]) == 0
    for name in OUTPUT_NAMES:
        assert (inputs_dir / name).exists() and (inputs_dir / f"{name}.provenance.json").exists()
    assert {p.name: p.read_bytes() for p in results_dir.glob("*walker*.provenance.json")} == sidecars_before


def test_the_cli_no_sidecars_flag_writes_no_sidecar(inputs_dir):
    assert wa.main(["--results-dir", str(inputs_dir), "--no-sidecars"]) == 0
    assert not list(inputs_dir.glob("*.provenance.json"))


# ---------------------------------------------------------------------------
# 9. Writers on awkward text
# ---------------------------------------------------------------------------
def test_gene_names_with_commas_quotes_and_unicode_round_trip_in_both_twins(tmp_path):
    awkward = ['Gene,with,commas', 'Gene "quoted"', "Genéµ", "A;B;C"]
    core = pd.DataFrame({
        ID: [f"C{i}" for i in range(4)], "Gene names": awkward,
        "log2FC": [1.0, -1.0, 2.0, -2.0], "regulated": ["UP", "DOWN", "UP", "DOWN"],
        "p_value": [0.1] * 4, "adj_p_value": [0.9] * 4})
    frame = wa.build_walker_frame(core, [])
    csv_path, txt_path = tmp_path / "w.csv", tmp_path / "w.txt"
    wa.write_walker(frame, csv_path, txt_path)
    for path, delim in ((csv_path, ","), (txt_path, "\t")):
        with open(path, newline="", encoding="utf-8") as fh:
            got = [r["Gene names"] for r in csv.DictReader(fh, delimiter=delim)]
        assert got == awkward


def test_write_walker_can_skip_the_txt_twin(tmp_path):
    frame = wa.build_walker_frame(pd.DataFrame({
        ID: ["C1"], "Gene names": ["a"], "log2FC": [1.0], "regulated": ["UP"],
        "p_value": [0.1], "adj_p_value": [0.9]}), [])
    written = wa.write_walker(frame, tmp_path / "w.csv", txt_path=False)
    assert [p.name for p in written] == ["w.csv"] and not (tmp_path / "w.txt").exists()
