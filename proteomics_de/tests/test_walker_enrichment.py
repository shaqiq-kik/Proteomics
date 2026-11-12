"""Contract for ``enrich/walker_enrichment.py`` (DECISIONS_LOG D21): ORA and STRING
over Dr. Walker's 721-protein list (D24: the first 968-protein version was superseded).

Offline by default. The stage REPLAYS the committed ``raw/`` API cache; every test
that runs it first makes ``requests.get/post`` raise, so a network call is a test
failure, not a silent fallback. Nothing here talks to g:Profiler or STRING. The one
live-network test is marked ``network`` (excluded by the default ``addopts``).

What these tests exist to catch:

1. **The 715 results being disturbed.** The walker analysis is additive. The old
   ``results/enrichment/`` files must be byte-identical before and after a run, and
   the new code must reproduce the old outputs when fed the old inputs (the
   equivalence test): if it does, differences between the 715 and the 721 runs are
   the list, not the code.
2. **A cache that answers a different question.** Replay must refuse a cached
   response whose request fingerprint differs, and must never go live by itself.
3. **Invented numbers.** The five assigned log2FC values (-6.643856) are placeholders.
   ORA and STRING use membership only: the query lists must not depend on them, the
   node table must leave log2FC blank on those rows, and no p-value reaches any output.
4. **Confounded comparisons.** The committed 715 run used older STRING / g:Profiler
   releases; the comparison must carry the same-version control and say which
   versions were used.
"""

from __future__ import annotations

import filecmp
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PDE = _REPO_ROOT / "proteomics_de"
for _p in (_PDE / "enrich", _PDE, _REPO_ROOT):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import walker_enrichment as we  # noqa: E402
import ora  # noqa: E402

from tools import freeze  # noqa: E402
import run_pipeline as rp  # noqa: E402

RESULTS = _PDE / "results"
OLD = RESULTS / "enrichment"
NEW = RESULTS / "enrichment_walker"
CONTROL = NEW / "control_715_same_version"
MANIFEST = _PDE / "tests" / "expected" / "outputs.sha256"

RUN_FILES = we.RUN_FILES
ASSIGNED = {"P01132", "Q61521", "O08542", "P55144", "Q9D1H9"}
NAMED = {"Egf": "P01132", "Ereg": "Q61521", "Efna4": "O08542", "Tyro3": "P55144",
         "Mfap4": "Q9D1H9", "Frzb": "P97401"}
FLAGS = {"kept", "gained", "lost", ""}


def _no_network(monkeypatch):
    import requests

    def boom(*a, **k):  # pragma: no cover - only fires on a regression
        raise AssertionError("network call attempted in an offline test")

    monkeypatch.setattr(requests, "post", boom)
    monkeypatch.setattr(requests, "get", boom)


@pytest.fixture(scope="module")
def walker_df():
    return we.load_walker_list()


@pytest.fixture(scope="module")
def replay(tmp_path_factory):
    """Replay the committed cache into a scratch dir; the real outputs stay put."""
    from _pytest.monkeypatch import MonkeyPatch

    mp = MonkeyPatch()
    _no_network(mp)
    out = tmp_path_factory.mktemp("walker_replay") / "enrichment_walker"
    for sub in (".", we.load_config()["control_dir_name"]):
        (out / sub / "raw").mkdir(parents=True, exist_ok=True)
        for f in (NEW / sub / "raw").iterdir():
            shutil.copy(f, out / sub / "raw" / f.name)
    before = we.tree_hashes(OLD)
    result = we.run(fetch=False, out_dir=out, sidecars=False)
    after = we.tree_hashes(OLD)
    yield {"out": out, "result": result, "old_before": before, "old_after": after}
    mp.undo()


# --------------------------------------------------------------------------
# The list that was analysed
# --------------------------------------------------------------------------
def test_config_section_is_load_bearing_and_names_the_six_proteins():
    cfg = we.load_config()
    assert cfg["output_dir"].endswith("results/enrichment_walker")
    assert {f["gene"]: f["accession"] for f in cfg["focus_proteins"]} == NAMED
    assert cfg["string"]["api_base"].startswith("https://string-db.org/api")
    assert cfg["comparison"]["top_n"] == 20


def test_walker_list_is_the_721_and_its_up_down_split_is_frozen(walker_df, frozen_counts):
    assert len(walker_df) == frozen_counts["n_ipa_walker"] == 721
    up, down, meta = we.build_queries(walker_df)
    assert (len(up), len(down)) == (frozen_counts["n_walker_ora_up"],
                                    frozen_counts["n_walker_ora_down"]) == (510, 211)
    assert len(up) + len(down) == 721
    # UP = core UP + FRZB ; DOWN = core DOWN + the 5 assigned (no recovered protein since D24)
    assert frozen_counts["n_walker_recovered"] == 0
    assert len(up) == frozen_counts["n_up"] + frozen_counts["n_walker_measured_added"]
    assert len(down) == frozen_counts["n_down"] + frozen_counts["n_walker_assigned"]
    assert meta["up_n"] == len(up) and meta["down_n"] == len(down)
    assert set(up).isdisjoint(down)


def test_the_five_assigned_proteins_are_in_the_down_query(walker_df):
    assigned = walker_df[walker_df["value_source"] == "assigned"]
    assert set(assigned[we.ACC_COL]) == ASSIGNED
    assert (assigned["regulated"] == "DOWN").all()
    _up, down, _ = we.build_queries(walker_df)
    assert {"Egf", "Ereg", "Efna4", "Tyro3", "Mfap4"} <= set(down)


def test_committed_gprofiler_responses_carry_the_721_queries(walker_df):
    """The cached responses must have been asked about THIS list, not the 715 or the old 968."""
    up, down, _ = we.build_queries(walker_df)
    for d, q in (("up", up), ("down", down)):
        resp = json.loads((NEW / "raw" / f"gprofiler_{d}.json").read_text())
        assert resp["meta"]["query_metadata"]["queries"]["query_1"] == q
        assert resp["meta"]["query_metadata"]["domain_scope"] == "custom"
        assert resp["meta"]["query_metadata"]["significance_threshold_method"] == "g_SCS"


def test_the_cached_string_response_answers_exactly_the_721_accessions(walker_df):
    """No response of the superseded 968 run may be left in the cache (it would be stale)."""
    ids = json.loads((NEW / "raw" / "string_get_ids.json").read_text())
    assert {r["queryItem"] for r in ids} <= set(walker_df[we.ACC_COL])
    manifest = json.loads((NEW / "raw" / we.MANIFEST_NAME).read_text())
    assert set(manifest) == {
        "gprofiler_up.json", "gprofiler_up_all.json.gz", "gprofiler_up_all_evidence.json.gz",
        "gprofiler_down.json", "gprofiler_down_all.json.gz", "gprofiler_down_all_evidence.json.gz",
        "string_get_ids.json", "string_network.tsv", "string_version.json"}
    up = json.loads((NEW / "raw" / "gprofiler_up.json").read_text())
    down = json.loads((NEW / "raw" / "gprofiler_down.json").read_text())
    assert len(up["meta"]["query_metadata"]["queries"]["query_1"]) == 510
    assert len(down["meta"]["query_metadata"]["queries"]["query_1"]) == 211


# --------------------------------------------------------------------------
# Background
# --------------------------------------------------------------------------
def test_every_query_symbol_is_inside_the_existing_background(walker_df):
    background, bg_meta = we.load_background()
    up, down, _ = we.build_queries(walker_df)
    assert set(up) | set(down) <= set(background)
    cover = we.assert_queries_in_background(up, down, background)
    assert cover == {"n_query_not_in_background": 0, "background_extended": False}
    assert bg_meta["background_row_union_n"] == 2552


def test_a_query_outside_the_background_is_refused_not_absorbed(walker_df):
    background, _ = we.load_background()
    up, down, _ = we.build_queries(walker_df)
    with pytest.raises(we.WalkerEnrichmentError, match="outside the custom background"):
        we.assert_queries_in_background(up + ["NotAGeneInTheUniverse"], down, background)


def test_all_721_accessions_are_in_the_ipa_reference_set(walker_df):
    ref = (RESULTS / "ipa_background_measured.txt").read_text().splitlines()[1:]
    assert set(walker_df[we.ACC_COL]) <= set(ref)
    assert len(ref) == 2552


# --------------------------------------------------------------------------
# Outputs on disk (committed run)
# --------------------------------------------------------------------------
def test_ora_outputs_record_the_list_the_versions_and_the_unmapped_symbols(walker_df, frozen_counts):
    meta = json.loads((NEW / "ora_meta.json").read_text())
    assert meta["list"] == "walker721"
    assert meta["gprofiler_version"]
    assert meta["domain_scope"] == "custom"
    assert meta["background_extended"] is False and meta["n_query_not_in_background"] == 0
    assert meta["n_query_submitted_up"] == frozen_counts["n_walker_ora_up"]
    assert meta["n_query_submitted_down"] == frozen_counts["n_walker_ora_down"]
    assert meta["n_failed_to_map_up"] == len(meta["failed_to_map_symbols_up"])
    for d in ("up", "down"):
        assert (meta[f"n_query_mapped_{d}"] + meta[f"n_failed_to_map_{d}"]
                + meta[f"n_ambiguous_{d}"]) == meta[f"n_query_submitted_{d}"]
    assert meta["n_terms_up"] == len(pd.read_csv(NEW / "ora_up.csv"))
    assert meta["n_terms_down"] == len(pd.read_csv(NEW / "ora_down.csv"))
    assert list(pd.read_csv(NEW / "ora_up.csv").columns) == ora.CSV_COLUMNS


def test_string_outputs_match_the_frozen_counts(walker_df, frozen_counts):
    meta = json.loads((NEW / "string_meta.json").read_text())
    nodes = pd.read_csv(NEW / "string_node_metrics.csv")
    edges = pd.read_csv(NEW / "string_edges.tsv", sep="\t")
    assert meta["n_seeds"] == 721 and meta["species"] == 10090
    assert meta["string_version"] and meta["string_stable_address"]
    assert meta["n_mapped"] == frozen_counts["n_walker_string_mapped"] == len(nodes)
    assert meta["n_nodes_in_network"] == frozen_counts["n_walker_string_nodes"] == len(nodes)
    assert meta["n_edges"] == frozen_counts["n_walker_string_edges"] == len(edges)
    assert meta["n_mapped"] + meta["n_unmapped"] == 721
    # nodes are a subset of the list; direction and tier agree with the list
    lst = walker_df.set_index(we.ACC_COL)
    assert set(nodes["accession"]) <= set(lst.index)
    assert (nodes["regulated"].values == lst.loc[nodes["accession"], "regulated"].values).all()
    assert (nodes["tier"].values == lst.loc[nodes["accession"], "tier"].values).all()
    assert set(nodes["accession"]).isdisjoint(meta["unmapped_accessions"])


def test_assigned_values_and_p_values_never_reach_the_outputs(walker_df):
    nodes = pd.read_csv(NEW / "string_node_metrics.csv")
    assigned = nodes[nodes["value_source"] == "assigned"]
    assert set(assigned["accession"]) == ASSIGNED
    assert assigned["log2FC"].isna().all(), "assigned placeholder leaked into the node table"
    measured = nodes[nodes["value_source"] == "measured"]
    assert measured["log2FC"].notna().all()
    for f in RUN_FILES:
        for root in (NEW, CONTROL):
            p = root / f
            if p.suffix in (".csv", ".tsv"):
                cols = pd.read_csv(p, sep="\t" if p.suffix == ".tsv" else ",").columns
                assert not {"p_value_input", "adj_p_value"} & set(cols), p
    assert "adj_p_value" not in nodes.columns
    assert "-6.643856" not in (NEW / "string_node_metrics.csv").read_text()
    assert "-6.64" not in (NEW / "string_meta.json").read_text()


def test_queries_and_graph_ignore_log2fc_and_p_values(walker_df):
    """Membership only: scrambling every number must not move a single query symbol."""
    scrambled = walker_df.copy()
    scrambled["log2FC"] = np.linspace(-9, 9, len(scrambled))
    assert we.build_queries(scrambled)[:2] == we.build_queries(walker_df)[:2]
    masked = we.mask_assigned(scrambled)
    assert masked.loc[masked["value_source"] == "assigned", "log2FC"].isna().all()
    assert masked.loc[masked["value_source"] == "measured", "log2FC"].notna().all()
    assert "p_value" not in walker_df.columns and "adj_p_value" not in walker_df.columns


def test_the_six_named_proteins_are_reported(frozen_counts):
    focus = pd.read_csv(NEW / "string_focus_proteins.csv")
    nodes = pd.read_csv(NEW / "string_node_metrics.csv").set_index("accession")
    assert dict(zip(focus["gene"], focus["accession"])) == NAMED
    for _, r in focus.iterrows():
        assert bool(r["in_network"]) and r["accession"] in nodes.index
        assert r["degree"] == nodes.loc[r["accession"], "degree"]
        assert r["degree"] == len(str(r["neighbours"]).split(";"))
    assert set(focus.loc[focus["gene"] != "Frzb", "value_source"]) == {"assigned"}
    assert focus.loc[focus["gene"] == "Frzb", "value_source"].iloc[0] == "measured"


# --------------------------------------------------------------------------
# The 715 run is untouched, and the new code reproduces it
# --------------------------------------------------------------------------
def test_old_enrichment_files_match_the_freeze_manifest():
    manifest = freeze.read_manifest(MANIFEST)
    old_entries = {k: v for k, v in manifest.items()
                   if k.startswith("proteomics_de/results/enrichment/")}
    assert len(old_entries) == 19
    for rel, (sha, _mode) in old_entries.items():
        assert freeze.digest(_REPO_ROOT / rel)[0] == sha, f"{rel} drifted"
    # and none of the walker files live there
    assert not any("walker" in k for k in old_entries)


def test_new_code_reproduces_the_committed_715_outputs_from_the_committed_cache(tmp_path, monkeypatch):
    """Equivalence: old inputs through the new code give byte-identical old outputs."""
    _no_network(monkeypatch)
    store = we.RawStore(OLD / "raw", require_manifest=False, writable=False)
    background, bg_meta = we.load_background()
    core = we.load_core_list()
    up, down, lm = we.build_queries(core)
    we.run_ora(up, down, background, store, tmp_path, label="eq", list_meta=lm,
               bg_meta=bg_meta, description="equivalence")
    we.run_string(core, store, tmp_path, api_base="https://string-db.org/api",
                  caller_identity="x", label="eq", description="eq", caveat="eq")
    for f in ("ora_up.csv", "ora_down.csv", "string_edges.tsv", "string_node_metrics.csv"):
        assert filecmp.cmp(OLD / f, tmp_path / f, shallow=False), f"{f} differs from the 715 run"
    old_detail = json.loads((OLD / "ora_top_terms_detail.json").read_text())["terms"]
    new_detail = json.loads((tmp_path / "ora_top_terms_detail.json").read_text())["terms"]
    assert old_detail == new_detail
    assert store.network_calls == 0


def test_gprofiler_payload_is_exactly_what_ora_sends(monkeypatch):
    import ora as ora_mod

    seen = {}

    class _R:
        def raise_for_status(self):
            pass

        def json(self):
            return {}

    def fake_post(url, json=None, timeout=None):
        seen["payload"], seen["url"] = json, url
        return _R()

    monkeypatch.setattr(ora_mod.requests, "post", fake_post)
    ora_mod._post(["A", "B"], ["A", "B", "C"], True, False)
    assert seen["payload"] == we.gprofiler_payload(["A", "B"], ["A", "B", "C"], True, False)
    assert seen["url"] == ora_mod.GPROFILER_URL


# --------------------------------------------------------------------------
# Replay: offline, deterministic, additive, and honest about a bad cache
# --------------------------------------------------------------------------
def test_replay_makes_no_network_call_and_is_byte_identical_to_the_committed_run(replay):
    assert replay["result"]["network_calls"] == 0
    out = replay["out"]
    n = 0
    for root_rel in (".", we.load_config()["control_dir_name"]):
        for f in RUN_FILES:
            assert filecmp.cmp(NEW / root_rel / f, out / root_rel / f, shallow=False), \
                f"{root_rel}/{f} is not reproduced by replay"
            n += 1
    for f in (we.COMPARISON_NAME, we.FOCUS_NAME, we.STRING_COMPARISON_NAME):
        assert filecmp.cmp(NEW / f, out / f, shallow=False), f
        n += 1
    assert n == 2 * len(RUN_FILES) + 3


def test_replay_leaves_the_715_outputs_byte_identical(replay):
    assert replay["old_before"] == replay["old_after"]
    assert replay["out"].resolve() != OLD.resolve() and OLD.resolve() not in replay["out"].resolve().parents


def test_replay_never_rewrites_the_raw_cache(replay):
    for sub in (".", we.load_config()["control_dir_name"]):
        for f in (NEW / sub / "raw").iterdir():
            assert filecmp.cmp(f, replay["out"] / sub / "raw" / f.name, shallow=False), f.name


def test_a_changed_list_is_refused_by_the_cache_not_silently_replayed(tmp_path, monkeypatch, walker_df):
    _no_network(monkeypatch)
    shutil.copytree(NEW / "raw", tmp_path / "raw")
    store = we.RawStore(tmp_path / "raw")
    background, _ = we.load_background()
    up, down, _ = we.build_queries(walker_df)
    with pytest.raises(we.CacheMismatch, match="different request"):
        we._gprofiler_get(store, "gprofiler_up.json", up[:-1], background, False, False, "t")
    # the unchanged request is accepted
    assert we._gprofiler_get(store, "gprofiler_up.json", up, background, False, False, "t")


def test_a_missing_cache_is_an_error_in_replay_mode_never_a_live_call(tmp_path, monkeypatch):
    _no_network(monkeypatch)
    store = we.RawStore(tmp_path / "raw", fetch=False)
    with pytest.raises(we.CacheMiss, match="--fetch"):
        store.get_text("x.json", "sha", "ep", lambda: "never", label="t")
    assert store.network_calls == 0


def test_a_cached_file_without_a_fingerprint_is_refused(tmp_path):
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "x.json").write_text("{}")
    with pytest.raises(we.CacheMismatch, match="no request fingerprint"):
        we.RawStore(raw).get_text("x.json", "sha", "ep", lambda: "", label="t")


def test_every_cached_response_has_a_fingerprint():
    for root in (NEW / "raw", CONTROL / "raw"):
        manifest = json.loads((root / we.MANIFEST_NAME).read_text())
        files = {p.name for p in root.iterdir() if p.name != we.MANIFEST_NAME}
        assert files == set(manifest), root
        assert all(len(v["request_sha256"]) == 64 for v in manifest.values())


def test_the_new_code_never_writes_into_the_715_directory(tmp_path):
    with pytest.raises(we.WalkerEnrichmentError, match="protected"):
        we.assert_not_inside(OLD / "walker", OLD)
    with pytest.raises(we.WalkerEnrichmentError):
        we.assert_not_inside(OLD, OLD)
    we.assert_not_inside(NEW, OLD)


# --------------------------------------------------------------------------
# comparison_vs_715.csv
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def comparison():
    return pd.read_csv(NEW / we.COMPARISON_NAME)


def test_comparison_schema(comparison, frozen_counts):
    cols = ["direction", "source", "term_id", "term_name"]
    for run in ("committed715", "control715", "walker721"):
        cols += [f"{run}_adj_p", f"{run}_hypergeom_p", f"{run}_rank",
                 f"{run}_intersection_size"]
    cols += ["sig_flag_vs_committed", "top20_flag_vs_committed",
             "sig_flag_vs_control", "top20_flag_vs_control"]
    assert list(comparison.columns) == cols
    assert len(comparison) == frozen_counts["n_walker_comparison_rows"]
    assert set(comparison["direction"]) == {"UP", "DOWN"}
    assert not comparison.duplicated(["direction", "term_id"]).any()
    for c in cols[-4:]:
        assert set(comparison[c].fillna("")) <= FLAGS
    assert ((comparison["walker721_adj_p"].dropna() >= 0) & (comparison["walker721_adj_p"].dropna() <= 1)).all()


def test_comparison_flags_follow_the_ranks(comparison):
    top = 20
    for tag, run in (("committed", "committed715"), ("control", "control715")):
        for _, r in comparison.iterrows():
            a = pd.notna(r[f"{run}_rank"]) and r[f"{run}_rank"] <= top
            b = pd.notna(r["walker721_rank"]) and r["walker721_rank"] <= top
            expected = we._flag(bool(a), bool(b))
            got = r[f"top{top}_flag_vs_{tag}"]
            assert (got if isinstance(got, str) else "") == expected, (r["term_id"], tag)


def test_comparison_lists_a_term_only_if_it_is_a_lead_in_some_run(comparison):
    for _, r in comparison.iterrows():
        ranks = [r[f"{x}_rank"] for x in ("committed715", "control715", "walker721")
                 if pd.notna(r[f"{x}_rank"])]
        assert (ranks and min(ranks) <= 20) or r["walker721_adj_p"] < 0.05


def test_comparison_ranks_agree_with_the_raw_responses(comparison):
    t = we.all_terms_table(OLD / "raw", "down").set_index("term_id")
    row = comparison[(comparison["direction"] == "DOWN")
                     & (comparison["term_id"] == t.index[0])].iloc[0]
    assert row["committed715_rank"] == 1
    assert row["committed715_adj_p"] == pytest.approx(t.iloc[0]["adj_p"])
    assert t.iloc[0]["adj_p"] == pytest.approx(0.7029801000701202)
    # the uncorrected p is a ranking aid that never exceeds the corrected p
    allt = we.all_terms_table(NEW / "raw", "down")
    assert (allt["raw_p"] <= allt["adj_p"] + 1e-12).all()
    assert allt["rank"].tolist() == list(range(1, len(allt) + 1))


def test_seven_of_the_top_twenty_down_leads_contain_an_assigned_protein(comparison):
    """Why the DOWN leads are fragile: the README's headline caveat, recomputed from the cache.

    (The first 968-protein run had 12 of 20; the 721 run has 7, and they are exactly the
    terms the comparison flags as gained against the same-version control.)
    """
    import gzip

    assigned = {"Egf", "Ereg", "Efna4", "Tyro3", "Mfap4"}
    resp = json.loads(gzip.open(NEW / "raw" / "gprofiler_down_all_evidence.json.gz", "rt").read())
    order = ora._query_order(resp)
    by_id = {t["native"]: t for t in resp["result"]}
    top = we.all_terms_table(NEW / "raw", "down").head(20)
    n, with_assigned = 0, set()
    for tid in top["term_id"]:
        genes = {g for g, e in zip(order, by_id[tid]["intersections"]) if e}
        n += bool(genes & assigned)
        if genes & assigned:
            with_assigned.add(tid)
    assert n == 7
    assert assigned <= set(order), "the five assigned proteins are in the mapped DOWN query"
    down = comparison[comparison["direction"] == "DOWN"]
    gained = set(down.loc[down["top20_flag_vs_control"] == "gained", "term_id"])
    assert gained == with_assigned, "every gained DOWN lead contains an assigned protein, and vice versa"


def test_leave_the_five_out_recomputation_reproduces_the_same_version_control():
    """Validates the README's approximate leave-five-out ranks.

    Removing the five assigned proteins from the 721 DOWN list leaves exactly the 715 control's
    DOWN list, so the approximate recomputation (uncorrected hypergeometric p with the intersection
    reduced by the assigned hits, query size by 5, background unchanged) must reproduce the
    control run's own ranks of the seven assigned-containing leads. It does, to the rank.
    """
    import gzip

    from scipy.stats import hypergeom

    assigned = {"Egf", "Ereg", "Efna4", "Tyro3", "Mfap4"}
    resp = json.loads(gzip.open(NEW / "raw" / "gprofiler_down_all_evidence.json.gz", "rt").read())
    order = ora._query_order(resp)
    n_assigned = len(assigned & set(order))
    recs = []
    for t in resp["result"]:
        genes = {g for g, e in zip(order, t["intersections"]) if e}
        k, K, M = t["intersection_size"] - len(genes & assigned), t["term_size"], t["effective_domain_size"]
        n = t["query_size"] - n_assigned
        recs.append((float(hypergeom.sf(k - 1, M, K, n)) if k > 0 else 1.0, t["native"], genes & assigned))
    recs.sort(key=lambda r: (r[0], r[1]))
    approx_rank = {tid: i + 1 for i, (_p, tid, _g) in enumerate(recs)}
    control = we.all_terms_table(CONTROL / "raw", "down").set_index("term_id")
    top = we.all_terms_table(NEW / "raw", "down").head(20)
    checked = 0
    for tid in top["term_id"]:
        hit = next(g for _p, t, g in recs if t == tid)
        if hit:
            assert approx_rank[tid] == control.loc[tid, "rank"], tid
            checked += 1
    assert checked == 7


def test_no_term_is_called_significant_when_none_passes(comparison):
    meta = json.loads((NEW / "ora_meta.json").read_text())
    if meta["n_terms_up"] == 0 and meta["n_terms_down"] == 0:
        assert comparison["sig_flag_vs_committed"].isna().all()
        assert comparison["sig_flag_vs_control"].isna().all()


# --------------------------------------------------------------------------
# Version drift and the STRING comparison
# --------------------------------------------------------------------------
def test_versions_are_recorded_and_the_drift_against_the_committed_run_is_stated():
    sc = json.loads((NEW / we.STRING_COMPARISON_NAME).read_text())
    committed_ora = json.loads((OLD / "ora_meta.json").read_text())
    committed_str = json.loads((OLD / "string_meta.json").read_text())
    walker_ora = json.loads((NEW / "ora_meta.json").read_text())
    walker_str = json.loads((NEW / "string_meta.json").read_text())
    control_ora = json.loads((CONTROL / "ora_meta.json").read_text())
    control_str = json.loads((CONTROL / "string_meta.json").read_text())
    assert sc["gprofiler_version_committed_715"] == committed_ora["gprofiler_version"]
    assert sc["gprofiler_version_walker_721"] == walker_ora["gprofiler_version"]
    assert sc["same_gprofiler_version_as_committed"] == (
        committed_ora["gprofiler_version"] == walker_ora["gprofiler_version"])
    assert sc["same_string_version_as_committed"] == (
        committed_str["string_version"] == walker_str["string_version"])
    # the control is genuinely same-version as the walker run; that is its whole point
    assert control_ora["gprofiler_version"] == walker_ora["gprofiler_version"]
    assert control_str["string_version"] == walker_str["string_version"]
    assert sc["same_string_version_as_control"] is True


def test_string_comparison_counts_agree_with_the_meta_files():
    sc = json.loads((NEW / we.STRING_COMPARISON_NAME).read_text())
    for key, path in (("committed_715", OLD), ("control_715_same_version", CONTROL),
                      ("walker_721", NEW)):
        m = json.loads((path / "string_meta.json").read_text())
        assert sc[key]["n_edges"] == m["n_edges"]
        assert sc[key]["n_nodes_in_network"] == m["n_nodes_in_network"]
    assert sc["core_core_edges_equal_control"] is True
    assert (sc["walker_edges_between_two_core_proteins"]
            + sc["walker_edges_touching_an_added_protein"]) == sc["walker_721"]["n_edges"]


def test_the_control_is_the_715_list_and_carries_no_assigned_rows():
    nodes = pd.read_csv(CONTROL / "string_node_metrics.csv")
    meta = json.loads((CONTROL / "string_meta.json").read_text())
    assert meta["n_seeds"] == 715
    assert not set(nodes["accession"]) & ASSIGNED
    assert "tier" not in nodes.columns  # same columns as the committed 715 table


# --------------------------------------------------------------------------
# Pipeline, freeze and docs
# --------------------------------------------------------------------------
def test_stage_is_registered_offline_and_after_its_inputs():
    stage = rp.STAGE_BY_ID["enrich_walker"]
    assert stage.network is False, "the stage replays the cache; live calls are --fetch only"
    assert stage.script == "proteomics_de/enrich/walker_enrichment.py"
    assert set(stage.depends_on) == {"ipa_export_walker", "enrich_string_ppi", "enrich_ora"}
    for out in stage.outputs:
        assert out.path.startswith("proteomics_de/results/enrichment_walker/")
        assert not (rp.ROOT / out.path).parent.samefile(OLD)
        assert (rp.ROOT / out.path).is_file()


def test_walker_enrichment_outputs_are_byte_frozen():
    manifest = freeze.read_manifest(MANIFEST)
    frozen = {p.relative_to(_REPO_ROOT).as_posix() for p in freeze.frozen_files(_REPO_ROOT)}
    walker = {k for k in frozen if k.startswith("proteomics_de/results/enrichment_walker/")}
    assert len(walker) >= 30
    assert walker <= set(manifest)
    for rel in walker:
        assert freeze.digest(_REPO_ROOT / rel)[0] == manifest[rel][0], rel
    assert not any(k.endswith(".provenance.json") for k in walker)


def test_ipa_analysis_walker_is_outside_the_freeze_by_policy():
    d = _PDE / "results" / "ipa_analysis_walker"
    assert (d / "UPLOAD_STEPS.md").is_file()
    excluded = dict(freeze.excluded_files(_REPO_ROOT))
    rel = "proteomics_de/results/ipa_analysis_walker/UPLOAD_STEPS.md"
    assert "external tool output" in excluded[rel]
    assert rel not in freeze.read_manifest(MANIFEST)


def test_readme_states_the_list_and_the_caveats():
    text = (NEW / "README.md").read_text(encoding="utf-8")
    low = text.lower()
    assert "721" in text and "ipa_input_walker" in text
    assert "510 UP / 211 DOWN" in text
    # D24: the 968-protein numbers must not be carried over into the live README.
    for stale in ("661", "937 ", "10,696", "618 UP", "12 of the 20", "66 terms"):
        assert stale not in text, stale
    assert "git history" in low and "D24" in text
    assert "assigned" in low and "not measurements" in low or "not a measurement" in low
    assert "membership" in low
    # D23: the walker README states the corrected design and the imputed-value
    # caveat, and must not claim a single/one biological replicate any more.
    assert "2 biological replicates per group" in low and "2 separate experiments" in low
    assert "single biological replicate" not in low and "one biological replicate" not in low
    assert "technical replicates only" not in low
    assert "imputed" in low and "fdr" in low
    meta = json.loads((NEW / "ora_meta.json").read_text())
    smeta = json.loads((NEW / "string_meta.json").read_text())
    assert meta["gprofiler_version"] in text
    assert f"STRING {smeta['string_version']}" in text


# --------------------------------------------------------------------------
# IPA preparation
# --------------------------------------------------------------------------
def test_walker_upload_has_the_extended_upload_format(walker_df):
    w = pd.read_csv(RESULTS / "ipa_input_walker.txt", sep="\t", dtype=str, keep_default_na=False)
    e = pd.read_csv(RESULTS / "ipa_input_extended.txt", sep="\t", dtype=str, keep_default_na=False)
    assert list(w.columns) == list(e.columns) + ["value_source"]
    raw = (RESULTS / "ipa_input_walker.txt").read_bytes()
    assert not raw.startswith(b"\xef\xbb\xbf") and b"\r" not in raw
    assert w.columns[0] == "UniProt Accession Number"
    # D24: walker minus the 5 assigned (716) is a subset of ipa_input_extended (963); the 247
    # proteins the extended file has beyond it are the recovered partial-tier proteins Dr. Walker
    # asked to drop. The sensitivity dataset is its own file, ipa_input_walker_without_assigned.txt.
    kept = w[w["value_source"] != "assigned"].drop(columns="value_source")
    assert len(kept) == 716
    assert set(kept["UniProt Accession Number"]) <= set(e["UniProt Accession Number"])
    assert len(set(e["UniProt Accession Number"]) - set(kept["UniProt Accession Number"])) == 247
    na = pd.read_csv(RESULTS / "ipa_input_walker_without_assigned.txt", sep="\t", dtype=str,
                     keep_default_na=False)
    assert list(na.columns) == list(w.columns) and len(na) == 716
    assert na.drop(columns="value_source").equals(kept.reset_index(drop=True))
    merged = kept.merge(e, on="UniProt Accession Number", suffixes=("_w", "_e"))
    for c in ("log2FC", "regulated", "p_value", "adj_p_value", "tier"):
        a, b = merged[f"{c}_w"], merged[f"{c}_e"]
        if c == "tier":
            continue  # recovered rows keep their source tier name; values are what matter
        assert (a == b).all(), c


def test_upload_steps_cover_what_the_user_must_do_and_warns_about_zscores():
    text = (_PDE / "results" / "ipa_analysis_walker" / "UPLOAD_STEPS.md").read_text()
    low = text.lower()
    for needle in ("ipa_input_walker.txt", "ipa_background_measured.txt", "reference set",
                   "canonical pathways", "upstream regulators", "diseases and functions",
                   "networks", "z-score", "-6.64", "sensitivity",
                   "ipa_input_walker_without_assigned.txt", "721", "716",
                   "log2fc", "ignore", "export"):
        assert needle in low, needle
    # D24: the old 968 / 963 row counts and the extended file must not be instructions any more.
    for stale in ("968", "963 rows"):
        assert stale not in text.replace("968-protein", "").replace("D24", ""), stale
    assert "assigned" in low and "not a measurement" in low or "not measurements" in low


# --------------------------------------------------------------------------
# Live (opt-in)
# --------------------------------------------------------------------------
@pytest.mark.network
def test_live_version_probes_answer():  # pragma: no cover - needs the network
    import requests

    assert requests.get("https://string-db.org/api/json/version", timeout=30).json()[0]["string_version"]
