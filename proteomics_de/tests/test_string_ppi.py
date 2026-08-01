"""Contract for ``enrich/string_ppi.py`` and its committed STRING outputs.

Offline by default. Everything here reads either the source module or the files
already on disk under ``results/enrichment/``; the one test that talks to
string-db.org is marked ``network`` and is excluded by the default ``addopts``.

Two things these tests exist to catch, both of which actually happened:

1. **Stale hardcoded counts.** ``load_seeds`` used to assert
   ``len(seeds) == 715 / n_up == 206 / n_down == 509`` as inline literals. The
   DECISIONS_LOG **D7** correction swapped UP and DOWN (206/509 -> 509/206) and
   two of those three literals silently became wrong -- the script would have
   aborted a legitimate re-run. The expectations now come from
   ``tests/expected/frozen_counts.json``;
   :func:`test_source_has_no_hardcoded_dataset_counts` keeps them there.
2. **An uncached response.** The ``/json/version`` probe was the only STRING
   call whose response was not written to ``results/enrichment/raw/``, which
   would leave a future offline-replay mode without a recorded version.
   :func:`test_get_string_version_caches_the_raw_response` pins the cache write
   using a stub, so it needs no network.

The topology itself (694 nodes / 5963 edges) is NOT re-derived here -- that is a
live STRING query and the freeze manifest's job. What is asserted is that the
node set is a subset of the seed set, that every node carries a direction, and
that the direction and log2FC on each node agree with ``foldchange_all.csv``.
"""

from __future__ import annotations

import ast
import json
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_PDE = _REPO_ROOT / "proteomics_de"
_RESULTS = _PDE / "results"
_ENRICHMENT = _RESULTS / "enrichment"
_RAW = _ENRICHMENT / "raw"

FOLDCHANGE_CSV = _RESULTS / "foldchange_all.csv"
NODE_METRICS_CSV = _ENRICHMENT / "string_node_metrics.csv"
EDGES_TSV = _ENRICHMENT / "string_edges.tsv"
META_JSON = _ENRICHMENT / "string_meta.json"
RAW_GET_IDS = _RAW / "string_get_ids.json"
RAW_NETWORK = _RAW / "string_network.tsv"
RAW_VERSION = _RAW / "string_version.json"

STRING_PPI_SRC = _PDE / "enrich" / "string_ppi.py"

#: Mus musculus. Confirmed by the user (DECISIONS_LOG D5); every STRING call in
#: the pipeline must carry it.
MOUSE_TAXID = 10090

ACCESSION_COL = "UniProt Accession Number"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def string_ppi():
    """The module under test. Imported, never executed -- ``main()`` is live."""
    sys.path.insert(0, str(_PDE / "enrich"))
    try:
        import string_ppi as mod
    finally:
        sys.path.pop(0)
    return mod


@pytest.fixture(scope="module")
def foldchange() -> pd.DataFrame:
    return pd.read_csv(FOLDCHANGE_CSV)


@pytest.fixture(scope="module")
def seeds(string_ppi) -> pd.DataFrame:
    return string_ppi.load_seeds()


@pytest.fixture(scope="module")
def node_metrics() -> pd.DataFrame:
    return pd.read_csv(NODE_METRICS_CSV)


@pytest.fixture(scope="module")
def meta() -> dict:
    return json.loads(META_JSON.read_text())


@pytest.fixture(scope="module")
def src() -> str:
    return STRING_PPI_SRC.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Seed set: exactly the regulated rows of foldchange_all.csv
# ---------------------------------------------------------------------------
def test_seed_set_is_exactly_the_regulated_rows(seeds, foldchange):
    """``load_seeds`` selects UP|DOWN and nothing else, by accession."""
    expected = foldchange[foldchange["regulated"].isin(["UP", "DOWN"])]
    assert set(seeds[ACCESSION_COL]) == set(expected[ACCESSION_COL])
    assert len(seeds) == len(expected)
    assert set(seeds["regulated"]) == {"UP", "DOWN"}


def test_no_nochange_or_onoff_protein_is_a_seed(seeds, foldchange):
    """The complement is genuinely excluded, not merely under-counted."""
    excluded = foldchange[~foldchange["regulated"].isin(["UP", "DOWN"])]
    assert set(seeds[ACCESSION_COL]).isdisjoint(set(excluded[ACCESSION_COL]))


def test_seed_accessions_are_unique(seeds):
    assert not seeds[ACCESSION_COL].duplicated().any()


def test_seed_counts_match_frozen_counts(seeds, frozen_counts):
    """715 seeds, 509 UP, 206 DOWN -- post-D7 (was 206 UP / 509 DOWN)."""
    n_up = int((seeds["regulated"] == "UP").sum())
    n_down = int((seeds["regulated"] == "DOWN").sum())
    assert len(seeds) == frozen_counts["ipa_input_rows"]
    assert n_up == frozen_counts["n_up"]
    assert n_down == frozen_counts["n_down"]
    assert n_up + n_down == len(seeds)


def test_frozen_counts_carry_the_post_d7_direction(frozen_counts):
    """Guards the direction itself: D7 made UP the larger class.

    If this ever fails, either the dataset changed (re-derive frozen_counts.json)
    or the control/treated assignment has been inverted again.
    """
    assert frozen_counts["n_up"] > frozen_counts["n_down"], (
        "post-D7 the UP class is the larger one (509 UP / 206 DOWN); "
        "n_up < n_down suggests the control/treated assignment is inverted again"
    )


def test_load_seeds_returns_the_columns_the_pipeline_uses(seeds):
    assert list(seeds.columns) == [ACCESSION_COL, "Gene names", "log2FC", "regulated"]


def test_seed_log2fc_signs_agree_with_the_direction_labels(seeds):
    """UP <=> log2FC >= +0.585, DOWN <=> log2FC <= -0.585."""
    up = seeds[seeds["regulated"] == "UP"]["log2FC"]
    down = seeds[seeds["regulated"] == "DOWN"]["log2FC"]
    assert (up >= 0.585).all()
    assert (down <= -0.585).all()


# ---------------------------------------------------------------------------
# The counts are read from frozen_counts.json, not typed into the source
# ---------------------------------------------------------------------------
def test_source_has_no_hardcoded_dataset_counts(src, frozen_counts):
    """No 715 / 509 / 206 numeric literal in executable code.

    This is the regression guard for the D7 breakage: the numbers belong in
    ``tests/expected/frozen_counts.json`` so a dataset change updates one file.
    Parsed with ``ast`` rather than grepped, so prose in docstrings and comments
    (which is allowed to name the counts) cannot trip it and, more importantly,
    cannot mask a real literal on the same line.
    """
    dataset_numbers = {frozen_counts[k] for k in ("ipa_input_rows", "n_up", "n_down")}
    tree = ast.parse(src)
    offenders = [
        f"line {node.lineno}: {node.value}"
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, int)
        and not isinstance(node.value, bool)
        and node.value in dataset_numbers
    ]
    assert not offenders, (
        "dataset counts are hardcoded in string_ppi.py; read them from "
        "tests/expected/frozen_counts.json instead:\n" + "\n".join(offenders)
    )


def test_load_seeds_reads_the_counts_from_frozen_counts_json(src):
    """...and reads them by the documented keys, not by index or guesswork."""
    assert "load_frozen_counts()" in src
    for key in ("ipa_input_rows", "n_up", "n_down"):
        assert f'"{key}"' in src or f"'{key}'" in src, f"{key} is never read"


def test_load_seeds_still_asserts_against_frozen_counts(string_ppi, monkeypatch):
    """De-hardcoding must not have quietly disabled the check.

    A wrong expectation has to raise. If this passes silently, ``load_seeds``
    is no longer validating anything.
    """
    real = string_ppi.core.load_frozen_counts

    def wrong(*a, **kw):
        counts = dict(real(*a, **kw))
        counts["n_up"] = 999999
        return counts

    monkeypatch.setattr(string_ppi.core, "load_frozen_counts", wrong)
    with pytest.raises(AssertionError, match="999999"):
        string_ppi.load_seeds()


def test_baseline_checks_are_on_by_default(string_ppi):
    assert string_ppi.core.baseline_checks_enabled({}) is True
    assert string_ppi.core.baseline_checks_enabled({"PDE_EXPECT_BASELINE": "0"}) is False


# ---------------------------------------------------------------------------
# Species: mouse, taxid 10090 (D5)
# ---------------------------------------------------------------------------
def test_module_species_is_mouse(string_ppi):
    assert string_ppi.SPECIES == MOUSE_TAXID


def test_meta_records_mouse(meta):
    assert meta["species"] == MOUSE_TAXID
    assert meta["species_name"] == "Mus musculus"


def test_every_string_call_carries_the_species(src):
    """No STRING request may go out without ``species``."""
    assert src.count('"species": SPECIES') >= 2


def test_raw_network_rows_are_all_mouse():
    """STRING's ``ncbiTaxonId`` column, checked on the cached response."""
    cols = pd.read_csv(RAW_NETWORK, sep="\t", header=None,
                       names=_network_columns())
    assert set(cols["ncbiTaxonId"].unique()) == {MOUSE_TAXID}


def _network_columns() -> list[str]:
    return [
        "stringId_A", "stringId_B", "preferredName_A", "preferredName_B",
        "ncbiTaxonId", "score", "nscore", "fscore", "pscore", "ascore",
        "escore", "dscore", "tscore",
    ]


# ---------------------------------------------------------------------------
# Rate limiting (STRING asks for >=1 s between calls)
# ---------------------------------------------------------------------------
def test_rate_limit_is_at_least_one_second(string_ppi):
    assert string_ppi.RATE_LIMIT_SECONDS >= 1.0


# ---------------------------------------------------------------------------
# Node metrics table
# ---------------------------------------------------------------------------
def test_node_metrics_columns(node_metrics):
    assert list(node_metrics.columns) == [
        "accession", "symbol", "string_id", "degree", "betweenness",
        "community", "log2FC", "regulated",
    ]


def test_node_metrics_row_count_matches_frozen_counts(node_metrics, frozen_counts):
    assert len(node_metrics) == frozen_counts["string_nodes"]


def test_every_node_has_a_direction(node_metrics):
    assert node_metrics["regulated"].notna().all()
    assert set(node_metrics["regulated"].unique()) <= {"UP", "DOWN"}
    assert set(node_metrics["regulated"].unique()) == {"UP", "DOWN"}


def test_node_directions_reflect_the_post_d7_orientation(node_metrics):
    """The mapped subset inherits the seed set's UP-heavy composition."""
    n_up = int((node_metrics["regulated"] == "UP").sum())
    n_down = int((node_metrics["regulated"] == "DOWN").sum())
    assert n_up + n_down == len(node_metrics)
    assert n_up > n_down, (
        "post-D7 the majority of network nodes are UP; a DOWN-heavy network "
        "means the node annotations predate the D7 correction"
    )


def test_node_direction_and_log2fc_agree_with_foldchange(node_metrics, foldchange):
    """Each node's direction/log2FC is the one its accession carries upstream.

    This is what actually breaks if the STRING outputs are regenerated from a
    stale seed table: identical topology, inverted annotations.
    """
    ref = foldchange.set_index(ACCESSION_COL)
    merged = node_metrics.join(
        ref[["regulated", "log2FC"]], on="accession", rsuffix="_src"
    )
    assert merged["regulated_src"].notna().all(), "a node accession is not in foldchange_all.csv"
    assert (merged["regulated"] == merged["regulated_src"]).all()
    assert (merged["log2FC"] - merged["log2FC_src"]).abs().max() < 1e-9


def test_node_accessions_are_unique_and_are_seeds(node_metrics, seeds):
    assert not node_metrics["accession"].duplicated().any()
    assert set(node_metrics["accession"]) <= set(seeds[ACCESSION_COL])


def test_node_symbols_and_string_ids_are_unique(node_metrics):
    """Duplicate labels would silently merge two proteins into one graph node."""
    assert not node_metrics["symbol"].duplicated().any()
    assert not node_metrics["string_id"].duplicated().any()
    assert node_metrics["symbol"].notna().all()


def test_node_metrics_are_in_range(node_metrics):
    assert (node_metrics["degree"] >= 0).all()
    assert node_metrics["betweenness"].between(0.0, 1.0).all()
    assert (node_metrics["community"] >= 0).all()


def test_node_metrics_sorted_by_degree_descending(node_metrics):
    assert node_metrics["degree"].is_monotonic_decreasing


# ---------------------------------------------------------------------------
# Edge list / meta consistency
# ---------------------------------------------------------------------------
def test_edges_parse_and_reference_known_nodes(node_metrics):
    edges = pd.read_csv(EDGES_TSV, sep="\t")
    assert list(edges.columns) == ["nodeA_symbol", "nodeB_symbol", "combined_score"]
    symbols = set(node_metrics["symbol"])
    assert set(edges["nodeA_symbol"]) <= symbols
    assert set(edges["nodeB_symbol"]) <= symbols
    assert (edges["nodeA_symbol"] != edges["nodeB_symbol"]).all(), "self-loop in edge list"


def test_edge_scores_respect_required_score(meta):
    edges = pd.read_csv(EDGES_TSV, sep="\t")
    assert (edges["combined_score"] * 1000 >= meta["required_score"] - 1).all()


def test_meta_matches_the_written_tables(meta, node_metrics):
    edges = pd.read_csv(EDGES_TSV, sep="\t")
    assert meta["n_nodes_in_network"] == len(node_metrics)
    assert meta["n_edges"] == len(edges)
    assert meta["n_mapped"] == len(node_metrics)
    assert 0 < meta["n_communities"] <= len(node_metrics)


def test_meta_seed_counts_match_frozen_counts(meta, frozen_counts):
    """The D7 swap has to have reached ``string_meta.json`` too."""
    assert meta["n_seeds"] == frozen_counts["ipa_input_rows"]
    assert meta["n_seeds_up"] == frozen_counts["n_up"]
    assert meta["n_seeds_down"] == frozen_counts["n_down"]
    assert meta["n_mapped"] <= meta["n_seeds"]


def test_meta_records_the_rate_limit_and_caveat(meta):
    assert meta["rate_limit_seconds_between_calls"] >= 1.0
    assert "hypothesis-generating" in meta["caveat"].lower()


# ---------------------------------------------------------------------------
# Cached raw STRING responses
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("path", [RAW_GET_IDS, RAW_NETWORK, RAW_VERSION])
def test_raw_response_exists_and_is_non_empty(path):
    assert path.exists(), f"missing cached STRING response: {path}"
    assert path.stat().st_size > 0


def test_raw_get_ids_is_valid_json_covering_the_nodes(node_metrics):
    records = json.loads(RAW_GET_IDS.read_text())
    assert isinstance(records, list) and records
    assert {"queryItem", "stringId"} <= set(records[0])
    assert set(node_metrics["string_id"]) <= {r["stringId"] for r in records}


def test_raw_version_is_valid_json_and_agrees_with_meta(meta):
    payload = json.loads(RAW_VERSION.read_text())
    assert isinstance(payload, list) and payload
    assert payload[0]["string_version"] == meta["string_version"]


def test_raw_network_is_valid_tsv_matching_the_edge_count():
    raw = pd.read_csv(RAW_NETWORK, sep="\t", header=None, names=_network_columns())
    edges = pd.read_csv(EDGES_TSV, sep="\t")
    deduped = raw.drop_duplicates(subset=["stringId_A", "stringId_B"])
    assert len(deduped) == len(edges)
    assert raw["score"].between(0.0, 1.0).all()


def test_get_string_version_caches_the_raw_response(string_ppi, tmp_path, monkeypatch):
    """The version probe writes ``raw/string_version.json`` -- no live call.

    Every other STRING response is cached; this one used to be discarded, so a
    future ``--offline`` replay would have had no version to replay.
    """
    payload = [{"string_version": "99.9", "stable_address": "https://example.invalid"}]

    class _Resp:
        def raise_for_status(self):
            return None

        def json(self):
            return payload

    called = {}

    def fake_get(url, **kw):
        called["url"] = url
        return _Resp()

    out = tmp_path / "string_version.json"
    monkeypatch.setattr(string_ppi.requests, "get", fake_get)
    monkeypatch.setattr(string_ppi, "RAW_VERSION_OUT", out)

    assert string_ppi.get_string_version() == "99.9"
    assert called["url"].endswith("/json/version")
    assert json.loads(out.read_text()) == payload


# ---------------------------------------------------------------------------
# Live STRING (excluded from the default run)
# ---------------------------------------------------------------------------
@pytest.mark.network
def test_live_string_version_matches_the_committed_meta(string_ppi, meta, tmp_path,
                                                        monkeypatch):
    """One live call to ``/json/version``; does not touch results/."""
    monkeypatch.setattr(string_ppi, "RAW_VERSION_OUT", tmp_path / "v.json")
    assert string_ppi.get_string_version() == meta["string_version"]
