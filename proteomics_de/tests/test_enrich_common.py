"""Contract tests for the enrichment layer's shared plumbing.

Three things are locked here, all of them offline. Nothing in this file makes a
network call; the live g:Profiler / Enrichr behaviour is exercised by running
``enrich/ora.py`` and ``enrich/gsea.py`` themselves, and the single
network-marked test below only checks that their *cached* responses on disk are
the ones the committed results were built from.

1. **The accession policy is shared, not duplicated.**
   ``enrich_common._first_token_or_accession`` used to inline
   ``str(acc).split(";")[0].strip()``. It now calls
   ``proteomics_de.etl.accessions.first_token``. The test asserts the two are
   the same function of their input -- over hand-picked edge cases *and* over
   every accession in the two committed tables -- rather than trusting a
   comment that says "behaviour-identical".

2. **The dataset counts come from one file.**
   ``tests/expected/frozen_counts.json`` is the only place ``n_up`` / ``n_down``
   / ``background_union`` are written down. DECISIONS_LOG D7 is exactly why:
   the control/treated assignment shipped inverted, 206 UP / 509 DOWN became
   509 UP / 206 DOWN, and inline literals in ``enrich_common.py`` were what made
   that a multi-file edit. If these ever disagree again, this test says so.

3. **Empty is the correct answer, and that is a contract.**
   ``ora_up.csv`` and ``ora_down.csv`` are header-only ON PURPOSE -- 0
   GO/KEGG/Reactome terms survive g:SCS correction against the honest
   detected-proteome background (DECISIONS_LOG D6, the report's headline
   finding). ``test_ora_direction_csv_is_header_only`` therefore asserts
   *exactly zero* data rows with the full 8-column header intact. A future
   change that starts emitting rows here has changed the science, not fixed a
   bug, and must be argued for in DECISIONS_LOG -- not absorbed by loosening
   this assertion.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPO_ROOT = _PKG_DIR.parent

# The pipeline scripts are flat modules; mirror the sys.path layout they expect
# (conftest.py does this too, but this file must also work standalone).
for _entry in (_PKG_DIR / "viz", _PKG_DIR / "enrich", _PKG_DIR, _REPO_ROOT):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))

import enrich_common as ec  # noqa: E402
from proteomics_de.etl import accessions  # noqa: E402

_RESULTS = _PKG_DIR / "results"
_ENRICHMENT = _RESULTS / "enrichment"
_RAW = _ENRICHMENT / "raw"

#: The eight columns of ora_{up,down}.csv, in order. Duplicated from ora.py on
#: purpose: a test that imports the constant it is checking cannot detect a
#: change to it.
ORA_CSV_COLUMNS = [
    "source",
    "term_id",
    "term_name",
    "p_value",
    "term_size",
    "query_size",
    "intersection_size",
    "intersecting_genes",
]


def _old_inline_first_token(acc):
    """The exact expression ``enrich_common.py`` carried before the refactor.

    Kept verbatim as the oracle for the equivalence test below.
    """
    return str(acc).split(";")[0].strip()


# ---------------------------------------------------------------------------
# 1. first_token equivalence with the inline logic it replaced
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "acc",
    [
        "P12345",
        "P08752;P20612;Q9DC51",
        "  P12345  ",
        "P12345;",
        ";P12345",
        "",
        ";",
        "P12345-2",
        "12345;678",  # junk index list -- first_token must not special-case it
        float("nan"),  # str(nan) == "nan"; both forms must agree on that
        123,
        None,
    ],
)
def test_first_token_matches_the_inline_logic_it_replaced(acc):
    assert accessions.first_token(acc) == _old_inline_first_token(acc)


def test_first_token_matches_inline_logic_on_every_committed_accession():
    """Equivalence over the real data, not just hand-picked edge cases."""
    fc = pd.read_csv(_RESULTS / "foldchange_all.csv")
    scp = pd.read_csv(_RESULTS / "single_condition_proteins.csv")
    values = list(fc["UniProt Accession Number"]) + list(scp["accession"])
    # Read, not typed: DECISIONS_LOG D11 moved this from 2554 to 2552 when the
    # 2 junk accessions were quarantined out of single_condition_proteins.csv.
    expected = ec.load_frozen_counts()["background_union"]
    assert len(values) == expected, (
        f"expected the full {expected}-row detected proteome, got {len(values)}"
    )
    mismatches = [
        v for v in values if accessions.first_token(v) != _old_inline_first_token(v)
    ]
    assert mismatches == []


def test_enrich_common_uses_the_shared_policy_object():
    """Not a copy of the helper -- the same function object."""
    assert ec.first_token is accessions.first_token


def test_accession_fallback_path_routes_through_first_token():
    """A blank gene name falls back to the accession's first token."""
    sym, was_split, used_fallback = ec._first_token_or_accession(
        "", "P08752;P20612;Q9DC51"
    )
    assert (sym, was_split, used_fallback) == ("P08752", False, True)


def test_gene_name_path_is_unaffected_by_the_refactor():
    assert ec._first_token_or_accession("Lama1;Lama2", "P12345") == (
        "Lama1",
        True,
        False,
    )
    assert ec._first_token_or_accession("Lama1", "P12345") == ("Lama1", False, False)


# ---------------------------------------------------------------------------
# 2. Counts come from frozen_counts.json, and the code agrees with it
# ---------------------------------------------------------------------------
def test_no_dataset_counts_are_hardcoded_in_enrich_common():
    """The D7 literals must not creep back in.

    206/509 were the pre-D7 UP/DOWN counts and 2554 is the background union;
    all three were once inline ``assert len(...) == N`` literals here. Prose
    may mention them (the docstrings explain the correction); executable code
    may not.

    Checked over the AST rather than over the text: comments never reach the
    AST and docstrings arrive as ``str`` constants, so only genuine numeric
    literals are inspected. A text scan is what a naive version of this test
    would do, and it is exactly what silently passes on ``_n = 206`` when the
    docstring-stripping is off by one.
    """
    import ast

    tree = ast.parse(
        (_PKG_DIR / "enrich" / "enrich_common.py").read_text(encoding="utf-8")
    )
    forbidden = {206, 509, 2554}
    found = [
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant)
        and isinstance(node.value, int)
        and not isinstance(node.value, bool)
        and node.value in forbidden
    ]
    assert found == [], (
        f"dataset counts {sorted(set(found))} are hardcoded as numeric literals "
        "in enrich_common.py; they belong in tests/expected/frozen_counts.json"
    )


def test_load_frozen_counts_reads_the_shared_file(frozen_counts):
    assert ec.load_frozen_counts() == frozen_counts


def test_query_and_background_sizes_match_frozen_counts(frozen_counts):
    background, up, down, list_meta = ec.load_background_and_queries()
    assert len(up) == frozen_counts["n_up"] == 509
    assert len(down) == frozen_counts["n_down"] == 206
    assert list_meta["background_row_union_n"] == frozen_counts["background_union"]
    # 2552 since DECISIONS_LOG D11 (was 2554): the 2 quarantined junk
    # accessions were never real proteins and so were never real background.
    assert list_meta["background_row_union_n"] == 2552
    # The background is a deduplicated symbol list, so it is necessarily
    # smaller than the row union it is built from.
    assert len(background) < list_meta["background_row_union_n"]
    assert len(background) == list_meta["background_unique_symbols_n"]


def test_up_and_down_queries_are_disjoint():
    """A protein cannot be both UP and DOWN in the same contrast."""
    _, up, down, _ = ec.load_background_and_queries()
    assert set(up).isdisjoint(set(down))


def test_queries_are_drawn_from_the_background():
    background, up, down, _ = ec.load_background_and_queries()
    bg = set(background)
    assert set(up) <= bg
    assert set(down) <= bg


def test_up_is_larger_than_down_after_d7():
    """Direction guard for the D7 correction.

    Pre-D7 the pipeline had 206 UP / 509 DOWN. Post-D7 it is 509 UP / 206 DOWN.
    A regression that re-inverts the assignment would restore the old ratio, and
    a count-only check against frozen_counts.json would still pass if both the
    data and the expectations were flipped together. This pins the direction.
    """
    _, up, down, _ = ec.load_background_and_queries()
    assert len(up) > len(down)


# ---------------------------------------------------------------------------
# 3. ora_meta.json records the organism and the custom background
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def ora_meta() -> dict:
    return json.loads((_ENRICHMENT / "ora_meta.json").read_text(encoding="utf-8"))


def test_ora_meta_records_mouse(ora_meta):
    """Organism is mouse -- confirmed by the user, DECISIONS_LOG D5 CLOSED."""
    assert ora_meta["organism"] == "mmusculus"
    assert ec.ORGANISM_GPROFILER == "mmusculus"


def test_ora_meta_records_the_custom_detected_proteome_background(
    ora_meta, frozen_counts
):
    """The custom background IS the finding (D6), not an implementation detail.

    Switching to g:Profiler's default whole-genome background turns this same
    UP query into 196 nominally "significant" terms -- a background-inflation
    artifact. ``domain_scope`` must stay ``"custom"``.
    """
    assert ora_meta["domain_scope"] == "custom"
    assert ora_meta["background_size_row_union"] == frozen_counts["background_union"]
    assert ora_meta["background_size_row_union"] == 2554
    assert ora_meta["background_size_unique_symbols"] < 2554


def test_ora_meta_query_sizes_match_frozen_counts(ora_meta, frozen_counts):
    assert ora_meta["up_n"] == frozen_counts["n_up"]
    assert ora_meta["down_n"] == frozen_counts["n_down"]


def test_ora_meta_thresholds(ora_meta):
    assert ora_meta["significance_threshold_method"] == "g_SCS"
    assert ora_meta["user_threshold"] == 0.05


def test_ora_meta_reports_zero_significant_terms(ora_meta):
    """D6: 0 terms pass g:SCS correction in either direction. See the module
    docstring -- this is the finding, not a failure."""
    assert ora_meta["n_terms_up"] == 0
    assert ora_meta["n_terms_down"] == 0


def test_ora_best_subthreshold_terms_do_not_pass_correction(ora_meta):
    """The best term in each direction is still above 0.05 by a wide margin."""
    for key in ("best_subthreshold_term_up", "best_subthreshold_term_down"):
        best = ora_meta[key]
        assert best is not None, f"{key} missing"
        assert best["p_value"] > 0.05, (
            f"{key} has g:SCS p={best['p_value']} -- something now PASSES "
            "correction where nothing did before. That is a scientific finding "
            "and must be reviewed against DECISIONS_LOG D6, not absorbed here."
        )


# ---------------------------------------------------------------------------
# 4. ora_{up,down}.csv are header-only BY DESIGN
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("direction", ["up", "down"])
def test_ora_direction_csv_is_header_only(direction):
    """Zero data rows with the full 8-column header. CORRECT per D6.

    Emptiness here is the result, not a missing file and not a failed run: no
    GO/KEGG/Reactome term survives g:SCS correction against the detected
    proteome. Asserted as a contract in both directions -- shape intact,
    content empty.
    """
    path = _ENRICHMENT / f"ora_{direction}.csv"
    assert path.exists(), f"{path} is missing -- ora.py must still write it"

    raw = path.read_text(encoding="utf-8")
    assert raw.strip() != "", "header-only is not the same as zero-byte"

    df = pd.read_csv(path)
    assert list(df.columns) == ORA_CSV_COLUMNS
    assert len(df.columns) == 8
    assert len(df) == 0, (
        f"ora_{direction}.csv has {len(df)} rows; D6 says it must have 0. A "
        "term now passing g:SCS<0.05 is a major finding -- report it, do not "
        "relax this test."
    )


# ---------------------------------------------------------------------------
# 5. Cached raw API responses exist and parse
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "filename",
    [
        "gprofiler_up.json",
        "gprofiler_down.json",
        "gprofiler_up_all.json",
        "gprofiler_down_all.json",
    ],
)
def test_cached_gprofiler_response_parses(filename):
    path = _RAW / filename
    assert path.exists(), f"{path} is missing -- every outbound call is cached"
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert "result" in payload
    assert payload["meta"]["query_metadata"]["organism"] == "mmusculus"
    assert payload["meta"]["query_metadata"]["domain_scope"] == "custom"


@pytest.mark.parametrize(
    "filename",
    ["gprofiler_up_all_evidence.json.gz", "gprofiler_down_all_evidence.json.gz"],
)
def test_cached_gprofiler_evidence_response_parses(filename):
    """The gene-level responses are gzipped (tens of MB uncompressed)."""
    path = _RAW / filename
    assert path.exists(), f"{path} is missing"
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        payload = json.load(fh)
    assert payload["meta"]["query_metadata"]["organism"] == "mmusculus"
    assert payload["result"], "the all_results=True response should carry terms"


@pytest.mark.parametrize("direction", ["up", "down"])
def test_cached_gprofiler_query_matches_the_current_query_set(direction):
    """The cached response is the one the current UP/DOWN lists produced.

    This is what proves the committed results are not stale after D7: if the
    query set on disk still had 206 genes under ``up``, the cache would predate
    the correction.
    """
    _, up, down, _ = ec.load_background_and_queries()
    expected = {"up": up, "down": down}[direction]
    payload = json.loads((_RAW / f"gprofiler_{direction}.json").read_text())
    queries = payload["meta"]["query_metadata"]["queries"]
    cached = queries.get("query_1") or next(iter(queries.values()))
    assert cached == expected


def test_cached_gprofiler_evidence_vector_aligns_with_mapped_genes():
    """Guards the intersections-alignment fix in ora.py.

    g:Profiler indexes each term's ``intersections`` vector by the genes it
    successfully MAPPED (509 submitted -> 473 mapped for UP), not by the
    submitted query. Zipping against the submitted list mislabels every gene
    after the first unmapped symbol.
    """
    with gzip.open(_RAW / "gprofiler_up_all_evidence.json.gz", "rt") as fh:
        payload = json.load(fh)
    entry = payload["meta"]["genes_metadata"]["query"]["query_1"]
    n_mapped = len(entry["ensgs"])
    n_submitted = len(payload["meta"]["query_metadata"]["queries"]["query_1"])
    assert n_mapped < n_submitted, (
        "some genes always fail to map; if they ever all map, the two orders "
        "coincide by luck and this guard stops being meaningful"
    )
    lengths = {len(t["intersections"]) for t in payload["result"] if "intersections" in t}
    assert lengths == {n_mapped}, (
        f"intersections vectors are length {lengths}, mapped genes = {n_mapped}"
    )


def test_cached_enrichr_libraries_parse(frozen_counts):
    """Both Enrichr libraries are cached gzipped, with their term counts."""
    path = _RAW / "enrichr_libraries.json.gz"
    assert path.exists(), f"{path} is missing"
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        payload = json.load(fh)
    libs = payload["libraries"]
    assert set(libs) == {"GO_Biological_Process_2021", "KEGG_2019_Mouse"}
    assert len(libs["GO_Biological_Process_2021"]) == 6034
    assert len(libs["KEGG_2019_Mouse"]) == 303
    assert payload["_n_terms_per_library"] == {
        name: len(terms) for name, terms in libs.items()
    }
    # every term maps to a non-empty gene list
    for name, terms in libs.items():
        assert all(isinstance(g, list) and g for g in terms.values()), name


# ---------------------------------------------------------------------------
# 6. GSEA is null too (D6's second half)
# ---------------------------------------------------------------------------
def test_gsea_meta_reports_zero_terms_passing_fdr(frozen_counts):
    meta = json.loads((_ENRICHMENT / "gsea_meta.json").read_text(encoding="utf-8"))
    assert meta["n_terms_tested_after_size_filter"] == frozen_counts["gsea_terms"]
    assert meta["n_fdr_lt_0.05"] == 0, (
        "a gene set now passes GSEA FDR<0.05 where none did before -- a major "
        "finding to report, not a test to relax"
    )
    assert meta["min_fdr_q_value"] > 0.05


def test_gsea_results_csv_agrees_with_its_meta(frozen_counts):
    df = pd.read_csv(_ENRICHMENT / "gsea_results.csv")
    assert len(df) == frozen_counts["gsea_terms"]
    assert (df["fdr_q_value"] >= 0.05).all()


def test_gsea_ranking_direction_follows_d7():
    """After D7 the score range is ~[-2.69, +3.80], the mirror of the pre-D7
    ~[-3.80, +2.69]. A re-inversion would flip it back."""
    meta = json.loads((_ENRICHMENT / "gsea_meta.json").read_text(encoding="utf-8"))
    lo, hi = meta["score_range"]
    assert lo < 0 < hi
    assert hi > abs(lo), (
        f"score range {meta['score_range']} is skewed the pre-D7 way; the "
        "control/treated assignment may have re-inverted"
    )


# ---------------------------------------------------------------------------
# Live re-query -- network-marked, excluded from the default run
# ---------------------------------------------------------------------------
@pytest.mark.network
def test_live_gprofiler_still_returns_nothing_significant():
    """Re-demonstrates D6 against the live API. Not run by default."""
    import ora  # noqa: PLC0415  -- imported lazily so the offline run stays offline

    background, up, _down, _ = ec.load_background_and_queries()
    resp = ora._post(up, background, all_results=False, no_evidences=True)
    significant = [t for t in resp["result"] if t.get("significant")]
    assert significant == [], (
        f"{len(significant)} terms now pass g:SCS<0.05 against the custom "
        "detected-proteome background; D6 needs revisiting"
    )
