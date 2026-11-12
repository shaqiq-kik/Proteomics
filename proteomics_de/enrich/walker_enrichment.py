"""
ORA (g:Profiler) and STRING over Dr. Walker's 721-protein list (DECISIONS_LOG D21, D24).

WHAT THIS IS
------------
``enrich/ora.py`` and ``enrich/string_ppi.py`` still run on the 715-protein core
list and write to ``results/enrichment/``. Dr. Walker's list
(``results/ipa_input_walker.csv``, D20/D24) has 721 rows: 715 core + FRZB + 5 control-only
proteins with an ASSIGNED fold change (the 247 recovered partial-tier proteins of the
first 968-protein version were dropped at Dr. Walker's request, D24; that version is in
git history, commits b9c75b9 and 9e345d9). This module runs the same two analyses on the 721,
ADDITIVELY:

* it imports the low-level pieces of ``ora.py`` / ``string_ppi.py`` (payload,
  term-row builder, g:Profiler gene-order resolver, node metric recipe) but
  never calls their ``main()`` and never writes into ``results/enrichment/``;
* everything goes to ``results/enrichment_walker/``;
* the 715 list is ALSO re-run, under the data versions current today, into
  ``results/enrichment_walker/control_715_same_version/``. The committed 715
  run used STRING 12.0 and g:Profiler e114_eg62_p19_27110d83; if those have moved
  since, a walker-vs-committed comparison would confound "the list changed" with
  "the databases changed". The same-version control separates the two.

WHICH INPUT, WHICH BACKGROUND
-----------------------------
Query: ``regulated`` (UP 510 / DOWN 211) of ``ipa_input_walker.csv``; one gene
symbol per row exactly as ``enrich_common.build_symbol_list`` resolves them (first
listed symbol, accession fallback). UP and DOWN are queried SEPARATELY, as in
ora.py.

Background: the SAME custom universe as the 715 run -- the measured proteome,
the 2552-row union of foldchange_all.csv and single_condition_proteins.csv
(2530 unique symbols; it is also the IPA Reference Set,
``ipa_background_measured.txt``, D19). That universe already contains every one
of the 721 query symbols, including the 5 control-only proteins (they are rows
of single_condition_proteins.csv) and FRZB (a row of foldchange_all.csv), so it
needs no extension. ``assert_queries_in_background``
fails loudly if that ever stops being true; it never silently extends or drops.

ASSIGNED VALUES ARE NEVER USED
------------------------------
ORA uses gene-list membership only. STRING uses membership only for the network;
node tables carry the UP/DOWN label and tier, and log2FC is left BLANK on the 5
``value_source == assigned`` rows (their -6.643856 is a placeholder, D20, not a
measurement), and no p-value of any kind is attached to any node. Nothing here
reads ``p_value`` / ``adj_p_value``.

OFFLINE REPLAY (how this stays testable)
----------------------------------------
Every API response is cached verbatim under ``<out>/raw/`` with a fingerprint of
the request that produced it (``raw/request_manifest.json``). The DEFAULT mode is
replay: read the cache, verify each response belongs to the request about to be
made, write the outputs, make NO network call. A missing or mismatched cache is an
error, never a silent fallback to the network. Live calls happen only with
``--fetch`` and only for responses not already cached. So the pipeline stage is
deterministic and offline, and the live run is an explicit, separate act.

Outputs (``results/enrichment_walker/``)
  ora_up.csv, ora_down.csv, ora_meta.json, ora_top_terms_detail.json,
  string_edges.tsv, string_node_metrics.csv, string_meta.json,
  string_focus_proteins.csv, string_comparison_vs_715.json,
  comparison_vs_715.csv, raw/, control_715_same_version/ (same files + raw/).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import requests
import yaml
from networkx.algorithms.community import greedy_modularity_communities
from scipy.stats import hypergeom

ENRICH_DIR = Path(__file__).resolve().parent
PDE_DIR = ENRICH_DIR.parent
REPO_ROOT = PDE_DIR.parent
for _p in (str(ENRICH_DIR), str(PDE_DIR), str(REPO_ROOT)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import enrich_common as ec  # noqa: E402
import ora  # noqa: E402
import string_ppi  # noqa: E402

RESULTS_DIR = PDE_DIR / "results"
CONFIG_PATH = PDE_DIR / "config" / "config.yaml"
OLD_ENRICHMENT_DIR = RESULTS_DIR / "enrichment"          # READ-ONLY here
WALKER_CSV = RESULTS_DIR / "ipa_input_walker.csv"

ACC_COL = "UniProt Accession Number"
GENE_COL = "Gene names"

VALUE_ASSIGNED = "assigned"
MANIFEST_NAME = "request_manifest.json"

COMPARISON_NAME = "comparison_vs_715.csv"
FOCUS_NAME = "string_focus_proteins.csv"
STRING_COMPARISON_NAME = "string_comparison_vs_715.json"

RUN_FILES = (
    "ora_up.csv", "ora_down.csv", "ora_meta.json", "ora_top_terms_detail.json",
    "string_edges.tsv", "string_node_metrics.csv", "string_meta.json",
)

WALKER_CAVEAT = (
    "n=2 biological replicates per group, from 2 separate experiments (low power); "
    "0/1938 proteins pass per-protein FDR<0.05 (limma). The list is Dr. Walker's "
    "721-protein regulated list (DECISIONS_LOG D20, D24): 715 core proteins and FRZB "
    "(which rests on the samples that were measured, with the missing value imputed, and "
    "is not FDR-significant), plus 5 control-only proteins whose 100-fold decrease is "
    "ASSIGNED, not measured. ORA and STRING use gene-list membership and direction only. "
    "Exploratory, hypothesis-generating."
)


class WalkerEnrichmentError(RuntimeError):
    """The inputs, the cache or the guard rails are not what the analysis needs."""


class CacheMiss(WalkerEnrichmentError):
    """A response is not cached and --fetch was not given."""


class CacheMismatch(WalkerEnrichmentError):
    """A cached response does not belong to the request about to be made."""


# --------------------------------------------------------------------------
# Config and inputs
# --------------------------------------------------------------------------
def load_config(path=None) -> dict:
    cfg = yaml.safe_load(Path(path or CONFIG_PATH).read_text(encoding="utf-8"))
    section = cfg.get("walker_enrichment")
    if not section:
        raise WalkerEnrichmentError("config.yaml has no walker_enrichment section")
    return section


def _frozen(counts_path=None) -> dict:
    return json.loads(
        Path(counts_path or ec.FROZEN_COUNTS_PATH).read_text(encoding="utf-8")
    )


def load_walker_list(results_dir=None, counts_path=None) -> pd.DataFrame:
    """The 721-row list, with the checks that make it safe to enrich."""
    results_dir = Path(results_dir or RESULTS_DIR)
    df = pd.read_csv(results_dir / "ipa_input_walker.csv", float_precision="round_trip")
    counts = _frozen(counts_path)
    if len(df) != counts["n_ipa_walker"]:
        raise WalkerEnrichmentError(
            f"expected {counts['n_ipa_walker']} rows in ipa_input_walker.csv "
            f"(frozen_counts.json:n_ipa_walker), got {len(df)}"
        )
    if df[ACC_COL].duplicated().any():
        raise WalkerEnrichmentError("duplicate accessions in ipa_input_walker.csv")
    if not df["regulated"].isin(["UP", "DOWN"]).all():
        raise WalkerEnrichmentError("ipa_input_walker.csv has a row that is not UP/DOWN")
    return df[[ACC_COL, GENE_COL, "log2FC", "regulated", "tier", "value_source"]].copy()


def load_core_list(results_dir=None, counts_path=None) -> pd.DataFrame:
    """The 715 core seeds exactly as string_ppi.load_seeds / ora select them."""
    results_dir = Path(results_dir or RESULTS_DIR)
    fc = pd.read_csv(results_dir / "foldchange_all.csv")
    df = fc[fc["regulated"].isin(["UP", "DOWN"])].reset_index(drop=True)
    counts = _frozen(counts_path)
    if len(df) != counts["ipa_input_rows"]:
        raise WalkerEnrichmentError(
            f"expected {counts['ipa_input_rows']} core rows, got {len(df)}"
        )
    return df[[ACC_COL, GENE_COL, "log2FC", "regulated"]].copy()


def mask_assigned(df: pd.DataFrame) -> pd.DataFrame:
    """Blank log2FC on assigned rows. ORA/STRING must not see a placeholder as data."""
    out = df.copy()
    if "value_source" in out.columns:
        out.loc[out["value_source"] == VALUE_ASSIGNED, "log2FC"] = np.nan
    return out


def build_queries(df: pd.DataFrame):
    """(up_symbols, down_symbols, list_meta) -- same resolver as enrich_common."""
    up_df = df[df["regulated"] == "UP"]
    down_df = df[df["regulated"] == "DOWN"]
    up, up_split, up_fb = ec.build_symbol_list(up_df, GENE_COL, ACC_COL)
    down, down_split, down_fb = ec.build_symbol_list(down_df, GENE_COL, ACC_COL)
    meta = {
        "up_n": len(up), "up_semicolon_split_n": up_split, "up_accession_fallback_n": up_fb,
        "down_n": len(down), "down_semicolon_split_n": down_split,
        "down_accession_fallback_n": down_fb,
    }
    return up, down, meta


def load_background():
    """The measured-proteome universe, exactly as the 715 ORA run builds it."""
    background, _up, _down, bg_meta = ec.load_background_and_queries()
    return background, {k: v for k, v in bg_meta.items() if k.startswith("background")}


def assert_queries_in_background(up, down, background) -> dict:
    """Every query symbol must be inside the custom background (g:Profiler's contract).

    Raises rather than extending: a silently grown universe would change the
    statistics without the reader being told.
    """
    bg = set(background)
    missing_up = sorted(set(up) - bg)
    missing_down = sorted(set(down) - bg)
    if missing_up or missing_down:
        raise WalkerEnrichmentError(
            f"query symbols outside the custom background: UP {missing_up[:10]}, "
            f"DOWN {missing_down[:10]} -- the universe must be extended explicitly "
            "and documented; refusing to continue"
        )
    return {"n_query_not_in_background": 0, "background_extended": False}


# --------------------------------------------------------------------------
# Raw response cache with request fingerprints
# --------------------------------------------------------------------------
def request_sha256(payload) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def gprofiler_payload(query, background, all_results, no_evidences) -> dict:
    """The exact JSON body ora._post sends (asserted equal by a test)."""
    return {
        "organism": ec.ORGANISM_GPROFILER,
        "query": query,
        "sources": ora.SOURCES,
        "user_threshold": ora.USER_THRESHOLD,
        "significance_threshold_method": ora.SIG_METHOD,
        "domain_scope": "custom",
        "background": background,
        "no_evidences": no_evidences,
        "all_results": all_results,
    }


def _gz_write_text(path: Path, text: str):
    # mtime=0: the gzip header otherwise carries the wall clock, and two writes of
    # the same payload would hash differently.
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
            gz.write(text.encode("utf-8"))


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return fh.read()
    return path.read_text(encoding="utf-8")


class RawStore:
    """``<dir>/raw`` cache. Replay by default; ``fetch=True`` fills only what is missing.

    ``require_manifest=False`` is for reading the COMMITTED 715 cache
    (``results/enrichment/raw``), which predates request fingerprints; it is only
    ever used read-only, for the equivalence test.
    """

    def __init__(self, raw_dir, *, fetch=False, require_manifest=True, writable=True):
        self.raw_dir = Path(raw_dir)
        self.fetch = fetch
        self.require_manifest = require_manifest
        self.writable = writable
        self.network_calls = 0
        self._manifest = {}
        mpath = self.raw_dir / MANIFEST_NAME
        if mpath.exists():
            self._manifest = json.loads(mpath.read_text(encoding="utf-8"))

    # -- lookup -----------------------------------------------------------
    def _candidates(self, name):
        base = name[:-3] if name.endswith(".gz") else name
        return [self.raw_dir / base, self.raw_dir / (base + ".gz")]

    def find(self, name):
        for p in self._candidates(name):
            if p.exists():
                return p
        return None

    def _check_fingerprint(self, path, req_sha, label):
        entry = self._manifest.get(path.name)
        if entry is None:
            if self.require_manifest:
                raise CacheMismatch(
                    f"{path.name} is cached but has no request fingerprint in "
                    f"{MANIFEST_NAME}; cannot prove it answers {label}"
                )
            return
        if entry.get("request_sha256") != req_sha:
            raise CacheMismatch(
                f"{path.name} was fetched for a different request than {label} "
                "(list, background, parameters or API base changed). Delete the "
                "raw/ cache deliberately and re-run with --fetch; this module "
                "never overwrites a cached response."
            )

    def _record(self, path, req_sha, endpoint):
        self._manifest[path.name] = {
            "endpoint": endpoint,
            "request_sha256": req_sha,
            "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        }
        (self.raw_dir / MANIFEST_NAME).write_text(
            json.dumps(self._manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    # -- generic get ------------------------------------------------------
    def get_text(self, name, req_sha, endpoint, fetch_fn, *, label):
        found = self.find(name)
        if found is not None:
            self._check_fingerprint(found, req_sha, label)
            return _read_text(found)
        if not self.fetch:
            raise CacheMiss(
                f"{self.raw_dir / name} is not cached for {label}; replay mode makes "
                "no network call. Re-run with --fetch to query the live API."
            )
        if not self.writable:
            raise WalkerEnrichmentError("read-only store cannot fetch")
        text = fetch_fn()
        self.network_calls += 1
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        target = self.raw_dir / name
        if name.endswith(".gz"):
            _gz_write_text(target, text)
        else:
            target.write_text(text, encoding="utf-8")
        self._record(target, req_sha, endpoint)
        print(f"[walker] cached raw response -> {target.name} "
              f"({target.stat().st_size / 1024:.0f} KB on disk)")
        return text

    def get_json(self, name, req_sha, endpoint, fetch_fn, *, label):
        text = self.get_text(
            name, req_sha, endpoint, lambda: json.dumps(fetch_fn(), indent=2), label=label
        )
        return json.loads(text)


# --------------------------------------------------------------------------
# ORA
# --------------------------------------------------------------------------
def _gprofiler_get(store, name, query, background, all_results, no_evidences, label):
    payload = gprofiler_payload(query, background, all_results, no_evidences)
    resp = store.get_json(
        name, request_sha256(payload), "gprofiler/gost/profile",
        lambda: ora._post(query, background, all_results, no_evidences),
        label=label,
    )
    sent = resp.get("meta", {}).get("query_metadata", {}).get("queries", {})
    got = sent.get("query_1")
    if got is not None and list(got) != list(query):
        raise CacheMismatch(
            f"{name}: the query recorded inside the g:Profiler response differs from "
            f"the {len(query)} symbols about to be analysed"
        )
    return resp


@dataclass
class OraResult:
    sig_rows: dict
    all_rows: dict
    version: str
    meta: dict
    detail: list


def run_ora(up, down, background, store, out_dir, *, label, list_meta, bg_meta,
            description) -> OraResult:
    """The g:Profiler stage of ora.py, over ``store`` instead of results/enrichment/raw."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    queries = {"up": up, "down": down}
    sig_rows, all_rows, primaries = {}, {}, {}
    for d, q in queries.items():
        lab = f"{label} {d.upper()} ORA"
        primary = _gprofiler_get(store, f"gprofiler_{d}.json", q, background,
                                 False, False, lab)
        qorder = ora._query_order(primary)
        sig = ora._rows_from_result(primary["result"], qorder, only_significant=True)
        df = pd.DataFrame(sig, columns=ora.CSV_COLUMNS + ["significant"])[ora.CSV_COLUMNS]
        df.to_csv(out_dir / f"ora_{d}.csv", index=False)
        full = _gprofiler_get(store, f"gprofiler_{d}_all.json.gz", q, background,
                              True, True, lab + " (all terms)")
        rows = ora._rows_from_result(full["result"], ora._query_order(full),
                                     only_significant=False)
        sig_rows[d], all_rows[d], primaries[d] = sig, rows, primary

    version = primaries["up"]["meta"].get("version")
    if primaries["down"]["meta"].get("version") != version:
        raise WalkerEnrichmentError("UP and DOWN responses carry different g:Profiler versions")

    # Top 5 per direction, selected exactly as ora.py does (stable sort by
    # corrected p; ties at 1.0 stay in API order -- NOT a meaningful ranking).
    top_pooled = []
    for d in ("up", "down"):
        for r in all_rows[d][:5]:
            top_pooled.append({**r, "direction": d.upper()})
    gene_maps = {}
    for d in ("up", "down"):
        ids = [r["term_id"] for r in top_pooled if r["direction"] == d.upper()]
        if not ids:
            gene_maps[d] = {}
            continue
        resp = _gprofiler_get(store, f"gprofiler_{d}_all_evidence.json.gz", queries[d],
                              background, True, False, f"{label} {d.upper()} ORA (evidence)")
        qorder = ora._query_order(resp)
        wanted, found = set(ids), {}
        for term in resp["result"]:
            if term["native"] not in wanted:
                continue
            inter = term.get("intersections") or []
            if len(inter) != len(qorder):
                raise WalkerEnrichmentError(
                    f"intersections vector for {term['native']} has {len(inter)} entries "
                    f"but {len(qorder)} mapped genes were resolved"
                )
            found[term["native"]] = [g for g, ev in zip(qorder, inter) if ev]
        gene_maps[d] = found

    detail = []
    for r in top_pooled:
        d = r["direction"].lower()
        detail.append({
            "direction": r["direction"], "source": r["source"], "term_id": r["term_id"],
            "term_name": r["term_name"], "p_value": r["p_value"],
            "term_size": r["term_size"], "query_size": r["query_size"],
            "intersection_size": r["intersection_size"],
            "intersecting_genes": gene_maps[d].get(r["term_id"], []),
            "significant_g_SCS_0.05": bool(r["significant"]),
        })
    n_sig = {d: len(sig_rows[d]) for d in sig_rows}
    (out_dir / "ora_top_terms_detail.json").write_text(json.dumps({
        "note": (
            f"{description}. Top 5 terms from EACH direction, selected exactly as "
            "enrich/ora.py selects them: stable sort by g:SCS-corrected p_value within "
            "direction. Many corrected p-values are exactly 1.0, and ties keep the API's "
            "order, so a block of p=1.0 terms is NOT a ranking -- use "
            f"{COMPARISON_NAME} (which breaks ties with the uncorrected hypergeometric "
            f"p) for an ordered view. Terms significant at g:SCS<0.05: UP {n_sig['up']}, "
            f"DOWN {n_sig['down']}; see significant_g_SCS_0.05 per term. Hypothesis-"
            "generating only."
        ),
        "terms": detail,
    }, indent=2))

    failed = {}
    ambiguous = {}
    queried = {}
    for d in ("up", "down"):
        gm = primaries[d]["meta"]["genes_metadata"]
        failed[d] = list(gm.get("failed", []))
        ambiguous[d] = len(gm.get("ambiguous", {}))
        queried[d] = len(gm["query"]["query_1"]["ensgs"])
    resp_ts = primaries["up"]["meta"].get("timestamp")
    meta = {
        "list": label,
        "description": description,
        "gprofiler_version": version,
        "gprofiler_response_timestamp_utc": resp_ts,
        "organism": ec.ORGANISM_GPROFILER,
        "sources": ora.SOURCES,
        "user_threshold": ora.USER_THRESHOLD,
        "significance_threshold_method": ora.SIG_METHOD,
        "domain_scope": "custom",
        "background_size_unique_symbols": len(background),
        "background_size_row_union": bg_meta["background_row_union_n"],
        "n_terms_up": n_sig["up"],
        "n_terms_down": n_sig["down"],
        "n_query_submitted_up": len(up),
        "n_query_submitted_down": len(down),
        "n_query_mapped_up": queried["up"],
        "n_query_mapped_down": queried["down"],
        "n_failed_to_map_up": len(failed["up"]),
        "n_failed_to_map_down": len(failed["down"]),
        "n_ambiguous_up": ambiguous["up"],
        "n_ambiguous_down": ambiguous["down"],
        "failed_to_map_symbols_up": failed["up"],
        "failed_to_map_symbols_down": failed["down"],
        "note": (
            "Same method as the 715 run (enrich/ora.py): g:Profiler g:GOSt, mmusculus, "
            "GO:BP/MF/CC+KEGG+REAC, g:SCS, user_threshold 0.05, domain_scope=custom with "
            "the measured proteome as background; UP and DOWN queried separately. "
            "Symbols that g:Profiler cannot map (renamed/obsolete mouse symbols, "
            "Excel-mangled names) drop out of the query exactly as in the 715 run; "
            "see failed_to_map_symbols_*. g:Profiler also discards symbols it reports as "
            "ambiguous (n_ambiguous_*); submitted = mapped + failed + ambiguous. ORA uses "
            "list membership only."
        ),
        "best_subthreshold_term_up": all_rows["up"][0] if all_rows["up"] else None,
        "best_subthreshold_term_down": all_rows["down"][0] if all_rows["down"] else None,
        "top5_subthreshold_terms_up": all_rows["up"][:5],
        "top5_subthreshold_terms_down": all_rows["down"][:5],
        **list_meta,
        **{k: v for k, v in bg_meta.items()},
    }
    (out_dir / "ora_meta.json").write_text(json.dumps(meta, indent=2, default=str))
    return OraResult(sig_rows, all_rows, version, meta, detail)


# --------------------------------------------------------------------------
# STRING
# --------------------------------------------------------------------------
@dataclass
class StringResult:
    graph: nx.Graph
    node_metrics: pd.DataFrame
    edges: pd.DataFrame
    meta: dict
    unmapped: list


def _string_post(api_base, endpoint, params, timeout=90):
    r = requests.post(f"{api_base}/{endpoint}", data=params, timeout=timeout)
    r.raise_for_status()
    return r


def _string_get_version(api_base):
    r = requests.get(f"{api_base}/json/version", timeout=30)
    r.raise_for_status()
    return r.json()


def run_string(seeds: pd.DataFrame, store, out_dir, *, api_base, caller_identity,
               label, description, caveat) -> StringResult:
    """The STRING stage of string_ppi.main(), over ``store``.

    The graph/metric recipe below is string_ppi.main()'s, line for line; a test
    replays the committed 715 cache through it and demands byte-identical
    string_edges.tsv / string_node_metrics.csv.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seeds = seeds.reset_index(drop=True)
    accessions = seeds[ACC_COL].tolist()
    lab = f"{label} STRING"

    ver_payload = {"endpoint": "json/version", "api_base": api_base}
    ver_text = store.get_text(
        "string_version.json", request_sha256(ver_payload), "string/json/version",
        lambda: json.dumps(_string_get_version(api_base), indent=2),
        label=lab,
    )
    version_payload = json.loads(ver_text)
    string_version = version_payload[0]["string_version"]
    stable_address = version_payload[0].get("stable_address")

    ids_params = {
        "identifiers": "\r".join(accessions), "species": string_ppi.SPECIES,
        "echo_query": 1, "limit": 1, "caller_identity": caller_identity,
    }

    def _fetch_ids():
        time.sleep(string_ppi.RATE_LIMIT_SECONDS)
        return _string_post(api_base, "json/get_string_ids", ids_params).json()

    id_records = store.get_json(
        "string_get_ids.json",
        request_sha256({"endpoint": "json/get_string_ids", "api_base": api_base,
                        "params": ids_params}),
        "string/json/get_string_ids", _fetch_ids, label=lab,
    )
    stray = {r.get("queryItem") for r in id_records} - set(accessions)
    if stray:
        raise CacheMismatch(
            f"string_get_ids.json answers accessions that are not in this list: "
            f"{sorted(stray)[:5]}"
        )

    acc_to_string_id, acc_to_preferred = {}, {}
    for rec in id_records:
        acc = rec.get("queryItem")
        if acc and acc not in acc_to_string_id:
            acc_to_string_id[acc] = rec["stringId"]
            acc_to_preferred[acc] = rec.get("preferredName", "")

    n_seeds = len(accessions)
    n_mapped = len(acc_to_string_id)
    mapped = seeds[seeds[ACC_COL].isin(acc_to_string_id)].copy()
    mapped["string_id"] = mapped[ACC_COL].map(acc_to_string_id)
    mapped["string_preferred_name"] = mapped[ACC_COL].map(acc_to_preferred)

    n_missing_gene_name = int(mapped[GENE_COL].isna().sum())
    mapped["symbol"] = mapped[GENE_COL]
    fb = mapped["symbol"].isna() | (mapped["symbol"].astype(str).str.strip() == "")
    mapped.loc[fb, "symbol"] = mapped.loc[fb, "string_preferred_name"]
    fb = mapped["symbol"].isna() | (mapped["symbol"].astype(str).str.strip() == "")
    mapped.loc[fb, "symbol"] = mapped.loc[fb, ACC_COL]
    dup_symbols = int(mapped["symbol"].duplicated().sum())
    if dup_symbols:
        mapped["symbol"] = np.where(
            mapped["symbol"].duplicated(keep=False),
            mapped["symbol"] + "_" + mapped[ACC_COL],
            mapped["symbol"],
        )

    string_ids = mapped["string_id"].tolist()
    net_params = {
        "identifiers": "\r".join(string_ids), "species": string_ppi.SPECIES,
        "required_score": string_ppi.REQUIRED_SCORE, "caller_identity": caller_identity,
    }

    def _fetch_net():
        time.sleep(string_ppi.RATE_LIMIT_SECONDS)
        return _string_post(api_base, "tsv-no-header/network", net_params).text

    net_text = store.get_text(
        "string_network.tsv",
        request_sha256({"endpoint": "tsv-no-header/network", "api_base": api_base,
                        "params": net_params}),
        "string/tsv-no-header/network", _fetch_net, label=lab,
    )
    if net_text.strip():
        from io import StringIO
        edges_raw = pd.read_csv(StringIO(net_text), sep="\t", header=None,
                                names=string_ppi.NETWORK_TSV_COLUMNS)
    else:
        edges_raw = pd.DataFrame(columns=string_ppi.NETWORK_TSV_COLUMNS)
    edges_raw = edges_raw.drop_duplicates(subset=["stringId_A", "stringId_B"])

    has_tier = "tier" in mapped.columns
    has_vs = "value_source" in mapped.columns
    G = nx.Graph()
    for _, row in mapped.iterrows():
        attrs = dict(
            accession=row[ACC_COL], symbol=row["symbol"],
            string_preferred_name=row["string_preferred_name"],
            log2FC=float(row["log2FC"]), regulated=row["regulated"],
        )
        if has_tier:
            attrs["tier"] = row["tier"]
        if has_vs:
            attrs["value_source"] = row["value_source"]
        G.add_node(row["string_id"], **attrs)
    for _, row in edges_raw.iterrows():
        a, b = row["stringId_A"], row["stringId_B"]
        if a == b or a not in G or b not in G:
            continue
        score = float(row["score"])
        G.add_edge(a, b, weight=score, combined_score=score)

    n_nodes, n_edges = G.number_of_nodes(), G.number_of_edges()
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G, weight=None, seed=string_ppi.SEED)
    if n_edges > 0:
        communities = list(greedy_modularity_communities(G, weight="weight"))
    else:
        communities = [{n} for n in G.nodes()]
    node_to_comm = {n: ci for ci, comm in enumerate(communities) for n in comm}

    rows = []
    for n, data in G.nodes(data=True):
        row = {
            "accession": data["accession"], "symbol": data["symbol"], "string_id": n,
            "degree": degree.get(n, 0), "betweenness": betweenness.get(n, 0.0),
            "community": node_to_comm.get(n, -1), "log2FC": data["log2FC"],
            "regulated": data["regulated"],
        }
        if has_tier:
            row["tier"] = data["tier"]
        if has_vs:
            row["value_source"] = data["value_source"]
        rows.append(row)
    node_metrics = pd.DataFrame(rows).sort_values(
        "degree", ascending=False).reset_index(drop=True)
    node_metrics.to_csv(out_dir / "string_node_metrics.csv", index=False)

    sym = {n: d["symbol"] for n, d in G.nodes(data=True)}
    edge_rows = [
        {"nodeA_symbol": sym.get(a, a), "nodeB_symbol": sym.get(b, b),
         "combined_score": d["combined_score"]}
        for a, b, d in G.edges(data=True)
    ]
    edges_df = pd.DataFrame(edge_rows, columns=["nodeA_symbol", "nodeB_symbol",
                                                "combined_score"]).sort_values(
        "combined_score", ascending=False).reset_index(drop=True)
    edges_df.to_csv(out_dir / "string_edges.tsv", sep="\t", index=False)

    unmapped = sorted(set(accessions) - set(acc_to_string_id))
    meta = {
        "list": label,
        "description": description,
        "string_version": string_version,
        "string_stable_address": stable_address,
        "string_api_base": api_base,
        "species": string_ppi.SPECIES,
        "species_name": "Mus musculus",
        "required_score": string_ppi.REQUIRED_SCORE,
        "n_seeds": n_seeds,
        "n_seeds_up": int((seeds["regulated"] == "UP").sum()),
        "n_seeds_down": int((seeds["regulated"] == "DOWN").sum()),
        "n_mapped": n_mapped,
        "mapping_rate": round(n_mapped / n_seeds, 4),
        "n_unmapped": len(unmapped),
        "unmapped_accessions": unmapped,
        "n_nodes_in_network": n_nodes,
        "n_edges": n_edges,
        "n_communities": len(communities),
        "n_seeds_missing_gene_name": n_missing_gene_name,
        "n_residual_duplicate_symbols": dup_symbols,
        "rate_limit_seconds_between_calls": string_ppi.RATE_LIMIT_SECONDS,
        "uses_values": (
            "list membership and UP/DOWN direction only; node log2FC is blank on "
            "value_source=assigned rows and no p-value is attached to any node"
        ),
        "caveat": caveat,
    }
    if has_vs:
        meta["n_seeds_assigned"] = int((seeds["value_source"] == VALUE_ASSIGNED).sum())
    (out_dir / "string_meta.json").write_text(json.dumps(meta, indent=2))
    return StringResult(G, node_metrics, edges_df, meta, unmapped)


# --------------------------------------------------------------------------
# Comparison tables
# --------------------------------------------------------------------------
def all_terms_table(store_or_dir, direction) -> pd.DataFrame:
    """Every term g:Profiler tested for ``direction``, ranked.

    Rank = position after sorting by (g:SCS-corrected p, uncorrected hypergeometric
    p, term_id). The corrected p is floored at exactly 1.0 for most terms, so the
    tie-break is what makes a rank meaningful; it is COMPUTED here from the
    response's own counts (effective_domain_size, term_size, query_size,
    intersection_size) as the upper-tail hypergeometric probability g:Profiler's
    test is built on. It is a ranking aid, never a significance claim.
    """
    raw_dir = Path(store_or_dir)
    base = f"gprofiler_{direction}_all.json"
    path = raw_dir / base if (raw_dir / base).exists() else raw_dir / (base + ".gz")
    resp = json.loads(_read_text(path))
    recs = []
    for t in resp["result"]:
        k, K = t["intersection_size"], t["term_size"]
        n, M = t["query_size"], t["effective_domain_size"]
        raw_p = float(hypergeom.sf(k - 1, M, K, n)) if k > 0 else 1.0
        recs.append({
            "source": t["source"], "term_id": t["native"], "term_name": t["name"],
            "adj_p": float(t["p_value"]), "raw_p": raw_p,
            "intersection_size": int(k), "term_size": int(K),
            "significant": bool(t.get("significant", False)),
        })
    df = pd.DataFrame(recs)
    df = df.sort_values(["adj_p", "raw_p", "term_id"], kind="mergesort").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)
    df.attrs["gprofiler_version"] = resp["meta"].get("version")
    return df


def _flag(in_a: bool, in_b: bool) -> str:
    if in_a and in_b:
        return "kept"
    if in_b:
        return "gained"
    if in_a:
        return "lost"
    return ""


def build_comparison(committed_raw, control_raw, walker_raw, top_n) -> pd.DataFrame:
    """Term-level old-vs-new table (see comparison_vs_715.csv in the README).

    committed715 = results/enrichment/raw (the committed 715 run, old data versions);
    control715   = the 715 list re-run today (same versions as walker721);
    walker721    = Dr. Walker's list.
    A term is listed if it is significant (g:SCS<0.05) or in the top ``top_n`` by rank
    in ANY of the three runs. ``sig_flag_*`` compares the significant sets;
    ``top{n}_flag_*`` compares top-n membership; both: kept / gained / lost, blank
    when the term is in neither set. "gained" = in the walker run only, "lost" = in
    the other run only.
    """
    runs = {"committed715": committed_raw, "control715": control_raw,
            "walker721": walker_raw}
    out = []
    for direction in ("up", "down"):
        tabs = {name: all_terms_table(path, direction).set_index("term_id")
                for name, path in runs.items()}
        ids = set()
        for t in tabs.values():
            ids |= set(t.index[t["significant"] | (t["rank"] <= top_n)])
        for tid in ids:
            rec = {"direction": direction.upper()}
            ref = next((t.loc[tid] for t in tabs.values() if tid in t.index), None)
            rec.update(source=ref["source"], term_id=tid, term_name=ref["term_name"])
            in_sig, in_top = {}, {}
            for name, t in tabs.items():
                if tid in t.index:
                    r = t.loc[tid]
                    rec[f"{name}_adj_p"] = r["adj_p"]
                    rec[f"{name}_hypergeom_p"] = r["raw_p"]
                    rec[f"{name}_rank"] = int(r["rank"])
                    rec[f"{name}_intersection_size"] = int(r["intersection_size"])
                    in_sig[name] = bool(r["significant"])
                    in_top[name] = bool(r["rank"] <= top_n)
                else:
                    for c in ("adj_p", "hypergeom_p", "rank", "intersection_size"):
                        rec[f"{name}_{c}"] = np.nan
                    in_sig[name] = in_top[name] = False
            for tag, other in (("committed", "committed715"), ("control", "control715")):
                rec[f"sig_flag_vs_{tag}"] = _flag(in_sig[other], in_sig["walker721"])
                rec[f"top{top_n}_flag_vs_{tag}"] = _flag(in_top[other], in_top["walker721"])
            rec["_sort"] = min(
                (rec[f"{n}_rank"] for n in runs if not math.isnan(rec[f"{n}_rank"])),
                default=10**9,
            )
            out.append(rec)
    df = pd.DataFrame(out)
    df = df.sort_values(["direction", "walker721_rank", "_sort", "term_id"],
                        ascending=[False, True, True, True], na_position="last",
                        kind="mergesort")
    df = df.drop(columns="_sort").reset_index(drop=True)
    for c in df.columns:
        if c.endswith(("_rank", "_intersection_size")):
            df[c] = df[c].astype("Int64")
    return df


COMPARISON_FLOAT_FORMAT = "%.10g"


def write_comparison(df: pd.DataFrame, path):
    df.to_csv(path, index=False, float_format=COMPARISON_FLOAT_FORMAT)


def focus_table(result: StringResult, focus: list) -> pd.DataFrame:
    """Where the named proteins sit in the network (degree, rank, neighbours)."""
    G = result.graph
    nm = result.node_metrics
    by_acc = {d["accession"]: n for n, d in G.nodes(data=True)}
    comm_size = nm.groupby("community").size().to_dict()
    deg_rank = nm["degree"].rank(ascending=False, method="min").astype(int)
    deg_rank.index = nm["accession"]
    recs = []
    for f in focus:
        acc, gene = f["accession"], f["gene"]
        sid = by_acc.get(acc)
        if sid is None:
            recs.append({"gene": gene, "accession": acc, "in_network": False,
                         "reason": "no STRING identifier for this accession"
                         if acc in result.unmapped else "not a node"})
            continue
        d = G.nodes[sid]
        row = nm[nm["accession"] == acc].iloc[0]
        nbrs = sorted(G.neighbors(sid),
                      key=lambda n: (-G[sid][n]["combined_score"], G.nodes[n]["symbol"]))
        recs.append({
            "gene": gene, "accession": acc, "in_network": True,
            "regulated": d["regulated"], "tier": d.get("tier", ""),
            "value_source": d.get("value_source", ""),
            "degree": int(row["degree"]),
            "degree_rank_of_n_nodes": f"{int(deg_rank[acc])}/{len(nm)}",
            "betweenness": float(row["betweenness"]),
            "community": int(row["community"]),
            "community_size": int(comm_size[int(row["community"])]),
            "n_neighbours_up": sum(G.nodes[n]["regulated"] == "UP" for n in nbrs),
            "n_neighbours_down": sum(G.nodes[n]["regulated"] == "DOWN" for n in nbrs),
            "n_neighbours_core": sum(G.nodes[n].get("tier") == "core" for n in nbrs),
            "n_neighbours_added": sum(G.nodes[n].get("tier") not in (None, "core")
                                      for n in nbrs),
            "neighbours": ";".join(
                f"{G.nodes[n]['symbol']}({G[sid][n]['combined_score']:g})" for n in nbrs),
        })
    return pd.DataFrame(recs)


def string_comparison(committed_dir, control: StringResult, walker: StringResult) -> dict:
    """Node/edge counts: committed 715, same-version 715, walker 721."""
    committed = json.loads((Path(committed_dir) / "string_meta.json").read_text())

    def pick(m):
        return {k: m[k] for k in ("string_version", "n_seeds", "n_mapped", "mapping_rate",
                                  "n_nodes_in_network", "n_edges", "n_communities")}

    core_nodes = set(walker.node_metrics.loc[walker.node_metrics["tier"] == "core", "symbol"])
    core_core = int(sum(1 for a, b in zip(walker.edges["nodeA_symbol"],
                                          walker.edges["nodeB_symbol"])
                        if a in core_nodes and b in core_nodes))
    added_nodes = int((walker.node_metrics["tier"] != "core").sum())
    return {
        "committed_715": pick(committed),
        "control_715_same_version": pick(control.meta),
        "walker_721": pick(walker.meta),
        "walker_nodes_added_tiers": added_nodes,
        "walker_edges_between_two_core_proteins": core_core,
        "walker_edges_touching_an_added_protein": int(walker.meta["n_edges"] - core_core),
        "control_715_edges": int(control.meta["n_edges"]),
        "core_core_edges_equal_control": core_core == int(control.meta["n_edges"]),
        "same_string_version_as_committed":
            committed["string_version"] == walker.meta["string_version"],
        "same_string_version_as_control":
            control.meta["string_version"] == walker.meta["string_version"],
    }


# --------------------------------------------------------------------------
# Guard rails
# --------------------------------------------------------------------------
def tree_hashes(directory) -> dict:
    directory = Path(directory)
    return {
        p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(directory.rglob("*")) if p.is_file()
    }


def assert_not_inside(out_dir, protected) -> None:
    out, prot = Path(out_dir).resolve(), Path(protected).resolve()
    if out == prot or prot in out.parents:
        raise WalkerEnrichmentError(f"{out} is inside protected directory {prot}")


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------
def run_list(frame, out_dir, *, background, bg_meta, fetch, api_base, caller_identity,
             label, description, caveat, results_dir=None):
    """ORA + STRING for one list into ``out_dir`` (+ ``out_dir/raw``)."""
    out_dir = Path(out_dir)
    store = RawStore(out_dir / "raw", fetch=fetch)
    frame = mask_assigned(frame)
    up, down, list_meta = build_queries(frame)
    cover = assert_queries_in_background(up, down, background)
    ora_res = run_ora(up, down, background, store, out_dir, label=label,
                      list_meta={**list_meta, **cover}, bg_meta=bg_meta,
                      description=description)
    str_res = run_string(frame, store, out_dir, api_base=api_base,
                         caller_identity=caller_identity, label=label,
                         description=description, caveat=caveat)
    return ora_res, str_res, store


def run(*, fetch=False, results_dir=None, cfg=None, control=None, sidecars=True,
        out_dir=None):
    results_dir = Path(results_dir or RESULTS_DIR)
    cfg = cfg or load_config()
    out_dir = Path(out_dir or (results_dir / Path(cfg["output_dir"]).name))
    old_dir = results_dir / "enrichment"
    assert_not_inside(out_dir, old_dir)
    control = cfg.get("control_715_rerun", True) if control is None else control
    control_dir = out_dir / cfg["control_dir_name"]

    before = tree_hashes(old_dir)

    background, bg_meta = load_background()
    walker_df = load_walker_list(results_dir)
    api_base = cfg["string"]["api_base"]
    who = cfg["string"]["caller_identity"]

    print(f"[walker] list: {len(walker_df)} rows "
          f"({(walker_df.regulated == 'UP').sum()} UP / {(walker_df.regulated == 'DOWN').sum()} DOWN); "
          f"mode: {'FETCH (live API for anything not cached)' if fetch else 'replay from cache'}")
    w_ora, w_str, w_store = run_list(
        walker_df, out_dir, background=background, bg_meta=bg_meta, fetch=fetch,
        api_base=api_base, caller_identity=who, label="walker721",
        description="Dr. Walker's 721-protein regulated list (D20, D24)", caveat=WALKER_CAVEAT,
    )

    c_ora = c_str = c_store = None
    if control:
        core_df = load_core_list(results_dir)
        c_ora, c_str, c_store = run_list(
            core_df, control_dir, background=background, bg_meta=bg_meta, fetch=fetch,
            api_base=api_base, caller_identity=who, label="control715_same_version",
            description=("the 715-protein core list re-run with the walker run's "
                         "API versions (control for database drift)"),
            caveat=(string_ppi_caveat()),
        )
        if c_ora.version != w_ora.version or c_str.meta["string_version"] != w_str.meta["string_version"]:
            print("[walker] WARNING: control and walker responses carry different "
                  "database versions; the same-version control is not same-version")

    comparison = None
    if control:
        comparison = build_comparison(
            old_dir / "raw", control_dir / "raw", out_dir / "raw", cfg["comparison"]["top_n"]
        )
        write_comparison(comparison, out_dir / COMPARISON_NAME)
        sc = string_comparison(old_dir, c_str, w_str)
        sc["gprofiler_version_committed_715"] = json.loads(
            (old_dir / "ora_meta.json").read_text())["gprofiler_version"]
        sc["gprofiler_version_control_715"] = c_ora.version
        sc["gprofiler_version_walker_721"] = w_ora.version
        sc["same_gprofiler_version_as_committed"] = (
            sc["gprofiler_version_committed_715"] == w_ora.version)
        (out_dir / STRING_COMPARISON_NAME).write_text(json.dumps(sc, indent=2))

    focus = focus_table(w_str, cfg["focus_proteins"])
    focus.to_csv(out_dir / FOCUS_NAME, index=False)

    after = tree_hashes(old_dir)
    if before != after:
        changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
        raise WalkerEnrichmentError(
            f"results/enrichment changed during the walker run: {changed[:5]}")

    if sidecars:
        write_sidecars(out_dir, w_ora, w_str, c_ora, c_str, comparison is not None)

    return {
        "out_dir": out_dir, "walker_ora": w_ora, "walker_string": w_str,
        "control_ora": c_ora, "control_string": c_str, "comparison": comparison,
        "focus": focus, "network_calls": w_store.network_calls + (c_store.network_calls if c_store else 0),
    }


def string_ppi_caveat() -> str:
    return (
        "n=2 biological replicates per group, from 2 separate experiments (low power); "
        "0/1938 proteins survive FDR<0.05 (limma). The 715 seed proteins are the "
        "|log2FC|>=0.585 HYPOTHESIS-GENERATING set. Same-version CONTROL run of the "
        "committed 715 list; exploratory."
    )


def write_sidecars(out_dir, w_ora, w_str, c_ora, c_str, has_comparison):
    """Provenance sidecars for the human-facing tables, written only when stale."""
    from proteomics_de.provenance import sha256_file, sidecar

    facts = {
        "list": "Dr. Walker's 721-protein list (DECISIONS_LOG D20/D21/D24)",
        "uses_values": "list membership and direction only; the 5 assigned log2FC "
                       "values (-6.643856, NOT measured) are never used",
        "assigned_values_are_not_measurements": True,
        "gprofiler_version": w_ora.version,
        "string_version": w_str.meta["string_version"],
        "enrichment_walker_caveat": WALKER_CAVEAT,
    }
    names = ["ora_up.csv", "ora_down.csv", "string_node_metrics.csv", "string_edges.tsv",
             FOCUS_NAME]
    if has_comparison:
        names.append(COMPARISON_NAME)
    for n in names:
        p = Path(out_dir) / n
        sc = p.with_name(p.name + ".provenance.json")
        if sc.exists():
            try:
                if json.loads(sc.read_text()).get("sha256") == sha256_file(p):
                    continue
            except json.JSONDecodeError:
                pass
        sidecar(p, **facts)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--fetch", action="store_true",
                    help="query g:Profiler / STRING live for responses that are not "
                         "cached (default: replay from raw/, no network)")
    ap.add_argument("--skip-control", action="store_true",
                    help="do not (re)build the same-version 715 control or the comparison")
    ap.add_argument("--no-sidecars", action="store_true")
    ap.add_argument("--results-dir", default=None)
    args = ap.parse_args(argv)

    out = run(fetch=args.fetch, results_dir=args.results_dir,
              control=False if args.skip_control else None, sidecars=not args.no_sidecars)
    w_ora, w_str = out["walker_ora"], out["walker_string"]
    print(f"\n[walker] g:Profiler {w_ora.version}; STRING {w_str.meta['string_version']}; "
          f"network calls this run: {out['network_calls']}")
    print(f"[walker] ORA significant terms: UP {w_ora.meta['n_terms_up']}, "
          f"DOWN {w_ora.meta['n_terms_down']}")
    print(f"[walker] STRING: {w_str.meta['n_mapped']}/{w_str.meta['n_seeds']} mapped, "
          f"{w_str.meta['n_nodes_in_network']} nodes, {w_str.meta['n_edges']} edges")
    print(f"[walker] outputs in {out['out_dir']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
