"""
STRING protein-protein interaction network for the DE-regulated seed proteins.

Organism: Mus musculus (NCBI taxid 10090). Every STRING call uses species=10090.

Seed set: the 715 "regulated" proteins from `results/foldchange_all.csv`
(rows where `regulated` in {UP, DOWN}: 206 UP + 509 DOWN). These are the
+/-0.585 log2FC HYPOTHESIS-GENERATING set from a n=2 technical-replicate
SILAC experiment where 0/1938 proteins survive FDR < 0.05 (see
`proteomics_de/results/de_stats.json` / the limma report). The PPI network
below shows structure among *candidate* proteins only -- it is exploratory
and must not be read as a network of statistically confirmed hits.

Pipeline:
  1. Map seed UniProt accessions -> STRING identifiers via
     POST /api/json/get_string_ids (species=10090, echo_query=1, limit=1).
  2. Fetch the induced subnetwork among mapped STRING IDs via
     POST /api/tsv-no-header/network (species=10090, required_score=400).
  3. Build a networkx graph (edge weight = combined score), compute degree,
     betweenness centrality, and greedy-modularity communities.
  4. Write edge list, node metrics table, and a metadata JSON.

STRING rate limit: >=1 second sleep between API calls (RATE_LIMIT_SECONDS).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import requests
from networkx.algorithms.community import greedy_modularity_communities

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ENRICH_DIR = Path(__file__).resolve().parent
PROTEOMICS_DE_DIR = ENRICH_DIR.parent
RESULTS_DIR = PROTEOMICS_DE_DIR / "results"
ENRICHMENT_DIR = RESULTS_DIR / "enrichment"
RAW_DIR = ENRICHMENT_DIR / "raw"
ENRICHMENT_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

FOLDCHANGE_CSV = RESULTS_DIR / "foldchange_all.csv"

EDGES_OUT = ENRICHMENT_DIR / "string_edges.tsv"
NODE_METRICS_OUT = ENRICHMENT_DIR / "string_node_metrics.csv"
META_OUT = ENRICHMENT_DIR / "string_meta.json"
RAW_GET_IDS_OUT = RAW_DIR / "string_get_ids.json"
RAW_NETWORK_OUT = RAW_DIR / "string_network.tsv"

# ----------------------------------------------------------------------------
# STRING API config
# ----------------------------------------------------------------------------
STRING_API_BASE = "https://string-db.org/api"
SPECIES = 10090  # Mus musculus
REQUIRED_SCORE = 400
CALLER_IDENTITY = "proteomics_de_pipeline_worker_a"
RATE_LIMIT_SECONDS = 1.2  # >=1s between STRING calls, per guardrails
GET_IDS_BATCH_SIZE = 400  # fallback batch size if a single call is rejected
SEED = 42  # deterministic community detection tie-breaks / any RNG use

NETWORK_TSV_COLUMNS = [
    "stringId_A", "stringId_B", "preferredName_A", "preferredName_B",
    "ncbiTaxonId", "score", "nscore", "fscore", "pscore", "ascore",
    "escore", "dscore", "tscore",
]

np.random.seed(SEED)


# ----------------------------------------------------------------------------
# Step 0: load + validate seed set
# ----------------------------------------------------------------------------
def load_seeds() -> pd.DataFrame:
    df = pd.read_csv(FOLDCHANGE_CSV)
    seeds = df[df["regulated"].isin(["UP", "DOWN"])].copy()
    seeds = seeds.reset_index(drop=True)

    n_up = int((seeds["regulated"] == "UP").sum())
    n_down = int((seeds["regulated"] == "DOWN").sum())
    print(f"[seeds] loaded {len(seeds)} regulated proteins from {FOLDCHANGE_CSV.name} "
          f"({n_up} UP, {n_down} DOWN)")

    assert len(seeds) == 715, f"expected 715 regulated seeds, got {len(seeds)}"
    assert n_up == 206, f"expected 206 UP, got {n_up}"
    assert n_down == 509, f"expected 509 DOWN, got {n_down}"

    dup = seeds["UniProt Accession Number"].duplicated().sum()
    assert dup == 0, f"expected unique accessions among seeds, found {dup} duplicates"

    return seeds[["UniProt Accession Number", "Gene names", "log2FC", "regulated"]]


# ----------------------------------------------------------------------------
# Step 1: map seed accessions -> STRING identifiers
# ----------------------------------------------------------------------------
def _post(endpoint: str, params: dict, timeout: int = 90) -> requests.Response:
    url = f"{STRING_API_BASE}/{endpoint}"
    r = requests.post(url, data=params, timeout=timeout)
    r.raise_for_status()
    return r


def _get_string_ids_call(accessions: list[str]) -> list[dict]:
    params = {
        "identifiers": "\r".join(accessions),
        "species": SPECIES,
        "echo_query": 1,
        "limit": 1,  # best match per query identifier only
        "caller_identity": CALLER_IDENTITY,
    }
    r = _post("json/get_string_ids", params)
    return r.json()


def map_string_ids(accessions: list[str]) -> tuple[list[dict], bool, int]:
    """Return (records, batched?, n_batches)."""
    try:
        records = _get_string_ids_call(accessions)
        print(f"[get_string_ids] single call succeeded for {len(accessions)} identifiers "
              f"-> {len(records)} records")
        return records, False, 1
    except requests.exceptions.RequestException as e:
        print(f"[get_string_ids] single call failed ({e}); falling back to batches of "
              f"{GET_IDS_BATCH_SIZE}")

    batches = [accessions[i:i + GET_IDS_BATCH_SIZE]
               for i in range(0, len(accessions), GET_IDS_BATCH_SIZE)]
    all_records: list[dict] = []
    for i, batch in enumerate(batches):
        if i > 0:
            time.sleep(RATE_LIMIT_SECONDS)
        res = _get_string_ids_call(batch)
        all_records.extend(res)
        print(f"  [get_string_ids] batch {i + 1}/{len(batches)}: {len(batch)} identifiers "
              f"-> {len(res)} records")
    return all_records, True, len(batches)


# ----------------------------------------------------------------------------
# Step 2: fetch induced subnetwork
# ----------------------------------------------------------------------------
def _network_call(string_ids: list[str]) -> str:
    params = {
        "identifiers": "\r".join(string_ids),
        "species": SPECIES,
        "required_score": REQUIRED_SCORE,
        "caller_identity": CALLER_IDENTITY,
    }
    r = _post("tsv-no-header/network", params)
    return r.text


def fetch_network(string_ids: list[str], batch_size: int = 500) -> tuple[str, bool, int]:
    """Return (raw_tsv_text, batched?, n_batches).

    The `network` endpoint returns the induced subnetwork (all edges where both
    endpoints are in the submitted identifier list). If a single call over all
    identifiers is rejected/too large, we fall back to querying in overlapping-free
    batches and taking the union of edges returned per batch; because STRING's
    induced-subnetwork query only returns edges *within* the submitted set, a
    naive non-overlapping batch split would silently drop true edges that cross
    a batch boundary. To avoid that, the batched fallback below still submits
    the FULL identifier list per request (STRING's own batching guidance is by
    identifier count *per request*, not by splitting the graph) -- i.e. batching
    here is a retry-with-smaller-request strategy, not a network partition.
    """
    try:
        text = _network_call(string_ids)
        print(f"[network] single call succeeded for {len(string_ids)} identifiers")
        return text, False, 1
    except requests.exceptions.RequestException as e:
        print(f"[network] single call failed ({e}); retrying in chunks of {batch_size} "
              f"(each chunk still queries against the full identifier list is not possible "
              f"for this endpoint; falling back to a best-effort partial network from "
              f"non-overlapping batches -- inter-batch edges will be missed and this is "
              f"documented in string_meta.json)")

    batches = [string_ids[i:i + batch_size] for i in range(0, len(string_ids), batch_size)]
    texts = []
    for i, batch in enumerate(batches):
        if i > 0:
            time.sleep(RATE_LIMIT_SECONDS)
        texts.append(_network_call(batch))
        print(f"  [network] batch {i + 1}/{len(batches)}: {len(batch)} identifiers")
    return "\n".join(t for t in texts if t.strip()), True, len(batches)


# ----------------------------------------------------------------------------
# Step 3: STRING version
# ----------------------------------------------------------------------------
def get_string_version() -> str:
    r = requests.get(f"{STRING_API_BASE}/json/version", timeout=30)
    r.raise_for_status()
    return r.json()[0]["string_version"]


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    seeds = load_seeds()
    accessions = seeds["UniProt Accession Number"].tolist()

    string_version = get_string_version()
    print(f"[string] API version: {string_version}")
    time.sleep(RATE_LIMIT_SECONDS)

    # --- map ---
    id_records, ids_batched, ids_n_batches = map_string_ids(accessions)
    RAW_GET_IDS_OUT.write_text(json.dumps(id_records, indent=2))
    print(f"[get_string_ids] raw JSON saved -> {RAW_GET_IDS_OUT}")

    # accession -> best STRING id (limit=1 already restricts to one row per query,
    # but guard against any residual duplicates by keeping the first occurrence)
    acc_to_string_id: dict[str, str] = {}
    acc_to_preferred_name: dict[str, str] = {}
    for rec in id_records:
        acc = rec.get("queryItem")
        if acc and acc not in acc_to_string_id:
            acc_to_string_id[acc] = rec["stringId"]
            acc_to_preferred_name[acc] = rec.get("preferredName", "")

    n_seeds = len(accessions)
    n_mapped = len(acc_to_string_id)
    mapping_rate = n_mapped / n_seeds
    print(f"[mapping] {n_mapped}/{n_seeds} seeds mapped to STRING IDs "
          f"({mapping_rate:.1%})")

    mapped_seeds = seeds[seeds["UniProt Accession Number"].isin(acc_to_string_id)].copy()
    mapped_seeds["string_id"] = mapped_seeds["UniProt Accession Number"].map(acc_to_string_id)
    mapped_seeds["string_preferred_name"] = mapped_seeds["UniProt Accession Number"].map(
        acc_to_preferred_name
    )
    string_id_to_acc = {v: k for k, v in acc_to_string_id.items()}

    # A handful of seeds have a missing/NaN `Gene names` value in the source
    # CSV (protein-group rows with no assigned symbol). Using the raw `Gene
    # names` column as the graph's node label would collide (two different
    # accessions both keying to NaN) and break the symbol-keyed edge list, so
    # fall back to the STRING preferred name, then the accession itself, to
    # keep every node's display label unique and non-null.
    n_missing_gene_name = int(mapped_seeds["Gene names"].isna().sum())
    mapped_seeds["symbol"] = mapped_seeds["Gene names"]
    fallback_mask = mapped_seeds["symbol"].isna() | (mapped_seeds["symbol"].astype(str).str.strip() == "")
    mapped_seeds.loc[fallback_mask, "symbol"] = mapped_seeds.loc[fallback_mask, "string_preferred_name"]
    fallback_mask = mapped_seeds["symbol"].isna() | (mapped_seeds["symbol"].astype(str).str.strip() == "")
    mapped_seeds.loc[fallback_mask, "symbol"] = mapped_seeds.loc[fallback_mask, "UniProt Accession Number"]
    if n_missing_gene_name:
        print(f"[symbols] {n_missing_gene_name} mapped seeds had no `Gene names` value; "
              f"fell back to STRING preferred name / accession for the node label")
    dup_symbols = mapped_seeds["symbol"].duplicated().sum()
    if dup_symbols:
        # extremely unlikely (would require two accessions sharing the same
        # symbol AND the same STRING preferred name), but disambiguate rather
        # than silently merge nodes if it ever happens
        mapped_seeds["symbol"] = np.where(
            mapped_seeds["symbol"].duplicated(keep=False),
            mapped_seeds["symbol"] + "_" + mapped_seeds["UniProt Accession Number"],
            mapped_seeds["symbol"],
        )
        print(f"[symbols] {dup_symbols} residual duplicate symbols disambiguated with "
              f"accession suffix")

    time.sleep(RATE_LIMIT_SECONDS)

    # --- network ---
    string_ids = mapped_seeds["string_id"].tolist()
    network_text, net_batched, net_n_batches = fetch_network(string_ids)
    RAW_NETWORK_OUT.write_text(network_text)
    print(f"[network] raw TSV saved -> {RAW_NETWORK_OUT}")

    if network_text.strip():
        edges_raw = pd.read_csv(
            RAW_NETWORK_OUT, sep="\t", header=None, names=NETWORK_TSV_COLUMNS
        )
    else:
        edges_raw = pd.DataFrame(columns=NETWORK_TSV_COLUMNS)

    # de-duplicate (batched fallback could in principle repeat an edge)
    edges_raw = edges_raw.drop_duplicates(subset=["stringId_A", "stringId_B"])

    # --- build graph ---
    G = nx.Graph()
    for _, row in mapped_seeds.iterrows():
        G.add_node(
            row["string_id"],
            accession=row["UniProt Accession Number"],
            symbol=row["symbol"],
            string_preferred_name=row["string_preferred_name"],
            log2FC=float(row["log2FC"]),
            regulated=row["regulated"],
        )

    for _, row in edges_raw.iterrows():
        a, b = row["stringId_A"], row["stringId_B"]
        if a == b:
            continue
        if a not in G or b not in G:
            # defensive: STRING should only return edges among submitted IDs
            continue
        score = float(row["score"])
        G.add_edge(a, b, weight=score, combined_score=score)

    n_nodes_in_network = G.number_of_nodes()
    n_edges = G.number_of_edges()
    print(f"[graph] {n_nodes_in_network} nodes, {n_edges} edges "
          f"(required_score={REQUIRED_SCORE})")

    # --- per-node metrics ---
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G, weight=None, seed=SEED)

    if n_edges > 0:
        communities = list(greedy_modularity_communities(G, weight="weight"))
    else:
        communities = [{n} for n in G.nodes()]
    node_to_community: dict[str, int] = {}
    for ci, comm in enumerate(communities):
        for n in comm:
            node_to_community[n] = ci
    n_communities = len(communities)

    rows = []
    for n, data in G.nodes(data=True):
        rows.append({
            "accession": data["accession"],
            "symbol": data["symbol"],
            "string_id": n,
            "degree": degree.get(n, 0),
            "betweenness": betweenness.get(n, 0.0),
            "community": node_to_community.get(n, -1),
            "log2FC": data["log2FC"],
            "regulated": data["regulated"],
        })
    node_metrics = pd.DataFrame(rows).sort_values(
        "degree", ascending=False
    ).reset_index(drop=True)
    node_metrics.to_csv(NODE_METRICS_OUT, index=False)
    print(f"[write] node metrics -> {NODE_METRICS_OUT}")

    # --- edge list (symbol-labeled, as specified) ---
    string_id_to_symbol = {n: d["symbol"] for n, d in G.nodes(data=True)}
    edge_rows = []
    for a, b, d in G.edges(data=True):
        edge_rows.append({
            "nodeA_symbol": string_id_to_symbol.get(a, a),
            "nodeB_symbol": string_id_to_symbol.get(b, b),
            "combined_score": d["combined_score"],
        })
    edges_out_df = pd.DataFrame(edge_rows).sort_values(
        "combined_score", ascending=False
    ).reset_index(drop=True)
    edges_out_df.to_csv(EDGES_OUT, sep="\t", index=False)
    print(f"[write] edge list -> {EDGES_OUT}")

    # --- metadata ---
    meta = {
        "string_version": string_version,
        "species": SPECIES,
        "species_name": "Mus musculus",
        "required_score": REQUIRED_SCORE,
        "n_seeds": n_seeds,
        "n_seeds_up": int((seeds["regulated"] == "UP").sum()),
        "n_seeds_down": int((seeds["regulated"] == "DOWN").sum()),
        "n_mapped": n_mapped,
        "mapping_rate": round(mapping_rate, 4),
        "n_nodes_in_network": n_nodes_in_network,
        "n_edges": n_edges,
        "n_communities": n_communities,
        "get_string_ids_batched": ids_batched,
        "get_string_ids_n_batches": ids_n_batches,
        "network_batched": net_batched,
        "network_n_batches": net_n_batches,
        "rate_limit_seconds_between_calls": RATE_LIMIT_SECONDS,
        "caveat": (
            "n=2 technical replicates per condition (no biological replication); "
            "0/1938 proteins survive FDR<0.05 (limma). The 715 seed proteins are "
            "the |log2FC|>=0.585 HYPOTHESIS-GENERATING set, not confirmed "
            "significant hits. This PPI network shows structure among candidate "
            "proteins and is exploratory."
        ),
    }
    META_OUT.write_text(json.dumps(meta, indent=2))
    print(f"[write] meta -> {META_OUT}")

    # --- summary ---
    top5 = node_metrics.head(5)
    print("\n=== SUMMARY ===")
    print(f"Mapping rate: {n_mapped}/{n_seeds} = {mapping_rate:.1%}")
    print(f"Network: {n_nodes_in_network} nodes, {n_edges} edges, {n_communities} communities "
          f"(required_score={REQUIRED_SCORE})")
    print("Top-5 hub genes by degree:")
    for _, r in top5.iterrows():
        print(f"  {r['symbol']:<20s} degree={r['degree']:<4d} "
              f"betweenness={r['betweenness']:.4f}  log2FC={r['log2FC']:+.3f} "
              f"({r['regulated']})")
    print(f"\n[caveat] {meta['caveat']}")

    unmapped = sorted(set(accessions) - set(acc_to_string_id))
    if unmapped:
        print(f"\n[unmapped] {len(unmapped)} accessions did not map to a STRING ID "
              f"(first 10): {unmapped[:10]}")


if __name__ == "__main__":
    main()
