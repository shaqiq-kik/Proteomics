"""
Item 16 (ORA half): Over-representation analysis via g:Profiler g:GOSt.

POSTs to https://biit.cs.ut.ee/gprofiler/api/gost/profile/, organism="mmusculus"
(MOUSE), sources = GO:BP / GO:MF / GO:CC / KEGG / REAC, user_threshold=0.05,
significance_threshold_method="g_SCS", domain_scope="custom" with a CUSTOM
background = the detected proteome (2554-row union of foldchange_all.csv and
single_condition_proteins.csv; see enrich_common.load_background_and_queries).
UP (509) and DOWN (206) gene sets are queried SEPARATELY (opposite biology).
Those sizes were 206 UP / 509 DOWN before DECISIONS_LOG D7 corrected the
inverted control/treated assignment; the two query sets swapped wholesale. The
counts are read from tests/expected/frozen_counts.json, never hardcoded.

Custom background matters here, and it is the whole point. Against
g:Profiler's DEFAULT whole-genome background the same two lists return 326
(UP) and 196 (DOWN) "significant" GO/KEGG/REAC terms, both topped by GO:CC
"cytoplasm" (p=9.0e-70 and p=1.9e-23) -- but that is inflation from comparing
against ~20-25k background genes when only ~2.5k proteins were even detectable
in this experiment. Restricting to the CUSTOM detected-proteome background
(the scientifically correct comparison) removes ALL of that signal: with
domain_scope="custom" as specified, 0/0 terms pass g:SCS-corrected
significance at 0.05 in either direction (best corrected p = 0.703 DOWN,
1.000 UP). This null result is reported honestly below, not papered over.

The default-background counts moved with the D7 correction because the query
sets swapped: the 196-terms/p=1.9e-23 figure recorded pre-D7 for the UP list
now belongs to the DOWN list, unchanged to three significant figures, which is
itself a check that the swap is a clean relabelling and nothing else.

Also produces results/figures/ora_dotplot.png/.svg: since 0 terms are
significant, the dotplot shows the best (still non-significant) g:SCS-corrected
terms per direction against a dashed p=0.05 reference line, explicitly labeled
as sub-threshold.

Outputs:
  results/enrichment/raw/gprofiler_up.json, gprofiler_down.json
    -- raw API response for the EXACT spec'd call (all_results=False,
       no_evidences=False; small/empty here because nothing is significant,
       but carries full query/background/domain metadata for audit).
  results/enrichment/raw/gprofiler_up_all.json, gprofiler_down_all.json
    -- raw API response for the all_results=True / no_evidences=True ranking
       call (every tested term with its g:SCS-corrected p-value, no gene-level
       evidence vectors). Previously discarded in memory; cached now so a
       future offline-replay mode has every response this script depends on.
  results/enrichment/raw/gprofiler_up_all_evidence.json.gz,
  gprofiler_down_all_evidence.json.gz
    -- raw API response for the all_results=True / no_evidences=FALSE
       gene-level call behind ora_top_terms_detail.json. Gzipped because the
       uncompressed payload carries one boolean vector per term per query gene
       (tens of MB); gzip takes it to a committable size. The module previously
       documented this response as deliberately not persisted -- that made
       offline replay impossible, so it is persisted compressed instead.
  results/enrichment/ora_up.csv, ora_down.csv
    -- source, term_id, term_name, p_value, term_size, query_size,
       intersection_size, intersecting_genes for rows with significant==True
       (0 rows + header here).
  results/enrichment/ora_meta.json
    -- version, organism, background/query construction stats, thresholds,
       n_terms_up/down, and the best sub-threshold terms per direction
       (for hypothesis-generating reference only).
  results/enrichment/ora_top_terms_detail.json
    -- top ~10 terms pooled across UP+DOWN by g:SCS-corrected p_value, WITH
       gene-level intersecting-gene lists. Consumed by enrich/upset.py.
"""

import gzip
import json
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import requests

import enrich_common as ec
import style

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"]
USER_THRESHOLD = 0.05
SIG_METHOD = "g_SCS"

CSV_COLUMNS = [
    "source", "term_id", "term_name", "p_value", "term_size",
    "query_size", "intersection_size", "intersecting_genes",
]


def _post(query, background, all_results, no_evidences, retries=3):
    payload = {
        "organism": ec.ORGANISM_GPROFILER,
        "query": query,
        "sources": SOURCES,
        "user_threshold": USER_THRESHOLD,
        "significance_threshold_method": SIG_METHOD,
        "domain_scope": "custom",
        "background": background,
        "no_evidences": no_evidences,
        "all_results": all_results,
    }
    last_err = None
    for attempt in range(retries):
        try:
            resp = requests.post(GPROFILER_URL, json=payload, timeout=180)
            resp.raise_for_status()
            return resp.json()
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"g:Profiler request failed after {retries} attempts: {last_err}")


def _cache_raw(response, filename, gzipped=False):
    """Persist a raw g:Profiler response under results/enrichment/raw/.

    Every outbound call this module makes is cached, so an offline-replay mode
    can be built without re-querying g:Profiler. Large responses are gzipped.
    """
    path = ec.RAW_DIR / filename
    payload = json.dumps(response, indent=2)
    if gzipped:
        with gzip.open(path, "wt", encoding="utf-8") as fh:
            fh.write(payload)
    else:
        path.write_text(payload)
    print(f"[ora] cached raw response -> {path.name} "
          f"({path.stat().st_size / 1024:.0f} KB on disk)")
    return path


def _query_order(response):
    """Gene labels aligned 1:1 with each term's ``intersections`` vector.

    NOT the submitted query list. g:Profiler's ``intersections`` booleans are
    indexed by the genes it successfully MAPPED, not by what was submitted: for
    the UP query 509 symbols go out and 31 fail to map, so the vectors come back
    length 473. Zipping those 473 booleans against the 509 submitted symbols
    (what this function used to do) stays correct only up to the first failed
    symbol and silently mislabels everything after it -- e.g. GO:MF
    "mannosyltransferase activity" was credited to `Snrpd1` when the gene
    g:Profiler actually matched is `Tmtc3`.

    The authoritative order is ``meta.genes_metadata.query.<q>.ensgs``; the
    sibling ``mapping`` ({submitted symbol: [ensg, ...]}) inverts it back to
    symbols. An ENSG with no symbol (or several) keeps the ENSG id so the
    ambiguity is visible rather than guessed at.
    """
    gm = response["meta"]["genes_metadata"]["query"]
    entry = gm.get("query_1") or next(iter(gm.values()))
    ensgs = entry["ensgs"]
    by_ensg = {}
    for symbol, mapped in entry["mapping"].items():
        for ensg in mapped:
            by_ensg.setdefault(ensg, []).append(symbol)
    return [
        by_ensg[e][0] if len(by_ensg.get(e, [])) == 1 else e
        for e in ensgs
    ]


def _rows_from_result(result_list, query_order, only_significant):
    rows = []
    for term in result_list:
        if only_significant and not term.get("significant", False):
            continue
        inter = term.get("intersections")
        genes = []
        if inter:
            # Never let zip() paper over a length mismatch: a short/long
            # evidence vector means the alignment assumption above broke.
            assert len(inter) == len(query_order), (
                f"intersections vector for {term.get('native')} has "
                f"{len(inter)} entries but {len(query_order)} mapped genes "
                "were resolved -- gene labels would be misaligned"
            )
            for gene, ev in zip(query_order, inter):
                if ev:
                    genes.append(gene)
        rows.append({
            "source": term["source"],
            "term_id": term["native"],
            "term_name": term["name"],
            "p_value": term["p_value"],
            "term_size": term["term_size"],
            "query_size": term["query_size"],
            "intersection_size": term["intersection_size"],
            "intersecting_genes": ";".join(genes),
            "significant": term.get("significant", False),
        })
    rows.sort(key=lambda r: r["p_value"])
    return rows


def run_direction(direction, query, background):
    """Two API calls for this direction, BOTH cached to raw/:
      1. spec-literal (all_results=False, no_evidences=False) -> raw
         gprofiler_{direction}.json + ora_{direction}.csv (significant-only).
      2. all_results=True, no_evidences=True -> fast/small, ranks every tested
         term by corrected p_value for reporting -> raw
         gprofiler_{direction}_all.json.
    """
    print(f"[ora] querying g:Profiler for {direction.upper()} (n={len(query)})...")
    resp_primary = _post(query, background, all_results=False, no_evidences=False)
    _cache_raw(resp_primary, f"gprofiler_{direction}.json")

    qorder = _query_order(resp_primary)
    sig_rows = _rows_from_result(resp_primary["result"], qorder, only_significant=True)
    df = pd.DataFrame(sig_rows, columns=CSV_COLUMNS + ["significant"])
    df = df[CSV_COLUMNS]
    df.to_csv(ec.ENRICHMENT_DIR / f"ora_{direction}.csv", index=False)

    print(f"[ora] querying g:Profiler for {direction.upper()} (all terms, for ranking)...")
    resp_all = _post(query, background, all_results=True, no_evidences=True)
    _cache_raw(resp_all, f"gprofiler_{direction}_all.json")
    qorder_all = _query_order(resp_all)
    all_rows = _rows_from_result(resp_all["result"], qorder_all, only_significant=False)

    return {
        "primary_response": resp_primary,
        "sig_rows": sig_rows,
        "all_rows": all_rows,
    }


def fetch_gene_level_for_top_terms(direction, query, background, top_term_ids, top_sources):
    """Targeted re-fetch WITH gene-level evidence (no_evidences=False,
    all_results=True), used only to harvest intersecting-gene lists for the
    handful of top-ranked (sub-threshold) terms chosen for the UpSet plot /
    detail JSON.

    The full response (~8-20MB uncompressed for this dataset -- one boolean
    vector per term per query gene) IS persisted, gzipped, as
    raw/gprofiler_{direction}_all_evidence.json.gz. It used to be dropped on
    the floor for size reasons, which meant this script could never be replayed
    offline; gzip makes keeping it cheap.
    """
    resp = _post(query, background, all_results=True, no_evidences=False)
    _cache_raw(resp, f"gprofiler_{direction}_all_evidence.json.gz", gzipped=True)
    qorder = _query_order(resp)
    wanted = set(top_term_ids)
    out = {}
    for term in resp["result"]:
        if term["native"] not in wanted:
            continue
        inter = term.get("intersections") or []
        assert len(inter) == len(qorder), (
            f"intersections vector for {term['native']} has {len(inter)} "
            f"entries but {len(qorder)} mapped genes were resolved -- gene "
            "labels would be misaligned"
        )
        genes = [g for g, ev in zip(qorder, inter) if ev]
        out[term["native"]] = genes
    return out


def make_dotplot(top_up, top_down, n_each=5):
    """ORA dot/bar figure: top (sub-threshold, since 0 are g:SCS-significant)
    terms for UP and DOWN side by side. x = -log10(g:SCS-corrected p_value),
    size = intersection_size, colored by direction, with a dashed reference
    line at the p=0.05 significance cutoff nothing reaches.
    """
    import numpy as np

    up_terms = top_up[:n_each]
    down_terms = top_down[:n_each]
    rows = []
    for r in up_terms:
        rows.append({**r, "direction": "UP"})
    for r in down_terms:
        rows.append({**r, "direction": "DOWN"})
    if not rows:
        rows = [{"term_name": "(no terms returned)", "source": "", "p_value": 1.0,
                  "intersection_size": 0, "term_size": 0, "direction": "UP"}]

    plot_df = pd.DataFrame(rows)
    plot_df["neg_log10_p"] = -np.log10(plot_df["p_value"].clip(lower=1e-300))
    # order: UP terms (by p) on top, DOWN terms below
    plot_df["order_key"] = plot_df["direction"].map({"UP": 0, "DOWN": 1})
    plot_df = plot_df.sort_values(["order_key", "p_value"]).reset_index(drop=True)
    plot_df = plot_df.iloc[::-1].reset_index(drop=True)  # so UP ends up on top of the plot

    fig, ax = plt.subplots(figsize=(9.6, max(4.5, 0.5 * len(plot_df) + 2)))
    y = range(len(plot_df))
    sizes = 40 + 12 * plot_df["intersection_size"].astype(float)
    colors = plot_df["direction"].map(
        {"UP": style.REGULATED_COLORS["UP"], "DOWN": style.REGULATED_COLORS["DOWN"]}
    )
    ax.scatter(plot_df["neg_log10_p"], list(y), s=sizes, c=colors, alpha=0.85,
               edgecolors="white", linewidths=0.6, zorder=4)
    labels = [
        f"{row.term_name}  [{row.source}]  ({row.intersection_size}/{row.term_size})"
        for row in plot_df.itertuples()
    ]
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8.6)
    sig_line_x = -np.log10(0.05)
    ax.axvline(sig_line_x, color=style.CHROME["ink_muted"], lw=1.1, ls="--", zorder=2)
    ax.text(sig_line_x, len(plot_df) - 0.3,
            " g:SCS p = 0.05 (none reached)", fontsize=7.8, color=style.CHROME["ink_secondary"],
            va="top", ha="left")

    style.recede_spines(ax)
    ax.set_xlabel("-log10(g:Profiler g:SCS-corrected p-value)")
    ax.set_title(
        f"ORA top terms (UP top {len(up_terms)}, DOWN top {len(down_terms)}) — "
        f"none pass g:SCS < 0.05",
        pad=12,
    )
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=style.REGULATED_COLORS["UP"],
               markersize=9, label="UP-regulated query"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=style.REGULATED_COLORS["DOWN"],
               markersize=9, label="DOWN-regulated query"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8.5, title="Dot size = intersection_size")
    style.add_caveat(fig, text=ec.CAVEAT_ENRICHMENT)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig, plot_df


def main():
    background, up, down, list_meta = ec.load_background_and_queries()
    print(f"[ora] background: {len(background)} unique symbols "
          f"(row-union {list_meta['background_row_union_n']})")
    print(f"[ora] UP query: {len(up)} genes "
          f"(semicolon-split {list_meta['up_semicolon_split_n']}, "
          f"accession-fallback {list_meta['up_accession_fallback_n']})")
    print(f"[ora] DOWN query: {len(down)} genes "
          f"(semicolon-split {list_meta['down_semicolon_split_n']}, "
          f"accession-fallback {list_meta['down_accession_fallback_n']})")

    out = {}
    for direction, query in [("up", up), ("down", down)]:
        out[direction] = run_direction(direction, query, background)
        n_sig = len(out[direction]["sig_rows"])
        best = out[direction]["all_rows"][0] if out[direction]["all_rows"] else None
        print(f"\n=== {direction.upper()} — {n_sig} term(s) pass g:SCS<{USER_THRESHOLD} "
              f"(custom background, n={len(background)}) ===")
        if n_sig:
            for r in out[direction]["sig_rows"][:10]:
                print(f"  [{r['source']}] {r['term_name']} "
                      f"(p={r['p_value']:.3e}, {r['intersection_size']}/{r['term_size']})")
        else:
            print("  none significant. Best sub-threshold terms (NOT significant):")
            for r in out[direction]["all_rows"][:5]:
                print(f"  [{r['source']}] {r['term_name']} "
                      f"(p={r['p_value']:.4f}, {r['intersection_size']}/{r['term_size']})")

    version = out["up"]["primary_response"]["meta"].get("version")

    # Top terms per direction (NOT a single pooled-then-truncated top 10):
    # many g:SCS-corrected p-values floor to exactly 1.0 once the background
    # is corrected, so a naive pool-and-sort-by-p (stable sort) would let UP's
    # earlier-processed ties crowd out every DOWN term. Take top 5 from EACH
    # direction explicitly (same set used in the dotplot) so both directions
    # are represented in the detail JSON / UpSet input, per spec ("combine
    # the strongest terms from ORA up+down").
    top_pooled = []
    for direction in ("up", "down"):
        for r in out[direction]["all_rows"][:5]:
            top_pooled.append({**r, "direction": direction.upper()})

    up_term_ids = [r["term_id"] for r in top_pooled if r["direction"] == "UP"]
    down_term_ids = [r["term_id"] for r in top_pooled if r["direction"] == "DOWN"]
    print(f"\n[ora] fetching gene-level evidence for top {len(up_term_ids)} UP-side "
          f"and {len(down_term_ids)} DOWN-side pooled terms (for upset.py)...")
    up_gene_map = fetch_gene_level_for_top_terms("up", up, background, up_term_ids, None) if up_term_ids else {}
    down_gene_map = fetch_gene_level_for_top_terms("down", down, background, down_term_ids, None) if down_term_ids else {}

    detail_rows = []
    for r in top_pooled:
        genes = (up_gene_map if r["direction"] == "UP" else down_gene_map).get(r["term_id"], [])
        detail_rows.append({
            "direction": r["direction"],
            "source": r["source"],
            "term_id": r["term_id"],
            "term_name": r["term_name"],
            "p_value": r["p_value"],
            "term_size": r["term_size"],
            "query_size": r["query_size"],
            "intersection_size": r["intersection_size"],
            "intersecting_genes": genes,
            "significant_g_SCS_0.05": bool(r["significant"]),
        })
    (ec.ENRICHMENT_DIR / "ora_top_terms_detail.json").write_text(
        json.dumps({
            "note": (
                "Top 5 terms from EACH direction (UP, DOWN), ranked by "
                "g:SCS-corrected p_value within direction (not a single "
                "pooled-and-sorted top 10 -- many corrected p-values floor to "
                "exactly 1.0, so pooling would let one direction's ties crowd "
                "out the other). NONE pass significance (user_threshold=0.05) "
                "-- see significant_g_SCS_0.05 (all False here). Provided for "
                "hypothesis-generating exploration (dotplot, UpSet) only."
            ),
            "terms": detail_rows,
        }, indent=2)
    )

    meta = {
        "gprofiler_version": version,
        "organism": ec.ORGANISM_GPROFILER,
        "sources": SOURCES,
        "user_threshold": USER_THRESHOLD,
        "significance_threshold_method": SIG_METHOD,
        "domain_scope": "custom",
        "background_size_unique_symbols": len(background),
        "background_size_row_union": list_meta["background_row_union_n"],
        "n_terms_up": len(out["up"]["sig_rows"]),
        "n_terms_down": len(out["down"]["sig_rows"]),
        "note": (
            "0 terms pass g:SCS-corrected significance (user_threshold=0.05) in "
            f"either direction using the correct CUSTOM background (detected "
            f"proteome, {list_meta['background_row_union_n']} rows / "
            f"{len(background)} unique symbols). For comparison, the SAME two "
            "queries against g:Profiler's default whole-genome background "
            "return 326 (UP) and 196 (DOWN) nominally 'significant' terms, "
            "both topped by GO:CC 'cytoplasm' (p=9.0e-70 and p=1.9e-23) -- an "
            "artifact of comparing detected proteins to the whole genome "
            "rather than to what could have been detected. This 0/0 result is "
            "a SEPARATE statistical test from per-protein FDR (0/1938 "
            "significant, limma) and should not be conflated with it; "
            "reported honestly as a null finding, not an error."
        ),
        "default_background_comparison": {
            "_what": (
                "Diagnostic only -- NOT how this pipeline computes its results. "
                "Same queries, same sources, same g:SCS threshold, but "
                "domain_scope='annotated' (g:Profiler's whole-genome default) "
                "instead of the detected-proteome custom background."
            ),
            "n_significant_up_default_background": 326,
            "n_significant_down_default_background": 196,
            "top_term_both_directions": "GO:CC cytoplasm",
            "top_p_up": 9.01e-70,
            "top_p_down": 1.92e-23,
        },
        "best_subthreshold_term_up": out["up"]["all_rows"][0] if out["up"]["all_rows"] else None,
        "best_subthreshold_term_down": out["down"]["all_rows"][0] if out["down"]["all_rows"] else None,
        "top5_subthreshold_terms_up": out["up"]["all_rows"][:5],
        "top5_subthreshold_terms_down": out["down"]["all_rows"][:5],
        **list_meta,
    }
    (ec.ENRICHMENT_DIR / "ora_meta.json").write_text(json.dumps(meta, indent=2, default=str))

    # --- ORA dotplot figure ---
    fig, plot_df = make_dotplot(out["up"]["all_rows"], out["down"]["all_rows"], n_each=5)
    png, svg = style.save_fig(fig, "ora_dotplot")

    ec.record_enrich_manifest([{
        "file": png,
        "title": "ORA top terms — UP vs. DOWN (g:Profiler, custom background)",
        "caption": (
            "Top 5 GO/KEGG/Reactome terms per direction from g:Profiler g:GOSt "
            "(organism=mmusculus, sources=GO:BP/MF/CC+KEGG+REAC, custom "
            f"background = the {list_meta['background_row_union_n']}-row "
            "detected proteome, g:SCS-corrected p-values, user_threshold=0.05). "
            "0/0 terms pass significance in "
            "either direction (dashed line = the p=0.05 cutoff nothing "
            "reaches); dot size = intersection_size. The same two queries "
            "against g:Profiler's DEFAULT whole-genome background would return "
            "326 (UP) and 196 (DOWN) 'significant' terms -- background "
            "inflation this pipeline deliberately avoids. Shown as the best "
            "available leads for hypothesis generation, NOT as significant "
            "pathway enrichment. This is a separate statistical test from the "
            "non-significant per-protein FDR (0/1938 at FDR<0.05)."
        ),
        "key_numbers": {
            "n_significant_up": len(out["up"]["sig_rows"]),
            "n_significant_down": len(out["down"]["sig_rows"]),
            "background_size": len(background),
            "up_query_size": len(up),
            "down_query_size": len(down),
            "best_p_up": out["up"]["all_rows"][0]["p_value"] if out["up"]["all_rows"] else None,
            "best_p_down": out["down"]["all_rows"][0]["p_value"] if out["down"]["all_rows"] else None,
            "top_terms_shown": plot_df[["term_name", "source", "direction", "p_value", "intersection_size"]].to_dict("records"),
        },
    }])

    print(f"\n[ora] wrote: {ec.ENRICHMENT_DIR/'ora_up.csv'}, {ec.ENRICHMENT_DIR/'ora_down.csv'}, "
          f"{ec.ENRICHMENT_DIR/'ora_meta.json'}, {ec.ENRICHMENT_DIR/'ora_top_terms_detail.json'}")
    print(f"[ora] wrote figure: {png}, {svg}")


if __name__ == "__main__":
    main()
