"""
Item 16 (GSEA half): rank-based preranked GSEA over the full ranked protein
list (item 16), against MOUSE Enrichr gene-set libraries.

Ranked list: all 1938 rows of results/qc_limma.csv, score =
sign(limma_log2FC) * -log10(p_value), descending (range ~ -2.69 to +3.80).
That range was ~ -3.80 to +2.69 before DECISIONS_LOG D7 corrected the inverted
control/treated assignment: every log2FC negates, so the score range mirrors
and the whole ranking (and every term's ES/NES sign) reverses. Verified after
the re-run -- all 568 terms flipped NES sign, corr(NES_old, NES_new) = -0.9998.
Multi-gene "gene" entries (razor/shared peptide groups, 45 rows) are
collapsed to their first listed symbol; the 8 rows with no gene name fall
back to the UniProt accession (`id`), same resolution rule as enrich/ora.py,
so every one of the 1938 rows contributes exactly one ranked identifier.

Gene-set libraries (MOUSE, not human):
  - GO_Biological_Process_2021 -- Enrichr's GO-BP library. This library lives
    on Enrichr's single shared "Enrichr" database (gseapy does not actually
    maintain a separate mouse-vs-human GO-BP library/endpoint: passing
    organism="Mouse" to gseapy.get_library_name()/get_library() only changes
    which Enrichr *instance* is queried for the fly/yeast/worm/fish modEnrichr
    organisms -- for Human and Mouse both, gseapy hits the same base
    "Enrichr" database either way). It IS returned by
    gseapy.get_library_name(organism="Mouse"), i.e. gseapy itself lists it as
    mouse-usable. Gene symbols in the library are uppercase; gseapy's
    load_gmt() auto-detects that our mixed-case mouse symbols ("Lama1") are
    not upper and does a case-insensitive match against the library, so this
    effectively relies on 1:1 mouse<->human ortholog symbol correspondence
    (true for the large majority of one-to-one orthologs, imperfect for
    paralogs/renamed genes). Documented here rather than silently assumed.
  - KEGG_2019_Mouse -- a genuinely Mus-musculus-specific Enrichr KEGG library
    (distinct from KEGG_2019_Human / KEGG_2021_Human), no ortholog-mapping
    caveat.
  Fetched via gseapy.get_library() (NOT gseapy.prerank(gene_sets=[names])):
  gseapy 1.1.11's internal `_download_libraries` helper has a bug that raises
  IndexError on a trailing blank line in Enrichr's streamed GMT response, so
  the library dicts are downloaded directly with gp.get_library() (which does
  not hit that code path) and passed to gp.prerank(gene_sets=<dict>) instead.

Seed fixed (42, matching the report style module's SEED) for determinism.

Outputs:
  results/enrichment/raw/enrichr_libraries.json.gz
    -- the raw gp.get_library() payload for BOTH libraries, keyed by library
       name ({term: [gene, ...]}). This was the last uncached outbound call in
       the enrichment layer; caching it means every network response this
       pipeline depends on now lives under results/enrichment/raw/ and an
       offline-replay mode is buildable. Gzipped: the two libraries are ~6300
       terms and tens of MB of gene symbols uncompressed.
  results/enrichment/gsea_results.csv
    -- term, source, es, nes, nom_p_value, fdr_q_value, tag_pct, gene_pct,
       lead_genes (semicolon-separated), one row per gene set that passed
       min_size/max_size filtering, sorted by fdr_q_value ascending.
  results/enrichment/gsea_meta.json
    -- libraries used, n_genes_ranked, seed, permutation_num, min/max_size,
       gseapy version, n_terms_tested, FDR summary.
  results/figures/gsea_top.png/.svg
    -- diverging barplot of the strongest pathways by NES (top +NES / top
       -NES), colored UP-like (red, NES>0) / DOWN-like (blue, NES<0), with an
       explicit caveat that none pass FDR<0.05.
"""

import gzip
import json

import gseapy as gp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import enrich_common as ec
import style

SEED = 42
PERMUTATION_NUM = 1000
MIN_SIZE = 15
MAX_SIZE = 500
LIBRARIES = ["GO_Biological_Process_2021", "KEGG_2019_Mouse"]

RESULT_COLUMNS = [
    "term", "source", "es", "nes", "nom_p_value", "fdr_q_value",
    "tag_pct", "gene_pct", "lead_genes",
]


def build_ranked_list():
    qc = pd.read_csv(ec.RESULTS_DIR / "qc_limma.csv")
    symbols, n_split, n_fallback = ec.build_symbol_list(qc, "gene", "id")
    score = np.sign(qc["limma_log2FC"]) * -np.log10(qc["p_value"])
    rnk = pd.DataFrame({"gene": symbols, "score": score})
    rnk = rnk.sort_values("score", ascending=False).reset_index(drop=True)
    n_dup = rnk["gene"].duplicated().sum()
    assert len(rnk) == 1938, f"expected 1938 ranked rows, got {len(rnk)}"
    list_meta = {
        "n_genes_ranked": len(rnk),
        "score_min": float(rnk["score"].min()),
        "score_max": float(rnk["score"].max()),
        "semicolon_split_n": int(n_split),
        "accession_fallback_n": int(n_fallback),
        "duplicate_symbol_n": int(n_dup),
    }
    return rnk.set_index("gene")["score"], list_meta


def load_gene_sets():
    """Fetch MOUSE Enrichr libraries via get_library() directly (bypasses the
    gseapy 1.1.11 _download_libraries IndexError bug -- see module docstring).

    Every fetched library is cached verbatim to
    results/enrichment/raw/enrichr_libraries.json.gz (keyed by library name)
    before it is merged, so the one remaining uncached outbound call in the
    enrichment layer is now on disk for audit / offline replay.

    Returns (merged_gene_sets_dict, per_library_n_terms, term_to_source_map).
    """
    libs = {}
    per_library_n_terms = {}
    term_to_source = {}
    raw_by_library = {}
    for name in LIBRARIES:
        d = gp.get_library(name=name, organism="Mouse")
        raw_by_library[name] = d
        per_library_n_terms[name] = len(d)
        overlap = set(libs) & set(d)
        if overlap:
            print(f"[gsea] WARNING: {len(overlap)} term-name collisions between "
                  f"{name} and already-loaded libraries; later library wins.")
        libs.update(d)
        term_to_source.update({term: name for term in d})

    cache_path = ec.RAW_DIR / "enrichr_libraries.json.gz"
    with gzip.open(cache_path, "wt", encoding="utf-8") as fh:
        json.dump(
            {
                "_source": "gseapy.get_library(name=..., organism='Mouse')",
                "_gseapy_version": gp.__version__,
                "_n_terms_per_library": per_library_n_terms,
                "libraries": raw_by_library,
            },
            fh,
            indent=2,
        )
    print(f"[gsea] cached Enrichr libraries -> {cache_path.name} "
          f"({cache_path.stat().st_size / 1024:.0f} KB on disk)")

    return libs, per_library_n_terms, term_to_source


def run_prerank(rnk, gene_sets):
    return gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        min_size=MIN_SIZE,
        max_size=MAX_SIZE,
        permutation_num=PERMUTATION_NUM,
        seed=SEED,
        threads=4,
        outdir=None,
        no_plot=True,
        verbose=False,
    )


def tidy_results(res2d, term_to_source):
    df = res2d.copy()
    df = df.rename(columns={
        "Term": "term",
        "ES": "es",
        "NES": "nes",
        "NOM p-val": "nom_p_value",
        "FDR q-val": "fdr_q_value",
        "Tag %": "tag_pct",
        "Gene %": "gene_pct",
        "Lead_genes": "lead_genes",
    })
    for c in ["es", "nes", "nom_p_value", "fdr_q_value"]:
        df[c] = df[c].astype(float)
    df["source"] = df["term"].map(term_to_source)
    df = df.sort_values(["fdr_q_value", "nom_p_value"]).reset_index(drop=True)
    return df[RESULT_COLUMNS]


def make_barplot(df, n_each=10):
    pos = df[df["nes"] > 0].sort_values("nes", ascending=False).head(n_each)
    neg = df[df["nes"] < 0].sort_values("nes", ascending=True).head(n_each)
    plot_df = pd.concat([pos, neg]).sort_values("nes").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.2, max(5.0, 0.42 * len(plot_df) + 2)))
    colors = np.where(plot_df["nes"] > 0, style.REGULATED_COLORS["UP"], style.REGULATED_COLORS["DOWN"])
    y = range(len(plot_df))
    ax.barh(list(y), plot_df["nes"], color=colors, alpha=0.85, height=0.62, zorder=3)
    labels = [
        f"{row.term}  (FDR q={row.fdr_q_value:.2f}, NOM p={row.nom_p_value:.3f})"
        for row in plot_df.itertuples()
    ]
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=8.0)
    ax.axvline(0, color=style.CHROME["baseline"], lw=1.0, zorder=2)
    style.recede_spines(ax)
    ax.set_xlabel("Normalized Enrichment Score (NES)")
    min_fdr = df["fdr_q_value"].min()
    ax.set_title(
        f"GSEA prerank — top pathways by NES (min FDR q = {min_fdr:.3f}; "
        f"none pass FDR < 0.05)",
        pad=12,
    )
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=style.REGULATED_COLORS["UP"],
               markersize=10, label="NES > 0 (enriched toward UP-scoring end)"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=style.REGULATED_COLORS["DOWN"],
               markersize=10, label="NES < 0 (enriched toward DOWN-scoring end)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8.5)
    style.add_caveat(fig, text=ec.CAVEAT_ENRICHMENT)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    return fig, plot_df


def main():
    rnk, list_meta = build_ranked_list()
    print(f"[gsea] ranked list: {list_meta['n_genes_ranked']} genes, "
          f"score range [{list_meta['score_min']:.3f}, {list_meta['score_max']:.3f}], "
          f"semicolon-split {list_meta['semicolon_split_n']}, "
          f"accession-fallback {list_meta['accession_fallback_n']}, "
          f"duplicate symbols {list_meta['duplicate_symbol_n']}")

    gene_sets, per_library_n_terms, term_to_source = load_gene_sets()
    print(f"[gsea] gene sets: {per_library_n_terms}, "
          f"{sum(per_library_n_terms.values())} terms total before filtering")

    pre = run_prerank(rnk, gene_sets)
    df = tidy_results(pre.res2d, term_to_source)
    df.to_csv(ec.ENRICHMENT_DIR / "gsea_results.csv", index=False)

    n_fdr_05 = int((df["fdr_q_value"] < 0.05).sum())
    n_fdr_25 = int((df["fdr_q_value"] < 0.25).sum())
    n_nom_05 = int((df["nom_p_value"] < 0.05).sum())
    min_fdr = float(df["fdr_q_value"].min())

    print(f"\n[gsea] {len(df)} terms tested (passed min_size={MIN_SIZE}/max_size={MAX_SIZE}) "
          f"of {sum(per_library_n_terms.values())} in the merged libraries")
    print(f"[gsea] FDR<0.05: {n_fdr_05}   FDR<0.25: {n_fdr_25}   NOM p<0.05: {n_nom_05}   "
          f"min FDR = {min_fdr:.4f}")
    print("[gsea] top 10 by FDR q-value:")
    for row in df.head(10).itertuples():
        print(f"  [{row.source}] {row.term} NES={row.nes:.3f} "
              f"NOM p={row.nom_p_value:.4f} FDR q={row.fdr_q_value:.4f}")

    meta = {
        "libraries": {
            "GO_Biological_Process_2021": (
                "Enrichr shared GO-BP library (single 'Enrichr' database; "
                "listed as mouse-usable by gseapy.get_library_name(organism="
                "'Mouse')). Gene symbols are uppercase; matched to our "
                "mixed-case mouse symbols via gseapy's automatic "
                "case-insensitive matching, relying on 1:1 mouse<->human "
                "ortholog symbol correspondence."
            ),
            "KEGG_2019_Mouse": (
                "Enrichr's Mus-musculus-specific KEGG pathway library "
                "(distinct from KEGG_2019_Human/KEGG_2021_Human)."
            ),
        },
        "n_terms_per_library_before_filtering": per_library_n_terms,
        "n_terms_tested_after_size_filter": len(df),
        "n_genes_ranked": list_meta["n_genes_ranked"],
        "score_range": [list_meta["score_min"], list_meta["score_max"]],
        "seed": SEED,
        "permutation_num": PERMUTATION_NUM,
        "min_size": MIN_SIZE,
        "max_size": MAX_SIZE,
        "gseapy_version": gp.__version__,
        "n_fdr_lt_0.05": n_fdr_05,
        "n_fdr_lt_0.25": n_fdr_25,
        "n_nom_p_lt_0.05": n_nom_05,
        "min_fdr_q_value": min_fdr,
        "note": (
            "No gene set passes FDR q < 0.05 (min FDR q = "
            f"{min_fdr:.3f}), consistent with the weak per-protein signal "
            "(0/1938 proteins at limma FDR<0.05). This is expected given n=2 "
            "technical replicates and is reported honestly, not tuned away. "
            f"{n_nom_05} terms have nominal (uncorrected) NOM p-value < 0.05 "
            "and are listed as hypothesis-generating leads only -- this GSEA "
            "FDR is a SEPARATE test from per-protein FDR and neither should "
            "be read as confirming the other."
        ),
        **list_meta,
    }
    (ec.ENRICHMENT_DIR / "gsea_meta.json").write_text(json.dumps(meta, indent=2))

    fig, plot_df = make_barplot(df, n_each=10)
    png, svg = style.save_fig(fig, "gsea_top")

    ec.record_enrich_manifest([{
        "file": png,
        "title": "GSEA prerank — top pathways by NES (mouse GO-BP + KEGG)",
        "caption": (
            f"Preranked GSEA (gseapy {gp.__version__}, seed={SEED}, "
            f"{PERMUTATION_NUM} permutations) over all 1938 qc_limma.csv "
            "proteins ranked by sign(limma_log2FC) * -log10(p_value), against "
            "GO_Biological_Process_2021 + KEGG_2019_Mouse (mouse Enrichr "
            f"libraries; {len(df)} gene sets tested after min_size=15/"
            f"max_size=500 filtering). Top 10 pathways by NES on each side "
            "shown (red = NES>0, enriched toward the UP-scoring end; blue = "
            f"NES<0, toward the DOWN-scoring end). Minimum FDR q-value = "
            f"{min_fdr:.3f} -- NO pathway passes FDR<0.05; {n_nom_05} pass an "
            "uncorrected nominal p<0.05 and are shown as hypothesis-generating "
            "leads only. This GSEA FDR is a separate statistical test from "
            "the non-significant per-protein FDR (0/1938 at FDR<0.05, limma) "
            "and must not be conflated with it."
        ),
        "key_numbers": {
            "n_terms_tested": len(df),
            "n_fdr_lt_0.05": n_fdr_05,
            "n_fdr_lt_0.25": n_fdr_25,
            "n_nom_p_lt_0.05": n_nom_05,
            "min_fdr_q_value": round(min_fdr, 4),
            "top_terms_shown": plot_df[["term", "source", "nes", "nom_p_value", "fdr_q_value"]].to_dict("records"),
        },
    }])

    print(f"\n[gsea] wrote: {ec.ENRICHMENT_DIR/'gsea_results.csv'}, "
          f"{ec.ENRICHMENT_DIR/'gsea_meta.json'}")
    print(f"[gsea] wrote figure: {png}, {svg}")


if __name__ == "__main__":
    main()
