"""
Item 18: UpSet plot -- gene-overlap structure across the top ORA pathways.

Design choice (documented, since none of these ORA terms are individually
significant -- see enrich/ora.py): with the CORRECT custom background, 0
GO/KEGG/REAC terms pass g:SCS-corrected significance in either direction
(min corrected p ~= 0.68-1.0), so there is no "top 6 SIGNIFICANT pathways" to
draw. Per the task's explicit fallback ("if a clean UpSet isn't meaningful
with the data, build the most defensible set-intersection view and document
the choice"), this plot uses the top-ranked (by g:SCS-corrected p_value)
terms per direction from results/enrichment/ora_top_terms_detail.json
(already carrying gene-level membership).

The spec's "~6" was tried literally first (top 3 UP + top 3 DOWN) but that
cutoff has ZERO gene overlap of any kind -- neither within a direction nor
across -- which makes an UpSet plot pointless (six disconnected bars). Widened
to top 5 UP + top 5 DOWN (10 terms, still the same small "best available
leads" pool already fetched by ora.py, nothing new queried) so the real
within-direction overlaps are visible: the three DOWN-side Reactome mRNA terms
share spliceosomal Sm-like genes -- "mRNA decay by 5' to 3' exoribonuclease"
and "mRNA Splicing - Major Pathway" share Lsm2 and Lsm6, and "mRNA Splicing -
Minor Pathway" and "Major Pathway" share Snrpg. This is the most defensible
set-intersection view available given the data: still hypothesis-generating
leads, not confirmed pathways, but at least shows genuine (if modest)
intersection structure instead of an empty grid.

Two corrections are folded into the current output. (1) DECISIONS_LOG D7
inverted the control/treated assignment, so the UP and DOWN query sets swapped
wholesale -- the Reactome mRNA cluster that used to sit on the UP side is now
on the DOWN side. (2) ora.py used to align g:Profiler's ``intersections``
booleans against the SUBMITTED query list, but that vector is indexed by the
genes g:Profiler successfully MAPPED (509 submitted -> 473 mapped for UP), so
every gene label after the first unmapped symbol was wrong. The gene names
below are the fixed ones; they are also the ones that make biological sense
(e.g. "succinate dehydrogenase activity" -> Sdha, not the old Ybx3).

A structural note baked into the result: UP and DOWN are, by construction,
mutually exclusive protein categories (a protein cannot be both up- and
down-regulated in this contrast), so no gene can ever appear in both a
UP-side pathway and a DOWN-side pathway here -- the UpSet necessarily shows
two disjoint intersection clusters (overlap only WITHIN each direction), not
cross-direction overlap. That is expected, not a display bug, and is called
out in the caption. Real within-direction overlap does show up, e.g. Lsm2 and
Lsm6 are shared between two of the DOWN-side Reactome mRNA terms.

Per-pathway total bars are colored by direction (red = UP-side term, blue =
DOWN-side term, matching the report palette) via UpSet.style_categories().

Output: results/figures/upset.png / .svg
"""

import json
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_contents

import enrich_common as ec
import style

N_PER_DIRECTION = 5  # 5 UP + 5 DOWN; see module docstring for why "~6" (3+3) was widened


def load_top_terms():
    detail = json.loads((ec.ENRICHMENT_DIR / "ora_top_terms_detail.json").read_text())
    terms = detail["terms"]
    up_terms = [t for t in terms if t["direction"] == "UP"][:N_PER_DIRECTION]
    down_terms = [t for t in terms if t["direction"] == "DOWN"][:N_PER_DIRECTION]
    return up_terms, down_terms


def label_for(term, width=42):
    short = f"{term['term_name']} [{term['source']}]"
    wrapped = "\n".join(textwrap.wrap(short, width=width))
    return wrapped


def main():
    up_terms, down_terms = load_top_terms()
    chosen = up_terms + down_terms
    if not chosen:
        raise RuntimeError("No ORA top terms available to build an UpSet plot from.")

    contents = {}
    label_direction = {}
    for t in chosen:
        label = label_for(t)
        genes = t["intersecting_genes"]
        if not genes:
            continue
        contents[label] = set(genes)
        label_direction[label] = t["direction"]

    print(f"[upset] using {len(up_terms)} UP-side + {len(down_terms)} DOWN-side terms:")
    for t in chosen:
        print(f"  [{t['direction']}] [{t['source']}] {t['term_name']} "
              f"-- {len(t['intersecting_genes'])} genes: {t['intersecting_genes']}")

    all_genes_in_sets = sorted(set().union(*contents.values())) if contents else []
    up_gene_set = set().union(*[set(t["intersecting_genes"]) for t in up_terms]) if up_terms else set()
    down_gene_set = set().union(*[set(t["intersecting_genes"]) for t in down_terms]) if down_terms else set()
    cross_overlap = up_gene_set & down_gene_set
    print(f"[upset] {len(all_genes_in_sets)} unique genes across the {len(contents)} "
          f"non-empty sets; UP/DOWN cross-overlap = {len(cross_overlap)} "
          "(expected 0 -- UP and DOWN are mutually exclusive by construction)")

    data = from_contents(contents)

    fig = plt.figure(figsize=(11.5, 7.2))
    upset = UpSet(
        data,
        subset_size="count",
        sort_by="cardinality",
        sort_categories_by="input",
        show_counts=True,
        facecolor=style.CHROME["ink_secondary"],
        element_size=42,
        min_subset_size=1,
    )
    up_labels = [lbl for lbl, d in label_direction.items() if d == "UP"]
    down_labels = [lbl for lbl, d in label_direction.items() if d == "DOWN"]
    if up_labels:
        upset.style_categories(up_labels, bar_facecolor=style.REGULATED_COLORS["UP"])
    if down_labels:
        upset.style_categories(down_labels, bar_facecolor=style.REGULATED_COLORS["DOWN"])

    axes = upset.plot(fig=fig)
    fig.suptitle(
        f"Gene overlap across top ORA leads ({len(up_terms)} UP-side + "
        f"{len(down_terms)} DOWN-side terms; NONE pass g:SCS < 0.05)",
        fontsize=12.5, fontweight="semibold", y=0.99,
    )
    axes["totals"].set_xlabel("genes in term")
    style.add_caveat(fig, text=(
        ec.CAVEAT_ENRICHMENT
        + " UP-side and DOWN-side terms never overlap here by construction "
        "(a protein cannot be both UP- and DOWN-regulated); shown are the "
        "best available (non-significant) ORA leads, not confirmed pathways."
    ), y=0.005)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))

    png, svg = style.save_fig(fig, "upset", tight=False)

    key_numbers = {
        "n_sets_shown": len(contents),
        "n_up_terms": len(up_terms),
        "n_down_terms": len(down_terms),
        "n_unique_genes": len(all_genes_in_sets),
        "up_down_cross_overlap_n": len(cross_overlap),
        "terms": [
            {"direction": t["direction"], "source": t["source"], "term_name": t["term_name"],
             "n_genes": len(t["intersecting_genes"]), "genes": t["intersecting_genes"],
             "p_value_g_SCS": t["p_value"]}
            for t in chosen
        ],
    }
    ec.record_enrich_manifest([{
        "file": png,
        "title": "UpSet plot — gene overlap across top ORA leads (UP + DOWN)",
        "caption": (
            f"UpSet plot of gene-level overlap across the top {len(up_terms)} "
            f"UP-side + top {len(down_terms)} DOWN-side g:Profiler ORA terms "
            "(by g:SCS-corrected p-value; NONE pass significance at 0.05 -- "
            "see enrich/ora.py). Chosen as the most defensible "
            "set-intersection view available since there are no significant "
            "pathways to rank by: these are the best (still non-significant) "
            "leads per direction, widened from a literal top-3+3 (which had "
            "zero overlap of any kind) so the one genuine within-direction "
            "overlap is visible. "
            "Per-pathway total bars are colored by direction (red = UP-side "
            "term, blue = DOWN-side term). UP-side and DOWN-side terms never "
            "share genes here by construction (UP/DOWN are mutually "
            "exclusive protein calls), so the plot shows two disjoint "
            "intersection clusters, not a display error; real overlap WITHIN "
            "a direction is visible (e.g. the DOWN-side Reactome mRNA terms "
            "share the Sm-like genes Lsm2/Lsm6, and Snrpg between the minor "
            "and major splicing pathways). Hypothesis-generating only."
        ),
        "key_numbers": key_numbers,
    }])

    print(f"\n[upset] wrote figure: {png}, {svg}")


if __name__ == "__main__":
    main()
