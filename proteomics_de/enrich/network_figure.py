"""
Static PPI network figure for the STRING network built by `enrich/string_ppi.py`.

Cytoscape/py4cytoscape is not available headless in this environment, so the
deliverable is a static networkx figure using the report's shared visual
system (`proteomics_de/viz/style.py`).

For legibility, only the LARGEST CONNECTED COMPONENT of the 694-node STRING
network is drawn (the remaining nodes are isolated/small satellite components
that would clutter the layout without adding readable structure); the number
of nodes/edges omitted is reported in the title, on the figure, and in the
manifest entry.

Encoding:
  - node color  = log2FC on the shared diverging blue(down)/grey/red(up) scale
  - node size   = proportional to STRING degree (hub proteins draw bigger)
  - node label  = top ~15 hub genes by degree only (avoid label clutter)
  - layout      = deterministic spring_layout(seed=42)

Same n=2 / hypothesis-generating caveat as every other report figure is
stamped on the figure via the shared `add_caveat` helper.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

import sys

ENRICH_DIR = Path(__file__).resolve().parent
PROTEOMICS_DE_DIR = ENRICH_DIR.parent
sys.path.insert(0, str(PROTEOMICS_DE_DIR))

from viz.style import (  # noqa: E402
    CHROME,
    DIVERGING_CMAP,
    FIGURES_DIR,
    RESULTS_DIR,
    CAVEAT_TEXT,
    add_caveat,
    save_fig,
)

ENRICHMENT_DIR = RESULTS_DIR / "enrichment"
NODE_METRICS_CSV = ENRICHMENT_DIR / "string_node_metrics.csv"
EDGES_TSV = ENRICHMENT_DIR / "string_edges.tsv"
META_JSON = ENRICHMENT_DIR / "string_meta.json"

MANIFEST_NETWORK_PATH = FIGURES_DIR / "figures_manifest_network.json"

SEED = 42
N_HUB_LABELS = 15
STEM = "ppi_network"


def load_graph() -> nx.Graph:
    node_metrics = pd.read_csv(NODE_METRICS_CSV)
    edges = pd.read_csv(EDGES_TSV, sep="\t")

    G = nx.Graph()
    for _, row in node_metrics.iterrows():
        G.add_node(
            row["symbol"],
            accession=row["accession"],
            degree=int(row["degree"]),
            betweenness=float(row["betweenness"]),
            community=int(row["community"]),
            log2FC=float(row["log2FC"]),
            regulated=row["regulated"],
        )
    for _, row in edges.iterrows():
        a, b = row["nodeA_symbol"], row["nodeB_symbol"]
        if a in G and b in G:
            G.add_edge(a, b, combined_score=float(row["combined_score"]))
    return G


def main():
    meta = json.loads(META_JSON.read_text())
    G = load_graph()
    n_nodes_full = G.number_of_nodes()
    n_edges_full = G.number_of_edges()

    components = sorted(nx.connected_components(G), key=len, reverse=True)
    lcc_nodes = components[0]
    LCC = G.subgraph(lcc_nodes).copy()

    n_nodes_plotted = LCC.number_of_nodes()
    n_edges_plotted = LCC.number_of_edges()
    n_nodes_omitted = n_nodes_full - n_nodes_plotted
    n_edges_omitted = n_edges_full - n_edges_plotted
    n_components_omitted = len(components) - 1

    print(f"[network_figure] full network: {n_nodes_full} nodes, {n_edges_full} edges")
    print(f"[network_figure] largest connected component: {n_nodes_plotted} nodes, "
          f"{n_edges_plotted} edges")
    print(f"[network_figure] omitted: {n_nodes_omitted} nodes across "
          f"{n_components_omitted} smaller components ({n_edges_omitted} edges)")

    # --- deterministic layout ---
    pos = nx.spring_layout(
        LCC, seed=SEED, weight="combined_score", k=1.1 / np.sqrt(n_nodes_plotted),
        iterations=80,
    )

    # --- encodings ---
    log2fc = np.array([LCC.nodes[n]["log2FC"] for n in LCC.nodes()])
    degree = np.array([LCC.nodes[n]["degree"] for n in LCC.nodes()])

    vmax = float(np.max(np.abs(log2fc))) if len(log2fc) else 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    # size: linearly scaled so hubs are visibly larger, floor so low-degree
    # nodes stay visible
    deg_min, deg_max = degree.min(), degree.max()
    if deg_max > deg_min:
        size_norm = (degree - deg_min) / (deg_max - deg_min)
    else:
        size_norm = np.zeros_like(degree, dtype=float)
    node_sizes = 12 + size_norm * 340

    node_list = list(LCC.nodes())
    color_list = [LCC.nodes[n]["log2FC"] for n in node_list]

    fig, ax = plt.subplots(figsize=(17, 15))
    ax.set_facecolor(CHROME["surface"])

    # edges: dense network -> thin, low-alpha lines drawn first
    nx.draw_networkx_edges(
        LCC, pos, ax=ax, edge_color=CHROME["baseline"], width=0.35, alpha=0.18,
    )

    nodes_artist = nx.draw_networkx_nodes(
        LCC, pos, ax=ax, nodelist=node_list,
        node_size=node_sizes, node_color=color_list, cmap=DIVERGING_CMAP,
        vmin=-vmax, vmax=vmax, linewidths=0.4, edgecolors=CHROME["ink_secondary"],
    )
    nodes_artist.set_zorder(3)

    # --- label top hub genes only ---
    # Hub proteins cluster tightly in the layout (that's the real structure --
    # they're the most interconnected nodes), so placing labels directly on
    # top of the nodes collides. Use adjustText to push labels apart and draw
    # thin leader lines back to the actual node position.
    top_hubs = sorted(LCC.nodes(data=True), key=lambda x: x[1]["degree"], reverse=True)[:N_HUB_LABELS]
    top_hub_names = [n for n, _ in top_hubs]

    # small dots marking the exact hub node center, so the leader line has a
    # precise anchor even after the label text is nudged away
    hub_x = [pos[n][0] for n in top_hub_names]
    hub_y = [pos[n][1] for n in top_hub_names]
    ax.scatter(hub_x, hub_y, s=6, color=CHROME["ink_primary"], zorder=5)

    texts = [
        ax.text(pos[n][0], pos[n][1], n, fontsize=10.5, fontweight="bold",
                 color=CHROME["ink_primary"], zorder=6)
        for n in top_hub_names
    ]
    adjust_text(
        texts, x=hub_x, y=hub_y, ax=ax,
        arrowprops=dict(arrowstyle="-", color=CHROME["ink_secondary"], lw=0.7, alpha=0.85),
        expand=(1.4, 1.7), force_text=(0.6, 0.8), force_static=(0.3, 0.3),
    )

    ax.set_axis_off()
    ax.set_title(
        f"STRING PPI network — DE-regulated seed proteins (mouse, taxid 10090)\n"
        f"largest connected component: {n_nodes_plotted} nodes / {n_edges_plotted} edges "
        f"({n_nodes_omitted} nodes in {n_components_omitted} smaller components omitted "
        f"for legibility)",
        fontsize=13, fontweight="semibold", loc="left",
    )

    # --- colorbar (log2FC) ---
    sm = ScalarMappable(norm=norm, cmap=DIVERGING_CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, shrink=0.55)
    cbar.set_label("log2 fold change (treated vs. control)", fontsize=9.5, color=CHROME["ink_secondary"])
    cbar.ax.tick_params(labelsize=8.5, color=CHROME["ink_muted"], labelcolor=CHROME["ink_secondary"])

    # --- node-size legend (degree) ---
    legend_degrees = sorted({int(deg_min), int(round((deg_min + deg_max) / 2)), int(deg_max)})
    legend_handles = []
    for d in legend_degrees:
        dn = (d - deg_min) / (deg_max - deg_min) if deg_max > deg_min else 0.0
        ms = np.sqrt(12 + dn * 340)
        legend_handles.append(
            Line2D([0], [0], marker="o", linestyle="", markerfacecolor=CHROME["baseline"],
                   markeredgecolor=CHROME["ink_secondary"], markersize=ms / 3.2, label=f"degree {d}")
        )
    size_legend = ax.legend(
        handles=legend_handles, title="node size ~\ndegree", loc="lower left",
        frameon=False, fontsize=8.5, title_fontsize=8.5, labelspacing=1.3,
        borderpad=0.8, handletextpad=1.0,
    )
    ax.add_artist(size_legend)

    fig.text(
        0.5, 0.965,
        "Exploratory network among candidate (not statistically confirmed) proteins — "
        "edges are STRING combined-score ≥ 0.4 associations, not experimental evidence "
        "from this dataset.",
        ha="center", va="top", fontsize=9, style="italic", color=CHROME["ink_muted"],
    )

    add_caveat(fig, text=CAVEAT_TEXT, y=0.01)

    png_rel, svg_rel = save_fig(fig, STEM, dpi=200, tight=True)
    print(f"[network_figure] saved -> {png_rel}, {svg_rel}")

    top_hub_summary = [
        {"symbol": n, "degree": int(d["degree"]), "log2FC": round(float(d["log2FC"]), 3),
         "regulated": d["regulated"]}
        for n, d in top_hubs
    ]

    manifest_entry = {
        "file": png_rel,
        "title": "STRING PPI network of DE-regulated proteins (mouse) — largest connected component",
        "caption": (
            f"Protein-protein interaction network from STRING v{meta['string_version']} "
            f"(species=10090, Mus musculus, required_score={meta['required_score']}) built "
            f"from the {meta['n_seeds']} DE-\"regulated\" seed proteins ({meta['n_seeds_up']} UP, "
            f"{meta['n_seeds_down']} DOWN; |log2FC|>=0.585). {meta['n_mapped']}/{meta['n_seeds']} "
            f"seeds ({meta['mapping_rate']:.1%}) mapped to a STRING identifier, yielding a "
            f"{meta['n_nodes_in_network']}-node / {meta['n_edges']}-edge network with "
            f"{meta['n_communities']} greedy-modularity communities. Only the largest connected "
            f"component ({n_nodes_plotted} nodes, {n_edges_plotted} edges) is drawn for "
            f"legibility; {n_nodes_omitted} nodes across {n_components_omitted} smaller "
            f"components are omitted. Node color = log2 fold change (blue = down, red = up); "
            f"node size is proportional to STRING degree; the top {N_HUB_LABELS} hub genes by "
            f"degree are labeled. IMPORTANT: this is an exploratory network over the "
            f"hypothesis-generating |log2FC|>=0.585 candidate set from an n=2 "
            f"technical-replicate SILAC design (no biological replication) in which 0/1938 "
            f"proteins pass FDR<0.05 (limma) — network structure reflects STRING's curated "
            f"evidence, not statistical significance in this dataset, and should not be read "
            f"as confirmed findings."
        ),
        "key_numbers": {
            "n_seeds": meta["n_seeds"],
            "n_seeds_up": meta["n_seeds_up"],
            "n_seeds_down": meta["n_seeds_down"],
            "n_mapped": meta["n_mapped"],
            "mapping_rate": meta["mapping_rate"],
            "required_score": meta["required_score"],
            "n_nodes_full_network": n_nodes_full,
            "n_edges_full_network": n_edges_full,
            "n_nodes_plotted": n_nodes_plotted,
            "n_edges_plotted": n_edges_plotted,
            "n_omitted": n_nodes_omitted,
            "n_edges_omitted": n_edges_omitted,
            "n_components_omitted": n_components_omitted,
            "n_communities": meta["n_communities"],
            "string_version": meta["string_version"],
            "top_hub_genes": top_hub_summary,
        },
    }

    existing = []
    if MANIFEST_NETWORK_PATH.exists():
        try:
            existing = json.loads(MANIFEST_NETWORK_PATH.read_text())
        except json.JSONDecodeError:
            existing = []
    by_file = {e["file"]: e for e in existing}
    by_file[manifest_entry["file"]] = manifest_entry
    MANIFEST_NETWORK_PATH.write_text(json.dumps(list(by_file.values()), indent=2))
    print(f"[network_figure] manifest written -> {MANIFEST_NETWORK_PATH}")

    print("\n=== SUMMARY ===")
    print(f"Plotted LCC: {n_nodes_plotted} nodes / {n_edges_plotted} edges "
          f"({n_nodes_omitted} nodes omitted in {n_components_omitted} components)")
    print(f"Top {N_HUB_LABELS} hub genes: {[h['symbol'] for h in top_hub_summary]}")


if __name__ == "__main__":
    main()
