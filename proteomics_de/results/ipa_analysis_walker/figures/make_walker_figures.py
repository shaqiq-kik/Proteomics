"""Two figures for the Dr. Walker update page (D22).

Reads only committed result files and writes PNG + SVG next to this script:
  walker_string_focus_network.png  the six focus proteins and their high-confidence STRING neighbours
  run2_top_pathways.png            the top IPA canonical pathways of run 2 (assigned values included)

Run:  python make_walker_figures.py   (needs matplotlib, networkx, pandas)
Nothing here is measured data: STRING links are prior knowledge, and the 5 assigned proteins carry a
stand-in log2FC of -6.64 (see DECISIONS_LOG D20-D22).
"""
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
RES = HERE.parent.parent  # proteomics_de/results
PAPER, INK, MUTED = "#f3eee2", "#1c2230", "#5d6272"
CORE, FRZB, ASSIGNED, RECOVERED = "#8b887c", "#1f6f6b", "#c4432b", "#b4801f"
ASSIGNED_GENES = {"Egf", "Ereg", "Efna4", "Tyro3", "Mfap4"}
FOCUS = ["Egf", "Ereg", "Efna4", "Tyro3", "Mfap4", "Frzb"]
MIN_SCORE = 0.7

plt.rcParams.update({
    "font.family": "DejaVu Serif", "text.color": INK, "axes.labelcolor": INK,
    "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": INK,
    "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    "svg.fonttype": "none",
})


def save(fig, name):
    fig.savefig(HERE / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(HERE / f"{name}.svg", bbox_inches="tight")
    plt.close(fig)


def network_figure():
    edges = pd.read_csv(RES / "enrichment_walker" / "string_edges.tsv", sep="\t")
    nodes = pd.read_csv(RES / "enrichment_walker" / "string_node_metrics.csv").set_index("symbol")
    hi = edges[edges.combined_score >= MIN_SCORE]
    keep = set(FOCUS)
    for _, r in hi.iterrows():
        if r.nodeA_symbol in FOCUS:
            keep.add(r.nodeB_symbol)
        if r.nodeB_symbol in FOCUS:
            keep.add(r.nodeA_symbol)
    g = nx.Graph()
    g.add_nodes_from(keep)
    for _, r in hi.iterrows():
        if r.nodeA_symbol in keep and r.nodeB_symbol in keep:
            g.add_edge(r.nodeA_symbol, r.nodeB_symbol, w=r.combined_score)
    lone = [n for n in FOCUS if g.degree(n) == 0]
    comps = [c for c in nx.connected_components(g) if len(c) > 1]
    comps.sort(key=len, reverse=True)
    main = g.subgraph(comps[0])
    others = [g.subgraph(c) for c in comps[1:]]
    tier_col = {"core": CORE, "partial": RECOVERED, "qualitative": ASSIGNED}

    def layout_main(sg):
        inner = [n for n in ("Egf",) if n in sg]
        mid = [n for n in ("Ereg", "Egfr", "Efna4") if n in sg]
        rest = sorted(n for n in sg if n not in inner + mid)
        return nx.shell_layout(sg, nlist=[inner, mid, rest], scale=1.0)

    def draw(ax, sg, pos, radial):
        ax.set_axis_off()
        nx.draw_networkx_edges(sg, pos, ax=ax, edge_color=MUTED, alpha=.55,
                               width=[0.9 + 3.0 * (sg[a][b]["w"] - MIN_SCORE) / (1 - MIN_SCORE) for a, b in sg.edges])
        cx = sum(x for x, _ in pos.values()) / len(pos)
        cy = sum(y for _, y in pos.values()) / len(pos)
        for n in sg:
            t = nodes.loc[n, "tier"] if n in nodes.index else "core"
            is_focus, assigned = n in FOCUS, n in ASSIGNED_GENES
            face = ASSIGNED if assigned else (FRZB if n == "Frzb" else tier_col.get(t, CORE))
            ax.scatter([pos[n][0]], [pos[n][1]], s=900 if is_focus else 380,
                       marker="s" if assigned else "o", c=face, edgecolors=INK if is_focus else PAPER,
                       linewidths=2.2 if is_focus else 1.0, hatch="////" if assigned else None, zorder=3)
            arrow = ""
            if n in nodes.index:
                arrow = "\u25B2 " if nodes.loc[n, "regulated"] == "UP" else "\u25BC "
            dx, dy = pos[n][0] - cx, pos[n][1] - cy
            norm = (dx * dx + dy * dy) ** 0.5 or 1.0
            if radial and not is_focus:
                off = (dx / norm * 30, dy / norm * 22)
                ha = "left" if dx > 0.15 else ("right" if dx < -0.15 else "center")
            else:
                off, ha = (0, 24 if is_focus else 17), "center"
            ax.annotate(arrow + n.upper(), pos[n], xytext=off, textcoords="offset points", ha=ha, va="center",
                        fontsize=11.5 if is_focus else 8.8, fontweight="bold" if is_focus else "normal", zorder=5,
                        color=INK, bbox=dict(boxstyle="round,pad=0.12", fc=PAPER, ec="none", alpha=.9))
        ax.margins(.2)

    fig = plt.figure(figsize=(11.2, 7.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.15, 1], height_ratios=[1, 1], wspace=.05, hspace=.12)
    ax_main = fig.add_subplot(gs[:, 0])
    draw(ax_main, main, layout_main(main), radial=True)
    for k, sg in enumerate(others[:2]):
        ax = fig.add_subplot(gs[k, 1])
        pos = nx.circular_layout(sg, scale=.6) if len(sg) > 2 else {n: (i * 1.4 - .7, 0) for i, n in enumerate(sorted(sg))}
        draw(ax, sg, pos, radial=False)
    if lone:
        best = {n: edges[(edges.nodeA_symbol == n) | (edges.nodeB_symbol == n)].combined_score.max() for n in lone}
        txt = "\n".join(f"{n.upper()}: no link \u2265 {MIN_SCORE}\n(its best link is {best[n]:.2f})" for n in lone)
        fig.text(0.78, 0.075, txt, fontsize=9, color=MUTED, va="bottom")
    fig.suptitle("EGF, EREG, EFNA4, TYRO3, MFAP4 and FRZB in the STRING network", x=0.01, ha="left",
                 fontsize=15, fontweight="bold", y=0.985)
    legend = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=ASSIGNED, markeredgecolor=INK, markersize=11,
               label="assigned stand-in value (not measured)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CORE, markersize=9, label="neighbour: core protein"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=RECOVERED, markersize=9,
               label="neighbour: recovered protein (imputed replicate)"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=INK, markeredgecolor=INK, markersize=8,
               label="higher after testosterone"),
        Line2D([0], [0], marker="v", color="none", markerfacecolor=INK, markeredgecolor=INK, markersize=8,
               label="lower after testosterone"),
    ]
    fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.01, 0.045), frameon=True, facecolor=PAPER,
               edgecolor=MUTED, fontsize=9)
    fig.text(0.01, 0.0,
             f"Links with STRING combined score \u2265 {MIN_SCORE} (thicker = stronger). \u25B2 up, \u25BC down after testosterone.\n"
             "STRING links are prior knowledge from databases and literature, not evidence from this data. "
             "Mouse proteins; 968-protein list.",
             fontsize=8.8, color=MUTED, va="top")
    save(fig, "walker_string_focus_network")
    print("network: main", main.number_of_nodes(), "nodes", main.number_of_edges(), "edges; others",
          [sorted(o.nodes) for o in others], "lone", lone)


def pathways_figure(top=15):
    fn = HERE.parent / "run2_assigned_included_canonical_pathways.txt"
    lines = fn.read_text(encoding="utf-8", errors="replace").splitlines()
    i = next(k for k, l in enumerate(lines) if "\t" in l)
    hdr = [h.strip() for h in lines[i].split("\t")]
    rows = [dict(zip(hdr, l.split("\t"))) for l in lines[i + 1:] if l.strip()]
    df = pd.DataFrame(rows)
    df["logp"] = pd.to_numeric(df["-log(p-value)"], errors="coerce")
    df["z"] = pd.to_numeric(df["z-score"], errors="coerce")
    df = df.dropna(subset=["logp"]).sort_values("logp", ascending=False).head(top).iloc[::-1]
    df["has_assigned"] = df["Molecules"].apply(
        lambda s: bool({m.strip().upper() for m in s.split(",")} & {g.upper() for g in ASSIGNED_GENES}))

    def col(z):
        if pd.isna(z) or abs(z) < 2:
            return "#c9c4b4"
        return "#d9772b" if z > 0 else "#2f5fa8"

    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    y = range(len(df))
    bars = ax.barh(list(y), df["logp"], color=[col(z) for z in df["z"]], edgecolor=INK, linewidth=0.6, height=.72)
    for b, h in zip(bars, df["has_assigned"]):
        if h:
            b.set_hatch("////")
    ax.set_yticks(list(y))
    ax.set_yticklabels(df["Ingenuity Canonical Pathways"].str.slice(0, 54), fontsize=9.5)
    for yi, (lp, z) in enumerate(zip(df["logp"], df["z"])):
        ax.text(lp + 0.25, yi, "z " + ("n/a" if pd.isna(z) else f"{z:+.2f}"), va="center", fontsize=9, color=INK)
    ax.set_xlabel("−log10 p (inflated: IPA default reference set; compare the order and z-scores, not the p-values)", fontsize=9.2)
    ax.set_xlim(0, df["logp"].max() * 1.13)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title("IPA run 2: top 15 canonical pathways", fontsize=14.5, fontweight="bold", loc="left", pad=30)
    sub = ("None of these 15 pathways lists an assigned stand-in protein, so this chart does not rest on the placeholder values."
           if not df["has_assigned"].any() else
           "Hatched pathways list an assigned stand-in protein; their z-score rests partly on the −6.64 placeholder.")
    ax.text(0, 1.025, sub, transform=ax.transAxes, fontsize=9.2, color=MUTED, va="bottom")
    legend = [Patch(fc="#d9772b", ec=INK, label="predicted activated (z ≥ 2)"),
              Patch(fc="#2f5fa8", ec=INK, label="predicted inhibited (z ≤ −2)"),
              Patch(fc="#c9c4b4", ec=INK, label="no clear direction (|z| < 2)")]
    if df["has_assigned"].any():
        legend.append(Patch(fc="#c9c4b4", ec=INK, hatch="////", label="lists an assigned stand-in protein"))
    ax.legend(handles=legend, loc="lower right", frameon=True, facecolor=PAPER, edgecolor=MUTED, fontsize=9)
    fig.text(0.0, -0.02,
             "Run 2 = the 968-protein list with the five assigned proteins analysed (934 molecules). "
             "Grey bars: IPA predicts no direction (|z| < 2).\n"
             "Nothing in this list is statistically significant (0 of 968 pass FDR < 0.05).",
             fontsize=8.8, color=MUTED, va="top")
    save(fig, "run2_top_pathways")
    print(df[["Ingenuity Canonical Pathways", "logp", "z", "has_assigned"]].iloc[::-1].to_string(index=False))


if __name__ == "__main__":
    network_figure()
    pathways_figure()
