"""Build report_v2_16nov.html (update page for Dr. Walker, 16 Nov 2026) and its three figures.

Run from the repo root:   python tools/build_report_v2.py      (needs pandas, numpy, matplotlib, networkx)

Every number on the page is read from the exports at build time:
  proteomics_de/results/ipa_input_walker.csv                       the 721 list (values, p, tiers)
  proteomics_de/results/enrichment_walker/ora_meta.json            g:Profiler
  proteomics_de/results/enrichment_walker/string_meta.json         STRING
  proteomics_de/results/enrichment_walker/string_node_metrics.csv  STRING degrees
  proteomics_de/results/enrichment_walker/string_edges.tsv         STRING links (figure 2)
  proteomics_de/results/ipa_analysis_walker/flagged_*.csv          flagged counts + table
  proteomics_de/results/ipa_analysis_walker/walker721_*/walker716_* IPA exports (figure 1, 3, androgen check)
  Pilot Project/CLEANED Silac Proteomics Soluble Factors.xlsx      sample IDs (pairing note)
The only hand-typed facts are the ones that live in no file (the 247 dropped proteins and the list of
five assigned gene names, both from Dr. Walker's 11 Nov reply) and the pinned commit used in the file links.
Writes three PNGs (walker721_*.png, figures/) and the self-contained HTML (base64 images, no scripts).
"""
from __future__ import annotations

import base64
import csv
import html
import json
import os
import re
import subprocess
import textwrap
import zipfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "proteomics_de" / "results"
IPA = RES / "ipa_analysis_walker"
ENR = RES / "enrichment_walker"
FIG = IPA / "figures"
OUT = ROOT / "report_v2_16nov.html"
SHA = "07461955bf428713e74de93291b961efbd32a9a6"  # pushed commit that holds every file linked from the page
GH = f"https://github.com/shaqiq-kik/Proteomics/blob/{SHA}/"

ASSIGNED_GENES = ["EGF", "EREG", "EFNA4", "TYRO3", "MFAP4"]  # Dr. Walker's reply, 11 Nov
N_DROPPED = 247  # the recovered proteins Dr. Walker asked to drop; lives in no export of this list
FOCUS = ["Egf", "Ereg", "Efna4", "Tyro3", "Mfap4", "Frzb"]

PAPER, INK, MUTED = "#f3eee2", "#1c2230", "#5d6272"
CORE, FRZB, ASSIGNED = "#8b887c", "#1f6f6b", "#c4432b"
plt.rcParams.update({
    "font.family": "DejaVu Serif", "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": MUTED,
    "xtick.color": MUTED, "ytick.color": INK, "figure.facecolor": PAPER, "axes.facecolor": PAPER,
    "savefig.facecolor": PAPER,
})


# ------------------------------------------------------------------ readers
def read_ipa(path: Path) -> pd.DataFrame:
    """IPA exports: a copyright line and a blank line, then a tab-separated table."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    i = next(k for k, ln in enumerate(lines) if "\t" in ln)
    hdr = [h.strip() for h in lines[i].split("\t")]
    rows = []
    for ln in lines[i + 1:]:
        if ln.strip():
            c = ln.split("\t")
            c += [""] * (len(hdr) - len(c))
            rows.append({h: v.strip() for h, v in zip(hdr, c)})
    return pd.DataFrame(rows)


def fnum(s):
    try:
        v = float(s)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def esc(s) -> str:
    return html.escape(str(s), quote=True)


def fmt_int(n) -> str:
    return str(int(n))  # plain digits (no thousands separator), e.g. 6870, 2277


def save_fig(fig, name: str) -> Path:
    """Write the figure unless it already exists (the pushed PNGs stay untouched and are what the page embeds).
    Set REFRESH_FIGURES=1 to redraw them."""
    p = FIG / name
    if os.environ.get("REFRESH_FIGURES") or not p.exists():
        fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


def b64(p: Path) -> str:
    return base64.b64encode(p.read_bytes()).decode("ascii")


# ------------------------------------------------------------------ data
inp = pd.read_csv(RES / "ipa_input_walker.csv")
inp["sym"] = inp["Gene names"].str.upper()
ora = json.loads((ENR / "ora_meta.json").read_text())
strm = json.loads((ENR / "string_meta.json").read_text())
nodes = pd.read_csv(ENR / "string_node_metrics.csv")
summ = {r["table"]: r for r in csv.DictReader(open(IPA / "flagged_summary.csv", newline=""))}
flag = list(csv.DictReader(open(IPA / "flagged_dependent_results.csv", newline="", encoding="utf-8")))

n_list = len(inp)
n_assigned = int((inp.value_source == "assigned").sum())
n_core = int((inp.tier == "core").sum())
n_frzb = n_list - n_core - n_assigned
assert n_frzb == 1 and set(inp[inp.value_source == "assigned"].sym) == set(ASSIGNED_GENES)
lfc_assigned = inp[inp.value_source == "assigned"]["log2FC"].round(6).unique()
assert len(lfc_assigned) == 1
LFC_ASSIGNED = f"{lfc_assigned[0]:.6f}"  # -6.643856
frzb = inp[inp.sym == "FRZB"].iloc[0]
LFC_FRZB, P_FRZB, ADJ_FRZB = float(frzb.log2FC), float(frzb.p_value), float(frzb.adj_p_value)
n_up, n_down = int((inp.regulated == "UP").sum()), int((inp.regulated == "DOWN").sum())
assert n_up == ora["n_query_submitted_up"] and n_down == ora["n_query_submitted_down"]
assert ora["n_terms_up"] == 0 and ora["n_terms_down"] == 0
n_gp_terms = ora["n_terms_up"] + ora["n_terms_down"]
n_nodes, n_edges = strm["n_nodes_in_network"], strm["n_edges"]
m = re.search(r"(\d+)/(\d+) proteins pass per-protein FDR<0\.05", strm["caveat"])
assert m, "FDR statement not found in string_meta.json"
FDR_PASS, FDR_TOTAL = int(m.group(1)), int(m.group(2))
min_adj_list = float(inp.adj_p_value.min())
deg = nodes.set_index("symbol")["degree"]
focus_deg = {g: int(deg[g]) for g in FOCUS}

flag_tot = {t: int(r["flagged_total"]) for t, r in summ.items()}
FLAG_ALL = sum(flag_tot.values())
assert FLAG_ALL == len(flag)
flips = {t: int(r["sign_flips"]) for t, r in summ.items()}
cross = {t: int(r["crossed_abs_z_2"]) for t, r in summ.items()}
FLIPS_ALL, CROSS_ALL = sum(flips.values()), sum(cross.values())
assert FLIPS_ALL == sum(r["sign_flip"] == "yes" for r in flag)
assert CROSS_ALL == sum(r["crosses_abs_z_2"] == "yes" for r in flag)
rows721 = {t: int(r["rows_721"]) for t, r in summ.items()}
rows716 = {t: int(r["rows_716"]) for t, r in summ.items()}
dis_rows = rows721["diseases_and_functions"]
assert dis_rows == 1000 == rows716["diseases_and_functions"]

TABLE_LABEL = {"canonical_pathways": "Canonical pathways", "upstream_regulators": "Upstream regulators",
               "diseases_and_functions": "Diseases and functions"}

# sample IDs from the cleaned sheet (xlsx is a zip; stdlib is enough)
xlsx = ROOT / "Pilot Project" / "CLEANED Silac Proteomics Soluble Factors.xlsx"
with zipfile.ZipFile(xlsx) as z:
    ss = z.read("xl/sharedStrings.xml").decode("utf-8")
samples = {(g, int(r)): sid for g, r, sid in re.findall(r"(Vehicle|Testosterone)_Rep(\d)_(\d+)", ss)}
PAIR = {k: (samples[("Vehicle", k)], samples[("Testosterone", k)]) for k in (1, 2)}
assert PAIR == {1: ("31579", "31578"), 2: ("31581", "31580")}, PAIR


# ------------------------------------------------------------------ androgen check
def upstream_z(run_file: str, name: str):
    df = read_ipa(IPA / run_file)
    r = df[df["Upstream Regulator"] == name]
    return fnum(r.iloc[0]["Activation z-score"]) if len(r) else None


def androgen_claim():
    flagged_names = {r["name"] for r in flag if r["table"] == "upstream_regulators"}
    z = {}
    for nm in ("AR", "metribolone"):
        a = upstream_z("walker721_with_assigned_upstream_regulators.txt", nm)
        b = upstream_z("walker716_without_assigned_upstream_regulators.txt", nm)
        if a is None or b is None or nm in flagged_names or abs(a - b) > 0.0005:
            return None
        z[nm] = (a, b)
    return z


ANDRO = androgen_claim()


# ------------------------------------------------------------------ figure 1: top pathways
def fig_pathways(top=15) -> Path:
    df = read_ipa(IPA / "walker721_with_assigned_canonical_pathways.txt")
    df["logp"] = pd.to_numeric(df["-log(p-value)"], errors="coerce")
    df["z"] = pd.to_numeric(df["z-score"], errors="coerce")
    df = df.dropna(subset=["logp"]).sort_values("logp", ascending=False).head(top).iloc[::-1]
    asg = set(ASSIGNED_GENES)
    df["has_assigned"] = df["Molecules"].apply(lambda s: bool({x.strip().upper() for x in s.split(",")} & asg))
    n_has = int(df.has_assigned.sum())
    global PATH_N_HAS
    PATH_N_HAS = n_has

    def col(z):
        if pd.isna(z) or abs(z) < 2:
            return "#c9c4b4"
        return "#d9772b" if z > 0 else "#2f5fa8"

    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    y = list(range(len(df)))
    bars = ax.barh(y, df["logp"], color=[col(z) for z in df["z"]], edgecolor=INK, linewidth=0.6, height=.72)
    for b, h in zip(bars, df["has_assigned"]):
        if h:
            b.set_hatch("////")
    ax.set_yticks(y)
    ax.set_yticklabels(["\n".join(textwrap.wrap(t, 48)) for t in df["Ingenuity Canonical Pathways"]], fontsize=9)
    for yi, (lp, z) in enumerate(zip(df["logp"], df["z"])):
        ax.text(lp + 0.25, yi, "z " + ("n/a" if pd.isna(z) else f"{z:+.2f}"), va="center", fontsize=9)
    ax.set_xlabel("−log10 p (inflated: IPA default reference set; compare the order and z, not the p-values)",
                  fontsize=9.2)
    ax.set_xlim(0, df["logp"].max() * 1.13)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_title(f"IPA run 2 ({n_list}-protein list): top {top} canonical pathways", fontsize=14, fontweight="bold",
                 loc="left", pad=30)
    sub = (f"{n_has} of these {top} pathways list an assigned stand-in protein (hatched); their z rests partly on the "
           f"{float(LFC_ASSIGNED):.2f} placeholder." if n_has else
           f"None of these {top} pathways lists an assigned stand-in protein.")
    ax.text(0, 1.025, sub, transform=ax.transAxes, fontsize=9.2, color=MUTED, va="bottom")
    legend = [Patch(fc="#d9772b", ec=INK, label="predicted activated (z ≥ 2)"),
              Patch(fc="#2f5fa8", ec=INK, label="predicted inhibited (z ≤ −2)"),
              Patch(fc="#c9c4b4", ec=INK, label="no clear direction (|z| < 2)")]
    if n_has:
        legend.append(Patch(fc="#c9c4b4", ec=INK, hatch="////", label="lists an assigned stand-in protein"))
    ax.legend(handles=legend, loc="lower right", frameon=True, facecolor=PAPER, edgecolor=MUTED, fontsize=9)
    fig.text(0.0, -0.02,
             f"Run 2 = the {n_list}-protein list with the five assigned proteins analysed. "
             f"None is statistically significant ({FDR_PASS} of {FDR_TOTAL} proteins pass FDR < 0.05).",
             fontsize=8.8, color=MUTED, va="top")
    return save_fig(fig, "walker721_top_pathways.png")


# ------------------------------------------------------------------ figure 2: STRING focus network
def fig_network(min_score=0.7) -> Path:
    edges = pd.read_csv(ENR / "string_edges.tsv", sep="\t")
    nd = nodes.set_index("symbol")
    hi = edges[edges.combined_score >= min_score]
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
    global NET_ASG_DRAWN, NET_FRZB_DRAWN
    comps = sorted((c for c in nx.connected_components(g) if len(c) > 1), key=len, reverse=True)
    main = g.subgraph(comps[0])
    others = [g.subgraph(c) for c in comps[1:]]
    _drawn = set(main.nodes)
    for _sg in others[:2]:
        _drawn |= set(_sg.nodes)
    NET_ASG_DRAWN = sorted(n for n in _drawn if n.upper() in ASSIGNED_GENES)
    NET_FRZB_DRAWN = "Frzb" in _drawn

    def layout_main(sg):
        inner = [n for n in ("Egf",) if n in sg]
        mid = [n for n in ("Ereg", "Egfr", "Efna4") if n in sg]
        rest = sorted(n for n in sg if n not in inner + mid)
        return nx.shell_layout(sg, nlist=[inner, mid, rest], scale=1.0)

    def draw(ax, sg, pos, radial):
        ax.set_axis_off()
        nx.draw_networkx_edges(sg, pos, ax=ax, edge_color=MUTED, alpha=.55,
                               width=[0.9 + 3.0 * (sg[a][b]["w"] - min_score) / (1 - min_score) for a, b in sg.edges])
        cx = sum(x for x, _ in pos.values()) / len(pos)
        cy = sum(yy for _, yy in pos.values()) / len(pos)
        for n in sg:
            is_focus, assigned = n in FOCUS, n.upper() in ASSIGNED_GENES
            face = ASSIGNED if assigned else (FRZB if n == "Frzb" else CORE)
            ax.scatter([pos[n][0]], [pos[n][1]], s=900 if is_focus else 380, marker="s" if assigned else "o",
                       c=face, edgecolors=INK if is_focus else PAPER, linewidths=2.2 if is_focus else 1.0,
                       hatch="////" if assigned else None, zorder=3)
            arrow = ""
            if n in nd.index:
                arrow = "▲ " if nd.loc[n, "regulated"] == "UP" else "▼ "
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
    draw(fig.add_subplot(gs[:, 0]), main, layout_main(main), radial=True)
    for k, sg in enumerate(others[:2]):
        ax = fig.add_subplot(gs[k, 1])
        pos = (nx.circular_layout(sg, scale=.6) if len(sg) > 2
               else {n: (i * 1.4 - .7, 0) for i, n in enumerate(sorted(sg))})
        draw(ax, sg, pos, radial=False)
    if lone:
        best = {n: edges[(edges.nodeA_symbol == n) | (edges.nodeB_symbol == n)].combined_score.max() for n in lone}
        txt = "\n".join(f"{n.upper()}: no link ≥ {min_score}\n(its best link is {best[n]:.2f})" for n in lone)
        fig.text(0.78, 0.075, txt, fontsize=9, color=MUTED, va="bottom")
    fig.suptitle("EGF, EREG, EFNA4, TYRO3, MFAP4 and FRZB in the STRING network", x=0.01, ha="left",
                 fontsize=15, fontweight="bold", y=0.985)
    legend = [
        Line2D([0], [0], marker="s", color="none", markerfacecolor=ASSIGNED, markeredgecolor=INK, markersize=11,
               label="assigned stand-in value (not measured)"),
    ]
    if NET_FRZB_DRAWN:
        legend.append(Line2D([0], [0], marker="o", color="none", markerfacecolor=FRZB, markeredgecolor=INK,
                             markersize=10, label="FRZB (measured, not FDR-significant)"))
    legend += [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CORE, markersize=9, label="neighbour: core protein"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=INK, markeredgecolor=INK, markersize=8,
               label="higher after testosterone"),
        Line2D([0], [0], marker="v", color="none", markerfacecolor=INK, markeredgecolor=INK, markersize=8,
               label="lower after testosterone"),
    ]
    fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.01, 0.045), frameon=True, facecolor=PAPER,
               edgecolor=MUTED, fontsize=9)
    fig.text(0.01, 0.0,
             f"Links with STRING combined score ≥ {min_score} (thicker = stronger). ▲ up, ▼ down after "
             "testosterone.\nSTRING links are prior knowledge from databases and literature, not evidence from this "
             f"data. Mouse proteins; {n_list}-protein list ({n_nodes} mapped to STRING).",
             fontsize=8.8, color=MUTED, va="top")
    return save_fig(fig, "walker721_string_focus_network.png")


# ------------------------------------------------------------------ figure 3: z scatter
def fig_scatter() -> tuple[Path, int, int]:
    d = pd.DataFrame([r for r in flag if r["table"] in ("upstream_regulators", "diseases_and_functions")])
    d["x"] = d["z_without_assigned_716"].map(fnum)
    d["y"] = d["z_with_assigned_721"].map(fnum)
    d = d.dropna(subset=["x", "y"])
    d["asg"] = d["lists_assigned_protein"].astype(str).str.len() > 0
    lim = float(max(d.x.abs().max(), d.y.abs().max())) * 1.06
    fig, ax = plt.subplots(figsize=(8.4, 7.0))
    ax.axhline(0, color=MUTED, lw=.5)
    ax.axvline(0, color=MUTED, lw=.5)
    for v in (-2, 2):
        ax.axhline(v, color=ASSIGNED, lw=.9, ls=":")
        ax.axvline(v, color=ASSIGNED, lw=.9, ls=":")
    ax.plot([-lim, lim], [-lim, lim], color=INK, lw=.9, ls="--")
    no, yes = d[~d.asg], d[d.asg]
    if len(no):
        ax.scatter(no.x, no.y, s=16, c="#9a978a", alpha=.55, edgecolors="none", label=f"does not list an assigned protein ({len(no)})")
    if len(yes):
        ax.scatter(yes.x, yes.y, s=26, c=ASSIGNED, alpha=.7, edgecolors="none", label=f"lists an assigned protein ({len(yes)})")
    reg = d[(d.table == "upstream_regulators") & (d.crosses_abs_z_2 == "yes")]
    offs = [(8, 10), (8, -14), (-8, 10), (-8, -14), (10, 0)]
    for k, (_, r) in enumerate(reg.iterrows()):
        dx, dy = offs[k % len(offs)]
        ax.annotate(r["name"], (r.x, r.y), xytext=(dx, dy), textcoords="offset points", fontsize=8.5,
                    ha="left" if dx > 0 else "right", arrowprops=dict(arrowstyle="-", color=MUTED, lw=.6),
                    bbox=dict(boxstyle="round,pad=0.1", fc=PAPER, ec="none", alpha=.85))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("z without the five assigned proteins (716 run)")
    ax.set_ylabel("z with the five assigned proteins (721 run)")
    ax.set_title("Flagged regulators and diseases: z with vs without the five", fontsize=13, fontweight="bold",
                 loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(loc="upper left", frameon=True, facecolor=PAPER, edgecolor=MUTED, fontsize=8.8)
    fig.text(0.0, -0.01,
             "Dashed line: no change. Dotted lines: |z| = 2. Labelled: regulators that cross |z| = 2.\n"
             f"Only results with a z-score in both runs are drawn ({len(d)} of "
             f"{flag_tot['upstream_regulators'] + flag_tot['diseases_and_functions']}).",
             fontsize=8.6, color=MUTED, va="top")
    return save_fig(fig, "walker721_z_scatter.png"), len(d), int(d.asg.sum())


P_PATH = fig_pathways()
P_NET = fig_network()
P_SCAT, N_SCAT, N_SCAT_ASG = fig_scatter()


# ------------------------------------------------------------------ flagged table
def zs(v: str) -> str:
    v = (v or "").strip()
    return "not in run" if fnum(v) is None else v


def pick_rows() -> list[dict]:
    chosen: list[dict] = []
    seen = set()

    def add(r, why):
        k = (r["table"], r["name"])
        if k not in seen:
            seen.add(k)
            chosen.append({**r, "why": why})

    def why_cross(r):
        a, b = fnum(r["z_with_assigned_721"]), fnum(r["z_without_assigned_716"])
        if a is None:
            return "|z| reaches 2 only without the five (no z with them)"
        if b is None:
            return "|z| reaches 2 only with the five (no z without them)"
        return ("|z| reaches 2 only without the five" if abs(b) >= 2 else "|z| reaches 2 only with the five")

    pr = [r for r in flag if r["table"] in ("canonical_pathways", "upstream_regulators")]
    for r in pr:
        if r["crosses_abs_z_2"] == "yes":
            add(r, why_cross(r))
    with_delta = sorted((r for r in pr if fnum(r["delta_z"]) is not None and (r["table"], r["name"]) not in seen),
                        key=lambda r: -abs(fnum(r["delta_z"])))
    for r in with_delta[:5]:
        add(r, f"large shift (Δz {fnum(r['delta_z']):+.3f})")
    dis = sorted((r for r in flag if r["table"] == "diseases_and_functions" and r["crosses_abs_z_2"] == "yes"
                  and fnum(r["delta_z"]) is not None), key=lambda r: -abs(fnum(r["delta_z"])))
    for r in dis[:5]:
        add(r, why_cross(r) + f" (Δz {fnum(r['delta_z']):+.3f})")
    assert 15 <= len(chosen) <= 25, len(chosen)
    return chosen


FLAG_ROWS = pick_rows()
N_PR_CROSS = sum(1 for r in FLAG_ROWS if r["crosses_abs_z_2"] == "yes" and r["table"] != "diseases_and_functions")
assert N_PR_CROSS == cross["canonical_pathways"] + cross["upstream_regulators"]


def flagged_table() -> str:
    out = ['<div class="tw"><table id="flagged-results"><thead><tr><th>Table</th><th>Name</th>'
           '<th>z with five (721)</th><th>z without (716)</th><th>Assigned protein(s)</th><th>Why listed</th></tr></thead><tbody>']
    for r in FLAG_ROWS:
        asg = r["lists_assigned_protein"].replace(",", ", ") or "none listed"
        out.append(
            f"<tr><td>{esc(TABLE_LABEL[r['table']])}</td><td class=\"nm\">{esc(r['name'])}</td>"
            f"<td class=\"num\">{esc(zs(r['z_with_assigned_721']))}</td><td class=\"num\">{esc(zs(r['z_without_assigned_716']))}</td>"
            f"<td>{esc(asg)}</td><td>{esc(r['why'])}</td></tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def counts_table() -> str:
    out = ['<div class="tw"><table class="counts"><thead><tr><th>Table</th><th>Rows, 721 run</th><th>Rows, 716 run</th>'
           '<th>List an assigned protein</th><th>z changed or new</th><th>Sign changes</th><th>|z| = 2 crossings</th>'
           '<th>Flagged</th></tr></thead><tbody>']
    for t, r in summ.items():
        out.append(f"<tr><td>{esc(TABLE_LABEL[t])}</td><td class=\"num\">{r['rows_721']}</td><td class=\"num\">{r['rows_716']}</td>"
                   f"<td class=\"num\">{r['lists_assigned_protein']}</td><td class=\"num\">{r['z_changed_or_new']}</td>"
                   f"<td class=\"num\">{r['sign_flips']}</td><td class=\"num\">{r['crossed_abs_z_2']}</td>"
                   f"<td class=\"num\"><b>{r['flagged_total']}</b></td></tr>")
    out.append(f"<tr><td><b>All three</b></td><td></td><td></td><td></td><td></td><td class=\"num\"><b>{FLIPS_ALL}</b></td>"
               f"<td class=\"num\"><b>{CROSS_ALL}</b></td><td class=\"num\"><b>{FLAG_ALL}</b></td></tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


# ------------------------------------------------------------------ paired-model sensitivity (pairing_sensitivity/*.csv)
PSD = IPA / "pairing_sensitivity"


def _ps(name):
    return list(csv.DictReader(open(PSD / f"{name}.csv", newline="", encoding="utf-8")))


_u, _p, _uv, _pv = (_ps(n) for n in ("unpaired_trend_robust", "paired_trend_robust", "unpaired_vanilla", "paired_vanilla"))
_nfdr = lambda rows: sum(float(r["adj"]) < 0.05 for r in rows)
_nraw = lambda rows: sum(float(r["p"]) < 0.05 for r in rows)
_minadj = lambda rows: min(float(r["adj"]) for r in rows)
_top10 = lambda rows: [r["id"] for r in sorted(rows, key=lambda r: float(r["p"]))[:10]]
_g1 = lambda r: r["gene"].split(";")[0]
PS_OVERLAP = len(set(_top10(_u)) & set(_top10(_p)))
PS_ENTERS = [_g1(r) for r in _p if r["id"] in set(_top10(_p)) - set(_top10(_u))]
PS_DROPS = [r["gene"] for r in _u if r["id"] in set(_top10(_u)) - set(_top10(_p))]
_fu = next(r for r in _u if r["gene"].upper() == "FRZB")
_fp = next(r for r in _p if r["gene"].upper() == "FRZB")
_list_ids = set(inp["UniProt Accession Number"].astype(str))
_list_acc = _list_ids | {a.strip() for x in _list_ids for a in x.split(";")}
PS_HITS = [r for r in _p if float(r["adj"]) < 0.05]


def _hit_text(r):
    inlist = r["id"] in _list_acc or any(a in _list_acc for a in r["id"].split(";"))
    extra = len(r["gene"].split(";")) - 1
    name = esc(_g1(r)) + (f" (shared by {extra + 1} near-identical histone entries)" if extra else "")
    note = "in the 721 list" if inlist else "<b>not in the 721 list</b>"
    if int(r["n_imputed"]) > 0:
        note += f"; {int(r['n_imputed'])} of its four values is imputed, so it leans on a filled-in number"
    return f"<li>{name}: {note}. Paired adjusted p {float(r['adj']):.4f}.</li>"


_inlist_genes = {_g1(r) for r in PS_HITS if r["id"] in _list_acc or any(a in _list_acc for a in r["id"].split(";"))}
assert _inlist_genes == {"Hist1h2bc", "Hist1h2ba", "Abhd14b"}, _inlist_genes
assert "Tuba1a" not in _inlist_genes
ps_hits_html = "".join(_hit_text(r) for r in PS_HITS)
ps_html = f"""
<h2><span class="n">5</span> Paired model: sensitivity check</h2>
<p class="sub">The main analysis stays unpaired: its list, fold changes, g:Profiler, STRING and IPA results are unchanged. I ran this only to ask what happens if each testosterone sample is paired with the control from its own experiment (model: experiment + group, limma trend with robust fitting).</p>
<div class="tw"><table class="counts"><thead><tr><th>Measure</th><th>Unpaired (as reported)</th><th>Paired</th></tr></thead><tbody>
<tr><td>Proteins with FDR &lt; 0.05</td><td class="num">{_nfdr(_u)}</td><td class="num">{_nfdr(_p)}</td></tr>
<tr><td>Proteins with raw p &lt; 0.05</td><td class="num">{_nraw(_u)}</td><td class="num">{_nraw(_p)}</td></tr>
<tr><td>Smallest adjusted p</td><td class="num">{_minadj(_u):.3f}</td><td class="num">{_minadj(_p):.4f}</td></tr>
<tr><td>FRZB p / adjusted p</td><td class="num">{float(_fu['p']):.4f} / {float(_fu['adj']):.3f}</td><td class="num">{float(_fp['p']):.4f} / {float(_fp['adj']):.3f}</td></tr>
<tr><td>Top 10 by raw p</td><td class="num" colspan="2">{PS_OVERLAP} of 10 shared; {esc(', '.join(PS_ENTERS))} enters, {esc(', '.join(PS_DROPS))} drops</td></tr>
<tr><td>Plain eBayes: FDR &lt; 0.05 / raw p &lt; 0.05</td><td class="num">{_nfdr(_uv)} / {_nraw(_uv)}</td><td class="num">{_nfdr(_pv)} / {_nraw(_pv)}</td></tr>
</tbody></table></div>
<h3>The {len(PS_HITS)} paired FDR hits</h3>
<ul class="tight">{ps_hits_html}</ul>
<p><b>Read with care.</b> With 2 pairs there is only 1 residual degree of freedom, so these hits are leads, not findings. They vanish under plain eBayes ({_nfdr(_pv)} hits). This is not a replacement for the main analysis. I would like your view on whether to adopt it.</p>
"""

# ------------------------------------------------------------------ files section
def link(path: str, name: str, desc: str, kind: str) -> str:
    """One Files row: a short human name as the link text (raw path only in the hover title), then a description."""
    assert (ROOT / path).exists(), path
    assert subprocess.run(["git", "-C", str(ROOT), "cat-file", "-e", f"{SHA}:{path}"]).returncode == 0, \
        f"{path} is not in pinned commit {SHA[:7]}"
    return (f'<li><a class="f" href="{GH}{path}" target="_blank" rel="noopener" title="{esc(path)}">{esc(name)}</a>'
            f'<span class="ty">{kind}</span><span class="d">{esc(desc)}</span></li>')


_R = "proteomics_de/results/"
_A = _R + "ipa_analysis_walker/"
_P = _A + "pairing_sensitivity/"
files_html = "<ul class=\"links\">" + "".join([
    link(_R + "ipa_input_walker.csv", "IPA input list (721)", "The 721 list. The last column, measured or assigned, marks the five stand-ins.", "spreadsheet"),
    link(_R + "enrichment_walker/ora_meta.json", "g:Profiler settings and counts", "g:Profiler settings and result counts.", "data"),
    link(_R + "enrichment_walker/string_meta.json", "STRING network settings", "STRING settings and network size.", "data"),
    link(_R + "enrichment_walker/string_focus_proteins.csv", "STRING focus proteins", "The six focus proteins and their STRING neighbours.", "spreadsheet"),
    link(_A + "walker721_with_assigned_canonical_pathways.txt", "IPA run 2: canonical pathways", "IPA run 2 (721, with the five).", "IPA export"),
    link(_A + "walker721_with_assigned_upstream_regulators.txt", "IPA run 2: upstream regulators", "IPA run 2 (721, with the five).", "IPA export"),
    link(_A + "walker721_with_assigned_diseases_and_functions.txt", "IPA run 2: diseases and functions", "IPA run 2 (721, with the five). First 1000 rows.", "IPA export"),
    link(_A + "walker721_with_assigned_networks.txt", "IPA run 2: networks", "IPA run 2 (721, with the five).", "IPA export"),
    link(_A + "walker716_without_assigned_canonical_pathways.txt", "IPA 716 comparison: canonical pathways", "Comparison run (716, without the five).", "IPA export"),
    link(_A + "walker716_without_assigned_upstream_regulators.txt", "IPA 716 comparison: upstream regulators", "Comparison run (716, without the five).", "IPA export"),
    link(_A + "walker716_without_assigned_diseases_and_functions.txt", "IPA 716 comparison: diseases and functions", "Comparison run (716, without the five). First 1000 rows.", "IPA export"),
    link(_A + "walker716_without_assigned_networks.txt", "IPA 716 comparison: networks", "Comparison run (716, without the five).", "IPA export"),
    link(_A + "flagged_dependent_results.csv", "Flagged results (all rows)", "Every flagged result with both z-scores.", "spreadsheet"),
    link(_A + "flagged_summary.csv", "Flagged results (counts)", "The flagged counts per table.", "spreadsheet"),
    link(_A + "figures/walker721_top_pathways.png", "Figure 1: top pathways", "Top canonical pathways of IPA run 2.", "image"),
    link(_A + "figures/walker721_string_focus_network.png", "Figure 2: STRING network", "The six focus proteins in the STRING network.", "image"),
    link(_A + "figures/walker721_z_scatter.png", "Figure 3: z-score scatter", "z with versus without the five, for flagged regulators and diseases.", "image"),
    link(_P + "refit.R", "Paired fit script (R)", "Refits the model with experiment in the design.", "script"),
    link(_P + "compare.R", "Paired vs unpaired comparison script (R)", "Compares the paired and unpaired fits.", "script"),
    link(_P + "unpaired_trend_robust.csv", "Unpaired fit, as reported", "Unpaired fit (trend, robust), the one used in the main analysis.", "spreadsheet"),
    link(_P + "paired_trend_robust.csv", "Paired fit", "Paired fit (trend, robust).", "spreadsheet"),
    link(_P + "unpaired_vanilla.csv", "Unpaired fit, plain eBayes", "Unpaired fit with plain eBayes.", "spreadsheet"),
    link(_P + "paired_vanilla.csv", "Paired fit, plain eBayes", "Paired fit with plain eBayes.", "spreadsheet"),
]) + "</ul>"

# ------------------------------------------------------------------ page
CSS = """
  :root{color-scheme:light;
    --paper:#f3eee2; --paper2:#ebe5d6; --card:#fbf8f0; --ink:#1c2230; --muted:#5d6272; --rule:#cfc7b4;
    --core:#8b887c; --frzb:#1f6f6b; --assigned:#c4432b; --recovered:#b4801f; --link:#1b4f8a;
    --shadow:0 1px 0 rgba(28,34,48,.06), 0 12px 28px -18px rgba(28,34,48,.35);
    --serif:"Iowan Old Style","Palatino Linotype",Palatino,"Book Antiqua",Georgia,serif;
    --mono:ui-monospace,"SF Mono","Cascadia Mono",Menlo,Consolas,monospace;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--serif);font-size:17px;line-height:1.55;
    background-image:radial-gradient(1200px 500px at 85% -10%,rgba(196,67,43,.07),transparent 60%),
                     radial-gradient(900px 500px at -10% 30%,rgba(31,111,107,.07),transparent 60%);}
  a{color:var(--link);text-underline-offset:3px;text-decoration-thickness:1px}
  .wrap{max-width:880px;margin:0 auto;padding:0 20px 80px}
  header.top{border-top:6px solid var(--ink);padding:22px 0 8px;display:flex;justify-content:space-between;gap:16px;flex-wrap:wrap;
    font-family:var(--mono);font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}
  h1{font-size:clamp(2rem,6vw,3.2rem);line-height:1.06;margin:.35em 0 .25em;font-weight:700;letter-spacing:-.01em}
  .lede{font-size:1.12rem;max-width:64ch;margin:0 0 16px}
  .lede b{background:linear-gradient(transparent 62%,rgba(180,128,31,.35) 0)}
  .read1{background:var(--card);border:1px solid var(--rule);border-left:6px solid var(--assigned);border-radius:8px;padding:12px 16px;margin:0 0 8px;box-shadow:var(--shadow)}
  .read1 b{color:var(--assigned)}
  h2{font-size:1.4rem;margin:44px 0 4px;display:flex;align-items:baseline;gap:12px}
  h2 .n{font-family:var(--mono);font-size:.78rem;color:var(--muted);border:1px solid var(--rule);border-radius:99px;padding:1px 9px}
  h3{font-size:1.05rem;margin:18px 0 6px}
  .sub{color:var(--muted);margin:0 0 14px;font-size:.98rem}
  .stats{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:14px 0 6px}
  .stat{background:var(--card);border:1px solid var(--rule);border-radius:10px;padding:14px 14px 12px;box-shadow:var(--shadow)}
  .stat b{display:block;font-family:var(--mono);font-size:1.8rem;line-height:1;letter-spacing:-.03em}
  .stat span{display:block;margin-top:7px;font-size:.84rem;color:var(--muted);line-height:1.3}
  @media(max-width:700px){.stats{grid-template-columns:repeat(2,1fr)}}
  .num{font-family:var(--mono)}
  ul.tight{margin:6px 0 10px;padding-left:20px}
  li + li{margin-top:.5em}
  .callout{border-left:4px solid var(--recovered);background:var(--paper2);padding:10px 14px;border-radius:0 8px 8px 0;margin:12px 0;font-size:.95rem}
  .callout.red{border-left-color:var(--assigned)}
  .cav{display:grid;gap:10px}
  details.c{background:var(--card);border:1px solid var(--rule);border-radius:9px;box-shadow:var(--shadow)}
  details.c summary{cursor:pointer;list-style:none;padding:13px 16px;display:flex;gap:14px;align-items:baseline;font-weight:700}
  details.c summary::-webkit-details-marker{display:none}
  details.c summary .k{font-family:var(--mono);color:var(--assigned);font-size:.9rem}
  details.c summary::after{content:"+";margin-left:auto;font-family:var(--mono);font-size:1.15rem;line-height:1;color:var(--muted)}
  details.c[open] summary::after{content:"\\2212"}
  details.c summary:focus-visible{outline:3px solid var(--link);outline-offset:2px;border-radius:9px}
  details.c > div{padding:0 16px 14px 46px;font-size:.97rem}
  .tw{overflow-x:auto;margin:10px 0;border-radius:10px;background:var(--card);box-shadow:var(--shadow)}
  table{border-collapse:collapse;width:100%;font-size:.84rem}
  th,td{padding:7px 10px;text-align:left;vertical-align:top;border:1px solid var(--rule)}
  th{font-family:var(--mono);font-size:.72rem;letter-spacing:.04em;text-transform:uppercase;color:var(--muted);background:var(--paper2);position:sticky;top:0}
  td.nm{min-width:230px;max-width:330px;overflow-wrap:anywhere}
  td.num{white-space:nowrap}
  table.counts td.num,table.counts th{text-align:right}
  table.counts td:first-child,table.counts th:first-child{text-align:left}
  .fig{margin:14px 0;background:var(--card);border:1px solid var(--rule);border-radius:10px;padding:10px;box-shadow:var(--shadow)}
  .fig img{display:block;width:100%;height:auto;border-radius:6px;border:1px solid var(--rule)}
  .fig figcaption{font-size:.9rem;margin-top:8px}
  .links{list-style:none;margin:0 0 8px;padding:0}
  .links li{border:1px solid var(--rule);border-radius:8px;background:var(--paper);padding:10px 14px;display:grid;grid-template-columns:1fr auto;gap:2px 14px;align-items:center}
  .links li:hover{background:var(--paper2)}
  .links .f{font-size:1rem;font-weight:700;color:var(--link)}
  .links .d{font-size:.9rem;color:var(--muted);grid-column:1/3}
  .links .ty{font-family:var(--mono);font-size:.68rem;letter-spacing:.05em;text-transform:uppercase;border:1px solid var(--rule);border-radius:4px;padding:0 6px;background:var(--paper2);color:var(--ink);white-space:nowrap}
"""


def figure(path: Path, alt: str, cap: str) -> str:
    return (f'<figure class="fig"><img src="data:image/png;base64,{b64(path)}" alt="{esc(alt)}">'
            f"<figcaption>{cap}</figcaption></figure>")


egf, ereg, efna4, tyro3, mfap4, frzb_d = (focus_deg[g] for g in FOCUS)
asg_list = ", ".join(ASSIGNED_GENES)
pair_rows = "".join(
    f"<li><b>Rep{k}</b> = control {PAIR[k][0]} with testosterone {PAIR[k][1]} (experiment {k})</li>" for k in (1, 2))

androgen_html = ""
if ANDRO:
    androgen_html = (
        f'<div class="callout">AR (z {ANDRO["AR"][0]:.3f}) and metribolone (z {ANDRO["metribolone"][0]:.3f}) have the '
        "same z-score in both runs, and neither is flagged. Metribolone is a synthetic androgen, so its score partly "
        "reflects the treatment itself.</div>")

body = f"""
<header class="top"><span>Proteomics &middot; testosterone vs control</span><span>Update &middot; 16 Nov 2026</span></header>

<h1>Update on the analysis</h1>
<p class="lede">Following your 11 Nov reply, I kept EFNA4, TYRO3 and MFAP4, dropped the {N_DROPPED} recovered proteins, and redid IPA run 2 on the {n_list}-protein list, flagging results that depend on the five assigned values.</p>

<div class="read1"><b>Read first.</b> Nothing in the main analysis is statistically significant ({FDR_PASS} of {FDR_TOTAL} proteins pass FDR &lt; 0.05), so everything is exploratory. The paired check in section 5 finds 4 leads, not findings. The design is 2 biological replicates per group, each control/testosterone pair from separate experiments (low power). The five assigned values are unconfirmed stand-ins, not measurements.</div>

<div class="stats">
  <div class="stat"><b>{n_list}</b><span>proteins in the list ({n_core} core + FRZB + {n_assigned} assigned)</span></div>
  <div class="stat"><b>{n_gp_terms}</b><span>g:Profiler terms significant, either direction</span></div>
  <div class="stat"><b>{n_nodes}</b><span>STRING proteins, {fmt_int(n_edges)} links</span></div>
  <div class="stat"><b>{FLAG_ALL}</b><span>IPA results flagged as sensitive to the five</span></div>
</div>

<h2><span class="n">1</span> What changed since 11 Nov</h2>
<ul class="tight">
  <li>Extra genes EFNA4, TYRO3 and MFAP4 are confirmed (alongside EGF and EREG, all control-only).</li>
  <li>The {N_DROPPED} recovered proteins are dropped. They are not in any result below.</li>
  <li>The replicate design is stated as 2 biological replicates per group, each control/testosterone pair from separate experiments.</li>
  <li>IPA run 2 uses the {n_list} list with the five assigned. A 716 run without them is the comparison.</li>
</ul>

<h2><span class="n">2</span> The list ({n_list})</h2>
<p class="sub">{n_core} core proteins, FRZB with its measured value, and the five assigned. {n_up} up and {n_down} down in total.</p>
<ul class="tight">
  <li><b>FRZB</b>: measured, log2FC {LFC_FRZB:.5f} (up). Raw p {P_FRZB:.3f}, adjusted p {ADJ_FRZB:.3f}, so not significant.</li>
  <li><b>{esc(asg_list)}</b>: control-only. I assigned each a 100-fold decrease, log2FC {LFC_ASSIGNED}. These are stand-ins, not measured.</li>
</ul>

<h2><span class="n">3</span> g:Profiler and STRING</h2>
<ul class="tight">
  <li><b>g:Profiler</b>: {n_gp_terms} significant terms ({ora['n_query_submitted_up']} up and {ora['n_query_submitted_down']} down proteins queried separately, measured proteome as background).</li>
  <li><b>STRING</b>: {n_nodes} proteins and {fmt_int(n_edges)} links in the network (STRING score 0.4 or higher; the figure below draws only 0.7 or higher). Links per protein: EGF {egf}, EREG {ereg}, EFNA4 {efna4}, TYRO3 {tyro3}, MFAP4 {mfap4}, FRZB {frzb_d}. STRING links are prior knowledge, not evidence from this data.</li>
</ul>
{figure(P_NET, "STRING network around the six focus proteins", f"Look at: how EGF sits in a large cluster while EFNA4 has almost no links. Do not conclude: that a well-connected protein matters more here, since the {len(NET_ASG_DRAWN)} red squares drawn ({', '.join(g.upper() for g in NET_ASG_DRAWN)}) carry assigned values, not measurements; links shown are STRING score 0.7 or higher, so proteins without one are listed in a note, not drawn.")}

<h2><span class="n">4</span> IPA run 2 ({n_list}) vs the 716 run</h2>
<p class="sub">Run 2 includes the five assigned proteins; the 716 run is the same settings without them. A result is <b>flagged</b> if it lists one of the five or its z-score differs between the runs (or it exists in only one). Flagged means sensitive to the assumption, not a finding.</p>
{counts_table()}
<p class="sub">Of the flagged results, {FLIPS_ALL} change sign ({flips['upstream_regulators']} regulators + {flips['diseases_and_functions']} diseases) and {CROSS_ALL} cross |z| = 2 ({cross['canonical_pathways']} pathways + {cross['upstream_regulators']} regulators + {cross['diseases_and_functions']} diseases).</p>
{figure(P_PATH, "Top 15 canonical pathways in IPA run 2", (("Look at: which pathways rank highest and which are hatched (list an assigned protein). " if PATH_N_HAS else "Look at: which pathways rank highest; none of the top 15 lists an assigned protein, so none is hatched. ") + "Do not conclude: that any is a finding; the p-values are inflated and nothing is FDR-significant."))}
<h3>Flagged results (table of {len(FLAG_ROWS)})</h3>
<p class="sub">Rule: all |z| = 2 crossings among canonical pathways and upstream regulators, then the largest |&Delta;z| shifts among pathways and regulators not already included, then a few disease crossings with the largest |&Delta;z|. &ldquo;not in run&rdquo; = no z-score in that run.</p>
{flagged_table()}
{figure(P_SCAT, "z with versus without the five for flagged regulators and diseases", f"Look at: points far from the dashed diagonal, or on opposite sides of a dotted |z| = 2 line. Do not conclude: that these are real effects; " + (f"every one of the {N_SCAT} drawn points lists an assigned protein, " if N_SCAT_ASG == N_SCAT else f"{N_SCAT_ASG} of the {N_SCAT} drawn points list an assigned protein, ") + "and the five values are stand-ins.")}
{androgen_html}
{ps_html}
<h2><span class="n">6</span> Caveats</h2>
<div class="cav">
  <details class="c" open><summary><span class="k">A</span> Nothing is statistically significant</summary>
    <div>{FDR_PASS} of {FDR_TOTAL} proteins pass FDR &lt; 0.05 and g:Profiler has {n_gp_terms} significant terms, so this is exploratory. With 2 biological replicates per group from separate experiments, power is low. (Smallest adjusted p on the list: {min_adj_list:.3f}.)</div></details>
  <details class="c"><summary><span class="k">B</span> The five assigned values are unconfirmed stand-ins</summary>
    <div>{esc(asg_list)} were seen in the controls only. I assigned the 100-fold decrease (log2FC {LFC_ASSIGNED}) so that IPA could accept these proteins. It is not measured and not an estimate of the true change.</div></details>
  <details class="c"><summary><span class="k">C</span> IPA used its default reference set</summary>
    <div>Both IPA runs used the default reference (the whole database), which makes IPA p-values inflated. Compare order and z-scores, not p-values.</div></details>
  <details class="c"><summary><span class="k">D</span> Row caps and different table sizes</summary>
    <div>The Diseases and functions export is capped at {dis_rows} rows in both runs, so a disease can be missing from one run only because of the cap. The upstream tables have {fmt_int(rows721['upstream_regulators'])} rows in run 2 ({n_list}) against {fmt_int(rows716['upstream_regulators'])} in the 716 run, so some regulators exist in only one run.</div></details>
</div>

<h2><span class="n">7</span> Sample pairing</h2>
<p class="sub">I took this from the cleaned data sheet (Pilot Project/CLEANED Silac Proteomics Soluble Factors.xlsx), which lays out the four samples in this order.</p>
<ul class="tight">{pair_rows}</ul>
<p>If the pairs were run differently, please reply with the right ones.</p>

<h2><span class="n">8</span> Files</h2>
{files_html}

"""

html_doc = (f'<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
            f'<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>Proteomics update, 16 Nov 2026</title><style>{CSS}</style></head><body><div class=\"wrap\">"
            f"{body}</div></body></html>\n")
assert "@@" not in html_doc, "unfilled placeholder in output"
OUT.write_text(html_doc, encoding="utf-8")
print("wrote", OUT, f"({len(html_doc) // 1024} KB)")
print("flagged rows:", len(FLAG_ROWS), "| androgen claim:", "kept" if ANDRO else "dropped")
