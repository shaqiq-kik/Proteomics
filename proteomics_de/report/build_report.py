#!/usr/bin/env python3
"""
build_report.py — reproducible assembler for the final interactive report.

Reads the hand-authored template `report_template.html` and produces a single,
fully self-contained `report.html`:

  * embeds the two display/body fonts (Fraunces + Newsreader, roman + italic
    variable woff2) as base64 `@font-face` sources, and
  * embeds all 13 figures as `data:` URIs — SVG (crisp vector) for every figure
    except ppi_network + rank_abundance, which are embedded as base64 PNG per
    the build spec (the ppi_network SVG is ~1.5 MB / 25k nodes).

Figure titles and captions are pulled verbatim from the verified digest
`report_facts.json` (the content source of truth) so no numbers are invented
here. The template carries `%%TOKEN%%` placeholders; this script substitutes
them and writes the result. It never modifies any pipeline script, result, or
figure — it only READS the committed figures to embed copies.

Usage:
    python3 build_report.py            # writes report.html next to the template
"""
from __future__ import annotations

import base64
import html
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent                       # proteomics_de/report
FIGDIR = HERE.parent / "results" / "figures"                 # proteomics_de/results/figures
FONTDIR = HERE / "assets" / "fonts"
TEMPLATE = HERE / "report_template.html"
OUTPUT = HERE / "report.html"

# report_facts.json (verified digest) lives in the session scratchpad; allow an
# override argument, and fall back to a copy beside this script if present.
FACTS_CANDIDATES = [
    Path(
        "/private/tmp/claude-501/-Users-abrarshakik-Documents-Proteomics/"
        "258a4089-6261-45f0-9686-91944327ce50/scratchpad/report_facts.json"
    ),
    HERE / "report_facts.json",
]


def load_facts() -> dict:
    for c in FACTS_CANDIDATES:
        if c.exists():
            return json.loads(c.read_text())
    raise SystemExit("report_facts.json not found; pass its path as argv[1].")


def b64(path: Path) -> str:
    return base64.b64encode(path.read_bytes()).decode("ascii")


# ---- font tokens ------------------------------------------------------------
FONTS = {
    "%%FONT_FRAUNCES_ROMAN%%": "fraunces-roman.woff2",
    "%%FONT_FRAUNCES_ITALIC%%": "fraunces-italic.woff2",
    "%%FONT_NEWSREADER_ROMAN%%": "newsreader-roman.woff2",
    "%%FONT_NEWSREADER_ITALIC%%": "newsreader-italic.woff2",
}

# ---- figure plate metadata --------------------------------------------------
# key -> (basename-without-ext, embed-format, figure-number, tag, hypothesis?)
FIGURES = [
    ("intensity_distributions", "svg", 1,  "QC · distributions",     False),
    ("missing_values",          "svg", 2,  "QC · missingness",       False),
    ("sample_correlation",      "svg", 3,  "QC · reproducibility",   False),
    ("rank_abundance",          "png", 4,  "QC · dynamic range",     False),
    ("volcano",                 "svg", 5,  "Track A · limma",        True),
    ("ma_plot",                 "svg", 6,  "Track A · limma",        True),
    ("ora_dotplot",             "svg", 7,  "Track A · enrichment",   True),
    ("pca_qc",                  "svg", 8,  "QC only · n = 4",        False),
    ("sample_dendrogram",       "svg", 9,  "QC only · n = 4",        False),
    ("heatmap_top_de",          "svg", 10, "Track B · leads",        True),
    ("ppi_network",             "png", 11, "Track B · leads",        True),
    ("gsea_top",                "svg", 12, "Track B · leads",        True),
    ("upset",                   "svg", 13, "Track B · leads",        True),
]

MIME = {"svg": "image/svg+xml", "png": "image/png"}


def facts_by_key(facts: dict) -> dict:
    out = {}
    for f in facts["figures"]:
        key = Path(f["file"]).stem            # "figures/volcano.png" -> "volcano"
        out[key] = f
    return out


def make_plate(key, fmt, num, tag, hypo, meta) -> str:
    src_path = FIGDIR / f"{key}.{fmt}"
    if not src_path.exists():
        raise SystemExit(f"missing figure file: {src_path}")
    data_uri = f"data:{MIME[fmt]};base64,{b64(src_path)}"
    title = html.escape(meta["title"])
    caption = html.escape(meta["caption"])
    warn = '<span class="fw">hypothesis-generating</span>' if hypo else ""
    return (
        '<div class="plate">\n'
        '  <figure>\n'
        f'    <div class="imgwrap"><img loading="lazy" decoding="async" alt="{title}" src="{data_uri}"></div>\n'
        '    <figcaption>\n'
        f'      <div class="ftag"><span>Figure&nbsp;{num:02d}</span><span>{html.escape(tag)}</span>{warn}</div>\n'
        f'      <div class="fttl">{title}</div>\n'
        f'      <p>{caption}</p>\n'
        '    </figcaption>\n'
        '  </figure>\n'
        '</div>'
    )


def main() -> None:
    if len(sys.argv) > 1:
        FACTS_CANDIDATES.insert(0, Path(sys.argv[1]))
    facts = load_facts()
    fmeta = facts_by_key(facts)

    html_text = TEMPLATE.read_text()

    # 1) fonts
    for token, fname in FONTS.items():
        fp = FONTDIR / fname
        if not fp.exists():
            raise SystemExit(f"missing font file: {fp}")
        html_text = html_text.replace(token, f"data:font/woff2;base64,{b64(fp)}")

    # 2) figure plates
    embedded = 0
    for key, fmt, num, tag, hypo in FIGURES:
        token = f"%%PLATE_{key}%%"
        if token not in html_text:
            raise SystemExit(f"template is missing placeholder {token}")
        html_text = html_text.replace(token, make_plate(key, fmt, num, tag, hypo, fmeta[key]))
        embedded += 1

    # 3) sanity: no unresolved tokens remain
    leftover = re.findall(r"%%[A-Z0-9_]+%%", html_text)
    if leftover:
        raise SystemExit(f"unresolved template tokens remain: {sorted(set(leftover))}")

    OUTPUT.write_text(html_text)

    # 4) self-check report
    size = OUTPUT.stat().st_size
    external = [
        u for u in re.findall(r'(?:src|href)\s*=\s*"([^"]+)"', html_text)
        if u.startswith("http://") or u.startswith("https://")
    ]
    print(f"wrote {OUTPUT}")
    print(f"  size: {size/1024/1024:.2f} MB ({size:,} bytes)")
    print(f"  figures embedded: {embedded}/13")
    print(f"  external src/href references (should be none): {len(external)}")
    for u in external:
        print(f"    ! {u}")


if __name__ == "__main__":
    main()
