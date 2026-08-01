"""
Shared visual system for the proteomics DE report figures (Section 6, items 11-14).

Palette, rcParams, and small helpers so every figure in `proteomics_de/viz/`
reads as one deliberate system rather than default matplotlib styling.

Palette source: the `dataviz` skill's validated reference palette
(references/palette.md). Colors are used exactly as documented there -
no eyeballed hex values.

Color jobs used in this report:
  - "regulated" (UP / DOWN / NO CHANGE) is treated as a DIVERGING / polarity
    encoding (up vs. down relative to a no-change baseline), per the skill's
    color-formula: two hues (blue <-> red) + a neutral gray midpoint.
    UP = red (categorical slot 6, #e34948), DOWN = blue (categorical slot 1,
    #2a78d6), NO CHANGE = neutral baseline gray, drawn recessively.
  - ON_OFF (present/absent between conditions) is a distinct nominal category,
    not part of the up/down polarity -> categorical slot 8 (orange).
  - Sample-level identity (control vs. treated) reuses the same blue/red pair
    so the whole report shares one color language: control = blue, treated = red.
  - Continuous magnitude (correlations, % present, rank-abundance) uses the
    documented single-hue sequential blue ramp (steps 100 -> 700).
  - The top-DE z-score heatmap uses the diverging blue/neutral/red triple as a
    3-stop colormap (documented endpoints only).
"""

import json
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------------------
# Deterministic behavior
# ----------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
VIZ_DIR = Path(__file__).resolve().parent
PROTEOMICS_DE_DIR = VIZ_DIR.parent
REPO_ROOT = PROTEOMICS_DE_DIR.parent
RESULTS_DIR = PROTEOMICS_DE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
MANIFEST_PATH = FIGURES_DIR / "figures_manifest.json"

# ----------------------------------------------------------------------------
# config package bootstrap
# ----------------------------------------------------------------------------
# `viz/*.py` scripts are run as file paths (`python proteomics_de/viz/volcano.py`),
# so sys.path[0] is *this* directory and the repo root is not importable. Tests,
# by contrast, get the repo root from tests/conftest.py. Try the normal import
# first (so pytest / `python -m` reuse the already-imported package object rather
# than a second copy) and only fall back to extending sys.path.
#
# The path is derived from `__file__`, never from the cwd, so this works from any
# working directory. It is appended, not prepended, so it can never shadow a
# module the caller already resolves.
try:  # pragma: no cover - exercised by both branches across entry points
    from proteomics_de.config import constants as _constants
    from proteomics_de.config import design as _design
except ImportError:
    if str(REPO_ROOT) not in sys.path:
        sys.path.append(str(REPO_ROOT))
    from proteomics_de.config import constants as _constants
    from proteomics_de.config import design as _design

# ----------------------------------------------------------------------------
# Palette (verbatim from the dataviz skill's references/palette.md)
# ----------------------------------------------------------------------------
CATEGORICAL = {
    "blue": "#2a78d6",
    "aqua": "#1baf7a",
    "yellow": "#eda100",
    "green": "#008300",
    "violet": "#4a3aa7",
    "red": "#e34948",
    "magenta": "#e87ba4",
    "orange": "#eb6834",
}

SEQUENTIAL_BLUE = {
    100: "#cde2fb", 150: "#b7d3f6", 200: "#9ec5f4", 250: "#86b6ef",
    300: "#6da7ec", 350: "#5598e7", 400: "#3987e5", 450: "#2a78d6",
    500: "#256abf", 550: "#1c5cab", 600: "#184f95", 650: "#104281",
    700: "#0d366b",
}

DIVERGING = {
    "neg": CATEGORICAL["blue"],     # #2a78d6
    "neg_dark": SEQUENTIAL_BLUE[700],  # #0d366b
    "mid": "#f0efec",               # neutral gray midpoint (light surface)
    "pos": CATEGORICAL["red"],      # #e34948
}

STATUS = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "serious": "#ec835a",
    "critical": "#d03b3b",
}

CHROME = {
    "surface": "#fcfcfb",
    "page": "#f9f9f7",
    "ink_primary": "#0b0b0b",
    "ink_secondary": "#52514e",
    "ink_muted": "#898781",
    "gridline": "#e1e0d9",
    "baseline": "#c3c2b7",
}

# ----------------------------------------------------------------------------
# Semantic color roles for this report
# ----------------------------------------------------------------------------
REGULATED_COLORS = {
    "UP": CATEGORICAL["red"],
    "DOWN": CATEGORICAL["blue"],
    "NO CHANGE": CHROME["baseline"],
    "ON_OFF": CATEGORICAL["orange"],
}
REGULATED_ALPHA = {
    "UP": 0.85,
    "DOWN": 0.85,
    "NO CHANGE": 0.45,
    "ON_OFF": 0.85,
}
REGULATED_ORDER = ["NO CHANGE", "DOWN", "UP", "ON_OFF"]  # draw order (background first)
REGULATED_LEGEND_ORDER = ["UP", "DOWN", "NO CHANGE", "ON_OFF"]

CONDITION_COLORS = {"control": CATEGORICAL["blue"], "treated": CATEGORICAL["red"]}

def derive_sample_maps(sheet=None):
    """Derive the three sample maps from ``config/sample_sheet.tsv``.

    Returns ``(cols, labels, condition)`` where

    * ``cols`` is the intensity column names in canonical order (control group
      first, then treated; ascending replicate within a group),
    * ``labels`` maps each intensity column to its display label
      ``f"{group}_{replicate}"`` (e.g. ``"control_1"``),
    * ``condition`` maps each intensity column to its group name.

    Both dicts are keyed in canonical order, so ``list(labels)`` == ``cols``.

    Everything is read from the sheet, nothing is assumed. In particular no
    sample id is hardwired to a group: DECISIONS_LOG D7 records that the
    control/treated assignment is inverted and will be flipped by editing the
    TSV, and that flip propagates through here (and thus through every figure)
    with no code change. Adding biological replicates likewise only means adding
    rows -- the maps grow to match.

    Parameters
    ----------
    sheet :
        Optional sheet / path override forwarded to ``config.design``. Defaults
        to the committed ``config/sample_sheet.tsv``.
    """
    df = _design.read_sample_sheet(sheet) if sheet is not None else _design.read_sample_sheet()
    cols = df["channel"].tolist()
    groups = df["group"].tolist()
    labels = {
        channel: f"{group}_{replicate}"
        for channel, group, replicate in zip(cols, groups, df["replicate"].tolist())
    }
    condition = dict(zip(cols, groups))
    return cols, labels, condition


#: For the committed (D7-corrected) sheet: ``["Intensity 31579", "Intensity
#: 31581", "Intensity 31578", "Intensity 31580"]`` with labels
#: control_1/control_2/treated_1/treated_2. Before D7 this module hardcoded the
#: literal ``["Intensity 31578", "Intensity 31580", "Intensity 31579",
#: "Intensity 31581"]`` with control_1/2 on the FIRST pair -- i.e. the
#: now-inverted assignment; deriving from the sheet is what let the D7 fix
#: reach every figure with no code edit here. `tests/test_style_samples.py`
#: asserts the current values.
SAMPLE_COLS, SAMPLE_LABELS, SAMPLE_CONDITION = derive_sample_maps()
SAMPLE_ORDER = SAMPLE_COLS  # canonical order: control group first, then treated

#: Re-exported so there is one source of truth. This exact string is stamped
#: onto 10+ committed figures via `add_caveat()` (and directly by
#: `viz/heatmap.py`), so the name must keep resolving here.
CAVEAT_TEXT = _constants.CAVEAT_TEXT

FC_THRESHOLD = 0.585  # |log2FC| >= 0.585 is the up/down call boundary used upstream
RAW_P_THRESHOLD = 0.05

# ----------------------------------------------------------------------------
# Diverging colormap for the heatmap (documented stops only: blue -> neutral -> red)
# ----------------------------------------------------------------------------
DIVERGING_CMAP = LinearSegmentedColormap.from_list(
    "report_diverging",
    [DIVERGING["neg"], DIVERGING["mid"], DIVERGING["pos"]],
    N=256,
)
# Exactly the documented diverging pair (blue <-> red) + documented neutral
# midpoint gray, verbatim from palette.md - no invented/eyeballed hex values.

SEQUENTIAL_CMAP = LinearSegmentedColormap.from_list(
    "report_sequential_blue",
    [SEQUENTIAL_BLUE[k] for k in sorted(SEQUENTIAL_BLUE)],
    N=256,
)

# ----------------------------------------------------------------------------
# rcParams / style application
# ----------------------------------------------------------------------------
def apply_style():
    plt.rcParams.update({
        "figure.facecolor": CHROME["surface"],
        "axes.facecolor": CHROME["surface"],
        "savefig.facecolor": CHROME["surface"],
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "text.color": CHROME["ink_primary"],
        "axes.edgecolor": CHROME["baseline"],
        "axes.labelcolor": CHROME["ink_secondary"],
        "axes.titlecolor": CHROME["ink_primary"],
        "axes.titlesize": 13,
        "axes.titleweight": "semibold",
        "axes.labelsize": 10.5,
        "xtick.color": CHROME["ink_muted"],
        "ytick.color": CHROME["ink_muted"],
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "grid.color": CHROME["gridline"],
        "grid.linewidth": 0.8,
        "grid.linestyle": "-",
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "figure.dpi": 100,
        "svg.fonttype": "none",
    })


apply_style()


def recede_spines(ax):
    """Thin, hairline, recessive axis spines (skill: 'recessive grid/axes')."""
    for side in ("left", "bottom"):
        ax.spines[side].set_color(CHROME["baseline"])
        ax.spines[side].set_linewidth(0.8)


def add_caveat(fig, text=CAVEAT_TEXT, y=0.005):
    """Small italic caveat footer, present on every statistical figure."""
    fig.text(
        0.5, y, text, ha="center", va="bottom",
        fontsize=7.8, style="italic", color=CHROME["ink_muted"], wrap=True,
    )


def save_fig(fig, stem, dpi=200, tight=True):
    """Save a matplotlib Figure as both PNG (dpi=200) and SVG under figures/.

    Returns the two paths written, relative to `RESULTS_DIR` (i.e. what
    the manifest records as `file`).
    """
    png_path = FIGURES_DIR / f"{stem}.png"
    svg_path = FIGURES_DIR / f"{stem}.svg"
    kwargs = {"facecolor": CHROME["surface"]}
    if tight:
        kwargs["bbox_inches"] = "tight"
    fig.savefig(png_path, dpi=dpi, **kwargs)
    fig.savefig(svg_path, **kwargs)
    plt.close(fig)
    return (
        str(png_path.relative_to(RESULTS_DIR)),
        str(svg_path.relative_to(RESULTS_DIR)),
    )


def record_manifest(entries):
    """Merge a list of manifest entry dicts into figures_manifest.json, keyed by `file`.

    Each entry: {file, title, caption, key_numbers}.
    Safe to call once per script; re-running a script overwrites only its own entries.
    """
    existing = []
    if MANIFEST_PATH.exists():
        try:
            existing = json.loads(MANIFEST_PATH.read_text())
        except json.JSONDecodeError:
            existing = []
    by_file = {e["file"]: e for e in existing}
    for e in entries:
        by_file[e["file"]] = e
    merged = list(by_file.values())
    MANIFEST_PATH.write_text(json.dumps(merged, indent=2, default=_json_default))
    return merged


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def safe_log2(series):
    """log2 transform treating values <= 0 (and NaN) as missing."""
    arr = series.astype(float)
    arr = arr.where(arr > 0)
    return np.log2(arr)


def fmt_p(p):
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.4f}"
