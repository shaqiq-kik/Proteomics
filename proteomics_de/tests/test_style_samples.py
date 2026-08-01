"""`viz/style.py`'s sample maps must be derived from the sample sheet.

Before this layer, ``style.py`` hardcoded ``SAMPLE_COLS`` / ``SAMPLE_LABELS`` /
``SAMPLE_CONDITION`` as literals. ``style`` is imported by all five viz scripts
plus ``gated/pca_cluster.py`` and four ``enrich/`` modules, so those literals
were the widest single copy of the design in the tree.

Three things are pinned here:

1. **No pixels moved.** The derived maps must equal the old literals *exactly*,
   including dict key order -- ``viz/qc_plots.py`` and ``viz/heatmap.py`` iterate
   ``SAMPLE_COLS`` to build DataFrame columns, so order is load-bearing for the
   committed figures' bytes.
2. **The derivation actually derives.** Point the loader at a synthetic sheet and
   the maps must follow it -- both when replicates are *added* (forward path) and
   when the control/treated assignment is *flipped* (DECISIONS_LOG D7). If these
   pass, the D7 flip is a one-line TSV edit rather than a code change.
3. **Import mechanics survive.** ``viz/*.py`` does a bare ``import style``
   relying on ``sys.path[0]``; ``enrich/network_figure.py`` does
   ``from viz.style import ...``. Adding a dependency on the
   ``proteomics_de.config`` package must not break either.
"""

from __future__ import annotations

import hashlib
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPO_ROOT = _PKG_DIR.parent
_VIZ_DIR = _PKG_DIR / "viz"

# ---------------------------------------------------------------------------
# The literals that `viz/style.py` carried before this layer, copied verbatim
# from the pre-refactor file (style.py:120-134 at commit 96d16fd). These are the
# values baked into the 27 committed PNG/SVG figures; nothing may change them.
# ---------------------------------------------------------------------------
LEGACY_SAMPLE_COLS = [
    "Intensity 31578", "Intensity 31580", "Intensity 31579", "Intensity 31581",
]
LEGACY_SAMPLE_LABELS = {
    "Intensity 31578": "control_1",
    "Intensity 31580": "control_2",
    "Intensity 31579": "treated_1",
    "Intensity 31581": "treated_2",
}
LEGACY_SAMPLE_CONDITION = {
    "Intensity 31578": "control",
    "Intensity 31580": "control",
    "Intensity 31579": "treated",
    "Intensity 31581": "treated",
}

#: sha256 of the caveat string as stamped on the committed figures (202 bytes,
#: em dash intact). Compared as a digest so an invisible character swap -- an
#: en dash for the em dash, say -- cannot pass review by eye.
CAVEAT_SHA256 = "02edd2eb390439ebef05126c0551c12de5d52092c1cdb6e5d122fd74e44f5446"


@pytest.fixture(scope="module")
def style():
    """The `style` module, imported the way `viz/*.py` imports it."""
    import style as _style

    return _style


def _write_sheet(path: Path, rows) -> Path:
    """Write a sample sheet TSV; `rows` are (sample, group, channel, replicate)."""
    lines = ["\t".join(["sample", "group", "channel", "replicate"])]
    lines += ["\t".join(str(v) for v in row) for row in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# 1. Derived == today's literals, exactly
# ---------------------------------------------------------------------------
def test_sample_cols_match_legacy_literals(style):
    assert style.SAMPLE_COLS == LEGACY_SAMPLE_COLS


def test_sample_labels_match_legacy_literals(style):
    assert style.SAMPLE_LABELS == LEGACY_SAMPLE_LABELS


def test_sample_condition_match_legacy_literals(style):
    assert style.SAMPLE_CONDITION == LEGACY_SAMPLE_CONDITION


def test_map_key_order_is_the_legacy_order(style):
    """Dict equality ignores order, but the figures do not."""
    assert list(style.SAMPLE_LABELS) == LEGACY_SAMPLE_COLS
    assert list(style.SAMPLE_CONDITION) == LEGACY_SAMPLE_COLS


def test_sample_order_still_aliases_sample_cols(style):
    """`gated/pca_cluster.py` and `qc_plots.py` both rely on this alias."""
    assert style.SAMPLE_ORDER == style.SAMPLE_COLS


def test_maps_agree_with_the_design_package(style):
    """Same facts, reached independently through `config.design`'s own API."""
    from proteomics_de.config import design

    assert style.SAMPLE_COLS == design.sample_columns()
    assert list(style.SAMPLE_CONDITION.values()) == design.group_vector()
    assert set(style.SAMPLE_CONDITION.values()) == set(design.group_names())


def test_no_sample_id_is_hardcoded_in_style_source():
    """D7 flips control/treated by editing the TSV; style.py must not pin ids.

    A literal `Intensity 31578` in the source would silently survive the flip.
    Sample ids may still appear in prose (comments/docstrings), so only code
    lines are inspected.
    """
    src = (_VIZ_DIR / "style.py").read_text(encoding="utf-8")
    offenders = []
    for lineno, line in enumerate(src.splitlines(), start=1):
        code = line.split("#", 1)[0]
        if any(sample_id in code for sample_id in ("31578", "31579", "31580", "31581")):
            offenders.append(f"{lineno}: {line.strip()}")
    assert not offenders, "sample ids hardcoded in style.py code:\n" + "\n".join(offenders)


# ---------------------------------------------------------------------------
# 2. CAVEAT_TEXT is re-exported, not re-typed
# ---------------------------------------------------------------------------
def test_caveat_text_is_the_constants_object(style):
    from proteomics_de.config import constants

    assert style.CAVEAT_TEXT is constants.CAVEAT_TEXT


def test_caveat_text_digest_is_unchanged(style):
    digest = hashlib.sha256(style.CAVEAT_TEXT.encode("utf-8")).hexdigest()
    assert digest == CAVEAT_SHA256
    assert len(style.CAVEAT_TEXT.encode("utf-8")) == 202
    assert "—" in style.CAVEAT_TEXT  # em dash, not en dash or hyphen


def test_add_caveat_still_defaults_to_the_caveat(style):
    """10+ call sites use the default arg; `heatmap.py` reads the name directly."""
    import inspect

    default = inspect.signature(style.add_caveat).parameters["text"].default
    assert default is style.CAVEAT_TEXT


# ---------------------------------------------------------------------------
# 3. Forward path -- more replicates
# ---------------------------------------------------------------------------
def test_derivation_grows_with_a_six_sample_sheet(style, tmp_path):
    sheet = _write_sheet(
        tmp_path / "sample_sheet.tsv",
        [
            ("31578", "control", "Intensity 31578", 1),
            ("31580", "control", "Intensity 31580", 2),
            ("31582", "control", "Intensity 31582", 3),
            ("31579", "treated", "Intensity 31579", 1),
            ("31581", "treated", "Intensity 31581", 2),
            ("31583", "treated", "Intensity 31583", 3),
        ],
    )
    cols, labels, condition = style.derive_sample_maps(sheet)

    assert len(cols) == len(labels) == len(condition) == 6
    assert cols == [
        "Intensity 31578", "Intensity 31580", "Intensity 31582",
        "Intensity 31579", "Intensity 31581", "Intensity 31583",
    ]
    assert labels == {
        "Intensity 31578": "control_1",
        "Intensity 31580": "control_2",
        "Intensity 31582": "control_3",
        "Intensity 31579": "treated_1",
        "Intensity 31581": "treated_2",
        "Intensity 31583": "treated_3",
    }
    assert list(condition.values()) == ["control"] * 3 + ["treated"] * 3


def test_row_order_in_the_tsv_does_not_permute_the_maps(style, tmp_path):
    """Canonical order comes from the `group` column, not from file row order."""
    shuffled = _write_sheet(
        tmp_path / "shuffled.tsv",
        [
            ("31581", "treated", "Intensity 31581", 2),
            ("31578", "control", "Intensity 31578", 1),
            ("31579", "treated", "Intensity 31579", 1),
            ("31580", "control", "Intensity 31580", 2),
        ],
    )
    cols, labels, condition = style.derive_sample_maps(shuffled)

    assert cols == LEGACY_SAMPLE_COLS
    assert labels == LEGACY_SAMPLE_LABELS
    assert condition == LEGACY_SAMPLE_CONDITION


# ---------------------------------------------------------------------------
# 4. The D7 flip -- control/treated inverted
# ---------------------------------------------------------------------------
def test_d7_flip_propagates_through_the_derivation(style, tmp_path):
    """DECISIONS_LOG D7: 31579/31581 are the controls (Vehicle), 31578/31580 the
    treated (Testosterone). The flip must reach the figures by editing the TSV
    alone. D7 also notes the replicate PAIRING is unaffected: Rep1 =
    31578/31579, Rep2 = 31580/31581.
    """
    flipped = _write_sheet(
        tmp_path / "flipped.tsv",
        [
            ("31578", "treated", "Intensity 31578", 1),
            ("31580", "treated", "Intensity 31580", 2),
            ("31579", "control", "Intensity 31579", 1),
            ("31581", "control", "Intensity 31581", 2),
        ],
    )
    cols, labels, condition = style.derive_sample_maps(flipped)

    # Condition follows the sheet, not the sample id.
    assert condition == {
        "Intensity 31579": "control",
        "Intensity 31581": "control",
        "Intensity 31578": "treated",
        "Intensity 31580": "treated",
    }
    # Controls still sort first, so the column order inverts with the labels.
    assert cols == [
        "Intensity 31579", "Intensity 31581", "Intensity 31578", "Intensity 31580",
    ]
    assert labels == {
        "Intensity 31579": "control_1",
        "Intensity 31581": "control_2",
        "Intensity 31578": "treated_1",
        "Intensity 31580": "treated_2",
    }
    # Every id that is control today is treated after the flip, and vice versa.
    for channel, group in condition.items():
        assert group != LEGACY_SAMPLE_CONDITION[channel]


def test_condition_colors_cover_every_derived_group(style, tmp_path):
    """`qc_plots.py:201` does CONDITION_COLORS[SAMPLE_CONDITION[c]] -- a group
    name the palette lacks would be a KeyError at figure time, including after
    the D7 flip."""
    flipped = _write_sheet(
        tmp_path / "flipped.tsv",
        [
            ("31578", "treated", "Intensity 31578", 1),
            ("31580", "treated", "Intensity 31580", 2),
            ("31579", "control", "Intensity 31579", 1),
            ("31581", "control", "Intensity 31581", 2),
        ],
    )
    _cols, _labels, condition = style.derive_sample_maps(flipped)

    for group in set(condition.values()) | set(style.SAMPLE_CONDITION.values()):
        assert group in style.CONDITION_COLORS


# ---------------------------------------------------------------------------
# 5. Import mechanics: bare `import style`, from any cwd
# ---------------------------------------------------------------------------
def test_import_style_works_under_pytest(style):
    """Trivially true if we got here -- the conftest sys.path layout resolves."""
    assert style.__file__ == str(_VIZ_DIR / "style.py")


@pytest.mark.parametrize(
    "cwd, label",
    [
        (str(_VIZ_DIR), "viz dir (sys.path[0], as viz scripts run)"),
        (str(_REPO_ROOT), "repo root"),
        (str(Path(_REPO_ROOT).anchor), "filesystem root (unrelated cwd)"),
    ],
)
def test_bare_import_style_from_various_cwds(cwd, label):
    """`python <script>` sets sys.path[0] to the script's dir, so a bare
    `import style` must still find `proteomics_de.config` from any cwd."""
    program = textwrap.dedent(
        f"""
        import sys
        sys.path.insert(0, {str(_VIZ_DIR)!r})   # what `python viz/volcano.py` does
        import style
        assert style.SAMPLE_COLS == {LEGACY_SAMPLE_COLS!r}, style.SAMPLE_COLS
        assert style.SAMPLE_LABELS == {LEGACY_SAMPLE_LABELS!r}, style.SAMPLE_LABELS
        assert style.SAMPLE_CONDITION == {LEGACY_SAMPLE_CONDITION!r}, style.SAMPLE_CONDITION
        assert len(style.CAVEAT_TEXT.encode()) == 202
        print("OK")
        """
    )
    proc = subprocess.run(
        [sys.executable, "-c", program],
        cwd=cwd, capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"bare `import style` failed from {label}:\n{proc.stderr}"
    assert "OK" in proc.stdout


def test_qualified_viz_style_import_still_works():
    """`enrich/network_figure.py` does `from viz.style import CAVEAT_TEXT, ...`."""
    program = textwrap.dedent(
        f"""
        import sys
        sys.path.insert(0, {str(_PKG_DIR)!r})    # what network_figure.py adds
        sys.path.insert(0, {str(_REPO_ROOT)!r})
        from viz.style import CAVEAT_TEXT, add_caveat, SAMPLE_COLS
        assert len(CAVEAT_TEXT.encode()) == 202
        assert SAMPLE_COLS == {LEGACY_SAMPLE_COLS!r}
        print("OK")
        """
    )
    proc = subprocess.run(
        [sys.executable, "-c", program],
        cwd=str(Path(_REPO_ROOT).anchor), capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"`from viz.style import ...` failed:\n{proc.stderr}"
    assert "OK" in proc.stdout


def test_style_does_not_leak_a_duplicate_config_module(style):
    """The bootstrap must reuse the package, not load a second copy of it.

    Two module objects would mean two `CAVEAT_TEXT` strings, and the identity
    re-export in `test_caveat_text_is_the_constants_object` would be a lie.
    """
    from proteomics_de.config import constants, design

    assert style._constants is constants
    assert style._design is design
