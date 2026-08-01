"""Tests for the byte-freeze gate itself.

The gate is safety-critical: every later package's acceptance test is "the
frozen outputs did not move." A gate that is silently always-green would let a
real regression through, so these tests exercise BOTH directions -- it must
stay green on legitimate regeneration noise and go red on genuine change.
"""

import shutil
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))
import freeze  # noqa: E402

MANIFEST = REPO_ROOT / "proteomics_de" / "tests" / "expected" / "outputs.sha256"


def test_manifest_exists_and_parses():
    entries = freeze.read_manifest(MANIFEST)
    # Deliberately an explicit number: adding a scientific artifact should be a
    # conscious act that updates this manifest, not something silently absorbed.
    # Grew during the D7 enrichment re-run: the previously-uncached STRING and
    # g:Profiler responses are now recorded under results/enrichment/raw/.
    assert len(entries) == 79, f"expected 79 frozen artifacts, got {len(entries)}"
    modes = {mode for _sha, mode in entries.values()}
    assert modes == {"raw", "svg-canon"}, modes


def test_manifest_freezes_outputs_not_sources():
    """The gate must not freeze source code.

    Freezing .py/.R files makes the gate fail by construction whenever anyone
    refactors -- which is the activity the gate exists to make safe.
    """
    for rel in freeze.read_manifest(MANIFEST):
        assert not rel.endswith((".py", ".R")), f"source file frozen: {rel}"
        assert rel.startswith("proteomics_de/results/") or "_limma_" in rel, rel


def test_tree_matches_manifest():
    ok, changed, missing, extra = freeze.verify(MANIFEST, REPO_ROOT)
    assert not changed, f"frozen outputs drifted: {changed}"
    assert not missing, f"frozen outputs missing: {missing}"
    assert not extra, f"untracked artifacts present: {extra}"
    assert len(ok) == 79


def test_svg_canonicalization_absorbs_regeneration_noise(tmp_path):
    """matplotlib's <dc:date> and salted ids must not count as drift."""
    src = REPO_ROOT / "proteomics_de/results/figures/volcano.svg"
    a = tmp_path / "a.svg"
    shutil.copy(src, a)
    raw = a.read_bytes()

    b = tmp_path / "b.svg"
    b.write_bytes(
        raw.replace(b"<dc:date>", b"<dc:date>2099-01-01T00:00:00.000000", 1)
    )
    assert a.read_bytes() != b.read_bytes(), "fixture did not actually differ"
    assert freeze.digest(a)[0] == freeze.digest(b)[0]


def test_svg_canonicalization_still_catches_real_change(tmp_path):
    """A genuine plot change must change the hash."""
    src = REPO_ROOT / "proteomics_de/results/figures/volcano.svg"
    a = tmp_path / "a.svg"
    shutil.copy(src, a)
    raw = a.read_bytes()

    marker = b'd="M '
    idx = raw.index(marker) + len(marker)
    b = tmp_path / "b.svg"
    b.write_bytes(raw[:idx] + b"9" + raw[idx + 1:])  # perturb one path coordinate
    assert freeze.digest(a)[0] != freeze.digest(b)[0]


def test_salted_ids_stay_distinct(tmp_path):
    """Renumbering must not collapse two distinct ids into one."""
    p = tmp_path / "x.svg"
    p.write_bytes(b'<g id="pabcdef0123"/><g id="p0123abcdef"/>')
    canon, mode = freeze.canonical_bytes(p)
    assert mode == "svg-canon"
    assert b"id0" in canon and b"id1" in canon


def test_raw_mode_is_byte_exact(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("a,b\n1,2\n")
    first = freeze.digest(p)[0]
    p.write_text("a,b\n1,3\n")
    assert freeze.digest(p)[0] != first
