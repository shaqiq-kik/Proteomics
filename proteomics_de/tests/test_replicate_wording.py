"""D23 guard: the live project must not call the samples technical replicates.

On 2026-11-11 the investigator, Dr. Walker, stated that the four samples are
2 BIOLOGICAL replicates per group from 2 SEPARATE experiments (each experiment
gave one control and one treated sample). Until then the whole repository said
the opposite -- "n = 2 technical replicates per group, one biological replicate"
-- in the config, the figure footers, every provenance sidecar, the report and the
READMEs. That claim is now false, and a false claim about the study design is
exactly the kind of text that comes back when someone copies an old docstring.

This test scans every live text file for the retired phrases and fails if any
reappears. It deliberately does NOT scan:

* ``DECISIONS_LOG.md`` / ``BUILD_LOG.md`` -- history is append-only, and D23 itself
  must quote what it supersedes;
* ``research1.md``, ``Pilot Project/``, ``report_11nov.html``, ``for_dr_walker.html``
  -- historical documents (the 11 Nov page was already emailed);
* ``results/**/raw/`` caches and the committed manual IPA exports ``run1_*`` /
  ``run2_*``;
* ``tests/`` -- the tests have to name the phrases in order to forbid them.
"""
from __future__ import annotations

import re
from pathlib import Path

import yaml

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PKG = _REPO_ROOT / "proteomics_de"

#: Retired claims. Matched case-insensitively with whitespace/NBSP collapsed.
RETIRED_PHRASES = (
    "technical replicates only",
    "one biological replicate",
    "single biological replicate",
    "n = 2 technical",
    "n=2 technical",
    "technical replicates per group",
    "technical replicates per condition",
    "no biological replication",
    "technical-replicates-only",
)

_TEXT_SUFFIXES = {
    ".py", ".md", ".yaml", ".yml", ".json", ".html", ".csv", ".tsv", ".txt",
    ".svg", ".toml", ".sh", ".r",
}
_SKIP_DIR_NAMES = {
    ".git", ".venv", "venv", "__pycache__", ".pytest_cache", "node_modules",
    "Pilot Project", "tests", "raw",
}
_SKIP_FILE_NAMES = {
    "DECISIONS_LOG.md", "BUILD_LOG.md", "research1.md",
    "report_11nov.html", "for_dr_walker.html",
}
#: Manual IPA exports (committed, historical).
_SKIP_FILE_PREFIXES = ("run1_", "run2_")


def _live_text_files():
    for path in _REPO_ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in _TEXT_SUFFIXES:
            continue
        rel_parts = path.relative_to(_REPO_ROOT).parts
        if any(part in _SKIP_DIR_NAMES for part in rel_parts[:-1]):
            continue
        if path.name in _SKIP_FILE_NAMES or path.name.startswith(_SKIP_FILE_PREFIXES):
            continue
        yield path


def _normalise(text: str) -> str:
    text = text.replace("&nbsp;", " ").replace(" ", " ")
    return re.sub(r"\s+", " ", text).lower()


def test_the_scan_actually_covers_the_live_files():
    """Guard the guard: a scan that finds nothing to read would always pass."""
    rels = {p.relative_to(_REPO_ROOT).as_posix() for p in _live_text_files()}
    for must in (
        "README.md",
        "proteomics_de/README.md",
        "proteomics_de/STATUS.md",
        "proteomics_de/config/config.yaml",
        "proteomics_de/config/constants.py",
        "proteomics_de/provenance.py",
        "proteomics_de/report/report.html",
        "proteomics_de/report/report_template.html",
        "proteomics_de/report/report_facts.json",
        "proteomics_de/results/enrichment_walker/README.md",
        "proteomics_de/results/ipa_analysis_walker/README.md",
        "proteomics_de/results/ipa_analysis_walker/UPLOAD_STEPS.md",
        "proteomics_de/results/foldchange_all.csv.provenance.json",
    ):
        assert must in rels, f"{must} is not being scanned"
    # ...and the historical files really are exempt.
    for exempt in ("proteomics_de/DECISIONS_LOG.md", "research1.md", "report_11nov.html"):
        assert exempt not in rels


def test_no_live_file_calls_the_samples_technical_replicates():
    hits = []
    for path in _live_text_files():
        text = _normalise(path.read_text(encoding="utf-8", errors="replace"))
        for phrase in RETIRED_PHRASES:
            if phrase in text:
                hits.append(f"{path.relative_to(_REPO_ROOT).as_posix()}: {phrase!r}")
    assert not hits, (
        "retired replicate wording is back (D23: 2 biological replicates from 2 "
        "separate experiments):\n  " + "\n  ".join(hits)
    )


def test_the_corrected_design_is_stated_in_the_single_source_files():
    from proteomics_de.config.constants import CAVEAT_TEXT

    assert "n = 2 biological replicates per group, from 2 separate experiments" in CAVEAT_TEXT
    assert "low power" in CAVEAT_TEXT
    assert "exploratory / hypothesis-generating" in CAVEAT_TEXT

    cfg = yaml.safe_load((_PKG / "config" / "config.yaml").read_text(encoding="utf-8"))
    design = cfg["design"]
    assert design["n_biological_replicates_per_group"] == 2
    # Not stated by the investigator; each channel is treated as one biological sample.
    assert design["n_technical_replicates_per_group"] == 1
    assert design["pairing_confirmed"] is False
    caveat = _normalise(cfg["notes"]["caveat"])
    assert "2 biological replicates per group" in caveat
    assert "2 separate experiments" in caveat
    assert "low statistical power" in caveat
    # the model is unchanged: a paired / blocked fit is an open question, not a change
    assert cfg["statistics"]["design_matrix"].startswith("~group")
