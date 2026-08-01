"""Tests for report/build_report.py and the generated report/report.html.

Two latent bugs are pinned here:

1. ``FACTS_CANDIDATES`` used to check a hardcoded absolute scratchpad path
   *before* the committed ``report_facts.json``. It worked only because that
   path happened not to exist -- a later session landing on the same UUID would
   have built the report from a foreign facts file, silently.
2. The unresolved-token check used ``%%[A-Z0-9_]+%%``, which cannot match a
   lowercase token. Every real token is upper-case, so the hole was invisible;
   a mistyped ``%%plate_volcano%%`` would have shipped as literal text.
"""

from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPORT_DIR = _PKG_DIR / "report"
_REPORT_HTML = _REPORT_DIR / "report.html"


def _load_build_report():
    spec = importlib.util.spec_from_file_location(
        "build_report", _REPORT_DIR / "build_report.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def build_report():
    return _load_build_report()


@pytest.fixture(scope="module")
def report_html() -> str:
    return _REPORT_HTML.read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# bug 1: the dead absolute scratchpad path
# ---------------------------------------------------------------------------
def test_facts_candidates_contains_only_the_committed_file(build_report):
    assert build_report.FACTS_CANDIDATES == [_REPORT_DIR / "report_facts.json"]


def test_no_scratchpad_path_is_consulted(build_report):
    """A session-scoped temp path must never outrank the committed digest."""
    for candidate in build_report.FACTS_CANDIDATES:
        text = str(candidate)
        assert "/private/tmp" not in text
        assert "scratchpad" not in text
        assert "claude-501" not in text


def test_build_report_source_has_no_hardcoded_absolute_temp_path():
    source = (_REPORT_DIR / "build_report.py").read_text(encoding="utf-8")
    # The comment explaining the removal may name it; the code must not use it.
    code = "\n".join(
        line for line in source.splitlines()
        if not line.lstrip().startswith("#")
    )
    assert "258a4089" not in code
    assert "/private/tmp/claude-501" not in code


def test_load_facts_reads_the_committed_digest(build_report):
    facts = build_report.load_facts()
    assert facts["n_figures"] == 13
    assert facts["counts"]["regulated_UP"] == 509


# ---------------------------------------------------------------------------
# bug 2: the case-blind unresolved-token check
# ---------------------------------------------------------------------------
def test_leftover_token_pattern_matches_lowercase():
    """This is the bug: [A-Z0-9_] silently ignored lowercase tokens."""
    source = (_REPORT_DIR / "build_report.py").read_text(encoding="utf-8")
    match = re.search(r're\.findall\(r"(%%\[[^\]]+\]\+%%)"', source)
    assert match, "could not locate the unresolved-token regex in build_report.py"
    pattern = match.group(1)

    assert re.findall(pattern, "%%PLATE_VOLCANO%%") == ["%%PLATE_VOLCANO%%"]
    assert re.findall(pattern, "%%plate_volcano%%") == ["%%plate_volcano%%"]
    assert re.findall(pattern, "%%Plate_Volcano%%") == ["%%Plate_Volcano%%"]
    assert re.findall(pattern, "%%font_fraunces_roman%%") == ["%%font_fraunces_roman%%"]


def test_old_uppercase_only_pattern_would_have_missed_it():
    """Pin the failure mode so the fix cannot be quietly reverted."""
    assert re.findall(r"%%[A-Z0-9_]+%%", "%%plate_volcano%%") == []
    assert re.findall(r"%%[A-Za-z0-9_]+%%", "%%plate_volcano%%") == ["%%plate_volcano%%"]


def test_builder_rejects_a_lowercase_unresolved_token(build_report, tmp_path, monkeypatch):
    """End-to-end: a lowercase leftover must abort the build, not ship."""
    template = "\n".join(
        [f"%%FONT_{n}%%" for n in
         ["FRAUNCES_ROMAN", "FRAUNCES_ITALIC", "NEWSREADER_ROMAN", "NEWSREADER_ITALIC"]]
        + [f"%%PLATE_{key}%%" for key, *_ in build_report.FIGURES]
        + ["<p>%%stray_lowercase_token%%</p>"]
    )
    src = tmp_path / "template.html"
    src.write_text(template, encoding="utf-8")
    out = tmp_path / "out.html"
    monkeypatch.setattr(build_report, "TEMPLATE", src)
    monkeypatch.setattr(build_report, "OUTPUT", out)

    with pytest.raises(SystemExit) as excinfo:
        build_report.main(argv=[])
    assert "%%stray_lowercase_token%%" in str(excinfo.value)
    assert not out.exists(), "a failed build must not leave an output file"


def test_main_does_not_mutate_the_module_level_candidates(
    build_report, tmp_path, monkeypatch
):
    """main() used to insert argv[1] into FACTS_CANDIDATES in place, so a second
    call in the same process inherited the first call's path."""
    before = list(build_report.FACTS_CANDIDATES)
    monkeypatch.setattr(build_report, "OUTPUT", tmp_path / "out.html")
    build_report.main(argv=[])
    assert build_report.FACTS_CANDIDATES == before


def test_a_named_but_missing_facts_file_is_an_error_not_a_fallback(
    build_report, tmp_path
):
    """Silently building from a different digest than the one asked for is the
    same failure mode as the removed scratchpad path."""
    with pytest.raises(SystemExit, match="facts file not found"):
        build_report.load_facts(tmp_path / "nope.json")


def test_load_facts_override_is_preferred_but_not_sticky(build_report, tmp_path):
    other = tmp_path / "other_facts.json"
    other.write_text('{"n_figures": 99, "figures": []}', encoding="utf-8")
    assert build_report.load_facts(other)["n_figures"] == 99
    # no override -> back to the committed digest
    assert build_report.load_facts()["n_figures"] == 13


# ---------------------------------------------------------------------------
# the generated report
# ---------------------------------------------------------------------------
def test_report_exists_and_is_self_contained(report_html):
    assert len(report_html) > 1_000_000, "report looks truncated"


def test_report_embeds_thirteen_data_uri_plates(report_html):
    plates = re.findall(r'<img[^>]+src="data:(image/[a-z+]+);base64,', report_html)
    assert len(plates) == 13, f"expected 13 embedded figures, found {len(plates)}"
    assert set(plates) == {"image/svg+xml", "image/png"}


def test_report_has_no_external_src_or_href(report_html):
    external = [
        u for u in re.findall(r'(?:src|href)\s*=\s*"([^"]+)"', report_html)
        if u.startswith("http://") or u.startswith("https://")
    ]
    assert external == [], f"report is not self-contained: {external}"


def test_report_has_no_unresolved_tokens_including_lowercase(report_html):
    leftover = re.findall(r"%%[A-Za-z0-9_]+%%", report_html)
    assert leftover == [], f"unresolved template tokens: {sorted(set(leftover))}"


def test_caveat_ribbon_is_present_and_prominent(report_html):
    assert 'class="caveat"' in report_html
    assert 'aria-label="Interpretation caveat"' in report_html
    assert "technical replicates" in report_html
    assert "hypothesis-generating" in report_html
    # The ribbon must still carry the headline null verbatim.
    assert re.search(r"0\s*of\s*1,938", report_html), "headline null missing from ribbon"


def test_report_carries_the_corrected_d7_orientation(report_html):
    """The numbers the hand-authored facts file had backwards."""
    assert re.search(r'class="val up num">509<', report_html)
    assert re.search(r'class="val down num">206<', report_html)
    assert "509 up" in report_html and "206 down" in report_html


def test_report_states_the_inversion_was_found_and_corrected(report_html):
    """A corrected result should not be presented as if nothing happened."""
    assert "inverted" in report_html
    assert "Vehicle_Rep1_31579" in report_html
    assert "sample_sheet.tsv" in report_html


def test_report_carries_the_d9_and_d11_numbers(report_html):
    assert "min adj.p = 0.116" in report_html          # trend/robust is default
    assert "Min adj.p, vanilla eBayes" in report_html  # vanilla is the comparison
    assert "2,552" in report_html                       # D11 background
    assert "604" in report_html                         # D11 single-condition
    # The old background may appear ONCE, in the sentence describing the change.
    assert report_html.count("2,554") == 1
    assert "2,554 → 2,552" in report_html


def test_stale_figure_captions_are_acknowledged_not_hidden(report_html):
    """Two QC figures predate the D11 quarantine and still say 2554/606.

    build_facts.py copies manifest captions verbatim, so those cannot be fixed
    by editing the digest -- only by re-running the figure layer. The report
    must therefore explain the mismatch rather than leave it looking wrong.
    """
    stale_captions = report_html.count("2554-protein identified universe")
    assert stale_captions == 2, f"expected 2 pre-D11 captions, found {stale_captions}"
    assert "were last drawn before the two junk accessions were quarantined" in report_html


def test_report_forward_path_claim_is_accurate(report_html):
    """The paragraph must separate what regenerates from what merely guards."""
    assert "Genuinely automatic" in report_html
    assert "Guarded, but not automatic" in report_html
    # the three stages that do NOT generalise must be named
    for stage in ("foldchange.py", "replicate_check.py", "qc/schema.py"):
        assert stage in report_html, f"{stage} not named in the forward-path section"
    assert "2-channel SILAC-specific" in report_html
    # and the claim that used to overreach must be gone
    assert "replicate-count-agnostic by construction" not in report_html
    assert "No code rewrite required" not in report_html


def test_report_matches_the_facts_digest_figure_titles(build_report, report_html):
    """Every plate caption in the HTML comes from report_facts.json."""
    import html as html_mod

    facts = build_report.load_facts()
    for entry in facts["figures"]:
        assert html_mod.escape(entry["title"]) in report_html
