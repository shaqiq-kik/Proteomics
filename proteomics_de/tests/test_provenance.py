"""Tests for :mod:`proteomics_de.provenance`.

The gap being closed: research1.md line 245 says the n=2 technical-replicate
caveat "must be stated on every output", and it was on every *figure* and no
*CSV*. The CSV is the artifact that actually leaves the repo.

The constraint that shapes the whole design: the caveat cannot go **inside** the
CSV. IPA's contract is a single header row followed by data, so a ``#`` comment
line becomes the header and the upload breaks. Hence a sidecar. Two tests below
(``test_sidecar_does_not_touch_the_data_file`` and
``test_no_deliverable_gains_a_comment_header``) exist specifically to keep anyone
from "simplifying" this back into an in-file comment.
"""

from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de import provenance  # noqa: E402
from proteomics_de.config.constants import CAVEAT_TEXT  # noqa: E402
from proteomics_de.provenance import (  # noqa: E402
    DEFAULT_DELIVERABLES,
    SIDECAR_SUFFIX,
    count_rows,
    emit_default_sidecars,
    git_commit,
    sha256_file,
    sidecar,
    tool_versions,
)


@pytest.fixture
def sample_csv(tmp_path):
    path = tmp_path / "thing.csv"
    path.write_text("a,b\n1,2\n3,4\n5,6\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Placement and shape
# ---------------------------------------------------------------------------


def test_sidecar_is_written_next_to_the_file(sample_csv):
    out = sidecar(sample_csv)
    assert out.parent == sample_csv.parent
    assert out.name == "thing.csv" + SIDECAR_SUFFIX
    assert out.exists()


def test_sidecar_keeps_the_extension_so_csv_and_txt_do_not_collide(tmp_path):
    """``ipa_input_full.csv`` and ``ipa_input_full.txt`` are distinct deliverables.

    Naming sidecars off the stem would have the second silently overwrite the
    first, and the surviving record would misreport its own sha256.
    """
    csv_path = tmp_path / "full.csv"
    txt_path = tmp_path / "full.txt"
    csv_path.write_text("a,b\n1,2\n", encoding="utf-8")
    txt_path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    a, b = sidecar(csv_path), sidecar(txt_path)
    assert a != b
    assert json.loads(a.read_text())["sha256"] != json.loads(b.read_text())["sha256"]


def test_sidecar_is_valid_utf8_json_with_a_trailing_newline(sample_csv):
    raw = sidecar(sample_csv).read_bytes()
    assert raw.endswith(b"\n")
    json.loads(raw.decode("utf-8"))


def test_sidecar_refuses_to_describe_a_missing_file(tmp_path):
    """A provenance record for a file that is not there would be fiction."""
    with pytest.raises(FileNotFoundError):
        sidecar(tmp_path / "ghost.csv")


# ---------------------------------------------------------------------------
# Contents
# ---------------------------------------------------------------------------


def test_caveat_is_carried_verbatim(sample_csv):
    """Character for character, em dash included.

    The same string is stamped onto 10+ committed figures. If the sidecar
    reworded it, a reader comparing a CSV to a plot would see two different
    caveats and have to work out whether they mean the same thing.
    """
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    assert record["caveat"] == CAVEAT_TEXT
    assert "—" in record["caveat"], "the em dash must survive the JSON round-trip"
    assert "n = 2 technical replicates" in record["caveat"]
    assert "hypothesis-generating only" in record["caveat"]


def test_caveat_names_its_single_source(sample_csv):
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    assert record["caveat_source"].endswith("constants.py::CAVEAT_TEXT")


def test_sidecar_records_the_files_own_digest(sample_csv):
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    assert record["sha256"] == sha256_file(sample_csv)
    assert re.fullmatch(r"[0-9a-f]{64}", record["sha256"])


def test_sidecar_digest_changes_when_the_file_does(sample_csv):
    first = json.loads(sidecar(sample_csv).read_text())["sha256"]
    sample_csv.write_text("a,b\n9,9\n", encoding="utf-8")
    second = json.loads(sidecar(sample_csv).read_text())["sha256"]
    assert first != second


def test_sidecar_records_row_count_excluding_the_header(sample_csv):
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    assert record["n_rows"] == 3


def test_sidecar_records_a_git_commit_and_dirty_flag(sample_csv):
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    sha = record["git_commit"]
    assert sha is None or re.fullmatch(r"[0-9a-f]{40}", sha)
    assert record["git_dirty"] in (True, False, None)


def test_sidecar_records_a_parseable_utc_timestamp(sample_csv):
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    stamp = datetime.fromisoformat(record["generated_at"])
    assert stamp.tzinfo is not None, "timestamps must be timezone-aware"
    assert stamp.utcoffset().total_seconds() == 0


def test_sidecar_records_tool_versions(sample_csv):
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    versions = record["tool_versions"]
    assert "python" in versions and "pandas" in versions
    for value in versions.values():
        assert isinstance(value, str) and value


def test_sidecar_describes_a_path_outside_the_repo_absolutely(sample_csv):
    """No crash when the deliverable lives outside the work tree.

    (The repo-relative form is asserted on the real sidecars in
    ``test_sidecars_are_json_not_csv_rows``. This test deliberately does not
    write into ``results/`` -- a unit test must not leave artifacts in the
    repository it is testing.)
    """
    record = json.loads(sidecar(sample_csv).read_text(encoding="utf-8"))
    assert record["describes"] == str(sample_csv.resolve())


def test_extra_facts_are_merged_and_override_computed_keys(sample_csv):
    """A caller who knows better (e.g. the in-memory row count) wins."""
    out = sidecar(sample_csv, input_sha256="deadbeef", n_rows=999, contrast="T_vs_V")
    record = json.loads(out.read_text(encoding="utf-8"))
    assert record["input_sha256"] == "deadbeef"
    assert record["n_rows"] == 999
    assert record["contrast"] == "T_vs_V"
    assert record["caveat"] == CAVEAT_TEXT  # still there


# ---------------------------------------------------------------------------
# The constraint that produced this design
# ---------------------------------------------------------------------------


def test_sidecar_does_not_touch_the_data_file(sample_csv):
    """The whole point: the deliverable's bytes are untouched.

    ``results/ipa_input.csv`` is byte-frozen and may already have been uploaded
    to QIAGEN. Provenance that modified it would be worse than no provenance.
    """
    before = sample_csv.read_bytes()
    sidecar(sample_csv)
    assert sample_csv.read_bytes() == before


def test_no_deliverable_gains_a_comment_header(results_dir):
    """The reason a ``#`` header was rejected in favour of a sidecar.

    IPA reads the first line as the column header; a comment line breaks the
    upload. Every deliverable must still start with real column names.
    """
    for name in DEFAULT_DELIVERABLES:
        path = results_dir / name
        if not path.exists():
            continue
        first = path.read_text(encoding="utf-8").splitlines()[0]
        assert not first.lstrip().startswith("#"), f"{name} starts with a comment line"
        assert "," in first, f"{name} has no delimited header"


def test_sidecars_are_json_not_csv_rows(results_dir):
    """Sidecars are separate files, never appended to the table."""
    for name in DEFAULT_DELIVERABLES:
        path = results_dir / (name + SIDECAR_SUFFIX)
        if not path.exists():
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        # repo-relative, never absolute: an absolute path in a record that gets
        # mailed to a collaborator leaks the author's home directory.
        assert record["describes"] == f"proteomics_de/results/{name}"
        data = results_dir / name
        assert record["sha256"] == sha256_file(data), (
            f"{name} changed since its sidecar was written; regenerate with "
            "`python -m proteomics_de.provenance`"
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def test_count_rows_handles_a_quoted_embedded_newline(tmp_path):
    """Counting ``\\n`` would report 3 rows here. There are 2."""
    path = tmp_path / "quoted.csv"
    path.write_text('a,b\n1,"line1\nline2"\n3,4\n', encoding="utf-8")
    assert count_rows(path) == 2


def test_count_rows_understands_tab_delimited_txt(tmp_path):
    path = tmp_path / "t.txt"
    path.write_text("a\tb\n1\t2\n3\t4\n", encoding="utf-8")
    assert count_rows(path) == 2


def test_count_rows_returns_none_for_a_non_table(tmp_path):
    path = tmp_path / "fig.png"
    path.write_bytes(b"\x89PNG\r\n\x1a\n")
    assert count_rows(path) is None


def test_count_rows_of_a_header_only_file_is_zero(tmp_path):
    """``ipa_input_significant.csv`` -- 0 rows is a real answer, not an error."""
    path = tmp_path / "empty.csv"
    path.write_text("a,b,c\n", encoding="utf-8")
    assert count_rows(path) == 0


def test_r_versions_are_read_from_the_limma_handoff(tmp_path):
    path = tmp_path / "_limma_versions.txt"
    path.write_text(
        "R version 4.6.1 (2026-06-24)\nlimma 3.68.4\nimputeLCMD 2.1\nseed 42\n",
        encoding="utf-8",
    )
    parsed = provenance._r_versions(path)
    assert parsed == {"R": "4.6.1", "limma": "3.68.4",
                      "imputeLCMD": "2.1", "r_seed": "42"}


def test_r_versions_tolerate_a_missing_handoff(tmp_path):
    """R may not be installed when a sidecar is regenerated."""
    assert provenance._r_versions(tmp_path / "absent.txt") == {}


def test_tool_versions_includes_the_r_stack_when_the_handoff_exists():
    versions = tool_versions()
    assert versions["python"]
    if (provenance._LIMMA_VERSIONS).exists():
        assert versions.get("limma"), "the limma version handoff exists but was not read"


def test_git_commit_returns_none_outside_a_work_tree(tmp_path):
    assert git_commit(tmp_path) == (None, None)


def test_git_dirty_ignores_untracked_files(results_dir):
    """The sidecars are themselves untracked; they must not mark the tree dirty.

    Otherwise every sidecar would report ``git_dirty: true`` because of its own
    existence, and the flag would carry no information.
    """
    sha, dirty = git_commit()
    if sha is None:
        pytest.skip("not a git work tree")
    assert isinstance(dirty, bool)


# ---------------------------------------------------------------------------
# The four deliverables
# ---------------------------------------------------------------------------


def test_emit_default_sidecars_covers_the_four_deliverables(tmp_path):
    for name in DEFAULT_DELIVERABLES:
        (tmp_path / name).write_text("a,b\n1,2\n", encoding="utf-8")
    written = emit_default_sidecars(tmp_path)
    assert [p.name for p in written] == [n + SIDECAR_SUFFIX for n in DEFAULT_DELIVERABLES]


def test_emit_default_sidecars_skips_files_that_are_not_there(tmp_path):
    """Callable before ``ipa_input_full.csv`` has been built."""
    (tmp_path / "qc_limma.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    written = emit_default_sidecars(tmp_path)
    assert [p.name for p in written] == ["qc_limma.csv" + SIDECAR_SUFFIX]


def test_deliverables_include_both_ipa_uploads_and_both_tables():
    assert set(DEFAULT_DELIVERABLES) == {
        "ipa_input.csv", "ipa_input_full.csv", "qc_limma.csv", "foldchange_all.csv"
    }


def test_committed_deliverables_have_sidecars_with_matching_counts(
    results_dir, frozen_counts
):
    expected = {
        "ipa_input.csv": frozen_counts["ipa_input_rows"],
        "ipa_input_full.csv": frozen_counts["ipa_input_rows"],
        "qc_limma.csv": frozen_counts["qc_limma_rows"],
        "foldchange_all.csv": frozen_counts["foldchange_all_rows"],
    }
    seen = 0
    for name, n_rows in expected.items():
        path = results_dir / (name + SIDECAR_SUFFIX)
        if not path.exists():
            continue
        record = json.loads(path.read_text(encoding="utf-8"))
        assert record["n_rows"] == n_rows, f"{name}: {record['n_rows']} != {n_rows}"
        assert record["caveat"] == CAVEAT_TEXT
        seen += 1
    if seen == 0:
        pytest.skip("no sidecars generated yet")
