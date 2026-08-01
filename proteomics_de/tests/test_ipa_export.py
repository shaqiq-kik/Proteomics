"""Tests for :mod:`proteomics_de.export.ipa_export`.

Three things are being defended here, in order of how expensive they are to get
wrong:

1. **``results/ipa_input.csv`` must never move.** It is a committed output that
   the user may already have uploaded to QIAGEN IPA. Every test that touches
   :func:`write_ipa` checks bytes, not values.
2. **``ipa_input_full.csv`` must agree with it textually.** The two files are
   uploaded to the same tool and cross-referenced by a human; ``log2FC`` reading
   ``1.0633144767505383`` in one and ``...385`` in the other would be a real
   support burden for a difference that is not real. This is what makes the
   ``float_precision="round_trip"`` reads load-bearing rather than decorative.
3. **``ipa_input_significant.csv`` stays empty.** 0/1938 pass FDR < 0.05
   (DECISIONS_LOG D2). The test asserts the emptiness as a *contract*, so that
   nobody later "fixes" the file by loosening the cutoff.
"""

from __future__ import annotations

import csv
import io
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # works with or without a rootdir conftest
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.export import ipa_export  # noqa: E402
from proteomics_de.export.ipa_export import (  # noqa: E402
    IPA_COLS,
    IPA_FULL_COLS,
    IPA_ID_COL,
    IPA_SIGNIFICANT_COLS,
    IpaValidationError,
    build_ipa_full,
    expected_ipa_rows,
    validate_ipa,
    validate_ipa_significant,
    write_ipa,
    write_ipa_full,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def ipa_csv(results_dir):
    path = results_dir / "ipa_input.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return path


@pytest.fixture(scope="module")
def committed_ipa(ipa_csv):
    """The committed 4-column IPA frame, read losslessly."""
    return pd.read_csv(ipa_csv, float_precision="round_trip")


@pytest.fixture(scope="module")
def committed_limma(results_dir):
    path = results_dir / "qc_limma.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return pd.read_csv(path, float_precision="round_trip")


@pytest.fixture
def tiny_ipa():
    """A 3-row stand-in with the same shape as the real file."""
    return pd.DataFrame(
        {
            IPA_ID_COL: ["P00001", "P00002;P00003", "Q00004"],
            "Gene names": ["Aaa", "Bbb;Ccc", None],
            "log2FC": [1.5, -0.9, 0.7],
            "regulated": ["UP", "DOWN", "UP"],
        }
    )


@pytest.fixture
def tiny_limma(tiny_ipa):
    return pd.DataFrame(
        {
            "id": list(tiny_ipa[IPA_ID_COL]) + ["Z99999"],
            "p_value": [0.01, 0.2, 0.3, 0.4],
            "adj_p_value": [0.9, 0.95, 0.99, 0.99],
        }
    )


def _write(path, header, rows, delim=","):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter=delim)
        w.writerow(header)
        w.writerows(rows)
    return path


# ---------------------------------------------------------------------------
# write_ipa is byte-frozen
# ---------------------------------------------------------------------------


def test_write_ipa_reproduces_the_committed_file_byte_for_byte(
    ipa_csv, committed_ipa, tmp_path
):
    """The frozen output regenerates exactly. This is the whole contract."""
    out = tmp_path / "ipa_input.csv"
    write_ipa(committed_ipa, out, IPA_COLS)
    assert out.read_bytes() == ipa_csv.read_bytes()


def test_write_ipa_emits_no_bom(tiny_ipa, tmp_path):
    """``encoding="utf-8"``, never ``"utf-8-sig"``.

    A BOM would make IPA read the first header cell as
    ``"\\ufeffUniProt Accession Number"`` -- the exact "IPA does not recognize the
    column headers" failure the earlier Excel->CSV fix was chasing.
    """
    out = tmp_path / "x.csv"
    write_ipa(tiny_ipa, out, IPA_COLS)
    assert out.read_bytes()[:3] != b"\xef\xbb\xbf"


def test_write_ipa_writes_no_index_column(tiny_ipa, tmp_path):
    """An index column would push the identifier out of position 0."""
    out = tmp_path / "x.csv"
    write_ipa(tiny_ipa.iloc[::-1], out, IPA_COLS)  # non-trivial index
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert header.split(",")[0] == IPA_ID_COL


def test_write_ipa_does_not_filter_or_reorder(tiny_ipa, tmp_path):
    """Row selection belongs to the caller; the writer is dumb on purpose."""
    out = tmp_path / "x.csv"
    write_ipa(tiny_ipa, out, IPA_COLS)
    back = pd.read_csv(out, float_precision="round_trip")
    assert list(back[IPA_ID_COL]) == list(tiny_ipa[IPA_ID_COL])


# ---------------------------------------------------------------------------
# The limma join
# ---------------------------------------------------------------------------


def test_build_ipa_full_shape_and_columns(tiny_ipa, tiny_limma):
    full = build_ipa_full(tiny_ipa, tiny_limma)
    assert list(full.columns) == IPA_FULL_COLS
    assert len(full) == len(tiny_ipa)
    assert list(full[IPA_ID_COL]) == list(tiny_ipa[IPA_ID_COL])


def test_build_ipa_full_on_the_real_data_has_no_missing_pvalues(
    committed_ipa, committed_limma, frozen_counts
):
    """The eligibility invariant, checked against the real files.

    Every IPA row is ``complete & regulated in {UP, DOWN}``, and limma tests
    every complete row (only ON_OFF proteins are excluded, and those never carry
    a directional label). So the left join must be total.
    """
    full = build_ipa_full(committed_ipa, committed_limma)
    assert len(full) == frozen_counts["ipa_input_rows"]
    assert full["p_value"].notna().all()
    assert full["adj_p_value"].notna().all()


def test_build_ipa_full_raises_when_a_row_has_no_limma_result(tiny_ipa, tiny_limma):
    """A NaN p-value means the eligibility assumption broke -- fail loudly.

    Shipping the file anyway would put blank measurement cells in front of IPA,
    which drops or misreads the row instead of erroring.
    """
    orphaned = tiny_limma[tiny_limma["id"] != "Q00004"]
    with pytest.raises(IpaValidationError, match="no limma"):
        build_ipa_full(tiny_ipa, orphaned)


def test_build_ipa_full_rejects_duplicate_accessions(tiny_ipa, tiny_limma):
    dupes = pd.concat([tiny_ipa, tiny_ipa.iloc[[0]]], ignore_index=True)
    with pytest.raises(IpaValidationError, match="duplicate accessions in the IPA frame"):
        build_ipa_full(dupes, tiny_limma)


def test_build_ipa_full_rejects_duplicate_limma_ids(tiny_ipa, tiny_limma):
    dupes = pd.concat([tiny_limma, tiny_limma.iloc[[0]]], ignore_index=True)
    with pytest.raises(IpaValidationError, match="duplicate accessions in qc_limma"):
        build_ipa_full(tiny_ipa, dupes)


def test_build_ipa_full_rejects_a_limma_frame_without_pvalues(tiny_ipa, tiny_limma):
    with pytest.raises(IpaValidationError, match="missing column 'adj_p_value'"):
        build_ipa_full(tiny_ipa, tiny_limma.drop(columns=["adj_p_value"]))


# ---------------------------------------------------------------------------
# The full file, and its textual agreement with the frozen one
# ---------------------------------------------------------------------------


def test_full_file_log2fc_text_matches_ipa_input_exactly(
    committed_ipa, committed_limma, ipa_csv, tmp_path
):
    """Byte-level agreement on the shared columns.

    This is the test that keeps ``float_precision="round_trip"`` in the source.
    pandas' default CSV float parser is not correctly rounded: drop the
    parameter and 87 of the 715 ``log2FC`` values shift by one ULP, so the two
    files a human uploads to IPA would disagree in their last digit for rows
    that are numerically identical. Asserting on the *text* catches that
    regardless of which pandas version introduced it.
    """
    full = build_ipa_full(committed_ipa, committed_limma)
    out = tmp_path / "ipa_input_full.csv"
    write_ipa_full(full, out, txt_path=False)

    frozen_lines = ipa_csv.read_text(encoding="utf-8").splitlines()
    full_lines = out.read_text(encoding="utf-8").splitlines()
    assert len(frozen_lines) == len(full_lines)
    for frozen, combined in zip(frozen_lines, full_lines):
        # strip the two appended measurement fields
        assert combined.rsplit(",", 2)[0] == frozen


def test_write_ipa_full_emits_csv_and_tab_delimited_txt(tiny_ipa, tiny_limma, tmp_path):
    """QIAGEN "strongly recommend ... tab-delimited text format (.txt ...)"."""
    full = build_ipa_full(tiny_ipa, tiny_limma)
    written = write_ipa_full(full, tmp_path / "full.csv")
    assert [p.name for p in written] == ["full.csv", "full.txt"]

    txt = written[1].read_bytes()
    assert txt[:3] != b"\xef\xbb\xbf"
    assert txt.splitlines()[0].count(b"\t") == len(IPA_FULL_COLS) - 1

    as_csv = pd.read_csv(written[0], float_precision="round_trip")
    as_txt = pd.read_csv(written[1], sep="\t", float_precision="round_trip")
    assert as_csv.equals(as_txt)


def test_committed_full_files_if_present(results_dir, frozen_counts):
    """The generated deliverables, when they exist on disk."""
    csv_path = results_dir / "ipa_input_full.csv"
    txt_path = results_dir / "ipa_input_full.txt"
    if not csv_path.exists():
        pytest.skip("ipa_input_full.csv not generated yet")

    for path, sep in ((csv_path, ","), (txt_path, "\t")):
        assert path.exists(), f"{path} missing"
        frame = pd.read_csv(path, sep=sep, float_precision="round_trip")
        assert list(frame.columns) == IPA_FULL_COLS
        assert len(frame) == frozen_counts["ipa_input_rows"]
        assert frame["p_value"].notna().all()
        assert frame["adj_p_value"].notna().all()
        counts = frame["regulated"].value_counts()
        assert int(counts["UP"]) == frozen_counts["n_up"]
        assert int(counts["DOWN"]) == frozen_counts["n_down"]


# ---------------------------------------------------------------------------
# validate_ipa -- research1.md line 187, for real
# ---------------------------------------------------------------------------


def test_validate_ipa_passes_on_the_committed_file(ipa_csv):
    report = validate_ipa(ipa_csv, required_columns=IPA_COLS)
    assert report["ok"], report["failures"]
    assert report["n_rows"] == expected_ipa_rows()
    assert all(c["ok"] for c in report["checks"].values())


def test_validate_ipa_passes_on_the_full_files(results_dir):
    for name, in (("ipa_input_full.csv",), ("ipa_input_full.txt",)):
        path = results_dir / name
        if not path.exists():
            pytest.skip(f"{name} not generated yet")
        report = validate_ipa(path, required_columns=IPA_FULL_COLS)
        assert report["ok"], report["failures"]


def test_validate_ipa_rejects_a_bom(committed_ipa, tmp_path):
    out = tmp_path / "bom.csv"
    write_ipa(committed_ipa, out, IPA_COLS)
    out.write_bytes(b"\xef\xbb\xbf" + out.read_bytes())
    with pytest.raises(IpaValidationError, match="BOM"):
        validate_ipa(out, required_columns=IPA_COLS)


def test_validate_ipa_rejects_identifier_not_leftmost(tiny_ipa, tmp_path):
    out = tmp_path / "swapped.csv"
    reordered = ["Gene names", IPA_ID_COL, "log2FC", "regulated"]
    tiny_ipa[reordered].to_csv(out, index=False, encoding="utf-8")
    report = validate_ipa(out, expected_rows=3, required_columns=IPA_COLS, strict=False)
    assert not report["checks"]["identifier_leftmost"]["ok"]


def test_validate_ipa_rejects_na_in_a_measurement_column(tmp_path):
    out = _write(
        tmp_path / "na.csv", IPA_COLS,
        [["P1", "A", "1.0", "UP"], ["P2", "B", "", "DOWN"]],
    )
    with pytest.raises(IpaValidationError, match="measurements_finite"):
        validate_ipa(out, expected_rows=2, required_columns=IPA_COLS)


def test_validate_ipa_rejects_inf_in_a_measurement_column(tmp_path):
    """Bug 3's failure mode, caught at the file boundary this time."""
    out = _write(
        tmp_path / "inf.csv", IPA_COLS,
        [["P1", "A", "inf", "UP"], ["P2", "B", "1.0", "DOWN"]],
    )
    with pytest.raises(IpaValidationError, match="measurements_finite"):
        validate_ipa(out, expected_rows=2, required_columns=IPA_COLS)


def test_validate_ipa_rejects_a_repeated_header_row(tmp_path):
    """Two files concatenated -- the classic way a "merged" export goes wrong."""
    out = _write(
        tmp_path / "double.csv", IPA_COLS,
        [["P1", "A", "1.0", "UP"], IPA_COLS, ["P2", "B", "1.0", "DOWN"]],
    )
    report = validate_ipa(out, expected_rows=3, required_columns=IPA_COLS, strict=False)
    assert not report["checks"]["single_header_row"]["ok"]


def test_validate_ipa_rejects_ragged_rows(tmp_path):
    out = tmp_path / "ragged.csv"
    out.write_text(
        ",".join(IPA_COLS) + "\nP1,A,1.0,UP\nP2,B,1.0\n", encoding="utf-8"
    )
    report = validate_ipa(out, expected_rows=2, required_columns=IPA_COLS, strict=False)
    assert not report["checks"]["single_header_row"]["ok"]


def test_validate_ipa_rejects_the_wrong_row_count(committed_ipa, tmp_path):
    out = tmp_path / "short.csv"
    write_ipa(committed_ipa.iloc[:-1], out, IPA_COLS)
    with pytest.raises(IpaValidationError, match="row_count"):
        validate_ipa(out, required_columns=IPA_COLS)


def test_validate_ipa_allows_protein_group_accessions(tmp_path):
    """``P05132;P68181`` is a legitimate protein group, not a stray character.

    21 of the 715 committed accessions look like this. A naive "identifier must
    match the UniProt regex" check would reject every one of them.
    """
    out = _write(
        tmp_path / "groups.csv", IPA_COLS,
        [["P05132;P68181", "Prkaca;Prkacb", "1.0", "UP"]],
    )
    report = validate_ipa(out, expected_rows=1, required_columns=IPA_COLS)
    assert report["checks"]["identifier_clean"]["ok"]


@pytest.mark.parametrize("bad", [" P00001", "P00001 ", "P0 0001", "", "P00001\t"])
def test_validate_ipa_rejects_stray_characters_in_identifiers(bad, tmp_path):
    out = _write(tmp_path / "dirty.csv", IPA_COLS, [[bad, "A", "1.0", "UP"]])
    report = validate_ipa(out, expected_rows=1, required_columns=IPA_COLS, strict=False)
    assert not report["checks"]["identifier_clean"]["ok"]


def test_validate_ipa_reports_a_missing_file(tmp_path):
    with pytest.raises(IpaValidationError, match="does not exist"):
        validate_ipa(tmp_path / "nope.csv")


def test_validate_ipa_non_strict_collects_every_failure(tmp_path):
    """The report lists all violations at once, not just the first."""
    out = tmp_path / "awful.csv"
    out.write_bytes(
        b"\xef\xbb\xbf" + b"Gene names,UniProt Accession Number,log2FC\nA,P1,\n"
    )
    report = validate_ipa(out, expected_rows=99, required_columns=IPA_COLS, strict=False)
    assert not report["ok"]
    assert len(report["failures"]) >= 4


# ---------------------------------------------------------------------------
# The header-only significant file -- a finding, not a bug
# ---------------------------------------------------------------------------


def test_significant_file_is_header_only_by_design(results_dir, frozen_counts):
    """0/1938 pass FDR < 0.05. That emptiness IS the result (DECISIONS_LOG D2).

    If this ever fails because rows appeared, the correct response is to check
    what changed in the data -- never to relax ADJ_P_THRESHOLD until the file
    fills up.
    """
    path = results_dir / "ipa_input_significant.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    report = validate_ipa_significant(path)
    assert report["ok"], report["failures"]
    assert report["columns"] == IPA_SIGNIFICANT_COLS
    assert report["n_rows"] == 0 == frozen_counts["n_significant_fdr05"]


def test_significant_file_contract_fails_if_a_row_appears(tmp_path):
    out = _write(
        tmp_path / "sig.csv", IPA_SIGNIFICANT_COLS,
        [["P00001", "Aaa", "1.5", "UP", "0.01"]],
    )
    with pytest.raises(IpaValidationError, match="row_count"):
        validate_ipa_significant(out)


def test_significant_file_contract_fails_on_wrong_columns(tmp_path):
    out = _write(tmp_path / "sig.csv", IPA_COLS, [])
    with pytest.raises(IpaValidationError, match="required_columns"):
        validate_ipa_significant(out)


# ---------------------------------------------------------------------------
# Expected counts come from the data / the JSON, never from literals
# ---------------------------------------------------------------------------


def test_expected_rows_derived_from_foldchange_matches_frozen_json(
    results_dir, frozen_counts
):
    """research1.md line 187 states the rule as ``(complete & UP|DOWN).sum()``.

    Recomputing it from ``foldchange_all.csv`` and comparing to the frozen JSON
    checks both directions at once: the JSON is not stale, and the derivation is
    not wrong.
    """
    fc = results_dir / "foldchange_all.csv"
    if not fc.exists():
        pytest.skip("foldchange_all.csv not present")
    assert expected_ipa_rows(fc) == frozen_counts["ipa_input_rows"]


def test_frozen_counts_are_internally_consistent(frozen_counts):
    """D7: 509 UP + 206 DOWN = the 715-row IPA total."""
    assert frozen_counts["n_up"] + frozen_counts["n_down"] == frozen_counts["ipa_input_rows"]


def test_committed_direction_counts_match_d7(committed_ipa, frozen_counts):
    """The control/treated inversion was corrected: 509 UP / 206 DOWN."""
    counts = committed_ipa["regulated"].value_counts()
    assert int(counts["UP"]) == frozen_counts["n_up"]
    assert int(counts["DOWN"]) == frozen_counts["n_down"]
    assert set(counts.index) == {"UP", "DOWN"}


def test_every_read_csv_asks_for_round_trip_float_precision():
    """A source-level policy check, because the data cannot fully test this.

    Dropping ``float_precision="round_trip"`` from the ``ipa_input.csv`` read is
    caught by ``test_build_full_from_results_preserves_every_source_digit`` (87
    of 715 ``log2FC`` values move). Dropping it from the ``qc_limma.csv`` read
    happens to be harmless *on this dataset* -- pandas' default parser reproduces
    all 715 p-values exactly -- so no data-driven test can catch that one. It is
    still wrong: the default parser is not correctly rounded, and the next
    dataset has no reason to be as lucky. Pin the policy at the source instead of
    pretending the coverage exists.
    """
    with open(ipa_export.__file__, encoding="utf-8") as fh:
        source = fh.read()
    calls = re.findall(r"pd\.read_csv\((?:[^()]|\([^()]*\))*\)", source)
    assert calls, "no pd.read_csv calls found -- did the module move?"
    missing = [c for c in calls if 'float_precision="round_trip"' not in c]
    assert not missing, (
        "these reads would silently perturb floats in the last ULP: " + str(missing)
    )


def test_module_hardcodes_no_dataset_specific_counts(frozen_counts):
    """The counts live in one place: tests/expected/frozen_counts.json.

    Eight pipeline files used to carry inline ``assert len(df) == 1948``-style
    literals. This test stops ``ipa_export.py`` from becoming the ninth.
    """
    source = ipa_export.__file__
    with open(source, encoding="utf-8") as fh:
        text = fh.read()
    # strip the prose, where "715" and "1938" legitimately appear as narrative
    code = re.sub(r'"""(?:.|\n)*?"""', "", text)
    code = re.sub(r"#.*", "", code)
    forbidden = {
        str(frozen_counts[k])
        for k in ("ipa_input_rows", "qc_limma_rows", "n_up", "n_down",
                  "foldchange_all_rows")
    }
    found = sorted(n for n in forbidden if re.search(rf"(?<![\w.]){n}(?![\w.])", code))
    assert not found, f"dataset-specific literals in ipa_export.py: {found}"


# ---------------------------------------------------------------------------
# End-to-end through the results directory
# ---------------------------------------------------------------------------


def test_build_full_from_results_is_idempotent(results_dir, tmp_path):
    """Running twice produces identical bytes -- no timestamps, no ordering drift."""
    for name in ("ipa_input.csv", "qc_limma.csv"):
        src = results_dir / name
        if not src.exists():
            pytest.skip(f"{name} not present")
        (tmp_path / name).write_bytes(src.read_bytes())

    first = ipa_export.build_full_from_results(tmp_path)
    snapshot = {p.name: p.read_bytes() for p in first}
    ipa_export.build_full_from_results(tmp_path)
    assert {p.name: p.read_bytes() for p in first} == snapshot


def test_build_full_from_results_leaves_ipa_input_untouched(results_dir, tmp_path):
    """The frozen file is read, never rewritten."""
    for name in ("ipa_input.csv", "qc_limma.csv"):
        src = results_dir / name
        if not src.exists():
            pytest.skip(f"{name} not present")
        (tmp_path / name).write_bytes(src.read_bytes())

    before = (tmp_path / "ipa_input.csv").read_bytes()
    ipa_export.build_full_from_results(tmp_path)
    assert (tmp_path / "ipa_input.csv").read_bytes() == before
    assert before == (results_dir / "ipa_input.csv").read_bytes()


def test_build_full_from_results_preserves_every_source_digit(results_dir, tmp_path):
    """The production path, checked textually against BOTH of its inputs.

    ``test_full_file_log2fc_text_matches_ipa_input_exactly`` proves the join and
    the writer are faithful, but it hands them frames the test itself read
    correctly. This one goes through :func:`build_full_from_results`, which does
    its own reading -- so it is the test that actually pins
    ``float_precision="round_trip"`` on those reads. Without it, pandas' default
    (not correctly rounded) parser shifts 87 of the 715 ``log2FC`` values by one
    ULP and ``ipa_input_full.csv`` stops agreeing with the ``ipa_input.csv``
    already sitting in QIAGEN.
    """
    for name in ("ipa_input.csv", "qc_limma.csv"):
        src = results_dir / name
        if not src.exists():
            pytest.skip(f"{name} not present")
        (tmp_path / name).write_bytes(src.read_bytes())

    ipa_export.build_full_from_results(tmp_path)
    full_lines = (tmp_path / "ipa_input_full.csv").read_text(
        encoding="utf-8").splitlines()

    # columns 0-3 must be the frozen file's text, character for character
    frozen_lines = (tmp_path / "ipa_input.csv").read_text(
        encoding="utf-8").splitlines()
    assert [line.rsplit(",", 2)[0] for line in full_lines] == frozen_lines

    # columns 4-5 must be qc_limma.csv's text, character for character
    limma_text = {}
    with open(tmp_path / "qc_limma.csv", newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            limma_text[row["id"]] = (row["p_value"], row["adj_p_value"])
    for line in full_lines[1:]:
        acc = line.split(",", 1)[0]
        p, adj = line.rsplit(",", 2)[1:]
        assert (p, adj) == limma_text[acc], f"{acc}: p-value text drifted"


def test_build_full_from_results_needs_limma_first(tmp_path, committed_ipa):
    write_ipa(committed_ipa, tmp_path / "ipa_input.csv", IPA_COLS)
    with pytest.raises(FileNotFoundError, match="qc_limma.csv"):
        ipa_export.build_full_from_results(tmp_path)


def test_no_ipa_output_contains_a_comment_header(results_dir):
    """A ``#`` line would become IPA's header row and break the upload.

    This is why the n=2 caveat ships as a sidecar rather than as a comment
    inside the CSV.
    """
    for name in ("ipa_input.csv", "ipa_input_full.csv", "ipa_input_full.txt",
                 "ipa_input_significant.csv"):
        path = results_dir / name
        if not path.exists():
            continue
        first = path.read_text(encoding="utf-8").splitlines()[0]
        assert not first.lstrip().startswith("#"), f"{name} starts with a comment"
        assert first.split("\t" if name.endswith(".txt") else ",")[0] == IPA_ID_COL


def test_delimiter_is_inferred_from_the_extension(tiny_ipa, tiny_limma, tmp_path):
    """A tab file validated as CSV would look like a single ragged column."""
    full = build_ipa_full(tiny_ipa, tiny_limma)
    written = write_ipa_full(full, tmp_path / "d.csv")
    for path in written:
        report = validate_ipa(path, expected_rows=3, required_columns=IPA_FULL_COLS)
        assert report["columns"] == IPA_FULL_COLS


def test_validate_ipa_reads_a_stringio_free_path_only(tmp_path):
    """Guard against someone passing a buffer: validation is about bytes on disk."""
    with pytest.raises((TypeError, IpaValidationError, AttributeError)):
        validate_ipa(io.StringIO("a,b\n1,2\n"))
