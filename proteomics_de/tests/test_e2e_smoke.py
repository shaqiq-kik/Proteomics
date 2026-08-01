"""End-to-end runs of the real pipeline — on new data, and on the old data.

WHAT THIS ADDS THAT THE OTHER 500 TESTS DO NOT
----------------------------------------------
Everything else in the suite either exercises a function in isolation or
asserts a frozen fact about the one committed workbook. Neither shape can
answer "does the pipeline still *work*?", because:

* a unit test passes over a hand-built three-row frame that never went through
  ``pd.read_excel``, the L/H merge, the dtype restore, the R subprocess, or the
  CSV round trip; and
* a frozen-count test passes on stale outputs — it reads files, it does not
  produce them — and it only ever describes one dataset, so any branch that
  dataset happens not to exercise is untested by construction. The committed
  workbook, for instance, encodes every absent measurement as a **zero**; the
  blank-cell form of the same thing had no coverage at all.

Two levels, run at different times:

**Level 1** (default, ~8 s) runs ``foldchange.py`` over
``fixtures/mini_sheet.xlsx`` — 31 hand-designed proteins, one per branch —
**including the R limma leg**, and checks every row against the expectation
designed into it. This is the regression net for *future* data.

**Level 2** (``slow`` + ``golden``, opt-in) runs the same entry point over the
real workbook and compares the outputs to the committed ones through
``tools/freeze.py``'s canonical digest. Deselected in CI on purpose — see
``.github/workflows/tests.yml``.

WHY THE DIGEST AND NOT ``shasum``
---------------------------------
Do not be tempted to swap ``freeze.digest`` for a raw byte hash. matplotlib
stamps a wall-clock ``<dc:date>`` into every SVG and salts element ids with
random hex, so two runs of *identical* code produce different SVG bytes.
``freeze.canonical_bytes`` strips exactly those two things and hashes
everything else verbatim. A raw hash would report drift on every single run and
would train everyone to ignore it.

A NOTE ON THE INTERMEDIATE GUARD
--------------------------------
``limma_test.py`` writes its Python<->R handoff files (``_limma_input.csv`` and
friends) next to itself in ``proteomics_de/``, **not** into ``--results-dir``.
Two of them are byte-frozen artifacts. So merely *running* this module would
otherwise corrupt the committed tree and turn ``test_freeze.py`` red. The
autouse :func:`_protect_limma_intermediates` fixture snapshots and restores
them, and then verifies the restore against the freeze manifest so a botched
restore fails loudly instead of silently rewriting the scientific record.
That leak is itself asserted on, deliberately —
see :func:`test_run_touches_only_the_known_handoff_files_outside_results`.

PARALLELISM
-----------
This module writes to shared paths in ``proteomics_de/`` (above) and must not
run concurrently with itself. The suite's default invocation is serial; if you
add ``-n`` for xdist, pass ``--dist loadfile`` so these tests stay on one
worker.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent                 # proteomics_de/
_REPO_ROOT = _PKG_DIR.parent                 # repo root
_FIXTURES_DIR = _TESTS_DIR / "fixtures"

sys.path.insert(0, str(_REPO_ROOT / "tools"))
sys.path.insert(0, str(_FIXTURES_DIR))
import freeze  # noqa: E402  (tools/freeze.py -- the canonical digest)
import make_mini_sheet as mini  # noqa: E402  (the fixture's design table)

FOLDCHANGE_SCRIPT = _PKG_DIR / "foldchange.py"
MINI_SHEET = _FIXTURES_DIR / "mini_sheet.xlsx"
REAL_WORKBOOK = _REPO_ROOT / "Copy of General Sheet.xlsx"
COMMITTED_RESULTS = _PKG_DIR / "results"
FREEZE_MANIFEST = _TESTS_DIR / "expected" / "outputs.sha256"

pytestmark = pytest.mark.e2e

#: Files ``limma_test.py`` writes into ``proteomics_de/`` regardless of
#: ``--results-dir``. See the module docstring; the leak is asserted, not
#: tolerated silently.
_LIMMA_HANDOFF_NAMES = frozenset({
    "_limma_input.csv",
    "_limma_design.tsv",
    "_limma_output.csv",
    "_limma_versions.txt",
    "_limma_output_vanilla.csv",
    "_limma_versions_vanilla.txt",
})

#: Every artifact ``foldchange.py`` is contracted to write under its
#: ``--results-dir``. Chained stages included: centering (Bug 5), replicate
#: correlation (Bug 6) and limma (Bug 7) all run inside the same invocation.
_EXPECTED_OUTPUTS = (
    "foldchange_all.csv",
    "ipa_input.csv",
    "single_condition_proteins.csv",
    "onoff_proteins.csv",
    "qc_centering.csv",
    "qc_replicate_correlation.csv",
    "replicate_correlation.png",
    "qc_limma.csv",
    "qc_limma_vanilla.csv",
    "ipa_input_significant.csv",
    "de/design.tsv",
    "de/intensity_matrix.tsv",
    "de/limma_results.tsv",
    "qc/qc_boundaries.json",
    "qc/quarantine_accessions.csv",
)


# ---------------------------------------------------------------------------
# Tree snapshots
# ---------------------------------------------------------------------------

def _snapshot(root: Path, pattern: str = "**/*") -> dict[str, tuple[int, int]]:
    """``{relpath: (size, mtime_ns)}`` for every file under `root`.

    mtime as well as size on purpose: "writes nothing" has to mean *nothing*,
    including a rewrite with identical content. A stage that reopens a committed
    file and writes the same bytes back is still reaching into the committed
    tree, and the next such stage might not be so lucky with its contents.
    """
    out: dict[str, tuple[int, int]] = {}
    for p in root.glob(pattern):
        if p.is_file():
            st = p.stat()
            out[p.relative_to(root).as_posix()] = (st.st_size, st.st_mtime_ns)
    return out


def _changed(before: dict, after: dict) -> list[str]:
    return sorted(
        set(before) ^ set(after)
        | {k for k in set(before) & set(after) if before[k] != after[k]}
    )


# ---------------------------------------------------------------------------
# Protecting the committed handoff intermediates
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module", autouse=True)
def _protect_limma_intermediates(tmp_path_factory):
    """Snapshot ``proteomics_de/_limma_*`` and put them back afterwards.

    Not defensive programming — running the pipeline in this module *does*
    overwrite two byte-frozen artifacts (``_limma_input.csv``,
    ``_limma_output.csv``), and without this every test session would leave the
    repository dirty and ``test_freeze.py`` red.

    The restore is verified against the freeze manifest, because a restore that
    quietly half-worked would be strictly worse than no restore: the tree would
    look clean and the scientific record would have moved.
    """
    backup = tmp_path_factory.mktemp("limma_intermediates")
    originals = sorted(_PKG_DIR.glob("_limma_*"))
    for p in originals:
        shutil.copy2(p, backup / p.name)
    original_names = {p.name for p in originals}
    try:
        yield
    finally:
        for p in sorted(_PKG_DIR.glob("_limma_*")):
            if p.name not in original_names:
                p.unlink()  # created by a test run; was never committed
        for p in sorted(backup.iterdir()):
            shutil.copy2(p, _PKG_DIR / p.name)

        manifest = freeze.read_manifest(FREEZE_MANIFEST)
        drifted = []
        for rel, (sha, _mode) in manifest.items():
            if "_limma_" not in rel:
                continue
            actual, _ = freeze.digest(_REPO_ROOT / rel)
            if actual != sha:
                drifted.append(rel)
        assert not drifted, (
            "FAILED TO RESTORE the byte-frozen limma intermediates after the "
            f"end-to-end run: {drifted}. The committed scientific record has "
            "moved -- restore them from git before doing anything else."
        )


# ---------------------------------------------------------------------------
# Running the pipeline
# ---------------------------------------------------------------------------

@dataclass
class PipelineRun:
    """One completed ``foldchange.py`` invocation, plus what it touched."""

    results_dir: Path
    stdout: str
    results_before: dict = field(repr=False, default_factory=dict)
    results_after: dict = field(repr=False, default_factory=dict)
    pkg_before: dict = field(repr=False, default_factory=dict)
    pkg_after: dict = field(repr=False, default_factory=dict)
    _cache: dict = field(repr=False, default_factory=dict)

    def read(self, relpath: str) -> pd.DataFrame:
        """Read one produced table, with the exact float parser (see below)."""
        if relpath not in self._cache:
            sep = "\t" if relpath.endswith(".tsv") else ","
            self._cache[relpath] = pd.read_csv(
                self.results_dir / relpath, sep=sep, float_precision="round_trip"
            )
        return self._cache[relpath]


def _run_foldchange(workbook: Path, results_dir: Path, *, expect_baseline: bool):
    """Invoke ``foldchange.py`` as a subprocess and return a :class:`PipelineRun`.

    A subprocess rather than ``foldchange.main([...])`` in-process, for three
    reasons that all bite in practice: the module mutates ``sys.path`` and
    process-global state in ``qc.boundaries`` at import time; the R worker is
    launched with ``cwd`` set to the package directory; and the CLI itself
    (``--input`` / ``--results-dir``) is part of what is under test. Calling
    ``main`` directly would test a slightly different program.
    """
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    if not expect_baseline:
        # The dataset-specific assertions inside foldchange.py (UP == 509, ...)
        # describe the committed workbook. A different dataset legitimately
        # produces different numbers, and this is the documented way to say so
        # -- the assertions are switched off, not deleted.
        env["PDE_EXPECT_BASELINE"] = "0"
    else:
        env.pop("PDE_EXPECT_BASELINE", None)

    results_before = _snapshot(COMMITTED_RESULTS)
    pkg_before = _snapshot(_PKG_DIR, "*")

    proc = subprocess.run(
        [sys.executable, str(FOLDCHANGE_SCRIPT),
         "--input", str(workbook),
         "--results-dir", str(results_dir)],
        cwd=str(_REPO_ROOT), env=env, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise AssertionError(
            f"foldchange.py failed (rc={proc.returncode}) on {workbook.name}\n"
            f"--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}"
        )

    return PipelineRun(
        results_dir=results_dir,
        stdout=proc.stdout,
        results_before=results_before,
        results_after=_snapshot(COMMITTED_RESULTS),
        pkg_before=pkg_before,
        pkg_after=_snapshot(_PKG_DIR, "*"),
    )


@pytest.fixture(scope="module")
def mini_run(tmp_path_factory, has_rscript) -> PipelineRun:
    """LEVEL 1: the whole pipeline over the 31-protein fixture. Runs once."""
    if not has_rscript:
        pytest.skip("Rscript not on PATH; the limma leg cannot run")
    results_dir = tmp_path_factory.mktemp("mini_results")
    return _run_foldchange(MINI_SHEET, results_dir, expect_baseline=False)


# ===========================================================================
# The fixture itself
# ===========================================================================

def test_committed_fixture_matches_its_generator():
    """``mini_sheet.xlsx`` must be exactly what ``make_mini_sheet.py`` produces.

    The workbook is committed so the tests need no generation step, which means
    it can silently fall behind the design table that every expectation below is
    read from. Comparing the frames (not the bytes) keeps this robust to
    openpyxl's own output changing between versions while still catching an
    edited-by-hand or stale fixture.
    """
    expected_L, expected_H = mini.build_frames()
    actual_L = pd.read_excel(MINI_SHEET, sheet_name=mini.SHEET_L)
    actual_H = pd.read_excel(MINI_SHEET, sheet_name=mini.SHEET_H)
    for name, expected, actual in (
        (mini.SHEET_L, expected_L, actual_L),
        (mini.SHEET_H, expected_H, actual_H),
    ):
        assert list(actual.columns) == list(expected.columns), name
        assert len(actual.columns) == 29, f"{name} is not the real 29-column layout"
        # The design table spells absence as ``None``; a blank cell read back
        # from .xlsx is ``NaN``. Normalise so the comparison is about content.
        pd.testing.assert_frame_equal(
            actual.reset_index(drop=True),
            expected.map(lambda v: np.nan if v is None else v).reset_index(drop=True),
            check_dtype=False,
            obj=f"{name} (regenerate with make_mini_sheet.py)",
        )


def test_committed_fixture_regenerates_byte_for_byte(tmp_path):
    """Regenerating the workbook must produce the committed bytes exactly.

    Stronger than the frame comparison above and it costs nothing: it proves the
    committed .xlsx really is the current output of the script, with no manual
    edit and no drift, and it keeps ``make_mini_sheet.write`` reproducible.

    Reproducibility here took work — an .xlsx hides two wall clocks (the zip
    member timestamps, and a ``dcterms:modified`` that openpyxl overwrites with
    ``now()`` on save, discarding the pinned document property). Both are
    flattened in ``_normalize_zip``. Without this test the pinning would rot
    silently, and the symptom would look like flakiness rather than a bug: the
    bytes only change when the wall clock crosses a second boundary.
    """
    regenerated = mini.write(tmp_path / "mini_sheet.xlsx")
    assert regenerated.read_bytes() == MINI_SHEET.read_bytes(), (
        "mini_sheet.xlsx is not what make_mini_sheet.py currently produces -- "
        "regenerate it (or, if the bytes differ while the frames match, a new "
        "openpyxl has introduced another unpinned timestamp)"
    )


def test_fixture_covers_every_branch():
    """The design table must actually reach every branch it claims to.

    A fixture that quietly lost its ON_OFF rows would leave the end-to-end test
    green and the branch untested -- the exact failure mode this whole module
    exists to prevent. Cheap insurance, and it needs no pipeline run.
    """
    both = mini.both_rows()
    regulated = {r["expect"]["regulated"] for r in both}
    assert regulated == {"UP", "DOWN", "NO CHANGE", "ON_OFF"}

    onoff = {r["expect"]["onoff"] for r in both}
    assert {"on_with_treatment", "off_with_treatment"} <= onoff, "ON_OFF is one-sided"

    assert any(r["expect"]["complete"] for r in both)
    assert any(not r["expect"]["complete"] for r in both)

    # absence encoded BOTH ways -- the committed workbook only has the zero form
    assert any(0 in (r["t1"], r["t2"], r["c1"], r["c2"]) for r in both), "no zero form"
    assert any(None in (r["t1"], r["t2"], r["c1"], r["c2"]) for r in both), "no blank form"

    # a zero denominator that the `complete` mask must keep away from a division
    assert any(
        r["c1"] == 0 and r["t1"] not in (0, None) and r["c2"] not in (0, None)
        for r in both
    ), "no zero-denominator row"

    # absent on both sides -> must stay NO CHANGE, never invented into ON_OFF
    assert any(
        (r["t1"], r["t2"], r["c1"], r["c2"]) == (0, 0, 0, 0)
        and r["expect"]["regulated"] == "NO CHANGE"
        for r in both
    ), "no fully-empty row"

    assert any(";" in r["acc"] for r in both), "no protein group on the merge path"
    assert any(";" in r["acc"] for r in mini.single_condition_rows()), \
        "no protein group on the single-condition path"
    assert any(r["gene"] is None for r in both), "no missing gene symbol"

    detected = {r["expect"]["detected_in"] for r in mini.single_condition_rows()}
    assert detected == {"control_only", "treated_only"}
    assert mini.quarantined_rows(), "no D11 junk accession"


# ===========================================================================
# LEVEL 1 -- the pipeline on new data
# ===========================================================================

@pytest.mark.r
def test_mini_run_writes_every_contracted_output(mini_run):
    missing = [rel for rel in _EXPECTED_OUTPUTS
               if not (mini_run.results_dir / rel).is_file()]
    assert not missing, f"foldchange.py did not write: {missing}"


@pytest.mark.r
def test_mini_class_counts_match_the_design(mini_run):
    """Every count, derived from the design table rather than typed in."""
    expected = mini.expected_counts()
    fc = mini_run.read("foldchange_all.csv")
    assert len(fc) == expected["foldchange_all_rows"]
    assert int(fc["complete"].sum()) == expected["complete_rows"]
    counts = fc["regulated"].value_counts()
    assert int(counts.get("UP", 0)) == expected["n_up"]
    assert int(counts.get("DOWN", 0)) == expected["n_down"]
    assert int(counts.get("NO CHANGE", 0)) == expected["n_nochange"]
    assert int(counts.get("ON_OFF", 0)) == expected["n_onoff"]
    assert len(mini_run.read("ipa_input.csv")) == expected["ipa_input_rows"]
    assert len(mini_run.read("qc_limma.csv")) == expected["qc_limma_rows"]
    assert len(mini_run.read("onoff_proteins.csv")) == expected["n_onoff"]
    assert (
        len(mini_run.read("single_condition_proteins.csv"))
        == expected["single_condition_rows"]
    )


@pytest.mark.r
@pytest.mark.parametrize(
    "row", mini.both_rows(), ids=[r["acc"] for r in mini.both_rows()]
)
def test_mini_row_lands_in_its_designed_branch(mini_run, row):
    """The heart of the regression net: one designed protein, one expectation.

    Parametrised per row so a failure names the protein and the branch rather
    than reporting "17 != 16" somewhere in an aggregate.
    """
    fc = mini_run.read("foldchange_all.csv").set_index("UniProt Accession Number")
    assert row["acc"] in fc.index, f"{row['acc']} missing from foldchange_all.csv"
    actual = fc.loc[row["acc"]]
    expect = row["expect"]
    why = f"{row['acc']} ({row['gene']}): {row['note']}"

    assert bool(actual["complete"]) == expect["complete"], f"complete -- {why}"
    assert actual["regulated"] == expect["regulated"], (
        f"regulated: expected {expect['regulated']!r}, got "
        f"{actual['regulated']!r} -- {why}"
    )
    onoff_actual = "" if pd.isna(actual["onoff"]) else actual["onoff"]
    assert onoff_actual == expect["onoff"], f"onoff -- {why}"

    if expect["log2FC"] is None:
        assert pd.isna(actual["log2FC"]), f"expected NaN log2FC -- {why}"
        return

    got = float(actual["log2FC"])
    want = expect["log2FC"]
    # Dyadic expectations (integer multiples of 0.5) are exactly representable
    # and every operation producing them is exact, so they are asserted to the
    # BIT. A tolerance there would hide exactly the systematic drift this file
    # exists to catch. The three log2(1.5) boundary rows are irrational and get
    # 1 ulp.
    if float(want * 2).is_integer():
        assert got == want, f"log2FC: expected exactly {want!r}, got {got!r} -- {why}"
    else:
        assert abs(got - want) <= 1e-15, f"log2FC: {got!r} vs {want!r} -- {why}"


@pytest.mark.r
def test_mini_direction_agrees_with_the_labels(mini_run):
    """Same D7 sign invariant the committed outputs are held to."""
    fc = mini_run.read("foldchange_all.csv")
    up = fc[fc["regulated"] == "UP"]["log2FC"]
    down = fc[fc["regulated"] == "DOWN"]["log2FC"]
    assert len(up) and len(down)
    assert (up > 0).all(), "an UP row has a non-positive log2FC"
    assert (down < 0).all(), "a DOWN row has a non-negative log2FC"


@pytest.mark.r
def test_mini_threshold_is_the_rounded_literal_not_log2_of_1_5(mini_run):
    """A 1.5-fold change is NOT called regulated, and that is on purpose.

    ``config/constants.LOG2_THRESHOLD`` is the rounded literal ``0.585``, which
    sits fractionally ABOVE ``log2(1.5) = 0.5849625``. So a protein up exactly
    1.5-fold falls just inside NO CHANGE. That is a real, load-bearing property
    of the classifier -- 'tidying' the literal into ``math.log2(1.5)`` would
    silently reclassify every protein sitting on the boundary. If this test
    fails, decide whether the change was intended and update the fixture's
    expectation deliberately; do not delete the row.
    """
    from proteomics_de.config.constants import LOG2_THRESHOLD

    assert LOG2_THRESHOLD > mini.L15, (
        "LOG2_THRESHOLD is no longer above log2(1.5); the boundary rows in "
        "mini_sheet.xlsx were designed around that and must be revisited."
    )
    fc = mini_run.read("foldchange_all.csv").set_index("UniProt Accession Number")
    # Bndry1: +1.5x in both replicates. Bndry3: the reciprocal. Both NO CHANGE.
    assert fc.loc["P10013", "regulated"] == "NO CHANGE"
    assert fc.loc["P10015", "regulated"] == "NO CHANGE"
    # Bndry2 is comfortably past the cutoff and must be called.
    assert fc.loc["P10014", "regulated"] == "UP"


@pytest.mark.r
def test_mini_ipa_input_is_the_complete_regulated_selection(mini_run):
    ipa = mini_run.read("ipa_input.csv")
    expected = [
        r["acc"] for r in mini.both_rows()
        if r["expect"]["complete"] and r["expect"]["regulated"] in ("UP", "DOWN")
    ]
    assert list(ipa["UniProt Accession Number"]) == expected
    assert set(ipa["regulated"]) <= {"UP", "DOWN"}
    assert ipa["log2FC"].notna().all()


@pytest.mark.r
def test_mini_onoff_file_carries_both_directions(mini_run):
    onoff = mini_run.read("onoff_proteins.csv").set_index("accession")
    expected = {
        r["acc"]: r["expect"]["onoff"]
        for r in mini.both_rows() if r["expect"]["regulated"] == "ON_OFF"
    }
    assert set(onoff.index) == set(expected)
    for acc, label in expected.items():
        assert onoff.loc[acc, "onoff"] == label, acc
    assert set(onoff["onoff"]) == {"on_with_treatment", "off_with_treatment"}


@pytest.mark.r
def test_mini_single_condition_labels_pin_the_d7_orientation(mini_run):
    """An L-sheet-only protein must come out ``treated_only``.

    This is the D7 correction seen from a completely different angle than the
    fold-change signs: 31578/31580 live in ``Protein Report L`` and are the
    TESTOSTERONE channels, so a protein detected only there was detected only in
    the treated condition. The label used to be hardcoded ``left_only ->
    control_only``, which was silently wrong for all 606 rescued proteins.
    """
    scp = mini_run.read("single_condition_proteins.csv").set_index("accession")
    for row in mini.single_condition_rows():
        assert row["acc"] in scp.index, row["acc"]
        assert scp.loc[row["acc"], "detected_in"] == row["expect"]["detected_in"], (
            f"{row['acc']}: {row['note']}"
        )
    # And the intensities land in the channels of the condition they came from.
    treated_only = scp[scp["detected_in"] == "treated_only"]
    assert treated_only[["Intensity 31578", "Intensity 31580"]].notna().all().all()
    assert treated_only[["Intensity 31579", "Intensity 31581"]].isna().all().all()


@pytest.mark.r
def test_mini_junk_accession_is_quarantined_and_the_group_survives(mini_run):
    """D11 on new data: token shape decides, never length."""
    quarantine = mini_run.read("qc/quarantine_accessions.csv")
    junk = [r["acc"] for r in mini.quarantined_rows()]
    assert list(quarantine["value"]) == junk
    assert set(quarantine["column"]) == {"accession"}

    scp = mini_run.read("single_condition_proteins.csv")
    assert not set(scp["accession"]) & set(junk), "a junk accession was shipped"
    # The legitimate semicolon group on the same path must NOT be swept up.
    groups = [r["acc"] for r in mini.single_condition_rows() if ";" in r["acc"]]
    assert groups, "the fixture lost its single-condition protein group"
    assert set(groups) <= set(scp["accession"]), (
        "a legitimate protein group was quarantined as junk -- the D11 "
        "discriminator has regressed from token shape to string length"
    )


@pytest.mark.r
def test_mini_limma_tests_exactly_the_eligible_rows(mini_run):
    """ON_OFF proteins are excluded; everything else in the table is tested."""
    qc = mini_run.read("qc_limma.csv")
    expected = [r["acc"] for r in mini.eligible_rows()]
    assert list(qc["id"]) == expected
    onoff_accs = {
        r["acc"] for r in mini.both_rows() if r["expect"]["regulated"] == "ON_OFF"
    }
    assert not set(qc["id"]) & onoff_accs


@pytest.mark.r
def test_mini_n_imputed_counts_the_values_r_invented(mini_run):
    """D10's load-bearing column, checked per row against the designed gaps.

    ``n_imputed`` is the only thing distinguishing a measured value from one
    MinProb drew at random, and at n=2 that distinction decides whether a
    protein's fold change is data or noise. It is counted in R on the log2
    matrix; the expectation here is counted by hand from the workbook.
    """
    qc = mini_run.read("qc_limma.csv").set_index("id")
    for row in mini.eligible_rows():
        expected = row["expect"]["n_imputed"]
        assert int(qc.loc[row["acc"], "n_imputed"]) == expected, (
            f"{row['acc']}: expected {expected} imputed value(s) -- {row['note']}"
        )
    assert set(qc["n_imputed"]) != {0}, "the fixture no longer exercises imputation"


@pytest.mark.r
def test_mini_pipeline_and_limma_log2fc_agree(mini_run):
    """The wiring check, on data the pipeline has never seen before.

    On the committed dataset this correlation is the evidence that the D7
    contrast is wired correctly. Re-running it on the fixture proves the wiring
    is a property of the *code*, not a coincidence of one workbook -- and here
    it is sharper, because every fully-observed row must match exactly (limma
    only diverges where it imputed).
    """
    fc = mini_run.read("foldchange_all.csv")
    qc = mini_run.read("qc_limma.csv")
    merged = fc[["UniProt Accession Number", "log2FC"]].merge(
        qc[["id", "limma_log2FC", "n_imputed"]],
        left_on="UniProt Accession Number", right_on="id", validate="one_to_one",
    )
    paired = merged.dropna(subset=["log2FC", "limma_log2FC"])
    assert len(paired) > 2
    r = float(np.corrcoef(paired["log2FC"], paired["limma_log2FC"])[0, 1])
    assert np.isfinite(r) and r > 0.9999, f"corr = {r!r} over {len(paired)} proteins"

    # Nothing was imputed on these rows, so the two computations are the same
    # arithmetic by two routes and must agree to qc_limma.csv's 6-dp rounding.
    fully_observed = paired[paired["n_imputed"] == 0]
    assert len(fully_observed) > 2
    assert np.allclose(
        fully_observed["log2FC"], fully_observed["limma_log2FC"], rtol=0, atol=5e-7
    ), "limma and the pipeline disagree on a protein with no imputed values"


@pytest.mark.r
def test_mini_ipa_significant_file_is_a_live_filter(mini_run):
    """Header-only again here, but for a reason that must be recomputable."""
    from proteomics_de.config.constants import ADJ_P_THRESHOLD

    sig = mini_run.read("ipa_input_significant.csv")
    qc = mini_run.read("qc_limma.csv")
    expected = qc[qc["regulated"].isin(["UP", "DOWN"]) & (qc["adj_p_value"] < ADJ_P_THRESHOLD)]
    assert len(sig) == len(expected)
    assert list(sig.columns) == [
        "UniProt Accession Number", "Gene names", "log2FC", "regulated", "adj_p_value",
    ]


@pytest.mark.r
def test_mini_boundary_record_follows_the_results_dir(mini_run):
    """``qc/qc_boundaries.json`` must land in --results-dir, not the package.

    This is the leak a sibling agent found: the boundary hooks take only
    ``(stage, df)``, so before ``boundaries.set_results_dir`` existed they wrote
    their QC record into the committed tree on every run, whatever
    ``--results-dir`` said. Standing test, not a one-off check.
    """
    import json

    record = json.loads(
        (mini_run.results_dir / "qc" / "qc_boundaries.json").read_text(encoding="utf-8")
    )
    stages = [r["stage"] for r in record["records"]]
    assert stages == ["after_load", "after_load", "after_merge", "after_foldchange"]
    # The permissive stages must have SEEN the junk accession and routed it.
    permissive = [r for r in record["records"] if r["policy"] == "permissive"]
    assert any(r["n_quarantinable_rows"] > 0 for r in permissive), (
        "no stage noticed the D11 junk accession"
    )
    # The strict stage must have passed: by then the junk rows are gone.
    strict = [r for r in record["records"] if r["policy"] == "strict"]
    assert strict and all(r["passed"] and not r["raised"] for r in strict)


@pytest.mark.r
def test_mini_de_contract_files_describe_the_same_run(mini_run):
    """``results/de/*`` must agree with ``qc_limma.csv`` row for row."""
    qc = mini_run.read("qc_limma.csv")
    matrix = mini_run.read("de/intensity_matrix.tsv")
    results = mini_run.read("de/limma_results.tsv")
    assert list(matrix["accession"]) == list(qc["id"])
    assert list(results["accession"]) == list(qc["id"])
    assert np.allclose(results["fold_change"], np.power(2.0, results["logFC"]), rtol=1e-12)

    design = mini_run.read("de/design.tsv")
    assert set(design["group"]) == {"control", "treated"}
    assert len(design) == 4


# ---------------------------------------------------------------------------
# The leak tests
# ---------------------------------------------------------------------------

@pytest.mark.r
def test_run_writes_nothing_into_the_committed_results_dir(mini_run):
    """``--results-dir <tmp>`` must leave ``proteomics_de/results/`` alone.

    A sibling agent found exactly this leak: the stage-boundary hooks wrote
    ``results/qc/qc_boundaries.json`` into the package's own results directory
    regardless of ``--results-dir``, so a test run silently edited the committed
    scientific record. It is fixed; this is the standing test that keeps it
    fixed. Size AND mtime, so a rewrite with identical bytes still counts as a
    write.
    """
    touched = _changed(mini_run.results_before, mini_run.results_after)
    assert not touched, (
        "running foldchange.py with --results-dir pointing elsewhere still "
        f"wrote into the committed proteomics_de/results/: {touched}"
    )


@pytest.mark.r
def test_run_touches_only_the_known_handoff_files_outside_results(mini_run):
    """The Python<->R handoff files leak into ``proteomics_de/`` — bounded here.

    ``limma_test.py`` resolves ``_limma_input.csv`` / ``_limma_output.csv`` /
    ``_limma_versions.txt`` / ``_limma_design.tsv`` against its OWN directory,
    not against ``--results-dir``, and the R worker is launched with ``cwd``
    there. So a run into a temp directory still overwrites files in the
    committed tree — and two of them (``_limma_input.csv``,
    ``_limma_output.csv``) are byte-frozen artifacts under
    ``tools/freeze.py``'s manifest. That is a real leak, of the same family as
    the ``qc_boundaries.json`` one above but one directory further out.

    The assertion is a SUBSET check on purpose:

    * fixing the leak (routing the handoff into ``--results-dir``) keeps this
      green — nothing changes, and the empty set is a subset;
    * the leak growing to any other file turns it red immediately.

    Not xfail: the current behaviour is bounded and understood, and an xfail
    would stop guarding the boundary while we wait for a fix.
    """
    touched = set(_changed(mini_run.pkg_before, mini_run.pkg_after))
    unexpected = sorted(touched - _LIMMA_HANDOFF_NAMES)
    assert not unexpected, (
        "the run wrote unexpected file(s) into the committed proteomics_de/ "
        f"directory: {unexpected}. Only the limma handoff intermediates "
        f"({sorted(_LIMMA_HANDOFF_NAMES)}) are known to leak there."
    )


# ===========================================================================
# LEVEL 2 -- byte-for-byte reproduction of the committed outputs
# ===========================================================================

@pytest.mark.slow
@pytest.mark.golden
@pytest.mark.r
def test_real_workbook_reproduces_the_committed_outputs(tmp_path, has_rscript):
    """Re-run the DE leg on the real data and canonically diff every artifact.

    ``test_freeze.py`` checks the files on disk against a stored manifest, which
    proves nobody edited them. This proves something else: that running the code
    *today* still produces them. A refactor can leave every committed file
    untouched (so the freeze gate stays green) while quietly changing what the
    pipeline would now generate; only re-running catches that.

    Scope is ``foldchange.py`` and the three stages it chains (centering,
    replicate correlation, limma). The viz and enrichment stages resolve their
    output paths from ``viz/style.py`` and cannot be redirected to a temp
    directory, and the enrichment stages make live network calls -- so they are
    out of scope here and stay covered by the freeze gate.

    LOCAL ONLY (``golden``). Different BLAS or R builds differ in the last
    float digit, so this is a this-machine guarantee. Do not "fix" a red CI by
    adding a tolerance; the correct response is to not run it in CI, which is
    what the workflow does.
    """
    if not has_rscript:
        pytest.skip("Rscript not on PATH; the limma leg cannot run")
    assert REAL_WORKBOOK.is_file(), f"input workbook missing: {REAL_WORKBOOK}"

    # Baseline assertions left ON: this IS the committed dataset, so
    # foldchange.py's own frozen-count checks should hold and are free evidence.
    run = _run_foldchange(REAL_WORKBOOK, tmp_path, expect_baseline=True)

    drifted, missing = [], []
    for rel in _EXPECTED_OUTPUTS:
        produced = run.results_dir / rel
        committed = COMMITTED_RESULTS / rel
        if not committed.is_file():
            missing.append(rel)
            continue
        if freeze.digest(produced)[0] != freeze.digest(committed)[0]:
            drifted.append(rel)

    assert not missing, f"no committed counterpart for: {missing}"
    assert not drifted, (
        "re-running the pipeline no longer reproduces the committed outputs: "
        f"{drifted}.\nThe committed files are unchanged (test_freeze.py covers "
        "that) -- what changed is what the code now produces. Investigate the "
        "diff; do not re-baseline until you know why."
    )

    # The WARN branch of the centering check writes an extra file on this
    # dataset (median log2FC = 0.2834, well outside the 0.10 tolerance). It is a
    # committed artifact, so its absence would be a silent regression.
    centered = run.results_dir / "foldchange_all_centered.csv"
    assert centered.is_file(), "the centering check no longer writes its WARN copy"
    assert (
        freeze.digest(centered)[0]
        == freeze.digest(COMMITTED_RESULTS / "foldchange_all_centered.csv")[0]
    )


@pytest.mark.slow
@pytest.mark.golden
@pytest.mark.r
def test_real_workbook_reproduces_the_frozen_limma_intermediates(tmp_path, has_rscript):
    """The R handoff files are frozen artifacts too, and must regenerate exactly.

    ``_limma_input.csv`` is what Python sent R and ``_limma_output.csv`` is what
    R sent back: together they are the audit trail for every p-value in the
    study. They are covered by the freeze manifest, so proving they *regenerate*
    identically is proving the Python/R boundary is deterministic — the same
    matrix, the same seed, the same MinProb draws.
    """
    if not has_rscript:
        pytest.skip("Rscript not on PATH; the limma leg cannot run")
    _run_foldchange(REAL_WORKBOOK, tmp_path, expect_baseline=True)

    manifest = freeze.read_manifest(FREEZE_MANIFEST)
    drifted = []
    for rel, (sha, _mode) in manifest.items():
        if "_limma_" not in rel:
            continue
        actual, _ = freeze.digest(_REPO_ROOT / rel)
        if actual != sha:
            drifted.append(rel)
    assert not drifted, (
        f"the regenerated limma intermediates differ from the manifest: {drifted}"
    )
