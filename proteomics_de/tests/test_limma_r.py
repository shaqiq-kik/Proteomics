"""Behavioural tests for the R worker, ``limma_test.R``.

Everything statistical runs against a **synthetic 40-protein fixture built in
tmp_path**, not the committed 1938-row data: it runs in milliseconds, and -- more
importantly -- it has a known ground truth, so these tests can assert what the
numbers should *be* rather than merely that they did not change.

The single most valuable assertion here is
:func:`test_contrast_direction_treated_minus_control`. Nothing else in the
repository checks the sign of the contrast, and coefficient 2 of
``model.matrix(~ factor(group, levels = c("control", "treated")))`` is
"treated - control" only because control sorts first. Invert that and every UP
becomes a DOWN, silently, in the fold-change table, the volcano plot, the
enrichment calls and the final report.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_REPO_ROOT = Path(__file__).resolve().parents[2]
_PKG_DIR = _REPO_ROOT / "proteomics_de"
for _p in (str(_REPO_ROOT), str(_PKG_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

R_SCRIPT = _PKG_DIR / "limma_test.R"
REAL_INPUT = _PKG_DIR / "_limma_input.csv"

# The fixture reuses the real handoff column names so that BOTH invocation forms
# (positional, which hardcodes them, and --design, which reads them) accept it.
HANDOFF_COLS = ["ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581"]
GROUPS = ["control", "control", "treated", "treated"]

N_PROTEINS = 40
N_PLANTED = 5           # proteins 0..4 carry a real, clean effect
PLANTED_LOG2FC = 3.0
BASE_LOG2 = 20.0        # a typical log2 intensity for this data
SD_PLANTED = 0.01       # tiny within-group variance -> unmistakable signal
SD_NOISE = 0.30

# Missing values at positions known to the test, all inside NOISE proteins so
# they cannot muddy the planted-effect recovery.
NA_POSITIONS = [(10, 0), (20, 3), (30, 1), (30, 2)]

_NO_RSCRIPT = shutil.which("Rscript") is None


def requires_r(fn):
    """Mark a test as needing the R toolchain, and skip it when absent."""
    fn = pytest.mark.skipif(_NO_RSCRIPT, reason="Rscript not on PATH")(fn)
    return pytest.mark.r(fn)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
def _synthetic_frame(with_nas=True):
    """40 proteins x 4 samples of RAW intensities with a known ground truth."""
    rng = np.random.default_rng(0)

    log2_vals = np.empty((N_PROTEINS, 4), dtype=float)
    for i in range(N_PROTEINS):
        if i < N_PLANTED:
            log2_vals[i, :2] = rng.normal(BASE_LOG2, SD_PLANTED, 2)
            log2_vals[i, 2:] = rng.normal(BASE_LOG2 + PLANTED_LOG2FC, SD_PLANTED, 2)
        else:
            log2_vals[i, :] = rng.normal(BASE_LOG2, SD_NOISE, 4)

    # The worker log2-transforms its input, so hand it 2**x.
    raw = np.power(2.0, log2_vals)
    if with_nas:
        for row, col in NA_POSITIONS:
            raw[row, col] = np.nan

    df = pd.DataFrame(
        {
            "id": [f"P{i:05d}" for i in range(N_PROTEINS)],
            "gene": [f"Gene{i}" for i in range(N_PROTEINS)],
        }
    )
    for j, col in enumerate(HANDOFF_COLS):
        df[col] = raw[:, j]
    return df


def _write_input(tmp_path, df=None, name="input.csv"):
    if df is None:
        df = _synthetic_frame()
    path = tmp_path / name
    df.to_csv(path, index=False, encoding="utf-8")
    return path


def _write_design(tmp_path, samples=HANDOFF_COLS, groups=GROUPS, name="design.tsv"):
    path = tmp_path / name
    pd.DataFrame({"sample": samples, "group": groups}).to_csv(
        path, sep="\t", index=False, lineterminator="\n"
    )
    return path


def _run_r(tmp_path, *args):
    return subprocess.run(
        ["Rscript", str(R_SCRIPT), *[str(a) for a in args]],
        capture_output=True, text=True, cwd=tmp_path,
    )


def _run_flags(tmp_path, inp, out, seed=42, mode="vanilla", design=None):
    args = ["--in", inp, "--out", out, "--seed", seed, "--mode", mode]
    if design is not None:
        args += ["--design", design]
    proc = _run_r(tmp_path, *args)
    assert proc.returncode == 0, f"R failed:\n{proc.stderr}"
    return pd.read_csv(out)


@pytest.fixture
def synthetic(tmp_path):
    """(input_csv, design_tsv) for the standard 40-protein fixture."""
    return _write_input(tmp_path), _write_design(tmp_path)


# ---------------------------------------------------------------------------
# THE contrast-direction test
# ---------------------------------------------------------------------------
@requires_r
def test_contrast_direction_treated_minus_control(tmp_path, synthetic):
    """treated >> control MUST produce a POSITIVE logFC.

    A sign flip here inverts every biological conclusion downstream, so this is
    asserted on the planted proteins (where treated is a clean +3 log2 above
    control) rather than inferred from the committed numbers.
    """
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    planted = res.iloc[:N_PLANTED]
    assert (planted["limma_log2FC"] > 0).all(), (
        "CONTRAST INVERTED: treated is 3 log2 units ABOVE control, but limma "
        f"reported negative logFC:\n{planted}"
    )
    assert planted["limma_log2FC"].min() > 2.0


@requires_r
def test_contrast_direction_flips_when_labels_are_swapped(tmp_path):
    """Mirror image: call the high samples 'control' and the sign must invert.

    Group labels stay in canonical order (control first), so the reference level
    is unchanged; only which samples wear which label changes.
    """
    inp = _write_input(tmp_path)
    relabelled = _write_design(
        tmp_path,
        samples=["trt_31579", "trt_31581", "ctrl_31578", "ctrl_31580"],
        groups=["control", "control", "treated", "treated"],
        name="relabelled.tsv",
    )
    res = _run_flags(tmp_path, inp, tmp_path / "relabelled_out.csv", design=relabelled)

    planted = res.iloc[:N_PLANTED]
    assert (planted["limma_log2FC"] < 0).all(), planted["limma_log2FC"].tolist()
    assert planted["limma_log2FC"].max() < -2.0


@requires_r
def test_reference_level_follows_design_row_order(tmp_path):
    """Listing treated FIRST makes it the reference, inverting every logFC.

    This is why ``config.design`` canonically sorts control-before-treated
    instead of trusting the sheet's row order: the reference level -- and
    therefore the sign of every fold change in the study -- is decided by
    whichever group appears first in the design file.
    """
    inp = _write_input(tmp_path)
    treated_first = _write_design(
        tmp_path,
        samples=["trt_31579", "trt_31581", "ctrl_31578", "ctrl_31580"],
        groups=["treated", "treated", "control", "control"],
        name="treated_first.tsv",
    )
    res = _run_flags(tmp_path, inp, tmp_path / "tf_out.csv", design=treated_first)

    planted = res.iloc[:N_PLANTED]
    assert (planted["limma_log2FC"] < 0).all(), planted["limma_log2FC"].tolist()


@requires_r
def test_relabelling_and_reordering_together_cancel_out(tmp_path):
    """The composition of the two flips above is a no-op -- a real invariant.

    Swapping the labels *and* the row order leaves the design matrix identical,
    so the committed numbers are reproduced exactly.
    """
    inp = _write_input(tmp_path)
    both = _write_design(
        tmp_path,
        groups=["treated", "treated", "control", "control"],
        name="both.tsv",
    )
    canonical = _write_design(tmp_path, name="canonical.tsv")

    a = tmp_path / "both_out.csv"
    b = tmp_path / "canonical_out.csv"
    _run_flags(tmp_path, inp, a, design=both)
    _run_flags(tmp_path, inp, b, design=canonical)

    assert a.read_bytes() == b.read_bytes()


# ---------------------------------------------------------------------------
# Recovery of the planted truth
# ---------------------------------------------------------------------------
@requires_r
def test_planted_proteins_recover_logfc(tmp_path, synthetic):
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    recovered = res.iloc[:N_PLANTED]["limma_log2FC"]
    assert np.allclose(recovered, PLANTED_LOG2FC, atol=0.1), recovered.tolist()


@requires_r
def test_planted_proteins_rank_top_five_by_pvalue(tmp_path, synthetic):
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    top5 = set(res["p_value"].nsmallest(N_PLANTED).index)
    assert top5 == set(range(N_PLANTED)), (
        f"expected the planted proteins to rank top-5; got rows {sorted(top5)}"
    )


@requires_r
def test_noise_proteins_have_no_systematic_shift(tmp_path, synthetic):
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    noise = res.iloc[N_PLANTED:]["limma_log2FC"]
    assert abs(noise.mean()) < 0.5, noise.mean()


@requires_r
def test_id_order_is_preserved(tmp_path, synthetic):
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)
    assert res["id"].tolist() == pd.read_csv(inp)["id"].tolist()


@requires_r
def test_adjusted_pvalues_in_range_and_bh_monotone(tmp_path, synthetic):
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    for col in ("p_value", "adj_p_value"):
        assert res[col].between(0.0, 1.0).all(), res[col]
    # BH never lowers a p-value, and is monotone in raw-p order.
    assert (res["adj_p_value"] >= res["p_value"] - 1e-12).all()
    ordered = res.sort_values("p_value")["adj_p_value"].to_numpy()
    assert np.all(np.diff(ordered) >= -1e-12), "BH adjusted p-values not monotone"


# ---------------------------------------------------------------------------
# Determinism / the seed
# ---------------------------------------------------------------------------
@requires_r
def test_same_seed_is_byte_identical(tmp_path, synthetic):
    inp, dsn = synthetic
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    _run_flags(tmp_path, inp, a, seed=42, design=dsn)
    _run_flags(tmp_path, inp, b, seed=42, design=dsn)
    assert a.read_bytes() == b.read_bytes()


@requires_r
def test_different_seeds_differ_when_values_are_missing(tmp_path, synthetic):
    """Proves MinProb is genuinely stochastic and set.seed() is load-bearing."""
    inp, dsn = synthetic
    a, b = tmp_path / "s42.csv", tmp_path / "s7.csv"
    _run_flags(tmp_path, inp, a, seed=42, design=dsn)
    _run_flags(tmp_path, inp, b, seed=7, design=dsn)
    assert a.read_bytes() != b.read_bytes(), (
        "seed 42 and seed 7 produced identical output on a matrix WITH missing "
        "values -- imputation is not actually being seeded"
    )


@requires_r
def test_seed_is_irrelevant_without_missing_values(tmp_path):
    """The control for the test above: no NAs, nothing to impute, no randomness."""
    inp = _write_input(tmp_path, _synthetic_frame(with_nas=False), name="complete.csv")
    dsn = _write_design(tmp_path)
    a, b = tmp_path / "c42.csv", tmp_path / "c7.csv"
    _run_flags(tmp_path, inp, a, seed=42, design=dsn)
    _run_flags(tmp_path, inp, b, seed=7, design=dsn)
    assert a.read_bytes() == b.read_bytes()


# ---------------------------------------------------------------------------
# --design equivalence: the core proof of this package
# ---------------------------------------------------------------------------
@requires_r
def test_design_flag_matches_positional_on_synthetic(tmp_path, synthetic):
    inp, dsn = synthetic
    pos, flag = tmp_path / "pos.csv", tmp_path / "flag.csv"

    proc = _run_r(tmp_path, inp, pos, 42, "vanilla")
    assert proc.returncode == 0, proc.stderr
    _run_flags(tmp_path, inp, flag, design=dsn)

    assert pos.read_bytes() == flag.read_bytes()


@requires_r
@pytest.mark.skipif(not REAL_INPUT.exists(), reason="_limma_input.csv not present")
def test_design_flag_matches_positional_on_real_input(tmp_path):
    """The committed 1938-protein handoff: --design must change nothing at all."""
    dsn = _write_design(tmp_path)
    pos, flag = tmp_path / "pos.csv", tmp_path / "flag.csv"

    proc = _run_r(tmp_path, REAL_INPUT, pos, 42, "vanilla")
    assert proc.returncode == 0, proc.stderr
    _run_flags(tmp_path, REAL_INPUT, flag, design=dsn)

    assert pos.read_bytes() == flag.read_bytes(), (
        "--design perturbed the design matrix on the real input"
    )


@requires_r
def test_design_flag_reproduces_committed_output(tmp_path):
    """And it still reproduces the frozen _limma_output.csv byte-for-byte."""
    committed = _PKG_DIR / "_limma_output.csv"
    if not (REAL_INPUT.exists() and committed.exists()):
        pytest.skip("committed limma intermediates not present")

    dsn = _write_design(tmp_path)
    out = tmp_path / "_limma_output.csv"
    _run_flags(tmp_path, REAL_INPUT, out, design=dsn)

    assert out.read_bytes() == committed.read_bytes()


# ---------------------------------------------------------------------------
# Fail-loud contract (no limma work needed -- these are fast)
# ---------------------------------------------------------------------------
def _assert_bug7_error(proc):
    assert proc.returncode != 0, f"expected nonzero exit, got 0\n{proc.stdout}"
    assert proc.stderr.startswith("BUG7 R ERROR"), (
        f"stderr must start with 'BUG7 R ERROR', got: {proc.stderr!r}"
    )


@requires_r
def test_missing_column_fails_loudly(tmp_path):
    df = _synthetic_frame().drop(columns=[HANDOFF_COLS[-1]])
    inp = _write_input(tmp_path, df, name="short.csv")
    _assert_bug7_error(_run_r(tmp_path, inp, tmp_path / "o.csv", 42))


@requires_r
def test_zero_rows_fails_loudly(tmp_path):
    inp = _write_input(tmp_path, _synthetic_frame().iloc[0:0], name="empty.csv")
    _assert_bug7_error(_run_r(tmp_path, inp, tmp_path / "o.csv", 42))


@requires_r
def test_non_integer_seed_fails_loudly(tmp_path, synthetic):
    inp, _ = synthetic
    _assert_bug7_error(_run_r(tmp_path, inp, tmp_path / "o.csv", "not-an-int"))


@requires_r
def test_unknown_mode_fails_loudly(tmp_path, synthetic):
    inp, _ = synthetic
    _assert_bug7_error(_run_r(tmp_path, inp, tmp_path / "o.csv", 42, "sideways"))


@requires_r
def test_unknown_flag_fails_loudly(tmp_path, synthetic):
    inp, _ = synthetic
    _assert_bug7_error(
        _run_r(tmp_path, "--in", inp, "--out", tmp_path / "o.csv",
               "--seed", 42, "--bogus", "x")
    )


@requires_r
def test_missing_design_file_fails_loudly(tmp_path, synthetic):
    inp, _ = synthetic
    _assert_bug7_error(
        _run_r(tmp_path, "--in", inp, "--out", tmp_path / "o.csv",
               "--seed", 42, "--design", tmp_path / "nope.tsv")
    )


@requires_r
def test_design_without_group_column_fails_loudly(tmp_path, synthetic):
    inp, _ = synthetic
    bad = tmp_path / "bad_design.tsv"
    pd.DataFrame({"sample": HANDOFF_COLS}).to_csv(bad, sep="\t", index=False)
    _assert_bug7_error(
        _run_r(tmp_path, "--in", inp, "--out", tmp_path / "o.csv",
               "--seed", 42, "--design", bad)
    )


@requires_r
def test_no_output_written_on_failure(tmp_path, synthetic):
    """A failed run must leave nothing behind that could look like a result."""
    inp, _ = synthetic
    out = tmp_path / "should_not_exist.csv"
    _assert_bug7_error(_run_r(tmp_path, inp, out, "not-an-int"))
    assert not out.exists()


# ---------------------------------------------------------------------------
# Python side of the fail-loud contract (no R needed)
# ---------------------------------------------------------------------------
def test_python_raises_and_writes_nothing_when_r_fails(tmp_path, monkeypatch):
    """run_limma_test must raise RuntimeError and write NO result files."""
    import limma_test

    fc_csv = _PKG_DIR / "results" / "foldchange_all.csv"
    if not fc_csv.exists():
        pytest.skip("foldchange_all.csv not present")

    work = tmp_path / "work"
    work.mkdir()
    outdir = tmp_path / "out"

    # Redirect the module's own directory so the real intermediates are untouched.
    monkeypatch.setattr(limma_test, "_HERE", str(work))

    def fake_run(*args, **kwargs):
        return subprocess.CompletedProcess(
            args=args[0] if args else [], returncode=1,
            stdout="", stderr="BUG7 R ERROR: simulated failure",
        )

    monkeypatch.setattr(limma_test.subprocess, "run", fake_run)

    with pytest.raises(RuntimeError, match="BUG7"):
        limma_test.run_limma_test(foldchange_csv=str(fc_csv), outdir=str(outdir))

    assert not (outdir / "qc_limma.csv").exists()
    assert not (outdir / "ipa_input_significant.csv").exists()


def test_python_raises_when_rscript_is_absent(tmp_path, monkeypatch):
    import limma_test

    fc_csv = _PKG_DIR / "results" / "foldchange_all.csv"
    if not fc_csv.exists():
        pytest.skip("foldchange_all.csv not present")

    work = tmp_path / "work"
    work.mkdir()
    outdir = tmp_path / "out"
    monkeypatch.setattr(limma_test, "_HERE", str(work))

    def boom(*args, **kwargs):
        raise FileNotFoundError("Rscript")

    monkeypatch.setattr(limma_test.subprocess, "run", boom)

    with pytest.raises(RuntimeError, match="Rscript"):
        limma_test.run_limma_test(foldchange_csv=str(fc_csv), outdir=str(outdir))

    assert not (outdir / "qc_limma.csv").exists()


def test_run_limma_test_signature_is_stable():
    """foldchange.py calls run_limma_test() with ZERO arguments -- keep it so."""
    import inspect

    import limma_test

    sig = inspect.signature(limma_test.run_limma_test)
    assert all(
        p.default is not inspect.Parameter.empty for p in sig.parameters.values()
    ), "run_limma_test() must remain callable with no arguments"
    assert list(sig.parameters) == [
        "foldchange_csv", "outdir", "ebayes_mode", "qc_filename",
        "output_name", "reuse_input", "write_ipa",
    ]
