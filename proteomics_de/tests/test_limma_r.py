"""Behavioural tests for the R worker, ``limma_test.R``.

Everything statistical runs against a **synthetic 40-protein fixture built in
tmp_path**, not the committed 1938-row data: it runs in milliseconds, and -- more
importantly -- it has a known ground truth, so these tests can assert what the
numbers should *be* rather than merely that they did not change.

The single most valuable assertion here is
:func:`test_contrast_direction_treated_minus_control`. Nothing else in the
repository checks the sign of the contrast, and coefficient 2 of
``model.matrix(~ factor(group, levels = group_levels))`` is "treated - control"
only because ``"control"`` is level 1. Invert that and every UP becomes a DOWN,
silently, in the fold-change table, the volcano plot, the enrichment calls and
the final report.

That reference level used to be ``unique(group)`` -- order of first appearance --
which made the design file's ROW ORDER decide the sign of every logFC in the
study with nothing to signal it. DECISIONS_LOG D7 removed that landmine: the
design file is now written in acquisition order (so MinProb imputation stays
invariant under relabelling), and ``"control"`` is pinned as the reference BY
NAME. Both halves of the new contract are asserted below --
:func:`test_reference_level_is_control_by_name_not_row_order`,
:func:`test_relabelling_negates_logfc_and_leaves_pvalues_invariant`, and
:func:`test_design_without_a_control_group_fails_loudly`.
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
REAL_DESIGN = _PKG_DIR / "_limma_design.tsv"
COMMITTED_OUTPUT = _PKG_DIR / "_limma_output.csv"

# The layout the POSITIONAL invocation hardcodes (limma_test.R:130). It is the
# original 2x2 form, kept working byte-for-byte, so it still carries the pre-D7
# naming: columns 1-2 named ctrl_*, columns 3-4 named trt_*.
#
# The synthetic fixture below is built with these names and this grouping, so
# BOTH invocation forms accept it -- positional (which hardcodes them) and
# --design (which reads them) -- which is what makes the equivalence tests
# possible. Since DECISIONS_LOG D7 the REAL handoff no longer looks like this:
# limma_test.py hands R the matrix in acquisition order, i.e.
# trt_31578, trt_31580, ctrl_31579, ctrl_31581. See REAL_HANDOFF_COLS.
HANDOFF_COLS = ["ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581"]
GROUPS = ["control", "control", "treated", "treated"]

# The real handoff, D7-corrected: acquisition order, with each column's group in
# its name. Must agree with the committed _limma_input.csv / _limma_design.tsv.
REAL_HANDOFF_COLS = ["trt_31578", "trt_31580", "ctrl_31579", "ctrl_31581"]
REAL_GROUPS = ["treated", "treated", "control", "control"]

N_PROTEINS = 40
N_PLANTED = 5           # proteins 0..4 carry a real, clean effect
PLANTED_LOG2FC = 3.0
BASE_LOG2 = 20.0        # a typical log2 intensity for this data
SD_PLANTED = 0.01       # tiny within-group variance -> unmistakable signal
SD_NOISE = 0.30

# Missing values at positions known to the test, all inside NOISE proteins so
# they cannot muddy the planted-effect recovery. Row 30 deliberately carries
# TWO, so n_imputed has to distinguish 0 / 1 / 2 rather than just "any".
NA_POSITIONS = [(10, 0), (20, 3), (30, 1), (30, 2)]

#: The eBayes flavour the worker uses when --mode is absent. DECISIONS_LOG D9
#: flipped this from "vanilla" to "trend_robust"; the tests below run the
#: DEFAULT unless they are specifically about the other flavour, so the path the
#: pipeline actually takes is the path under test.
DEFAULT_MODE = "trend_robust"

#: The R worker's output schema. The first four are the original contract and
#: are position-pinned; D10 appends the rest (research1.md line 169).
R_OUTPUT_ORIGINAL_COLUMNS = ["id", "limma_log2FC", "p_value", "adj_p_value"]
R_OUTPUT_D10_COLUMNS = ["n_imputed", "AveExpr", "t", "B"]

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


def _run_flags(tmp_path, inp, out, seed=42, mode=DEFAULT_MODE, design=None):
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
def test_reference_level_is_control_by_name_not_row_order(tmp_path):
    """Design-file ROW ORDER must not touch the sign of a single logFC. (D7)

    The worker used to derive its factor levels as ``unique(group)``, so the
    group that happened to appear first in the design file became the reference
    -- a silent landmine, since the sign of every fold change in the study rode
    on the order of four lines in a TSV. That behaviour was deliberately
    REMOVED: ``"control"`` is now level 1 by name.

    It matters because the design file is no longer written control-first. It is
    written in acquisition order, which for today's D7-corrected sheet genuinely
    puts *treated* first (REAL_GROUPS above).

    Run with no missing values so imputation cannot contribute any difference:
    reordering the rows permutes the matrix columns, and MinProb is stochastic
    per column.
    """
    inp = _write_input(tmp_path, _synthetic_frame(with_nas=False), name="complete.csv")
    control_first = _write_design(tmp_path, name="control_first.tsv")
    treated_first = _write_design(
        tmp_path,
        samples=["trt_31579", "trt_31581", "ctrl_31578", "ctrl_31580"],
        groups=["treated", "treated", "control", "control"],
        name="treated_first.tsv",
    )

    cf = _run_flags(tmp_path, inp, tmp_path / "cf_out.csv", design=control_first)
    tf = _run_flags(tmp_path, inp, tmp_path / "tf_out.csv", design=treated_first)

    # Same sample->group mapping, treated listed first: still treated - control.
    planted = tf.iloc[:N_PLANTED]["limma_log2FC"]
    assert (planted > 0).all(), (
        "reference level followed design-file row order: listing treated first "
        f"inverted the contrast. planted logFC = {planted.tolist()}"
    )
    assert planted.min() > 2.0

    # And not merely the sign -- the whole result is unmoved by row order.
    assert np.array_equal(tf["id"], cf["id"])
    assert np.allclose(tf["limma_log2FC"], cf["limma_log2FC"], rtol=0, atol=1e-12)
    assert np.allclose(tf["p_value"], cf["p_value"], rtol=0, atol=1e-12)
    assert np.allclose(tf["adj_p_value"], cf["adj_p_value"], rtol=0, atol=1e-12)


@requires_r
def test_relabelling_negates_logfc_and_leaves_pvalues_invariant(tmp_path):
    """Swapping the two labels negates logFC and moves no p-value at all. (D7)

    This is the property the D7 correction rests on: swapping the levels of a
    two-group contrast negates logFC and t, but leaves |t| -- and therefore p
    and adj-p -- untouched. Verified on the real data when D7 landed (limma
    logFC: max |new + old| = 0.000e+00; p_value invariant to ~1e-15).

    It only holds because the matrix handed to R stays in acquisition order:
    both runs below pass the SAME sample order and differ only in the `group`
    column, so MinProb draws the same values for the same samples. Build the
    input WITH missing values, so that is actually exercised.
    """
    inp = _write_input(tmp_path)
    canonical = _write_design(tmp_path, name="canonical.tsv")
    relabelled = _write_design(
        tmp_path,
        groups=["treated", "treated", "control", "control"],
        name="relabelled.tsv",
    )

    a = _run_flags(tmp_path, inp, tmp_path / "canonical_out.csv", design=canonical)
    b = _run_flags(tmp_path, inp, tmp_path / "relabelled_out.csv", design=relabelled)

    assert np.array_equal(b["id"], a["id"])
    assert np.allclose(b["limma_log2FC"], -a["limma_log2FC"], rtol=0, atol=1e-12)
    assert np.allclose(b["p_value"], a["p_value"], rtol=0, atol=1e-12)
    assert np.allclose(b["adj_p_value"], a["adj_p_value"], rtol=0, atol=1e-12)
    # Sanity: the planted effect really did change sign, it was not ~0 already.
    assert (a["limma_log2FC"].iloc[:N_PLANTED] > 2.0).all()


@requires_r
def test_design_without_a_control_group_fails_loudly(tmp_path, synthetic):
    """No `control` group means no reference level to pin -- refuse to guess.

    With the levels no longer taken from order of appearance, a design file that
    never says "control" has nothing to anchor the sign of the contrast. Falling
    back to alphabetical (or first-seen) order is exactly the silent behaviour
    D7 removed, so the worker must abort instead.
    """
    inp, _ = synthetic
    no_control = _write_design(
        tmp_path,
        groups=["vehicle", "vehicle", "testosterone", "testosterone"],
        name="no_control.tsv",
    )
    proc = _run_r(tmp_path, "--in", inp, "--out", tmp_path / "o.csv",
                  "--seed", 42, "--design", no_control)
    _assert_bug7_error(proc)
    assert "control" in proc.stderr
    assert not (tmp_path / "o.csv").exists()


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
def test_positional_invocation_is_refused(tmp_path, synthetic):
    """The design must come from the sample sheet or not at all.

    limma_test.R used to fall back to an inline
    ``ctrl_31578/ctrl_31580/trt_31579/trt_31581`` layout with groups
    ``control, control, treated, treated`` when --design was absent. Per
    DECISIONS_LOG D7 that assignment is INVERTED: 31578/31580 are testosterone,
    31579/31581 are vehicle. Fed matching data the fallback would have silently
    negated every logFC in the study.

    It was removed rather than corrected, because a second hardcoded source of
    truth for the design is precisely how the original error survived. This test
    pins that: no --design is a loud, immediate failure.
    """
    inp, _dsn = synthetic
    proc = _run_r(tmp_path, inp, tmp_path / "pos.csv", 42, "vanilla")
    assert proc.returncode != 0, "positional invocation must fail, not guess a design"
    assert "BUG7 R ERROR" in proc.stderr, proc.stderr
    assert "--design is required" in proc.stderr, proc.stderr
    assert not (tmp_path / "pos.csv").exists(), "no output may be written on refusal"


@requires_r
@pytest.mark.skipif(not REAL_INPUT.exists(), reason="_limma_input.csv not present")
def test_design_drives_the_real_input_end_to_end(tmp_path):
    """1938 real proteins, real missingness, design read from the committed TSV.

    Replaces the old positional-vs---design equivalence check: with the fallback
    gone there is exactly one code path, so what matters now is that the
    committed design file is the one actually in force and that it reproduces
    the committed output byte-for-byte.
    """
    real = pd.read_csv(REAL_INPUT)
    assert list(real.columns) == ["id", "gene"] + REAL_HANDOFF_COLS, (
        f"committed handoff columns changed: {list(real.columns)}"
    )
    out = tmp_path / "out.csv"
    _run_flags(tmp_path, REAL_INPUT, out, design=REAL_DESIGN)
    assert out.read_bytes() == COMMITTED_OUTPUT.read_bytes(), (
        "re-running the committed design did not reproduce _limma_output.csv"
    )



@requires_r
def test_design_flag_reproduces_committed_output(tmp_path):
    """The committed _limma_design.tsv + _limma_input.csv reproduce _limma_output.csv.

    Byte-for-byte, so this also pins the D7-corrected sign: the committed output
    is the post-flip one (limma logFC negated relative to the pre-D7 artifact,
    p-values unchanged).
    """
    committed = _PKG_DIR / "_limma_output.csv"
    if not (REAL_INPUT.exists() and committed.exists() and REAL_DESIGN.exists()):
        pytest.skip("committed limma intermediates not present")

    # The design as the pipeline actually writes it: acquisition order, treated
    # first. Read here rather than reconstructed, so a drift shows up as a diff.
    dsn = pd.read_csv(REAL_DESIGN, sep="\t", dtype=str)
    assert dsn["sample"].tolist() == REAL_HANDOFF_COLS
    assert dsn["group"].tolist() == REAL_GROUPS

    out = tmp_path / "_limma_output.csv"
    _run_flags(tmp_path, REAL_INPUT, out, design=REAL_DESIGN)

    assert out.read_bytes() == committed.read_bytes()


# ---------------------------------------------------------------------------
# D9 -- eBayes(trend=TRUE, robust=TRUE) is the default flavour
# ---------------------------------------------------------------------------
@requires_r
def test_omitting_mode_gives_the_trend_robust_default(tmp_path, synthetic):
    """No --mode must mean trend/robust, not vanilla. (D9)

    The worker shipped with "vanilla" as its default because that was the
    byte-reproducible baseline. research1.md line 124 specifies
    ``eBayes(trend=TRUE, robust=TRUE)`` and it is the proteomics field standard,
    so D9 flipped the default. Pinned here because the default is what
    ``limma_test.py`` -- and therefore the whole pipeline -- actually gets.
    """
    inp, dsn = synthetic
    args = ["--in", inp, "--out", tmp_path / "implicit.csv", "--seed", 42,
            "--design", dsn]
    proc = _run_r(tmp_path, *args)
    assert proc.returncode == 0, f"R failed:\n{proc.stderr}"
    implicit = pd.read_csv(tmp_path / "implicit.csv")

    explicit_tr = _run_flags(tmp_path, inp, tmp_path / "tr.csv",
                             mode="trend_robust", design=dsn)
    explicit_van = _run_flags(tmp_path, inp, tmp_path / "van.csv",
                              mode="vanilla", design=dsn)

    assert np.allclose(implicit["p_value"], explicit_tr["p_value"], rtol=0, atol=0), (
        "omitting --mode did not run trend/robust"
    )
    assert not np.array_equal(implicit["p_value"], explicit_van["p_value"]), (
        "trend/robust and vanilla produced identical p-values -- the mode "
        "toggle is not reaching eBayes"
    )


@requires_r
def test_the_two_flavours_differ_only_in_the_variance_model(tmp_path, synthetic):
    """logFC bit-identical, p-values not. (D9)

    ``eBayes`` moderates residual variances; it does not refit the linear model,
    so the fitted coefficients -- and therefore every fold change -- are
    untouched. This is the entire argument for changing the default without
    re-litigating any biological conclusion, so it is asserted rather than
    assumed. Run on the synthetic fixture, where the ground truth is known, in
    addition to the real-data check in test_limma_contract.py.
    """
    inp, dsn = synthetic
    van = _run_flags(tmp_path, inp, tmp_path / "v.csv", mode="vanilla", design=dsn)
    tr = _run_flags(tmp_path, inp, tmp_path / "t.csv", mode="trend_robust", design=dsn)

    assert van["id"].tolist() == tr["id"].tolist()
    assert np.array_equal(
        van["limma_log2FC"].to_numpy(), tr["limma_log2FC"].to_numpy()
    ), (
        "eBayes moved logFC between flavours; max |diff| = "
        f"{np.abs(van['limma_log2FC'] - tr['limma_log2FC']).max():.3e}"
    )
    # AveExpr is a property of the data, not of the moderation, so it must not
    # move either. B and t are moderation outputs and are expected to.
    assert np.array_equal(van["AveExpr"].to_numpy(), tr["AveExpr"].to_numpy())
    assert np.array_equal(van["n_imputed"].to_numpy(), tr["n_imputed"].to_numpy())
    assert not np.array_equal(van["p_value"].to_numpy(), tr["p_value"].to_numpy())


@requires_r
def test_both_flavours_still_recover_the_planted_effect(tmp_path, synthetic):
    """Changing the default must not cost the worker its known-truth recovery."""
    inp, dsn = synthetic
    for mode in ("vanilla", "trend_robust"):
        res = _run_flags(tmp_path, inp, tmp_path / f"{mode}.csv", mode=mode,
                         design=dsn)
        planted = res.iloc[:N_PLANTED]
        assert np.allclose(planted["limma_log2FC"], PLANTED_LOG2FC, atol=0.1), (
            f"{mode}: {planted['limma_log2FC'].tolist()}"
        )
        top5 = set(res["p_value"].nsmallest(N_PLANTED).index)
        assert top5 == set(range(N_PLANTED)), f"{mode}: ranked rows {sorted(top5)}"


@requires_r
@pytest.mark.skipif(not REAL_INPUT.exists(), reason="_limma_input.csv not present")
def test_vanilla_companion_reproduces_its_committed_intermediate(tmp_path):
    """results/qc_limma_vanilla.csv's worker output is reproducible. (D9)

    D9 keeps vanilla alive so both flavours stay comparable. A preserved file
    nobody can regenerate is not a baseline, so the committed
    ``_limma_output_vanilla.csv`` is re-derived here from the same input, seed
    and design -- byte-for-byte.
    """
    committed = _PKG_DIR / "_limma_output_vanilla.csv"
    if not (committed.exists() and REAL_DESIGN.exists()):
        pytest.skip("committed vanilla intermediate not present")

    out = tmp_path / "_limma_output_vanilla.csv"
    _run_flags(tmp_path, REAL_INPUT, out, mode="vanilla", design=REAL_DESIGN)
    assert out.read_bytes() == committed.read_bytes()


# ---------------------------------------------------------------------------
# D10 -- the restored output columns
# ---------------------------------------------------------------------------
@requires_r
def test_output_columns_are_the_original_four_plus_the_d10_four(tmp_path, synthetic):
    """D10 APPENDS to the worker's contract; it must not reorder or rename."""
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)
    assert list(res.columns) == R_OUTPUT_ORIGINAL_COLUMNS + R_OUTPUT_D10_COLUMNS


@requires_r
def test_n_imputed_counts_the_planted_missing_values(tmp_path, synthetic):
    """Ground truth: the fixture knows exactly which cells it blanked.

    ``n_imputed`` must be counted on the log2 matrix BEFORE ``impute.MinProb``
    runs -- one line later every cell is finite and the distinction between a
    measured value and an invented one no longer exists anywhere. This is the
    column's whole reason to exist at n=2 with a stochastic imputer, so it is
    checked against the fixture's own NA_POSITIONS rather than range-checked.
    """
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    expected = [0] * N_PROTEINS
    for row, _col in NA_POSITIONS:
        expected[row] += 1
    assert res["n_imputed"].tolist() == expected
    # The fixture plants a 2-NA row on purpose; if that ever stops being true
    # the test degenerates into "0 or 1" and stops proving anything.
    assert max(expected) == 2
    assert res["n_imputed"].between(0, 4).all()


@requires_r
def test_n_imputed_is_all_zero_when_nothing_is_missing(tmp_path):
    """The control: no NAs in, no imputation claimed."""
    inp = _write_input(tmp_path, _synthetic_frame(with_nas=False), name="complete.csv")
    dsn = _write_design(tmp_path)
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)
    assert (res["n_imputed"] == 0).all()


@requires_r
def test_n_imputed_survives_a_zero_intensity_not_just_a_blank(tmp_path):
    """<= 0 is missing too, and must be counted as such.

    ``etl.build_matrix.intensity_series`` treats a 0 intensity as "below the
    detection limit" (MNAR), identically to a blank cell; the R worker maps
    ``<= 0 -> NA`` for the same reason. A zero that slipped through as a
    measurement would log2 to -Inf and be counted as observed.
    """
    df = _synthetic_frame(with_nas=False)
    df.loc[3, HANDOFF_COLS[1]] = 0.0
    df.loc[7, HANDOFF_COLS[0]] = 0.0
    df.loc[7, HANDOFF_COLS[2]] = 0.0
    inp = _write_input(tmp_path, df, name="zeros.csv")
    dsn = _write_design(tmp_path)
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    assert res["n_imputed"].iloc[3] == 1
    assert res["n_imputed"].iloc[7] == 2
    assert res["n_imputed"].drop(index=[3, 7]).eq(0).all()


@requires_r
def test_avexpr_t_and_b_are_finite_and_consistent(tmp_path, synthetic):
    inp, dsn = synthetic
    res = _run_flags(tmp_path, inp, tmp_path / "out.csv", design=dsn)

    for col in ("AveExpr", "t", "B"):
        assert np.isfinite(res[col]).all(), f"{col} has non-finite entries"
    # AveExpr is the row mean of the imputed log2 matrix; the fixture sits at
    # ~BASE_LOG2, with the planted rows pulled up by half the planted effect.
    assert np.allclose(
        res["AveExpr"].iloc[N_PLANTED:], BASE_LOG2, atol=1.0
    ), res["AveExpr"].iloc[N_PLANTED:].tolist()
    # t and logFC are the same contrast: the signs cannot disagree.
    same_sign = np.sign(res["t"]) == np.sign(res["limma_log2FC"])
    assert same_sign[res["limma_log2FC"] != 0].all()
    # B (log-odds of differential expression) must rank the planted effect top.
    assert set(res["B"].nlargest(N_PLANTED).index) == set(range(N_PLANTED))


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
        # D9/D10 additions; both default to the pipeline's behaviour, so
        # foldchange.py's zero-argument call keeps producing the full output set.
        "write_contract", "vanilla_companion",
    ]
    assert sig.parameters["ebayes_mode"].default == "trend_robust", (
        "DECISIONS_LOG D9: the default eBayes flavour must be trend/robust"
    )
    assert sig.parameters["vanilla_companion"].default is True


def test_default_ebayes_mode_agrees_between_python_and_r():
    """One default, two languages -- assert they cannot drift apart.

    ``limma_test.py`` forwards its default as ``--mode``, so a divergence here
    would be invisible in the pipeline but would change what a bare
    ``Rscript limma_test.R`` produces.
    """
    import limma_test

    r_source = R_SCRIPT.read_text(encoding="utf-8")
    assert f'DEFAULT_EBAYES_MODE <- "{limma_test.DEFAULT_EBAYES_MODE}"' in r_source, (
        "limma_test.R's DEFAULT_EBAYES_MODE does not match limma_test.py's "
        f"{limma_test.DEFAULT_EBAYES_MODE!r}"
    )
    assert limma_test.DEFAULT_EBAYES_MODE == DEFAULT_MODE
