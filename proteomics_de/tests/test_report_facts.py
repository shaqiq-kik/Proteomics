"""Tests for report/build_facts.py and the committed report/report_facts.json.

``report_facts.json`` used to be hand-authored, which meant a full pipeline
re-run refreshed the report's figures but left its numbers behind. These tests
pin the file to its sources: every block must equal what the committed
artifacts say, and the file on disk must be exactly what the generator emits.

The strict-JSON tests exist because the hand-authored file shipped a bare
``NaN`` literal. ``json.load`` accepts that (it is a Python extension, not
JSON); JavaScript's ``JSON.parse``, Go's ``encoding/json`` and ``serde_json``
all reject it outright, and jq silently coerces it to ``null``. Either way the
digest was not portable, and the value was wrong besides.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPORT_DIR = _PKG_DIR / "report"
_FACTS_PATH = _REPORT_DIR / "report_facts.json"


def _load_build_facts():
    """Import report/build_facts.py by path (report/ is not a package)."""
    spec = importlib.util.spec_from_file_location(
        "build_facts", _REPORT_DIR / "build_facts.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def build_facts():
    return _load_build_facts()


@pytest.fixture(scope="module")
def generated(build_facts) -> dict:
    """The digest as build_facts.py would write it right now."""
    return build_facts.build()


@pytest.fixture(scope="module")
def committed() -> dict:
    return json.loads(_FACTS_PATH.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def figures_dir() -> Path:
    return _PKG_DIR / "results" / "figures"


@pytest.fixture(scope="module")
def enrich_dir() -> Path:
    return _PKG_DIR / "results" / "enrichment"


# ---------------------------------------------------------------------------
# strict JSON
# ---------------------------------------------------------------------------
def _reject_constant(token: str):
    raise AssertionError(
        f"non-strict JSON constant {token!r} in report_facts.json. NaN/Infinity "
        f"are Python extensions; JS/Go/serde reject them and jq silently turns "
        f"them into null."
    )


def test_committed_facts_parse_under_strict_json():
    """The committed file must contain no NaN/Infinity literals."""
    text = _FACTS_PATH.read_text(encoding="utf-8")
    json.loads(text, parse_constant=_reject_constant)


def test_generated_facts_round_trip_under_strict_json(build_facts, generated):
    """What the generator emits must survive a strict round-trip unchanged."""
    text = build_facts.serialise(generated)
    reparsed = json.loads(text, parse_constant=_reject_constant)
    assert reparsed == generated
    # and serialising the reparsed copy is byte-stable
    assert build_facts.serialise(reparsed) == text


def test_generator_refuses_to_write_non_finite_values(build_facts):
    """allow_nan=False is what makes the strict guarantee real, not a habit."""
    with pytest.raises(ValueError):
        build_facts.serialise({"bad": float("nan")})
    with pytest.raises(ValueError):
        build_facts.serialise({"bad": float("inf")})


def test_no_bare_nan_token_anywhere_in_the_file():
    """Belt and braces: the literal that broke every non-Python parser."""
    text = _FACTS_PATH.read_text(encoding="utf-8")
    for token in (": NaN", ": Infinity", ": -Infinity"):
        assert token not in text, f"{token!r} present in report_facts.json"


def test_committed_file_is_exactly_what_the_generator_produces(build_facts, generated):
    """No hand-edits. If this fails, re-run build_facts.py rather than patching."""
    assert _FACTS_PATH.read_text(encoding="utf-8") == build_facts.serialise(generated)


# ---------------------------------------------------------------------------
# figures: the verbatim union of the four manifests
# ---------------------------------------------------------------------------
def test_figures_are_the_union_of_the_four_manifests(committed, build_facts, figures_dir):
    union = []
    for name in build_facts.MANIFESTS:
        union.extend(json.loads((figures_dir / name).read_text(encoding="utf-8")))

    assert len(union) == 13, "expected 7 + 3 + 2 + 1 manifest entries"
    assert committed["n_figures"] == 13
    assert len(committed["figures"]) == 13
    assert [f["file"] for f in committed["figures"]] == [f["file"] for f in union]


@pytest.mark.parametrize("idx", range(13))
def test_each_figure_entry_matches_its_manifest_source_byte_for_byte(
    idx, committed, build_facts, figures_dir
):
    union = []
    for name in build_facts.MANIFESTS:
        union.extend(json.loads((figures_dir / name).read_text(encoding="utf-8")))

    got, want = committed["figures"][idx], union[idx]
    assert got["file"] == want["file"]
    assert got["title"] == want["title"], f"title drift for {want['file']}"
    assert got["caption"] == want["caption"], f"caption drift for {want['file']}"
    assert got.get("key_numbers") == want.get("key_numbers"), (
        f"key_numbers drift for {want['file']}"
    )


def test_every_figure_file_referenced_actually_exists(committed, figures_dir):
    for entry in committed["figures"]:
        stem = Path(entry["file"]).stem
        assert (figures_dir / f"{stem}.png").exists() or (
            figures_dir / f"{stem}.svg"
        ).exists(), f"no figure file on disk for {entry['file']}"


# ---------------------------------------------------------------------------
# numeric blocks vs their source artifacts
# ---------------------------------------------------------------------------
def test_counts_match_the_committed_csvs(committed, results_dir):
    fa = pd.read_csv(results_dir / "foldchange_all.csv")
    sc = pd.read_csv(results_dir / "single_condition_proteins.csv")
    oo = pd.read_csv(results_dir / "onoff_proteins.csv")
    vc = fa["regulated"].value_counts()
    c = committed["counts"]

    assert c["both_condition"] == len(fa)
    assert c["single_condition"] == len(sc)
    assert c["on_off"] == len(oo)
    assert c["detected_universe"] == len(fa) + len(sc)
    assert c["regulated_UP"] == int(vc.get("UP", 0))
    assert c["regulated_DOWN"] == int(vc.get("DOWN", 0))
    assert c["no_change"] == int(vc.get("NO CHANGE", 0))
    assert c["regulated_total"] == c["regulated_UP"] + c["regulated_DOWN"]


def test_counts_match_frozen_counts(committed, frozen_counts):
    """The same numbers the rest of the suite asserts against (D7 + D11)."""
    c = committed["counts"]
    assert c["both_condition"] == frozen_counts["foldchange_all_rows"]
    assert c["single_condition"] == frozen_counts["single_condition_rows"]
    assert c["on_off"] == frozen_counts["onoff_rows"]
    assert c["detected_universe"] == frozen_counts["background_union"]
    assert c["regulated_UP"] == frozen_counts["n_up"]
    assert c["regulated_DOWN"] == frozen_counts["n_down"]
    assert c["regulated_total"] == frozen_counts["ipa_input_rows"]


def test_d7_orientation_is_the_corrected_one(committed):
    """509 UP / 206 DOWN, not the 206/509 the hand-authored file carried."""
    assert committed["counts"]["regulated_UP"] == 509
    assert committed["counts"]["regulated_DOWN"] == 206
    assert committed["string"]["n_seeds_up"] == 509
    assert committed["string"]["n_seeds_down"] == 206
    assert committed["ora"]["up_n"] == 509
    assert committed["ora"]["down_n"] == 206


def test_d11_quarantine_is_reflected(committed):
    assert committed["counts"]["single_condition"] == 604
    assert committed["counts"]["detected_universe"] == 2552
    assert committed["ora"]["background_size_row_union"] == 2552


def test_qc_centering_matches_its_csv(committed, results_dir):
    row = pd.read_csv(results_dir / "qc_centering.csv").iloc[0]
    got = committed["qc"]["centering"]
    for key in row.index:
        assert got[key] == pytest.approx(row[key]) if isinstance(
            row[key], (int, float, np.floating, np.integer)
        ) else got[key] == row[key], f"qc.centering.{key} drifted"


def test_qc_replicate_matches_its_csv_one_to_one(committed, results_dir):
    row = pd.read_csv(results_dir / "qc_replicate_correlation.csv").iloc[0]
    got = committed["qc"]["replicate"]
    assert set(got) == set(row.index)
    for key in row.index:
        if isinstance(row[key], (float, np.floating)):
            assert got[key] == pytest.approx(row[key])
        else:
            assert got[key] == row[key]


def test_replicate_raw_r_is_recomputed_against_the_corrected_sheet(
    committed, results_dir
):
    """The recomputed values must be keyed by the sample sheet, not by
    replicate_check.py's hardcoded pre-D7 channel pairs."""
    from proteomics_de.config import design as design_cfg

    fa = pd.read_csv(results_dir / "foldchange_all.csv")
    by_group = committed["qc"]["replicate_raw_r_by_group"]

    for group in ("control", "treated"):
        cols = design_cfg.columns_for_group(group)
        assert by_group[group]["channels"] == cols
        x = pd.to_numeric(fa[cols[0]], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(fa[cols[1]], errors="coerce").to_numpy(dtype=float)
        m = (x > 0) & (y > 0)
        r = float(np.corrcoef(np.log10(x[m]), np.log10(y[m]))[0, 1])
        assert by_group[group]["raw_log10_pearson_r"] == pytest.approx(r, abs=5e-5)
        assert by_group[group]["n"] == int(m.sum())


def test_limma_block_matches_recomputation(committed, results_dir):
    trend = pd.read_csv(results_dir / "qc_limma.csv")
    vanilla = pd.read_csv(results_dir / "qc_limma_vanilla.csv")
    lim = committed["limma"]

    assert lim["n_tested"] == len(trend)
    assert lim["n_sig_fdr005"] == int((trend["adj_p_value"] < 0.05).sum())
    assert lim["min_raw_p"] == pytest.approx(float(trend["p_value"].min()))
    assert lim["n_raw_p_lt05"] == int((trend["p_value"] < 0.05).sum())
    assert lim["min_adj_p"] == pytest.approx(float(trend["adj_p_value"].min()))

    assert lim["vanilla_n_sig_fdr005"] == int((vanilla["adj_p_value"] < 0.05).sum())
    assert lim["vanilla_min_raw_p"] == pytest.approx(float(vanilla["p_value"].min()))
    assert lim["vanilla_n_raw_p_lt05"] == int((vanilla["p_value"] < 0.05).sum())
    assert lim["vanilla_min_adj_p"] == pytest.approx(
        float(vanilla["adj_p_value"].min())
    )


def test_d9_default_is_trend_robust_and_headline_null_holds(committed):
    """Trend/robust is the default; vanilla is now the comparison run."""
    lim = committed["limma"]
    assert "trend" in lim["ebayes_mode"] and "robust" in lim["ebayes_mode"]
    assert lim["min_adj_p"] == pytest.approx(0.116, abs=5e-4)
    assert lim["vanilla_min_adj_p"] == pytest.approx(0.305, abs=5e-4)
    # The headline is unchanged by D7/D9/D11: still nothing at FDR < 0.05.
    assert lim["n_tested"] == 1938
    assert lim["n_sig_fdr005"] == 0
    assert lim["vanilla_n_sig_fdr005"] == 0
    # ...and D6 still holds: no enrichment term survives the honest background.
    assert committed["ora"]["n_terms_up"] == 0
    assert committed["ora"]["n_terms_down"] == 0
    assert committed["gsea"]["n_fdr_lt_0.05"] == 0


def test_corr_limma_vs_pipeline_is_one(committed):
    """The value that used to be a bare NaN. It is +1.0, not missing."""
    corr = committed["limma"]["corr_limma_vs_pipeline_log2FC"]
    assert isinstance(corr, float)
    assert math.isfinite(corr)
    assert corr > 0.9999


def test_corr_is_computed_on_the_paired_non_null_subset(committed, results_dir):
    """Recompute independently, and prove the naive version really is NaN."""
    trend = pd.read_csv(results_dir / "qc_limma.csv")
    fa = pd.read_csv(results_dir / "foldchange_all.csv")
    merged = trend.merge(
        fa[["UniProt Accession Number", "log2FC"]],
        left_on="id",
        right_on="UniProt Accession Number",
        how="inner",
    )
    a = pd.to_numeric(merged["limma_log2FC"], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(merged["log2FC"], errors="coerce").to_numpy(dtype=float)
    paired = np.isfinite(a) & np.isfinite(b)

    assert committed["limma"]["corr_n_paired"] == int(paired.sum())
    assert committed["limma"]["corr_n_unpaired_partial_rows"] == int((~paired).sum())
    assert committed["limma"]["corr_limma_vs_pipeline_log2FC"] == pytest.approx(
        float(np.corrcoef(a[paired], b[paired])[0, 1])
    )
    # The bug being guarded against: unmasked, numpy propagates the NaNs.
    assert not np.isfinite(np.corrcoef(a, b)[0, 1])


def test_top_candidates_are_the_ten_smallest_p_from_the_default_run(
    committed, results_dir
):
    top = pd.read_csv(results_dir / "qc_limma.csv").nsmallest(10, "p_value")
    got = committed["top_candidates"]

    assert len(got) == 10
    assert [g["acc"] for g in got] == top["id"].tolist()
    assert [g["regulated"] for g in got] == top["regulated"].tolist()
    for g, (_, row) in zip(got, top.iterrows()):
        assert g["p"] == pytest.approx(float(row["p_value"]))
        assert g["log2FC"] == pytest.approx(round(float(row["limma_log2FC"]), 3))
        assert g["n_imputed"] == int(row["n_imputed"])  # D10


@pytest.mark.parametrize(
    "block,filename",
    [
        ("string", "string_meta.json"),
        ("ora", "ora_meta.json"),
        ("gsea", "gsea_meta.json"),
    ],
)
def test_enrichment_blocks_are_verbatim_copies_of_their_meta(
    block, filename, committed, enrich_dir
):
    assert committed[block] == json.loads(
        (enrich_dir / filename).read_text(encoding="utf-8")
    )


def test_pca_variance_matches_its_csv(committed, results_dir):
    df = pd.read_csv(results_dir / "gated" / "pca_variance.csv")
    got = committed["pca_variance"]
    assert len(got) == len(df)
    for rec, (_, row) in zip(got, df.iterrows()):
        assert rec["PC"] == row["PC"]
        assert rec["variance_explained"] == pytest.approx(row["variance_explained"])


def test_gate_skips_are_the_skip_log_as_records(committed, results_dir):
    df = pd.read_csv(results_dir / "gated" / "skip_log.csv")
    assert committed["gate_skips"] == df.to_dict(orient="records")


def test_dataset_block_is_read_from_the_sample_sheet(committed):
    from proteomics_de.config import design as design_cfg

    sheet = design_cfg.read_sample_sheet()
    ds = committed["dataset"]
    assert ds["control_samples"] == sheet.loc[
        sheet["group"] == "control", "sample"
    ].tolist()
    assert ds["treated_samples"] == sheet.loc[
        sheet["group"] == "treated", "sample"
    ].tolist()
    assert ds["n_samples"] == design_cfg.n_samples()
    assert ds["replicates_per_group"] == design_cfg.replicates_per_group()
    # D7: the corrected assignment, not the inherited one.
    assert ds["control_samples"] == ["31579", "31581"]
    assert ds["treated_samples"] == ["31578", "31580"]
