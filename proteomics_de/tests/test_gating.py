"""The forward-path gate, actually exercised.

``gated/pca_cluster.py`` has always claimed that its gated analyses "switch
themselves on" as the sample sheet grows. Until this file existed the claim was
untested in both directions: the dispatcher was only ever evaluated at the
committed n=4, so *no* gate had ever been observed to flip, and the sample count
was a hardcoded ``N_SAMPLES = 4`` that no sheet could move anyway.

These tests pin the behaviour at three sample counts:

    n=4   the committed sheet -- PCA and clustering QC_ONLY, everything else SKIP
    n=6   research1.md SS5's trustworthy gate -- PCA and clustering flip, nothing else
    n=20  the WGCNA FAQ's "at least 20 samples" -- WGCNA flips too

n=6 and n=20 are synthetic sheets written into ``tmp_path`` and read back
through ``config.design``, so what is under test is the real path from TSV to
gate decision, not a hand-passed integer.

Nothing here is allowed to be a loose smoke test. If one of these fails, the
forward path is broken and the correct response is to fix the dispatcher -- not
to relax the expectation.
"""

from __future__ import annotations

import importlib
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

# Reproduce the sys.path layout the pipeline scripts expect. Module-level, because
# conftest's autouse fixture has not run by the time these imports execute.
os.environ.setdefault("MPLBACKEND", "Agg")

_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPO_ROOT = _PKG_DIR.parent
for _entry in (_PKG_DIR / "viz", _PKG_DIR, _REPO_ROOT):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))

from gated import pca_cluster  # noqa: E402
from proteomics_de.config import design  # noqa: E402

GATED_RESULTS_DIR = _PKG_DIR / "results" / "gated"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def write_sheet(tmp_path: Path, replicates_per_group: int, groups=("control", "treated")) -> Path:
    """Write a synthetic balanced sample sheet and return its path.

    Same four columns as the committed sheet, so ``config.design`` reads it with
    no special-casing. Sample ids are arbitrary but unique.
    """
    lines = ["\t".join(design.REQUIRED_COLUMNS)]
    sample_id = 40000
    for group in groups:
        for replicate in range(1, replicates_per_group + 1):
            sample_id += 1
            lines.append(f"{sample_id}\t{group}\tIntensity {sample_id}\t{replicate}")
    path = tmp_path / "sample_sheet.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def statuses_for_sheet(sheet_path: Path) -> dict[str, str]:
    """Gate statuses obtained by reading `sheet_path` -- the full TSV -> gate path."""
    n_samples = design.n_samples(sheet_path)
    n_replicates = design.replicates_per_group(sheet_path)
    registry = pca_cluster.build_registry(n_groups=design.n_groups(sheet_path))
    rows = pca_cluster.evaluate_registry(registry, n_samples, n_replicates)
    return {row["analysis"]: row["status"] for row in rows}


def rows_for(n_samples: int, n_replicates_per_group: int) -> list[dict]:
    return pca_cluster.evaluate_registry(
        pca_cluster.REGISTRY, n_samples, n_replicates_per_group
    )


ALL_ANALYSES = {
    "pca",
    "hierarchical_clustering",
    "wgcna",
    "umap_tsne",
    "mixomics_splsda",
    "vae_gnn",
}


# ---------------------------------------------------------------------------
# The counts come from the sample sheet, not from a literal
# ---------------------------------------------------------------------------
def test_module_counts_are_read_from_the_sample_sheet():
    """The whole point of the refactor: no hardcoded 4 anywhere in the module."""
    assert pca_cluster.N_SAMPLES == design.n_samples()
    assert pca_cluster.N_GROUPS == design.n_groups()
    assert pca_cluster.N_REPLICATES_PER_GROUP == design.replicates_per_group()


def test_module_counts_follow_a_different_sample_sheet(tmp_path, monkeypatch):
    """The decisive test for the refactor.

    ``N_SAMPLES == design.n_samples()`` alone proves nothing while both are 4 --
    a hardcoded ``N_SAMPLES = 4`` would satisfy it. Point the design layer at a
    10-sample sheet, re-import the module, and the gates must move with it.
    """
    sheet = write_sheet(tmp_path, replicates_per_group=5)  # 10 samples, 2 groups
    monkeypatch.setattr(design, "DEFAULT_SAMPLE_SHEET", sheet)
    try:
        reloaded = importlib.reload(pca_cluster)
        assert reloaded.N_SAMPLES == 10
        assert reloaded.N_GROUPS == 2
        assert reloaded.N_REPLICATES_PER_GROUP == 5

        statuses = {row["analysis"]: row["status"] for row in reloaded.dispatch()}
        assert statuses["pca"] == "RUN"          # was QC_ONLY at n=4
        assert statuses["hierarchical_clustering"] == "RUN"
        assert statuses["wgcna"] == "SKIP"       # 10 < 15, still off
    finally:
        monkeypatch.undo()
        importlib.reload(pca_cluster)

    # ...and the module is back on the committed sheet for every other test.
    assert pca_cluster.N_SAMPLES == design.n_samples() == 4


def test_dispatch_is_a_thin_wrapper_over_evaluate_registry():
    """dispatch() must add nothing but the design read."""
    assert pca_cluster.dispatch() == rows_for(
        design.n_samples(), design.replicates_per_group()
    )


def test_registry_covers_exactly_the_declared_analyses():
    assert {e["analysis"] for e in pca_cluster.REGISTRY} == ALL_ANALYSES


def test_statuses_come_from_the_declared_vocabulary():
    for n, r in ((4, 2), (6, 3), (20, 10), (1000, 500)):
        for row in rows_for(n, r):
            assert row["status"] in pca_cluster.STATUSES


# ---------------------------------------------------------------------------
# n=4 -- the committed sheet
# ---------------------------------------------------------------------------
EXPECTED_AT_N4 = {
    "pca": "QC_ONLY",
    "hierarchical_clustering": "QC_ONLY",
    "wgcna": "SKIP",
    "umap_tsne": "SKIP",
    "mixomics_splsda": "SKIP",
    "vae_gnn": "SKIP",
}


def test_gate_table_at_the_committed_design():
    """PCA/clustering run as QC only; every multivariate/ML analysis is off."""
    assert design.n_samples() == 4 and design.replicates_per_group() == 2, (
        "this test pins the n=4 committed design; if the sheet grew, the "
        "expectations below must be re-derived, not deleted"
    )
    statuses = {row["analysis"]: row["status"] for row in pca_cluster.dispatch()}
    assert statuses == EXPECTED_AT_N4


def test_every_row_at_n4_carries_a_written_reason():
    """A gate decision with no stated reason is not auditable."""
    for row in pca_cluster.dispatch():
        reason = row["reason"]
        assert isinstance(reason, str) and reason.strip(), row["analysis"]
        # The reason must name the gate it is reporting against...
        assert str(row["min_samples"]) in reason, row["analysis"]
        # ...and the count we actually have.
        assert str(design.n_samples()) in reason, row["analysis"]


def test_skip_reasons_state_the_shortfall():
    """Each SKIP explains what is missing, in the units of the sheet."""
    n = design.n_samples()
    for row in pca_cluster.dispatch():
        if row["status"] != "SKIP":
            continue
        assert f"have {n}." in row["reason"], row["analysis"]


def test_qc_only_reasons_name_the_trustworthy_gate():
    gate = pca_cluster.TRUSTWORTHY_MIN_SAMPLES
    for row in pca_cluster.dispatch():
        if row["status"] != "QC_ONLY":
            continue
        assert f"trustworthy gate n>={gate}" in row["reason"], row["analysis"]
        assert "research1.md" in row["reason"], row["analysis"]


# ---------------------------------------------------------------------------
# n=6 -- research1.md SS5's trustworthy gate
# ---------------------------------------------------------------------------
def test_n6_flips_exactly_pca_and_clustering(tmp_path):
    """At the n>=6 gate, PCA and clustering become trustworthy -- and nothing else.

    "Nothing else" is the load-bearing half: a dispatcher that turned everything
    on at once would pass a test that only checked pca/clustering.
    """
    sheet = write_sheet(tmp_path, replicates_per_group=3)
    assert design.n_samples(sheet) == 6
    assert design.replicates_per_group(sheet) == 3

    statuses = statuses_for_sheet(sheet)
    assert statuses == {
        "pca": "RUN",
        "hierarchical_clustering": "RUN",
        "wgcna": "SKIP",
        "umap_tsne": "SKIP",
        "mixomics_splsda": "SKIP",
        "vae_gnn": "SKIP",
    }

    flipped = {a for a, s in statuses.items() if s != EXPECTED_AT_N4[a]}
    assert flipped == {"pca", "hierarchical_clustering"}


# ---------------------------------------------------------------------------
# n=20 -- the WGCNA FAQ's "at least 20 samples"
# ---------------------------------------------------------------------------
def test_n20_switches_wgcna_on(tmp_path):
    """WGCNA's gate is 15 to run / 20 to trust, per the WGCNA FAQ quoted in
    research1.md:235. UMAP still needs dozens; VAE/GNN needs thousands."""
    sheet = write_sheet(tmp_path, replicates_per_group=10)
    assert design.n_samples(sheet) == 20
    assert design.replicates_per_group(sheet) == 10

    statuses = statuses_for_sheet(sheet)
    assert statuses["pca"] == "RUN"
    assert statuses["hierarchical_clustering"] == "RUN"
    assert statuses["wgcna"] == "RUN"
    assert statuses["vae_gnn"] == "SKIP"
    # Full table, so an over-eager gate cannot hide in an unasserted row.
    assert statuses == {
        "pca": "RUN",
        "hierarchical_clustering": "RUN",
        "wgcna": "RUN",
        "umap_tsne": "SKIP",
        "mixomics_splsda": "RUN",
        "vae_gnn": "SKIP",
    }


def test_wgcna_is_qc_only_between_its_run_and_trustworthy_gates(tmp_path):
    """15 samples runs it, 20 makes it trustworthy -- the band between is QC_ONLY."""
    sheet = write_sheet(tmp_path, replicates_per_group=8)  # 16 samples
    assert design.n_samples(sheet) == 16
    assert statuses_for_sheet(sheet)["wgcna"] == "QC_ONLY"


def test_vae_gnn_stays_off_at_every_realistic_sample_count():
    """Thousands of samples is not a proteomics-core sample sheet."""
    for n, r in ((4, 2), (6, 3), (20, 10), (100, 50), (500, 250)):
        statuses = {row["analysis"]: row["status"] for row in rows_for(n, r)}
        assert statuses["vae_gnn"] == "SKIP", n


def test_reasons_are_rendered_at_the_evaluated_n():
    """A reason may never go stale by reporting a sample count it was not
    evaluated at -- that was the bug the hardcoded 'have 4.' strings invited."""
    for row in rows_for(20, 10):
        assert "have 4." not in row["reason"], row["analysis"]
        if row["status"] == "SKIP":
            assert "have 20." in row["reason"], row["analysis"]


def test_build_registry_scales_the_splsda_gate_with_group_count():
    """sPLS-DA's gate is samples >> classes, so it depends on the design's
    group count -- the one registry number that is not a fixed threshold."""
    two = {e["analysis"]: e for e in pca_cluster.build_registry(n_groups=2)}
    four = {e["analysis"]: e for e in pca_cluster.build_registry(n_groups=4)}
    assert two["mixomics_splsda"]["min_samples_to_run"] == 20
    assert four["mixomics_splsda"]["min_samples_to_run"] == 40
    # The pure research1.md thresholds must NOT move with group count.
    assert two["wgcna"]["min_samples_to_run"] == four["wgcna"]["min_samples_to_run"] == 15


# ---------------------------------------------------------------------------
# The committed artifacts agree with the dispatcher
# ---------------------------------------------------------------------------
def test_skip_log_has_one_row_per_registered_analysis():
    """Row count is derived from the registry, never a literal 6: adding a
    candidate analysis must show up in the log without editing this test."""
    df = pd.read_csv(GATED_RESULTS_DIR / "skip_log.csv")
    assert len(df) == len(pca_cluster.REGISTRY)
    assert list(df["analysis"]) == [e["analysis"] for e in pca_cluster.REGISTRY]


def test_committed_skip_log_matches_the_current_dispatcher():
    """The committed artifact and the code that writes it must not drift."""
    df = pd.read_csv(GATED_RESULTS_DIR / "skip_log.csv")
    expected = pd.DataFrame(
        pca_cluster.dispatch(),
        columns=["analysis", "min_samples", "have_samples", "status", "reason"],
    )
    pd.testing.assert_frame_equal(df, expected)


def test_skip_log_reports_the_sheets_sample_count():
    df = pd.read_csv(GATED_RESULTS_DIR / "skip_log.csv")
    assert set(df["have_samples"]) == {design.n_samples()}


# ---------------------------------------------------------------------------
# PCA's rank degeneracy, as research1.md:246 describes it
# ---------------------------------------------------------------------------
def test_pca_variance_has_exactly_n_minus_one_components():
    """rank <= n-1 (research1.md:246): centering costs a degree of freedom, so
    at n=4 exactly 3 PCs exist. Derived from the sheet, not written as '3'."""
    variance = pd.read_csv(GATED_RESULTS_DIR / "pca_variance.csv")
    assert len(variance) == design.n_samples() - 1
    assert list(variance["PC"]) == [
        f"PC{i + 1}" for i in range(design.n_samples() - 1)
    ]


def test_pca_variance_sums_to_one_which_is_the_degeneracy():
    """All variance lands in n-1 PCs by construction -- that IS the rank
    deficiency, and it is why the output is QC-only."""
    variance = pd.read_csv(GATED_RESULTS_DIR / "pca_variance.csv")
    assert variance["variance_explained"].sum() == pytest.approx(1.0, abs=1e-9)


def test_pca_coords_has_one_row_per_sample():
    coords = pd.read_csv(GATED_RESULTS_DIR / "pca_coords.csv")
    assert len(coords) == design.n_samples()
    assert list(coords.columns) == ["sample"] + [
        f"PC{i + 1}" for i in range(design.n_samples() - 1)
    ]


def test_gated_layer_reads_the_frozen_qc_limma_table(frozen_counts):
    """PCA/clustering must be built from the same table the rest of the
    pipeline froze; expected row count comes from frozen_counts.json."""
    _df, log2_complete, n_total, n_complete = pca_cluster.load_complete_matrix()
    assert n_total == frozen_counts["qc_limma_rows"]
    assert 0 < n_complete <= n_total
    assert list(log2_complete.columns) == pca_cluster.SAMPLE_ORDER
    assert len(log2_complete.columns) == design.n_samples()
