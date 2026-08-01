# STATUS — what actually exists in `proteomics_de/`

> **This file is GENERATED. Do not hand-edit.**
> Regenerate with `.venv/bin/python tools/status.py` (source: `tools/status.py`). Every number below is read off the filesystem at generation time, so it cannot silently go stale.

* Generated: **2026-08-01 02:08:01 EDT**
* Git HEAD: **`924199b`**
* Repo root: `/Users/abrarshakik/Documents/Proteomics`

Companions: `DECISIONS_LOG.md` (human decisions D1–D6) · `BUILD_LOG.md` (per-work-package history) · `../research1.md` §Build Log (per-bug narrative).

---

## a) research1.md Section 6 build list (items 1-20)

**8 implemented as proposed · 12 implemented elsewhere · 0 missing** (of 20).

`↔ implemented elsewhere` is not a gap: the build kept the flat, individually-verified script layout instead of the doc's `io/ etl/ de/` tree (DECISIONS_LOG D1), and substituted Python for R in three places.

| # | Doc proposes | Actual file(s) on disk | Status | Note |
|---|---|---|---|---|
| 1 | `config/sample_sheet.tsv + config.yaml` | `config/config.yaml`, `config/sample_sheet.tsv` | ✅ implemented | Descriptive only -- the frozen scripts do not read it (DECISIONS_LOG D1). |
| 2 | `io/load.py` | `foldchange.py` | ↔ implemented elsewhere | Inlined: foldchange.py reads 'Protein Report L'/'H' directly via pd.read_excel. |
| 3 | `qc/schema.py (pandera)` | `qc/schema.py`<br>+ `qc/validate.py`, `qc/__init__.py` | ✅ implemented | Built as specified; validate.py is the runner that emits results/qc/qc_report.*. |
| 4 | `etl/merge.py` | `foldchange.py` | ↔ implemented elsewhere | Inlined: outer merge with indicator; exclusives -> single_condition_proteins.csv. |
| 5 | `etl/foldchange.py` | `foldchange.py` | ↔ implemented elsewhere | Same module, flat layout instead of etl/ (DECISIONS_LOG D1). |
| 6 | `qc/replicate_qc.py` | `replicate_check.py`, `centering_check.py` | ↔ implemented elsewhere | Split into two flat modules: replicate correlation (bug 6) + centering check. |
| 7 | `etl/build_matrix.py` | `limma_test.py` | ↔ implemented elsewhere | Inlined: limma_test.py writes the handoff CSV (_limma_input.csv) directly. |
| 8 | `de/run_limma.R` | `limma_test.R` | ↔ implemented elsewhere | Same role, flat layout; emits _limma_output.csv + _limma_versions.txt. |
| 9 | `de/run_limma.py` | `limma_test.py` | ↔ implemented elsewhere | Same role, flat layout. |
| 10 | `export/ipa_export.py` | `foldchange.py`, `limma_test.py` | ↔ implemented elsewhere | Inlined: foldchange.py writes ipa_input.csv; limma_test.py writes ipa_input_significant.csv. |
| 11 | `viz/volcano.py / .R` | `viz/volcano.py` | ✅ implemented | Python variant built (the doc allowed either). |
| 12 | `viz/ma_plot.py` | `viz/ma_plot.py` | ✅ implemented |  |
| 13 | `viz/heatmap.R (pheatmap)` | `viz/heatmap.py` | ↔ implemented elsewhere | Built in Python (matplotlib/scipy) instead of R/pheatmap -- no extra R dependency. |
| 14 | `viz/qc_plots.py` | `viz/qc_plots.py` | ✅ implemented |  |
| 15 | `enrich/string_ppi.py` | `enrich/string_ppi.py`<br>+ `enrich/enrich_common.py` | ✅ implemented | networkx used for graph metrics instead of igraph; STRING species=10090 (D5). |
| 16 | `enrich/ora_gsea.R` | `enrich/ora.py`, `enrich/gsea.py` | ↔ implemented elsewhere | Python + web APIs (g:Profiler, gseapy/Enrichr) instead of R Bioconductor (D3). |
| 17 | `viz/network_cytoscape.py` | `enrich/network_figure.py` | ↔ implemented elsewhere | Static networkx figure -- py4cytoscape needs a running Cytoscape desktop (D3). |
| 18 | `viz/upset.py` | `enrich/upset.py` | ↔ implemented elsewhere | Lives under enrich/ rather than viz/ (it consumes enrichment sets). |
| 19 | `gated/pca_cluster.py` | `gated/pca_cluster.py` | ✅ implemented | Gate fires: n=4 < 6, so output is QC-only + results/gated/skip_log.csv. |
| 20 | `report/build_report.py` | `report/build_report.py`<br>+ `report/report_template.html`, `report/report.html`, `report/report_facts.json` | ✅ implemented | Self-contained single-file HTML report. |

---

## b) Artifact inventory (`proteomics_de/results/`)

**76 files**, 20.7 MB total, 24 tabular (`.csv`/`.tsv`).

Row counts EXCLUDE the header. Three files are header-only **by design** — they are the honest scientific result, not a failure. See `DECISIONS_LOG.md` D2 and D6.

| File | Size | Data rows | Note |
|---|---|---|---|
| `results/de/design.tsv` | 69 B | 4 |  |
| `results/de/intensity_matrix.tsv` | 110.8 KB | 1,938 |  |
| `results/de/limma_results.tsv` | 272.5 KB | 1,938 |  |
| `results/enrichment/gsea_meta.json` | 1.5 KB | — |  |
| `results/enrichment/gsea_results.csv` | 114.1 KB | 568 | matches expected 568 |
| `results/enrichment/ora_down.csv` | 91 B | **0 rows (expected -- 0 GO/KEGG/Reactome terms survive the honest detected-proteome background (DECISIONS_LOG D6))** |  |
| `results/enrichment/ora_meta.json` | 5.4 KB | — |  |
| `results/enrichment/ora_top_terms_detail.json` | 4.1 KB | — |  |
| `results/enrichment/ora_up.csv` | 91 B | **0 rows (expected -- 0 GO/KEGG/Reactome terms survive the honest detected-proteome background (DECISIONS_LOG D6))** |  |
| `results/enrichment/raw/enrichr_libraries.json.gz` | 686.7 KB | — |  |
| `results/enrichment/raw/gprofiler_down.json` | 27.4 KB | — |  |
| `results/enrichment/raw/gprofiler_down_all.json` | 3.7 MB | — |  |
| `results/enrichment/raw/gprofiler_down_all_evidence.json.gz` | 650.2 KB | — |  |
| `results/enrichment/raw/gprofiler_up.json` | 63.2 KB | — |  |
| `results/enrichment/raw/gprofiler_up_all.json` | 5.5 MB | — |  |
| `results/enrichment/raw/gprofiler_up_all_evidence.json.gz` | 1.1 MB | — |  |
| `results/enrichment/raw/string_get_ids.json` | 386.7 KB | — |  |
| `results/enrichment/raw/string_network.tsv` | 574.3 KB | 5,962 |  |
| `results/enrichment/raw/string_version.json` | 100 B | — |  |
| `results/enrichment/string_edges.tsv` | 103.9 KB | 5,963 |  |
| `results/enrichment/string_meta.json` | 760 B | — |  |
| `results/enrichment/string_node_metrics.csv` | 57.5 KB | 694 | matches expected 694 |
| `results/figures/figures_manifest.json` | 8.2 KB | — |  |
| `results/figures/figures_manifest_enrich.json` | 13.2 KB | — |  |
| `results/figures/figures_manifest_gated.json` | 2.9 KB | — |  |
| `results/figures/figures_manifest_network.json` | 3.5 KB | — |  |
| `results/figures/gsea_top.png` | 376.9 KB | — |  |
| `results/figures/gsea_top.svg` | 27.3 KB | — |  |
| `results/figures/heatmap_top_de.png` | 163.9 KB | — |  |
| `results/figures/heatmap_top_de.svg` | 89.4 KB | — |  |
| `results/figures/intensity_distributions.png` | 99.1 KB | — |  |
| `results/figures/intensity_distributions.svg` | 34.3 KB | — |  |
| `results/figures/ma_plot.png` | 248.6 KB | — |  |
| `results/figures/ma_plot.svg` | 226.7 KB | — |  |
| `results/figures/missing_values.png` | 204.1 KB | — |  |
| `results/figures/missing_values.svg` | 32.3 KB | — |  |
| `results/figures/ora_dotplot.png` | 214.4 KB | — |  |
| `results/figures/ora_dotplot.svg` | 25.6 KB | — |  |
| `results/figures/pca_qc.png` | 134.4 KB | — |  |
| `results/figures/pca_qc.svg` | 18.3 KB | — |  |
| `results/figures/ppi_network.png` | 1.2 MB | — |  |
| `results/figures/ppi_network.svg` | 1.4 MB | — |  |
| `results/figures/rank_abundance.png` | 101.0 KB | — |  |
| `results/figures/rank_abundance.svg` | 155.7 KB | — |  |
| `results/figures/sample_correlation.png` | 121.0 KB | — |  |
| `results/figures/sample_correlation.svg` | 20.3 KB | — |  |
| `results/figures/sample_dendrogram.png` | 112.3 KB | — |  |
| `results/figures/sample_dendrogram.svg` | 12.4 KB | — |  |
| `results/figures/upset.png` | 271.5 KB | — |  |
| `results/figures/upset.svg` | 43.6 KB | — |  |
| `results/figures/volcano.png` | 208.4 KB | — |  |
| `results/figures/volcano.svg` | 230.4 KB | — |  |
| `results/foldchange_all.csv` | 281.6 KB | 1,948 | matches expected 1,948 |
| `results/foldchange_all.csv.provenance.json` | 806 B | — |  |
| `results/foldchange_all_centered.csv` | 327.7 KB | 1,948 |  |
| `results/gated/pca_coords.csv` | 287 B | 4 |  |
| `results/gated/pca_variance.csv` | 93 B | 3 |  |
| `results/gated/skip_log.csv` | 1.2 KB | 6 |  |
| `results/ipa_input.csv` | 25.3 KB | 715 | matches expected 715 |
| `results/ipa_input.csv.provenance.json` | 800 B | — |  |
| `results/ipa_input_full.csv` | 50.5 KB | 715 |  |
| `results/ipa_input_full.csv.provenance.json` | 805 B | — |  |
| `results/ipa_input_full.txt` | 50.5 KB | — |  |
| `results/ipa_input_significant.csv` | 65 B | **0 rows (expected -- 0/1938 proteins pass FDR<0.05 at n=2 technical replicates (DECISIONS_LOG D2))** |  |
| `results/onoff_proteins.csv` | 620 B | 10 | matches expected 10 |
| `results/qc/qc_boundaries.json` | 4.6 KB | — |  |
| `results/qc/qc_report.json` | 5.7 KB | — |  |
| `results/qc/qc_report.md` | 3.5 KB | — |  |
| `results/qc/quarantine_accessions.csv` | 33.3 KB | 2 |  |
| `results/qc_centering.csv` | 110 B | 1 |  |
| `results/qc_limma.csv` | 319.8 KB | 1,938 | matches expected 1,938 |
| `results/qc_limma.csv.provenance.json` | 800 B | — |  |
| `results/qc_limma_vanilla.csv` | 317.9 KB | 1,938 |  |
| `results/qc_replicate_correlation.csv` | 179 B | 1 |  |
| `results/replicate_correlation.png` | 99.7 KB | — |  |
| `results/single_condition_proteins.csv` | 28.8 KB | 604 | matches expected 604 |

All 7 headline row counts match the contract in `config/config.yaml`.

---

## c) Byte-freeze drift

Manifest: `proteomics_de/tests/expected/outputs.sha256` (79 files).

**79 OK · 0 CHANGED · 0 MISSING**

✅ **No drift.** Every frozen file is byte-identical to its baseline.

---

## d) Test coverage

| Test file | Test functions |
|---|---|
| `tests/test_accessions.py` | 14 |
| `tests/test_boundaries.py` | 36 |
| `tests/test_build_matrix.py` | 11 |
| `tests/test_build_report.py` | 21 |
| `tests/test_design.py` | 23 |
| `tests/test_e2e_smoke.py` | 22 |
| `tests/test_enrich_common.py` | 27 |
| `tests/test_foldchange_core.py` | 29 |
| `tests/test_freeze.py` | 7 |
| `tests/test_gating.py` | 22 |
| `tests/test_golden_outputs.py` | 27 |
| `tests/test_ipa_export.py` | 41 |
| `tests/test_limma_contract.py` | 28 |
| `tests/test_limma_r.py` | 37 |
| `tests/test_merge_guard.py` | 16 |
| `tests/test_provenance.py` | 30 |
| `tests/test_report_facts.py` | 24 |
| `tests/test_run_pipeline.py` | 40 |
| `tests/test_schemas.py` | 18 |
| `tests/test_string_ppi.py` | 36 |
| `tests/test_style_samples.py` | 18 |

**21 test files · 527 test functions.**

Module coverage — a module counts as covered if a test file is named after it (`test_<module>.py`) or mentions it in its filename:

| Pipeline module | Test file? |
|---|---|
| `centering_check.py` | ❌ none |
| `enrich/enrich_common.py` | ✅ yes |
| `enrich/gsea.py` | ❌ none |
| `enrich/network_figure.py` | ❌ none |
| `enrich/ora.py` | ❌ none |
| `enrich/string_ppi.py` | ✅ yes |
| `enrich/upset.py` | ❌ none |
| `etl/accessions.py` | ✅ yes |
| `etl/build_matrix.py` | ✅ yes |
| `etl/foldchange_core.py` | ✅ yes |
| `etl/merge_guard.py` | ✅ yes |
| `export/ipa_export.py` | ✅ yes |
| `foldchange.py` | ✅ yes |
| `gated/pca_cluster.py` | ❌ none |
| `limma_test.py` | ❌ none |
| `provenance.py` | ✅ yes |
| `qc/boundaries.py` | ✅ yes |
| `qc/schema.py` | ✅ yes |
| `qc/validate.py` | ❌ none |
| `replicate_check.py` | ❌ none |
| `report/build_facts.py` | ❌ none |
| `report/build_report.py` | ✅ yes |
| `viz/heatmap.py` | ❌ none |
| `viz/ma_plot.py` | ❌ none |
| `viz/qc_plots.py` | ❌ none |
| `viz/style.py` | ✅ yes |
| `viz/volcano.py` | ❌ none |

**13/27 pipeline modules have a matching test file.**

---

## e) Environment

* Python: `3.13.7` (`/Users/abrarshakik/Documents/Proteomics/.venv/bin/python`)
* pandas: `2.3.3`
* `.venv/`: ✅ present (`/Users/abrarshakik/Documents/Proteomics/.venv`)
* `pyproject.toml`: ✅ present
* `requirements-dev.txt`: ✅ present
* `requirements-lock.txt`: ✅ present
* Rscript: `Rscript (R) version 4.6.1 (2026-06-24)`
* R package `limma`: ✅ 3.68.4
* R package `imputeLCMD`: ✅ 2.1

