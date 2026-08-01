# `proteomics_de/` — module map

The SILAC differential-expression pipeline: 13 executable stages plus shared
support modules. This file is the map — what each script does, what it reads,
what it writes, and in what order. For the experiment itself and the n=2 caveat,
read the [repo README](../README.md) first.

**Companion documents.** Read them together; each answers a different question.

| File | Answers |
|---|---|
| [`DECISIONS_LOG.md`](DECISIONS_LOG.md) | *Why* — human decisions D1–D11, open questions, what still needs a person. |
| [`BUILD_LOG.md`](BUILD_LOG.md) | *How it got built* — append-only, per work package, with verification evidence. |
| [`STATUS.md`](STATUS.md) | *What exists right now* — generated from the filesystem by `tools/status.py`. Never hand-edited, always re-runnable. |
| [`../research1.md`](../research1.md) | *The design* — the original technical doc and its Section 6 build list (items 1–20), which every module traces back to. |

---

## Running it

```sh
.venv/bin/python run_pipeline.py --all      # from the repo root
.venv/bin/python run_pipeline.py --list     # the same table, as the runner sees it
```

`run_pipeline.py` holds the stage order as executable data, so this README and
the runner cannot drift apart in what they *claim* — but the runner is the
source of truth. Execution is serial by design (see its docstring; `--jobs` is a
deliberate non-goal).

---

## Run order

Stage 1 does the heavy lifting and chains three modules in-process. Stages 3–6
are mutually independent. 8 must follow 7, 10 must follow 9, 13 must be last.

```
 1  foldchange.py ─────────────────┐   (chains centering_check, replicate_check,
    │                              │    limma_test → limma_test.R in-process)
    ├──> 2  qc/validate.py         │
    │                              │
    ├──> 3  viz/qc_plots.py     ┐  │
    ├──> 4  viz/volcano.py      │  │  3–6 mutually independent
    ├──> 5  viz/ma_plot.py      │  │
    ├──> 6  viz/heatmap.py      ┘  │
    │                              │
    ├──> 7  enrich/string_ppi.py      [network]
    │        └──> 8  enrich/network_figure.py
    │
    ├──> 9  enrich/ora.py             [network]
    │        └──> 10 enrich/upset.py
    │
    ├──> 11 enrich/gsea.py            [network] [slow]
    │
    └──> 12 gated/pca_cluster.py
                                      │
    13 report/build_report.py  <──────┘  (reads every manifest; must be last)
```

`[network]` = live outbound API call (STRING, g:Profiler, Enrichr) with mouse
gene identifiers only, never intensities. Skip with `--skip-network`.

---

## The stages

### 1. `foldchange.py` — fold change + all of the DE core

The one stage that reads raw data. Also the only stage still **cwd-locked to the
repo root** (its input path and `RESULTS_DIR` are relative), which is why the
runner always invokes stages with `cwd=<repo root>`.

Fixes three correctness bugs in the legacy `proteomics_analysis.py`: log2 each
replicate ratio *before* averaging; symmetric ±log2(1.5) = ±0.585 thresholds;
never emit inf/NaN. It also replaces the legacy inner join with an outer merge,
so proteins seen in only one condition are kept rather than silently dropped.

* **In:** `Copy of General Sheet.xlsx`, sheets `Protein Report L` (control,
  intensities 31578 / 31580) and `Protein Report H` (treated, 31579 / 31581).
* **Out:** `results/foldchange_all.csv` (1948) · `single_condition_proteins.csv`
  (606) · `onoff_proteins.csv` (10) · `ipa_input.csv` (715 regulated).

It then chains three modules in-process, each of which only ever *adds* files:

| Chained module | Writes |
|---|---|
| `centering_check.py` | `results/qc_centering.csv`; on WARN also `foldchange_all_centered.csv` |
| `replicate_check.py` | `results/qc_replicate_correlation.csv`, `results/replicate_correlation.png` |
| `limma_test.py` → `limma_test.R` | `results/qc_limma.csv` (1938), `results/ipa_input_significant.csv` (**header-only by design**) |

`limma_test.py` owns the Python↔R subprocess boundary; the statistics live in
`limma_test.R` (MinProb imputation via `imputeLCMD`, then a limma moderated
t-test, per Peng et al. 2024). They talk over CSV files, not pipes, so the exact
matrix R saw stays auditable on disk: `_limma_input.csv`, `_limma_output.csv`,
`_limma_versions.txt`. If the R worker fails, it raises and writes nothing —
never partial p-values.

Only the 1938 genuine both-condition proteins are tested. The 10 ON/OFF
proteins are excluded on purpose: testing them would require inventing an
entire absent condition.

### 2. `qc/validate.py` (+ `qc/schema.py`)

pandera schema validation over the DE outputs — types, UniProt accession regex,
uniqueness, cardinality.

* **In:** the `results/*.csv` written by stage 1.
* **Out:** `results/qc/qc_report.json`, `results/qc/qc_report.md`.

### 3. `viz/qc_plots.py`

* **In:** `results/foldchange_all.csv`, `results/single_condition_proteins.csv`.
* **Out:** `results/figures/` — `sample_correlation`, `intensity_distributions`,
  `missing_values`, `rank_abundance` (each `.png` + `.svg`), and
  `figures_manifest.json`.

### 4. `viz/volcano.py`

log2FC vs −log10(raw p), coloured UP / DOWN / NO CHANGE at |log2FC| ≥ 0.585.
No FDR line is drawn, because nothing crosses it; the raw p = 0.05 line is shown
for orientation only.

* **In:** `results/qc_limma.csv`. · **Out:** `figures/volcano.{png,svg}`.

### 5. `viz/ma_plot.py`

Mean log2 intensity vs limma log2 fold change.

* **In:** `results/qc_limma.csv`. · **Out:** `figures/ma_plot.{png,svg}`.

### 6. `viz/heatmap.py`

Row-z-scored log2 intensities for the top 40 proteins **by smallest raw
p-value** (not by largest |log2FC| — see the script docstring), clustered by
protein only; sample columns stay in fixed control/treated order.

* **In:** `results/qc_limma.csv`. · **Out:** `figures/heatmap_top_de.{png,svg}`.

### 7. `enrich/string_ppi.py` `[network]`

STRING v12 REST, **species 10090 (mouse)**, `required_score=400`, seeded with
the 715 regulated proteins. 694/715 (97.1%) map. networkx supplies the graph
metrics (igraph was not needed).

* **In:** `results/foldchange_all.csv`.
* **Out:** `results/enrichment/string_node_metrics.csv` (694) ·
  `string_edges.tsv` · `string_meta.json` · raw API responses under
  `results/enrichment/raw/`.

### 8. `enrich/network_figure.py`

The publication network figure. Delivered as a static networkx render rather
than via py4cytoscape, which needs a running Cytoscape desktop (D3).

* **In:** `string_node_metrics.csv`, `string_edges.tsv`.
* **Out:** `figures/ppi_network.{png,svg}`, `figures/figures_manifest_network.json`.

### 9. `enrich/ora.py` `[network]`

g:Profiler g:GOSt over-representation, organism `mmusculus`, sources
GO:BP/MF/CC + KEGG + REAC, g:SCS correction — against a **custom background of
the 2554 detected proteins**, not the whole genome. That choice is the whole
point; see the header note below and DECISIONS_LOG D6.

* **In:** `results/foldchange_all.csv` + `single_condition_proteins.csv` (their
  union is the background).
* **Out:** `results/enrichment/ora_up.csv` and `ora_down.csv` (**both header-only
  by design**) · `ora_meta.json` · `ora_top_terms_detail.json` ·
  `raw/gprofiler_{up,down}.json` · `figures/ora_dotplot.{png,svg}`.

### 10. `enrich/upset.py`

Gene-level overlap across the top 5 UP-side and top 5 DOWN-side ORA leads.
Widened from a literal top-3+3, which had zero overlap of any kind and so
produced an empty plot.

* **In:** `results/enrichment/ora_top_terms_detail.json`.
* **Out:** `figures/upset.{png,svg}`.

### 11. `enrich/gsea.py` `[network] [slow]`

gseapy prerank (seed 42, 1000 permutations) over all 1938 `qc_limma.csv`
proteins ranked by `sign(log2FC) × −log10(p)`, against mouse Enrichr libraries
(GO_Biological_Process_2021 + KEGG_2019_Mouse). 568 gene sets survive the
min 15 / max 500 size filter. Minimum FDR q = 0.125 — nothing passes 0.05.

* **In:** `results/qc_limma.csv`.
* **Out:** `results/enrichment/gsea_results.csv` (568) · `gsea_meta.json` ·
  `figures/gsea_top.{png,svg}` · `figures/figures_manifest_enrich.json`.

### 12. `gated/pca_cluster.py`

PCA and hierarchical clustering behind an n ≥ 6 gate. **The gate fires**: this
dataset has 4 samples, so the output is explicitly QC-only, and WGCNA, UMAP/tSNE,
mixOmics sPLS-DA and VAE/GNN are skipped outright with reasons recorded. At n=4
the sample covariance matrix is rank-deficient (exactly 3 PCs exist by
construction) and no bootstrap support is possible for the dendrogram — so
neither is evidence of biological structure, and both say so on the figure.

* **In:** `results/qc_limma.csv`.
* **Out:** `results/gated/skip_log.csv`, `pca_coords.csv`, `pca_variance.csv` ·
  `figures/pca_qc.{png,svg}`, `figures/sample_dendrogram.{png,svg}` ·
  `figures/figures_manifest_gated.json`.

### 13. `report/build_report.py`

Assembles `report/report.html`: one self-contained file with every figure
inlined as a data URI and the fonts embedded, so it opens from disk with no
server and no sibling assets. Numbers come from `report_facts.json`, never from
prose, so nothing in the report is invented.

* **In:** `report/report_template.html`, `report/report_facts.json`, and all four
  `figures/figures_manifest*.json` plus the figures themselves.
* **Out:** `report/report.html`.

---

## Support modules (not stages)

| Module | Role |
|---|---|
| `viz/style.py` | Shared palette, `apply_style()`, `save_fig()`, and `record_manifest()`. Every figure goes through it. |
| `enrich/enrich_common.py` | Shared background/query construction and the enrichment figure manifest. |
| `config/config.yaml`, `config/sample_sheet.tsv` | The declarative design. **Descriptive only today** — the frozen scripts do not read them yet (D1). |
| `config/design.py`, `config/constants.py` | Readers for the above, and one home for values currently duplicated across nine files. Wiring them in is a later wave. |
| `etl/accessions.py` | One documented policy for UniProt accession fields, which today have three implicit and conflicting ones. |
| `export/ipa_export.py`, `provenance.py`, `qc/boundaries.py` | Wave-0 shims with real docstrings and no-op bodies. They write no files, deliberately, so the sha256 freeze stays green. |
| `tests/` | pytest suite. `tests/expected/frozen_counts.json` is the single source of expected row counts; `tests/expected/protected.sha256` is the byte-freeze manifest. |

---

## The three header-only outputs — do not "fix" these

Three files in `results/` contain a header and **zero data rows**. All three are
correct. A naive "the file is empty, something broke" reflex would destroy the
most defensible finding in this project.

| File | Rows | Why |
|---|---|---|
| `results/ipa_input_significant.csv` | 0 | 0 of 1938 proteins pass FDR < 0.05; the minimum adjusted p-value in the experiment is 0.305. The expected ceiling of n=2 **technical** replicates. **D2.** |
| `results/enrichment/ora_up.csv` | 0 | 0 GO/KEGG/Reactome terms survive g:SCS correction against the detected-proteome background. **D6.** |
| `results/enrichment/ora_down.csv` | 0 | Same. **D6.** |

The ORA pair carries the sharper point. The identical UP query against
g:Profiler's *default whole-genome* background returns **196 "significant"**
terms, top p = 1.9e-23 ("cytoplasm") — a textbook background-inflation artifact.
The honest detected-proteome background returns zero. Both were reproduced by a
fresh live re-query. The zero is the result.

`run_pipeline.py` encodes this directly: each output declares its expected shape,
and a table contracted to have 0 rows that has 0 rows reports
**`PASS_EMPTY_EXPECTED`**, printed distinctly with its reason. Emptiness is a
failure only where rows were expected. The inverse is also enforced: if one of
these ever becomes non-empty the runner FAILs and says so, because that is a
scientific event a human must sign off on, not a silent improvement.

---

## Reproducing the trend/robust limma variant

The committed baseline uses vanilla `eBayes()`. research1.md specifies
`trend=TRUE, robust=TRUE`, and D9 makes that the intended default; the vanilla
run is kept as the byte-reproducible baseline so the two stay comparable.

Since **D9** (2026-08-01) `trend_robust` is the DEFAULT, so the field-standard
variant is simply what `qc_limma.csv` now contains. What was the experiment is
now the shipped result, and the roles are reversed: it is *vanilla* that is
carried as the comparison.

`run_limma_test()` runs both on every invocation. Vanilla is written to
`results/qc_limma_vanilla.csv`, and the two are cross-checked in-process:

    logFC must be BIT-IDENTICAL between the flavours; only p-values may move.

That is D9's whole justification — `eBayes` moderates variances, it does not
refit coefficients — so it is asserted on every run rather than trusted. If
logFC ever moves, the run raises.

The conclusion does not change under either flavour: still **0/1938 significant**
at FDR<0.05. The minimum adjusted p-value is **0.116** (trend/robust, shipped)
versus **0.305** (vanilla), and the raw-p<0.05 count 55 versus 63 — the
refinement sharpens the extremes without crossing the line.

> The former `results/qc_limma_trend.csv` was **deleted**, not renamed. It
> predated the **D7** correction: its logFC is the exact negation of current
> values, so it described the inverted experiment. Under D9 the trend/robust
> result *is* `qc_limma.csv`.

---

## Two things to know before trusting a number

Both are open items in `DECISIONS_LOG.md`, and neither is a code defect.

* **D7 — the control/treated direction was inverted.** `proteomics_de/` inherited
  `31578, 31580 = control` from research1.md line 10. The lab's own earlier
  Pilot Project (`Pilot Project/Analysis/General Analysis/step1_data_cleaning.py`)
  says the opposite, and the 30 proteins shared between the two tables confirm
  it: correlation −0.82, sign agreement 2/30. The direction has been decided —
  the pilot is right and the pipeline is being flipped — but **the artifacts
  currently in `results/` are still the pre-flip orientation**. When the flip
  lands, every log2FC negates and the UP/DOWN counts swap to 509 UP / 206 DOWN.
  The 715-row IPA total is unchanged, and **p-values are unchanged**: swapping
  group labels negates logFC and t but leaves |t| and p invariant. So the
  headline "0 significant at FDR<0.05" survives the flip untouched. Read any
  directional claim in the committed results with D7 in hand.
* **D8 — the SILAC quantity is an open question.** Each sample carries its own
  `Intensity L`, `Intensity H` and `Ratio H/L` (verified: `Intensity L +
  Intensity H == Intensity` exactly). The pipeline uses combined `Intensity` —
  the sum of both isotope channels — and never touches the ratios; the two
  correlate at r = 0.066. All four samples show median log2(H/L) between −1.05
  and −1.53, the same direction every run, so this is not a reciprocal
  label-swap design. If Heavy is a common spike-in reference (super-SILAC), the
  ratio is the correct abundance measure and every result changes. **Unresolved;
  needs the lab.**

---

## Tests

```sh
.venv/bin/pytest                  # default: offline and fast
.venv/bin/pytest -m network       # live STRING / g:Profiler / Enrichr calls
.venv/bin/pytest -m r             # limma / imputeLCMD round-trip
.venv/bin/pytest -m golden        # sha256 byte-identity; LOCAL ONLY, never CI
```

The default `addopts` is `-m 'not network and not slow'`. `golden` must be
deselected in CI: a different BLAS or R build changes the last float digit,
which changes the CSV bytes, which changes the hash — an environment difference,
not a regression. Full detail in [`../ENVIRONMENT.md`](../ENVIRONMENT.md) §5.

Byte-freeze status on demand:

```sh
.venv/bin/python run_pipeline.py --verify-frozen   # per-file OK/CHANGED/MISSING
.venv/bin/python tools/status.py --check           # same verdict, inside STATUS.md
```

Both read `tests/expected/protected.sha256` through the same `freeze_check()`
implementation in `tools/status.py`, so they cannot disagree.
