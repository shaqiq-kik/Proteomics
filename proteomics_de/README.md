# `proteomics_de/` — module map

The SILAC differential-expression pipeline: 13 executable stages plus shared
support modules. This file is the map — what each script does, what it reads,
what it writes, and in what order. For the experiment itself and the n=2 caveat
(2 biological replicates per group from 2 separate experiments, low power; D23),
read the [repo README](../README.md) first.

**Companion documents.** Read them together; each answers a different question.

| File | Answers |
|---|---|
| [`DECISIONS_LOG.md`](DECISIONS_LOG.md) | *Why* — human decisions D1–D23, open questions, what still needs a person. |
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

* **In:** `Copy of General Sheet.xlsx`, sheets `Protein Report L` (treated,
  intensities 31578 / 31580) and `Protein Report H` (control, 31579 / 31581).
  Sheet residency is not condition — which sheet is which is set by
  `config/sample_sheet.tsv` alone, per **D7**.
* **Out:** `results/foldchange_all.csv` (1948) · `single_condition_proteins.csv`
  (604) · `onoff_proteins.csv` (10) · `ipa_input.csv` (715 regulated).

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
min 15 / max 500 size filter. Minimum FDR q = 0.127 — nothing passes 0.05.

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
| `viz/style.py` | Shared palette, `apply_style()`, `save_fig()`, and `record_manifest()`. Every figure goes through it; its sample maps are derived from `config/sample_sheet.tsv`. |
| `enrich/enrich_common.py` | Shared background/query construction and the enrichment figure manifest. |
| `config/sample_sheet.tsv` | The design. **Load-bearing**, not descriptive — `config/design.py` reads it and drives the limma design matrix, `viz/style.py`'s sample maps, and `gated/pca_cluster.py`'s dispatcher. `foldchange.py`, `replicate_check.py` and `qc/schema.py` assert-match it and fail loudly on drift rather than deriving from it (D1: their logic is inherently 2-channel-SILAC-specific). `config/config.yaml` stays descriptive documentation only. |
| `config/design.py`, `config/constants.py` | The design reader, and one home for thresholds/seed/caveat text previously duplicated across several files. |
| `etl/accessions.py` | One documented policy for UniProt accession fields, replacing three implicit and conflicting ones. |
| `export/ipa_export.py` | Builds `ipa_input_full.csv`/`.txt` (limma p-value + FDR joined onto the regulated set) and validates all four IPA deliverables. Pipeline stage 2. |
| `export/regulated_lists.py` | Splits `ipa_input_full.csv` into `regulated_up.csv`/`regulated_down.csv` with a linear `fold_change` column, sorted by descending magnitude of change. Pipeline stage 3. |
| `export/supplementary_lists.py` | Two more files closing the gap `regulated_up/down.csv` structurally cannot: `regulated_up_partial.csv`/`regulated_down_partial.csv` (152/96 tier-3 proteins reclassified from limma's imputed log2FC) and `qualitative_changes.csv` (614 tier-1/2 proteins never identified on one side at all, with a `direction` call but no fold change). Pipeline stage 4. **D17.** |
| `export/ipa_extended.py` | Carries D17's two supplementary tiers into the QIAGEN uploads, which the core export structurally cannot see: `ipa_input_extended.csv`/`.txt` (963 rows, `tier`-tagged) plus `ipa_qualitative_up.txt`/`ipa_qualitative_down.txt` (372/242 bare accessions). Reads only; re-hashes the six byte-frozen IPA/regulated files around every build. Pipeline stage 5. **D18.** |
| `export/walker_additions.py` | Dr. Walker's list (**D20**): the 715 core rows + FRZB (real values) + 5 control-only proteins with an **assigned** 100-fold decrease (`value_source=assigned`, p-values blank): `ipa_input_walker.csv`/`.txt` (**721 rows**, D24; the 247 recovered proteins of the first 968-row version were dropped at Dr. Walker's request, that version is in git history) and `walker_additions_manifest.csv` (6 rows), plus `ipa_input_walker_without_assigned.csv`/`.txt` (716 rows, the IPA sensitivity dataset). Driven by `config.yaml: walker_additions`. Pipeline stage `ipa_export_walker`. |
| `enrich/walker_enrichment.py` | ORA (g:Profiler) + STRING over Dr. Walker's 721-protein list (**D21**, list size **D24**), additive: writes `results/enrichment_walker/` (and a same-version re-run of the 715 list, `control_715_same_version/`, plus `comparison_vs_715.csv`), reusing `ora.py`/`string_ppi.py` pieces and never touching `results/enrichment/`. **Replays its committed `raw/` API cache offline**; live calls only with `--fetch`. Driven by `config.yaml: walker_enrichment`. Pipeline stage `enrich_walker`. |
| `provenance.py` | Writes `<name>.provenance.json` sidecars (caveat text, git commit, sha256, row count, tool versions) for the CSV deliverables. |
| `qc/boundaries.py` | Real pandera validation at `foldchange.py`'s load/merge/foldchange boundaries; `after_load` is permissive (junk accessions route to quarantine), later stages are strict. Appends to `results/qc/qc_boundaries.json`. |
| `tests/` | pytest suite. `tests/expected/frozen_counts.json` is the single source of expected row counts; `tests/expected/outputs.sha256` is the byte-freeze manifest (scientific outputs only — see `tools/freeze.py`). |

---

## The three header-only outputs — do not "fix" these

Three files in `results/` contain a header and **zero data rows**. All three are
correct. A naive "the file is empty, something broke" reflex would destroy the
most defensible finding in this project.

| File | Rows | Why |
|---|---|---|
| `results/ipa_input_significant.csv` | 0 | 0 of 1938 proteins pass FDR < 0.05; the minimum adjusted p-value in the experiment is 0.116. The expected ceiling at n=2 per group (low power; 2 biological replicates from 2 separate experiments, **D23**). **D2.** |
| `results/enrichment/ora_up.csv` | 0 | 0 GO/KEGG/Reactome terms survive g:SCS correction against the detected-proteome background. **D6.** |
| `results/enrichment/ora_down.csv` | 0 | Same. **D6.** |

The ORA pair carries the sharper point. The identical UP query against
g:Profiler's *default whole-genome* background returns **326 "significant"**
terms, top p = 9.0e-70 ("cytoplasm"); the DOWN query returns **196**, top
p = 1.9e-23 (same term) — a textbook background-inflation artifact either way.
The honest detected-proteome background returns zero for both. Both were
reproduced by a fresh live re-query. The zero is the result.

`run_pipeline.py` encodes this directly: each output declares its expected shape,
and a table contracted to have 0 rows that has 0 rows reports
**`PASS_EMPTY_EXPECTED`**, printed distinctly with its reason. Emptiness is a
failure only where rows were expected. The inverse is also enforced: if one of
these ever becomes non-empty the runner FAILs and says so, because that is a
scientific event a human must sign off on, not a silent improvement.

---

## Two supplementary files, and why they exist

`regulated_up.csv`/`regulated_down.csv` are gated on `complete=True` --
every one of the 4 raw replicate intensities present and non-zero -- decided
before limma ever runs. That gate is correct for the frozen QIAGEN
deliverable, but it silently drops two tiers of real signal: 360 proteins
with 1-2 missing replicates that limma still tests after MinProb imputation
(FRZB, GAS6, LUM, SLIT3 among them), and 614 proteins never identified at
all on one side (EGF, EREG among them). Both tiers were the subject of a
direct professor report ("EGF, EREG, FRZB are missing"). D17 adds two files
rather than loosening the gate, keeping `ipa_input.csv`/`ipa_input_full.csv`/
`regulated_up.csv`/`regulated_down.csv` byte-for-byte untouched:

| File | Rows | Source | What it is NOT |
|---|---|---|---|
| `results/regulated_up_partial.csv` / `..._down_partial.csv` | 152 / 96 | `qc_limma.csv` filtered to `n_imputed > 0`, reclassified at `LOG2_THRESHOLD` | Not complete-case data -- every row carries `n_imputed >= 1`. |
| `results/qualitative_changes.csv` | 614 (372 UP / 242 DOWN) | `single_condition_proteins.csv` + `onoff_proteins.csv`, unioned | Has no fold-change or p-value column, on purpose -- neither is valid when one whole condition has zero data points or zero variance. |

See DECISIONS_LOG D17 for the full evidence and the pilot-panel reconciliation.

Dr. Walker's follow-up decision on these gaps (FRZB added with its real values, EGF/EREG plus three more control-only proteins added with an *assigned* 100-fold decrease) is in the section "Dr. Walker's list (D20)" below.

---

## What to upload to QIAGEN IPA

D17 wrote the two supplementary files but nothing carried them to QIAGEN: every
IPA file is built through the same `complete=True` gate, so the 862 proteins
D17 recovered had **zero** overlap with the 715 in `ipa_input.csv` —
Frzb/Egf/Ereg included. D18 adds three uploads rather than editing the frozen
ones, and registers `ipa_extended.py` as a pipeline stage so they cannot go
stale the way `ipa_input_full.csv` once did (**D13**).

| Upload | Rows | What it is |
|---|---|---|
| `results/ipa_input_extended.txt` (`.csv` twin) | 963 | The quantitative dataset: the 715 core rows **verbatim** from `ipa_input_full.csv`, plus the 248 tier-3 partial rows, with `log2FC`, `p_value`, `adj_p_value` and a `tier` column. |
| `results/ipa_qualitative_up.txt` | 372 | ID-only list — gained after testosterone. Run as its own Core Analysis. |
| `results/ipa_qualitative_down.txt` | 242 | ID-only list — lost after testosterone. Same. |

The three are pairwise disjoint and cover 1577 distinct accessions. The older
`ipa_input.csv`/`ipa_input_full.csv` remain valid and byte-identical; the
extended file is a superset of the latter, not a replacement for it.

> **The `tier` column is not cosmetic.** In `ipa_input_extended`, `log2FC`
> carries **two different quantities**: for `tier=core` it is the raw
> mean-of-log2-ratios on unimputed data, for `tier=partial` it is limma's
> estimate after MinProb imputation of 1–2 of the 4 replicates. Both are
> legitimate IPA input; reading a partial row as a complete-case measurement is
> not. Say so in any methods write-up, and split on `tier` when it matters.

The two qualitative lists carry no fold change and no p-value on purpose —
neither is defined when one whole condition has zero data points. IPA accepts
an identifier-only list, so over-representation still runs, but **no activation
z-scores come back**, since those are driven by fold-change direction. Inventing
a sentinel ±1 log2FC to force them was considered and rejected (**D18**).

Still open: ORA and STRING are seeded with the same 715 and would need widening
to the 963 to match these uploads. GSEA already covers all 1938. Actually
running the analyses in IPA is manual — **D4**.

---

## Dr. Walker's list (D20): measured vs ASSIGNED values

Dr. Walker reviewed the 715-protein list (the professor's report that EGF, EREG
and FRZB were missing, D17/D18) and decided what to add. The result is a **new**
upload, `results/ipa_input_walker.txt` (`.csv` twin), alongside the unchanged
older files:

| Step | Rows | What it is |
|---|---|---|
| core regulated list | 715 | `ipa_input_full.csv`, verbatim (`value_source=measured`) |
| + FRZB P97401 | 716 | its **real** values: log2FC +4.407 (~21x UP), p 0.0203, adj.P 0.998, one imputed replicate (`value_source=measured`) |
| + 5 control-only proteins | 721 | EGF, EREG, Efna4, Tyro3, Mfap4: **DOWN, assigned 100-fold decrease** (`value_source=assigned`) |
| *(D24)* the 247 recovered partial-tier proteins of the first 968-row version | dropped | Dr. Walker's final answer: drop them. The 721 is the final list. |

`results/walker_additions_manifest.csv` (6 rows: FRZB and the five assigned) lists every addition with its
source, reason, `value_source` and `n_imputed`. **D24:** the list is now **721 rows (510 UP / 211 DOWN)**; the
same list without the five assigned rows (716 rows, 510 UP / 206 DOWN) is
`results/ipa_input_walker_without_assigned.txt`/`.csv`, the dataset for the IPA run "without the assigned
values". The 968-row version (with the 247 recovered proteins) is in git history (commits `b9c75b9`, `9e345d9`).

> **Methods note -- the assigned values are not measurements.** EGF, EREG, Efna4,
> Tyro3 and Mfap4 were detected in both control readings and in neither
> testosterone reading, so no fold change can be computed from the data. By
> Dr. Walker's decision each is **assigned** a 100-fold decrease (fold change
> 0.01, log2FC = -6.643856). `value_source = assigned` marks these rows, and
> `p_value`/`adj_p_value` are blank on purpose. "Not detected" may mean *below the
> detection limit*, so the true change is unknown and may be smaller or larger
> than 100-fold. Report them as "lost or strongly reduced after testosterone",
> never as a measured 100-fold change. Unlike the ID-only lists above, these rows
> carry a fold-change direction, so IPA can include them in activation z-scores;
> five values at -6.64 can dominate a z-score built from mostly-small real ones.

**The "~860 recovered proteins"** are the 862 that D17/D18 recovered (248 partial
+ 614 qualitative). Dr. Walker asked for only those altered by at least 50% with
"no problems". The strict rule (`|log2FC| >= 0.585`, seen in both replicates of
both groups, `n_imputed == 0`) admits **none**: every partial-tier protein is in
that tier *because* a reading is missing and imputed, and no qualitative protein has
a fold change (862 -> 856 not already listed -> 247 clear the cutoff -> 0 -> 0). D20
had widened the rule to admit those 247 (each rests on the measured samples with an
imputed value; none FDR-significant). **Dr. Walker's final answer (D24) was to drop
them**, so both "no problems" switches
(`walker_additions.recovered_rule.require_all_replicates_seen` /
`require_no_imputed_values`) are `true` again = the strict rule = his original
wording, and zero recovered proteins are added. The rule machinery is kept: setting
both switches to `false` re-admits the 247 (the tests cover that path). See DECISIONS_LOG D20, D24.

Everything is generated by `run_pipeline.py --only ipa_export_walker`; no existing
result file is modified, and the builder re-hashes them around every run.

---

## ORA and STRING on Dr. Walker's list (D21)

Dr. Walker wants his list (721 proteins since D24; D20) run through IPA, g:Profiler and STRING. IPA is a manual
upload (steps in `results/ipa_analysis_walker/UPLOAD_STEPS.md`). g:Profiler and STRING are done, **additively**,
in `results/enrichment_walker/` (see its `README.md` for the numbers); the 715-based `results/enrichment/`
is untouched and still what `enrich_ora` / `enrich_string_ppi` produce.

* **List and background.** UP 510 / DOWN 211, queried separately (g:Profiler mapped 474 and 201). The background is
  the same measured proteome as the 715 run (2530 unique symbols); all 721 query symbols were verified to be inside
  it, so it was not extended, and the code refuses to run if that stops holding.
* **Membership only.** The 5 assigned proteins (Egf, Ereg, Efna4, Tyro3, Mfap4) enter ORA/STRING as DOWN
  gene-list members; their assigned -6.64 is never used, STRING leaves their `log2FC` blank, and no p-value
  reaches any output. FRZB rests on the measured samples with an imputed value; none of the 721 is FDR-significant.
* **Result.** 0 g:SCS-significant terms in either direction (the D6 null again). STRING: 700/721 mapped,
  700 nodes, 6,870 edges (the control's 6,796 core-core edges plus 74 touching FRZB or an assigned protein).
  `comparison_vs_715.csv` shows which leads came and went; 7 of the top 20 DOWN leads contain an assigned
  protein, and they are exactly the 7 leads that are new against the same-version control (approximate
  leave-the-five-out ranks equal the control's). Nothing is significant (the only corrected p below 1.0 is 0.897),
  so ranks are tie-break aids only.
* **Version drift.** The committed 715 run used STRING 12.0 and g:Profiler `e114_eg62_p19_27110d83`; this run
  uses STRING 12.5 and `e114_eg62_p19_9dc92ac8` (unchanged since the first, 968-protein run, so the control was not
  refreshed). STRING 12.5 alone adds 833 edges to the 715 network (5,963 to 6,796), so the 715 list was re-run
  under these versions (`control_715_same_version/`) and **the 721 run must be compared with that control**, not
  with the committed files.
* **Superseded 968 run.** The first version of this analysis covered 968 proteins; it is in git history
  (commits `b9c75b9`, `9e345d9`). See DECISIONS_LOG D24.
* **Offline by default.** `python run_pipeline.py --only enrich_walker` replays the committed `raw/` cache
  and makes no network call (a missing/mismatched cache is an error). `python
  proteomics_de/enrich/walker_enrichment.py --fetch` queries the live APIs for anything not cached.

---

## IPA on Dr. Walker's list (D22)

**D24: the two runs below were made on the earlier 968-protein list and are SUPERSEDED by the 721-protein list
(the 247 recovered proteins were dropped at Dr. Walker's request). New IPA runs on the 721 list, and on the 716-row
`ipa_input_walker_without_assigned.txt`, are pending (manual).** The old exports and figures are kept for the
record; their numbers describe the 968 list, not the current one.

IPA is a manual upload; the exports are in `results/ipa_analysis_walker/` (outside the byte-freeze) and its
`README.md` has the numbers. Two runs were made on 2026-11-10 (968 list):

* **Run 1, `assigned_dropped`:** `p_value` / `adj_p_value` mapped. IPA drops rows that are blank in a mapped column,
  so the 5 assigned proteins (blank p-values) were mapped but never analysed: 929 analysis-ready molecules, the
  same as the 963-protein dataset without them. It serves as the "without the assigned values" analysis.
* **Run 2, `assigned_included`:** only `log2FC` mapped, no cutoffs: 934 analysis-ready molecules, the 5 included.
  This is the run Dr. Walker asked for.
* **Reference set.** IPA refused the measured-proteome reference (overlap too large), so both runs use IPA's
  default Knowledge Base reference. IPA p-values are inflated by that; compare rankings and z-scores, never
  p-values, and not with the in-house ORA.
* **What the 5 assigned values change.** Adding them changes 60 pathway, 156 regulator and 393 function z-scores;
  6 pathways, 16 regulators and 79 functions with abs z >= 2 list an assigned protein (e.g. Organismal death, z -6.08).
  Those z-scores rest partly on the assigned -6.64 placeholder and must not be quoted as findings. The androgen
  entries barely move (AR z 2.08 and metribolone z 3.38 in both runs).
* 0 of the proteins pass FDR < 0.05 (as for the 968 and the 721 list); every IPA number is exploratory. See DECISIONS_LOG D22, D24.

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

## Two things we caught, and had to check with the lab

Both are corrected-and-verified items now. Both are in `DECISIONS_LOG.md`,
and neither was a code defect.

* **D7 — the control/treated direction was inverted, and has been corrected.**
  `proteomics_de/` originally inherited `31578, 31580 = control` from
  research1.md line 10. The lab's own earlier Pilot Project
  (`Pilot Project/Analysis/General Analysis/step1_data_cleaning.py`) says the
  opposite, and the 30 proteins shared between the two tables confirmed it:
  correlation −0.82, sign agreement 2/30 pre-fix. The pilot was right;
  `config/sample_sheet.tsv` now carries the corrected orientation, and the fix
  was verified as a pure sign inversion (foldchange and limma logFC negate to
  within float noise; p-values are exactly unchanged, as they must be —
  swapping group labels negates logFC and t but leaves |t| and p invariant).
  Every artifact in `results/` reflects the corrected orientation: UP/DOWN is
  **509 UP / 206 DOWN** (post-fix correlation with the pilot: +0.82, sign
  agreement 93.3%), the 715-row IPA total is unchanged, and the headline
  "0 significant at FDR<0.05" is untouched. Full evidence in
  `DECISIONS_LOG.md` D7.
* **D8 — the SILAC quantity question, resolved.** Each sample carries its own
  `Intensity L`, `Intensity H` and `Ratio H/L` (verified: `Intensity L +
  Intensity H == Intensity` exactly). The pipeline uses combined `Intensity` —
  the sum of both isotope channels — and never touches the ratios; the two
  correlate at r = 0.066. The professor confirmed the SILAC metabolic-labelling
  step wasn't actually completed during the experiment, so `Ratio H/L` carries
  no real signal and was never a valid alternative. Since the pipeline has
  never used it, **nothing about the results changes.** Full evidence in
  `DECISIONS_LOG.md` D8.

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

Both read `tests/expected/outputs.sha256` through the same `freeze_check()`
implementation in `tools/status.py`, so they cannot disagree.
