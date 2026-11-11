# ORA (g:Profiler) and STRING on Dr. Walker's 968-protein list

DECISIONS_LOG **D21**. Produced by `proteomics_de/enrich/walker_enrichment.py`
(config: `config.yaml: walker_enrichment`; pipeline stage `enrich_walker`).
Nothing under `results/enrichment/` (the 715-protein run) was changed; it is read-only here.

## Which list, exactly

`results/ipa_input_walker.csv` (D20), **968 proteins**: 715 core + FRZB + 5 control-only proteins + 247
recovered partial-tier proteins. Split as the file's `regulated` column has it:
**661 UP / 307 DOWN** (UP = 509 core + FRZB + 151 recovered; DOWN = 206 core + 5 assigned + 96 recovered).
UP and DOWN were queried **separately**, as in the 715 run. One gene symbol per protein, resolved by
`enrich_common.build_symbol_list` (first listed symbol, accession fallback).

Three things about the rows you must carry with any result:

* **The 5 assigned proteins (Egf, Ereg, Efna4, Tyro3, Mfap4) are not measured.** Their log2FC of -6.643856
  is an assignment (D20), not a measurement. **ORA and STRING use gene-list membership and direction
  only**: no fold change and no p-value of any kind enters either analysis. In the STRING node table
  their `log2FC` cell is blank. They sit in the DOWN list because they were detected in both control
  readings and neither testosterone reading; that is all the data says.
* **FRZB and the 247 recovered proteins rest on a single biological replicate with an imputed
  replicate** (limma MinProb; `n_imputed` 1 or 2). They are not the "no problems" set Dr. Walker asked for
  (D20), and **none is FDR-significant** (0 of 968 pass FDR < 0.05; n = 2 technical replicates, D2).
* Everything below is **hypothesis-generating**.

## Background (universe)

The same custom universe as the 715 run: the **measured proteome**, 2552 rows / **2530 unique symbols**
(union of `foldchange_all.csv` and `single_condition_proteins.csv`; also the IPA Reference Set,
`ipa_background_measured.txt`, D19). **It was checked, not assumed:** all 968 query symbols are inside it
(the 5 assigned proteins are rows of `single_condition_proteins.csv`, the 248 partial-tier proteins rows of
`foldchange_all.csv`), so it was **not extended** (`background_extended: false`,
`n_query_not_in_background: 0` in `ora_meta.json`). The code refuses to run if that stops being true.
The 968 accessions are also all in `ipa_background_measured.txt`.

## Versions, and why there is a same-version control

| | g:Profiler | STRING |
|---|---|---|
| Committed 715 run (`results/enrichment/`) | `e114_eg62_p19_27110d83` | 12.0 |
| **This run (968)** and the control | `e114_eg62_p19_9dc92ac8` | **STRING 12.5** |

**Both databases moved after the committed 715 run.** Comparing the 968 run to the committed 715 run
directly would confound "the list changed" with "the database changed". So the 715 list was **re-run in this
folder, with the same versions as the 968 run**: `control_715_same_version/`. Compare the 968 run to the
control, not to the committed files, when you want the effect of the additions. (g:Profiler's drift did not
change any rank of the terms compared here; STRING's did: see below.) A different STRING release can be
pinned by setting `walker_enrichment.string.api_base` to its stable address, then re-running with `--fetch`.

## ORA (g:Profiler g:GOSt)

Same method as `enrich/ora.py`: organism `mmusculus`, GO:BP / GO:MF / GO:CC / KEGG / REAC, g:SCS correction,
`user_threshold` 0.05, `domain_scope=custom` with the background above. Sizes sent: **UP 661, DOWN 307**.
g:Profiler mapped **618 UP and 294 DOWN**; **37 UP and 9 DOWN symbols did not map** (renamed or obsolete
mouse symbols such as `Gnb2l1`, and Excel-mangled names such as `2018-09-09 00:00:00`) and 6 UP and 4 DOWN
were reported as ambiguous and discarded (661 = 618 + 37 + 6; 307 = 294 + 9 + 4; the failed symbols are listed, the ambiguous ones counted, in `ora_meta.json`). The
715 run lost 31 + 5 of 509 and 7 + 3 of 206 the same way. This is inherited method, not changed here.

**Result: 0 terms pass g:SCS < 0.05 in either direction** (`ora_up.csv`, `ora_down.csv` are header-only),
the same null as the 715 run (D6). Best corrected p-values: UP 1.0, DOWN 1.0 (the 715 run's best, DOWN
"dephosphorylation", was 0.703; it went to 1.0 here because 294 mapped DOWN genes instead of 196 dilute it).
Nearly every corrected p-value is floored at exactly 1.0, so the order of terms in
`ora_top_terms_detail.json` (selected exactly as `ora.py` selects, ties in API order) and in `ora_meta.json`'s
"top5" lists is **not a ranking**. `ora_top_terms_detail.json` lists the top 5 per direction for the UpSet-style use only.

### `comparison_vs_715.csv`: what the additions changed

One row per term that is significant, or in the **top 20 by rank, in any of three runs**: the committed 715
run, the same-version control, and the 968 run (66 terms: 30 UP, 36 DOWN). For each run: corrected p
(`*_adj_p`), the **uncorrected hypergeometric p** (`*_hypergeom_p`), the rank, the intersection size.

* **Rank** = position after sorting by (corrected p, uncorrected p, term id). The uncorrected p is
  **computed here** from each response's own counts (effective domain size, term size, query size,
  intersection size); it is a tie-breaker that makes a rank mean something when the corrected p is 1.0, not
  a significance claim.
* **Flags**: `sig_flag_*` compares the significant sets (blank for every row: nothing is significant
  anywhere); `top20_flag_*` compares top-20 membership. `gained` = in the 968 run only, `lost` = in the
  other run only, `kept` = in both. `*_vs_committed` compares with the committed 715 run, `*_vs_control` with the
  same-version control; **they are identical on every row here**, i.e. none of this movement is database drift.
* Counts of `top20` flags (vs control): UP 10 kept / 10 gained / 10 lost; DOWN 4 kept / 16 gained / 16 lost.
  UP leads stay translation initiation / ribosome / COPII / spliceosome terms. DOWN leads change almost completely.

**Read the DOWN change with this in mind.** 12 of the 20 top DOWN leads of the 968 run contain at least one of
the five **assigned** proteins, but they do not all depend on them equally:

* **Four leads exist only because the assigned proteins were added:** `ERBB2-EGFR signaling pathway`, `SHC1
  events in ERBB2 signaling`, `EGFR interacts with phospholipase C-gamma` and `ERBB2 Activates PTK6 Signaling`
  (ranks 11, 12, 13, 15). Each is a 3-gene term whose hits are exactly **Egf, Ereg, Egfr**, two of the three being
  assigned. Recomputing the ranking with the five assigned proteins removed from the DOWN query (an
  **approximate** recomputation done offline from the cached responses: uncorrected hypergeometric p with the
  intersection reduced by the assigned hits and the query size by 5, background unchanged) sends them to about
  rank 2,540-2,730 of 6,417 tested terms.
* **The other eight leads that contain an assigned protein move, but are also driven by core and recovered
  proteins.** `rhythmic process` (rank 1; 15 hits, 2 of them assigned, Ereg and Tyro3) would be about rank 4 without
  the five; `female sex differentiation` (2, Ereg and Tyro3) about 11; `positive regulation of Wnt signaling
  pathway` (3) about 10; `ovarian follicle development` (7) about 15; `regulation of nuclear division` (6; 8
  hits, two of them assigned) falls to about 141; `regulation of fibroblast proliferation` (16),
  `neuron cellular homeostasis` (17) and `neuron migration` (19), one assigned hit each, fall to about
  41, 63 and 36. In every case the drop is a change of rank among terms that are all at corrected p = 1.0.
* **The 96 recovered DOWN proteins also move the DOWN ranks.** The 968 DOWN list is the 206 core proteins plus
  the 5 assigned plus these 96; the change from the 715 run is the sum of all of them, not the assigned five alone.

None of this is significant: every corrected p is 1.0 in the 968 run (the 715 runs' only value below 1.0 was
0.703), so the ranks are tie-break aids for reading `comparison_vs_715.csv`, not findings.

## STRING (v12.5, species 10090, required score 400)

Same recipe as `enrich/string_ppi.py` (STRING ids via `get_string_ids`, induced network, networkx degree /
betweenness / greedy-modularity communities); the 715 committed cache replayed through this code reproduces
`results/enrichment/string_edges.tsv` and `string_node_metrics.csv` byte for byte (a test).

| | seeds | mapped = nodes | edges | communities | isolated nodes |
|---|---|---|---|---|---|
| committed 715 (STRING 12.0) | 715 | 694 | 5,963 | 57 | n/a |
| control 715 (STRING 12.5) | 715 | 694 | 6,796 | 35 | 25 |
| **walker 968 (STRING 12.5)** | 968 | **937** | **10,696** | 39 | 31 |

* The database update alone added 833 edges to the 715 network (5,963 to 6,796) and merged communities (57 to
  35). So **do not compare the 968 network with the committed network**; compare with the control.
* The 968 network contains **exactly the control's 6,796 core-core edges** plus **3,900 edges that touch an
  added protein**. 937 of 968 map (31 do not: 21 core, 10 partial; the control lost 21 core). Nodes: 694 core,
  238 partial (incl. FRZB), 5 assigned.
* STRING edges are database knowledge (co-expression, literature, experiments), not evidence from this data;
  a protein's degree reflects how well studied it is.

### Where the six proteins sit (`string_focus_proteins.csv`)

| Protein | List status | Degree (rank of 937) | Community (size) | Neighbours |
|---|---|---|---|---|
| **Egf** | assigned | **66 (68th)** | 0 (299) | Egfr .999, Ereg .874, Tgfb1 .826, Adam17 .817, Pcsk9 .812 ...; 44 UP / 22 DOWN, 49 core / 17 added |
| **Ereg** | assigned | 11 (548th) | 0 (299) | Egfr .995, Egf .874, Adam17 .699, Hsp90aa1 .550, Memo1 .528 ... |
| **Efna4** | assigned | 2 (841st) | 5 (7) | Egfr .908, Efnb1 .843 |
| **Tyro3** | assigned | 8 (633rd) | 0 (299) | Gas6 .994, Adam17 .646, Actb .499, Anxa5 .444, Stat1 .442 ... |
| **Mfap4** | assigned | 16 (422nd) | 0 (299) | Mfap2 .847, Fbln5 .796, Efemp2 .652, Fbln1 .612, Col1a2 .586 ... |
| **Frzb** | measured (partial) | 3 (801st) | 0 (299) | Man2c1 .431, Nid2 .408, Ogn .400 |

All six are in the network; five are in the giant community 0. Egf is a mid-size hub (degree 66; the largest
hubs are Actb 187, Hsp90aa1 135, Rps3 133), Efna4 is nearly isolated (its only neighbours are Egfr and Efnb1). Egf's and
Ereg's links to Egfr and each other are STRING's prior knowledge about the genes; they say nothing about their
measured behaviour here, and the DOWN direction of the five is an assignment.

## Files

| File | Contents |
|---|---|
| `ora_up.csv`, `ora_down.csv` | significant g:SCS terms (header-only: none) |
| `ora_meta.json` | versions, background and query sizes, unmapped symbols, best sub-threshold terms |
| `ora_top_terms_detail.json` | top 5 per direction with intersecting genes (selected like `ora.py`) |
| `string_edges.tsv`, `string_node_metrics.csv`, `string_meta.json` | network, per-node metrics (+`tier`, `value_source`), run facts |
| `string_focus_proteins.csv` | the six named proteins: degree, rank, community, neighbours with scores |
| `string_comparison_vs_715.json` | node/edge counts: committed 715, same-version 715, walker 968; both versions |
| `comparison_vs_715.csv` | term-level comparison described above |
| `control_715_same_version/` | the committed 715 list re-run today (same files, same `raw/`) |
| `raw/` | every API response verbatim (`gprofiler_*_all*.json.gz` gzipped) + `request_manifest.json` (request fingerprints) |

## Reproducing

`python run_pipeline.py --only enrich_walker` **replays `raw/` offline** and makes no network call; a missing
or mismatched cache is an error, not a silent live query. To query g:Profiler/STRING live for whatever is not
cached: `python proteomics_de/enrich/walker_enrichment.py --fetch` (it never overwrites a cached response; delete
`raw/` deliberately to refresh). These files were fetched live on 2026-11-09 with the user's approval, then
replayed twice, byte-identical.
