# ORA (g:Profiler) and STRING on Dr. Walker's 721-protein list

DECISIONS_LOG **D21**, list size updated by **D24**. Produced by `proteomics_de/enrich/walker_enrichment.py`
(config: `config.yaml: walker_enrichment`; pipeline stage `enrich_walker`).
Nothing under `results/enrichment/` (the 715-protein run) was changed; it is read-only here.

**History.** This folder first held a run on a 968-protein list (the 721 plus 247 "recovered" partial-tier
proteins). Dr. Walker then answered that the 247 should be dropped (D24), so the list is now 721 and
everything below was **re-run and re-computed** for it; no number was carried over from the 968 run. The 968
version (outputs and its raw API cache) is in git history, commits `b9c75b9` and `9e345d9`.

## Which list, exactly

`results/ipa_input_walker.csv` (D20, D24), **721 proteins**: 715 core + FRZB + 5 control-only proteins. Split as
the file's `regulated` column has it: **510 UP / 211 DOWN** (UP = 509 core + FRZB; DOWN = 206 core + 5 assigned).
UP and DOWN were queried **separately**, as in the 715 run. One gene symbol per protein, resolved by
`enrich_common.build_symbol_list` (first listed symbol, accession fallback).

Three things about the rows you must carry with any result:

* **The 5 assigned proteins (Egf, Ereg, Efna4, Tyro3, Mfap4) are not measured.** Their log2FC of -6.643856
  is an assignment (D20), not a measurement. **ORA and STRING use gene-list membership and direction
  only**: no fold change and no p-value of any kind enters either analysis. In the STRING node table
  their `log2FC` cell is blank. They sit in the DOWN list because they were detected in both control
  samples and in neither testosterone-treated sample. Since the two samples per group come from 2
  separate experiments (D23), that absence is reproduced across both experiments; it can still mean
  below the detection limit, and that is all the data says.
* **FRZB rests on the samples that were measured, with the missing value imputed** (limma MinProb;
  `n_imputed` 1). It is not FDR-significant, and **none of the 721 is** (0 of 1938 measured proteins pass
  FDR < 0.05; n = 2 biological replicates per group from 2 separate experiments, so low power; D2, D23).
* Everything below is **exploratory / hypothesis-generating**.

## Background (universe)

The same custom universe as the 715 run: the **measured proteome**, 2552 rows / **2530 unique symbols**
(union of `foldchange_all.csv` and `single_condition_proteins.csv`; also the IPA Reference Set,
`ipa_background_measured.txt`, D19). **It was checked, not assumed:** all 721 query symbols are inside it
(the 5 assigned proteins are rows of `single_condition_proteins.csv`, FRZB a row of `foldchange_all.csv`), so
it was **not extended** (`background_extended: false`, `n_query_not_in_background: 0` in `ora_meta.json`).
The code refuses to run if that stops being true. All 721 accessions (and the 716 of
`ipa_input_walker_without_assigned.csv`) are also in `ipa_background_measured.txt`.

## Versions, and why there is a same-version control

| | g:Profiler | STRING |
|---|---|---|
| Committed 715 run (`results/enrichment/`) | `e114_eg62_p19_27110d83` | 12.0 |
| **This run (721)** and the control | `e114_eg62_p19_9dc92ac8` | **STRING 12.5** |

**Both databases moved after the committed 715 run.** Comparing the 721 run to the committed 715 run
directly would confound "the list changed" with "the database changed". So the 715 list was **re-run in this
folder, with the same versions as the 721 run**: `control_715_same_version/`. Compare the 721 run to the
control, not to the committed files, when you want the effect of the five assigned proteins and FRZB. The
live queries for the 721 list on 2026-11-11 returned the **same two versions** as the control's
(2026-11-09) responses, so the control did not need refreshing and was kept as it was. A different STRING
release can be pinned by setting `walker_enrichment.string.api_base` to its stable address, then re-running
with `--fetch`.

## ORA (g:Profiler g:GOSt)

Same method as `enrich/ora.py`: organism `mmusculus`, GO:BP / GO:MF / GO:CC / KEGG / REAC, g:SCS correction,
`user_threshold` 0.05, `domain_scope=custom` with the background above. Sizes sent: **UP 510, DOWN 211**.
g:Profiler mapped **474 UP and 201 DOWN**; **31 UP and 7 DOWN symbols did not map** (renamed or obsolete
mouse symbols such as `Gnb2l1`, and Excel-mangled names such as `2018-09-09 00:00:00`) and 5 UP and 3 DOWN
were reported as ambiguous and discarded (510 = 474 + 31 + 5; 211 = 201 + 7 + 3; the failed symbols are listed,
the ambiguous ones counted, in `ora_meta.json`). The 715 control lost 31 + 5 of 509 and 7 + 3 of 206 the same
way, so FRZB and the five assigned proteins all mapped. This is inherited method, not changed here.

**Result: 0 terms pass g:SCS < 0.05 in either direction** (`ora_up.csv`, `ora_down.csv` are header-only),
the same null as the 715 run (D6). Best corrected p-values: UP 1.0; DOWN 0.897 ("dephosphorylation", the
same term and the same 11 hits as in the 715 runs, where it was 0.703; the DOWN query grew from 196 to 201
mapped genes, which dilutes it). Apart from that one DOWN term every corrected p-value in both directions is
floored at exactly 1.0, so the order of terms in `ora_top_terms_detail.json` (selected exactly as `ora.py`
selects, ties in API order) and in `ora_meta.json`'s "top5" lists is **not a ranking**.
`ora_top_terms_detail.json` lists the top 5 per direction for the UpSet-style use only.

### `comparison_vs_715.csv`: what the additions changed

One row per term that is significant, or in the **top 20 by rank, in any of three runs**: the committed 715
run, the same-version control, and the 721 run (47 terms: 20 UP, 27 DOWN). For each run: corrected p
(`*_adj_p`), the **uncorrected hypergeometric p** (`*_hypergeom_p`), the rank, the intersection size.

* **Rank** = position after sorting by (corrected p, uncorrected p, term id). The uncorrected p is
  **computed here** from each response's own counts (effective domain size, term size, query size,
  intersection size); it is a tie-breaker that makes a rank mean something when the corrected p is 1.0, not
  a significance claim.
* **Flags**: `sig_flag_*` compares the significant sets (blank for every row: nothing is significant
  anywhere); `top20_flag_*` compares top-20 membership. `gained` = in the 721 run only, `lost` = in the
  other run only, `kept` = in both. `*_vs_committed` compares with the committed 715 run, `*_vs_control` with the
  same-version control; **they are identical on every row here**, i.e. none of this movement is database drift.
* Counts of `top20` flags (vs control): UP 20 kept / 0 gained / 0 lost; DOWN 13 kept / 7 gained / 7 lost.
  The UP leads are the control's (translation initiation / ribosome / COPII / spliceosome terms; only the
  order of the last few shifts by one place). Only the DOWN leads move.

**Read the DOWN change with this in mind.** Seven of the 20 top DOWN leads of the 721 run contain at least one
of the five **assigned** proteins, and they are **exactly the seven terms that were "gained"** against the
control; the seven that were "lost" (ranks 13 to 20 in the control, for example `chemical homeostasis` and
`Chemokine signaling pathway`) were pushed just below rank 20 by them, to ranks 23 to 32.

| Rank (721) | Term | Assigned hits | Rank in the control |
|---|---|---|---|
| 5 | female sex differentiation | Ereg, Tyro3 | 152 |
| 8 | Colorectal cancer | Egf, Ereg | 215 |
| 9 | neuron migration | Tyro3 | 75 |
| 12 | Pancreatic cancer | Egf | 91 |
| 16 | regulation of cellular component organization | Egf, Ereg | 40 |
| 18 | rhythmic process | Ereg, Tyro3 | 173 |
| 20 | export from cell | Tyro3 | 45 |

* **Leave-the-five-out check (approximate).** Removing the five assigned proteins from the DOWN query leaves
  exactly the 715 control's DOWN list, so the recomputation can be checked against the control's own ranks. It
  was done offline from the cached evidence response (`gprofiler_down_all_evidence.json.gz`): uncorrected
  hypergeometric p with the intersection reduced by the assigned hits and the query size reduced by 5,
  background unchanged. It puts the seven terms at ranks 152, 215, 75, 91, 40, 173 and 45, which **equal the
  control's ranks**, so the seven leads exist only because the assigned values were added (a test pins this).
* **Four DOWN leads have no assigned hit** but still moved by a place or two (for example `female sex
  differentiation` pushed `inorganic ion homeostasis` from rank 5 to 6), because the assigned proteins changed
  the query size and the ranks of the terms above them.
* No ERBB/EGFR term is among the compared DOWN leads (Egfr is not on the 721 list).

None of this is significant: the only corrected p below 1.0 in the 721 run is 0.897 (DOWN
"dephosphorylation", which is not an assigned-driven term), so the ranks are tie-break aids for reading
`comparison_vs_715.csv`, not findings.

## STRING (v12.5, species 10090, required score 400)

Same recipe as `enrich/string_ppi.py` (STRING ids via `get_string_ids`, induced network, networkx degree /
betweenness / greedy-modularity communities); the 715 committed cache replayed through this code reproduces
`results/enrichment/string_edges.tsv` and `string_node_metrics.csv` byte for byte (a test).

| | seeds | mapped = nodes | edges | communities | isolated nodes |
|---|---|---|---|---|---|
| committed 715 (STRING 12.0) | 715 | 694 | 5,963 | 57 | n/a |
| control 715 (STRING 12.5) | 715 | 694 | 6,796 | 35 | 25 |
| **walker 721 (STRING 12.5)** | 721 | **700** | **6,870** | 35 | 25 |

* The database update alone added 833 edges to the 715 network (5,963 to 6,796) and merged communities (57 to
  35). So **do not compare the 721 network with the committed network**; compare with the control.
* The 721 network contains **exactly the control's 6,796 core-core edges** plus **74 edges that touch an added
  protein** (the five assigned proteins and FRZB). 700 of 721 map; the 21 that do not are all core
  protein-group accessions (listed in `string_meta.json`), the same 21 the control loses. Nodes: 694 core,
  5 assigned, 1 FRZB (partial tier). The largest community has 276 nodes (272 in the control).
* STRING edges are database knowledge (co-expression, literature, experiments), not evidence from this data;
  a protein's degree reflects how well studied it is.

### Where the six proteins sit (`string_focus_proteins.csv`)

| Protein | List status | Degree (rank of 700) | Community (size) | Neighbours |
|---|---|---|---|---|
| **Egf** | assigned | **51 (61st)** | 0 (276) | Ereg .874, Tgfb1 .826, Adam17 .817, Pcsk9 .812, Nes .803 ...; 33 UP / 18 DOWN, 49 core / 2 added |
| **Ereg** | assigned | 8 (441st) | 0 (276) | Egf .874, Adam17 .699, Hsp90aa1 .550, Memo1 .528, Lrp2 .515 ... |
| **Efna4** | assigned | 1 (653rd) | 6 (6) | Efnb1 .843 |
| **Tyro3** | assigned | 6 (496th) | 0 (276) | Adam17 .646, Actb .499, Anxa5 .444, Stat1 .442, Casp3 .416, Egf .403 |
| **Mfap4** | assigned | 8 (441st) | 0 (276) | Mfap2 .847, Efemp2 .652, Fbln1 .612, Col1a2 .586, Col6a2 .429 ... |
| **Frzb** | measured (partial) | 2 (624th) | 4 (25) | Man2c1 .431, Nid2 .408 |

All six are in the network; four (Egf, Ereg, Tyro3, Mfap4) are in the giant community 0. Egf is a mid-size hub
(degree 51; the largest hubs are Actb 153, Hsp90aa1 112, Rps3 111); Efna4 is nearly isolated (its only
neighbour is Efnb1). Egf's links to Ereg and Tyro3 are STRING's prior knowledge about the genes; they say nothing
about their measured behaviour here, and the DOWN direction of the five is an assignment.

## Files

| File | Contents |
|---|---|
| `ora_up.csv`, `ora_down.csv` | significant g:SCS terms (header-only: none) |
| `ora_meta.json` | versions, background and query sizes, unmapped symbols, best sub-threshold terms |
| `ora_top_terms_detail.json` | top 5 per direction with intersecting genes (selected like `ora.py`) |
| `string_edges.tsv`, `string_node_metrics.csv`, `string_meta.json` | network, per-node metrics (+`tier`, `value_source`), run facts |
| `string_focus_proteins.csv` | the six named proteins: degree, rank, community, neighbours with scores |
| `string_comparison_vs_715.json` | node/edge counts: committed 715, same-version 715, walker 721; all versions |
| `comparison_vs_715.csv` | term-level comparison described above |
| `control_715_same_version/` | the committed 715 list re-run under the same versions (same files, same `raw/`) |
| `raw/` | every API response verbatim (`gprofiler_*_all*.json.gz` gzipped) + `request_manifest.json` (request fingerprints) |

## Reproducing

`python run_pipeline.py --only enrich_walker` **replays `raw/` offline** and makes no network call; a missing
or mismatched cache is an error, not a silent live query. To query g:Profiler/STRING live for whatever is not
cached: `python proteomics_de/enrich/walker_enrichment.py --fetch` (it never overwrites a cached response; delete
`raw/` deliberately to refresh). The 721 responses were fetched live on 2026-11-11 with the user's approval (the
superseded 968 responses were removed from `raw/` first, since the fingerprints no longer matched), then replayed
twice, byte-identical.
