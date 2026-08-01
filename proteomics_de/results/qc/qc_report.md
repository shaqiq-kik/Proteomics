# proteomics_de QC validation report

Generated: 2026-08-01T06:10:20.050610+00:00
pandera 0.32.1 / pandas 2.3.3 / Python 3.13.7

**Overall: PASS** (9/9 files passed)

This report validates the frozen, committed pipeline outputs under
`proteomics_de/results/` against the pandera schemas in
`proteomics_de/qc/schema.py`. It is READ-ONLY: no pipeline output was
modified to produce this report.

## Per-file results

| File | Rows | Expected | Row count OK | Schema | Failures |
|---|---|---|---|---|---|
| foldchange_all.csv | 1948 | 1948 | True | PASS | 0 |
| qc_limma.csv | 1938 | 1938 | True | PASS | 0 |
| ipa_input.csv | 715 | 715 | True | PASS | 0 |
| single_condition_proteins.csv | 604 | 604 | True | PASS | 0 |
| onoff_proteins.csv | 10 | 10 | True | PASS | 0 |
| de/intensity_matrix.tsv | 1938 | 1938 | True | PASS | 0 |
| de/design.tsv | 4 | 4 | True | PASS | 0 |
| de/limma_results.tsv | 1938 | 1938 | True | PASS | 0 |
| ipa_input_full.csv | 715 | 715 | True | PASS | 0 |

## Cross-file consistency checks

| Check | Detail | Passed |
|---|---|---|
| foldchange_all.csv rows minus onoff_proteins.csv rows equals qc_limma.csv rows | 1948 - 10 == 1938 | True |
| foldchange_all.csv rows plus single_condition_proteins.csv rows equals the detected-proteome background size used by enrichment (tests/expected/frozen_counts.json background_union) | 1948 + 604 == 2552 (expected 2552) | True |
| every quarantined accession carries a reason (results/qc/quarantine_accessions.csv, DECISIONS_LOG D11) | 2 quarantined row(s), all with a reason: True | True |

## Documented data quirks accommodated by these schemas

- **UniProt accession regex extended for O/P/Q-prefixed accessions.** The literal task-spec regex `^[A-NR-Z0-9]([A-Z0-9]{5}|[A-Z0-9]{9})(-\d+)?$` excludes accessions starting with O, P, or Q, which are common and legitimate in this dataset (e.g. P19137, Q9JHU4, O08528). `schema.py`'s `UNIPROT_TOKEN_RE` adds the missing `[OPQ][A-Z0-9]{5}` branch.
- **Multi-protein-group accessions** (';'-joined UniProt IDs, e.g. `P05132;P68181`) occur in foldchange_all.csv (48 rows), qc_limma.csv (48), ipa_input.csv (21), single_condition_proteins.csv (67). Each semicolon-delimited token is validated independently.
- **NaN gene names** occur in foldchange_all.csv (8), qc_limma.csv (8), ipa_input.csv (2), single_condition_proteins.csv (15) -- unannotated or ambiguous UniProt entries. Gene columns are nullable.
- **2 numeric-junk 'accession' rows are QUARANTINED, not excepted** (DECISIONS_LOG D11). They are ';'-joined lists of bare MaxQuant row indices (32,759 and 681 characters), an upstream raw-data artifact. An earlier version of `schema.py` carved an exception so they would PASS validation; that exception is deleted. They are now removed from `single_condition_proteins.csv` (606 -> 604) and recorded in full, with a reason, in `results/qc/quarantine_accessions.csv`. Discrimination is by TOKEN SHAPE (every ';'-token is a bare integer) via `etl/accessions.py::is_junk_index_list`, never by length -- the two 69-character `P08752;P20612;...` rows are legitimate protein groups and are retained.

## Stage-boundary validation

`qc/boundaries.py` validates the frames crossing each pipeline stage and records the outcome in `results/qc/qc_boundaries.json`. `after_load` and `after_merge` are PERMISSIVE (they see raw MaxQuant data, where a malformed accession is a fact about the input and is routed to quarantine); `after_foldchange` and `before_ipa_export` are STRICT (they see frames this repo produced, where a schema failure is a bug and stops the run).
