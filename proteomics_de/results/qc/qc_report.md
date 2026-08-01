# proteomics_de QC validation report

Generated: 2026-08-01T04:23:17.468685+00:00
pandera 0.32.1 / pandas 2.3.3 / Python 3.13.7

**Overall: PASS** (5/5 files passed)

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
| single_condition_proteins.csv | 606 | 606 | True | PASS | 0 |
| onoff_proteins.csv | 10 | 10 | True | PASS | 0 |

## Cross-file consistency checks

| Check | Detail | Passed |
|---|---|---|
| foldchange_all.csv rows minus onoff_proteins.csv rows equals qc_limma.csv rows | 1948 - 10 == 1938 | True |
| foldchange_all.csv rows plus single_condition_proteins.csv rows equals the detected-proteome background size used by enrichment (config.yaml enrichment.gprofiler.background_n) | 1948 + 606 == 2554 (expected 2554) | True |

## Documented data quirks accommodated by these schemas

- **UniProt accession regex extended for O/P/Q-prefixed accessions.** The literal task-spec regex `^[A-NR-Z0-9]([A-Z0-9]{5}|[A-Z0-9]{9})(-\d+)?$` excludes accessions starting with O, P, or Q, which are common and legitimate in this dataset (e.g. P19137, Q9JHU4, O08528). `schema.py`'s `UNIPROT_TOKEN_RE` adds the missing `[OPQ][A-Z0-9]{5}` branch.
- **Multi-protein-group accessions** (';'-joined UniProt IDs, e.g. `P05132;P68181`) occur in foldchange_all.csv (48 rows), qc_limma.csv (48), ipa_input.csv (21), single_condition_proteins.csv (67). Each semicolon-delimited token is validated independently.
- **NaN gene names** occur in foldchange_all.csv (8), qc_limma.csv (8), ipa_input.csv (2), single_condition_proteins.csv (15) -- unannotated or ambiguous UniProt entries. Gene columns are nullable.
- **2 numeric-junk 'accession' rows in single_condition_proteins.csv** (long ';'-joined lists of bare integers, paired with NaN gene and NaN in all 4 intensity columns) are an upstream raw-data artifact passed through unmodified by the frozen pipeline. The schema accepts them via a narrowly-scoped exception requiring ALL of: accession matches `^[0-9;]+$`, gene is NaN, and all 4 intensities are NaN -- see `schema.py::_single_condition_accession_ok`.
