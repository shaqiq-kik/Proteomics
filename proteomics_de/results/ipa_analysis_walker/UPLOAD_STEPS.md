# IPA upload steps for Dr. Walker's 968-protein list (DECISIONS_LOG D20 / D21)

IPA is a licensed web app, so this is a manual job (D4). This page is the checklist. Drop IPA's
exports into this folder when you are done (it is outside the byte-freeze, like `results/ipa_analysis/`).

## Files you need (all under `proteomics_de/results/`)

| File | Rows | What it is |
|---|---|---|
| `ipa_input_walker.txt` | 968 | The dataset to analyse (tab-delimited; `ipa_input_walker.csv` is the identical CSV twin). |
| `ipa_background_measured.txt` | 2,552 | The Reference Set: every protein the experiment measured. It contains all 968 accessions (checked). |
| `ipa_input_extended.txt` | 963 | The sensitivity-check dataset (step 7). Already used for the earlier `core963` run. |

`ipa_input_walker.txt` has the same layout as `ipa_input_extended.txt` (tab-delimited, one header
row, accession leftmost, UTF-8, no BOM) with one extra trailing column, `value_source`. Columns:
`UniProt Accession Number, Gene names, log2FC, regulated, p_value, adj_p_value, tier, value_source`.

## Steps

1. In IPA, start a new project (or reuse "Testosterone Proteomics 2026-08") and choose **Upload data**.
   Pick `ipa_input_walker.txt`. Tick "my file has a header row".
2. **Map the columns.**
   * `UniProt Accession Number` -> **Identifier**, type **UniProt/Swiss-Prot Accession** only
     (the earlier runs mapped identically with GenPept removed).
   * `log2FC` -> **Expr Log Ratio** (observation 1).
   * `p_value` and `adj_p_value` -> **Ignore**. **Do not map them.** They are blank on the 5 assigned
     rows on purpose (do not fill them in), and IPA silently drops every row that is blank in a mapped
     column: with them mapped, 935 IDs mapped but only 929 were analysis-ready and the 5 assigned
     proteins were never analysed (run 1). With them ignored, 934 are analysis-ready (run 2). With
     no cutoffs set, the analysis does not use them anyway. (Correction made 2026-11-10; the first
     version of this step said to map them.)
   * **Ignore** (set to "Ignore"/"Not used"): `Gene names`, `regulated`, `tier`, `value_source`.
3. Choose **Core Analysis -> Expression Analysis** on Expr Log Ratio. Species: all. Leave IPA's
   relationship and confidence filters at default as before.
4. **Set NO expression cutoff.** The file is already filtered to |log2FC| >= 0.585, and a p-value /
   FDR cutoff would drop the 5 assigned rows (no p-value) and nearly everything else (0 of 968 pass FDR < 0.05).
5. **Reference Set (note: IPA refused this on 2026-11-10, so the runs used IPA's default reference;**
   **see `README.md`).** Upload `ipa_background_measured.txt` as its own dataset (identifier only,
   UniProt/Swiss-Prot Accession) and select it as the **Reference Set** for the Core Analysis. Do not
   leave the default (the whole Ingenuity Knowledge Base): that is the background inflation D6 rejected,
   and `results/ipa_analysis/README.md` says the earlier IPA p-values carry it. If you run both, label
   them "default reference" and "measured reference" and compare the *rankings*, not the p-values.
6. Run the analyses and export each as tab-delimited `.txt` (and PNG where useful):
   * **Canonical Pathways** (-log p, ratio, z-score, molecules)
   * **Upstream Regulators** (activation z-score, p-value of overlap; look at AR and the androgen entries)
   * **Diseases and Functions** (activation state, z-score)
   * **Networks** (score, focus molecules, top functions)
7. **Sensitivity check (do this; it is the one that matters for the 5 assigned rows).**
   **Read this first: IPA uses the direction of the fold change to compute activation z-scores.**
   The 5 assigned rows (Egf, Ereg, Efna4, Tyro3, Mfap4) carry log2FC = **-6.643856**, an *assigned*
   placeholder, not a measurement ("not detected" may mean "below the detection limit"). Next to 963
   real values that are mostly within +/-1.5, five rows at -6.64 can move or even decide a z-score
   (EGF and EREG alone feed many growth-factor regulators and pathways). So any z-score that depends on
   them is partly an artefact of the assignment.
   **Status (2026-11-10): done without a new upload.** Run 1 (p-values mapped) never analysed the 5
   assigned rows, so it already is the "without assigned" analysis; run 2 (p-values ignored) includes
   them. Comparison and numbers are in `README.md` in this folder. The instructions below are kept for
   repeating the check from scratch.
   Run the **same analysis, same Reference Set, without the 5 assigned rows**:
   * Upload `ipa_input_extended.txt` (963 rows). It is exactly `ipa_input_walker.txt` minus the 5
     assigned rows (same 963 accessions and values; verified), so no new file is needed. The row *order* differs
     for 3 rows (FRZB is row 716 in the walker file and row 718 in the extended file); that does not matter to IPA.
   * Or, equivalently, delete the 5 `value_source = assigned` lines from a copy of
     `ipa_input_walker.txt`.
   * The other option, "the existing qualitative handling" (D18), is already done: the 5 are on
     `ipa_qualitative_down.txt` (242 accessions, identifiers only, no z-scores). It tests whether they
     enrich anywhere, but it cannot say anything about direction.
   Then compare the 968 and 963 runs: for every pathway, regulator or function you intend to quote,
   check that its z-score, its sign and its "activated/inhibited" call survive. Quote only the ones that do;
   report the others as "depends on the assigned values".
8. **Export back.** Save every IPA table/figure into `proteomics_de/results/ipa_analysis_walker/`
   with names that say which run it is, for example `walker968_canonical_pathways.txt`,
   `walker968_upstream_regulators.txt`, `walker968_diseases_and_functions.txt`,
   `walker968_networks.txt`, and `sensitivity963_...` for step 7. Add a short README beside them
   (rows uploaded, rows IPA mapped, reference set used) in the style of `results/ipa_analysis/README.md`.

## Caveats to repeat with any IPA number from this list

* 0 of 968 proteins pass FDR < 0.05; n = 2 technical replicates, one biological replicate (D2).
* 247 + FRZB rest on one biological replicate with an imputed replicate (D20); they are not the
  "no problems" set Dr. Walker asked for.
* The 5 assigned values are not measurements. Never read a fold change of -6.64 for Egf, Ereg, Efna4,
  Tyro3 or Mfap4 as a result.
* The in-house ORA (g:Profiler) and STRING runs on the same 968 list are in `results/enrichment_walker/`.
