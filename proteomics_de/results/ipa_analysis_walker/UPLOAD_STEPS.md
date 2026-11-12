# IPA upload steps for Dr. Walker's 721-protein list (DECISIONS_LOG D20 / D21 / D24)

> **D24:** the list is now 721 proteins (715 core + FRZB + the 5 assigned). The 247 recovered proteins of the earlier
> 968-protein list were dropped at Dr. Walker's request, and the two earlier IPA runs in this folder are superseded.
> **New IPA runs are pending (manual):** one on `ipa_input_walker.txt` (721, the run to report, with the five assigned
> proteins included; every result that depends on them must be flagged) and one on
> `ipa_input_walker_without_assigned.txt` (716, the sensitivity run).

IPA is a licensed web app, so this is a manual job (D4). This page is the checklist. Drop IPA's
exports into this folder when you are done (it is outside the byte-freeze, like `results/ipa_analysis/`).

## Files you need (all under `proteomics_de/results/`)

| File | Rows | What it is |
|---|---|---|
| `ipa_input_walker.txt` | 721 | The dataset to analyse (tab-delimited; `ipa_input_walker.csv` is the identical CSV twin). 510 UP / 211 DOWN. |
| `ipa_input_walker_without_assigned.txt` | 716 | The sensitivity dataset (step 7): the 721 minus the 5 `value_source=assigned` rows (715 core + FRZB; 510 UP / 206 DOWN). `.csv` twin alongside. |
| `ipa_background_measured.txt` | 2,552 | The Reference Set: every protein the experiment measured. It contains all 721 accessions and all 716 (checked). |

`ipa_input_walker.txt` has the same layout as `ipa_input_extended.txt` (tab-delimited, one header
row, accession leftmost, UTF-8, no BOM) with one extra trailing column, `value_source`;
`ipa_input_walker_without_assigned.txt` has exactly the same columns. Columns:
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
     column: in the earlier 968-protein runs, with them mapped the 5 assigned proteins were never
     analysed (run 1); with them ignored they were (run 2). With no cutoffs set, the analysis does not
     use them anyway. (For the 716-row file this does not matter: it has no blank p-values.) (Correction made 2026-11-10; the first
     version of this step said to map them.)
   * **Ignore** (set to "Ignore"/"Not used"): `Gene names`, `regulated`, `tier`, `value_source`.
3. Choose **Core Analysis -> Expression Analysis** on Expr Log Ratio. Species: all. Leave IPA's
   relationship and confidence filters at default as before.
4. **Set NO expression cutoff.** The file is already filtered to |log2FC| >= 0.585, and a p-value /
   FDR cutoff would drop the 5 assigned rows (no p-value) and nearly everything else (0 of the proteins pass FDR < 0.05).
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
   placeholder, not a measurement ("not detected" may mean "below the detection limit"). Next to 716
   real values that are mostly within +/-1.5, five rows at -6.64 can move or even decide a z-score
   (EGF and EREG alone feed many growth-factor regulators and pathways). So any z-score that depends on
   them is partly an artefact of the assignment.
   **Status: pending.** The earlier runs on the 968-protein list are superseded (D24) and cannot be reused.
   Run the **same analysis, same settings, same Reference Set, on the dataset without the 5 assigned rows**:
   * Upload `ipa_input_walker_without_assigned.txt` (716 rows). It is `ipa_input_walker.txt` minus the 5
     assigned rows, in the same order with identical values (verified by a test), so no manual editing is needed.
     (Do not use `ipa_input_extended.txt` for this: it also carries the 247 proteins Dr. Walker dropped.)
   * The other option, "the existing qualitative handling" (D18), is already done: the 5 are on
     `ipa_qualitative_down.txt` (242 accessions, identifiers only, no z-scores). It tests whether they
     enrich anywhere, but it cannot say anything about direction.
   Then compare the 721 and 716 runs: for every pathway, regulator or function you intend to quote,
   check that its z-score, its sign and its "activated/inhibited" call survive. Quote only the ones that do;
   report the others as "depends on the assigned values".
8. **Export back.** Save every IPA table/figure into `proteomics_de/results/ipa_analysis_walker/`
   with names that say which run it is, for example `walker721_canonical_pathways.txt`,
   `walker721_upstream_regulators.txt`, `walker721_diseases_and_functions.txt`,
   `walker721_networks.txt`, and `sensitivity716_...` for step 7 (the earlier `run1_*` / `run2_*` files are the
   superseded 968-protein runs). Add a short README beside them
   (rows uploaded, rows IPA mapped, reference set used) in the style of `results/ipa_analysis/README.md`.

## Caveats to repeat with any IPA number from this list

* 0 of the 721 proteins pass FDR < 0.05; n = 2 biological replicates per group from 2 separate experiments (D23), which is low power (D2).
* FRZB rests on the measured samples with an imputed value (D20); it is not a "no problems" protein.
  (The 247 recovered proteins that also were not are no longer on the list, D24.)
* The 5 assigned values are not measurements. Never read a fold change of -6.64 for Egf, Ereg, Efna4,
  Tyro3 or Mfap4 as a result.
* The in-house ORA (g:Profiler) and STRING runs on the same 721 list are in `results/enrichment_walker/`.
