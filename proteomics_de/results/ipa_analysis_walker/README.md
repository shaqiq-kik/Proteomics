> **These exports and the two figures are from the earlier 968-protein list (superseded by the 721-protein list, D24; the 247 recovered proteins were dropped at Dr. Walker's request). New IPA runs on the 721 list are pending.**

# IPA runs on Dr. Walker's 968-protein list (D20 / D21)

> Replicate design corrected on 2026-11-11 (D23): the samples are 2 biological replicates from 2 separate experiments, not technical replicates; the numbers below are unchanged.

This folder holds IPA exports for `results/ipa_input_walker.txt`. Upload steps are in
`UPLOAD_STEPS.md`. This folder is outside the byte-freeze, like `results/ipa_analysis/`.

## Two runs, and why

IPA drops any row that has a blank value in a column it is told to read. The 5 assigned rows (Egf,
Ereg, Efna4, Tyro3, Mfap4; log2FC -6.643856) have blank `p_value` and `adj_p_value`. So:

| | Run 1: `assigned_dropped` | Run 2: `assigned_included` |
|---|---|---|
| IPA analysis name | `walker968_default_reference` | `walker968_assigned_included` |
| Dataset in IPA | `ipa_input_walker` | `walker968_no_pvalues` |
| Columns read | log2FC, p_value, adj_p_value | log2FC only (p_value, adj_p_value set to Ignore) |
| Mapped IDs | 935 | 935 |
| Analysis-ready molecules | **929** (289 down, 640 up) | **934** (294 down, 640 up) |
| The 5 assigned proteins analysed? | **No** (mapped but not analysis-ready) | **Yes** |

Run 1 was meant to include the 5 assigned proteins and did not; the 963-protein dataset
(`ipa_input_extended`, the list minus the 5) also gives exactly 929 analysis-ready molecules, which is
how this was found. Run 1 is therefore effectively the 963-protein run, and doubles as the
"without the assigned values" comparison. **Run 2 is the run that does what Dr. Walker asked.**
Both runs use no expression cutoff and IPA's default reference set (see below), run 2026-11-10.

Files (each run has four IPA tab-delimited exports):

* `run1_assigned_dropped_{canonical_pathways,upstream_regulators,diseases_and_functions,networks}.txt`
* `run2_assigned_included_{canonical_pathways,upstream_regulators,diseases_and_functions,networks}.txt`

## Reference set: the measured-proteome reference was refused

`UPLOAD_STEPS.md` asks for the 2,552-protein measured proteome as the Reference Set. IPA refused it:
"The set of molecules you are analyzing has too large of an overlap with the reference set you have
chosen. Analysis-ready molecules should be a subset of your dataset. Therefore, please filter the
analysis-ready molecules more stringently and try again." IPA's only offered fix is a stricter cutoff,
which this project does not use, so both runs use IPA's default reference (the whole Ingenuity
Knowledge Base), as the earlier `results/ipa_analysis/` runs did. That means IPA p-values here carry
the background inflation `results/ipa_analysis/README.md` describes (for example -log10 p of 20 or
more for ribosome terms). **Compare rankings and z-scores, not p-values, and do not compare these
p-values with the in-house ORA** (`results/enrichment_walker/`, which uses the measured proteome).

## What adding the 5 assigned proteins changed (run 2 against run 1)

The only difference in input is the 5 rows, so every difference below comes from them.

| Table | Rows in both | z-scores that changed | Sign flips | Crossed the abs z = 2 line | Rows listing an assigned protein | of those, abs z >= 2 |
|---|---|---|---|---|---|---|
| Canonical pathways | 1,155 | 60 | 1 | 3 | 81 | 6 |
| Upstream regulators | 2,373 | 156 | 9 | 4 | 237 | 16 |
| Diseases and functions | 970 | 393 | 20 | 14 | 727 | 79 |

* Most changes are small. What moved most are results that contain Egf, Ereg and Tyro3:
  "Organismal death" (z -6.08), "Cell death of tumor cell lines" (-4.89), "Apoptosis of carcinoma cell
  lines" (-4.75), "Cell movement of carcinoma cell lines" (4.72), "Cell proliferation of carcinoma cell
  lines" (4.60). Among regulators, the entries with an assigned target include lipopolysaccharide (z
  2.96, Ereg), GABA (-2.89, Egf), INSR (2.79, Ereg), IL4 (2.77, Tyro3). Among pathways, ID1 Signaling
  (-2.67, Egf), CLEAR Signaling (2.34, Egf), EPH-Ephrin signaling (2.14, Efna4).
* **Any of these z-scores depends partly on the assigned placeholder (-6.64), not on a measurement.**
  Do not quote an abs z >= 2 result that contains an assigned protein as a finding. The same protein
  lists are in each file's Molecules / Target Molecules column.
* EGF as a regulator: z 1.088 (run 1) to 0.825 (run 2), no predicted state in either run.
* Results that do not contain an assigned protein can still move slightly, because the z-score of a
  regulator or function is computed from all analysis-ready molecules in it.
* The androgen entries barely moved: AR z 2.080 in both runs (Activated), metribolone z 3.379 in both
  (Activated), testosterone z 1.954 to 1.662 (no state), EGFR z 1.695 to 1.672.

## Headline results (exploratory)

* Top canonical pathways by abs z are the same type as before (translation initiation, neutrophil
  degranulation, rRNA processing); see the exports.
* AR is predicted Activated (z 2.08) and metribolone Activated (z 3.38). Metribolone is a synthetic
  androgen; see `results/ipa_analysis/README.md` for the bias caveat.

## Caveats that apply to every number here

* 0 of 968 proteins pass FDR < 0.05; n = 2 biological replicates per group from 2 separate experiments (D23), which is low power (D2).
* 247 + FRZB rest on the measured samples with an imputed value (D20); they are not the "no
  problems" set Dr. Walker asked for.
* The 5 assigned values are not measurements. Never read a fold change of -6.64 for Egf, Ereg, Efna4,
  Tyro3 or Mfap4 as a result.
* The default reference set inflates IPA p-values (see above).
* These are single IPA runs on the 2026-11-10 knowledge base; IPA's knowledge base changes between
  releases, so rerunning later can give slightly different numbers.
