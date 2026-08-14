# QIAGEN IPA results — 2026-08-13

Exported from IPA (QIAGEN Ingenuity Pathway Analysis), project **Testosterone Proteomics 2026-08**.
These are the first IPA analyses to include the proteins recovered by **D17** — every earlier
upload contained 715 core proteins and **zero** of the 862 recovered ones (see D18).

## Read this before quoting any p-value

**IPA ran against its default Reference Set (the full Ingenuity Knowledge Base), not against the
proteins we actually measured.** That is the same naive-background choice this project explicitly
rejected for its in-house ORA in **D6**:

| Analysis | Background | Result |
|---|---|---|
| In-house ORA (`results/enrichment/`) | 2,530 measured proteins (custom) | **0 significant terms** |
| In-house ORA, diagnostic re-run | g:Profiler whole genome | 326 UP / 196 DOWN "significant", topped by GO:CC "cytoplasm" |
| IPA, below | full Ingenuity Knowledge Base | -log(p) up to 21.7, topped by translation/ribosome terms |

The IPA numbers here carry the same background inflation D6 documented. Translation, ribosome and
proteasome terms dominating a whole-proteome background is the classic signature of it. **Treat the
*ranking* as informative and the *absolute p-values* as optimistic**, and do not present these
alongside the ORA null result as if they disagree about biology — they disagree about background.

To make them comparable, upload the full measured proteome (the 2,552-row background used in
`enrichment/ora_meta.json`) as a dataset and set it as the IPA Reference Set, then re-run. That has
**not** been done here.

Second caveat, from D2: **nothing in this dataset is statistically significant.** 0 of 1,938
proteins pass FDR < 0.05 (minimum adjusted p = 0.116) with n = 2 technical replicates. The uploaded
963 are fold-change calls (|log2FC| > 0.585), not significance calls.

## Files

| File | Analysis | Rows | Contents |
|---|---|---|---|
| `core963_canonical_pathways.txt` | Quantitative 963 | 1,119 | Canonical pathways, -log(p), ratio, z-score, molecules |
| `core963_upstream_regulators.txt` | Quantitative 963 | 2,370 | Upstream regulators, activation z-score, p-value of overlap |
| `core963_diseases_and_functions.txt` | Quantitative 963 | 1,000 | Diseases & bio functions, activation state, z-score |
| `qualitative_up372_canonical_pathways.txt` | Qualitative UP 372 | 875 | Canonical pathways (no z-scores — see below) |
| `qualitative_up372_summary.txt` | Qualitative UP 372 | — | Full analysis summary, all sections |
| `qualitative_down242_canonical_pathways.txt` | Qualitative DOWN 242 | 750 | Canonical pathways (no z-scores — see below) |
| `core963_networks.txt` | Quantitative 963 | 25 | Networks, score, focus molecules, top diseases/functions |

Tab-delimited, UTF-8. Each file carries a QIAGEN copyright line then a header row.

### Figures (`figures/`)

| File | Analysis | Size |
|---|---|---|
| `core963_canonical_pathways.png` | Quantitative 963 | 2030 × 828 |
| `core963_diseases_and_functions_heatmap.png` | Quantitative 963 | 2293 × 800 |
| `core963_graphical_summary.png` | Quantitative 963 | 1823 × 1724 |
| `qualitative_up372_canonical_pathways.png` | Qualitative UP 372 | 2030 × 867 |
| `qualitative_down242_canonical_pathways.png` | Qualitative DOWN 242 | 2030 × 828 |

PNG, 150 dpi, exported with analysis details embedded.

**The bar charts are truncated to roughly the top 20 pathways**, so they are readable rather than
complete — the full ranking is in the matching `.txt`. IPA has no "top N" control, so the cap was
applied through its `-log(p-value)` score cutoff, set from the 20th-ranked value in each exported
table: **12.05** for core963 (exactly 20 pathways), **3.13** for DOWN (exactly 20), and **3.04** for
UP (**21** pathways — ranks 20 and 21 are tied at 3.04, so no cutoff separates them). Bar colour is
z-score; grey means no activity pattern is available, which is every bar on the two qualitative
charts.

Only the quantitative analysis has a Graphical Summary — see the z-score note below.

These are **external tool outputs**, not pipeline artifacts: they are not in `outputs.sha256` and
`run_pipeline.py --verify-frozen` does not check them. Re-running IPA will not reproduce them
byte-for-byte.

## What was uploaded

| Dataset | Uploaded | Mapped by IPA | Analysis-ready |
|---|---|---|---|
| `ipa_input_extended.txt` (963 = 715 core + 248 partial) | 963 | 930 | 929 (640 up / 289 down) |
| `ipa_qualitative_up.txt` | 372 | 332 | 332 |
| `ipa_qualitative_down.txt` | 242 | 216 | 216 |

Settings: Expression Analysis on **Expr Log Ratio**; `p_value` → Expr p-value, `adj_p_value` → Expr
False Discovery Rate; identifier type **UniProt/Swiss-Prot Accession only** (IPA also auto-detected
GenPept — removing it left the mapped count unchanged at 930, so it was contributing nothing).
**No expression cutoff was applied in IPA**, because the uploaded file is already filtered to
|log2FC| > 0.585; a second cutoff would double-filter. Biological filters were left at IPA defaults
(direct + indirect relationships, experimentally observed confidence, species All).

### 33 proteins failed to map

The 33 unmapped rows in the quantitative file are almost all semicolon-joined MaxQuant protein
groups (e.g. `Q6ZWY9;Q64525;Q64478;P10854`, log2FC 10.10) that IPA cannot parse as one accession.
**This is pre-existing and affects the original `ipa_input.csv` too**, not something D18 introduced.
Splitting to the leading accession would recover most of them; that is not done here.

At 930/963 (96.6%) the mapping rate is unremarkable, and no result quoted here depends on the
missing 33. Note the largest fold change in `ipa_input_extended` is **Tuba1a (P68369) at 10.55**, a
`tier=partial` row that maps to IPA without trouble — the unmapped histone group at 10.10 is second,
and a four-gene protein group is not a row to draw biology from either way.

### The two qualitative lists have no z-scores, by design

`ipa_qualitative_up.txt` and `ipa_qualitative_down.txt` carry identifiers only — no fold change and
no p-value, because neither is valid when one whole condition has zero data points. IPA therefore
runs over-representation on them but produces **no activation z-scores** (the `z-score` column reads
`NaN`), since z-scores are driven by fold-change direction. Fabricating a sentinel ±1 log2FC to
force z-scores was considered and rejected (D18).

IPA confirms this downstream: **neither qualitative analysis can produce a Graphical Summary**, and
says so explicitly — *"Unable to create Graphical Summary for this analysis. This occurs when the
analysis does not contain enough connectable entities with sufficiently high z-scores."* That is the
predicted consequence of an identifier-only upload, not an error.

### The DOWN analysis is flagged "failed" in IPA — read this before assuming its results are bad

IPA's status bar reports the Qualitative DOWN 242 analysis as `failed`. **The failure is confined to
one step.** On its Summary tab every stage reads `Done` — Molecules, Canonical Pathways, My
Pathways, ML Disease Pathways, My Lists, Tox Lists, Diseases & Functions, Upstream Regulators,
Regulator Effects, Graphical Summary, Networks — except **Causal Networks**, which never left
`Coming soon...`. That single incomplete step is what marks the whole analysis failed.

The exported `qualitative_down242_canonical_pathways.txt` is complete and populated (750 pathways,
with EGF appearing as a member molecule of three of them). No Causal Networks results exist for this
analysis, and none are included here.

## Results worth noting

**Quantitative 963** — top canonical pathways are Signaling by ROBO receptors (-log p 21.7,
z 2.94), GAIT Translation Signaling (21.3, z **-3.24**), Eukaryotic Translation Initiation (20.8,
z 4.32), Microautophagy (20.0), Neutrophil degranulation (19.1, z 4.57). Subject to the background
caveat above.

**The androgen signal is present and was not put there by us.** Among 2,370 upstream regulators,
**metribolone (R1881), a synthetic androgen-receptor agonist, is predicted Activated with
z = 3.379**, alongside beta-estradiol (Activated, z 2.157) and ESR2. For a testosterone-treated
experiment this is the sanity check you want, and it is derived from IPA's own knowledge base rather
than from anything this pipeline asserts.

**The D17 fix is visible end-to-end.** EGF — one of the two proteins reported missing from the
down-regulated list — now appears as a member molecule of Ephrin A Signaling, Lysosome Positioning
Signaling and Glucocorticoid Receptor Signaling in the qualitative DOWN analysis. FRZB reaches IPA
through `ipa_input_extended` as an UP/partial row (log2FC +4.41). Neither was reachable by any
earlier IPA upload.

**Glucocorticoid Receptor Signaling** appearing among the DOWN-only pathways is worth a look given
the steroid context.
