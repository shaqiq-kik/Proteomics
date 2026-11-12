# Human Decisions Log

Running log of decisions that need **your** (the human's) input, or that I made on
your behalf and you may want to override. Maintained by the orchestrator across
build layers. Newest section at the bottom.

Legend: 🔵 OPEN (needs your call) · 🟢 DECIDED (my default, override if you disagree) · ⚪ FYI

---

## Pre-build assessment (2026-07-12)

**🟢 D1 — Architecture: continue forward vs. retrofit the Section-6 module tree.**
research1.md Section 6 envisions a formal package (`config/ io/ qc/ etl/ de/ …`)
with `config.yaml` + `sample_sheet.tsv` (item 1) and a **pandera** schema module
(item 3). The actual build grew as flat, individually-verified bug-fix scripts
(`foldchange.py`, `replicate_check.py`, `centering_check.py`, `limma_test.py/.R`).
The DE core they produce is complete and its outputs are frozen/verified.
*My default:* build the remaining layers (viz → enrichment → gated PCA → report)
**forward on top of the existing verified scripts** rather than refactoring them
into the formal tree first. Rationale: refactoring risks the byte-frozen outputs
for no gain toward the professor-facing report; the roadmap's *analytical* value
is in the unbuilt layers, not in restructuring working code.
*Override if:* you want the formal `config.yaml`/pandera/`io`/`etl` split as a
deliverable in its own right (e.g., for reproducibility grading).

**⚪ D2 — The headline scientific result is "0 significant at FDR<0.05".**
limma reports 0/1938 proteins significant (min adj.p = 0.305; 63 at raw p<0.05
— **the parenthesis is superseded; see D2 CORRECTED immediately below**).
This is the *expected* n=2 technical-replicate ceiling (research1.md §5), not a
bug — corr(limma log2FC, pipeline log2FC) = 1.0000. Every report figure/output
must carry the "technical replicates → hypothesis-generating only" caveat. No
action needed; flagging so it is a conscious, owned framing in the final report.

**⚪ D2 CORRECTED (2026-08-13) — those two numbers are the VANILLA run's, not
the shipped run's. D2's headline claim is unaffected.** `min adj.p = 0.305; 63
at raw p<0.05` was written pre-**D9**, when plain `eBayes()` was what shipped;
both figures match `results/qc_limma_vanilla.csv` today (0.304713459848272 and
63 rows), not the deliverable. Since D9 made `trend=TRUE, robust=TRUE` the
default, the shipped `results/qc_limma.csv` gives **min adj p =
0.116075908668444 and 55 rows at raw p<0.05**. What D2 actually asserts —
**0/1938 significant at FDR<0.05** — is true in BOTH flavours and does not
move, which is exactly why the stale pair survived unnoticed: it supported a
conclusion that was right either way. The original sentence is annotated rather
than rewritten, as with D5/D8. D9 already recorded the 0.305 → 0.116 move, and
`README.md` § "Reproducing the trend/robust limma variant" already quotes
0.116/55 correctly, so this log was the last place still carrying the old pair.

**⚪ D2 — SUPERSEDED in part by D23 (2026-11-11):** the replicate design is 2 biological replicates per group from 2 separate experiments, not technical replicates; the original text above is left as written.

**🔵 D3 — Enrichment layer (items 15–18) needs a data-source choice.**
GO/KEGG/Reactome ORA+GSEA and STRING PPI are not yet built. Two routes:
(a) **R Bioconductor** (clusterProfiler/fgsea/ReactomePA + org.Hs.eg.db) — all
currently MISSING; a heavy offline install. (b) **Web APIs** (g:Profiler / Enrichr
/ STRING REST) — light, but makes **external network calls** with your gene list.
Deferred until we reach that layer; noting now that it needs your OK on external
API calls and/or the package install. (STRING/g:Profiler send only UniProt/gene
IDs, no raw intensities.)

**🔵 D4 — QIAGEN IPA is a manual, human-run external tool.**
The pipeline produces `ipa_input.csv` (715 regulated) and header-only
`ipa_input_significant.csv` (0 significant, by design). Actually uploading to IPA
and running Core/Pathway/Upstream analyses is a step only you can perform (licensed
Salesforce app; blocks automation). No code can close this loop.

---

## Enrichment layer (2026-07-12)

**⚪ D5 — CORRECTION: the organism is MOUSE, not human.** research1.md §"Details"
and Section 4 (STRING) assume *Homo sapiens* / taxid **9606**. The actual data is
**Mus musculus / taxid 10090**: 1924/1930 gene symbols (99.7%) are title-case MGI
nomenclature (Lama1, Sptan1, Cryab…), and `H2-K1` is a mouse MHC class-I gene.
All enrichment is organism-specific, so I am building the layer for **mouse**:
STRING species=10090, g:Profiler organism=`mmusculus`, Enrichr mouse libraries /
R `org.Mm.eg.db` (NOT `org.Hs.eg.db`). Using human parameters would return
meaningless mappings. Flagging so the professor sees this was caught and corrected;
please confirm the experiment is indeed a mouse system (testosterone vs vehicle).

**🟢 D3 (DECIDED 2026-07-12) — enrichment data source = external web APIs.**
You chose web APIs over the offline R Bioconductor install. Built with STRING REST
(species 10090), g:Profiler g:GOSt (organism `mmusculus`, custom background), and
gseapy prerank against Enrichr mouse libraries. Only gene/protein IDs are sent —
never raw intensities. Item 17 (py4cytoscape live network) needs a running
Cytoscape desktop app (unavailable headless), so the network is delivered as a
static networkx figure instead. Raw API responses are cached under
`results/enrichment/raw/` for auditability.

**⚪ D6 — FINDING: enrichment is null with the correct background (and the report
should feature *why*).** With a CUSTOM background of the 2554 detected proteins,
**0 GO/KEGG/Reactome terms pass g:Profiler g:SCS significance** in either direction
(best corrected p ≈ 0.70), and **0/568 GSEA terms pass FDR<0.05** — consistent with
the 0/1938 per-protein result. Critically, the same UP query against g:Profiler's
DEFAULT whole-genome background returns **196 "significant" terms** (top p=1.9e-23,
"cytoplasm") — a textbook background-inflation artifact. The pipeline uses the
detected-proteome background and therefore does NOT fall into that trap. This is a
strong, defensible teaching point for the professor: the honest answer at n=2 is
"no enrichment survives correction," and the naive background would have manufactured
a false pathway story. (Independently reproduced by a fresh live re-query.)

---

## Build complete (2026-07-12)

All Section-6 roadmap layers are built, independently verified, and merged to `main`:
viz (11–14), enrichment (15–18), gated PCA (19), config+pandera (1,3 / D1 resolved),
and the final interactive HTML report (20). Every layer: an author worker + a
separate fresh-context verifier; protected/frozen outputs proven byte-identical
(sha256) after every layer. The report (`proteomics_de/report/report.html`) is a
single self-contained file, passed a dedicated correctness pass (all numbers trace
to source) and a dedicated visual-quality pass, and is mobile-hardened.

**🔵 Still needs YOU (nothing code can close):**
- **Confirm the organism** is mouse (D5) — I inferred it from the gene symbols; please
  verify the experiment is a mouse system so the enrichment parameters are right.
- **D4 — run QIAGEN IPA** yourself if desired: upload `results/ipa_input.csv` (715
  regulated leads). `ipa_input_significant.csv` is header-only by design (0 pass FDR).
- **Optional:** publish the report as a shareable Artifact / present it to the professor.

**🟢 Owned defaults you may override:** D1 (kept the flat verified scripts, added
config/pandera additively rather than refactoring), D3 (web APIs), the balanced
dual-track report framing, and delivering item 17 as a static network figure
(no headless Cytoscape).

---

## Test & execution hardening (2026-07-31)

**🔴 D7 — CORRECTION: the control/treated assignment was INVERTED. Now flipped.**
`proteomics_de/` had `31578, 31580 = control (Light)` and `31579, 31581 = treated
(Heavy)`, inherited from research1.md line 10. The lab's own earlier Pilot Project
says the opposite: `Pilot Project/CLEANED Silac Proteomics Soluble Factors.xlsx`
and `.../General Analysis/step1_data_cleaning.py:63-66` name the columns
**`Vehicle_Rep1_31579`, `Vehicle_Rep2_31581`, `Testosterone_Rep1_31578`,
`Testosterone_Rep2_31580`**.
Empirical confirmation on the 30 proteins shared between the pilot's cleaned table
and `results/foldchange_all.csv`: corr(pilot log2FC, pipeline log2FC) = **-0.82**,
sign agreement **2/30 (6.7%)** — mirror images.
The replicate PAIRING is unaffected (Rep1 = 31578/31579, Rep2 = 31580/31581,
exactly as already implemented); only the DIRECTION inverts.
Consequences: every log2FC negates; **UP/DOWN swap → 509 UP, 206 DOWN** (was 206
UP / 509 DOWN); the 715-row IPA total is unchanged; limma **p-values are
unchanged** (swapping group labels negates logFC and t but leaves |t| and p
invariant); ORA up/down term lists swap; the GSEA ranking reverses.
*Decided by you (2026-07-31): the pilot is right — flip the pipeline.* The fix is a
one-line change to `config/sample_sheet.tsv` once the config-driven refactor lands,
which is precisely the payoff of making that sheet load-bearing.

**🔵 D8 — OPEN, needs the professor: is the pipeline analysing the right quantity?**
Each of the four samples is a complete SILAC experiment carrying its own
`Intensity L`, `Intensity H` and `Ratio H/L` columns — verified that
`Intensity L + Intensity H == Intensity` exactly (median ratio 1.0000). The
pipeline uses the combined `Intensity`, i.e. the **sum of both isotope channels**,
and never touches the SILAC ratios.
All four samples show median log2(H/L normalized) between -1.05 and -1.53 — the
same direction in every run, so this is **not** a reciprocal label-swap design.
The pipeline's log2FC correlates with the mean native log2(H/L) at **r = 0.066**:
the two approaches answer different questions.
If Heavy is a common spike-in reference (super-SILAC), the correct per-sample
abundance is the ratio to that reference, not summed intensity, and every result
would change. If the SILAC labelling is incidental to this comparison, the current
approach is fine.
*Status: flagged, not acted on.* The pipeline continues to use `Intensity` as
built. **This needs a human answer from the lab before the results are presented.**

**🟢 D9 — eBayes default becomes `trend=TRUE, robust=TRUE`.** research1.md line 124
specifies it; the build shipped vanilla as a byte-reproducible baseline. The
experiment already proved the conclusion is unchanged (still 0/1938 significant;
min adj.p 0.305 → 0.116). Vanilla results are preserved as `qc_limma_vanilla.csv`
so both flavours stay comparable. *Decided by you (2026-07-31).*

**🟢 D10 — limma output regains `n_imputed`, `AveExpr`, `t`, `B`** (research1.md
line 169). `n_imputed` matters most: with stochastic MinProb at n=2, nothing
currently tells a consumer which values were measured and which were invented.
*Decided by you (2026-07-31).*

**🟢 D11 — the 2 junk accessions are quarantined.** Two rows of
`single_condition_proteins.csv` carry a bare-integer MaxQuant row-index list as an
accession (32,759 and 681 characters). `qc/schema.py` currently carves an
exception to let them PASS validation. They move to a quarantine file instead.
Two OTHER long accessions (69 chars, e.g. `P08752;P20612;…`) are legitimate protein
groups and stay. *Decided by you (2026-07-31).*

**⚪ D5 CLOSED — organism confirmed as mouse** by you (2026-07-31). All enrichment
parameters (STRING 10090, g:Profiler `mmusculus`, Enrichr mouse libraries) stand.

---

## Effort complete (2026-08-01)

All waves are built, tested and merged to `main`. **640 tests pass; 75 frozen
artifacts, zero drift; a full end-to-end run completes in ~12 s in a single
directory** — something this repo had never demonstrated before.

**⚪ D12 — FINDING: `replicate_check.py` had the pre-D7 mapping.** Caught during
the final pass. Because that module only *labels* two correlations and computes
no contrast, the error was silent — the numbers were right and the names were
swapped, so `qc_replicate_correlation.csv` reported the **treated** pair's
r=0.8624 (n=1723) as `control_raw_r`, and the report repeated it. Corrected:
control r=0.8407 (n=1656), treated r=0.8624 (n=1723). Now derived from the
sample sheet. This is the same failure mode as D7 — a hardcoded condition
assignment producing plausible output — and it is why the config-driven work
mattered more than the individual fixes.

**⚪ D13 — FINDING: `ipa_input_full.csv` shipped stale for one commit.** Its
p-values were the vanilla eBayes numbers while the report quoted trend/robust.
The code was right; `export/ipa_export.py` simply was not a `run_pipeline.py`
stage, so `--all` refreshed `qc_limma.csv` and never rebuilt the file quoting
it. Now stage 2 of 14. Found by the cross-file invariant suite, **not** by the
byte-freeze — a manifest will happily freeze a self-consistent but stale file.

**⚪ D14 — the byte-freeze gate was rebuilt twice, both times because it was
wrong in a way that would have trained people to ignore it.** First it froze
source files, so any refactor failed it by construction. Then it could not
survive a re-run of its own pipeline (matplotlib salts SVG ids and stamps
timestamps). It now covers 75 scientific outputs, is idempotent across
consecutive full runs, and is verified in BOTH directions — it catches a single
changed path coordinate and ignores pure regeneration noise.

---

## 🔵 STILL NEEDS YOU — nothing code can close

1. **D7 — confirm the corrected orientation** with whoever ran the experiment.
   The evidence is strong (the lab's own Pilot Project labelling, plus r = +0.82
   after correction versus -0.82 before), but it deserves a human yes.
2. **Run QIAGEN IPA.** Upload `results/ipa_input_full.csv` — it now carries
   `p_value` and `adj_p_value`, which unlock IPA's expression cutoffs.
   `ipa_input_significant.csv` is header-only by design (0 pass FDR).
3. **Biological replicates.** The only thing that lifts the n=2 ceiling. Adding
   them is now genuinely a sample-sheet edit for the design matrix, viz labels
   and the gated dispatcher (verified at n=6 and n=20); `foldchange.py`,
   `replicate_check.py` and `qc/schema.py` will fail loudly as 2-channel-SILAC
   specific and need regenerating.
4. **Decide on CI.** `.github/workflows/tests.yml` exists but only runs if you
   push to GitHub.

(D8 — is the pipeline analysing the right quantity — closed 2026-08-02, see
"SILAC quantity resolved" above.)

---

## Documentation cleanup (2026-08-01)

**⚪ D15 — stale pre-D7 examples in docstrings and comments, corrected.** D7
fixed the code; a sweep for `31578`/`31579`/`31580`/`31581` across every
`.py`/`.R` file found nine remaining spots that still illustrated the OLD
(inverted) assignment in prose — `config/design.py`'s worked example and sheet
table, `foldchange.py`'s `SINGLE_COLS`/`ONOFF_COLS` inline comments,
`qc/schema.py`'s `INTENSITY_COLUMNS` comments, `limma_test.R`/`limma_test.py`'s
handoff-name examples, `etl/build_matrix.py`'s contract note, and
`viz/style.py`'s "used to hardcode" comment (which had drifted from historical
framing into a present-tense claim that was simply false post-D7). No code, no
values, no committed output changed — comments and docstrings only, verified by
running the full suite before and after. Found independently by two exploratory
sessions; folded into one pass here rather than landing duplicate fixes.

**⚪ D16 — deleted `tests/expected/protected.sha256`; nothing consumed it.**
An exploratory session found the file — kept in D15's predecessor commit as an
"honest sources inventory" — had caused a second, live bug: `run_pipeline.py
--verify-frozen`'s banner still printed its path while the hashing it actually
delegated to had moved to `outputs.sha256` (fixed since `226dffe`, but the
banner string was a separate literal that fix missed). Grepped for every
consumer of the file, `constants.PROTECTED_FILES`, and the `protected_sha256`
fixture: zero. Deleted the file and all three dead readers rather than keep
patching references to a manifest nothing reads — git already tracks source
history better than a hand-copied hash list, and its mere presence had already
caused two separate investigations to mistake it for the gate.

---

## SILAC quantity resolved (2026-08-02)

**⚪ D8 CLOSED — confirmed by the lab: the SILAC labelling step was not
actually completed, so `Ratio H/L` carries no real signal.** The professor
replied to the report (`for_dr_walker.html`) confirming that the SILAC
metabolic-labelling step wasn't correctly executed during the wet-lab
experiment — an intended addition was omitted — so the native `Intensity L`
/ `Intensity H` / `Ratio H/L` channels do not reflect a real heavy/light
label split and should not be used. He asked whether this changes the
pipeline's interpretation.

It does not, because the premise of D8 already assumed this might be true and
the pipeline was built defensively: it has **never** read, computed from, or
exported anything derived from the native `Ratio H/L` columns. Every
consumer of intensity data — `foldchange.py` (the fold-change/regulation
call), `limma_test.R` (the moderated t-test), `enrich/gsea.py` (the GSEA
prerank ranking metric), and `enrich/ora.py` (the UP/DOWN query gene sets) —
has exclusively used the combined/summed `Intensity` (H+L) column from the
start. No result in `results/` changes, no re-analysis is needed, and no
number in the report moves as a consequence of this confirmation.

D8's original open question — "is Heavy a spike-in reference, making the
ratio the correct quantity?" — is now answered "no": there is no valid H/L
ratio to have used in the first place. The `r = 0.066` correlation noted in
the original D8 entry was, in hindsight, exactly what you'd expect from
comparing a real intensity signal against a ratio computed from a labelling
step that didn't run as designed.

---

## Supplementary regulated exports (2026-08-03)

**⚪ D17 — Three tiers of missing-data proteins never reached
`regulated_up.csv`/`regulated_down.csv`; two new additive files close the
gap.** A professor review of the 509 UP / 206 DOWN lists noticed EGF, EREG
missing from DOWN and FRZB missing from UP, despite FRZB being highly
induced in the earlier Pilot Project (log2FC +4.42 there). He also asked how
many of the ~40 proteins identified in that earlier analysis still show up
in the current lists.

Root cause: `regulated_up.csv`/`regulated_down.csv` are keyed on
`complete=True` — every one of the 4 raw replicate intensities present and
non-zero (`etl/foldchange_core.mark_complete`) — decided BEFORE limma ever
runs. That correctly excludes two tiers of otherwise-real signal:

* **360 "tier-3" partial-missingness proteins** (1-2 of 4 replicates
  zero/missing, e.g. FRZB P97401, GAS6 Q61592, LUM P51885, SLIT3 Q9WVB4)
  never get a raw log2FC or regulated call (parked in the default "NO
  CHANGE"), but ARE tested by limma after standard MinProb imputation —
  `qc_limma.csv` already carries a legitimate `limma_log2FC`/`p_value`/
  `adj_p_value`/`n_imputed` for every one of them. Reclassifying at the
  existing `LOG2_THRESHOLD` (0.585) newly calls 248 of the 360 (**152 UP,
  96 DOWN**); 0/360 pass FDR<0.05, consistent with the dataset-wide 0/1938
  finding (D2).
* **614 "tier-1/2" fully-undetected proteins** (EGF P01132, EREG Q61521
  among them) split across `single_condition_proteins.csv` (604 — absent
  from one whole raw MaxQuant sheet; EGF/EREG are `detected_in=control_only`)
  and `onoff_proteins.csv` (10 — present in both sheets structurally, every
  replicate on one side null/zero). Both are correctly excluded from limma
  (testing them would invent an entire absent condition) and so never reach
  any regulated file, complete or otherwise.

*Decided with the professor's evidence in hand, not by loosening the raw
classifier*: `ipa_input.csv`, `ipa_input_full.csv`, `regulated_up.csv` and
`regulated_down.csv` stay byte-for-byte untouched (all four ripple into
STRING/GSEA/ORA/the report, and the first is already uploaded to QIAGEN).
Two new files are added instead, both purely additive, built by the new
`export/supplementary_lists.py` module:

* `results/regulated_up_partial.csv` / `results/regulated_down_partial.csv`
  — the 152/96 tier-3 proteins, reclassified from `qc_limma.csv`'s
  `limma_log2FC` at the existing threshold, carrying `n_imputed` on every
  row so a reader can never mistake an imputed call for a complete-case one.
* `results/qualitative_changes.csv` — the 614 tier-1/2 proteins (372 UP /
  242 DOWN: UP = detected only in treated / `on_with_treatment`, i.e. gained
  after testosterone; DOWN = detected only in control / `off_with_treatment`,
  i.e. lost), with NO fold-change or p-value column — neither is
  statistically valid when one whole condition has zero data points or zero
  variance — and an explicit `note` column stating why.

New `run_pipeline.py` stage `regulated_lists_supplementary` (depends only on
`foldchange`, since all three source files are `foldchange` outputs); new
keys in `frozen_counts.json`; new regression tests in `test_golden_outputs.py`
(`test_frzb_gas6_lum_slit3_are_in_regulated_up_partial`,
`test_egf_ereg_are_in_qualitative_changes_as_down`,
`test_supplementary_files_are_disjoint_from_the_core_regulated_set`) pinning
FRZB/GAS6/LUM/SLIT3 into `regulated_up_partial.csv` and EGF/EREG into
`qualitative_changes.csv` with `direction=DOWN`, so this exact gap can never
silently regress again.

Separately, a one-off reconciliation (`tools/reconcile_pilot_panel.py`, not
a pipeline stage — it answers one specific email, and its input lives
entirely outside `proteomics_de/`) against the Pilot Project's 42-gene
panel, matching by UniProt accession first and gene symbol as fallback,
finds **33/42** pilot proteins present somewhere in the union of all five
current deliverables. Accession matching resolves the pilot's `C4A/C4B`
(human dual-gene nomenclature — the mouse carries one gene, `C4b`, at the
very accession the pilot cites, P01029) and three gene-symbol mismatches
(`LYZ`→`Lyz2`, `NAXE`→`Apoa1bp`, `CYRIB`→`Fam49b`) that a symbol-only match
would silently miss. The 9 still unmatched (AOC1, BGN, AGT, AMH, TGFBI,
CST3, CST12, CRISP1/CRISP3, BMP1) are all genuinely detected in the current
dataset — 8 `complete=True, regulated="NO CHANGE"` with real log2FC between
0.20 and 0.58, 1 tier-3-tested (AOC1, Q8JZQ5) at 0.174 — simply below the
±0.585 threshold in this larger, more complete run. Not a pipeline gap.

---

## Extended IPA uploads (2026-08-13)

**⚪ D18 — every IPA file QIAGEN has ever seen holds 715 proteins and NONE of
the 862 D17 recovered; three new upload files close that, additively.** The
uploads were stale in *scope*, not in content. `ipa_input.csv` and
`ipa_input_full.csv` carry the same 715 core proteins, and their overlap with
D17's 862 recovered proteins (248 tier-3 partial + 614 tier-1/2 qualitative) is
**exactly 0** — verified, not assumed. Frzb (P97401), Egf (P01132) and Ereg
(Q61521), the three proteins the professor named, are absent from **all six**
frozen deliverables: `ipa_input.csv`, `ipa_input_full.csv`,
`ipa_input_full.txt`, `ipa_input_significant.csv`, `regulated_up.csv`,
`regulated_down.csv`.

Root cause, and why re-running fixes nothing: `etl/foldchange_core.py:313`
(`build_ipa_frame`, line 307) selects `df[complete] & (regulated != "NO
CHANGE")` — the same `complete == True` gate D17 diagnosed, decided BEFORE
limma ever runs. D17 was deliberately additive and said so in as many words
(lines 317–320 above: "`ipa_input.csv`, `ipa_input_full.csv`,
`regulated_up.csv` and `regulated_down.csv` stay byte-for-byte untouched"), so
it left the IPA path alone. `--all` therefore reproduces the same 715 rows for
ever; re-running `export/ipa_export.py` was verified to emit byte-identical
output.

**This is D13's failure mode recurring** — a self-consistent but stale
deliverable. Every value in `ipa_input_full.csv` traces correctly to
`qc_limma.csv`, every file passes `validate_ipa`, and the byte-freeze is
perfectly content, because *a manifest will happily freeze a stale file* (D13).
D13 was stale in its **values**; D18 is stale in its **scope**. Hashing cannot
catch either.

*Decision: additive again, on D17's precedent, not by loosening the gate.* New
module `export/ipa_extended.py`; new `run_pipeline.py` stage
`ipa_export_extended` (stage 5 of 17, registered for precisely the reason D13
gives — an export that is not a stage goes stale the moment an upstream stage
re-runs); new filenames beside the old ones. The six files above stay
byte-for-byte untouched: `ipa_input.csv` is **already uploaded to QIAGEN**, and
the rest ripple into STRING/GSEA/ORA/the report. That is not a promise in prose
— `ipa_extended.assert_frozen_inputs_untouched` re-hashes all six around every
build and raises if one moved.

Three new files, all read-only over existing results — nothing is recomputed,
imputed or invented:

* `results/ipa_input_extended.csv` + `.txt` — **963 rows = 715 core + 248
  tier-3 partial**, columns `UniProt Accession Number, Gene names, log2FC,
  regulated, p_value, adj_p_value, tier`. The core block is lifted verbatim
  from `ipa_input_full.csv`, same rows in the same order, so a reviewer diffing
  the two sees 715 identical lines followed by 248 additions.
  `assert_core_rows_are_verbatim` compares the raw CSV *cells* rather than the
  parsed floats: dropping `float_precision="round_trip"` on any read in the
  module would move **87 of the 715** `log2FC` values in the last ULP, with no
  other symptom.
* `results/ipa_qualitative_up.txt` (**372** accessions) and
  `results/ipa_qualitative_down.txt` (**242**) — one bare identifier column, no
  measurements, for a SEPARATE ID-only Core Analysis.

The three partition the recovered proteome: 963 + 372 + 242 = **1577 distinct
accessions, pairwise disjoint**, with the 1577 derived from `frozen_counts.json`
rather than typed.

**The scientific caveat, which must be repeated in any methods write-up:** in
`ipa_input_extended`, the `log2FC` column carries **two different quantities**.
For `tier=core` it is the raw mean-of-log2-ratios on unimputed data; for
`tier=partial` it is limma's estimate after MinProb imputation of 1–2 of the 4
replicates. Both are legitimate IPA input; conflating them is not. The `tier`
column is mandatory on every row so that the distinction is auditable and a
reader can split the file on it.

**Why the 614 qualitative proteins ship as identifiers only:** they carry no
fold change and no p-value *deliberately* — neither quantity is valid when one
whole condition has zero data points. QIAGEN treats an ID-only list as a
first-class upload (research1.md SECTION 3), so nothing is lost structurally,
but one consequence must be stated plainly: **IPA will run over-representation
on these two lists and produce NO activation z-scores**, because z-scores are
driven by fold-change direction. Fabricating a sentinel ±1 log2FC to force
z-scores out of them was considered and **REJECTED** — it invents a magnitude
the data does not contain, and would put a made-up number in the same column as
two real ones.

Statistics are unchanged, and D18 improves none of them: **0/963 pass
FDR<0.05** (min adj p 0.116), consistent with **D2** and the dataset-wide
0/1938.

Pinned by six new tests in `tests/test_golden_outputs.py`
(`test_ipa_extended_is_the_core_plus_the_partial_tier`,
`test_ipa_extended_core_rows_are_ipa_input_full_unchanged`,
`test_ipa_extended_and_the_qualitative_lists_partition_1577_proteins`,
`test_frzb_reaches_ipa_as_an_up_regulated_partial_row`,
`test_egf_and_ereg_reach_ipa_on_the_qualitative_down_list`,
`test_the_qualitative_id_lists_carry_no_measurement_columns`), the
`n_ipa_extended` / `n_ipa_qualitative_up` / `n_ipa_qualitative_down` keys in
`frozen_counts.json`, and 4 new byte-freeze entries (80 → 84).

**🔵 Still needs YOU — D18 does not close these:**

1. **Running IPA is manual** (**D4**): licensed Salesforce app, no automation
   possible. Upload `results/ipa_input_extended.txt` as the quantitative
   dataset, then `ipa_qualitative_up.txt` and `ipa_qualitative_down.txt` as two
   separate ID-only Core Analyses. Expect pathways/functions from those two but
   no z-scores, per the paragraph above.
2. **ORA and STRING still run on the 715.** `enrich/ora.py`'s UP/DOWN queries
   (509/206) and `enrich/string_ppi.py`'s seed set both come from
   `foldchange_all.csv`'s `regulated` column — the same complete-gated set.
   Widening them to the 963 would make the enrichment layer agree with what
   QIAGEN is now being given; **D18 does not do that, and it stays open.** GSEA
   is unaffected: it already ranks all 1938 `qc_limma.csv` proteins, so the 248
   tier-3 proteins are in its ranking (the 614 never enter limma at all, and by
   D17's reasoning cannot).

---

## Dr. Walker's additions (2026-11-09)

**⚪ D20 — the regulated list for Dr. Walker: 715 core + FRZB (measured) + 5
control-only proteins with an ASSIGNED 100-fold decrease + 247 partial-tier
proteins admitted by a rule the user widened. Additive; assigned values are
labelled everywhere.** Dr. Walker reviewed the 715-protein list and decided (the
user approved each point):

1. **EGF (P01132), EREG (Q61521)**, plus three control-only proteins the user
   chose — **Efna4 (O08542), Tyro3 (P55144), Mfap4 (Q9D1H9)** — are added as
   **DOWN with an assigned fold change of 0.01** (a 100-fold decrease),
   `log2FC = log2(0.01) = -6.643856`. All five were **verified from the data,
   not assumed**: each is `detected_in=control_only` in
   `single_condition_proteins.csv`, detected (> 0) in **both** control readings
   (Intensity 31579, 31581 per `sample_sheet.tsv`, D7) and absent in **both**
   testosterone readings (31578, 31580), and each is a `DOWN` row of
   `qualitative_changes.csv`. The builder re-checks this on every run and refuses
   to assign a value to a protein that fails it (`verify_control_only`). It also
   refuses an assigned `|log2FC|` below the pipeline's regulation cutoff (0.585),
   so a factor such as 1.2 cannot produce a "regulated DOWN" row the pipeline
   itself would call NO CHANGE.
2. **FRZB (P97401)** is added with its **real measured values**: log2FC 4.40696
   (~21.2x UP), p 0.0203, adj.P 0.998, `n_imputed = 1`, copied unchanged from
   `regulated_up_partial.csv`. The 100-fold rule is **not** applied to it.
3. **The "~860 recovered proteins".** Dr. Walker wanted only those altered by at
   least 50% with "no problems".

**What the "~860" is.** It is the **862** proteins D17/D18 recovered
(`ipa_export_extended`'s stage description and D18 name the number): **248**
tier-3 partial rows (152 UP + 96 DOWN; 1–2 of 4 readings missing, limma-imputed)
plus **614** tier-1/2 qualitative rows (never identified on one side, no fold
change). That is unambiguous in the repo.

**The rule, and how it was widened.** The user first chose a *strict* rule:
`|log2FC| >= 0.585` AND seen in both replicates of both groups AND
`n_imputed == 0`. Rows surviving each filter were 862 → 856 (not already listed)
→ 247 (`|log2FC| >= 0.585`; only partial rows have a log2FC) → **0** (all four
readings seen) → **0** (`n_imputed == 0`). That emptiness is structural: every
tier-3 row is in that tier *because* a reading is missing, so it has
`n_imputed >= 1` (179 with 1, 69 with 2, FRZB included); none of the 614
qualitative rows has a fold change. I did not loosen the rule myself; I reported
it, and **the user then decided to widen it**:

> add every **partial-tier** protein with a real log2FC and `|log2FC| >= 0.585`
> that is not already on the list, **even though one replicate was imputed**.

`config.yaml: walker_additions.recovered_rule` now reads `min_abs_log2fc: 0.585`,
`require_partial_tier: true`, `require_all_replicates_seen: false`,
`require_no_imputed_values: false`. The last two are the old "no problems"
conditions kept as switches: setting both back to `true` restores the strict rule
(and the build then adds 0). Qualitative rows still never qualify under either
rule. Rows surviving each active filter: **862 → 856 → 247 → 247**.

**The 247 do have a problem, and the files say so.** Dr. Walker said "no
problems"; each of these proteins has one: a missing replicate imputed by the
pipeline (MinProb), `n_imputed` **1 (178 proteins) or 2 (69)**. They are **151
UP and 96 DOWN**. Each rests on **one biological replicate** (the design has two
*technical* replicates, D2) in which one technical replicate was imputed rather
than measured. Their `log2FC`, `p_value` and `adj_p_value` are the real limma
values, copied unchanged from `regulated_up_partial.csv` /
`regulated_down_partial.csv`, `value_source=measured`, `tier=partial` (as in the
source files). The manifest carries `n_imputed` per protein and a `reason` that
says the rule was widened by the user's decision and names the imputed replicate
count. `ipa_input_walker.csv`/`.txt` carry no `n_imputed` column (their column
set is unchanged); it lives in the manifest only.

**The statistics are unchanged and still do not support confident calls:** 0 of
the 247 (and of FRZB) pass FDR < 0.05. Minimum adj.P is 0.116 (Tuba1a); 246 of
the 247 sit at adj.P = 0.998, as FRZB does, and 4 have a raw p < 0.05. This is
the same picture as the core list (0 of the 715 pass FDR < 0.05, D2): n=2
technical replicates, hypothesis-generating only.

**Row counts:** 715 core → **716** (+FRZB) → **721** (+5 assigned) → **968**
(+247 recovered; 715 + 1 + 5 + 247 = 968, as expected).

**ASSIGNED IS NOT MEASURED — and this is the one place the repo knowingly writes
such a number into a fold-change column.** D18 considered and REJECTED a
sentinel ±1 log2FC for the qualitative lists because "it invents a magnitude the
data does not contain". D20 supersedes that **for these five proteins only**, on
the professor's explicit decision, and pays for it with labelling that cannot be
filtered away by accident:

* `value_source` column — `assigned` on the five rows, `measured` on the other
  963;
* `p_value` and `adj_p_value` are **blank** on assigned rows (a test statistic
  for an unmeasured value would be invented), in `ipa_input_walker.csv` and
  `.txt`. The manifest has no p-value columns; there, assigned rows are marked by
  `value_source=assigned`, a blank `n_imputed`, and a `reason` that says "NOT
  measured";
* the provenance sidecars (`assigned_values_are_not_measurements: true`, the
  assignment note, the rule's step counts), the manifest `reason` column, the
  README methods note, and this entry all say the 100-fold decrease is assigned;
* **"not detected" may mean below the detection limit.** The true change is
  unknown and could be much smaller or larger than 100-fold. Treat these five as
  "lost or strongly reduced after testosterone", not as a 100x measurement.

Consequence for IPA: unlike D18's ID-only lists, these rows *do* carry a
fold-change direction, so IPA can include them in activation z-scores. Five
proteins at −6.64 can dominate a z-score computed from 963 real values that are
mostly within ±1.5; say so in any methods write-up.

**Files (all NEW; nothing existing is rewritten).** `export/walker_additions.py`;
`results/ipa_input_walker.csv` + `.txt` (**968** rows; columns `UniProt
Accession Number, Gene names, log2FC, regulated, p_value, adj_p_value, tier,
value_source`; core block verbatim from `ipa_input_full.csv` in the same order,
then FRZB, the five assigned in config order, then the 247 by descending
`|log2FC|`); `results/walker_additions_manifest.csv` (253 rows: source, reason,
value_source and n_imputed of each addition); three provenance sidecars;
`config.yaml: walker_additions` (the 100, the five accessions, FRZB and the rule
— no magic numbers in code); `run_pipeline.py` stage `ipa_export_walker`;
frozen-count keys `n_walker_measured_added` (1) / `n_walker_assigned` (5) /
`n_walker_recovered` (247) / `n_walker_additions_listed` (253) / `n_ipa_walker`
(968); 3 new byte-freeze entries (85 → 88). `ipa_input.csv`, `ipa_input_full.*`,
`ipa_input_extended.*` and the other existing results stay byte-for-byte
untouched, and the builder re-hashes them around every run. The builder stages
its three files in a scratch directory, validates them there, and only then
renames them into place, so a failed validation leaves no deliverable behind.

**Two existing assertions changed, deliberately, not weakened:**
`test_freeze.py` pinned 85 frozen artifacts (now 88, +3 new files) and
`test_run_pipeline.py` pinned 17 stages and their order (now 18, with
`ipa_export_walker` after `ipa_export_extended`). Both are explicit counts whose
purpose is to force a conscious update when something is added.

**🔵 Still needs YOU:** (1) run IPA manually on `ipa_input_walker.txt` (**D4**);
(2) decide whether Dr. Walker should be told plainly that the 247 are *not* the
"no problems" set he asked for — they were included at the user's decision, and
the manifest says so per protein; (3) the same n=2 technical-replicate caveat
applies to everything on this list: 0 of the core 715, of FRZB and of the 247
pass FDR < 0.05, and the five assigned rows have no p-value at all.

**⚪ D20 — SUPERSEDED in part by D23 (2026-11-11):** the replicate design is 2 biological replicates per group from 2 separate experiments, not technical replicates; the original text above is left as written.

**⚪ D20 — SUPERSEDED in part by D24 (2026-11-11):** the walker list is now 721 proteins (the 247 recovered proteins were dropped at Dr. Walker's request); the 968-protein numbers above are history, kept as written.

---

## ORA and STRING on Dr. Walker's list (2026-11-09)

**⚪ D21 — g:Profiler ORA and STRING now also run on the 968-protein walker list, additively; 0 terms are
significant; the database versions moved since the 715 run, so a same-version 715 control was added.**
Dr. Walker wants the D20 list run through IPA, g:Profiler and STRING. `enrich/ora.py` and `enrich/string_ppi.py`
still run only on the 715 (D18 left that open); nothing there was changed.

*What was built (all NEW).* `enrich/walker_enrichment.py` (reuses `ora._post`, `_query_order`,
`_rows_from_result` and the `string_ppi` graph/metric recipe; never calls their `main()`); `config.yaml:
walker_enrichment`; `run_pipeline.py` stage `enrich_walker` (stage 16 of 19); `results/enrichment_walker/`
(`ora_up/down.csv`, `ora_meta.json`, `ora_top_terms_detail.json`, `string_edges.tsv`, `string_node_metrics.csv`,
`string_meta.json`, `string_focus_proteins.csv`, `string_comparison_vs_715.json`, `comparison_vs_715.csv`,
`raw/`, `README.md`, and `control_715_same_version/` with the same files); `results/ipa_analysis_walker/
UPLOAD_STEPS.md`; `tools/freeze.py` excludes `results/ipa_analysis_walker/` (future manual IPA exports land
there; the prefix test does not match `ipa_analysis`). `results/enrichment/*` and every existing hash are
byte-identical (the run hashes the old tree before and after and raises on any change; the 88 prior freeze
entries are unchanged, 38 new ones added: 88 to 126).

*Input and background.* Query = `regulated` of `ipa_input_walker.csv`: **UP 661 / DOWN 307**, queried separately,
one symbol per row as `enrich_common` resolves them. Background = the same measured proteome as the 715 run
(2552 rows / 2530 unique symbols; also the IPA Reference Set, D19). **Checked:** all 968 query symbols are inside
it (the 5 assigned proteins are rows of `single_condition_proteins.csv`; the partial tier is in
`foldchange_all.csv`), and all 968 accessions are in `ipa_background_measured.txt`, so the universe was **not
extended** (a missing symbol makes the code refuse, never silently extend or drop).

*Assigned values are not used.* ORA and STRING use membership and direction only. The five assigned
log2FC = -6.643856 never reach either analysis; the STRING node table leaves their `log2FC` blank; no p-value
is attached to anything. A test scrambles every number in the list and demands identical query symbols.

*Results.* **ORA: 0 terms pass g:SCS < 0.05 in either direction** (same null as D6). g:Profiler mapped 618/661
UP and 294/307 DOWN symbols (37 and 9 obsolete/mangled symbols failed to map and 6 and 4 were ambiguous and
dropped; the 715 run lost 31+5 and 7+3).
`comparison_vs_715.csv` (66 terms: significant or top 20 by rank in any of the committed 715, same-version 715
and 968 runs; rank = corrected p, then the uncorrected hypergeometric p computed from the response's own counts)
shows UP leads mostly kept (10 kept / 10 gained / 10 lost; translation initiation, ribosome, COPII, spliceosome)
and DOWN leads almost wholly replaced (4 kept / 16 / 16). **12 of the top 20 DOWN leads contain an assigned
protein**, and only four of them (the ERBB/EGFR 3-gene terms, hits exactly Egf, Ereg, Egfr) are artefacts of the
additions; the other eight contain an assigned protein but are also driven by core and recovered proteins (an
approximate leave-the-five-out recomputation moves them rank, it does not remove them), and the 96 recovered DOWN
additions move the DOWN ranks too. Nothing is significant, so ranks are tie-break aids only. The 715 run's best DOWN term ("dephosphorylation", corrected p 0.703) is diluted to 1.0 by the
larger DOWN list. **STRING: 937/968 map, 937 nodes, 10,696 edges** (control 694 nodes / 6,796 edges); the
walker network is exactly the control's 6,796 core-core edges plus 3,900 edges touching an added protein.
Egf degree 66 (68th of 937, community 0 of 299), Ereg 11, Tyro3 8, Mfap4 16, Efna4 2 (neighbours Egfr, Efnb1),
Frzb 3. That is STRING's prior knowledge of these genes, not evidence from this data.

**⚠ Version drift, and what was done about it.** The committed 715 run used **STRING 12.0** and g:Profiler
**e114_eg62_p19_27110d83**; today's API answers **STRING 12.5** and **e114_eg62_p19_9dc92ac8**. A walker-vs-committed
comparison would therefore confound the list with the databases. STRING 12.5 alone changes the 715 network from
5,963 to 6,796 edges and from 57 to 35 communities; the g:Profiler change moved no term rank in the comparison
(committed and same-version ranks are identical on every row). The 715 list was re-run under today's versions
into `control_715_same_version/`, and the comparison file carries both `*_vs_committed` and `*_vs_control`
flags. **Compare the 968 run with the control, not the committed files.** The committed 715 outputs (and the
`ppi_network` figure, the `ora_dotplot` and the report built on them) remain at STRING 12.0 / the older
g:Profiler release; refreshing them is a separate decision (the new code can pin a release via
`string.api_base`, e.g. `version-12-0.string-db.org`).

*Offline replay and caching.* Every response is cached verbatim in `raw/` with a request fingerprint
(`request_manifest.json`; the g:Profiler responses also carry the query they answered). The default mode REPLAYS
the cache and makes no network call; a missing, unfingerprinted or mismatched response is an error, never a silent
live query; `--fetch` queries only what is missing and never overwrites. The stage is therefore `network=False`,
deterministic, and run twice gave byte-identical frozen files. The same code replays the *committed 715 cache*
into byte-identical `ora_up/down.csv`, `string_edges.tsv` and `string_node_metrics.csv` (a test), which is what
licenses attributing differences to the list rather than the code. `raw/` is committed (about 10.5 MB across both
runs; the large `*_all` responses are gzipped, mtime 0).

*Tests.* `tests/test_walker_enrichment.py` (list, UP/DOWN split from the file, background coverage, no invented
numbers, old files byte-identical, equivalence on the 715, replay determinism offline, cache mismatch/miss,
comparison schema and flags, version drift, stage, freeze, README, IPA files). Existing assertions changed only
as forced counts: `test_freeze.py` 88 to 126; `test_run_pipeline.py` 18 to 19 stages.

*IPA.* `ipa_input_walker.txt` has the same layout as `ipa_input_extended.txt` plus one trailing `value_source`
column; every accession is in `ipa_background_measured.txt`, so no new reference set was needed.
`UPLOAD_STEPS.md` covers the column mapping, the reference set, the four analyses and the warning that **IPA
computes activation z-scores from fold-change direction, so the five assigned rows at -6.64 can move them**. The
sensitivity check needs no new file: `ipa_input_extended.txt` is exactly the walker list minus the 5 assigned
rows (verified), so uploading it as a second run, with the same reference set, shows which z-scores survive.

**🔵 Still needs YOU:** (1) upload to IPA by hand (D4) and run the sensitivity check; (2) decide whether to
refresh the committed 715 STRING/ORA outputs to the current database versions or leave them at 12.0;
(3) `Gene names` carry Excel-mangled values (e.g. `2018-09-09 00:00:00`) and obsolete symbols that g:Profiler cannot
map (46 of the 968 here, plus 10 ambiguous); querying by UniProt accession instead would recover them, but it would change the 715
method too and was not done; (4) all caveats of D20 and D2 apply: 0 of 968 pass FDR < 0.05.


**⚪ D21 — SUPERSEDED in part by D23 (2026-11-11):** the replicate design is 2 biological replicates per group from 2 separate experiments, not technical replicates; the original text above is left as written.

**⚪ D21 — SUPERSEDED in part by D24 (2026-11-11):** the walker list is now 721 proteins (the 247 recovered proteins were dropped at Dr. Walker's request); the 968-protein numbers above are history, kept as written.

---

## IPA on Dr. Walker's list (2026-11-10)

**⚪ D22 — IPA was run on the 968-protein walker list; the first run silently dropped the 5 assigned proteins
(blank p-values), so it was rerun with the p-value columns ignored; the 5 assigned values move many z-scores;
IPA refused the measured-proteome reference, so both runs use the default reference.**
All IPA work is manual (D4) and its exports live in `results/ipa_analysis_walker/` (outside the byte-freeze).
Nothing in the code, the frozen results or the tests changed.

*What happened.* Run 1 (`walker968_default_reference`) used the dataset `ipa_input_walker` with `p_value` and
`adj_p_value` mapped to Expr p-value / Expr False Discovery Rate, as `UPLOAD_STEPS.md` first said. IPA mapped 935 of
968 IDs but only **929** were analysis-ready (289 down, 640 up). The 963-protein dataset
(`ipa_input_extended`, the list minus the 5) also gives exactly 929, so the 5 assigned proteins (Egf, Ereg, Efna4,
Tyro3, Mfap4) were mapped but never analysed: **IPA drops every row that is blank in a mapped column**, and the
assigned rows have blank p-values by design (D20). Run 2 (`walker968_assigned_included`, dataset
`walker968_no_pvalues`) maps only `log2FC`; with no cutoffs the p-value columns are not needed, and **934**
molecules are analysis-ready (294 down, 640 up). Run 1 is therefore effectively the 963-protein run and doubles as
the "without the assigned values" comparison. `UPLOAD_STEPS.md` was corrected (Ignore the p-value columns).

*Reference set.* IPA refused the 2,552-protein measured proteome as the Reference Set ("too large of an overlap
with the reference set you have chosen"; the only fix IPA offers is a stricter cutoff, which this project does
not use). **Both runs use IPA's default reference (the whole Knowledge Base)**, as the D18 runs did, so the IPA
p-values carry the same background inflation D6 rejected for ORA. Compare rankings and z-scores, not p-values,
and do not compare them with the in-house ORA.

*What the 5 assigned proteins change (run 2 against run 1; the only input difference is those 5 rows).* 60
canonical-pathway, 156 upstream-regulator and 393 diseases-and-functions z-scores change; 1, 9 and 20 flip sign;
3, 4 and 14 cross the abs z = 2 line. Results that list an assigned protein: 81, 237 and 727, of which 6, 16 and
79 have abs z >= 2 (for example Organismal death z -6.08, Cell death of tumor cell lines -4.89, Cell
proliferation of carcinoma cell lines 4.60; regulators lipopolysaccharide 2.96 with Ereg, GABA -2.89 with Egf).
**Their z-scores rest partly on the assigned placeholder (-6.64), so they are not findings.** The androgen
entries barely moved: AR z 2.080 in both runs (Activated), metribolone 3.379 in both (Activated), testosterone
1.954 to 1.662 (no state). EGF itself as a regulator: z 1.088 to 0.825, no state.

*Files.* `run1_assigned_dropped_*` and `run2_assigned_included_*`, four IPA tab-delimited tables each (canonical
pathways, upstream regulators, diseases and functions, networks), and `README.md` (settings, the comparison table,
caveats).

**🔵 Still needs YOU:** (1) whether to report IPA numbers from run 2 at all, given that the 5 assigned values
shape 6 + 16 + 79 of the z-scored results (report run 1, the without-assigned analysis, as the defensible
version, with run 2 as a sensitivity result, or the reverse); (2) whether to rerun IPA with a log2FC cutoff so
the measured-proteome reference is accepted (this changes the analysed set and was not done); (3) all caveats of
D20, D21 and D2 apply: 0 of 968 pass FDR < 0.05, n = 2 technical replicates, one biological replicate.

**⚪ D22 — SUPERSEDED in part by D23 (2026-11-11):** the replicate design is 2 biological replicates per group from 2 separate experiments, not technical replicates; the original text above is left as written.

**⚪ D22 — SUPERSEDED in part by D24 (2026-11-11):** the walker list is now 721 proteins (the 247 recovered proteins were dropped at Dr. Walker's request); the 968-protein numbers above are history, kept as written.

---

## Replicate design corrected (2026-11-11)

**🔴 D23 — Dr. Walker corrected the replicate design: 2 biological replicates from 2 separate experiments (2026-11-11).**
*What he wrote:* "there were 2 biological replicates. Each pair of control and testosterone treated were from separate
experiments. Two replicates may not be enough statistical power, but at least more than one experiment was tested."
So control 31579/31581 and testosterone 31578/31580 are 2 **biological** replicates per group from 2 **separate
experiments**; each experiment gave one control and one treated sample. He did not say which treated sample goes with
which control; by replicate index in `config/sample_sheet.tsv` the assumed (UNCONFIRMED) pairs are 31579 with 31578 and
31581 with 31580. The investigator's statement is final.

*What it supersedes.* D2's "n=2 technical-replicate ceiling" and "technical replicates -> hypothesis-generating only"
caveat, and every later "technical replicates only / one biological replicate / no biological replication" statement
(D20, D21 and D22 included, and the same words in the config, docs, report, figure footers and sidecars). The original
log entries are left as written and carry a one-line pointer to this entry.

*What changed (words only).* `config.yaml` (design block: `n_biological_replicates_per_group` 1 -> 2,
`n_technical_replicates_per_group` set to 1 because the readings per sample were not stated, new `experiments: 2` and
`pairing_confirmed: false`; the caveat note and the design-matrix note), `config/constants.py` `CAVEAT_TEXT` (the single
source for every figure footer and every `.provenance.json` sidecar), docstrings and strings in `provenance.py`,
`config/design.py`, `viz/*`, `gated/pca_cluster.py`, `enrich/*`, `export/walker_additions.py`, `export/ipa_export.py`,
`replicate_check.py`, `run_pipeline.py`, `tools/status.py`, the report template, `build_facts.py`, both READMEs and
the walker READMEs; the figures' footers and captions were regenerated, the sidecars, the report and STATUS.md were
rebuilt, and the frozen hashes of exactly the text-bearing outputs were refrozen. A test now fails if the old
phrases return in the live files.

*What did NOT change.* Every number: the fold changes, the limma p-values and FDR, the thresholds, the regulated lists,
the IPA inputs, the enrichment tables. The limma design stays `~group`; the 715/968 lists and the 0 of 1938 result are
untouched.

*What is better.* A change or an absence seen in **both** experiments is a reproduced biological observation, not just a
repeated instrument reading. In particular the control-only "zero" proteins (EGF, EREG, EFNA4, TYRO3, MFAP4 and the other
control-only ones) are present in both control experiments and absent in both treated ones (absence can still mean below
the detection limit).

*What stays true.* n = 2 per group is low statistical power. 0 proteins pass FDR < 0.05 (min adj.p 0.116), so every
result stays exploratory / hypothesis-generating at this n.

**🔵 Still needs YOU / him:** (1) confirm which treated sample pairs with which control (31579 with 31578 and
31581 with 31580 is assumed); (2) decide whether to refit limma with an experiment (paired / blocked) term, which would
change the p-values and FDR but not the fold-change based regulated lists; (3) the number of readings per sample (technical
replicates) was not stated and is unknown; the pipeline treats each channel as one biological sample.

**🔴 D24 — Dr. Walker's answers: keep EFNA4/TYRO3/MFAP4, drop the 247, report the IPA run with the five assigned proteins included; the list is now 721 (2026-11-11).**
*What he answered (final):* (1) EFNA4, TYRO3 and MFAP4 are the genes he meant, so they stay with EGF and EREG as the 5 assigned
proteins. (2) The 247 recovered proteins: **drop them**. (3) IPA: report the run that includes the five assigned proteins, with
every result that depends on them flagged. The IPA reruns are manual and later; this entry only prepares the files.

*Decision.* The list shrinks from 968 to **721 = 715 core + FRZB (real values, partial tier) + 5 assigned control-only DOWN
proteins** (UP 510 = 509 core + FRZB; DOWN 211 = 206 core + 5 assigned). A new file holds the same list without the five
assigned rows (**716**: 715 core + FRZB; 510 UP / 206 DOWN) for the IPA sensitivity run "without the assigned values".

*Alternatives rejected.* (a) Deleting the recovered-protein rule from the code: rejected, because the answer could change
again and the widened path stays tested. The two switches are simply `true` again (the strict rule = his original "no
problems" wording, which admits zero). (b) Reusing the 968 ORA/STRING outputs and subtracting the 247: rejected, because
ORA and STRING are not additive; the analysis was re-run on the new list. (c) Building the without-assigned dataset by hand
or reusing `ipa_input_extended.txt`: rejected, because that file also carries the dropped 247; a build function plus
validation cannot drift from the 721 file.

*What changed.* `config.yaml: walker_additions.recovered_rule` `require_all_replicates_seen` and `require_no_imputed_values`
-> `true` (comment cites this answer; the rule machinery is untouched). `export/walker_additions.py`: new `drop_assigned()`
and `assert_without_assigned_is_the_list_minus_assigned()`, new outputs `results/ipa_input_walker_without_assigned.csv`/`.txt`
(+ provenance sidecars) written atomically with the other walker files and validated like them, sidecar notes now say STRICT
or WIDENED according to the switches. `run_pipeline.py` stage `ipa_export_walker` registers the two new outputs;
`frozen_counts.json`: `n_ipa_walker` 968 -> 721, `n_walker_recovered` 247 -> 0, `n_walker_additions_listed` 253 -> 6, new
`n_ipa_walker_without_assigned` 716, `n_walker_ora_up/down` 661/307 -> 510/211, `n_walker_string_mapped/nodes/edges`
937/937/10,696 -> 700/700/6,870, `n_walker_comparison_rows` 66 -> 47 (and a `_d24_note`). `walker_enrichment.py` labels
`walker968` -> `walker721` (comparison columns and JSON keys follow). The manifest is 6 rows (FRZB + the five). The recovered
pool is 862, 856 not already listed, 247 clear the cutoff, 0 pass the strict rule. Rows 1-715 of the 721 equal `ipa_input.csv`,
row 716 is FRZB with its real values, rows 717-721 are the five assigned (checked independently of the builder); all 721 and
all 716 accessions are in `ipa_background_measured.txt`. Tests: the widened-rule tests now set the switches to `false` in a
tmp config (the 247/968 path is still exercised) and the default config must yield 0 recovered; new tests cover the 716 file,
`drop_assigned`, the 721 layout and the new count key. Full default suite: 944 passed (885 before; 59 more, new and split
tests), 1 failed (`test_e2e_smoke.py::test_committed_fixture_regenerates_byte_for_byte`, which also fails at HEAD), 5
deselected (network/slow). Freeze: 126 -> 128 entries; only the walker/ORA/STRING files
that depend on the list changed hash (plus the two new files); every other hash is unchanged.

*ORA and STRING on the 721 (re-run live 2026-11-11 with the user's approval; recomputed, nothing copied from the 968 run).*
g:Profiler `e114_eg62_p19_9dc92ac8` and STRING 12.5, **both unchanged** from the 968 run, so the same-version 715 control
(fetched 2026-11-09) was kept, not refreshed. Queries UP 510 / DOWN 211, all inside the measured-proteome background (2530
symbols, not extended). g:Profiler mapped 474 and 201 (31 + 5 and 7 + 3 lost, as in the 715 control). **0 significant terms in
either direction** (best corrected p: UP 1.0, DOWN 0.897 "dephosphorylation", 0.703 in the 715 runs). STRING: 700 of 721 map
(the 21 unmapped are core protein-group accessions), 700 nodes, 6,870 edges (6,796 core-core, identical to the control, + 74
touching an added protein), 35 communities, 25 isolated nodes. Where the six sit (degree, rank of 700): Egf 51 (61st), Ereg 8
(441st), Efna4 1 (653rd, only neighbour Efnb1), Tyro3 6 (496th), Mfap4 8 (441st), Frzb 2 (624th). Seven of the top 20 DOWN
leads contain an assigned protein and they are exactly the seven that are new against the control; the approximate
leave-the-five-out ranks (uncorrected hypergeometric p, intersection reduced by the assigned hits, query size by 5) equal the
control's ranks to the rank (this is checkable because removing the five leaves the control's DOWN list). No ERBB/EGFR term is
among the DOWN leads any more (Egfr was among the dropped proteins).

*Against the 968 run.* The 968 run had UP 661 / DOWN 307, 937 nodes and 10,696 edges (3,900 touching an added protein),
66 comparison rows, 12 of the top 20 DOWN leads with an assigned protein (four ERBB/EGFR 3-gene leads that existed only
through the assigned and recovered proteins) and an all-1.0 corrected-p table. Dropping the 247 removes most of the
network growth (74 instead of 3,900 added edges), most of the DOWN lead movement (7 instead of 12) and the ERBB/EGFR
leads; the UP leads stay the control's, and nothing is significant in either version.

*Superseded but kept.* The 968-protein outputs and cache (commits `b9c75b9`, `9e345d9` hold them), the two IPA runs
(`run1_*`, `run2_*` exports) and the two figures in `results/ipa_analysis_walker/figures/` are SUPERSEDED by the 721 list.
The old files stay as the record, with a banner in `ipa_analysis_walker/README.md`; the 968-specific raw cache was removed
from `results/enrichment_walker/raw/` (it would have been stale: the request fingerprints no longer match) and the file
`results/ipa_analysis_walker/figures/make_walker_figures.py` would now draw the 721 network if re-run.

*What stays true.* n = 2 biological replicates per group from 2 separate experiments (D23); 0 of 1938 pass FDR < 0.05; the
five assigned values are assignments, not measurements; everything is exploratory.

**🔵 Still needs YOU:** (1) new IPA runs (manual): the 721 list (`ipa_input_walker.txt`, the run to report, five assigned
included) and the 716 list (`ipa_input_walker_without_assigned.txt`), then flag every z-score that depends on the assigned
values (steps in `results/ipa_analysis_walker/UPLOAD_STEPS.md`); (2) a v2 of the emailable report page (the 11 Nov page
still describes the 968 list and the first IPA runs); (3) D23's open question: which treated sample pairs with which
control (and whether to refit limma with an experiment term).
