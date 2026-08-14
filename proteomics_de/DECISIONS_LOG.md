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
