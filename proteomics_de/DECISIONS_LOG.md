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
limma reports 0/1938 proteins significant (min adj.p = 0.305; 63 at raw p<0.05).
This is the *expected* n=2 technical-replicate ceiling (research1.md §5), not a
bug — corr(limma log2FC, pipeline log2FC) = 1.0000. Every report figure/output
must carry the "technical replicates → hypothesis-generating only" caveat. No
action needed; flagging so it is a conscious, owned framing in the final report.

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

1. **D8 — is the pipeline analysing the right quantity?** THE OPEN SCIENTIFIC
   QUESTION. Each sample carries its own `Intensity L`, `Intensity H` and
   `Ratio H/L`; the pipeline uses the **sum of both isotope channels** and never
   touches the SILAC ratios. The two readouts correlate at **r = 0.066**. If
   Heavy is a spike-in reference, the correct per-sample quantity is the ratio to
   that reference and every result changes. Ask the lab before presenting.
2. **D7 — confirm the corrected orientation** with whoever ran the experiment.
   The evidence is strong (the lab's own Pilot Project labelling, plus r = +0.82
   after correction versus -0.82 before), but it deserves a human yes.
3. **Run QIAGEN IPA.** Upload `results/ipa_input_full.csv` — it now carries
   `p_value` and `adj_p_value`, which unlock IPA's expression cutoffs.
   `ipa_input_significant.csv` is header-only by design (0 pass FDR).
4. **Biological replicates.** The only thing that lifts the n=2 ceiling. Adding
   them is now genuinely a sample-sheet edit for the design matrix, viz labels
   and the gated dispatcher (verified at n=6 and n=20); `foldchange.py`,
   `replicate_check.py` and `qc/schema.py` will fail loudly as 2-channel-SILAC
   specific and need regenerating.
5. **Decide on CI.** `.github/workflows/tests.yml` exists but only runs if you
   push to GitHub.
