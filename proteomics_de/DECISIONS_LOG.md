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
