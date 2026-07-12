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
