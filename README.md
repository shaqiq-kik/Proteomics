# Proteomics — SILAC differential expression (testosterone vs. vehicle, mouse)

A quantitative proteomics analysis of a **mouse** SILAC experiment comparing
**testosterone-treated** against **vehicle-control** samples, with **n = 2
technical replicates per condition** (no biological replication).

The analysis lives in [`proteomics_de/`](proteomics_de/README.md). It takes the
two SILAC channel sheets in `Copy of General Sheet.xlsx`, computes per-protein
log2 fold changes, runs a moderated t-test (limma, with MinProb imputation),
and produces QC figures, a PPI network, pathway enrichment, and a single
self-contained HTML report.

---

## Read this before you read any result

**Nothing in this dataset is statistically significant, and that is the
finding.**

* **0 of 1938** tested proteins pass FDR < 0.05. The smallest adjusted p-value
  in the entire experiment is **0.116**. (55 proteins reach an *uncorrected*
  p < 0.05 — fewer than the ~97 you would expect from 1938 tests by chance
  alone at that threshold.)
* **0** GO / KEGG / Reactome terms survive over-representation testing against
  the detected proteome, and **0 of 568** GSEA gene sets pass FDR < 0.05.
* This is the expected ceiling of the design, not a bug and not a processing
  error. Two **technical** replicates measure instrument reproducibility. They
  do not measure biological variability, so there is no denominator from which
  a per-protein p-value could reach significance.

Everything the pipeline emits — every fold change, every "UP"/"DOWN" call,
every network hub, every pathway — is therefore **hypothesis-generating only**.
Treat the outputs as a ranked list of candidates worth a properly replicated
follow-up, never as confirmed differential expression.

There is a second, sharper lesson buried in the enrichment layer. Running the
same gene list against g:Profiler's *default whole-genome* background returns
**196 "significant" terms** (top p = 1.9e-23). Running it against the honest
**detected-proteome** background returns **zero**. The pipeline uses the honest
background. The 196 terms are a textbook background-inflation artifact — the
false pathway story this analysis specifically declines to tell. See
[`proteomics_de/DECISIONS_LOG.md`](proteomics_de/DECISIONS_LOG.md) D6.

---

## Quickstart

Set up the environment once — see **[ENVIRONMENT.md](ENVIRONMENT.md)** for the
full instructions and the two traps that will bite you if you skip it (the venv
*must* be created with `--system-site-packages`, and a bare `python3` is a coin
flip between four interpreters, only one of which has pandas).

```sh
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -r requirements-dev.txt
```

Then the whole pipeline is one command:

```sh
.venv/bin/python run_pipeline.py --all
```

That runs all 13 stages in the one order that works, verifies each stage's
declared outputs, and prints a summary table. Useful variants:

```sh
.venv/bin/python run_pipeline.py --list             # what runs, in what order
.venv/bin/python run_pipeline.py --all --dry-run    # print the plan, run nothing
.venv/bin/python run_pipeline.py --all --skip-network   # no outbound API calls
.venv/bin/python run_pipeline.py --from viz_volcano  # that stage and everything after
.venv/bin/python run_pipeline.py --only viz_volcano,viz_heatmap
.venv/bin/python run_pipeline.py --verify-frozen    # sha256 the byte-frozen outputs
```

Three stages (`enrich_string_ppi`, `enrich_ora`, `enrich_gsea`) make **live
outbound API calls** to STRING, g:Profiler and Enrichr with the mouse gene list
— identifiers only, never raw intensities. `--skip-network` skips them and
reports them as SKIPPED; their committed artifacts stay in place so everything
downstream still runs.

A full run needs `Rscript` with `limma` and `imputeLCMD` (stage 1 shells out to
R). Expect a few minutes end to end; execution is deliberately serial.

Run the tests with:

```sh
.venv/bin/pytest
```

---

## Where things live

| Path | What |
|---|---|
| `run_pipeline.py` | The runner. The stage order, and what "correct output" means, as executable data. |
| [`proteomics_de/`](proteomics_de/README.md) | The pipeline itself — start with its README for the module map. |
| `proteomics_de/results/` | Every generated artifact: CSVs, figures, enrichment tables. |
| `proteomics_de/results/figures/` | All PNG + SVG figures, plus per-figure JSON manifests with captions and key numbers. |
| **`proteomics_de/report/report.html`** | **The deliverable.** One self-contained interactive HTML file — open it in a browser, no server, no assets. |
| `proteomics_de/results/ipa_input.csv` | 715 regulated proteins, formatted for QIAGEN IPA upload. |
| `proteomics_de/results/regulated_up.csv`, `regulated_down.csv` | The same 715 proteins split UP/DOWN (509/206) with a linear `fold_change` column, sorted by magnitude of change — for handing to a person rather than QIAGEN. |
| `proteomics_de/DECISIONS_LOG.md` | Human decisions D1–D16: what was chosen, why, and what still needs a person. |
| `proteomics_de/STATUS.md` | Generated inventory of what exists. `.venv/bin/python tools/status.py` to refresh. |
| `proteomics_de/BUILD_LOG.md` | Append-only build history, per work package. |
| `research1.md` | The original technical design doc (the Section 6 build list everything traces back to). |
| `proteomics_analysis.py` | The legacy 88-line script. Frozen baseline, deliberately never modified. |
| `Pilot Project/` | Earlier lab work. Legacy scripts; pytest is configured never to collect from here. |

---

## What still needs a human

Both open scientific questions are now resolved. Full detail is in
`proteomics_de/DECISIONS_LOG.md`.

Resolved: **D8 — the SILAC labelling step was not actually completed.** Each
sample carries its own `Intensity L`, `Intensity H` and `Ratio H/L`, and the
pipeline was flagged as possibly analysing the wrong quantity if those
channels were meaningful (D8, correlate at only r = 0.066 with the combined
`Intensity`). The professor confirmed the SILAC metabolic-labelling step
wasn't correctly executed during the experiment, so `Ratio H/L` carries no
real signal and only total intensity should be used. **This changes nothing:**
the pipeline has, from the start, used only the combined `Intensity` (H+L)
across `foldchange.py`, `limma_test.R`, `enrich/gsea.py` and `enrich/ora.py`,
and has never read or computed from the native SILAC ratio.

Resolved: **D7 — the control/treated direction was inverted.** The lab's own
earlier Pilot Project assigns the samples the other way round, and the 30
shared proteins confirmed it (correlation −0.82, sign agreement 2/30). The
flip has landed (`config/sample_sheet.tsv` now carries the corrected
orientation) and was verified as a pure sign inversion; every artifact in
`results/` reflects it — UP/DOWN is **509 UP / 206 DOWN**. p-values were
unaffected throughout — swapping group labels negates the fold change but
leaves |t| and p invariant.

Optional: `results/ipa_input.csv` is ready to upload to QIAGEN IPA, which is a
licensed manual tool no script can drive (D4).
