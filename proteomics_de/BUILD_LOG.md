# Build Log — per work package

Append-only record of **work**: what each work package changed, why, and what was
actually run to prove it. One entry per package, newest at the bottom.

## How this relates to the other three logs

| Document | Granularity | Voice | Maintained by |
|---|---|---|---|
| `DECISIONS_LOG.md` | per **decision** (D1–D6) | "here is a fork in the road; here is my default; override if you disagree" | hand-written |
| `../research1.md` § Build Log | per **bug** (Bug 6, Bug 7, …) | narrative — What / Design / Files / Result / **Lesson** / Env / Open follow-up | hand-written |
| **`BUILD_LOG.md` (this file)** | per **work package** | operational — what changed, what was run, what the counts were before and after | hand-written |
| `STATUS.md` | none — a live snapshot | generated from the filesystem by `tools/status.py`; never hand-edited | generated |

Read them as: *why we chose it* (DECISIONS_LOG) → *what we learned* (research1.md)
→ *what we did* (this file) → *what exists right now* (STATUS.md).

The overlap with `research1.md`'s Build Log is deliberate and small. That log is
scientific: it exists so a reader understands why the pipeline gives the answer it
gives. This log is operational: it exists so a reader can reconstruct the sequence
of changes, find the commit, and see the evidence each change was verified with.
A bug fix may appear in both — as a *lesson* there and as a *package* here.

## Rules

1. **Append, never rewrite.** If an entry turns out to be wrong, add a new entry
   that corrects it and mark the original `❌ reverted` / `⚠️ partial`. Do not edit
   history — a log you can rewrite is a log you cannot trust.
2. **Verification means something that was run.** Record the command and its
   actual output. "Should work", "looks correct", and "verified by inspection" are
   not verification. If nothing was run, write `none — not verified`.
3. **Never claim a file exists without looking.** In a parallel wave, describe only
   what is on disk at the time of writing, and say so.
4. **Counts before → after** exists to make silent drift visible. If a count moved
   that should not have, that is the finding, and it goes in the entry.

## Entry template

```
## <PKG-ID> — <title> (<date>)
**Wave:** N · **Commit:** <sha> · **Status:** ✅ done / ⚠️ partial / ❌ reverted
**Closes:** <gap ids / research1.md refs>
**What changed:** ...
**Why:** ...
**Files touched:** ...
**Verification:** <what was actually run, and its result — not "should work">
**Counts before → after:** ...
**Open follow-ups:** ...
```

---

## State at the start of this effort (baseline: HEAD `421814c`, 2026-07-31)

Recorded so the before/after is checkable rather than remembered. Everything below
was verified against the tree at `421814c`, not recalled.

**What existed — a complete, working, verified analysis.**

* **93 tracked files under `proteomics_de/`**, spanning the whole roadmap: the DE
  core (`foldchange.py`, `replicate_check.py`, `centering_check.py`,
  `limma_test.py` + `limma_test.R`), five viz modules, five enrichment modules, the
  gated PCA, the pandera schema layer, and the final self-contained HTML report.
* **All of research1.md's Section-6 items 1–20 are functionally present.** Eight sit
  at the path the doc proposed; twelve were delivered elsewhere — four are *inlined*
  into the flat scripts rather than existing as their own modules (item 2 `io/load.py`,
  4 `etl/merge.py`, 7 `etl/build_matrix.py`, 10 `export/ipa_export.py`), and the rest
  differ in directory or language (Python instead of R for items 13 and 16, networkx
  instead of py4cytoscape for 17). Nothing is missing. See `STATUS.md` §a for the
  per-item mapping, regenerated from disk.
* **Every artifact is on disk and internally consistent**: `foldchange_all.csv` 1948
  rows, `single_condition_proteins.csv` 606, `onoff_proteins.csv` 10, `qc_limma.csv`
  1938, `ipa_input.csv` 715, `gsea_results.csv` 568, `string_node_metrics.csv` 694 —
  all matching the contract written into `config/config.yaml`. Three files are
  header-only **by design** (`ipa_input_significant.csv`, `enrichment/ora_up.csv`,
  `enrichment/ora_down.csv`); that is the honest scientific result at n=2, per
  DECISIONS_LOG D2 and D6, not a broken step.
* Both build logs already existed: `DECISIONS_LOG.md` (D1–D6) and the per-bug
  narrative at the end of `research1.md`.

**What did not exist — everything around the analysis.**

* **Zero tests.** No `tests/` directory, no `conftest.py`, no `pytest` configuration
  anywhere in the tracked tree. The pipeline's correctness rested entirely on
  one-off manual verification recorded in prose.
* **No runner.** No `Makefile`, no `run_all.py`, no `__main__`, no shell driver. The
  order in which the scripts must run existed only in the authors' heads and in
  comments. Reproducing the analysis meant knowing to run `foldchange.py` before
  `limma_test.py` before `viz/*` before `enrich/*` before `report/build_report.py`.
* **No dependency manifest.** No `pyproject.toml`, no `requirements*.txt`, no
  environment documentation. Nothing recorded that the R side needs limma 3.68.4 and
  imputeLCMD 2.1, or which Python interpreter the scripts were verified against.
* **`config/config.yaml` was read by nothing.** Confirmed by grep: no `.py` in the
  tree loads it, and no `.py` parses `config/sample_sheet.tsv` (the two mentions in
  `gated/pca_cluster.py` are comments). The config layer was accurate *documentation*
  of a design that lived, duplicated, as literals in nine scripts. Consciously chosen
  — DECISIONS_LOG D1 — but it means the "forward path" promise in research1.md ("adding
  biological replicates is a config change, not a rewrite") was not yet true.
* **No freeze.** Outputs had been checked byte-identical during the build, but no
  manifest of expected hashes was committed, so nothing could detect drift after the
  fact.

In one line: **the science was done and the engineering around it was not.** This
effort adds the engineering without touching a byte of the science.

---

## W0 — Foundation (2026-07-31)

**Wave:** 0 · **Commit:** _pending — the orchestrator commits at the end of the wave;
baseline is HEAD `421814c`_ · **Status:** ✅ done

**Closes:** the "no tests / no manifest / no runner / config-read-by-nothing" gaps
listed in the baseline section above; research1.md Recommendation 6 ("make
reproducibility first-class: pin versions, set seeds, run pytest + pandera in CI"),
partially — pinning and the pytest harness land here, CI does not.

**What changed.** Wave 0 adds scaffolding only. **No frozen script was modified and
no committed output was regenerated** — that is the wave's defining constraint, and
the freeze manifest exists to prove it. Five strands, written in parallel by three
agents:

1. **Test scaffolding.** `pyproject.toml` (pytest configuration only — deliberately
   no `[build-system]`/`[project]`, because this repo is an analysis project and not
   a distributable package), with `testpaths` scoped to `proteomics_de/tests` so
   pytest can never wander into the legacy `Pilot Project/` scripts, a default
   `addopts` of `-m 'not network and not slow'`, and five declared markers
   (`network`, `slow`, `r`, `golden`, `e2e`). Plus `requirements-dev.txt` (pytest
   9.1.1, pytest-xdist 3.8.0), `requirements-lock.txt`, `ENVIRONMENT.md`,
   `proteomics_de/tests/{__init__.py,conftest.py}`, and the first tests.
2. **Config / design layer.** `config/design.py` reads `config/sample_sheet.tsv` and
   becomes the single source of truth for the design; `config/constants.py` collects
   the values currently duplicated as literals across nine scripts. Nothing is wired
   into the frozen scripts yet — by design, so the wiring is a separate, reviewable
   change.
3. **Accession policy.** `etl/accessions.py` states, in one place, the three
   previously-implicit and never-reconciled policies for the semicolon-separated
   `UniProt Accession Number` field (merge on the whole string; take token 0 for
   external queries; validate every token).
4. **Shims.** `export/ipa_export.py`, `provenance.py`, `qc/boundaries.py` are
   deliberately no-op / pass-through today. They exist so Wave-1 call sites can be
   installed now and filled in later without the call sites changing. Critically,
   `provenance.sidecar()` writes **no file** — an extra file on disk would break the
   "all committed outputs stay sha256-identical" acceptance test.
5. **Logging and status system** (this strand): `BUILD_LOG.md` (this file),
   `tools/status.py`, and its generated output `proteomics_de/STATUS.md`.
   `tools/status.py` reports Section-6 coverage, a full artifact inventory with data
   row counts, byte-freeze drift, test coverage (parsed with `ast`, never executed),
   and the Python/R environment. It resolves every path from `Path(__file__)`, so it
   is cwd-independent, and `--check` makes it usable as a CI gate.

**Why.** Three problems, one root cause. (a) The analysis was verified but not
*re-verifiable*: with no tests and no manifest, the only defence against a future
edit silently changing a result was human memory. (b) The design was documented but
not *authoritative*: `config.yaml` described literals it did not control, so
research1.md's forward-path promise was aspirational. (c) "What exists" was
described in prose that goes stale the moment a file is added — hence a generated
inventory instead of a written one. Wave 0 lays the floor for all three without
risking the frozen outputs; the wiring happens in later waves, on top of a freeze
that can prove nothing moved.

**Files touched.** All new; nothing modified except `.gitignore`. Observed on disk at
the time of writing (a parallel wave — sibling agents may have added more since):

* Mine: `tools/status.py`, `proteomics_de/BUILD_LOG.md`, `proteomics_de/STATUS.md`
  (generated).
* Siblings': `pyproject.toml`, `requirements-dev.txt`, `requirements-lock.txt`,
  `ENVIRONMENT.md`,
  `proteomics_de/tests/{__init__.py,conftest.py,test_design.py,test_accessions.py}`,
  `proteomics_de/tests/expected/{protected.sha256,frozen_counts.json}`,
  `proteomics_de/config/{__init__.py,constants.py,design.py}`,
  `proteomics_de/etl/{__init__.py,accessions.py}`,
  `proteomics_de/export/{__init__.py,ipa_export.py}`,
  `proteomics_de/provenance.py`, `proteomics_de/qc/boundaries.py`.

**Verification.** What I ran, and what it printed:

* `.venv/bin/python tools/status.py` → wrote `proteomics_de/STATUS.md`; reported
  `Section 6: 8 implemented, 12 elsewhere, 0 missing` and
  `freeze: 93 OK, 0 CHANGED, 0 MISSING`.
* `.venv/bin/python tools/status.py --check` → **exit 0**. No byte-frozen output has
  drifted; Wave 0 is additive, as intended.
* Ran again from `/tmp` with absolute paths → byte-identical report. Path resolution
  is genuinely cwd-independent, not accidentally correct because of where it was run.
* Row counts re-derived independently with pandas (`len(pd.read_csv(...))`, so quoted
  embedded newlines cannot inflate a `wc -l`): 1948 / 606 / 10 / 1938 / 715 / 568 /
  694. All seven match. The three header-only files parse to exactly 0 data rows,
  as expected.
* Cross-check on the freeze manifest: its 93 entries are *exactly* the 93 files
  tracked under `proteomics_de/` at `421814c` (`diff` of the two sorted lists is
  empty). The freeze has no blind spots and no phantom entries.
* R probe: `Rscript` 4.6.1, limma 3.68.4, imputeLCMD 2.1 — all present, matching what
  research1.md's Bug-7 entry recorded.
* **Not run by me:** `pytest`. The test suite belongs to a sibling agent in this same
  wave and was still being written; its verification is that agent's to record.

**Counts before → after.**

| | before (`421814c`) | after (W0) |
|---|---|---|
| Section-6 items covered | 20/20 (8 at proposed path, 12 elsewhere) | 20/20 — unchanged, no pipeline code touched |
| Committed artifacts | 93 files, 0 drifted | 93 files, **0 drifted** (sha256-verified) |
| Test files | 0 | 2 files / 34 test functions **as of this writing, mid-wave** — live count in `STATUS.md` §d |
| Dependency manifests | 0 | 3 (`pyproject.toml`, `requirements-dev.txt`, `requirements-lock.txt`) |
| Generated status report | none | `STATUS.md`, regenerable in one command |
| Scripts reading `config/` | 0 | 0 — the reader (`config/design.py`) exists, wiring is Wave 1 |

**Open follow-ups.**

* Wire the frozen scripts to `config/design.py` and `config/constants.py` (Wave 1).
  Until then the duplicated literals are still the operative source of truth, and the
  forward-path promise remains unfulfilled.
* Fill in the three shims (`export/ipa_export.py`, `provenance.py`,
  `qc/boundaries.py`) — packages P6/P7 per their own docstrings.
* Still no runner: the stage order lives in comments, not in code.
* No CI. `tools/status.py --check` is designed as the gate; nothing invokes it
  automatically yet. Note that the `golden` pytest marker is explicitly local-only
  (different BLAS/R builds differ in the last float digit), so a CI job must exclude
  it or it will fail for reasons that are not bugs.
* `STATUS.md` is a snapshot. Anyone reading it should re-run `tools/status.py` first;
  the numbers in this entry are from the moment of writing, mid-wave.
* Discrepancy worth flagging: this effort was briefed as "18/20 Section-6 items
  implemented". Direct inspection of the tree says **20/20 are functionally covered,
  0 missing** — the four "gaps" (items 2, 4, 7, 10) exist as inlined code inside
  `foldchange.py` and `limma_test.py` rather than as standalone modules. Whether
  "inlined" counts as "implemented" is a judgement call; `STATUS.md` reports the
  mapping so a reader can make it themselves.
