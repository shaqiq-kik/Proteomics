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

---

## P1 — foldchange.py made testable and cwd-independent (2026-08-01)

**Wave:** 1 · **Commit:** _this commit; branch `worktree-agent-a37cd2116177abdff`,
branched from `96d16fd`_ · **Status:** ✅ done

**Closes:** research1.md Bug 4 / Section 2 "Duplicate / cardinality detection" (lines
52 and 183 — the merge cardinality guard, which did not exist); Section 6 items 4
(`etl/merge.py`) and 10 (`export/ipa_export.py` call site); the Wave-0 follow-up
"wire the frozen scripts to `config/design.py` and `config/constants.py`", for
`foldchange.py` only.

**What changed.**

1. **Logic extracted to `etl/foldchange_core.py`** (13 stage functions: `read_sheets`,
   `merge_with_indicator`, `split_both_single`, `restore_left_order`,
   `mark_complete`, `compute_ratios`, `compute_log2fc`, `classify_regulated`,
   `detect_onoff`, `summarize`, `build_ipa_frame`, `build_single_condition_frame`,
   `build_onoff_frame`). Every one is a plain function over DataFrames. Before this,
   all of it lived inside `if __name__ == "__main__"`, so none of it could be
   imported, and the only way to exercise the Bug 1–4b fixes was to run the full
   2,315-protein pipeline and read the console.
   `foldchange.py` is now wiring: it reads the workbook, calls the functions in
   order, and writes the CSVs. Print strings, assert order, `to_csv` call order and
   arguments, and column order are unchanged — deliberately, statement for statement.
2. **cwd lock fixed.** `INPUT_FILE` / `RESULTS_DIR` were cwd-relative, so the script
   only ran from the repo root. Both now resolve from `Path(__file__)`. Added
   `--input` / `--results-dir` / `--sample-sheet`, defaulting to today's values, and
   an explicit `sys.path.insert` for `_HERE` and `_ROOT` so the bare
   `from centering_check import ...` also resolves under `-m`.
3. **Merge cardinality guard added** (`etl/merge_guard.py`). `assert_merge_safe`
   counts duplicate accessions per input sheet, asserts
   `len(both) <= min(len(df_L), len(df_H))`, and asserts
   `len(merged) <= len(df_L) + len(df_H)`. `assert_classification_partition` asserts
   `n_up + n_down + n_nochange + n_onoff == len(df)`. Every bound is derived from the
   frames in hand — no dataset number appears in the module (there is a test that
   parses the module's AST and proves it).
4. **Frozen literals moved out of source.** `n_up == 206`, `n_down == 509`,
   `n_nc + n_onoff == 1233` and `len(df_ipa) == 715` now come from
   `tests/expected/frozen_counts.json`. The checks stay **active by default**;
   `PDE_EXPECT_BASELINE=0` disables them for a future dataset.
5. **Shim call sites installed.** `qc.boundaries.check("after_load"|"after_merge"|
   "after_foldchange", df)` and `export.ipa_export.write_ipa(...)`. The IPA writer
   receives the **live in-memory frame**, never a re-read of the CSV just written.
6. **`config.design.assert_matches(INTENSITY_COLS)`** guards the hardcoded column
   list against the sample sheet. Assert-match only: this stage is not made
   design-driven, because the L/H sheet split and the left-order/Heavy-dtype restore
   are inherently 2-channel SILAC. `LOG2_THRESHOLD` now comes from
   `config/constants.py` and is re-exported for `centering_check.py`.

**Why.** The extraction is the precondition for everything else: a stage that cannot
be called cannot be tested, and a stage that cannot be tested accumulates exactly the
kind of silent arithmetic error research1.md's Bugs 1–3 already were. The cardinality
guard is the one research1.md asks for twice and the pipeline never had; it is a
no-op on today's clean sheets, which is precisely when it is worth installing,
because the failure it catches (a many-to-many row explosion on a duplicated
accession) is invisible in every downstream count.

**Files touched.** Modified `proteomics_de/foldchange.py`. Added
`proteomics_de/etl/foldchange_core.py`, `proteomics_de/etl/merge_guard.py`,
`proteomics_de/tests/test_foldchange_core.py`, `proteomics_de/tests/test_merge_guard.py`.
Appended this entry. Nothing else.

**Verification** (all run, with output).

* Full pipeline re-run from the repo root: **all 13 files `foldchange.py` writes are
  sha256-identical to the manifest** — `foldchange_all.csv`, `ipa_input.csv`,
  `single_condition_proteins.csv`, `onoff_proteins.csv`, `qc_centering.csv`,
  `foldchange_all_centered.csv`, `qc_replicate_correlation.csv`,
  `replicate_correlation.png`, `_limma_input.csv`, `_limma_output.csv`,
  `_limma_versions.txt`, `qc_limma.csv`, `ipa_input_significant.csv`. The R leg
  (Rscript 4.6.1 / limma 3.68.4 / imputeLCMD 2.1, seed 42) reproduces bit for bit.
* **Console output diffed against the pre-refactor script** (`96d16fd`'s
  `foldchange.py`, run into the same tree): `diff` is empty. Every print string, in
  order, is preserved.
* Re-run **from `/tmp`** and **as `python -m proteomics_de.foldchange`**: both
  byte-identical. Path resolution is genuinely cwd-free.
* Re-run with `--results-dir <tmpdir>`: all 10 result files byte-identical to the
  committed ones, and `results/` untouched. The CLI is real, not decorative.
* Frozen-count assertions proven **live**, not merely present: with
  `load_frozen_counts` monkeypatched to `n_up=999` the run raises
  `AssertionError: UP changed to 206, expected 999`; with `PDE_EXPECT_BASELINE=0` the
  same doctored run completes. Design guard proven live too: a sample sheet naming
  `Intensity WRONG` fails at step 0 with `ValueError: design drift in foldchange.py`,
  before any file is read.
* `pytest -q`: **106 passed** (51 before, +55 from the two new files). New tests are
  hand-built 1–4 row DataFrames with no file I/O, except `test_merge_guard.py`'s two
  positive tests, which read the real workbook (~0.4 s) to prove today's sheets pass.
* `tools/status.py --check`: **`freeze: 91 OK, 2 CHANGED, 0 MISSING`, exit 1.**
  Both CHANGED entries are explained below — neither is an output.

**Counts before → after.**

| | before (`96d16fd`) | after (P1) |
|---|---|---|
| Pipeline outputs drifted | 0 / 13 | **0 / 13** |
| UP / DOWN / NO CHANGE / ON_OFF | 206 / 509 / 1223 / 10 | 206 / 509 / 1223 / 10 |
| `ipa_input.csv` rows | 715 | 715 |
| `foldchange.py` lines under `__main__` | 195 | 0 |
| Importable fold-change functions | 0 | 19 (15 in `etl/foldchange_core.py`, 4 in `etl/merge_guard.py`) |
| Dataset literals in `foldchange.py` source | 4 (206, 509, 1233, 715) | 0 |
| Tests | 51 | 106 |
| Working directories the script runs from | 1 (repo root) | any |

**Open follow-ups.**

* **The freeze manifest covers source files, not just outputs.** All 93 tracked files
  at `421814c` are in `protected.sha256`, `foldchange.py` among them, so *any* source
  edit makes `tools/status.py --check` report drift and exit 1. It cannot report
  `93 OK` for a wave whose whole purpose is to edit sources. Two entries are CHANGED
  after this package: `proteomics_de/foldchange.py` (this package's intended edit) and
  `proteomics_de/DECISIONS_LOG.md` (**pre-existing** — commit `96d16fd` added D7–D11
  to a file the manifest had frozen at `421814c`; it was already CHANGED before I
  touched anything). The orchestrator should decide whether to split the manifest into
  an *outputs* freeze (the real gate — a 13-file subset, currently 13/13 OK) and a
  *sources* baseline that is re-cut per wave. **I did not re-baseline anything.**
* `LOG2_THRESHOLD` is `0.585`, the 3-decimal rounding of `log2(1.5) = 0.5849625…`, so
  a ratio of *exactly* 1.5 (or 1/1.5) is called NO CHANGE by 7.2e-5 of log2. The
  symmetry that Bug 2 is about is unaffected — the same rounded value bounds both
  sides — but the cutoff is a hair stricter than the "1.5-fold" it is named for. This
  is recorded, not fixed: changing the constant would move real proteins across the
  boundary and rewrite every committed output. Test
  `test_bug2_threshold_is_a_rounded_log2_of_1_point_5` pins the current behaviour.
* `qc.boundaries.check("before_ipa_export", ...)` is defined in `STAGES` but has no
  call site yet; it belongs with P7, who owns `export/ipa_export.py`.
* `centering_check.py`, `replicate_check.py` and `limma_test.py` still carry their own
  copies of the intensity-column literals and are not guarded by
  `design.assert_matches`. They were out of scope here (not my files).
* `run_limma_test` is now called with `foldchange_csv=` / `outdir=` so `--results-dir`
  flows through, but it still writes `_limma_input.csv` / `_limma_output.csv` /
  `_limma_versions.txt` into `proteomics_de/` regardless of that flag. Harmless today;
  worth tidying when `limma_test.py` is next opened.
