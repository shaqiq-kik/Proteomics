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

## P5 — `viz/style.py` sample maps derived from the sample sheet (2026-07-31)
**Wave:** 2 (Tier 2) · **Commit:** see below · **Status:** ✅ done
**Closes:** the widest single copy of the design literals; unblocks D7's one-line flip.

**What changed.** `viz/style.py` no longer hardcodes `SAMPLE_COLS`,
`SAMPLE_LABELS`, `SAMPLE_CONDITION`. A new `derive_sample_maps(sheet=None)` reads
`config/sample_sheet.tsv` through `config.design.read_sample_sheet()` and returns
the three maps; module scope calls it once. Labels are built as
`f"{group}_{replicate}"` from the sheet's own columns, so nothing pins a sample id
to a group. `CAVEAT_TEXT` is now `_constants.CAVEAT_TEXT` — a re-export, not a
second copy — and `style.CAVEAT_TEXT` keeps resolving for the 10+ `add_caveat()`
call sites and `viz/heatmap.py`'s direct reference. `SAMPLE_ORDER` still aliases
`SAMPLE_COLS`.

**Why.** `style` is imported by all five `viz/*` scripts plus `gated/pca_cluster.py`
and four `enrich/*` modules, so it was the broadest hardcoded copy of the design in
the tree. Deriving it once propagates to every figure. Per D7 the control/treated
assignment is inverted and will be flipped by editing the TSV; after this change
that flip reaches the figures with no code edit.

**Import mechanics.** `viz/*.py` runs as a file path, so `sys.path[0]` is `viz/`
and the repo root is not importable. `style.py` tries the plain
`from proteomics_de.config import ...` first (so pytest and `python -m` reuse the
already-imported package rather than loading a second copy — a duplicate module
would mean two `CAVEAT_TEXT` objects and a false identity re-export) and only on
`ImportError` appends `REPO_ROOT`, derived from `__file__`, never from the cwd. It
appends rather than prepends so it cannot shadow a caller's resolution.

**Files touched:** `proteomics_de/viz/style.py` (modified),
`proteomics_de/tests/test_style_samples.py` (new, 20 tests). No viz/enrich/gated
script was touched.

**Verification.**
* Derived maps equal the pre-refactor literals exactly, including dict key order
  (`qc_plots.py`/`heatmap.py` build DataFrame columns by iterating `SAMPLE_COLS`,
  so order is load-bearing for figure bytes).
* `CAVEAT_TEXT` identity holds against `config.constants`; digest pinned at
  sha256 `02edd2eb…5446`, 202 bytes, em dash asserted present.
* Forward path: a synthetic 6-sample sheet yields 6 entries with `control_3` /
  `treated_3`. Row order in the TSV does not permute the maps.
* D7 flip: a synthetic sheet with `31579/31581 = control` inverts
  `SAMPLE_CONDITION`, the column order, and the labels. Every channel's group is
  asserted to differ from today's.
* Import mechanics exercised by subprocess from three cwds (viz dir, repo root,
  `/`) plus the `from viz.style import ...` form that `enrich/network_figure.py`
  uses.
* Regenerated all four viz scripts. **All 13 PNGs byte-identical; all 4 JSON
  manifests byte-identical.** The 13 SVGs are content-identical but not
  byte-identical — see below.
* `pytest proteomics_de/tests/` → **71 passed** (51 pre-existing + 20 new).

**🔴 Finding — the SVG byte-freeze is not achievable, and never was.**
`tools/status.py --check` reports `91 OK, 2 CHANGED, 0 MISSING` (exit 1). Neither
CHANGED file is a figure:
1. `proteomics_de/DECISIONS_LOG.md` — **pre-existing**, drifted before this package
   began. The manifest was frozen at `0c41b20`; `96d16fd` then added D7–D11.
2. `proteomics_de/viz/style.py` — the file this package exists to modify.
   `protected.sha256` freezes 20 `.py` sources alongside the outputs, so *any*
   refactor of a frozen script shows as CHANGED by construction. The "93 OK"
   target is unreachable for any package that edits a frozen script.

Separately, re-running a viz script rewrites its SVG with different bytes **even
with completely unmodified code**. Proven by control experiment: `git stash`ed the
change, re-ran `volcano.py` on the pristine committed `style.py`, and
`volcano.svg` still drifted while `volcano.png` stayed byte-identical. Two
consecutive runs of the *same* pristine code also differ from each other. Causes,
both matplotlib-internal and 100% of the diff:
* `<dc:date>` — a wall-clock timestamp embedded in every SVG.
* Randomly salted element ids (`p<10hex>`, `m<10hex>`, `image<10hex>`) regenerated
  per run.
Normalizing only those two, every regenerated SVG is **identical** to its committed
counterpart. The embedded base64 rasters match exactly. So this change moved zero
pixels; the SVG rows of the freeze are simply not a valid byte-level gate.
**No figure was re-baselined** — the working tree's figures are the committed bytes.

*Remedy (needs a human decision, not done here):* set
`plt.rcParams["svg.hashsalt"]` to a fixed string in `apply_style()` and pass
`metadata={"Date": None}` in `save_fig`'s SVG call. That makes SVG output
deterministic, but it changes every committed SVG's bytes once, so it requires a
deliberate re-baseline of the 13 SVG entries in `protected.sha256`.
## P2 — Config-driven limma design matrix + file contract (2026-07-31)
**Wave:** 1 · **Commit:** (this commit) · **Status:** ✅ done
**Closes:** research1.md Section 6 item 7 (`etl/build_matrix.py`); the §1 file
contract (`intensity_matrix.tsv` + `design.tsv`); W0's follow-up "wire the frozen
scripts to `config/design.py`" — for `limma_test.{py,R}` only.

**What changed:**

* **`etl/build_matrix.py` (new).** `build(eligible_df, outdir, sheet=None)` emits
  `results/de/intensity_matrix.tsv` (`accession, gene, <sample_1..n>`) and
  `results/de/design.tsv` (`sample, group`), both TAB-separated, both shaped
  entirely by `config/sample_sheet.tsv`. Called as a *library* from
  `limma_test.py`, so the existing `foldchange.py -> run_limma_test()` chain is
  untouched; the `__main__` block is for standalone inspection only.
* **`limma_test.R`.** Hand-rolled parsing for `--in/--design/--out/--seed/--mode`,
  alongside the original positional form. With `--design`, the group vector and
  the required column list are read from the design file instead of the two
  hardcoded literals (old lines 59 and 95). Both paths now feed one shared
  `model.matrix(~ factor(group, levels = group_levels))` call, which is what makes
  the equivalence structural rather than coincidental. No `optparse` — not
  installed, and not worth a dependency for five flags.
* **`limma_test.R` error contract.** Added `bug7_abort()`. Previously every failure
  reached stderr behind R's own banner (`Error: BUG7 R ERROR: ...`,
  or `Error in parse_flags(args) :`), so the documented "stderr starts with
  BUG7 R ERROR" contract was not actually true. It is now, and is tested.
* **`limma_test.py`.** `_CTRL_COLS` / `_TRT_COLS` / `_INTENSITY_COLS` /
  `_HANDOFF_NAMES` (old lines 52-57) and the four literal intensity columns in the
  `qc_limma.csv` writer now come from `config.design`. A guard raises if the sheet
  ever stops ordering controls before treated. Writes `_limma_design.tsv` and
  passes `--design` to R. `_missing_to_blank` delegates to
  `build_matrix.intensity_series` so the CSV handoff and the TSV matrix cannot
  drift apart.

**Why:** the sample sheet drove nothing. The design matrix — the thing that
decides the sign and meaning of every fold change — was two literals in two
languages. Adding a biological replicate meant editing R.

**Deliberately NOT changed:** the eBayes mode, and the output column set. Those
are the next commit, so that the first numeric change is isolated and auditable.

**Files touched:** `limma_test.py`, `limma_test.R`, `etl/build_matrix.py` (new),
`tests/test_build_matrix.py` (new), `tests/test_limma_r.py` (new),
`tests/test_limma_contract.py` (new), `BUILD_LOG.md`. `STATUS.md` is regenerated
output of the gate command, not a hand edit.

**Verification:**

* `.venv/bin/python tools/status.py --check` →
  `freeze: 90 OK, 3 CHANGED, 0 MISSING`, exit 1. The three CHANGED are
  `limma_test.R` and `limma_test.py` — the two files this package was
  commissioned to modify, both of which the manifest freezes as *source* — plus
  `DECISIONS_LOG.md`, which was **already drifted before this work started**
  (main's own commit `96d16fd` edited it after W0 froze the manifest at
  `0c41b20`). Verified at baseline: `92 OK, 1 CHANGED`. See "Open follow-ups".
* **Every data artifact is byte-identical.** All 19 csv / 14 png / 13 svg /
  13 json / 3 tsv / 1 txt manifest entries hash to their frozen values, including
  `results/qc_limma.csv`, `_limma_output.csv`, `_limma_versions.txt` and
  `_limma_input.csv`. Zero numbers changed, which was the whole constraint.
* **Control run, before any edit:** re-running the untouched worker reproduced
  `_limma_output.csv` byte-for-byte, so the environment is stable and later
  identity results are meaningful (R 4.6.1 / limma 3.68.4 / imputeLCMD 2.1,
  matching `_limma_versions.txt`).
* **`--design` equivalence, on the real 1938-protein input:** positional and
  `--design` invocations both produce sha256
  `132039e12e802af79112f992df7d4455b8ae96f1b8531ac970b8beda05121146` — equal to
  each other and to the committed `_limma_output.csv`.
* `.venv/bin/pytest proteomics_de/tests -q` → **99 passed in 16.9s**.

**Counts before → after.**

| | before | after |
|---|---|---|
| Design literals in `style.py` | 12 (4 cols + 4 labels + 4 conditions) | 0 |
| Copies of `CAVEAT_TEXT` | 2 (`style.py`, `constants.py`) | 1 (`constants.py`) |
| Modules reached by the derivation | 0 | 10 (5 `viz/`, 4 `enrich/`, `gated/pca_cluster.py`) |
| Tests | 51 | 71 |
| Figures content-changed | — | **0** |

**Open follow-ups.**
* `style.FC_THRESHOLD` (0.585) and `style.RAW_P_THRESHOLD` (0.05) are still
  literals, duplicating `constants.LOG2_THRESHOLD` / `constants.RAW_P_THRESHOLD`.
  Left alone deliberately — out of this package's scope and `foldchange.py` was
  being edited concurrently.
* `protected.sha256` needs a decision on its two weak rows: the 13 SVGs (not
  byte-reproducible, see above) and the 20 `.py` sources (guaranteed to drift as
  the refactor proceeds). Suggest splitting it into an *outputs* freeze that gates
  CI and a *sources* inventory that does not.
* `DECISIONS_LOG.md`'s manifest row is stale and will keep failing the gate until
  someone re-baselines that single entry.
## P4 — one-command pipeline runner + the missing READMEs (2026-07-31)

**Wave:** 1 · **Commit:** see this entry's commit · **Status:** ✅ done
**Closes:** the Wave-0 open follow-up *"Still no runner: the stage order lives in
comments, not in code."* Also closes the documentation gap — the repo had no
`README.md` at any level.

**What changed.** Four new files, no existing pipeline file touched.

* **`run_pipeline.py`** (repo root) — the runner. `STAGES` is a declarative table
  of 13 stages: id, description, script path, declared outputs *with their
  expected shapes*, `network`/`slow`/`requires_r` flags, and dependencies. CLI:
  `--all`, `--from <stage>`, `--only a,b,c`, `--list`, `--dry-run`,
  `--skip-network`, `--verify-frozen` (+ `--allow-drift`). Per-stage timing, a
  summary table, nonzero exit on any FAIL.
* **`README.md`** (repo root) — what the experiment is, the one-command
  quickstart, where results live, the n=2 caveat stated plainly, and the two
  open human questions (D7, D8). Points at `ENVIRONMENT.md` for setup.
* **`proteomics_de/README.md`** — the module map: every script's inputs and
  outputs, the dependency graph, the run order, the support modules, how to
  reproduce the trend/robust limma variant, and an explicit section on the three
  by-design header-only outputs.
* **`proteomics_de/tests/test_run_pipeline.py`** — 40 tests, all offline, none
  of which executes a pipeline stage.

**Why.** Three design points are worth stating, because each was a fork.

1. **Emptiness is not failure.** Three outputs are header-only *on purpose*
   (`ipa_input_significant.csv` per D2; `ora_up.csv` / `ora_down.csv` per D6),
   and a naive "non-empty output = success" check false-positives on all three —
   it would report a failure on every single run and train everyone to ignore
   the runner. So every output declares its own expected shape and each stage
   resolves to `PASS` / `PASS_EMPTY_EXPECTED` / `FAIL`. A table contracted to 0
   rows that has 0 rows is `PASS_EMPTY_EXPECTED`, printed distinctly with its
   reason. The inverse is enforced too: 0 → non-zero is a **FAIL**, because that
   is a scientific event a human must sign off on, not a silent improvement.
2. **Serial on purpose.** Stages 3–6 are mutually independent and would
   parallelise cleanly, saving ~20 s on a multi-minute run. Execution is
   nevertheless strictly serial: a linear log where a traceback belongs to
   obviously one stage, and where stage N demonstrably wrote before stage N+1
   read, is worth more than 20 s on a pipeline a human runs a handful of times.
   `--jobs` is documented in the module docstring as a deliberate non-goal.
3. **Row counts are read, never hardcoded.** All seven headline counts resolve
   from `tests/expected/frozen_counts.json` by key. A test
   (`test_row_counts_come_from_frozen_counts_not_hardcoded_literals`) fails the
   build if anyone inlines a non-zero literal. Zero is the one permitted literal
   — the ORA term counts have no key in `frozen_counts.json`, and that 0 is the
   D6 assertion itself.

`--verify-frozen` **imports `tools/status.py` and calls its `freeze_check()`**
rather than reimplementing the manifest parsing and hashing, so the runner and
`STATUS.md` cannot drift into disagreeing. (`tools/` is not a package, so it is
loaded via `importlib.util.spec_from_file_location`.) Verdict semantics match
`tools/status.py --check` exactly: CHANGED *or* MISSING is a failure, downgraded
to a warning by `--allow-drift`.

Two smaller decisions worth recording. A **FAILED** dependency marks downstream
stages `BLOCKED` and refuses to run them (their inputs may be wrong); a
deliberately **SKIPPED** one (`--skip-network`) does *not* block downstream,
because the committed artifacts are still on disk and still valid. And stages are
launched with `sys.executable` and `cwd=<repo root>` — never bare `python3`
(four exist on PATH, one has pandas) and never the caller's cwd, because
`foldchange.py` is still cwd-locked to the repo root.

**Files touched.** All new except this log. Observed on disk at the time of
writing (parallel wave — siblings are concurrently editing `foldchange.py`,
`limma_test.py/.R`, `viz/style.py` and `enrich/*.py` in their own worktrees;
nothing of theirs was touched here):

* `run_pipeline.py`, `README.md`, `proteomics_de/README.md`,
  `proteomics_de/tests/test_run_pipeline.py`, and this entry.

**Verification.** What I ran, and what it printed:

* `.venv/bin/pytest` (whole default suite) → **91 passed in 0.35s**. My file
  alone: **40 passed**. No test spawns a subprocess, makes a network call, or
  writes into `results/`.
* `.venv/bin/python run_pipeline.py --list` → all 13 stages, numbered, with 34
  declared outputs and exactly **3** marked `0 rows EXPECTED`.
* `.venv/bin/python run_pipeline.py --from viz_volcano --dry-run --skip-network`
  → `10 stage(s) selected`, the three network stages shown as SKIPPED, and
  `No subprocess was spawned`. A test snapshots every mtime+size under
  `results/` around a dry run and asserts nothing moved.
* Contract check against the **committed** artifacts — `check_output()` run over
  all 34 declared outputs without executing any stage: **0 FAIL**, and exactly
  the three expected `PASS_EMPTY_EXPECTED`
  (`ipa_input_significant.csv`, `ora_up.csv`, `ora_down.csv`). The declared
  expectations therefore match reality as built, not as remembered.
* `.venv/bin/python run_pipeline.py --verify-frozen` → `92 OK · 1 CHANGED ·
  0 MISSING (of 93 frozen files)`, **exit 1**. With `--allow-drift`, exit 0.
  `.venv/bin/python tools/status.py --check` independently printed `freeze: 92
  OK, 1 CHANGED, 0 MISSING` and exit 1 — the two agree exactly, which is the
  point of sharing `freeze_check()`.
* cwd-independence: `--list` run from an unrelated directory with an absolute
  path produces identical output; a test asserts it under `monkeypatch.chdir`.
* **Not run by me:** any live pipeline stage. Stages 7/9/11 make live outbound
  calls to STRING / g:Profiler / Enrichr and the user has not approved a live
  run, so the runner's logic was validated entirely through `--list`,
  `--dry-run`, and direct calls into the verification layer.

**Counts before → after.**

| | before | after (P4) |
|---|---|---|
| Documented stage order | prose + comments only | `run_pipeline.py:STAGES`, executable, 13 stages |
| Commands to run the pipeline | 13, by hand, order written down nowhere | **1** |
| `README.md` files | 0 | 2 (repo root, `proteomics_de/`) |
| Test files / collected tests | 2 / 51 | 3 / 91 (mine: +40) |
| Pipeline files modified | — | **0** |
| Frozen outputs drifted by this work | — | **0** (the 1 CHANGED predates it; see below) |

**Open follow-ups.**

* ⚠️ **Pre-existing freeze drift, not caused by this work.**
  `proteomics_de/DECISIONS_LOG.md` reports **CHANGED** against
  `tests/expected/protected.sha256`. The manifest was generated at `421814c`;
  commit `96d16fd` then added D7–D11 to `DECISIONS_LOG.md` without regenerating
  it. So `--verify-frozen` and `tools/status.py --check` both exit 1 on a clean
  tree today. Either regenerate the manifest entry or drop hand-written docs
  from the freeze — a doc that is *expected* to grow does not belong in a
  byte-freeze whose whole job is to make change alarming.
* The runner cannot verify a stage that no-ops while exiting 0. It records a
  staleness NOTE when a declared output's mtime predates the stage that
  supposedly wrote it, but deliberately does not fail on it. A real
  content-dependency model would be needed to do better.
* `--from` is the manual answer to "resume after a crash". Automatic staleness
  resolution is not implemented and is not obviously worth it.
* No CI still invokes anything. `run_pipeline.py --verify-frozen` and
  `tools/status.py --check` are both usable as gates; nothing calls them.
* When D7's flip lands, the UP/DOWN counts swap (206/509 → 509/206). Nothing in
  the runner's contracts asserts direction — `foldchange_all.csv` is still 1948
  rows and `ipa_input.csv` still 715 — so the flip will not trip a false FAIL.
  Both READMEs describe the current state and defer to D7 rather than asserting
  the present orientation as settled.
| Test files / tests | 2 / 51 | 5 / **99** |
| Frozen data artifacts drifted | 0 | **0** |
| Frozen *source* files drifted | 0 | 2 (this package's mandate) |
| Hardcoded design literals in `limma_test.py` | 4 lists + 4 inline column names | 0 |
| Hardcoded design literals in `limma_test.R` | 2 (`required_cols`, `group`) | 0 when `--design` is passed |
| Scripts reading `config/` | 0 | 1 (`limma_test.py`, and `etl/build_matrix.py`) |
| Section-6 item 7 | inlined in `limma_test.py` | implemented at the proposed path |

**Findings worth recording.**

* **The reference level is decided by design-file row order.** `group_levels` is
  `unique(group)`, so whichever group appears first becomes the denominator.
  Listing treated first silently inverts every logFC in the study. This is
  precisely what `config/design.py`'s canonical sort prevents, and it is now
  pinned by `test_reference_level_follows_design_row_order`. Corollary, also
  pinned: swapping the labels *and* the row order cancels out and reproduces the
  committed numbers exactly.
* **Nothing checked the contrast direction before now.** A sign error would have
  been invisible to the entire suite and would have inverted every conclusion in
  the report. `test_contrast_direction_treated_minus_control` asserts it against
  a planted +3.0 log2 ground truth rather than against the committed numbers.
* **MinProb's stochasticity is confined to imputation.** Seeds 42 and 7 differ on
  a matrix *with* NAs and are identical on the same matrix *without* them, so
  `set.seed` is load-bearing exactly where the code claims.

**Open follow-ups.**

* **The gate cannot reach `0 CHANGED` while the freeze manifest covers source
  files.** `protected.sha256` freezes 21 `.py`/`.R` files, so *any* sanctioned
  code change drifts it. `tests/expected/protected.sha256` was outside this
  package's file scope and was deliberately **not** re-baselined. A human needs to
  decide whether the manifest should split "outputs" (never change without
  review) from "source" (changes with every commit), or whether re-baselining the
  source rows is part of every commit's ritual. The same decision covers
  `DECISIONS_LOG.md`, which `96d16fd` already drifted.
* **Contract wart in `design.tsv`.** research1.md's R sketch assumes
  `design_df$sample` names the matrix columns. It does not: `design.tsv` carries
  bare sample ids (`31578`) while `intensity_matrix.tsv`'s columns are channels
  (`Intensity 31578`), so the two contract files join only through the sample
  sheet. `limma_test.py` therefore hands R a separate `_limma_design.tsv` whose
  `sample` column holds the handoff names (`ctrl_31578`). Fixing this properly
  means changing `design.write_design_tsv`, which is W0's file, not this package's.
* Nothing reads `results/de/*.tsv` yet — the R worker is still fed by the older
  `_limma_input.csv`. Migrating the handoff onto the contract files is a separate,
  numerically-neutral step.
* `tools/status.py`'s Section-6 table still maps item 7 to `limma_test.py`; it is
  a static table inside `tools/status.py` (out of scope here) and now understates
  what is on disk.
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

---

## CHAIN-B — enrichment layer re-run after the D7 correction (2026-08-01)

Files: `enrich/enrich_common.py`, `enrich/ora.py`, `enrich/gsea.py`, `enrich/upset.py`,
new `tests/test_enrich_common.py`. Outputs regenerated live against g:Profiler and
Enrichr (user-approved; gene/protein identifiers only, per D3). Nothing committed.

**D6 re-demonstrated, not assumed. It survives the swap intact.**

| | UP | DOWN |
|---|---|---|
| query size | 509 (was 206) | 206 (was 509) |
| terms passing g:SCS < 0.05, custom background | **0** | **0** |
| best corrected p | 1.000 — GO:MF mannosyltransferase activity | 0.703 — GO:BP dephosphorylation |
| terms passing against g:Profiler's DEFAULT background | 326 (GO:CC cytoplasm, p=9.0e-70) | 196 (GO:CC cytoplasm, p=1.9e-23) |

GSEA: **0 / 568** terms pass FDR<0.05 (min FDR q = 0.1157; 15 at q<0.25, 37 at nominal
p<0.05). Term count and library sizes unchanged (6034 GO-BP + 303 KEGG).

**The swap is confirmed exactly.** The new UP query list is *identical, element for
element and in the same order*, to the list cached as DOWN at `83b22a1`, and vice
versa (`raw/gprofiler_*.json` `meta.query_metadata.queries.query_1`); the new and old
UP sets share zero genes. Three independent corroborations: the two ORA best-terms
traded places at unchanged p-values (0.702980 and 1.000000 to 6 dp); the GSEA score
range mirrored exactly, `[-3.8034621343, +2.6912806382]` → `[-2.6912806382,
+3.8034621343]`, with all 568 NES signs flipped and corr(NES_old, NES_new) = −0.9998;
and the pre-D7 "196 significant terms, cytoplasm p=1.9e-23" default-background figure
reappears unchanged — now attached to DOWN. Swapping two label sets manufactured no
enrichment, exactly as D6 predicts.

**Work done**

* **De-hardcoded the frozen literals.** `enrich_common.load_background_and_queries`
  read `== 206` / `== 509` / `== 2554` inline; the first two were wrong after D7. All
  three now come from `tests/expected/frozen_counts.json` and the assertions stay
  active (and now name the key they came from in the failure message). Docstrings that
  stated the old counts are corrected in `enrich_common.py`, `ora.py` and `gsea.py`.
  `test_no_dataset_counts_are_hardcoded_in_enrich_common` walks the module's **AST**
  (not its text — comments and docstrings never reach the AST) so the literals cannot
  creep back in.
* **Routed accession handling through the shared policy.** The inline
  `str(acc).split(";")[0].strip()` is now `etl.accessions.first_token`. Equivalence is
  asserted, not assumed: over 12 edge cases (NaN, empty, leading/trailing separators,
  non-str) and over all 2554 committed accessions.
* **Every outbound call is now cached** — offline replay is buildable.
  `raw/gprofiler_{up,down}_all.json` (the ranking call), `raw/enrichr_libraries.json.gz`
  (both Enrichr libraries, 687 KB gzipped). Also cached the gene-level evidence call as
  `raw/gprofiler_{up,down}_all_evidence.json.gz` — `ora.py` had documented it as
  deliberately *not* persisted because it is ~8–20 MB, which is precisely what made
  replay impossible; gzip takes it to 650 KB–1.1 MB.
* `upset.py` unchanged in behaviour: `ora.py` still writes the top-5 **sub-threshold**
  terms per direction, which is the only reason `upset.py` does not hit its
  `RuntimeError`. Preserved deliberately.

**🔴 Bug found and fixed in `ora.py`: every ORA gene label was misaligned.**

Not caused by D7 — it was there from the original build. g:Profiler indexes each term's
`intersections` boolean vector by the genes it successfully **mapped**, not by the
genes submitted: 509 UP symbols go out, 31 fail to map, the vectors come back length
473. `_query_order()` returned the *submitted* list, and `zip()` silently truncated —
so labels stayed correct only up to the first unmapped symbol and were wrong after it.
GO:MF "mannosyltransferase activity" was credited to `Snrpd1`; the gene g:Profiler
actually matched is `Tmtc3`.

Fixed by deriving the order from `meta.genes_metadata.query.<q>.ensgs`, inverted back
to symbols via the sibling `mapping`, plus a hard length assertion so a mismatch can
never again be absorbed by `zip()`. The corrected labels are self-evidently right:
succinate dehydrogenase activity → `Sdha` (was `Ybx3`), microfilament motor activity →
`Myh9` (was `Hadh`), SNARE binding → `Napa` (was `Clybl`), mRNA splicing → `Lsm2/Lsm6/
Snrpg/Prpf19/Ddx46` (was `Prosc/Tpd52l2/Fabp4`).

Affects `ora_top_terms_detail.json` and the UpSet figure only. `ora_up.csv` /
`ora_down.csv` are header-only, so no shipped table carried a wrong gene. **No p-value,
term, or count changes** — D6 is untouched by this. `upset.py`'s docstring and caption
referenced the old mislabelled overlap (`Prosc`/`Tpd52l2`) and now name the real one
(`Lsm2`/`Lsm6`, plus `Snrpg` between the minor and major splicing pathways).

**Also corrected:** `ora_meta.json`'s note claimed "the SAME UP query against the
default background returns 196 terms". Post-D7 that query is the DOWN query. Both
directions were re-measured against `domain_scope="annotated"` and the real numbers
(326 UP / 196 DOWN) are now recorded under a new `default_background_comparison` key,
flagged as a diagnostic that is *not* how the pipeline computes its results.

**Tests:** `tests/test_enrich_common.py`, 43 offline + 1 `@pytest.mark.network`.
Full suite: **306 passed**. The only failure is `test_freeze.py::test_tree_matches_manifest`
— expected drift from the regenerated enrichment outputs (and CHAIN-A's STRING/network
files). Not re-baselined; that is the orchestrator's call.

The header-only contract is asserted as a **contract**, in both directions: exactly 0
data rows with all 8 columns intact. Emptiness here is D6's result, not a failure, and
the test message says so — if a term ever passes, that is a finding to report, not a
test to relax. Same guard on `gsea_meta.json` (`n_fdr_lt_0.05 == 0`) and on both ORA
best-subthreshold p-values (`> 0.05`). No assertion was weakened anywhere.

**Not done / for the orchestrator:** `results/enrichment/STALE_AFTER_D7.md` left in
place as instructed. `tools/freeze.py --write` not run. Nothing committed.

---

## CHAIN-A — STRING PPI network re-run after the D7 correction (2026-08-01)
**Wave:** 2 · **Commit:** see below · **Status:** ✅ done
**Closes:** the STRING half of the stale-enrichment problem D7 created.

**What changed.** Re-ran `enrich/string_ppi.py` and `enrich/network_figure.py`
live against STRING v12.0 (species 10090, `required_score=400`) so the network's
direction annotations reflect the corrected 509 UP / 206 DOWN assignment.
De-hardcoded the three dataset literals in `load_seeds` (`715`/`206`/`509`) to
read from `tests/expected/frozen_counts.json`; two of them had silently become
wrong at D7 and would have aborted a legitimate re-run. Cached the `/json/version`
probe to `raw/string_version.json` — it was the only STRING call whose response
was never recorded, which would have left an offline-replay mode without a
version to pin.

**Why.** `results/enrichment/` described the inverted experiment. The seed SET is
unchanged (the same 715 regulated proteins), so this is a relabelling, not a
re-analysis — but every node's UP/DOWN annotation and the figure's colour mapping
were backwards.

**Verification (run, not assumed).**
- Node set identical: 694 → 694 rows, same accessions, verified by set equality.
- Direction swap is exactly clean: 496 DOWN→UP and 198 UP→DOWN, **zero
  off-diagonal** in the transition table. Nothing moved except the labels.
- `string_edges.tsv` **byte-unchanged**, which is the expected consequence of an
  identical seed set — a change there would have meant something else moved.
- `string_meta.json`: 715 seeds (509 UP / 206 DOWN), 694 mapped (97.1%),
  5963 edges, 57 communities, species 10090.
- `tests/test_string_ppi.py`: **37 passing**, 1 network-marked and deselected by
  default. Includes an AST guard that fails if anyone re-inlines a dataset count.

**Counts before → after:** UP 198 → 496, DOWN 496 → 198 (network nodes);
seeds 715 → 715; edges 5963 → 5963.

**Open follow-ups.** The agent was interrupted before writing its tests; the
orchestrator finished `test_string_ppi.py` (it needed one missing `import ast`)
and wrote this entry. No further live STRING calls were made — the re-run itself
was already complete and independently verified before the interruption.

---

## P8 — gated dispatcher reads the sample sheet (2026-08-01)
**Wave:** 3 (Tier 3, config-driven refactor) · **Commit:** see below · **Status:** ✅ done
**Closes:** the last hardcoded sample count in `gated/pca_cluster.py`; makes
research1.md's FORWARD-PATH promise executable rather than aspirational.

**What changed.**
- `N_SAMPLES = 4` / `N_GROUPS = 2` / `N_REPLICATES_PER_GROUP` now come from
  `config.design.n_samples()` / `n_groups()` / `replicates_per_group()`. The stale
  comment telling a future reader to "read from `config/sample_sheet.tsv` once it
  exists" is gone — Wave 0 built it, and this module now uses it.
- Extracted `evaluate_registry(registry, n_samples, n_replicates_per_group)` as a
  pure function; `dispatch()` is now a three-line wrapper whose only contribution
  is the design read. The gating logic had never been exercised at any n except 4,
  so the "switches itself on" behaviour had literally never been observed.
- Registry `reason_*` fields became format templates. They previously hardcoded
  `"; have 4."` and `"rank<=n-1=3"`, so at n=20 the log would have emitted a
  confidently wrong sentence. Reasons are now rendered against the n actually
  evaluated, and `_render_reason` guarantees no row is ever emitted with a blank
  reason.
- `main()` is genuinely dispatcher-driven: it executes exactly the analyses whose
  status is RUN or QC_ONLY instead of running PCA and clustering unconditionally.
  Figure banners, manifest captions, titles and `gate_status` are all built from
  the dispatched status, so a figure can no longer carry a "QC-ONLY / NOT
  INTERPRETABLE" caveat that its own gate contradicts.
- `n_components = min(3, n-1)` → `min(n_samples - 1, n_features)`. The literal 3
  was the rank bound written out longhand (research1.md:246).
- Status vocabulary `TRUSTWORTHY` → `RUN`. No committed artifact contained the old
  value (no row reaches it at n=4) and nothing else in the repo referenced it.

**Docstring.** The module claimed gated analyses "switch themselves on" as
thresholds are met. Only `pca` and `hierarchical_clustering` are implemented; the
other four registry entries have no code behind them, so growing the sheet moves
them out of SKIP in `skip_log.csv` — which is the signal to go build them — but
does not conjure the analysis. The docstring now says exactly that.

**Verification (run, not assumed).**
- **Control first:** re-ran the *unmodified* script and confirmed `freeze.py
  --check` reported 70 OK, so any later drift would be attributable to the edit.
- After the refactor: `tools/freeze.py --check` → **70 OK, 0 CHANGED**. All five
  gated artifacts (`skip_log.csv`, `pca_coords.csv`, `pca_variance.csv`,
  `pca_qc.png/.svg`, `sample_dendrogram.png/.svg`, `figures_manifest_gated.json`)
  byte-identical. Script **stdout is also byte-identical** to the baseline.
- The two SVGs show as modified in `git status` purely from matplotlib's
  `<dc:date>` + salted element ids; they hash OK under `freeze.py`'s canonical
  mode and were restored rather than committed as noise.
- `tests/test_gating.py`: **22 passing**. Full suite: **329 passing, 2 deselected**.
- **Mutation-checked, so the suite is known to bite:** re-hardcoding
  `N_SAMPLES = 4` fails 1 test; pinning the status so the gate never flips to RUN
  fails 3. Neither mutant passed silently.

**Gating table (the deliverable).**

| analysis | n=4 (committed) | n=6, 3/group | n=20, 10/group |
|---|---|---|---|
| pca | QC_ONLY | RUN | RUN |
| hierarchical_clustering | QC_ONLY | RUN | RUN |
| wgcna | SKIP | SKIP | RUN |
| umap_tsne | SKIP | SKIP | SKIP |
| mixomics_splsda | SKIP | SKIP | RUN |
| vae_gnn | SKIP | SKIP | SKIP |

n=6 and n=20 are synthetic sheets written to `tmp_path` and read back through
`config.design`, so what is under test is the real TSV → gate path. At n=6 exactly
two analyses flip and no others; WGCNA's 15-to-run / 20-to-trust band is covered
by a separate n=16 → QC_ONLY case.

**D7 check.** PCA and clustering are direction-agnostic, and the outputs were
confirmed byte-identical, so nothing here depends on the control/treated
assignment — as expected.

**Not done / for the orchestrator:** `tools/freeze.py --write` not run (no
re-baseline needed — nothing changed). Only `gated/pca_cluster.py` and the new
`tests/test_gating.py` were touched.
