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
