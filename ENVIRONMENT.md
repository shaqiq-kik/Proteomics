# Environment

How to get a working environment for `proteomics_de/` on this machine (macOS
arm64), and the two traps that will bite you if you skip it.

---

## 1. Setup

```sh
python3 -m venv --system-site-packages .venv
.venv/bin/pip install -r requirements-dev.txt
```

That is the whole setup. Then always call the tools by path:

```sh
.venv/bin/python   ...
.venv/bin/pytest   ...
```

`--system-site-packages` is **essential, not optional.** The scientific stack
(pandas, numpy, scipy, scikit-learn, pandera, matplotlib, seaborn, gseapy,
networkx, upsetplot, adjustText, openpyxl, PyYAML) is installed into the
Homebrew interpreter, not into the venv. A plain `python3 -m venv .venv` creates
an isolated environment that cannot see any of it, and every pipeline import
fails immediately. The venv exists only to hold `pytest` + `pytest-xdist`; it
borrows everything else.

Exact versions the committed results were produced with are recorded in
`requirements-lock.txt`.

---

## 2. Trap: Homebrew Python is externally managed (PEP 668)

`/opt/homebrew/bin/python3` is marked externally-managed, so installing into it
directly **fails**:

```
error: externally-managed-environment
```

Do not work around this with `--break-system-packages` or `--user`; it puts the
interpreter into a state Homebrew will fight with on the next `brew upgrade`.
The repo-local venv above is the fix. `pip install` targets `.venv/bin/pip`
only, and only for test tooling.

---

## 3. Trap: four different `python3` on PATH

```
/opt/homebrew/bin/python3                                   3.13.7   <-- has the stack
/Library/Frameworks/Python.framework/Versions/3.13/bin/python3  3.13.3   no pandas
/usr/local/bin/python3                                      3.13.3   no pandas
/usr/bin/python3                                            3.9.6    no pandas
```

Only the Homebrew one has the scientific packages. Which one a bare `python3`
resolves to depends entirely on PATH order in the shell that happens to be
running, so a bare `python3` is a coin flip that usually lands on
`ModuleNotFoundError: No module named 'pandas'`.

Rules:

- In the shell: call `.venv/bin/python`, never bare `python3`.
- In code that shells out to Python (runners, subprocess calls, test helpers):
  use `sys.executable`, never the string `"python3"`.

The same applies to `pytest`: use `.venv/bin/pytest`.

---

## 4. R side

The limma step (`proteomics_de/limma_test.py` -> `proteomics_de/limma_test.R`)
shells out to `Rscript`.

```
Rscript      /opt/homebrew/bin/Rscript
R            4.6.1 (2026-06-24)
limma        3.68.4
imputeLCMD   2.1
```

Install:

```r
install.packages("BiocManager")
BiocManager::install("limma")
install.packages("imputeLCMD")
```

These versions are also written next to the results in
`proteomics_de/_limma_versions.txt`, which is what the provenance checks read.
Tests that need R are marked `r`; they should skip cleanly when `Rscript` is not
on PATH rather than fail.

---

## 5. Running tests

```sh
.venv/bin/pytest                 # default: everything except network + slow
.venv/bin/pytest -m golden       # byte-identity / sha256 freeze checks
.venv/bin/pytest -m network      # live STRING / g:Profiler / Enrichr calls
.venv/bin/pytest -m "r"          # limma / imputeLCMD round-trip
.venv/bin/pytest -n auto         # parallel (pytest-xdist)
```

The default `addopts` in `pyproject.toml` is `-m 'not network and not slow'`, so
a bare `pytest` run is offline and fast. Deselecting happens by marker, so to
add a marked group back in you must name it explicitly as above.

Markers (declared in `pyproject.toml`):

| marker    | meaning                                                        |
|-----------|----------------------------------------------------------------|
| `network` | makes live outbound HTTP calls (STRING / g:Profiler / Enrichr) |
| `slow`    | takes more than ~10 seconds                                     |
| `r`       | requires Rscript + limma + imputeLCMD on PATH                   |
| `golden`  | byte-identity / sha256 freeze checks; **local only**            |
| `e2e`     | end-to-end pipeline run                                         |

**`golden` must not run in CI.** Those tests assert sha256 identity of committed
result files. A different BLAS/LAPACK build behind numpy/scipy, or a different R
build, changes the last float digit, which changes the CSV bytes, which changes
the hash — a real difference in the environment, not a regression in the code.
On a CI runner, deselect them:

```sh
pytest -m "not network and not slow and not golden"
```

`testpaths` is pinned to `proteomics_de/tests` so pytest never wanders into
`Pilot Project/`, which contains loose legacy scripts that error on import.
