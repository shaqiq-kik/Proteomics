"""Shared pytest fixtures for the proteomics_de test suite.

Two things worth knowing before adding tests here:

1. The production scripts use *flat*, sys.path-relative imports. ``viz/*.py``
   does a bare ``import style``; ``enrich/*.py`` does ``import enrich_common``;
   ``qc/validate.py`` does ``import schema``. Meanwhile
   ``enrich/network_figure.py`` does ``from viz.style import ...`` and
   ``centering_check.py`` does ``from proteomics_de.foldchange import ...``.
   Nothing is installed as a package. The autouse ``_pipeline_sys_path`` fixture
   below reproduces the sys.path layout those scripts expect, so tests can
   import any pipeline module directly.

2. Every path is resolved from ``Path(__file__)``, never from the current
   working directory. pytest may be invoked from anywhere.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import pytest

# Force a non-interactive matplotlib backend before any pipeline module (most of
# viz/ and enrich/ import pyplot at module scope) can pick a GUI one.
os.environ.setdefault("MPLBACKEND", "Agg")

# proteomics_de/tests/conftest.py -> proteomics_de/tests -> proteomics_de -> repo root
_TESTS_DIR = Path(__file__).resolve().parent
_PKG_DIR = _TESTS_DIR.parent
_REPO_ROOT = _PKG_DIR.parent

# Order matters: most specific first. `proteomics_de/` must be present for
# `from viz.style import ...`, and the repo root for `from proteomics_de.foldchange
# import ...`; the leaf dirs make the bare `import style` / `import enrich_common`
# / `import schema` forms resolve.
_SYS_PATH_ENTRIES = (
    _PKG_DIR / "viz",
    _PKG_DIR / "enrich",
    _PKG_DIR / "qc",
    _PKG_DIR,
    _REPO_ROOT,
)


@pytest.fixture(scope="session", autouse=True)
def _pipeline_sys_path():
    """Make every pipeline module importable, the way the scripts expect."""
    added = []
    for entry in _SYS_PATH_ENTRIES:
        s = str(entry)
        if s not in sys.path:
            sys.path.insert(0, s)
            added.append(s)
    yield
    for s in added:
        try:
            sys.path.remove(s)
        except ValueError:
            pass


@pytest.fixture(scope="session")
def repo_root() -> Path:
    """Absolute path to the repository root (resolved from __file__, not cwd)."""
    return _REPO_ROOT


@pytest.fixture(scope="session")
def results_dir(repo_root: Path) -> Path:
    """proteomics_de/results/ -- the committed pipeline outputs."""
    return repo_root / "proteomics_de" / "results"


@pytest.fixture(scope="session")
def frozen_counts() -> dict:
    """Expected row/class counts for the committed dataset.

    See tests/expected/frozen_counts.json. Keys starting with "_" are commentary.
    """
    path = _TESTS_DIR / "expected" / "frozen_counts.json"
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


# NOTE: this module used to also provide a `protected_sha256` fixture reading
# tests/expected/protected.sha256. No test ever used it -- the real byte-freeze
# gate is tools/freeze.py against tests/expected/outputs.sha256, exercised
# directly by test_freeze.py. Removed along with the manifest itself.


@pytest.fixture(scope="session")
def has_rscript() -> bool:
    """True when Rscript is on PATH (limma/imputeLCMD steps can run)."""
    return shutil.which("Rscript") is not None
