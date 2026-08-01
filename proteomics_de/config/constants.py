"""Single source of truth for values currently duplicated across the pipeline.

Every value below already exists as a literal in at least one frozen script. The
source of each is named in its comment. This module does not change any of them
-- it exists so that a future edit happens in exactly one place instead of being
copy-pasted into nine files.

Nothing here is wired into the frozen scripts yet; that happens in a later wave.
"""

from __future__ import annotations

from pathlib import Path

# ---------------------------------------------------------------------------
# Statistical thresholds and seeds
# ---------------------------------------------------------------------------

#: |log2FC| call boundary for UP/DOWN. = log2(1.5); the DOWN side is the
#: symmetric -0.585 (Bug 2 fix). Source: ``foldchange.py:24`` (also duplicated as
#: ``FC_THRESHOLD`` in ``viz/style.py:143``).
LOG2_THRESHOLD = 0.585

#: FDR cutoff for calling a protein "significant". Source: ``limma_test.py:39``.
ADJ_P_THRESHOLD = 0.05

#: Raw (unadjusted) p-value cutoff used only for figure shading, never for
#: significance calls. Source: ``viz/style.py:144``.
RAW_P_THRESHOLD = 0.05

#: Seed passed through to R so ``imputeLCMD::impute.MinProb`` is reproducible.
#: Source: ``limma_test.py:40``, consumed at ``limma_test.R`` ``set.seed(seed)``.
R_SEED = 42

# --- MinProb imputation parameters (limma_test.R:87) ------------------------
# impute.MinProb(mat, q = 0.01, tune.sigma = 1)

#: Quantile of the observed distribution that the imputed values are drawn
#: around. Source: ``limma_test.R:87``.
MINPROB_Q = 0.01

#: Width multiplier for the imputation draw. Source: ``limma_test.R:87``.
MINPROB_TUNE_SIGMA = 1

# ---------------------------------------------------------------------------
# Figure caveat
# ---------------------------------------------------------------------------
# VERBATIM copy of ``viz/style.py:137-141``. This string is stamped onto 10+
# committed figures; changing a single character (including the em dash) changes
# those figures' bytes and breaks the protected-file check below. If the design
# ever gains biological replicates, this text is regenerated, not hand-edited.
CAVEAT_TEXT = (
    "SILAC design: n = 2 technical replicates per condition (no biological "
    "replication) — 0 / 1938 proteins pass FDR < 0.05 (limma); results below "
    "are hypothesis-generating only, not confirmed significant."
)

# ---------------------------------------------------------------------------
# Protected files
# ---------------------------------------------------------------------------

_PROTECTED_MANIFEST = (
    Path(__file__).resolve().parent.parent / "tests" / "expected" / "protected.sha256"
)


def _read_protected(manifest: Path = _PROTECTED_MANIFEST) -> list[str]:
    """Parse the ``sha256sum``-format manifest into a list of repo-relative paths.

    The manifest is authored by the refactor's baseline step and is NEVER
    regenerated from here -- this function only reads it, so that the list of
    protected paths cannot silently drift by being recomputed against a
    already-modified tree.
    """
    if not manifest.exists():
        return []
    paths: list[str] = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        # "<64-hex>  <path>"  (two spaces, per sha256sum's text mode)
        _digest, _sep, path = line.partition("  ")
        path = path.strip()
        if path:
            paths.append(path)
    return paths


#: Repo-relative paths that must stay **byte-identical** through the refactor.
#:
#: "Protected" means: the file is either a frozen, individually-verified pipeline
#: script or one of its committed outputs (CSV, JSON, PNG/SVG figure, report).
#: The refactor is allowed to add new modules and to change *how* these files are
#: produced, but re-running the pipeline must reproduce every one of them with an
#: unchanged sha256. Any diff in this list is a regression, not an improvement,
#: and must be explained before it is accepted.
#:
#: Derived by parsing ``proteomics_de/tests/expected/protected.sha256``.
PROTECTED_FILES = _read_protected()

__all__ = [
    "LOG2_THRESHOLD",
    "ADJ_P_THRESHOLD",
    "RAW_P_THRESHOLD",
    "R_SEED",
    "MINPROB_Q",
    "MINPROB_TUNE_SIGMA",
    "CAVEAT_TEXT",
    "PROTECTED_FILES",
]
