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
# Single source of the caveat text. It is stamped onto 10+ committed figures and
# written into every output's .provenance.json sidecar; changing a single
# character (including the em dash) changes those figures' bytes, so the frozen
# hashes in tests/expected/outputs.sha256 must be refrozen with it.
# Corrected by DECISIONS_LOG D23 (2026-11-11): Dr. Walker stated the four samples
# are 2 BIOLOGICAL replicates per group from 2 separate experiments, replacing the
# earlier technical-replicate wording. If the design changes again, this
# text is regenerated, not hand-edited.
CAVEAT_TEXT = (
    "SILAC design: n = 2 biological replicates per group, from 2 separate "
    "experiments — 0 / 1938 proteins pass FDR < 0.05 (limma); n = 2 is low "
    "power, so results below are exploratory / hypothesis-generating, not "
    "confirmed significant."
)

# NOTE: this module used to also expose PROTECTED_FILES, parsed from a
# tests/expected/protected.sha256 that froze all 93 tracked source AND output
# files. Nothing in the codebase ever consumed the parsed list -- the actual
# byte-freeze gate is tools/freeze.py + tests/expected/outputs.sha256, scoped to
# scientific outputs only (see that module's docstring for why). The unused
# file caused two separate incidents of confusion before it was deleted: a
# runner banner that named it while hashing something else, and two separate
# investigations mistaking it for the gate. Removed rather than kept "for
# reference" -- git already tracks source history better than a hand-copied
# hash list, and a file nothing reads cannot help anyone despite looking load-
# bearing.

__all__ = [
    "LOG2_THRESHOLD",
    "ADJ_P_THRESHOLD",
    "RAW_P_THRESHOLD",
    "R_SEED",
    "MINPROB_Q",
    "MINPROB_TUNE_SIGMA",
    "CAVEAT_TEXT",
]
