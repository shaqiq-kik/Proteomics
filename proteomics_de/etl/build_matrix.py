"""Emit the limma file contract: ``intensity_matrix.tsv`` + ``design.tsv``.

research1.md Section 6 item 7, against the Section 1 "Output schemas / naming
conventions" contract::

    intensity_matrix.tsv : accession, gene, <sample_1..n>   (one row per protein)
    design.tsv           : sample, group                    (group in {control, treated})

Both files are TAB-separated, and **both are derived entirely from
``config/sample_sheet.tsv``** -- the sample columns are named and ordered by
:func:`proteomics_de.config.design.sample_columns`, so adding biological
replicates is an edit to the sheet and nothing else. That is the whole point of
the module: the file contract is replicate-count-agnostic by construction.

This is a *library*, not a pipeline stage. ``limma_test.py`` imports it and calls
:func:`build` inline, which keeps the existing single-entry-point chain
(``foldchange.py`` -> ``limma_test.run_limma_test``) untouched. The ``__main__``
block exists only so the two files can be regenerated and eyeballed on their own.

**These are additive files.** Nothing reads them yet -- the R worker is fed by
the older ``_limma_input.csv`` handoff, which stays byte-frozen. They exist so
the contract is on disk and testable now, ahead of the handoff migrating onto it.

**Known contract wart.** research1.md's R sketch does
``intensity_cols <- design_df$sample`` -- i.e. it assumes ``design.tsv``'s
``sample`` column *names the matrix columns*. It does not, here:
``design.tsv`` carries bare sample ids (``31578``) while the matrix columns are
channels (``Intensity 31578``). The two contract files are therefore joinable
only through the sample sheet, not directly. ``limma_test.py`` consequently
hands the R worker its own ``_limma_design.tsv``, whose ``sample`` column holds
the handoff names (``ctrl_31578``) that really do match its matrix. Resolving
the wart means changing ``design.write_design_tsv``, which is out of scope here.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:  # importable however the caller was launched
    sys.path.insert(0, str(_REPO_ROOT))

from proteomics_de.config import design  # noqa: E402

#: Subdirectory of ``results/`` that holds the limma file contract.
DE_SUBDIR = "de"

INTENSITY_MATRIX_NAME = "intensity_matrix.tsv"
DESIGN_NAME = "design.tsv"

#: Source columns in ``foldchange_all.csv`` for the contract's first two columns.
ACCESSION_COL = "UniProt Accession Number"
GENE_COL = "Gene names"


def intensity_series(series) -> pd.Series:
    """Raw intensity -> float, with blank/non-numeric/``<=0`` turned into NaN.

    The single missing-value rule for everything crossing the limma boundary.
    NaN is written by ``to_csv`` as an empty cell (the default ``na_rep``), which
    R reads back as ``NA`` so MinProb fills it -- an intensity of 0 means "below
    the detection limit" (MNAR), not "measured as zero", so it must not survive
    the ``log2`` as ``-Inf``.

    ``limma_test._missing_to_blank`` delegates here, so the ``_limma_input.csv``
    handoff and ``intensity_matrix.tsv`` cannot drift apart.
    """
    num = pd.to_numeric(series, errors="coerce")
    return num.where(num > 0)  # keep > 0; everything else (<=0, NaN) -> NaN


def build(eligible_df, outdir, sheet=None) -> dict[str, Path]:
    """Write ``<outdir>/de/intensity_matrix.tsv`` and ``<outdir>/de/design.tsv``.

    Parameters
    ----------
    eligible_df :
        The proteins to test -- the both-condition, non-ON_OFF rows of
        ``foldchange_all.csv``. Must carry :data:`ACCESSION_COL`,
        :data:`GENE_COL` and every channel named by
        :func:`~proteomics_de.config.design.sample_columns`. Row order is
        preserved verbatim; this function never sorts or filters.
    outdir :
        The results directory (e.g. ``proteomics_de/results``). The ``de/``
        subdirectory is created if needed.
    sheet :
        Optional sample-sheet override (a DataFrame or a path), forwarded to
        ``config.design``. Defaults to the committed sheet.

    Returns
    -------
    dict
        ``{"intensity_matrix": Path, "design": Path}``.

    Raises
    ------
    ValueError
        If `eligible_df` is empty or is missing a required column. Emitting a
        truncated matrix would silently under-test the experiment, so this fails
        loudly instead.
    """
    sample_cols = design.sample_columns(sheet)

    if len(eligible_df) == 0:
        raise ValueError(
            "build_matrix: eligible_df has zero rows -- refusing to write an "
            "empty intensity_matrix.tsv."
        )
    required = [ACCESSION_COL, GENE_COL] + sample_cols
    missing = [c for c in required if c not in eligible_df.columns]
    if missing:
        raise ValueError(
            f"build_matrix: eligible_df is missing column(s): {missing}.\n"
            f"  expected: {required}\n"
            f"  present : {list(eligible_df.columns)}\n"
            "The sample columns come from config/sample_sheet.tsv; if the sheet "
            "changed, the upstream fold-change table must be regenerated too."
        )

    de_dir = Path(outdir) / DE_SUBDIR
    de_dir.mkdir(parents=True, exist_ok=True)

    matrix = pd.DataFrame(
        {
            "accession": eligible_df[ACCESSION_COL].values,
            "gene": eligible_df[GENE_COL].values,
        }
    )
    for col in sample_cols:
        matrix[col] = intensity_series(eligible_df[col]).values

    matrix_path = de_dir / INTENSITY_MATRIX_NAME
    matrix.to_csv(
        matrix_path, sep="\t", index=False, encoding="utf-8", lineterminator="\n"
    )

    design_path = design.write_design_tsv(de_dir / DESIGN_NAME, sheet)

    return {"intensity_matrix": matrix_path, "design": Path(design_path)}


if __name__ == "__main__":  # pragma: no cover - standalone inspection only
    # Regenerate both contract files from the committed fold-change table, using
    # the same eligibility rule as limma_test (ON_OFF proteins are not testable:
    # they are absent from one condition entirely).
    _RESULTS = Path(__file__).resolve().parents[1] / "results"
    _FC = _RESULTS / "foldchange_all.csv"

    fc = pd.read_csv(_FC, dtype=str, keep_default_na=False)
    eligible = fc[fc["onoff"].str.strip() == ""].reset_index(drop=True)

    written = build(eligible, _RESULTS)
    print(f"[build_matrix] {len(eligible)} eligible proteins")
    for key, path in written.items():
        print(f"  {key:17s}: {path.relative_to(_REPO_ROOT)}")
