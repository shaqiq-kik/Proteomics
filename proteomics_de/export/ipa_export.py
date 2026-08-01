"""IPA export writer.

Wave-0 shim. P7 replaces the body with validated export + p-value/FDR merge.
Call sites installed in Wave 1 must not change when the body is filled in.
"""

from __future__ import annotations


def write_ipa(df, path, columns):
    """Write `df`'s `columns` to `path` as the IPA input CSV.

    Wave-0 shim: this is literally today's behaviour at ``foldchange.py:171``
    (``df[cols].to_csv(path, index=False, encoding="utf-8")``), byte for byte.
    It exists so Wave 1 can install the call site without changing any output;
    the acceptance test for that wave is that all 13 committed output files stay
    sha256-identical.

    Parameters
    ----------
    df :
        Already-filtered frame. This function does NOT filter, sort, or
        reorder -- the caller owns row selection.
    path :
        Destination CSV path.
    columns :
        Column names to emit, in output order.

    Returns
    -------
    None
    """
    df[columns].to_csv(path, index=False, encoding="utf-8")


__all__ = ["write_ipa"]
