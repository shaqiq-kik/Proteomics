"""Stage-boundary validation hooks.

Wave-0 shim. P6 replaces the body with real pandera validation at each named
boundary. Call sites installed in Wave 1 must not change when the body is filled
in: :func:`check` returns its input frame unchanged and performs no validation,
no copying, and no I/O.
"""

from __future__ import annotations

#: The boundaries a stage may validate at. P6 will attach a schema to each.
#:
#: * ``after_load``        -- raw sheets read from the input workbook
#: * ``after_merge``       -- L/H outer merge, before fold changes exist
#: * ``after_foldchange``  -- log2FC / regulated / onoff columns populated
#: * ``before_ipa_export`` -- the filtered frame handed to the IPA writer
STAGES = ("after_load", "after_merge", "after_foldchange", "before_ipa_export")


def check(stage, df, schema=None):
    """Validate `df` at the named `stage` and return it.

    Wave-0 shim: returns `df` unchanged, with no validation and no side effects.
    Deliberately not even a defensive copy -- a copy would change identity and
    could perturb downstream ``.loc`` assignment behaviour, and this shim must be
    genuinely behaviour-neutral.

    Parameters
    ----------
    stage :
        One of :data:`STAGES`.
    df :
        The frame crossing the boundary.
    schema :
        Optional explicit schema override; ignored by the shim.

    Returns
    -------
    The same object that was passed in as `df`.
    """
    return df


__all__ = ["STAGES", "check"]
