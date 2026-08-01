"""Provenance sidecars for pipeline outputs.

Wave-0 shim. P7 replaces the body with the real writer. Call sites installed in
Wave 1 must not change when the body is filled in: :func:`sidecar` currently does
nothing at all and, critically, writes NO file -- an extra file on disk would
break the "all committed outputs stay sha256-identical" acceptance test.

Intent (P7): for an output at ``<path>``, write ``<name>.provenance.json``
alongside it carrying

* ``caveat``       -- ``config.constants.CAVEAT_TEXT``, so the n=2 technical-
  replicate limitation travels with the data and not only with the figures
* ``git_commit``   -- HEAD sha of the tree that produced the file
* ``input_sha256`` -- digest of the upstream input(s)
* ``n_rows``       -- row count of the emitted table
* ``tool_versions`` -- python / pandas / numpy / R / limma / imputeLCMD versions
"""

from __future__ import annotations


def sidecar(path, **facts):
    """Record provenance for the output at `path`.

    Wave-0 shim: a no-op returning ``None``. Writes nothing, reads nothing.

    Parameters
    ----------
    path :
        The output file the sidecar will describe.
    **facts :
        Arbitrary provenance key/values (row counts, input digests, versions).
        Accepted and discarded by the shim so that Wave-1 call sites can already
        pass the real facts.

    Returns
    -------
    None
    """
    return None


__all__ = ["sidecar"]
