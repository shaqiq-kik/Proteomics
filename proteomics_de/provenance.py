"""Provenance sidecars for pipeline outputs.

The problem this solves
-----------------------
research1.md line 245 is unambiguous about the n=2 limitation: *"This is the
single most important interpretive caveat and must be stated on every output."*
The figures comply -- :data:`config.constants.CAVEAT_TEXT` is stamped onto 10+
committed PNG/SVGs. **The CSVs did not.** ``results/ipa_input.csv`` is the file
that actually leaves this repo (it gets uploaded to QIAGEN IPA and mailed to
collaborators) and it carried no indication that its 715 "regulated" proteins
are hypothesis-generating leads from two technical replicates, none of which
survive FDR correction.

The obvious fix -- a ``#`` comment header inside the CSV -- is unavailable. IPA's
file contract (research1.md SECTION 3) is a single header row followed by data;
a comment line becomes the header and the upload fails to recognize the columns.
That is precisely the failure mode the earlier Excel->CSV fix was chasing.

So the caveat travels in a **sidecar**: for an output at ``<path>``, a sibling
``<path>.provenance.json``. The data file stays byte-identical and contract-clean;
anyone who receives the pair gets the caveat, the exact commit that produced the
numbers, a checksum to prove the file was not edited in transit, and the tool
versions needed to reproduce it.

What a sidecar carries
----------------------
``caveat``
    :data:`config.constants.CAVEAT_TEXT`, **verbatim**. Deliberately not
    reworded here -- it is the same string stamped on the figures, so a reader
    comparing a CSV to a plot sees identical wording, and there is exactly one
    place to edit if the design ever gains biological replicates.
``git_commit`` / ``git_dirty``
    HEAD of the tree that produced the file, and whether that tree had
    uncommitted changes (a clean sha alone would be a false promise otherwise).
``sha256``
    Digest of the file itself, so tampering or truncation is detectable.
``n_rows``
    Data rows, excluding the header. Cross-checks the file against the counts in
    ``tests/expected/frozen_counts.json``.
``generated_at``
    UTC ISO-8601 timestamp.
``tool_versions``
    python / pandas / numpy from the live interpreter, plus R / limma /
    imputeLCMD read from the ``_limma_versions.txt`` handoff that the R step
    already writes.

Sidecars are NEW files under ``results/``. They are expected to show up as
untracked in the freeze check; they never modify a frozen output.

Usage::

    python -m proteomics_de.provenance            # emit for the four deliverables
    python -m proteomics_de.provenance --path results/foo.csv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_HERE = Path(__file__).resolve().parent          # proteomics_de/
_ROOT = _HERE.parent                             # repo root

if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from proteomics_de.config.constants import CAVEAT_TEXT  # noqa: E402

#: Suffix appended to the full filename (not the stem): ``ipa_input.csv`` ->
#: ``ipa_input.csv.provenance.json``. Keeping the extension means
#: ``ipa_input_full.csv`` and ``ipa_input_full.txt`` get distinct sidecars
#: instead of silently overwriting each other.
SIDECAR_SUFFIX = ".provenance.json"

#: The deliverables that leave the repo, and therefore must carry the caveat.
#: ``ipa_input.csv`` and ``ipa_input_full.csv`` go to QIAGEN; ``qc_limma.csv``
#: and ``foldchange_all.csv`` are the tables collaborators actually open.
DEFAULT_DELIVERABLES = (
    "ipa_input.csv",
    "ipa_input_full.csv",
    "qc_limma.csv",
    "foldchange_all.csv",
)

#: Where the R step records its versions (written by ``limma_test.py``).
_LIMMA_VERSIONS = _HERE / "_limma_versions.txt"

_DELIMITERS = {".csv": ",", ".txt": "\t", ".tsv": "\t"}


def sha256_file(path, chunk_size=1 << 20):
    """sha256 of `path`, streamed so a large figure or matrix is not slurped."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def count_rows(path):
    """Data-row count for a delimited text file, or ``None`` if not one.

    Uses :mod:`csv` rather than counting newlines: a quoted field may legally
    contain a newline, and miscounting rows in a provenance record is worse than
    declining to count.
    """
    suffix = Path(path).suffix.lower()
    if suffix not in _DELIMITERS:
        return None
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter=_DELIMITERS[suffix])
        n = sum(1 for _ in reader)
    return max(n - 1, 0)  # minus the single header row


def git_commit(repo_dir=None):
    """``(sha, dirty)`` for the tree containing the outputs.

    ``dirty`` counts **tracked** modifications only (``--untracked-files=no``).
    The sidecars themselves are new untracked files under ``results/``, so
    including untracked paths would make every sidecar report ``dirty: true``
    because of its own existence. The question worth answering is "did the
    committed sources or outputs differ from HEAD when this ran", and that is
    what the tracked-only form asks.

    Returns ``(None, None)`` when git is unavailable or this is not a work tree,
    rather than raising -- a sidecar with an honest ``null`` commit is more
    useful than no sidecar at all.
    """
    repo_dir = str(repo_dir or _ROOT)
    try:
        sha = subprocess.run(
            ["git", "-C", repo_dir, "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True, timeout=15,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", repo_dir, "status", "--porcelain", "--untracked-files=no"],
            capture_output=True, text=True, check=True, timeout=30,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None, None
    return sha or None, bool(status.strip())


def _r_versions(path=_LIMMA_VERSIONS):
    """Parse ``_limma_versions.txt`` -> ``{"R": ..., "limma": ..., ...}``.

    The file is four lines written by the R step::

        R version 4.6.1 (2026-06-24)
        limma 3.68.4
        imputeLCMD 2.1
        seed 42
    """
    out = {}
    path = Path(path)
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "R" and len(parts) >= 3 and parts[1] == "version":
            out["R"] = parts[2]
        elif parts[0] in ("limma", "imputeLCMD") and len(parts) >= 2:
            out[parts[0]] = parts[1]
        elif parts[0] == "seed" and len(parts) >= 2:
            out["r_seed"] = parts[1]
    return out


def tool_versions():
    """Versions of everything that touched the numbers.

    Python-side versions come from the live interpreter (authoritative); the R
    side is read from the handoff file the R step writes, because R is a
    subprocess and may not even be installed when a sidecar is regenerated.
    """
    versions = {"python": platform.python_version()}
    for name in ("pandas", "numpy", "scipy", "matplotlib"):
        try:
            module = __import__(name)
        except ImportError:
            continue
        version = getattr(module, "__version__", None)
        if version:
            versions[name] = version
    versions.update(_r_versions())
    return versions


def sidecar(path, **facts):
    """Write ``<path><SIDECAR_SUFFIX>`` describing the output at `path`.

    Parameters
    ----------
    path :
        The deliverable to describe. Must exist -- a provenance record for a
        file that is not there would be fiction.
    **facts :
        Extra key/values merged into the record (``input_sha256``, ``contrast``,
        ``notes``, ...). Caller-supplied keys **override** the computed ones, so
        a caller who knows better (e.g. the exact in-memory row count before a
        write) can say so.

    Returns
    -------
    pathlib.Path
        The sidecar written.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"cannot write provenance for missing file: {path}")

    sha, dirty = git_commit()
    try:
        described = str(path.resolve().relative_to(_ROOT))
    except ValueError:
        described = str(path.resolve())

    record = {
        "describes": described,
        # VERBATIM from config/constants.py -- the same string stamped on every
        # figure. Do not reword it here; edit it in exactly one place.
        "caveat": CAVEAT_TEXT,
        "caveat_source": "proteomics_de/config/constants.py::CAVEAT_TEXT",
        "git_commit": sha,
        "git_dirty": dirty,
        "sha256": sha256_file(path),
        "n_rows": count_rows(path),
        "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tool_versions": tool_versions(),
    }
    record.update(facts)

    out = path.with_name(path.name + SIDECAR_SUFFIX)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(record, fh, indent=2, ensure_ascii=False, sort_keys=False)
        fh.write("\n")
    return out


def emit_default_sidecars(results_dir=None, names=DEFAULT_DELIVERABLES, **facts):
    """Write sidecars for :data:`DEFAULT_DELIVERABLES` that exist on disk.

    Missing files are skipped with a warning rather than raising: a caller may
    legitimately run this before ``ipa_input_full.csv`` has been built.

    Returns
    -------
    list[pathlib.Path]
    """
    results_dir = Path(results_dir) if results_dir else _HERE / "results"
    written = []
    for name in names:
        target = results_dir / name
        if not target.exists():
            print(f"provenance: skipping {target} (not found)", file=sys.stderr)
            continue
        written.append(sidecar(target, **facts))
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", default=None,
                    help="directory holding the deliverables "
                         "(default: proteomics_de/results)")
    ap.add_argument("--path", action="append", default=None,
                    help="write a sidecar for this file instead of the defaults "
                         "(repeatable)")
    args = ap.parse_args(argv)

    if args.path:
        written = [sidecar(p) for p in args.path]
    else:
        written = emit_default_sidecars(args.results_dir)

    for p in written:
        print(f"Saved {p}")
    return 0


__all__ = [
    "CAVEAT_TEXT",
    "SIDECAR_SUFFIX",
    "DEFAULT_DELIVERABLES",
    "sha256_file",
    "count_rows",
    "git_commit",
    "tool_versions",
    "sidecar",
    "emit_default_sidecars",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
