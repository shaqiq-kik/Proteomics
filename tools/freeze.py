"""Canonical hashing for the byte-freeze gate.

WHY THIS EXISTS
---------------
The pipeline's scientific artifacts must not change when we refactor the code
that produces them. A sha256 manifest is how we prove that. Two lessons from
the first attempt are baked in here:

1. **Freeze outputs, not sources.** The original manifest hashed all 93 tracked
   files, including 21 .py/.R source files. That makes the gate fail by
   construction the moment anyone refactors a script -- which is exactly the
   activity the gate is supposed to make safe. Git already versions sources.
   Only `proteomics_de/results/**` and the `_limma_*` handoff intermediates are
   frozen here.

2. **SVGs are not byte-reproducible, and never were.** matplotlib stamps a
   wall-clock `<dc:date>` into every SVG and salts element ids with random hex.
   Two consecutive runs of *unmodified* code produce different SVG bytes. This
   was verified by control experiment (stash the change, re-run, still drifts),
   so it is a property of matplotlib, not of any edit we made. PNGs, by
   contrast, are byte-identical run to run.

   Rather than exclude SVGs from the gate (which would leave 13 figures
   unguarded) we hash their *canonical* form: the date is dropped and each
   distinct salted id is renumbered by first-appearance order. That removes the
   volatile parts while preserving structure -- a genuine change to a plot still
   changes the hash, because the path data, transforms and embedded rasters are
   all still hashed verbatim.

The manifest records which mode each file used, so the normalization is never
silent.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

#: Files whose bytes are the scientific record. Refactors must not change these.
FROZEN_GLOBS = (
    "proteomics_de/results/**/*",
    "proteomics_de/_limma_input.csv",
    "proteomics_de/_limma_output.csv",
    "proteomics_de/_limma_versions.txt",
)

#: Deliberately NOT frozen: provenance sidecars.
#
# A sidecar records WHEN a deliverable was generated and FROM WHAT STATE
# (git_commit, git_dirty, generated_at). Those are the point of the file, and
# they legitimately differ between runs -- freezing them would mean the gate
# reported drift after every single run, which is how a gate gets ignored.
#
# They are not left unverified. `test_provenance.py` enforces the property that
# actually matters: each sidecar's recorded sha256 must equal the current hash
# of the file it describes. That is strictly stronger than byte-freezing the
# sidecar, because it fails when the DESCRIBED file changes -- which is the
# failure anyone cares about, and it caught exactly that during D9.
FROZEN_EXCLUDE_SUFFIXES = (".provenance.json",)

# matplotlib's SVG volatility.
_SVG_DATE = re.compile(rb"<dc:date>[^<]*</dc:date>")

#: Every id matplotlib DECLARES, whatever its prefix.
#
# This used to enumerate prefixes (`p`, `m`, `image`, the font names). That was
# a guess dressed as a rule, and it silently missed upsetplot's marker
# collections, which declare ids like `C0_0_cc31541d2c` -- so upset.svg reported
# as drifted after every run and the gate could not verify it at all.
#
# Reading the declarations instead of guessing their shape is exact: an id that
# matplotlib salts must appear in an id="..." attribute, so the rename map is
# built from what the file actually declares. The 10-hex-suffix requirement is
# what distinguishes a salted id from a stable one (e.g. the `id1` clip paths).
_SVG_DECLARED_ID = re.compile(rb'id="([^"]*[0-9a-f]{10})"')

#: Fields whose value is a wall-clock stamp: real provenance, but not content.
_TIMESTAMP_KEYS = ("generated_at", "generated", "timestamp", "created_at", "run_at")
_JSON_TIMESTAMP = re.compile(
    rb'("(?:' + b"|".join(k.encode() for k in _TIMESTAMP_KEYS) + rb')"\s*:\s*)"[^"]*"'
)
_MD_TIMESTAMP = re.compile(
    rb"(?im)^(.*(?:generated|timestamp|run at)[^\n]*?)\d{4}-\d{2}-\d{2}[T ][0-9:.+Z-]*"
)


def canonical_bytes(path: Path) -> tuple[bytes, str]:
    """Return (bytes_to_hash, mode) for *path*.

    Three modes. ``raw`` hashes the bytes. ``svg-canon`` drops matplotlib's
    ``<dc:date>`` and renumbers its salted element ids. ``ts-canon`` blanks
    wall-clock stamps in files that record when they were generated -- those
    stamps are genuine provenance and are kept IN the file, but hashing them
    would mean every re-run reports drift, which trains people to ignore the
    gate. Everything else in those files is still hashed verbatim.
    """
    data = path.read_bytes()
    suffix = path.suffix.lower()

    if suffix == ".svg":
        data = _SVG_DATE.sub(b"<dc:date></dc:date>", data)
        # Renumber by first appearance so distinct ids stay distinct: a blanket
        # substitution would collapse them and could mask a real change.
        mapping: dict[bytes, bytes] = {}
        for m in _SVG_DECLARED_ID.finditer(data):
            token = m.group(1)
            if token not in mapping:
                mapping[token] = b"canonid%d" % len(mapping)
        # Longest-first so a shorter id cannot corrupt a longer one containing it.
        for token in sorted(mapping, key=len, reverse=True):
            data = data.replace(token, mapping[token])
        return data, "svg-canon"

    if suffix == ".json" and _JSON_TIMESTAMP.search(data):
        return _JSON_TIMESTAMP.sub(rb'\1"<ts>"', data), "ts-canon"

    if suffix == ".md" and _MD_TIMESTAMP.search(data):
        return _MD_TIMESTAMP.sub(rb"\1<ts>", data), "ts-canon"

    return data, "raw"


def digest(path: Path) -> tuple[str, str]:
    """Return (sha256_hex, mode) of *path*'s canonical form."""
    data, mode = canonical_bytes(path)
    return hashlib.sha256(data).hexdigest(), mode


def frozen_files(root: Path | None = None) -> list[Path]:
    """Every frozen artifact, sorted, as absolute paths."""
    root = root or REPO_ROOT
    found: set[Path] = set()
    for pattern in FROZEN_GLOBS:
        for p in root.glob(pattern):
            if p.is_file() and not p.name.endswith(FROZEN_EXCLUDE_SUFFIXES):
                found.add(p)
    return sorted(found)


def build_manifest(root: Path | None = None) -> list[tuple[str, str, str]]:
    """Return [(relpath, sha256, mode)] for every frozen artifact."""
    root = root or REPO_ROOT
    rows = []
    for p in frozen_files(root):
        sha, mode = digest(p)
        rows.append((p.relative_to(root).as_posix(), sha, mode))
    return rows


def write_manifest(path: Path, root: Path | None = None) -> int:
    """Write the manifest to *path*. Returns the number of entries."""
    rows = build_manifest(root)
    lines = [
        "# Byte-freeze manifest for proteomics_de scientific outputs.",
        "# Format: <sha256>  <mode>  <relpath>",
        "# mode=raw       -> sha256 of the file's bytes",
        "# mode=svg-canon -> sha256 after dropping matplotlib's <dc:date> and",
        "#                   renumbering its randomly-salted element ids, which",
        "#                   differ between runs of identical code.",
        "# Source files are deliberately NOT frozen; git versions those.",
        f"# entries: {len(rows)}",
    ]
    lines += [f"{sha}  {mode}  {rel}" for rel, sha, mode in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(rows)


def read_manifest(path: Path) -> dict[str, tuple[str, str]]:
    """Parse a manifest into {relpath: (sha256, mode)}."""
    out: dict[str, tuple[str, str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        sha, mode, rel = line.split(None, 2)
        out[rel.strip()] = (sha, mode)
    return out


def verify(manifest_path: Path, root: Path | None = None):
    """Compare the tree against *manifest_path*.

    Returns (ok, changed, missing, extra) as lists of relpaths.
    """
    root = root or REPO_ROOT
    expected = read_manifest(manifest_path)
    ok, changed, missing = [], [], []
    for rel, (sha, _mode) in sorted(expected.items()):
        p = root / rel
        if not p.is_file():
            missing.append(rel)
            continue
        actual, _ = digest(p)
        (ok if actual == sha else changed).append(rel)
    present = {p.relative_to(root).as_posix() for p in frozen_files(root)}
    extra = sorted(present - set(expected))
    return ok, changed, missing, extra


if __name__ == "__main__":
    import argparse
    import sys

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--write", action="store_true", help="regenerate the manifest")
    ap.add_argument("--check", action="store_true", help="verify, nonzero exit on drift")
    ap.add_argument(
        "--manifest",
        default=str(REPO_ROOT / "proteomics_de/tests/expected/outputs.sha256"),
    )
    args = ap.parse_args()
    mpath = Path(args.manifest)

    if args.write:
        n = write_manifest(mpath)
        print(f"wrote {mpath.relative_to(REPO_ROOT)} ({n} artifacts)")
    if args.check:
        ok, changed, missing, extra = verify(mpath)
        print(f"freeze: {len(ok)} OK, {len(changed)} CHANGED, {len(missing)} MISSING, {len(extra)} UNTRACKED")
        for rel in changed:
            print(f"  CHANGED  {rel}")
        for rel in missing:
            print(f"  MISSING  {rel}")
        for rel in extra:
            print(f"  UNTRACKED {rel}")
        sys.exit(1 if (changed or missing) else 0)
