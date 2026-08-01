"""One documented policy for UniProt accession fields.

Today three implicit policies coexist in the codebase and nobody states them:

1. ``foldchange.py`` merges the L and H sheets on the **whole** accession string
   (``pd.merge(..., on="UniProt Accession Number")``), so ``"P05132;P68181"`` and
   ``"P05132"`` are different proteins.
2. ``enrich/enrich_common.py:86`` silently takes ``str(acc).split(";")[0]`` when
   it needs one query identifier.
3. ``qc/schema.py`` validates **every** semicolon token against a UniProt regex,
   and carves out an exception for two junk rows.

Those three are individually defensible but were never reconciled. This module
states the policy explicitly and gives every stage the same helpers.

POLICY
------
* **Merging** is on the *whole* accession string. A protein group is an identity,
  not a set -- collapsing groups would silently fuse distinct MaxQuant rows.
* **Enrichment** uses the *first token* only (:func:`first_token`), because
  g:Profiler / STRING accept a single identifier per query. This is a lossy
  projection and is documented as such at the call site.
* **Validation** checks *every* token (:func:`is_valid_uniprot` applied across
  :func:`split_group`), not just the first. A group is valid only if all of its
  members are.
* **Isoforms** (a trailing ``-N``, e.g. ``P12345-2``) are **detected and
  reported**, never silently stripped. Stripping would merge an isoform into its
  canonical entry and change protein counts. Today's data carries **zero**
  isoform suffixes, so :func:`detect_isoforms` returning a non-empty list is a
  signal that the input changed. :func:`strip_isoform` exists for callers that
  have made an explicit, logged decision -- it is not applied anywhere by default.
* **Junk index lists** (:func:`is_junk_index_list`) are **quarantined**, not
  parsed. MaxQuant occasionally emits a semicolon-joined list of *row indices*
  where an accession belongs; two such rows exist in
  ``results/single_condition_proteins.csv`` (32,759 and 681 characters). They are
  not accessions at all and must not be fed to enrichment or validation. Long
  *legitimate* protein groups exist too (e.g. the two 69-character
  ``P08752;P20612;...`` rows), so the test is "all tokens are bare integers",
  never "the string is long".

The UniProt regex is imported-by-value from ``qc/schema.py``'s settled pattern:
the task spec's regex rejects real accessions such as ``P19137``, ``Q9JHU4`` and
``O08528`` (its first character class excludes ``O``, ``P`` and ``Q``), so a
second branch ``[OPQ][A-Z0-9]{5}`` was added there. That exact pattern is reused
here so the two modules cannot disagree.
"""

from __future__ import annotations

import re

import pandas as pd

# Identical to ``qc/schema.py``'s ``UNIPROT_TOKEN_RE``. Branch 1 is the task
# spec's literal regex (non-[OPQ] accessions, 6 or 10 chars, optional isoform
# suffix); branch 2 extends it to [OPQ]-prefixed accessions, which are always 6
# characters and which the spec's version wrongly rejects on real data.
UNIPROT_TOKEN_RE = re.compile(
    r"^(?:[A-NR-Z0-9](?:[A-Z0-9]{5}|[A-Z0-9]{9})|[OPQ][A-Z0-9]{5})(?:-\d+)?$"
)

#: Trailing UniProt isoform suffix, e.g. the ``-2`` of ``P12345-2``.
ISOFORM_SUFFIX_RE = re.compile(r"-\d+$")

_INTEGER_RE = re.compile(r"^\d+$")


def _as_text(acc) -> str:
    """Normalise any cell value to the text ``str(acc)`` would produce.

    NaN becomes ``"nan"``, matching the existing ``str(acc).split(";")[0]``
    behaviour in ``enrich_common.py`` exactly. Callers that care about missing
    values must check for them before calling; this module does not invent a
    different missing-value convention than the one already shipping.
    """
    return str(acc)


def first_token(acc) -> str:
    """First semicolon token, stripped.

    Behaviour-identical to ``enrich_common.py``'s
    ``str(acc).split(";")[0].strip()``, including its treatment of NaN (which
    yields the string ``"nan"``) and of the empty string (which yields ``""``).
    """
    return _as_text(acc).split(";")[0].strip()


def split_group(acc) -> list[str]:
    """All semicolon tokens, each stripped of surrounding whitespace.

    A single accession yields a one-element list. Empty input yields ``[""]``,
    mirroring ``str.split``'s behaviour rather than special-casing it.
    """
    return [tok.strip() for tok in _as_text(acc).split(";")]


def is_valid_uniprot(tok: str) -> bool:
    """True if `tok` is a single, well-formed UniProt accession.

    Accepts an optional ``-N`` isoform suffix. Does NOT accept a semicolon-joined
    group -- pass one token at a time (see :func:`is_valid_group`).
    """
    return bool(UNIPROT_TOKEN_RE.match(_as_text(tok).strip()))


def is_valid_group(acc) -> bool:
    """True if EVERY token of `acc` is a valid UniProt accession.

    This is the validation policy: a protein group is valid only if all of its
    members are. Mirrors ``qc/schema.py``'s ``_accession_field_ok`` for non-NaN
    input; NaN is rejected here too.
    """
    if _is_missing(acc):
        return False
    return all(is_valid_uniprot(tok) for tok in split_group(acc))


def _is_missing(acc) -> bool:
    try:
        return bool(pd.isna(acc))
    except (TypeError, ValueError):
        return False


def strip_isoform(tok: str) -> str:
    """Remove a trailing ``-N`` isoform suffix from one token.

    Not applied anywhere by default -- see the module POLICY. Use only after an
    explicit, logged decision to collapse isoforms into canonical entries.
    """
    return ISOFORM_SUFFIX_RE.sub("", _as_text(tok).strip())


def has_isoform(tok: str) -> bool:
    """True if one token carries an isoform suffix."""
    return bool(ISOFORM_SUFFIX_RE.search(_as_text(tok).strip()))


def detect_isoforms(accs) -> list[str]:
    """Accession values (as given) that carry an isoform suffix on any token.

    Detection only -- nothing is stripped. Today's committed data yields an empty
    list; a non-empty result means the input changed and the isoform decision has
    to be made explicitly before those rows flow downstream.

    Order is preserved and duplicates are kept, so the result length is the
    number of affected rows.
    """
    found: list[str] = []
    for acc in accs:
        if _is_missing(acc):
            continue
        text = _as_text(acc)
        if any(has_isoform(tok) for tok in split_group(text)):
            found.append(text)
    return found


def is_junk_index_list(acc) -> bool:
    """True for a MaxQuant row-index list masquerading as an accession.

    The signature is: a semicolon-joined list whose tokens are ALL bare integers.
    Length is irrelevant -- ``P08752;P20612;Q9DC51;...`` is 69 characters and is a
    perfectly legitimate protein group, while the two junk rows in
    ``results/single_condition_proteins.csv`` are 32,759 and 681 characters of
    nothing but digits and semicolons.

    A single bare integer with no semicolon is NOT flagged: it is a malformed
    accession, and the accession validator should report it as such rather than
    have it quietly quarantined here.
    """
    if _is_missing(acc):
        return False
    text = _as_text(acc)
    if ";" not in text:
        return False
    tokens = split_group(text)
    return all(bool(_INTEGER_RE.match(tok)) for tok in tokens)


__all__ = [
    "UNIPROT_TOKEN_RE",
    "ISOFORM_SUFFIX_RE",
    "first_token",
    "split_group",
    "is_valid_uniprot",
    "is_valid_group",
    "strip_isoform",
    "has_isoform",
    "detect_isoforms",
    "is_junk_index_list",
]
