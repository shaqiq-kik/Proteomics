"""IPA upload files that also carry the 862 proteins the core export cannot see.

The problem this module exists to fix
-------------------------------------
``ipa_input.csv`` / ``ipa_input_full.csv`` are built from
``etl.foldchange_core.build_ipa_frame``, whose row selection is
``df[complete] & (regulated != "NO CHANGE")`` -- a ``complete == True`` gate
decided BEFORE limma ever runs. That gate is correct for what those files are,
and they are byte-frozen because ``ipa_input.csv`` has already been uploaded to
QIAGEN. But it means the 862 proteins DECISIONS_LOG D17 recovered are simply
absent from every file QIAGEN has ever seen:

* **248 tier-3 "partial missingness" proteins** (``regulated_up_partial.csv`` +
  ``regulated_down_partial.csv``) -- 1-2 of 4 raw replicates were zero/missing,
  so the raw classifier never ran, but limma DID test them after MinProb
  imputation and produced a real ``log2FC``/``p_value``/``adj_p_value``. FRZB
  P97401, GAS6 Q61592, LUM P51885 and SLIT3 Q9WVB4 are in this set -- the exact
  proteins a professor review reported missing.
* **614 tier-1/2 "fully undetected" proteins** (``qualitative_changes.csv``) --
  never identified at all on one side, so they carry no fold change and no
  p-value, deliberately (see below).

Re-running the pipeline does not fix this. ``build_ipa_frame`` still gates on
``complete``, so ``--all`` faithfully regenerates the same 715 rows. New export
code was the only honest option, and this module is it.

What it writes, and why the split
---------------------------------
``ipa_input_extended.csv`` / ``ipa_input_extended.txt``
    963 rows = the 715 core rows, values taken verbatim from
    ``ipa_input_full.csv``, plus the 248 tier-3 rows. A ``tier`` column
    (``core`` / ``partial``) is mandatory on every row: the two log2FC values
    are DIFFERENT QUANTITIES -- the core rows carry the raw mean-of-log2-ratios
    on unimputed data, the partial rows carry limma's estimate after MinProb
    imputation of 1-2 replicates -- and a reader who cannot tell them apart
    would be reading an imputed call as a complete-case measurement. Both are
    legitimate IPA input; conflating them is not. The ``.txt`` twin exists
    because QIAGEN "strongly recommend that the dataset files be in
    tab-delimited text format (.txt ...) for faster upload" (research1.md
    line 195); the ``.csv`` because it is the format the rest of this pipeline
    speaks.

``ipa_background_measured.txt``
    2552 bare accessions, one column: the **Reference Set** -- every protein
    this experiment actually measured, for IPA to test the foreground against
    instead of its full knowledge base. Without it a Core Analysis compares the
    963/372/242 uploads to every gene QIAGEN knows, which silently overstates
    enrichment for anything the assay was simply never able to see; a
    mass-spec run is not an unbiased sample of the proteome. This is the same
    correction the in-house ORA already applies via g:Profiler's custom
    ``domain_scope`` (DECISIONS_LOG D6), reusing the *same two source files*, so
    the two analyses cannot end up testing against two different universes --
    see :data:`BACKGROUND_SOURCES`.

``ipa_qualitative_up.txt`` / ``ipa_qualitative_down.txt``
    372 + 242 bare accessions, one column, no fold change and no p-value. That
    omission is the whole point, not an incomplete implementation: neither
    quantity is statistically valid when one entire condition has zero data
    points (tier 1) or zero variance (tier 2). QIAGEN's own contract makes this
    a first-class upload -- "the identifier column is the only strictly required
    element ... you can upload a list of molecules without [any measurement
    values]" (research1.md SECTION 3) -- so these are meant for a SEPARATE
    ID-only Core Analysis, not to be merged into the quantitative file. Merging
    them would require inventing a fold change, which this repo does not do.

Additive only
-------------
Following the precedent D17 set (``export/supplementary_lists.py``), nothing
here rewrites ``ipa_input.csv``, ``ipa_input_full.csv``, ``ipa_input_full.txt``,
``ipa_input_significant.csv``, ``regulated_up.csv``, ``regulated_down.csv``, or
the reference set's two sources ``foldchange_all.csv`` and
``single_condition_proteins.csv``. This module only READS them and writes new
filenames beside them -- see :data:`BYTE_FROZEN_INPUTS`, which is re-hashed
around every build so the promise is checked rather than asserted in prose.

The float-precision landmine, again
-----------------------------------
Every read here passes ``float_precision="round_trip"``, for the reason
``ipa_export.py`` documents at length: pandas' default parser is fast rather
than correctly rounded, and re-reading ``ipa_input_full.csv`` with it perturbs
log2FC values in the last ULP. The core rows of ``ipa_input_extended.csv`` must
be the *same text* as the file the user already uploaded, so that a reviewer
diffing the two sees 715 identical rows and 248 new ones. That property is
asserted, not hoped for -- see :func:`assert_core_rows_are_verbatim`, which
compares the raw CSV cells character by character rather than the parsed floats.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# proteomics_de/export/ipa_extended.py -> export -> proteomics_de -> repo root
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_ROOT = _PKG_DIR.parent

if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from proteomics_de.etl.accessions import (  # noqa: E402
    is_junk_index_list,
    is_valid_group,
)
from proteomics_de.export.ipa_export import (  # noqa: E402
    IPA_FULL_COLS,
    IPA_ID_COL,
    MEASUREMENT_COLS,
    REGULATED_LABELS,
    validate_ipa,
)

#: ``tier`` values. Not cosmetic -- see the module docstring: the two tiers'
#: ``log2FC`` columns are different quantities computed on different data.
TIER_CORE = "core"
TIER_PARTIAL = "partial"

#: Column order for ``ipa_input_extended.*``. Identifier leftmost, per QIAGEN's
#: contract (research1.md SECTION 3); ``IPA_FULL_COLS`` unchanged in the middle
#: so the 715 core rows are a literal prefix-projection of ``ipa_input_full.csv``;
#: ``tier`` appended last so a diff against the older file shows one added
#: column rather than a reshuffle.
IPA_EXTENDED_COLS = list(IPA_FULL_COLS) + ["tier"]

#: The tier-3 source files, and the directional label each one implies. The
#: label comes from WHICH FILE a row came from, never from re-thresholding its
#: log2FC here -- ``supplementary_lists.select_partial_regulated`` already made
#: that call at ``config.constants.LOG2_THRESHOLD`` and this module must not
#: quietly make a second, possibly different one.
PARTIAL_SOURCES = (
    ("regulated_up_partial.csv", "UP"),
    ("regulated_down_partial.csv", "DOWN"),
)

IPA_EXTENDED_CSV_FILENAME = "ipa_input_extended.csv"
IPA_EXTENDED_TXT_FILENAME = "ipa_input_extended.txt"
IPA_QUALITATIVE_UP_FILENAME = "ipa_qualitative_up.txt"
IPA_QUALITATIVE_DOWN_FILENAME = "ipa_qualitative_down.txt"

IPA_BACKGROUND_FILENAME = "ipa_background_measured.txt"

#: ``qualitative_changes.csv``'s direction label -> output filename.
QUALITATIVE_SOURCES = (
    ("UP", IPA_QUALITATIVE_UP_FILENAME, "n_ipa_qualitative_up"),
    ("DOWN", IPA_QUALITATIVE_DOWN_FILENAME, "n_ipa_qualitative_down"),
)

#: The three quantitative/qualitative uploads, with the frozen-count key each
#: one's size is pinned to. One list, used both by :func:`assert_partition` (the
#: three must partition 1577 proteins) and by
#: :func:`assert_background_covers_foreground` (the background must strictly
#: contain all three). Keeping it in one place is what stops a fourth upload
#: being added to one check and forgotten by the other.
FOREGROUND_UPLOADS = (
    (IPA_EXTENDED_CSV_FILENAME, "n_ipa_extended"),
    (IPA_QUALITATIVE_UP_FILENAME, "n_ipa_qualitative_up"),
    (IPA_QUALITATIVE_DOWN_FILENAME, "n_ipa_qualitative_down"),
)

#: ``(filename, accession column)`` of the two files that between them enumerate
#: every protein the experiment measured -- the Reference Set's definition.
#:
#: This is DELIBERATELY the same pair ``enrich/enrich_common.py``'s
#: ``load_background_and_queries`` unions to build the ORA ``domain_scope``
#: background (DECISIONS_LOG D6/D11): ``foldchange_all.csv`` is every protein
#: quantified in BOTH conditions, ``single_condition_proteins.csv`` every protein
#: quantified in exactly one, and nothing measured falls outside those two. The
#: two files are disjoint by accession, so the union is ``1948 + 604`` rows --
#: the ``background_union`` key of ``frozen_counts.json``, asserted below rather
#: than restated here.
#:
#: The ONE difference from the ORA background is the key. ORA projects each row
#: onto a single gene symbol (falling back to the accession's first token) and
#: dedupes to 2530 unique symbols, because g:Profiler takes one identifier per
#: query. This file is keyed on the accession instead -- IPA's identifier type is
#: set to UniProt/Swiss-Prot Accession for the other uploads, and a background in
#: a different namespace than the foreground would not match. See
#: :func:`build_measured_background` for why the WHOLE group string is kept.
BACKGROUND_SOURCES = (
    ("foldchange_all.csv", "UniProt Accession Number"),
    ("single_condition_proteins.csv", "accession"),
)

#: Files this module reads and must never write. Enforced by
#: :func:`assert_frozen_inputs_untouched` around every build, because "additive
#: only" is a promise about bytes and deserves a runtime check rather than a
#: comment (DECISIONS_LOG D17; ``ipa_input.csv`` is already at QIAGEN).
BYTE_FROZEN_INPUTS = (
    "ipa_input.csv",
    "ipa_input_full.csv",
    "ipa_input_full.txt",
    "ipa_input_significant.csv",
    "regulated_up.csv",
    "regulated_down.csv",
    # The two Reference Set sources. Read by build_measured_background; they are
    # the fold-change stage's frozen outputs and the ORA background's inputs, so
    # an accidental write here would silently move the enrichment universe too.
    "foldchange_all.csv",
    "single_condition_proteins.csv",
)


class IpaExtendedError(AssertionError):
    """An extended IPA output violates its contract.

    Subclasses ``AssertionError`` on purpose, matching
    :class:`proteomics_de.export.ipa_export.IpaValidationError` and
    :class:`proteomics_de.export.supplementary_lists.SupplementaryListsError`:
    these checks replace bare ``assert`` statements, and callers (and pytest)
    should be able to treat them the same way.
    """


def _frozen_counts(path=None):
    path = Path(path) if path else _PKG_DIR / "tests" / "expected" / "frozen_counts.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _read(path):
    """Read a committed CSV at maximal float precision. No exceptions to this."""
    return pd.read_csv(path, float_precision="round_trip")


# ---------------------------------------------------------------------------
# ipa_input_extended.{csv,txt}
# ---------------------------------------------------------------------------
def build_extended_frame(df_full, partial_frames):
    """Stack the 715 core rows and the 248 tier-3 rows into one IPA frame.

    Row order is the deliverable, so it is spelled out rather than left to a
    sort: **core rows first, in exactly the order ``ipa_input_full.csv`` has
    them**, then the partial rows by descending ``|log2FC|``. Keeping the core
    block untouched means a reviewer can diff the new file against the uploaded
    one and see 715 unchanged lines followed by 248 additions, instead of an
    interleaving that looks like a wholesale rewrite. The sort of the partial
    block uses ``kind="mergesort"`` for the same reason
    ``regulated_lists.sort_by_magnitude`` does: the byte-freeze gate
    (DECISIONS_LOG D14) requires a re-run to reproduce identical bytes, and the
    default quicksort is not stable across ties.

    Parameters
    ----------
    df_full :
        ``ipa_input_full.csv`` as a frame; needs :data:`IPA_FULL_COLS`.
    partial_frames :
        Iterable of ``(frame, label)`` pairs -- typically the two files named in
        :data:`PARTIAL_SOURCES`. `label` (``UP``/``DOWN``) becomes the row's
        ``regulated`` value; see :data:`PARTIAL_SOURCES` for why it is not
        re-derived from ``log2FC`` here.

    Returns
    -------
    pandas.DataFrame
        :data:`IPA_EXTENDED_COLS`, ``RangeIndex``.
    """
    missing = [c for c in IPA_FULL_COLS if c not in df_full.columns]
    if missing:
        raise IpaExtendedError(f"ipa_input_full frame is missing columns: {missing}")

    core = df_full[list(IPA_FULL_COLS)].copy()
    core["tier"] = TIER_CORE

    blocks = []
    for frame, label in partial_frames:
        if label not in REGULATED_LABELS:
            raise IpaExtendedError(
                f"partial block label {label!r} is not one of {list(REGULATED_LABELS)}"
            )
        needed = [IPA_ID_COL, "Gene names", "log2FC", "p_value", "adj_p_value"]
        absent = [c for c in needed if c not in frame.columns]
        if absent:
            raise IpaExtendedError(f"partial frame ({label}) is missing columns: {absent}")
        block = frame[needed].copy()
        block["regulated"] = label
        block["tier"] = TIER_PARTIAL
        blocks.append(block[IPA_EXTENDED_COLS])

    if not blocks:
        raise IpaExtendedError("no partial frames supplied; the extended file would be "
                               "identical to ipa_input_full.csv and pointless")

    partial = pd.concat(blocks, ignore_index=True)
    partial = (
        partial.assign(_abs=partial["log2FC"].abs())
        .sort_values("_abs", ascending=False, kind="mergesort")
        .drop(columns="_abs")
        .reset_index(drop=True)
    )

    extended = pd.concat(
        [core[IPA_EXTENDED_COLS], partial[IPA_EXTENDED_COLS]], ignore_index=True
    )

    # Structural invariants of the frame itself. Checked here, before anything
    # reaches disk, because a duplicate accession or a NaN measurement is a
    # SILENT corruption at the QIAGEN end -- IPA drops or misreads the row
    # rather than erroring (ipa_export.py, MEASUREMENT_COLS).
    if not extended[IPA_ID_COL].is_unique:
        dupes = extended.loc[extended[IPA_ID_COL].duplicated(), IPA_ID_COL].tolist()
        raise IpaExtendedError(
            "duplicate accessions in the extended frame -- the core and tier-3 "
            f"sets are supposed to be disjoint by construction: {dupes[:10]}"
        )
    for col in MEASUREMENT_COLS:
        series = pd.to_numeric(extended[col], errors="coerce")
        n_bad = int((~np.isfinite(series.to_numpy(dtype=float))).sum())
        if n_bad:
            raise IpaExtendedError(
                f"{n_bad} non-finite value(s) in {col!r}. IPA does not error on "
                "these; it silently drops or misreads the row."
            )
    labels = set(extended["regulated"].unique())
    if not labels <= set(REGULATED_LABELS):
        raise IpaExtendedError(
            f"unexpected 'regulated' values: {sorted(labels - set(REGULATED_LABELS))}"
        )
    return extended[IPA_EXTENDED_COLS]


def write_ipa_extended(df_ext, csv_path, txt_path=None):
    """Write ``ipa_input_extended.csv`` and its tab-delimited ``.txt`` twin.

    Mirrors :func:`proteomics_de.export.ipa_export.write_ipa_full`: identifier
    leftmost, UTF-8 **without** a BOM (``"utf-8"``, never ``"utf-8-sig"`` -- the
    latter makes IPA read ``\\ufeffUniProt Accession Number`` as the first column
    name), a single header row. ``txt_path=False`` writes only the CSV;
    ``None`` derives the ``.txt`` sibling from `csv_path`.

    Returns
    -------
    list[pathlib.Path]
        The paths written, CSV first.
    """
    csv_path = Path(csv_path)
    if txt_path is None:
        txt_path = csv_path.with_suffix(".txt")

    written = [csv_path]
    df_ext[IPA_EXTENDED_COLS].to_csv(csv_path, index=False, encoding="utf-8")
    if txt_path is not False:
        txt_path = Path(txt_path)
        df_ext[IPA_EXTENDED_COLS].to_csv(txt_path, index=False, sep="\t", encoding="utf-8")
        written.append(txt_path)
    return written


# ---------------------------------------------------------------------------
# ipa_qualitative_{up,down}.txt
# ---------------------------------------------------------------------------
def select_qualitative_ids(qualitative, direction):
    """The accessions of ``qualitative_changes.csv`` rows with `direction`.

    Row ORDER is inherited from the source file, never re-sorted here.
    ``supplementary_lists.build_qualitative_changes`` already sorts
    deterministically (UP before DOWN, then case-folded gene name, mergesort),
    so inheriting it keeps these files byte-reproducible without this module
    owning a second, divergent ordering rule.

    Returns
    -------
    pandas.Series
        Accessions, ``RangeIndex``.
    """
    for col in (IPA_ID_COL, "direction"):
        if col not in qualitative.columns:
            raise IpaExtendedError(f"qualitative frame is missing column {col!r}")

    ids = qualitative.loc[qualitative["direction"] == direction, IPA_ID_COL]
    ids = ids.reset_index(drop=True)
    if ids.isna().any():
        raise IpaExtendedError(f"{direction}: blank accession in qualitative_changes.csv")
    if not ids.is_unique:
        dupes = ids[ids.duplicated()].tolist()
        raise IpaExtendedError(f"{direction}: duplicate accessions: {dupes[:10]}")
    return ids


def write_id_list(ids, path, header=IPA_ID_COL):
    """Write a one-column, header + N-line identifier file for IPA.

    Deliberately the whole file format: an ID-only upload is a first-class IPA
    input (research1.md SECTION 3, quoting QIAGEN: *"you can upload a list of
    molecules without [any measurement values]"*), and adding a placeholder
    fold-change column would be fabricating the number this tier does not have.

    ``.txt`` because that is QIAGEN's recommended upload extension; with one
    column there is no delimiter to choose, so the file is simultaneously valid
    tab- and comma-delimited.
    """
    path = Path(path)
    pd.Series(list(ids), name=header).to_frame().to_csv(
        path, index=False, encoding="utf-8"
    )
    return path


# ---------------------------------------------------------------------------
# ipa_background_measured.txt -- the IPA Reference Set
# ---------------------------------------------------------------------------
def build_measured_background(source_frames):
    """Every accession the experiment measured, deduplicated and sorted.

    `source_frames` is an iterable of ``(frame, column, label)`` -- typically the
    two files named in :data:`BACKGROUND_SOURCES`, whose docstring gives the
    derivation. `label` is only used in error messages.

    The WHOLE accession string is kept, semicolons and all
    -------------------------------------------------------
    95 of the 1577 foreground accessions are protein groups like
    ``P05132;P68181``, and ``etl/accessions.py``'s stated policy is that merging
    is on the whole string -- "a protein group is an identity, not a set". So the
    foreground files carry group strings verbatim, and a background built from
    ``first_token`` would be a set of 2531 single accessions that does NOT
    contain those 95 rows. It would look plausible, upload cleanly, and quietly
    exclude 95 measured proteins from the reference set. Keeping the whole string
    makes the background live in the same namespace as the foreground, which is
    what lets :func:`assert_background_covers_foreground` be a real check rather
    than a formality.

    Sorted ascending by the accession string, purely for determinism: this file
    is a SET, its order carries no meaning, and the byte-freeze gate
    (DECISIONS_LOG D14) requires a re-run to reproduce identical bytes.

    Returns
    -------
    pandas.Series
        Accessions, ``RangeIndex``, named :data:`IPA_ID_COL`.
    """
    seen = {}
    for frame, column, label in source_frames:
        if column not in frame.columns:
            raise IpaExtendedError(f"{label} is missing its accession column {column!r}")
        values = frame[column]
        if values.isna().any():
            raise IpaExtendedError(
                f"{label}: {int(values.isna().sum())} blank accession(s); a "
                "measured protein with no identifier cannot go in the reference set"
            )
        ids = values.astype(str).str.strip()
        if (ids == "").any():
            raise IpaExtendedError(f"{label}: empty-string accession(s)")
        if not ids.is_unique:
            dupes = ids[ids.duplicated()].tolist()
            raise IpaExtendedError(f"{label}: duplicate accessions: {dupes[:10]}")

        # The two sources are disjoint by construction -- a protein is quantified
        # in both conditions or in exactly one, never both -- and that is what
        # makes background_union a plain 1948 + 604 sum. An overlap means one of
        # the two files was rebuilt from the wrong mask, so it is an error here
        # rather than something the set union quietly absorbs.
        overlap = sorted(set(ids) & set(seen))
        if overlap:
            raise IpaExtendedError(
                f"{label} shares {len(overlap)} accession(s) with "
                f"{seen[overlap[0]]}; the background sources must be disjoint: "
                f"{overlap[:10]}"
            )
        for acc in ids:
            seen[acc] = label

    # Junk-index quarantine, re-checked at the point of upload. DECISIONS_LOG D11
    # found two rows of single_condition_proteins.csv carrying a ';'-joined list
    # of bare MaxQuant row indices where an accession belongs; they are dropped
    # upstream now, but this file is the one place where such a row would become
    # a reference-set member that IPA cannot map and nobody would notice.
    bad = sorted(acc for acc in seen if not is_valid_group(acc))
    if bad:
        junk = [acc for acc in bad if is_junk_index_list(acc)]
        raise IpaExtendedError(
            f"{len(bad)} accession(s) are not well-formed UniProt groups "
            f"({len(junk)} look like MaxQuant row-index lists, cf. D11): "
            f"{[acc[:40] for acc in bad[:5]]}"
        )

    return pd.Series(sorted(seen), name=IPA_ID_COL)


# ---------------------------------------------------------------------------
# Cross-file validation
# ---------------------------------------------------------------------------
def _measurement_cells(path, id_col=IPA_ID_COL):
    """``{accession: {column: raw_cell_text}}`` for the measurement columns.

    Reads the file as TEXT, not as floats. Comparing parsed floats would prove
    the two files agree numerically; comparing the cells proves the bytes agree,
    which is the stronger statement and the one the module docstring promises.
    """
    delim = "\t" if Path(path).suffix.lower() in (".txt", ".tsv") else ","
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter=delim))
    header, data = rows[0], rows[1:]
    idx = {c: header.index(c) for c in MEASUREMENT_COLS if c in header}
    key = header.index(id_col)
    return {r[key]: {c: r[i] for c, i in idx.items()} for r in data if r}


def assert_core_rows_are_verbatim(extended_csv, full_csv):
    """The 715 core rows must be textually identical to ``ipa_input_full.csv``.

    Not a nicety. ``ipa_input.csv`` is already uploaded to QIAGEN, and the whole
    argument for shipping a second file rather than editing the first is that
    the overlap is unchanged. If ``float_precision="round_trip"`` were ever
    dropped from a read in this module, 87 of the 715 ``log2FC`` values would
    move in the last ULP and that argument would quietly become false -- with no
    other symptom. This check is what makes the regression loud.
    """
    ext_cells = _measurement_cells(extended_csv)
    full_cells = _measurement_cells(full_csv)

    missing = sorted(set(full_cells) - set(ext_cells))
    if missing:
        raise IpaExtendedError(
            f"{len(missing)} accession(s) from {Path(full_csv).name} are absent "
            f"from {Path(extended_csv).name}: {missing[:10]}"
        )

    drifted = [
        (acc, col, full_cells[acc][col], ext_cells[acc][col])
        for acc in full_cells
        for col in full_cells[acc]
        if full_cells[acc][col] != ext_cells[acc][col]
    ]
    if drifted:
        raise IpaExtendedError(
            f"{len(drifted)} core measurement cell(s) differ textually from "
            f"{Path(full_csv).name} -- a float_precision regression, not a data "
            f"change. First: {drifted[:5]}"
        )


def assert_twins_agree(csv_path, txt_path):
    """The ``.csv`` and its tab-delimited ``.txt`` must carry identical data.

    Two files go to QIAGEN and only one of them tends to get eyeballed. Parsing
    both and comparing cell-by-cell is the only way to know the unexamined one
    is not a stale copy from an earlier run.
    """
    def cells(path):
        delim = "\t" if Path(path).suffix.lower() in (".txt", ".tsv") else ","
        with open(path, newline="", encoding="utf-8") as fh:
            return list(csv.reader(fh, delimiter=delim))

    a, b = cells(csv_path), cells(txt_path)
    if a != b:
        n = sum(1 for x, y in zip(a, b) if x != y)
        raise IpaExtendedError(
            f"{Path(csv_path).name} and {Path(txt_path).name} disagree: "
            f"{len(a)} vs {len(b)} lines, {n} differing rows"
        )


def assert_partition(results_dir, counts_path=None):
    """The three new uploads must be pairwise disjoint and cover 1577 proteins.

    This is the check that says the deliverable is *complete*: every protein
    D17 recovered reaches QIAGEN exactly once, in exactly the file whose
    statistical contract it satisfies. 963 + 372 + 242 = 1577 distinct
    accessions, and the 1577 is DERIVED from the three frozen counts rather than
    typed, so it cannot drift out of agreement with them.

    Returns
    -------
    dict
        ``{"n_extended", "n_up", "n_down", "n_union"}``.
    """
    results_dir = Path(results_dir)
    counts = _frozen_counts(counts_path)

    ext = set(_read(results_dir / IPA_EXTENDED_CSV_FILENAME)[IPA_ID_COL])
    up = set(_read(results_dir / IPA_QUALITATIVE_UP_FILENAME)[IPA_ID_COL])
    down = set(_read(results_dir / IPA_QUALITATIVE_DOWN_FILENAME)[IPA_ID_COL])

    by_name = {
        IPA_EXTENDED_CSV_FILENAME: ext,
        IPA_QUALITATIVE_UP_FILENAME: up,
        IPA_QUALITATIVE_DOWN_FILENAME: down,
    }
    for label, key in FOREGROUND_UPLOADS:
        ids = by_name[label]
        if len(ids) != counts[key]:
            raise IpaExtendedError(
                f"{label}: {len(ids)} distinct accessions, expected "
                f"{counts[key]} (frozen_counts {key})"
            )

    for (an, a), (bn, b) in (
        ((IPA_EXTENDED_CSV_FILENAME, ext), (IPA_QUALITATIVE_UP_FILENAME, up)),
        ((IPA_EXTENDED_CSV_FILENAME, ext), (IPA_QUALITATIVE_DOWN_FILENAME, down)),
        ((IPA_QUALITATIVE_UP_FILENAME, up), (IPA_QUALITATIVE_DOWN_FILENAME, down)),
    ):
        overlap = a & b
        if overlap:
            raise IpaExtendedError(
                f"{an} and {bn} share {len(overlap)} accession(s); the three "
                f"uploads must partition the protein set: {sorted(overlap)[:10]}"
            )

    union = ext | up | down
    expected_union = (
        counts["n_ipa_extended"]
        + counts["n_ipa_qualitative_up"]
        + counts["n_ipa_qualitative_down"]
    )
    if len(union) != expected_union:
        raise IpaExtendedError(
            f"the three uploads cover {len(union)} distinct accessions, expected "
            f"{expected_union}"
        )
    return {"n_extended": len(ext), "n_up": len(up), "n_down": len(down),
            "n_union": len(union)}


def assert_background_covers_foreground(results_dir, counts_path=None):
    """The Reference Set must be a STRICT superset of all three uploads.

    A background that does not contain its own foreground is not a background:
    IPA would be told that some tested molecules are outside the universe they
    were drawn from, and every enrichment p-value computed from that pair is
    arithmetic on an incoherent contingency table. So this is an error, not a
    warning -- and *strict* is asserted too, because a background exactly equal
    to the foreground would mean 975 measured-but-unregulated proteins had
    silently gone missing from the file, leaving a "background" that cannot
    discriminate anything.

    Returns
    -------
    dict
        ``{"n_background", "n_foreground", "n_background_only"}``.
    """
    results_dir = Path(results_dir)
    counts = _frozen_counts(counts_path)

    background = set(_read(results_dir / IPA_BACKGROUND_FILENAME)[IPA_ID_COL])
    expected = int(counts["n_ipa_background_measured"])
    if len(background) != expected:
        raise IpaExtendedError(
            f"{IPA_BACKGROUND_FILENAME}: {len(background)} distinct accessions, "
            f"expected {expected} (frozen_counts n_ipa_background_measured)"
        )

    foreground = set()
    for label, _key in FOREGROUND_UPLOADS:
        ids = set(_read(results_dir / label)[IPA_ID_COL])
        orphans = sorted(ids - background)
        if orphans:
            raise IpaExtendedError(
                f"{len(orphans)} accession(s) in {label} are absent from "
                f"{IPA_BACKGROUND_FILENAME}. The reference set must contain "
                f"every uploaded molecule: {orphans[:10]}"
            )
        foreground |= ids

    if not background > foreground:
        raise IpaExtendedError(
            f"{IPA_BACKGROUND_FILENAME} ({len(background)}) is not a STRICT "
            f"superset of the {len(foreground)} foreground accessions -- the "
            "measured-but-unregulated proteins are what make it a background"
        )
    return {
        "n_background": len(background),
        "n_foreground": len(foreground),
        "n_background_only": len(background - foreground),
    }


def assert_frozen_inputs_untouched(results_dir, before):
    """Re-hash :data:`BYTE_FROZEN_INPUTS` and refuse to have moved them.

    `before` is the mapping :func:`_hash_frozen_inputs` returned before the
    build. Cheap insurance on the one promise this module cannot be allowed to
    break.
    """
    after = _hash_frozen_inputs(results_dir)
    moved = sorted(k for k in before if before[k] != after.get(k))
    if moved:
        raise IpaExtendedError(
            "this module modified a byte-frozen input, which it must never do: "
            f"{moved}"
        )


def _hash_frozen_inputs(results_dir):
    import hashlib

    results_dir = Path(results_dir)
    out = {}
    for name in BYTE_FROZEN_INPUTS:
        p = results_dir / name
        if p.exists():
            out[name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_extended_from_results(results_dir=None, *, validate=True, counts_path=None):
    """Build all five new IPA files from the committed results directory.

    Reads ``ipa_input_full.csv``, ``regulated_up_partial.csv``,
    ``regulated_down_partial.csv``, ``qualitative_changes.csv`` and the two
    :data:`BACKGROUND_SOURCES` -- every number is lifted from a file that already
    exists; nothing is recomputed, imputed or invented here.

    Returns
    -------
    list[pathlib.Path]
        ``[extended.csv, extended.txt, qualitative_up.txt, qualitative_down.txt,
        background_measured.txt]``.
    """
    results_dir = Path(results_dir) if results_dir else _PKG_DIR / "results"

    full_path = results_dir / "ipa_input_full.csv"
    qual_path = results_dir / "qualitative_changes.csv"
    partial_paths = [(results_dir / name, label) for name, label in PARTIAL_SOURCES]
    background_paths = [(results_dir / name, col) for name, col in BACKGROUND_SOURCES]
    for p in ([full_path, qual_path]
              + [pp for pp, _ in partial_paths]
              + [bp for bp, _ in background_paths]):
        if not p.exists():
            raise FileNotFoundError(
                f"{p} not found -- run run_pipeline.py --only "
                "foldchange,ipa_export,regulated_lists_supplementary first."
            )

    before = _hash_frozen_inputs(results_dir)

    df_full = _read(full_path)
    partial_frames = [(_read(p), label) for p, label in partial_paths]
    extended = build_extended_frame(df_full, partial_frames)

    csv_path = results_dir / IPA_EXTENDED_CSV_FILENAME
    txt_path = results_dir / IPA_EXTENDED_TXT_FILENAME
    written = write_ipa_extended(extended, csv_path, txt_path)

    qualitative = _read(qual_path)
    for direction, filename, _key in QUALITATIVE_SOURCES:
        ids = select_qualitative_ids(qualitative, direction)
        written.append(write_id_list(ids, results_dir / filename))

    background_frames = [(_read(p), col, p.name) for p, col in background_paths]
    background = build_measured_background(background_frames)
    background_path = results_dir / IPA_BACKGROUND_FILENAME
    written.append(write_id_list(background, background_path))

    if validate:
        counts = _frozen_counts(counts_path)
        n_extended = int(counts["n_ipa_extended"])

        # Reuse ipa_export's file-contract validator verbatim: no BOM, decodes as
        # UTF-8, exactly one header row with uniform field counts, identifier
        # leftmost and free of stray characters, no NA/inf in any measurement
        # column, expected row count. Passing `required_columns` explicitly
        # because that function infers 4- vs 6-column layouts and knows nothing
        # about the 7-column extended one.
        for p in (csv_path, txt_path):
            validate_ipa(p, expected_rows=n_extended, required_columns=IPA_EXTENDED_COLS)
        for _direction, filename, key in QUALITATIVE_SOURCES:
            validate_ipa(
                results_dir / filename,
                expected_rows=int(counts[key]),
                required_columns=[IPA_ID_COL],
                measurement_columns=[],
            )
        validate_ipa(
            background_path,
            expected_rows=int(counts["n_ipa_background_measured"]),
            required_columns=[IPA_ID_COL],
            measurement_columns=[],
        )

        # The Reference Set's size is the two sources' row counts, and that has
        # to be the SAME number the ORA background is built from -- the whole
        # point of reusing BACKGROUND_SOURCES is that IPA and g:Profiler test
        # against one universe. Checked against background_union rather than
        # only against its own key, so the two cannot drift apart silently.
        n_background = len(background)
        n_source_rows = sum(len(f) for f, _col, _name in background_frames)
        if n_background != n_source_rows:
            raise IpaExtendedError(
                f"{n_background} background accessions from {n_source_rows} source "
                "rows -- the two sources are supposed to be disjoint and unique"
            )
        if n_background != int(counts["background_union"]):
            raise IpaExtendedError(
                f"{n_background} background accessions, but frozen_counts "
                f"background_union (the ORA domain_scope background) is "
                f"{counts['background_union']}; IPA and g:Profiler would be "
                "testing against different universes"
            )

        n_core = len(df_full)
        n_partial = sum(len(f) for f, _ in partial_frames)
        if n_core + n_partial != n_extended:
            raise IpaExtendedError(
                f"{n_core} core + {n_partial} partial = {n_core + n_partial} rows, "
                f"expected {n_extended} (frozen_counts n_ipa_extended)"
            )
        tiers = extended["tier"].value_counts().to_dict()
        if tiers.get(TIER_CORE) != n_core or tiers.get(TIER_PARTIAL) != n_partial:
            raise IpaExtendedError(f"tier split is {tiers}, expected "
                                   f"{{{TIER_CORE!r}: {n_core}, {TIER_PARTIAL!r}: {n_partial}}}")

        assert_core_rows_are_verbatim(csv_path, full_path)
        assert_twins_agree(csv_path, txt_path)
        assert_partition(results_dir, counts_path)
        assert_background_covers_foreground(results_dir, counts_path)
        assert_frozen_inputs_untouched(results_dir, before)

    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", default=None,
                    help="directory holding ipa_input_full.csv / "
                         "regulated_{up,down}_partial.csv / qualitative_changes.csv "
                         "/ foldchange_all.csv / single_condition_proteins.csv "
                         "(default: proteomics_de/results)")
    ap.add_argument("--no-sidecars", action="store_true",
                    help="skip writing the .provenance.json sidecars")
    args = ap.parse_args(argv)

    results_dir = Path(args.results_dir) if args.results_dir else _PKG_DIR / "results"

    written = build_extended_from_results(results_dir)

    # Re-validate from disk and report. `expected_rows` is passed explicitly for
    # every file: ipa_export.validate_ipa defaults it to the 715-row
    # complete&UP|DOWN count, which is precisely the number these files exist to
    # stop being the whole story.
    counts = _frozen_counts()
    contracts = {
        IPA_EXTENDED_CSV_FILENAME: (counts["n_ipa_extended"], IPA_EXTENDED_COLS, None),
        IPA_EXTENDED_TXT_FILENAME: (counts["n_ipa_extended"], IPA_EXTENDED_COLS, None),
        IPA_QUALITATIVE_UP_FILENAME: (counts["n_ipa_qualitative_up"], [IPA_ID_COL], []),
        IPA_QUALITATIVE_DOWN_FILENAME: (counts["n_ipa_qualitative_down"], [IPA_ID_COL], []),
        IPA_BACKGROUND_FILENAME: (counts["n_ipa_background_measured"], [IPA_ID_COL], []),
    }
    for p in written:
        expected, columns, measurements = contracts[p.name]
        rep = validate_ipa(p, expected_rows=expected, required_columns=columns,
                           measurement_columns=measurements)
        print(f"Saved {p} ({rep['n_rows']} rows, {len(rep['columns'])} columns, "
              f"{len(rep['checks'])} checks passed)")

    summary = assert_partition(results_dir)
    print(
        f"Partition OK: {summary['n_extended']} quantitative + {summary['n_up']} "
        f"qualitative UP + {summary['n_down']} qualitative DOWN = "
        f"{summary['n_union']} distinct accessions, pairwise disjoint"
    )

    coverage = assert_background_covers_foreground(results_dir)
    print(
        f"Reference Set OK: {coverage['n_background']} measured accessions "
        f"strictly contain all {coverage['n_foreground']} foreground accessions "
        f"({coverage['n_background_only']} measured but not uploaded)"
    )

    if not args.no_sidecars:
        from proteomics_de.provenance import sidecar

        for p in written:
            print(f"Saved {sidecar(p)}")

    return 0


__all__ = [
    "TIER_CORE",
    "TIER_PARTIAL",
    "IPA_EXTENDED_COLS",
    "PARTIAL_SOURCES",
    "QUALITATIVE_SOURCES",
    "FOREGROUND_UPLOADS",
    "BACKGROUND_SOURCES",
    "IPA_EXTENDED_CSV_FILENAME",
    "IPA_EXTENDED_TXT_FILENAME",
    "IPA_QUALITATIVE_UP_FILENAME",
    "IPA_QUALITATIVE_DOWN_FILENAME",
    "IPA_BACKGROUND_FILENAME",
    "BYTE_FROZEN_INPUTS",
    "IpaExtendedError",
    "build_extended_frame",
    "write_ipa_extended",
    "select_qualitative_ids",
    "write_id_list",
    "build_measured_background",
    "assert_core_rows_are_verbatim",
    "assert_twins_agree",
    "assert_partition",
    "assert_background_covers_foreground",
    "assert_frozen_inputs_untouched",
    "build_extended_from_results",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
