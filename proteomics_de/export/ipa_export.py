"""IPA export writers and the IPA file-contract validator.

What IPA actually needs (research1.md SECTION 3, quoting QIAGEN's own docs):

* the **identifier column is the only strictly required element**, and it must be
  leftmost;
* fold change alone runs every IPA module (all z-scores are driven by the
  *direction* of the fold change), but shipping the FDR alongside it is what
  "unlocks **cutoff-based filtering**" -- IPA's Standard/Bonferroni/FDR expression
  cutoffs (research1.md line 210: *"ship fold change + adj.P.Val once limma is in
  place"*);
* QIAGEN "strongly recommend that the dataset files be in tab-delimited text
  format (.txt, rather than Excel-based formats) for faster upload"
  (research1.md line 195);
* UTF-8 with **no BOM**, a single header row, no NA/inf in any measurement column.

So this module emits three things, deliberately kept apart:

``ipa_input.csv``
    The historical four-column file (accession, gene, log2FC, regulated). It is
    a **frozen output** -- the user may already have uploaded it to QIAGEN, so
    its bytes must never move. :func:`write_ipa` is therefore left exactly as it
    was and is the only writer that touches it.

``ipa_input_full.csv`` / ``ipa_input_full.txt``
    New. The same 715 rows plus ``p_value`` and ``adj_p_value`` joined from
    ``qc_limma.csv``, in both CSV and QIAGEN's preferred tab-delimited form.
    This is the file to upload when you want expression cutoffs.

``ipa_input_significant.csv``
    Written by ``limma_test.py``, not here, and **header-only by design**: 0 of
    1938 proteins pass FDR < 0.05 (DECISIONS_LOG D2). That emptiness is the
    headline finding, not a bug, so :func:`validate_ipa_significant` asserts it
    as a contract rather than trying to repair it.

A float-precision landmine, learned the hard way
------------------------------------------------
``pandas.read_csv``'s default float parser (``float_precision="high"``) is not
correctly rounded: re-reading ``ipa_input.csv`` and writing it straight back out
perturbs **87 of the 715** ``log2FC`` values in the last ULP. That would leave
``ipa_input_full.csv`` textually disagreeing with the ``ipa_input.csv`` the user
already uploaded, for rows that are numerically identical. Every read in this
module therefore passes ``float_precision="round_trip"``, which was verified to
reproduce ``ipa_input.csv`` byte for byte. Do not remove it.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import pandas as pd

# proteomics_de/export/ipa_export.py -> export -> proteomics_de -> repo root
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_ROOT = _PKG_DIR.parent

#: Leftmost column in every IPA file. QIAGEN keys the upload off this.
IPA_ID_COL = "UniProt Accession Number"

#: The frozen four-column layout of ``ipa_input.csv``.
IPA_COLS = [IPA_ID_COL, "Gene names", "log2FC", "regulated"]

#: ``ipa_input_full.*`` -- the same rows with limma's raw p and BH-adjusted p.
IPA_FULL_COLS = IPA_COLS + ["p_value", "adj_p_value"]

#: ``ipa_input_significant.csv`` (written by limma_test.py; validated here).
IPA_SIGNIFICANT_COLS = IPA_COLS + ["adj_p_value"]

#: Columns IPA treats as "measurement value types" (research1.md line 193). These
#: are the ones that must be finite -- an NA or inf here is a silent upload
#: corruption, because IPA drops or misreads the row rather than erroring.
MEASUREMENT_COLS = ("log2FC", "p_value", "adj_p_value")

#: The two directional labels that make a protein IPA-eligible.
REGULATED_LABELS = ("UP", "DOWN")

_UTF8_BOM = b"\xef\xbb\xbf"

#: Delimiter per extension. QIAGEN accepts ``.csv`` and ``.txt`` (tab-delimited).
_DELIMITERS = {".csv": ",", ".txt": "\t", ".tsv": "\t"}


class IpaValidationError(AssertionError):
    """An IPA output violates the upload contract.

    Subclasses ``AssertionError`` on purpose: these checks replace bare
    ``assert`` statements in the pipeline, and callers (and pytest) should be
    able to treat them the same way.
    """


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def write_ipa(df, path, columns):
    """Write `df`'s `columns` to `path` as the IPA input CSV.

    **Byte-frozen.** This is literally the original ``foldchange.py:171``
    behaviour (``df[cols].to_csv(path, index=False, encoding="utf-8")``) and it
    must stay that way: ``results/ipa_input.csv`` is a committed output that the
    user may already have uploaded to QIAGEN, so its sha256 is part of the
    freeze gate. Anything new goes in :func:`write_ipa_full` instead.

    Note ``encoding="utf-8"`` and not ``"utf-8-sig"`` -- the latter would prepend
    a BOM, which IPA reads as part of the first column name.

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


def build_ipa_full(df_ipa, limma, *, id_col=IPA_ID_COL, limma_id_col="id"):
    """Join limma's p-values onto the IPA frame.

    The join is a LEFT join on accession and every row must find a match. That
    is not an optimistic assumption, it is a checked one: an IPA row is by
    construction ``complete & regulated in {UP, DOWN}``, and limma tests every
    complete row (only the ON_OFF proteins are excluded, and those never carry a
    directional label). So a NaN p-value here means the eligibility invariant
    broke somewhere upstream -- e.g. an ON_OFF protein leaked into the IPA set,
    or limma silently dropped rows -- and we raise rather than ship a file with
    blank measurement cells that IPA would quietly ignore.

    Parameters
    ----------
    df_ipa :
        Frame with :data:`IPA_COLS`.
    limma :
        ``qc_limma.csv`` as a frame; needs `limma_id_col`, ``p_value`` and
        ``adj_p_value``.

    Returns
    -------
    pandas.DataFrame
        :data:`IPA_FULL_COLS`, same rows and same order as `df_ipa`.
    """
    missing = [c for c in IPA_COLS if c not in df_ipa.columns]
    if missing:
        raise IpaValidationError(f"IPA frame is missing columns: {missing}")
    for col in (limma_id_col, "p_value", "adj_p_value"):
        if col not in limma.columns:
            raise IpaValidationError(f"limma frame is missing column {col!r}")

    if not df_ipa[id_col].is_unique:
        dupes = df_ipa.loc[df_ipa[id_col].duplicated(), id_col].tolist()
        raise IpaValidationError(f"duplicate accessions in the IPA frame: {dupes[:10]}")
    if not limma[limma_id_col].is_unique:
        dupes = limma.loc[limma[limma_id_col].duplicated(), limma_id_col].tolist()
        raise IpaValidationError(f"duplicate accessions in qc_limma: {dupes[:10]}")

    n_before = len(df_ipa)
    merged = df_ipa[IPA_COLS].merge(
        limma[[limma_id_col, "p_value", "adj_p_value"]],
        how="left",
        left_on=id_col,
        right_on=limma_id_col,
        validate="one_to_one",
    )
    if len(merged) != n_before:
        raise IpaValidationError(
            f"the limma join changed the row count: {n_before} -> {len(merged)}"
        )

    for col in ("p_value", "adj_p_value"):
        n_nan = int(merged[col].isna().sum())
        if n_nan:
            lost = merged.loc[merged[col].isna(), id_col].tolist()[:10]
            raise IpaValidationError(
                f"{n_nan} IPA row(s) have no limma {col}. Every IPA row is "
                "'complete & regulated in {UP, DOWN}' and therefore limma-eligible, "
                "so this means the eligibility assumption broke (an ON_OFF or "
                f"incomplete protein reached the IPA set, or limma dropped rows). "
                f"First offenders: {lost}"
            )

    return merged[IPA_FULL_COLS]


def write_ipa_full(df_full, csv_path, txt_path=None):
    """Write ``ipa_input_full.csv`` and its tab-delimited ``.txt`` twin.

    Both carry :data:`IPA_FULL_COLS` with the identifier leftmost, UTF-8 without
    a BOM and a single header row. The ``.txt`` exists because QIAGEN
    "strongly recommend ... tab-delimited text format (.txt ...) for faster
    upload" (research1.md line 195); the ``.csv`` stays because it is the format
    the rest of this pipeline speaks and the one the earlier Excel->CSV fix
    landed on.

    Passing ``txt_path=False`` writes only the CSV; passing ``None`` derives the
    ``.txt`` sibling from `csv_path`.

    Returns
    -------
    list[pathlib.Path]
        The paths written, CSV first.
    """
    csv_path = Path(csv_path)
    if txt_path is None:
        txt_path = csv_path.with_suffix(".txt")

    written = [csv_path]
    df_full[IPA_FULL_COLS].to_csv(csv_path, index=False, encoding="utf-8")
    if txt_path is not False:
        txt_path = Path(txt_path)
        df_full[IPA_FULL_COLS].to_csv(
            txt_path, index=False, sep="\t", encoding="utf-8"
        )
        written.append(txt_path)
    return written


# ---------------------------------------------------------------------------
# Expected row count -- derived, never inlined
# ---------------------------------------------------------------------------


def frozen_counts(path=None):
    """Load ``tests/expected/frozen_counts.json``.

    The dataset-specific expectations live in exactly one place. Nothing in this
    module types ``715`` (or ``0``, or ``1938``) as a literal.
    """
    path = Path(path) if path else _PKG_DIR / "tests" / "expected" / "frozen_counts.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def expected_ipa_rows(foldchange_csv=None, counts_path=None):
    """How many rows ``ipa_input.csv`` must have.

    research1.md line 187 states the rule as *"assert IPA file row count == count
    of (complete & UP|DOWN)"*, so when ``foldchange_all.csv`` is available the
    number is **recomputed from the data** -- that is the check the roadmap
    actually asks for. The frozen JSON is the fallback for callers that only
    have the IPA file in hand; both agree at 715 for the committed dataset.
    """
    if foldchange_csv is None:
        default = _PKG_DIR / "results" / "foldchange_all.csv"
        foldchange_csv = default if default.exists() else None
    if foldchange_csv is not None and Path(foldchange_csv).exists():
        fc = pd.read_csv(foldchange_csv, float_precision="round_trip")
        mask = fc["complete"].astype(bool) & fc["regulated"].isin(REGULATED_LABELS)
        return int(mask.sum())
    return int(frozen_counts(counts_path)["ipa_input_rows"])


# ---------------------------------------------------------------------------
# Validator
# ---------------------------------------------------------------------------


def _delimiter_for(path):
    return _DELIMITERS.get(Path(path).suffix.lower(), ",")


def _check(report, failures, name, ok, detail):
    report[name] = {"ok": bool(ok), "detail": detail}
    if not ok:
        failures.append(f"{name}: {detail}")


def validate_ipa(
    path,
    *,
    expected_rows=None,
    required_columns=None,
    measurement_columns=None,
    id_col=IPA_ID_COL,
    strict=True,
):
    """Check an IPA upload file against QIAGEN's contract.

    This is research1.md line 187 ("Export validation") implemented for real. The
    only export check the pipeline had before was a bare
    ``assert len(df_ipa) == 715`` on an in-memory frame, which says nothing about
    the bytes that actually reach QIAGEN.

    Checks, in order:

    ``exists``
        The file is there and non-empty.
    ``no_bom``
        First three bytes are not ``EF BB BF``. A BOM turns the first header
        cell into ``"\\ufeffUniProt Accession Number"``, which is exactly the
        class of "IPA does not recognize the column headers" failure the
        Excel->CSV fix was chasing.
    ``utf8``
        The whole file decodes as UTF-8.
    ``single_header_row``
        Exactly one header line, every data line has the same field count, and no
        data line repeats the header (the classic symptom of two files
        concatenated).
    ``identifier_leftmost``
        Column 0 is the identifier column (research1.md: *"Identifier column
        leftmost"*).
    ``identifier_clean``
        No identifier is empty, whitespace-bearing or control-character-bearing.
        ``;`` is explicitly allowed -- 21 of the 715 accessions are legitimate
        protein groups like ``P05132;P68181``.
    ``required_columns``
        The expected column set is present, in order.
    ``measurements_finite``
        No NA and no inf in any measurement column present in the file.
    ``row_count``
        Equals ``(complete & regulated in {UP, DOWN}).sum()``.

    Parameters
    ----------
    path :
        File to validate (``.csv``, or ``.txt``/``.tsv`` for tab-delimited).
    expected_rows :
        Row count to require. Defaults to :func:`expected_ipa_rows`. Pass an
        explicit number for files with a different contract (e.g. the header-only
        significant file, which :func:`validate_ipa_significant` wraps).
    required_columns :
        Defaults to :data:`IPA_COLS` for a 4-column file, :data:`IPA_FULL_COLS`
        for a 6-column one -- inferred from the header so both files validate
        with the same call.
    strict :
        Raise :class:`IpaValidationError` listing every failure (default). With
        ``strict=False`` the report is returned and the caller decides.

    Returns
    -------
    dict
        ``{"path", "n_rows", "columns", "checks", "ok", "failures"}``.
    """
    path = Path(path)
    report = {"path": str(path), "n_rows": None, "columns": None}
    checks = {}
    failures = []
    report["checks"] = checks

    if not path.exists():
        _check(checks, failures, "exists", False, f"{path} does not exist")
        return _finish(report, failures, strict)

    raw = path.read_bytes()
    _check(checks, failures, "exists", len(raw) > 0, f"{len(raw)} bytes")
    if not raw:
        return _finish(report, failures, strict)

    _check(
        checks, failures, "no_bom",
        raw[:3] != _UTF8_BOM,
        "no UTF-8 BOM" if raw[:3] != _UTF8_BOM else
        "file starts with a UTF-8 BOM; IPA will read it as part of the first column name",
    )

    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        _check(checks, failures, "utf8", False, f"not valid UTF-8: {exc}")
        return _finish(report, failures, strict)
    _check(checks, failures, "utf8", True, "decodes as UTF-8")

    delim = _delimiter_for(path)

    # --- one header row, uniform field counts -------------------------------
    rows = list(csv.reader(text.splitlines(), delimiter=delim))
    if not rows:
        _check(checks, failures, "single_header_row", False, "no rows at all")
        return _finish(report, failures, strict)

    header, data_rows = rows[0], rows[1:]
    report["columns"] = header
    report["n_rows"] = len(data_rows)

    ragged = [i + 2 for i, r in enumerate(data_rows) if len(r) != len(header)]
    repeats = [i + 2 for i, r in enumerate(data_rows) if r == header]
    header_ok = not ragged and not repeats
    _check(
        checks, failures, "single_header_row", header_ok,
        f"1 header + {len(data_rows)} data rows, all {len(header)} fields wide"
        if header_ok else
        f"ragged lines at {ragged[:5]}; header repeated at lines {repeats[:5]}",
    )

    _check(
        checks, failures, "identifier_leftmost",
        header[:1] == [id_col],
        f"column 0 is {header[0]!r}" if header else "empty header",
    )

    if required_columns is None:
        required_columns = IPA_FULL_COLS if len(header) == len(IPA_FULL_COLS) else IPA_COLS
    _check(
        checks, failures, "required_columns",
        header == list(required_columns),
        f"header {header} == {list(required_columns)}" if header == list(required_columns)
        else f"header {header} != expected {list(required_columns)}",
    )

    # --- identifiers free of stray characters -------------------------------
    if header and header[0] == id_col:
        dirty = []
        for i, r in enumerate(data_rows):
            if not r:
                continue
            v = r[0]
            if v == "" or v != v.strip() or any(ch.isspace() or ord(ch) < 32 for ch in v):
                dirty.append((i + 2, v))
        _check(
            checks, failures, "identifier_clean",
            not dirty,
            f"{len(data_rows)} identifiers, none empty or whitespace/control bearing"
            if not dirty else f"stray characters at {dirty[:5]}",
        )

    # --- no NA / inf in measurement columns ---------------------------------
    if measurement_columns is None:
        measurement_columns = [c for c in MEASUREMENT_COLS if c in header]
    frame = pd.read_csv(path, sep=delim, float_precision="round_trip")
    bad = {}
    for col in measurement_columns:
        if col not in frame.columns:
            bad[col] = "column absent"
            continue
        series = pd.to_numeric(frame[col], errors="coerce")
        n_na = int(series.isna().sum())
        n_inf = int(series.apply(lambda x: isinstance(x, float) and math.isinf(x)).sum())
        if n_na or n_inf:
            bad[col] = f"{n_na} NA, {n_inf} inf"
    _check(
        checks, failures, "measurements_finite",
        not bad,
        f"{list(measurement_columns)} all finite" if not bad else f"{bad}",
    )

    # --- row count ----------------------------------------------------------
    if expected_rows is None:
        expected_rows = expected_ipa_rows()
    _check(
        checks, failures, "row_count",
        len(data_rows) == expected_rows,
        f"{len(data_rows)} rows == expected {expected_rows}"
        if len(data_rows) == expected_rows
        else f"{len(data_rows)} rows, expected {expected_rows} "
             "(complete & regulated in {UP, DOWN})",
    )

    return _finish(report, failures, strict)


def _finish(report, failures, strict):
    report["failures"] = failures
    report["ok"] = not failures
    if failures and strict:
        raise IpaValidationError(
            f"{report['path']} violates the IPA upload contract:\n  - "
            + "\n  - ".join(failures)
        )
    return report


def validate_ipa_significant(path, *, counts_path=None, strict=True):
    """Assert the header-only contract of ``ipa_input_significant.csv``.

    **This file being empty is the result, not a defect.** 0 of 1938 proteins
    pass FDR < 0.05 -- the expected n=2 technical-replicate ceiling
    (DECISIONS_LOG D2, research1.md SECTION 5). So the contract is: correct
    header, exactly ``frozen_counts["n_significant_fdr05"]`` data rows. If that
    number is ever non-zero the JSON changes with the data, not this code, and
    the check keeps holding. Never "fix" this file by relaxing the FDR cutoff.
    """
    expected = int(frozen_counts(counts_path)["n_significant_fdr05"])
    return validate_ipa(
        path,
        expected_rows=expected,
        required_columns=IPA_SIGNIFICANT_COLS,
        strict=strict,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_full_from_results(results_dir=None, *, validate=True):
    """Build ``ipa_input_full.{csv,txt}`` from the committed results directory.

    Reads the two files rather than recomputing them, deliberately: the log2FC
    that ships in the full file is then *the same text* as in the
    ``ipa_input.csv`` the user uploaded, which is only true because every read
    here uses ``float_precision="round_trip"`` (see the module docstring).

    Returns
    -------
    list[pathlib.Path]
    """
    results_dir = Path(results_dir) if results_dir else _PKG_DIR / "results"
    ipa_path = results_dir / "ipa_input.csv"
    limma_path = results_dir / "qc_limma.csv"
    for p in (ipa_path, limma_path):
        if not p.exists():
            raise FileNotFoundError(
                f"{p} not found -- run foldchange.py (which runs limma) first."
            )

    df_ipa = pd.read_csv(ipa_path, float_precision="round_trip")
    limma = pd.read_csv(limma_path, float_precision="round_trip")
    df_full = build_ipa_full(df_ipa, limma)

    written = write_ipa_full(df_full, results_dir / "ipa_input_full.csv")
    if validate:
        for p in written:
            validate_ipa(p, required_columns=IPA_FULL_COLS)
    return written


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", default=None,
                    help="directory holding ipa_input.csv / qc_limma.csv "
                         "(default: proteomics_de/results)")
    ap.add_argument("--no-sidecars", action="store_true",
                    help="skip writing the .provenance.json sidecars")
    args = ap.parse_args(argv)

    results_dir = Path(args.results_dir) if args.results_dir else _PKG_DIR / "results"

    written = build_full_from_results(results_dir)
    for p in written:
        rep = validate_ipa(p, required_columns=IPA_FULL_COLS)
        print(f"Saved {p} ({rep['n_rows']} rows, {len(rep['columns'])} columns)")

    # Re-validate the frozen four-column file too. It is not rewritten here; the
    # point is to prove the file on disk still satisfies the upload contract.
    for name, validator in (
        ("ipa_input.csv", lambda p: validate_ipa(p, required_columns=IPA_COLS)),
        ("ipa_input_significant.csv", validate_ipa_significant),
    ):
        p = results_dir / name
        if p.exists():
            rep = validator(p)
            print(f"Validated {p} ({rep['n_rows']} rows, "
                  f"{len(rep['checks'])} checks passed)")

    if not args.no_sidecars:
        if str(_ROOT) not in sys.path:
            sys.path.insert(0, str(_ROOT))
        from proteomics_de.provenance import emit_default_sidecars

        for p in emit_default_sidecars(results_dir):
            print(f"Saved {p}")

    return 0


__all__ = [
    "IPA_ID_COL",
    "IPA_COLS",
    "IPA_FULL_COLS",
    "IPA_SIGNIFICANT_COLS",
    "MEASUREMENT_COLS",
    "IpaValidationError",
    "write_ipa",
    "build_ipa_full",
    "write_ipa_full",
    "frozen_counts",
    "expected_ipa_rows",
    "validate_ipa",
    "validate_ipa_significant",
    "build_full_from_results",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
