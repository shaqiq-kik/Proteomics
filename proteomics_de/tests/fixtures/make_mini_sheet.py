"""Generate ``mini_sheet.xlsx`` — a hand-designed 31-protein SILAC workbook.

WHY THIS EXISTS
---------------
The suite has 500+ tests, and every one of them either exercises a unit in
isolation or asserts a frozen fact about the **one** committed dataset. Nothing
proved the pipeline still works on *different* data. A refactor that broke, say,
the ON_OFF detection for blank-cell (rather than zero) absences would sail past
the whole suite, because the committed workbook happens to encode absence as a
zero.

This module is the other half: a tiny workbook in the real 29-column MaxQuant
layout whose every row was designed to land in a **known** branch of
``etl/foldchange_core.py``, together with the expected outcome for that row.
``tests/test_e2e_smoke.py`` runs the real ``foldchange.py`` (including the R
limma leg) over it and checks every row against :data:`ROWS`.

The design table below is the single source of truth: the workbook and the
expectations are generated from the same rows, so they cannot drift apart. The
.xlsx is committed so the test needs no generation step; re-run this file to
regenerate it::

    .venv/bin/python proteomics_de/tests/fixtures/make_mini_sheet.py

D7 ORIENTATION — READ THIS BEFORE CHANGING A NUMBER
---------------------------------------------------
``config/sample_sheet.tsv`` (and DECISIONS_LOG **D7**) say:

======================  ==========  ================  =========
channel                 sheet       condition         replicate
======================  ==========  ================  =========
``Intensity 31578``     Report L    **treated**       1
``Intensity 31580``     Report L    **treated**       2
``Intensity 31579``     Report H    **control**       1
``Intensity 31581``     Report H    **control**       2
======================  ==========  ================  =========

So ``ratio_rep1 = 31578 / 31579`` and ``ratio_rep2 = 31580 / 31581`` — treated
over control — and a protein with more signal in the *Light* sheet is UP. The
pipeline shipped this inverted once; the fixture pins the corrected orientation
in two independent ways (per-row ``log2FC`` signs, and the ``detected_in``
labels of the single-condition rescue file, where an L-only protein must come
out as ``treated_only``).

WHY THE INTENSITIES ARE POWERS OF TWO
-------------------------------------
Every ratio is an exact power of two (or exactly 1.5), so every expected
``log2FC`` is an exact binary float and the test can assert equality rather than
"close enough". A tolerance would hide precisely the kind of small systematic
drift these tests exist to catch. The three ``Bndry`` rows are the exception and
are the point of the exercise — see :data:`ROWS`.

WHAT EACH ROW COVERS
--------------------
* complete rows in all three classes (UP / DOWN / NO CHANGE), positive and
  negative, single- and double-replicate agreement;
* the **Bug 1** witness (``Nochg1``): mean-of-log2 says NO CHANGE where the
  legacy mean-of-ratios would have said UP;
* the **threshold razor edge**: ``config/constants.LOG2_THRESHOLD`` is the
  *rounded literal* ``0.585``, which is very slightly **above**
  ``log2(1.5) = 0.5849625``. A protein up exactly 1.5-fold therefore lands in
  NO CHANGE. Three rows pin that (``Bndry1``/``Bndry3`` just inside,
  ``Bndry2`` outside), so replacing the literal with ``math.log2(1.5)`` — a
  change that looks like a harmless tidy-up — flips a label and fails loudly;
* ON_OFF in **both** directions and in **both** absence encodings (a zero and
  an empty cell);
* the "absent on both sides" row, which must stay NO CHANGE rather than being
  invented into an ON_OFF signal;
* incomplete-but-not-ON_OFF rows, including the **zero denominator** that Bug 3
  exists to prevent from ever reaching a division;
* semicolon protein groups, in both the fold-change and single-condition paths;
* a missing gene symbol;
* single-condition rescue in both directions (L-only and H-only);
* a **D11** junk accession (a MaxQuant row-index list) that must be quarantined
  out of ``single_condition_proteins.csv`` rather than shipped.

DELIBERATE LAYOUT CHOICES
-------------------------
* The two ``Quant Report`` sheets of the real workbook are **not** reproduced.
  ``foldchange.py`` reads only ``Protein Report L`` / ``Protein Report H``, and
  a sheet nothing reads is a sheet nothing tests.
* The Heavy sheet carries **no blank intensities**, so its intensity columns
  load as ``int64`` — which is what makes
  ``foldchange_core.restore_left_order``'s dtype restore do real work. The Light
  sheet does carry blanks (it must, to encode absence both ways), so its columns
  load as ``float64``. That int/float asymmetry mirrors the real workbook, and
  it is why ``foldchange_all.csv`` shows ``2.0e9`` for 31578/31580 but bare
  integers for 31579/31581.
* Document properties are pinned to a fixed timestamp so regenerating the
  workbook from unchanged inputs produces identical bytes.
"""

from __future__ import annotations

import datetime as _dt
import math
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
OUTPUT_PATH = HERE / "mini_sheet.xlsx"

SHEET_L = "Protein Report L"
SHEET_H = "Protein Report H"

#: Light-sheet channels, in workbook order. TREATED under D7.
L_SAMPLES = ("31578", "31580")
#: Heavy-sheet channels, in workbook order. CONTROL under D7.
H_SAMPLES = ("31579", "31581")

GB = 1_000_000_000

#: ``log2(1.5) = 0.5849625...`` -- the value ``config/constants.LOG2_THRESHOLD``
#: is a rounded-up literal of (``0.585``). Written as an expression, not as
#: digits, so the three ``Bndry`` rows below state *why* they sit where they do.
#: These are the only expected fold changes that are not exact binary fractions;
#: the test asserts them to 1 ulp and everything else to the bit.
L15 = math.log2(1.5)

#: openpyxl stamps a wall clock into docProps/core.xml; pin it so the committed
#: .xlsx is reproducible from this script.
_FIXED_TIMESTAMP = _dt.datetime(2026, 1, 1, 0, 0, 0)

# ---------------------------------------------------------------------------
# The design table
# ---------------------------------------------------------------------------
# Keys:
#   acc          accession (the merge key; whole string, groups included)
#   gene         gene symbol, or None for a blank cell
#   t1, t2       Intensity 31578 / 31580   (Light sheet, TREATED)
#   c1, c2       Intensity 31579 / 31581   (Heavy sheet, CONTROL)
#                a number, 0 (measured-as-absent), or None (blank cell).
#                None on a side means "omit from that sheet" only when the
#                corresponding `sheets` entry says so -- see `sheets`.
#   sheets       "both" | "L" | "H".  Which protein report(s) the row appears
#                in. "L"/"H" rows are the single-condition rescue path.
#   note         why the row exists
#   expect       what the pipeline must produce, as a dict:
#       complete    bool
#       regulated   UP | DOWN | NO CHANGE | ON_OFF
#       onoff       "" | on_with_treatment | off_with_treatment
#       log2FC      exact expected value, or None for NaN
#       n_imputed   how many of the 4 channels limma must impute (both-rows)
#       detected_in control_only | treated_only          (single-condition rows)
#       quarantined True                                  (D11 rows)
# ---------------------------------------------------------------------------

ROWS: tuple[dict, ...] = (
    # --- complete, UP -----------------------------------------------------
    dict(acc="P10001", gene="Upreg1", t1=8 * GB, t2=16 * GB, c1=1 * GB, c2=2 * GB,
         sheets="both", note="both replicates agree, 8x up",
         expect=dict(complete=True, regulated="UP", onoff="", log2FC=3.0, n_imputed=0)),
    dict(acc="P10002", gene="Upreg2", t1=4 * GB, t2=2 * GB, c1=1 * GB, c2=2 * GB,
         sheets="both", note="one replicate 4x up, the other flat -> mean 1.0",
         expect=dict(complete=True, regulated="UP", onoff="", log2FC=1.0, n_imputed=0)),
    dict(acc="O10003", gene="Upreg3", t1=2 * GB, t2=8 * GB, c1=1 * GB, c2=4 * GB,
         sheets="both", note="O-prefixed accession: the [OPQ] regex branch",
         expect=dict(complete=True, regulated="UP", onoff="", log2FC=1.0, n_imputed=0)),
    dict(acc="P10004", gene="Upreg4", t1=16 * GB, t2=2 * GB, c1=1 * GB, c2=2 * GB,
         sheets="both", note="16x / flat -> mean 2.0",
         expect=dict(complete=True, regulated="UP", onoff="", log2FC=2.0, n_imputed=0)),

    # --- complete, DOWN ---------------------------------------------------
    dict(acc="P10005", gene="Dnreg1", t1=1 * GB, t2=2 * GB, c1=8 * GB, c2=16 * GB,
         sheets="both", note="both replicates agree, 8x down",
         expect=dict(complete=True, regulated="DOWN", onoff="", log2FC=-3.0, n_imputed=0)),
    dict(acc="P10006", gene="Dnreg2", t1=1 * GB, t2=2 * GB, c1=4 * GB, c2=2 * GB,
         sheets="both", note="one replicate 4x down, the other flat",
         expect=dict(complete=True, regulated="DOWN", onoff="", log2FC=-1.0, n_imputed=0)),
    dict(acc="Q10007", gene="Dnreg3", t1=1 * GB, t2=4 * GB, c1=2 * GB, c2=8 * GB,
         sheets="both", note="Q-prefixed accession, 2x down in both replicates",
         expect=dict(complete=True, regulated="DOWN", onoff="", log2FC=-1.0, n_imputed=0)),

    # --- complete, NO CHANGE ---------------------------------------------
    dict(acc="P10008", gene="Nochg1", t1=4 * GB, t2=1 * GB, c1=1 * GB, c2=4 * GB,
         sheets="both",
         note="BUG 1 WITNESS. ratios are 4 and 0.25: log2-then-mean = 0.0 (NO "
              "CHANGE, correct); the legacy mean-then-log2 gives "
              "log2(mean(4, 0.25)) = log2(2.125) = +1.09 and would call this UP.",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="", log2FC=0.0, n_imputed=0)),
    dict(acc="P10009", gene="Nochg2", t1=2 * GB, t2=4 * GB, c1=1 * GB, c2=4 * GB,
         sheets="both", note="+0.5 -- inside the +0.585 cutoff",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="", log2FC=0.5, n_imputed=0)),
    dict(acc="P10010", gene="Nochg3", t1=1 * GB, t2=4 * GB, c1=2 * GB, c2=4 * GB,
         sheets="both", note="-0.5 -- inside the -0.585 cutoff (Bug 2 symmetry)",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="", log2FC=-0.5, n_imputed=0)),
    dict(acc="P10011", gene="Nochg4", t1=1 * GB, t2=4 * GB, c1=1 * GB, c2=4 * GB,
         sheets="both", note="flat in both replicates",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="", log2FC=0.0, n_imputed=0)),
    dict(acc="Q10012", gene="Nochg5", t1=1 * GB, t2=16 * GB, c1=4 * GB, c2=4 * GB,
         sheets="both", note="replicates disagree violently (-2, +2) and cancel",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="", log2FC=0.0, n_imputed=0)),

    # --- the threshold razor edge ----------------------------------------
    dict(acc="P10013", gene="Bndry1", t1=1_500_000_000, t2=3 * GB, c1=1 * GB, c2=2 * GB,
         sheets="both",
         note="up EXACTLY 1.5-fold in both replicates -> log2FC = 0.5849625, "
              "which is BELOW the rounded 0.585 literal in config/constants.py. "
              "NO CHANGE. Replacing that literal with math.log2(1.5) flips this "
              "row to UP -- which is the point of the row.",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="",
                     log2FC=L15, n_imputed=0)),
    dict(acc="P10014", gene="Bndry2", t1=1_500_000_000, t2=4 * GB, c1=1 * GB, c2=2 * GB,
         sheets="both", note="1.5x and 2x -> 0.7925, clear of the cutoff",
         expect=dict(complete=True, regulated="UP", onoff="",
                     log2FC=(L15 + 1.0) / 2, n_imputed=0)),
    dict(acc="P10015", gene="Bndry3", t1=1 * GB, t2=2 * GB, c1=1_500_000_000, c2=3 * GB,
         sheets="both", note="the negative razor edge: -0.5849625, still NO CHANGE",
         expect=dict(complete=True, regulated="NO CHANGE", onoff="",
                     log2FC=-L15, n_imputed=0)),

    # --- protein group and missing gene symbol ---------------------------
    dict(acc="P10016;Q10017", gene="Grpaa;Grpbb", t1=4 * GB, t2=8 * GB, c1=1 * GB, c2=2 * GB,
         sheets="both",
         note="semicolon protein group -- the accession is an IDENTITY, never "
              "split for the merge (etl/accessions.py POLICY)",
         expect=dict(complete=True, regulated="UP", onoff="", log2FC=2.0, n_imputed=0)),
    dict(acc="P10018", gene=None, t1=1 * GB, t2=2 * GB, c1=4 * GB, c2=8 * GB,
         sheets="both", note="no resolved gene symbol -- gene columns are nullable",
         expect=dict(complete=True, regulated="DOWN", onoff="", log2FC=-2.0, n_imputed=0)),

    # --- ON_OFF, both directions, both absence encodings ------------------
    dict(acc="P10019", gene="Onoff1", t1=500_000_000, t2=600_000_000, c1=0, c2=0,
         sheets="both", note="control side absent as ZEROS -> on_with_treatment",
         expect=dict(complete=False, regulated="ON_OFF", onoff="on_with_treatment",
                     log2FC=None, n_imputed=None)),
    dict(acc="P10020", gene="Onoff2", t1=900_000_000, t2=1_100_000_000, c1=0, c2=0,
         sheets="both", note="second on_with_treatment, different magnitudes",
         expect=dict(complete=False, regulated="ON_OFF", onoff="on_with_treatment",
                     log2FC=None, n_imputed=None)),
    dict(acc="P10021", gene="Onoff3", t1=0, t2=0, c1=700_000_000, c2=800_000_000,
         sheets="both", note="treated side absent as ZEROS -> off_with_treatment",
         expect=dict(complete=False, regulated="ON_OFF", onoff="off_with_treatment",
                     log2FC=None, n_imputed=None)),
    dict(acc="P10022", gene="Onoff4", t1=None, t2=None, c1=600_000_000, c2=500_000_000,
         sheets="both",
         note="treated side absent as BLANK CELLS -> off_with_treatment. The "
              "committed workbook only ever encodes absence as a zero, so this "
              "is the branch no frozen-count test can reach.",
         expect=dict(complete=False, regulated="ON_OFF", onoff="off_with_treatment",
                     log2FC=None, n_imputed=None)),

    # --- incomplete but NOT ON_OFF ---------------------------------------
    dict(acc="P10023", gene="Zerodn", t1=3 * GB, t2=2 * GB, c1=0, c2=1 * GB,
         sheets="both",
         note="ZERO DENOMINATOR. Without the `complete` mask, 3e9/0 = inf would "
              "reach log2FC (Bug 3). The control side is not entirely absent, "
              "so this is NOT ON_OFF -- it is an incomplete NO CHANGE.",
         expect=dict(complete=False, regulated="NO CHANGE", onoff="", log2FC=None,
                     n_imputed=1)),
    dict(acc="P10024", gene="Partia", t1=3 * GB, t2=None, c1=1 * GB, c2=2 * GB,
         sheets="both", note="one blank numerator -> incomplete, not ON_OFF",
         expect=dict(complete=False, regulated="NO CHANGE", onoff="", log2FC=None,
                     n_imputed=1)),
    dict(acc="P10025", gene="Bothpt", t1=0, t2=2 * GB, c1=1 * GB, c2=0,
         sheets="both", note="one value missing on EACH side; neither side fully absent",
         expect=dict(complete=False, regulated="NO CHANGE", onoff="", log2FC=None,
                     n_imputed=2)),
    dict(acc="P10026", gene="Empty1", t1=0, t2=0, c1=0, c2=0,
         sheets="both",
         note="absent on BOTH sides. An empty row is not evidence of anything, "
              "so it must stay NO CHANGE -- calling it ON_OFF would invent a "
              "signal out of nothing.",
         expect=dict(complete=False, regulated="NO CHANGE", onoff="", log2FC=None,
                     n_imputed=4)),

    # --- single condition: Light sheet only -> treated_only (D7) ----------
    dict(acc="P10027", gene="Tonly1", t1=2 * GB, t2=3 * GB, c1=None, c2=None,
         sheets="L", note="detected only in the treated channels",
         expect=dict(detected_in="treated_only")),
    dict(acc="P10028;Q10029", gene="Tgrpa;Tgrpb", t1=1 * GB, t2=1_500_000_000,
         c1=None, c2=None, sheets="L",
         note="protein group on the single-condition path",
         expect=dict(detected_in="treated_only")),
    dict(acc="12345;67890;24680", gene=None, t1=None, t2=None, c1=None, c2=None,
         sheets="L",
         note="DECISIONS_LOG D11: a MaxQuant ROW-INDEX LIST where an accession "
              "belongs. Every token is a bare integer, so it is quarantined to "
              "results/qc/quarantine_accessions.csv and never reaches "
              "single_condition_proteins.csv. Note the row above is a 13-char "
              "protein GROUP and must survive -- the discriminator is token "
              "shape, never length.",
         expect=dict(quarantined=True)),

    # --- single condition: Heavy sheet only -> control_only (D7) ----------
    dict(acc="P10030", gene="Conly1", t1=None, t2=None, c1=400_000_000, c2=500_000_000,
         sheets="H", note="detected only in the control channels",
         expect=dict(detected_in="control_only")),
    dict(acc="Q10031", gene="Conly2", t1=None, t2=None, c1=1_200_000_000, c2=900_000_000,
         sheets="H", note="second control-only protein",
         expect=dict(detected_in="control_only")),
    dict(acc="P10032", gene="Conly3", t1=None, t2=None, c1=0, c2=300_000_000,
         sheets="H", note="control-only with one channel measured as zero",
         expect=dict(detected_in="control_only")),
)


# ---------------------------------------------------------------------------
# Derived views of the design table (imported by tests/test_e2e_smoke.py)
# ---------------------------------------------------------------------------

def both_rows() -> list[dict]:
    """Rows present in both protein reports -- the fold-change population."""
    return [r for r in ROWS if r["sheets"] == "both"]


def light_only_rows() -> list[dict]:
    return [r for r in ROWS if r["sheets"] == "L"]


def heavy_only_rows() -> list[dict]:
    return [r for r in ROWS if r["sheets"] == "H"]


def single_condition_rows() -> list[dict]:
    """Single-condition rows that SURVIVE quarantine (i.e. reach the CSV)."""
    return [r for r in ROWS
            if r["sheets"] in ("L", "H") and not r["expect"].get("quarantined")]


def quarantined_rows() -> list[dict]:
    return [r for r in ROWS if r["expect"].get("quarantined")]


def eligible_rows() -> list[dict]:
    """Both-sheet rows that limma tests: everything except the ON_OFF ones."""
    return [r for r in both_rows() if r["expect"]["regulated"] != "ON_OFF"]


def expected_counts() -> dict:
    """Every count the end-to-end test asserts, derived from :data:`ROWS`.

    Mirrors the shape of ``tests/expected/frozen_counts.json`` but for the
    fixture, so the same reasoning applies to both datasets. Nothing here is a
    literal: change a row above and these follow.
    """
    both = both_rows()
    reg = [r["expect"]["regulated"] for r in both]
    onoff = [r["expect"]["onoff"] for r in both]
    return {
        "sheet_L_rows": len(both) + len(light_only_rows()),
        "sheet_H_rows": len(both) + len(heavy_only_rows()),
        "foldchange_all_rows": len(both),
        "complete_rows": sum(1 for r in both if r["expect"]["complete"]),
        "n_up": reg.count("UP"),
        "n_down": reg.count("DOWN"),
        "n_nochange": reg.count("NO CHANGE"),
        "n_onoff": reg.count("ON_OFF"),
        "n_on": onoff.count("on_with_treatment"),
        "n_off": onoff.count("off_with_treatment"),
        "ipa_input_rows": sum(
            1 for r in both
            if r["expect"]["complete"] and r["expect"]["regulated"] in ("UP", "DOWN")
        ),
        "qc_limma_rows": len(eligible_rows()),
        "single_condition_rows": len(single_condition_rows()),
        "single_condition_rows_pre_quarantine": (
            len(light_only_rows()) + len(heavy_only_rows())
        ),
        "quarantined_accessions": len(quarantined_rows()),
        "background_union": len(both) + len(single_condition_rows()),
    }


# ---------------------------------------------------------------------------
# Workbook construction
# ---------------------------------------------------------------------------

def _protein_report_columns(samples: tuple[str, str]) -> list[str]:
    """The real 29-column ``Protein Report`` header for a two-channel sheet.

    Verbatim from ``Copy of General Sheet.xlsx`` (both report sheets share this
    layout; only the sample ids in the per-channel column names differ).
    """
    a, b = samples
    return [
        "UniProt Accession Number",
        "Protein names",
        "Gene names",
        "Mol. weight [kDa]",
        "Number of proteins",
        f"Peptides {a}", f"Peptides {b}",
        f"Razor + unique peptides {a}", f"Razor + unique peptides {b}",
        f"Unique peptides {a}", f"Unique peptides {b}",
        f"Sequence coverage {a} [%]", f"Sequence coverage {b} [%]",
        "Q-value",
        "Score",
        f"Identification type {a}", f"Identification type {b}",
        f"Ratio H/L {a}", f"Ratio H/L {b}",
        f"Ratio H/L normalized {a}", f"Ratio H/L normalized {b}",
        f"Ratio H/L count {a}", f"Ratio H/L count {b}",
        f"Intensity {a}", f"Intensity L {a}", f"Intensity H {a}",
        f"Intensity {b}", f"Intensity L {b}", f"Intensity H {b}",
    ]


def _split_isotopes(total):
    """Split a combined intensity into (Light, Heavy) halves.

    DECISIONS_LOG **D8** records the verified property of the real workbook that
    ``Intensity L + Intensity H == Intensity`` exactly. The fixture reproduces
    it so anyone acting on D8 (i.e. switching the pipeline onto the SILAC
    ratios) finds a fixture that already carries the columns they need.
    """
    if total is None:
        return None, None
    light = int(total) // 2
    return light, int(total) - light


def _sheet_frame(rows: list[dict], samples: tuple[str, str], side: str) -> pd.DataFrame:
    """Build one protein-report sheet from `rows`.

    `side` is "L" or "H" and selects which pair of intensities each row
    contributes (``t1/t2`` for Light, ``c1/c2`` for Heavy).
    """
    a, b = samples
    columns = _protein_report_columns(samples)
    records = []
    for i, row in enumerate(rows):
        va, vb = (row["t1"], row["t2"]) if side == "L" else (row["c1"], row["c2"])
        la, ha = _split_isotopes(va)
        lb, hb = _split_isotopes(vb)
        # Filler columns are deterministic functions of the row position, so the
        # workbook regenerates identically and nothing downstream can depend on
        # a value that "looks real" but is arbitrary.
        n_pep = 2 + (i % 7)
        records.append({
            "UniProt Accession Number": row["acc"],
            "Protein names": f"Synthetic fixture protein {i + 1}",
            "Gene names": row["gene"],
            "Mol. weight [kDa]": round(20.0 + i * 3.5, 3),
            "Number of proteins": row["acc"].count(";") + 1,
            f"Peptides {a}": n_pep, f"Peptides {b}": n_pep + 1,
            f"Razor + unique peptides {a}": n_pep, f"Razor + unique peptides {b}": n_pep + 1,
            f"Unique peptides {a}": max(1, n_pep - 1),
            f"Unique peptides {b}": max(1, n_pep),
            f"Sequence coverage {a} [%]": round(1.0 + i * 0.7, 3),
            f"Sequence coverage {b} [%]": round(1.2 + i * 0.7, 3),
            "Q-value": 0,
            "Score": round(30.0 + i * 4.25, 3),
            f"Identification type {a}": "By MS/MS",
            f"Identification type {b}": "By MS/MS",
            f"Ratio H/L {a}": None, f"Ratio H/L {b}": None,
            f"Ratio H/L normalized {a}": None, f"Ratio H/L normalized {b}": None,
            f"Ratio H/L count {a}": n_pep, f"Ratio H/L count {b}": n_pep + 1,
            f"Intensity {a}": va, f"Intensity L {a}": la, f"Intensity H {a}": ha,
            f"Intensity {b}": vb, f"Intensity L {b}": lb, f"Intensity H {b}": hb,
        })
    return pd.DataFrame(records, columns=columns)


def build_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(Protein Report L, Protein Report H)`` as DataFrames."""
    l_rows = [r for r in ROWS if r["sheets"] in ("both", "L")]
    h_rows = [r for r in ROWS if r["sheets"] in ("both", "H")]
    return _sheet_frame(l_rows, L_SAMPLES, "L"), _sheet_frame(h_rows, H_SAMPLES, "H")


def _validate_design() -> None:
    """Fail loudly on a design table that could not possibly be right.

    Cheap, and it runs on import, so a typo in :data:`ROWS` surfaces as a clear
    message at collection time rather than as a baffling assertion failure deep
    inside the end-to-end test.
    """
    accs = [r["acc"] for r in ROWS]
    dupes = {a for a in accs if accs.count(a) > 1}
    assert not dupes, f"duplicate accessions in the design table: {sorted(dupes)}"
    for r in ROWS:
        if r["sheets"] == "both":
            exp = r["expect"]
            missing = {"complete", "regulated", "onoff", "log2FC"} - set(exp)
            assert not missing, f"{r['acc']}: expectation missing {sorted(missing)}"
            assert exp["regulated"] in ("UP", "DOWN", "NO CHANGE", "ON_OFF"), r["acc"]
            # A label and a fold change must agree: this is the same direction
            # invariant tests/test_golden_outputs.py asserts on the real data.
            if exp["regulated"] == "UP":
                assert exp["log2FC"] is not None and exp["log2FC"] > 0, r["acc"]
            if exp["regulated"] == "DOWN":
                assert exp["log2FC"] is not None and exp["log2FC"] < 0, r["acc"]
            if exp["regulated"] == "ON_OFF":
                assert exp["log2FC"] is None, r["acc"]
                assert exp["onoff"] in ("on_with_treatment", "off_with_treatment"), r["acc"]
            else:
                assert exp["onoff"] == "", r["acc"]
            assert exp["complete"] == (exp["log2FC"] is not None), r["acc"]
        else:
            assert r["expect"].get("detected_in") or r["expect"].get("quarantined"), r["acc"]


_validate_design()


_CORE_PROPS_NAME = "docProps/core.xml"
_MODIFIED_RE = re.compile(
    rb"(<dcterms:modified[^>]*>)[^<]*(</dcterms:modified>)"
)


def _normalize_zip(path: Path) -> None:
    """Make the .xlsx byte-reproducible: two wall clocks are hidden inside it.

    1. An .xlsx is a zip, and ``zipfile`` stamps every member with the time it
       was written. Flattened to :data:`_FIXED_TIMESTAMP`.
    2. openpyxl **overwrites** ``docProps/core.xml``'s ``dcterms:modified`` with
       ``datetime.now()`` as it saves, discarding whatever ``book.properties``
       was set to. So the pinned property alone is not enough, and the symptom
       is deceptive: consecutive regenerations agree until the wall clock ticks
       over a second, which reads as "flaky" rather than "unpinned".

    Member order and every member's content are otherwise untouched. Same idea
    as ``tools/freeze.py``'s SVG canonicalisation -- strip the clock, keep the
    data -- except that here it can be done at write time, so the committed file
    itself is stable rather than merely comparable.
    """
    import zipfile

    stamp = _FIXED_TIMESTAMP.strftime("%Y-%m-%dT%H:%M:%SZ").encode()
    with zipfile.ZipFile(path) as src:
        members = [(info, src.read(info.filename)) for info in src.infolist()]

    fixed = (
        _FIXED_TIMESTAMP.year, _FIXED_TIMESTAMP.month, _FIXED_TIMESTAMP.day,
        _FIXED_TIMESTAMP.hour, _FIXED_TIMESTAMP.minute, _FIXED_TIMESTAMP.second,
    )
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as out:
        for info, data in members:
            if info.filename == _CORE_PROPS_NAME:
                data, n = _MODIFIED_RE.subn(rb"\g<1>" + stamp + rb"\g<2>", data)
                if n != 1:
                    raise RuntimeError(
                        f"expected exactly one <dcterms:modified> in "
                        f"{_CORE_PROPS_NAME}, found {n} -- openpyxl's core "
                        "properties layout changed and the fixture is no longer "
                        "byte-reproducible"
                    )
            entry = zipfile.ZipInfo(info.filename, date_time=fixed)
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = info.external_attr
            entry.create_system = info.create_system
            out.writestr(entry, data)


def write(path: Path = OUTPUT_PATH) -> Path:
    """Write the two protein-report sheets to `path` and return it.

    Byte-reproducible: rerunning this over an unchanged :data:`ROWS` produces an
    identical file (document properties pinned, zip timestamps flattened).
    """
    df_L, df_H = build_frames()
    path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        df_L.to_excel(writer, sheet_name=SHEET_L, index=False)
        df_H.to_excel(writer, sheet_name=SHEET_H, index=False)
        props = writer.book.properties
        props.created = _FIXED_TIMESTAMP
        props.modified = _FIXED_TIMESTAMP
        props.creator = "proteomics_de/tests/fixtures/make_mini_sheet.py"
        props.lastModifiedBy = props.creator
    _normalize_zip(path)
    return path


def main() -> int:
    path = write()
    counts = expected_counts()
    print(f"wrote {path}")
    print(f"  {SHEET_L}: {counts['sheet_L_rows']} rows")
    print(f"  {SHEET_H}: {counts['sheet_H_rows']} rows")
    width = max(len(k) for k in counts)
    for key, value in counts.items():
        print(f"  {key:<{width}} : {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
