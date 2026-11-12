"""Dr. Walker's additions to the regulated list (DECISIONS_LOG D20).

What this module builds
-----------------------
``ipa_input_walker.csv`` / ``ipa_input_walker.txt``
    The core regulated list (every row of ``ipa_input_full.csv``, verbatim) plus
    three kinds of addition, in this order::

        core (measured)  ->  measured additions  ->  assigned  ->  recovered

    Columns: ``UniProt Accession Number, Gene names, log2FC, regulated,
    p_value, adj_p_value, tier, value_source``.

``walker_additions_manifest.csv``
    One row per protein ADDED to the core list: where it came from, why it was
    added, and whether its value is ``measured`` or ``assigned``.

Everything is driven by the ``walker_additions`` section of
``config/config.yaml``; this file contains no thresholds, factors or
accessions.

ASSIGNED values are not measurements
------------------------------------
Five control-only proteins (detected in both control readings, absent from both
testosterone readings) cannot have a fold change computed from the data. By
decision of Dr. Walker each is given a fixed **100-fold decrease**
(``fold_change = 0.01``, ``log2FC = log2(0.01) = -6.643856``). That number is a
labelled placeholder, not a result:

* ``value_source`` is ``assigned`` on these rows (``measured`` everywhere else),
  so the distinction survives any later filtering of the file;
* ``p_value`` and ``adj_p_value`` are left BLANK -- a test statistic for a value
  nobody measured would be invented;
* "not detected" may mean "below the detection limit", so the true change is
  unknown and may be smaller or larger than 100-fold.

This is the one place the repo knowingly writes a non-measured number into a
fold-change column. D18 rejected a sentinel +/-1 log2FC for the qualitative
lists; D20 supersedes that for these five proteins only, on the professor's
explicit decision, with the labelling above as the price.

The "recovered" rule (strict again since D24)
---------------------------------------------
D17/D18 recovered 862 proteins (248 tier-3 partial + 614 tier-1/2 qualitative).
Dr. Walker asked for only those altered by at least 50% with "no problems". The
user first chose a strict rule (``|log2FC| >= 0.585`` AND seen in all four
readings AND ``n_imputed == 0``), which admits NOTHING: every tier-3 row is in
that tier because a reading is missing, and no qualitative row has a fold change.
D20 WIDENED it (a partial-tier row with a real log2FC, ``|log2FC| >= 0.585``, not
already listed), which admitted 247 proteins that each carry an imputed replicate.
Dr. Walker's final answer (D24) is to DROP those 247, so the two "no problems"
switches (``require_all_replicates_seen``, ``require_no_imputed_values``) are
``true`` again: the strict rule, which admits zero. The machinery is kept -- set
both switches to ``false`` and the 247 return; the recovered rows would then rest
on the samples that were measured, with the missing value imputed, and carry their
real limma ``p_value``/``adj_p_value`` unchanged. :func:`apply_recovered_rule`
reports the rows surviving each active filter.

``ipa_input_walker_without_assigned.csv`` / ``.txt``
    The same list MINUS the ``value_source == assigned`` rows (D24): the dataset
    for the IPA run "without the five assigned values", a sensitivity check on how
    much the assigned -6.643856 values drive the IPA result. Same columns, same row
    order, every remaining row byte-identical to its row in ``ipa_input_walker.csv``.

Additive only
-------------
Nothing here rewrites an existing result; :data:`BYTE_FROZEN_INPUTS` is re-hashed
around every build, as ``ipa_extended.py`` does.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

# proteomics_de/export/walker_additions.py -> export -> proteomics_de -> repo root
_HERE = Path(__file__).resolve().parent
_PKG_DIR = _HERE.parent
_ROOT = _PKG_DIR.parent

if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from proteomics_de.config import constants, design  # noqa: E402
from proteomics_de.export.ipa_export import (  # noqa: E402
    IPA_FULL_COLS,
    IPA_ID_COL,
    REGULATED_LABELS,
    validate_ipa,
)
from proteomics_de.export.ipa_extended import (  # noqa: E402
    BYTE_FROZEN_INPUTS as _EXTENDED_FROZEN_INPUTS,
    PARTIAL_SOURCES,
    assert_core_rows_are_verbatim,
    assert_twins_agree,
)

# ---------------------------------------------------------------------------
# Names and vocabulary
# ---------------------------------------------------------------------------
WALKER_CSV_FILENAME = "ipa_input_walker.csv"
WALKER_TXT_FILENAME = "ipa_input_walker.txt"
WALKER_MANIFEST_FILENAME = "walker_additions_manifest.csv"
#: The list without the ASSIGNED rows (D24): the IPA "without the five" run.
WALKER_NOASSIGNED_CSV_FILENAME = "ipa_input_walker_without_assigned.csv"
WALKER_NOASSIGNED_TXT_FILENAME = "ipa_input_walker_without_assigned.txt"

#: ``value_source`` vocabulary. ``assigned`` is the label that must never be
#: dropped: it is the only thing telling a reader the number was not measured.
VALUE_MEASURED = "measured"
VALUE_ASSIGNED = "assigned"
VALUE_SOURCES = (VALUE_MEASURED, VALUE_ASSIGNED)

#: ``tier`` vocabulary. ``core`` / ``partial`` match ipa_input_extended.csv;
#: ``qualitative`` is the tier-1/2 (control-only) source of the assigned rows.
#: A recovered protein keeps the tier of the file it came from (``partial``);
#: the manifest's ``source``/``reason`` columns say it was added by the rule.
TIER_CORE = "core"
TIER_PARTIAL = "partial"
TIER_QUALITATIVE = "qualitative"

#: Columns of ``ipa_input_walker.*``. ``IPA_FULL_COLS`` stays a literal prefix so
#: the core block is a prefix-projection of ``ipa_input_full.csv``.
WALKER_COLS = list(IPA_FULL_COLS) + ["tier", "value_source"]

#: Columns of ``walker_additions_manifest.csv``.
MANIFEST_COLS = [
    IPA_ID_COL, "Gene names", "log2FC", "regulated", "tier", "value_source",
    "n_imputed", "source", "reason",
]

#: Source tier -> output ``tier`` for a recovered protein. Only ``partial`` rows
#: can ever pass (qualitative rows have no fold change), but the mapping is
#: explicit so a row never silently takes the wrong label.
_SOURCE_TO_TIER = {"partial": TIER_PARTIAL, "qualitative": TIER_QUALITATIVE}

#: Only ``log2FC`` is a measurement column for IPA's NA check: assigned rows
#: carry no p-value by design, so the default (log2FC, p, adj p) check would
#: reject them.
IPA_MEASUREMENT_COLS = ["log2FC"]

#: Frozen-count keys this module's outputs are pinned to.
COUNT_KEYS = {
    "core": "ipa_input_rows",
    "measured": "n_walker_measured_added",
    "assigned": "n_walker_assigned",
    "recovered": "n_walker_recovered",
    "total": "n_ipa_walker",
    "additions_listed": "n_walker_additions_listed",
    "without_assigned": "n_ipa_walker_without_assigned",
}

#: Files read by this module and never written. Re-hashed around every build.
BYTE_FROZEN_INPUTS = tuple(_EXTENDED_FROZEN_INPUTS) + (
    "ipa_input_extended.csv",
    "ipa_input_extended.txt",
    "qc_limma.csv",
    "qualitative_changes.csv",
    "regulated_up_partial.csv",
    "regulated_down_partial.csv",
    "onoff_proteins.csv",
)

_REQUIRED_CONFIG = {
    "assigned_fold_change": ("decrease_factor", "log2fc_decimals", "direction",
                             "expected_detected_in"),
    "recovered_rule": ("min_abs_log2fc", "require_partial_tier",
                       "require_all_replicates_seen", "require_no_imputed_values"),
}


class WalkerAdditionsError(AssertionError):
    """A Walker-additions input or output violates its contract.

    Subclasses ``AssertionError`` on purpose, matching ``IpaExtendedError`` and
    ``RegulatedListsError``.
    """


# ---------------------------------------------------------------------------
# Config / IO
# ---------------------------------------------------------------------------
def _frozen_counts(path=None):
    path = Path(path) if path else _PKG_DIR / "tests" / "expected" / "frozen_counts.json"
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _read(path):
    """Read a committed CSV at maximal float precision. No exceptions to this."""
    return pd.read_csv(path, float_precision="round_trip")


def load_config(path=None):
    """The ``walker_additions`` section of ``config.yaml``, checked for shape."""
    path = Path(path) if path else _PKG_DIR / "config" / "config.yaml"
    with open(path, encoding="utf-8") as fh:
        full = yaml.safe_load(fh)
    cfg = (full or {}).get("walker_additions")
    if not isinstance(cfg, dict):
        raise WalkerAdditionsError(f"{path}: no 'walker_additions' section")
    for section, keys in _REQUIRED_CONFIG.items():
        block = cfg.get(section)
        if not isinstance(block, dict):
            raise WalkerAdditionsError(f"walker_additions.{section} is missing")
        absent = [k for k in keys if k not in block]
        if absent:
            raise WalkerAdditionsError(f"walker_additions.{section} lacks {absent}")
    for section in ("assigned_proteins", "measured_additions"):
        entries = cfg.get(section)
        if not isinstance(entries, list) or not entries:
            raise WalkerAdditionsError(f"walker_additions.{section} must be a non-empty list")
        for e in entries:
            if not (isinstance(e, dict) and e.get("accession") and e.get("gene")):
                raise WalkerAdditionsError(
                    f"walker_additions.{section}: each entry needs accession + gene, got {e!r}")
    factor = cfg["assigned_fold_change"]["decrease_factor"]
    if not (isinstance(factor, (int, float)) and factor > 1):
        raise WalkerAdditionsError(f"decrease_factor must be a number > 1, got {factor!r}")
    if cfg["assigned_fold_change"]["direction"] not in REGULATED_LABELS:
        raise WalkerAdditionsError("assigned_fold_change.direction must be UP or DOWN")
    return cfg


def assigned_log2fc(cfg):
    """``log2(1/decrease_factor)`` rounded to the configured decimals.

    The assigned fold change is ``1/decrease_factor`` (a decrease), so its log2
    is ``-log2(decrease_factor)``: -6.643856 for a 100-fold decrease.
    """
    block = cfg["assigned_fold_change"]
    return round(-math.log2(block["decrease_factor"]), int(block["log2fc_decimals"]))


def assigned_fold_change(cfg):
    """The linear fold change that is assigned: ``1/decrease_factor`` (0.01)."""
    return 1.0 / cfg["assigned_fold_change"]["decrease_factor"]


# ---------------------------------------------------------------------------
# Measured additions (FRZB)
# ---------------------------------------------------------------------------
def build_measured_rows(cfg, partial_frames):
    """Rows for ``measured_additions``, copied from the partial-tier files.

    `partial_frames` is an iterable of ``(frame, label)`` as for
    :data:`ipa_extended.PARTIAL_SOURCES`. The label comes from WHICH FILE the row
    is in, never from re-thresholding here. No value is changed: this is a lookup.
    """
    rows = []
    for entry in cfg["measured_additions"]:
        acc, gene = entry["accession"], entry["gene"]
        hits = [(f.loc[f[IPA_ID_COL] == acc].iloc[0], label)
                for f, label in partial_frames if (f[IPA_ID_COL] == acc).any()]
        if len(hits) != 1:
            raise WalkerAdditionsError(
                f"{acc} ({gene}) is expected in exactly one partial-tier file, found {len(hits)}")
        src, label = hits[0]
        if str(src["Gene names"]) != gene:
            raise WalkerAdditionsError(
                f"{acc}: config says {gene!r} but the data says {src['Gene names']!r}")
        rows.append({
            IPA_ID_COL: acc,
            "Gene names": src["Gene names"],
            "log2FC": src["log2FC"],
            "regulated": label,
            "p_value": src["p_value"],
            "adj_p_value": src["adj_p_value"],
            "tier": TIER_PARTIAL,
            "value_source": VALUE_MEASURED,
            "n_imputed": int(src["n_imputed"]),
            "source": "regulated_up_partial.csv" if label == "UP" else "regulated_down_partial.csv",
            "reason": (f"Dr. Walker: add {gene}. Real limma estimate "
                       f"(n_imputed={int(src['n_imputed'])}); the assigned-value rule "
                       "is NOT applied."),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Assigned rows (EGF, EREG, ...)
# ---------------------------------------------------------------------------
def verify_control_only(accession, gene, single, qualitative, control_cols, treated_cols,
                        expected_detected_in):
    """Prove one protein is control-only in the data, or raise.

    Control-only = detected (finite, > 0) in BOTH control readings and in NEITHER
    testosterone reading (blank, or zero). Checked against the data rather than
    trusted from the config, because the assigned value is only defensible for
    proteins that genuinely have no testosterone-side measurement.

    Returns the ``single_condition_proteins.csv`` row.
    """
    hit = single.loc[single["accession"] == accession]
    if len(hit) != 1:
        raise WalkerAdditionsError(
            f"{accession} ({gene}): expected one row in single_condition_proteins.csv, got {len(hit)}")
    row = hit.iloc[0]
    if str(row["gene"]) != gene:
        raise WalkerAdditionsError(f"{accession}: config says {gene!r}, data says {row['gene']!r}")
    if row["detected_in"] != expected_detected_in:
        raise WalkerAdditionsError(
            f"{accession} ({gene}) is detected_in={row['detected_in']!r}, not "
            f"{expected_detected_in!r}; it is NOT control-only")
    ctrl = [row[c] for c in control_cols]
    trt = [row[c] for c in treated_cols]
    if not all(pd.notna(v) and v > 0 for v in ctrl):
        raise WalkerAdditionsError(
            f"{accession} ({gene}) is not detected in both control readings: {ctrl}")
    if any(pd.notna(v) and v > 0 for v in trt):
        raise WalkerAdditionsError(
            f"{accession} ({gene}) has a testosterone-side reading: {trt}")
    q = qualitative.loc[qualitative[IPA_ID_COL] == accession]
    if len(q) != 1 or q.iloc[0]["direction"] != "DOWN":
        raise WalkerAdditionsError(
            f"{accession} ({gene}) is not a DOWN row of qualitative_changes.csv")
    return row


def build_assigned_rows(cfg, single, qualitative, control_cols, treated_cols):
    """Rows for ``assigned_proteins``: fixed fold change, blank p-values."""
    block = cfg["assigned_fold_change"]
    lfc = assigned_log2fc(cfg)
    direction = block["direction"]
    if (direction == "DOWN") != (lfc < 0):
        raise WalkerAdditionsError(
            f"assigned log2FC {lfc} is inconsistent with direction {direction}")
    # A row is only a "regulated" call if it clears the pipeline's own cutoff;
    # a small factor (say 1.2, log2FC -0.26) would otherwise be written as a
    # "regulated DOWN" row that the pipeline itself would call NO CHANGE.
    if abs(lfc) < constants.LOG2_THRESHOLD:
        raise WalkerAdditionsError(
            f"assigned |log2FC| {abs(lfc)} is below the regulation cutoff "
            f"{constants.LOG2_THRESHOLD}; decrease_factor "
            f"{block['decrease_factor']!r} cannot yield a 'regulated' row")
    factor = block["decrease_factor"]
    rows = []
    for entry in cfg["assigned_proteins"]:
        acc, gene = entry["accession"], entry["gene"]
        src = verify_control_only(acc, gene, single, qualitative, control_cols,
                                  treated_cols, block["expected_detected_in"])
        rows.append({
            IPA_ID_COL: acc,
            "Gene names": gene,
            "log2FC": lfc,
            "regulated": direction,
            "p_value": np.nan,
            "adj_p_value": np.nan,
            "tier": TIER_QUALITATIVE,
            "value_source": VALUE_ASSIGNED,
            "n_imputed": pd.NA,
            "source": "single_condition_proteins.csv",
            "reason": (f"Dr. Walker: control-only protein ({src['detected_in']}); ASSIGNED "
                       f"{factor:g}-fold decrease, NOT measured. 'Not detected' may mean "
                       "below the detection limit."),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Recovered proteins and the strict rule
# ---------------------------------------------------------------------------
POOL_COLS = [IPA_ID_COL, "Gene names", "log2FC", "regulated", "p_value",
             "adj_p_value", "n_imputed", "n_seen", "source"]


def build_recovered_pool(partial_frames, limma, qualitative, replicate_cols):
    """The D17/D18 "recovered" set as one frame with the columns the rule reads.

    * tier-3 partial rows: ``log2FC``/``p``/``adj p``/``n_imputed`` from the partial
      files; ``n_seen`` = how many of the 4 raw intensities (from ``qc_limma.csv``)
      are present and > 0.
    * tier-1/2 qualitative rows: no fold change, no ``n_imputed`` (both NaN); only
      ``n_seen`` is computed, from the row's own intensity columns.

    The qualitative rows therefore fail the rule at the first filter that needs a
    value, which is the honest outcome: they have none.
    """
    blocks = []
    limma_idx = limma.set_index("id")
    for frame, label in partial_frames:
        ids = frame[IPA_ID_COL]
        missing = [i for i in ids if i not in limma_idx.index]
        if missing:
            raise WalkerAdditionsError(f"partial rows absent from qc_limma.csv: {missing[:5]}")
        raw = limma_idx.loc[ids, replicate_cols]
        seen = (raw.fillna(0).to_numpy() > 0).sum(axis=1)
        blocks.append(pd.DataFrame({
            IPA_ID_COL: ids.to_numpy(),
            "Gene names": frame["Gene names"].to_numpy(),
            "log2FC": frame["log2FC"].to_numpy(),
            "regulated": label,
            "p_value": frame["p_value"].to_numpy(),
            "adj_p_value": frame["adj_p_value"].to_numpy(),
            "n_imputed": frame["n_imputed"].to_numpy(),
            "n_seen": seen,
            "source": "partial",
        }))
    q_seen = (qualitative[replicate_cols].fillna(0).to_numpy() > 0).sum(axis=1)
    blocks.append(pd.DataFrame({
        IPA_ID_COL: qualitative[IPA_ID_COL].to_numpy(),
        "Gene names": qualitative["Gene names"].to_numpy(),
        "log2FC": np.nan,
        "regulated": qualitative["direction"].to_numpy(),
        "p_value": np.nan,
        "adj_p_value": np.nan,
        "n_imputed": np.nan,
        "n_seen": q_seen,
        "source": "qualitative",
    }))
    pool = pd.concat(blocks, ignore_index=True)
    if not pool[IPA_ID_COL].is_unique:
        dupes = pool.loc[pool[IPA_ID_COL].duplicated(), IPA_ID_COL].tolist()
        raise WalkerAdditionsError(f"recovered pool has duplicate accessions: {dupes[:5]}")
    return pool[POOL_COLS]


def apply_recovered_rule(pool, rule, already_listed, n_replicates):
    """Filter the recovered `pool` by the configured rule, counting every step.

    The rule (``walker_additions.recovered_rule``) has one threshold and three
    switches. As widened by the user's decision (D20) it is::

        not already listed  AND  partial tier with a real log2FC
        AND |log2FC| >= min_abs_log2fc

    with ``require_all_replicates_seen`` and ``require_no_imputed_values`` both
    ``false``. Setting either to ``true`` restores the original strict rule
    (all four readings present and > 0; ``n_imputed == 0``), which admits nothing
    on the committed data.

    Parameters
    ----------
    pool :
        :func:`build_recovered_pool` output (or a frame with :data:`POOL_COLS`).
    rule :
        ``walker_additions.recovered_rule`` from the config.
    already_listed :
        Accessions already on the list (core, measured additions, assigned); a
        recovered protein is never added twice.
    n_replicates :
        Number of raw replicates per protein (4): "seen in both replicates of both
        groups" means ``n_seen == n_replicates``.

    Returns
    -------
    (pandas.DataFrame, list[tuple[str, int]])
        The survivors, and ``[(step label, rows remaining)]`` starting with the
        pool. Only the steps that are switched on appear.
    """
    steps = [("recovered pool (partial + qualitative)", len(pool))]
    df = pool[~pool[IPA_ID_COL].isin(set(already_listed))]
    steps.append(("not already on the list", len(df)))

    if rule["require_partial_tier"]:
        df = df[(df["source"] == "partial") & df["log2FC"].notna()]
        steps.append(("partial tier with a real log2FC", len(df)))

    thr = float(rule["min_abs_log2fc"])
    lfc = df["log2FC"]
    df = df[lfc.notna() & (lfc.abs() >= thr)]
    steps.append((f"|log2FC| >= {thr:g}", len(df)))

    if rule["require_all_replicates_seen"]:
        df = df[df["n_seen"] == n_replicates]
        steps.append((f"seen in all {n_replicates} readings (both replicates of both groups)",
                      len(df)))

    if rule["require_no_imputed_values"]:
        df = df[df["n_imputed"].notna() & (df["n_imputed"] == 0)]
        steps.append(("n_imputed == 0", len(df)))
    return df.reset_index(drop=True), steps


_RECOVERED_COLS = [IPA_ID_COL, "Gene names", "log2FC", "regulated", "p_value",
                   "adj_p_value", "tier", "value_source", "n_imputed", "source", "reason"]

#: Why a recovered row is on the list. The wording matters: Dr. Walker asked for
#: proteins with "no problems", and these have one, so the reason says so.
_RECOVERED_REASON = (
    "Added by the recovered-protein rule, WIDENED by the user's decision "
    "(Dr. Walker said 'no problems'; this protein has one: a missing replicate "
    "imputed by the pipeline, n_imputed={n}). Real limma log2FC, p_value and "
    "adj_p_value copied from {src}; rests on the samples that were measured, "
    "with the missing value imputed by the pipeline.")


def build_recovered_rows(passed):
    """Turn rule survivors into walker rows (measured, tier of the source file).

    Sorted by descending ``|log2FC|`` with a stable sort, so ties keep their
    source-file order and a re-run reproduces identical bytes.
    """
    if passed.empty:
        return pd.DataFrame(columns=_RECOVERED_COLS)
    out = passed.copy()
    unknown = set(out["source"]) - set(_SOURCE_TO_TIER)
    if unknown:
        raise WalkerAdditionsError(f"unknown recovered source tier(s): {sorted(unknown)}")
    out["tier"] = out["source"].map(_SOURCE_TO_TIER)
    out["value_source"] = VALUE_MEASURED
    file_names = dict(zip(("UP", "DOWN"), (n for n, _l in PARTIAL_SOURCES)))
    srcs = [file_names.get(r, "qualitative_changes.csv") for r in out["regulated"]]
    out["reason"] = [
        _RECOVERED_REASON.format(n=int(n), src=src) if n > 0 else
        f"Added by the recovered-protein rule (strict conditions satisfied: no imputed "
        f"values). Real limma values copied from {src}."
        for n, src in zip(out["n_imputed"], srcs)]
    out["source"] = srcs
    out["n_imputed"] = out["n_imputed"].astype("Int64")
    out = (out.assign(_abs=out["log2FC"].abs())
              .sort_values("_abs", ascending=False, kind="mergesort")
              .drop(columns="_abs").reset_index(drop=True))
    return out[_RECOVERED_COLS]


# ---------------------------------------------------------------------------
# The assembled frame
# ---------------------------------------------------------------------------
def build_walker_frame(df_full, additions):
    """Stack the core rows and the ordered `additions` into the Walker frame.

    Core rows come first, in exactly the order ``ipa_input_full.csv`` has them
    (so a reviewer diffing sees the core block unchanged, then the additions).
    Invariants are checked here, before anything is written.
    """
    missing = [c for c in IPA_FULL_COLS if c not in df_full.columns]
    if missing:
        raise WalkerAdditionsError(f"ipa_input_full frame is missing columns: {missing}")
    core = df_full[list(IPA_FULL_COLS)].copy()
    core["tier"] = TIER_CORE
    core["value_source"] = VALUE_MEASURED

    blocks = [core[WALKER_COLS]]
    for add in additions:
        if len(add):
            blocks.append(add[WALKER_COLS])
    frame = pd.concat(blocks, ignore_index=True)
    check_walker_frame(frame)
    return frame[WALKER_COLS]


def check_walker_frame(frame):
    """Structural invariants of a Walker frame; raises on the first violation."""
    if list(frame.columns) != WALKER_COLS:
        raise WalkerAdditionsError(f"columns are {list(frame.columns)}, expected {WALKER_COLS}")
    if not frame[IPA_ID_COL].is_unique:
        dupes = frame.loc[frame[IPA_ID_COL].duplicated(), IPA_ID_COL].tolist()
        raise WalkerAdditionsError(f"duplicate accessions: {dupes[:10]}")
    bad_labels = set(frame["regulated"].unique()) - set(REGULATED_LABELS)
    if bad_labels:
        raise WalkerAdditionsError(f"unexpected 'regulated' values: {sorted(bad_labels)}")
    bad_src = set(frame["value_source"].unique()) - set(VALUE_SOURCES)
    if bad_src:
        raise WalkerAdditionsError(f"unexpected 'value_source' values: {sorted(bad_src)}")
    lfc = pd.to_numeric(frame["log2FC"], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(lfc).all():
        raise WalkerAdditionsError("non-finite log2FC; IPA silently drops or misreads the row")
    up = frame["regulated"].to_numpy() == "UP"
    if ((lfc > 0) != up).any() or (lfc == 0).any():
        bad = frame.loc[(lfc > 0) != up, IPA_ID_COL].tolist()
        raise WalkerAdditionsError(f"'regulated' disagrees with the sign of log2FC: {bad[:10]}")
    assigned = frame["value_source"] == VALUE_ASSIGNED
    if frame.loc[assigned, ["p_value", "adj_p_value"]].notna().any().any():
        raise WalkerAdditionsError("an assigned row carries a p-value; none may be invented")
    measured = ~assigned
    if frame.loc[measured, ["p_value", "adj_p_value"]].isna().any().any():
        raise WalkerAdditionsError("a measured row is missing its p-value / adj p-value")


def drop_assigned(frame):
    """The Walker frame minus every ``value_source == assigned`` row (D24).

    Used for the IPA "without the five assigned values" run. It is a pure row
    filter: the surviving rows keep their order and values, so each one is
    byte-identical to its row in the full list, and nothing is re-derived. The
    result is checked as a Walker frame in its own right (and must contain no
    assigned row, otherwise the file would not be what its name says).
    """
    check_walker_frame(frame)
    keep = frame["value_source"] != VALUE_ASSIGNED
    out = frame.loc[keep].reset_index(drop=True)
    if (out["value_source"] == VALUE_ASSIGNED).any():   # pragma: no cover - defensive
        raise WalkerAdditionsError("an assigned row survived drop_assigned")
    if len(out) != len(frame) - int((~keep).sum()):      # pragma: no cover - defensive
        raise WalkerAdditionsError("drop_assigned lost or invented rows")
    check_walker_frame(out)
    return out


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------
def write_walker(frame, csv_path, txt_path=None):
    """Write ``ipa_input_walker.csv`` and its tab-delimited ``.txt`` twin.

    UTF-8 without a BOM, single header row, identifier leftmost, blank cells for
    the assigned rows' p-values -- the same conventions as ``write_ipa_extended``.
    ``txt_path=None`` derives the sibling; ``False`` writes only the CSV.
    """
    csv_path = Path(csv_path)
    if txt_path is None:
        txt_path = csv_path.with_suffix(".txt")
    frame[WALKER_COLS].to_csv(csv_path, index=False, encoding="utf-8")
    written = [csv_path]
    if txt_path is not False:
        txt_path = Path(txt_path)
        frame[WALKER_COLS].to_csv(txt_path, index=False, sep="\t", encoding="utf-8")
        written.append(txt_path)
    return written


def write_manifest(manifest, path):
    manifest[MANIFEST_COLS].to_csv(path, index=False, encoding="utf-8")
    return Path(path)


# ---------------------------------------------------------------------------
# Byte-freeze guard
# ---------------------------------------------------------------------------
def _hash_frozen_inputs(results_dir):
    results_dir = Path(results_dir)
    return {name: hashlib.sha256((results_dir / name).read_bytes()).hexdigest()
            for name in BYTE_FROZEN_INPUTS if (results_dir / name).exists()}


def assert_frozen_inputs_untouched(results_dir, before):
    after = _hash_frozen_inputs(results_dir)
    moved = sorted(k for k in before if before[k] != after.get(k))
    if moved:
        raise WalkerAdditionsError(
            f"this module modified a byte-frozen input, which it must never do: {moved}")


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
@dataclass
class WalkerResult:
    frame: pd.DataFrame
    manifest: pd.DataFrame
    frame_without_assigned: pd.DataFrame = None
    steps: list = field(default_factory=list)
    counts: dict = field(default_factory=dict)
    paths: list = field(default_factory=list)


_INPUTS = ("ipa_input_full.csv", "qc_limma.csv", "qualitative_changes.csv",
           "single_condition_proteins.csv", "regulated_up_partial.csv",
           "regulated_down_partial.csv")


def compute_walker(results_dir, cfg):
    """Build the Walker frame and manifest from `results_dir`; writes nothing."""
    results_dir = Path(results_dir)
    for name in _INPUTS:
        if not (results_dir / name).exists():
            raise FileNotFoundError(
                f"{results_dir / name} not found -- run run_pipeline.py --only "
                "foldchange,ipa_export,regulated_lists_supplementary first.")

    df_full = _read(results_dir / "ipa_input_full.csv")
    limma = _read(results_dir / "qc_limma.csv")
    qualitative = _read(results_dir / "qualitative_changes.csv")
    single = _read(results_dir / "single_condition_proteins.csv")
    partial_frames = [(_read(results_dir / name), label) for name, label in PARTIAL_SOURCES]

    control_cols = design.control_columns()
    treated_cols = design.treated_columns()
    replicate_cols = control_cols + treated_cols

    measured = build_measured_rows(cfg, partial_frames)
    assigned = build_assigned_rows(cfg, single, qualitative, control_cols, treated_cols)

    already = (set(df_full[IPA_ID_COL]) | set(measured[IPA_ID_COL])
               | set(assigned[IPA_ID_COL]))
    overlap = (set(measured[IPA_ID_COL]) | set(assigned[IPA_ID_COL])) & set(df_full[IPA_ID_COL])
    if overlap:
        raise WalkerAdditionsError(
            f"additions already in the core list (would duplicate): {sorted(overlap)}")
    if set(measured[IPA_ID_COL]) & set(assigned[IPA_ID_COL]):
        raise WalkerAdditionsError("a protein is both measured-added and assigned")

    pool = build_recovered_pool(partial_frames, limma, qualitative, replicate_cols)
    passed, steps = apply_recovered_rule(pool, cfg["recovered_rule"], already,
                                         n_replicates=len(replicate_cols))
    recovered = build_recovered_rows(passed)

    frame = build_walker_frame(df_full, [measured, assigned, recovered])
    frame_without_assigned = drop_assigned(frame)
    manifest = pd.concat(
        [f[MANIFEST_COLS].astype({"n_imputed": "Int64"})
         for f in (measured, assigned, recovered) if len(f)],
        ignore_index=True)
    manifest["n_imputed"] = manifest["n_imputed"].astype("Int64")
    counts = {
        "core": len(df_full),
        "measured": len(measured),
        "assigned": len(assigned),
        "recovered": len(recovered),
        "total": len(frame),
        "additions_listed": len(manifest),
        "without_assigned": len(frame_without_assigned),
        "recovered_pool": len(pool),
    }
    return WalkerResult(frame=frame, manifest=manifest,
                        frame_without_assigned=frame_without_assigned,
                        steps=steps, counts=counts)


def build_walker_from_results(results_dir=None, *, validate=True, counts_path=None,
                              config_path=None):
    """Build the Walker files from the committed results directory.

    Returns
    -------
    WalkerResult
        With ``.paths`` = ``[walker.csv, walker.txt, manifest.csv,
        without_assigned.csv, without_assigned.txt]``.
    """
    results_dir = Path(results_dir) if results_dir else _PKG_DIR / "results"
    cfg = load_config(config_path)

    before = _hash_frozen_inputs(results_dir)
    result = compute_walker(results_dir, cfg)

    # Count checks need no files, so they run before anything is written.
    counts = _frozen_counts(counts_path) if validate else None
    if validate:
        _check_counts(result, counts)

    csv_path = results_dir / WALKER_CSV_FILENAME
    txt_path = results_dir / WALKER_TXT_FILENAME
    manifest_path = results_dir / WALKER_MANIFEST_FILENAME
    na_csv_path = results_dir / WALKER_NOASSIGNED_CSV_FILENAME
    na_txt_path = results_dir / WALKER_NOASSIGNED_TXT_FILENAME

    # Stage in a scratch directory INSIDE results_dir (same filesystem, so the
    # final os.replace is an atomic rename), validate the staged files under
    # their real names (validators pick the delimiter from the extension), and
    # only then move them into place. A failed validation therefore leaves no
    # deliverable behind -- nothing is renamed -- and the scratch directory is
    # always removed.
    stage = Path(tempfile.mkdtemp(prefix=".walker_stage_", dir=results_dir))
    try:
        staged_csv, staged_txt = write_walker(
            result.frame, stage / WALKER_CSV_FILENAME, stage / WALKER_TXT_FILENAME)
        staged_manifest = write_manifest(result.manifest, stage / WALKER_MANIFEST_FILENAME)
        staged_na_csv, staged_na_txt = write_walker(
            result.frame_without_assigned, stage / WALKER_NOASSIGNED_CSV_FILENAME,
            stage / WALKER_NOASSIGNED_TXT_FILENAME)

        if validate:
            for p in (staged_csv, staged_txt):
                validate_ipa(p, expected_rows=int(counts[COUNT_KEYS["total"]]),
                             required_columns=WALKER_COLS,
                             measurement_columns=IPA_MEASUREMENT_COLS)
            assert_core_rows_are_verbatim(staged_csv, results_dir / "ipa_input_full.csv")
            assert_twins_agree(staged_csv, staged_txt)
            for p in (staged_na_csv, staged_na_txt):
                validate_ipa(p, expected_rows=int(counts[COUNT_KEYS["without_assigned"]]),
                             required_columns=WALKER_COLS,
                             measurement_columns=IPA_MEASUREMENT_COLS)
            assert_twins_agree(staged_na_csv, staged_na_txt)
            assert_without_assigned_is_the_list_minus_assigned(staged_csv, staged_na_csv)
            assert_frozen_inputs_untouched(results_dir, before)

        for staged, final in ((staged_csv, csv_path), (staged_txt, txt_path),
                              (staged_manifest, manifest_path),
                              (staged_na_csv, na_csv_path), (staged_na_txt, na_txt_path)):
            os.replace(staged, final)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
    result.paths = [csv_path, txt_path, manifest_path, na_csv_path, na_txt_path]
    return result


def assert_without_assigned_is_the_list_minus_assigned(full_csv, without_csv):
    """Re-read both written CSVs and prove the second is the first minus assigned rows.

    Compares the FILES (as text, so a formatting difference counts), not the
    in-memory frames: the guarantee is about what a reader of the two files sees.
    """
    full = pd.read_csv(full_csv, dtype=str, keep_default_na=False, float_precision="round_trip")
    wout = pd.read_csv(without_csv, dtype=str, keep_default_na=False, float_precision="round_trip")
    expected = full.loc[full["value_source"] != VALUE_ASSIGNED].reset_index(drop=True)
    if not expected.equals(wout):
        raise WalkerAdditionsError(
            f"{Path(without_csv).name} is not {Path(full_csv).name} minus its "
            "value_source=assigned rows")
    if (wout["value_source"] == VALUE_ASSIGNED).any():
        raise WalkerAdditionsError(f"{Path(without_csv).name} still has an assigned row")


def _check_counts(result, counts):
    """Pin the result to the frozen counts (raises before anything is written)."""
    for what, key in COUNT_KEYS.items():
        if result.counts[what] != int(counts[key]):
            raise WalkerAdditionsError(
                f"{what}: {result.counts[what]} rows, expected {counts[key]} "
                f"(frozen_counts {key})")
    if int(counts[COUNT_KEYS["total"]]) != sum(
            int(counts[COUNT_KEYS[k]]) for k in ("core", "measured", "assigned", "recovered")):
        raise WalkerAdditionsError("frozen_counts n_ipa_walker is not core + additions")
    if int(counts[COUNT_KEYS["additions_listed"]]) != sum(
            int(counts[COUNT_KEYS[k]]) for k in ("measured", "assigned", "recovered")):
        raise WalkerAdditionsError(
            "frozen_counts n_walker_additions_listed is not measured + assigned + recovered")
    if int(counts[COUNT_KEYS["without_assigned"]]) != int(counts[COUNT_KEYS["total"]]) - int(
            counts[COUNT_KEYS["assigned"]]):
        raise WalkerAdditionsError(
            "frozen_counts n_ipa_walker_without_assigned is not n_ipa_walker - n_walker_assigned")


def _n_imputed_counts(rows):
    """``{"1": n, "2": n, ...}`` for a manifest slice (blank n_imputed ignored)."""
    return {str(int(k)): int(v) for k, v in rows.groupby("n_imputed").size().items()}


def sidecar_facts(result, cfg, which="list"):
    """Extra provenance fields: what is assigned, and how the rule narrowed.

    `which` scopes the per-file facts to the file the sidecar describes:
    ``"list"`` (``ipa_input_walker.csv``/``.txt``: every row of the list),
    ``"manifest"`` (``walker_additions_manifest.csv``: only the additions) or
    ``"without_assigned"`` (``ipa_input_walker_without_assigned.*``: the list minus
    the assigned rows). A count taken from the wrong frame would describe a
    different file.
    """
    if which not in ("list", "manifest", "without_assigned"):
        raise WalkerAdditionsError(
            f"sidecar scope must be 'list', 'manifest' or 'without_assigned', got {which!r}")
    block = cfg["assigned_fold_change"]
    scoped = {"list": result.frame, "manifest": result.manifest,
              "without_assigned": result.frame_without_assigned}[which]
    scope_text = {
        "list": "all rows of ipa_input_walker.csv/.txt",
        "manifest": "all rows of walker_additions_manifest.csv (the additions only)",
        "without_assigned": ("all rows of ipa_input_walker_without_assigned.csv/.txt "
                             "(the list minus its value_source=assigned rows)"),
    }[which]
    measured_acc = {e["accession"] for e in cfg["measured_additions"]}
    man = result.manifest
    is_measured = man["value_source"] == VALUE_MEASURED
    is_listed_measured = man[IPA_ID_COL].isin(measured_acc) & is_measured
    rule = cfg["recovered_rule"]
    strict = bool(rule["require_all_replicates_seen"] or rule["require_no_imputed_values"])
    n_recovered = int(result.counts["recovered"])
    if strict:
        rule_note = (
            "The recovered-protein rule is STRICT (Dr. Walker's original 'no problems' "
            "wording, restored by his final answer, DECISIONS_LOG D24): a protein must be "
            "seen in all four readings with n_imputed == 0. That admits "
            f"{n_recovered} proteins on the committed data, so the recovered pool is empty; "
            "the 247 proteins the widened rule had admitted (D20) were dropped at his request.")
    else:
        rule_note = (
            "The recovered-protein rule is WIDENED (both 'no problems' switches false): "
            "partial-tier proteins with |log2FC| >= the cutoff are added although one "
            "replicate was imputed by the pipeline. They rest on the samples that were "
            "measured, with the missing value imputed. Real limma p_value / adj_p_value are "
            "kept; see the manifest's n_imputed and reason columns.")
    facts = {
        "value_source_counts_scope": scope_text,
        "value_source_counts": {k: int(v) for k, v in
                                scoped["value_source"].value_counts().items()},
        "assigned_values_are_not_measurements": True,
        "assigned_fold_change": assigned_fold_change(cfg),
        "assigned_log2FC": assigned_log2fc(cfg),
        "assigned_accessions": [e["accession"] for e in cfg["assigned_proteins"]],
        "assignment_note": (
            f"The {block['decrease_factor']:g}-fold decrease on value_source=assigned rows is "
            "ASSIGNED by Dr. Walker's decision, not measured. These proteins were detected in "
            "both control readings and in neither testosterone reading; 'not detected' may "
            "mean below the detection limit. p_value / adj_p_value are blank on those rows "
            "on purpose (DECISIONS_LOG D20)."),
        "recovered_rule": rule,
        "recovered_rule_steps": [{"step": s, "n": n} for s, n in result.steps],
        "recovered_rule_note": rule_note,
        # FRZB (a measured addition) and the rule-recovered proteins are reported
        # separately so neither count silently includes the other.
        "measured_additions_n_imputed": {
            str(a): int(n) for a, n in
            man.loc[is_listed_measured].set_index(IPA_ID_COL)["n_imputed"].items()},
        "recovered_n_imputed_counts": _n_imputed_counts(
            man.loc[is_measured & ~is_listed_measured]),
        "decision": "DECISIONS_LOG D20; D24 (Dr. Walker's final answers)",
    }
    if which == "without_assigned":
        facts["without_assigned_note"] = (
            "This file is ipa_input_walker.csv minus its value_source=assigned rows (the "
            f"{int(result.counts['assigned'])} control-only proteins whose "
            f"{assigned_log2fc(cfg)} is assigned, not measured). It exists for the IPA sensitivity run 'without the "
            "assigned values' (D24); every row is a measured row, unchanged.")
    return facts


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results-dir", default=None,
                    help="directory holding ipa_input_full.csv / qc_limma.csv / "
                         "qualitative_changes.csv / ... (default: proteomics_de/results)")
    ap.add_argument("--no-sidecars", action="store_true",
                    help="skip writing the .provenance.json sidecars")
    args = ap.parse_args(argv)

    results_dir = Path(args.results_dir) if args.results_dir else _PKG_DIR / "results"
    cfg = load_config()
    result = build_walker_from_results(results_dir)

    c = result.counts
    print(f"core {c['core']} -> +{c['measured']} measured -> +{c['assigned']} assigned "
          f"(value_source=assigned, NOT measured) -> +{c['recovered']} recovered "
          f"= {c['total']} rows")
    print(f"without the {c['assigned']} assigned rows: {c['without_assigned']} rows")
    print("recovered-protein rule, rows surviving each step:")
    for label, n in result.steps:
        print(f"  {n:>5}  {label}")
    for p in result.paths:
        print(f"Saved {p}")

    if not args.no_sidecars:
        from proteomics_de.provenance import sidecar

        for p in result.paths:
            if p.name == WALKER_MANIFEST_FILENAME:
                which = "manifest"
            elif p.name.startswith("ipa_input_walker_without_assigned"):
                which = "without_assigned"
            else:
                which = "list"
            print(f"Saved {sidecar(p, **sidecar_facts(result, cfg, which))}")
    return 0


__all__ = [
    "WALKER_CSV_FILENAME", "WALKER_TXT_FILENAME", "WALKER_MANIFEST_FILENAME",
    "WALKER_NOASSIGNED_CSV_FILENAME", "WALKER_NOASSIGNED_TXT_FILENAME",
    "VALUE_MEASURED", "VALUE_ASSIGNED", "VALUE_SOURCES",
    "TIER_CORE", "TIER_PARTIAL", "TIER_QUALITATIVE",
    "WALKER_COLS", "MANIFEST_COLS", "COUNT_KEYS", "BYTE_FROZEN_INPUTS",
    "WalkerAdditionsError", "WalkerResult",
    "load_config", "assigned_log2fc", "assigned_fold_change",
    "build_measured_rows", "verify_control_only", "build_assigned_rows",
    "build_recovered_pool", "apply_recovered_rule", "build_recovered_rows",
    "build_walker_frame", "check_walker_frame", "drop_assigned",
    "assert_without_assigned_is_the_list_minus_assigned",
    "write_walker", "write_manifest",
    "assert_frozen_inputs_untouched",
    "compute_walker", "build_walker_from_results", "sidecar_facts", "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
