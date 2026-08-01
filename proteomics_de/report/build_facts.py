#!/usr/bin/env python3
"""build_facts.py — regenerate ``report/report_facts.json`` from committed artifacts.

WHY THIS EXISTS
---------------
``report_facts.json`` used to be hand-authored. Nothing produced it, so a full
pipeline re-run refreshed the report's *figures* while its *numbers* silently
rotted: the committed file still carried the pre-D7 orientation (206 UP / 509
DOWN), 606 single-condition rows, a 2554-protein background and vanilla
p-values, long after DECISIONS_LOG D7/D9/D10/D11 had moved every one of them.

This script closes that gap. Every value it emits is read or recomputed from a
committed artifact under ``proteomics_de/results/`` (or from
``config/sample_sheet.tsv``). There are no numeric literals describing the
dataset anywhere below -- if the pipeline is re-run, re-running this script
updates the report's numbers to match.

STRICT JSON
-----------
The previous hand-authored file contained a bare ``NaN`` literal on the
``corr_limma_vs_pipeline_log2FC`` line. Python's ``json.load`` tolerates that
extension; ``jq``, JavaScript ``JSON.parse``, Go ``encoding/json`` and Rust
``serde_json`` all reject it, so the digest was unreadable outside Python. This
script writes with ``allow_nan=False``, which turns any stray non-finite value
into a hard failure at build time rather than a corrupt artifact.

That NaN was not a real result. It came from ``np.corrcoef`` propagating NaNs
from the 360 both-condition rows that are *partial* (limma imputes them, so they
have a ``limma_log2FC``, but they have no complete-case pipeline ``log2FC``).
Computed on the paired, fully-observed subset the correlation is +1.0 -- see
:func:`limma_block`.

Usage::

    python3 build_facts.py            # rewrites report_facts.json in place
    python3 build_facts.py --check    # verify the committed file is current
    python3 build_facts.py -o out.json
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent                  # proteomics_de/report
PKG = HERE.parent                                       # proteomics_de
REPO = PKG.parent                                       # repo root
RESULTS = PKG / "results"
FIGDIR = RESULTS / "figures"
ENRICH = RESULTS / "enrichment"
GATED = RESULTS / "gated"
OUTPUT = HERE / "report_facts.json"

# config.design is the single source of truth for the design; import it the way
# the rest of the pipeline does (repo root on sys.path).
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from proteomics_de.config import design as design_cfg  # noqa: E402

#: The four per-figure manifests, in the order their entries are concatenated.
#: The union of these IS the ``figures`` array -- titles, captions and
#: key_numbers are copied verbatim, never re-derived here, so a figure's caption
#: can only ever be changed by re-running the layer that draws it.
MANIFESTS = [
    "figures_manifest.json",          # viz core          -> 7 entries
    "figures_manifest_enrich.json",   # enrichment layer  -> 3 entries
    "figures_manifest_gated.json",    # gated PCA layer   -> 2 entries
    "figures_manifest_network.json",  # STRING network    -> 1 entry
]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def jsonable(obj):
    """Convert numpy/pandas scalars to plain Python; NaN -> None.

    ``allow_nan=False`` at dump time rejects non-finite floats outright, so any
    genuinely-missing cell must become ``null`` here (a real JSON value) rather
    than riding through as ``NaN``.
    """
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        obj = float(obj)
    if isinstance(obj, float):
        return None if not math.isfinite(obj) else obj
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if obj is None or isinstance(obj, (str, int, bool)):
        return obj
    if isinstance(obj, np.ndarray):
        return jsonable(obj.tolist())
    if obj is pd.NaT or (isinstance(obj, float) and pd.isna(obj)):
        return None
    return obj


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def one_row(path: Path) -> dict:
    """Read a single-row CSV into a dict (the QC check-record format)."""
    df = pd.read_csv(path)
    if len(df) != 1:
        raise SystemExit(f"{path} should hold exactly 1 row, has {len(df)}")
    return jsonable(df.iloc[0].to_dict())


def records(path: Path) -> list[dict]:
    return jsonable(pd.read_csv(path).to_dict(orient="records"))


# ---------------------------------------------------------------------------
# blocks
# ---------------------------------------------------------------------------
def dataset_block() -> dict:
    """Design metadata, resolved from the sample sheet and the STRING meta.

    The control/treated assignment is read from ``config/sample_sheet.tsv``,
    which is the file DECISIONS_LOG D7 corrected. Nothing here restates it.

    The isotope channel is deliberately NOT presented as the contrast. Per
    ``foldchange.py`` (``SHEET_L_COLS`` / ``SHEET_H_COLS``, and the explicit
    "sheet-H channels: dtype, not condition" comment) the workbook's "Protein
    Report L" / "Protein Report H" sheets are where a channel's row LIVES, not
    which condition it belongs to; and DECISIONS_LOG D8 records that all four
    samples are complete SILAC experiments carrying their own L/H columns.
    """
    sheet = design_cfg.read_sample_sheet()
    string_meta = read_json(ENRICH / "string_meta.json")

    control_ids = sheet.loc[sheet["group"] == "control", "sample"].tolist()
    treated_ids = sheet.loc[sheet["group"] == "treated", "sample"].tolist()
    n_rep = design_cfg.replicates_per_group()

    return {
        "organism": f"{string_meta['species_name']} (taxid {string_meta['species']})",
        "assay": "SILAC (Stable Isotope Labeling by Amino acids in Cell culture)",
        "design": (
            f"testosterone-treated ({', '.join(treated_ids)}) vs "
            f"vehicle control ({', '.join(control_ids)})"
        ),
        "control_samples": control_ids,
        "treated_samples": treated_ids,
        "n_samples": design_cfg.n_samples(),
        "replicates_per_group": n_rep,
        "replicates": (
            f"n={n_rep} TECHNICAL replicates per group (no biological replication)"
        ),
        "design_source": "config/sample_sheet.tsv",
        "channel_note": (
            "The workbook's 'Protein Report L' / 'Protein Report H' sheets are "
            "row residency, not condition: 31578/31580 live in sheet L and "
            "31579/31581 in sheet H, but which of those is control is set by "
            "the sample sheet alone. Every sample is itself a complete SILAC "
            "run; the pipeline analyses the combined Intensity and never the "
            "native H/L ratio (DECISIONS_LOG D8, open)."
        ),
    }


def counts_block() -> dict:
    """Row and class counts, straight off the committed CSVs."""
    fa = pd.read_csv(RESULTS / "foldchange_all.csv")
    sc = pd.read_csv(RESULTS / "single_condition_proteins.csv")
    oo = pd.read_csv(RESULTS / "onoff_proteins.csv")
    vc = fa["regulated"].value_counts()

    n_up = int(vc.get("UP", 0))
    n_down = int(vc.get("DOWN", 0))
    return {
        "both_condition": len(fa),
        "single_condition": len(sc),
        "on_off": len(oo),
        # The enrichment background: every protein seen in at least one
        # condition. D11 quarantined 2 junk accessions out of single_condition,
        # so this shrank 2554 -> 2552 without any other count moving.
        "detected_universe": len(fa) + len(sc),
        "regulated_UP": n_up,
        "regulated_DOWN": n_down,
        "no_change": int(vc.get("NO CHANGE", 0)),
        "regulated_total": n_up + n_down,
    }


def qc_block() -> dict:
    """The two standalone QC check-records.

    ``replicate`` is copied 1:1 from ``qc_replicate_correlation.csv``.
    ``replicate_raw_r_by_group`` is recomputed here from the sample sheet
    because that CSV's control/treated LABELS are a known pre-D7 residue:
    ``replicate_check.py`` hardcodes its channel pairs and never reads the
    sheet, so its two raw-r values are correct numbers under swapped names.
    Both are emitted so the discrepancy is visible rather than papered over.
    """
    centering = one_row(RESULTS / "qc_centering.csv")
    replicate = one_row(RESULTS / "qc_replicate_correlation.csv")

    fa = pd.read_csv(RESULTS / "foldchange_all.csv")
    by_group = {}
    for group in ("control", "treated"):
        cols = design_cfg.columns_for_group(group)
        if len(cols) != 2:
            raise SystemExit(
                f"qc_block assumes 2 replicates per group, got {len(cols)} for {group}"
            )
        x = pd.to_numeric(fa[cols[0]], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(fa[cols[1]], errors="coerce").to_numpy(dtype=float)
        mask = (x > 0) & (y > 0)
        by_group[group] = {
            "channels": cols,
            "raw_log10_pearson_r": round(
                float(np.corrcoef(np.log10(x[mask]), np.log10(y[mask]))[0, 1]), 4
            ),
            "n": int(mask.sum()),
        }

    return {
        "centering": centering,
        "replicate": replicate,
        "replicate_raw_r_by_group": by_group,
        "replicate_label_drift": (
            "qc_replicate_correlation.csv labels control_raw_r/treated_raw_r "
            "using replicate_check.py's hardcoded pre-D7 channel pairs, which "
            "are inverted relative to config/sample_sheet.tsv. The VALUES are "
            "right, the two names are swapped; replicate_raw_r_by_group above "
            "is recomputed against the corrected sheet."
        ),
    }


def limma_block() -> dict:
    """limma summary. Under D9 the DEFAULT is eBayes(trend=TRUE, robust=TRUE).

    ``qc_limma.csv`` is that default run; ``qc_limma_vanilla.csv`` is the
    plain-eBayes comparison kept for continuity. The unprefixed keys therefore
    describe trend/robust, and the vanilla numbers carry a ``vanilla_`` prefix
    -- the reverse of the pre-D9 file, where vanilla was the baseline.
    """
    trend = pd.read_csv(RESULTS / "qc_limma.csv")
    vanilla = pd.read_csv(RESULTS / "qc_limma_vanilla.csv")

    if len(trend) != len(vanilla):
        raise SystemExit(
            f"limma flavours disagree on row count: trend={len(trend)} "
            f"vanilla={len(vanilla)}"
        )

    def summarise(df):
        return {
            "n_sig_fdr005": int((df["adj_p_value"] < 0.05).sum()),
            "min_raw_p": float(df["p_value"].min()),
            "n_raw_p_lt05": int((df["p_value"] < 0.05).sum()),
            "min_adj_p": float(df["adj_p_value"].min()),
        }

    t, v = summarise(trend), summarise(vanilla)

    # corr(limma log2FC, pipeline log2FC) on the paired, fully-observed subset.
    # The pipeline log2FC is complete-case only, so the 360 partial rows are
    # NaN there; np.corrcoef over the raw columns propagates that to NaN, which
    # is where the old hand-authored `NaN` literal came from. It is not a
    # result -- restricted to rows where both are finite the answer is +1.0.
    fa = pd.read_csv(RESULTS / "foldchange_all.csv")
    merged = trend.merge(
        fa[["UniProt Accession Number", "log2FC"]],
        left_on="id",
        right_on="UniProt Accession Number",
        how="inner",
    )
    a = pd.to_numeric(merged["limma_log2FC"], errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(merged["log2FC"], errors="coerce").to_numpy(dtype=float)
    paired = np.isfinite(a) & np.isfinite(b)
    if paired.sum() < 2:
        raise SystemExit("fewer than 2 paired finite log2FC rows; cannot correlate")
    corr = float(np.corrcoef(a[paired], b[paired])[0, 1])
    if not math.isfinite(corr):
        raise SystemExit("correlation is non-finite after masking; refusing to write")

    return {
        "ebayes_mode": "trend=TRUE, robust=TRUE (DECISIONS_LOG D9 default)",
        "n_tested": len(trend),
        "n_sig_fdr005": t["n_sig_fdr005"],
        "min_raw_p": t["min_raw_p"],
        "n_raw_p_lt05": t["n_raw_p_lt05"],
        "min_adj_p": t["min_adj_p"],
        "vanilla_n_sig_fdr005": v["n_sig_fdr005"],
        "vanilla_min_raw_p": v["min_raw_p"],
        "vanilla_n_raw_p_lt05": v["n_raw_p_lt05"],
        "vanilla_min_adj_p": v["min_adj_p"],
        "corr_limma_vs_pipeline_log2FC": corr,
        "corr_n_paired": int(paired.sum()),
        "corr_n_unpaired_partial_rows": int((~paired).sum()),
    }


def top_candidates_block(n: int = 10) -> list[dict]:
    """The n smallest raw-p proteins from the DEFAULT (trend/robust) limma run.

    This is the same ranking the volcano labels and the heatmap selects on, so
    the table in the report and the figures agree by construction.
    """
    df = pd.read_csv(RESULTS / "qc_limma.csv").nsmallest(n, "p_value")
    return [
        {
            "gene": jsonable(r["gene"]),
            "acc": jsonable(r["id"]),
            "log2FC": round(float(r["limma_log2FC"]), 3),
            "p": float(r["p_value"]),
            "regulated": jsonable(r["regulated"]),
            "n_imputed": jsonable(r["n_imputed"]),  # D10
        }
        for _, r in df.iterrows()
    ]


def figures_block() -> list[dict]:
    """The 13 figure entries: the verbatim union of the four manifests."""
    out = []
    for name in MANIFESTS:
        entries = read_json(FIGDIR / name)
        if not isinstance(entries, list):
            raise SystemExit(f"{name} should be a JSON list, got {type(entries).__name__}")
        out.extend(entries)
    return out


def build() -> dict:
    facts = {
        "dataset": dataset_block(),
        "counts": counts_block(),
        "qc": qc_block(),
        "limma": limma_block(),
        "top_candidates": top_candidates_block(),
        # The three enrichment metas are copied whole: they are already the
        # layer's own machine-written record, so restating them would just add
        # a second place to go stale.
        "string": read_json(ENRICH / "string_meta.json"),
        "ora": read_json(ENRICH / "ora_meta.json"),
        "gsea": read_json(ENRICH / "gsea_meta.json"),
        "pca_variance": records(GATED / "pca_variance.csv"),
        "gate_skips": records(GATED / "skip_log.csv"),
        "figures": figures_block(),
    }
    facts["n_figures"] = len(facts["figures"])
    return jsonable(facts)


def serialise(facts: dict) -> str:
    # allow_nan=False is the whole point: strict JSON, parseable by jq/JS/Go/serde.
    return json.dumps(facts, indent=2, allow_nan=False, ensure_ascii=False) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--out", type=Path, default=OUTPUT,
                    help="output path (default: report/report_facts.json)")
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if the committed file is not current")
    args = ap.parse_args()

    text = serialise(build())

    if args.check:
        if not args.out.exists():
            raise SystemExit(f"{args.out} does not exist")
        if args.out.read_text(encoding="utf-8") != text:
            raise SystemExit(
                f"{args.out} is STALE -- re-run: python3 {Path(__file__).name}"
            )
        print(f"{args.out} is current")
        return

    args.out.write_text(text, encoding="utf-8")
    facts = json.loads(text)
    print(f"wrote {args.out}")
    print(f"  figures:            {facts['n_figures']}")
    print(f"  detected universe:  {facts['counts']['detected_universe']}")
    print(f"  UP / DOWN:          {facts['counts']['regulated_UP']} / "
          f"{facts['counts']['regulated_DOWN']}")
    print(f"  limma min adj.p:    {facts['limma']['min_adj_p']:.4f} "
          f"(vanilla {facts['limma']['vanilla_min_adj_p']:.4f})")
    print(f"  corr limma vs pipe: {facts['limma']['corr_limma_vs_pipeline_log2FC']:.6f} "
          f"(n={facts['limma']['corr_n_paired']})")


if __name__ == "__main__":
    main()
