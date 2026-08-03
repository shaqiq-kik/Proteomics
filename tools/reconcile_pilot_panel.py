#!/usr/bin/env python3
"""Reconciles the Pilot Project's 42-gene panel against every protein in the
current pipeline's regulated-adjacent deliverables, answering "of the ~40
proteins flagged in the earlier pilot analysis, how many still show up?"

One-off diagnostic, not a pipeline stage (DECISIONS_LOG D17): its input lives
entirely outside proteomics_de/ and outside the frozen-count machinery, and
the question it answers is a specific reply to a professor, not a recurring
pipeline concern. Kept as a script rather than a one-off computation so it is
reproducible if asked again after a data refresh.

Matching is by UniProt accession first (case-insensitive), falling back to
gene symbol (case-insensitive, both sides split on ';' for a protein group
and '/' for the pilot's human-nomenclature duplicate-gene rows like
"CRISP1/CRISP3"). Accession is primary because gene-symbol matching alone
produces known false negatives: the pilot's "C4A/C4B" is human nomenclature
for two paralogous genes with no separate mouse equivalent -- the mouse
genome carries one gene, C4b, at the exact UniProt accession (P01029) the
pilot cites -- and three further genes differ only by a gene-symbol
convention the pilot's human-style symbol does not use (LYZ -> Lyz2,
NAXE -> Apoa1bp, CYRIB -> Fam49b). All four resolve correctly by accession
and would be undercounted by symbol matching alone.

Usage:
    python tools/reconcile_pilot_panel.py
    python tools/reconcile_pilot_panel.py --write results/qc/pilot_panel_reconciliation.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
PDE = REPO_ROOT / "proteomics_de"
DEFAULT_PILOT_CSV = (
    REPO_ROOT / "Pilot Project" / "Analysis" / "General Analysis"
    / "cleaned_proteomics_data_with_QC_flags.csv"
)

#: The union that counts as "still in the current pipeline's output".
CURRENT_FILES = (
    "regulated_up.csv",
    "regulated_down.csv",
    "regulated_up_partial.csv",
    "regulated_down_partial.csv",
    "qualitative_changes.csv",
)


def _accession_tokens(value) -> set:
    if pd.isna(value):
        return set()
    return {tok.strip().upper() for tok in str(value).split(";") if tok.strip()}


def _gene_tokens(value) -> set:
    if pd.isna(value):
        return set()
    out = set()
    for tok in str(value).split(";"):
        for sub in tok.split("/"):
            sub = sub.strip().upper()
            if sub:
                out.add(sub)
    return out


def load_current_universe(results_dir):
    """``(accessions, genes)`` across every file in CURRENT_FILES that exists."""
    accs, genes = set(), set()
    for name in CURRENT_FILES:
        path = results_dir / name
        if not path.exists():
            continue
        df = pd.read_csv(path, float_precision="round_trip")
        for v in df["UniProt Accession Number"]:
            accs |= _accession_tokens(v)
        for v in df["Gene names"]:
            genes |= _gene_tokens(v)
    return accs, genes


def _current_status(accession, foldchange_all):
    """Best-effort context for an unmatched pilot protein: is it detected in
    the current data at all, and if so why didn't it reach a regulated file?
    """
    row = foldchange_all[foldchange_all["UniProt Accession Number"] == accession]
    if row.empty:
        return "not present in foldchange_all.csv (may be single-condition/onoff)"
    r = row.iloc[0]
    return f"complete={r['complete']}, regulated={r['regulated']!r}, log2FC={r.get('log2FC')}"


def reconcile(pilot_csv=DEFAULT_PILOT_CSV, results_dir=None):
    results_dir = Path(results_dir) if results_dir else PDE / "results"
    pilot = pd.read_csv(pilot_csv, float_precision="round_trip")
    current_accs, current_genes = load_current_universe(results_dir)

    fc_path = results_dir / "foldchange_all.csv"
    foldchange_all = (
        pd.read_csv(fc_path, float_precision="round_trip") if fc_path.exists() else None
    )

    rows = []
    for _, r in pilot.iterrows():
        acc_hit = bool(_accession_tokens(r["UniProt_Accession"]) & current_accs)
        gene_hit = bool(_gene_tokens(r["Gene"]) & current_genes)
        found = acc_hit or gene_hit
        status = None
        if not found and foldchange_all is not None:
            status = _current_status(str(r["UniProt_Accession"]).strip(), foldchange_all)
        rows.append({
            "Gene": r["Gene"],
            "UniProt_Accession": r["UniProt_Accession"],
            "pilot_log2FC": r["log_2_fold_change"],
            "matched_by_accession": acc_hit,
            "matched_by_gene_symbol": gene_hit,
            "found_in_current_deliverables": found,
            "current_status_if_not_found": status,
        })
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--pilot-csv", default=str(DEFAULT_PILOT_CSV))
    ap.add_argument("--results-dir", default=None)
    ap.add_argument("--write", default=None,
                     help="also write the full comparison table to this CSV path")
    args = ap.parse_args(argv)

    table = reconcile(args.pilot_csv, args.results_dir)
    n_found = int(table["found_in_current_deliverables"].sum())
    print(f"{n_found}/{len(table)} pilot-panel proteins found across "
          f"{', '.join(CURRENT_FILES)}")

    missing = table[~table["found_in_current_deliverables"]]
    if len(missing):
        print("\nNot found:")
        for _, r in missing.iterrows():
            print(f"  {r['Gene']} ({r['UniProt_Accession']}): {r['current_status_if_not_found']}")

    if args.write:
        table.to_csv(args.write, index=False, encoding="utf-8")
        print(f"\nSaved {args.write}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
