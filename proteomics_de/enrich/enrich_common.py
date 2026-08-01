"""
Shared utilities for the enrichment layer (Section 6, items 16 & 18):
over-representation analysis (ORA via g:Profiler), rank-based GSEA (gseapy),
and UpSet plots.

Reuses the committed report style (proteomics_de/viz/style.py) for palette /
save_fig, but writes to a SEPARATE manifest (figures_manifest_enrich.json) so
the existing report manifest (figures_manifest.json) is never touched.

Organism: MOUSE (Mus musculus) throughout -- g:Profiler organism="mmusculus",
gseapy Enrichr libraries are mouse libraries, gene symbols are mouse-style
(single leading capital, e.g. "Lama1"), NOT human uppercase.
"""

import json
import sys
from pathlib import Path

import pandas as pd

# ----------------------------------------------------------------------------
# Import the shared report style module (proteomics_de/viz/style.py)
# ----------------------------------------------------------------------------
ENRICH_DIR_PATH = Path(__file__).resolve().parent
PROTEOMICS_DE_DIR = ENRICH_DIR_PATH.parent
REPO_ROOT = PROTEOMICS_DE_DIR.parent
VIZ_DIR = PROTEOMICS_DE_DIR / "viz"
sys.path.insert(0, str(VIZ_DIR))
# The repo root has to be importable for the `proteomics_de.etl.*` package form
# (`proteomics_de/` is a namespace package -- no __init__.py). These scripts are
# run as bare files (`python proteomics_de/enrich/ora.py`), so only the script's
# own directory lands on sys.path automatically; the repo root does not.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
import style  # noqa: E402  (palette, apply_style(), save_fig, CHROME, etc.)

from proteomics_de.etl.accessions import first_token  # noqa: E402

RESULTS_DIR = style.RESULTS_DIR
FIGURES_DIR = style.FIGURES_DIR
ENRICHMENT_DIR = RESULTS_DIR / "enrichment"
RAW_DIR = ENRICHMENT_DIR / "raw"
ENRICHMENT_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

# NEW manifest -- never touches the existing figures_manifest.json
ENRICH_MANIFEST_PATH = FIGURES_DIR / "figures_manifest_enrich.json"

# Dataset-specific row/class counts live in exactly one place. They are NOT
# inlined here: the D7 correction (control/treated were inverted) swapped
# n_up/n_down from 206/509 to 509/206, and inline literals in this file were
# precisely what made that a multi-file edit. See tests/expected/frozen_counts.json.
FROZEN_COUNTS_PATH = PROTEOMICS_DE_DIR / "tests" / "expected" / "frozen_counts.json"


def load_frozen_counts():
    """Expected row/class counts for the committed dataset (single source of truth)."""
    return json.loads(FROZEN_COUNTS_PATH.read_text(encoding="utf-8"))

ORGANISM_GPROFILER = "mmusculus"
CAVEAT_ENRICHMENT = (
    "n=2 technical replicates; 0/1938 proteins survive per-protein FDR<0.05 "
    "(limma). ORA/GSEA below are HYPOTHESIS-GENERATING: pathway-level enrichment "
    "p-values (g:Profiler g:SCS / GSEA FDR) are a SEPARATE statistical test from "
    "the non-significant per-protein FDR -- a pathway can show up here among "
    "nominally-regulated proteins even though no single protein passes "
    "protein-level FDR. This is not a claim that any protein is significant."
)


def record_enrich_manifest(entries):
    """Merge manifest entries into figures_manifest_enrich.json (NEW file,
    keyed by `file`). Mirrors style.record_manifest but targets the new path
    so the legacy figures_manifest.json is never written to."""
    existing = []
    if ENRICH_MANIFEST_PATH.exists():
        try:
            existing = json.loads(ENRICH_MANIFEST_PATH.read_text())
        except json.JSONDecodeError:
            existing = []
    by_file = {e["file"]: e for e in existing}
    for e in entries:
        by_file[e["file"]] = e
    merged = list(by_file.values())
    ENRICH_MANIFEST_PATH.write_text(
        json.dumps(merged, indent=2, default=style._json_default)
    )
    return merged


# ----------------------------------------------------------------------------
# Gene symbol / background construction
# ----------------------------------------------------------------------------
def _first_token_or_accession(gene_val, acc_val):
    """Resolve one row to a single query identifier.

    Returns (symbol, was_semicolon_split, used_accession_fallback).
    Multi-gene "Gene names" entries (razor-peptide shared groups) are
    collapsed to their first listed gene symbol, per task spec. Rows with a
    missing/blank gene name fall back to the (first token of the) UniProt
    accession number so every row still contributes exactly one identifier.

    The accession fallback goes through
    :func:`proteomics_de.etl.accessions.first_token`, the one documented
    accession policy, rather than re-implementing ``str(acc).split(";")[0]``
    inline. That helper is behaviour-identical to the inline form it replaces
    (including its NaN -> ``"nan"`` handling); tests/test_enrich_common.py
    asserts the equivalence directly rather than trusting the comment.
    """
    if isinstance(gene_val, str) and gene_val.strip():
        parts = gene_val.split(";")
        return parts[0].strip(), len(parts) > 1, False
    return first_token(acc_val), False, True


def build_symbol_list(df, gene_col, acc_col):
    """One resolved identifier per row of df; returns (symbols, n_split, n_fallback)."""
    symbols, n_split, n_fallback = [], 0, 0
    for _, row in df.iterrows():
        sym, split, fb = _first_token_or_accession(row[gene_col], row[acc_col])
        symbols.append(sym)
        n_split += int(split)
        n_fallback += int(fb)
    return symbols, n_split, n_fallback


def load_background_and_queries():
    """Build the ORA background (detected proteome) and the UP/DOWN query lists.

    Background = union of foldchange_all.csv (1948 rows) and
    single_condition_proteins.csv (606 rows) -> 1948+606 = 2554 rows (the two
    files are disjoint by UniProt accession, verified during development), one
    resolved gene-symbol identifier per row (falling back to accession for the
    23 rows across both files with a missing gene name). After deduplicating
    identical symbols the background has fewer unique strings than 2554 rows
    (both counts are reported).

    UP / DOWN queries = foldchange_all.csv rows with regulated == UP (509) /
    DOWN (206), same per-row resolution. Those two counts were 206 / 509 until
    DECISIONS_LOG D7 corrected the inverted control/treated assignment; the
    query SETS swap wholesale (what was queried as UP is now the DOWN set and
    vice versa). Nothing here hardcodes them -- every expected count is read
    from tests/expected/frozen_counts.json.
    """
    counts = load_frozen_counts()
    fc = pd.read_csv(RESULTS_DIR / "foldchange_all.csv")
    scp = pd.read_csv(RESULTS_DIR / "single_condition_proteins.csv")

    fc_symbols, fc_split, fc_fb = build_symbol_list(
        fc, "Gene names", "UniProt Accession Number"
    )
    scp_symbols, scp_split, scp_fb = build_symbol_list(scp, "gene", "accession")
    background = sorted(set(fc_symbols) | set(scp_symbols))

    up_df = fc[fc["regulated"] == "UP"]
    down_df = fc[fc["regulated"] == "DOWN"]
    up_symbols, up_split, up_fb = build_symbol_list(
        up_df, "Gene names", "UniProt Accession Number"
    )
    down_symbols, down_split, down_fb = build_symbol_list(
        down_df, "Gene names", "UniProt Accession Number"
    )

    n_up_expected = counts["n_up"]
    n_down_expected = counts["n_down"]
    background_expected = counts["background_union"]
    assert len(up_df) == n_up_expected, (
        f"expected {n_up_expected} UP rows (frozen_counts.json:n_up), got {len(up_df)}"
    )
    assert len(down_df) == n_down_expected, (
        f"expected {n_down_expected} DOWN rows (frozen_counts.json:n_down), "
        f"got {len(down_df)}"
    )
    assert len(fc_symbols) + len(scp_symbols) == background_expected, (
        f"expected {background_expected}-row background union "
        f"(frozen_counts.json:background_union), got "
        f"{len(fc_symbols) + len(scp_symbols)}"
    )

    list_meta = {
        "background_rows_foldchange_all": len(fc_symbols),
        "background_rows_single_condition": len(scp_symbols),
        "background_row_union_n": len(fc_symbols) + len(scp_symbols),
        "background_unique_symbols_n": len(background),
        "foldchange_all_semicolon_split_n": fc_split,
        "foldchange_all_accession_fallback_n": fc_fb,
        "single_condition_semicolon_split_n": scp_split,
        "single_condition_accession_fallback_n": scp_fb,
        "up_n": len(up_symbols),
        "up_semicolon_split_n": up_split,
        "up_accession_fallback_n": up_fb,
        "down_n": len(down_symbols),
        "down_semicolon_split_n": down_split,
        "down_accession_fallback_n": down_fb,
    }
    return background, up_symbols, down_symbols, list_meta
