"""Flag the IPA results that depend on the five assigned stand-in proteins (Dr. Walker, 2026-11-11).

Run from the repo root:  python3 tools/make_flagged_results.py

Compares the 721-protein IPA run (with the five assigned proteins) against the 716-protein run
(same settings, without them) table by table, and writes
proteomics_de/results/ipa_analysis_walker/flagged_dependent_results.csv  (one row per result that
depends on them) and flagged_summary.csv (counts per table). Standard library only.

A result is flagged (depends_on_assigned = yes) if ANY of these is true:
  lists_assigned_protein : one of EGF, EREG, EFNA4, TYRO3, MFAP4 is in its molecule list
                           (or is the regulator itself)
  z_changed_or_new       : its z-score differs between the two runs, or it exists in only one run
Nothing here is a finding about biology: the five values are assigned stand-ins, not measurements.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
IPA = ROOT / "proteomics_de" / "results" / "ipa_analysis_walker"
ASSIGNED = {"EGF", "EREG", "EFNA4", "TYRO3", "MFAP4"}
EPS = 0.0005  # IPA prints z-scores to 3 decimals
# table -> (key columns, z column, molecule column)
TABLES = {
    "canonical_pathways": (["Ingenuity Canonical Pathways"], "z-score", "Molecules"),
    "upstream_regulators": (["Upstream Regulator"], "Activation z-score", "Target Molecules in Dataset"),
    "diseases_and_functions": (["Categories", "Diseases or Functions Annotation"], "Activation z-score", "Molecules"),
}


def read_ipa(path):
    """IPA exports start with a copyright line and a blank line, then a tab-delimited table."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    i = next(k for k, l in enumerate(lines) if "\t" in l)
    hdr = [h.strip() for h in lines[i].split("\t")]
    rows = []
    for l in lines[i + 1:]:
        if l.strip():
            cells = l.split("\t")
            cells += [""] * (len(hdr) - len(cells))
            rows.append({h: c.strip() for h, c in zip(hdr, cells)})
    return rows


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def sign(x):
    return 0 if x is None else (1 if x > 0 else -1 if x < 0 else 0)


def main():
    out, summary = [], []
    for table, (keys, zcol, mcol) in TABLES.items():
        a = {tuple(r[k] for k in keys): r for r in read_ipa(IPA / f"walker721_with_assigned_{table}.txt")}
        b = {tuple(r[k] for k in keys): r for r in read_ipa(IPA / f"walker716_without_assigned_{table}.txt")}
        n_list = n_chg = n_flip = n_cross = n_flag = 0
        for key in sorted(set(a) | set(b)):
            ra, rb = a.get(key), b.get(key)
            z721 = num(ra[zcol]) if ra else None
            z716 = num(rb[zcol]) if rb else None
            mols = {m.strip().upper() for m in (ra[mcol] if ra else "").split(",")}
            hit = sorted((mols | {key[-1].upper()}) & ASSIGNED) if table == "upstream_regulators" else sorted(mols & ASSIGNED)
            if table == "upstream_regulators":
                hit = sorted((mols & ASSIGNED) | ({key[-1].upper()} & ASSIGNED))
            only_one = (ra is None) or (rb is None)
            changed = only_one or (z721 is None) != (z716 is None) or (z721 is not None and abs(z721 - z716) > EPS)
            flip = z721 is not None and z716 is not None and sign(z721) != sign(z716) and sign(z721) != 0 and sign(z716) != 0
            cross = (abs(z721 or 0) >= 2) != (abs(z716 or 0) >= 2)
            dep = bool(hit) or changed
            n_list += bool(hit); n_chg += changed; n_flip += flip; n_cross += cross; n_flag += dep
            if dep:
                out.append({
                    "table": table, "name": " | ".join(key),
                    "z_with_assigned_721": "" if z721 is None else z721,
                    "z_without_assigned_716": "" if z716 is None else z716,
                    "delta_z": "" if (z721 is None or z716 is None) else round(z721 - z716, 3),
                    "lists_assigned_protein": ",".join(hit),
                    "z_changed_or_new": "yes" if changed else "no",
                    "crosses_abs_z_2": "yes" if cross else "no",
                    "sign_flip": "yes" if flip else "no",
                    "depends_on_assigned": "yes",
                })
        summary.append({"table": table, "rows_721": len(a), "rows_716": len(b), "lists_assigned_protein": n_list,
                        "z_changed_or_new": n_chg, "sign_flips": n_flip, "crossed_abs_z_2": n_cross,
                        "flagged_total": n_flag})
    cols = ["table", "name", "z_with_assigned_721", "z_without_assigned_716", "delta_z", "lists_assigned_protein",
            "z_changed_or_new", "crosses_abs_z_2", "sign_flip", "depends_on_assigned"]
    with open(IPA / "flagged_dependent_results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, cols); w.writeheader(); w.writerows(out)
    with open(IPA / "flagged_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, list(summary[0])); w.writeheader(); w.writerows(summary)
    for s in summary:
        print(s)
    print("flagged rows written:", len(out))


if __name__ == "__main__":
    main()
