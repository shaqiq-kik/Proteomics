"""
proteomics_de/qc/validate.py

Runnable validator for the config + data-validation reproducibility layer
(research1.md Section 6, items 1 & 3; DECISIONS_LOG.md D1).

Loads each of the five committed, frozen pipeline output CSVs under
proteomics_de/results/, runs the pandera schema for it from schema.py, and
writes:
  - proteomics_de/results/qc/qc_report.json  (machine-readable)
  - proteomics_de/results/qc/qc_report.md    (human-readable)

READ-ONLY over the pipeline outputs: this script only ever opens the CSVs for
reading (pandas.read_csv) and only ever writes into results/qc/. It never
modifies foldchange_all.csv, qc_limma.csv, ipa_input.csv,
single_condition_proteins.csv, onoff_proteins.csv, or any other existing file.

Usage:
    python proteomics_de/qc/validate.py
    (equivalently: python proteomics_de/qc/schema.py, which delegates here)

Exit code is 0 if every file passes its schema, 1 otherwise.
"""

from __future__ import annotations

import datetime
import json
import platform
import sys
from pathlib import Path

import pandas as pd
import pandera
import pandera.pandas as pa

# Make this runnable both as `python proteomics_de/qc/validate.py` (cwd
# anywhere) and as `python -m` from within qc/, without requiring
# proteomics_de/ to be turned into a package (it currently is not; the
# existing flat scripts under proteomics_de/ are frozen and out of scope).
sys.path.insert(0, str(Path(__file__).resolve().parent))
import schema as qc_schema  # noqa: E402

QC_DIR = Path(__file__).resolve().parent
RESULTS_DIR = QC_DIR.parent / "results"
QC_OUT_DIR = RESULTS_DIR / "qc"

# Validation order mirrors the task's file-contract list.
FILE_ORDER = [
    "foldchange_all.csv",
    "qc_limma.csv",
    "ipa_input.csv",
    "single_condition_proteins.csv",
    "onoff_proteins.csv",
]

MAX_FAILURE_CASES_REPORTED = 50


def _failure_cases_to_records(failure_cases: pd.DataFrame) -> list[dict]:
    records = failure_cases.head(MAX_FAILURE_CASES_REPORTED).copy()
    # JSON can't serialize NaN/NaT/numpy scalars directly via json.dumps in
    # all cases; go through pandas' own JSON-safe conversion.
    records = json.loads(records.to_json(orient="records", date_format="iso"))
    return records


def validate_file(filename: str) -> dict:
    path = RESULTS_DIR / filename
    entry: dict = {
        "file": filename,
        "path": str(path),
        "expected_rows": qc_schema.EXPECTED_ROWS[filename],
    }

    if not path.exists():
        entry.update(
            {
                "exists": False,
                "passed": False,
                "rows": None,
                "columns": None,
                "n_failures": None,
                "failure_cases": [],
                "error": f"file not found: {path}",
            }
        )
        return entry

    df = pd.read_csv(path)
    entry["exists"] = True
    entry["rows"] = int(len(df))
    entry["columns"] = list(df.columns)
    entry["row_count_matches_expected"] = int(len(df)) == entry["expected_rows"]

    file_schema = qc_schema.FILE_SCHEMAS[filename]
    try:
        file_schema.validate(df, lazy=True)
    except pa.errors.SchemaErrors as exc:
        failure_cases = exc.failure_cases
        entry.update(
            {
                "passed": False,
                "n_failures": int(len(failure_cases)),
                "failure_cases": _failure_cases_to_records(failure_cases),
            }
        )
        return entry

    entry.update({"passed": True, "n_failures": 0, "failure_cases": []})
    return entry


def cross_file_checks(results_by_file: dict[str, dict]) -> list[dict]:
    """Checks that span more than one file -- not expressible in a single
    pandera DataFrameSchema, but real contract checks on the data lineage
    documented in config/config.yaml (enrichment.gprofiler.background_n).
    """
    checks = []

    fc_rows = results_by_file.get("foldchange_all.csv", {}).get("rows")
    qc_rows = results_by_file.get("qc_limma.csv", {}).get("rows")
    onoff_rows = results_by_file.get("onoff_proteins.csv", {}).get("rows")
    scp_rows = results_by_file.get("single_condition_proteins.csv", {}).get("rows")

    if fc_rows is not None and qc_rows is not None and onoff_rows is not None:
        ok = (fc_rows - onoff_rows) == qc_rows
        checks.append(
            {
                "name": "foldchange_all.csv rows minus onoff_proteins.csv rows equals qc_limma.csv rows",
                "detail": f"{fc_rows} - {onoff_rows} == {qc_rows}",
                "passed": bool(ok),
            }
        )

    if fc_rows is not None and scp_rows is not None:
        total = fc_rows + scp_rows
        ok = total == 2554
        checks.append(
            {
                "name": "foldchange_all.csv rows plus single_condition_proteins.csv rows equals the "
                "detected-proteome background size used by enrichment (config.yaml "
                "enrichment.gprofiler.background_n)",
                "detail": f"{fc_rows} + {scp_rows} == {total} (expected 2554)",
                "passed": bool(ok),
            }
        )

    return checks


def build_markdown_report(report: dict) -> str:
    lines = []
    lines.append("# proteomics_de QC validation report")
    lines.append("")
    lines.append(f"Generated: {report['generated_at']}")
    lines.append(
        f"pandera {report['environment']['pandera_version']} / "
        f"pandas {report['environment']['pandas_version']} / "
        f"Python {report['environment']['python_version']}"
    )
    lines.append("")
    lines.append(
        f"**Overall: {'PASS' if report['overall_passed'] else 'FAIL'}** "
        f"({sum(1 for f in report['files'] if f['passed'])}/{len(report['files'])} files passed)"
    )
    lines.append("")
    lines.append("This report validates the frozen, committed pipeline outputs under")
    lines.append("`proteomics_de/results/` against the pandera schemas in")
    lines.append("`proteomics_de/qc/schema.py`. It is READ-ONLY: no pipeline output was")
    lines.append("modified to produce this report.")
    lines.append("")
    lines.append("## Per-file results")
    lines.append("")
    lines.append("| File | Rows | Expected | Row count OK | Schema | Failures |")
    lines.append("|---|---|---|---|---|---|")
    for f in report["files"]:
        status = "PASS" if f["passed"] else "FAIL"
        rows = f.get("rows")
        expected = f.get("expected_rows")
        row_ok = f.get("row_count_matches_expected")
        n_fail = f.get("n_failures")
        lines.append(
            f"| {f['file']} | {rows} | {expected} | {row_ok} | {status} | {n_fail} |"
        )
    lines.append("")

    for f in report["files"]:
        if f["passed"]:
            continue
        lines.append(f"### FAILURES: {f['file']}")
        lines.append("")
        lines.append(
            f"{f['n_failures']} failure case(s) "
            f"(showing up to {MAX_FAILURE_CASES_REPORTED}):"
        )
        lines.append("")
        lines.append("```")
        lines.append(json.dumps(f["failure_cases"], indent=2))
        lines.append("```")
        lines.append("")

    lines.append("## Cross-file consistency checks")
    lines.append("")
    lines.append("| Check | Detail | Passed |")
    lines.append("|---|---|---|")
    for c in report["cross_file_checks"]:
        lines.append(f"| {c['name']} | {c['detail']} | {c['passed']} |")
    lines.append("")

    lines.append("## Documented data quirks accommodated by these schemas")
    lines.append("")
    lines.append(
        "- **UniProt accession regex extended for O/P/Q-prefixed accessions.** "
        "The literal task-spec regex `^[A-NR-Z0-9]([A-Z0-9]{5}|[A-Z0-9]{9})(-\\d+)?$` "
        "excludes accessions starting with O, P, or Q, which are common and "
        "legitimate in this dataset (e.g. P19137, Q9JHU4, O08528). "
        "`schema.py`'s `UNIPROT_TOKEN_RE` adds the missing `[OPQ][A-Z0-9]{5}` branch."
    )
    lines.append(
        "- **Multi-protein-group accessions** (';'-joined UniProt IDs, e.g. "
        "`P05132;P68181`) occur in foldchange_all.csv (48 rows), qc_limma.csv "
        "(48), ipa_input.csv (21), single_condition_proteins.csv (67). Each "
        "semicolon-delimited token is validated independently."
    )
    lines.append(
        "- **NaN gene names** occur in foldchange_all.csv (8), qc_limma.csv (8), "
        "ipa_input.csv (2), single_condition_proteins.csv (15) -- unannotated "
        "or ambiguous UniProt entries. Gene columns are nullable."
    )
    lines.append(
        "- **2 numeric-junk 'accession' rows in single_condition_proteins.csv** "
        "(long ';'-joined lists of bare integers, paired with NaN gene and "
        "NaN in all 4 intensity columns) are an upstream raw-data artifact "
        "passed through unmodified by the frozen pipeline. The schema accepts "
        "them via a narrowly-scoped exception requiring ALL of: accession "
        "matches `^[0-9;]+$`, gene is NaN, and all 4 intensities are NaN -- "
        "see `schema.py::_single_condition_accession_ok`."
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    QC_OUT_DIR.mkdir(parents=True, exist_ok=True)

    files_report = [validate_file(name) for name in FILE_ORDER]
    results_by_file = {f["file"]: f for f in files_report}
    cross_checks = cross_file_checks(results_by_file)

    overall_passed = all(f["passed"] for f in files_report) and all(
        c["passed"] for c in cross_checks
    )

    report = {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "environment": {
            "pandera_version": pandera.__version__,
            "pandas_version": pd.__version__,
            "python_version": platform.python_version(),
        },
        "results_dir": str(RESULTS_DIR),
        "overall_passed": overall_passed,
        "files": files_report,
        "cross_file_checks": cross_checks,
    }

    json_path = QC_OUT_DIR / "qc_report.json"
    json_path.write_text(json.dumps(report, indent=2))

    md_path = QC_OUT_DIR / "qc_report.md"
    md_path.write_text(build_markdown_report(report))

    print("=" * 72)
    print("proteomics_de QC validation")
    print("=" * 72)
    for f in files_report:
        status = "PASS" if f["passed"] else "FAIL"
        print(
            f"  [{status}] {f['file']:<32} rows={f.get('rows')!s:<6} "
            f"expected={f.get('expected_rows')!s:<6} failures={f.get('n_failures')}"
        )
    print("-" * 72)
    for c in cross_checks:
        status = "PASS" if c["passed"] else "FAIL"
        print(f"  [{status}] {c['name']}: {c['detail']}")
    print("-" * 72)
    print(f"OVERALL: {'PASS' if overall_passed else 'FAIL'}")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")

    return 0 if overall_passed else 1


if __name__ == "__main__":
    sys.exit(main())
