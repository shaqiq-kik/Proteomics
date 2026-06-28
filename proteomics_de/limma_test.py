"""Bug 7 — per-protein statistical significance via limma (R) + MinProb imputation.

The earlier bugs gave every both-condition protein a fold-change and a +/-0.585
"regulated" call, but a fold-change alone says nothing about whether the change is
*reproducible* — a protein measured twice can look "UP" by chance. Bug 7 adds the
missing half: a per-protein p-value and a Benjamini-Hochberg adjusted p-value
(FDR), computed with **limma** (a moderated t-test that borrows variance across
proteins) after **MinProb imputation** of missing intensities, following Peng et
al. 2024 (Nat Commun, DOI 10.1038/s41467-024-47899-w).

The statistics live in R (``limma_test.R``); this module owns the subprocess
boundary and the CSV file handoff. The split is deliberate:

  * Fork 1 (handoff)  — Python <-> R talk over CSV files, not stdin pipes, so the
                        exact matrix R saw is auditable on disk.
  * Fork 2 (eligibility) — only the ~1938 genuine both-condition proteins are
                        tested; the 10 ON_OFF proteins are excluded (testing them
                        would invent an entire absent condition) and the 606
                        single-condition proteins are not in foldchange_all.csv.
  * Fork 3 (transform) — log2 + imputation happen in R, close to limma, per Peng.

It only ADDS files (qc_limma.csv, ipa_input_significant.csv) and the auditable
intermediates; it overwrites NONE of the existing pipeline outputs. It reuses the
existing ``regulated`` column for the +/-0.585 part (it does not recompute the
threshold) and imports nothing from ``foldchange`` (the module graph stays acyclic;
it reads ``foldchange_all.csv`` from disk).

If the R worker fails for any reason, this module raises loudly and writes NO
result files — it never proceeds on empty or partial p-values.
"""

import os
import subprocess

import numpy as np
import pandas as pd

# --- Constants (single source of truth) -------------------------------------
ADJ_P_THRESHOLD = 0.05          # FDR cutoff for "significant"
R_SEED          = 42            # passed to R so MinProb imputation is reproducible
LIMMA_INPUT     = "_limma_input.csv"
LIMMA_OUTPUT    = "_limma_output.csv"
LIMMA_VERSIONS  = "_limma_versions.txt"

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
_R_SCRIPT = "limma_test.R"  # resolved relative to cwd=_HERE in the subprocess call

DEFAULT_RESULTS_DIR = os.path.join(_HERE, "results")
DEFAULT_FOLDCHANGE_CSV = os.path.join(DEFAULT_RESULTS_DIR, "foldchange_all.csv")

# Raw intensity columns -> fixed handoff order: control, control, treated, treated.
_CTRL_COLS = ["Intensity 31578", "Intensity 31580"]
_TRT_COLS = ["Intensity 31579", "Intensity 31581"]
_INTENSITY_COLS = _CTRL_COLS + _TRT_COLS
_HANDOFF_NAMES = ["ctrl_31578", "ctrl_31580", "trt_31579", "trt_31581"]
_HANDOFF_COLS = ["id", "gene"] + _HANDOFF_NAMES


def _missing_to_blank(series):
    """Raw intensity -> float, with blank/non-numeric/<=0 turned into NaN.

    NaN is written by ``to_csv`` as an empty cell (its default ``na_rep``), which
    R's ``read.csv`` reads back as ``NA`` so MinProb fills it. This is the single
    missing-value rule for the handoff.
    """
    num = pd.to_numeric(series, errors="coerce")
    return num.where(num > 0)  # keep > 0; everything else (<=0, NaN) -> NaN


def run_limma_test(foldchange_csv=DEFAULT_FOLDCHANGE_CSV, outdir=DEFAULT_RESULTS_DIR,
                   ebayes_mode="vanilla", qc_filename="qc_limma.csv",
                   output_name=LIMMA_OUTPUT, reuse_input=False, write_ipa=True):
    """Add limma p-values / FDR to the both-condition proteins.

    Reads ``foldchange_csv`` (the both-condition group), tests the eligible rows
    (``onoff`` empty) in R via MinProb + limma, and writes two NEW files into
    ``outdir``: ``qc_limma.csv`` (full audit) and ``ipa_input_significant.csv``
    (the honest "regulated AND significant" IPA list). Returns a summary-counts
    dict. Raises ``RuntimeError`` if the R worker fails (no partial outputs).

    The default call reproduces the committed baseline exactly. The extra args
    exist for the trend/robust side experiment and never change the default path:
      * ``ebayes_mode``  — "vanilla" (default) or "trend_robust"; forwarded to R.
      * ``qc_filename``  — name of the qc CSV under ``outdir`` (so the experiment
                           writes ``qc_limma_trend.csv`` instead of ``qc_limma.csv``).
      * ``output_name``  — name of the R worker's output intermediate in
                           ``proteomics_de/`` (so the experiment does not clobber
                           the committed ``_limma_output.csv``).
      * ``reuse_input``  — if True, reuse the existing ``_limma_input.csv`` instead
                           of regenerating it (proves identical data going in).
      * ``write_ipa``    — if False, skip the IPA file entirely (side experiment).
    """
    if ebayes_mode not in ("vanilla", "trend_robust"):
        raise ValueError(
            f"ebayes_mode must be 'vanilla' or 'trend_robust', got: {ebayes_mode!r}"
        )

    # 1) Eligibility (Fork 2). Read everything as strings + keep empty cells empty
    #    so the raw intensities pass through untouched and no NaN literal sneaks in.
    fc = pd.read_csv(foldchange_csv, dtype=str, keep_default_na=False)
    onoff_nonempty = fc["onoff"].str.strip() != ""
    n_onoff = int(onoff_nonempty.sum())
    eligible = fc[~onoff_nonempty].reset_index(drop=True)
    eligible_count = len(eligible)

    assert eligible_count == len(fc) - n_onoff, (
        f"eligible_count {eligible_count} != both-group {len(fc)} - ON_OFF {n_onoff}"
    )
    print("[Bug 7] limma statistical testing (R + MinProb imputation)")
    print(f"  eligible_count : {eligible_count}  "
          f"(both-group {len(fc)} - ON_OFF {n_onoff})")

    input_path = os.path.join(_HERE, LIMMA_INPUT)
    output_path = os.path.join(_HERE, output_name)
    # Versions filename is derived from the output name (mirrors the R worker), so a
    # non-default output writes its own versions file rather than the committed one.
    versions_name = os.path.splitext(output_name)[0].replace(
        "_limma_output", "_limma_versions") + ".txt"
    versions_path = os.path.join(_HERE, versions_name)

    # 2) Write the handoff CSV (Fork 1) in the EXACT column order — unless reusing
    #    the existing intermediate (so the data going in is provably the same).
    if reuse_input:
        if not os.path.exists(input_path):
            raise RuntimeError(
                f"BUG7: reuse_input=True but {LIMMA_INPUT} is missing — run the "
                "baseline first to produce it."
            )
        print(f"  reusing existing {LIMMA_INPUT} (not regenerated)")
    else:
        inp = pd.DataFrame(
            {
                "id": eligible["UniProt Accession Number"].values,
                "gene": eligible["Gene names"].values,
            }
        )
        for handoff_name, raw_col in zip(_HANDOFF_NAMES, _INTENSITY_COLS):
            inp[handoff_name] = _missing_to_blank(eligible[raw_col])
        inp[_HANDOFF_COLS].to_csv(input_path, index=False, encoding="utf-8")

    # Remove any stale worker artifacts so a leftover from a previous run can never
    # be mistaken for a fresh success (critical for the fail-loud guarantee).
    for stale in (output_path, versions_path):
        if os.path.exists(stale):
            os.remove(stale)

    # 3) Call the R worker. cwd=_HERE so the script + relative file names resolve
    #    and the intermediates land in proteomics_de/.
    try:
        result = subprocess.run(
            ["Rscript", _R_SCRIPT, LIMMA_INPUT, output_name, str(R_SEED), ebayes_mode],
            capture_output=True, text=True, cwd=_HERE,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "BUG7: 'Rscript' not found on PATH — cannot run the limma worker. "
            "Install R (and the limma + imputeLCMD packages). No outputs written."
        ) from exc

    if result.stdout:
        print(result.stdout.rstrip())
    if result.returncode != 0 or not os.path.exists(output_path):
        raise RuntimeError(
            "BUG7: limma R worker failed — refusing to continue with no/partial "
            f"p-values. returncode={result.returncode}\n"
            f"--- R stderr ---\n{result.stderr}\n--- R stdout ---\n{result.stdout}"
        )

    # 4) Read results and validate before trusting them.
    res = pd.read_csv(output_path)
    for col in ("p_value", "adj_p_value", "limma_log2FC"):
        vals = pd.to_numeric(res[col], errors="coerce")
        assert vals.notna().all() and np.isfinite(vals).all(), (
            f"BUG7: non-finite values found in '{col}' of {output_name}"
        )
    assert len(res) == eligible_count, (
        f"BUG7: R returned {len(res)} rows, expected eligible_count {eligible_count}"
    )

    # 5) Merge back onto the eligible rows by id (ids are unique here, so 1:1).
    merged = eligible.merge(
        res[["id", "limma_log2FC", "p_value", "adj_p_value"]],
        left_on="UniProt Accession Number", right_on="id", how="left",
    )
    assert merged["adj_p_value"].notna().all(), (
        "BUG7: an eligible protein has no limma result after the merge."
    )
    merged["significant"] = merged["adj_p_value"] < ADJ_P_THRESHOLD
    is_regulated = merged["regulated"].isin(["UP", "DOWN"])
    merged["regulated_significant"] = is_regulated & merged["significant"]

    os.makedirs(outdir, exist_ok=True)

    # 6a) qc_limma.csv — full audit of every eligible protein.
    qc = pd.DataFrame(
        {
            "id": merged["UniProt Accession Number"],
            "gene": merged["Gene names"],
            "Intensity 31578": merged["Intensity 31578"],
            "Intensity 31580": merged["Intensity 31580"],
            "Intensity 31579": merged["Intensity 31579"],
            "Intensity 31581": merged["Intensity 31581"],
            "limma_log2FC": merged["limma_log2FC"].round(6),
            "p_value": merged["p_value"],
            "adj_p_value": merged["adj_p_value"],
            "significant": merged["significant"],
            "regulated": merged["regulated"],
        }
    )
    qc_path = os.path.join(outdir, qc_filename)
    qc.to_csv(qc_path, index=False, encoding="utf-8")

    # 6b) ipa_input_significant.csv — only "regulated AND significant" rows, in the
    #     existing ipa_input.csv column layout (UniProt / Gene / log2FC / regulated)
    #     plus adj_p_value, for drop-in use. log2FC here is the pipeline's existing
    #     fold-change (UP/DOWN rows are always complete, so it is finite). Skipped
    #     entirely when write_ipa is False (the side experiment never touches IPA).
    ipa_sig_path = None
    if write_ipa:
        ipa_sig_rows = merged[merged["regulated_significant"]]
        ipa_sig = pd.DataFrame(
            {
                "UniProt Accession Number": ipa_sig_rows["UniProt Accession Number"],
                "Gene names": ipa_sig_rows["Gene names"],
                "log2FC": ipa_sig_rows["log2FC"],
                "regulated": ipa_sig_rows["regulated"],
                "adj_p_value": ipa_sig_rows["adj_p_value"],
            }
        )
        ipa_sig_path = os.path.join(outdir, "ipa_input_significant.csv")
        ipa_sig.to_csv(ipa_sig_path, index=False, encoding="utf-8")

    # 7) Summary: of the existing regulated proteins, how many are also significant
    #    (trustworthy) vs not (probably noise)?
    n_significant = int(merged["significant"].sum())
    n_regulated = int(is_regulated.sum())
    n_reg_sig = int(merged["regulated_significant"].sum())
    n_reg_not_sig = n_regulated - n_reg_sig

    print(f"  eBayes mode                : {ebayes_mode}")
    print(f"  proteins tested            : {eligible_count}")
    print(f"  significant (adj p < {ADJ_P_THRESHOLD}) : {n_significant}")
    print(f"  of {n_regulated} regulated (+/-0.585) proteins:")
    print(f"    also significant (trustworthy)   : {n_reg_sig}")
    print(f"    not significant (probably noise) : {n_reg_not_sig}")
    written = os.path.relpath(qc_path, _ROOT)
    if ipa_sig_path is not None:
        written += f", {os.path.relpath(ipa_sig_path, _ROOT)}"
    print(f"  wrote: {written}")
    if os.path.exists(versions_path):
        print(f"  versions: {os.path.relpath(versions_path, _ROOT)}")

    return {
        "eligible_count": eligible_count,
        "n_tested": eligible_count,
        "ebayes_mode": ebayes_mode,
        "n_significant": n_significant,
        "n_regulated": n_regulated,
        "n_regulated_significant": n_reg_sig,
        "n_regulated_not_significant": n_reg_not_sig,
    }


if __name__ == "__main__":
    run_limma_test()
