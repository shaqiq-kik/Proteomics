"""Acceptance checklist: has everything Dr. Walker asked for been done?

Run from the repo root:  python tools/walker_acceptance.py [--online]

Each requirement prints PASS, FAIL or PENDING (not produced yet). Exit code 0 only when every
line is PASS. --online also checks that every GitHub link in the v2 report returns HTTP 200.
Requirements come from his messages (2026-11-09 and 2026-11-11), see DECISIONS_LOG D20-D24.
"""
import csv, glob, json, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "proteomics_de" / "results"
ACC = "UniProt Accession Number"
ASSIGNED = {"Egf": "P01132", "Ereg": "Q61521", "Efna4": "O08542", "Tyro3": "P55144", "Mfap4": "Q9D1H9"}
OUT = []


def rec(rid, text, status, detail=""):
    OUT.append((rid, text, status, detail))


def rows(p):
    with open(p, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def run(*cmd):
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)


# ---- the list ---------------------------------------------------------------------------------
try:
    w = rows(R / "ipa_input_walker.csv")
    core = rows(R / "ipa_input.csv")
    by = {r["Gene names"]: r for r in w}
    rec("R1", "EGF and EREG are in the list as DOWN with the assigned 100-fold decrease (log2FC -6.643856), labelled assigned",
        "PASS" if all(by.get(g, {}).get("regulated") == "DOWN" and by[g]["log2FC"] == "-6.643856"
                      and by[g]["value_source"] == "assigned" and by[g]["p_value"] == "" for g in ("Egf", "Ereg")) else "FAIL")
    f = by.get("Frzb", {})
    rec("R2", "FRZB is in the list with its real measured value (UP, log2FC 4.40696, not the assigned value)",
        "PASS" if f.get("regulated") == "UP" and f.get("value_source") == "measured" and abs(float(f.get("log2FC", 0)) - 4.40696) < 1e-4 else "FAIL")
    sc = {r["gene"]: r for r in rows(R / "single_condition_proteins.csv")}
    ok = all(g in by and by[g]["regulated"] == "DOWN" and by[g]["value_source"] == "assigned"
             and sc.get(g, {}).get("detected_in") == "control_only" for g in ASSIGNED)
    rec("R3", "At least 2 more genes like EGF/EREG (control-only): EFNA4, TYRO3, MFAP4 added (and confirmed by him); all 5 assigned proteins are control-only in the data",
        "PASS" if ok and len(ASSIGNED) >= 5 else "FAIL")
    m = rows(R / "walker_additions_manifest.csv")
    first_ok = [(r[ACC], r["log2FC"], r["regulated"]) for r in w[:715]] == [(r[ACC], r["log2FC"], r["regulated"]) for r in core]
    rec("R4", "The 247 recovered proteins are dropped: list = 715 core + FRZB + 5 assigned = 721 rows (510 UP / 211 DOWN), manifest has 6 rows",
        "PASS" if len(w) == 721 and first_ok and sum(r["regulated"] == "UP" for r in w) == 510 and len(m) == 6
        and sum(r["tier"] == "partial" for r in w) == 1 else "FAIL", f"rows={len(w)} manifest={len(m)}")
except Exception as e:  # noqa: BLE001
    rec("R1-4", "walker list files readable", "FAIL", repr(e))

# ---- g:Profiler and STRING ----------------------------------------------------------------------
try:
    meta = json.loads((R / "enrichment_walker" / "ora_meta.json").read_text())
    vals = [meta.get("n_query_submitted_up"), meta.get("n_query_submitted_down")]
    nodes = {r["symbol"] for r in rows(R / "enrichment_walker" / "string_node_metrics.csv")}
    rec("R5", "g:Profiler and STRING were run on the 721 list (510 up / 211 down queried; all five assigned proteins are STRING nodes)",
        "PASS" if vals == [510, 211] and set(ASSIGNED) <= nodes else "FAIL", f"n_query_submitted up/down={vals}")
except Exception as e:  # noqa: BLE001
    rec("R5", "enrichment_walker outputs readable", "FAIL", repr(e))

# ---- IPA ----------------------------------------------------------------------------------------
tables = ["canonical_pathways", "upstream_regulators", "diseases_and_functions", "networks"]
IPA = R / "ipa_analysis_walker"
for tag, nm in (("walker721_with_assigned", "R6"), ("walker716_without_assigned", "R7")):
    have = [t for t in tables if (IPA / f"{tag}_{t}.txt").exists() and (IPA / f"{tag}_{t}.txt").stat().st_size > 1000]
    desc = ("IPA exports for the 721 list with the five assigned proteins" if nm == "R6"
            else "IPA exports for the 716 list without the five assigned proteins (the comparison run)")
    rec(nm, f"{desc} ({tag}_*.txt, 4 tables)", "PASS" if len(have) == 4 else "PENDING", f"{len(have)}/4 tables")
flag = IPA / "flagged_dependent_results.csv"
okf = False
if flag.exists():
    try:
        fr = rows(flag)
        okf = bool(fr) and "depends_on_assigned" in fr[0]
    except Exception:  # noqa: BLE001
        okf = False
rec("R8", "IPA results that depend on the five assigned proteins are flagged (ipa_analysis_walker/flagged_dependent_results.csv with a depends_on_assigned column)",
    "PASS" if okf else "PENDING")

# ---- replicate design correction ----------------------------------------------------------------
cfg = (ROOT / "proteomics_de" / "config" / "config.yaml").read_text()
t = run(sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "proteomics_de/tests/test_replicate_wording.py")
rec("R9", "Replicate design corrected everywhere: config says 2 biological replicates; no live file still carries the retired design claim (see test_replicate_wording.py)",
    "PASS" if re.search(r"n_biological_replicates_per_group:\s*2\b", cfg) and t.returncode == 0 else "FAIL",
    t.stdout.strip().splitlines()[-1] if t.stdout.strip() else t.stderr[-200:])

# ---- the v2 report ------------------------------------------------------------------------------
rep = sorted(glob.glob(str(ROOT / "report_v2*.html")))
if not rep:
    rec("R10", "A v2 report page for Dr. Walker exists (report_v2*.html) and reflects the 721 list, the corrected replicates and the flagged IPA run", "PENDING")
else:
    h = Path(rep[-1]).read_text(encoding="utf-8")
    # built from pieces so this file does not itself trip tests/test_replicate_wording.py
    retired = ("technical" + " replicates only", "one biological" + " replicate", "single biological" + " replicate", "n = 2 " + "technical")
    bad = [p for p in retired if p in h.lower()]
    good = "721" in h and re.search(r"2 biological replicates", h) is not None
    shas = set(re.findall(r"github\.com/shaqiq-kik/Proteomics/(?:blob|tree)/([0-9a-f]{40})/", h))
    pushed = all(run("git", "merge-base", "--is-ancestor", s, "origin/main").returncode == 0 for s in shas) and bool(shas)
    ok200 = True
    if "--online" in sys.argv:
        for u in sorted(set(re.findall(r'https://github\.com/shaqiq-kik/Proteomics/(?:blob|tree)/[0-9a-f]{40}/[^"\' <]+', h))):
            code = subprocess.run(["curl", "-s", "-o", "/dev/null", "-w", "%{http_code}", "-L", u], capture_output=True, text=True).stdout
            if code != "200":
                ok200 = False
                print("   link not 200:", code, u)
    rec("R10", f"v2 report ({Path(rep[-1]).name}) reflects 721, 2 biological replicates, no retired wording, links pinned to pushed commits" + (" and all links 200" if "--online" in sys.argv else ""),
        "PASS" if good and not bad and pushed and ok200 else "FAIL", f"bad_phrases={bad} pushed={pushed}")

# ---- v2 report content: every point of his 2026-11-11 reply is actually in the page ---------------
if rep:
    _fr = rows(IPA / "flagged_dependent_results.csv")
    _flip_cross = (sum(r["sign_flip"] == "yes" for r in _fr), sum(r["crosses_abs_z_2"] == "yes" for r in _fr))
    low = h.lower()
    text_only = re.sub(r"<[^>]+>", " ", re.sub(r"(?is)<(script|style).*?</\1>", " ", h))
    flat = re.sub(r"\s+", " ", text_only.replace("&nbsp;", " ")).lower()
    need = {
        "721 proteins (715 core + FRZB + 5 assigned)": "721" in flat and "715" in flat,
        "all five assigned genes named": all(g in flat for g in ("egf", "ereg", "efna4", "tyro3", "mfap4")),
        "assigned = stand-in, not measured": "stand-in" in flat and "not measurement" in flat or "not measured" in flat,
        "log2FC -6.643856": "-6.643856" in flat,
        "247 recovered proteins dropped": "247" in flat and "dropped" in flat,
        "2 separate experiments + low power": "separate experiments" in flat and "power" in flat,
        "0 proteins pass FDR 0.05 / exploratory": "fdr" in flat and "exploratory" in flat,
        "g:Profiler 0 significant terms": "g:profiler" in flat and "0 significant" in flat,
        "STRING 700 nodes / 6870 edges": "700" in flat and "6870" in flat or "6,870" in flat,
        "IPA run 2 with assigned vs 716 comparison": "run 2" in flat and "716" in flat,
        "flagged counts 81 / 411 / 739": all(n in flat for n in ("81", "411", "739")),
        "sign-flip and |z|=2-crossing counts (derived from flagged_dependent_results.csv)": all(str(n) in flat for n in _flip_cross),
        "flagged results table present": "flagged" in flat and "<table" in low,
        "caveat: default reference": "reference" in flat and "inflated" in flat,
        "caveat: 1000-row cap": "1000" in flat or "1,000" in flat,
        "caveat: 2277 vs 2431": ("2277" in flat or "2,277" in flat) and ("2431" in flat or "2,431" in flat),
        "pairing unconfirmed + question (31579/31578, 31581/31580)": "unconfirmed" in flat and all(s in flat for s in ("31579", "31578", "31581", "31580")),
        "links-go-live-after-push note": "links go live after the push" in flat,
        "no external scripts / stylesheets (self-contained)": not re.search(r'(?i)<script[^>]+src=|<link[^>]+href=["\']https?:', h),
    }
    miss = [k for k, v in need.items() if not v]
    # pinned GitHub links must point to files that exist at that (pushed) commit
    broken = []
    for s, pth in set(re.findall(r"github\.com/shaqiq-kik/Proteomics/(?:blob|tree)/([0-9a-f]{40})/([^\"' <#?]+)", h)):
        if run("git", "cat-file", "-e", f"{s}:{pth.rstrip('/')}").returncode != 0:
            broken.append(pth)
    rec("R13", "v2 report covers every point of his reply (721 / 247 dropped / 2 experiments / g:Profiler / STRING / IPA run 2 + flagged table / caveats / pairing question / go-live note) and pinned links resolve at their commit",
        "PASS" if not miss and not broken else "FAIL", f"missing={miss} broken_links={broken[:3]}")
else:
    rec("R13", "v2 report covers every point of his reply and pinned links resolve at their commit", "PENDING")

# ---- repo state ---------------------------------------------------------------------------------
dirty = [l for l in run("git", "status", "--porcelain").stdout.splitlines() if l[3:] not in (".DS_Store", "for_dr_walker.html")]
head = run("git", "rev-parse", "HEAD").stdout.strip()
remote = run("git", "ls-remote", "origin", "main").stdout.split()[:1]
rec("R11", "Everything is committed and pushed (no uncommitted changes besides .DS_Store / the old page; local HEAD equals origin/main)",
    "PASS" if not dirty and remote and remote[0] == head else "FAIL", f"dirty={len(dirty)} head={head[:8]}")
fr = run(sys.executable, "tools/freeze.py", "--check")
line = [l for l in fr.stdout.splitlines() if l.startswith("freeze:")]
rec("R12", "Freeze check is clean (no changed, missing or untracked frozen outputs)",
    "PASS" if line and "0 CHANGED, 0 MISSING, 0 UNTRACKED" in line[0] else "FAIL", line[0] if line else "")

# ---- report -------------------------------------------------------------------------------------
width = max(len(r[0]) for r in OUT)
for rid, text, st, det in OUT:
    print(f"[{st:7}] {rid:<{width}}  {text}" + (f"   ({det})" if det and st != "PASS" else ""))
n_ok = sum(1 for r in OUT if r[2] == "PASS")
print(f"\n{n_ok}/{len(OUT)} requirements PASS")
sys.exit(0 if n_ok == len(OUT) else 1)
