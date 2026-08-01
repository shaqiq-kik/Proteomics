# ⚠️ These enrichment results predate the D7 correction — DO NOT USE

**Status: STALE. Regeneration requires network approval.**

Every file in this directory was computed **before** `DECISIONS_LOG.md` **D7**
corrected the control/treated assignment (2026-08-01). They describe the
inverted experiment.

## What D7 changed

The pipeline had `31578/31580 = control` and `31579/31581 = treated`. The lab's
own Pilot Project labels them the other way round — `Vehicle = 31579/31581`,
`Testosterone = 31578/31580` — and the data agrees (pipeline vs. pilot log2FC
correlated at **r = −0.82** before the fix, **+0.82** after).

The fix is a pure sign inversion: **206 UP / 509 DOWN → 509 UP / 206 DOWN**.

## Why that invalidates each file here

| File | Why it is wrong now |
|---|---|
| `ora_up.csv`, `raw/gprofiler_up.json` | Queried the 206 proteins now known to be **DOWN** |
| `ora_down.csv`, `raw/gprofiler_down.json` | Queried the 509 proteins now known to be **UP** |
| `ora_meta.json`, `ora_top_terms_detail.json` | Same swap |
| `gsea_results.csv`, `gsea_meta.json` | Ranked by `sign(logFC) × …`; the whole ranking **reverses** |
| `string_edges.tsv`, `string_node_metrics.csv`, `string_meta.json` | Seed SET is the same 715 proteins, so the topology is unaffected — but every node's UP/DOWN annotation is inverted |

Note the *headline finding* is expected to survive: 0 terms passed correction in
either direction against the detected-proteome background (D6), and swapping two
label sets does not create enrichment. But that must be **re-demonstrated**, not
assumed.

## To regenerate

Requires outbound calls to STRING, g:Profiler and Enrichr (gene/protein
identifiers only — never intensities), per D3:

```bash
.venv/bin/python run_pipeline.py --only string_ppi,network_fig,ora,upset,gsea
```
