"""Tests for ``run_pipeline.py`` -- the one-command pipeline runner.

Nothing here executes a pipeline stage, makes a network call, or writes into
``proteomics_de/results/``. The runner's whole job is *deciding what to run and
whether the result is right*, and that decision layer is testable in isolation:
selection is pure, the stage table is data, and output verification takes a
path. The only place a subprocess could appear is stubbed out and asserted
against.

The load-bearing test is ``test_pass_empty_expected_is_not_a_failure`` and its
neighbours. Three of this pipeline's outputs are header-only **by design**
(DECISIONS_LOG D2 and D6), so a runner that treats "0 rows" as "broken" would
report a false failure on every single run and train everyone to ignore it.
"""

from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import pytest

# run_pipeline.py lives at the repo root, which is not importable by default:
# pytest's prepend import mode inserts `proteomics_de/` (the first ancestor
# without an __init__.py), and conftest's sys.path fixture runs only at session
# setup -- after this module is imported at collection time. Resolve it from
# __file__ so `.venv/bin/pytest` works from any cwd, per the repo convention.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import run_pipeline as rp  # noqa: E402


def _freeze_module():
    """Import tools/freeze.py, the single definition of the freeze gate."""
    tools = _REPO_ROOT / "tools"
    if str(tools) not in sys.path:
        sys.path.insert(0, str(tools))
    import freeze

    return freeze


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _run_cli(argv: list[str], capsys) -> tuple[int, str]:
    """Invoke ``main(argv)`` and return ``(exit_code, stdout)``."""
    code = rp.main(argv)
    return code, capsys.readouterr().out


@pytest.fixture()
def frozen() -> dict:
    return rp.load_frozen_counts()


@pytest.fixture()
def no_subprocess(monkeypatch):
    """Make ANY subprocess launch an immediate, loud test failure."""

    def _boom(*args, **kwargs):  # pragma: no cover - only runs when a test fails
        raise AssertionError(f"a subprocess was spawned: args={args!r} kwargs={kwargs!r}")

    monkeypatch.setattr(rp.subprocess, "run", _boom)
    monkeypatch.setattr(rp.subprocess, "Popen", _boom)
    monkeypatch.setattr(rp.subprocess, "call", _boom)
    monkeypatch.setattr(rp.subprocess, "check_output", _boom)
    return _boom


class _FakeCompleted:
    def __init__(self, returncode: int = 0):
        self.returncode = returncode


# --------------------------------------------------------------------------
# Stage table integrity
# --------------------------------------------------------------------------
def test_thirteen_stages_in_the_documented_order():
    """The order is the deliverable; pin it explicitly so a reorder is visible."""
    assert [s.id for s in rp.STAGES] == [
        "foldchange",
        "qc_validate",
        "viz_qc_plots",
        "viz_volcano",
        "viz_ma_plot",
        "viz_heatmap",
        "enrich_string_ppi",
        "enrich_network_figure",
        "enrich_ora",
        "enrich_upset",
        "enrich_gsea",
        "gated_pca_cluster",
        "report",
    ]


def test_stage_table_validates_at_import():
    """``validate_stage_table`` runs at import; assert it is genuinely clean."""
    rp.validate_stage_table()  # must not raise


def test_dependencies_form_a_dag_with_no_forward_references():
    """Every dependency must appear EARLIER in the list than its dependant.

    Execution is serial and in list order, so a forward reference would mean a
    stage reading an artifact that has not been written yet. This is checked
    independently of ``validate_stage_table`` rather than by calling it, so a
    bug in the validator cannot hide a bug in the table.
    """
    position = {s.id: i for i, s in enumerate(rp.STAGES)}
    ids = set(position)
    for stage in rp.STAGES:
        for dep in stage.depends_on:
            assert dep in ids, f"{stage.id} depends on unknown stage {dep!r}"
            assert position[dep] < position[stage.id], (
                f"{stage.id} (position {position[stage.id] + 1}) has a FORWARD "
                f"dependency on {dep} (position {position[dep] + 1})"
            )

    # Reachability closure: no cycles can exist given the above, but assert it
    # directly so the property is stated, not merely implied.
    def ancestors(sid: str, seen: frozenset[str] = frozenset()) -> set[str]:
        assert sid not in seen, f"dependency cycle through {sid}"
        acc: set[str] = set()
        for dep in rp.STAGE_BY_ID[sid].depends_on:
            acc.add(dep)
            acc |= ancestors(dep, seen | {sid})
        return acc

    for stage in rp.STAGES:
        assert stage.id not in ancestors(stage.id)


def test_documented_ordering_constraints_are_encoded():
    """The three hard ordering rules from the build notes, as assertions."""
    anc = {}
    for stage in rp.STAGES:
        seen = set()
        stack = list(stage.depends_on)
        while stack:
            d = stack.pop()
            if d in seen:
                continue
            seen.add(d)
            stack.extend(rp.STAGE_BY_ID[d].depends_on)
        anc[stage.id] = seen

    assert "enrich_string_ppi" in anc["enrich_network_figure"], "8 must follow 7"
    assert "enrich_ora" in anc["enrich_upset"], "10 must follow 9"
    assert rp.STAGES[-1].id == "report", "13 must be last"
    # Everything is downstream of stage 1.
    for stage in rp.STAGES[1:]:
        assert "foldchange" in anc[stage.id], f"{stage.id} must be downstream of foldchange"


def test_every_stage_script_exists_on_disk():
    for stage in rp.STAGES:
        assert (rp.ROOT / stage.script).is_file(), f"{stage.id}: missing {stage.script}"


def test_network_and_slow_flags_match_the_known_stages():
    assert {s.id for s in rp.STAGES if s.network} == {
        "enrich_string_ppi", "enrich_ora", "enrich_gsea",
    }
    assert "enrich_gsea" in {s.id for s in rp.STAGES if s.slow}


# --------------------------------------------------------------------------
# frozen_counts.json is the single source of expected counts
# --------------------------------------------------------------------------
def test_row_counts_come_from_frozen_counts_not_hardcoded_literals(frozen):
    """Every non-zero expected row count must be a KEY, never an inline int.

    Zero is the one allowed literal: ``frozen_counts.json`` has no entry for
    the ORA term counts, and a hardcoded 0 there is a deliberate D6 assertion.
    """
    for stage in rp.STAGES:
        for out in stage.outputs:
            if isinstance(out.expected_rows, int) and out.expected_rows != 0:
                pytest.fail(
                    f"{stage.id}/{out.path} hardcodes {out.expected_rows}; use a "
                    f"frozen_counts.json key instead"
                )
            if isinstance(out.expected_rows, str):
                assert out.expected_rows in frozen, (
                    f"{stage.id}/{out.path} references missing key "
                    f"{out.expected_rows!r}"
                )


def test_headline_counts_resolve_to_the_verified_values(frozen):
    """Spot-check the resolved contract against the known-good build."""
    resolved = {
        out.path: rp.resolve_expected_rows(out, frozen)
        for stage in rp.STAGES for out in stage.outputs
    }
    assert resolved["proteomics_de/results/foldchange_all.csv"] == 1948
    # 604, not 606: DECISIONS_LOG D11 quarantines the 2 junk MaxQuant
    # row-index-list accessions to results/qc/quarantine_accessions.csv.
    assert resolved["proteomics_de/results/single_condition_proteins.csv"] == 604
    assert resolved["proteomics_de/results/onoff_proteins.csv"] == 10
    assert resolved["proteomics_de/results/qc_limma.csv"] == 1938
    assert resolved["proteomics_de/results/ipa_input.csv"] == 715
    assert resolved["proteomics_de/results/enrichment/gsea_results.csv"] == 568
    assert resolved["proteomics_de/results/enrichment/string_node_metrics.csv"] == 694


def test_the_three_by_design_empty_outputs_are_declared_with_reasons(frozen):
    empty = {
        out.path: out
        for stage in rp.STAGES for out in stage.outputs
        if rp.resolve_expected_rows(out, frozen) == 0
    }
    assert set(empty) == {
        "proteomics_de/results/ipa_input_significant.csv",
        "proteomics_de/results/enrichment/ora_up.csv",
        "proteomics_de/results/enrichment/ora_down.csv",
    }
    for path, out in empty.items():
        assert out.empty_reason, f"{path} expects 0 rows but gives no reason"
        assert "DECISIONS_LOG" in out.empty_reason


def test_validate_stage_table_rejects_a_reasonless_empty_output():
    bad = (
        rp.Stage(
            id="x", description="d", script="run_pipeline.py",
            outputs=(rp.Output("results/whatever.csv", expected_rows=0),),
        ),
    )
    with pytest.raises(ValueError, match="empty_reason"):
        rp.validate_stage_table(bad)


def test_validate_stage_table_rejects_a_forward_reference():
    bad = (
        rp.Stage(id="a", description="d", script="run_pipeline.py", depends_on=("b",)),
        rp.Stage(id="b", description="d", script="run_pipeline.py"),
    )
    with pytest.raises(ValueError, match="FORWARD"):
        rp.validate_stage_table(bad)


# --------------------------------------------------------------------------
# Output verification: emptiness is not failure
# --------------------------------------------------------------------------
def test_pass_empty_expected_is_not_a_failure(tmp_path, frozen):
    """A header-only file that was CONTRACTED to be header-only is a success."""
    csv = tmp_path / "ora_up.csv"
    csv.write_text("term_name,source,p_value\n", encoding="utf-8")

    out = rp.Output(str(csv), expected_rows=0, empty_reason="D6: nothing survives the background")
    verdict, detail = rp.check_output(out, frozen)

    assert verdict == rp.PASS_EMPTY_EXPECTED
    assert verdict != rp.FAIL
    assert "EXPECTED" in detail and "D6" in detail
    # And it must not be the worst outcome in a fold.
    assert rp._SEVERITY[rp.PASS_EMPTY_EXPECTED] < rp._SEVERITY[rp.FAIL]


def test_emptiness_IS_a_failure_when_rows_were_expected(tmp_path, frozen):
    csv = tmp_path / "gsea_results.csv"
    csv.write_text("Term,NES,FDR\n", encoding="utf-8")

    verdict, detail = rp.check_output(rp.Output(str(csv), expected_rows=568), frozen)
    assert verdict == rp.FAIL
    assert "568" in detail and "0" in detail


def test_unexpected_rows_in_a_by_design_empty_file_is_a_failure(tmp_path, frozen):
    """0 -> non-zero is a scientific event, and must stop the run for a human."""
    csv = tmp_path / "ora_up.csv"
    csv.write_text("term_name\nsome term\n", encoding="utf-8")

    out = rp.Output(str(csv), expected_rows=0, empty_reason="D6")
    verdict, detail = rp.check_output(out, frozen)
    assert verdict == rp.FAIL
    assert "SCIENTIFIC" in detail


def test_row_count_mismatch_and_missing_file_are_failures(tmp_path, frozen):
    csv = tmp_path / "t.csv"
    csv.write_text("a\n1\n2\n", encoding="utf-8")
    assert rp.check_output(rp.Output(str(csv), expected_rows=2), frozen)[0] == rp.PASS
    assert rp.check_output(rp.Output(str(csv), expected_rows=99), frozen)[0] == rp.FAIL
    assert rp.check_output(rp.Output(str(tmp_path / "nope.csv")), frozen)[0] == rp.FAIL


def test_zero_byte_file_fails_even_when_zero_rows_are_expected(tmp_path, frozen):
    """A header-only file still has a HEADER. Zero bytes is always broken."""
    csv = tmp_path / "empty.csv"
    csv.write_bytes(b"")
    out = rp.Output(str(csv), expected_rows=0, empty_reason="D6")
    verdict, detail = rp.check_output(out, frozen)
    assert verdict == rp.FAIL
    assert "zero bytes" in detail


def test_unasserted_row_count_passes_with_whatever_it_finds(tmp_path, frozen):
    csv = tmp_path / "skip_log.csv"
    csv.write_text("analysis,reason\npca,gate\n", encoding="utf-8")
    verdict, detail = rp.check_output(rp.Output(str(csv)), frozen)
    assert verdict == rp.PASS
    assert "1 rows" in detail


def test_unparsable_json_output_fails(tmp_path, frozen):
    bad = tmp_path / "meta.json"
    bad.write_text("{not json", encoding="utf-8")
    verdict, detail = rp.check_output(rp.Output(str(bad), kind="json"), frozen)
    assert verdict == rp.FAIL
    assert "JSON" in detail


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------
def test_from_selects_the_correct_suffix():
    selected = [s.id for s in rp.select_stages(from_="viz_volcano")]
    assert selected[0] == "viz_volcano"
    assert selected == [s.id for s in rp.STAGES][3:]
    assert len(selected) == 10
    assert "foldchange" not in selected and selected[-1] == "report"


def test_from_the_first_stage_is_equivalent_to_all():
    assert [s.id for s in rp.select_stages(from_="foldchange")] == \
           [s.id for s in rp.select_stages(all_=True)]


def test_only_selects_exactly_those_stages():
    assert [s.id for s in rp.select_stages(only="viz_volcano,viz_ma_plot")] == \
        ["viz_volcano", "viz_ma_plot"]
    assert [s.id for s in rp.select_stages(only="report")] == ["report"]


def test_only_normalises_order_whitespace_and_duplicates():
    """`--only` never lets the CLI reorder the pipeline."""
    assert [s.id for s in rp.select_stages(only=" report , foldchange ,report")] == \
        ["foldchange", "report"]


def test_unknown_stage_ids_are_rejected():
    with pytest.raises(SystemExit, match="viz_volcano_typo"):
        rp.select_stages(from_="viz_volcano_typo")
    with pytest.raises(SystemExit, match="nope"):
        rp.select_stages(only="viz_volcano,nope")


def test_no_selection_returns_nothing():
    assert rp.select_stages() == []


# --------------------------------------------------------------------------
# --list
# --------------------------------------------------------------------------
def test_list_enumerates_all_thirteen_stages(capsys, no_subprocess):
    code, out = _run_cli(["--list"], capsys)
    assert code == 0
    for stage in rp.STAGES:
        assert stage.id in out, f"--list omitted {stage.id}"
        assert stage.script in out
    assert f"{len(rp.STAGES)} stages" in out
    assert "13 stages" in out
    # Numbered 1..13.
    for n in range(1, len(rp.STAGES) + 1):
        assert f"{n:>2}. " in out


def test_list_marks_the_by_design_empty_outputs(capsys, no_subprocess):
    _, out = _run_cli(["--list"], capsys)
    # The bracketed form is the per-output marker; the legend at the foot of
    # --list quotes it instead, so this counts the outputs, not the prose.
    assert out.count("[0 rows EXPECTED") == 3
    assert "DECISIONS_LOG D2" in out
    assert "DECISIONS_LOG D6" in out


def test_list_marks_network_and_slow_stages(capsys, no_subprocess):
    _, out = _run_cli(["--list"], capsys)
    assert "[network]" in out
    assert "[network,slow]" in out


# --------------------------------------------------------------------------
# --dry-run executes nothing
# --------------------------------------------------------------------------
def test_dry_run_spawns_no_subprocess(capsys, no_subprocess):
    """The `no_subprocess` fixture turns any launch into an AssertionError."""
    code, out = _run_cli(["--all", "--dry-run"], capsys)
    assert code == 0
    assert "DRY RUN" in out
    assert "NOTHING WILL EXECUTE" in out
    assert "No subprocess was spawned" in out
    for stage in rp.STAGES:
        assert stage.id in out


def test_dry_run_writes_no_files(tmp_path, capsys, no_subprocess, monkeypatch):
    """Snapshot results/ mtimes+sizes around a dry run; nothing may move."""
    results = rp.PDE / "results"
    before = {p: (p.stat().st_mtime_ns, p.stat().st_size)
              for p in results.rglob("*") if p.is_file()}
    _run_cli(["--all", "--dry-run"], capsys)
    after = {p: (p.stat().st_mtime_ns, p.stat().st_size)
             for p in results.rglob("*") if p.is_file()}
    assert before == after


def test_dry_run_honours_from_and_skip_network(capsys, no_subprocess):
    code, out = _run_cli(["--from", "viz_volcano", "--dry-run", "--skip-network"], capsys)
    assert code == 0
    assert "10 stage(s) selected" in out
    assert "enrich_string_ppi  ->  SKIPPED" in out
    assert "enrich_ora  ->  SKIPPED" in out
    assert "enrich_gsea  ->  SKIPPED" in out
    # Non-network downstream stages are still planned.
    assert "would run: " in out
    assert "proteomics_de/enrich/network_figure.py" in out


def test_no_mode_flag_prints_help_and_exits_nonzero(capsys, no_subprocess):
    code, _ = _run_cli([], capsys)
    assert code == 2


# --------------------------------------------------------------------------
# Execution semantics (subprocess stubbed -- no stage ever really runs)
# --------------------------------------------------------------------------
def _fake_stage(stage_id: str, outputs=(), depends_on=(), network=False) -> rp.Stage:
    return rp.Stage(
        id=stage_id, description="fake", script="run_pipeline.py",
        outputs=tuple(outputs), depends_on=tuple(depends_on), network=network,
    )


def test_execute_returns_zero_when_a_stage_is_pass_empty_expected(
    tmp_path, frozen, monkeypatch
):
    """End-to-end proof: a header-only-by-design stage does NOT fail the run."""
    csv = tmp_path / "ora_up.csv"
    csv.write_text("term_name\n", encoding="utf-8")

    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return _FakeCompleted(0)

    monkeypatch.setattr(rp.subprocess, "run", fake_run)

    stage = _fake_stage(
        "fake_ora",
        outputs=[rp.Output(str(csv), expected_rows=0, empty_reason="D6: nothing survives")],
    )
    buf = io.StringIO()
    code = rp.execute([stage], skip_network=False, frozen=frozen, stream=buf)

    assert code == 0, "PASS_EMPTY_EXPECTED must not fail the run"
    assert len(calls) == 1
    text = buf.getvalue()
    assert "PASS_EMPTY_EXPECTED" in text
    assert "1 PASS_EMPTY_EXPECTED" in text
    assert "0 FAIL" in text
    assert "D6" in text


def test_execute_returns_nonzero_and_blocks_downstream_on_failure(frozen, monkeypatch):
    monkeypatch.setattr(rp.subprocess, "run", lambda cmd, **kw: _FakeCompleted(1))

    upstream = _fake_stage("up")
    downstream = _fake_stage("down", depends_on=["up"])
    buf = io.StringIO()
    code = rp.execute([upstream, downstream], skip_network=False, frozen=frozen, stream=buf)

    assert code == 1
    text = buf.getvalue()
    assert "FAIL" in text
    assert "BLOCKED" in text
    assert "upstream failure: up" in text


def test_skip_network_reports_skipped_not_failed_and_does_not_block(frozen, monkeypatch):
    """A deliberately skipped stage leaves committed artifacts in place."""
    calls = []

    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        return _FakeCompleted(0)

    monkeypatch.setattr(rp.subprocess, "run", fake_run)

    net = _fake_stage("net_stage", network=True)
    downstream = _fake_stage("after_net", depends_on=["net_stage"])
    buf = io.StringIO()
    code = rp.execute([net, downstream], skip_network=True, frozen=frozen, stream=buf)

    text = buf.getvalue()
    assert code == 0, "a skipped network stage must not fail the run"
    assert "SKIPPED  net_stage" in text
    assert "0 BLOCKED" in text, "a SKIPPED dependency must not block downstream"
    assert "0 FAIL" in text
    assert len(calls) == 1, "only the non-network downstream stage should run"


def test_nonzero_exit_code_skips_output_verification(tmp_path, frozen, monkeypatch):
    """A crashed stage is FAIL on its exit code; its stale outputs prove nothing."""
    monkeypatch.setattr(rp.subprocess, "run", lambda cmd, **kw: _FakeCompleted(3))
    csv = tmp_path / "good.csv"
    csv.write_text("a\n1\n", encoding="utf-8")

    stage = _fake_stage("boom", outputs=[rp.Output(str(csv), expected_rows=1)])
    result = rp.run_stage(stage, frozen, stream=io.StringIO())
    assert result.outcome == rp.FAIL
    assert "exit code 3" in result.detail
    assert "NOT verified" in result.detail


def test_stage_uses_sys_executable_and_repo_root_cwd(frozen, monkeypatch):
    """Never bare `python3`: only sys.executable has the scientific stack."""
    seen = {}

    def fake_run(cmd, **kwargs):
        seen["cmd"] = cmd
        seen["cwd"] = kwargs.get("cwd")
        return _FakeCompleted(0)

    monkeypatch.setattr(rp.subprocess, "run", fake_run)
    rp.run_stage(_fake_stage("s"), frozen, stream=io.StringIO())

    assert seen["cmd"][0] == sys.executable
    assert seen["cmd"][1] == str(rp.ROOT / "run_pipeline.py")
    assert seen["cwd"] == str(rp.ROOT)


# --------------------------------------------------------------------------
# --verify-frozen reuses tools/status.py
# --------------------------------------------------------------------------
def test_verify_frozen_delegates_to_tools_status(monkeypatch):
    """The runner must not reimplement the sha256 bookkeeping."""
    module = rp._load_status_tool()
    assert hasattr(module, "freeze_check")
    # Assert against tools/freeze.py's own default rather than a literal
    # filename, so renaming the manifest can't silently break the delegation.
    # The gate deliberately covers scientific OUTPUTS only -- freezing source
    # files made it fail on any refactor. See tools/freeze.py.
    assert module.PROTECTED_MANIFEST.name == "outputs.sha256"
    assert module.PROTECTED_MANIFEST.parent == rp.PDE / "tests" / "expected"
    entries = _freeze_module().read_manifest(module.PROTECTED_MANIFEST)
    assert entries, "outputs manifest is empty"
    assert not any(rel.endswith((".py", ".R")) for rel in entries), (
        "the freeze gate must not cover source files"
    )

    fake_rows = [("a.csv", "OK"), ("b.csv", "CHANGED"), ("c.csv", "MISSING")]
    fake_counts = {"OK": 1, "CHANGED": 1, "MISSING": 1}
    monkeypatch.setattr(rp, "_load_status_tool",
                        lambda: type("M", (), {"freeze_check": staticmethod(
                            lambda: (fake_rows, fake_counts))})())

    buf = io.StringIO()
    assert rp.verify_frozen(allow_drift=False, stream=buf) == 1
    text = buf.getvalue()
    assert "CHANGED  b.csv" in text
    assert "MISSING  c.csv" in text
    assert "1 OK · 1 CHANGED · 1 MISSING" in text

    buf = io.StringIO()
    assert rp.verify_frozen(allow_drift=True, stream=buf) == 0
    assert "--allow-drift" in buf.getvalue()


def test_verify_frozen_exits_zero_when_clean(monkeypatch):
    monkeypatch.setattr(rp, "_load_status_tool",
                        lambda: type("M", (), {"freeze_check": staticmethod(
                            lambda: ([("a.csv", "OK")], {"OK": 1, "CHANGED": 0, "MISSING": 0}))})())
    buf = io.StringIO()
    assert rp.verify_frozen(stream=buf) == 0
    assert "No drift" in buf.getvalue()


# --------------------------------------------------------------------------
# cwd independence
# --------------------------------------------------------------------------
def test_all_paths_resolve_from_file_not_cwd(tmp_path, monkeypatch, capsys, no_subprocess):
    """`cd /tmp && python /abs/run_pipeline.py --list` must behave identically."""
    monkeypatch.chdir(tmp_path)
    code, out = _run_cli(["--list"], capsys)
    assert code == 0
    assert str(rp.ROOT) in out
    assert rp.ROOT == Path(rp.__file__).resolve().parent
    assert rp.FROZEN_COUNTS_PATH.is_file()
    # And the frozen counts still load with an unrelated cwd.
    assert rp.load_frozen_counts()["foldchange_all_rows"] == 1948


def test_frozen_counts_file_is_valid_json_without_prose_keys(frozen):
    raw = json.loads(rp.FROZEN_COUNTS_PATH.read_text(encoding="utf-8"))
    assert any(k.startswith("_") for k in raw), "expected a _comment key in the raw file"
    assert not any(k.startswith("_") for k in frozen), "prose keys must be stripped"
