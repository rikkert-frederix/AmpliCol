from __future__ import annotations

import importlib.util
from io import StringIO
from pathlib import Path


def _load_validation_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "validate_ufo_sm_fixture.py"
    )
    spec = importlib.util.spec_from_file_location(
        "pyamplicol_validate_ufo_sm_fixture",
        path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _option(command: list[str], name: str) -> str:
    return command[command.index(name) + 1]


def test_non_lc_fixture_generation_matches_uniform_matrix_chunking(tmp_path) -> None:
    validation = _load_validation_module()
    command = validation._generation_command(
        {
            "process": "g g > t t~ g g g",
            "color_accuracy": "full",
            "base_key": "gg_tt_jets",
            "n_final": 5,
        },
        cache_payloads={"lc": {}, "nlc": {}, "full": {}},
        process_output=tmp_path / "process",
        model_cache=tmp_path / "model-cache",
        n_cores=5,
    )

    assert _option(command, "--batch-size") == "64"
    assert _option(command, "--symbolica-output-chunk-size") == "128"
    assert _option(command, "--symbolica-output-chunk-strategy") == "uniform"
    assert _option(command, "--max-coupling-order") == "QED=0"
    assert "--monitor" in command


def test_lc_fixture_generation_retains_selected_flow_auto_chunking(tmp_path) -> None:
    validation = _load_validation_module()
    command = validation._generation_command(
        {
            "process": "d d~ > z g g g g",
            "color_accuracy": "lc",
            "base_key": "dd_z_jets",
            "n_final": 5,
        },
        cache_payloads={
            "lc": {
                "entries": {
                    "dd_z_jets": {
                        "5": {"amplicol": {"reference_color_order": [1, 3, 4, 5, 6, 2]}}
                    }
                }
            },
            "nlc": {},
            "full": {},
        },
        process_output=tmp_path / "process",
        model_cache=tmp_path / "model-cache",
        n_cores=5,
    )

    assert _option(command, "--batch-size") == "128"
    assert _option(command, "--symbolica-output-chunk-size") == "128"
    assert _option(command, "--symbolica-output-chunk-strategy") == "auto"
    assert _option(command, "--lc-sector-ids") == "0"
    assert _option(command, "--reference-color-order") == "1,3,4,5,6,2"
    assert "--monitor" in command


def test_generation_command_streams_and_retains_monitor_output(
    monkeypatch,
    tmp_path,
    capsys,
) -> None:
    validation = _load_validation_module()
    calls = []

    class FakeProcess:
        stdout = StringIO("phase one\nphase two\n")

        @staticmethod
        def wait() -> int:
            return 0

    def fake_popen(command, **kwargs):
        calls.append((command, kwargs))
        return FakeProcess()

    monkeypatch.setattr(validation.subprocess, "Popen", fake_popen)
    log_path = tmp_path / "logs" / "case.log"

    validation._run_generation_command(
        ["python", "-m", "pyamplicol", "generate-process"],
        cwd=tmp_path,
        env={"PYTHONPATH": "src"},
        log_path=log_path,
    )

    assert len(calls) == 1
    assert calls[0][1]["stderr"] is validation.subprocess.STDOUT
    assert log_path.read_text(encoding="utf-8") == (
        "$ python -m pyamplicol generate-process\nphase one\nphase two\n"
    )
    assert "phase one\n" in capsys.readouterr().out


def test_archive_process_output_preserves_incomplete_artifact(tmp_path) -> None:
    validation = _load_validation_module()
    output_root = tmp_path / "validation"
    process_output = output_root / "processes" / "full_case_n5"
    process_output.mkdir(parents=True)
    (process_output / "partial.evaluator.bin").write_bytes(b"partial")

    archive = validation._archive_process_output(output_root, process_output)

    assert not process_output.exists()
    assert archive.parent.parent == output_root / "archive"
    assert (archive / "partial.evaluator.bin").read_bytes() == b"partial"
