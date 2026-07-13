from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path


def _load_z_runner_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "docs"
        / "run_z_performance_table.py"
    )
    spec = importlib.util.spec_from_file_location("run_z_performance_table_docs", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_z_table_module():
    path = Path(__file__).resolve().parents[2] / "docs" / "z_performance_table.py"
    spec = importlib.util.spec_from_file_location("z_performance_table_docs", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _args() -> argparse.Namespace:
    return argparse.Namespace(
        batch_size=64,
        all_flow_batch_size=64,
        all_flow_output_chunk_size=8192,
    )


def test_z_table_reproduction_uses_latest_completed_jit_o1_row() -> None:
    table = _load_z_table_module()
    rendered = table.render_table(
        {
            "model_profile": "external-sm",
            "model_label": "external UFO/JSON SM",
            "entries": {
                "3": {"modes": {"jit_o1": {"status": "ok"}}},
                "4": {"modes": {"jit_o1": {"status": "missing"}}},
            },
        }
    )

    assert r"Reproduce the \(n=3\)" in rendered
    assert "'d d~ > z g g'" in rendered
    assert "manual/n3/jit_o1" in rendered
    assert "--reference-color-order 2,4,5,1,3" in rendered
    assert "manual/n7/jit_o1" not in rendered


def test_external_z_table_compares_wall_times_to_matching_builtin_row() -> None:
    table = _load_z_table_module()
    external = {
        "model_profile": "external-sm",
        "model_label": "external UFO/JSON SM",
        "entries": {
            "3": {
                "modes": {
                    "jit_o1": {
                        "status": "ok",
                        "wall_us_per_point": 3.0,
                        "all_flow_status": "ok",
                        "all_flow_wall_us_per_point": 5.0,
                    }
                }
            }
        },
    }
    built_in = {
        "entries": {
            "3": {
                "modes": {
                    "jit_o1": {
                        "status": "ok",
                        "wall_us_per_point": 2.0,
                        "all_flow_status": "ok",
                        "all_flow_wall_us_per_point": 10.0,
                    }
                }
            }
        }
    }

    rendered = table.render_table(external, built_in_data=built_in)

    assert rendered.count(r"\textbf{vs built-in}") == 4
    assert r"\textcolor{speedorange}{\texttt{x1.50}}" in rendered
    assert r"\textcolor{speedgreen}{\texttt{x0.50}}" in rendered


def _write_lc_cache(
    path: Path,
    *,
    selected_batch_size: int,
    selected_output_dir: str = "",
    all_flow_output_dir: str = "",
) -> None:
    path.write_text(
        json.dumps(
            {
                "entries": {
                    "dd_z_jets": {
                        "9": {
                            "pyamplicol_jit": {
                                "status": "ok",
                                "selected_output_dir": selected_output_dir,
                                "all_flow_output_dir": all_flow_output_dir,
                                "selected_generation_s": 2.0,
                                "wall_us_per_point": 3.0,
                                "runtime_us_per_point": 2.5,
                                "all_flow_generation_s": 20.0,
                                "all_flow_wall_us_per_point": 30.0,
                                "all_flow_runtime_us_per_point": 25.0,
                                "matrix_settings": {
                                    "selected_runtime_batch_size": (
                                        selected_batch_size
                                    ),
                                    "all_flow_runtime_batch_size": 64,
                                    "all_flow_symbolica_output_chunk_size": 8192,
                                    "all_flow_symbolica_jit_optimization_level": 1,
                                },
                            }
                        }
                    }
                }
            }
        )
    )


def test_external_z_generation_uses_matched_minimal_coupling_policy(
    monkeypatch,
    tmp_path: Path,
) -> None:
    runner = _load_z_runner_module()
    output_dir = tmp_path / "process"
    monkeypatch.setattr(runner, "MODEL_PROFILE", "external-sm")
    monkeypatch.setattr(runner, "MODEL_SOURCE", "compiled-model.json")

    command = runner._generate_command(
        "d d~ > Z",
        output_dir,
        mode="jit_o1",
        n_final=1,
        n_cores=5,
        batch_size=64,
        output_chunk_size=128,
        all_flows=False,
    )

    assert command[command.index("--coupling-order-policy") + 1] == "minimal"
    assert command[command.index("--model-cache-dir") + 1] == str(
        (output_dir / ".model_cache").resolve()
    )
    assert "--skip-generic-plan" in command
    assert "--no-runtime-lc-sector-selector" in command


def test_jit_o1_seed_rejects_mismatched_selected_batch(
    monkeypatch,
    tmp_path: Path,
) -> None:
    runner = _load_z_runner_module()
    cache = tmp_path / "lc.json"
    _write_lc_cache(cache, selected_batch_size=128)
    records: list[tuple[object, ...]] = []
    monkeypatch.setattr(runner, "LC_MATRIX_DATA", cache)
    monkeypatch.setattr(
        runner,
        "_record_row",
        lambda *args, **kwargs: records.append((*args, kwargs)),
    )

    assert runner._seed_from_lc_cache(9, "jit_o1", args=_args()) is False
    assert records == []


def test_jit_o1_reuses_matching_lc_all_flow_record(
    monkeypatch,
    tmp_path: Path,
) -> None:
    runner = _load_z_runner_module()
    cache = tmp_path / "lc.json"
    _write_lc_cache(cache, selected_batch_size=128)
    monkeypatch.setattr(runner, "LC_MATRIX_DATA", cache)

    record = runner._jit_o1_all_flow_record_from_lc(9, args=_args())

    assert record is not None
    assert record["all_flow_status"] == "ok"
    assert record["all_flow_generation_s"] == 20.0
    assert record["all_flow_time_batch_size"] == 64
    assert record["all_flow_symbolica_output_chunk_size"] == 8192


def test_retime_recovers_matching_lc_all_flow_artifact(
    monkeypatch,
    tmp_path: Path,
) -> None:
    runner = _load_z_runner_module()
    cache = tmp_path / "lc.json"
    lc_artifact = tmp_path / "lc-all-flow"
    fallback = tmp_path / "z-outputs"
    _write_lc_cache(
        cache,
        selected_batch_size=128,
        all_flow_output_dir=str(lc_artifact),
    )
    monkeypatch.setattr(runner, "LC_MATRIX_DATA", cache)
    monkeypatch.setattr(runner, "OUTPUT_ROOT", fallback)

    recovered = runner._retime_all_flow_artifact_path(
        9,
        "jit_o1",
        existing={"all_flow_generation_s": 20.0},
        args=_args(),
    )
    rejected = runner._retime_all_flow_artifact_path(
        9,
        "jit_o1",
        existing={"all_flow_generation_s": 21.0},
        args=_args(),
    )

    assert recovered == lc_artifact
    assert rejected == fallback / "n9" / "jit_o1_all_flows"


def test_z_runner_preserves_artifact_and_logs(tmp_path: Path) -> None:
    runner = _load_z_runner_module()
    output_dir = tmp_path / "jit_o1"
    output_dir.mkdir()
    (output_dir / "process_manifest.json").write_text("{}\n")
    log = tmp_path / "jit_o1.generate.log"
    log.write_text("old\n")

    preserved = runner._preserve_generated_output(output_dir)

    assert not output_dir.exists()
    preserved_dirs = list(tmp_path.glob("jit_o1_before_*"))
    preserved_logs = list(tmp_path.glob("jit_o1.generate.log.before_*"))
    assert len(preserved_dirs) == 1
    assert preserved == preserved_dirs[0]
    assert len(preserved_logs) == 1
    assert (preserved_dirs[0] / "process_manifest.json").is_file()
    assert preserved_logs[0].read_text() == "old\n"


def test_z_generation_uses_median_of_three_and_preserves_every_artifact(
    monkeypatch,
    tmp_path: Path,
) -> None:
    runner = _load_z_runner_module()
    output_dir = tmp_path / "jit_o1"
    elapsed = iter((3.0, 1.0, 2.0))
    preserved = iter((tmp_path / "sample-1", tmp_path / "sample-2"))

    monkeypatch.setattr(
        runner,
        "_run_json_command",
        lambda *args, **kwargs: {
            "_command_elapsed_s": next(elapsed),
            "generation_s": 0.5,
            "jit_compile_s": 0.4,
        },
    )
    monkeypatch.setattr(
        runner,
        "_preserve_generated_output",
        lambda path: next(preserved),
    )

    payload = runner._run_generation_json_command(
        ["generate"],
        timeout=None,
        log_path=tmp_path / "generate.log",
        output_dir=output_dir,
        median_samples=3,
    )

    assert payload["_command_elapsed_s"] == 2.0
    assert payload["_generation_sample_count"] == 3
    assert payload["_generation_median_sample_index"] == 3
    assert payload["_artifact_sample_index"] == 3
    assert [sample["artifact_dir"] for sample in payload["_generation_samples"]] == [
        str(tmp_path / "sample-1"),
        str(tmp_path / "sample-2"),
        str(output_dir),
    ]


def test_force_regeneration_keeps_amplicol_seed_but_rebuilds_jit(
    monkeypatch,
) -> None:
    runner = _load_z_runner_module()
    seeded: list[str] = []
    regenerated: list[str] = []

    def fake_seed(n_final, mode, *, args):
        assert n_final == 1
        seeded.append(mode)
        return True

    monkeypatch.setattr(runner, "_seed_from_lc_cache", fake_seed)
    monkeypatch.setattr(
        runner,
        "_run_pyamplicol_mode",
        lambda n_final, mode, args: regenerated.append(mode),
    )
    monkeypatch.setattr(runner, "_render_table", lambda: None)

    assert (
        runner.main(
            [
                "--n",
                "1",
                "--modes",
                "amplicol",
                "jit_o1",
                "--force-pyamplicol-regeneration",
            ]
        )
        == 0
    )
    assert seeded == ["amplicol"]
    assert regenerated == ["jit_o1"]


def test_retime_existing_preserves_generation_times(
    monkeypatch,
    tmp_path: Path,
) -> None:
    runner = _load_z_runner_module()
    output_root = tmp_path / "outputs"
    selected = output_root / "n9" / "jit_o3"
    all_flow = output_root / "n9" / "jit_o3_all_flows"
    selected.mkdir(parents=True)
    all_flow.mkdir(parents=True)
    (selected / "process_manifest.json").write_text("{}\n")
    (all_flow / "process_manifest.json").write_text("{}\n")
    data = tmp_path / "z.json"
    data.write_text(
        json.dumps(
            {
                "entries": {
                    "9": {
                        "modes": {
                            "jit_o3": {
                                "generation_s": 12.0,
                                "generation_measurement": {"sample_count": 3},
                                "all_flow_generation_s": 34.0,
                                "all_flow_generation_measurement": {
                                    "sample_count": 3
                                },
                                "all_flow_status": "ok",
                            }
                        }
                    }
                }
            }
        )
    )
    timed = iter(
        (
            {
                "profile": {
                    "wall_us_per_point": 5.0,
                    "core_evaluator_us_per_point": 4.0,
                }
            },
            {
                "profile": {
                    "wall_us_per_point": 7.0,
                    "core_evaluator_us_per_point": 6.0,
                }
            },
        )
    )
    records: list[tuple[tuple[object, ...], dict[str, object]]] = []
    monkeypatch.setattr(runner, "OUTPUT_ROOT", output_root)
    monkeypatch.setattr(runner, "Z_DATA", data)
    monkeypatch.setattr(
        runner,
        "_time_existing_artifact",
        lambda *args, **kwargs: next(timed),
    )
    monkeypatch.setattr(
        runner,
        "_record_row",
        lambda *args, **kwargs: records.append((args, kwargs)),
    )
    args = argparse.Namespace(
        target_runtime=10.0,
        batch_size=64,
        all_flow_batch_size=64,
        all_flow_output_chunk_size=8192,
    )

    runner._retime_pyamplicol_mode(9, "jit_o3", args)

    assert len(records) == 1
    _, record = records[0]
    assert record["generation_s"] == 12.0
    assert record["generation_measurement"] == {"sample_count": 3}
    assert record["wall_us_per_point"] == 5.0
    assert record["runtime_us_per_point"] == 4.0
    assert record["all_flow_generation_s"] == 34.0
    assert record["all_flow_generation_measurement"] == {"sample_count": 3}
    assert record["all_flow_wall_us_per_point"] == 7.0
    assert record["all_flow_runtime_us_per_point"] == 6.0
