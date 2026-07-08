from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import pyamplicol
import pyamplicol.__main__ as cli
from pyamplicol.__main__ import (
    _attach_native_probe_comparison,
    _attach_rusticol_probe_comparison_from_runtime,
    _generation_build_kwargs,
    _runtime_backend,
    _runtime_evaluator_kwargs,
    _z_gluon_family_passed,
    main,
    parse_args,
)
from pyamplicol.native import ExternalMomentum
from pyamplicol.processes import ProcessOptions
from pyamplicol.reference import AmplicolProbePoint, AmplicolWorkflowResult


def _assert_legacy_native_command_payload(
    capsys,
    *,
    command: str,
    process: str,
    runtime_backend: str | None = None,
) -> dict[str, object]:
    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert payload["command"] == command
    assert payload["process"] == process
    assert "legacy native-kernel command" in payload["error"]
    assert "generate-process PROCESS OUTPUT_DIR" in payload["error"]
    assert "time-process OUTPUT_DIR" in payload["error"]
    if runtime_backend is not None:
        assert payload["runtime_backend"] == runtime_backend
    return payload


def test_cli_processes_json_and_legacy_export(
    capsys,
    tmp_path: Path,
) -> None:
    export = tmp_path / "processes.txt"

    assert (
        main(
            [
                "processes",
                "d d~ > z g g",
                "--legacy-output",
                str(export),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["n_unique_processes"] == 1
    assert payload["n_groups"] == 1
    assert payload["n_records"] == 1
    assert export.exists()


def test_cli_help_hides_legacy_non_dag_commands(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        parse_args(["--help"])

    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert "generate-process" in help_text
    assert "time-process" in help_text
    assert "process-plan" in help_text
    assert "validate-z-gluon-family" not in help_text
    assert "benchmark-z-gluon-modes" not in help_text
    assert "profile-tensor-evaluator" not in help_text
    assert "profile-dag-evaluator" not in help_text


def test_cli_generate_process_help_hides_legacy_backend_options(capsys) -> None:
    with pytest.raises(SystemExit) as exc:
        parse_args(["generate-process", "--help"])

    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert "--compiled-dag-evaluator" not in help_text
    assert "--runtime-backend" not in help_text
    assert "numeric-tensor-network" not in help_text
    assert "--symbolica-evaluator-backend" in help_text


def test_cli_compare_amplicol_help_hides_legacy_me_test_and_native_backends(
    capsys,
) -> None:
    with pytest.raises(SystemExit) as exc:
        parse_args(["compare-amplicol", "--help"])

    assert exc.value.code == 0
    help_text = capsys.readouterr().out
    assert "--runtime-backend {rusticol}" in help_text
    assert "--me-test" not in help_text
    assert "--compiled-dag-evaluator" not in help_text
    assert "numeric-tensor-network" not in help_text


def test_cli_process_plan_supports_line_pairing_sector_strategy(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "process-plan",
                "--lc-sector-strategy",
                "line-pairing-representatives",
                "--coupling-order-policy",
                "minimal",
                "d d~ > u u~ s s~ c c~",
                str(tmp_path / "plan"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    status = payload["lowering_status"]
    assert status["current_color_sector_count"] == 14
    assert status["current_color_sectors"] == [
        0,
        1,
        2,
        3,
        5,
        6,
        7,
        8,
        9,
        11,
        14,
        15,
        21,
        23,
    ]
    assert status["closure_count"] > 0
    assert status["truncated"] is False


def test_cli_process_parser_accepts_uppercase_bosons(capsys) -> None:
    assert main(["processes", "d d~ > Z g", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["request"]["rest"] == ["z", "g"]
    assert payload["n_entries"] == 1


def test_cli_processes_json_reports_process_sets(capsys) -> None:
    assert main(["processes", "d d~ > Z g | u u~ > Z g", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["n_entries"] == 2
    assert payload["default_key"] == "d_dbar_to_z_g"
    assert [entry["key"] for entry in payload["entries"]] == [
        "d_dbar_to_z_g",
        "u_ubar_to_z_g",
    ]


def test_cli_process_plan_writes_single_generic_manifest(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "process-plan",
                "d d~ > Z g",
                str(tmp_path / "plan"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest_path = Path(payload["manifest"])

    assert payload["available"] is True
    assert payload["kind"] == "pyamplicol-generic-process-plan"
    assert payload["process"] == "d d~ > z g"
    assert payload["key"] == "d_dbar_to_z_g"
    assert payload["planning_status"]["color_ready"] is True
    assert payload["planning_status"]["color_sector_count"] == 1
    assert payload["lowering_status"]["closure_count"] > 0
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["process_ir"]["outgoing_particles"] == ["d~", "d", "z", "g"]


def test_cli_process_plan_defaults_to_safe_topology_representatives(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "process-plan",
                "d d~ > z g g",
                str(tmp_path / "safe"),
                "--json",
            ]
        )
        == 0
    )
    safe = json.loads(capsys.readouterr().out)

    assert safe["lowering_status"]["current_color_sector_count"] == 1
    assert safe["lowering_status"]["amplitude_root_count"] == 24

    assert (
        main(
            [
                "process-plan",
                "g g > g g",
                str(tmp_path / "unsafe"),
                "--json",
            ]
        )
        == 0
    )
    unsafe = json.loads(capsys.readouterr().out)

    assert unsafe["lowering_status"]["current_color_sector_count"] == 3
    assert unsafe["lowering_status"]["amplitude_root_count"] == 9


def test_cli_process_plan_writes_process_set_with_unsupported_diagnostics(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "process-plan",
                "--color-accuracy",
                "full",
                "--max-qed-order",
                "2",
                "--no-species-reachability-pruning",
                "d d~ > z g | d d~ > z z g",
                str(tmp_path / "plans"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    root_manifest_path = Path(payload["manifest"])

    assert payload["available"] is True
    assert payload["kind"] == "pyamplicol-generic-process-set-plan"
    assert payload["default_process_key"] == "d_dbar_to_z_g"
    assert payload["generic_generation"]["pruning"]["max_coupling_orders"] == {
        "QED": 2
    }
    assert payload["generic_generation"]["pruning"]["species_reachability_pruning"] is False
    assert payload["n_processes"] == 2
    assert payload["processes"][1]["lowering_status"][
        "unimplemented_vertex_kinds"
    ] == []
    assert payload["processes"][0]["planning_status"]["color_ready"] is True
    assert payload["processes"][0]["planning_status"]["idenso_required"] is False
    assert (
        payload["processes"][0]["planning_status"]["generic_evaluator_ready"]
        is True
    )
    assert root_manifest_path.exists()
    root_payload = json.loads(root_manifest_path.read_text(encoding="utf-8"))
    assert root_payload["generic_generation"] == payload["generic_generation"]
    assert (
        root_manifest_path.parent
        / "subprocesses"
        / "d_dbar_to_z_z_g"
        / "generic_process_manifest.json"
    ).exists()


def test_cli_rejects_compiled_dag_shortcut_on_legacy_native_command() -> None:
    with pytest.raises(SystemExit) as exc:
        parse_args(
            [
                "evaluate",
                "--compiled-dag-evaluator",
                "d d~ > z g",
            ]
        )

    assert exc.value.code == 2


def test_cli_evaluate_photon_gluon_is_legacy_native_command(capsys) -> None:
    assert (
        main(
            [
                "evaluate",
                "d d~ > a g",
                "--runtime-backend",
                "python",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > a g",
        runtime_backend="python",
    )


def test_cli_evaluate_w_gluon_is_legacy_native_command(capsys) -> None:
    assert (
        main(
            [
                "evaluate",
                "u d~ > w+ g",
                "--runtime-backend",
                "python",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="u d~ > w+ g",
        runtime_backend="python",
    )


def test_cli_evaluate_w_gluon_dag_backend_is_legacy_native_command(capsys) -> None:
    assert (
        main(
            [
                "evaluate",
                "u d~ > w+ g",
                "--runtime-backend",
                "dag",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="u d~ > w+ g",
        runtime_backend="dag",
    )


def test_cli_evaluate_neutral_dilepton_gluon_is_legacy_native_command(capsys) -> None:
    assert (
        main(
            [
                "evaluate",
                "d d~ > e+ e- g",
                "--runtime-backend",
                "python",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > e+ e- g",
        runtime_backend="python",
    )


def test_cli_evaluate_zero_gluon_neutral_dilepton_is_legacy_native_command(
    capsys,
) -> None:
    assert (
        main(
            [
                "evaluate",
                "d d~ > e+ e-",
                "--runtime-backend",
                "python",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > e+ e-",
        runtime_backend="python",
    )


def test_cli_evaluate_charged_leptonic_w_gluon_is_legacy_native_command(
    capsys,
) -> None:
    assert (
        main(
            [
                "evaluate",
                "u d~ > e+ ve g",
                "--runtime-backend",
                "python",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="u d~ > e+ ve g",
        runtime_backend="python",
    )


def test_cli_evaluate_zero_gluon_charged_leptonic_w_is_legacy_native_command(
    capsys,
) -> None:
    assert (
        main(
            [
                "evaluate",
                "u d~ > e+ ve",
                "--runtime-backend",
                "python",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="u d~ > e+ ve",
        runtime_backend="python",
    )


def test_cli_rusticol_availability_allows_charged_leptonic_w_with_gluon() -> None:
    assert cli._rusticol_artifact_unavailable_message("u d~ > e+ ve g") is None
    assert cli._rusticol_artifact_unavailable_message("u d~ > e+ ve") is None


def test_cli_rusticol_availability_allows_neutral_dilepton_with_gluon() -> None:
    assert cli._rusticol_artifact_unavailable_message("d d~ > e+ e- g") is None
    assert cli._rusticol_artifact_unavailable_message("d d~ > e+ e-") is None


def test_cli_rusticol_runtime_backend_is_available() -> None:
    args = parse_args(
        [
            "profile-dag-evaluator",
            "--runtime-backend",
            "rusticol",
            "d d~ > z g",
        ]
    )

    assert _runtime_backend(args) == "rusticol"


def test_cli_compare_amplicol_defaults_to_production_rusticol_backend() -> None:
    args = parse_args(["compare-amplicol", "d d~ > z g"])

    assert _runtime_backend(args) == "rusticol"


def test_cli_compare_generic_pruning_keeps_reference_sector_default() -> None:
    args = parse_args(
        [
            "compare-amplicol",
            "--max-quark-pairs",
            "3",
            "d d~ > z g",
        ]
    )

    kwargs = cli._compare_generic_dag_pruning_kwargs(args, process=args.process)

    assert kwargs["max_quark_pairs"] == 3
    assert "selected_color_sector_ids" not in kwargs


def test_cli_compare_generic_pruning_accepts_explicit_sector_ids() -> None:
    args = parse_args(
        [
            "compare-amplicol",
            "--lc-sector-ids",
            "0,24",
            "d d~ > u u~",
        ]
    )

    kwargs = cli._compare_generic_dag_pruning_kwargs(args, process=args.process)

    assert kwargs["selected_color_sector_ids"] == {0, 24}


def test_cli_compare_amplicol_rejects_legacy_native_backends() -> None:
    with pytest.raises(SystemExit) as exc:
        parse_args(
            [
                "compare-amplicol",
                "--runtime-backend",
                "python",
                "d d~ > z g",
            ]
        )

    assert exc.value.code == 2


def test_cli_profile_dag_defaults_to_fast_rusticol_cxx_o3() -> None:
    args = parse_args(["profile-dag-evaluator", "d d~ > z g"])
    kwargs = _runtime_evaluator_kwargs(args)

    assert _runtime_backend(args) == "rusticol"
    assert kwargs["symbolica_evaluator_backend"] == "compiled-complex"
    assert kwargs["symbolica_compiled_preset"] == "runtime-o3"
    assert kwargs["symbolica_n_cores"] == 10
    assert kwargs["symbolica_compiled_chunk_compile_workers"] == 10
    assert kwargs["batch_size"] == 64


def test_cli_profile_dag_generate_only_flag_is_available(tmp_path: Path) -> None:
    args = parse_args(
        [
            "profile-dag-evaluator",
            "--runtime-backend",
            "rusticol",
            "--generate-only",
            "--save-evaluator-dir",
            str(tmp_path / "process"),
            "d d~ > z g",
        ]
    )

    assert _runtime_backend(args) == "rusticol"
    assert args.generate_only is True
    assert args.save_evaluator_dir == tmp_path / "process"


def test_cli_generate_process_minimal_command_uses_fast_rusticol_defaults(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "process"
    args = parse_args(["generate-process", "d d~ > z g", str(output_dir)])
    kwargs = _runtime_evaluator_kwargs(args)

    assert args.process == "d d~ > z g"
    assert args.output_dir == output_dir
    assert _runtime_backend(args) == "rusticol"
    assert kwargs["symbolica_evaluator_backend"] == "compiled-complex"
    assert kwargs["symbolica_compiled_preset"] == "runtime-o3"
    assert kwargs["symbolica_n_cores"] == 10
    assert kwargs["symbolica_compiled_chunk_compile_workers"] == 10
    assert kwargs["batch_size"] == 64
    assert args.color_accuracy == "lc"
    assert args.append is False
    assert args.replace is False
    pruning = cli._generic_dag_pruning_kwargs(args)
    assert args.numerical_filter_current is None
    assert args.numerical_current_merging is None
    assert pruning["numerical_filter_current"] is True
    assert pruning["numerical_current_merging"] is True
    assert args.numerical_current_samples == 10
    assert args.numerical_current_seed == 12345


def test_cli_generate_process_disables_numerical_current_passes_for_nlc_default(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "process"
    args = parse_args(
        [
            "generate-process",
            "--color-accuracy",
            "nlc",
            "g g > g g",
            str(output_dir),
        ]
    )
    pruning = cli._generic_dag_pruning_kwargs(args)

    assert args.numerical_filter_current is None
    assert args.numerical_current_merging is None
    assert pruning["numerical_filter_current"] is False
    assert pruning["numerical_current_merging"] is False


def test_cli_generate_process_can_disable_numerical_current_passes(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "process"
    args = parse_args(
        [
            "generate-process",
            "--no-numerical-filter-current",
            "--no-numerical-current-merging",
            "--numerical-current-samples",
            "3",
            "d d~ > z g",
            str(output_dir),
        ]
    )

    assert args.numerical_filter_current is False
    assert args.numerical_current_merging is False
    assert args.numerical_current_samples == 3


def test_cli_generate_process_runtime_o3_resolves_to_plain_cpp_o3(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "process"
    args = parse_args(["generate-process", "d d~ > z g g g g g", str(output_dir)])
    kwargs = _runtime_evaluator_kwargs(args)
    settings = cli._symbolica_settings_from_runtime_kwargs(
        kwargs,
        process=args.process,
    )

    assert settings.compiled_inline_asm == "none"
    assert settings.compiled_optimization_level == 3
    assert settings.compiled_output_chunk_size == 64


def test_cli_generate_process_explicit_jit_backend_overrides_fast_default(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "process"
    args = parse_args(
        [
            "generate-process",
            "--symbolica-evaluator-backend",
            "jit",
            "--symbolica-jit-opt-level",
            "2",
            "d d~ > z g",
            str(output_dir),
        ]
    )
    kwargs = _generation_build_kwargs(args, _runtime_backend(args), output_dir)

    assert kwargs["symbolica_evaluator_backend"] == "jit"
    assert kwargs["symbolica_jit_optimization_level"] == 2
    assert kwargs["symbolica_compiled_output_dir"] == str(output_dir / "compiled")


def test_cli_generate_process_set_writes_root_manifest(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "set"
    launched: list[list[str]] = []
    launched_kwargs: list[dict[str, object]] = []

    class FakePopen:
        next_pid = 5000

        def __init__(self, cmd, **kwargs):
            self.cmd = list(cmd)
            launched.append(self.cmd)
            launched_kwargs.append(dict(kwargs))
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 0
            if len(self.cmd) > 5 and self.cmd[3] == "generate-process":
                Path(self.cmd[5]).mkdir(parents=True, exist_ok=True)

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.125}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "--max-qed-order",
                "1",
                "--max-quark-pairs",
                "1",
                "--lc-sector-ids",
                "0",
                "--no-species-reachability-pruning",
                "--symbolica-jit-opt-level",
                "3",
                "d d~ > z g | u u~ > z g",
                str(output_dir),
                "--n_cores",
                "2",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    assert payload["available"] is True
    assert [cmd[4] for cmd in launched] == ["d d~ > z g", "u u~ > z g"]
    assert all("--n_cores" in cmd and cmd[cmd.index("--n_cores") + 1] == "1" for cmd in launched)
    assert all(
        "--symbolica-jit-optimization-level" in cmd
        and cmd[cmd.index("--symbolica-jit-optimization-level") + 1] == "3"
        for cmd in launched
    )
    assert all(kwargs["start_new_session"] is True for kwargs in launched_kwargs)
    assert manifest["default_process_key"] == "d_dbar_to_z_g"
    assert [entry["key"] for entry in manifest["processes"]] == [
        "d_dbar_to_z_g",
        "u_ubar_to_z_g",
    ]
    generation = manifest["generic_generation"]
    assert generation["n_cores"] == 2
    assert generation["lc_sector_strategy"] == "topology-representatives"
    assert generation["pruning"]["max_quark_pairs"] == 1
    assert generation["pruning"]["max_coupling_orders"] == {"QED": 1}
    assert generation["pruning"]["selected_color_sector_ids"] == [0]
    assert generation["pruning"]["species_reachability_pruning"] is False
    assert [entry["generation_request"] for entry in manifest["processes"]] == [
        generation,
        generation,
    ]
    assert (output_dir / "check_standalone.py").exists()


def test_cli_generate_process_set_reuses_crossing_equivalent_artifact(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "set"
    launched: list[list[str]] = []

    class FakePopen:
        next_pid = 5200

        def __init__(self, cmd, **kwargs):
            del kwargs
            self.cmd = list(cmd)
            launched.append(self.cmd)
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 0
            if len(self.cmd) > 5 and self.cmd[3] == "generate-process":
                Path(self.cmd[5]).mkdir(parents=True, exist_ok=True)

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.125}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "d d~ > u u~ g | u u~ > d~ d g",
                str(output_dir),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    assert [cmd[4] for cmd in launched] == ["d d~ > u u~ g"]
    assert [entry["key"] for entry in manifest["processes"]] == [
        "d_dbar_to_u_ubar_g",
        "u_ubar_to_dbar_d_g",
    ]
    first, second = manifest["processes"]
    assert second["path"] == first["path"]
    assert second["crossing_alias_of"] == "d_dbar_to_u_ubar_g"
    assert second["generation_request"] == manifest["generic_generation"]
    assert second["input_crossing_map"] == [
        {"sign": -1.0, "source_index": 2, "target_index": 0},
        {"sign": -1.0, "source_index": 3, "target_index": 1},
        {"sign": -1.0, "source_index": 1, "target_index": 2},
        {"sign": -1.0, "source_index": 0, "target_index": 3},
        {"sign": 1.0, "source_index": 4, "target_index": 4},
    ]
    assert payload["crossing_aliases"][0]["key"] == "u_ubar_to_dbar_d_g"


def test_process_input_crossing_map_handles_self_conjugate_initial_crossing() -> None:
    assert cli._process_input_crossing_map(
        "d d~ > g g z",
        "g d~ > d~ g z",
        color_accuracy="lc",
        options=ProcessOptions(),
    ) == [
        {"sign": -1.0, "source_index": 2, "target_index": 0},
        {"sign": 1.0, "source_index": 1, "target_index": 1},
        {"sign": -1.0, "source_index": 0, "target_index": 2},
        {"sign": 1.0, "source_index": 3, "target_index": 3},
        {"sign": 1.0, "source_index": 4, "target_index": 4},
    ]


def test_cli_generate_process_set_append_uses_real_crossing_representative(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "set"
    output_dir.joinpath("subprocesses", "d_dbar_to_u_ubar_g").mkdir(
        parents=True,
        exist_ok=True,
    )
    output_dir.joinpath("process_set_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "kind": "pyamplicol-generic-dag-process-set",
                "default_process_key": "u_ubar_to_dbar_d_g",
                "processes": [
                    {
                        "key": "u_ubar_to_dbar_d_g",
                        "process": "u u~ > d~ d g",
                        "path": "subprocesses/d_dbar_to_u_ubar_g",
                        "crossing_alias_of": "d_dbar_to_u_ubar_g",
                    },
                    {
                        "key": "d_dbar_to_u_ubar_g",
                        "process": "d d~ > u u~ g",
                        "path": "subprocesses/d_dbar_to_u_ubar_g",
                        "runtime_available": True,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    class UnexpectedPopen:
        def __init__(self, *args, **kwargs):
            del args, kwargs
            raise AssertionError("crossing-equivalent append should not launch a child")

    monkeypatch.setattr(cli.subprocess, "Popen", UnexpectedPopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "d d~ > u~ u g | d d~ > u u~ g",
                str(output_dir),
                "--append",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    appended = next(
        entry
        for entry in manifest["processes"]
        if entry["key"] == "d_dbar_to_ubar_u_g"
    )
    assert appended["path"] == "subprocesses/d_dbar_to_u_ubar_g"
    assert appended["crossing_alias_of"] == "d_dbar_to_u_ubar_g"
    assert payload["crossing_aliases"][0]["crossing_alias_of"] == (
        "d_dbar_to_u_ubar_g"
    )


def test_cli_generate_process_set_respects_parallel_worker_limit_and_reports_ram(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "set"
    live_processes: list[object] = []
    max_live = 0
    ram_active_sizes: list[int] = []

    class FakePopen:
        next_pid = 5050

        def __init__(self, cmd, **kwargs):
            nonlocal max_live
            del kwargs
            self.cmd = list(cmd)
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = None
            self.poll_count = 0
            live_processes.append(self)
            max_live = max(max_live, len(live_processes))
            Path(self.cmd[5]).mkdir(parents=True, exist_ok=True)

        def poll(self):
            self.poll_count += 1
            if self.poll_count >= 2:
                self.returncode = 0
                if self in live_processes:
                    live_processes.remove(self)
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.01}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    def fake_rss(active):
        ram_active_sizes.append(len(active))
        return f"{len(active)} children"

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", fake_rss)

    assert (
        main(
            [
                "generate-process",
                "d d~ > z g | u u~ > z g | s s~ > z g",
                str(output_dir),
                "--n_cores",
                "2",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is True
    assert max_live <= 2
    assert max(ram_active_sizes) <= 2
    assert any(size > 0 for size in ram_active_sizes)


def test_cli_process_set_child_progress_events_update_parent_progress(
    monkeypatch,
    tmp_path: Path,
) -> None:
    events = cli.queue.SimpleQueue()
    updates: list[dict[str, object]] = []

    class Entry:
        key = "d_dbar_to_z_g"

    class Progress:
        def update(self, **kwargs):
            updates.append(kwargs)

    class Process:
        pid = 1234

    line = (
        cli._CHILD_PROGRESS_PREFIX
        + json.dumps({"stage": "jit initialize", "item": "amp out=8 p=64"})
    )
    payload = cli._parse_child_progress_event(line)
    assert payload == {"stage": "jit initialize", "item": "amp out=8 p=64"}
    events.put((1234, payload))
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "1.5GB")

    cli._drain_child_progress_events(
        events,
        {1234: (Process(), Entry(), tmp_path, 0.0)},
        Progress(),
    )

    assert updates == [
        {
            "stage": "child jit initialize",
            "item": "d_dbar_to_z_g:amp out=8 p=64",
            "ram": "1.5GB",
        }
    ]
    assert cli._child_generation_environment()[cli._CHILD_PROGRESS_ENV] == "1"


def test_cli_process_set_child_progress_events_are_throttled(
    monkeypatch,
    tmp_path: Path,
) -> None:
    events = cli.queue.SimpleQueue()
    updates: list[dict[str, object]] = []

    class Entry:
        key = "d_dbar_to_z_g"

    class Progress:
        def update(self, **kwargs):
            updates.append(kwargs)

    class Process:
        pid = 1234

    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "1.5GB")
    monkeypatch.setattr(cli.time, "perf_counter", lambda: 10.0)
    events.put((1234, {"stage": "jit wait", "item": "elapsed=1.0s"}))
    events.put((1234, {"stage": "jit wait", "item": "elapsed=1.1s"}))
    seen: dict[int, tuple[str, str, float]] = {}

    cli._drain_child_progress_events(
        events,
        {1234: (Process(), Entry(), tmp_path, 0.0)},
        Progress(),
        seen,
        min_interval_s=0.5,
    )

    assert updates == [
        {
            "stage": "child jit wait",
            "item": "d_dbar_to_z_g:elapsed=1.0s",
            "ram": "1.5GB",
        }
    ]

    monkeypatch.setattr(cli.time, "perf_counter", lambda: 10.6)
    events.put((1234, {"stage": "jit wait", "item": "elapsed=1.6s"}))
    cli._drain_child_progress_events(
        events,
        {1234: (Process(), Entry(), tmp_path, 0.0)},
        Progress(),
        seen,
        min_interval_s=0.5,
    )

    assert updates[-1] == {
        "stage": "child jit wait",
        "item": "d_dbar_to_z_g:elapsed=1.6s",
        "ram": "1.5GB",
    }


def test_cli_generate_process_expands_builtin_p_to_child_artifacts(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "pp-set"
    launched: list[list[str]] = []

    class FakePopen:
        next_pid = 5100

        def __init__(self, cmd, **kwargs):
            del kwargs
            self.cmd = list(cmd)
            launched.append(self.cmd)
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 0
            Path(self.cmd[5]).mkdir(parents=True, exist_ok=True)

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.01}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "p p > z g",
                str(output_dir),
                "--n_cores",
                "3",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    assert payload["available"] is True
    assert len(launched) == 5
    assert [cmd[4] for cmd in launched] == [
        "d d~ > g z",
        "u u~ > g z",
        "s s~ > g z",
        "c c~ > g z",
        "b b~ > g z",
    ]
    assert [entry["key"] for entry in manifest["processes"][:5]] == [
        "d_dbar_to_g_z",
        "u_ubar_to_g_z",
        "s_sbar_to_g_z",
        "c_cbar_to_g_z",
        "b_bbar_to_g_z",
    ]
    assert [entry["crossing_alias_of"] for entry in manifest["processes"][5:]] == [
        "d_dbar_to_g_z",
        "u_ubar_to_g_z",
        "s_sbar_to_g_z",
        "c_cbar_to_g_z",
        "b_bbar_to_g_z",
    ]
    assert manifest["default_process_key"] == "d_dbar_to_g_z"
    assert all("--n_cores" in cmd and cmd[cmd.index("--n_cores") + 1] == "1" for cmd in launched)


def test_cli_generate_zero_gluon_process_expands_builtin_p(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "pp-z-set"
    launched: list[list[str]] = []

    class FakePopen:
        next_pid = 5200

        def __init__(self, cmd, **kwargs):
            del kwargs
            self.cmd = list(cmd)
            launched.append(self.cmd)
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 0
            Path(self.cmd[5]).mkdir(parents=True, exist_ok=True)

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.01}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "p p > z",
                str(output_dir),
                "--n_cores",
                "2",
                "--json",
            ]
        )
        == 0
    )

    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    assert len(launched) == 5
    assert [cmd[4] for cmd in launched[:3]] == [
        "d d~ > z",
        "u u~ > z",
        "s s~ > z",
    ]
    assert [entry["crossing_alias_of"] for entry in manifest["processes"][5:8]] == [
        "d_dbar_to_z",
        "u_ubar_to_z",
        "s_sbar_to_z",
    ]
    assert manifest["default_process_key"] == "d_dbar_to_z"


def test_cli_generate_charged_current_process_expands_builtin_p(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "pp-w-set"
    launched: list[list[str]] = []

    class FakePopen:
        next_pid = 5300

        def __init__(self, cmd, **kwargs):
            del kwargs
            self.cmd = list(cmd)
            launched.append(self.cmd)
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 0
            Path(self.cmd[5]).mkdir(parents=True, exist_ok=True)

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.01}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "p p > w+ g",
                str(output_dir),
                "--n_cores",
                "2",
                "--json",
            ]
        )
        == 0
    )

    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    assert [cmd[4] for cmd in launched] == [
        "u d~ > g w+",
        "c s~ > g w+",
    ]
    assert [entry["key"] for entry in manifest["processes"]] == [
        "u_dbar_to_g_wplus",
        "c_sbar_to_g_wplus",
        "dbar_u_to_g_wplus",
        "sbar_c_to_g_wplus",
    ]
    assert [entry.get("crossing_alias_of") for entry in manifest["processes"]] == [
        None,
        None,
        "u_dbar_to_g_wplus",
        "c_sbar_to_g_wplus",
    ]
    assert manifest["default_process_key"] == "u_dbar_to_g_wplus"


def test_generate_process_child_command_forwards_process_options(tmp_path: Path) -> None:
    args = parse_args(
        [
            "generate-process",
            "--flavour-scheme",
            "4",
            "--include-3qqbar",
            "--include-cc",
            "--include-resonance",
            "--color-accuracy",
            "full",
            "--max-coupling-order",
            "NP=0",
            "--max-qcd-order",
            "3",
            "--max-qed-order",
            "1",
            "--coupling-order-policy",
            "minimal",
            "--max-lc-current-line-groups",
            "2",
            "--max-quark-pairs",
            "3",
            "--ignore-particles",
            "h,a",
            "--ignore-vertex-kinds",
            "16,20",
            "--no-color-order-mask-pruning",
            "--no-species-reachability-pruning",
            "--no-numerical-filter-current",
            "--no-numerical-current-merging",
            "--numerical-current-samples",
            "6",
            "--numerical-current-seed",
            "321",
            "--numerical-current-relative-tolerance",
            "1e-13",
            "--numerical-current-zero-tolerance",
            "1e-250",
            "--lc-sector-strategy",
            "topology-representatives",
            "--lc-sector-ids",
            "0,24",
            "--symbolica-split-vertex-current-stages",
            "d d~ > z g",
            str(tmp_path / "set"),
        ]
    )

    command = cli._generate_process_child_command(
        args,
        "u d~ > e+ ve g",
        tmp_path / "subprocess",
    )

    assert command[3:5] == ["generate-process", "u d~ > e+ ve g"]
    assert command[command.index("--flavour-scheme") + 1] == "4"
    assert "--include-3qqbar" in command
    assert "--include-cc" in command
    assert "--include-resonance" in command
    assert command[command.index("--color-accuracy") + 1] == "full"
    assert command[command.index("--max-coupling-order") + 1] == "NP=0"
    assert command[command.index("--max-qcd-order") + 1] == "3"
    assert command[command.index("--max-qed-order") + 1] == "1"
    assert command[command.index("--coupling-order-policy") + 1] == "minimal"
    assert command[command.index("--max-lc-current-line-groups") + 1] == "2"
    assert command[command.index("--max-quark-pairs") + 1] == "3"
    assert command[command.index("--ignore-particles") + 1] == "h,a"
    assert command[command.index("--ignore-vertex-kinds") + 1] == "16,20"
    assert "--no-color-order-mask-pruning" in command
    assert "--no-species-reachability-pruning" in command
    assert "--no-numerical-filter-current" in command
    assert "--no-numerical-current-merging" in command
    assert command[command.index("--numerical-current-samples") + 1] == "6"
    assert command[command.index("--numerical-current-seed") + 1] == "321"
    assert command[command.index("--numerical-current-relative-tolerance") + 1] == "1e-13"
    assert command[command.index("--numerical-current-zero-tolerance") + 1] == "1e-250"
    assert (
        command[command.index("--lc-sector-strategy") + 1]
        == "topology-representatives"
    )
    assert command[command.index("--lc-sector-ids") + 1] == "0,24"
    assert "--symbolica-split-vertex-current-stages" in command
    assert command[command.index("--n_cores") + 1] == "1"


def test_cli_generate_process_set_append_and_replace(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "set"
    calls: list[tuple[str, Path]] = []

    class FakePopen:
        next_pid = 8000

        def __init__(self, cmd, **kwargs):
            del kwargs
            self.cmd = list(cmd)
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 0
            if len(self.cmd) > 5 and self.cmd[3] == "generate-process":
                process = self.cmd[4]
                output = Path(self.cmd[5])
                calls.append((process, output))
                output.mkdir(parents=True, exist_ok=True)
                (output / "marker.txt").write_text(process, encoding="utf-8")

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": True, "generation_s": 0.125}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
            [
                "generate-process",
                "d d~ > z g | u u~ > z g",
                str(output_dir),
                "--json",
            ]
        )
        == 0
    )
    capsys.readouterr()
    assert calls == [
        ("d d~ > z g", output_dir / "subprocesses" / "d_dbar_to_z_g"),
        ("u u~ > z g", output_dir / "subprocesses" / "u_ubar_to_z_g"),
    ]

    assert (
        main(
            [
                "generate-process",
                "d d~ > z g | u u~ > z g | s s~ > z g",
                str(output_dir),
                "--append",
                "--json",
            ]
        )
        == 0
    )
    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((output_dir / "process_set_manifest.json").read_text())
    assert payload["skipped"] == ["d_dbar_to_z_g", "u_ubar_to_z_g"]
    assert calls[-1] == ("s s~ > z g", output_dir / "subprocesses" / "s_sbar_to_z_g")
    assert len(calls) == 3
    assert manifest["default_process_key"] == "d_dbar_to_z_g"
    assert [entry["key"] for entry in manifest["processes"]] == [
        "d_dbar_to_z_g",
        "u_ubar_to_z_g",
        "s_sbar_to_z_g",
    ]

    assert (
        main(
            [
                "generate-process",
                "d d~ > z g | u u~ > z g | s s~ > z g",
                str(output_dir),
                "--replace",
                "--json",
            ]
        )
        == 0
    )
    capsys.readouterr()
    replaced_manifest = json.loads(
        (output_dir / "process_set_manifest.json").read_text()
    )
    assert calls[-3:] == [
        ("d d~ > z g", output_dir / "subprocesses" / "d_dbar_to_z_g"),
        ("u u~ > z g", output_dir / "subprocesses" / "u_ubar_to_z_g"),
        ("s s~ > z g", output_dir / "subprocesses" / "s_sbar_to_z_g"),
    ]
    assert replaced_manifest["default_process_key"] == "d_dbar_to_z_g"


def test_process_set_standalone_checker_selects_subprocess(
    tmp_path: Path,
) -> None:
    root = tmp_path / "set"
    first = root / "subprocesses" / "d_dbar_to_z_g"
    second = root / "subprocesses" / "u_ubar_to_z_g"
    first.mkdir(parents=True)
    second.mkdir(parents=True)
    root.joinpath("process_set_manifest.json").write_text(
        json.dumps(
            {
                "kind": "pyamplicol-rusticol-process-set",
                "default_process_key": "d_dbar_to_z_g",
                "processes": [
                    {
                        "key": "d_dbar_to_z_g",
                        "process": "d d~ > z g",
                        "path": "subprocesses/d_dbar_to_z_g",
                    },
                    {
                        "key": "u_ubar_to_z_g",
                        "process": "u u~ > z g",
                        "path": "subprocesses/u_ubar_to_z_g",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )
    root.joinpath("check_standalone.py").write_text(
        cli._PROCESS_SET_STANDALONE_CHECK_SCRIPT,
        encoding="utf-8",
    )
    nested_script = (
        "import json, sys\n"
        "from pathlib import Path\n"
        "Path('selected.json').write_text(json.dumps(sys.argv))\n"
    )
    (first / "check_standalone.py").write_text(nested_script, encoding="utf-8")
    (second / "check_standalone.py").write_text(nested_script, encoding="utf-8")

    result = subprocess.run(
        [
            sys.executable,
            str(root / "check_standalone.py"),
            "--process",
            "u u~ > z g",
            "--precision",
            "32",
            "--profile",
            "--target-runtime",
            "0.5",
            "--rusticol-folder",
            str(tmp_path / "site-packages"),
        ],
        cwd=second,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    forwarded = json.loads((second / "selected.json").read_text(encoding="utf-8"))
    assert forwarded[0] == str(second / "check_standalone.py")
    assert "--precision" in forwarded
    assert forwarded[forwarded.index("--precision") + 1] == "32"
    assert "--profile" in forwarded
    assert "--target-runtime" in forwarded
    assert forwarded[forwarded.index("--target-runtime") + 1] == "0.5"
    assert "--rusticol-folder" in forwarded


def test_cli_generate_process_set_surfaces_child_json_error(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    class FakePopen:
        next_pid = 6000

        def __init__(self, cmd, **kwargs):
            del cmd, kwargs
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = 1

        def poll(self):
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return json.dumps({"available": False, "error": "missing lowering"}), ""

        def terminate(self):
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")

    assert (
        main(
                [
                    "generate-process",
                    "d d~ > z g | u u~ > z g",
                    str(tmp_path / "set"),
                    "--json",
                ]
        )
        == 1
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert payload["failures"][0]["error"] == "missing lowering"
    assert payload["failures"][0]["payload"]["error"] == "missing lowering"


def test_cli_generate_process_set_handles_keyboard_interrupt(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    instances: list[object] = []

    class FakePopen:
        next_pid = 7000
        interrupt_once = True

        def __init__(self, cmd, **kwargs):
            del cmd, kwargs
            self.pid = FakePopen.next_pid
            FakePopen.next_pid += 1
            self.returncode = None
            self.terminated = False
            self.killed = False
            instances.append(self)

        def poll(self):
            if FakePopen.interrupt_once:
                FakePopen.interrupt_once = False
                raise KeyboardInterrupt
            return self.returncode

        def communicate(self, timeout=None):
            del timeout
            return "", ""

        def terminate(self):
            self.terminated = True
            self.returncode = -15

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def kill(self):
            self.killed = True
            self.returncode = -9

    monkeypatch.setattr(cli.subprocess, "Popen", FakePopen)
    monkeypatch.setattr(cli, "_rss_text_for_active_processes", lambda active: "0B")
    monkeypatch.setattr(cli, "_process_tree_pids", lambda root_pids: set())

    assert (
        main(
            [
                "generate-process",
                "d d~ > z g | u u~ > z g",
                str(tmp_path / "set"),
                "--n_cores",
                "2",
                "--json",
            ]
        )
        == 130
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert "interrupted" in payload["error"]
    assert instances
    assert all(getattr(instance, "terminated") for instance in instances)


def test_process_set_child_cleanup_is_idempotent_and_clears_active(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    class VanishingProcess:
        pid = 9000
        returncode = None

        def poll(self):
            return self.returncode

        def terminate(self):
            raise ProcessLookupError

        def wait(self, timeout=None):
            del timeout
            self.returncode = -15
            return self.returncode

        def kill(self):
            raise ProcessLookupError

        def communicate(self, timeout=None):
            del timeout
            return "", ""

    active = {9000: (VanishingProcess(), object(), tmp_path, 0.0)}
    monkeypatch.setattr(cli, "_process_tree_pids", lambda root_pids: {9000})
    monkeypatch.setattr(cli.os, "killpg", lambda pid, sig: (_ for _ in ()).throw(ProcessLookupError))
    monkeypatch.setattr(cli.os, "kill", lambda pid, sig: (_ for _ in ()).throw(ProcessLookupError))

    cli._terminate_generation_children(active)
    cli._terminate_generation_children(active)

    assert active == {}


def test_process_set_rss_helper_reports_parent_and_active_trees(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen_root_pids: list[tuple[int, ...]] = []

    def fake_tree_rss(root_pids: tuple[int, ...]) -> int:
        seen_root_pids.append(root_pids)
        return 42

    class FakeProcess:
        pass

    monkeypatch.setattr(cli, "_process_tree_rss_bytes", fake_tree_rss)
    assert cli._rss_text_for_active_processes({1234: (FakeProcess(), object(), Path("."), 0.0)}) == "42B"

    assert seen_root_pids
    assert seen_root_pids[0][0] == os.getpid()
    assert 1234 in seen_root_pids[0]


def test_cli_time_process_defaults_to_double_precision_and_ten_seconds(
    tmp_path: Path,
) -> None:
    process_dir = tmp_path / "process"
    args = parse_args(["time-process", str(process_dir)])

    assert args.process_dir == process_dir
    assert args.precision == 16
    assert args.target_runtime == 10.0
    assert args.process_key is None


def test_time_process_precision_upcasts_double_literals_with_trailing_zeros() -> None:
    assert cli._upcast_decimal_literal("1.23456789012345", 16) == "1.23456789012345"
    assert cli._upcast_decimal_literal("-1.25e+3", 8) == "-1.25e+3"
    assert cli._upcast_decimal_literal("1.25e+3", 18) == f"1.25{'0' * 15}e+3"
    assert cli._upcast_decimal_literal("-0.00125E-4", 18) == f"-0.00125{'0' * 15}E-4"
    assert cli._upcast_decimal_literal("500", 18) == f"500.{'0' * 15}"
    assert cli._upcast_decimal_literal("0", 18) == f"0.{'0' * 17}"


def test_cli_time_process_preserves_crossing_alias_selection(
    capsys,
    monkeypatch,
    tmp_path: Path,
) -> None:
    root = tmp_path / "set"
    representative = root / "subprocesses" / "d_dbar_to_u_ubar_g"
    representative.mkdir(parents=True)
    representative.joinpath("process_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "kind": "pyamplicol-generic-dag-process",
                "artifact_class": "generic-dag-schema-v2",
                "process": "d d~ > u u~ g",
                "compiled": {"runtime_available": True, "stage_evaluators": {}},
                "metadata": {"batch_size": 1},
            }
        ),
        encoding="utf-8",
    )
    representative.joinpath("validation_momenta.json").write_text(
        json.dumps(
            {
                "available": True,
                "points": [
                    [
                        {"momentum": ["1", "2", "3", "4"]},
                        {"momentum": ["5", "6", "7", "8"]},
                        {"momentum": ["9", "10", "11", "12"]},
                        {"momentum": ["13", "14", "15", "16"]},
                        {"momentum": ["17", "18", "19", "20"]},
                    ]
                ],
            }
        ),
        encoding="utf-8",
    )
    crossing_map = [
        {"sign": -1.0, "source_index": 2, "target_index": 0},
        {"sign": -1.0, "source_index": 3, "target_index": 1},
        {"sign": -1.0, "source_index": 1, "target_index": 2},
        {"sign": -1.0, "source_index": 0, "target_index": 3},
        {"sign": 1.0, "source_index": 4, "target_index": 4},
    ]
    root.joinpath("process_set_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 2,
                "kind": "pyamplicol-generic-dag-process-set",
                "default_process_key": "d_dbar_to_u_ubar_g",
                "processes": [
                    {
                        "key": "d_dbar_to_u_ubar_g",
                        "process": "d d~ > u u~ g",
                        "path": "subprocesses/d_dbar_to_u_ubar_g",
                    },
                    {
                        "key": "u_ubar_to_dbar_d_g",
                        "process": "u u~ > d~ d g",
                        "path": "subprocesses/d_dbar_to_u_ubar_g",
                        "crossing_alias_of": "d_dbar_to_u_ubar_g",
                        "input_crossing_map": crossing_map,
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    load_calls: list[tuple[str, str | None]] = []
    evaluated_points: list[object] = []

    class FakeRuntime:
        def metadata(self):
            return {
                "process": "u u~ > d~ d g",
                "schema_version": 2,
                "artifact_class": "generic-dag-schema-v2",
            }

        def evaluate(self, points):
            evaluated_points.append(points.tolist())
            return [7.0]

    class FakeRuntimeFactory:
        @staticmethod
        def load(process_dir, process_key=None):
            load_calls.append((process_dir, process_key))
            return FakeRuntime()

    class FakeRusticolModule:
        Runtime = FakeRuntimeFactory

        @staticmethod
        def build_profile():
            return "release"

        @staticmethod
        def build_target():
            return "test-target"

    monkeypatch.setitem(sys.modules, "rusticol", FakeRusticolModule())
    monkeypatch.setattr(
        cli,
        "_profile_rusticol_process",
        lambda runtime, points, **kwargs: {
            "wall_us_per_point": 1.0,
            "core_evaluator_us_per_point": 0.5,
            "block_size": 1,
            "samples": 1,
            "block_count": 1,
        },
    )
    args = parse_args(
        [
            "time-process",
            str(root),
            "--process",
            "u_ubar_to_dbar_d_g",
            "--target-runtime",
            "0",
            "--json",
        ]
    )
    args._suppress_rusticol_exit = True

    assert cli._cmd_time_process(args) == 0

    payload = json.loads(capsys.readouterr().out)
    assert load_calls == [(str(root), "u_ubar_to_dbar_d_g")]
    assert payload["process"] == "u u~ > d~ d g"
    assert evaluated_points == [
        [
            [
                [-13.0, -14.0, -15.0, -16.0],
                [-9.0, -10.0, -11.0, -12.0],
                [-1.0, -2.0, -3.0, -4.0],
                [-5.0, -6.0, -7.0, -8.0],
                [17.0, 18.0, 19.0, 20.0],
            ]
        ]
    ]


def test_cli_compacts_rusticol_profile_values_for_json_payload() -> None:
    compact = cli._compact_rusticol_profile(
        {
            "points": 4,
            "stage_evaluator_time_s": 1.0,
            "values": [1.0, 2.0, 3.0, 4.0],
        }
    )

    assert compact == {
        "points": 4,
        "stage_evaluator_time_s": 1.0,
        "value_count": 4,
    }


def test_cli_rejects_compiled_dag_flags_on_hidden_profile_command() -> None:
    with pytest.raises(SystemExit) as exc:
        parse_args(
            [
                "profile-dag-evaluator",
                "--compiled-dag-evaluator",
                "d d~ > z g",
            ]
        )

    assert exc.value.code == 2


def test_cli_generic_stage_jit_direct_translation_default() -> None:
    staged_args = parse_args(
        [
            "profile-dag-evaluator",
            "d d~ > z g",
        ]
    )
    disabled_args = parse_args(
        [
            "profile-dag-evaluator",
            "--symbolica-no-jit-direct-translation",
            "d d~ > z g",
        ]
    )

    assert _runtime_evaluator_kwargs(staged_args)[
        "symbolica_jit_direct_translation"
    ] is False
    assert _runtime_evaluator_kwargs(disabled_args)[
        "symbolica_jit_direct_translation"
    ] is False


def test_cli_generic_stage_uses_tuned_common_pair_defaults() -> None:
    staged_args = parse_args(
        [
            "profile-dag-evaluator",
            "d d~ > z g",
        ]
    )
    explicit_args = parse_args(
        [
            "profile-dag-evaluator",
            "--symbolica-cpe-iterations",
            "5",
            "--symbolica-max-common-pair-distance",
            "75",
            "d d~ > z g",
        ]
    )

    staged_kwargs = _runtime_evaluator_kwargs(staged_args)
    explicit_kwargs = _runtime_evaluator_kwargs(explicit_args)

    assert staged_kwargs["symbolica_cpe_iterations"] is None
    assert staged_kwargs["symbolica_max_common_pair_distance"] == 100
    assert explicit_kwargs["symbolica_cpe_iterations"] == 5
    assert explicit_kwargs["symbolica_max_common_pair_distance"] == 75


def test_cli_generate_writes_metadata_cache(capsys, tmp_path: Path) -> None:
    assert (
        main(
            [
                "generate",
                "d d~ > z g g",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 1
    )

    _assert_legacy_native_command_payload(
        capsys,
        command="generate",
        process="d d~ > z g g",
    )


def test_cli_generate_zero_gluon_writes_symbolica_scalar_artifact(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate",
                "d d~ > z",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 1
    )

    _assert_legacy_native_command_payload(
        capsys,
        command="generate",
        process="d d~ > z",
    )


def test_cli_generate_three_gluon_graph_records_auxiliary_tensor_route(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate",
                "d d~ > z g g g",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 1
    )

    _assert_legacy_native_command_payload(
        capsys,
        command="generate",
        process="d d~ > z g g g",
    )


def test_cli_profile_tensor_evaluator_is_legacy_command(capsys) -> None:
    assert (
        main(
            [
                "profile-tensor-evaluator",
                "d d~ > z g",
                "--repetitions",
                "2",
                "--evaluator-repetitions",
                "2",
                "--json",
            ]
        )
        == 1
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert payload["command"] == "profile-tensor-evaluator"
    assert "legacy Z+gluon-only command" in payload["error"]
    assert "generate-process" in payload["error"]


def test_cli_compare_amplicol_dry_run_defaults_to_library_momenta_probe(
    capsys,
) -> None:
    assert (
        main(
            [
                "compare-amplicol",
                "d d~ > z g g",
                "--points",
                "11",
                "--dry-run",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["mode"] == "amplicol_momenta_probe_library"
    assert payload["commands"][-1][:3] == [
        "./amplicol_generate",
        "--amplicol_momenta_probe=11",
        "--timing=11",
    ]
    assert payload["commands"][2][:2] == ["./amplicol_generate", "--library=create"]
    assert payload["commands"][3][:2] == ["make", "-j8"]
    assert "--library=use" in payload["commands"][-1]
    assert payload["pyamplicol_cache_dir"]


def test_cli_compare_amplicol_dry_run_keeps_explicit_legacy_me_test(
    capsys,
) -> None:
    assert (
        main(
            [
                "compare-amplicol",
                "d d~ > z g g",
                "--points",
                "11",
                "--me-test",
                "--dry-run",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["mode"] == "me_test"
    assert payload["commands"][-1][:3] == [
        "./amplicol_generate",
        "--me_test=11",
        "--timing=11",
    ]
    assert all("--library=create" not in command for command in payload["commands"])
    assert all("--library=use" not in command for command in payload["commands"])


def test_cli_compare_amplicol_dry_run_uses_library_backed_fixed_probe(
    capsys,
) -> None:
    assert (
        main(
            [
                "compare-amplicol",
                "d d~ > z",
                "--points",
                "5",
                "--amplicol-probe",
                "--dry-run",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["mode"] == "amplicol_fixed_probe_library"
    assert payload["commands"][-1][:3] == [
        "./amplicol_generate",
        "--amplicol_fixed_probe=5",
        "--timing=5",
    ]
    assert payload["commands"][2][:2] == ["./amplicol_generate", "--library=create"]
    assert payload["commands"][3][:2] == ["make", "-j8"]
    assert "--library=use" in payload["commands"][-1]


def test_cli_compare_amplicol_dry_run_uses_library_backed_probe_for_generic_process(
    capsys,
) -> None:
    assert (
        main(
            [
                "compare-amplicol",
                "d d~ > e+ e-",
                "--points",
                "5",
                "--amplicol-probe",
                "--dry-run",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["mode"] == "amplicol_probe_library"
    assert payload["commands"][:2] == [
        ["make", "cleanlib"],
        ["make", "-j8", "amplicol_generate"],
    ]
    assert payload["commands"][2][:2] == ["./amplicol_generate", "--library=create"]
    assert payload["commands"][3][:2] == ["make", "-j8"]
    assert payload["commands"][-1][:3] == [
        "./amplicol_generate",
        "--amplicol_probe=5",
        "--timing=5",
    ]
    assert "--library=use" in payload["commands"][-1]


def test_cli_validate_z_gluon_family_is_legacy_command(
    capsys,
) -> None:
    assert (
        main(
            [
                "validate-z-gluon-family",
                "--max-gluons",
                "2",
                "--points",
                "7",
                "--runtime-backend",
                "dag",
                "--dry-run",
                "--json",
            ]
        )
        == 1
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert payload["command"] == "validate-z-gluon-family"
    assert "legacy Z+gluon-only command" in payload["error"]
    assert "compare-amplicol --runtime-backend rusticol" in payload["error"]


def test_z_gluon_family_passed_enforces_complete_rows_and_tolerance() -> None:
    assert _z_gluon_family_passed(
        {
            "requested_rows": 7,
            "validated_rows": 7,
            "max_relative_difference": 2.0e-12,
        },
        rel_tol=1.0e-8,
    )
    assert not _z_gluon_family_passed(
        {
            "requested_rows": 7,
            "validated_rows": 6,
            "max_relative_difference": 2.0e-12,
        },
        rel_tol=1.0e-8,
    )
    assert not _z_gluon_family_passed(
        {
            "requested_rows": 7,
            "validated_rows": 7,
            "max_relative_difference": 2.0e-7,
        },
        rel_tol=1.0e-8,
    )


def test_native_probe_comparison_reports_runtime_metrics() -> None:
    mass = 91.188
    run = AmplicolWorkflowResult(
        commands=(),
        process_file=Path("processes.txt"),
        probe_points=(
            AmplicolProbePoint(
                point=1,
                group=1,
                integral=1,
                particles=(
                    ExternalMomentum(1, (mass / 2.0, 0.0, 0.0, mass / 2.0)),
                    ExternalMomentum(-1, (mass / 2.0, 0.0, 0.0, -mass / 2.0)),
                    ExternalMomentum(23, (mass, 0.0, 0.0, 0.0)),
                ),
                matrix_element=142.10653372872991,
            ),
        ),
    )
    payload: dict[str, object] = {}

    _attach_native_probe_comparison("d d~ > z", run, payload)

    points = payload["native_probe_points"]
    assert isinstance(points, list)
    assert points[0]["relative_difference"] == 0.0
    assert points[0]["native_runtime_s"] >= 0.0
    runtime = payload["native_runtime"]
    assert isinstance(runtime, dict)
    assert runtime["evaluated_points"] == 1
    assert runtime["total_s"] >= 0.0
    assert runtime["mean_per_point_s"] >= 0.0


def test_rusticol_probe_comparison_uses_probe_momenta() -> None:
    class FakeRusticolRuntime:
        def __init__(self) -> None:
            self.seen_shape: tuple[int, ...] | None = None

        def evaluate(self, momenta):
            self.seen_shape = tuple(momenta.shape)
            return [12.0]

        def metadata(self):
            return {
                "process": "d d~ > z",
                "schema_version": 2,
                "runtime_available": True,
            }

    run = AmplicolWorkflowResult(
        commands=(),
        process_file=Path("processes.txt"),
        probe_points=(
            AmplicolProbePoint(
                point=3,
                group=1,
                integral=1,
                particles=(
                    ExternalMomentum(1, (50.0, 0.0, 0.0, 50.0)),
                    ExternalMomentum(-1, (50.0, 0.0, 0.0, -50.0)),
                    ExternalMomentum(23, (100.0, 0.0, 0.0, 0.0)),
                ),
                matrix_element=10.0,
            ),
        ),
    )
    payload: dict[str, object] = {}
    runtime = FakeRusticolRuntime()

    _attach_rusticol_probe_comparison_from_runtime(run, payload, runtime)

    assert runtime.seen_shape == (1, 3, 4)
    points = payload["native_probe_points"]
    assert isinstance(points, list)
    assert points[0]["point"] == 3
    assert points[0]["native_matrix_element"] == 12.0
    assert points[0]["reference_matrix_element"] == 10.0
    assert points[0]["relative_difference"] == pytest.approx(2.0 / 12.0)
    backend = payload["native_runtime_backend"]
    assert isinstance(backend, dict)
    assert backend["backend"] == "rusticol-generic-schema-v2"
    summary = payload["native_probe_summary"]
    assert isinstance(summary, dict)
    assert summary["available_points"] == 1


def test_rusticol_probe_comparison_reorders_to_artifact_pdg_order() -> None:
    class FakeRusticolRuntime:
        def __init__(self) -> None:
            self.seen_momenta = None

        def evaluate(self, momenta):
            self.seen_momenta = momenta.copy()
            return [2.0]

        def metadata(self):
            return {
                "process": "d d~ > z g",
                "schema_version": 2,
                "runtime_available": True,
            }

    run = AmplicolWorkflowResult(
        commands=(),
        process_file=Path("processes.txt"),
        probe_points=(
            AmplicolProbePoint(
                point=1,
                group=1,
                integral=1,
                particles=(
                    ExternalMomentum(1, (50.0, 0.0, 0.0, 50.0)),
                    ExternalMomentum(-1, (50.0, 0.0, 0.0, -50.0)),
                    ExternalMomentum(21, (30.0, 1.0, 2.0, 3.0)),
                    ExternalMomentum(23, (70.0, -1.0, -2.0, -3.0)),
                ),
                matrix_element=2.0,
            ),
        ),
    )
    payload: dict[str, object] = {}
    runtime = FakeRusticolRuntime()

    _attach_rusticol_probe_comparison_from_runtime(
        run,
        payload,
        runtime,
        external_pdg_order=(1, -1, 23, 21),
    )

    assert runtime.seen_momenta is not None
    assert tuple(runtime.seen_momenta.shape) == (1, 4, 4)
    assert runtime.seen_momenta[0, 2, 0] == pytest.approx(70.0)
    assert runtime.seen_momenta[0, 3, 0] == pytest.approx(30.0)
    points = payload["native_probe_points"]
    assert isinstance(points, list)
    assert points[0]["relative_difference"] == pytest.approx(0.0)


def test_cli_evaluate_zero_gluon_process_uses_native_kernel(capsys) -> None:
    assert main(["evaluate", "d d~ > z", "--json"]) == 1
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > z",
        runtime_backend="auto",
    )


def test_cli_evaluate_one_gluon_process_uses_native_kernel(capsys) -> None:
    assert main(["evaluate", "d d~ > z g", "--sqrt-s", "1000", "--json"]) == 1
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > z g",
        runtime_backend="auto",
    )


def test_cli_evaluate_two_gluon_process_uses_native_dag(capsys) -> None:
    assert main(["evaluate", "d d~ > z g g", "--sqrt-s", "1000", "--json"]) == 1
    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > z g g",
        runtime_backend="auto",
    )


def test_cli_evaluate_one_gluon_process_uses_numeric_tensor_network(capsys) -> None:
    assert (
        main(
            [
                "evaluate",
                "d d~ > z g",
                "--sqrt-s",
                "1000",
                "--runtime-backend",
                "numeric-tensor-network",
                "--json",
            ]
        )
        == 1
    )

    _assert_legacy_native_command_payload(
        capsys,
        command="evaluate",
        process="d d~ > z g",
        runtime_backend="numeric-tensor-network",
    )


def test_cli_evaluate_reports_cached_evaluator_artifact(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate",
                "d d~ > z g",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="generate",
        process="d d~ > z g",
    )


def test_cli_evaluate_reports_symbolica_artifact_cross_check(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate",
                "d d~ > z",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 1
    )
    _assert_legacy_native_command_payload(
        capsys,
        command="generate",
        process="d d~ > z",
    )


def test_cli_profile_reports_artifact_and_native_runtime(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "profile",
                "d d~ > z g",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 1
    )

    _assert_legacy_native_command_payload(
        capsys,
        command="profile",
        process="d d~ > z g",
        runtime_backend="auto",
    )


def test_cli_profile_dag_evaluator_deprecates_compiled_dag_backend(capsys) -> None:
    assert (
        main(
            [
                "profile-dag-evaluator",
                "d d~ > z g",
                "--points",
                "16",
                "--repetitions",
                "1",
                "--batch-size",
                "16",
                "--symbolica-evaluator-backend",
                "jit",
                "--json",
            ]
        )
        == 1
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert payload["runtime_backend"] == "rusticol"
    assert "legacy Z+gluon profiler" in payload["error"]
    assert "generate-process PROCESS OUTPUT_DIR" in payload["error"]


def test_cli_profile_dag_evaluator_reports_unsupported_process_json(capsys) -> None:
    for process, expected_error in (
        (
            "d d~ > z z g",
            "generate-process PROCESS OUTPUT_DIR",
        ),
    ):
        assert (
            main(
                [
                    "profile-dag-evaluator",
                    process,
                    "--runtime-backend",
                    "rusticol",
                    "--json",
                ]
            )
            == 1
        )

        payload = json.loads(capsys.readouterr().out)
        assert payload["available"] is False
        assert payload["process"] == process
        assert payload["runtime_backend"] == "rusticol"
        assert expected_error in payload["error"]


def test_cli_rusticol_generate_only_delegates_to_generic_artifact_json(
    capsys,
    tmp_path: Path,
) -> None:
    for index, process in enumerate(
        (
            "d d~ > z z g",
        )
    ):
        assert (
            main(
                [
                    "profile-dag-evaluator",
                    process,
                    "--runtime-backend",
                    "rusticol",
                    "--generate-only",
                    "--save-evaluator-dir",
                    str(tmp_path / f"process-{index}"),
                    "--symbolica-evaluator-backend",
                    "jit",
                    "--json",
                ]
            )
            == 0
        )

        payload = json.loads(capsys.readouterr().out)
        assert payload["available"] is True
        assert payload["process"] == process
        assert payload["runtime_backend"] == "rusticol"
        assert payload["kind"] == "pyamplicol-generic-dag-process"
        assert payload["artifact_class"] == "generic-dag-schema-v2"
        assert payload["lowering_status"]["full_tensor_network_ready"] is True
        assert (tmp_path / f"process-{index}" / "process_manifest.json").exists()


def test_cli_generate_process_writes_multi_singlet_schema_v2_artifact_json(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate-process",
                "d d~ > z z g",
                str(tmp_path / "process"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is True
    assert payload["process"] == "d d~ > z z g"
    assert payload["runtime_backend"] == "rusticol"
    assert payload["artifact_class"] == "generic-dag-schema-v2"
    assert payload["runtime_available"] is True
    assert payload["lowering_status"]["unimplemented_vertex_kinds"] == []
    assert payload["lowering_status"]["closure_count"] > 0


def test_cli_generate_process_writes_charged_leptonic_schema_v2_json(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate-process",
                "u d~ > e+ ve g",
                str(tmp_path / "process"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is True
    assert payload["process"] == "u d~ > e+ ve g"
    assert payload["runtime_backend"] == "rusticol"
    assert payload["artifact_class"] == "generic-dag-schema-v2"
    assert payload["runtime_available"] is True
    assert payload["lowering_status"]["unimplemented_vertex_kinds"] == []
    assert payload["lowering_status"]["closure_count"] > 0


def test_cli_generate_process_set_writes_multi_singlet_schema_v2_json(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate-process",
                "d d~ > z g | d d~ > z z g",
                str(tmp_path / "process-set"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((tmp_path / "process-set" / "process_set_manifest.json").read_text())
    assert payload["available"] is True
    assert manifest["artifact_class"] == "generic-dag-schema-v2"
    assert [entry["key"] for entry in payload["generated"]] == [
        "d_dbar_to_z_g",
        "d_dbar_to_z_z_g",
    ]
    assert all(entry["runtime_available"] is True for entry in payload["generated"])
    assert all(
        entry["artifact_class"] == "generic-dag-schema-v2"
        for entry in manifest["processes"]
    )


def test_cli_generate_process_set_writes_non_family_schema_v2_json(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate-process",
                "--symbolica-evaluator-backend",
                "jit",
                "d d~ > e+ e- g | u d~ > w+ z",
                str(tmp_path / "process-set"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    manifest = json.loads((tmp_path / "process-set" / "process_set_manifest.json").read_text())
    assert payload["available"] is True
    assert manifest["artifact_class"] == "generic-dag-schema-v2"
    assert manifest["default_process_key"] == "d_dbar_to_eplus_eminus_g"
    assert [entry["key"] for entry in payload["generated"]] == [
        "d_dbar_to_eplus_eminus_g",
        "u_dbar_to_wplus_z",
    ]
    assert all(entry["runtime_available"] is True for entry in payload["generated"])
    assert all(
        entry["artifact_class"] == "generic-dag-schema-v2"
        for entry in manifest["processes"]
    )


def test_cli_generate_process_accepts_supported_full_colour_json(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate-process",
                "--color-accuracy",
                "full",
                "--symbolica-evaluator-backend",
                "jit",
                "d d~ > z g",
                str(tmp_path / "process"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is True
    assert payload["runtime_available"] is True
    assert payload["planning_status"]["color_ready"] is True
    manifest = json.loads((tmp_path / "process" / "process_manifest.json").read_text())
    contraction = manifest["runtime_schema"]["amplitude_stage"]["color_contraction"]
    assert contraction["supported"] is True
    assert contraction["includes_color_factor"] is True


def test_cli_generate_process_accepts_multi_quark_nlc_colour_class_locally(
    capsys,
    tmp_path: Path,
) -> None:
    assert (
        main(
            [
                "generate-process",
                "--color-accuracy",
                "nlc",
                "--symbolica-evaluator-backend",
                "jit",
                "d d~ > u u~ s s~",
                str(tmp_path / "process"),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is True
    assert payload["runtime_available"] is True
    assert payload["runtime_unavailable_message"] is None
    manifest = json.loads((tmp_path / "process" / "process_manifest.json").read_text())
    contraction = manifest["runtime_schema"]["amplitude_stage"]["color_contraction"]
    assert contraction["supported"] is True


def test_cli_rusticol_generate_only_supports_zero_gluon_process(
    tmp_path: Path,
) -> None:
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(pyamplicol.__file__).resolve().parents[1]),
    )
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pyamplicol",
            "profile-dag-evaluator",
            "d d~ > z",
            "--runtime-backend",
            "rusticol",
            "--generate-only",
            "--save-evaluator-dir",
            str(tmp_path / "zero"),
            "--json",
        ],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stderr + result.stdout
    payload = json.loads(result.stdout)
    assert payload["available"] is True
    assert payload["runtime_backend"] == "rusticol"
    assert payload["kind"] == "pyamplicol-generic-dag-process"
    assert payload["artifact_class"] == "generic-dag-schema-v2"
    assert payload["lowering_status"]["full_tensor_network_ready"] is True
    assert (tmp_path / "zero" / "process_manifest.json").exists()
    assert (tmp_path / "zero" / "check_standalone.py").exists()


def test_cli_benchmark_z_gluon_modes_is_legacy_command(capsys) -> None:
    assert (
        main(
            [
                "benchmark-z-gluon-modes",
                "--only-legacy-shared",
                "--min-gluons",
                "1",
                "--max-gluons",
                "1",
                "--json",
            ]
        )
        == 1
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["available"] is False
    assert payload["command"] == "benchmark-z-gluon-modes"
    assert "legacy Z+gluon-only command" in payload["error"]
    assert "time-process" in payload["error"]


def test_memory_watchdog_allows_small_command() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    script = repo_root / "pyAmpliCol" / "scripts" / "run_with_memory_watch.py"

    completed = subprocess.run(
        [
            str(script),
            "--limit-gb",
            "30",
            "--",
            "python3",
            "-c",
            "print('ok')",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert completed.stdout.strip() == "ok"


def test_memory_watchdog_writes_peak_rss_report(tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[3]
    script = repo_root / "pyAmpliCol" / "scripts" / "run_with_memory_watch.py"
    report_path = tmp_path / "watchdog.json"

    completed = subprocess.run(
        [
            str(script),
            "--limit-gb",
            "30",
            "--poll-s",
            "0.01",
            "--report-json",
            str(report_path),
            "--",
            "python3",
            "-c",
            "print('ok')",
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    assert completed.stdout.strip() == "ok"
    report = json.loads(report_path.read_text())
    assert report["returncode"] == 0
    assert report["limit_gb"] == 30
    assert report["limit_bytes"] == 30 * 1024**3
    assert isinstance(report["rss_polling_available"], bool)
    if report["peak_rss_bytes"] is not None:
        assert report["peak_rss_bytes"] > 0
        assert report["peak_rss_gb"] > 0.0


def test_cli_version_mode_does_not_import_native_dependencies(capsys) -> None:
    assert main(["--version"]) == 0

    captured = capsys.readouterr()
    assert captured.out.strip() == pyamplicol.__version__
