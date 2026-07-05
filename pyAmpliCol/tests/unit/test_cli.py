from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pyamplicol
import pyamplicol.__main__ as cli
from pyamplicol.__main__ import (
    _attach_native_probe_comparison,
    _generation_build_kwargs,
    _native_generation_profile,
    _runtime_backend,
    _runtime_evaluator_kwargs,
    _z_gluon_family_passed,
    main,
    parse_args,
)
from pyamplicol.native import ExternalMomentum
from pyamplicol.processes import ProcessOptions
from pyamplicol.reference import AmplicolProbePoint, AmplicolWorkflowResult


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


def test_cli_process_parser_accepts_uppercase_bosons(capsys) -> None:
    assert main(["processes", "d d~ > Z g", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["request"]["rest"] == ["z", "g"]


def test_cli_compiled_dag_shortcut_selects_runtime_backend() -> None:
    args = parse_args(
        [
            "evaluate",
            "--compiled-dag-evaluator",
            "d d~ > z g",
        ]
    )

    assert _runtime_backend(args) == "compiled-dag"
    assert args.compiled_dag_lowering == "spenso"


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


def test_cli_profile_dag_defaults_to_fast_rusticol_cxx_o3() -> None:
    args = parse_args(["profile-dag-evaluator", "d d~ > z g"])
    kwargs = _runtime_evaluator_kwargs(args)

    assert _runtime_backend(args) == "rusticol"
    assert kwargs["symbolica_evaluator_backend"] == "compiled-complex"
    assert kwargs["symbolica_compiled_preset"] == "runtime-o3"
    assert kwargs["symbolica_n_cores"] == 10
    assert kwargs["symbolica_compiled_chunk_compile_workers"] == 10
    assert kwargs["batch_size"] == 128


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
    assert kwargs["batch_size"] == 128


def test_cli_generate_process_explicit_jit_backend_overrides_fast_default(
    tmp_path: Path,
) -> None:
    output_dir = tmp_path / "process"
    args = parse_args(
        [
            "generate-process",
            "--symbolica-evaluator-backend",
            "jit",
            "d d~ > z g",
            str(output_dir),
        ]
    )
    kwargs = _generation_build_kwargs(args, _runtime_backend(args), output_dir)

    assert kwargs["symbolica_evaluator_backend"] == "jit"
    assert kwargs["symbolica_compiled_output_dir"] == str(output_dir / "compiled")


def test_cli_time_process_defaults_to_double_precision_and_ten_seconds(
    tmp_path: Path,
) -> None:
    process_dir = tmp_path / "process"
    args = parse_args(["time-process", str(process_dir)])

    assert args.process_dir == process_dir
    assert args.precision == 16
    assert args.target_runtime == 10.0


def test_cli_compiled_dag_helicity_filter_flags() -> None:
    args = parse_args(
        [
            "generate",
            "--compiled-dag-evaluator",
            "--no-compiled-dag-helicity-filter",
            "--compiled-dag-helicity-filter-samples",
            "3",
            "--compiled-dag-helicity-filter-seed",
            "91",
            "--compiled-dag-helicity-filter-phase-space",
            "canonical",
            "d d~ > z g",
        ]
    )
    kwargs = _runtime_evaluator_kwargs(args)

    assert _runtime_backend(args) == "compiled-dag"
    assert kwargs["compiled_dag_helicity_filter"] is False
    assert kwargs["compiled_dag_helicity_filter_samples"] == 3
    assert kwargs["compiled_dag_helicity_filter_seed"] == 91
    assert kwargs["compiled_dag_helicity_filter_phase_space"] == "canonical"


def test_cli_compiled_dag_inline_external_wavefunction_flags() -> None:
    default_args = parse_args(
        [
            "profile-dag-evaluator",
            "--compiled-dag-evaluator",
            "d d~ > z g",
        ]
    )
    disabled_args = parse_args(
        [
            "profile-dag-evaluator",
            "--compiled-dag-evaluator",
            "--no-compiled-dag-inline-external-wavefunctions",
            "d d~ > z g",
        ]
    )

    assert _runtime_evaluator_kwargs(default_args)[
        "compiled_dag_inline_external_wavefunctions"
    ] is True
    assert _runtime_evaluator_kwargs(disabled_args)[
        "compiled_dag_inline_external_wavefunctions"
    ] is False


def test_cli_compiled_dag_defers_jit_direct_translation_default() -> None:
    compiled_args = parse_args(
        [
            "profile-dag-evaluator",
            "--compiled-dag-evaluator",
            "d d~ > z g",
        ]
    )
    staged_args = parse_args(
        [
            "profile-dag-evaluator",
            "d d~ > z g",
        ]
    )
    disabled_args = parse_args(
        [
            "profile-dag-evaluator",
            "--compiled-dag-evaluator",
            "--symbolica-no-jit-direct-translation",
            "d d~ > z g",
        ]
    )

    assert _runtime_evaluator_kwargs(compiled_args)[
        "symbolica_jit_direct_translation"
    ] is None
    assert _runtime_evaluator_kwargs(staged_args)[
        "symbolica_jit_direct_translation"
    ] is False
    assert _runtime_evaluator_kwargs(disabled_args)[
        "symbolica_jit_direct_translation"
    ] is False


def test_cli_compiled_dag_uses_tuned_common_pair_defaults() -> None:
    compiled_args = parse_args(
        [
            "profile-dag-evaluator",
            "--compiled-dag-evaluator",
            "d d~ > z g",
        ]
    )
    staged_args = parse_args(
        [
            "profile-dag-evaluator",
            "d d~ > z g",
        ]
    )
    explicit_args = parse_args(
        [
            "profile-dag-evaluator",
            "--compiled-dag-evaluator",
            "--symbolica-cpe-iterations",
            "5",
            "--symbolica-max-common-pair-distance",
            "75",
            "d d~ > z g",
        ]
    )

    compiled_kwargs = _runtime_evaluator_kwargs(compiled_args)
    staged_kwargs = _runtime_evaluator_kwargs(staged_args)
    explicit_kwargs = _runtime_evaluator_kwargs(explicit_args)

    assert compiled_kwargs["symbolica_cpe_iterations"] is None
    assert compiled_kwargs["symbolica_max_common_pair_distance"] == 250
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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["supported_native_target"] is True
    assert payload["backend"] == "native-python-recursion-staged"
    assert payload["graph"]["amplitudes"]
    assert payload["symbolic_lowering"]["tensor_network_probe"]["engine"] == "spenso"
    assert payload["symbolic_lowering"]["color_algebra_probe"]["engine"] == "idenso"
    assert payload["symbolic_lowering"]["recursion_plan"]["engine"] == "symbolica"
    assert payload["symbolic_lowering"]["vertex_lowering"][
        "full_tensor_network_ready"
    ] is True
    assert payload["symbolic_lowering"]["tensor_network_blueprint"][
        "full_me_tensor_network_ready"
    ] is True
    assert payload["symbolic_lowering"]["tensor_network_blueprint"][
        "expression_executed"
    ] is True
    assert payload["symbolic_lowering"]["tensor_network_blueprint"][
        "propagator_lowering_ready"
    ] is True
    cache_file = Path(payload["cache_file"])
    artifact_file = Path(payload["artifact_file"])
    artifact = json.loads(artifact_file.read_text())
    assert cache_file.exists()
    assert artifact_file.exists()
    assert payload["artifact_fingerprint"]
    assert payload["artifact_cache_hit"] is False
    assert artifact["kernel"] == "symbolica-z-gluon-tensor-network"
    assert artifact["tensor_network_scalar_evaluator_ready"] is True

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
        == 0
    )
    cached_payload = json.loads(capsys.readouterr().out)
    assert cached_payload["artifact_cache_hit"] is True
    assert cached_payload["artifact_fingerprint"] == payload["artifact_fingerprint"]


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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    artifact = json.loads(Path(payload["artifact_file"]).read_text())
    symbolic = artifact["symbolic_scalar_evaluator"]
    assert artifact["kernel"] == "symbolica-zero-gluon"
    assert artifact["symbolic_scalar_evaluator_ready"] is True
    assert symbolic["engine"] == "symbolica"
    assert symbolic["kernel"] == "symbolica-zero-gluon"
    assert symbolic["parameter_count"] == len(symbolic["parameter_names"])
    assert symbolic["evaluator_state_b64"]
    artifact_text = Path(payload["artifact_file"]).read_text()

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
        == 0
    )
    cached_payload = json.loads(capsys.readouterr().out)
    assert cached_payload["artifact_cache_hit"] is True
    assert Path(cached_payload["artifact_file"]).read_text() == artifact_text


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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    vertex_kinds = {
        interaction["vertex_kind"] for interaction in payload["graph"]["interactions"]
    }
    current_pdgs = {current["pdg"] for current in payload["graph"]["currents"]}
    assert payload["backend"] == "native-python-recursion-staged"
    assert {-21, 21}.issubset(current_pdgs)
    assert {1, 2, 3}.issubset(vertex_kinds)
    assert payload["symbolic_lowering"]["recursion_plan"][
        "tensor_route_vertex_kinds"
    ] == [1, 2, 3]


def test_cli_profile_tensor_evaluator_reports_hot_path_timing(capsys) -> None:
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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["process"] == "d d~ > z g"
    assert payload["gluon_count"] == 1
    assert payload["helicity_calls"] == 12
    assert payload["pure_evaluator_call_s"] > 0.0
    assert payload["full_evaluate_s"] > payload["pure_evaluator_call_s"]


def test_cli_compare_amplicol_dry_run_includes_timing_matching_me_test(capsys) -> None:
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
    assert payload["commands"][-1][:3] == [
        "./amplicol_generate",
        "--me_test=11",
        "--timing=11",
    ]
    assert payload["pyamplicol_cache_dir"]


def test_cli_compare_amplicol_dry_run_can_use_direct_probe(capsys) -> None:
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
    assert payload["mode"] == "amplicol_fixed_probe"
    assert payload["commands"][-1][:3] == [
        "./amplicol_generate",
        "--amplicol_fixed_probe=5",
        "--timing=5",
    ]
    assert all("--library=create" not in command for command in payload["commands"])


def test_cli_validate_z_gluon_family_dry_run_summarizes_probe_commands(
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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["family"] == "d d~ -> Z + n g"
    assert payload["max_gluons"] == 2
    assert payload["points"] == 7
    assert payload["rel_tol"] == 1e-08
    assert payload["pyamplicol_runtime_backend"] == "dag"
    assert [row["gluon_count"] for row in payload["rows"]] == [0, 1, 2]
    assert payload["rows"][0]["mode"] == "amplicol_fixed_probe"
    assert payload["rows"][0]["runtime_backend"] == "python"
    assert payload["rows"][0]["commands"][-1][:3] == [
        "./amplicol_generate",
        "--amplicol_fixed_probe=7",
        "--timing=7",
    ]
    assert payload["rows"][1]["mode"] == "amplicol_probe"
    assert payload["rows"][1]["runtime_backend"] == "dag"
    assert payload["rows"][1]["commands"][-1][:3] == [
        "./amplicol_generate",
        "--amplicol_probe=7",
        "--timing=7",
    ]


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


def test_native_generation_profile_records_artifact_performance(
    tmp_path: Path,
) -> None:
    profile = _native_generation_profile(
        "d d~ > z g",
        tmp_path,
        ProcessOptions(),
    )

    assert profile["backend"] == "native-python-recursion-staged"
    assert profile["kernel"] == "symbolica-one-gluon-tensor-network"
    assert profile["artifact_cache_hit"] is False
    assert isinstance(profile["generation_time_s"], float)
    assert profile["generation_time_s"] >= 0.0
    assert profile["artifact_load_s"] is not None
    assert isinstance(profile["artifact_file"], str)
    assert profile["full_me_tensor_network_ready"] is True
    assert profile["graph_counts"] == {
        "currents": 12,
        "interactions": 8,
        "amplitudes": 2,
    }
    blueprint = profile["tensor_network_blueprint"]
    assert isinstance(blueprint, dict)
    assert blueprint["engine"] == "spenso"
    assert blueprint["expression_built"] is True
    assert blueprint["expression_executed"] is True
    assert blueprint["full_me_tensor_network_ready"] is True
    assert Path(profile["artifact_file"]).exists()

    cached = _native_generation_profile(
        "d d~ > z g",
        tmp_path,
        ProcessOptions(),
    )
    assert cached["artifact_cache_hit"] is True
    assert cached["artifact_fingerprint"] == profile["artifact_fingerprint"]


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


def test_cli_evaluate_zero_gluon_process_uses_native_kernel(capsys) -> None:
    assert main(["evaluate", "d d~ > z", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["evaluation_available"] is True
    assert payload["backend"] == "native-python-zero-gluon"
    assert payload["matrix_element"] == 142.10653372872991


def test_cli_evaluate_one_gluon_process_uses_native_kernel(capsys) -> None:
    assert main(["evaluate", "d d~ > z g", "--sqrt-s", "1000", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["evaluation_available"] is True
    assert payload["backend"] == "native-spenso-symbolica-shared-helicity-current-dag"
    assert payload["kernel"] == "spenso-symbolica-shared-helicity-current-dag"
    assert payload["matrix_element"] > 0.0


def test_cli_evaluate_two_gluon_process_uses_native_dag(capsys) -> None:
    assert main(["evaluate", "d d~ > z g g", "--sqrt-s", "1000", "--json"]) == 0

    payload = json.loads(capsys.readouterr().out)
    assert payload["evaluation_available"] is True
    assert payload["backend"] == "native-spenso-symbolica-shared-helicity-current-dag"
    assert payload["kernel"] == "spenso-symbolica-shared-helicity-current-dag"
    assert payload["identical_factor"] == 2
    assert payload["matrix_element"] > 0.0


def test_cli_evaluate_six_gluon_process_uses_native_dag(capsys) -> None:
    assert (
        main(
            [
                "evaluate",
                "d d~ > z g g g g g g",
                "--sqrt-s",
                "1000",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["evaluation_available"] is True
    assert payload["backend"] == "native-spenso-symbolica-shared-helicity-current-dag"
    assert payload["kernel"] == "spenso-symbolica-shared-helicity-current-dag"
    assert payload["native_runtime_backend"]["evaluator_metadata"]["gluon_count"] == 6
    assert payload["matrix_element"] > 0.0


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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["evaluation_available"] is True
    assert payload["backend"] == "native-spenso-numeric-tensor-network"
    assert payload["kernel"] == "spenso-numeric-tensor-network"
    assert payload["native_runtime_backend"]["evaluator_metadata"]["gluon_count"] == 1
    assert payload["matrix_element"] > 0.0


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
        == 0
    )
    capsys.readouterr()

    assert (
        main(
            [
                "evaluate",
                "d d~ > z g",
                "--sqrt-s",
                "1000",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    artifact = payload["evaluator_artifact"]
    assert artifact["available"] is True
    assert artifact["kernel"] == "symbolica-one-gluon-tensor-network"
    assert artifact["full_me_tensor_network_ready"] is True
    assert artifact["tensor_network_scalar_evaluator_ready"] is True
    assert Path(artifact["file"]).exists()
    assert artifact["load_time_s"] >= 0.0
    symbolic = payload["symbolic_evaluation"]
    assert symbolic["available"] is True
    assert symbolic["kernel"] == "symbolica-one-gluon-tensor-network"
    assert symbolic["relative_difference"] < 1e-12
    assert abs(symbolic["matrix_element"] - payload["matrix_element"]) < 1e-10


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
        == 0
    )
    capsys.readouterr()

    assert (
        main(
            [
                "evaluate",
                "d d~ > z",
                "--cache-dir",
                str(tmp_path),
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    artifact = payload["evaluator_artifact"]
    symbolic = payload["symbolic_evaluation"]
    assert artifact["kernel"] == "symbolica-zero-gluon"
    assert artifact["symbolic_scalar_evaluator_ready"] is True
    assert symbolic["available"] is True
    assert symbolic["kernel"] == "symbolica-zero-gluon"
    assert symbolic["relative_difference"] < 1e-12
    assert abs(symbolic["matrix_element"] - payload["matrix_element"]) < 1e-10


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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    profile = payload["profile"]
    assert payload["artifact_file"] is not None
    assert profile["artifact_load_s"] is not None
    assert profile["per_point_runtime_s"] is not None
    assert profile["native_runtime_backend"]["backend"] == (
        "native-spenso-symbolica-shared-helicity-current-dag"
    )
    assert profile["tensor_network_reduction_s"] is not None
    assert profile["symbolica_optimization_jit_s"] is not None
    assert profile["per_point_runtime_s"] >= 0.0


def test_cli_profile_dag_evaluator_supports_compiled_dag_backend(capsys) -> None:
    assert (
        main(
            [
                "profile-dag-evaluator",
                "--compiled-dag-evaluator",
                "--no-compiled-dag-helicity-filter",
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
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["runtime_backend"] == "compiled-dag"
    assert payload["metadata"]["effective_lowering"] == "spenso"
    assert payload["metadata"]["spenso_body_count"] > 0
    assert payload["runtime_evaluator_only_us_per_point"] >= 0.0


def test_cli_profile_dag_evaluator_reports_unsupported_process_json(capsys) -> None:
    for process in (
        "u d~ > w+ g g",
        "d d~ > e+ e- g g",
        "d d~ > z z g",
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
        assert "currently supports q q~ -> Z" in payload["error"]


def test_cli_rusticol_generate_only_reports_unsupported_process_json(
    capsys,
    tmp_path: Path,
) -> None:
    for index, process in enumerate(
        (
            "u d~ > w+ g g",
            "d d~ > e+ e- g g",
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
                    "--json",
                ]
            )
            == 1
        )

        payload = json.loads(capsys.readouterr().out)
        assert payload["available"] is False
        assert payload["process"] == process
        assert payload["runtime_backend"] == "rusticol"
        assert payload["error"] == f"no native graph available for {process}"


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
    assert payload["metadata"]["kernel"] == "symbolica-zero-gluon"
    assert (tmp_path / "zero" / "process_manifest.json").exists()
    assert (tmp_path / "zero" / "check_standalone.py").exists()


def test_cli_benchmark_routes_compiled_dag_backend(monkeypatch, capsys) -> None:
    captured: dict[str, object] = {}

    def fake_benchmark_z_gluon_modes(**kwargs):
        captured.update(kwargs)
        return {
            "rows": [],
            "summary": {
                "mode_success_counts": {
                    "legacy": 0,
                    "python": 0,
                    "numeric_tn": 0,
                    "parametric_tn": 0,
                    "shared_dag": 0,
                },
                "max_relative_difference_to_legacy": None,
                "all_four_modes_match_for_all_rows": False,
            },
        }

    monkeypatch.setattr(
        cli,
        "benchmark_z_gluon_modes",
        fake_benchmark_z_gluon_modes,
    )

    assert (
        main(
            [
                "benchmark-z-gluon-modes",
                "--compiled-dag-evaluator",
                "--only-legacy-shared",
                "--min-gluons",
                "1",
                "--max-gluons",
                "1",
                "--json",
            ]
        )
        == 0
    )

    payload = json.loads(capsys.readouterr().out)
    assert payload["rows"] == []
    evaluator_kwargs = captured["evaluator_build_kwargs"]
    assert isinstance(evaluator_kwargs, dict)
    assert evaluator_kwargs["runtime_backend"] == "compiled-dag"


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
