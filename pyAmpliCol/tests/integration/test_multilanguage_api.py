from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence
import warnings

import numpy as np
import pytest

import pyamplicol


PROJECT_ROOT = Path(__file__).resolve().parents[2]
WATCHDOG = PROJECT_ROOT / "scripts" / "run_with_memory_watch.py"
PROCESS_KEYS = ("d_dbar_to_z_g", "d_dbar_to_z_g_g")
LANGUAGES = ("python", "cpp", "fortran")


def _environment() -> dict[str, str]:
    source_root = Path(pyamplicol.__file__).resolve().parents[1]
    executable_dir = Path(sys.executable).parent
    return {
        **os.environ,
        "PYTHONPATH": str(source_root),
        "PATH": str(executable_dir) + os.pathsep + os.environ.get("PATH", ""),
    }


def _run(
    command: Sequence[str | Path],
    *,
    cwd: Path | None = None,
    timeout: float = 900.0,
) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        [str(part) for part in command],
        cwd=cwd,
        env=_environment(),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if completed.returncode != 0:
        pytest.fail(
            "command failed\n"
            f"command: {' '.join(str(part) for part in command)}\n"
            f"exit code: {completed.returncode}\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
    return completed


def _native_config() -> Path:
    executable_dir = "Scripts" if os.name == "nt" else "bin"
    candidate = Path(sys.prefix) / executable_dir / "rusticol-config"
    if candidate.exists():
        return candidate
    discovered = shutil.which("rusticol-config")
    if discovered is None:
        pytest.skip("rusticol-config is not installed")
    return Path(discovered)


def _require_native_toolchain() -> None:
    missing = [
        name
        for name in ("make", "c++", "gfortran")
        if shutil.which(name) is None
    ]
    if missing:
        pytest.skip(f"native Rusticol examples require: {', '.join(missing)}")


def _generate_process(
    output: Path,
    process: str,
    *,
    color_accuracy: str = "lc",
    selected_lc_sector: int | None = None,
    append: bool = False,
) -> None:
    command: list[str | Path] = [
        WATCHDOG,
        "--limit-gb",
        "30",
        "--",
        sys.executable,
        "-m",
        "pyamplicol",
        "generate-process",
        process,
        output,
        "--color-accuracy",
        color_accuracy,
        "--n-cores",
        "1",
        "--symbolica-evaluator-backend",
        "jit",
        "--symbolica-jit-optimization-level",
        "1",
        "--symbolica-n-cores",
        "1",
        "--batch-size",
        "16",
        "--lc-sector-strategy",
        "all",
        "--runtime-lc-sector-selector",
        "--json",
    ]
    if append:
        command.append("--append")
    if selected_lc_sector is not None:
        command.extend(("--lc-sector-ids", str(selected_lc_sector)))
    completed = _run(command)
    payload = json.loads(completed.stdout)
    assert payload["available"] is True


def _build_native_runners(root: Path) -> None:
    _require_native_toolchain()
    config = _native_config()
    for language in ("cpp", "fortran"):
        _run(
            [
                "make",
                "-C",
                root / "API" / language,
                f"RUSTICOL_CONFIG={config}",
            ]
        )


def _runner_command(root: Path, language: str) -> list[str | Path]:
    if language == "python":
        return [sys.executable, root / "API" / "python" / "check_standalone.py"]
    return [root / "API" / language / "check_standalone"]


def _run_runner(
    root: Path,
    language: str,
    *,
    process_key: str | None = None,
    extra_arguments: Sequence[str | Path] = (),
) -> dict[str, Any]:
    command = _runner_command(root, language)
    if process_key is not None:
        command.extend(("--process", process_key))
    command.extend(extra_arguments)
    command.append("--json")
    completed = _run(command)
    try:
        payload = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        pytest.fail(
            f"{language} runner emitted invalid JSON: {error}\n"
            f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
        )
    assert payload["language"] == language
    return payload


def _assert_same_result(reference: dict[str, Any], candidate: dict[str, Any]) -> None:
    for key in (
        "available",
        "process",
        "process_key",
        "color_accuracy",
        "external_particles",
        "helicities",
        "colors",
        "shape",
    ):
        assert candidate[key] == reference[key]
    for key in ("values", "resolved_sum", "compatibility_total"):
        np.testing.assert_allclose(
            np.asarray(candidate[key], dtype=np.float64),
            np.asarray(reference[key], dtype=np.float64),
            rtol=1.0e-12,
            atol=1.0e-15,
        )


def _run_all_languages(
    root: Path,
    *,
    process_key: str | None = None,
    extra_arguments: Sequence[str | Path] = (),
) -> dict[str, dict[str, Any]]:
    results = {
        language: _run_runner(
            root,
            language,
            process_key=process_key,
            extra_arguments=extra_arguments,
        )
        for language in LANGUAGES
    }
    for language in ("cpp", "fortran"):
        _assert_same_result(results["python"], results[language])
    return results


@pytest.fixture(scope="module")
def mixed_process_set(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("multilanguage-process-set") / "process-set"
    _generate_process(root, "d d~ > z g | d d~ > z g g")
    _build_native_runners(root)
    assert (root / "API" / "validation_points.dat").exists()
    assert not (root / "check_standalone.py").exists()
    for process_key in PROCESS_KEYS:
        subprocess_root = root / "subprocesses" / process_key
        assert subprocess_root.exists()
        assert not (subprocess_root / "API").exists()
    return root


def test_mixed_process_set_runners_agree_for_all_model_input_paths(
    mixed_process_set: Path,
    tmp_path: Path,
) -> None:
    model_card = tmp_path / "model-parameters.json"
    model_card.write_text(
        json.dumps({"normalization.alpha_s_me_check": [0.101, 0.0]}),
        encoding="utf-8",
    )
    scenarios: dict[str, Sequence[str | Path]] = {
        "default": (),
        "model-card": ("--model-parameters", model_card),
        "direct": (
            "--set-parameter",
            "normalization.alpha_s_me_check",
            "0.097",
            "0.0",
        ),
    }

    for process_key in PROCESS_KEYS:
        totals: dict[str, float] = {}
        for scenario, arguments in scenarios.items():
            results = _run_all_languages(
                mixed_process_set,
                process_key=process_key,
                extra_arguments=arguments,
            )
            payload = results["python"]
            assert payload["available"] is True
            assert payload["shape"][0] == 1
            assert len(payload["values"]) == int(np.prod(payload["shape"]))
            np.testing.assert_allclose(
                payload["resolved_sum"],
                payload["compatibility_total"],
                rtol=1.0e-12,
                atol=1.0e-15,
            )
            totals[scenario] = float(payload["compatibility_total"][0])
        assert totals["model-card"] != pytest.approx(totals["default"], rel=1.0e-8)
        assert totals["direct"] != pytest.approx(totals["default"], rel=1.0e-8)


def test_python_typed_selectors_and_arbitrary_precision(
    mixed_process_set: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    runtime = rusticol.Runtime.load(
        str(mixed_process_set),
        process_key="d_dbar_to_z_g_g",
    )
    physics = runtime.physics
    assert type(physics).__name__ == "ProcessPhysics"
    assert all(
        type(item).__name__ == "ExternalParticle"
        for item in physics.external_particles
    )
    assert all(
        type(item).__name__ == "HelicityConfiguration"
        for item in physics.helicities
    )
    assert all(type(item).__name__ == "ColorFlow" for item in physics.color_flows)
    assert all(
        type(item).__name__ == "ModelParameter"
        for item in physics.model_parameters
    )

    point = _validation_point(mixed_process_set, physics.process_key)
    batch = np.asarray([point], dtype=np.float64)
    helicity = next(item for item in physics.helicities if item.computed)
    color_flow = next(item for item in physics.color_flows if item.computed)
    selected = runtime.evaluate_resolved(
        batch,
        helicities=helicity,
        color_flows=color_flow,
    )
    assert type(selected).__name__ == "ResolvedEvaluation"
    assert selected.shape == (1, 1, 1)
    assert selected.helicities[0].id == helicity.id
    assert selected.color_flows[0].id == color_flow.id

    high_precision = runtime.evaluate_resolved_with_prec(
        [point],
        50,
        helicities=[helicity.id],
        color_flows=[color_flow.id],
    )
    assert high_precision.shape == (1, 1, 1)
    assert type(high_precision.values[0][0][0]).__name__ == "Decimal"

    runner_payload = _run_runner(
        mixed_process_set,
        "python",
        process_key="d_dbar_to_z_g",
        extra_arguments=("--precision", "32"),
    )
    assert runner_payload["precision"] == 32
    assert isinstance(runner_payload["values"][0], str)


def test_model_parameter_batch_is_atomic(mixed_process_set: Path) -> None:
    rusticol = pytest.importorskip("rusticol")
    runtime = rusticol.Runtime.load(
        str(mixed_process_set),
        process_key="d_dbar_to_z_g",
    )
    point = np.asarray([_validation_point(mixed_process_set, runtime.physics.process_key)])
    original = float(runtime.evaluate(point)[0])

    with pytest.raises(ValueError, match="does[.]not[.]exist.*not used"):
        runtime.set_model_parameters(
            {
                "normalization.alpha_s_me_check": (0.101, 0.0),
                "does.not.exist": (1.0, 0.0),
            }
        )
    assert float(runtime.evaluate(point)[0]) == pytest.approx(original, rel=1.0e-15)

    runtime.set_model_parameter("normalization.alpha_s_me_check", 0.101)
    assert float(runtime.evaluate(point)[0]) != pytest.approx(original, rel=1.0e-8)


def test_resolved_warnings_are_once_per_handle_and_mutable(
    mixed_process_set: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    process_key = "d_dbar_to_z_g"
    point = np.asarray([_validation_point(mixed_process_set, process_key)])

    runtime = rusticol.Runtime.load(str(mixed_process_set), process_key=process_key)
    folded = next(item for item in runtime.physics.helicities if not item.computed)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        runtime.evaluate_resolved(point, helicities=folded)
        runtime.evaluate_resolved(point, helicities=folded)
    assert len(caught) == 1
    assert "exact symmetry representative" in str(caught[0].message)

    muted_runtime = rusticol.Runtime.load(
        str(mixed_process_set),
        process_key=process_key,
    )
    muted_runtime.mute_warnings()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        muted_runtime.evaluate_resolved(point, helicities=folded.id)
    assert caught == []
    muted_runtime.unmute_warnings()
    with pytest.warns(UserWarning, match="exact symmetry representative"):
        muted_runtime.evaluate_resolved(point, helicities=folded.id)


def test_old_schema_v2_keeps_total_but_requires_regeneration_for_resolved(
    mixed_process_set: Path,
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    source = mixed_process_set / "subprocesses" / "d_dbar_to_z_g"
    root = tmp_path / "old-schema-v2"
    shutil.copytree(source, root)
    manifest_path = root / "process_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["runtime_schema"].pop("physics")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    point = _manifest_validation_batch(root)

    runtime = rusticol.Runtime.load(str(root))
    assert np.isfinite(float(runtime.evaluate(point)[0]))
    with pytest.raises(ValueError, match="no resolved physics metadata"):
        _ = runtime.physics
    with pytest.raises(ValueError, match="resolved evaluation is unavailable"):
        runtime.evaluate_resolved(point)


def test_rusticol_explicitly_rejects_schema_v1_artifact(
    mixed_process_set: Path,
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    source = mixed_process_set / "subprocesses" / "d_dbar_to_z_g"
    root = tmp_path / "schema-v1"
    shutil.copytree(source, root)
    manifest_path = root / "process_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["schema_version"] = 1
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ValueError, match="Schema-v1 execution has been removed"):
        rusticol.Runtime.load(str(root))


def test_process_set_crossing_alias_matches_direct_artifact(
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    process_set = tmp_path / "crossing-set"
    direct = tmp_path / "direct-alias"
    alias_process = "d d~ > u~ u g"
    alias_key = "d_dbar_to_ubar_u_g"

    _generate_process(
        process_set,
        "d d~ > u u~ g | d d~ > z g",
    )
    _generate_process(
        process_set,
        alias_process,
        append=True,
    )
    _generate_process(direct, alias_process)
    _build_native_runners(process_set)

    native_results = _run_all_languages(process_set, process_key=alias_key)
    assert native_results["python"]["process"] == alias_process

    point = np.asarray([_validation_point(process_set, alias_key)])
    alias_runtime = rusticol.Runtime.load(str(process_set), process_key=alias_key)
    direct_runtime = rusticol.Runtime.load(str(direct))
    alias_physics = alias_runtime.physics
    direct_physics = direct_runtime.physics
    assert alias_physics.process == direct_physics.process == alias_process
    assert [item.pdg for item in alias_physics.external_particles] == [
        item.pdg for item in direct_physics.external_particles
    ]
    assert [item.helicities for item in alias_physics.helicities] == [
        item.helicities for item in direct_physics.helicities
    ]
    assert [item.word for item in alias_physics.color_flows] == [
        item.word for item in direct_physics.color_flows
    ]

    alias_resolved = alias_runtime.evaluate_resolved(point)
    direct_resolved = direct_runtime.evaluate_resolved(point)
    assert alias_resolved.shape == direct_resolved.shape
    np.testing.assert_allclose(
        alias_resolved.values,
        direct_resolved.values,
        rtol=1.0e-12,
        atol=1.0e-15,
    )
    np.testing.assert_allclose(
        alias_runtime.evaluate(point),
        direct_runtime.evaluate(point),
        rtol=1.0e-12,
        atol=1.0e-15,
    )


def test_runners_load_metadata_without_a_validation_point(
    mixed_process_set: Path,
    tmp_path: Path,
) -> None:
    root = tmp_path / "no-validation-point"
    shutil.copytree(mixed_process_set, root)
    (root / "API" / "validation_points.dat").write_text(
        "RUSTICOL_VALIDATION_POINTS_V1\n"
        "# unavailable\td_dbar_to_z_g\tremoved by integration test\n",
        encoding="utf-8",
    )

    for language in LANGUAGES:
        payload = _run_runner(
            root,
            language,
            process_key="d_dbar_to_z_g",
        )
        assert payload["available"] is False
        assert payload["process"] == "d d~ > z g"
        assert payload["external_particles"]
        assert payload["helicities"]
        assert payload["colors"]


def test_native_runners_reject_non_f64_precision(
    mixed_process_set: Path,
) -> None:
    for language in ("cpp", "fortran"):
        command = _runner_command(mixed_process_set, language)
        command.extend(
            (
                "--process",
                "d_dbar_to_z_g",
                "--precision",
                "32",
                "--json",
            )
        )
        completed = subprocess.run(
            [str(part) for part in command],
            env=_environment(),
            capture_output=True,
            text=True,
        )
        assert completed.returncode != 0
        assert "double precision" in completed.stderr


def test_single_process_bundle_runs_without_process_key(
    tmp_path: Path,
) -> None:
    root = tmp_path / "single-process"
    _generate_process(root, "d d~ > z g")
    _build_native_runners(root)

    results = _run_all_languages(root)
    assert results["python"]["process_key"] == "d_dbar_to_z_g"
    assert not (root / "check_standalone.py").exists()


@pytest.mark.parametrize("color_accuracy", ("nlc", "full"))
def test_native_runners_expose_one_contracted_color_component(
    tmp_path: Path,
    color_accuracy: str,
) -> None:
    root = tmp_path / color_accuracy
    _generate_process(root, "d d~ > z g", color_accuracy=color_accuracy)
    _build_native_runners(root)

    results = _run_all_languages(root)
    payload = results["python"]
    assert payload["color_accuracy"] == color_accuracy
    assert payload["shape"][2] == 1
    assert payload["colors"] == [
        {"id": "contracted", "kind": "contracted", "word": []}
    ]

    rusticol = pytest.importorskip("rusticol")
    runtime = rusticol.Runtime.load(str(root))
    point = np.asarray([_validation_point(root, runtime.physics.process_key)])
    with pytest.raises(ValueError, match="LC color-flow selection is unavailable"):
        runtime.evaluate_resolved(point, color_flows=["flow:missing"])


def _validation_point(root: Path, process_key: str) -> list[list[float]]:
    path = root / "API" / "validation_points.dat"
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[0] != process_key:
            continue
        count = int(fields[1])
        values = [float(value) for value in fields[2:]]
        assert len(values) == count * 4
        return [values[index : index + 4] for index in range(0, len(values), 4)]
    raise AssertionError(f"missing validation point for {process_key}")


def _manifest_validation_batch(root: Path) -> np.ndarray:
    payload = json.loads(
        (root / "validation_momenta.json").read_text(encoding="utf-8")
    )
    return np.asarray(
        [
            [
                [float(component) for component in particle["momentum"]]
                for particle in point
            ]
            for point in payload["points"]
        ],
        dtype=np.float64,
    )
