from __future__ import annotations

import json
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest

import pyamplicol.dag_runtime as dag_runtime_module
from pyamplicol.dag_runtime import (
    _RUSTICOL_STANDALONE_CHECK_SCRIPT,
    _build_shared_global_parameter_layout,
    _build_shared_helicity_current_table,
    _rusticol_layout_manifest,
    _rusticol_model_manifest,
    _rusticol_table_manifest,
    _write_zero_gluon_rusticol_process_artifacts,
    _z_gluon_graph,
    ZGluonDAGEvaluator,
)
from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.native import LeadingColorZJetsNativeEvaluator, NativeEvaluationError
from pyamplicol.process_runtime import load_process, load_process_manifest


def test_rusticol_process_manifest_helpers_are_process_shaped() -> None:
    model = AmplicolSMLeadingColorModel()
    graph = _z_gluon_graph("d d~ > z g", model)
    table = _build_shared_helicity_current_table(graph, model, gluon_count=1)
    layout = _build_shared_global_parameter_layout(table)

    layout_payload = _rusticol_layout_manifest(table, layout)
    table_payload = _rusticol_table_manifest(table)

    json.dumps(layout_payload)
    json.dumps(table_payload)
    model_payload = _rusticol_model_manifest(model)
    json.dumps(model_payload)

    assert layout_payload["parameter_count"] == layout.parameter_count
    assert len(layout_payload["current_offsets"]) == len(table.currents)
    assert table_payload["currents"][0]["id"] == 0
    assert table_payload["sources"]
    assert table_payload["amplitudes"]
    particles_by_pdg = {
        particle["pdg"]: particle
        for particle in model_payload["particles"]
    }
    assert {1, 21, 23}.issubset(particles_by_pdg)
    assert particles_by_pdg[1]["anti_pdg"] == -1
    assert model_payload["vertices"]
    assert "import rusticol" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "rusticol.Runtime.load(str(root))" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "--precision" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "--profile" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "--rusticol-folder" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "runtime.evaluate_with_prec(points, precision)" in (
        _RUSTICOL_STANDALONE_CHECK_SCRIPT
    )
    assert "runtime.profile(batch, precision=precision)" in (
        _RUSTICOL_STANDALONE_CHECK_SCRIPT
    )


def test_load_process_manifest_rejects_wrong_kind(tmp_path: Path) -> None:
    process_dir = tmp_path / "process"
    process_dir.mkdir()
    (process_dir / "process_manifest.json").write_text(
        json.dumps({"kind": "not-a-pyamplicol-process"}),
        encoding="utf-8",
    )

    with pytest.raises(NativeEvaluationError, match="unsupported process artifact kind"):
        load_process_manifest(process_dir)


def test_rusticol_load_rejects_incomplete_manifest(tmp_path: Path) -> None:
    _skip_if_rusticol_unavailable()
    process_dir = tmp_path / "process"
    process_dir.mkdir()
    (process_dir / "process_manifest.json").write_text("{}", encoding="utf-8")

    result = _run_rusticol_load_subprocess(process_dir)
    assert result.returncode == 0, result.stderr
    assert "could not parse process manifest" in result.stdout


def test_rusticol_load_rejects_wrong_kind_and_schema(tmp_path: Path) -> None:
    _skip_if_rusticol_unavailable()
    for name, patch, expected in (
        (
            "wrong-kind",
            {"kind": "not-a-pyamplicol-process", "schema_version": 1},
            "unsupported process artifact kind",
        ),
        (
            "wrong-schema",
            {"kind": "pyamplicol-rusticol-process", "schema_version": 999},
            "unsupported process manifest schema_version",
        ),
    ):
        process_dir = tmp_path / name
        process_dir.mkdir()
        payload = _minimal_process_manifest_payload()
        payload.update(patch)
        (process_dir / "process_manifest.json").write_text(
            json.dumps(payload),
            encoding="utf-8",
        )

        result = _run_rusticol_load_subprocess(process_dir)
        assert result.returncode == 0, result.stderr
        assert expected in result.stdout


def test_rusticol_load_rejects_unsupported_family(tmp_path: Path) -> None:
    _skip_if_rusticol_unavailable()
    process_dir = tmp_path / "process"
    process_dir.mkdir()
    payload = _minimal_process_manifest_payload()
    payload.update(
        {
            "process": "u d~ > w+ g g",
            "family": "unsupported-leading-color-family",
            "gluon_count": 2,
            "external_pdg_order": [2, -1, 24, 21, 21],
        }
    )
    (process_dir / "process_manifest.json").write_text(
        json.dumps(payload),
        encoding="utf-8",
    )

    result = _run_rusticol_load_subprocess(process_dir)
    assert result.returncode == 0, result.stderr
    assert "currently supports q-qbar-z-gluons" in result.stdout


def test_zero_gluon_process_artifact_loads_in_python_and_rusticol(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()
    model = AmplicolSMLeadingColorModel()
    process_dir = tmp_path / "zero"
    manifest_path = _write_zero_gluon_rusticol_process_artifacts(
        process_dir,
        process="d d~ > z",
        model=model,
    )
    manifest = load_process_manifest(process_dir)
    assert manifest_path.name == "process_manifest.json"
    assert manifest.gluon_count == 0
    assert manifest.compiled_kind == "zero-gluon-symbolic-scalar"
    manifest_payload = json.loads(manifest_path.read_text())
    assert manifest_payload["backend_settings"]["runtime"] == "rusticol"
    assert (
        manifest_payload["backend_settings"]["compiled_kind"]
        == "zero-gluon-symbolic-scalar"
    )
    particles_by_pdg = {
        particle["pdg"]: particle
        for particle in manifest_payload["model"]["particles"]
    }
    assert {1, 23}.issubset(particles_by_pdg)
    assert particles_by_pdg[1]["anti_pdg"] == -1
    assert manifest_payload["model"]["vertices"]
    assert "dependency_fingerprint" in manifest_payload
    assert "symbolica" in manifest_payload["dependency_fingerprint"]
    assert manifest_payload["helicity_filter"]["retained_amplitude_count"] == 1
    assert manifest_payload["helicity_filter"]["raw_sum_weights"] == [1.0]
    assert [particle["pdg"] for particle in manifest_payload["external_particles"]] == [
        1,
        -1,
        23,
    ]
    assert [particle["role"] for particle in manifest_payload["external_particles"]] == [
        "initial",
        "initial",
        "final",
    ]
    assert manifest_payload["momentum_conventions"]["input_shape"] == [
        "batch",
        3,
        4,
    ]
    assert manifest_payload["momentum_conventions"]["incoming_labels"] == [1, 2]
    assert manifest_payload["momentum_conventions"]["final_state_labels"] == [3]

    payload = json.loads((process_dir / "validation_momenta.json").read_text())
    momenta = [
        [
            [float(component) for component in particle["momentum"]]
            for particle in point
        ]
        for point in payload["points"]
    ]
    native = LeadingColorZJetsNativeEvaluator(model)
    reference = native.evaluate(
        "d d~ > z",
        particles=native.canonical_zero_gluon_point("d d~ > z", sqrt_s=1000.0),
    ).matrix_element
    python_runtime = load_process(process_dir, runtime="python")
    python_value = python_runtime.evaluate(momenta)[0]
    assert abs(python_value - reference) / abs(reference) < 1.0e-15
    assert python_runtime.metadata["runtime"] == "python-zero-gluon-symbolic-scalar"
    python_profile = python_runtime.profile(momenta)
    assert python_profile["points"] == 1
    assert python_profile["values"] == [python_value]
    assert python_profile["timing"] is None
    python_diagnostics = python_runtime.stage_diagnostics(momenta)
    assert python_diagnostics["stages"][0]["stage"] == "zero_gluon_symbolic_scalar"
    with pytest.raises(NativeEvaluationError, match=r"\(batch, 3, 4\)"):
        python_runtime.evaluate([momenta[0][:2]])

    code = f"""
import json
from decimal import Decimal
from pathlib import Path
import numpy as np
import pyamplicol
import rusticol

root = Path({str(process_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)
decimal_momenta = [
    [[Decimal(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
]
runtime = rusticol.Runtime.load(str(root))
runtime_from_pyamplicol = pyamplicol.load_process(root, runtime="rusticol")
value = float(runtime.evaluate(momenta)[0])
value_from_pyamplicol = float(runtime_from_pyamplicol.evaluate(momenta)[0])
value16 = float(runtime.evaluate_with_prec(momenta, 16)[0])
value32 = runtime.evaluate_with_prec(decimal_momenta, 32)[0]
value50 = runtime.evaluate_with_prec(decimal_momenta, 50)[0]
profile = dict(runtime.profile(momenta, precision=16))
batched_momenta = np.repeat(momenta, 3, axis=0)
batched_values = runtime.evaluate(batched_momenta)
batched_profile = dict(runtime.profile(batched_momenta, precision=16))
assert abs(value - {reference!r}) / abs({reference!r}) < 1.0e-14
assert value_from_pyamplicol == value
assert value == value16
assert value32.__class__.__name__ == "Decimal"
assert value50.__class__.__name__ == "Decimal"
assert profile["points"] == 1
assert profile["core_evaluator_time_s"] >= 0.0
assert len(batched_values) == 3
assert all(abs(float(item) - value) / abs(value) < 1.0e-14 for item in batched_values)
assert batched_profile["points"] == 3
assert batched_profile["batch_size"] == 3

bad_momenta = momenta[:, :2, :]
bad_decimal_momenta = [point[:2] for point in decimal_momenta]
for call in (
    lambda: runtime.evaluate(bad_momenta),
    lambda: runtime.evaluate_with_prec(bad_momenta, 16),
    lambda: runtime.profile(bad_momenta, precision=16),
    lambda: runtime.evaluate_with_prec(bad_decimal_momenta, 32),
    lambda: runtime.evaluate_with_prec(bad_decimal_momenta, 50),
):
    try:
        call()
    except ValueError as exc:
        assert "expected 3" in str(exc) or "(batch, 3, 4)" in str(exc)
    else:
        raise AssertionError("rusticol accepted a momentum batch with the wrong leg count")
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout

    rusticol_spec = importlib.util.find_spec("rusticol")
    assert rusticol_spec is not None and rusticol_spec.origin is not None
    rusticol_folder = Path(rusticol_spec.origin).parent
    for args in (
        ("--precision", "16"),
        ("--precision", "32"),
        ("--precision", "50"),
        ("--precision", "16", "--profile", "--target-runtime", "0.001"),
        ("--precision", "16", "--rusticol-folder", str(rusticol_folder)),
    ):
        standalone = subprocess.run(
            [sys.executable, str(process_dir / "check_standalone.py"), *args],
            check=False,
            env=env,
            text=True,
            capture_output=True,
        )
        assert standalone.returncode == 0, standalone.stderr + standalone.stdout
        assert "d d~ > z" in standalone.stdout
        assert "family" in standalone.stdout
        assert "gluons" in standalone.stdout
        assert "backend" in standalone.stdout
        assert "retained amplitudes" in standalone.stdout


def test_python_process_runtime_loads_nonzero_manifest_without_regenerating_graph(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _skip_if_rusticol_unavailable()
    model = AmplicolSMLeadingColorModel()
    process = "d d~ > z g"
    process_dir = tmp_path / "n1"
    stale_root_evaluator = process_dir / "shared_stage_stale.evaluator.bin"
    process_dir.mkdir()
    stale_root_evaluator.write_bytes(b"stale")
    generator = ZGluonDAGEvaluator(
        process,
        model=model,
        batch_size=16,
        symbolica_evaluator_backend="jit",
    )
    generator.save_evaluator_artifact(process_dir)
    assert list(process_dir.glob("*.evaluator.bin")) == []
    evaluator_files = sorted((process_dir / "evaluators").glob("*.evaluator.bin"))
    assert evaluator_files
    manifest_payload = json.loads((process_dir / "process_manifest.json").read_text())
    assert manifest_payload["backend_settings"]["runtime"] == "rusticol"
    assert manifest_payload["backend_settings"]["compiled_kind"] == "shared-compiled-sweep"
    assert (
        manifest_payload["backend_settings"]["evaluator_manifest"]
        == "manifest.json"
    )
    assert "symbolica_evaluator_settings" in manifest_payload["backend_settings"]
    assert "dependency_fingerprint" in manifest_payload
    assert "install_manifest" in manifest_payload["dependency_fingerprint"]
    particles_by_pdg = {
        particle["pdg"]: particle
        for particle in manifest_payload["model"]["particles"]
    }
    assert {1, 21, 23}.issubset(particles_by_pdg)
    assert particles_by_pdg[1]["anti_pdg"] == -1
    assert manifest_payload["model"]["vertices"]
    assert [particle["label"] for particle in manifest_payload["external_particles"]] == [
        1,
        2,
        3,
        4,
    ]
    assert [particle["role"] for particle in manifest_payload["external_particles"]] == [
        "initial",
        "initial",
        "final",
        "final",
    ]
    assert manifest_payload["momentum_conventions"]["input_shape"] == [
        "batch",
        4,
        4,
    ]
    assert manifest_payload["momentum_conventions"]["source_current_crossing"][
        "crossed_incoming_labels"
    ] == [1, 2]
    helicity_filter = manifest_payload["helicity_filter"]
    assert helicity_filter["strategy"] == "shared-current-retained-amplitudes"
    assert helicity_filter["retained_amplitude_count"] == len(
        generator.table.amplitudes
    )
    assert len(helicity_filter["amplitudes"]) == len(generator.table.amplitudes)
    point = LeadingColorZJetsNativeEvaluator(model).canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    reference = generator.evaluate(point).matrix_element
    momenta = [
        [
            [float(component) for component in particle.momentum]
            for particle in point
        ]
    ]

    def fail_graph_rebuild(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        raise AssertionError("process runtime rebuilt the graph from the process string")

    monkeypatch.setattr(dag_runtime_module, "_z_gluon_graph", fail_graph_rebuild)
    runtime = load_process(process_dir, runtime="python")
    metadata = runtime.metadata
    assert metadata["loaded_from_process_manifest"] is True
    assert metadata["current_count"] == generator.metadata.shared_current_count
    value = runtime.evaluate(momenta)[0]
    assert abs(value - reference) / abs(reference) < 1.0e-12
    profile = runtime.profile(momenta)
    assert profile["points"] == 1
    assert profile["values"] == [value]
    assert isinstance(profile["timing"], dict)
    diagnostics = runtime.stage_diagnostics(momenta)
    assert diagnostics["points"] == 1
    assert diagnostics["stages"][0]["stage"] == "sources"

    code = f"""
import json
import numpy as np
import pyamplicol
import rusticol
from decimal import Decimal

runtime = rusticol.Runtime.load({str(process_dir)!r})
runtime_from_pyamplicol = pyamplicol.load_process({str(process_dir)!r}, runtime="rusticol")
momenta = np.asarray({momenta!r}, dtype=np.float64)
decimal_momenta = [
    [[Decimal(str(component)) for component in particle] for particle in point]
    for point in {momenta!r}
]
single = float(runtime.evaluate(momenta)[0])
single_from_pyamplicol = float(runtime_from_pyamplicol.evaluate(momenta)[0])
value32 = runtime.evaluate_with_prec(decimal_momenta, 32)[0]
value50 = runtime.evaluate_with_prec(decimal_momenta, 50)[0]
batch = np.repeat(momenta, 4, axis=0)
values = runtime.evaluate(batch)
profile = dict(runtime.profile(batch, precision=16))
assert single_from_pyamplicol == single
assert abs(single - {reference!r}) / abs({reference!r}) < 1.0e-12
assert value32.__class__.__name__ == "Decimal"
assert value50.__class__.__name__ == "Decimal"
assert abs(float(value32) - {reference!r}) / abs({reference!r}) < 1.0e-12
assert abs(float(value50) - {reference!r}) / abs({reference!r}) < 1.0e-12
assert values.shape == (4,)
assert all(abs(float(value) - single) / abs(single) < 1.0e-12 for value in values)
assert profile["points"] == 4
assert profile["batch_size"] == 4
assert "current_rss_bytes" in profile
assert "peak_rss_bytes" in profile
diagnostics = dict(runtime.stage_diagnostics(momenta))
diagnostics["stages"] = [dict(stage) for stage in diagnostics["stages"]]
print(json.dumps(diagnostics, sort_keys=True))
"""
    batch_env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    batch_result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=batch_env,
        text=True,
        capture_output=True,
    )
    assert batch_result.returncode == 0, batch_result.stderr + batch_result.stdout
    rust_diagnostics = json.loads(batch_result.stdout)
    assert rust_diagnostics["points"] == diagnostics["points"]
    assert len(rust_diagnostics["stages"]) == len(diagnostics["stages"])
    for rust_stage, python_stage in zip(
        rust_diagnostics["stages"],
        diagnostics["stages"],
        strict=True,
    ):
        assert rust_stage["stage"] == python_stage["stage"]
        assert rust_stage["output_len"] == python_stage["output_len"]
        for key in ("sum_re", "sum_im", "sum_abs2", "max_abs"):
            assert abs(float(rust_stage[key]) - float(python_stage[key])) < 1.0e-9

    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    for args in (
        ("--precision", "16"),
        ("--precision", "32"),
        ("--precision", "16", "--profile", "--target-runtime", "0.001"),
    ):
        standalone = subprocess.run(
            [sys.executable, str(process_dir / "check_standalone.py"), *args],
            check=False,
            env=env,
            text=True,
            capture_output=True,
        )
        assert standalone.returncode == 0, standalone.stderr + standalone.stdout
        assert process in standalone.stdout
        assert "family" in standalone.stdout
        assert "retained amplitudes" in standalone.stdout


def test_python_and_rusticol_stage_diagnostics_match_when_artifact_exists() -> None:
    _skip_if_rusticol_unavailable()
    artifact_dir = _local_rusticol_smoke_artifact()
    if not (artifact_dir / "process_manifest.json").exists():
        pytest.skip("local rusticol smoke process artifact is not available")

    code = f"""
import json
from pathlib import Path
import numpy as np
import rusticol
import pyamplicol

root = Path({str(artifact_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)

rust_runtime = rusticol.Runtime.load(str(root))
python_runtime = pyamplicol.load_process(root, runtime="python")
rust_diag = dict(rust_runtime.stage_diagnostics(momenta))
python_diag = python_runtime.stage_diagnostics(momenta)

assert rust_diag["points"] == python_diag["points"]
assert len(rust_diag["stages"]) == len(python_diag["stages"])
for rust_stage, python_stage in zip(rust_diag["stages"], python_diag["stages"], strict=True):
    rust_stage = dict(rust_stage)
    assert rust_stage["stage"] == python_stage["stage"]
    assert rust_stage["output_len"] == python_stage["output_len"]
    for key in ("sum_re", "sum_im", "sum_abs2", "max_abs"):
        assert abs(float(rust_stage[key]) - float(python_stage[key])) < 1.0e-9
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


def test_rusticol_api_smoke_when_artifact_exists() -> None:
    _skip_if_rusticol_unavailable()
    artifact_dir = _local_rusticol_smoke_artifact()
    if not (artifact_dir / "process_manifest.json").exists():
        pytest.skip("local rusticol smoke process artifact is not available")

    code = f"""
import json
from decimal import Decimal
from pathlib import Path
import numpy as np
import rusticol

root = Path({str(artifact_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
decimal_momenta = [
    [[Decimal(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
]
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)

runtime = rusticol.Runtime.load(str(root))
values_a = runtime.evaluate(momenta)
values_b = runtime.evaluate(momenta)
values_16 = runtime.evaluate_with_prec(momenta, 16)
values_32 = runtime.evaluate_with_prec(decimal_momenta, 32)
profile = dict(runtime.profile(momenta, precision=16))

assert isinstance(values_a, np.ndarray)
assert values_a.dtype == np.float64
assert np.array_equal(values_a, values_b)
assert np.array_equal(values_a, values_16)
assert values_32 and values_32[0].__class__.__name__ == "Decimal"
for key in (
    "points",
    "batch_size",
    "values",
    "source_fill_time_s",
    "momentum_setup_time_s",
    "parameter_pack_time_s",
    "stage_evaluator_time_s",
    "amplitude_evaluator_time_s",
    "result_reduction_time_s",
    "output_transfer_time_s",
    "total_time_s",
    "core_evaluator_time_s",
    "peak_rss_bytes",
):
    assert key in profile
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


def _local_rusticol_smoke_artifact() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "pyAmpliCol"
        / "outputs"
        / "rusticol_high_precision_smoke"
        / "n1"
    )


def _skip_if_rusticol_unavailable() -> None:
    if importlib.util.find_spec("rusticol") is None:
        pytest.skip("rusticol extension is not installed")


def _run_rusticol_load_subprocess(process_dir: Path) -> subprocess.CompletedProcess[str]:
    code = (
        "import rusticol\n"
        "try:\n"
        f"    rusticol.Runtime.load({str(process_dir)!r})\n"
        "except Exception as exc:\n"
        "    print(type(exc).__name__ + ': ' + str(exc))\n"
        "else:\n"
        "    raise SystemExit('rusticol load unexpectedly succeeded')\n"
    )
    env = dict(os.environ)
    return subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )


def _minimal_process_manifest_payload() -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "pyamplicol-rusticol-process",
        "process": "d d~ > z g",
        "family": "q-qbar-z-gluons-leading-color",
        "gluon_count": 1,
        "external_pdg_order": [1, -1, 23, 21],
        "model": {
            "alpha_s_me_check": 0.13,
            "alpha_ew": 0.0075,
            "mass_z": 91.188,
        },
        "normalization": {
            "color_factor": 1.0,
            "average_factor": 1.0,
            "identical_factor": 1.0,
            "coupling_factor": 1.0,
        },
        "layout": {
            "parameter_count": 0,
            "current_offsets": [],
            "momentum_offsets_and_labels": [],
        },
        "table": {
            "currents": [],
            "sources": [],
            "amplitudes": [],
        },
        "compiled": {
            "kind": "shared-compiled-sweep",
            "stages": [],
            "amplitude_stage": {
                "output_length": 0,
                "raw_sum_weights": [],
                "amplitude_evaluator": None,
                "raw_sum_evaluator": None,
            },
        },
    }
