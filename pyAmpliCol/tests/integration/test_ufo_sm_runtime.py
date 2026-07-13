from __future__ import annotations

import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

import pyamplicol

from pyamplicol.model_assets import bundled_model_path


def _generate_process(
    output: Path,
    *,
    model: Path | None = None,
    process: str | None = None,
    color_accuracy: str = "lc",
    lc_sector_id: int = 0,
    max_coupling_orders: tuple[str, ...] = (),
    model_restriction: str = "default",
) -> dict[str, object]:
    selected_process = process or (
        "d d~ > Z g" if model is not None else "d d~ > z g"
    )
    arguments = [
        "generate-process",
        selected_process,
        str(output),
        "--color-accuracy",
        color_accuracy,
        "--symbolica-evaluator-backend",
        "jit",
        "--symbolica-jit-optimization-level",
        "1",
        "--symbolica-n-cores",
        "1",
        "--batch-size",
        "16",
        "--model-cache-dir",
        str(output.parent / "model-cache"),
        "--json",
    ]
    if color_accuracy == "lc":
        arguments[3:3] = ["--lc-sector-ids", str(lc_sector_id)]
    for coupling_order in max_coupling_orders:
        arguments.extend(("--max-coupling-order", coupling_order))
    if model is not None:
        arguments[1:1] = ["--model", str(model)]
        if model_restriction != "default":
            arguments[3:3] = ["--model-restriction", model_restriction]

    environment = dict(
        os.environ,
        PYTHONPATH=str(Path(pyamplicol.__file__).resolve().parents[1]),
    )
    completed = subprocess.run(
        [sys.executable, "-m", "pyamplicol", *arguments],
        check=True,
        capture_output=True,
        env=environment,
        text=True,
    )
    return json.loads(completed.stdout)


@pytest.mark.parametrize("color_accuracy", ("nlc", "full"))
def test_external_model_generation_preserves_requested_color_accuracy(
    tmp_path: Path,
    color_accuracy: str,
) -> None:
    output = tmp_path / color_accuracy
    _generate_process(
        output,
        model=bundled_model_path("sm", "json"),
        process="d d~ > t t~",
        color_accuracy=color_accuracy,
        max_coupling_orders=("QED=0",),
    )

    manifest = json.loads(
        (output / "process_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["process_ir"]["color_accuracy"] == color_accuracy
    assert manifest["normalization"]["color_accuracy"] == color_accuracy


@pytest.mark.parametrize(
    ("process", "color_accuracy"),
    (
        ("d d~ > t t~", "nlc"),
        ("d d~ > t t~", "full"),
        ("g g > t t~", "nlc"),
        ("g g > t t~", "full"),
        ("g g > g g", "nlc"),
        ("g g > g g", "full"),
    ),
)
def test_ufo_sm_qcd_color_modes_match_builtin(
    tmp_path: Path,
    process: str,
    color_accuracy: str,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime parity requires a release Rusticol extension")

    builtin_dir = tmp_path / "builtin"
    ufo_dir = tmp_path / "ufo"
    options = {
        "process": process,
        "color_accuracy": color_accuracy,
        "max_coupling_orders": ("QED=0",),
    }
    _generate_process(builtin_dir, **options)
    _generate_process(
        ufo_dir,
        model=bundled_model_path("sm", "json"),
        **options,
    )

    momenta = _validation_momenta(ufo_dir)
    builtin_value = float(rusticol.Runtime.load(str(builtin_dir)).evaluate(momenta)[0])
    ufo_value = float(rusticol.Runtime.load(str(ufo_dir)).evaluate(momenta)[0])
    assert abs(ufo_value - builtin_value) / max(
        abs(ufo_value),
        abs(builtin_value),
        1.0e-300,
    ) <= 1.0e-10


def _validation_momenta(process_dir: Path) -> np.ndarray:
    payload = json.loads(
        (process_dir / "validation_momenta.json").read_text(encoding="utf-8")
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


def test_ufo_sm_z_gluon_matches_builtin_and_refreshes_derived_parameters(
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime parity requires a release Rusticol extension")

    builtin_dir = tmp_path / "builtin"
    ufo_dir = tmp_path / "ufo"
    _generate_process(builtin_dir)
    _generate_process(
        ufo_dir,
        model=bundled_model_path("sm", "json"),
    )

    builtin_manifest = json.loads(
        (builtin_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    ufo_manifest = json.loads(
        (ufo_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    for field in (
        "amplitude_root_count",
        "current_count",
        "interaction_count",
        "source_count",
        "truncated",
    ):
        assert ufo_manifest["dag_summary"][field] == builtin_manifest["dag_summary"][
            field
        ]
    assert ufo_manifest["generation_filters"]["zero_current"]["skipped"] is False
    assert ufo_manifest["generation_filters"]["current_merging"]["skipped"] is False
    for stage in ufo_manifest["compiled"]["stage_evaluators"]["stages"]:
        assert stage["evaluator"]["settings"]["collect_factors"] is True

    momenta = _validation_momenta(ufo_dir)
    builtin_runtime = rusticol.Runtime.load(str(builtin_dir))
    ufo_runtime = rusticol.Runtime.load(str(ufo_dir))
    builtin_value = float(builtin_runtime.evaluate(momenta)[0])
    ufo_value = float(ufo_runtime.evaluate(momenta)[0])
    assert abs(ufo_value - builtin_value) / abs(builtin_value) <= 1.0e-10

    parameter_names = set(ufo_runtime.metadata()["model_parameter_names"])
    assert "aS" in parameter_names
    assert not any(name.startswith("derived_coupling_") for name in parameter_names)

    override = tmp_path / "half-alpha-s.json"
    override.write_text('{"aS": [0.059, 0.0]}\n', encoding="utf-8")
    overridden_runtime = rusticol.Runtime.load(
        str(ufo_dir),
        model_parameters=str(override),
    )
    overridden_value = float(overridden_runtime.evaluate(momenta)[0])
    assert overridden_value == pytest.approx(0.5 * ufo_value, rel=1.0e-12)


def test_ufo_sm_higgs_interference_matches_builtin(tmp_path: Path) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime parity requires a release Rusticol extension")

    builtin_dir = tmp_path / "builtin"
    ufo_dir = tmp_path / "ufo"
    options = {
        "process": "d d~ > z z z",
        "color_accuracy": "full",
    }
    _generate_process(builtin_dir, **options)
    _generate_process(
        ufo_dir,
        model=bundled_model_path("sm", "json"),
        **options,
    )
    matched_alpha_ew = 0.007546771114
    parameter_card = tmp_path / "matched-built-in-parameters.json"
    parameter_card.write_text(
        json.dumps(
            {
                "aEWM1": [1.0 / matched_alpha_ew, 0.0],
                "Gf": [
                    1.16639e-5 * matched_alpha_ew / (1.0 / 132.507),
                    0.0,
                ],
                "WH": [0.006382339, 0.0],
            }
        ),
        encoding="utf-8",
    )

    momenta = _validation_momenta(ufo_dir)
    builtin_value = float(
        rusticol.Runtime.load(str(builtin_dir)).evaluate(momenta)[0]
    )
    ufo_value = float(
        rusticol.Runtime.load(
            str(ufo_dir),
            model_parameters=str(parameter_card),
        ).evaluate(momenta)[0]
    )
    assert ufo_value == pytest.approx(builtin_value, rel=1.0e-10)


def test_ufo_sm_unitary_gauge_excludes_internal_goldstone_double_count(
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime parity requires a release Rusticol extension")

    builtin_dir = tmp_path / "builtin"
    ufo_dir = tmp_path / "ufo"
    options = {
        "process": "d d~ > t t~ z h",
        "color_accuracy": "full",
    }
    _generate_process(builtin_dir, **options)
    _generate_process(
        ufo_dir,
        model=bundled_model_path("sm", "json"),
        **options,
    )
    builtin_manifest = json.loads(
        (builtin_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    ufo_manifest = json.loads(
        (ufo_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    for field in ("current_count", "interaction_count", "amplitude_root_count"):
        assert ufo_manifest["dag_summary"][field] == builtin_manifest["dag_summary"][
            field
        ]
    generic_manifest = json.loads(
        (ufo_dir / "generic_process_manifest.json").read_text(encoding="utf-8")
    )
    assert all(
        250 not in interaction["vertex_particles"]
        for interaction in generic_manifest["interactions"]
    )

    matched_alpha_ew = 0.007546771114
    parameter_card = tmp_path / "matched-built-in-parameters.json"
    parameter_card.write_text(
        json.dumps(
            {
                "aEWM1": [1.0 / matched_alpha_ew, 0.0],
                "Gf": [
                    1.16639e-5 * matched_alpha_ew / (1.0 / 132.507),
                    0.0,
                ],
                "WH": [0.006382339, 0.0],
            }
        ),
        encoding="utf-8",
    )
    momenta = _validation_momenta(ufo_dir)
    builtin_value = float(
        rusticol.Runtime.load(str(builtin_dir)).evaluate(momenta)[0]
    )
    ufo_value = float(
        rusticol.Runtime.load(
            str(ufo_dir),
            model_parameters=str(parameter_card),
        ).evaluate(momenta)[0]
    )

    assert ufo_value == pytest.approx(builtin_value, rel=1.0e-10)


@pytest.mark.parametrize(
    ("process", "lc_sector_id"),
    (
        ("d d~ > t t~", 0),
        ("d d~ > t t~", 1),
        ("g g > t t~", 0),
        ("g g > t t~", 1),
        ("g g > g g", 0),
        ("g g > g g", 1),
        ("g g > g g", 2),
    ),
)
def test_ufo_sm_qcd_color_flows_match_builtin(
    tmp_path: Path,
    process: str,
    lc_sector_id: int,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime parity requires a release Rusticol extension")

    builtin_dir = tmp_path / "builtin"
    ufo_dir = tmp_path / "ufo"
    generation_options = {
        "process": process,
        "lc_sector_id": lc_sector_id,
        "max_coupling_orders": ("QED=0",),
    }
    _generate_process(builtin_dir, **generation_options)
    _generate_process(
        ufo_dir,
        model=bundled_model_path("sm", "json"),
        **generation_options,
    )

    builtin_manifest = json.loads(
        (builtin_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    ufo_manifest = json.loads(
        (ufo_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    if process != "g g > g g":
        for field in (
            "amplitude_root_count",
            "current_count",
            "interaction_count",
            "source_count",
            "truncated",
        ):
            assert ufo_manifest["dag_summary"][field] == builtin_manifest[
                "dag_summary"
            ][field]
    else:
        assert ufo_manifest["lowering_status"]["full_tensor_network_ready"] is True
        assert "ufo_contact_auxiliary_no_propagator" in {
            name
            for name, _count in ufo_manifest["lowering_status"][
                "required_propagator_kernel_counts"
            ]
        }

    momenta = _validation_momenta(ufo_dir)
    builtin_value = float(rusticol.Runtime.load(str(builtin_dir)).evaluate(momenta)[0])
    ufo_value = float(rusticol.Runtime.load(str(ufo_dir)).evaluate(momenta)[0])
    assert abs(ufo_value - builtin_value) / max(
        abs(ufo_value),
        abs(builtin_value),
        1.0e-300,
    ) <= 1.0e-10


def test_five_scalar_contact_tree_matches_direct_ufo_normalization(
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime contact validation requires a release Rusticol extension")

    output = tmp_path / "five-scalar-contact"
    _generate_process(
        output,
        model=bundled_model_path("scalars", "json"),
        process="scalar_0 scalar_0 > scalar_0 scalar_0 scalar_0",
        max_coupling_orders=("QCD=1",),
    )

    manifest = json.loads(
        (output / "process_manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["dag_summary"]["current_count"] == 18
    assert manifest["dag_summary"]["interaction_count"] == 18
    propagator_kernels = {
        name
        for name, _count in manifest["lowering_status"][
            "required_propagator_kernel_counts"
        ]
    }
    assert propagator_kernels == {"ufo_contact_auxiliary_no_propagator"}

    momenta = _validation_momenta(output)
    value = float(rusticol.Runtime.load(str(output)).evaluate(momenta)[0])
    assert value == pytest.approx(1.0 / math.factorial(3), rel=1.0e-12)


def test_custom_scalar_propagator_runs_through_generated_stage(
    tmp_path: Path,
) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime custom-propagator validation requires release Rusticol")

    payload = json.loads(
        bundled_model_path("scalars", "json").read_text(encoding="utf-8")
    )
    particle = next(item for item in payload["particles"] if item["name"] == "scalar_1")
    propagator = next(
        item for item in payload["propagators"] if item["particle"] == "scalar_1"
    )
    propagator["name"] += "_custom"
    particle["propagator"] = propagator["name"]
    custom_model = tmp_path / "custom-scalars.json"
    custom_model.write_text(json.dumps(payload), encoding="utf-8")

    reference_dir = tmp_path / "reference"
    custom_dir = tmp_path / "custom"
    options = {
        "process": "scalar_1 scalar_1 > scalar_1 scalar_1",
        "max_coupling_orders": ("QCD=2",),
    }
    _generate_process(
        reference_dir,
        model=bundled_model_path("scalars", "json"),
        **options,
    )
    _generate_process(
        custom_dir,
        model=custom_model,
        model_restriction="none",
        **options,
    )

    manifest = json.loads(
        (custom_dir / "process_manifest.json").read_text(encoding="utf-8")
    )
    assert "ufo_custom_propagator" in {
        name
        for name, _count in manifest["lowering_status"][
            "required_propagator_kernel_counts"
        ]
    }
    momenta = _validation_momenta(custom_dir)
    reference_value = float(
        rusticol.Runtime.load(str(reference_dir)).evaluate(momenta)[0]
    )
    custom_value = float(rusticol.Runtime.load(str(custom_dir)).evaluate(momenta)[0])
    assert custom_value == pytest.approx(reference_value, rel=1.0e-12)


def test_scalar_gravity_external_spin2_runs_at_high_precision(tmp_path: Path) -> None:
    rusticol = pytest.importorskip("rusticol")
    if rusticol.build_profile() != "release":
        pytest.skip("runtime spin-2 validation requires a release Rusticol extension")

    output = tmp_path / "scalar-gravity"
    _generate_process(
        output,
        model=bundled_model_path("scalar_gravity", "json"),
        process="scalar_0 scalar_0 > graviton graviton",
    )

    momenta = _validation_momenta(output)
    runtime = rusticol.Runtime.load(str(output))
    double_value = float(runtime.evaluate(momenta)[0])
    high_precision_value = float(runtime.evaluate_with_prec(momenta, 50)[0])
    assert high_precision_value == pytest.approx(double_value, rel=1.0e-12)
