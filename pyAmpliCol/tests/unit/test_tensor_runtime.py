from __future__ import annotations

import pytest

pytestmark = pytest.mark.skip(
    reason=(
        "legacy Z+gluon tensor-network evaluators are retired from production; "
        "schema-v2 generic DAG/Rusticol tests cover the active path"
    )
)

from pyamplicol.native import LeadingColorZJetsNativeEvaluator
from pyamplicol.lowering import (
    _current_momentum_currents,
    _current_momentum_parameter_head,
)
from pyamplicol.tensor_runtime import (
    ZGluonNumericTensorNetworkEvaluator,
    ZGluonTensorNetworkEvaluator,
)


@pytest.mark.parametrize("gluon_count", [1, 2, 3])
def test_z_gluon_tensor_network_evaluator_matches_native(
    gluon_count: int,
) -> None:
    process = "d d~ > z " + " ".join(["g"] * gluon_count)
    native_evaluator = LeadingColorZJetsNativeEvaluator()
    particles = native_evaluator.canonical_z_gluon_point(
        process,
        gluon_count=gluon_count,
        sqrt_s=1000.0,
    )
    native = native_evaluator.evaluate(process, particles=particles)

    evaluator = ZGluonTensorNetworkEvaluator(process)
    tensor_result = evaluator.evaluate(particles)
    assert _relative_difference(
        tensor_result.matrix_element,
        native.matrix_element,
    ) < 1e-11
    assert _relative_difference(
        tensor_result.raw_helicity_sum,
        native.raw_helicity_sum,
    ) < 1e-11
    assert _max_helicity_amplitude_difference(tensor_result, native) < 1e-12


def test_z_gluon_tensor_network_evaluator_roundtrips_artifact_payload() -> None:
    process = "d d~ > z g g"
    native_evaluator = LeadingColorZJetsNativeEvaluator()
    particles = native_evaluator.canonical_z_gluon_point(
        process,
        gluon_count=2,
        sqrt_s=1000.0,
    )
    native = native_evaluator.evaluate(process, particles=particles)

    evaluator = ZGluonTensorNetworkEvaluator(process)
    tensor_result = evaluator.evaluate(particles)
    assert _relative_difference(
        tensor_result.matrix_element,
        native.matrix_element,
    ) < 1e-11

    payload = evaluator.to_artifact_payload()
    assert payload["kind"] == "pyamplicol-z-gluon-tensor-network-evaluator"
    metadata = payload["runtime_metadata"]
    assert isinstance(metadata, dict)
    assert metadata["kernel"] == "symbolica-z-gluon-tensor-network"
    assert metadata["gluon_count"] == 2
    assert metadata["parameter_count"] == 62

    restored = ZGluonTensorNetworkEvaluator.from_artifact_payload(process, payload)
    restored_result = restored.evaluate(particles)
    assert _relative_difference(
        restored_result.matrix_element,
        native.matrix_element,
    ) < 1e-11
    assert _max_helicity_amplitude_difference(restored_result, native) < 1e-12


def test_z_gluon_tensor_network_evaluator_payload_can_store_helicity_filter() -> None:
    process = "d d~ > z g"
    native_evaluator = LeadingColorZJetsNativeEvaluator()
    particles = native_evaluator.canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )

    evaluator = ZGluonTensorNetworkEvaluator(process, build_helicity_filter=True)
    payload = evaluator.to_artifact_payload()
    helicity_filter = payload["helicity_filter"]
    assert isinstance(helicity_filter, dict)
    assert helicity_filter["original_count"] == 12
    assert helicity_filter["kept_count"] <= 12

    restored = ZGluonTensorNetworkEvaluator.from_artifact_payload(process, payload)
    filtered = restored.evaluate(particles)
    full = ZGluonTensorNetworkEvaluator(process).evaluate(particles)
    assert _relative_difference(
        filtered.matrix_element,
        full.matrix_element,
    ) < 1e-11


def test_z_gluon_tensor_network_marks_momentum_parameters_real() -> None:
    evaluator = ZGluonTensorNetworkEvaluator("d d~ > z g g")
    momentum_indices: list[int] = []
    for current in _current_momentum_currents(evaluator.graph):
        start, stop = evaluator.bundle.param_builder.positions[
            _current_momentum_parameter_head(current)
        ]
        momentum_indices.extend(range(start, stop))

    assert set(momentum_indices) <= set(evaluator.bundle.param_builder.real_valued_inputs)
    assert not [
        parameter_range
        for parameter_range in evaluator.bundle.param_builder.ranges.values()
        if "coupl" in parameter_range.role.lower()
        or any("coupl" in part.lower() for part in parameter_range.head)
    ]


@pytest.mark.parametrize("gluon_count", [1, 2])
def test_z_gluon_numeric_tensor_network_evaluator_matches_native(
    gluon_count: int,
) -> None:
    process = "d d~ > z " + " ".join(["g"] * gluon_count)
    native_evaluator = LeadingColorZJetsNativeEvaluator()
    particles = native_evaluator.canonical_z_gluon_point(
        process,
        gluon_count=gluon_count,
        sqrt_s=1000.0,
    )
    native = native_evaluator.evaluate(process, particles=particles)

    evaluator = ZGluonNumericTensorNetworkEvaluator(process)
    numeric = evaluator.evaluate(particles)
    assert evaluator.metadata.kernel == "spenso-numeric-tensor-network"
    assert _relative_difference(
        numeric.matrix_element,
        native.matrix_element,
    ) < 1e-11
    assert _relative_difference(
        numeric.raw_helicity_sum,
        native.raw_helicity_sum,
    ) < 1e-11
    assert _max_helicity_amplitude_difference(numeric, native) < 1e-12


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1e-300)


def _max_helicity_amplitude_difference(left: object, right: object) -> float:
    left_contributions = getattr(left, "helicity_contributions")
    right_contributions = getattr(right, "helicity_contributions")
    assert len(left_contributions) == len(right_contributions)
    return max(
        abs(left_item.amplitude - right_item.amplitude)
        for left_item, right_item in zip(left_contributions, right_contributions)
    )
