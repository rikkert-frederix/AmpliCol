from __future__ import annotations

import json

import pytest

from pyamplicol.dag_runtime import (
    SymbolicaEvaluatorSettings,
    ZGluonDAGEvaluator,
    _resolve_compiled_preset,
)
from pyamplicol.native import LeadingColorZJetsNativeEvaluator


@pytest.mark.parametrize("gluon_count", [1, 2, 3, 4, 5, 6])
def test_z_gluon_dag_evaluator_matches_staged_recursion(
    gluon_count: int,
) -> None:
    process = "d d~ > z " + " ".join(["g"] * gluon_count)
    reference = LeadingColorZJetsNativeEvaluator()
    particles = reference.canonical_z_gluon_point(
        process,
        gluon_count=gluon_count,
        sqrt_s=1000.0,
    )
    expected = reference.evaluate(process, particles=particles)

    evaluator = ZGluonDAGEvaluator(process)
    actual = evaluator.evaluate(particles)

    metadata = evaluator.metadata
    assert metadata.kernel == "spenso-symbolica-shared-helicity-current-dag"
    assert metadata.gluon_count == gluon_count
    assert metadata.shared_amplitude_count == 2 * 3 * (2**gluon_count)
    assert metadata.shared_source_current_count == 2 * gluon_count + 7
    assert metadata.symbolica_output_count == metadata.shared_amplitude_count
    assert metadata.merge_evaluators_strategy is False
    assert metadata.symbolica_evaluator_settings["backend"] == "jit"
    assert _relative_difference(
        actual.matrix_element,
        expected.matrix_element,
    ) < 1.0e-11
    assert _relative_difference(
        actual.raw_helicity_sum,
        expected.raw_helicity_sum,
    ) < 1.0e-11
    assert _max_helicity_amplitude_difference(actual, expected) < 1.0e-12


def test_z_gluon_dag_marks_stage_momentum_parameters_real() -> None:
    evaluator = ZGluonDAGEvaluator(
        "d d~ > z g g",
        merge_evaluators_strategy=False,
    )
    momentum_indices: list[int] = []
    for offset in evaluator.compiled.layout.momentum_offsets.values():
        momentum_indices.extend(range(offset, offset + 4))
    assert set(momentum_indices) <= set(evaluator.compiled.layout.real_valued_inputs)


def test_z_gluon_dag_raw_sum_final_stage_is_opt_in() -> None:
    default_evaluator = ZGluonDAGEvaluator("d d~ > z g")
    assert (
        default_evaluator.metadata.symbolica_evaluator_settings[
            "raw_sum_final_stage"
        ]
        is False
    )

    raw_sum_evaluator = ZGluonDAGEvaluator(
        "d d~ > z g",
        symbolica_evaluator_backend="compiled-complex",
        symbolica_raw_sum_final_stage=True,
    )
    assert (
        raw_sum_evaluator.metadata.symbolica_evaluator_settings[
            "raw_sum_final_stage"
        ]
        is True
    )


def test_z_gluon_dag_compiled_balanced_preset_sets_generic_cpp_o1() -> None:
    evaluator = ZGluonDAGEvaluator(
        "d d~ > z g",
        symbolica_evaluator_backend="compiled-complex",
        symbolica_compiled_preset="balanced",
    )
    settings = evaluator.metadata.symbolica_evaluator_settings
    assert settings["compiled_preset"] == "balanced"
    assert settings["compiled_inline_asm"] == "none"
    assert settings["compiled_optimization_level"] == 1


def test_z_gluon_dag_compiled_adaptive_preset_switches_by_multiplicity() -> None:
    assert (
        _resolve_compiled_preset(
            "adaptive",
            gluon_count=5,
            inline_asm="default",
            optimization_level=3,
            output_chunk_size=None,
        )
        == ("none", 1, 64)
    )
    assert (
        _resolve_compiled_preset(
            "adaptive",
            gluon_count=6,
            inline_asm="none",
            optimization_level=1,
            output_chunk_size=None,
        )
        == ("none", 1, 64)
    )
    assert (
        _resolve_compiled_preset(
            "adaptive",
            gluon_count=7,
            inline_asm="none",
            optimization_level=1,
            output_chunk_size=32,
        )
        == ("none", 1, 32)
    )
    assert (
        _resolve_compiled_preset(
            "adaptive",
            gluon_count=8,
            inline_asm="none",
            optimization_level=1,
            output_chunk_size=None,
        )
        == ("default", 3, None)
    )


def test_z_gluon_dag_compiled_runtime_preset_uses_measured_fast_path() -> None:
    assert (
        _resolve_compiled_preset(
            "runtime",
            gluon_count=5,
            inline_asm="default",
            optimization_level=3,
            output_chunk_size=None,
        )
        == ("none", 2, 64)
    )
    assert (
        _resolve_compiled_preset(
            "runtime",
            gluon_count=6,
            inline_asm="default",
            optimization_level=3,
            output_chunk_size=None,
        )
        == ("none", 2, 64)
    )
    assert (
        _resolve_compiled_preset(
            "runtime",
            gluon_count=7,
            inline_asm="default",
            optimization_level=3,
            output_chunk_size=None,
        )
        == ("none", 1, 96)
    )
    assert (
        _resolve_compiled_preset(
            "runtime",
            gluon_count=8,
            inline_asm="default",
            optimization_level=3,
            output_chunk_size=None,
        )
        == ("none", 1, 96)
    )


def test_z_gluon_dag_compiled_runtime_o3_preset_forces_o3() -> None:
    assert (
        _resolve_compiled_preset(
            "runtime-o3",
            gluon_count=5,
            inline_asm="default",
            optimization_level=1,
            output_chunk_size=None,
        )
        == ("none", 3, 64)
    )
    assert (
        _resolve_compiled_preset(
            "runtime-o3",
            gluon_count=7,
            inline_asm="default",
            optimization_level=1,
            output_chunk_size=None,
        )
        == ("none", 3, 96)
    )
    assert (
        _resolve_compiled_preset(
            "runtime-o3",
            gluon_count=9,
            inline_asm="default",
            optimization_level=1,
            output_chunk_size=None,
        )
        == ("none", 3, 96)
    )


def test_symbolica_compiled_output_dir_is_reported_in_metadata() -> None:
    settings = SymbolicaEvaluatorSettings(
        backend="compiled-complex",
        compiled_output_dir="pyAmpliCol/outputs/unit-test",
    )

    assert (
        settings.to_json_dict()["compiled_output_dir"]
        == "pyAmpliCol/outputs/unit-test"
    )


def test_z_gluon_dag_compiled_evaluator_artifact_roundtrip(tmp_path) -> None:
    process = "d d~ > z g"
    reference = LeadingColorZJetsNativeEvaluator()
    particles = reference.canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    artifact_dir = tmp_path / "compiled-artifact"
    evaluator = ZGluonDAGEvaluator(
        process,
        symbolica_evaluator_backend="compiled-complex",
        symbolica_compiled_preset="runtime-o3",
        symbolica_compiled_output_dir=artifact_dir / "compiled",
        symbolica_n_cores=1,
    )
    expected = evaluator.evaluate(particles)

    manifest = evaluator.save_evaluator_artifact(artifact_dir)
    loaded = ZGluonDAGEvaluator(
        process,
        symbolica_load_evaluator_dir=artifact_dir,
    )
    actual = loaded.evaluate(particles)

    assert manifest.exists()
    assert (
        loaded.metadata.symbolica_evaluator_settings["loaded_from_artifact"]
        is True
    )
    assert _relative_difference(
        actual.matrix_element,
        expected.matrix_element,
    ) < 1.0e-12


def test_z_gluon_dag_jit_evaluator_artifact_roundtrip(tmp_path) -> None:
    process = "d d~ > z g"
    reference = LeadingColorZJetsNativeEvaluator()
    particles = reference.canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    artifact_dir = tmp_path / "jit-artifact"
    evaluator = ZGluonDAGEvaluator(
        process,
        symbolica_evaluator_backend="jit",
        symbolica_n_cores=1,
    )
    expected = evaluator.evaluate(particles)

    manifest_path = evaluator.save_evaluator_artifact(artifact_dir)
    manifest = json.loads(manifest_path.read_text())
    loaded = ZGluonDAGEvaluator(
        process,
        symbolica_load_evaluator_dir=artifact_dir,
    )
    actual = loaded.evaluate(particles)

    assert "jit-symbolica-evaluator" in _artifact_kinds(manifest["compiled"])
    assert (
        loaded.metadata.symbolica_evaluator_settings["loaded_from_artifact"]
        is True
    )
    assert _relative_difference(
        actual.matrix_element,
        expected.matrix_element,
    ) < 1.0e-12


def test_z_gluon_dag_split_vertex_current_stages_match_default() -> None:
    process = "d d~ > z g g"
    reference = LeadingColorZJetsNativeEvaluator()
    particles = reference.canonical_z_gluon_point(
        process,
        gluon_count=2,
        sqrt_s=1000.0,
    )
    default_evaluator = ZGluonDAGEvaluator(process)
    split_evaluator = ZGluonDAGEvaluator(
        process,
        split_vertex_current_stages=True,
    )

    default_result = default_evaluator.evaluate(particles)
    split_result = split_evaluator.evaluate(particles)

    assert split_evaluator.metadata.split_vertex_current_stages is True
    assert default_evaluator.metadata.split_vertex_current_stages is False
    assert _relative_difference(
        split_result.matrix_element,
        default_result.matrix_element,
    ) < 1.0e-11


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1e-300)


def _artifact_kinds(value: object) -> set[str]:
    kinds: set[str] = set()
    if isinstance(value, dict):
        kind = value.get("kind")
        if isinstance(kind, str):
            kinds.add(kind)
        for child in value.values():
            kinds.update(_artifact_kinds(child))
    elif isinstance(value, list):
        for child in value:
            kinds.update(_artifact_kinds(child))
    return kinds


def _max_helicity_amplitude_difference(left: object, right: object) -> float:
    left_contributions = getattr(left, "helicity_contributions")
    right_contributions = getattr(right, "helicity_contributions")
    assert len(left_contributions) == len(right_contributions)
    return max(
        abs(left_item.amplitude - right_item.amplitude)
        for left_item, right_item in zip(left_contributions, right_contributions)
    )
