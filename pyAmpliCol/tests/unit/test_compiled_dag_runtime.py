from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from pyamplicol.compiled_dag_runtime import (
    ZGluonCompiledDAGEvaluator,
    _compiled_dag_artifact_fingerprint,
    _external_momentum_array,
    _fill_compiled_dag_source_parameters,
    _fill_compiled_dag_source_parameters_batch,
    _resolve_compiled_dag_compiled_preset,
)
from pyamplicol.evaluation import NativeRuntimeEvaluator
from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.native import LeadingColorZJetsNativeEvaluator, NativeEvaluationError
from pyamplicol.phase_space import rambo_z_gluon_point


def test_symbolica_evaluator_alias_hook_is_available() -> None:
    from symbolica import Expression, S

    x = S("x")
    y = S("y")
    alias = S("pyamplicol_alias_test")
    evaluator = Expression.evaluator_multiple(
        [alias + 1],
        [x, y],
        aliases=[(alias, x * y + 2)],
        jit_compile=False,
    )

    assert evaluator.evaluate([[3.0, 5.0]]) == [[18.0]]


def test_symbolica_expression_alias_helper_feeds_nested_evaluator_aliases() -> None:
    from symbolica import Expression, S

    x = S("x")
    y = S("y")
    a = S("pyamplicol_alias_helper_a")
    b = S("pyamplicol_alias_helper_b")
    inner = (x + 1).alias(a)
    outer = (inner[0] * inner[0] + y).alias(b)
    evaluator = Expression.evaluator_multiple(
        [outer[0] + inner[0]],
        [x, y],
        aliases=[inner, outer],
        jit_compile=False,
    )

    assert hasattr(Expression, "alias")
    assert evaluator.evaluate([[2.0, 5.0]]) == [[17.0]]


def test_z_gluon_compiled_dag_matches_staged_recursion() -> None:
    process = "d d~ > z g"
    reference = LeadingColorZJetsNativeEvaluator()
    particles = reference.canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    expected = reference.evaluate(process, particles=particles)

    evaluator = ZGluonCompiledDAGEvaluator(
        process,
        lowering="symbolic",
        symbolica_evaluator_backend="jit",
        compiled_dag_helicity_filter=False,
    )
    actual = evaluator.evaluate(particles)

    assert evaluator.metadata.kernel == "symbolica-compiled-shared-current-alias-dag"
    assert evaluator.metadata.symbolica_alias_available is True
    assert evaluator.metadata.opaque_alias_component_count == evaluator.metadata.alias_component_count
    assert evaluator.metadata.symbolica_output_count == 12
    assert evaluator.metadata.symbolica_parameter_count < 6 * len(evaluator.table.currents)
    assert _relative_difference(actual.matrix_element, expected.matrix_element) < 1.0e-12
    assert _relative_difference(actual.raw_helicity_sum, expected.raw_helicity_sum) < 1.0e-12


def test_z_gluon_compiled_dag_spenso_lowering_matches_symbolic() -> None:
    reference = LeadingColorZJetsNativeEvaluator()
    for gluon_count, process in ((1, "d d~ > z g"), (2, "d d~ > z g g")):
        particles = reference.canonical_z_gluon_point(
            process,
            gluon_count=gluon_count,
            sqrt_s=1000.0,
        )
        symbolic = ZGluonCompiledDAGEvaluator(
            process,
            lowering="symbolic",
            compiled_dag_helicity_filter=False,
        )
        spenso = ZGluonCompiledDAGEvaluator(
            process,
            lowering="spenso",
            compiled_dag_helicity_filter=False,
        )

        symbolic_value = symbolic.evaluate(particles)
        spenso_value = spenso.evaluate(particles)

        assert spenso.metadata.effective_lowering == "spenso"
        assert spenso.metadata.lowering_note is None
        assert spenso.metadata.spenso_body_count > 0
        assert spenso.metadata.spenso_body_execution_time_s > 0.0
        assert (
            _relative_difference(
                spenso_value.matrix_element,
                symbolic_value.matrix_element,
            )
            < 1.0e-12
        )


def test_z_gluon_compiled_dag_cross_checks_spenso_against_symbolic() -> None:
    process = "d d~ > z g"
    particles = LeadingColorZJetsNativeEvaluator().canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    evaluator = ZGluonCompiledDAGEvaluator(
        process,
        lowering="spenso",
        cross_check_lowering=True,
        compiled_dag_helicity_filter=False,
    )

    evaluator.evaluate(particles)

    assert evaluator.metadata.cross_check_lowering is True
    assert evaluator.metadata.cross_check_max_relative_difference is not None
    assert evaluator.metadata.cross_check_max_relative_difference < 1.0e-12


def test_z_gluon_compiled_dag_marks_momenta_real() -> None:
    evaluator = ZGluonCompiledDAGEvaluator(
        "d d~ > z g",
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
    )
    momentum_indices: list[int] = []
    for offset in evaluator.compiled.layout.momentum_offsets.values():
        momentum_indices.extend(range(offset, offset + 4))

    assert set(momentum_indices) <= set(evaluator.compiled.layout.real_valued_inputs)


def test_z_gluon_compiled_dag_uses_adaptive_cpe_defaults() -> None:
    low = ZGluonCompiledDAGEvaluator(
        "d d~ > z g",
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
    )
    high = ZGluonCompiledDAGEvaluator(
        "d d~ > z g g g g g g",
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
    )

    assert low.metadata.symbolica_evaluator_settings["cpe_iterations"] == 2
    assert high.metadata.symbolica_evaluator_settings["cpe_iterations"] is None
    assert low.metadata.symbolica_evaluator_settings[
        "max_common_pair_distance"
    ] == 250
    assert high.metadata.symbolica_evaluator_settings[
        "max_common_pair_distance"
    ] == 250


def test_z_gluon_compiled_dag_vectorized_source_fill_matches_scalar() -> None:
    process = "d d~ > z g g"
    model = AmplicolSMLeadingColorModel()
    reference = LeadingColorZJetsNativeEvaluator(model)
    points = (
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1000.0),
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1325.0),
    )
    evaluator = ZGluonCompiledDAGEvaluator(
        process,
        model=model,
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
        compiled_dag_inline_external_wavefunctions=False,
    )

    vectorized = np.zeros(
        (len(points), evaluator.compiled.layout.parameter_count),
        dtype=np.complex128,
    )
    _fill_compiled_dag_source_parameters_batch(
        vectorized,
        evaluator.compiled._source_parameter_specs,
        _external_momentum_array(points, evaluator.compiled._external_label_count),
        model,
        gluon_count=2,
    )
    scalar = np.zeros_like(vectorized)
    for index, point in enumerate(points):
        _fill_compiled_dag_source_parameters(
            scalar[index],
            evaluator.compiled._source_parameter_specs,
            point,
            model,
            gluon_count=2,
        )

    assert np.allclose(vectorized, scalar, rtol=0.0, atol=1.0e-15)


def test_z_gluon_compiled_dag_inline_sources_match_source_parameters() -> None:
    process = "d d~ > z g g"
    model = AmplicolSMLeadingColorModel()
    reference = LeadingColorZJetsNativeEvaluator(model)
    points = (
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1000.0),
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1325.0),
    )
    inline = ZGluonCompiledDAGEvaluator(
        process,
        model=model,
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
    )
    source_parameter = ZGluonCompiledDAGEvaluator(
        process,
        model=model,
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
        compiled_dag_inline_external_wavefunctions=False,
    )

    inline_values = inline.evaluate_matrix_elements_array(points)
    source_values = source_parameter.evaluate_matrix_elements_array(points)

    assert np.allclose(inline_values, source_values, rtol=1.0e-12, atol=0.0)


def test_z_gluon_compiled_dag_real_arithmetic_matches_complex_path() -> None:
    process = "d d~ > z g g"
    reference = LeadingColorZJetsNativeEvaluator()
    points = (
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1000.0),
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1325.0),
    )
    complex_path = ZGluonCompiledDAGEvaluator(
        process,
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
        compiled_dag_real_arithmetic=False,
    )
    real_path = ZGluonCompiledDAGEvaluator(
        process,
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
        compiled_dag_real_arithmetic=True,
    )

    assert real_path.metadata.real_arithmetic is True
    assert real_path.compiled.evaluator_output_length == 2 * complex_path.compiled.output_length
    assert np.allclose(
        real_path.evaluate_matrix_elements_array(points),
        complex_path.evaluate_matrix_elements_array(points),
        rtol=1.0e-12,
        atol=0.0,
    )


def test_z_gluon_compiled_dag_matrix_element_array_matches_tuple_api() -> None:
    process = "d d~ > z g g"
    reference = LeadingColorZJetsNativeEvaluator()
    points = (
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1000.0),
        reference.canonical_z_gluon_point(process, gluon_count=2, sqrt_s=1375.0),
    )
    evaluator = ZGluonCompiledDAGEvaluator(
        process,
        lowering="symbolic",
        batch_size=16,
        compiled_dag_helicity_filter=False,
    )

    array_values = evaluator.evaluate_matrix_elements_array(points)
    tuple_values = evaluator.evaluate_matrix_elements_many(points)

    assert isinstance(array_values, np.ndarray)
    assert array_values.dtype == np.float64
    assert np.allclose(array_values, tuple_values, rtol=1.0e-14, atol=0.0)


def test_z_gluon_compiled_dag_save_load_roundtrip(tmp_path: Path) -> None:
    process = "d d~ > z g"
    generated = NativeRuntimeEvaluator(
        process,
        runtime_backend="compiled-dag",
        compiled_dag_lowering="symbolic",
        compiled_dag_helicity_filter_phase_space="canonical",
    )
    manifest_path = generated.save_evaluator_artifact(tmp_path)
    manifest = json.loads(manifest_path.read_text())

    loaded = NativeRuntimeEvaluator(
        process,
        runtime_backend="compiled-dag",
        compiled_dag_lowering="symbolic",
        symbolica_load_evaluator_dir=str(tmp_path),
    )
    point = LeadingColorZJetsNativeEvaluator().canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    before = generated.evaluate(particles=point)
    after = loaded.evaluate(particles=point)

    assert manifest["kind"] == "pyamplicol-compiled-dag-evaluator"
    assert manifest["symbolica"]["available"] is True
    assert manifest["symbolica"]["evaluator_aliases_hook"] is True
    assert "symbolica_commit" in manifest["symbolica"]
    assert manifest["metadata"]["symbolica_evaluator_settings"]["backend"] == "jit"
    assert manifest["compiled"]["evaluator"]["kind"] == "jit-symbolica-evaluator"
    state_path = tmp_path / manifest["compiled"]["evaluator"]["evaluator_state_path"]
    assert state_path.exists()
    assert state_path.stat().st_size > 0
    assert isinstance(manifest["current_metadata"], list)
    assert len(manifest["current_metadata"]) == manifest["metadata"]["shared_current_count"]
    assert isinstance(manifest["root_metadata"], list)
    assert len(manifest["root_metadata"]) == manifest["metadata"]["symbolica_output_count"]
    assert manifest["root_metadata"][0]["output_index"] == 0
    assert manifest["artifact_fingerprint"] == _compiled_dag_artifact_fingerprint(
        manifest
    )
    assert len(manifest["artifact_fingerprint"]) == 64
    assert manifest["metadata"]["helicity_filter"]["original_count"] == 12
    assert manifest["metadata"]["helicity_filter"]["kept_count"] <= 12
    assert manifest["metadata"]["helicity_filter"]["phase_space_mode"] == "canonical"
    assert manifest["compiled"]["kind"] == "compiled-dag-multi-output"
    assert _relative_difference(after.matrix_element, before.matrix_element) < 1.0e-14


def test_z_gluon_compiled_dag_chunked_artifact_roundtrip(tmp_path: Path) -> None:
    process = "d d~ > z g"
    generated = NativeRuntimeEvaluator(
        process,
        runtime_backend="compiled-dag",
        compiled_dag_lowering="symbolic",
        compiled_dag_helicity_filter=False,
        symbolica_compiled_output_chunk_size=5,
    )
    manifest_path = generated.save_evaluator_artifact(tmp_path)
    manifest = json.loads(manifest_path.read_text())
    point = LeadingColorZJetsNativeEvaluator().canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    loaded = NativeRuntimeEvaluator(
        process,
        runtime_backend="compiled-dag",
        compiled_dag_lowering="symbolic",
        symbolica_load_evaluator_dir=str(tmp_path),
    )
    before = generated.evaluate(particles=point)
    after = loaded.evaluate(particles=point)

    assert manifest["compiled"]["evaluator"]["kind"] == "chunked-symbolica-evaluator"
    assert len(manifest["compiled"]["evaluator"]["chunks"]) == 3
    assert _relative_difference(after.matrix_element, before.matrix_element) < 1.0e-14


def test_z_gluon_compiled_dag_raw_sum_final_stage_reports_backend_limit() -> None:
    try:
        ZGluonCompiledDAGEvaluator(
            "d d~ > z g",
            lowering="symbolic",
            compiled_dag_helicity_filter=False,
            symbolica_raw_sum_final_stage=True,
        )
    except NativeEvaluationError as exc:
        message = str(exc)
    else:
        raise AssertionError("compiled DAG unexpectedly accepted raw-sum final stage")

    assert "conj()" in message


def test_z_gluon_compiled_dag_default_filter_preserves_matrix_element() -> None:
    process = "d d~ > z g"
    point = LeadingColorZJetsNativeEvaluator().canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    unfiltered = ZGluonCompiledDAGEvaluator(
        process,
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
    )
    filtered = ZGluonCompiledDAGEvaluator(
        process,
        lowering="symbolic",
        compiled_dag_helicity_filter_phase_space="canonical",
    )

    unfiltered_value = unfiltered.evaluate(point)
    filtered_value = filtered.evaluate(point)

    assert filtered.metadata.helicity_filter is not None
    assert filtered.metadata.helicity_filter["original_count"] == len(unfiltered.table.amplitudes)
    assert filtered.metadata.symbolica_output_count == len(filtered.table.amplitudes)
    assert filtered.metadata.symbolica_output_count <= unfiltered.metadata.symbolica_output_count
    assert _relative_difference(
        filtered_value.matrix_element,
        unfiltered_value.matrix_element,
    ) < 1.0e-12


def test_rambo_z_gluon_point_conserves_momentum_and_masses() -> None:
    model = LeadingColorZJetsNativeEvaluator().model
    point = rambo_z_gluon_point(
        "d d~ > z g g",
        model,
        gluon_count=2,
        sqrt_s=1000.0,
        seed=17,
    )
    final_total = [0.0, 0.0, 0.0, 0.0]
    for particle in point[2:]:
        for index, component in enumerate(particle.momentum):
            final_total[index] += component

    assert abs(final_total[0] - 1000.0) < 1.0e-10
    assert max(abs(component) for component in final_total[1:]) < 1.0e-10
    for particle in point[2:-1]:
        energy, px, py, pz = particle.momentum
        assert abs(energy * energy - px * px - py * py - pz * pz) < 1.0e-8
    z = point[-1]
    energy, px, py, pz = z.momentum
    assert abs(
        energy * energy - px * px - py * py - pz * pz - model.mass(23) ** 2
    ) < 1.0e-8


def test_native_runtime_compiled_dag_backend_routes() -> None:
    evaluator = NativeRuntimeEvaluator(
        "d d~ > z g",
        runtime_backend="compiled-dag",
        compiled_dag_lowering="symbolic",
        compiled_dag_helicity_filter_phase_space="canonical",
    )

    assert evaluator.metadata.backend == "native-symbolica-compiled-shared-current-alias-dag"
    assert evaluator.metadata.evaluator_metadata is not None
    assert evaluator.metadata.evaluator_metadata["effective_lowering"] == "symbolic"
    assert evaluator.metadata.evaluator_metadata["helicity_filter"]["sample_count"] == 10


def test_native_runtime_compiled_dag_default_uses_spenso_lowering() -> None:
    evaluator = NativeRuntimeEvaluator(
        "d d~ > z g",
        runtime_backend="compiled-dag",
        compiled_dag_helicity_filter_phase_space="canonical",
    )

    assert evaluator.metadata.evaluator_metadata is not None
    assert evaluator.metadata.evaluator_metadata["lowering"] == "spenso"
    assert evaluator.metadata.evaluator_metadata["effective_lowering"] == "spenso"
    assert evaluator.metadata.evaluator_metadata["spenso_body_count"] > 0


def test_compiled_dag_rejects_unsupported_processes_with_clear_diagnostic() -> None:
    for process in (
        "d d~ > e+ e- g g",
        "u d~ > w+ g g",
        "d d~ > z z g",
    ):
        try:
            ZGluonCompiledDAGEvaluator(process)
        except NativeEvaluationError as exc:
            message = str(exc)
        else:
            raise AssertionError(
                f"compiled DAG unexpectedly accepted unsupported process {process!r}"
            )

        assert (
            "q q~ -> Z plus ordered gluons" in message
            or "no native graph" in message
            or "unsupported" in message.lower()
        )


def test_compiled_dag_presets_do_not_auto_chunk_outputs() -> None:
    assert _resolve_compiled_dag_compiled_preset(
        "runtime-o3",
        gluon_count=5,
        inline_asm="default",
        optimization_level=3,
        output_chunk_size=None,
    ) == ("none", 3, None)
    assert _resolve_compiled_dag_compiled_preset(
        "runtime-o3",
        gluon_count=5,
        inline_asm="default",
        optimization_level=3,
        output_chunk_size=64,
    ) == ("none", 3, 64)


def test_compiled_dag_jit_direct_translation_default_is_multiplicity_adaptive() -> None:
    low = ZGluonCompiledDAGEvaluator(
        "d d~ > z g",
        compiled_dag_helicity_filter=False,
    )
    high = ZGluonCompiledDAGEvaluator(
        "d d~ > z g g g g",
        compiled_dag_helicity_filter=False,
    )

    assert low.metadata.symbolica_evaluator_settings[
        "jit_direct_translation"
    ] is False
    assert high.metadata.symbolica_evaluator_settings[
        "jit_direct_translation"
    ] is True


def test_compiled_dag_inline_external_wavefunctions_default_and_opt_out() -> None:
    inline = ZGluonCompiledDAGEvaluator(
        "d d~ > z g",
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
    )
    source_parameter = ZGluonCompiledDAGEvaluator(
        "d d~ > z g",
        lowering="symbolic",
        compiled_dag_helicity_filter=False,
        compiled_dag_inline_external_wavefunctions=False,
    )

    assert inline.metadata.inline_external_wavefunctions is True
    assert source_parameter.metadata.inline_external_wavefunctions is False
    assert inline.compiled.parameter_count < source_parameter.compiled.parameter_count


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)
