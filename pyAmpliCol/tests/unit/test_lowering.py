from __future__ import annotations

from symbolica.community.idenso import simplify_color

from pyamplicol.lowering import (
    _GraphTensorExpressionBuilder,
    _clean_symbolica_string,
    build_symbolic_lowering_report,
)
from pyamplicol.legacy_matrix import NativeMatrixElementGenerator
from pyamplicol.model import AmplicolSMLeadingColorModel


def test_symbolic_lowering_report_exercises_spenso_and_idenso_hooks() -> None:
    report = build_symbolic_lowering_report(AmplicolSMLeadingColorModel())

    tensor_probe = report.tensor_network_probe
    assert report.tensor_library == "TensorLibrary.hep_lib_atom"
    assert tensor_probe.engine == "spenso"
    assert tensor_probe.output_rank == 2
    assert tensor_probe.output_size == 16
    assert tensor_probe.nonzero_entries == 4
    assert tensor_probe.max_abs_entry == 1.5
    assert "two_gluon_to_tensor" in tensor_probe.expression
    assert "tensor_gluon_to_gluon" in tensor_probe.expression

    color_probe = report.color_algebra_probe
    assert color_probe.engine == "idenso"
    assert "f(" in color_probe.input_expression
    assert "\x1b" not in tensor_probe.expression
    assert "\x1b" not in color_probe.input_expression
    assert color_probe.simplified_expression == "-8*CA"
    assert report.full_me_tensor_network_ready is False


def test_symbolic_lowering_plan_is_derived_from_recursion_graph() -> None:
    result = NativeMatrixElementGenerator().generate(
        "d d~ > z g g g",
        write_cache_metadata=False,
    )

    assert result.graph is not None
    assert result.symbolic_lowering is not None
    plan = result.symbolic_lowering.recursion_plan
    vertex_lowering = result.symbolic_lowering.vertex_lowering
    blueprint = result.symbolic_lowering.tensor_network_blueprint
    assert plan is not None
    assert vertex_lowering is not None
    assert blueprint is not None
    assert plan.engine == "symbolica"
    assert plan.current_count == len(result.graph.currents)
    assert plan.interaction_count == len(result.graph.interactions)
    assert plan.amplitude_count == len(result.graph.amplitudes)
    assert plan.color_order == result.graph.color_order
    assert plan.tensor_route_vertex_kinds == (1, 2, 3)
    assert {kind for kind, _ in plan.vertex_kind_counts}.issuperset({1, 2, 3, 6, 10})
    assert "matrix_element_plan" in plan.expression
    assert plan.expression_length >= len(plan.expression)
    assert "assign(current" in plan.first_assignments[0]

    assert vertex_lowering.total_interactions == len(result.graph.interactions)
    assert vertex_lowering.full_tensor_network_ready is True
    assert vertex_lowering.tensor_names == (
        "current_momentum",
        "g",
        "gluon_tensor_to_gluon",
        "quark_vector_weyl_minus",
        "quark_vector_weyl_plus",
        "tensor_gluon_to_gluon",
        "two_gluon_to_tensor",
    )
    assert set(vertex_lowering.ready_vertex_kind_counts) == {
        (0, 4),
        (1, 2),
        (2, 1),
        (3, 1),
        (6, 24),
        (10, 8),
    }
    assert vertex_lowering.pending_vertex_kind_counts == ()
    assert vertex_lowering.first_steps[0].backend == "spenso"
    assert vertex_lowering.first_steps[0].full_tensor_network_ready is True

    assert blueprint.engine == "spenso"
    assert blueprint.current_count == len(result.graph.currents)
    assert blueprint.interaction_count == len(result.graph.interactions)
    assert blueprint.amplitude_count == len(result.graph.amplitudes)
    assert blueprint.expression_built is True
    assert blueprint.expression_executed is True
    assert blueprint.status == "scalar-skeleton"
    assert blueprint.full_me_tensor_network_ready is True
    assert blueprint.propagator_lowering_ready is True
    assert blueprint.ready_interactions == 40
    assert blueprint.pending_interactions == 0
    assert blueprint.placeholder_vertex_kinds == ()
    assert blueprint.parametric_external_current_count == 7
    assert blueprint.parametric_source_current_parameter_count == 22
    assert blueprint.parametric_current_momentum_count == 18
    assert blueprint.parametric_momentum_parameter_count == 72
    assert blueprint.parametric_parameter_count == 94
    assert blueprint.registered_tensor_names == (
        "current_momentum",
        "g",
        "gluon_tensor_to_gluon",
        "quark_vector_weyl_minus",
        "quark_vector_weyl_plus",
        "tensor_gluon_to_gluon",
        "two_gluon_to_tensor",
    )
    assert blueprint.expression is not None
    assert "quark_vector_weyl" in blueprint.expression
    assert "current_momentum" in blueprint.expression
    assert "propagator" in blueprint.expression
    assert "vertex_kind_0" not in blueprint.expression
    assert "vertex_kind_6" not in blueprint.expression
    assert blueprint.executed_expression is not None
    assert blueprint.execution_time_s is not None
    assert blueprint.execution_time_s >= 0.0


def test_one_gluon_blueprint_executes_with_parametric_source_currents() -> None:
    result = NativeMatrixElementGenerator().generate(
        "d d~ > z g",
        write_cache_metadata=False,
    )

    assert result.symbolic_lowering is not None
    blueprint = result.symbolic_lowering.tensor_network_blueprint
    assert blueprint is not None
    assert blueprint.status == "scalar-skeleton"
    assert blueprint.expression_built is True
    assert blueprint.expression_executed is True
    assert blueprint.full_me_tensor_network_ready is True
    assert blueprint.propagator_lowering_ready is True
    assert blueprint.placeholder_vertex_kinds == ()
    assert blueprint.parametric_external_current_count == 5
    assert blueprint.parametric_source_current_parameter_count == 14
    assert blueprint.parametric_current_momentum_count == 4
    assert blueprint.parametric_momentum_parameter_count == 16
    assert blueprint.parametric_parameter_count == 30
    assert blueprint.executed_expression is not None
    assert "propagator" not in blueprint.executed_expression
    assert "vertex_kind_6" not in blueprint.executed_expression
    assert "vertex_kind_10" not in blueprint.executed_expression


def test_photon_gluon_process_builds_native_recursion_graph_without_z_artifact() -> None:
    result = NativeMatrixElementGenerator().generate(
        "d d~ > a g",
        write_cache_metadata=False,
    )

    assert result.supported_native_target is True
    assert result.graph is not None
    assert result.graph.process == (1, -1, 21, 22)
    assert result.symbolic_lowering is not None
    assert "gamma" in " ".join(result.notes)
    assert result.backend == "native-python-recursion-staged"


def test_w_gluon_process_builds_flavour_changing_recursion_graph() -> None:
    result = NativeMatrixElementGenerator().generate(
        "u d~ > w+ g",
        write_cache_metadata=False,
    )

    assert result.supported_native_target is True
    assert result.graph is not None
    assert result.graph.process == (2, -1, 21, 24)
    assert result.backend == "native-python-recursion-staged"
    assert "W+" in " ".join(result.notes)
    charged_vertices = [
        interaction
        for interaction in result.graph.interactions
        if interaction.vertex_kind == 10
    ]
    assert charged_vertices
    assert {vertex.left.pdg for vertex in charged_vertices} == {1}
    assert {vertex.right.pdg for vertex in charged_vertices} == {24}
    assert {vertex.result.pdg for vertex in charged_vertices} == {2}
    assert {amplitude[0].pdg for amplitude in result.graph.amplitudes} == {2}
    assert {amplitude[1].pdg for amplitude in result.graph.amplitudes} == {-2}


def test_neutral_dilepton_gluon_process_reports_native_recursion_support() -> None:
    result = NativeMatrixElementGenerator().generate(
        "d d~ > e+ e- g",
        write_cache_metadata=False,
    )

    assert result.supported_native_target is True
    assert result.graph is None
    assert result.symbolic_lowering is None
    assert result.backend == "native-python-recursion-staged"
    assert "neutral dilepton" in " ".join(result.notes)


def test_zero_gluon_neutral_dilepton_process_reports_native_recursion_support() -> None:
    result = NativeMatrixElementGenerator().generate(
        "d d~ > e+ e-",
        write_cache_metadata=False,
    )

    assert result.supported_native_target is True
    assert result.graph is None
    assert result.symbolic_lowering is None
    assert result.backend == "native-python-recursion-staged"
    assert "neutral dilepton" in " ".join(result.notes)


def test_charged_leptonic_w_gluon_process_reports_native_recursion_support() -> None:
    result = NativeMatrixElementGenerator().generate(
        "u d~ > e+ ve g",
        write_cache_metadata=False,
    )

    assert result.supported_native_target is True
    assert result.graph is None
    assert result.symbolic_lowering is None
    assert result.backend == "native-python-recursion-staged"
    notes = " ".join(result.notes)
    assert "charged-current leptonic W" in notes


def test_zero_gluon_charged_leptonic_w_process_reports_native_recursion_support() -> None:
    result = NativeMatrixElementGenerator().generate(
        "u d~ > e+ ve",
        write_cache_metadata=False,
    )

    assert result.supported_native_target is True
    assert result.graph is None
    assert result.symbolic_lowering is None
    assert result.backend == "native-python-recursion-staged"
    notes = " ".join(result.notes)
    assert "charged-current leptonic W" in notes


def test_one_gluon_tensor_network_input_remains_factorized() -> None:
    model = AmplicolSMLeadingColorModel()
    result = NativeMatrixElementGenerator(model=model).generate(
        "d d~ > z g",
        write_cache_metadata=False,
    )

    assert result.graph is not None
    raw_expression = _GraphTensorExpressionBuilder(
        model,
        result.graph,
    ).matrix_element_skeleton()
    expression = _clean_symbolica_string(str(simplify_color(raw_expression)))

    assert expression.count(")*current_m1_1_p0(") == 2
    assert expression.count("quark_vector_weyl_plus(") == 4
    assert expression.count("quark_vector_weyl_minus(") == 4
    assert expression.count("quark_weyl_propagator_") == 4
    assert "current_p1_2_p1(" in expression
    assert "current_p1_2_m1(" in expression
    assert "pyamplicol::vertex_kind" not in expression
    assert "current_momentum" not in expression


def test_four_gluon_blueprint_builds_but_guards_execution() -> None:
    result = NativeMatrixElementGenerator().generate(
        "d d~ > z g g g g",
        write_cache_metadata=False,
    )

    assert result.symbolic_lowering is not None
    blueprint = result.symbolic_lowering.tensor_network_blueprint
    assert blueprint is not None
    assert blueprint.status == "execution-size-guarded"
    assert blueprint.expression_built is True
    assert blueprint.expression_executed is False
    assert blueprint.full_me_tensor_network_ready is False
    assert blueprint.propagator_lowering_ready is True
    assert blueprint.current_count == 38
    assert blueprint.interaction_count == 75
    assert blueprint.parametric_parameter_count == 130
    assert blueprint.expression is not None
    assert blueprint.executed_expression is None
