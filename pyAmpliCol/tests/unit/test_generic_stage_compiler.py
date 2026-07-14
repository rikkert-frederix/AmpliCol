from __future__ import annotations

from types import SimpleNamespace

import pytest
from symbolica import S

from pyamplicol.generic_artifact import build_generic_process_manifest
from pyamplicol.generic_dag import (
    ColorState,
    CurrentIndex,
    CurrentNode,
    InteractionNode,
)
from pyamplicol.generic_stage_compiler import (
    _compact_interaction_contribution,
    _fanout_aware_current_order,
    _interaction_contribution,
    _select_measured_chunk_candidate,
    _specialize_stage_symbolica_functions,
    _stage_symbolica_settings,
    build_generic_stage_compiler_blueprint,
)
from pyamplicol.symbolica_evaluator import SymbolicaEvaluatorSettings


def _current(current_id: int, *, particle_id: int = 21, chirality: int = 0) -> CurrentNode:
    return CurrentNode(
        id=current_id,
        index=CurrentIndex(
            particle_id=particle_id,
            external_mask=1 << current_id,
            external_labels=(current_id + 1,),
            ordered_external_labels=(current_id + 1,),
            helicity_ancestry=1 << current_id,
            chirality=chirality,
            spin_state=1,
            flavour_flow=(particle_id,),
            charge_flow=0,
            color_state=ColorState(accuracy="lc", sector_id=0),
            momentum_mask=1 << current_id,
        ),
        dimension=2,
        is_source=current_id < 2,
    )


class _ToyModel:
    def __init__(self) -> None:
        self.evaluation_count = 0

    def vertex_component_expression(self, _kind, left, right, **_kwargs):
        self.evaluation_count += 1
        return (left[0] + 2 * right[0], left[1] - right[1])


class _EquivalentToyModel(_ToyModel):
    def vertex_component_expression(self, kind, left, right, **_kwargs):
        self.evaluation_count += 1
        if kind == 7:
            return (-(right[0] + 2 * left[0]), -(right[1] - left[1]))
        return (left[0] + 2 * right[0], left[1] - right[1])


def test_interaction_contribution_applies_current_stage_color_weight() -> None:
    dag = SimpleNamespace(currents=(_current(0), _current(1), _current(2)))
    interaction = {
        "left_value_slot": {"value_slot_id": 10},
        "right_value_slot": {"value_slot_id": 11},
        "momentum_slots": {"left": 20, "right": 21},
        "left_current_id": 0,
        "right_current_id": 1,
        "result_current_id": 2,
        "vertex_kind": 6,
        "coupling": [1.0, 0.0],
        "color_weight": [-1.0, 0.0],
    }

    contribution = _interaction_contribution(
        dag,
        _ToyModel(),
        interaction,
        value_components_by_slot_id={
            10: (3.0, 5.0),
            11: (7.0, 11.0),
        },
        momentum_components_by_slot_id={
            20: (0.0, 0.0, 0.0, 0.0),
            21: (0.0, 0.0, 0.0, 0.0),
        },
        model_parameter_symbols={},
    )

    assert contribution == (-17.0, 6.0)


def test_interaction_contribution_applies_complex_current_stage_color_weight() -> None:
    dag = SimpleNamespace(currents=(_current(0), _current(1), _current(2)))
    interaction = {
        "left_value_slot": {"value_slot_id": 10},
        "right_value_slot": {"value_slot_id": 11},
        "momentum_slots": {"left": 20, "right": 21},
        "left_current_id": 0,
        "right_current_id": 1,
        "result_current_id": 2,
        "vertex_kind": 6,
        "coupling": [1.0, 0.0],
        "color_weight": [1.0, 1.0],
    }

    contribution = _interaction_contribution(
        dag,
        _ToyModel(),
        interaction,
        value_components_by_slot_id={
            10: (3.0, 5.0),
            11: (7.0, 11.0),
        },
        momentum_components_by_slot_id={
            20: (0.0, 0.0, 0.0, 0.0),
            21: (0.0, 0.0, 0.0, 0.0),
        },
        model_parameter_symbols={},
    )

    assert contribution == pytest.approx((17.0 + 17.0j, -6.0 - 6.0j))


def test_compact_interaction_contribution_evaluates_signed_fanout_once() -> None:
    currents = (_current(0), _current(1), _current(2), _current(3))
    interactions = tuple(
        InteractionNode(
            id=index,
            vertex_kind=6,
            vertex_particles=(21, 21, 21),
            left_id=0,
            right_id=1,
            result_id=2 + index,
            coupling=(1.0, 0.0),
            color_weight=((1.0, 0.0), (-1.0, 0.0))[index],
            lowering_backend="toy",
            full_tensor_network_ready=True,
            evaluation_group_id=7,
            evaluation_factor=(1.0, 0.0),
        )
        for index in range(2)
    )
    dag = SimpleNamespace(currents=currents, interactions=interactions)
    model = _ToyModel()
    evaluation_cache = {}
    common = {
        "value_components_by_slot_id": {
            10: (3.0, 5.0),
            11: (7.0, 11.0),
        },
        "input_value_slot_by_current_id": {0: 10, 1: 11},
        "momentum_components_by_slot_id": {
            20: (0.0, 0.0, 0.0, 0.0),
            21: (0.0, 0.0, 0.0, 0.0),
        },
        "momentum_slot_by_mask": {1: 20, 2: 21},
        "model_parameter_symbols": {},
        "coupling_cache": {},
        "evaluation_cache": evaluation_cache,
    }

    positive = _compact_interaction_contribution(dag, model, 0, **common)
    negative = _compact_interaction_contribution(dag, model, 1, **common)

    assert positive == (17.0, -6.0)
    assert negative == (-17.0, 6.0)
    assert model.evaluation_count == 1
    assert len(evaluation_cache) == 1


def test_compact_interaction_contribution_reuses_verified_cross_kind_relation() -> None:
    currents = (_current(0), _current(1), _current(2), _current(3))
    interactions = (
        InteractionNode(
            id=0,
            vertex_kind=6,
            vertex_particles=(21, 21, 21),
            left_id=0,
            right_id=1,
            result_id=2,
            coupling=(1.0, 0.0),
            color_weight=(1.0, 0.0),
            lowering_backend="toy",
            full_tensor_network_ready=True,
            evaluation_group_id=7,
            evaluation_factor=(1.0, 0.0),
        ),
        InteractionNode(
            id=1,
            vertex_kind=7,
            vertex_particles=(21, 21, 21),
            left_id=1,
            right_id=0,
            result_id=3,
            coupling=(1.0, 0.0),
            color_weight=(1.0, 0.0),
            lowering_backend="toy",
            full_tensor_network_ready=True,
            evaluation_group_id=7,
            evaluation_factor=(-1.0, 0.0),
        ),
    )
    dag = SimpleNamespace(currents=currents, interactions=interactions)
    model = _EquivalentToyModel()
    common = {
        "value_components_by_slot_id": {
            10: (3.0, 5.0),
            11: (7.0, 11.0),
        },
        "input_value_slot_by_current_id": {0: 10, 1: 11},
        "momentum_components_by_slot_id": {
            20: (0.0, 0.0, 0.0, 0.0),
            21: (0.0, 0.0, 0.0, 0.0),
        },
        "momentum_slot_by_mask": {1: 20, 2: 21},
        "model_parameter_symbols": {},
        "coupling_cache": {},
        "evaluation_cache": {},
    }

    positive = _compact_interaction_contribution(dag, model, 0, **common)
    negative = _compact_interaction_contribution(dag, model, 1, **common)

    assert positive == (17.0, -6.0)
    assert negative == (-17.0, 6.0)
    assert model.evaluation_count == 1


def test_fanout_aware_current_order_reduces_cross_chunk_evaluations() -> None:
    order, before, after = _fanout_aware_current_order(
        (0, 1, 2, 3),
        output_size_by_current={current_id: 2 for current_id in range(4)},
        evaluation_groups_by_current={
            0: frozenset((1, 10)),
            1: frozenset((2, 20)),
            2: frozenset((1, 30)),
            3: frozenset((2, 40)),
        },
        chunk_size=4,
    )

    assert order in {(0, 2, 1, 3), (1, 3, 0, 2)}
    assert before == 8
    assert after == 6


def test_fanout_aware_current_order_scales_to_large_stages() -> None:
    current_count = 4096
    group_count = current_count // 2
    current_ids = tuple(range(current_count))

    order, before, after = _fanout_aware_current_order(
        current_ids,
        output_size_by_current={current_id: 1 for current_id in current_ids},
        evaluation_groups_by_current={
            current_id: frozenset((current_id % group_count,))
            for current_id in current_ids
        },
        chunk_size=128,
    )

    assert len(order) == current_count
    assert set(order) == set(current_ids)
    assert after < before


def test_stage_function_specialization_inlines_single_use_calls() -> None:
    function, left, right, formal_left, formal_right = S(
        "stage_function",
        "left",
        "right",
        "formal_left",
        "formal_right",
    )

    outputs, definitions = _specialize_stage_symbolica_functions(
        (function(left, right),),
        ((function, (formal_left, formal_right), formal_left**2 + formal_right),),
    )

    assert outputs == (left**2 + right,)
    assert definitions == ()


def test_stage_function_specialization_inlines_distinct_calls_of_same_head() -> None:
    repeated, unused, left, right, formal_left, formal_right = S(
        "repeated_stage_function",
        "unused_stage_function",
        "left",
        "right",
        "formal_left",
        "formal_right",
    )
    repeated_definition = (
        repeated,
        (formal_left, formal_right),
        formal_left**2 + formal_right,
    )
    unused_definition = (unused, (formal_left,), formal_left + 1)
    original_outputs = (repeated(left, right), repeated(right, left))

    outputs, definitions = _specialize_stage_symbolica_functions(
        original_outputs,
        (repeated_definition, unused_definition),
    )

    assert outputs == (left**2 + right, left + right**2)
    assert definitions == ()


def test_stage_function_specialization_inlines_identical_repeated_calls() -> None:
    function, left, right, formal_left, formal_right = S(
        "repeated_concrete_stage_function",
        "left",
        "right",
        "formal_left",
        "formal_right",
    )
    definition = (
        function,
        (formal_left, formal_right),
        formal_left**2 + formal_right,
    )
    original_outputs = (function(left, right), function(left, right))

    outputs, definitions = _specialize_stage_symbolica_functions(
        original_outputs,
        (definition,),
    )

    assert outputs == (left**2 + right, left**2 + right)
    assert definitions == ()


def test_stage_function_specialization_inlines_mixed_concrete_calls() -> None:
    function, left, right, formal_left, formal_right = S(
        "partially_repeated_stage_function",
        "left",
        "right",
        "formal_left",
        "formal_right",
    )
    definition = (
        function,
        (formal_left, formal_right),
        formal_left**2 + formal_right,
    )

    outputs, definitions = _specialize_stage_symbolica_functions(
        (
            function(left, right),
            function(left, right) + function(right, left),
        ),
        (definition,),
    )

    assert outputs == (
        left**2 + right,
        left + left**2 + right + right**2,
    )
    assert definitions == ()


def test_stage_function_specialization_inlines_nested_calls() -> None:
    outer, inner, value, formal_outer, formal_inner = S(
        "outer_stage_function",
        "inner_stage_function",
        "value",
        "formal_outer",
        "formal_inner",
    )

    outputs, definitions = _specialize_stage_symbolica_functions(
        (outer(value),),
        (
            (outer, (formal_outer,), inner(formal_outer) + 1),
            (inner, (formal_inner,), formal_inner**2),
        ),
    )

    assert outputs == (value**2 + 1,)
    assert definitions == ()


def test_stage_function_specialization_retains_large_kernel_abbreviation() -> None:
    function, value, formal_value = S(
        "large_stage_function",
        "value",
        "formal_value",
    )
    body = sum(formal_value**power for power in range(1, 400))
    definition = (function, (formal_value,), body)

    outputs, definitions = _specialize_stage_symbolica_functions(
        (function(value),),
        (definition,),
    )

    assert body.get_byte_size() > 1024
    assert outputs == (function(value),)
    assert definitions == (definition,)


def test_stage_blueprint_reports_each_current_and_amplitude_stage() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    events: list[tuple[str, int, int]] = []

    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        stage_local_parameter_layout=True,
        progress_callback=lambda label, index, total: events.append(
            (label, index, total)
        ),
    )

    assert len(events) == blueprint.stage_count
    assert events[-1] == (
        "amplitude stage",
        blueprint.stage_count,
        blueprint.stage_count,
    )
    assert events[:-1] == [
        ("current stage", index, blueprint.stage_count)
        for index in range(1, blueprint.stage_count)
    ]


def test_tapered_stage_chunking_uses_recursion_position() -> None:
    current_stages = tuple(
        SimpleNamespace(stage_index=index, stage_kind="current")
        for index in range(1, 9)
    )
    amplitude_stage = SimpleNamespace(stage_index=9, stage_kind="amplitude-roots")
    blueprint = SimpleNamespace(stages=current_stages)
    settings = SymbolicaEvaluatorSettings(
        compiled_output_chunk_size=128,
        output_chunk_strategy="tapered-stage",
    )

    sizes = [
        _stage_symbolica_settings(item, blueprint, settings).compiled_output_chunk_size
        for item in (*current_stages, amplitude_stage)
    ]

    assert sizes == [None, None, 256, 256, 256, 128, 128, 64, None]


def test_measured_chunk_selection_requires_minimum_gain() -> None:
    assert (
        _select_measured_chunk_candidate(
            {64: 9.0, 128: 10.0, 256: 9.6, None: 11.0},
            baseline_size=128,
            minimum_gain=0.05,
        )
        == 64
    )
    assert (
        _select_measured_chunk_candidate(
            {64: 9.6, 128: 10.0, 256: 9.8, None: 11.0},
            baseline_size=128,
            minimum_gain=0.05,
        )
        == 128
    )


def test_auto_chunk_strategy_tunes_normal_jit_but_not_large_all_flow_chunks() -> None:
    stage = SimpleNamespace(
        stage_index=1,
        stage_kind="current-combine",
        output_length=256,
    )
    blueprint = SimpleNamespace(stages=(stage,))

    selected = _stage_symbolica_settings(
        stage,
        blueprint,
        SymbolicaEvaluatorSettings(
            backend="jit",
            compiled_output_chunk_size=128,
            output_chunk_strategy="auto",
        ),
    )
    all_flow = _stage_symbolica_settings(
        stage,
        blueprint,
        SymbolicaEvaluatorSettings(
            backend="jit",
            compiled_output_chunk_size=8192,
            output_chunk_strategy="auto",
        ),
    )

    assert selected.output_chunk_strategy == "measured-stage"
    assert all_flow.output_chunk_strategy == "uniform"


def test_auto_chunk_strategy_skips_tuning_when_stage_already_fits() -> None:
    stage = SimpleNamespace(
        stage_index=1,
        stage_kind="current-combine",
        output_length=128,
    )

    selected = _stage_symbolica_settings(
        stage,
        SimpleNamespace(stages=(stage,)),
        SymbolicaEvaluatorSettings(
            backend="jit",
            compiled_output_chunk_size=128,
            output_chunk_strategy="auto",
        ),
    )

    assert selected.output_chunk_strategy == "uniform"
