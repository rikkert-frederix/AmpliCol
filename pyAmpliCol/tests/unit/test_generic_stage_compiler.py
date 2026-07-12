from __future__ import annotations

from types import SimpleNamespace

import pytest

from pyamplicol.generic_artifact import build_generic_process_manifest
from pyamplicol.generic_dag import ColorState, CurrentIndex, CurrentNode
from pyamplicol.generic_stage_compiler import (
    _interaction_contribution,
    _select_measured_chunk_candidate,
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
    def vertex_component_expression(self, _kind, left, right, **_kwargs):
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


def test_interaction_contribution_rejects_complex_current_stage_color_weight() -> None:
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

    with pytest.raises(ValueError, match="complex color weights"):
        _interaction_contribution(
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
    stage = SimpleNamespace(stage_index=1, stage_kind="current-combine")
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
