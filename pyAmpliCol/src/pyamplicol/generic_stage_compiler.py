from __future__ import annotations

import heapq
import time
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .generic_artifact import (
    GenericProcessManifest,
    LC_SECTOR_SELECTOR_PARAMETER,
    _generic_runtime_schema_payload,
    _runtime_coupling_parameter_names,
)
from .generic_dag import GenericDAG
from .model import AmplicolSMLeadingColorModel, Model
from .params import ParamBuilder

_EXPRESSION_PREVIEW_LIMIT = 512
_MODEL_FUNCTION_INLINE_MAX_BYTES = 1024


@dataclass(frozen=True)
class GenericStageInputComponent:
    kind: str
    source_id: int
    component: int
    global_component: int
    parameter_index: int
    real_valued: bool = False

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "source_id": self.source_id,
            "component": self.component,
            "global_component": self.global_component,
            "parameter_index": self.parameter_index,
            "real_valued": self.real_valued,
        }


@dataclass(frozen=True)
class GenericStageOutputSlot:
    value_slot_id: int
    current_id: int
    variant: str
    component_start: int
    component_stop: int
    output_start: int
    output_stop: int

    def to_json_dict(self) -> dict[str, object]:
        return {
            "value_slot_id": self.value_slot_id,
            "current_id": self.current_id,
            "variant": self.variant,
            "component_start": self.component_start,
            "component_stop": self.component_stop,
            "output_start": self.output_start,
            "output_stop": self.output_stop,
        }


@dataclass(frozen=True)
class GenericCompiledStageBlueprint:
    stage_index: int
    stage_kind: str
    subset_size: int | None
    evaluator_label: str
    parameter_layout: str
    output_length: int
    output_slots: tuple[GenericStageOutputSlot, ...]
    input_value_slot_ids: tuple[int, ...]
    output_value_slot_ids: tuple[int, ...]
    interaction_ids: tuple[int, ...]
    input_components: tuple[GenericStageInputComponent, ...]
    parameter_count: int
    value_parameter_count: int
    momentum_parameter_count: int
    model_parameter_count: int
    real_valued_inputs: tuple[int, ...]
    expression_ready: bool
    blockers: tuple[str, ...]
    first_output_previews: tuple[str, ...]
    evaluation_groups_by_current: tuple[tuple[int, tuple[int, ...]], ...] = field(
        default=(),
        repr=False,
        compare=False,
    )
    fanout_chunk_size: int | None = None
    fanout_evaluation_occurrences_before: int | None = None
    fanout_evaluation_occurrences_after: int | None = None
    parameter_symbols: tuple[Any, ...] = field(
        default=(),
        repr=False,
        compare=False,
    )
    output_expressions: tuple[Any, ...] = field(
        default=(),
        repr=False,
        compare=False,
    )
    symbolica_functions: tuple[
        tuple[Any, tuple[Any, ...], Any], ...
    ] = field(
        default=(),
        repr=False,
        compare=False,
    )

    def to_json_dict(self) -> dict[str, object]:
        return {
            "stage_index": self.stage_index,
            "stage_kind": self.stage_kind,
            "subset_size": self.subset_size,
            "evaluator_label": self.evaluator_label,
            "parameter_layout": self.parameter_layout,
            "output_length": self.output_length,
            "output_slots": [slot.to_json_dict() for slot in self.output_slots],
            "input_value_slot_ids": list(self.input_value_slot_ids),
            "output_value_slot_ids": list(self.output_value_slot_ids),
            "interaction_ids": list(self.interaction_ids),
            "input_components": [
                component.to_json_dict() for component in self.input_components
            ],
            "parameter_count": self.parameter_count,
            "value_parameter_count": self.value_parameter_count,
            "momentum_parameter_count": self.momentum_parameter_count,
            "model_parameter_count": self.model_parameter_count,
            "real_valued_inputs": list(self.real_valued_inputs),
            "expression_ready": self.expression_ready,
            "blockers": list(self.blockers),
            "first_output_previews": list(self.first_output_previews),
            "symbolica_function_count": len(self.symbolica_functions),
            "evaluation_fanout": {
                "strategy": (
                    "shared-evaluation-affinity"
                    if self.fanout_chunk_size is not None
                    else "natural-current-order"
                ),
                "chunk_size": self.fanout_chunk_size,
                "evaluation_occurrences_before": (
                    self.fanout_evaluation_occurrences_before
                ),
                "evaluation_occurrences_after": (
                    self.fanout_evaluation_occurrences_after
                ),
            },
        }


@dataclass(frozen=True)
class GenericStageCompilerBlueprint:
    kind: str
    runtime_available: bool
    parameter_count: int
    value_parameter_count: int
    momentum_parameter_count: int
    model_parameter_count: int
    real_valued_inputs: tuple[int, ...]
    stage_count: int
    stages: tuple[GenericCompiledStageBlueprint, ...]
    amplitude_stage: GenericCompiledStageBlueprint
    expression_ready: bool
    blockers: tuple[str, ...]
    parameter_symbols: tuple[Any, ...] = field(
        default=(),
        repr=False,
        compare=False,
    )

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "runtime_available": self.runtime_available,
            "parameter_count": self.parameter_count,
            "value_parameter_count": self.value_parameter_count,
            "momentum_parameter_count": self.momentum_parameter_count,
            "model_parameter_count": self.model_parameter_count,
            "real_valued_inputs": list(self.real_valued_inputs),
            "stage_count": self.stage_count,
            "stages": [stage.to_json_dict() for stage in self.stages],
            "amplitude_stage": self.amplitude_stage.to_json_dict(),
            "expression_ready": self.expression_ready,
            "blockers": list(self.blockers),
        }


@dataclass(frozen=True)
class _StageLocalInputs:
    parameter_symbols: tuple[Any, ...]
    input_components: tuple[GenericStageInputComponent, ...]
    value_symbols: Mapping[int, tuple[Any, ...]]
    momentum_symbols: Mapping[int, tuple[Any, ...]]
    model_parameter_symbols: Mapping[str, Any]
    value_parameter_count: int
    momentum_parameter_count: int
    model_parameter_count: int
    real_valued_inputs: tuple[int, ...]


StageEvaluatorCompiler = Callable[
    [GenericCompiledStageBlueprint, tuple[Any, ...], tuple[int, ...]],
    dict[str, object],
]
StageBlueprintProgress = Callable[[str, int, int], None]
StageBlueprintConsumer = Callable[
    [GenericCompiledStageBlueprint, int, int],
    None,
]


class _RuntimeParameterizedModel(AmplicolSMLeadingColorModel):
    """Use runtime symbols for numeric model constants without changing topology."""

    def __init__(self, base: Model, parameters: Mapping[str, Any]) -> None:
        self._base_model = base
        self._runtime_parameters = parameters
        self.name = base.name
        self.particles = base.particles
        self.vertices = base.vertices
        for attribute in (
            "alpha_s_mz",
            "alpha_s_me_check",
            "alpha_ew",
            "sin_weak",
            "sqrt_s",
        ):
            if hasattr(base, attribute):
                setattr(self, attribute, getattr(base, attribute))

    def mass(self, pdg: int) -> Any:
        particle = self._base_model.particle(pdg)
        name = f"particle.{int(particle.pdg)}.mass"
        return self._runtime_parameters.get(name, self._base_model.mass(pdg))

    def width(self, pdg: int) -> Any:
        particle = self._base_model.particle(pdg)
        name = f"particle.{int(particle.pdg)}.width"
        return self._runtime_parameters.get(name, self._base_model.width(pdg))

    def propagator_lowering_rule(self, particle_id: int, chirality: int = 0) -> Any:
        return self._base_model.propagator_lowering_rule(particle_id, chirality)

    def is_chiral_eligible(self, pdg: int) -> bool:
        return self._base_model.is_chiral_eligible(pdg)

    def with_runtime_parameters(
        self,
        parameters: Mapping[str, Any],
    ) -> "_RuntimeParameterizedModel":
        return _RuntimeParameterizedModel(self._base_model, parameters)


def build_generic_stage_compiler_blueprint(
    manifest: GenericProcessManifest | GenericDAG,
    *,
    model: Model | None = None,
    selected_color_sector_ids: set[int] | None = None,
    enable_lc_sector_runtime_selector: bool | None = None,
    runtime_schema: Mapping[str, object] | None = None,
    stage_local_parameter_layout: bool = False,
    progress_callback: StageBlueprintProgress | None = None,
    stage_consumer: StageBlueprintConsumer | None = None,
    release_consumed_expressions: bool = False,
) -> GenericStageCompilerBlueprint:
    """Build evaluator-ready symbolic stage metadata for schema-v2 DAGs.

    This is intentionally separate from the legacy shared-current table path.
    It consumes the process-generic current DAG and runtime schema, asks the
    model for local vertex and propagator component expressions, and records
    stage output slots in terms of schema-v2 value slots.  The output is the
    narrow bridge toward serialized Symbolica evaluators and Rusticol schema-v2
    execution.
    """

    selected_model = model or _manifest_model(manifest)
    generic_manifest = (
        manifest
        if isinstance(manifest, GenericProcessManifest)
        else GenericProcessManifest(
            dag=manifest,
            model=selected_model,
            color_plan=manifest.color_plan,
        )
    )
    if enable_lc_sector_runtime_selector is None:
        enable_lc_sector_runtime_selector = selected_color_sector_ids is None
    schema = _dict(
        runtime_schema
        if runtime_schema is not None
        else _generic_runtime_schema_payload(
            generic_manifest.dag,
            generic_manifest.model,
            selected_color_sector_ids=selected_color_sector_ids,
            enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
        )
    )
    parameter_layout = _dict(schema["parameter_layout"])
    global_value_component_count = int(parameter_layout["value_component_count"])
    global_momentum_parameter_count = int(
        parameter_layout["momentum_parameter_count"]
    )
    global_model_parameter_count = int(
        parameter_layout.get("model_parameter_count", 0)
    )
    global_parameter_count = (
        global_value_component_count
        + global_momentum_parameter_count
        + global_model_parameter_count
    )
    if stage_local_parameter_layout:
        # Every compiled stage constructs its own compact symbols below. Avoid
        # materializing the large, otherwise-unused global Symbolica input set.
        parameter_symbols: tuple[Any, ...] = ()
        value_symbols: tuple[Any, ...] = ()
        momentum_symbols: tuple[Any, ...] = ()
        model_parameter_symbols: tuple[Any, ...] = ()
        global_real_valued_inputs = tuple(
            range(global_value_component_count, global_parameter_count)
        )
    else:
        builder = _parameter_builder(schema)
        parameter_symbols = tuple(builder.parameter_symbols())
        value_symbols = parameter_symbols[:global_value_component_count]
        momentum_start = global_value_component_count
        momentum_stop = momentum_start + global_momentum_parameter_count
        momentum_symbols = parameter_symbols[momentum_start:momentum_stop]
        model_parameter_symbols = parameter_symbols[momentum_stop:]
        global_real_valued_inputs = tuple(
            int(index) for index in builder.real_valued_inputs
        )
    model_parameter_records = tuple(
        _dict(item) for item in _list(schema.get("model_parameters", []))
    )
    model_parameter_symbols_by_name = (
        {}
        if stage_local_parameter_layout
        else _logical_model_parameter_symbols(
            model_parameter_records,
            {
                str(record["name"]): model_parameter_symbols[
                    int(record["parameter_index"])
                ]
                for record in model_parameter_records
            },
        )
    )
    expression_model = (
        selected_model.with_runtime_parameters(model_parameter_symbols_by_name)
        if hasattr(selected_model, "with_runtime_parameters")
        else _RuntimeParameterizedModel(
            selected_model,
            model_parameter_symbols_by_name,
        )
    )
    value_slots = _value_slots_by_id(schema)
    current_slots = _current_slots_by_id(schema)
    momentum_slots = _momentum_slots_by_id(schema)
    stage_records = tuple(_dict(stage) for stage in _list(schema["stages"]))
    stage_total = len(stage_records) + 1
    compiled_stages: list[GenericCompiledStageBlueprint] = []
    for stage_index, stage in enumerate(stage_records, start=1):
        if progress_callback is not None:
            progress_callback("current stage", stage_index, stage_total)
        compiled_stage = _compile_current_stage_blueprint(
            generic_manifest.dag,
            expression_model,
            stage,
            value_slots=value_slots,
            current_slots=current_slots,
            momentum_slots=momentum_slots,
            global_value_component_count=global_value_component_count,
            global_momentum_parameter_count=global_momentum_parameter_count,
            model_parameter_records=model_parameter_records,
            global_parameter_symbols=parameter_symbols,
            global_value_symbols=value_symbols,
            global_momentum_symbols=momentum_symbols,
            global_model_parameter_symbols=model_parameter_symbols_by_name,
            global_real_valued_inputs=global_real_valued_inputs,
            stage_local_parameter_layout=stage_local_parameter_layout,
        )
        if stage_consumer is not None:
            stage_consumer(compiled_stage, stage_index - 1, len(stage_records))
        if release_consumed_expressions and stage_consumer is not None:
            compiled_stage = replace(
                compiled_stage,
                parameter_symbols=(),
                output_expressions=(),
                symbolica_functions=(),
            )
        compiled_stages.append(compiled_stage)
    stages = tuple(compiled_stages)
    if progress_callback is not None:
        progress_callback("amplitude stage", stage_total, stage_total)
    amplitude_stage = _compile_amplitude_stage_blueprint(
        expression_model,
        _dict(schema["amplitude_stage"]),
        value_slots=value_slots,
        global_value_component_count=global_value_component_count,
        global_momentum_parameter_count=global_momentum_parameter_count,
        model_parameter_records=model_parameter_records,
        global_parameter_symbols=parameter_symbols,
        global_value_symbols=value_symbols,
        global_model_parameter_symbols=model_parameter_symbols_by_name,
        global_real_valued_inputs=global_real_valued_inputs,
        stage_local_parameter_layout=stage_local_parameter_layout,
    )
    if stage_consumer is not None:
        stage_consumer(amplitude_stage, len(stage_records), len(stage_records))
    if release_consumed_expressions and stage_consumer is not None:
        amplitude_stage = replace(
            amplitude_stage,
            parameter_symbols=(),
            output_expressions=(),
            symbolica_functions=(),
        )
    blockers = tuple(
        blocker
        for stage in (*stages, amplitude_stage)
        for blocker in stage.blockers
    )
    return GenericStageCompilerBlueprint(
        kind="pyamplicol-generic-stage-compiler-blueprint",
        runtime_available=False,
        parameter_count=global_parameter_count,
        value_parameter_count=int(schema["parameter_layout"]["value_component_count"]),
        momentum_parameter_count=global_momentum_parameter_count,
        model_parameter_count=global_model_parameter_count,
        real_valued_inputs=global_real_valued_inputs,
        stage_count=len(stages) + 1,
        stages=stages,
        amplitude_stage=amplitude_stage,
        expression_ready=not blockers,
        blockers=blockers,
        parameter_symbols=parameter_symbols,
    )


def _chunk_evaluation_occurrence_count(
    current_order: Sequence[int],
    *,
    output_size_by_current: Mapping[int, int],
    evaluation_groups_by_current: Mapping[int, frozenset[int]],
    chunk_size: int,
) -> int:
    chunk_groups: list[set[int]] = []
    output_offset = 0
    for current_id in current_order:
        output_size = int(output_size_by_current[current_id])
        if output_size <= 0:
            continue
        first_chunk = output_offset // chunk_size
        last_chunk = (output_offset + output_size - 1) // chunk_size
        while len(chunk_groups) <= last_chunk:
            chunk_groups.append(set())
        groups = evaluation_groups_by_current.get(current_id, frozenset())
        for chunk_index in range(first_chunk, last_chunk + 1):
            chunk_groups[chunk_index].update(groups)
        output_offset += output_size
    return sum(len(groups) for groups in chunk_groups)


def _fanout_aware_current_order(
    current_ids: Sequence[int],
    *,
    output_size_by_current: Mapping[int, int],
    evaluation_groups_by_current: Mapping[int, frozenset[int]],
    chunk_size: int,
) -> tuple[tuple[int, ...], int, int]:
    """Cluster current outputs whose kernel evaluations have shared fan-out.

    Shared evaluation groups form a sparse hypergraph over result currents.
    Indexed heaps reproduce the overlap/benefit greedy choice without scanning
    every unplaced current.  Pathologically large fan-outs use a deterministic
    anchor-group sort, keeping the construction bounded for very large stages.
    """

    natural_order = tuple(int(current_id) for current_id in current_ids)
    before = _chunk_evaluation_occurrence_count(
        natural_order,
        output_size_by_current=output_size_by_current,
        evaluation_groups_by_current=evaluation_groups_by_current,
        chunk_size=chunk_size,
    )
    frequencies = Counter(
        group_id
        for current_id in natural_order
        for group_id in evaluation_groups_by_current.get(current_id, frozenset())
    )
    shared_groups = {
        group_id for group_id, frequency in frequencies.items() if frequency > 1
    }
    if not shared_groups or len(natural_order) < 2:
        return natural_order, before, before

    shared_by_current = {
        current_id: tuple(
            sorted(
                group_id
                for group_id in evaluation_groups_by_current.get(
                    current_id,
                    frozenset(),
                )
                if group_id in shared_groups
            )
        )
        for current_id in natural_order
    }
    members_by_group: dict[int, list[int]] = {
        group_id: [] for group_id in shared_groups
    }
    for current_id, group_ids in shared_by_current.items():
        for group_id in group_ids:
            members_by_group[group_id].append(current_id)

    benefit_by_current = {
        current_id: sum(
            frequencies[group_id] - 1
            for group_id in shared_by_current[current_id]
        )
        for current_id in natural_order
    }
    large_fanout_limit = max(1024, 8 * chunk_size)
    if max(frequencies.values()) > large_fanout_limit:

        def anchor_key(current_id: int) -> tuple[int, int, int, int, int]:
            group_ids = shared_by_current[current_id]
            if not group_ids:
                return (1, 0, 0, 0, current_id)
            anchor = min(
                group_ids,
                key=lambda group_id: (-frequencies[group_id], group_id),
            )
            return (
                0,
                -frequencies[anchor],
                anchor,
                -benefit_by_current[current_id],
                current_id,
            )

        candidate_order = tuple(sorted(natural_order, key=anchor_key))
    else:
        remaining = set(natural_order)
        seed_heap = [
            (
                -benefit_by_current[current_id],
                -len(
                    evaluation_groups_by_current.get(
                        current_id,
                        frozenset(),
                    )
                ),
                current_id,
            )
            for current_id in natural_order
        ]
        heapq.heapify(seed_heap)
        fitting_heaps: dict[int, list[tuple[int, int]]] = {}
        for current_id in natural_order:
            size = int(output_size_by_current[current_id])
            fitting_heaps.setdefault(size, []).append(
                (-benefit_by_current[current_id], current_id)
            )
        for heap in fitting_heaps.values():
            heapq.heapify(heap)
        output_sizes = tuple(sorted(fitting_heaps))

        def pop_seed() -> int:
            while seed_heap:
                _benefit, _group_count, current_id = heapq.heappop(seed_heap)
                if current_id in remaining:
                    return current_id
            raise ValueError("fan-out ordering lost an unplaced current")

        def pop_fitting(capacity: int) -> int | None:
            selected_size: int | None = None
            selected_key: tuple[int, int, int] | None = None
            for size in output_sizes:
                if size > capacity:
                    break
                heap = fitting_heaps[size]
                while heap and heap[0][1] not in remaining:
                    heapq.heappop(heap)
                if not heap:
                    continue
                benefit, current_id = heap[0]
                key = (benefit, -size, current_id)
                if selected_key is None or key < selected_key:
                    selected_size = size
                    selected_key = key
            if selected_size is None:
                return None
            _benefit, current_id = heapq.heappop(
                fitting_heaps[selected_size]
            )
            return current_id

        bins: list[list[int]] = []
        while remaining:
            seed = pop_seed()
            current_bin: list[int] = []
            groups_in_bin: set[int] = set()
            overlap_by_current: dict[int, int] = {}
            candidate_heap: list[tuple[int, int, int, int]] = []
            used = 0

            def add_current(current_id: int) -> None:
                nonlocal used
                remaining.remove(current_id)
                current_bin.append(current_id)
                used += int(output_size_by_current[current_id])
                for group_id in shared_by_current[current_id]:
                    if group_id in groups_in_bin:
                        continue
                    groups_in_bin.add(group_id)
                    for candidate_id in members_by_group[group_id]:
                        if candidate_id not in remaining:
                            continue
                        overlap = overlap_by_current.get(candidate_id, 0) + 1
                        overlap_by_current[candidate_id] = overlap
                        heapq.heappush(
                            candidate_heap,
                            (
                                -overlap,
                                -benefit_by_current[candidate_id],
                                -int(output_size_by_current[candidate_id]),
                                candidate_id,
                            ),
                        )

            def pop_overlapping(capacity: int) -> int | None:
                postponed: list[tuple[int, int, int, int]] = []
                selected: int | None = None
                while candidate_heap:
                    item = heapq.heappop(candidate_heap)
                    overlap, _benefit, _size, current_id = item
                    if current_id not in remaining:
                        continue
                    if -overlap != overlap_by_current.get(current_id, 0):
                        continue
                    if int(output_size_by_current[current_id]) > capacity:
                        postponed.append(item)
                        continue
                    selected = current_id
                    break
                for item in postponed:
                    heapq.heappush(candidate_heap, item)
                return selected

            add_current(seed)
            while remaining:
                capacity = chunk_size - used
                if capacity < 0:
                    break
                selected = pop_overlapping(capacity)
                if selected is None:
                    selected = pop_fitting(capacity)
                if selected is None:
                    break
                add_current(selected)
            bins.append(current_bin)

        candidate_order = tuple(
            current_id for current_bin in bins for current_id in current_bin
        )
    after = _chunk_evaluation_occurrence_count(
        candidate_order,
        output_size_by_current=output_size_by_current,
        evaluation_groups_by_current=evaluation_groups_by_current,
        chunk_size=chunk_size,
    )
    if after >= before:
        return natural_order, before, before
    return candidate_order, before, after


def _stage_with_fanout_aware_output_order(
    stage: GenericCompiledStageBlueprint,
    *,
    chunk_size: int | None,
) -> GenericCompiledStageBlueprint:
    if (
        chunk_size is None
        or int(chunk_size) < 1
        or not stage.evaluation_groups_by_current
        or not stage.output_slots
        or str(stage.stage_kind).startswith("amplitude")
    ):
        return stage

    slots_by_current: dict[int, list[GenericStageOutputSlot]] = {}
    for slot in stage.output_slots:
        slots_by_current.setdefault(slot.current_id, []).append(slot)
    natural_order = tuple(slots_by_current)
    output_size_by_current = {
        current_id: sum(slot.output_stop - slot.output_start for slot in slots)
        for current_id, slots in slots_by_current.items()
    }
    groups_by_current = {
        int(current_id): frozenset(int(group_id) for group_id in group_ids)
        for current_id, group_ids in stage.evaluation_groups_by_current
        if current_id in slots_by_current
    }
    current_order, before, after = _fanout_aware_current_order(
        natural_order,
        output_size_by_current=output_size_by_current,
        evaluation_groups_by_current=groups_by_current,
        chunk_size=int(chunk_size),
    )
    if current_order == natural_order:
        return replace(
            stage,
            fanout_chunk_size=int(chunk_size),
            fanout_evaluation_occurrences_before=before,
            fanout_evaluation_occurrences_after=after,
        )

    outputs: list[Any] = []
    output_slots: list[GenericStageOutputSlot] = []
    for current_id in current_order:
        for slot in slots_by_current[current_id]:
            components = stage.output_expressions[slot.output_start : slot.output_stop]
            start = len(outputs)
            outputs.extend(components)
            output_slots.append(
                replace(
                    slot,
                    output_start=start,
                    output_stop=len(outputs),
                )
            )
    if len(outputs) != len(stage.output_expressions):
        raise ValueError("fan-out output ordering lost stage expressions")
    return replace(
        stage,
        output_slots=tuple(output_slots),
        first_output_previews=_expression_previews(outputs),
        output_expressions=tuple(outputs),
        fanout_chunk_size=int(chunk_size),
        fanout_evaluation_occurrences_before=before,
        fanout_evaluation_occurrences_after=after,
    )


def _prepare_stage_for_output_chunking(
    stage: GenericCompiledStageBlueprint,
    *,
    blueprint: GenericStageCompilerBlueprint | None,
    symbolica_settings: Any | None,
    current_stage_position: int | None = None,
    current_stage_count: int | None = None,
) -> GenericCompiledStageBlueprint:
    if symbolica_settings is None:
        return stage
    settings = _stage_symbolica_settings(
        stage,
        blueprint,
        symbolica_settings,
        current_stage_position=current_stage_position,
        current_stage_count=current_stage_count,
    )
    return _stage_with_fanout_aware_output_order(
        stage,
        chunk_size=getattr(settings, "compiled_output_chunk_size", None),
    )


def write_generic_stage_evaluator_artifacts(
    blueprint: GenericStageCompilerBlueprint,
    artifact_dir: str | Path,
    *,
    compiler: StageEvaluatorCompiler | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
) -> dict[str, object]:
    """Serialize evaluator artifacts for a generic schema-v2 stage blueprint.

    The function is intentionally opt-in.  Normal schema-v2 process manifest
    generation remains cheap, while this path is the bridge from the
    process-generic current DAG to concrete Symbolica evaluator artifacts.
    Rusticol can validate, load, and execute the resulting schema-v2 metadata
    through its generic staged runtime.
    """

    if not blueprint.expression_ready:
        raise ValueError(
            "cannot write generic evaluator artifacts with lowering blockers: "
            + "; ".join(blueprint.blockers)
        )

    output_dir = Path(artifact_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    build_started = time.perf_counter()
    if progress_callback is not None:
        progress_callback(
            {
                "stage": "stage compile",
                "item": "start",
                "total": blueprint.stage_count,
            }
        )

    def compile_stage(stage: GenericCompiledStageBlueprint) -> dict[str, object]:
        return _compile_stage_evaluator_artifact(
            stage,
            output_dir,
            compiler=compiler,
            blueprint=blueprint,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
        )

    stage_payloads = []
    stage_timings: list[dict[str, object]] = []
    for stage in blueprint.stages:
        prepared_stage = _prepare_stage_for_output_chunking(
            stage,
            blueprint=blueprint,
            symbolica_settings=symbolica_settings,
        )
        payload = prepared_stage.to_json_dict()
        payload["evaluator"] = compile_stage(prepared_stage)
        stage_timings.append(
            _stage_build_timing_record(
                prepared_stage.evaluator_label,
                payload["evaluator"],
            )
        )
        stage_payloads.append(payload)
        if progress_callback is not None:
            timing = stage_timings[-1]
            progress_callback(
                {
                    "stage": "stage complete",
                    "item": stage.evaluator_label,
                    "increment": 1,
                    "total": blueprint.stage_count,
                    "duration_s": timing["stage_evaluator_build_s"],
                }
            )

    prepared_amplitude_stage = _prepare_stage_for_output_chunking(
        blueprint.amplitude_stage,
        blueprint=blueprint,
        symbolica_settings=symbolica_settings,
    )
    amplitude_payload = prepared_amplitude_stage.to_json_dict()
    amplitude_payload["evaluator"] = compile_stage(prepared_amplitude_stage)
    stage_timings.append(
        _stage_build_timing_record(
            blueprint.amplitude_stage.evaluator_label,
            amplitude_payload["evaluator"],
        )
    )
    if progress_callback is not None:
        timing = stage_timings[-1]
        progress_callback(
            {
                "stage": "stage complete",
                "item": blueprint.amplitude_stage.evaluator_label,
                "increment": 1,
                "total": blueprint.stage_count,
                "duration_s": timing["stage_evaluator_build_s"],
            }
        )

    return _finalize_stage_evaluator_payload(
        blueprint,
        stage_payloads=stage_payloads,
        amplitude_payload=amplitude_payload,
        stage_timings=stage_timings,
        build_started=build_started,
    )


def write_model_parameter_evaluator_artifact(
    model: Model,
    runtime_schema: Mapping[str, object],
    artifact_dir: str | Path,
    *,
    symbolica_settings: Any | None = None,
    jit_compile: bool = True,
) -> dict[str, object] | None:
    schema = _dict(runtime_schema)
    records = tuple(
        sorted(
            (_dict(item) for item in _list(schema.get("model_parameters", []))),
            key=lambda item: int(item["parameter_index"]),
        )
    )
    input_records = tuple(
        record
        for record in records
        if str(record.get("kind"))
        in {"external_parameter", "external_parameter_component"}
    )
    derived_components: dict[str, dict[str, int]] = {}
    for record in records:
        if str(record.get("kind")) != "derived_parameter_component":
            continue
        runtime_name = record.get("runtime_name")
        component = record.get("complex_component")
        if isinstance(runtime_name, str) and component in {"real", "imag"}:
            derived_components.setdefault(runtime_name, {})[str(component)] = int(
                record["parameter_index"]
            )
    requested_output_names = tuple(
        name
        for name, components in sorted(
            derived_components.items(),
            key=lambda item: min(item[1].values()),
        )
        if set(components) == {"real", "imag"}
    )
    if not requested_output_names:
        return None

    definitions_provider = getattr(
        model,
        "runtime_derived_parameter_definitions",
        None,
    )
    if not callable(definitions_provider):
        return None
    definitions_subset_provider = getattr(
        model,
        "runtime_derived_parameter_definitions_for",
        None,
    )
    definition_values = (
        definitions_subset_provider(requested_output_names)
        if callable(definitions_subset_provider)
        else definitions_provider()
    )
    definitions = {
        str(name): str(expression)
        for name, expression in definition_values.items()
        if str(name) in requested_output_names
    }
    output_names = tuple(
        name for name in requested_output_names if name in definitions
    )
    if not output_names:
        return None

    from symbolica import E, S

    builder = ParamBuilder()
    parameter_symbols = tuple(
        builder.add_parameter_list(
            ("generic_schema_v2", "external_model_parameters"),
            len(input_records),
            role="generic_external_model_parameters",
            real_valued=True,
        )
    )
    slot_symbols = {
        str(record["name"]): parameter_symbols[index]
        for index, record in enumerate(input_records)
    }
    logical_symbols = _logical_model_parameter_symbols(
        input_records,
        slot_symbols,
    )
    outputs = []
    for name in output_names:
        expression = E(definitions[name])
        for parameter_name, symbol in logical_symbols.items():
            expression = expression.replace(S(f"UFO::{parameter_name}"), symbol)
        outputs.append(expression)

    stage = GenericCompiledStageBlueprint(
        stage_index=0,
        stage_kind="model-parameter-derivation",
        subset_size=None,
        evaluator_label="generic_model_parameter_derivation",
        parameter_layout="external-model-parameters",
        output_length=len(outputs),
        output_slots=(),
        input_value_slot_ids=(),
        output_value_slot_ids=(),
        interaction_ids=(),
        input_components=(),
        parameter_count=len(parameter_symbols),
        value_parameter_count=0,
        momentum_parameter_count=0,
        model_parameter_count=len(parameter_symbols),
        real_valued_inputs=tuple(range(len(parameter_symbols))),
        expression_ready=True,
        blockers=(),
        first_output_previews=tuple(
            expression.to_canonical_string()[:_EXPRESSION_PREVIEW_LIMIT]
            for expression in outputs[:3]
        ),
        parameter_symbols=parameter_symbols,
        output_expressions=tuple(outputs),
    )
    parameter_evaluator_settings = (
        None
        if symbolica_settings is None
        else replace(
            symbolica_settings,
            compiled_output_chunk_size=None,
            output_chunk_strategy="uniform",
        )
    )
    evaluator = _compile_stage_evaluator_artifact(
        stage,
        Path(artifact_dir).expanduser(),
        compiler=None,
        blueprint=None,
        symbolica_settings=parameter_evaluator_settings,
        merge_evaluators_strategy=False,
        verbose_evaluator_build=False,
        jit_compile=jit_compile,
        progress_callback=None,
    )
    return {
        "kind": "generic-model-parameter-evaluator",
        "input_parameter_indices": [
            int(record["parameter_index"]) for record in input_records
        ],
        "outputs": [
            {
                "runtime_name": name,
                "output_index": output_index,
                "real_parameter_index": derived_components[name]["real"],
                "imag_parameter_index": derived_components[name]["imag"],
            }
            for output_index, name in enumerate(output_names)
        ],
        "evaluator": evaluator,
    }


def _finalize_stage_evaluator_payload(
    blueprint: GenericStageCompilerBlueprint,
    *,
    stage_payloads: list[dict[str, object]],
    amplitude_payload: dict[str, object],
    stage_timings: list[dict[str, object]],
    build_started: float,
    total_build_s_override: float | None = None,
) -> dict[str, object]:
    stage_local_layout = (
        blueprint.amplitude_stage.parameter_layout == "stage-local-value-momentum"
        and all(
            stage.parameter_layout == "stage-local-value-momentum"
            for stage in blueprint.stages
        )
    )
    total_build_s = (
        time.perf_counter() - build_started
        if total_build_s_override is None
        else float(total_build_s_override)
    )
    jit_compile_s = sum(
        float(record.get("jit_compile_s") or 0.0)
        for record in stage_timings
    )
    timing_totals: dict[str, object] = {}
    non_additive_stage_timing_keys = {
        "output_chunk_autotune_batch_size",
        "output_chunk_autotune_baseline_size",
        "output_chunk_autotune_selected_size",
        "output_chunk_autotune_baseline_us",
        "output_chunk_autotune_selected_us",
        "output_chunk_autotune_gain",
        "output_chunk_autotune_shared_pack_scoring",
    }
    for record in stage_timings:
        for key, value in record.items():
            if (
                key == "evaluator_label"
                or key == "jit_compile_s"
                or key in non_additive_stage_timing_keys
            ):
                continue
            if isinstance(value, (float, int)):
                timing_totals[key] = float(timing_totals.get(key, 0.0)) + float(value)
    timing_totals["stage_evaluator_build_s"] = total_build_s
    timing_totals["jit_compile_s"] = jit_compile_s
    timing_totals["jit_fraction_of_stage_evaluator_build"] = (
        None if total_build_s <= 0.0 else jit_compile_s / total_build_s
    )
    timing_totals["stages"] = stage_timings
    return {
        "kind": "generic-dag-stage-evaluator-artifacts",
        "runtime_available": True,
        "runtime_unavailable_message": None,
        "parameter_count": 0 if stage_local_layout else blueprint.parameter_count,
        "value_parameter_count": (
            0 if stage_local_layout else blueprint.value_parameter_count
        ),
        "momentum_parameter_count": (
            0 if stage_local_layout else blueprint.momentum_parameter_count
        ),
        "model_parameter_count": (
            0 if stage_local_layout else blueprint.model_parameter_count
        ),
        "real_valued_inputs": (
            [] if stage_local_layout else list(blueprint.real_valued_inputs)
        ),
        "parameter_layout": (
            "stage-local-value-momentum"
            if stage_local_layout
            else "global-value-momentum"
        ),
        "stage_count": blueprint.stage_count,
        "build_timing": timing_totals,
        "stages": stage_payloads,
        "amplitude_stage": amplitude_payload,
    }


def _compile_stage_evaluator_artifact(
    stage: GenericCompiledStageBlueprint,
    artifact_dir: Path,
    *,
    compiler: StageEvaluatorCompiler | None,
    blueprint: GenericStageCompilerBlueprint | None,
    symbolica_settings: Any | None,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    jit_compile: bool,
    progress_callback: Any | None,
    current_stage_position: int | None = None,
    current_stage_count: int | None = None,
) -> dict[str, object]:
    if not stage.output_expressions:
        raise ValueError(
            f"generic stage {stage.evaluator_label!r} has no output expressions"
        )
    started = time.perf_counter()
    if compiler is not None:
        manifest = compiler(
            stage,
            stage.parameter_symbols,
            stage.real_valued_inputs,
        )
    else:
        manifest = _compile_default_stage_evaluator(
            stage,
            blueprint,
            artifact_dir,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
            current_stage_position=current_stage_position,
            current_stage_count=current_stage_count,
        )
    if not isinstance(manifest, dict):
        raise TypeError(
            f"generic stage compiler for {stage.evaluator_label!r} "
            "did not return a manifest dictionary"
        )
    build_s = time.perf_counter() - started
    manifest.setdefault("build_timing", {})
    timing = manifest["build_timing"]
    if isinstance(timing, dict):
        previous_stage_build_s = timing.get("stage_evaluator_build_s")
        timing["stage_evaluator_build_s"] = build_s
        if previous_stage_build_s is not None:
            timing["stage_compiler_wrapper_s"] = build_s - float(
                previous_stage_build_s
            )
        timing.setdefault("symbolica_evaluator_build_s", build_s)
        if _manifest_uses_jit_evaluator(manifest):
            timing.setdefault("jit_compile_s", build_s)
    return manifest


def build_and_write_generic_stage_evaluator_artifacts(
    manifest: GenericProcessManifest | GenericDAG,
    runtime_schema: Mapping[str, object],
    artifact_dir: str | Path,
    *,
    model: Model | None = None,
    enable_lc_sector_runtime_selector: bool | None = None,
    stage_local_parameter_layout: bool = False,
    compiler: StageEvaluatorCompiler | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    blueprint_progress_callback: StageBlueprintProgress | None = None,
    evaluator_progress_callback: Any | None = None,
) -> tuple[GenericStageCompilerBlueprint, dict[str, object]]:
    """Lower, compile, and release one recursion stage at a time."""

    output_dir = Path(artifact_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    schema = _dict(runtime_schema)
    current_stage_count = len(_list(schema["stages"]))
    stage_count = current_stage_count + 1
    build_started = time.perf_counter()
    if evaluator_progress_callback is not None:
        evaluator_progress_callback(
            {
                "stage": "stage compile",
                "item": "start",
                "total": stage_count,
            }
        )

    stage_payloads: list[dict[str, object]] = []
    amplitude_payload: dict[str, object] | None = None
    stage_timings: list[dict[str, object]] = []

    def consume_stage(
        stage: GenericCompiledStageBlueprint,
        position: int,
        reported_current_stage_count: int,
    ) -> None:
        nonlocal amplitude_payload
        if reported_current_stage_count != current_stage_count:
            raise ValueError("streamed stage count changed during blueprint lowering")
        if not stage.expression_ready:
            raise ValueError(
                "cannot write generic evaluator artifact with lowering blockers: "
                + "; ".join(stage.blockers)
            )
        prepared_stage = _prepare_stage_for_output_chunking(
            stage,
            blueprint=None,
            symbolica_settings=symbolica_settings,
            current_stage_position=position,
            current_stage_count=current_stage_count,
        )
        payload = prepared_stage.to_json_dict()
        payload["evaluator"] = _compile_stage_evaluator_artifact(
            prepared_stage,
            output_dir,
            compiler=compiler,
            blueprint=None,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=evaluator_progress_callback,
            current_stage_position=position,
            current_stage_count=current_stage_count,
        )
        timing = _stage_build_timing_record(
            prepared_stage.evaluator_label,
            payload["evaluator"],
        )
        stage_timings.append(timing)
        if str(prepared_stage.stage_kind).startswith("amplitude"):
            amplitude_payload = payload
        else:
            stage_payloads.append(payload)
        if evaluator_progress_callback is not None:
            evaluator_progress_callback(
                {
                    "stage": "stage complete",
                    "item": prepared_stage.evaluator_label,
                    "increment": 1,
                    "total": stage_count,
                    "duration_s": timing["stage_evaluator_build_s"],
                }
            )

    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        model=model,
        enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
        runtime_schema=schema,
        stage_local_parameter_layout=stage_local_parameter_layout,
        progress_callback=blueprint_progress_callback,
        stage_consumer=consume_stage,
        release_consumed_expressions=True,
    )
    if amplitude_payload is None or len(stage_payloads) != current_stage_count:
        raise ValueError("streamed stage compilation produced incomplete metadata")
    return blueprint, _finalize_stage_evaluator_payload(
        blueprint,
        stage_payloads=stage_payloads,
        amplitude_payload=amplitude_payload,
        stage_timings=stage_timings,
        build_started=build_started,
        total_build_s_override=sum(
            float(timing["stage_evaluator_build_s"])
            for timing in stage_timings
        ),
    )


def _stage_build_timing_record(
    evaluator_label: str,
    evaluator_manifest: object,
) -> dict[str, object]:
    manifest = evaluator_manifest if isinstance(evaluator_manifest, dict) else {}
    raw_timing = manifest.get("build_timing") if isinstance(manifest, dict) else None
    timing = raw_timing if isinstance(raw_timing, dict) else {}
    record: dict[str, object] = {
        "evaluator_label": evaluator_label,
        "stage_evaluator_build_s": float(
            timing.get("stage_evaluator_build_s") or 0.0
        ),
        "symbolica_evaluator_build_s": float(
            timing.get("symbolica_evaluator_build_s") or 0.0
        ),
        "jit_compile_s": (
            None
            if timing.get("jit_compile_s") is None
            else float(timing.get("jit_compile_s") or 0.0)
        ),
    }
    for key, value in timing.items():
        if key in record:
            continue
        if isinstance(value, (float, int)):
            record[str(key)] = float(value)
    return record


def _manifest_uses_jit_evaluator(manifest: Mapping[str, object]) -> bool:
    if str(manifest.get("kind", "")) == "jit-symbolica-evaluator":
        return True
    if str(manifest.get("kind", "")) == "chunked-symbolica-evaluator":
        chunks = manifest.get("chunks")
        if isinstance(chunks, Sequence) and chunks:
            return all(
                isinstance(chunk, Mapping)
                and _manifest_uses_jit_evaluator(chunk)
                for chunk in chunks
            )
    settings = manifest.get("settings")
    if isinstance(settings, Mapping):
        return str(settings.get("backend", "")) == "jit"
    return False


def _compile_default_stage_evaluator(
    stage: GenericCompiledStageBlueprint,
    blueprint: GenericStageCompilerBlueprint | None,
    artifact_dir: Path,
    *,
    symbolica_settings: Any | None,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    jit_compile: bool,
    progress_callback: Any | None,
    current_stage_position: int | None = None,
    current_stage_count: int | None = None,
) -> dict[str, object]:
    from .symbolica_evaluator import (
        SymbolicaEvaluatorSettings,
        _compile_symbolica_outputs,
        _symbolica_evaluator_artifact_manifest,
    )

    settings = _stage_symbolica_settings(
        stage,
        blueprint,
        symbolica_settings or SymbolicaEvaluatorSettings(),
        current_stage_position=current_stage_position,
        current_stage_count=current_stage_count,
    )
    symbolica_started = time.perf_counter()
    params = list(stage.parameter_symbols)

    def compile_with(candidate_settings: Any, candidate_label: str) -> Any:
        return _compile_symbolica_outputs(
            stage.output_expressions,
            params,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            real_params=stage.real_valued_inputs,
            symbolica_settings=candidate_settings,
            jit_compile=jit_compile,
            label=candidate_label,
            progress_callback=progress_callback,
            functions={
                (function, arguments): body
                for function, arguments, body in stage.symbolica_functions
            },
        )

    autotune_timing: dict[str, float] = {}
    if getattr(settings, "output_chunk_strategy", "uniform") == "measured-stage":
        evaluator, autotune_timing = _compile_measured_stage_output_chunks(
            settings=settings,
            output_count=len(stage.output_expressions),
            parameter_count=len(params),
            real_params=stage.real_valued_inputs,
            label=stage.evaluator_label,
            compile_with=compile_with,
            progress_callback=progress_callback,
            jit_compile=jit_compile,
        )
    else:
        evaluator = compile_with(settings, stage.evaluator_label)
    symbolica_build_s = time.perf_counter() - symbolica_started
    artifact_started = time.perf_counter()
    manifest = _symbolica_evaluator_artifact_manifest(evaluator, artifact_dir)
    artifact_manifest_s = time.perf_counter() - artifact_started
    timing = manifest.setdefault("build_timing", {})
    if isinstance(timing, dict):
        timing.update(autotune_timing)
        timing["symbolica_evaluator_build_s"] = symbolica_build_s
        timing["artifact_manifest_s"] = artifact_manifest_s
        timing["stage_evaluator_build_s"] = symbolica_build_s + artifact_manifest_s
    return manifest


def _compile_measured_stage_output_chunks(
    *,
    settings: Any,
    output_count: int,
    parameter_count: int,
    real_params: Sequence[int],
    label: str,
    compile_with: Callable[[Any, str], Any],
    progress_callback: Any | None,
    jit_compile: bool,
) -> tuple[Any, dict[str, float]]:
    base = getattr(settings, "compiled_output_chunk_size", None)
    if base is None or getattr(settings, "backend", None) != "jit" or not jit_compile:
        uniform = replace(settings, output_chunk_strategy="uniform")
        return compile_with(uniform, label), {}

    requested_sizes = (
        int(base),
        max(1, int(base) // 2),
        max(1, 3 * int(base) // 4),
        max(1, 3 * int(base) // 2),
        int(base) * 2,
        None,
    )
    effective_sizes: list[int | None] = []
    for requested in requested_sizes:
        effective = (
            None
            if requested is None or output_count <= int(requested)
            else int(requested)
        )
        if effective not in effective_sizes:
            effective_sizes.append(effective)
    baseline_size = (
        None if output_count <= int(base) else int(base)
    )

    started = time.perf_counter()
    candidates: dict[int | None, Any] = {}
    for chunk_size in effective_sizes:
        suffix = "none" if chunk_size is None else str(chunk_size)
        candidate_settings = replace(
            settings,
            compiled_output_chunk_size=chunk_size,
            output_chunk_strategy="uniform",
        )
        candidate_label = (
            label
            if chunk_size == baseline_size
            else f"{label}_autotune_chunk_{suffix}"
        )
        candidates[chunk_size] = compile_with(candidate_settings, candidate_label)

    autotune_batch_size = int(
        getattr(settings, "output_chunk_autotune_batch_size", 128)
    )
    rows = np.full(
        (autotune_batch_size, parameter_count),
        complex(0.75, 0.125),
        dtype=np.complex128,
    )
    if real_params:
        rows[:, list(real_params)] = 0.75

    materialize_started = time.perf_counter()
    for evaluator in candidates.values():
        evaluator.evaluate_complex(rows)
    materialize_s = time.perf_counter() - materialize_started

    shared_pack_scoring = all(
        bool(getattr(evaluator, "supports_complex_profiled", lambda: False)())
        for evaluator in candidates.values()
    )

    def score_candidate(evaluator: Any) -> float:
        if shared_pack_scoring:
            _output, profile = evaluator.evaluate_complex_profiled(rows)
            return max(sum(float(value) for value in profile), 1.0e-9)
        probe_started = time.perf_counter()
        evaluator.evaluate_complex(rows)
        return max(time.perf_counter() - probe_started, 1.0e-6)

    benchmark_started = time.perf_counter()
    scores: dict[int | None, float] = {}
    repeats: dict[int | None, int] = {}
    for chunk_size, evaluator in candidates.items():
        probe_s = max(score_candidate(evaluator), 1.0e-6)
        repeats[chunk_size] = max(1, min(256, int(0.01 / probe_s)))
        scores[chunk_size] = probe_s

    samples: dict[int | None, list[float]] = {
        chunk_size: [] for chunk_size in candidates
    }
    ordered_sizes = list(candidates)
    for round_index in range(5):
        rotated = ordered_sizes[round_index:] + ordered_sizes[:round_index]
        for chunk_size in rotated:
            evaluator = candidates[chunk_size]
            count = repeats[chunk_size]
            if shared_pack_scoring:
                samples[chunk_size].append(
                    sum(score_candidate(evaluator) for _ in range(count)) / count
                )
            else:
                sample_started = time.perf_counter()
                for _ in range(count):
                    evaluator.evaluate_complex(rows)
                samples[chunk_size].append(
                    (time.perf_counter() - sample_started) / count
                )
    scores = {
        chunk_size: sorted(values)[len(values) // 2]
        for chunk_size, values in samples.items()
    }
    benchmark_s = time.perf_counter() - benchmark_started
    selected_size = _select_measured_chunk_candidate(
        scores,
        baseline_size=baseline_size,
        minimum_gain=0.05,
    )
    selected = candidates[selected_size]
    baseline_score = scores[baseline_size]
    selected_score = scores[selected_size]
    autotune_s = time.perf_counter() - started
    if progress_callback is not None:
        progress_callback(
            {
                "stage": "chunk autotune",
                "item": (
                    f"{label} base={baseline_size or 'none'} "
                    f"selected={selected_size or 'none'} "
                    f"gain={(1.0 - selected_score / baseline_score):.1%}"
                    + (" shared-pack" if shared_pack_scoring else "")
                ),
            }
        )
    evaluator_timing = getattr(selected, "build_timing", None)
    if isinstance(evaluator_timing, dict):
        evaluator_timing.update(
            {
                "output_chunk_autotune_s": autotune_s,
                "output_chunk_autotune_materialize_s": materialize_s,
                "output_chunk_autotune_benchmark_s": benchmark_s,
                "output_chunk_autotune_candidate_count": float(len(candidates)),
                "output_chunk_autotune_batch_size": float(autotune_batch_size),
                "output_chunk_autotune_baseline_size": float(baseline_size or 0),
                "output_chunk_autotune_selected_size": float(selected_size or 0),
                "output_chunk_autotune_baseline_us": baseline_score * 1.0e6,
                "output_chunk_autotune_selected_us": selected_score * 1.0e6,
                "output_chunk_autotune_gain": 1.0 - selected_score / baseline_score,
                "output_chunk_autotune_shared_pack_scoring": float(
                    shared_pack_scoring
                ),
            }
        )
    return selected, {
        "output_chunk_autotune_s": autotune_s,
        "output_chunk_autotune_materialize_s": materialize_s,
        "output_chunk_autotune_benchmark_s": benchmark_s,
        "output_chunk_autotune_candidate_count": float(len(candidates)),
        "output_chunk_autotune_batch_size": float(autotune_batch_size),
        "output_chunk_autotune_baseline_size": float(baseline_size or 0),
        "output_chunk_autotune_selected_size": float(selected_size or 0),
        "output_chunk_autotune_baseline_us": baseline_score * 1.0e6,
        "output_chunk_autotune_selected_us": selected_score * 1.0e6,
        "output_chunk_autotune_gain": 1.0 - selected_score / baseline_score,
        "output_chunk_autotune_shared_pack_scoring": float(shared_pack_scoring),
    }


def _select_measured_chunk_candidate(
    scores: Mapping[int | None, float],
    *,
    baseline_size: int | None,
    minimum_gain: float,
) -> int | None:
    if baseline_size not in scores:
        raise ValueError("measured chunk scores do not include the baseline")
    best_size = min(scores, key=scores.__getitem__)
    if scores[best_size] <= scores[baseline_size] * (1.0 - minimum_gain):
        return best_size
    return baseline_size


def _stage_symbolica_settings(
    stage: GenericCompiledStageBlueprint,
    blueprint: GenericStageCompilerBlueprint | None,
    settings: Any,
    *,
    current_stage_position: int | None = None,
    current_stage_count: int | None = None,
) -> Any:
    """Apply the optional recursion-stage output-chunk taper."""

    strategy = getattr(settings, "output_chunk_strategy", "uniform")
    if strategy == "auto":
        base = getattr(settings, "compiled_output_chunk_size", None)
        output_count = int(
            getattr(
                stage,
                "output_length",
                0 if base is None else int(base) + 1,
            )
        )
        strategy = (
            "measured-stage"
            if getattr(settings, "backend", None) == "jit"
            and base is not None
            and int(base) <= 256
            and output_count > int(base)
            else "uniform"
        )
        settings = replace(settings, output_chunk_strategy=strategy)
    if strategy != "tapered-stage":
        return settings
    base = getattr(settings, "compiled_output_chunk_size", None)
    if base is None:
        return settings
    if str(stage.stage_kind).startswith("amplitude"):
        return replace(settings, compiled_output_chunk_size=None)

    if blueprint is not None:
        try:
            position = next(
                index
                for index, current_stage in enumerate(blueprint.stages)
                if current_stage.stage_index == stage.stage_index
            )
        except StopIteration:
            return settings
        stage_count = len(blueprint.stages)
    else:
        if current_stage_position is None or current_stage_count is None:
            return settings
        position = int(current_stage_position)
        stage_count = int(current_stage_count)
    remaining = stage_count - position - 1
    if position < 2:
        chunk_size = None
    elif remaining == 0:
        chunk_size = max(1, int(base) // 2)
    elif remaining <= 2:
        chunk_size = int(base)
    else:
        chunk_size = int(base) * 2
    return replace(settings, compiled_output_chunk_size=chunk_size)


def _compile_current_stage_blueprint(
    dag: GenericDAG,
    model: Model,
    stage: dict[str, Any],
    *,
    value_slots: dict[int, dict[str, Any]],
    current_slots: dict[int, dict[str, Any]],
    momentum_slots: dict[int, dict[str, Any]],
    global_value_component_count: int,
    global_momentum_parameter_count: int,
    model_parameter_records: Sequence[dict[str, Any]],
    global_parameter_symbols: Sequence[Any],
    global_value_symbols: Sequence[Any],
    global_momentum_symbols: Sequence[Any],
    global_model_parameter_symbols: Mapping[str, Any],
    global_real_valued_inputs: Sequence[int],
    stage_local_parameter_layout: bool,
) -> GenericCompiledStageBlueprint:
    blockers: list[str] = []
    outputs: list[Any] = []
    output_slots: list[GenericStageOutputSlot] = []
    interactions_compacted = bool(stage.get("interactions_compacted", False))
    interactions = (
        []
        if interactions_compacted
        else [_dict(item) for item in _list(stage["interactions"])]
    )
    interaction_ids = (
        tuple(int(value) for value in _list(stage.get("interaction_ids", [])))
        if interactions_compacted
        else tuple(int(interaction["interaction_id"]) for interaction in interactions)
    )
    input_value_slot_ids = tuple(
        int(value) for value in _list(stage["input_value_slot_ids"])
    )
    input_momentum_slot_ids = (
        tuple(
            int(value)
            for value in _list(stage.get("input_momentum_slot_ids", []))
        )
        if interactions_compacted
        else _stage_input_momentum_slot_ids(interactions)
    )
    output_slot_ids = tuple(int(value) for value in _list(stage["output_value_slot_ids"]))
    output_slots_by_current: dict[int, list[dict[str, Any]]] = {}
    for slot_id in output_slot_ids:
        slot = value_slots[int(slot_id)]
        output_slots_by_current.setdefault(int(slot["current_id"]), []).append(slot)
    stage_model_parameter_records = (
        _current_stage_model_parameter_records(
            model,
            model_parameter_records,
            dag=dag,
            interactions=interactions,
            interaction_ids=interaction_ids,
            output_slots_by_current=output_slots_by_current,
            current_slots=current_slots,
        )
        if stage_local_parameter_layout
        else model_parameter_records
    )
    local_inputs = (
        _stage_local_inputs(
            value_slot_ids=input_value_slot_ids,
            momentum_slot_ids=input_momentum_slot_ids,
            value_slots=value_slots,
            momentum_slots=momentum_slots,
            global_value_component_count=global_value_component_count,
            global_momentum_parameter_count=global_momentum_parameter_count,
            model_parameter_records=stage_model_parameter_records,
        )
        if stage_local_parameter_layout
        else _global_stage_inputs(
            parameter_symbols=global_parameter_symbols,
            value_symbols=global_value_symbols,
            momentum_symbols=global_momentum_symbols,
            model_parameter_symbols=global_model_parameter_symbols,
            value_parameter_count=global_value_component_count,
            momentum_parameter_count=len(global_momentum_symbols),
            model_parameter_count=len(model_parameter_records),
            real_valued_inputs=global_real_valued_inputs,
        )
    )
    stage_model = (
        model.with_runtime_parameters(local_inputs.model_parameter_symbols)
        if hasattr(model, "with_runtime_parameters")
        else _RuntimeParameterizedModel(model, local_inputs.model_parameter_symbols)
    )
    interactions_by_result: dict[int, list[dict[str, Any] | int]] = {}
    if interactions_compacted:
        for interaction_id in interaction_ids:
            interactions_by_result.setdefault(
                int(dag.interactions[interaction_id].result_id),
                [],
            ).append(interaction_id)
    else:
        for interaction in interactions:
            interactions_by_result.setdefault(
                int(interaction["result_current_id"]),
                [],
            ).append(interaction)
    value_components_by_slot_id = {
        int(slot_id): _value_components(
            value_slots[int(slot_id)],
            local_inputs.value_symbols,
        )
        for slot_id in input_value_slot_ids
    }
    momentum_components_by_slot_id = {
        int(slot_id): _momentum_components(
            int(slot_id),
            local_inputs.momentum_symbols,
            momentum_slots,
            by_slot_id=True,
        )
        for slot_id in input_momentum_slot_ids
    }
    momentum_components_by_mask = {
        int(slot["momentum_mask"]): momentum_components_by_slot_id[int(slot_id)]
        for slot_id, slot in momentum_slots.items()
        if int(slot_id) in momentum_components_by_slot_id
    }
    input_value_slot_by_current_id = {
        int(value_slots[slot_id]["current_id"]): int(slot_id)
        for slot_id in input_value_slot_ids
    }
    momentum_slot_by_mask = {
        int(slot["momentum_mask"]): int(slot_id)
        for slot_id, slot in momentum_slots.items()
    }
    compact_coupling_cache: dict[
        tuple[int, tuple[int, ...], tuple[float, ...]],
        tuple[Any, Any],
    ] = {}
    compact_evaluation_cache: dict[int, tuple[Any, ...]] = {}
    evaluation_groups_by_current = tuple(
        (
            current_id,
            tuple(
                sorted(
                    {
                        (
                            int(interaction.evaluation_group_id)
                            if interaction.evaluation_group_id is not None
                            else -(interaction.id + 1)
                        )
                        for interaction_item in interactions_by_result[current_id]
                        for interaction in (
                            dag.interactions[
                                int(interaction_item)
                                if interactions_compacted
                                else int(_dict(interaction_item)["interaction_id"])
                            ],
                        )
                    }
                )
            ),
        )
        for current_id in sorted(interactions_by_result)
    )

    for current_id in sorted(interactions_by_result):
        current_slot = current_slots[current_id]
        dimension = int(current_slot["dimension"])
        total = tuple(0j for _ in range(dimension))
        for interaction_item in interactions_by_result[current_id]:
            interaction_id = (
                int(interaction_item)
                if interactions_compacted
                else int(_dict(interaction_item)["interaction_id"])
            )
            try:
                contribution = _compact_interaction_contribution(
                    dag,
                    stage_model,
                    interaction_id,
                    value_components_by_slot_id=value_components_by_slot_id,
                    input_value_slot_by_current_id=input_value_slot_by_current_id,
                    momentum_components_by_slot_id=momentum_components_by_slot_id,
                    momentum_slot_by_mask=momentum_slot_by_mask,
                    model_parameter_symbols=local_inputs.model_parameter_symbols,
                    coupling_cache=compact_coupling_cache,
                    evaluation_cache=compact_evaluation_cache,
                )
            except ValueError as error:
                blockers.append(
                    f"interaction {interaction_id}: {error}"
                )
                continue
            total = _sum_components(total, contribution)
        result_slots = output_slots_by_current.get(current_id, ())
        for slot in result_slots:
            variant = str(slot["variant"])
            try:
                components = (
                    stage_model.propagator_component_expression(
                        int(current_slot["particle_id"]),
                        total,
                        momentum_components_by_mask[int(current_slot["momentum_mask"])],
                        chirality=int(current_slot["chirality"]),
                    )
                    if variant == "propagated"
                    else total
                )
            except ValueError as error:
                blockers.append(f"value slot {slot['value_slot_id']}: {error}")
                continue
            start = len(outputs)
            outputs.extend(components)
            output_slots.append(
                GenericStageOutputSlot(
                    value_slot_id=int(slot["value_slot_id"]),
                    current_id=current_id,
                    variant=variant,
                    component_start=int(slot["component_start"]),
                    component_stop=int(slot["component_stop"]),
                    output_start=start,
                    output_stop=len(outputs),
                )
            )

    specialized_outputs, symbolica_functions = (
        _specialize_stage_symbolica_functions(
            outputs,
            _model_symbolica_functions(stage_model),
        )
    )
    return GenericCompiledStageBlueprint(
        stage_index=int(stage["stage_index"]),
        stage_kind=str(stage["stage_kind"]),
        subset_size=int(stage["subset_size"]),
        evaluator_label=(
            f"generic_stage_{int(stage['stage_index'])}_subset_{int(stage['subset_size'])}"
        ),
        parameter_layout=(
            "stage-local-value-momentum"
            if stage_local_parameter_layout
            else "global-value-momentum"
        ),
        output_length=len(outputs),
        output_slots=tuple(output_slots),
        input_value_slot_ids=input_value_slot_ids,
        output_value_slot_ids=output_slot_ids,
        interaction_ids=interaction_ids,
        input_components=local_inputs.input_components,
        parameter_count=len(local_inputs.parameter_symbols),
        value_parameter_count=local_inputs.value_parameter_count,
        momentum_parameter_count=local_inputs.momentum_parameter_count,
        model_parameter_count=local_inputs.model_parameter_count,
        real_valued_inputs=local_inputs.real_valued_inputs,
        expression_ready=not blockers,
        blockers=tuple(blockers),
        first_output_previews=_expression_previews(specialized_outputs),
        evaluation_groups_by_current=evaluation_groups_by_current,
        parameter_symbols=local_inputs.parameter_symbols,
        output_expressions=specialized_outputs,
        symbolica_functions=symbolica_functions,
    )


def _compile_amplitude_stage_blueprint(
    model: Model,
    stage: dict[str, Any],
    *,
    value_slots: dict[int, dict[str, Any]],
    global_value_component_count: int,
    global_momentum_parameter_count: int,
    model_parameter_records: Sequence[dict[str, Any]],
    global_parameter_symbols: Sequence[Any],
    global_value_symbols: Sequence[Any],
    global_model_parameter_symbols: Mapping[str, Any],
    global_real_valued_inputs: Sequence[int],
    stage_local_parameter_layout: bool,
) -> GenericCompiledStageBlueprint:
    blockers: list[str] = []
    outputs: list[Any] = []
    output_slots: list[GenericStageOutputSlot] = []
    input_value_slot_ids = tuple(
        sorted(
            {
                int(root[side]["value_slot_id"])
                for root in (_dict(item) for item in _list(stage["roots"]))
                for side in ("left_value_slot", "right_value_slot")
            }
        )
    )
    local_inputs = (
        _stage_local_inputs(
            value_slot_ids=input_value_slot_ids,
            momentum_slot_ids=(),
            value_slots=value_slots,
            momentum_slots={},
            global_value_component_count=global_value_component_count,
            global_momentum_parameter_count=global_momentum_parameter_count,
            model_parameter_records=(
                _amplitude_stage_model_parameter_records(
                    model_parameter_records,
                    roots=tuple(_dict(item) for item in _list(stage["roots"])),
                )
                if stage_local_parameter_layout
                else model_parameter_records
            ),
        )
        if stage_local_parameter_layout
        else _global_stage_inputs(
            parameter_symbols=global_parameter_symbols,
            value_symbols=global_value_symbols,
            momentum_symbols=(),
            model_parameter_symbols=global_model_parameter_symbols,
            value_parameter_count=global_value_component_count,
            momentum_parameter_count=0,
            model_parameter_count=len(model_parameter_records),
            real_valued_inputs=global_real_valued_inputs,
        )
    )
    stage_model = (
        model.with_runtime_parameters(local_inputs.model_parameter_symbols)
        if hasattr(model, "with_runtime_parameters")
        else _RuntimeParameterizedModel(model, local_inputs.model_parameter_symbols)
    )
    lc_sector_selector = local_inputs.model_parameter_symbols.get(
        LC_SECTOR_SELECTOR_PARAMETER
    )
    for root in (_dict(item) for item in _list(stage["roots"])):
        try:
            output = _amplitude_root_expression(
                stage_model,
                root,
                value_symbols=local_inputs.value_symbols,
                model_parameter_symbols=local_inputs.model_parameter_symbols,
                value_slots=value_slots,
            )
        except ValueError as error:
            blockers.append(f"amplitude root {root['root_id']}: {error}")
            continue
        if lc_sector_selector is not None:
            output = _lc_sector_guard(
                output,
                selector=lc_sector_selector,
                sector_id=int(root["color_sector_id"]),
            )
        start = len(outputs)
        outputs.append(output)
        output_slots.append(
            GenericStageOutputSlot(
                value_slot_id=-1,
                current_id=-1,
                variant="amplitude-root",
                component_start=int(root["output_index"]),
                component_stop=int(root["output_index"]) + 1,
                output_start=start,
                output_stop=len(outputs),
            )
        )
    return GenericCompiledStageBlueprint(
        stage_index=0,
        stage_kind=str(stage["stage_kind"]),
        subset_size=None,
        evaluator_label="generic_amplitude_stage",
        parameter_layout=(
            "stage-local-value-momentum"
            if stage_local_parameter_layout
            else "global-value-momentum"
        ),
        output_length=len(outputs),
        output_slots=tuple(output_slots),
        input_value_slot_ids=input_value_slot_ids,
        output_value_slot_ids=(),
        interaction_ids=(),
        input_components=local_inputs.input_components,
        parameter_count=len(local_inputs.parameter_symbols),
        value_parameter_count=local_inputs.value_parameter_count,
        momentum_parameter_count=local_inputs.momentum_parameter_count,
        model_parameter_count=local_inputs.model_parameter_count,
        real_valued_inputs=local_inputs.real_valued_inputs,
        expression_ready=not blockers,
        blockers=tuple(blockers),
        first_output_previews=_expression_previews(outputs),
        parameter_symbols=local_inputs.parameter_symbols,
        output_expressions=tuple(outputs),
    )


def _interaction_contribution(
    dag: GenericDAG,
    model: Model,
    interaction: dict[str, Any],
    *,
    value_components_by_slot_id: Mapping[int, tuple[Any, ...]],
    momentum_components_by_slot_id: Mapping[int, tuple[Any, ...]],
    model_parameter_symbols: Mapping[str, Any],
) -> tuple[Any, ...]:
    left_slot = _dict(interaction["left_value_slot"])
    right_slot = _dict(interaction["right_value_slot"])
    left = value_components_by_slot_id[int(left_slot["value_slot_id"])]
    right = value_components_by_slot_id[int(right_slot["value_slot_id"])]
    momenta = _dict(interaction["momentum_slots"])
    left_current = dag.currents[int(interaction["left_current_id"])]
    right_current = dag.currents[int(interaction["right_current_id"])]
    result_current = dag.currents[int(interaction["result_current_id"])]
    components = model.vertex_component_expression(
        int(interaction["vertex_kind"]),
        left,
        right,
        result_particle_id=int(result_current.index.particle_id),
        result_chirality=int(result_current.index.chirality),
        left_chirality=int(left_current.index.chirality),
        right_chirality=int(right_current.index.chirality),
        coupling=_runtime_coupling(interaction, model_parameter_symbols),
        left_momentum=momentum_components_by_slot_id[int(momenta["left"])],
        right_momentum=momentum_components_by_slot_id[int(momenta["right"])],
    )
    color_weight = _coupling(interaction.get("color_weight"))
    if color_weight == (1.0, 0.0):
        return components
    weight = color_weight[0] + 1j * color_weight[1]
    return tuple(weight * component for component in components)


def _compact_interaction_contribution(
    dag: GenericDAG,
    model: Model,
    interaction_id: int,
    *,
    value_components_by_slot_id: Mapping[int, tuple[Any, ...]],
    input_value_slot_by_current_id: Mapping[int, int],
    momentum_components_by_slot_id: Mapping[int, tuple[Any, ...]],
    momentum_slot_by_mask: Mapping[int, int],
    model_parameter_symbols: Mapping[str, Any],
    coupling_cache: dict[
        tuple[int, tuple[int, ...], tuple[float, ...]],
        tuple[Any, Any],
    ],
    evaluation_cache: dict[int, tuple[Any, ...]],
) -> tuple[Any, ...]:
    interaction = dag.interactions[interaction_id]
    evaluation_group_id = interaction.evaluation_group_id
    canonical_components = (
        None
        if evaluation_group_id is None
        else evaluation_cache.get(evaluation_group_id)
    )
    evaluation_factor = complex(*interaction.evaluation_factor)
    if evaluation_factor == 0j:
        raise ValueError("interaction evaluation factor must be nonzero")
    if canonical_components is None:
        left_current = dag.currents[interaction.left_id]
        right_current = dag.currents[interaction.right_id]
        result_current = dag.currents[interaction.result_id]
        left = value_components_by_slot_id[
            input_value_slot_by_current_id[interaction.left_id]
        ]
        right = value_components_by_slot_id[
            input_value_slot_by_current_id[interaction.right_id]
        ]
        coupling_key = (
            int(interaction.vertex_kind),
            interaction.vertex_particles,
            interaction.coupling,
        )
        coupling = coupling_cache.get(coupling_key)
        if coupling is None:
            resolved_coupling = list(interaction.coupling)
            names = _runtime_coupling_parameter_names(
                interaction.vertex_kind,
                interaction.vertex_particles,
                interaction.coupling,
                model=model,
            )
            for index, name in enumerate(names):
                if isinstance(name, str) and name in model_parameter_symbols:
                    resolved_coupling[index] = model_parameter_symbols[name]
            coupling = (resolved_coupling[0], resolved_coupling[1])
            coupling_cache[coupling_key] = coupling
        components = model.vertex_component_expression(
            int(interaction.vertex_kind),
            left,
            right,
            result_particle_id=int(result_current.index.particle_id),
            result_chirality=int(result_current.index.chirality),
            left_chirality=int(left_current.index.chirality),
            right_chirality=int(right_current.index.chirality),
            coupling=coupling,
            left_momentum=momentum_components_by_slot_id[
                momentum_slot_by_mask[left_current.index.momentum_mask]
            ],
            right_momentum=momentum_components_by_slot_id[
                momentum_slot_by_mask[right_current.index.momentum_mask]
            ],
        )
        canonical_components = (
            components
            if evaluation_factor == 1.0 + 0.0j
            else tuple(component / evaluation_factor for component in components)
        )
        if evaluation_group_id is not None:
            evaluation_cache[evaluation_group_id] = canonical_components
    color_weight = complex(*interaction.color_weight)
    attachment_weight = color_weight * evaluation_factor
    if attachment_weight == 1.0 + 0.0j:
        return canonical_components
    return tuple(
        attachment_weight * component for component in canonical_components
    )


def _amplitude_root_expression(
    model: Model,
    root: dict[str, Any],
    *,
    value_symbols: Sequence[Any],
    model_parameter_symbols: Mapping[str, Any],
    value_slots: dict[int, dict[str, Any]],
) -> Any:
    left = _value_components(_dict(root["left_value_slot"]), value_symbols)
    right = _value_components(_dict(root["right_value_slot"]), value_symbols)
    kind = str(root["kind"])
    contraction = str(root.get("contraction", ""))
    coupling = _runtime_coupling(root, model_parameter_symbols)
    color_weight = _coupling(root.get("color_weight"))
    weight = color_weight[0] + 1j * color_weight[1]
    if kind == "direct-contraction":
        return weight * _contract_components(contraction, left, right)
    if kind == "vertex-closure":
        vertex_kind = root.get("vertex_kind")
        particles = root.get("vertex_particles")
        if vertex_kind is None or not isinstance(particles, list) or len(particles) != 3:
            raise ValueError("vertex closure is missing vertex metadata")
        components = model.vertex_component_expression(
            int(vertex_kind),
            left,
            right,
            result_particle_id=int(particles[2]),
            result_chirality=0,
            coupling=coupling,
        )
        if contraction == "scalar" and len(components) == 1:
            return weight * components[0]
        raise ValueError(
            f"vertex closure contraction {contraction!r} is not scalar-lowered"
        )
    raise ValueError(f"unsupported amplitude root kind {kind!r}")


def _contract_components(contraction: str, left: Sequence[Any], right: Sequence[Any]) -> Any:
    if contraction == "scalar":
        return left[0] * right[0]
    if contraction == "weyl":
        return left[0] * right[0] + left[1] * right[1]
    if contraction == "dirac":
        return sum(left[index] * right[index] for index in range(min(len(left), len(right))))
    if contraction == "lorentz":
        return left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3]
    if contraction == "antisymmetric-tensor":
        return sum(left[index] * right[index] for index in range(min(len(left), len(right))))
    raise ValueError(f"unsupported direct contraction {contraction!r}")


def _stage_input_momentum_slot_ids(
    interactions: Sequence[dict[str, Any]],
) -> tuple[int, ...]:
    slot_ids: set[int] = set()
    for interaction in interactions:
        momentum_slots = _dict(interaction["momentum_slots"])
        for key in ("left", "right", "result"):
            slot_ids.add(int(momentum_slots[key]))
    return tuple(sorted(slot_ids))


def _current_stage_model_parameter_records(
    model: Model,
    model_parameter_records: Sequence[dict[str, Any]],
    *,
    dag: GenericDAG,
    interactions: Sequence[dict[str, Any]],
    interaction_ids: Sequence[int],
    output_slots_by_current: Mapping[int, Sequence[dict[str, Any]]],
    current_slots: Mapping[int, dict[str, Any]],
) -> tuple[dict[str, Any], ...]:
    used_names = _coupling_parameter_names_used_by_records(interactions)
    if not interactions:
        processed_signatures: set[
            tuple[int, tuple[int, ...], tuple[float, ...]]
        ] = set()
        for interaction_id in interaction_ids:
            interaction = dag.interactions[int(interaction_id)]
            signature = (
                int(interaction.vertex_kind),
                interaction.vertex_particles,
                interaction.coupling,
            )
            if signature in processed_signatures:
                continue
            processed_signatures.add(signature)
            used_names.update(
                name
                for name in _runtime_coupling_parameter_names(
                    interaction.vertex_kind,
                    interaction.vertex_particles,
                    interaction.coupling,
                    model=model,
                )
                if isinstance(name, str)
            )
    if _has_model_parameter_record(
        model_parameter_records,
        LC_SECTOR_SELECTOR_PARAMETER,
    ):
        used_names.add(LC_SECTOR_SELECTOR_PARAMETER)
    for current_id, slots in output_slots_by_current.items():
        if not any(str(slot["variant"]) == "propagated" for slot in slots):
            continue
        current_slot = current_slots[int(current_id)]
        used_names.update(
            _particle_model_parameter_names(
                model,
                int(current_slot["particle_id"]),
            )
        )
    return _filter_model_parameter_records(model_parameter_records, used_names)


def _amplitude_stage_model_parameter_records(
    model_parameter_records: Sequence[dict[str, Any]],
    *,
    roots: Sequence[dict[str, Any]],
) -> tuple[dict[str, Any], ...]:
    used_names = _coupling_parameter_names_used_by_records(roots)
    if _has_model_parameter_record(
        model_parameter_records,
        LC_SECTOR_SELECTOR_PARAMETER,
    ):
        used_names.add(LC_SECTOR_SELECTOR_PARAMETER)
    return _filter_model_parameter_records(model_parameter_records, used_names)


def _has_model_parameter_record(
    model_parameter_records: Sequence[dict[str, Any]],
    name: str,
) -> bool:
    return any(str(record.get("name")) == name for record in model_parameter_records)


def _lc_sector_guard(
    expression: Any,
    *,
    selector: Any,
    sector_id: int,
) -> Any:
    from symbolica import Expression

    selected = Expression.IF(selector - int(sector_id), 0, expression)
    return Expression.IF(selector + 1, selected, expression)


def _coupling_parameter_names_used_by_records(
    records: Sequence[dict[str, Any]],
) -> set[str]:
    used_names: set[str] = set()
    for record in records:
        names = record.get("coupling_parameter_names")
        if not isinstance(names, list):
            continue
        for name in names:
            if isinstance(name, str):
                used_names.add(name)
    return used_names


def _particle_model_parameter_names(model: Model, pdg: int) -> tuple[str, ...]:
    runtime_names = getattr(model, "runtime_parameter_names_for_particle", None)
    if callable(runtime_names):
        return tuple(str(name) for name in runtime_names(int(pdg)))
    try:
        particle = model.particle(pdg)
    except KeyError:
        return ()
    return (
        f"particle.{int(particle.pdg)}.mass",
        f"particle.{int(particle.pdg)}.width",
    )


def _filter_model_parameter_records(
    model_parameter_records: Sequence[dict[str, Any]],
    used_names: set[str],
) -> tuple[dict[str, Any], ...]:
    return tuple(
        record
        for record in sorted(
            model_parameter_records,
            key=lambda item: int(item["parameter_index"]),
        )
        if str(record.get("runtime_name", record["name"])) in used_names
        and not (
            record.get("complex_domain") == "real"
            and record.get("complex_component") == "imag"
        )
        and not (
            record.get("complex_domain") == "imaginary"
            and record.get("complex_component") == "real"
        )
    )


def _logical_model_parameter_symbols(
    model_parameter_records: Sequence[dict[str, Any]],
    slot_symbols: Mapping[str, Any],
) -> dict[str, Any]:
    logical_symbols: dict[str, Any] = {}
    complex_components: dict[str, dict[str, Any]] = {}
    complex_domains: dict[str, str] = {}
    for record in model_parameter_records:
        slot_name = str(record["name"])
        symbol = slot_symbols[slot_name]
        runtime_name = record.get("runtime_name")
        component = record.get("complex_component")
        if isinstance(runtime_name, str) and component in {"real", "imag"}:
            complex_components.setdefault(runtime_name, {})[str(component)] = symbol
            domain = str(record.get("complex_domain", "complex"))
            previous = complex_domains.setdefault(runtime_name, domain)
            if previous != domain:
                raise ValueError(
                    f"runtime model parameter {runtime_name!r} has conflicting domains"
                )
        else:
            logical_symbols[slot_name] = symbol
    for runtime_name, components in complex_components.items():
        domain = complex_domains[runtime_name]
        if "real" not in components and domain == "imaginary":
            components["real"] = 0.0
        if "imag" not in components and domain == "real":
            components["imag"] = 0.0
        if set(components) != {"real", "imag"}:
            raise ValueError(
                f"runtime model parameter {runtime_name!r} is missing a real or imaginary slot"
            )
        logical_symbols[runtime_name] = (
            components["real"] + 1j * components["imag"]
        )
    return logical_symbols


def _stage_local_inputs(
    *,
    value_slot_ids: Sequence[int],
    momentum_slot_ids: Sequence[int],
    value_slots: Mapping[int, dict[str, Any]],
    momentum_slots: Mapping[int, dict[str, Any]],
    global_value_component_count: int,
    global_momentum_parameter_count: int,
    model_parameter_records: Sequence[dict[str, Any]],
) -> _StageLocalInputs:
    builder = ParamBuilder()
    input_components: list[GenericStageInputComponent] = []
    value_symbols: dict[int, tuple[Any, ...]] = {}
    momentum_symbols: dict[int, tuple[Any, ...]] = {}
    model_parameter_slot_symbols: dict[str, Any] = {}

    value_spans = tuple(
        (
            int(value_slot_id),
            int(value_slots[int(value_slot_id)]["component_start"]),
            int(value_slots[int(value_slot_id)]["component_stop"]),
        )
        for value_slot_id in value_slot_ids
    )
    momentum_spans = tuple(
        (
            int(momentum_slot_id),
            int(momentum_slots[int(momentum_slot_id)]["component_start"]),
            int(momentum_slots[int(momentum_slot_id)]["component_stop"]),
        )
        for momentum_slot_id in momentum_slot_ids
    )
    sorted_model_parameter_records = tuple(
        sorted(
            model_parameter_records,
            key=lambda item: int(item["parameter_index"]),
        )
    )
    value_parameter_count = sum(stop - start for _, start, stop in value_spans)
    momentum_parameter_count = sum(
        stop - start for _, start, stop in momentum_spans
    )
    model_parameter_count = len(sorted_model_parameter_records)
    parameter_count = (
        value_parameter_count + momentum_parameter_count + model_parameter_count
    )
    parameter_symbols = (
        builder.add_parameter_list(
            ("generic_schema_v2_stage", "inputs"),
            parameter_count,
            role="generic_stage_input_storage",
        )
        if parameter_count
        else ()
    )
    real_parameter_start = value_parameter_count
    builder.real_valued_inputs.extend(range(real_parameter_start, parameter_count))
    parameter_cursor = 0

    for value_slot_id, start, stop in value_spans:
        length = stop - start
        local_symbols = parameter_symbols[parameter_cursor : parameter_cursor + length]
        value_symbols[value_slot_id] = local_symbols
        parameter_start = len(input_components)
        for component, global_component in enumerate(range(start, stop)):
            input_components.append(
                GenericStageInputComponent(
                    kind="value",
                    source_id=int(value_slot_id),
                    component=component,
                    global_component=global_component,
                    parameter_index=parameter_start + component,
                    real_valued=False,
                )
            )
        parameter_cursor += length

    for momentum_slot_id, start, stop in momentum_spans:
        length = stop - start
        local_symbols = parameter_symbols[parameter_cursor : parameter_cursor + length]
        momentum_symbols[momentum_slot_id] = local_symbols
        parameter_start = len(input_components)
        for component, local_component in enumerate(range(start, stop)):
            input_components.append(
                GenericStageInputComponent(
                    kind="momentum",
                    source_id=int(momentum_slot_id),
                    component=component,
                    global_component=global_value_component_count + local_component,
                    parameter_index=parameter_start + component,
                    real_valued=True,
                )
            )
        parameter_cursor += length

    model_parameter_global_start = (
        global_value_component_count + global_momentum_parameter_count
    )
    for record in sorted_model_parameter_records:
        name = str(record["name"])
        parameter_index = int(record["parameter_index"])
        symbol = parameter_symbols[parameter_cursor]
        model_parameter_slot_symbols[name] = symbol
        input_components.append(
            GenericStageInputComponent(
                kind="model_parameter",
                source_id=parameter_index,
                component=0,
                global_component=model_parameter_global_start + parameter_index,
                parameter_index=len(input_components),
                real_valued=True,
            )
        )
        parameter_cursor += 1

    if parameter_cursor != parameter_count:
        raise RuntimeError("stage-local parameter layout cursor mismatch")

    model_parameter_symbols = _logical_model_parameter_symbols(
        sorted_model_parameter_records,
        model_parameter_slot_symbols,
    )

    return _StageLocalInputs(
        parameter_symbols=tuple(parameter_symbols),
        input_components=tuple(input_components),
        value_symbols=value_symbols,
        momentum_symbols=momentum_symbols,
        model_parameter_symbols=model_parameter_symbols,
        value_parameter_count=value_parameter_count,
        momentum_parameter_count=momentum_parameter_count,
        model_parameter_count=model_parameter_count,
        real_valued_inputs=tuple(int(index) for index in builder.real_valued_inputs),
    )


def _global_stage_inputs(
    *,
    parameter_symbols: Sequence[Any],
    value_symbols: Sequence[Any],
    momentum_symbols: Sequence[Any],
    model_parameter_symbols: Mapping[str, Any],
    value_parameter_count: int,
    momentum_parameter_count: int,
    model_parameter_count: int,
    real_valued_inputs: Sequence[int],
) -> _StageLocalInputs:
    return _StageLocalInputs(
        parameter_symbols=tuple(parameter_symbols),
        input_components=(),
        value_symbols=tuple(value_symbols),
        momentum_symbols=tuple(momentum_symbols),
        model_parameter_symbols=dict(model_parameter_symbols),
        value_parameter_count=int(value_parameter_count),
        momentum_parameter_count=int(momentum_parameter_count),
        model_parameter_count=int(model_parameter_count),
        real_valued_inputs=tuple(int(index) for index in real_valued_inputs),
    )


def _parameter_builder(schema: dict[str, Any]) -> ParamBuilder:
    layout = _dict(schema["parameter_layout"])
    builder = ParamBuilder()
    value_component_count = int(layout["value_component_count"])
    if value_component_count:
        builder.add_parameter_list(
            ("generic_schema_v2", "values"),
            value_component_count,
            role="generic_value_storage",
        )
    momentum_parameter_count = int(layout["momentum_parameter_count"])
    if momentum_parameter_count:
        builder.add_parameter_list(
            ("generic_schema_v2", "momenta"),
            momentum_parameter_count,
            role="generic_momentum_storage",
            real_valued=True,
        )
    model_parameter_count = int(layout.get("model_parameter_count", 0))
    if model_parameter_count:
        builder.add_parameter_list(
            ("generic_schema_v2", "model_parameters"),
            model_parameter_count,
            role="generic_model_parameters",
            real_valued=True,
        )
    return builder


def _manifest_model(manifest: GenericProcessManifest | GenericDAG) -> Model:
    if isinstance(manifest, GenericProcessManifest):
        return manifest.model
    return AmplicolSMLeadingColorModel()


def _model_symbolica_functions(
    model: Model,
) -> tuple[tuple[Any, tuple[Any, ...], Any], ...]:
    getter = getattr(model, "symbolica_function_definitions", None)
    if not callable(getter):
        return ()
    definitions = getter()
    if not isinstance(definitions, Mapping):
        raise TypeError("model Symbolica function definitions must be a mapping")
    return tuple(
        (function, tuple(arguments), body)
        for (function, arguments), body in definitions.items()
    )


def _specialize_stage_symbolica_functions(
    outputs: Sequence[Any],
    definitions: Sequence[tuple[Any, tuple[Any, ...], Any]],
) -> tuple[
    tuple[Any, ...],
    tuple[tuple[Any, tuple[Any, ...], Any], ...],
]:
    """Inline model-kernel calls before constructing a stage evaluator."""

    output_expressions = tuple(outputs)
    function_definitions = tuple(definitions)
    if not output_expressions or not function_definitions:
        return output_expressions, ()

    from symbolica import Replacement, S

    replacements: list[Any] = []
    for definition_index, (function, arguments, body) in enumerate(
        function_definitions
    ):
        if int(body.get_byte_size()) > _MODEL_FUNCTION_INLINE_MAX_BYTES:
            continue
        wildcards = tuple(
            S(
                "pyamplicol_inline_function_"
                f"{definition_index}_argument_{argument_index}_"
            )
            for argument_index in range(len(arguments))
        )
        pattern = function(*wildcards)
        replacement = body
        for argument, wildcard in zip(arguments, wildcards, strict=True):
            replacement = replacement.replace(
                argument,
                wildcard,
                allow_new_wildcards_on_rhs=True,
            )
        replacements.append(
            Replacement(
                pattern,
                replacement,
                allow_new_wildcards_on_rhs=True,
            )
        )

    rewritten = tuple(
        expression.replace_multiple(replacements, repeat=True)
        for expression in output_expressions
    )

    function_names = {function for function, _arguments, _body in function_definitions}
    required = {
        symbol
        for expression in rewritten
        for symbol in _expression_symbols(expression, enter_functions=True)
        if symbol in function_names
    }
    while True:
        dependencies = {
            symbol
            for function, _arguments, body in function_definitions
            if function in required
            for symbol in _expression_symbols(body, enter_functions=True)
            if symbol in function_names
        }
        expanded = required | dependencies
        if expanded == required:
            break
        required = expanded
    retained = tuple(
        definition
        for definition in function_definitions
        if definition[0] in required
    )
    return rewritten, retained


def _expression_symbols(
    expression: Any,
    *,
    enter_functions: bool,
) -> set[Any]:
    getter = getattr(expression, "get_all_symbols", None)
    if not callable(getter):
        return set()
    return set(getter(enter_functions))


def _value_components(
    slot: dict[str, Any],
    value_symbols: Sequence[Any] | Mapping[int, tuple[Any, ...]],
) -> tuple[Any, ...]:
    if isinstance(value_symbols, Mapping):
        return tuple(value_symbols[int(slot["value_slot_id"])])
    start = int(slot["component_start"])
    stop = int(slot["component_stop"])
    return tuple(value_symbols[index] for index in range(start, stop))


def _momentum_components(
    key: int,
    momentum_symbols: Sequence[Any] | Mapping[int, tuple[Any, ...]],
    momentum_slots: dict[int, dict[str, Any]],
    *,
    by_slot_id: bool = False,
) -> tuple[Any, ...]:
    slot = momentum_slots[key] if by_slot_id else _momentum_slot_by_mask(momentum_slots, key)
    if isinstance(momentum_symbols, Mapping):
        return tuple(momentum_symbols[int(slot["momentum_slot_id"])])
    start = int(slot["component_start"])
    stop = int(slot["component_stop"])
    return tuple(momentum_symbols[index] for index in range(start, stop))


def _momentum_slot_by_mask(
    momentum_slots: dict[int, dict[str, Any]],
    momentum_mask: int,
) -> dict[str, Any]:
    for slot in momentum_slots.values():
        if int(slot["momentum_mask"]) == momentum_mask:
            return slot
    raise ValueError(f"no momentum slot for mask {momentum_mask}")


def _value_slots_by_id(schema: dict[str, Any]) -> dict[int, dict[str, Any]]:
    storage = _dict(schema["value_storage"])
    return {
        int(slot["value_slot_id"]): _dict(slot)
        for slot in _list(storage["value_slots"])
    }


def _current_slots_by_id(schema: dict[str, Any]) -> dict[int, dict[str, Any]]:
    storage = _dict(schema["current_storage"])
    return {
        int(slot["current_id"]): _dict(slot)
        for slot in _list(storage["current_slots"])
    }


def _momentum_slots_by_id(schema: dict[str, Any]) -> dict[int, dict[str, Any]]:
    return {
        int(slot["momentum_slot_id"]): _dict(slot)
        for slot in _list(schema["momentum_slots"])
    }


def _sum_components(left: Sequence[Any], right: Sequence[Any]) -> tuple[Any, ...]:
    if len(left) != len(right):
        raise ValueError(f"component dimensions differ: {len(left)} != {len(right)}")
    return tuple(left[index] + right[index] for index in range(len(left)))


def _coupling(value: object) -> tuple[Any, Any]:
    if not isinstance(value, list | tuple) or len(value) != 2:
        raise ValueError("coupling metadata must have two entries")
    return value[0], value[1]


def _runtime_coupling(
    record: Mapping[str, Any],
    model_parameter_symbols: Mapping[str, Any],
) -> tuple[Any, Any]:
    values = list(_coupling(record.get("coupling")))
    names = record.get("coupling_parameter_names")
    if isinstance(names, list):
        for index, name in enumerate(names[: len(values)]):
            if isinstance(name, str) and name in model_parameter_symbols:
                values[index] = model_parameter_symbols[name]
    return values[0], values[1]


def _expression_previews(expressions: Sequence[Any]) -> tuple[str, ...]:
    return tuple(_preview_expression(expression) for expression in expressions[:4])


def _preview_expression(expression: Any) -> str:
    text = str(expression)
    if len(text) <= _EXPRESSION_PREVIEW_LIMIT:
        return text
    return text[:_EXPRESSION_PREVIEW_LIMIT] + "...<truncated>"


def _dict(value: object) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise TypeError("expected JSON object")
    return value


def _list(value: object) -> list[Any]:
    if not isinstance(value, list):
        raise TypeError("expected JSON array")
    return value


__all__ = [
    "GenericCompiledStageBlueprint",
    "GenericStageCompilerBlueprint",
    "GenericStageOutputSlot",
    "build_generic_stage_compiler_blueprint",
    "write_generic_stage_evaluator_artifacts",
]
