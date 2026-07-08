from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from .generic_artifact import GenericProcessManifest
from .generic_dag import GenericDAG
from .model import AmplicolSMLeadingColorModel, Model
from .params import ParamBuilder

_EXPRESSION_PREVIEW_LIMIT = 512


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
    real_valued_inputs: tuple[int, ...]
    expression_ready: bool
    blockers: tuple[str, ...]
    first_output_previews: tuple[str, ...]
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
            "real_valued_inputs": list(self.real_valued_inputs),
            "expression_ready": self.expression_ready,
            "blockers": list(self.blockers),
            "first_output_previews": list(self.first_output_previews),
        }


@dataclass(frozen=True)
class GenericStageCompilerBlueprint:
    kind: str
    runtime_available: bool
    parameter_count: int
    value_parameter_count: int
    momentum_parameter_count: int
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
    value_parameter_count: int
    momentum_parameter_count: int
    real_valued_inputs: tuple[int, ...]


StageEvaluatorCompiler = Callable[
    [GenericCompiledStageBlueprint, tuple[Any, ...], tuple[int, ...]],
    dict[str, object],
]


def build_generic_stage_compiler_blueprint(
    manifest: GenericProcessManifest | GenericDAG,
    *,
    model: Model | None = None,
    selected_color_sector_ids: set[int] | None = None,
    stage_local_parameter_layout: bool = False,
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
    payload = generic_manifest.to_json_dict(
        selected_color_sector_ids=selected_color_sector_ids,
    )
    schema = _dict(payload["runtime_schema"])
    builder = _parameter_builder(schema)
    parameter_symbols = tuple(builder.parameter_symbols())
    global_value_component_count = int(
        schema["parameter_layout"]["value_component_count"]
    )
    value_symbols = parameter_symbols[:global_value_component_count]
    momentum_symbols = parameter_symbols[global_value_component_count:]
    global_real_valued_inputs = tuple(int(index) for index in builder.real_valued_inputs)
    value_slots = _value_slots_by_id(schema)
    current_slots = _current_slots_by_id(schema)
    momentum_slots = _momentum_slots_by_id(schema)
    stages = tuple(
        _compile_current_stage_blueprint(
            generic_manifest.dag,
            selected_model,
            _dict(stage),
            value_slots=value_slots,
            current_slots=current_slots,
            momentum_slots=momentum_slots,
            global_value_component_count=global_value_component_count,
            global_parameter_symbols=parameter_symbols,
            global_value_symbols=value_symbols,
            global_momentum_symbols=momentum_symbols,
            global_real_valued_inputs=global_real_valued_inputs,
            stage_local_parameter_layout=stage_local_parameter_layout,
        )
        for stage in _list(schema["stages"])
    )
    amplitude_stage = _compile_amplitude_stage_blueprint(
        selected_model,
        _dict(schema["amplitude_stage"]),
        value_slots=value_slots,
        global_value_component_count=global_value_component_count,
        global_parameter_symbols=parameter_symbols,
        global_value_symbols=value_symbols,
        global_real_valued_inputs=global_real_valued_inputs,
        stage_local_parameter_layout=stage_local_parameter_layout,
    )
    blockers = tuple(
        blocker
        for stage in (*stages, amplitude_stage)
        for blocker in stage.blockers
    )
    return GenericStageCompilerBlueprint(
        kind="pyamplicol-generic-stage-compiler-blueprint",
        runtime_available=False,
        parameter_count=len(parameter_symbols),
        value_parameter_count=int(schema["parameter_layout"]["value_component_count"]),
        momentum_parameter_count=int(schema["parameter_layout"]["momentum_parameter_count"]),
        real_valued_inputs=global_real_valued_inputs,
        stage_count=len(stages) + 1,
        stages=stages,
        amplitude_stage=amplitude_stage,
        expression_ready=not blockers,
        blockers=blockers,
        parameter_symbols=parameter_symbols,
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

    def compile_stage(stage: GenericCompiledStageBlueprint) -> dict[str, object]:
        if not stage.output_expressions:
            raise ValueError(
                f"generic stage {stage.evaluator_label!r} has no output expressions"
            )
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
                output_dir,
                symbolica_settings=symbolica_settings,
                merge_evaluators_strategy=merge_evaluators_strategy,
                verbose_evaluator_build=verbose_evaluator_build,
                jit_compile=jit_compile,
                progress_callback=progress_callback,
            )
        if not isinstance(manifest, dict):
            raise TypeError(
                f"generic stage compiler for {stage.evaluator_label!r} "
                "did not return a manifest dictionary"
            )
        return manifest

    stage_payloads = []
    for stage in blueprint.stages:
        payload = stage.to_json_dict()
        payload["evaluator"] = compile_stage(stage)
        stage_payloads.append(payload)

    amplitude_payload = blueprint.amplitude_stage.to_json_dict()
    amplitude_payload["evaluator"] = compile_stage(blueprint.amplitude_stage)

    stage_local_layout = (
        blueprint.amplitude_stage.parameter_layout == "stage-local-value-momentum"
        and all(stage.parameter_layout == "stage-local-value-momentum" for stage in blueprint.stages)
    )
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
        "real_valued_inputs": (
            [] if stage_local_layout else list(blueprint.real_valued_inputs)
        ),
        "parameter_layout": (
            "stage-local-value-momentum"
            if stage_local_layout
            else "global-value-momentum"
        ),
        "stage_count": blueprint.stage_count,
        "stages": stage_payloads,
        "amplitude_stage": amplitude_payload,
    }


def _compile_default_stage_evaluator(
    stage: GenericCompiledStageBlueprint,
    blueprint: GenericStageCompilerBlueprint,
    artifact_dir: Path,
    *,
    symbolica_settings: Any | None,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    jit_compile: bool,
    progress_callback: Any | None,
) -> dict[str, object]:
    from .symbolica_evaluator import (
        SymbolicaEvaluatorSettings,
        _compile_symbolica_outputs,
        _symbolica_evaluator_artifact_manifest,
    )

    settings = symbolica_settings or SymbolicaEvaluatorSettings()
    evaluator = _compile_symbolica_outputs(
        stage.output_expressions,
        list(stage.parameter_symbols),
        merge_evaluators_strategy=merge_evaluators_strategy,
        verbose_evaluator_build=verbose_evaluator_build,
        real_params=stage.real_valued_inputs,
        symbolica_settings=settings,
        jit_compile=jit_compile,
        label=stage.evaluator_label,
        progress_callback=progress_callback,
    )
    return _symbolica_evaluator_artifact_manifest(evaluator, artifact_dir)


def _compile_current_stage_blueprint(
    dag: GenericDAG,
    model: Model,
    stage: dict[str, Any],
    *,
    value_slots: dict[int, dict[str, Any]],
    current_slots: dict[int, dict[str, Any]],
    momentum_slots: dict[int, dict[str, Any]],
    global_value_component_count: int,
    global_parameter_symbols: Sequence[Any],
    global_value_symbols: Sequence[Any],
    global_momentum_symbols: Sequence[Any],
    global_real_valued_inputs: Sequence[int],
    stage_local_parameter_layout: bool,
) -> GenericCompiledStageBlueprint:
    blockers: list[str] = []
    outputs: list[Any] = []
    output_slots: list[GenericStageOutputSlot] = []
    interactions = [_dict(item) for item in _list(stage["interactions"])]
    input_value_slot_ids = tuple(
        int(value) for value in _list(stage["input_value_slot_ids"])
    )
    input_momentum_slot_ids = _stage_input_momentum_slot_ids(interactions)
    local_inputs = (
        _stage_local_inputs(
            value_slot_ids=input_value_slot_ids,
            momentum_slot_ids=input_momentum_slot_ids,
            value_slots=value_slots,
            momentum_slots=momentum_slots,
            global_value_component_count=global_value_component_count,
        )
        if stage_local_parameter_layout
        else _global_stage_inputs(
            parameter_symbols=global_parameter_symbols,
            value_symbols=global_value_symbols,
            momentum_symbols=global_momentum_symbols,
            value_parameter_count=global_value_component_count,
            momentum_parameter_count=len(global_momentum_symbols),
            real_valued_inputs=global_real_valued_inputs,
        )
    )
    interactions_by_result: dict[int, list[dict[str, Any]]] = {}
    for interaction in interactions:
        interactions_by_result.setdefault(
            int(interaction["result_current_id"]),
            [],
        ).append(interaction)
    output_slots_by_current: dict[int, list[dict[str, Any]]] = {}
    for slot_id in _list(stage["output_value_slot_ids"]):
        slot = value_slots[int(slot_id)]
        output_slots_by_current.setdefault(int(slot["current_id"]), []).append(slot)

    for current_id in sorted(interactions_by_result):
        current_slot = current_slots[current_id]
        dimension = int(current_slot["dimension"])
        total = tuple(0j for _ in range(dimension))
        for interaction in interactions_by_result[current_id]:
            try:
                contribution = _interaction_contribution(
                    dag,
                    model,
                    interaction,
                    value_symbols=local_inputs.value_symbols,
                    momentum_symbols=local_inputs.momentum_symbols,
                    value_slots=value_slots,
                    momentum_slots=momentum_slots,
                )
            except ValueError as error:
                blockers.append(
                    f"interaction {interaction['interaction_id']}: {error}"
                )
                continue
            total = _sum_components(total, contribution)
        result_slots = output_slots_by_current.get(current_id, ())
        for slot in result_slots:
            variant = str(slot["variant"])
            try:
                components = (
                    model.propagator_component_expression(
                        int(current_slot["particle_id"]),
                        total,
                        _momentum_components(
                            int(current_slot["momentum_mask"]),
                            local_inputs.momentum_symbols,
                            momentum_slots,
                        ),
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
        output_value_slot_ids=tuple(int(value) for value in _list(stage["output_value_slot_ids"])),
        interaction_ids=tuple(int(interaction["interaction_id"]) for interaction in interactions),
        input_components=local_inputs.input_components,
        parameter_count=len(local_inputs.parameter_symbols),
        value_parameter_count=local_inputs.value_parameter_count,
        momentum_parameter_count=local_inputs.momentum_parameter_count,
        real_valued_inputs=local_inputs.real_valued_inputs,
        expression_ready=not blockers,
        blockers=tuple(blockers),
        first_output_previews=_expression_previews(outputs),
        parameter_symbols=local_inputs.parameter_symbols,
        output_expressions=tuple(outputs),
    )


def _compile_amplitude_stage_blueprint(
    model: Model,
    stage: dict[str, Any],
    *,
    value_slots: dict[int, dict[str, Any]],
    global_value_component_count: int,
    global_parameter_symbols: Sequence[Any],
    global_value_symbols: Sequence[Any],
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
        )
        if stage_local_parameter_layout
        else _global_stage_inputs(
            parameter_symbols=global_parameter_symbols,
            value_symbols=global_value_symbols,
            momentum_symbols=(),
            value_parameter_count=global_value_component_count,
            momentum_parameter_count=0,
            real_valued_inputs=global_real_valued_inputs,
        )
    )
    for root in (_dict(item) for item in _list(stage["roots"])):
        try:
            output = _amplitude_root_expression(
                model,
                root,
                value_symbols=local_inputs.value_symbols,
                value_slots=value_slots,
            )
        except ValueError as error:
            blockers.append(f"amplitude root {root['root_id']}: {error}")
            continue
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
    value_symbols: Sequence[Any],
    momentum_symbols: Sequence[Any],
    value_slots: dict[int, dict[str, Any]],
    momentum_slots: dict[int, dict[str, Any]],
) -> tuple[Any, ...]:
    left_slot = _dict(interaction["left_value_slot"])
    right_slot = _dict(interaction["right_value_slot"])
    left = _value_components(left_slot, value_symbols)
    right = _value_components(right_slot, value_symbols)
    momenta = _dict(interaction["momentum_slots"])
    left_current = dag.currents[int(interaction["left_current_id"])]
    right_current = dag.currents[int(interaction["right_current_id"])]
    result_current = dag.currents[int(interaction["result_current_id"])]
    return model.vertex_component_expression(
        int(interaction["vertex_kind"]),
        left,
        right,
        result_particle_id=int(result_current.index.particle_id),
        result_chirality=int(result_current.index.chirality),
        left_chirality=int(left_current.index.chirality),
        right_chirality=int(right_current.index.chirality),
        coupling=_coupling(interaction.get("coupling")),
        left_momentum=_momentum_components(
            int(momenta["left"]),
            momentum_symbols,
            momentum_slots,
            by_slot_id=True,
        ),
        right_momentum=_momentum_components(
            int(momenta["right"]),
            momentum_symbols,
            momentum_slots,
            by_slot_id=True,
        ),
    )


def _amplitude_root_expression(
    model: Model,
    root: dict[str, Any],
    *,
    value_symbols: Sequence[Any],
    value_slots: dict[int, dict[str, Any]],
) -> Any:
    left = _value_components(_dict(root["left_value_slot"]), value_symbols)
    right = _value_components(_dict(root["right_value_slot"]), value_symbols)
    kind = str(root["kind"])
    contraction = str(root.get("contraction", ""))
    coupling = _coupling(root.get("coupling"))
    color_weight = _coupling(root.get("color_weight"))
    weight = color_weight[0]
    if color_weight[1] != 0.0:
        raise ValueError("complex color weights are not lowered in LC stage blueprint")
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


def _stage_local_inputs(
    *,
    value_slot_ids: Sequence[int],
    momentum_slot_ids: Sequence[int],
    value_slots: Mapping[int, dict[str, Any]],
    momentum_slots: Mapping[int, dict[str, Any]],
    global_value_component_count: int,
) -> _StageLocalInputs:
    builder = ParamBuilder()
    input_components: list[GenericStageInputComponent] = []
    value_symbols: dict[int, tuple[Any, ...]] = {}
    momentum_symbols: dict[int, tuple[Any, ...]] = {}

    for value_slot_id in value_slot_ids:
        slot = value_slots[int(value_slot_id)]
        start = int(slot["component_start"])
        stop = int(slot["component_stop"])
        symbols = builder.add_parameter_list(
            ("generic_schema_v2_stage", "value", str(value_slot_id)),
            stop - start,
            role="generic_stage_value_storage",
        )
        value_symbols[int(value_slot_id)] = symbols
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

    for momentum_slot_id in momentum_slot_ids:
        slot = momentum_slots[int(momentum_slot_id)]
        start = int(slot["component_start"])
        stop = int(slot["component_stop"])
        symbols = builder.add_parameter_list(
            ("generic_schema_v2_stage", "momentum", str(momentum_slot_id)),
            stop - start,
            role="generic_stage_momentum_storage",
            real_valued=True,
        )
        momentum_symbols[int(momentum_slot_id)] = symbols
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

    return _StageLocalInputs(
        parameter_symbols=tuple(builder.parameter_symbols()),
        input_components=tuple(input_components),
        value_symbols=value_symbols,
        momentum_symbols=momentum_symbols,
        value_parameter_count=sum(
            int(value_slots[int(value_slot_id)]["component_stop"])
            - int(value_slots[int(value_slot_id)]["component_start"])
            for value_slot_id in value_slot_ids
        ),
        momentum_parameter_count=sum(
            int(momentum_slots[int(momentum_slot_id)]["component_stop"])
            - int(momentum_slots[int(momentum_slot_id)]["component_start"])
            for momentum_slot_id in momentum_slot_ids
        ),
        real_valued_inputs=tuple(int(index) for index in builder.real_valued_inputs),
    )


def _global_stage_inputs(
    *,
    parameter_symbols: Sequence[Any],
    value_symbols: Sequence[Any],
    momentum_symbols: Sequence[Any],
    value_parameter_count: int,
    momentum_parameter_count: int,
    real_valued_inputs: Sequence[int],
) -> _StageLocalInputs:
    return _StageLocalInputs(
        parameter_symbols=tuple(parameter_symbols),
        input_components=(),
        value_symbols=tuple(value_symbols),
        momentum_symbols=tuple(momentum_symbols),
        value_parameter_count=int(value_parameter_count),
        momentum_parameter_count=int(momentum_parameter_count),
        real_valued_inputs=tuple(int(index) for index in real_valued_inputs),
    )


def _parameter_builder(schema: dict[str, Any]) -> ParamBuilder:
    layout = _dict(schema["parameter_layout"])
    builder = ParamBuilder()
    builder.add_parameter_list(
        ("generic_schema_v2", "values"),
        int(layout["value_component_count"]),
        role="generic_value_storage",
    )
    builder.add_parameter_list(
        ("generic_schema_v2", "momenta"),
        int(layout["momentum_parameter_count"]),
        role="generic_momentum_storage",
        real_valued=True,
    )
    return builder


def _manifest_model(manifest: GenericProcessManifest | GenericDAG) -> Model:
    if isinstance(manifest, GenericProcessManifest):
        return manifest.model
    return AmplicolSMLeadingColorModel()


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
