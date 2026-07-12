from __future__ import annotations

import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .generic_artifact import (
    GenericProcessManifest,
    LC_SECTOR_SELECTOR_PARAMETER,
    _generic_runtime_schema_payload,
)
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
    model_parameter_count: int
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
            "model_parameter_count": self.model_parameter_count,
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
        else {
            str(record["name"]): model_parameter_symbols[
                int(record["parameter_index"])
            ]
            for record in model_parameter_records
        }
    )
    expression_model = _RuntimeParameterizedModel(
        selected_model,
        model_parameter_symbols_by_name,
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
        compiled_stages.append(
            _compile_current_stage_blueprint(
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
        )
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

    stage_payloads = []
    stage_timings: list[dict[str, object]] = []
    for stage in blueprint.stages:
        payload = stage.to_json_dict()
        payload["evaluator"] = compile_stage(stage)
        stage_timings.append(
            _stage_build_timing_record(stage.evaluator_label, payload["evaluator"])
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

    amplitude_payload = blueprint.amplitude_stage.to_json_dict()
    amplitude_payload["evaluator"] = compile_stage(blueprint.amplitude_stage)
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

    stage_local_layout = (
        blueprint.amplitude_stage.parameter_layout == "stage-local-value-momentum"
        and all(stage.parameter_layout == "stage-local-value-momentum" for stage in blueprint.stages)
    )
    total_build_s = time.perf_counter() - build_started
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

    settings = _stage_symbolica_settings(
        stage,
        blueprint,
        symbolica_settings or SymbolicaEvaluatorSettings(),
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
    blueprint: GenericStageCompilerBlueprint,
    settings: Any,
) -> Any:
    """Apply the optional recursion-stage output-chunk taper."""

    strategy = getattr(settings, "output_chunk_strategy", "uniform")
    if strategy == "auto":
        base = getattr(settings, "compiled_output_chunk_size", None)
        strategy = (
            "measured-stage"
            if getattr(settings, "backend", None) == "jit"
            and base is not None
            and int(base) <= 256
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

    try:
        position = next(
            index
            for index, current_stage in enumerate(blueprint.stages)
            if current_stage.stage_index == stage.stage_index
        )
    except StopIteration:
        return settings
    remaining = len(blueprint.stages) - position - 1
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
    interactions = [_dict(item) for item in _list(stage["interactions"])]
    input_value_slot_ids = tuple(
        int(value) for value in _list(stage["input_value_slot_ids"])
    )
    input_momentum_slot_ids = _stage_input_momentum_slot_ids(interactions)
    output_slot_ids = tuple(int(value) for value in _list(stage["output_value_slot_ids"]))
    output_slots_by_current: dict[int, list[dict[str, Any]]] = {}
    for slot_id in output_slot_ids:
        slot = value_slots[int(slot_id)]
        output_slots_by_current.setdefault(int(slot["current_id"]), []).append(slot)
    stage_model_parameter_records = (
        _current_stage_model_parameter_records(
            model,
            model_parameter_records,
            interactions=interactions,
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
            model_parameter_count=len(global_model_parameter_symbols),
            real_valued_inputs=global_real_valued_inputs,
        )
    )
    stage_model = (
        model.with_runtime_parameters(local_inputs.model_parameter_symbols)
        if isinstance(model, _RuntimeParameterizedModel)
        else _RuntimeParameterizedModel(model, local_inputs.model_parameter_symbols)
    )
    interactions_by_result: dict[int, list[dict[str, Any]]] = {}
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

    for current_id in sorted(interactions_by_result):
        current_slot = current_slots[current_id]
        dimension = int(current_slot["dimension"])
        total = tuple(0j for _ in range(dimension))
        for interaction in interactions_by_result[current_id]:
            try:
                contribution = _interaction_contribution(
                    dag,
                    stage_model,
                    interaction,
                    value_components_by_slot_id=value_components_by_slot_id,
                    momentum_components_by_slot_id=momentum_components_by_slot_id,
                    model_parameter_symbols=local_inputs.model_parameter_symbols,
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
        interaction_ids=tuple(int(interaction["interaction_id"]) for interaction in interactions),
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
            model_parameter_count=len(global_model_parameter_symbols),
            real_valued_inputs=global_real_valued_inputs,
        )
    )
    stage_model = (
        model.with_runtime_parameters(local_inputs.model_parameter_symbols)
        if isinstance(model, _RuntimeParameterizedModel)
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
    if color_weight[1] != 0.0:
        raise ValueError("complex color weights are not lowered in current stages")
    weight = color_weight[0]
    if weight == 1.0:
        return components
    return tuple(weight * component for component in components)


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


def _current_stage_model_parameter_records(
    model: Model,
    model_parameter_records: Sequence[dict[str, Any]],
    *,
    interactions: Sequence[dict[str, Any]],
    output_slots_by_current: Mapping[int, Sequence[dict[str, Any]]],
    current_slots: Mapping[int, dict[str, Any]],
) -> tuple[dict[str, Any], ...]:
    used_names = _coupling_parameter_names_used_by_records(interactions)
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
        if str(record["name"]) in used_names
    )


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
    model_parameter_symbols: dict[str, Any] = {}

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
        model_parameter_symbols[name] = symbol
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
