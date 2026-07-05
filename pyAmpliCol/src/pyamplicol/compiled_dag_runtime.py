from __future__ import annotations

import hashlib
import json
import math
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Literal, Sequence

import numpy as np

from .dag_runtime import (
    DAGEvaluationTiming,
    SharedCurrentTable,
    SymbolicaEvaluatorSettings,
    _ChunkedSymbolicaEvaluator,
    _artifact_path_for_manifest,
    _artifact_path_from_manifest,
    _build_shared_helicity_current_table,
    _compile_symbolica_outputs,
    _concatenate_complex_outputs,
    _evaluate_complex_outputs,
    _expr_dot_weyl,
    _expr_minkowski_dot,
    _expr_minkowski_dot_momentum,
    _expr_minkowski_square,
    _expr_propagate,
    _expr_vertex,
    _expr_vector_slash_terms,
    _final_state_identical_factor,
    _initial_state_average_factor,
    _load_symbolica_evaluator_artifact,
    _momentum_labels_for_current,
    _prepare_process_output_directory,
    _read_evaluator_artifact_manifest,
    _register_shared_momentum_parameters,
    _resolve_compiled_preset,
    _shared_current_stage_ids,
    _shared_momentum_current_ids,
    _symbolica_evaluator_artifact_manifest,
    _validate_z_gluon_graph,
    _validate_z_gluon_point,
    _weighted_abs2_sums_array,
    _z_gluon_graph,
)
from .lowering import _quark_vector_weyl_tensor_name, _weyl_coupling_for_chirality
from .matrix import RecursionGraph
from .model import AmplicolSMLeadingColorModel
from .native import (
    ExternalMomentum,
    HelicityContribution,
    LeadingColorZJetsNativeEvaluator,
    MatrixElementEvaluation,
    NativeEvaluationError,
    _ext_antiquark_weyl,
    _ext_gluon_cmplx,
    _ext_massive_vector,
    _ext_quark_weyl,
    _negate_momentum,
)
from .params import ParamBuilder
from .phase_space import rambo_z_gluon_point
from .symbols import symbols
from .tensor_runtime import HelicityFilter, HelicityFilterEntry

CompiledDAGLowering = Literal["spenso", "symbolic"]
CompiledDAGHelicityFilterMode = Literal["rambo", "canonical"]


@dataclass(frozen=True)
class CompiledDAGMetadata:
    process: str
    kernel: str
    lowering: str
    effective_lowering: str
    lowering_note: str | None
    gluon_count: int
    graph_current_count: int
    graph_interaction_count: int
    graph_amplitude_count: int
    shared_current_count: int
    shared_source_current_count: int
    shared_interaction_count: int
    shared_amplitude_count: int
    current_table_build_time_s: float
    helicity_filter_build_time_s: float
    alias_construction_time_s: float
    symbolica_evaluator_build_time_s: float
    symbolica_parameter_count: int
    symbolica_output_count: int
    alias_component_count: int
    opaque_alias_component_count: int
    symbolica_alias_available: bool
    spenso_body_execution_time_s: float
    spenso_body_count: int
    inline_external_wavefunctions: bool
    real_arithmetic: bool
    source_alias_component_count: int
    current_alias_component_count: int
    zero_current_component_count: int
    momentum_alias_component_count: int
    propagator_denominator_alias_count: int
    vertex_temporary_alias_count: int
    cross_check_lowering: bool
    cross_check_max_relative_difference: float | None
    batch_size: int
    symbolica_evaluator_settings: dict[str, object]
    helicity_filter: dict[str, object] | None = None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "kernel": self.kernel,
            "lowering": self.lowering,
            "effective_lowering": self.effective_lowering,
            "lowering_note": self.lowering_note,
            "gluon_count": self.gluon_count,
            "graph_current_count": self.graph_current_count,
            "graph_interaction_count": self.graph_interaction_count,
            "graph_amplitude_count": self.graph_amplitude_count,
            "shared_current_count": self.shared_current_count,
            "shared_source_current_count": self.shared_source_current_count,
            "shared_interaction_count": self.shared_interaction_count,
            "shared_amplitude_count": self.shared_amplitude_count,
            "current_table_build_time_s": self.current_table_build_time_s,
            "helicity_filter_build_time_s": self.helicity_filter_build_time_s,
            "alias_construction_time_s": self.alias_construction_time_s,
            "symbolica_evaluator_build_time_s": self.symbolica_evaluator_build_time_s,
            "symbolica_parameter_count": self.symbolica_parameter_count,
            "symbolica_output_count": self.symbolica_output_count,
            "alias_component_count": self.alias_component_count,
            "opaque_alias_component_count": self.opaque_alias_component_count,
            "symbolica_alias_available": self.symbolica_alias_available,
            "spenso_body_execution_time_s": self.spenso_body_execution_time_s,
            "spenso_body_count": self.spenso_body_count,
            "inline_external_wavefunctions": self.inline_external_wavefunctions,
            "real_arithmetic": self.real_arithmetic,
            "source_alias_component_count": self.source_alias_component_count,
            "current_alias_component_count": self.current_alias_component_count,
            "zero_current_component_count": self.zero_current_component_count,
            "momentum_alias_component_count": self.momentum_alias_component_count,
            "propagator_denominator_alias_count": (
                self.propagator_denominator_alias_count
            ),
            "vertex_temporary_alias_count": self.vertex_temporary_alias_count,
            "cross_check_lowering": self.cross_check_lowering,
            "cross_check_max_relative_difference": self.cross_check_max_relative_difference,
            "batch_size": self.batch_size,
            "symbolica_evaluator_settings": self.symbolica_evaluator_settings,
            "helicity_filter": self.helicity_filter,
        }


@dataclass(frozen=True)
class _CompiledDAGParameterLayout:
    parameter_symbols: tuple[Any, ...]
    parameter_count: int
    source_offsets: tuple[int, ...]
    source_component_offsets: tuple[tuple[int, ...], ...]
    momentum_offsets: dict[int, int]
    external_momentum_offsets: tuple[tuple[int, ...], ...]
    real_valued_inputs: tuple[int, ...]
    inline_external_wavefunctions: bool

    def to_json_dict(self) -> dict[str, object]:
        return {
            "parameter_count": self.parameter_count,
            "source_offsets": list(self.source_offsets),
            "source_component_offsets": [
                list(offsets) for offsets in self.source_component_offsets
            ],
            "momentum_offsets": {
                str(current_id): offset
                for current_id, offset in self.momentum_offsets.items()
            },
            "external_momentum_offsets": [
                list(offsets) for offsets in self.external_momentum_offsets
            ],
            "real_valued_inputs": list(self.real_valued_inputs),
            "inline_external_wavefunctions": self.inline_external_wavefunctions,
        }


@dataclass(frozen=True)
class _CompiledDAGSourceParameterSpec:
    current_id: int
    offset: int
    component_offsets: tuple[int, ...]
    dimension: int
    leg_label: int
    helicity: int
    physical_helicity: int
    chirality: int


@dataclass(frozen=True)
class _CompiledDAGSourceFillPlan:
    anti_offsets: dict[int, tuple[int, ...]]
    quark_offsets: dict[int, tuple[int, ...]]
    gluon_offsets: dict[tuple[int, int], tuple[int, ...]]
    z_offsets: dict[int, tuple[int, ...]]
    z_label: int


@dataclass(frozen=True)
class _RealPairExpression:
    real: Any
    imag: Any


class ZGluonCompiledDAGEvaluator:
    """Single multi-output Symbolica evaluator built from shared current aliases."""

    def __init__(
        self,
        process: str,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        graph: RecursionGraph | None = None,
        batch_size: int = 16,
        lowering: CompiledDAGLowering = "spenso",
        cross_check_lowering: bool = False,
        verbose_evaluator_build: bool = False,
        symbolica_evaluator_backend: str = "jit",
        symbolica_iterations: int = 1,
        symbolica_cpe_iterations: int | None = None,
        symbolica_n_cores: int = 4,
        symbolica_direct_translation: bool = True,
        symbolica_jit_direct_translation: bool | None = None,
        symbolica_jit_optimization_level: int = 3,
        symbolica_max_horner_scheme_variables: int = 500,
        symbolica_max_common_pair_cache_entries: int = 1000000,
        symbolica_max_common_pair_distance: int = 250,
        symbolica_collect_factors: bool = False,
        symbolica_compiled_preset: str = "adaptive",
        symbolica_compiled_inline_asm: str = "default",
        symbolica_compiled_optimization_level: int = 3,
        symbolica_compiled_native: bool = True,
        symbolica_compiler_path: str | None = None,
        symbolica_compiler_flags: Sequence[str] = (),
        symbolica_compiled_output_chunk_size: int | None = None,
        symbolica_compiled_chunk_compile_workers: int = 1,
        symbolica_compiled_output_dir: str | Path | None = None,
        symbolica_load_evaluator_dir: str | Path | None = None,
        symbolica_raw_sum_final_stage: bool = False,
        symbolica_real_param_sqrt_real: bool = True,
        symbolica_real_param_log_real: bool = False,
        symbolica_real_param_powf_real: bool = True,
        symbolica_real_param_real_if_args_real: bool = True,
        compiled_dag_real_arithmetic: bool | None = None,
        compiled_dag_helicity_filter: bool = True,
        compiled_dag_helicity_filter_samples: int = 10,
        compiled_dag_helicity_filter_seed: int = 12345,
        compiled_dag_helicity_filter_relative_tolerance: float = 1.0e-12,
        compiled_dag_helicity_filter_zero_tolerance: float = 1.0e-300,
        compiled_dag_helicity_filter_phase_space: CompiledDAGHelicityFilterMode = "rambo",
        compiled_dag_inline_external_wavefunctions: bool = True,
    ) -> None:
        if batch_size < 1:
            raise NativeEvaluationError("batch_size must be positive")
        if lowering not in ("spenso", "symbolic"):
            raise NativeEvaluationError(
                "compiled DAG lowering must be 'spenso' or 'symbolic'"
            )
        if cross_check_lowering and lowering == "symbolic":
            raise NativeEvaluationError(
                "compiled DAG cross-check lowering is only meaningful from the spenso route"
            )
        if compiled_dag_helicity_filter_samples < 1:
            raise NativeEvaluationError(
                "compiled DAG helicity filter sample count must be positive"
            )
        if compiled_dag_helicity_filter_phase_space not in ("rambo", "canonical"):
            raise NativeEvaluationError(
                "compiled DAG helicity filter phase-space mode must be 'rambo' or 'canonical'"
            )
        if symbolica_raw_sum_final_stage:
            raise NativeEvaluationError(
                "compiled-DAG raw-sum final stage is disabled: Symbolica "
                "JIT/compiled evaluators do not currently support the conj() "
                "operation needed for |amplitude|^2"
            )
        self.process = process
        self.model = model or AmplicolSMLeadingColorModel()
        self.graph = graph or _z_gluon_graph(process, self.model)
        self.gluon_count = _validate_z_gluon_graph(self.graph)
        self._process_pdgs = tuple(int(pdg) for pdg in self.graph.process)
        self._matrix_element_normalization = _z_gluon_matrix_element_normalization(
            self.model,
            self._process_pdgs,
            gluon_count=self.gluon_count,
        )
        self.batch_size = int(batch_size)
        self.lowering = lowering
        effective_jit_direct_translation = (
            self.gluon_count >= 3
            if symbolica_jit_direct_translation is None
            else bool(symbolica_jit_direct_translation)
        )
        self.cross_check_lowering = bool(cross_check_lowering)
        self._cross_check_compiled: _CompiledDAGMultiOutputEvaluator | None = None
        self.last_cross_check_max_relative_difference: float | None = None
        self.last_runtime_timing = DAGEvaluationTiming()
        self.loaded_evaluator_artifact_metadata: dict[str, Any] | None = None
        self.helicity_filter: HelicityFilter | None = None
        self.helicity_filter_build_time_s = 0.0

        (
            effective_compiled_inline_asm,
            effective_compiled_optimization_level,
            effective_compiled_output_chunk_size,
        ) = _resolve_compiled_dag_compiled_preset(
            symbolica_compiled_preset,
            gluon_count=self.gluon_count,
            inline_asm=symbolica_compiled_inline_asm,
            optimization_level=symbolica_compiled_optimization_level,
            output_chunk_size=symbolica_compiled_output_chunk_size,
        )
        effective_cpe_iterations = (
            2
            if symbolica_cpe_iterations is None and self.gluon_count <= 5
            else symbolica_cpe_iterations
        )
        self.symbolica_settings = SymbolicaEvaluatorSettings(
            backend=symbolica_evaluator_backend,
            iterations=symbolica_iterations,
            cpe_iterations=effective_cpe_iterations,
            n_cores=symbolica_n_cores,
            direct_translation=symbolica_direct_translation,
            jit_direct_translation=effective_jit_direct_translation,
            jit_optimization_level=symbolica_jit_optimization_level,
            max_horner_scheme_variables=symbolica_max_horner_scheme_variables,
            max_common_pair_cache_entries=symbolica_max_common_pair_cache_entries,
            max_common_pair_distance=symbolica_max_common_pair_distance,
            collect_factors=symbolica_collect_factors,
            compiled_preset=symbolica_compiled_preset,
            compiled_inline_asm=effective_compiled_inline_asm,
            compiled_optimization_level=effective_compiled_optimization_level,
            compiled_native=symbolica_compiled_native,
            compiler_path=symbolica_compiler_path,
            compiler_flags=tuple(symbolica_compiler_flags),
            compiled_output_chunk_size=effective_compiled_output_chunk_size,
            compiled_chunk_compile_workers=symbolica_compiled_chunk_compile_workers,
            compiled_output_dir=(
                None
                if symbolica_compiled_output_dir is None
                else str(symbolica_compiled_output_dir)
            ),
            real_param_sqrt_real=symbolica_real_param_sqrt_real,
            real_param_log_real=symbolica_real_param_log_real,
            real_param_powf_real=symbolica_real_param_powf_real,
            real_param_real_if_args_real=symbolica_real_param_real_if_args_real,
        )
        if compiled_dag_real_arithmetic is None:
            effective_real_arithmetic = False
        else:
            effective_real_arithmetic = bool(compiled_dag_real_arithmetic)
        if effective_real_arithmetic and (
            lowering != "symbolic"
            or not compiled_dag_inline_external_wavefunctions
            or self.symbolica_settings.backend != "jit"
            or self.symbolica_settings.compiled_output_chunk_size is not None
        ):
            raise NativeEvaluationError(
                "compiled-DAG real arithmetic currently requires symbolic "
                "lowering, inline external wavefunctions, the JIT backend, "
                "and no output chunking"
            )
        self.real_arithmetic = effective_real_arithmetic

        table_start = time.perf_counter()
        full_table = _build_shared_helicity_current_table(
            self.graph,
            self.model,
            gluon_count=self.gluon_count,
        )
        self.current_table_build_time_s = time.perf_counter() - table_start
        self.table = full_table

        if symbolica_load_evaluator_dir is not None:
            artifact_dir = Path(symbolica_load_evaluator_dir).expanduser()
            loaded_manifest = _read_evaluator_artifact_manifest(artifact_dir)
            if loaded_manifest.get("kind") != "pyamplicol-compiled-dag-evaluator":
                raise NativeEvaluationError(
                    f"unsupported compiled DAG artifact kind: {loaded_manifest.get('kind')!r}"
                )
            loaded_metadata = loaded_manifest.get("metadata")
            if isinstance(loaded_metadata, dict):
                self.loaded_evaluator_artifact_metadata = loaded_metadata
                self.helicity_filter_build_time_s = float(
                    loaded_metadata.get("helicity_filter_build_time_s", 0.0)
                )
                loaded_filter = loaded_metadata.get("helicity_filter")
                if isinstance(loaded_filter, dict):
                    self.helicity_filter = HelicityFilter.from_json_dict(loaded_filter)
                    self.table = _apply_helicity_filter_to_shared_table(
                        full_table,
                        self.helicity_filter,
                    )
            compiled_manifest = loaded_manifest.get("compiled")
            if not isinstance(compiled_manifest, dict):
                raise NativeEvaluationError("compiled DAG artifact is missing evaluator data")
            self.compiled = _CompiledDAGMultiOutputEvaluator.from_artifact(
                self.table,
                artifact_dir,
                manifest=compiled_manifest,
                symbolica_settings=self.symbolica_settings,
            )
            self.alias_construction_time_s = 0.0
            self.symbolica_evaluator_build_time_s = 0.0
            self.effective_lowering = str(
                self.loaded_evaluator_artifact_metadata.get("effective_lowering", lowering)
                if self.loaded_evaluator_artifact_metadata is not None
                else lowering
            )
            self.lowering_note = (
                None
                if self.loaded_evaluator_artifact_metadata is None
                else self.loaded_evaluator_artifact_metadata.get("lowering_note")
            )
        else:
            if compiled_dag_helicity_filter:
                filter_start = time.perf_counter()
                self.helicity_filter = build_compiled_dag_helicity_filter(
                    self.process,
                    full_table,
                    self.model,
                    gluon_count=self.gluon_count,
                    sample_count=compiled_dag_helicity_filter_samples,
                    seed=compiled_dag_helicity_filter_seed,
                    relative_tolerance=compiled_dag_helicity_filter_relative_tolerance,
                    zero_tolerance=compiled_dag_helicity_filter_zero_tolerance,
                    phase_space_mode=compiled_dag_helicity_filter_phase_space,
                    progress_enabled=verbose_evaluator_build,
                )
                self.helicity_filter_build_time_s = time.perf_counter() - filter_start
                self.table = _apply_helicity_filter_to_shared_table(
                    full_table,
                    self.helicity_filter,
                )
            compiled_start = time.perf_counter()
            self.compiled = _CompiledDAGMultiOutputEvaluator(
                self.table,
                model=self.model,
                lowering=lowering,
                inline_external_wavefunctions=compiled_dag_inline_external_wavefunctions,
                real_arithmetic=self.real_arithmetic,
                verbose_evaluator_build=verbose_evaluator_build,
                symbolica_settings=self.symbolica_settings,
            )
            self.alias_construction_time_s = self.compiled.alias_construction_time_s
            self.symbolica_evaluator_build_time_s = time.perf_counter() - compiled_start
            self.effective_lowering = self.compiled.effective_lowering
            self.lowering_note = self.compiled.lowering_note

        if self.cross_check_lowering:
            self._cross_check_compiled = _CompiledDAGMultiOutputEvaluator(
                self.table,
                model=self.model,
                lowering="symbolic",
                inline_external_wavefunctions=compiled_dag_inline_external_wavefunctions,
                real_arithmetic=self.real_arithmetic,
                verbose_evaluator_build=verbose_evaluator_build,
                symbolica_settings=self.symbolica_settings,
            )

    def materialize_compiled_evaluators(self) -> None:
        materialize = getattr(self.compiled, "materialize", None)
        if callable(materialize):
            materialize()

    def save_evaluator_artifact(self, output_dir: str | Path) -> Path:
        self.materialize_compiled_evaluators()
        output_path = Path(output_dir).expanduser()
        _prepare_process_output_directory(output_path)
        artifact = {
            "schema_version": 1,
            "kind": "pyamplicol-compiled-dag-evaluator",
            "process": self.process,
            "gluon_count": self.gluon_count,
            "symbolica": _symbolica_provenance(),
            "metadata": self.metadata.to_json_dict(),
            "parameter_layout": self.compiled.layout.to_json_dict(),
            "current_metadata": _compiled_dag_current_metadata(self.table),
            "root_metadata": _compiled_dag_root_metadata(self.table),
            "compiled": self.compiled.artifact_manifest(output_path),
        }
        artifact["artifact_fingerprint"] = _compiled_dag_artifact_fingerprint(
            artifact
        )
        manifest_path = output_path / "manifest.json"
        manifest_path.write_text(
            json.dumps(artifact, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        return manifest_path

    @property
    def metadata(self) -> CompiledDAGMetadata:
        settings = self.symbolica_settings.to_json_dict()
        if self.loaded_evaluator_artifact_metadata is not None:
            loaded_settings = self.loaded_evaluator_artifact_metadata.get(
                "symbolica_evaluator_settings"
            )
            if isinstance(loaded_settings, dict):
                settings = dict(loaded_settings)
                settings["loaded_from_artifact"] = True
        return CompiledDAGMetadata(
            process=self.process,
            kernel="symbolica-compiled-shared-current-alias-dag",
            lowering=self.lowering,
            effective_lowering=self.effective_lowering,
            lowering_note=self.lowering_note,
            gluon_count=self.gluon_count,
            graph_current_count=len(self.graph.currents),
            graph_interaction_count=len(self.graph.interactions),
            graph_amplitude_count=len(self.graph.amplitudes),
            shared_current_count=len(self.table.currents),
            shared_source_current_count=len(self.table.sources),
            shared_interaction_count=len(self.table.interactions),
            shared_amplitude_count=len(self.table.amplitudes),
            current_table_build_time_s=self.current_table_build_time_s,
            helicity_filter_build_time_s=self.helicity_filter_build_time_s,
            alias_construction_time_s=self.alias_construction_time_s,
            symbolica_evaluator_build_time_s=self.symbolica_evaluator_build_time_s,
            symbolica_parameter_count=self.compiled.parameter_count,
            symbolica_output_count=self.compiled.output_length,
            alias_component_count=self.compiled.alias_component_count,
            opaque_alias_component_count=self.compiled.opaque_alias_component_count,
            symbolica_alias_available=self.compiled.symbolica_alias_available,
            spenso_body_execution_time_s=self.compiled.spenso_body_execution_time_s,
            spenso_body_count=self.compiled.spenso_body_count,
            inline_external_wavefunctions=(
                self.compiled.inline_external_wavefunctions
            ),
            real_arithmetic=self.compiled.real_arithmetic,
            source_alias_component_count=self.compiled.source_alias_component_count,
            current_alias_component_count=self.compiled.current_alias_component_count,
            zero_current_component_count=(
                self.compiled.zero_current_component_count
            ),
            momentum_alias_component_count=self.compiled.momentum_alias_component_count,
            propagator_denominator_alias_count=(
                self.compiled.propagator_denominator_alias_count
            ),
            vertex_temporary_alias_count=self.compiled.vertex_temporary_alias_count,
            cross_check_lowering=self.cross_check_lowering,
            cross_check_max_relative_difference=(
                self.last_cross_check_max_relative_difference
            ),
            batch_size=self.batch_size,
            symbolica_evaluator_settings=settings,
            helicity_filter=(
                None
                if self.helicity_filter is None
                else self.helicity_filter.to_json_dict()
            ),
        )

    def evaluate(
        self,
        particles: Sequence[ExternalMomentum],
    ) -> MatrixElementEvaluation:
        point = _validate_z_gluon_point(particles, gluon_count=self.gluon_count)
        return self.evaluate_many((point,))[0]

    def evaluate_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> tuple[MatrixElementEvaluation, ...]:
        evaluations: list[MatrixElementEvaluation] = []
        timing = DAGEvaluationTiming()
        for start in range(0, len(points), self.batch_size):
            batch = tuple(
                _validate_z_gluon_point(particles, gluon_count=self.gluon_count)
                for particles in points[start : start + self.batch_size]
            )
            amplitude_rows = self.compiled.evaluate_amplitude_rows(
                batch,
                self.model,
                gluon_count=self.gluon_count,
            )
            self._cross_check_amplitude_rows(batch, amplitude_rows)
            reduction_start = time.perf_counter()
            evaluations.extend(
                self._evaluation_from_amplitudes(point, amplitudes)
                for point, amplitudes in zip(batch, amplitude_rows, strict=True)
            )
            timing += self.compiled.last_timing.with_result_reduction(
                time.perf_counter() - reduction_start
            )
        self.last_runtime_timing = timing
        return tuple(evaluations)

    def evaluate_matrix_elements_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> tuple[float, ...]:
        return tuple(float(value) for value in self.evaluate_matrix_elements_array(points))

    def evaluate_matrix_elements_array(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
        *,
        validate: bool = True,
    ) -> np.ndarray:
        matrix_element_chunks: list[np.ndarray] = []
        timing = DAGEvaluationTiming()
        for start in range(0, len(points), self.batch_size):
            batch = _compiled_dag_point_batch(
                points[start : start + self.batch_size],
                gluon_count=self.gluon_count,
                validate=validate,
            )
            raw_sums = self.compiled.evaluate_raw_sum_array(
                batch,
                self.model,
                gluon_count=self.gluon_count,
            )
            if self._cross_check_compiled is not None:
                self._cross_check_raw_sums(
                    batch,
                    tuple(float(value) for value in raw_sums),
                )
            reduction_start = time.perf_counter()
            matrix_element_chunks.append(raw_sums * self._matrix_element_normalization)
            timing += self.compiled.last_timing.with_result_reduction(
                time.perf_counter() - reduction_start
            )
        if not matrix_element_chunks:
            self.last_runtime_timing = timing
            return np.empty(0, dtype=np.float64)
        reduction_start = time.perf_counter()
        matrix_elements = (
            matrix_element_chunks[0]
            if len(matrix_element_chunks) == 1
            else np.concatenate(matrix_element_chunks)
        )
        timing += DAGEvaluationTiming(
            result_reduction_time_s=time.perf_counter() - reduction_start
        )
        self.last_runtime_timing = timing
        return matrix_elements

    def _evaluation_from_amplitudes(
        self,
        point: tuple[ExternalMomentum, ...],
        amplitudes: tuple[complex, ...],
    ) -> MatrixElementEvaluation:
        raw_sum = 0.0
        helicity_contributions: list[HelicityContribution] = []
        for amplitude_record, amplitude in zip(
            self.table.amplitudes,
            amplitudes,
            strict=True,
        ):
            squared = float((amplitude * amplitude.conjugate()).real)
            weighted_squared = squared * amplitude_record.multiplicity
            raw_sum += weighted_squared
            helicity_contributions.append(
                HelicityContribution(
                    helicities=amplitude_record.helicities,
                    amplitude=amplitude,
                    squared=weighted_squared,
                )
            )
        pdgs = tuple(particle.pdg for particle in point)
        return MatrixElementEvaluation(
            process=self.process,
            particles=point,
            matrix_element=self._matrix_element_from_raw_sum(point, raw_sum),
            raw_helicity_sum=raw_sum,
            color_factor=self.model.leading_color_factor(pdgs),
            average_factor=_initial_state_average_factor(pdgs[:2]),
            coupling_factor=(
                (4.0 * math.pi * self.model.alpha_s_me_check) ** self.gluon_count
                * (2.0 * 4.0 * math.pi * self.model.alpha_ew)
            ),
            helicity_contributions=tuple(helicity_contributions),
            identical_factor=_final_state_identical_factor(
                particle.pdg for particle in point[2:]
            ),
        )

    def _matrix_element_from_raw_sum(
        self,
        point: tuple[ExternalMomentum, ...],
        raw_sum: float,
    ) -> float:
        pdgs = tuple(particle.pdg for particle in point)
        if pdgs == self._process_pdgs:
            return raw_sum * self._matrix_element_normalization
        return raw_sum * _z_gluon_matrix_element_normalization(
            self.model,
            pdgs,
            gluon_count=self.gluon_count,
        )

    def _cross_check_amplitude_rows(
        self,
        points: tuple[tuple[ExternalMomentum, ...], ...],
        amplitude_rows: tuple[tuple[complex, ...], ...],
    ) -> None:
        if self._cross_check_compiled is None:
            return
        reference_rows = self._cross_check_compiled.evaluate_amplitude_rows(
            points,
            self.model,
            gluon_count=self.gluon_count,
        )
        max_relative = _max_relative_complex_rows_difference(
            amplitude_rows,
            reference_rows,
        )
        self.last_cross_check_max_relative_difference = max_relative
        if max_relative > 1.0e-8:
            raise NativeEvaluationError(
                "compiled DAG Spenso/symbolic lowering cross-check failed: "
                f"max relative difference {max_relative:.3e}"
            )

    def _cross_check_raw_sums(
        self,
        points: tuple[tuple[ExternalMomentum, ...], ...],
        raw_sums: tuple[float, ...],
    ) -> None:
        if self._cross_check_compiled is None:
            return
        reference_raw_sums = self._cross_check_compiled.evaluate_raw_sum_rows(
            points,
            self.model,
            gluon_count=self.gluon_count,
        )
        max_relative = _max_relative_real_rows_difference(raw_sums, reference_raw_sums)
        self.last_cross_check_max_relative_difference = max_relative
        if max_relative > 1.0e-8:
            raise NativeEvaluationError(
                "compiled DAG Spenso/symbolic raw-sum cross-check failed: "
                f"max relative difference {max_relative:.3e}"
            )


class _CompiledDAGMultiOutputEvaluator:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        model: AmplicolSMLeadingColorModel,
        lowering: CompiledDAGLowering,
        inline_external_wavefunctions: bool,
        real_arithmetic: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> None:
        self.table = table
        self.inline_external_wavefunctions = bool(inline_external_wavefunctions)
        self.real_arithmetic = bool(real_arithmetic)
        self.layout = _build_compiled_dag_parameter_layout(
            table,
            inline_external_wavefunctions=self.inline_external_wavefunctions,
        )
        self._init_common(table)
        self.lowering = lowering
        self.effective_lowering = lowering
        self.lowering_note = None
        self.symbolica_alias_available = _symbolica_evaluator_alias_available()
        self.progress_enabled = bool(verbose_evaluator_build)
        self.alias_definitions: tuple[tuple[Any, Any], ...] = ()
        self.spenso_body_execution_time_s = 0.0
        self.spenso_body_count = 0
        self.source_alias_component_count = 0
        self.current_alias_component_count = 0
        self.zero_current_component_count = 0
        self.momentum_alias_component_count = 0
        self.propagator_denominator_alias_count = 0
        self.vertex_temporary_alias_count = 0
        self.alias_vertex_temporaries = False
        self.alias_momentum_components = False
        self.alias_single_use_currents = True
        self.alias_propagator_denominators = True
        self.prune_zero_current_components = False
        alias_start = time.perf_counter()
        self.outputs = self._build_alias_outputs(table, model=model)
        self.alias_construction_time_s = time.perf_counter() - alias_start
        self.output_length = len(table.amplitudes)
        self.evaluator_output_length = len(self.outputs)
        self.alias_component_count = (
            self.source_alias_component_count
            + self.current_alias_component_count
            + self.momentum_alias_component_count
            + self.propagator_denominator_alias_count
            + self.vertex_temporary_alias_count
        )
        self.opaque_alias_component_count = (
            self.alias_component_count if self.symbolica_alias_available else 0
        )
        self.raw_sum_weights = np.asarray(
            [amplitude.multiplicity for amplitude in table.amplitudes],
            dtype=np.float64,
        )
        self.evaluator = _compile_symbolica_outputs(
            self.outputs,
            list(self.layout.parameter_symbols),
            merge_evaluators_strategy=False,
            verbose_evaluator_build=verbose_evaluator_build,
            aliases=self.alias_definitions,
            real_params=self.layout.real_valued_inputs,
            symbolica_settings=symbolica_settings,
            label="compiled_dag_amplitudes",
        )
        self.parameter_count = self.layout.parameter_count
        self.last_timing = DAGEvaluationTiming()

    @classmethod
    def from_artifact(
        cls,
        table: SharedCurrentTable,
        artifact_dir: Path,
        *,
        manifest: dict[str, Any],
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> "_CompiledDAGMultiOutputEvaluator":
        instance = cls.__new__(cls)
        instance.table = table
        instance.inline_external_wavefunctions = bool(
            manifest.get("inline_external_wavefunctions", False)
        )
        instance.real_arithmetic = bool(manifest.get("real_arithmetic", False))
        instance.layout = _build_compiled_dag_parameter_layout(
            table,
            inline_external_wavefunctions=instance.inline_external_wavefunctions,
        )
        instance._init_common(table)
        instance.lowering = str(manifest.get("lowering", "spenso"))
        instance.effective_lowering = str(manifest.get("effective_lowering", "symbolic"))
        lowering_note = manifest.get("lowering_note")
        instance.lowering_note = lowering_note if isinstance(lowering_note, str) else None
        instance.outputs = ()
        instance.output_length = int(manifest["output_length"])
        instance.evaluator_output_length = int(
            manifest.get("evaluator_output_length", instance.output_length)
        )
        instance.alias_component_count = int(manifest.get("alias_component_count", 0))
        instance.opaque_alias_component_count = int(
            manifest.get("opaque_alias_component_count", 0)
        )
        instance.source_alias_component_count = int(
            manifest.get("source_alias_component_count", 0)
        )
        instance.current_alias_component_count = int(
            manifest.get("current_alias_component_count", 0)
        )
        instance.zero_current_component_count = int(
            manifest.get("zero_current_component_count", 0)
        )
        instance.momentum_alias_component_count = int(
            manifest.get("momentum_alias_component_count", 0)
        )
        instance.propagator_denominator_alias_count = int(
            manifest.get("propagator_denominator_alias_count", 0)
        )
        instance.vertex_temporary_alias_count = int(
            manifest.get("vertex_temporary_alias_count", 0)
        )
        instance.alias_vertex_temporaries = False
        instance.alias_momentum_components = False
        instance.alias_single_use_currents = True
        instance.symbolica_alias_available = bool(
            manifest.get("symbolica_alias_available", False)
        )
        instance.spenso_body_execution_time_s = float(
            manifest.get("spenso_body_execution_time_s", 0.0)
        )
        instance.spenso_body_count = int(manifest.get("spenso_body_count", 0))
        instance.raw_sum_weights = np.asarray(
            manifest["raw_sum_weights"],
            dtype=np.float64,
        )
        instance.evaluator = _load_compiled_dag_evaluator_artifact(
            manifest["evaluator"],
            artifact_dir,
        )
        instance.parameter_count = instance.layout.parameter_count
        instance.alias_construction_time_s = 0.0
        instance.last_timing = DAGEvaluationTiming()
        _ = symbolica_settings
        if len(table.amplitudes) != instance.output_length:
            raise NativeEvaluationError(
                "compiled DAG artifact amplitude count does not match process"
            )
        return instance

    def materialize(self) -> None:
        return None

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        return {
            "kind": "compiled-dag-multi-output",
            "lowering": self.lowering,
            "effective_lowering": self.effective_lowering,
            "lowering_note": self.lowering_note,
            "parameter_count": self.parameter_count,
            "output_length": self.output_length,
            "evaluator_output_length": self.evaluator_output_length,
            "alias_component_count": self.alias_component_count,
            "opaque_alias_component_count": self.opaque_alias_component_count,
            "source_alias_component_count": self.source_alias_component_count,
            "current_alias_component_count": self.current_alias_component_count,
            "zero_current_component_count": self.zero_current_component_count,
            "momentum_alias_component_count": self.momentum_alias_component_count,
            "propagator_denominator_alias_count": (
                self.propagator_denominator_alias_count
            ),
            "vertex_temporary_alias_count": self.vertex_temporary_alias_count,
            "symbolica_alias_available": self.symbolica_alias_available,
            "spenso_body_execution_time_s": self.spenso_body_execution_time_s,
            "spenso_body_count": self.spenso_body_count,
            "inline_external_wavefunctions": self.inline_external_wavefunctions,
            "real_arithmetic": self.real_arithmetic,
            "raw_sum_weights": self.raw_sum_weights.tolist(),
            "evaluator": _compiled_dag_evaluator_artifact_manifest(
                self.evaluator,
                artifact_dir,
                label="compiled_dag_amplitudes",
            ),
        }

    def evaluate_amplitude_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[tuple[complex, ...], ...]:
        if not points:
            return ()
        parameter_rows, timing = self._parameter_rows(points, model, gluon_count=gluon_count)
        start = time.perf_counter()
        if self.real_arithmetic:
            evaluated = _evaluate_real_outputs(self.evaluator, parameter_rows)
        else:
            evaluated = _evaluate_complex_outputs(self.evaluator, parameter_rows)
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        if self.real_arithmetic:
            rows = tuple(
                tuple(
                    complex(row[2 * index], row[2 * index + 1])
                    for index in range(self.output_length)
                )
                for row in evaluated
            )
        else:
            evaluated_array = _concatenate_complex_outputs(evaluated)
            rows = tuple(
                tuple(row[: self.output_length].tolist())
                for row in evaluated_array
            )
        output_transfer_time_s = time.perf_counter() - transfer_start
        self.last_timing = timing + DAGEvaluationTiming(
            evaluator_time_s=evaluator_time_s,
            output_transfer_time_s=output_transfer_time_s,
        )
        return rows

    def evaluate_raw_sum_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[float, ...]:
        return tuple(
            float(value)
            for value in self.evaluate_raw_sum_array(
                points,
                model,
                gluon_count=gluon_count,
            )
        )

    def evaluate_raw_sum_array(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> np.ndarray:
        if not points:
            return np.empty(0, dtype=np.float64)
        parameter_rows, timing = self._parameter_rows(points, model, gluon_count=gluon_count)
        start = time.perf_counter()
        if self.real_arithmetic:
            evaluated = _evaluate_real_outputs(self.evaluator, parameter_rows)
        else:
            evaluated = _evaluate_complex_outputs(self.evaluator, parameter_rows)
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        if self.real_arithmetic:
            raw_sums = _weighted_real_pair_abs2_sums_array(
                evaluated,
                self.raw_sum_weights,
                self.output_length,
            )
        else:
            raw_sums = _weighted_abs2_sums_array(
                evaluated,
                self.raw_sum_weights,
                self.output_length,
            )
        output_transfer_time_s = time.perf_counter() - transfer_start
        self.last_timing = timing + DAGEvaluationTiming(
            evaluator_time_s=evaluator_time_s,
            output_transfer_time_s=output_transfer_time_s,
        )
        return raw_sums

    def _build_alias_outputs(
        self,
        table: SharedCurrentTable,
        *,
        model: AmplicolSMLeadingColorModel,
    ) -> tuple[Any, ...]:
        if self.real_arithmetic:
            return self._build_real_pair_alias_outputs(table, model=model)

        from symbolica import Expression

        parameter_symbols = self.layout.parameter_symbols
        zero = Expression.num(0)
        alias_definitions: list[tuple[Any, Any]] = []
        if self.layout.inline_external_wavefunctions:
            external_momenta = _compiled_dag_external_momentum_expressions(
                self.layout,
            )
            current_expressions = _compiled_dag_inline_source_expressions(
                table,
                external_momenta,
                model,
            )
            momentum_expressions = _compiled_dag_inline_momentum_expressions(
                table,
                external_momenta,
            )
            if self.symbolica_alias_available:
                self.source_alias_component_count = _alias_inline_source_expressions(
                    table,
                    current_expressions,
                    alias_definitions,
                )
                if self.alias_momentum_components:
                    self.momentum_alias_component_count = (
                        _alias_inline_momentum_expressions(
                            momentum_expressions,
                            alias_definitions,
                        )
                    )
        else:
            current_expressions = {}
            for current in table.currents:
                if not current.is_source:
                    continue
                component_offsets = self.layout.source_component_offsets[current.id]
                current_expressions[current.id] = tuple(
                    (
                        zero
                        if offset < 0
                        else parameter_symbols[offset]
                    )
                    for offset in component_offsets
                )

            momentum_expressions = {}
            for current_id, start in self.layout.momentum_offsets.items():
                momentum_expressions[current_id] = tuple(
                    parameter_symbols[index] for index in range(start, start + 4)
                )

        spenso_lowerer = (
            _SpensoCompiledDAGCurrentBodyLowerer(model)
            if self.lowering == "spenso"
            else None
        )
        current_use_counts = _compiled_dag_current_use_counts(table)
        current_work = tuple(
            (stage_index, current_id)
            for stage_index, current_ids in enumerate(_shared_current_stage_ids(table), 1)
            for current_id in current_ids
        )
        for _stage_index, current_id in _compiled_dag_progress(
            current_work,
            enabled=self.progress_enabled,
            label="compiled DAG aliases",
            metadata=f"curr={len(current_work):<5} out={len(table.amplitudes):<5}",
        ):
            current = table.currents[current_id]
            if spenso_lowerer is None:
                total = tuple(0j for _ in range(current.dimension))
                for interaction_id in table.interactions_by_result[current_id]:
                    interaction = table.interactions[interaction_id]
                    if (
                        self.symbolica_alias_available
                        and self.alias_vertex_temporaries
                    ):
                        contribution, alias_count = _expr_vertex_with_aliases(
                            interaction,
                            current_expressions[interaction.left_id],
                            current_expressions[interaction.right_id],
                            momentum_expressions,
                            alias_definitions,
                        )
                        self.vertex_temporary_alias_count += alias_count
                    else:
                        contribution = _expr_vertex(
                            interaction,
                            current_expressions[interaction.left_id],
                            current_expressions[interaction.right_id],
                            momentum_expressions,
                        )
                    total = tuple(
                        left + right
                        for left, right in zip(total, contribution, strict=True)
                    )
                if current.needs_propagator:
                    if (
                        self.symbolica_alias_available
                        and self.alias_propagator_denominators
                    ):
                        denominator_alias = _propagator_denominator_alias(current.id)
                        alias_definitions.append(
                            (
                                denominator_alias,
                                _expr_minkowski_square(
                                    momentum_expressions[current.id]
                                ),
                            )
                        )
                        self.propagator_denominator_alias_count += 1
                        total = _expr_propagate_with_denominator(
                            current.key,
                            total,
                            momentum_expressions[current.id],
                            denominator_alias,
                        )
                    else:
                        total = _expr_propagate(
                            current.key,
                            total,
                            momentum_expressions[current.id],
                        )
            else:
                total = spenso_lowerer.lower_current_body(
                    table,
                    current_id=current_id,
                    current_expressions=current_expressions,
                    momentum_expressions=momentum_expressions,
                )
            aliases = tuple(
                _current_component_alias(
                    current_id,
                    component_index,
                )
                for component_index in range(current.dimension)
            )
            if self.symbolica_alias_available:
                should_alias_current = (
                    spenso_lowerer is not None
                    or self.alias_single_use_currents
                    or current_use_counts[current_id] > 1
                )
            else:
                should_alias_current = False
            if should_alias_current:
                aliased_total: list[Any] = []
                for alias, body in zip(aliases, total, strict=True):
                    if (
                        self.prune_zero_current_components
                        and _is_symbolica_zero(body)
                    ):
                        aliased_total.append(zero)
                        self.zero_current_component_count += 1
                        continue
                    alias_definitions.append((alias, body))
                    aliased_total.append(alias)
                    self.current_alias_component_count += 1
                current_expressions[current_id] = tuple(aliased_total)
            else:
                current_expressions[current_id] = total

        if spenso_lowerer is not None:
            self.spenso_body_execution_time_s = spenso_lowerer.execution_time_s
            self.spenso_body_count = spenso_lowerer.body_count
        self.alias_definitions = tuple(alias_definitions)
        return tuple(
            _expr_dot_weyl(
                current_expressions[amplitude.left_id],
                current_expressions[amplitude.right_id],
            )
            for amplitude in table.amplitudes
        )

    def _build_real_pair_alias_outputs(
        self,
        table: SharedCurrentTable,
        *,
        model: AmplicolSMLeadingColorModel,
    ) -> tuple[Any, ...]:
        from symbolica import Expression

        external_momenta = _compiled_dag_external_momentum_expressions(
            self.layout,
        )
        current_expressions = _compiled_dag_inline_source_pair_expressions(
            table,
            external_momenta,
            model,
        )
        momentum_expressions = _compiled_dag_inline_momentum_expressions(
            table,
            external_momenta,
        )
        alias_definitions: list[tuple[Any, Any]] = []
        self.source_alias_component_count = _alias_inline_source_pair_expressions(
            table,
            current_expressions,
            alias_definitions,
        )
        zero = Expression.num(0)
        current_work = tuple(
            (stage_index, current_id)
            for stage_index, current_ids in enumerate(_shared_current_stage_ids(table), 1)
            for current_id in current_ids
        )
        for _stage_index, current_id in _compiled_dag_progress(
            current_work,
            enabled=self.progress_enabled,
            label="compiled DAG real aliases",
            metadata=f"curr={len(current_work):<5} out={len(table.amplitudes):<5}",
        ):
            current = table.currents[current_id]
            total = tuple(_rp_zero(zero) for _ in range(current.dimension))
            for interaction_id in table.interactions_by_result[current_id]:
                interaction = table.interactions[interaction_id]
                contribution = _expr_vertex_pair(
                    interaction,
                    current_expressions[interaction.left_id],
                    current_expressions[interaction.right_id],
                    momentum_expressions,
                )
                total = tuple(
                    _rp_add(left, right)
                    for left, right in zip(total, contribution, strict=True)
                )
            if current.needs_propagator:
                if self.alias_propagator_denominators:
                    denominator_alias = _propagator_denominator_alias(current.id)
                    alias_definitions.append(
                        (
                            denominator_alias,
                            _expr_minkowski_square(momentum_expressions[current.id]),
                        )
                    )
                    self.propagator_denominator_alias_count += 1
                    total = _expr_propagate_pair_with_denominator(
                        current.key,
                        total,
                        momentum_expressions[current.id],
                        denominator_alias,
                    )
                else:
                    total = _expr_propagate_pair_with_denominator(
                        current.key,
                        total,
                        momentum_expressions[current.id],
                        _expr_minkowski_square(momentum_expressions[current.id]),
                    )
            aliases = tuple(
                _RealPairExpression(
                    _current_component_real_alias(current_id, component_index),
                    _current_component_imag_alias(current_id, component_index),
                )
                for component_index in range(current.dimension)
            )
            for alias, body in zip(aliases, total, strict=True):
                alias_definitions.append((alias.real, body.real))
                alias_definitions.append((alias.imag, body.imag))
            current_expressions[current_id] = aliases

        self.alias_definitions = tuple(alias_definitions)
        outputs: list[Any] = []
        for amplitude in table.amplitudes:
            value = _expr_dot_weyl_pair(
                current_expressions[amplitude.left_id],
                current_expressions[amplitude.right_id],
            )
            outputs.extend((value.real, value.imag))
        return tuple(outputs)

    def _init_common(self, table: SharedCurrentTable) -> None:
        self._source_parameter_specs = _compiled_dag_source_parameter_specs(
            table,
            self.layout,
        )
        self._source_fill_plan = _compiled_dag_source_fill_plan(
            self._source_parameter_specs
        )
        self._momentum_offsets_and_labels = _compiled_dag_momentum_offsets_and_labels(
            table,
            self.layout,
        )
        self._external_label_count = max(
            label
            for current in table.currents
            for label in current.external_labels
        )
        self._external_momentum_signs = np.ones(
            (self._external_label_count, 1),
            dtype=np.float64,
        )
        self._external_momentum_signs[:2, 0] = -1.0
        self._momentum_label_matrix = np.zeros(
            (len(self._momentum_offsets_and_labels), self._external_label_count),
            dtype=np.float64,
        )
        self._momentum_flat_columns = np.asarray(
            [
                column
                for offset, _ in self._momentum_offsets_and_labels
                for column in range(offset, offset + 4)
            ],
            dtype=np.intp,
        )
        for row, (_offset, labels) in enumerate(self._momentum_offsets_and_labels):
            for label in labels:
                self._momentum_label_matrix[row, label - 1] = 1.0

    def _parameter_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[np.ndarray, DAGEvaluationTiming]:
        pack_start = time.perf_counter()
        rows = np.empty(
            (len(points), self.layout.parameter_count),
            dtype=(np.float64 if self.real_arithmetic else np.complex128),
        )
        external_momenta = _external_momentum_array(
            points,
            self._external_label_count,
        )
        if self.layout.inline_external_wavefunctions:
            rows[:, :] = external_momenta.reshape(len(points), -1)
            parameter_pack_time_s = time.perf_counter() - pack_start
            return (
                rows,
                DAGEvaluationTiming(
                    parameter_pack_time_s=parameter_pack_time_s,
                ),
            )
        parameter_pack_time_s = time.perf_counter() - pack_start
        source_start = time.perf_counter()
        if not _fill_compiled_dag_source_parameters_batch_fast(
            rows,
            self._source_fill_plan,
            external_momenta,
            model,
            gluon_count=gluon_count,
        ):
            _fill_compiled_dag_source_parameters_batch(
                rows,
                self._source_parameter_specs,
                external_momenta,
                model,
                gluon_count=gluon_count,
            )
        source_fill_time_s = time.perf_counter() - source_start

        momentum_start = time.perf_counter()
        if self._momentum_offsets_and_labels:
            signed_momenta = external_momenta * self._external_momentum_signs
            momentum_values = np.einsum(
                "ml,blc->bmc",
                self._momentum_label_matrix,
                signed_momenta,
                optimize=True,
            )
            flat_momentum_values = momentum_values.reshape(len(points), -1)
            rows[:, self._momentum_flat_columns] = flat_momentum_values
        momentum_setup_time_s = time.perf_counter() - momentum_start
        return (
            rows,
            DAGEvaluationTiming(
                source_fill_time_s=source_fill_time_s,
                momentum_setup_time_s=momentum_setup_time_s,
                parameter_pack_time_s=parameter_pack_time_s,
            ),
        )


class _SpensoCompiledDAGCurrentBodyLowerer:
    def __init__(self, model: AmplicolSMLeadingColorModel) -> None:
        from symbolica.community.spenso import Representation, TensorName

        self.model = model
        self._mink = Representation.mink(4)
        self._aux6 = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
        self._weyl = Representation("pyamplicol::weyl_spinor", 2)
        self._tensor_name = TensorName
        self.execution_time_s = 0.0
        self.body_count = 0

    def lower_current_body(
        self,
        table: SharedCurrentTable,
        *,
        current_id: int,
        current_expressions: dict[int, tuple[Any, ...]],
        momentum_expressions: dict[int, tuple[Any, ...]],
    ) -> tuple[Any, ...]:
        from symbolica import Expression
        from symbolica.community.spenso import LibraryTensor, TensorNetwork

        current = table.currents[current_id]
        library = self.model.build_tensor_library()
        registered_currents: set[int] = set()
        registered_momenta: set[int] = set()

        def register_current_tensor(parent_id: int) -> None:
            if parent_id in registered_currents:
                return
            parent = table.currents[parent_id]
            library.register(
                LibraryTensor.dense(
                    self._tensor_name(self._current_tensor_name(parent_id))(
                        self._representation_for_current(parent)
                    ),
                    current_expressions[parent_id],
                )
            )
            registered_currents.add(parent_id)

        def register_momentum_tensor(momentum_current_id: int) -> None:
            if momentum_current_id in registered_momenta:
                return
            library.register(
                LibraryTensor.dense(
                    self._tensor_name(self._momentum_tensor_name(momentum_current_id))(
                        self._mink
                    ),
                    momentum_expressions[momentum_current_id],
                )
            )
            registered_momenta.add(momentum_current_id)

        output_slots = self._slots_for_current(current, f"j{current_id}_out")
        result_slots = (
            self._slots_for_current(current, f"j{current_id}_preprop")
            if current.needs_propagator
            else output_slots
        )
        total = Expression.num(0)
        for interaction_id in table.interactions_by_result[current_id]:
            interaction = table.interactions[interaction_id]
            register_current_tensor(interaction.left_id)
            register_current_tensor(interaction.right_id)
            if int(interaction.vertex_kind) == 0:
                register_momentum_tensor(interaction.left_id)
                register_momentum_tensor(interaction.right_id)
            left_slots = self._slots_for_current(
                table.currents[interaction.left_id],
                f"j{current_id}_i{interaction_id}_left",
            )
            right_slots = self._slots_for_current(
                table.currents[interaction.right_id],
                f"j{current_id}_i{interaction_id}_right",
            )
            total = total + (
                self._vertex_expression(
                    interaction,
                    left_slots=left_slots,
                    right_slots=right_slots,
                    output_slots=result_slots,
                )
                * self._tensor_name(self._current_tensor_name(interaction.left_id))(
                    *left_slots
                ).to_expression()
                * self._tensor_name(self._current_tensor_name(interaction.right_id))(
                    *right_slots
                ).to_expression()
            )

        if current.needs_propagator:
            register_momentum_tensor(current.id)
            library.register(
                LibraryTensor.dense(
                    self._tensor_name(self._propagator_tensor_name(current.id))(
                        self._representation_for_current(current),
                        self._representation_for_current(current),
                    ),
                    self._propagator_tensor_data(
                        current,
                        momentum_expressions[current.id],
                    ),
                )
            )
            total = (
                self._tensor_name(self._propagator_tensor_name(current.id))(
                    result_slots[0],
                    output_slots[0],
                ).to_expression()
                * total
            )

        start = time.perf_counter()
        network = TensorNetwork(total, library)
        network.execute(library=library)
        result = network.result_tensor(library)
        self.execution_time_s += time.perf_counter() - start
        self.body_count += 1
        return tuple(result[index] for index in range(current.dimension))

    def _vertex_expression(
        self,
        interaction: Any,
        *,
        left_slots: tuple[Any, ...],
        right_slots: tuple[Any, ...],
        output_slots: tuple[Any, ...],
    ) -> Any:
        from symbolica import Expression

        kind = int(interaction.vertex_kind)
        if kind == 0:
            return self.model.three_gluon_current_expression(
                left_slot=left_slots[0],
                right_slot=right_slots[0],
                output_slot=output_slots[0],
                left_momentum_tensor_name=self._momentum_tensor_name(interaction.left_id),
                right_momentum_tensor_name=self._momentum_tensor_name(interaction.right_id),
                dummy_prefix=f"compiled_dag_j{interaction.result_id}_i{interaction.id}",
            )
        if kind == 1:
            return self._tensor_name(str(symbols.two_gluon_to_tensor))(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 2:
            return self._tensor_name(str(symbols.tensor_gluon_to_gluon))(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 3:
            return self._tensor_name(str(symbols.gluon_tensor_to_gluon))(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 6:
            return self._tensor_name(
                _quark_vector_weyl_tensor_name(int(interaction.result.chirality))
            )(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 10:
            chirality = int(interaction.result.chirality)
            coupling = _weyl_coupling_for_chirality(chirality, interaction.coupling)
            return Expression.num(coupling) * self._tensor_name(
                _quark_vector_weyl_tensor_name(chirality)
            )(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        raise NativeEvaluationError(f"unsupported shared DAG vertex kind: {kind}")

    def _propagator_tensor_data(
        self,
        current: Any,
        momentum: tuple[Any, ...],
    ) -> list[Any]:
        pdg = int(current.pdg)
        if pdg == 21:
            return self.model.gluon_propagator_tensor_data(momentum)
        if 1 <= abs(pdg) <= 6 and int(current.chirality) != 0:
            return self.model.quark_weyl_propagator_tensor_data(
                momentum,
                chirality=int(current.chirality),
            )
        raise NativeEvaluationError(
            f"no Spenso propagator tensor for current {current.key}"
        )

    def _slots_for_current(self, current: Any, prefix: str) -> tuple[Any, ...]:
        pdg = int(current.pdg)
        if pdg == -21:
            return (self._aux6(f"{prefix}_A"),)
        if pdg in (21, 22, 23):
            return (self._mink(f"{prefix}_mu"),)
        if 1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16:
            return (self._weyl(f"{prefix}_alpha"),)
        raise NativeEvaluationError(f"unsupported current type for Spenso lowering: {pdg}")

    def _representation_for_current(self, current: Any) -> Any:
        pdg = int(current.pdg)
        if pdg == -21:
            return self._aux6
        if pdg in (21, 22, 23):
            return self._mink
        if 1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16:
            return self._weyl
        raise NativeEvaluationError(f"unsupported current type for Spenso lowering: {pdg}")

    @staticmethod
    def _current_tensor_name(current_id: int) -> str:
        return f"pyamplicol::compiled_dag_parent_current_{current_id}"

    @staticmethod
    def _momentum_tensor_name(current_id: int) -> str:
        return f"pyamplicol::compiled_dag_momentum_{current_id}"

    @staticmethod
    def _propagator_tensor_name(current_id: int) -> str:
        return f"pyamplicol::compiled_dag_propagator_{current_id}"


def _build_compiled_dag_parameter_layout(
    table: SharedCurrentTable,
    *,
    inline_external_wavefunctions: bool = False,
) -> _CompiledDAGParameterLayout:
    builder = ParamBuilder()
    source_offsets = [-1 for _ in table.currents]
    source_component_offsets: list[tuple[int, ...]] = [
        () for _ in table.currents
    ]
    if inline_external_wavefunctions:
        external_label_count = max(
            label
            for current in table.currents
            for label in current.external_labels
        )
        external_momentum_offsets: list[tuple[int, ...]] = []
        for label in range(1, external_label_count + 1):
            head = ("compiled_dag", "external_momentum", str(label))
            builder.add_parameter_list(
                head,
                4,
                role="compiled_dag_external_momentum",
                real_valued=True,
            )
            start = builder.positions[head][0]
            external_momentum_offsets.append(tuple(range(start, start + 4)))
        return _CompiledDAGParameterLayout(
            parameter_symbols=tuple(builder.parameter_symbols()),
            parameter_count=len(builder.values),
            source_offsets=tuple(source_offsets),
            source_component_offsets=tuple(source_component_offsets),
            momentum_offsets={},
            external_momentum_offsets=tuple(external_momentum_offsets),
            real_valued_inputs=tuple(builder.real_valued_inputs),
            inline_external_wavefunctions=True,
        )

    source_by_current_id = {source.current_id: source for source in table.sources}
    for current in table.currents:
        if not current.is_source:
            continue
        source = source_by_current_id[current.id]
        offsets: list[int] = []
        for component_index in range(current.dimension):
            if _compiled_dag_source_component_is_structural_zero(
                current,
                source,
                component_index,
            ):
                offsets.append(-1)
                continue
            head = (
                "compiled_dag",
                "source",
                str(current.id),
                str(component_index),
            )
            builder.add_parameter_list(
                head,
                1,
                role="compiled_dag_source_current",
                real_valued=_compiled_dag_source_component_is_structural_real(
                    current,
                    source,
                    component_index,
                ),
            )
            offsets.append(builder.positions[head][0])
        source_component_offsets[current.id] = tuple(offsets)
        source_offsets[current.id] = next(
            (offset for offset in offsets if offset >= 0),
            -1,
        )
    momentum_offsets = _register_shared_momentum_parameters(builder, table)
    return _CompiledDAGParameterLayout(
        parameter_symbols=tuple(builder.parameter_symbols()),
        parameter_count=len(builder.values),
        source_offsets=tuple(source_offsets),
        source_component_offsets=tuple(source_component_offsets),
        momentum_offsets=momentum_offsets,
        external_momentum_offsets=(),
        real_valued_inputs=tuple(builder.real_valued_inputs),
        inline_external_wavefunctions=False,
    )


def _compiled_dag_source_component_is_structural_zero(
    current: Any,
    source: SharedSourceCurrent,
    component_index: int,
) -> bool:
    pdg = int(current.pdg)
    if pdg == 21 and component_index == 0:
        return True
    if pdg == 23 and abs(int(source.helicity)) == 1 and component_index == 0:
        return True
    return False


def _compiled_dag_source_component_is_structural_real(
    current: Any,
    source: SharedSourceCurrent,
    component_index: int,
) -> bool:
    pdg = int(current.pdg)
    if pdg == 21 and component_index == 3:
        return True
    if pdg == 23 and component_index in (0, 3):
        return True
    if pdg == 23 and int(source.helicity) == 0:
        return True
    return False


def _compiled_dag_external_momentum_expressions(
    layout: _CompiledDAGParameterLayout,
) -> tuple[tuple[Any, ...], ...]:
    parameter_symbols = layout.parameter_symbols
    return tuple(
        tuple(parameter_symbols[offset] for offset in offsets)
        for offsets in layout.external_momentum_offsets
    )


def _compiled_dag_inline_momentum_expressions(
    table: SharedCurrentTable,
    external_momenta: tuple[tuple[Any, ...], ...],
) -> dict[int, tuple[Any, ...]]:
    from symbolica import Expression

    zero = Expression.num(0)
    momentum_expressions: dict[int, tuple[Any, ...]] = {}
    for current_id in _shared_momentum_current_ids(table):
        labels = _momentum_labels_for_current(table.currents[current_id])
        components: list[Any] = []
        for component_index in range(4):
            total = zero
            for label in labels:
                sign = -1.0 if label <= 2 else 1.0
                total = total + Expression.num(sign) * external_momenta[
                    label - 1
                ][component_index]
            components.append(total)
        momentum_expressions[current_id] = tuple(components)
    return momentum_expressions


def _compiled_dag_inline_source_expressions(
    table: SharedCurrentTable,
    external_momenta: tuple[tuple[Any, ...], ...],
    model: AmplicolSMLeadingColorModel,
) -> dict[int, tuple[Any, ...]]:
    source_expressions: dict[int, tuple[Any, ...]] = {}
    for source in table.sources:
        current = table.currents[source.current_id]
        pdg = int(current.pdg)
        if source.leg_label == 1:
            source_expressions[source.current_id] = (
                _inline_antiquark_source_expression(
                    external_momenta[0],
                    physical_helicity=source.physical_helicity,
                )
            )
        elif source.leg_label == 2:
            source_expressions[source.current_id] = (
                _inline_quark_source_expression(
                    external_momenta[1],
                    chirality=source.chirality,
                )
            )
        elif pdg == 21:
            source_expressions[source.current_id] = (
                _inline_gluon_source_expression(
                    external_momenta[source.leg_label - 1],
                    helicity=source.helicity,
                )
            )
        elif pdg == 23:
            source_expressions[source.current_id] = (
                _inline_z_source_expression(
                    external_momenta[source.leg_label - 1],
                    helicity=source.helicity,
                    mass=model.mass(23),
                )
            )
        else:
            raise NativeEvaluationError(
                f"cannot inline source wavefunction for current {current.key}"
            )
    return source_expressions


def _compiled_dag_current_use_counts(table: SharedCurrentTable) -> tuple[int, ...]:
    counts = [0 for _ in table.currents]
    for interaction in table.interactions:
        counts[interaction.left_id] += 1
        counts[interaction.right_id] += 1
    for amplitude in table.amplitudes:
        counts[amplitude.left_id] += 1
        counts[amplitude.right_id] += 1
    return tuple(counts)


def _is_symbolica_zero(value: Any) -> bool:
    if isinstance(value, (int, float, complex)):
        return value == 0
    return str(value) == "0"


def _alias_inline_source_expressions(
    table: SharedCurrentTable,
    source_expressions: dict[int, tuple[Any, ...]],
    alias_definitions: list[tuple[Any, Any]],
) -> int:
    alias_count = 0
    for source in table.sources:
        current = table.currents[source.current_id]
        aliased_components: list[Any] = []
        for component_index, expression in enumerate(
            source_expressions[source.current_id]
        ):
            if _compiled_dag_inline_source_component_is_structural_zero(
                current,
                source,
                component_index,
            ):
                aliased_components.append(expression)
                continue
            alias = _source_component_alias(
                current.id,
                component_index,
                real=_compiled_dag_inline_source_component_is_structural_real(
                    current,
                    source,
                    component_index,
                ),
            )
            alias_definitions.append((alias, expression))
            aliased_components.append(alias)
            alias_count += 1
        source_expressions[source.current_id] = tuple(aliased_components)
    return alias_count


def _alias_inline_momentum_expressions(
    momentum_expressions: dict[int, tuple[Any, ...]],
    alias_definitions: list[tuple[Any, Any]],
) -> int:
    alias_count = 0
    for current_id, expressions in list(momentum_expressions.items()):
        aliases = tuple(
            _momentum_component_alias(current_id, component_index)
            for component_index in range(4)
        )
        alias_definitions.extend(zip(aliases, expressions, strict=True))
        momentum_expressions[current_id] = aliases
        alias_count += 4
    return alias_count


def _compiled_dag_inline_source_component_is_structural_zero(
    current: Any,
    source: SharedSourceCurrent,
    component_index: int,
) -> bool:
    if source.leg_label == 1:
        if int(source.physical_helicity) == 1 and component_index == 1:
            return True
        if int(source.physical_helicity) == -1 and component_index == 0:
            return True
    if source.leg_label == 2:
        if int(source.chirality) == 1 and component_index == 0:
            return True
        if int(source.chirality) == -1 and component_index == 1:
            return True
    return _compiled_dag_source_component_is_structural_zero(
        current,
        source,
        component_index,
    )


def _compiled_dag_inline_source_component_is_structural_real(
    current: Any,
    source: SharedSourceCurrent,
    component_index: int,
) -> bool:
    if source.leg_label in (1, 2):
        return True
    return _compiled_dag_source_component_is_structural_real(
        current,
        source,
        component_index,
    )


def _compiled_dag_inline_source_pair_expressions(
    table: SharedCurrentTable,
    external_momenta: tuple[tuple[Any, ...], ...],
    model: AmplicolSMLeadingColorModel,
) -> dict[int, tuple[_RealPairExpression, ...]]:
    source_expressions: dict[int, tuple[_RealPairExpression, ...]] = {}
    for source in table.sources:
        current = table.currents[source.current_id]
        pdg = int(current.pdg)
        if source.leg_label == 1:
            source_expressions[source.current_id] = _inline_antiquark_source_pair(
                external_momenta[0],
                physical_helicity=source.physical_helicity,
            )
        elif source.leg_label == 2:
            source_expressions[source.current_id] = _inline_quark_source_pair(
                external_momenta[1],
                chirality=source.chirality,
            )
        elif pdg == 21:
            source_expressions[source.current_id] = _inline_gluon_source_pair(
                external_momenta[source.leg_label - 1],
                helicity=source.helicity,
            )
        elif pdg == 23:
            source_expressions[source.current_id] = _inline_z_source_pair(
                external_momenta[source.leg_label - 1],
                helicity=source.helicity,
                mass=model.mass(23),
            )
        else:
            raise NativeEvaluationError(
                f"cannot inline source wavefunction for current {current.key}"
            )
    return source_expressions


def _alias_inline_source_pair_expressions(
    table: SharedCurrentTable,
    source_expressions: dict[int, tuple[_RealPairExpression, ...]],
    alias_definitions: list[tuple[Any, Any]],
) -> int:
    from symbolica import Expression

    zero = Expression.num(0)
    alias_count = 0
    for source in table.sources:
        current = table.currents[source.current_id]
        aliased_components: list[_RealPairExpression] = []
        for component_index, expression in enumerate(
            source_expressions[source.current_id]
        ):
            if _compiled_dag_inline_source_component_is_structural_zero(
                current,
                source,
                component_index,
            ):
                aliased_components.append(expression)
                continue
            real_alias = _source_component_real_alias(current.id, component_index)
            imag_alias = _source_component_imag_alias(current.id, component_index)
            alias_definitions.append((real_alias, expression.real))
            alias_count += 1
            if _compiled_dag_inline_source_pair_imag_is_structural_zero(
                current,
                source,
                component_index,
            ):
                aliased_components.append(
                    _RealPairExpression(real_alias, zero)
                )
            else:
                alias_definitions.append((imag_alias, expression.imag))
                alias_count += 1
                aliased_components.append(
                    _RealPairExpression(real_alias, imag_alias)
                )
        source_expressions[source.current_id] = tuple(aliased_components)
    return alias_count


def _compiled_dag_inline_source_pair_imag_is_structural_zero(
    current: Any,
    source: SharedSourceCurrent,
    component_index: int,
) -> bool:
    pdg = int(current.pdg)
    if source.leg_label in (1, 2):
        return True
    if pdg == 21 and component_index in (0, 3):
        return True
    if pdg == 23 and int(source.helicity) == 0:
        return True
    if pdg == 23 and component_index in (0, 3):
        return True
    return False


def _rp_zero(zero: Any) -> _RealPairExpression:
    return _RealPairExpression(zero, zero)


def _rp_real(value: Any, zero: Any) -> _RealPairExpression:
    return _RealPairExpression(value, zero)


def _rp_add(left: _RealPairExpression, right: _RealPairExpression) -> _RealPairExpression:
    return _RealPairExpression(left.real + right.real, left.imag + right.imag)


def _rp_sub(left: _RealPairExpression, right: _RealPairExpression) -> _RealPairExpression:
    return _RealPairExpression(left.real - right.real, left.imag - right.imag)


def _rp_mul(left: _RealPairExpression, right: _RealPairExpression) -> _RealPairExpression:
    return _RealPairExpression(
        left.real * right.real - left.imag * right.imag,
        left.real * right.imag + left.imag * right.real,
    )


def _rp_mul_real(value: _RealPairExpression, factor: Any) -> _RealPairExpression:
    return _RealPairExpression(value.real * factor, value.imag * factor)


def _rp_i_times(value: _RealPairExpression) -> _RealPairExpression:
    return _RealPairExpression(-value.imag, value.real)


def _rp_minus_i_times(value: _RealPairExpression) -> _RealPairExpression:
    return _RealPairExpression(value.imag, -value.real)


def _rp_i_times_real_factor(
    value: _RealPairExpression,
    factor: Any,
) -> _RealPairExpression:
    return _rp_mul_real(_rp_i_times(value), factor)


def _rp_minus_i_times_real_factor(
    value: _RealPairExpression,
    factor: Any,
) -> _RealPairExpression:
    return _rp_mul_real(_rp_minus_i_times(value), factor)


def _rp_minkowski_dot(
    left: tuple[_RealPairExpression, ...],
    right: tuple[_RealPairExpression, ...],
) -> _RealPairExpression:
    total = _rp_mul(left[0], right[0])
    for index in range(1, 4):
        total = _rp_sub(total, _rp_mul(left[index], right[index]))
    return total


def _rp_minkowski_dot_momentum(
    vector: tuple[_RealPairExpression, ...],
    momentum: tuple[Any, ...],
) -> _RealPairExpression:
    total = _rp_mul_real(vector[0], momentum[0])
    for index in range(1, 4):
        total = _rp_sub(total, _rp_mul_real(vector[index], momentum[index]))
    return total


def _expr_vertex_pair(
    interaction: Any,
    left: tuple[_RealPairExpression, ...],
    right: tuple[_RealPairExpression, ...],
    momentum_expressions: dict[int, tuple[Any, ...]],
) -> tuple[_RealPairExpression, ...]:
    kind = int(interaction.vertex_kind)
    if kind == 0:
        return _expr_three_gluon_pair(
            left,
            momentum_expressions[interaction.left_id],
            right,
            momentum_expressions[interaction.right_id],
        )
    if kind == 1:
        return _expr_two_gluon_to_tensor_pair(left, right)
    if kind == 2:
        return _expr_tensor_gluon_to_gluon_pair(left, right)
    if kind == 3:
        return _expr_gluon_tensor_to_gluon_pair(left, right)
    if kind == 6:
        return _expr_quark_vector_weyl_pair(
            left,
            right,
            int(interaction.result.chirality),
        )
    if kind == 10:
        chirality = int(interaction.result.chirality)
        coupling = _weyl_coupling_for_chirality(chirality, interaction.coupling)
        return tuple(
            _rp_mul_real(component, coupling)
            for component in _expr_quark_vector_weyl_pair(left, right, chirality)
        )
    raise NativeEvaluationError(f"unsupported shared DAG vertex kind: {kind}")


def _expr_three_gluon_pair(
    left: tuple[_RealPairExpression, ...],
    left_momentum: tuple[Any, ...],
    right: tuple[_RealPairExpression, ...],
    right_momentum: tuple[Any, ...],
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    tmp1 = _rp_minkowski_dot(left, right)
    tmp2 = _rp_minkowski_dot_momentum(left, right_momentum)
    tmp3 = _rp_minkowski_dot_momentum(right, left_momentum)
    inv_sqrt2 = Expression.num(1.0 / math.sqrt(2.0))
    two = Expression.num(2.0)
    output: list[_RealPairExpression] = []
    for index in range(4):
        bracket = _rp_add(
            _rp_mul_real(tmp1, left_momentum[index] - right_momentum[index]),
            _rp_mul_real(
                _rp_sub(
                    _rp_mul(tmp2, right[index]),
                    _rp_mul(tmp3, left[index]),
                ),
                two,
            ),
        )
        output.append(_rp_i_times_real_factor(bracket, inv_sqrt2))
    return tuple(output)


def _expr_two_gluon_to_tensor_pair(
    left: tuple[_RealPairExpression, ...],
    right: tuple[_RealPairExpression, ...],
) -> tuple[_RealPairExpression, ...]:
    return (
        _rp_sub(_rp_mul(left[0], right[1]), _rp_mul(left[1], right[0])),
        _rp_sub(_rp_mul(left[0], right[2]), _rp_mul(left[2], right[0])),
        _rp_sub(_rp_mul(left[0], right[3]), _rp_mul(left[3], right[0])),
        _rp_sub(_rp_mul(left[1], right[2]), _rp_mul(left[2], right[1])),
        _rp_sub(_rp_mul(left[1], right[3]), _rp_mul(left[3], right[1])),
        _rp_sub(_rp_mul(left[2], right[3]), _rp_mul(left[3], right[2])),
    )


def _expr_tensor_gluon_to_gluon_pair(
    tensor: tuple[_RealPairExpression, ...],
    gluon: tuple[_RealPairExpression, ...],
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    half = Expression.num(0.5)
    return tuple(
        _rp_i_times_real_factor(value, half)
        for value in (
            _rp_add(
                _rp_add(_rp_mul(tensor[0], gluon[1]), _rp_mul(tensor[1], gluon[2])),
                _rp_mul(tensor[2], gluon[3]),
            ),
            _rp_add(
                _rp_add(_rp_mul(tensor[0], gluon[0]), _rp_mul(tensor[3], gluon[2])),
                _rp_mul(tensor[4], gluon[3]),
            ),
            _rp_add(
                _rp_sub(_rp_mul(tensor[1], gluon[0]), _rp_mul(tensor[3], gluon[1])),
                _rp_mul(tensor[5], gluon[3]),
            ),
            _rp_sub(
                _rp_sub(_rp_mul(tensor[2], gluon[0]), _rp_mul(tensor[4], gluon[1])),
                _rp_mul(tensor[5], gluon[2]),
            ),
        )
    )


def _expr_gluon_tensor_to_gluon_pair(
    gluon: tuple[_RealPairExpression, ...],
    tensor: tuple[_RealPairExpression, ...],
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    half = Expression.num(0.5)
    return tuple(
        _rp_i_times_real_factor(value, half)
        for value in (
            _rp_sub(
                _rp_sub(
                    _rp_sub(_rp_zero(Expression.num(0)), _rp_mul(gluon[1], tensor[0])),
                    _rp_mul(gluon[2], tensor[1]),
                ),
                _rp_mul(gluon[3], tensor[2]),
            ),
            _rp_sub(
                _rp_sub(
                    _rp_sub(_rp_zero(Expression.num(0)), _rp_mul(gluon[0], tensor[0])),
                    _rp_mul(gluon[2], tensor[3]),
                ),
                _rp_mul(gluon[3], tensor[4]),
            ),
            _rp_sub(
                _rp_add(
                    _rp_sub(_rp_zero(Expression.num(0)), _rp_mul(gluon[0], tensor[1])),
                    _rp_mul(gluon[1], tensor[3]),
                ),
                _rp_mul(gluon[3], tensor[5]),
            ),
            _rp_add(
                _rp_add(
                    _rp_sub(_rp_zero(Expression.num(0)), _rp_mul(gluon[0], tensor[2])),
                    _rp_mul(gluon[1], tensor[4]),
                ),
                _rp_mul(gluon[2], tensor[5]),
            ),
        )
    )


def _expr_vector_slash_terms_pair(
    vector: tuple[_RealPairExpression, ...],
) -> tuple[_RealPairExpression, _RealPairExpression, _RealPairExpression, _RealPairExpression]:
    return (
        _rp_add(vector[0], vector[3]),
        _rp_sub(vector[0], vector[3]),
        _rp_add(vector[1], _rp_i_times(vector[2])),
        _rp_sub(vector[1], _rp_i_times(vector[2])),
    )


def _expr_quark_vector_weyl_pair(
    quark: tuple[_RealPairExpression, ...],
    vector: tuple[_RealPairExpression, ...],
    chirality: int,
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    tmp1, tmp2, tmp3, tmp4 = _expr_vector_slash_terms_pair(vector)
    inv_sqrt2 = Expression.num(1.0 / math.sqrt(2.0))
    q1, q2 = quark
    if chirality == 1:
        return (
            _rp_i_times_real_factor(
                _rp_sub(_rp_mul(tmp2, q1), _rp_mul(tmp3, q2)),
                inv_sqrt2,
            ),
            _rp_i_times_real_factor(
                _rp_sub(_rp_mul(tmp1, q2), _rp_mul(tmp4, q1)),
                inv_sqrt2,
            ),
        )
    if chirality == -1:
        return (
            _rp_i_times_real_factor(
                _rp_add(_rp_mul(tmp1, q1), _rp_mul(tmp3, q2)),
                inv_sqrt2,
            ),
            _rp_i_times_real_factor(
                _rp_add(_rp_mul(tmp2, q2), _rp_mul(tmp4, q1)),
                inv_sqrt2,
            ),
        )
    raise NativeEvaluationError("Weyl quark-vector expression needs nonzero chirality")


def _expr_propagate_pair_with_denominator(
    current: Any,
    value: tuple[_RealPairExpression, ...],
    momentum: tuple[Any, ...],
    denominator: Any,
) -> tuple[_RealPairExpression, ...]:
    pdg = int(current.pdg)
    if pdg == 21:
        return tuple(_rp_minus_i_times_real_factor(component, 1.0 / denominator) for component in value)
    if 1 <= abs(pdg) <= 6 and int(current.chirality) != 0:
        return _expr_quark_propagator_weyl_pair_with_denominator(
            value,
            momentum,
            int(current.chirality),
            denominator,
        )
    return value


def _expr_quark_propagator_weyl_pair_with_denominator(
    quark: tuple[_RealPairExpression, ...],
    momentum: tuple[Any, ...],
    chirality: int,
    denominator: Any,
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    energy, px, py, pz = momentum
    zero = Expression.num(0)
    tmp1 = _rp_real(energy + pz, zero)
    tmp2 = _rp_real(energy - pz, zero)
    tmp3 = _RealPairExpression(px, py)
    tmp4 = _RealPairExpression(px, -py)
    q1, q2 = quark
    inv_den = 1.0 / denominator
    if chirality == 1:
        return (
            _rp_i_times_real_factor(
                _rp_add(_rp_mul(tmp1, q1), _rp_mul(tmp3, q2)),
                inv_den,
            ),
            _rp_i_times_real_factor(
                _rp_add(_rp_mul(tmp2, q2), _rp_mul(tmp4, q1)),
                inv_den,
            ),
        )
    if chirality == -1:
        return (
            _rp_i_times_real_factor(
                _rp_sub(_rp_mul(tmp2, q1), _rp_mul(tmp3, q2)),
                inv_den,
            ),
            _rp_i_times_real_factor(
                _rp_sub(_rp_mul(tmp1, q2), _rp_mul(tmp4, q1)),
                inv_den,
            ),
        )
    raise NativeEvaluationError("Weyl quark propagator expression needs nonzero chirality")


def _expr_dot_weyl_pair(
    left: tuple[_RealPairExpression, ...],
    right: tuple[_RealPairExpression, ...],
) -> _RealPairExpression:
    return _rp_add(_rp_mul(left[0], right[0]), _rp_mul(left[1], right[1]))


def _expr_vertex_with_aliases(
    interaction: Any,
    left: tuple[Any, ...],
    right: tuple[Any, ...],
    momentum_expressions: dict[int, tuple[Any, ...]],
    alias_definitions: list[tuple[Any, Any]],
) -> tuple[tuple[Any, ...], int]:
    kind = int(interaction.vertex_kind)
    if kind == 0:
        return _expr_three_gluon_with_aliases(
            interaction,
            left,
            momentum_expressions[interaction.left_id],
            right,
            momentum_expressions[interaction.right_id],
            alias_definitions,
        )
    if kind == 6:
        return _expr_quark_vector_weyl_with_aliases(
            interaction,
            left,
            right,
            int(interaction.result.chirality),
            alias_definitions,
            coupling=None,
        )
    if kind == 10:
        chirality = int(interaction.result.chirality)
        return _expr_quark_vector_weyl_with_aliases(
            interaction,
            left,
            right,
            chirality,
            alias_definitions,
            coupling=_weyl_coupling_for_chirality(chirality, interaction.coupling),
        )
    return _expr_vertex(interaction, left, right, momentum_expressions), 0


def _expr_three_gluon_with_aliases(
    interaction: Any,
    left: tuple[Any, ...],
    left_momentum: tuple[Any, ...],
    right: tuple[Any, ...],
    right_momentum: tuple[Any, ...],
    alias_definitions: list[tuple[Any, Any]],
) -> tuple[tuple[Any, ...], int]:
    tmp1 = _interaction_temporary_alias(interaction.id, 0)
    tmp2 = _interaction_temporary_alias(interaction.id, 1)
    tmp3 = _interaction_temporary_alias(interaction.id, 2)
    alias_definitions.extend(
        (
            (tmp1, _expr_minkowski_dot(left, right)),
            (tmp2, _expr_minkowski_dot_momentum(left, right_momentum)),
            (tmp3, _expr_minkowski_dot_momentum(right, left_momentum)),
        )
    )
    prefactor = 1j / math.sqrt(2.0)
    return (
        tuple(
            prefactor
            * (
                tmp1 * (left_momentum[index] - right_momentum[index])
                + 2.0 * (tmp2 * right[index] - tmp3 * left[index])
            )
            for index in range(4)
        ),
        3,
    )


def _expr_quark_vector_weyl_with_aliases(
    interaction: Any,
    quark: tuple[Any, ...],
    vector: tuple[Any, ...],
    chirality: int,
    alias_definitions: list[tuple[Any, Any]],
    *,
    coupling: float | None,
) -> tuple[tuple[Any, ...], int]:
    tmp_aliases = tuple(
        _interaction_temporary_alias(interaction.id, index)
        for index in range(4)
    )
    alias_definitions.extend(
        zip(tmp_aliases, _expr_vector_slash_terms(vector), strict=True)
    )
    prefactor = 1j / math.sqrt(2.0)
    if coupling is not None:
        prefactor *= coupling
    tmp1, tmp2, tmp3, tmp4 = tmp_aliases
    q1, q2 = quark
    if chirality == 1:
        return (
            (
                prefactor * (tmp2 * q1 - tmp3 * q2),
                prefactor * (tmp1 * q2 - tmp4 * q1),
            ),
            4,
        )
    if chirality == -1:
        return (
            (
                prefactor * (tmp1 * q1 + tmp3 * q2),
                prefactor * (tmp2 * q2 + tmp4 * q1),
            ),
            4,
        )
    raise NativeEvaluationError("Weyl quark-vector expression needs nonzero chirality")


def _expr_propagate_with_denominator(
    current: Any,
    value: tuple[Any, ...],
    momentum: tuple[Any, ...],
    denominator: Any,
) -> tuple[Any, ...]:
    pdg = int(current.pdg)
    if pdg == 21:
        prefactor = -1j / denominator
        return tuple(component * prefactor for component in value)
    if 1 <= abs(pdg) <= 6 and int(current.chirality) != 0:
        return _expr_quark_propagator_weyl_with_denominator(
            value,
            momentum,
            int(current.chirality),
            denominator,
        )
    return value


def _expr_quark_propagator_weyl_with_denominator(
    quark: tuple[Any, ...],
    momentum: tuple[Any, ...],
    chirality: int,
    denominator: Any,
) -> tuple[Any, ...]:
    energy, px, py, pz = momentum
    prefactor = 1j / denominator
    tmp1 = energy + pz
    tmp2 = energy - pz
    tmp3 = px + 1j * py
    tmp4 = px - 1j * py
    q1, q2 = quark
    if chirality == 1:
        return (
            (tmp1 * q1 + tmp3 * q2) * prefactor,
            (tmp2 * q2 + tmp4 * q1) * prefactor,
        )
    if chirality == -1:
        return (
            (tmp2 * q1 - tmp3 * q2) * prefactor,
            (tmp1 * q2 - tmp4 * q1) * prefactor,
        )
    raise NativeEvaluationError("Weyl quark propagator expression needs nonzero chirality")


def _inline_antiquark_source_expression(
    momentum: tuple[Any, ...],
    *,
    physical_helicity: int,
) -> tuple[Any, ...]:
    from symbolica import Expression

    zero = Expression.num(0)
    norm = (Expression.num(2.0) * momentum[0]).sqrt()
    if physical_helicity == 1:
        return (norm, zero)
    if physical_helicity == -1:
        return (zero, norm)
    raise NativeEvaluationError(
        f"unexpected antiquark physical helicity: {physical_helicity}"
    )


def _inline_quark_source_expression(
    momentum: tuple[Any, ...],
    *,
    chirality: int,
) -> tuple[Any, ...]:
    from symbolica import Expression

    zero = Expression.num(0)
    norm = (Expression.num(2.0) * momentum[0]).sqrt()
    if chirality == 1:
        return (zero, norm)
    if chirality == -1:
        return (-norm, zero)
    raise NativeEvaluationError(f"unexpected quark chirality: {chirality}")


def _inline_antiquark_source_pair(
    momentum: tuple[Any, ...],
    *,
    physical_helicity: int,
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    zero = Expression.num(0)
    norm = (Expression.num(2.0) * momentum[0]).sqrt()
    if physical_helicity == 1:
        return (_rp_real(norm, zero), _rp_zero(zero))
    if physical_helicity == -1:
        return (_rp_zero(zero), _rp_real(norm, zero))
    raise NativeEvaluationError(
        f"unexpected antiquark physical helicity: {physical_helicity}"
    )


def _inline_quark_source_pair(
    momentum: tuple[Any, ...],
    *,
    chirality: int,
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    zero = Expression.num(0)
    norm = (Expression.num(2.0) * momentum[0]).sqrt()
    if chirality == 1:
        return (_rp_zero(zero), _rp_real(norm, zero))
    if chirality == -1:
        return (_rp_real(-norm, zero), _rp_zero(zero))
    raise NativeEvaluationError(f"unexpected quark chirality: {chirality}")


def _inline_gluon_source_expression(
    momentum: tuple[Any, ...],
    *,
    helicity: int,
) -> tuple[Any, ...]:
    from symbolica import Expression

    energy, px, py, pz = momentum
    zero = Expression.num(0)
    imag = Expression.num(1j)
    sqh = Expression.num(math.sqrt(0.5))
    hel = Expression.num(float(helicity))
    pt = (px * px + py * py).sqrt()
    pzpt = pz / (energy * pt) * sqh * hel
    return (
        zero,
        -px * pzpt - imag * py / pt * sqh,
        -py * pzpt + imag * px / pt * sqh,
        hel * pt / energy * sqh,
    )


def _inline_gluon_source_pair(
    momentum: tuple[Any, ...],
    *,
    helicity: int,
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    energy, px, py, pz = momentum
    zero = Expression.num(0)
    sqh = Expression.num(math.sqrt(0.5))
    hel = Expression.num(float(helicity))
    pt = (px * px + py * py).sqrt()
    pzpt = pz / (energy * pt) * sqh * hel
    return (
        _rp_zero(zero),
        _RealPairExpression(-px * pzpt, -py / pt * sqh),
        _RealPairExpression(-py * pzpt, px / pt * sqh),
        _rp_real(hel * pt / energy * sqh, zero),
    )


def _inline_z_source_expression(
    momentum: tuple[Any, ...],
    *,
    helicity: int,
    mass: float,
) -> tuple[Any, ...]:
    from symbolica import Expression

    if mass == 0.0:
        raise NativeEvaluationError("massless vector wavefunction not expected for Z")
    energy, px, py, pz = momentum
    zero = Expression.num(0)
    imag = Expression.num(1j)
    sqh = Expression.num(math.sqrt(0.5))
    mass_expr = Expression.num(mass)
    pt2 = px * px + py * py
    pt = pt2.sqrt()
    pp = (pt2 + pz * pz).sqrt()
    emp = energy / (mass_expr * pp)
    if helicity == 0:
        return (
            pp / mass_expr,
            px * emp,
            py * emp,
            pz * emp,
        )
    if helicity not in (-1, 1):
        raise NativeEvaluationError(f"unexpected Z helicity: {helicity}")
    hel = Expression.num(float(helicity))
    pzpt = pz / (pp * pt) * sqh * hel
    return (
        zero,
        -px * pzpt - imag * py / pt * sqh,
        -py * pzpt + imag * px / pt * sqh,
        hel * pt / pp * sqh,
    )


def _inline_z_source_pair(
    momentum: tuple[Any, ...],
    *,
    helicity: int,
    mass: float,
) -> tuple[_RealPairExpression, ...]:
    from symbolica import Expression

    if mass == 0.0:
        raise NativeEvaluationError("massless vector wavefunction not expected for Z")
    energy, px, py, pz = momentum
    zero = Expression.num(0)
    sqh = Expression.num(math.sqrt(0.5))
    mass_expr = Expression.num(mass)
    pt2 = px * px + py * py
    pt = pt2.sqrt()
    pp = (pt2 + pz * pz).sqrt()
    emp = energy / (mass_expr * pp)
    if helicity == 0:
        return (
            _rp_real(pp / mass_expr, zero),
            _rp_real(px * emp, zero),
            _rp_real(py * emp, zero),
            _rp_real(pz * emp, zero),
        )
    if helicity not in (-1, 1):
        raise NativeEvaluationError(f"unexpected Z helicity: {helicity}")
    hel = Expression.num(float(helicity))
    pzpt = pz / (pp * pt) * sqh * hel
    return (
        _rp_zero(zero),
        _RealPairExpression(-px * pzpt, -py / pt * sqh),
        _RealPairExpression(-py * pzpt, px / pt * sqh),
        _rp_real(hel * pt / pp * sqh, zero),
    )


def _compiled_dag_source_parameter_specs(
    table: SharedCurrentTable,
    layout: _CompiledDAGParameterLayout,
) -> tuple[_CompiledDAGSourceParameterSpec, ...]:
    specs: list[_CompiledDAGSourceParameterSpec] = []
    for source in table.sources:
        current = table.currents[source.current_id]
        specs.append(
            _CompiledDAGSourceParameterSpec(
                current_id=source.current_id,
                offset=layout.source_offsets[source.current_id],
                component_offsets=layout.source_component_offsets[source.current_id],
                dimension=current.dimension,
                leg_label=source.leg_label,
                helicity=source.helicity,
                physical_helicity=source.physical_helicity,
                chirality=source.chirality,
            )
        )
    return tuple(specs)


def _compiled_dag_source_fill_plan(
    specs: Sequence[_CompiledDAGSourceParameterSpec],
) -> _CompiledDAGSourceFillPlan:
    anti_offsets: dict[int, tuple[int, ...]] = {}
    quark_offsets: dict[int, tuple[int, ...]] = {}
    gluon_offsets: dict[tuple[int, int], tuple[int, ...]] = {}
    z_offsets: dict[int, tuple[int, ...]] = {}
    z_label = max((spec.leg_label for spec in specs), default=0)
    for spec in specs:
        if spec.leg_label == 1:
            anti_offsets[spec.physical_helicity] = spec.component_offsets
        elif spec.leg_label == 2:
            quark_offsets[spec.chirality] = spec.component_offsets
        elif spec.leg_label == z_label:
            z_offsets[spec.helicity] = spec.component_offsets
        elif 3 <= spec.leg_label < z_label:
            gluon_offsets[(spec.leg_label - 3, spec.helicity)] = (
                spec.component_offsets
            )
    return _CompiledDAGSourceFillPlan(
        anti_offsets=anti_offsets,
        quark_offsets=quark_offsets,
        gluon_offsets=gluon_offsets,
        z_offsets=z_offsets,
        z_label=z_label,
    )


def _external_momentum_array(
    points: Sequence[tuple[ExternalMomentum, ...]],
    external_label_count: int,
) -> np.ndarray:
    if not points:
        return np.empty((0, external_label_count, 4), dtype=np.float64)
    return np.fromiter(
        (
            component
            for point in points
            for particle in point
            for component in particle.momentum
        ),
        dtype=np.float64,
        count=len(points) * external_label_count * 4,
    ).reshape(len(points), external_label_count, 4)


def _compiled_dag_point_batch(
    points: Sequence[Sequence[ExternalMomentum]],
    *,
    gluon_count: int,
    validate: bool,
) -> tuple[tuple[ExternalMomentum, ...], ...]:
    if validate:
        return tuple(
            _validate_z_gluon_point(particles, gluon_count=gluon_count)
            for particles in points
        )
    return tuple(tuple(particles) for particles in points)


def _fill_compiled_dag_source_parameters_batch_fast(
    rows: np.ndarray,
    plan: _CompiledDAGSourceFillPlan,
    external_momenta: np.ndarray,
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
) -> bool:
    """Fill source currents for physical q q~ -> Z + gluons batches.

    The generic batch filler mirrors the scalar wavefunction routines and handles
    both energy signs. In the compiled-DAG hot path we evaluate physical
    incoming beams and positive-energy final states, so the sign branches are
    known and several helicities share the same kinematic auxiliaries.
    """

    if rows.shape[0] == 0:
        return True
    if plan.z_label != gluon_count + 3:
        return False
    if (
        set(plan.anti_offsets) != {-1, 1}
        or set(plan.quark_offsets) != {-1, 1}
        or set(plan.z_offsets) != {-1, 0, 1}
    ):
        return False
    for gluon_index in range(gluon_count):
        if (gluon_index, -1) not in plan.gluon_offsets:
            return False
        if (gluon_index, 1) not in plan.gluon_offsets:
            return False
    if np.any(external_momenta[:, :, 0] <= 0.0):
        return False

    anti_closure_momentum = -external_momenta[:, 0, :]
    quark_start_momentum = -external_momenta[:, 1, :]
    if np.any(anti_closure_momentum[:, 0] >= 0.0):
        return False
    if np.any(quark_start_momentum[:, 0] >= 0.0):
        return False

    _fill_antiquark_sources_negative(
        rows,
        plan.anti_offsets,
        anti_closure_momentum,
    )
    _fill_quark_sources_negative(
        rows,
        plan.quark_offsets,
        quark_start_momentum,
    )
    _fill_gluon_sources_positive(
        rows,
        plan.gluon_offsets,
        external_momenta[:, 2 : 2 + gluon_count, :],
    )
    _fill_z_sources_positive(
        rows,
        plan.z_offsets,
        external_momenta[:, plan.z_label - 1, :],
        mass=model.mass(23),
    )
    return True


def _assign_component_column(
    rows: np.ndarray,
    offsets: tuple[int, ...],
    component_index: int,
    values: np.ndarray | complex | float,
) -> None:
    offset = offsets[component_index]
    if offset >= 0:
        rows[:, offset] = values


def _fill_quark_sources_negative(
    rows: np.ndarray,
    offsets_by_chirality: dict[int, tuple[int, ...]],
    momenta: np.ndarray,
) -> None:
    energy = momenta[:, 0]
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]

    sqp0p3 = -np.sqrt(np.maximum(-(energy + pz), 0.0))
    zero = sqp0p3 == 0.0
    safe = ~zero

    def chi2_for(helicity: int) -> np.ndarray:
        chi2 = np.empty_like(sqp0p3, dtype=np.complex128)
        if np.any(zero):
            chi2[zero] = complex(-helicity) * np.sqrt(2.0 * np.abs(energy[zero]))
        if np.any(safe):
            chi2[safe] = (
                helicity * px[safe] + 1j * py[safe]
            ) / sqp0p3[safe]
        return chi2

    # chirality +1 uses helicity -1 and maps to (chi1, chi2).
    offsets = offsets_by_chirality[1]
    _assign_component_column(rows, offsets, 0, sqp0p3)
    _assign_component_column(rows, offsets, 1, chi2_for(-1))

    # chirality -1 uses helicity +1 and maps to (chi2, chi1).
    offsets = offsets_by_chirality[-1]
    _assign_component_column(rows, offsets, 0, chi2_for(1))
    _assign_component_column(rows, offsets, 1, sqp0p3)


def _fill_antiquark_sources_negative(
    rows: np.ndarray,
    offsets_by_physical_helicity: dict[int, tuple[int, ...]],
    momenta: np.ndarray,
) -> None:
    energy = momenta[:, 0]
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]

    sqp0p3 = np.sqrt(np.maximum(-(energy + pz), 0.0))
    zero = sqp0p3 == 0.0
    safe = ~zero

    def chi2_for(helicity: int) -> np.ndarray:
        chi2 = np.empty_like(sqp0p3, dtype=np.complex128)
        if np.any(zero):
            chi2[zero] = complex(-helicity) * np.sqrt(2.0 * np.abs(energy[zero]))
        if np.any(safe):
            chi2[safe] = (
                -helicity * px[safe] - 1j * py[safe]
            ) / sqp0p3[safe]
        return chi2

    # physical helicity +1 calls (helicity=+1, chirality=-1) and maps to
    # (chi1, chi2) in the negative-energy branch.
    offsets = offsets_by_physical_helicity[1]
    _assign_component_column(rows, offsets, 0, sqp0p3)
    _assign_component_column(rows, offsets, 1, chi2_for(1))

    # physical helicity -1 calls (helicity=-1, chirality=+1) and maps to
    # (chi2, chi1) in the negative-energy branch.
    offsets = offsets_by_physical_helicity[-1]
    _assign_component_column(rows, offsets, 0, chi2_for(-1))
    _assign_component_column(rows, offsets, 1, sqp0p3)


def _fill_gluon_sources_positive(
    rows: np.ndarray,
    offsets_by_gluon_helicity: dict[tuple[int, int], tuple[int, ...]],
    momenta: np.ndarray,
) -> None:
    if momenta.shape[1] == 0:
        return
    sqh = math.sqrt(0.5)
    for gluon_index in range(momenta.shape[1]):
        momentum = momenta[:, gluon_index, :]
        energy = momentum[:, 0]
        px = momentum[:, 1]
        py = momentum[:, 2]
        pz = momentum[:, 3]
        pt = np.sqrt(px * px + py * py)
        nonzero_pt = pt != 0.0
        for helicity in (-1, 1):
            offsets = offsets_by_gluon_helicity[(gluon_index, helicity)]
            hel = float(helicity)
            _assign_component_column(rows, offsets, 0, 0.0)
            _assign_component_column(rows, offsets, 3, hel * pt / energy * sqh)
            if np.any(nonzero_pt):
                pzpt = np.zeros_like(pt)
                pzpt[nonzero_pt] = (
                    pz[nonzero_pt]
                    / (energy[nonzero_pt] * pt[nonzero_pt])
                    * sqh
                    * hel
                )
                values1 = np.empty_like(pzpt, dtype=np.complex128)
                values2 = np.empty_like(pzpt, dtype=np.complex128)
                values1[nonzero_pt] = (
                    -px[nonzero_pt] * pzpt[nonzero_pt]
                    - 1j * py[nonzero_pt] / pt[nonzero_pt] * sqh
                )
                values2[nonzero_pt] = (
                    -py[nonzero_pt] * pzpt[nonzero_pt]
                    + 1j * px[nonzero_pt] / pt[nonzero_pt] * sqh
                )
            else:
                values1 = np.empty_like(pt, dtype=np.complex128)
                values2 = np.empty_like(pt, dtype=np.complex128)
            if np.any(~nonzero_pt):
                zero = ~nonzero_pt
                values1[zero] = -hel * sqh
                values2[zero] = 1j * np.copysign(sqh, pz[zero])
            _assign_component_column(rows, offsets, 1, values1)
            _assign_component_column(rows, offsets, 2, values2)


def _fill_z_sources_positive(
    rows: np.ndarray,
    offsets_by_helicity: dict[int, tuple[int, ...]],
    momenta: np.ndarray,
    *,
    mass: float,
) -> None:
    if mass == 0.0:
        raise NativeEvaluationError("massless vector wavefunction not expected for Z")
    energy = momenta[:, 0]
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]
    sqh = math.sqrt(0.5)
    pt2 = px * px + py * py
    pp = np.minimum(energy, np.sqrt(pt2 + pz * pz))
    pt = np.minimum(pp, np.sqrt(pt2))
    zero_pp = pp == 0.0
    nonzero_pp = ~zero_pp
    nonzero_pt = pt != 0.0
    safe_pt = nonzero_pp & nonzero_pt
    zero_pt = nonzero_pp & ~nonzero_pt

    emp = np.zeros_like(energy)
    emp[nonzero_pp] = energy[nonzero_pp] / (mass * pp[nonzero_pp])
    for helicity in (-1, 0, 1):
        offsets = offsets_by_helicity[helicity]
        hel = float(helicity)
        nsvahl = abs(helicity)
        hel0 = 1.0 - abs(hel)

        v0 = hel0 * pp / mass
        v1 = np.empty_like(energy, dtype=np.complex128)
        v2 = np.empty_like(energy, dtype=np.complex128)
        v3 = hel0 * pz * emp
        v3[nonzero_pp] += hel * pt[nonzero_pp] / pp[nonzero_pp] * sqh

        if np.any(zero_pp):
            v1[zero_pp] = -hel * sqh
            v2[zero_pp] = 1j * nsvahl * sqh
            v3[zero_pp] = hel0
        if np.any(safe_pt):
            pzpt = (
                pz[safe_pt] / (pp[safe_pt] * pt[safe_pt]) * sqh * hel
            )
            v1[safe_pt] = (
                hel0 * px[safe_pt] * emp[safe_pt]
                - px[safe_pt] * pzpt
                - 1j * nsvahl * py[safe_pt] / pt[safe_pt] * sqh
            )
            v2[safe_pt] = (
                hel0 * py[safe_pt] * emp[safe_pt]
                - py[safe_pt] * pzpt
                + 1j * nsvahl * px[safe_pt] / pt[safe_pt] * sqh
            )
        if np.any(zero_pt):
            v1[zero_pt] = -hel * sqh
            v2[zero_pt] = 1j * nsvahl * np.copysign(sqh, pz[zero_pt])

        _assign_component_column(rows, offsets, 0, v0)
        _assign_component_column(rows, offsets, 1, v1)
        _assign_component_column(rows, offsets, 2, v2)
        _assign_component_column(rows, offsets, 3, v3)


def _fill_compiled_dag_source_parameters_batch(
    rows: np.ndarray,
    specs: Sequence[_CompiledDAGSourceParameterSpec],
    external_momenta: np.ndarray,
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
) -> None:
    if rows.shape[0] == 0:
        return
    z_mass = model.mass(23)
    anti_closure_momentum = -external_momenta[:, 0, :]
    quark_start_momentum = -external_momenta[:, 1, :]

    anti_cache: dict[int, np.ndarray] = {}
    quark_cache: dict[int, np.ndarray] = {}
    gluon_cache: dict[tuple[int, int], np.ndarray] = {}
    z_cache: dict[int, np.ndarray] = {}

    for spec in specs:
        if spec.leg_label == 1:
            wavefunctions = anti_cache.get(spec.physical_helicity)
            if wavefunctions is None:
                if spec.physical_helicity == 1:
                    wavefunctions = _ext_antiquark_weyl_batch(
                        anti_closure_momentum,
                        helicity=1,
                        chirality=-1,
                    )
                else:
                    wavefunctions = _ext_antiquark_weyl_batch(
                        anti_closure_momentum,
                        helicity=-1,
                        chirality=1,
                    )
                anti_cache[spec.physical_helicity] = wavefunctions
        elif spec.leg_label == 2:
            wavefunctions = quark_cache.get(spec.chirality)
            if wavefunctions is None:
                if spec.chirality == 1:
                    wavefunctions = _ext_quark_weyl_batch(
                        quark_start_momentum,
                        helicity=-1,
                        chirality=1,
                    )
                else:
                    wavefunctions = _ext_quark_weyl_batch(
                        quark_start_momentum,
                        helicity=1,
                        chirality=-1,
                    )
                quark_cache[spec.chirality] = wavefunctions
        elif 3 <= spec.leg_label < gluon_count + 3:
            gluon_index = spec.leg_label - 3
            cache_key = (gluon_index, spec.helicity)
            wavefunctions = gluon_cache.get(cache_key)
            if wavefunctions is None:
                wavefunctions = _ext_gluon_cmplx_batch(
                    external_momenta[:, spec.leg_label - 1, :],
                    spec.helicity,
                )
                gluon_cache[cache_key] = wavefunctions
        elif spec.leg_label == gluon_count + 3:
            wavefunctions = z_cache.get(spec.helicity)
            if wavefunctions is None:
                wavefunctions = _ext_massive_vector_batch(
                    external_momenta[:, spec.leg_label - 1, :],
                    spec.helicity,
                    z_mass,
                )
                z_cache[spec.helicity] = wavefunctions
        else:
            raise NativeEvaluationError(
                f"unexpected source leg label: {spec.leg_label}"
            )
        for component_index, offset in enumerate(spec.component_offsets):
            if offset >= 0:
                rows[:, offset] = wavefunctions[:, component_index]


def _ext_quark_weyl_batch(
    momenta: np.ndarray,
    *,
    helicity: int,
    chirality: int,
) -> np.ndarray:
    energy = momenta[:, 0]
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]
    output = np.zeros((momenta.shape[0], 2), dtype=np.complex128)

    positive = energy > 0.0
    if np.any(positive):
        idx = np.nonzero(positive)[0]
        sqp0p3 = np.sqrt(np.maximum(energy[idx] + pz[idx], 0.0))
        chi1 = sqp0p3.astype(np.complex128)
        chi2 = np.empty(idx.shape[0], dtype=np.complex128)
        zero = sqp0p3 == 0.0
        if np.any(zero):
            chi2[zero] = complex(-helicity) * np.sqrt(2.0 * energy[idx][zero])
        if np.any(~zero):
            chi2[~zero] = (
                helicity * px[idx][~zero] - 1j * py[idx][~zero]
            ) / sqp0p3[~zero]
        if helicity == 1 and chirality == 1:
            output[idx, 0] = chi1
            output[idx, 1] = chi2
        elif helicity == -1 and chirality == -1:
            output[idx, 0] = chi2
            output[idx, 1] = chi1

    negative = ~positive
    if np.any(negative):
        idx = np.nonzero(negative)[0]
        sqp0p3 = -np.sqrt(np.maximum(-(energy[idx] + pz[idx]), 0.0))
        chi1 = sqp0p3.astype(np.complex128)
        chi2 = np.empty(idx.shape[0], dtype=np.complex128)
        zero = sqp0p3 == 0.0
        if np.any(zero):
            chi2[zero] = complex(-helicity) * np.sqrt(2.0 * np.abs(energy[idx][zero]))
        if np.any(~zero):
            chi2[~zero] = (
                helicity * px[idx][~zero] + 1j * py[idx][~zero]
            ) / sqp0p3[~zero]
        if helicity == -1 and chirality == 1:
            output[idx, 0] = chi1
            output[idx, 1] = chi2
        elif helicity == 1 and chirality == -1:
            output[idx, 0] = chi2
            output[idx, 1] = chi1

    return output


def _ext_antiquark_weyl_batch(
    momenta: np.ndarray,
    *,
    helicity: int,
    chirality: int,
) -> np.ndarray:
    energy = momenta[:, 0]
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]
    output = np.zeros((momenta.shape[0], 2), dtype=np.complex128)

    positive = energy > 0.0
    if np.any(positive):
        idx = np.nonzero(positive)[0]
        sqp0p3 = -np.sqrt(np.maximum(energy[idx] + pz[idx], 0.0))
        chi1 = sqp0p3.astype(np.complex128)
        chi2 = np.empty(idx.shape[0], dtype=np.complex128)
        zero = sqp0p3 == 0.0
        if np.any(zero):
            chi2[zero] = complex(-helicity) * np.sqrt(2.0 * energy[idx][zero])
        if np.any(~zero):
            chi2[~zero] = (
                -helicity * px[idx][~zero] + 1j * py[idx][~zero]
            ) / sqp0p3[~zero]
        if helicity == 1 and chirality == 1:
            output[idx, 0] = chi2
            output[idx, 1] = chi1
        elif helicity == -1 and chirality == -1:
            output[idx, 0] = chi1
            output[idx, 1] = chi2

    negative = ~positive
    if np.any(negative):
        idx = np.nonzero(negative)[0]
        sqp0p3 = np.sqrt(np.maximum(-(energy[idx] + pz[idx]), 0.0))
        chi1 = sqp0p3.astype(np.complex128)
        chi2 = np.empty(idx.shape[0], dtype=np.complex128)
        zero = sqp0p3 == 0.0
        if np.any(zero):
            chi2[zero] = complex(-helicity) * np.sqrt(2.0 * np.abs(energy[idx][zero]))
        if np.any(~zero):
            chi2[~zero] = (
                -helicity * px[idx][~zero] - 1j * py[idx][~zero]
            ) / sqp0p3[~zero]
        if helicity == -1 and chirality == 1:
            output[idx, 0] = chi2
            output[idx, 1] = chi1
        elif helicity == 1 and chirality == -1:
            output[idx, 0] = chi1
            output[idx, 1] = chi2

    return output


def _ext_gluon_cmplx_batch(momenta: np.ndarray, helicity: int) -> np.ndarray:
    energy = momenta[:, 0]
    if np.any(energy == 0.0):
        raise NativeEvaluationError("cannot generate external gluon with zero energy")
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]
    sqh = math.sqrt(0.5)
    output = np.zeros((momenta.shape[0], 4), dtype=np.complex128)

    positive = energy > 0.0
    if np.any(positive):
        idx = np.nonzero(positive)[0]
        hel = float(helicity)
        pp = energy[idx]
        pt = np.sqrt(px[idx] ** 2 + py[idx] ** 2)
        output[idx, 3] = hel * pt / pp * sqh
        nonzero_pt = pt != 0.0
        if np.any(nonzero_pt):
            nz = idx[nonzero_pt]
            pt_nz = pt[nonzero_pt]
            pzpt = pz[nz] / (pp[nonzero_pt] * pt_nz) * sqh * hel
            output[nz, 1] = -px[nz] * pzpt - 1j * py[nz] / pt_nz * sqh
            output[nz, 2] = -py[nz] * pzpt + 1j * px[nz] / pt_nz * sqh
        if np.any(~nonzero_pt):
            z = idx[~nonzero_pt]
            output[z, 1] = -hel * sqh
            output[z, 2] = 1j * np.copysign(sqh, pz[z])

    negative = ~positive
    if np.any(negative):
        idx = np.nonzero(negative)[0]
        hel = float(-helicity)
        pp = -energy[idx]
        pt = np.sqrt(px[idx] ** 2 + py[idx] ** 2)
        output[idx, 3] = hel * pt / pp * sqh
        nonzero_pt = pt != 0.0
        if np.any(nonzero_pt):
            nz = idx[nonzero_pt]
            pt_nz = pt[nonzero_pt]
            pzpt = -pz[nz] / (pp[nonzero_pt] * pt_nz) * sqh * hel
            output[nz, 1] = px[nz] * pzpt + 1j * py[nz] / pt_nz * sqh
            output[nz, 2] = py[nz] * pzpt - 1j * px[nz] / pt_nz * sqh
        if np.any(~nonzero_pt):
            z = idx[~nonzero_pt]
            output[z, 1] = -hel * sqh
            output[z, 2] = -1j * np.copysign(sqh, pz[z])

    return output


def _ext_massive_vector_batch(
    momenta: np.ndarray,
    helicity: int,
    mass: float,
) -> np.ndarray:
    if mass == 0.0:
        raise NativeEvaluationError("massless vector wavefunction not expected for Z")
    energy = momenta[:, 0]
    px = momenta[:, 1]
    py = momenta[:, 2]
    pz = momenta[:, 3]
    sqh = math.sqrt(0.5)
    hel = float(helicity)
    nsvahl = abs(helicity)
    hel0 = 1.0 - abs(hel)
    pt2 = px**2 + py**2
    pp = np.minimum(energy, np.sqrt(pt2 + pz**2))
    pt = np.minimum(pp, np.sqrt(pt2))
    output = np.zeros((momenta.shape[0], 4), dtype=np.complex128)

    zero_pp = pp == 0.0
    if np.any(zero_pp):
        idx = np.nonzero(zero_pp)[0]
        output[idx, 1] = -hel * sqh
        output[idx, 2] = 1j * nsvahl * sqh
        output[idx, 3] = hel0

    nonzero_pp = ~zero_pp
    if np.any(nonzero_pp):
        idx = np.nonzero(nonzero_pp)[0]
        emp = energy[idx] / (mass * pp[idx])
        output[idx, 0] = hel0 * pp[idx] / mass
        output[idx, 3] = hel0 * pz[idx] * emp + hel * pt[idx] / pp[idx] * sqh
        nonzero_pt = pt[idx] != 0.0
        if np.any(nonzero_pt):
            nz = idx[nonzero_pt]
            pt_nz = pt[nz]
            pp_nz = pp[nz]
            emp_nz = energy[nz] / (mass * pp_nz)
            pzpt = pz[nz] / (pp_nz * pt_nz) * sqh * hel
            output[nz, 1] = (
                hel0 * px[nz] * emp_nz
                - px[nz] * pzpt
                - 1j * nsvahl * py[nz] / pt_nz * sqh
            )
            output[nz, 2] = (
                hel0 * py[nz] * emp_nz
                - py[nz] * pzpt
                + 1j * nsvahl * px[nz] / pt_nz * sqh
            )
        if np.any(~nonzero_pt):
            z = idx[~nonzero_pt]
            output[z, 1] = -hel * sqh
            output[z, 2] = 1j * nsvahl * np.copysign(sqh, pz[z])

    return output


def _fill_compiled_dag_source_parameters(
    row: np.ndarray,
    specs: Sequence[_CompiledDAGSourceParameterSpec],
    particles: tuple[ExternalMomentum, ...],
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
) -> None:
    incoming_quark = particles[0]
    incoming_antiquark = particles[1]
    anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
    quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
    gluons = particles[2 : 2 + gluon_count]
    z_boson = particles[-1]

    anti_cache: dict[int, tuple[complex, ...]] = {}
    quark_cache: dict[int, tuple[complex, ...]] = {}
    gluon_cache: dict[tuple[int, int], tuple[complex, ...]] = {}
    z_cache: dict[int, tuple[complex, ...]] = {}
    z_mass = model.mass(23)

    for spec in specs:
        if spec.leg_label == 1:
            wavefunction = anti_cache.get(spec.physical_helicity)
            if wavefunction is None:
                if spec.physical_helicity == 1:
                    wavefunction = _ext_antiquark_weyl(
                        anti_closure_momentum,
                        1,
                        -1,
                    )
                else:
                    wavefunction = _ext_antiquark_weyl(
                        anti_closure_momentum,
                        -1,
                        1,
                    )
                anti_cache[spec.physical_helicity] = wavefunction
        elif spec.leg_label == 2:
            wavefunction = quark_cache.get(spec.chirality)
            if wavefunction is None:
                if spec.chirality == 1:
                    wavefunction = _ext_quark_weyl(
                        quark_start_momentum,
                        -1,
                        1,
                    )
                else:
                    wavefunction = _ext_quark_weyl(
                        quark_start_momentum,
                        1,
                        -1,
                    )
                quark_cache[spec.chirality] = wavefunction
        elif 3 <= spec.leg_label < gluon_count + 3:
            gluon_index = spec.leg_label - 3
            cache_key = (gluon_index, spec.helicity)
            wavefunction = gluon_cache.get(cache_key)
            if wavefunction is None:
                wavefunction = _ext_gluon_cmplx(
                    gluons[gluon_index].momentum,
                    spec.helicity,
                )
                gluon_cache[cache_key] = wavefunction
        elif spec.leg_label == gluon_count + 3:
            wavefunction = z_cache.get(spec.helicity)
            if wavefunction is None:
                wavefunction = _ext_massive_vector(
                    z_boson.momentum,
                    spec.helicity,
                    z_mass,
                )
                z_cache[spec.helicity] = wavefunction
        else:
            raise NativeEvaluationError(
                f"unexpected source leg label: {spec.leg_label}"
            )
        for component_index, offset in enumerate(spec.component_offsets):
            if offset >= 0:
                row[offset] = wavefunction[component_index]


def _z_gluon_matrix_element_normalization(
    model: AmplicolSMLeadingColorModel,
    pdgs: Sequence[int],
    *,
    gluon_count: int,
) -> float:
    pdg_tuple = tuple(int(pdg) for pdg in pdgs)
    color_factor = model.leading_color_factor(pdg_tuple)
    average_factor = _initial_state_average_factor(pdg_tuple[:2])
    identical_factor = _final_state_identical_factor(pdg_tuple[2:])
    coupling_factor = (
        (4.0 * math.pi * model.alpha_s_me_check) ** gluon_count
        * (2.0 * 4.0 * math.pi * model.alpha_ew)
    )
    return color_factor * coupling_factor / (average_factor * identical_factor)


def build_compiled_dag_helicity_filter(
    process: str,
    table: SharedCurrentTable,
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
    sample_count: int = 10,
    seed: int = 12345,
    relative_tolerance: float = 1.0e-12,
    zero_tolerance: float = 1.0e-300,
    phase_space_mode: CompiledDAGHelicityFilterMode = "rambo",
    progress_enabled: bool = False,
) -> HelicityFilter:
    if sample_count < 1:
        raise NativeEvaluationError("helicity filter sample count must be positive")
    points = _compiled_dag_helicity_filter_sample_points(
        process,
        model,
        gluon_count=gluon_count,
        sample_count=sample_count,
        seed=seed,
        phase_space_mode=phase_space_mode,
    )
    helicities = tuple(amplitude.helicities for amplitude in table.amplitudes)
    squared_by_helicity: dict[tuple[int, ...], list[float]] = {
        entry: [] for entry in helicities
    }
    native = LeadingColorZJetsNativeEvaluator(model)
    for point in _compiled_dag_progress(
        points,
        enabled=progress_enabled,
        label="helicity warmup",
        metadata=f"samples={len(points):<4} hel={len(helicities):<5}",
    ):
        evaluation = native.evaluate_z_gluons(
            process,
            point,
            gluon_count=gluon_count,
        )
        amplitudes = {
            contribution.helicities: contribution.amplitude
            for contribution in evaluation.helicity_contributions
        }
        for entry in helicities:
            amplitude = amplitudes.get(entry, 0j)
            squared_by_helicity[entry].append(
                float((amplitude * amplitude.conjugate()).real)
            )

    active = [
        entry
        for entry in helicities
        if any(value > zero_tolerance for value in squared_by_helicity[entry])
    ]

    assigned: set[tuple[int, ...]] = set()
    filter_entries: list[HelicityFilterEntry] = []
    for entry in active:
        if entry in assigned:
            continue
        equivalents = [entry]
        assigned.add(entry)
        for candidate in active:
            if candidate in assigned:
                continue
            if _same_squared_helicity_signature(
                squared_by_helicity[entry],
                squared_by_helicity[candidate],
                relative_tolerance=relative_tolerance,
                zero_tolerance=zero_tolerance,
            ):
                equivalents.append(candidate)
                assigned.add(candidate)
        filter_entries.append(
            HelicityFilterEntry(
                helicities=entry,
                multiplicity=len(equivalents),
                equivalents=tuple(equivalents),
            )
        )

    return HelicityFilter(
        entries=tuple(filter_entries),
        original_count=len(helicities),
        sample_count=len(points),
        relative_tolerance=relative_tolerance,
        zero_tolerance=zero_tolerance,
        sample_sqrt_s=tuple(point[0].momentum[0] + point[1].momentum[0] for point in points),
        phase_space_mode=phase_space_mode,
        seed=seed,
    )


def _compiled_dag_helicity_filter_sample_points(
    process: str,
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
    sample_count: int,
    seed: int,
    phase_space_mode: CompiledDAGHelicityFilterMode,
) -> tuple[tuple[ExternalMomentum, ...], ...]:
    base_sqrt_s = max(model.sqrt_s, 1000.0)
    if phase_space_mode == "canonical":
        native = LeadingColorZJetsNativeEvaluator(model)
        return tuple(
            native.canonical_z_gluon_point(
                process,
                gluon_count=gluon_count,
                sqrt_s=base_sqrt_s * (1.0 + 0.037 * index),
            )
            for index in range(sample_count)
        )
    if phase_space_mode == "rambo":
        return tuple(
            rambo_z_gluon_point(
                process,
                model,
                gluon_count=gluon_count,
                sqrt_s=base_sqrt_s * (1.0 + 0.037 * index),
                seed=seed + index,
            )
            for index in range(sample_count)
        )
    raise NativeEvaluationError(
        f"unsupported compiled DAG helicity filter phase-space mode: {phase_space_mode!r}"
    )


def _apply_helicity_filter_to_shared_table(
    table: SharedCurrentTable,
    helicity_filter: HelicityFilter,
) -> SharedCurrentTable:
    amplitudes_by_helicity = {
        amplitude.helicities: amplitude
        for amplitude in table.amplitudes
    }
    filtered_amplitudes: list[Any] = []
    for entry in helicity_filter.entries:
        amplitude = amplitudes_by_helicity.get(entry.helicities)
        if amplitude is None:
            raise NativeEvaluationError(
                f"helicity filter entry is not present in shared table: {entry.helicities}"
            )
        filtered_amplitudes.append(
            replace(
                amplitude,
                multiplicity=entry.multiplicity,
            )
        )
    return replace(table, amplitudes=tuple(filtered_amplitudes))


def _same_squared_helicity_signature(
    left: Sequence[float],
    right: Sequence[float],
    *,
    relative_tolerance: float,
    zero_tolerance: float,
) -> bool:
    if len(left) != len(right):
        return False
    for left_value, right_value in zip(left, right, strict=True):
        denominator = max(abs(left_value) + abs(right_value), zero_tolerance)
        if abs(left_value - right_value) / denominator > relative_tolerance:
            return False
    return True


def _max_relative_complex_rows_difference(
    left_rows: Sequence[Sequence[complex]],
    right_rows: Sequence[Sequence[complex]],
) -> float:
    max_relative = 0.0
    for left_row, right_row in zip(left_rows, right_rows, strict=True):
        for left, right in zip(left_row, right_row, strict=True):
            denominator = max(abs(left), abs(right), 1.0e-300)
            max_relative = max(max_relative, abs(left - right) / denominator)
    return max_relative


def _max_relative_real_rows_difference(
    left_values: Sequence[float],
    right_values: Sequence[float],
) -> float:
    max_relative = 0.0
    for left, right in zip(left_values, right_values, strict=True):
        denominator = max(abs(left), abs(right), 1.0e-300)
        max_relative = max(max_relative, abs(left - right) / denominator)
    return max_relative


def _compiled_dag_momentum_offsets_and_labels(
    table: SharedCurrentTable,
    layout: _CompiledDAGParameterLayout,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    if layout.inline_external_wavefunctions:
        return ()
    labels_by_offset: dict[int, tuple[int, ...]] = {}
    for current_id in _shared_momentum_current_ids(table):
        offset = layout.momentum_offsets[current_id]
        labels_by_offset.setdefault(
            offset,
            _momentum_labels_for_current(table.currents[current_id]),
        )
    return tuple(sorted(labels_by_offset.items()))


def _compiled_dag_current_metadata(
    table: SharedCurrentTable,
) -> list[dict[str, object]]:
    return [
        {
            "id": current.id,
            "pdg": current.pdg,
            "chirality": current.chirality,
            "external_labels": list(current.external_labels),
            "source_ids": list(current.source_ids),
            "ext_source_bits": current.ext_source_bits,
            "is_source": current.is_source,
            "needs_propagator": current.needs_propagator,
            "dimension": current.dimension,
        }
        for current in table.currents
    ]


def _compiled_dag_root_metadata(
    table: SharedCurrentTable,
) -> list[dict[str, object]]:
    return [
        {
            "output_index": index,
            "left_current_id": amplitude.left_id,
            "right_current_id": amplitude.right_id,
            "helicities": list(amplitude.helicities),
            "multiplicity": amplitude.multiplicity,
        }
        for index, amplitude in enumerate(table.amplitudes)
    ]


def _symbolica_provenance() -> dict[str, object]:
    try:
        import symbolica
        from symbolica import Expression
    except Exception as exc:
        return {
            "available": False,
            "error": str(exc),
        }
    local_versions = getattr(symbolica, "LOCAL_VERSIONS", None)
    return {
        "available": True,
        "version": getattr(symbolica, "__version__", None),
        "local_versions": (
            dict(local_versions) if isinstance(local_versions, dict) else None
        ),
        "symbolica_commit": (
            local_versions.get("symbolica")
            if isinstance(local_versions, dict)
            else None
        ),
        "spenso_commit": (
            local_versions.get("spenso") if isinstance(local_versions, dict) else None
        ),
        "idenso_commit": (
            local_versions.get("idenso") if isinstance(local_versions, dict) else None
        ),
        "evaluator_aliases_hook": _symbolica_evaluator_alias_available(),
        "expression_alias_method": hasattr(Expression, "alias"),
    }


def _compiled_dag_artifact_fingerprint(artifact: dict[str, object]) -> str:
    payload = dict(artifact)
    payload.pop("artifact_fingerprint", None)
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _resolve_compiled_dag_compiled_preset(
    preset: str,
    *,
    gluon_count: int | None,
    inline_asm: str,
    optimization_level: int,
    output_chunk_size: int | None,
) -> tuple[str, int, int | None]:
    effective_inline_asm, effective_optimization_level, effective_chunk_size = (
        _resolve_compiled_preset(
            preset,
            gluon_count=gluon_count,
            inline_asm=inline_asm,
            optimization_level=optimization_level,
            output_chunk_size=output_chunk_size,
        )
    )
    if output_chunk_size is None:
        effective_chunk_size = None
    return effective_inline_asm, effective_optimization_level, effective_chunk_size


def _compiled_dag_progress(
    items: Sequence[Any],
    *,
    enabled: bool,
    label: str,
    metadata: str,
) -> Iterable[Any]:
    if not enabled:
        return items
    try:
        from colorama import Fore, Style  # type: ignore[import-untyped]
        from progressbar import Bar, ETA, Percentage, ProgressBar  # type: ignore[import-not-found]
    except ImportError:
        return items
    widgets = [
        Fore.CYAN,
        f" {label:<23}",
        Style.RESET_ALL,
        " ",
        Fore.YELLOW,
        f"{metadata:<22}",
        Style.RESET_ALL,
        " ",
        Percentage(),
        " ",
        Bar(),
        " ",
        ETA(),
    ]
    return ProgressBar(max_value=len(items), widgets=widgets)(items)


def _current_component_alias(current_id: int, component_index: int) -> Any:
    return symbols.symbol(f"compiled_dag_alias__j{current_id}__c{component_index}")


def _current_component_real_alias(current_id: int, component_index: int) -> Any:
    return symbols.real_symbol(
        f"compiled_dag_alias_re__j{current_id}__c{component_index}"
    )


def _current_component_imag_alias(current_id: int, component_index: int) -> Any:
    return symbols.real_symbol(
        f"compiled_dag_alias_im__j{current_id}__c{component_index}"
    )


def _source_component_alias(
    current_id: int,
    component_index: int,
    *,
    real: bool = False,
) -> Any:
    if real:
        return symbols.real_symbol(
            f"compiled_dag_source_alias_real__j{current_id}__c{component_index}"
        )
    return symbols.symbol(f"compiled_dag_source_alias__j{current_id}__c{component_index}")


def _source_component_real_alias(current_id: int, component_index: int) -> Any:
    return symbols.real_symbol(
        f"compiled_dag_source_alias_re__j{current_id}__c{component_index}"
    )


def _source_component_imag_alias(current_id: int, component_index: int) -> Any:
    return symbols.real_symbol(
        f"compiled_dag_source_alias_im__j{current_id}__c{component_index}"
    )


def _momentum_component_alias(current_id: int, component_index: int) -> Any:
    return symbols.real_symbol(
        f"compiled_dag_momentum_alias_real__j{current_id}__c{component_index}"
    )


def _propagator_denominator_alias(current_id: int) -> Any:
    return symbols.real_symbol(f"compiled_dag_prop_den_real__j{current_id}")


def _interaction_temporary_alias(interaction_id: int, temporary_index: int) -> Any:
    return symbols.symbol(
        f"compiled_dag_vertex_tmp__i{interaction_id}__t{temporary_index}"
    )


def _symbolica_evaluator_alias_available() -> bool:
    try:
        from symbolica import Expression, S

        x = S("pyamplicol_alias_probe_x")
        Expression.evaluator_multiple(
            [x],
            [x],
            aliases=[],
            jit_compile=False,
        )
        return True
    except TypeError:
        return False
    except Exception:
        return False


def _compiled_dag_evaluator_artifact_manifest(
    evaluator: Any,
    artifact_dir: Path,
    *,
    label: str,
) -> dict[str, Any]:
    try:
        return _symbolica_evaluator_artifact_manifest(evaluator, artifact_dir)
    except NativeEvaluationError:
        chunks = getattr(evaluator, "_evaluators", None)
        if isinstance(chunks, tuple):
            return {
                "kind": "chunked-symbolica-evaluator",
                "chunks": [
                    _compiled_dag_evaluator_artifact_manifest(
                        chunk,
                        artifact_dir,
                        label=f"{label}_chunk_{index}",
                    )
                    for index, chunk in enumerate(chunks)
                ],
            }
        save = getattr(evaluator, "save", None)
        if not callable(save):
            raise
        state_dir = artifact_dir / "evaluators"
        state_dir.mkdir(parents=True, exist_ok=True)
        state_path = state_dir / f"{label}.evaluator.bin"
        state_path.write_bytes(save())
        return {
            "kind": "symbolica-evaluator-state",
            "state_path": _artifact_path_for_manifest(state_path, artifact_dir),
        }


def _load_compiled_dag_evaluator_artifact(
    manifest: Any,
    artifact_dir: Path,
) -> Any:
    if not isinstance(manifest, dict):
        raise NativeEvaluationError("compiled DAG evaluator artifact is invalid")
    if manifest.get("kind") == "symbolica-evaluator-state":
        from symbolica import Evaluator

        state_path = _artifact_path_from_manifest(
            str(manifest["state_path"]),
            artifact_dir,
        )
        state = state_path.read_bytes()
        try:
            return Evaluator.load(state, {})
        except TypeError:
            return Evaluator.load(state)
    if manifest.get("kind") == "chunked-symbolica-evaluator":
        chunks = manifest.get("chunks")
        if not isinstance(chunks, list):
            raise NativeEvaluationError("chunked evaluator artifact is missing chunks")
        return _ChunkedSymbolicaEvaluator(
            tuple(
                _load_compiled_dag_evaluator_artifact(chunk, artifact_dir)
                for chunk in chunks
            )
        )
    return _load_symbolica_evaluator_artifact(manifest, artifact_dir)


def _evaluate_real_outputs(evaluator: Any, parameter_rows: Any) -> np.ndarray:
    return np.asarray(evaluator.evaluate(parameter_rows), dtype=np.float64)


def _weighted_real_pair_abs2_sums_array(
    evaluated: np.ndarray,
    weights: np.ndarray,
    output_length: int,
) -> np.ndarray:
    real = evaluated[:, : 2 * output_length : 2]
    imag = evaluated[:, 1 : 2 * output_length : 2]
    squared = real * real + imag * imag
    return np.dot(squared, weights)


__all__ = [
    "CompiledDAGLowering",
    "CompiledDAGMetadata",
    "ZGluonCompiledDAGEvaluator",
    "build_compiled_dag_helicity_filter",
]
