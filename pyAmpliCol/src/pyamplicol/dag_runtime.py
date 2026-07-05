from __future__ import annotations

import base64
from importlib import metadata as importlib_metadata
import json
import math
import os
import platform
import shutil
import tempfile
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field, replace
from itertools import product
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import numpy as np

from .lowering import (
    _current_dimension,
    _current_key_tuple,
    _current_needs_propagator,
    _source_currents,
    _weyl_coupling_for_chirality,
)
from .matrix import CurrentKey, NativeMatrixElementGenerator, RecursionGraph
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
    _final_state_identical_factor,
    _initial_state_average_factor,
)
from .params import ParamBuilder
from .symbols import symbols

ProgressCallback = Callable[[dict[str, object]], None]


def _report_progress(
    callback: ProgressCallback | None,
    *,
    stage: str,
    item: str,
    increment: int = 0,
    total: int | None = None,
) -> None:
    if callback is None:
        return
    payload: dict[str, object] = {
        "stage": stage,
        "item": item,
    }
    if increment:
        payload["increment"] = int(increment)
    if total is not None:
        payload["total"] = int(total)
    callback(payload)


def generation_progress_total(*, gluon_count: int, split_vertex_current_stages: bool = False) -> int:
    """Return the number of coarse generation steps exposed to the CLI."""

    source_total = 4 + 2 * gluon_count + 3
    gluon_total = sum(
        (gluon_count - length + 1) * (2**length)
        for length in range(2, gluon_count + 1)
    )
    quark_no_z_total = 2 * sum(2**end for end in range(1, gluon_count + 1))
    quark_with_z_total = 2 * 3 * sum(2**end for end in range(0, gluon_count + 1))
    amplitude_total = 2 * 3 * (2**gluon_count)
    stage_total = gluon_count + 1
    compile_total = 2 * stage_total if split_vertex_current_stages else stage_total
    compile_total += 1  # amplitude stage
    return (
        source_total
        + gluon_total
        + quark_no_z_total
        + quark_with_z_total
        + amplitude_total
        + compile_total
        + 2  # compiled materialization plus artifact/process save
    )


@dataclass(frozen=True)
class SymbolicaEvaluatorSettings:
    backend: str = "jit"
    iterations: int = 1
    cpe_iterations: int | None = None
    n_cores: int = 4
    direct_translation: bool = True
    jit_direct_translation: bool = False
    jit_optimization_level: int = 3
    max_horner_scheme_variables: int = 500
    max_common_pair_cache_entries: int = 1000000
    max_common_pair_distance: int = 100
    collect_factors: bool = False
    compiled_preset: str = "adaptive"
    compiled_inline_asm: str = "default"
    compiled_optimization_level: int = 3
    compiled_native: bool = True
    compiler_path: str | None = None
    compiler_flags: tuple[str, ...] = ()
    compiled_output_chunk_size: int | None = None
    compiled_chunk_compile_workers: int = 1
    compiled_output_dir: str | None = None
    raw_sum_final_stage: bool = False
    real_param_sqrt_real: bool = False
    real_param_log_real: bool = False
    real_param_powf_real: bool = False
    real_param_real_if_args_real: bool = False

    def __post_init__(self) -> None:
        if self.backend not in ("jit", "compiled-complex", "compiled-complex-4x"):
            raise NativeEvaluationError(
                "symbolica evaluator backend must be 'jit', "
                "'compiled-complex', or 'compiled-complex-4x'"
            )
        if self.compiled_preset not in (
            "manual",
            "adaptive",
            "generation",
            "balanced",
            "runtime",
            "runtime-o3",
        ):
            raise NativeEvaluationError(
                "symbolica compiled preset must be 'manual', 'adaptive', "
                "'generation', 'balanced', 'runtime', or 'runtime-o3'"
            )
        if self.iterations < 1:
            raise NativeEvaluationError("symbolica iterations must be positive")
        if self.cpe_iterations is not None and self.cpe_iterations < 0:
            raise NativeEvaluationError("symbolica cpe iterations must be non-negative")
        if self.n_cores < 1:
            raise NativeEvaluationError("symbolica n_cores must be positive")
        if self.jit_optimization_level not in (0, 1, 2, 3):
            raise NativeEvaluationError(
                "symbolica jit optimization level must be 0, 1, 2, or 3"
            )
        if self.compiled_optimization_level not in (0, 1, 2, 3):
            raise NativeEvaluationError(
                "symbolica compiled optimization level must be 0, 1, 2, or 3"
            )
        if (
            self.compiled_output_chunk_size is not None
            and self.compiled_output_chunk_size < 1
        ):
            raise NativeEvaluationError(
                "symbolica compiled output chunk size must be positive"
            )
        if self.compiled_chunk_compile_workers < 1:
            raise NativeEvaluationError(
                "symbolica compiled chunk compile workers must be positive"
            )
        if self.max_horner_scheme_variables < 1:
            raise NativeEvaluationError(
                "symbolica max_horner_scheme_variables must be positive"
            )
        if self.max_common_pair_cache_entries < 1:
            raise NativeEvaluationError(
                "symbolica max_common_pair_cache_entries must be positive"
            )
        if self.max_common_pair_distance < 1:
            raise NativeEvaluationError(
                "symbolica max_common_pair_distance must be positive"
            )

    def to_json_dict(self) -> dict[str, object]:
        return {
            "backend": self.backend,
            "iterations": self.iterations,
            "cpe_iterations": self.cpe_iterations,
            "n_cores": self.n_cores,
            "direct_translation": self.direct_translation,
            "jit_direct_translation": self.jit_direct_translation,
            "jit_optimization_level": self.jit_optimization_level,
            "max_horner_scheme_variables": self.max_horner_scheme_variables,
            "max_common_pair_cache_entries": self.max_common_pair_cache_entries,
            "max_common_pair_distance": self.max_common_pair_distance,
            "collect_factors": self.collect_factors,
            "compiled_preset": self.compiled_preset,
            "compiled_inline_asm": self.compiled_inline_asm,
            "compiled_optimization_level": self.compiled_optimization_level,
            "compiled_native": self.compiled_native,
            "compiler_path": self.compiler_path,
            "compiler_flags": list(self.compiler_flags),
            "effective_compiler_flags": list(_compiled_compiler_flags(self)),
            "compiled_output_chunk_size": self.compiled_output_chunk_size,
            "compiled_chunk_compile_workers": self.compiled_chunk_compile_workers,
            "compiled_output_dir": self.compiled_output_dir,
            "raw_sum_final_stage": self.raw_sum_final_stage,
            "real_param_sqrt_real": self.real_param_sqrt_real,
            "real_param_log_real": self.real_param_log_real,
            "real_param_powf_real": self.real_param_powf_real,
            "real_param_real_if_args_real": self.real_param_real_if_args_real,
        }


@dataclass(frozen=True)
class DAGEvaluationTiming:
    source_fill_time_s: float = 0.0
    momentum_setup_time_s: float = 0.0
    parameter_pack_time_s: float = 0.0
    evaluator_time_s: float = 0.0
    output_transfer_time_s: float = 0.0
    result_reduction_time_s: float = 0.0

    def __add__(self, other: "DAGEvaluationTiming") -> "DAGEvaluationTiming":
        return DAGEvaluationTiming(
            source_fill_time_s=self.source_fill_time_s + other.source_fill_time_s,
            momentum_setup_time_s=(
                self.momentum_setup_time_s + other.momentum_setup_time_s
            ),
            parameter_pack_time_s=(
                self.parameter_pack_time_s + other.parameter_pack_time_s
            ),
            evaluator_time_s=self.evaluator_time_s + other.evaluator_time_s,
            output_transfer_time_s=(
                self.output_transfer_time_s + other.output_transfer_time_s
            ),
            result_reduction_time_s=(
                self.result_reduction_time_s + other.result_reduction_time_s
            ),
        )

    @property
    def python_overhead_time_s(self) -> float:
        return (
            self.source_fill_time_s
            + self.momentum_setup_time_s
            + self.parameter_pack_time_s
            + self.output_transfer_time_s
            + self.result_reduction_time_s
        )

    @property
    def measured_time_s(self) -> float:
        return self.python_overhead_time_s + self.evaluator_time_s

    def with_result_reduction(self, elapsed_s: float) -> "DAGEvaluationTiming":
        return DAGEvaluationTiming(
            source_fill_time_s=self.source_fill_time_s,
            momentum_setup_time_s=self.momentum_setup_time_s,
            parameter_pack_time_s=self.parameter_pack_time_s,
            evaluator_time_s=self.evaluator_time_s,
            output_transfer_time_s=self.output_transfer_time_s,
            result_reduction_time_s=self.result_reduction_time_s + elapsed_s,
        )

    def to_json_dict(self) -> dict[str, float]:
        return {
            "source_fill_time_s": self.source_fill_time_s,
            "momentum_setup_time_s": self.momentum_setup_time_s,
            "parameter_pack_time_s": self.parameter_pack_time_s,
            "evaluator_time_s": self.evaluator_time_s,
            "output_transfer_time_s": self.output_transfer_time_s,
            "result_reduction_time_s": self.result_reduction_time_s,
            "python_overhead_time_s": self.python_overhead_time_s,
            "measured_time_s": self.measured_time_s,
        }


@dataclass(frozen=True)
class DAGEvaluatorMetadata:
    process: str
    kernel: str
    gluon_count: int
    graph_current_count: int
    graph_interaction_count: int
    graph_amplitude_count: int
    shared_current_count: int
    shared_source_current_count: int
    shared_interaction_count: int
    shared_amplitude_count: int
    current_table_build_time_s: float
    block_build_time_s: float
    symbolica_evaluator_build_time_s: float
    symbolica_parameter_count: int
    symbolica_output_count: int
    batch_size: int
    merge_evaluators_strategy: bool
    split_vertex_current_stages: bool
    verbose_evaluator_build: bool
    symbolica_evaluator_settings: dict[str, object]

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "kernel": self.kernel,
            "gluon_count": self.gluon_count,
            "graph_current_count": self.graph_current_count,
            "graph_interaction_count": self.graph_interaction_count,
            "graph_amplitude_count": self.graph_amplitude_count,
            "shared_current_count": self.shared_current_count,
            "shared_source_current_count": self.shared_source_current_count,
            "shared_interaction_count": self.shared_interaction_count,
            "shared_amplitude_count": self.shared_amplitude_count,
            "current_table_build_time_s": self.current_table_build_time_s,
            "block_build_time_s": self.block_build_time_s,
            "symbolica_evaluator_build_time_s": self.symbolica_evaluator_build_time_s,
            "symbolica_parameter_count": self.symbolica_parameter_count,
            "symbolica_output_count": self.symbolica_output_count,
            "batch_size": self.batch_size,
            "merge_evaluators_strategy": self.merge_evaluators_strategy,
            "split_vertex_current_stages": self.split_vertex_current_stages,
            "verbose_evaluator_build": self.verbose_evaluator_build,
            "symbolica_evaluator_settings": self.symbolica_evaluator_settings,
        }


@dataclass(frozen=True)
class DAGEvaluatorBuildOptions:
    merge_evaluators_strategy: bool = False
    split_vertex_current_stages: bool = False
    verbose_evaluator_build: bool = False
    symbolica_settings: SymbolicaEvaluatorSettings = field(
        default_factory=SymbolicaEvaluatorSettings
    )


@dataclass(frozen=True)
class SharedCurrentNode:
    id: int
    key: CurrentKey
    ext_source_bits: int
    source_ids: tuple[int, ...]
    is_source: bool
    needs_propagator: bool
    dimension: int

    @property
    def pdg(self) -> int:
        return self.key.pdg

    @property
    def external_labels(self) -> tuple[int, ...]:
        return self.key.external_labels

    @property
    def chirality(self) -> int:
        return self.key.chirality


@dataclass(frozen=True)
class SharedSourceCurrent:
    current_id: int
    leg_label: int
    helicity: int
    physical_helicity: int
    chirality: int
    source_bit: int


@dataclass(frozen=True)
class SharedInteractionNode:
    id: int
    vertex_kind: int
    left_id: int
    right_id: int
    result_id: int
    left: CurrentKey
    right: CurrentKey
    result: CurrentKey
    coupling: tuple[float, float]


@dataclass(frozen=True)
class SharedAmplitudeRecord:
    left_id: int
    right_id: int
    helicities: tuple[int, ...]
    multiplicity: int = 1


@dataclass(frozen=True)
class SharedCurrentTable:
    currents: tuple[SharedCurrentNode, ...]
    sources: tuple[SharedSourceCurrent, ...]
    interactions: tuple[SharedInteractionNode, ...]
    interactions_by_result: tuple[tuple[int, ...], ...]
    amplitudes: tuple[SharedAmplitudeRecord, ...]


@dataclass(frozen=True)
class _SharedGlobalParameterLayout:
    parameter_symbols: tuple[Any, ...]
    parameter_count: int
    current_offsets: tuple[int, ...]
    momentum_offsets: dict[int, int]
    real_valued_inputs: tuple[int, ...]
    interaction_offsets: tuple[int, ...] = ()


class ZGluonDAGEvaluator:
    """Shared-current DAG evaluator using spenso/Symbolica block kernels."""

    def __init__(
        self,
        process: str,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        graph: RecursionGraph | None = None,
        batch_size: int = 16,
        merge_evaluators_strategy: bool = False,
        split_vertex_current_stages: bool = False,
        verbose_evaluator_build: bool = False,
        symbolica_evaluator_backend: str = "jit",
        symbolica_iterations: int = 1,
        symbolica_cpe_iterations: int | None = None,
        symbolica_n_cores: int = 4,
        symbolica_direct_translation: bool = True,
        symbolica_jit_direct_translation: bool = False,
        symbolica_jit_optimization_level: int = 3,
        symbolica_max_horner_scheme_variables: int = 500,
        symbolica_max_common_pair_cache_entries: int = 1000000,
        symbolica_max_common_pair_distance: int = 100,
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
        progress_callback: ProgressCallback | None = None,
    ) -> None:
        if batch_size < 1:
            raise NativeEvaluationError("batch_size must be positive")
        self.process = process
        self.model = model or AmplicolSMLeadingColorModel()
        self.graph = graph or _z_gluon_graph(process, self.model)
        self.gluon_count = _validate_z_gluon_graph(self.graph)
        self.batch_size = int(batch_size)
        self.progress_callback = progress_callback
        self.last_runtime_timing = DAGEvaluationTiming()
        self.loaded_evaluator_artifact_metadata: dict[str, Any] | None = None
        (
            effective_compiled_inline_asm,
            effective_compiled_optimization_level,
            effective_compiled_output_chunk_size,
        ) = _resolve_compiled_preset(
            symbolica_compiled_preset,
            gluon_count=self.gluon_count,
            inline_asm=symbolica_compiled_inline_asm,
            optimization_level=symbolica_compiled_optimization_level,
            output_chunk_size=symbolica_compiled_output_chunk_size,
        )
        symbolica_settings = SymbolicaEvaluatorSettings(
            backend=symbolica_evaluator_backend,
            iterations=symbolica_iterations,
            cpe_iterations=symbolica_cpe_iterations,
            n_cores=symbolica_n_cores,
            direct_translation=symbolica_direct_translation,
            jit_direct_translation=symbolica_jit_direct_translation,
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
            raw_sum_final_stage=symbolica_raw_sum_final_stage,
        )
        self.build_options = DAGEvaluatorBuildOptions(
            merge_evaluators_strategy=merge_evaluators_strategy,
            split_vertex_current_stages=split_vertex_current_stages,
            verbose_evaluator_build=verbose_evaluator_build,
            symbolica_settings=symbolica_settings,
        )
        table_start = time.perf_counter()
        self.table = _build_shared_helicity_current_table(
            self.graph,
            self.model,
            gluon_count=self.gluon_count,
            progress_callback=progress_callback,
        )
        self.current_table_build_time_s = time.perf_counter() - table_start
        compiled_start = time.perf_counter()
        self.compiled: Any
        if symbolica_load_evaluator_dir is not None:
            if split_vertex_current_stages:
                raise NativeEvaluationError(
                    "compiled evaluator loading is currently implemented only for "
                    "the default shared-current D-mode, not split vertex/current stages"
                )
            loaded_manifest = _read_evaluator_artifact_manifest(
                Path(symbolica_load_evaluator_dir)
            )
            loaded_metadata = loaded_manifest.get("metadata")
            if isinstance(loaded_metadata, dict):
                self.loaded_evaluator_artifact_metadata = loaded_metadata
            self.compiled = _SharedCompiledSweepEvaluator.from_artifact(
                self.table,
                Path(symbolica_load_evaluator_dir),
                symbolica_settings=symbolica_settings,
                progress_callback=progress_callback,
            )
        elif split_vertex_current_stages:
            self.compiled = _SharedSplitCompiledSweepEvaluator(
                self.table,
                merge_evaluators_strategy=merge_evaluators_strategy,
                verbose_evaluator_build=verbose_evaluator_build,
                symbolica_settings=symbolica_settings,
                progress_callback=progress_callback,
            )
        else:
            self.compiled = _SharedCompiledSweepEvaluator(
                self.table,
                merge_evaluators_strategy=merge_evaluators_strategy,
                verbose_evaluator_build=verbose_evaluator_build,
                symbolica_settings=symbolica_settings,
                progress_callback=progress_callback,
            )
        self.symbolica_evaluator_build_time_s = time.perf_counter() - compiled_start
        self.blocks: _DAGBlockEvaluators | None = None
        self.block_build_time_s = 0.0

    def materialize_compiled_evaluators(self) -> None:
        _report_progress(
            self.progress_callback,
            stage="materialize",
            item="compiled evaluators",
        )
        materialize = getattr(self.compiled, "materialize", None)
        if callable(materialize):
            materialize()
        _report_progress(
            self.progress_callback,
            stage="materialize",
            item="compiled evaluators",
            increment=1,
        )

    def save_evaluator_artifact(self, output_dir: str | Path) -> Path:
        self.materialize_compiled_evaluators()
        output_path = Path(output_dir).expanduser()
        _prepare_process_output_directory(output_path)
        artifact = {
            "schema_version": 1,
            "kind": "pyamplicol-zgluon-shared-dag-compiled",
            "process": self.process,
            "gluon_count": self.gluon_count,
            "table_counts": {
                "currents": len(self.table.currents),
                "sources": len(self.table.sources),
                "interactions": len(self.table.interactions),
                "amplitudes": len(self.table.amplitudes),
            },
            "metadata": self.metadata.to_json_dict(),
            "compiled": self.compiled.artifact_manifest(output_path),
        }
        manifest_path = output_path / "manifest.json"
        manifest_path.write_text(
            json.dumps(artifact, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        _write_rusticol_process_artifacts(
            output_path,
            evaluator_manifest=artifact,
            evaluator_manifest_name=manifest_path.name,
            process=self.process,
            gluon_count=self.gluon_count,
            table=self.table,
            layout=self.compiled.layout,
            model=self.model,
        )
        _report_progress(
            self.progress_callback,
            stage="save",
            item="process artifacts",
            increment=1,
        )
        return manifest_path

    @property
    def metadata(self) -> DAGEvaluatorMetadata:
        symbolica_evaluator_settings = (
            self.build_options.symbolica_settings.to_json_dict()
        )
        if self.loaded_evaluator_artifact_metadata is not None:
            loaded_settings = self.loaded_evaluator_artifact_metadata.get(
                "symbolica_evaluator_settings"
            )
            if isinstance(loaded_settings, dict):
                symbolica_evaluator_settings = dict(loaded_settings)
                symbolica_evaluator_settings["loaded_from_artifact"] = True
        return DAGEvaluatorMetadata(
            process=self.process,
            kernel="spenso-symbolica-shared-helicity-current-dag",
            gluon_count=self.gluon_count,
            graph_current_count=len(self.graph.currents),
            graph_interaction_count=len(self.graph.interactions),
            graph_amplitude_count=len(self.graph.amplitudes),
            shared_current_count=len(self.table.currents),
            shared_source_current_count=len(self.table.sources),
            shared_interaction_count=len(self.table.interactions),
            shared_amplitude_count=len(self.table.amplitudes),
            current_table_build_time_s=self.current_table_build_time_s,
            block_build_time_s=self.block_build_time_s,
            symbolica_evaluator_build_time_s=self.symbolica_evaluator_build_time_s,
            symbolica_parameter_count=self.compiled.parameter_count,
            symbolica_output_count=len(self.table.amplitudes),
            batch_size=self.batch_size,
            merge_evaluators_strategy=self.build_options.merge_evaluators_strategy,
            split_vertex_current_stages=(
                self.build_options.split_vertex_current_stages
            ),
            verbose_evaluator_build=self.build_options.verbose_evaluator_build,
            symbolica_evaluator_settings=symbolica_evaluator_settings,
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
                _validate_z_gluon_point(
                    particles,
                    gluon_count=self.gluon_count,
                )
                for particles in points[start : start + self.batch_size]
            )
            evaluations.extend(self._evaluate_shared_batch(batch))
            timing += self.last_runtime_timing
        self.last_runtime_timing = timing
        return tuple(evaluations)

    def _evaluate_shared_point(
        self,
        point: tuple[ExternalMomentum, ...],
    ) -> MatrixElementEvaluation:
        return self._evaluate_shared_batch((point,))[0]

    def _evaluate_shared_batch(
        self,
        points: tuple[tuple[ExternalMomentum, ...], ...],
    ) -> tuple[MatrixElementEvaluation, ...]:
        amplitude_rows = self.compiled.evaluate_amplitude_rows(
            points,
            self.model,
            gluon_count=self.gluon_count,
        )
        reduction_start = time.perf_counter()
        evaluations = tuple(
            self._evaluation_from_amplitudes(point, amplitudes)
            for point, amplitudes in zip(points, amplitude_rows, strict=True)
        )
        self.last_runtime_timing = self.compiled.last_timing.with_result_reduction(
            time.perf_counter() - reduction_start
        )
        return evaluations

    def evaluate_matrix_elements_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> tuple[float, ...]:
        matrix_elements: list[float] = []
        timing = DAGEvaluationTiming()
        for start in range(0, len(points), self.batch_size):
            batch = tuple(
                _validate_z_gluon_point(
                    particles,
                    gluon_count=self.gluon_count,
                )
                for particles in points[start : start + self.batch_size]
            )
            matrix_elements.extend(self._evaluate_matrix_element_batch(batch))
            timing += self.last_runtime_timing
        self.last_runtime_timing = timing
        return tuple(matrix_elements)

    def _evaluate_matrix_element_batch(
        self,
        points: tuple[tuple[ExternalMomentum, ...], ...],
    ) -> tuple[float, ...]:
        raw_sums = self.compiled.evaluate_raw_sum_rows(
            points,
            self.model,
            gluon_count=self.gluon_count,
        )
        reduction_start = time.perf_counter()
        matrix_elements = tuple(
            self._matrix_element_from_raw_sum(point, raw_sum)
            for point, raw_sum in zip(points, raw_sums, strict=True)
        )
        self.last_runtime_timing = self.compiled.last_timing.with_result_reduction(
            time.perf_counter() - reduction_start
        )
        return matrix_elements

    def stage_diagnostics_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> dict[str, object]:
        batch = tuple(
            _validate_z_gluon_point(
                particles,
                gluon_count=self.gluon_count,
            )
            for particles in points
        )
        diagnostics = self.compiled.stage_diagnostics(
            batch,
            self.model,
            gluon_count=self.gluon_count,
        )
        self.last_runtime_timing = self.compiled.last_timing
        return diagnostics

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
        color_factor = self.model.leading_color_factor(pdgs)
        average_factor = _initial_state_average_factor(pdgs[:2])
        identical_factor = _final_state_identical_factor(
            particle.pdg for particle in point[2:]
        )
        coupling_factor = (
            (4.0 * math.pi * self.model.alpha_s_me_check) ** self.gluon_count
            * (2.0 * 4.0 * math.pi * self.model.alpha_ew)
        )
        matrix_element = (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
        )
        return MatrixElementEvaluation(
            process=self.process,
            particles=point,
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=color_factor,
            average_factor=average_factor,
            coupling_factor=coupling_factor,
            helicity_contributions=tuple(helicity_contributions),
            identical_factor=identical_factor,
        )

    def _matrix_element_from_raw_sum(
        self,
        point: tuple[ExternalMomentum, ...],
        raw_sum: float,
    ) -> float:
        pdgs = tuple(particle.pdg for particle in point)
        color_factor = self.model.leading_color_factor(pdgs)
        average_factor = _initial_state_average_factor(pdgs[:2])
        identical_factor = _final_state_identical_factor(
            particle.pdg for particle in point[2:]
        )
        coupling_factor = (
            (4.0 * math.pi * self.model.alpha_s_me_check) ** self.gluon_count
            * (2.0 * 4.0 * math.pi * self.model.alpha_ew)
        )
        return (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
        )

    def _helicity_amplitudes(
        self,
        particles: tuple[ExternalMomentum, ...],
    ) -> tuple[tuple[tuple[int, ...], complex], ...]:
        incoming_quark = particles[0]
        incoming_antiquark = particles[1]
        gluons = particles[2 : 2 + self.gluon_count]
        z_boson = particles[-1]
        anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
        quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
        gluon_momenta = tuple(gluon.momentum for gluon in gluons)
        z_momentum = z_boson.momentum

        anti_plus = _ext_antiquark_weyl(anti_closure_momentum, 1, -1)
        anti_minus = _ext_antiquark_weyl(anti_closure_momentum, -1, 1)
        quark_minus = _ext_quark_weyl(quark_start_momentum, -1, 1)
        quark_plus = _ext_quark_weyl(quark_start_momentum, 1, -1)
        z_vectors = {
            helicity: _ext_massive_vector(z_momentum, helicity, self.model.mass(23))
            for helicity in (-1, 0, 1)
        }
        amplitudes: list[tuple[tuple[int, ...], complex]] = []
        for (
            physical_quark_helicity,
            physical_antiquark_helicity,
            chirality,
            quark_wf,
            anti_wf,
        ) in (
            (1, -1, 1, quark_minus, anti_plus),
            (-1, 1, -1, quark_plus, anti_minus),
        ):
            for z_helicity in (-1, 0, 1):
                for gluon_helicities in product(
                    (-1, 1),
                    repeat=self.gluon_count,
                ):
                    gluon_wfs = tuple(
                        _ext_gluon_cmplx(momentum, helicity)
                        for momentum, helicity in zip(
                            gluon_momenta,
                            gluon_helicities,
                            strict=True,
                        )
                    )
                    currents = self._evaluate_currents(
                        particles,
                        chirality=chirality,
                        quark_wf=quark_wf,
                        anti_wf=anti_wf,
                        gluon_wfs=gluon_wfs,
                        z_wf=z_vectors[z_helicity],
                    )
                    amplitudes.extend(
                        (
                            (
                                physical_quark_helicity,
                                physical_antiquark_helicity,
                                *gluon_helicities,
                                z_helicity,
                            ),
                            _dot_weyl(
                                currents[_current_key_tuple(left)],
                                currents[_current_key_tuple(right)],
                            ),
                        )
                        for left, right in self.graph.amplitudes
                        if int(left.chirality) == chirality
                    )
        return tuple(amplitudes)

    def _evaluate_currents(
        self,
        particles: tuple[ExternalMomentum, ...],
        *,
        chirality: int,
        quark_wf: tuple[complex, complex],
        anti_wf: tuple[complex, complex],
        gluon_wfs: tuple[tuple[complex, complex, complex, complex], ...],
        z_wf: tuple[complex, complex, complex, complex],
    ) -> dict[tuple[int, tuple[int, ...], int], tuple[complex, ...]]:
        currents: dict[tuple[int, tuple[int, ...], int], tuple[complex, ...]] = {}
        for current in _source_currents(self.graph):
            currents[_current_key_tuple(current)] = tuple(
                0j for _ in range(_current_dimension(current))
            )
        if self.blocks is None:
            self.blocks = _DAGBlockEvaluators(self.model)

        incoming_quark = particles[0]
        incoming_antiquark = particles[1]
        currents[_current_key_tuple(CurrentKey(-incoming_quark.pdg, (1,), 0))] = anti_wf
        currents[
            _current_key_tuple(CurrentKey(-incoming_antiquark.pdg, (2,), chirality))
        ] = quark_wf
        for offset, gluon_wf in enumerate(gluon_wfs, start=3):
            currents[_current_key_tuple(CurrentKey(21, (offset,), 0))] = gluon_wf
        currents[
            _current_key_tuple(CurrentKey(23, (self.gluon_count + 3,), 0))
        ] = z_wf

        momenta_by_label = _current_momenta_by_label(particles)
        interactions_by_result: dict[tuple[int, tuple[int, ...], int], list[Any]] = {}
        for interaction in self.graph.interactions:
            interactions_by_result.setdefault(
                _current_key_tuple(interaction.result),
                [],
            ).append(interaction)

        for current in self.graph.currents:
            current_key = _current_key_tuple(current)
            interactions = interactions_by_result.get(current_key)
            if not interactions:
                continue
            if int(current.chirality) not in (0, chirality):
                continue
            total = tuple(0j for _ in range(_current_dimension(current)))
            for interaction in interactions:
                left = currents[_current_key_tuple(interaction.left)]
                right = currents[_current_key_tuple(interaction.right)]
                contribution = self.blocks.vertex(
                    interaction,
                    left,
                    right,
                    momenta_by_label,
                )
                total = _sum_components(total, contribution)
            if _current_needs_propagator(self.graph, current):
                total = self.blocks.propagate(
                    current,
                    total,
                    _sum_momenta(
                        momenta_by_label[label] for label in current.external_labels
                    ),
                )
            currents[current_key] = total
        return currents


class _SharedCompiledSweepEvaluator:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
        progress_callback: ProgressCallback | None = None,
    ) -> None:
        self.table = table
        self.layout = _build_shared_global_parameter_layout(table)
        self._init_common(table)
        stage_id_groups = _shared_current_stage_ids(table)
        stages: list[_CompiledCurrentStage] = []
        for stage_index, current_ids in enumerate(stage_id_groups, start=1):
            _report_progress(
                progress_callback,
                stage="compile",
                item=f"stage {stage_index}/{len(stage_id_groups)}",
            )
            stages.append(
                _CompiledCurrentStage(
                    table,
                    layout=self.layout,
                    current_ids=current_ids,
                    stage_index=stage_index,
                    merge_evaluators_strategy=merge_evaluators_strategy,
                    verbose_evaluator_build=verbose_evaluator_build,
                    symbolica_settings=symbolica_settings,
                )
            )
            _report_progress(
                progress_callback,
                stage="compile",
                item=f"stage {stage_index}/{len(stage_id_groups)}",
                increment=1,
            )
        self.stages = tuple(stages)
        _report_progress(
            progress_callback,
            stage="compile",
            item="amplitude stage",
        )
        self.amplitude_stage = _CompiledAmplitudeStage(
            table,
            layout=self.layout,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            symbolica_settings=symbolica_settings,
        )
        _report_progress(
            progress_callback,
            stage="compile",
            item="amplitude stage",
            increment=1,
        )
        self.output_length = len(table.amplitudes)
        self.parameter_count = self.layout.parameter_count
        self.last_evaluator_time_s = 0.0
        self.last_timing = DAGEvaluationTiming()

    def _init_common(self, table: SharedCurrentTable) -> None:
        self._momentum_offsets_and_labels = _unique_momentum_offsets_and_labels(
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
        self._momentum_column_slice = _contiguous_column_slice(
            self._momentum_flat_columns,
        )
        for row, (_, labels) in enumerate(self._momentum_offsets_and_labels):
            for label in labels:
                self._momentum_label_matrix[row, label - 1] = 1.0
        self.zero_current_values = np.zeros(
            (len(table.currents), 6),
            dtype=np.complex128,
        )

    @classmethod
    def from_artifact(
        cls,
        table: SharedCurrentTable,
        artifact_dir: Path,
        *,
        symbolica_settings: SymbolicaEvaluatorSettings,
        progress_callback: ProgressCallback | None = None,
    ) -> "_SharedCompiledSweepEvaluator":
        manifest = _read_evaluator_artifact_manifest(artifact_dir)
        if manifest.get("kind") != "pyamplicol-zgluon-shared-dag-compiled":
            raise NativeEvaluationError(
                f"unsupported evaluator artifact kind: {manifest.get('kind')!r}"
            )
        compiled_manifest = manifest.get("compiled")
        if not isinstance(compiled_manifest, dict):
            raise NativeEvaluationError("compiled evaluator artifact is missing")
        if compiled_manifest.get("kind") != "shared-compiled-sweep":
            raise NativeEvaluationError(
                "compiled evaluator artifact does not contain a default "
                "shared-current sweep"
            )
        instance = cls.__new__(cls)
        instance.table = table
        instance.layout = _build_shared_global_parameter_layout(table)
        instance._init_common(table)
        stage_manifests = compiled_manifest.get("stages")
        if not isinstance(stage_manifests, list):
            raise NativeEvaluationError("compiled evaluator artifact has no stages")
        stages: list[_CompiledCurrentStage] = []
        for index, stage_manifest in enumerate(stage_manifests, start=1):
            _report_progress(
                progress_callback,
                stage="load",
                item=f"stage {index}/{len(stage_manifests)}",
            )
            stages.append(
                _CompiledCurrentStage.from_artifact(
                    table,
                    layout=instance.layout,
                    manifest=stage_manifest,
                    artifact_dir=artifact_dir,
                )
            )
            _report_progress(
                progress_callback,
                stage="load",
                item=f"stage {index}/{len(stage_manifests)}",
                increment=1,
            )
        instance.stages = tuple(stages)
        amplitude_manifest = compiled_manifest.get("amplitude_stage")
        if not isinstance(amplitude_manifest, dict):
            raise NativeEvaluationError(
                "compiled evaluator artifact has no amplitude stage"
            )
        _report_progress(progress_callback, stage="load", item="amplitude stage")
        instance.amplitude_stage = _CompiledAmplitudeStage.from_artifact(
            table,
            layout=instance.layout,
            manifest=amplitude_manifest,
            artifact_dir=artifact_dir,
            symbolica_settings=symbolica_settings,
        )
        _report_progress(
            progress_callback,
            stage="load",
            item="amplitude stage",
            increment=1,
        )
        instance.output_length = len(table.amplitudes)
        instance.parameter_count = instance.layout.parameter_count
        instance.last_evaluator_time_s = 0.0
        instance.last_timing = DAGEvaluationTiming()
        return instance

    def materialize(self) -> None:
        self.amplitude_stage.materialize()

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        return {
            "kind": "shared-compiled-sweep",
            "parameter_count": self.parameter_count,
            "output_length": self.output_length,
            "stages": [
                stage.artifact_manifest(artifact_dir)
                for stage in self.stages
            ],
            "amplitude_stage": self.amplitude_stage.artifact_manifest(artifact_dir),
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
        state_rows, current_rows, timing = self._source_current_rows(
            points,
            model,
            gluon_count=gluon_count,
        )
        for stage in self.stages:
            timing += stage.evaluate_and_assign(
                state_rows,
                current_rows,
                self.table,
            )
        amplitude_rows, amplitude_timing = self.amplitude_stage.evaluate_rows(
            state_rows,
            current_rows
        )
        timing += amplitude_timing
        self.last_evaluator_time_s = timing.evaluator_time_s
        self.last_timing = timing
        return amplitude_rows

    def evaluate_raw_sum_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[float, ...]:
        if not points:
            return ()
        state_rows, current_rows, timing = self._source_current_rows(
            points,
            model,
            gluon_count=gluon_count,
        )
        for stage in self.stages:
            timing += stage.evaluate_and_assign(
                state_rows,
                current_rows,
                self.table,
            )
        raw_sums, raw_sum_timing = self.amplitude_stage.evaluate_raw_sums(
            state_rows,
            current_rows
        )
        timing += raw_sum_timing
        self.last_evaluator_time_s = timing.evaluator_time_s
        self.last_timing = timing
        return raw_sums

    def stage_diagnostics(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> dict[str, object]:
        state_rows, current_rows, timing = self._source_current_rows(
            points,
            model,
            gluon_count=gluon_count,
        )
        stages: list[dict[str, object]] = [
            {
                "stage": "sources",
                "output_len": int(current_rows.size),
                **_complex_array_checksum(current_rows),
            }
        ]
        for stage_index, stage in enumerate(self.stages, start=1):
            timing += stage.evaluate_and_assign(
                state_rows,
                current_rows,
                self.table,
            )
            values = _selected_current_components(
                current_rows,
                stage.output_current_ids,
                stage.output_current_components,
            )
            stages.append(
                {
                    "stage": stage_index,
                    "output_len": int(values.size),
                    **_complex_array_checksum(values),
                }
            )
        self.last_evaluator_time_s = timing.evaluator_time_s
        self.last_timing = timing
        return {
            "points": len(points),
            "stages": stages,
            "timing": timing.to_json_dict(),
        }

    def _source_current_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        DAGEvaluationTiming,
    ]:
        state_rows = np.zeros(
            (len(points), self.layout.parameter_count),
            dtype=np.complex128,
        )
        current_rows = state_rows[
            :,
            : 6 * len(self.table.currents),
        ].reshape((len(points), len(self.table.currents), 6))
        source_start = time.perf_counter()
        for row_index, point in enumerate(points):
            _fill_shared_source_currents(
                current_rows[row_index],
                self.table,
                point,
                model,
                gluon_count=gluon_count,
                zero_current_values=self.zero_current_values,
                clear_values=False,
            )
        source_fill_time_s = time.perf_counter() - source_start
        momentum_start = time.perf_counter()
        if self._momentum_offsets_and_labels:
            external_momenta = np.empty(
                (len(points), self._external_label_count, 4),
                dtype=np.float64,
            )
            for row_index, point in enumerate(points):
                for label_index, particle in enumerate(point):
                    external_momenta[row_index, label_index, :] = particle.momentum
            external_momenta *= self._external_momentum_signs
            momentum_values = np.einsum(
                "ml,blc->bmc",
                self._momentum_label_matrix,
                external_momenta,
                optimize=True,
            )
            flat_momentum_values = momentum_values.reshape(len(points), -1)
            if self._momentum_column_slice is None:
                state_rows[:, self._momentum_flat_columns] = flat_momentum_values
            else:
                state_rows[:, self._momentum_column_slice] = flat_momentum_values
        momentum_setup_time_s = time.perf_counter() - momentum_start
        return (
            state_rows,
            current_rows,
            DAGEvaluationTiming(
                source_fill_time_s=source_fill_time_s,
                momentum_setup_time_s=momentum_setup_time_s,
            ),
        )

    def evaluate_amplitudes(
        self,
        point: tuple[ExternalMomentum, ...],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[complex, ...]:
        return self.evaluate_amplitude_rows(
            (point,),
            model,
            gluon_count=gluon_count,
        )[0]


class _SharedSplitCompiledSweepEvaluator:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
        progress_callback: ProgressCallback | None = None,
    ) -> None:
        self.table = table
        self.layout = _build_shared_split_global_parameter_layout(table)
        self._momentum_offsets_and_labels = _unique_momentum_offsets_and_labels(
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
        self._momentum_column_slice = _contiguous_column_slice(
            self._momentum_flat_columns,
        )
        for row, (_, labels) in enumerate(self._momentum_offsets_and_labels):
            for label in labels:
                self._momentum_label_matrix[row, label - 1] = 1.0
        self.zero_current_values = np.zeros(
            (len(table.currents), 6),
            dtype=np.complex128,
        )
        stage_id_groups = _shared_current_stage_ids(table)
        stages: list[_CompiledSplitCurrentStage] = []
        for stage_index, current_ids in enumerate(stage_id_groups, start=1):
            _report_progress(
                progress_callback,
                stage="compile",
                item=f"vertices {stage_index}/{len(stage_id_groups)}",
            )
            stage = _CompiledSplitCurrentStage(
                table,
                layout=self.layout,
                current_ids=current_ids,
                stage_index=stage_index,
                stage_count=len(stage_id_groups),
                merge_evaluators_strategy=merge_evaluators_strategy,
                verbose_evaluator_build=verbose_evaluator_build,
                symbolica_settings=symbolica_settings,
                progress_callback=progress_callback,
            )
            stages.append(stage)
        self.stages = tuple(stages)
        _report_progress(progress_callback, stage="compile", item="amplitude stage")
        self.amplitude_stage = _CompiledAmplitudeStage(
            table,
            layout=self.layout,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            symbolica_settings=symbolica_settings,
        )
        _report_progress(
            progress_callback,
            stage="compile",
            item="amplitude stage",
            increment=1,
        )
        self.output_length = len(table.amplitudes)
        self.parameter_count = self.layout.parameter_count
        self.last_evaluator_time_s = 0.0
        self.last_timing = DAGEvaluationTiming()

    def evaluate_amplitude_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[tuple[complex, ...], ...]:
        if not points:
            return ()
        state_rows, current_rows, interaction_rows, timing = self._source_current_rows(
            points,
            model,
            gluon_count=gluon_count,
        )
        for stage in self.stages:
            timing += stage.evaluate_and_assign(
                state_rows,
                current_rows,
                interaction_rows,
                self.table,
            )
        amplitude_rows, amplitude_timing = self.amplitude_stage.evaluate_rows(
            state_rows,
            current_rows,
        )
        timing += amplitude_timing
        self.last_evaluator_time_s = timing.evaluator_time_s
        self.last_timing = timing
        return amplitude_rows

    def evaluate_raw_sum_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[float, ...]:
        if not points:
            return ()
        state_rows, current_rows, interaction_rows, timing = self._source_current_rows(
            points,
            model,
            gluon_count=gluon_count,
        )
        for stage in self.stages:
            timing += stage.evaluate_and_assign(
                state_rows,
                current_rows,
                interaction_rows,
                self.table,
            )
        raw_sums, raw_sum_timing = self.amplitude_stage.evaluate_raw_sums(
            state_rows,
            current_rows,
        )
        timing += raw_sum_timing
        self.last_evaluator_time_s = timing.evaluator_time_s
        self.last_timing = timing
        return raw_sums

    def _source_current_rows(
        self,
        points: Sequence[tuple[ExternalMomentum, ...]],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        DAGEvaluationTiming,
    ]:
        state_rows = np.zeros(
            (len(points), self.layout.parameter_count),
            dtype=np.complex128,
        )
        current_rows = state_rows[
            :,
            : 6 * len(self.table.currents),
        ].reshape((len(points), len(self.table.currents), 6))
        interaction_start = self.layout.interaction_offsets[0]
        interaction_rows = state_rows[
            :,
            interaction_start : interaction_start + 6 * len(self.table.interactions),
        ].reshape((len(points), len(self.table.interactions), 6))
        source_start = time.perf_counter()
        for row_index, point in enumerate(points):
            _fill_shared_source_currents(
                current_rows[row_index],
                self.table,
                point,
                model,
                gluon_count=gluon_count,
                zero_current_values=self.zero_current_values,
                clear_values=False,
            )
        source_fill_time_s = time.perf_counter() - source_start
        momentum_start = time.perf_counter()
        if self._momentum_offsets_and_labels:
            external_momenta = np.empty(
                (len(points), self._external_label_count, 4),
                dtype=np.float64,
            )
            for row_index, point in enumerate(points):
                for label_index, particle in enumerate(point):
                    external_momenta[row_index, label_index, :] = particle.momentum
            external_momenta *= self._external_momentum_signs
            momentum_values = np.einsum(
                "ml,blc->bmc",
                self._momentum_label_matrix,
                external_momenta,
                optimize=True,
            )
            flat_momentum_values = momentum_values.reshape(len(points), -1)
            if self._momentum_column_slice is None:
                state_rows[:, self._momentum_flat_columns] = flat_momentum_values
            else:
                state_rows[:, self._momentum_column_slice] = flat_momentum_values
        momentum_setup_time_s = time.perf_counter() - momentum_start
        return (
            state_rows,
            current_rows,
            interaction_rows,
            DAGEvaluationTiming(
                source_fill_time_s=source_fill_time_s,
                momentum_setup_time_s=momentum_setup_time_s,
            ),
        )

    def evaluate_amplitudes(
        self,
        point: tuple[ExternalMomentum, ...],
        model: AmplicolSMLeadingColorModel,
        *,
        gluon_count: int,
    ) -> tuple[complex, ...]:
        return self.evaluate_amplitude_rows(
            (point,),
            model,
            gluon_count=gluon_count,
        )[0]


def _current_slot_index_arrays(
    slots: Sequence[tuple[int, int, int]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    columns: list[int] = []
    current_ids: list[int] = []
    components: list[int] = []
    for current_id, start, stop in slots:
        for component, column in enumerate(range(start, stop)):
            columns.append(column)
            current_ids.append(current_id)
            components.append(component)
    return (
        np.asarray(columns, dtype=np.intp),
        np.asarray(current_ids, dtype=np.intp),
        np.asarray(components, dtype=np.intp),
    )


def _output_slot_index_arrays(
    slots: Sequence[tuple[int, tuple[int, int]]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    columns: list[int] = []
    current_ids: list[int] = []
    components: list[int] = []
    for current_id, (start, stop) in slots:
        for component, column in enumerate(range(start, stop)):
            columns.append(column)
            current_ids.append(current_id)
            components.append(component)
    return (
        np.asarray(columns, dtype=np.intp),
        np.asarray(current_ids, dtype=np.intp),
        np.asarray(components, dtype=np.intp),
    )


def _output_slots_to_manifest(
    slots: Sequence[tuple[int, tuple[int, int]]],
) -> list[dict[str, int]]:
    return [
        {
            "id": int(output_id),
            "start": int(bounds[0]),
            "stop": int(bounds[1]),
        }
        for output_id, bounds in slots
    ]


def _output_slices_from_manifest(
    manifest: dict[str, Any],
) -> dict[int, tuple[int, int]]:
    slots = manifest.get("output_slots")
    if not isinstance(slots, list):
        raise NativeEvaluationError("compiled evaluator stage is missing output slots")
    return {
        int(slot["id"]): (int(slot["start"]), int(slot["stop"]))
        for slot in slots
    }


def _read_evaluator_artifact_manifest(artifact_dir: Path) -> dict[str, Any]:
    manifest_path = artifact_dir.expanduser() / "manifest.json"
    if not manifest_path.exists():
        raise NativeEvaluationError(
            f"compiled evaluator artifact manifest does not exist: {manifest_path}"
        )
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def _prepare_process_output_directory(output_path: Path) -> None:
    output_path.mkdir(parents=True, exist_ok=True)
    for stale_file in output_path.glob("*.evaluator.bin"):
        stale_file.unlink()
    for stale_pattern in ("pyamplicol_*.cpp", "libpyamplicol_*"):
        for stale_file in output_path.glob(stale_pattern):
            if stale_file.is_file() or stale_file.is_symlink():
                stale_file.unlink()
    pycache = output_path / "__pycache__"
    if pycache.exists():
        shutil.rmtree(pycache)


def _rusticol_backend_settings(
    evaluator_manifest: dict[str, Any],
    *,
    evaluator_manifest_name: str | None,
) -> dict[str, Any]:
    metadata = evaluator_manifest.get("metadata")
    evaluator_settings = (
        metadata.get("symbolica_evaluator_settings")
        if isinstance(metadata, dict)
        else None
    )
    compiled = evaluator_manifest.get("compiled")
    return {
        "runtime": "rusticol",
        "compiled_kind": (
            compiled.get("kind") if isinstance(compiled, dict) else None
        ),
        "evaluator_manifest": evaluator_manifest_name,
        "symbolica_evaluator_settings": (
            evaluator_settings if isinstance(evaluator_settings, dict) else {}
        ),
    }


def _rusticol_dependency_fingerprint() -> dict[str, Any]:
    fingerprint: dict[str, Any] = {
        "pyamplicol_version": _installed_package_version("pyamplicol"),
        "rusticol_version": _installed_package_version("rusticol"),
        "symbolica": _symbolica_runtime_fingerprint(),
    }
    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "dependencies"
        / "install_manifest.json"
    )
    if manifest_path.exists():
        try:
            install_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            fingerprint["install_manifest_error"] = str(exc)
        else:
            if isinstance(install_manifest, dict):
                fingerprint["install_manifest"] = {
                    key: install_manifest[key]
                    for key in (
                        "schema_version",
                        "dependency_patches",
                        "rusticol",
                        "symbolica",
                        "symbolica_community",
                    )
                    if key in install_manifest
                }
    return fingerprint


def _installed_package_version(name: str) -> str | None:
    try:
        return importlib_metadata.version(name)
    except importlib_metadata.PackageNotFoundError:
        return None


def _symbolica_runtime_fingerprint() -> dict[str, Any]:
    try:
        import symbolica
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
    }


def _rusticol_helicity_filter_manifest(table: SharedCurrentTable) -> dict[str, Any]:
    return {
        "strategy": "shared-current-retained-amplitudes",
        "retained_amplitude_count": len(table.amplitudes),
        "amplitudes": [
            {
                "index": index,
                "left_current_id": amplitude.left_id,
                "right_current_id": amplitude.right_id,
                "helicities": list(amplitude.helicities),
                "multiplicity": amplitude.multiplicity,
                "raw_sum_weight": float(amplitude.multiplicity),
            }
            for index, amplitude in enumerate(table.amplitudes)
        ],
        "raw_sum_weights": [
            float(amplitude.multiplicity)
            for amplitude in table.amplitudes
        ],
    }


def _rusticol_external_particles_manifest(
    point: Sequence[ExternalMomentum],
) -> list[dict[str, Any]]:
    return [
        {
            "label": index + 1,
            "index": index,
            "pdg": int(particle.pdg),
            "role": "initial" if index < 2 else "final",
            "momentum_slot": index,
            "momentum_components": ["E", "px", "py", "pz"],
        }
        for index, particle in enumerate(point)
    ]


def _rusticol_momentum_conventions_manifest(
    point: Sequence[ExternalMomentum],
) -> dict[str, Any]:
    incoming_labels = [
        index + 1
        for index in range(min(2, len(point)))
    ]
    final_labels = [
        index + 1
        for index in range(2, len(point))
    ]
    return {
        "input_shape": ["batch", len(point), 4],
        "component_order": ["E", "px", "py", "pz"],
        "input_momenta": "physical external four-momenta in process order",
        "incoming_labels": incoming_labels,
        "final_state_labels": final_labels,
        "source_current_crossing": {
            "crossed_incoming_labels": incoming_labels,
            "operation": "negate four-vector when filling source currents",
        },
        "metric": "mostly-minus",
    }


def _rusticol_model_manifest(model: AmplicolSMLeadingColorModel) -> dict[str, Any]:
    return {
        "name": model.name,
        "alpha_s_me_check": model.alpha_s_me_check,
        "alpha_ew": model.alpha_ew,
        "mass_z": model.mass(23),
        "sqrt_s_default": model.sqrt_s,
        "particles": [
            {
                "pdg": particle.pdg,
                "anti_pdg": particle.anti_pdg,
                "spin": particle.spin,
                "dimension": particle.dimension,
                "color_rep": particle.color_rep,
                "mass": particle.mass,
                "width": particle.width,
                "charge": particle.charge,
                "weak_isospin": list(particle.weak_isospin),
                "weak_hypercharge": list(particle.weak_hypercharge),
            }
            for particle in sorted(
                model.particles.values(),
                key=lambda item: item.pdg,
            )
        ],
        "vertices": [
            {
                "kind": vertex.kind,
                "particles": list(vertex.particles),
                "coupling": list(vertex.coupling),
            }
            for vertex in model.vertices
        ],
    }


def _write_rusticol_process_artifacts(
    output_path: Path,
    *,
    evaluator_manifest: dict[str, Any],
    evaluator_manifest_name: str,
    process: str,
    gluon_count: int,
    table: SharedCurrentTable,
    layout: _SharedGlobalParameterLayout,
    model: AmplicolSMLeadingColorModel,
) -> None:
    validation_points = _rusticol_validation_points(
        process,
        gluon_count=gluon_count,
        model=model,
    )
    process_manifest = {
        "schema_version": 1,
        "kind": "pyamplicol-rusticol-process",
        "process": process,
        "family": "q-qbar-z-gluons-leading-color",
        "gluon_count": gluon_count,
        "external_pdg_order": [
            int(particle.pdg)
            for particle in validation_points[0]
        ],
        "external_particles": _rusticol_external_particles_manifest(
            validation_points[0]
        ),
        "momentum_conventions": _rusticol_momentum_conventions_manifest(
            validation_points[0]
        ),
        "model": _rusticol_model_manifest(model),
        "normalization": {
            "color_factor": model.leading_color_factor(
                particle.pdg for particle in validation_points[0]
            ),
            "average_factor": _initial_state_average_factor(
                particle.pdg for particle in validation_points[0][:2]
            ),
            "identical_factor": _final_state_identical_factor(
                particle.pdg for particle in validation_points[0][2:]
            ),
            "coupling_factor": (
                (4.0 * math.pi * model.alpha_s_me_check) ** gluon_count
                * (2.0 * 4.0 * math.pi * model.alpha_ew)
            ),
        },
        "helicity_filter": _rusticol_helicity_filter_manifest(table),
        "layout": _rusticol_layout_manifest(table, layout),
        "table": _rusticol_table_manifest(table),
        "evaluator_manifest": evaluator_manifest_name,
        "compiled": evaluator_manifest["compiled"],
        "backend_settings": _rusticol_backend_settings(
            evaluator_manifest,
            evaluator_manifest_name=evaluator_manifest_name,
        ),
        "dependency_fingerprint": _rusticol_dependency_fingerprint(),
        "metadata": evaluator_manifest.get("metadata", {}),
    }
    (output_path / "process_manifest.json").write_text(
        json.dumps(process_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (output_path / "validation_momenta.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pyamplicol-rusticol-validation-momenta",
                "process": process,
                "points": [
                    [
                        {
                            "pdg": int(particle.pdg),
                            "momentum": [
                                _decimal_string(component)
                                for component in particle.momentum
                            ],
                        }
                        for particle in point
                    ]
                    for point in validation_points
                ],
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    (output_path / "check_standalone.py").write_text(
        _RUSTICOL_STANDALONE_CHECK_SCRIPT,
        encoding="utf-8",
    )


def _write_zero_gluon_rusticol_process_artifacts(
    output_path: Path,
    *,
    process: str,
    model: AmplicolSMLeadingColorModel,
) -> Path:
    from .symbolic import ZeroGluonSymbolicEvaluator

    output_path.mkdir(parents=True, exist_ok=True)
    validation_points = _rusticol_validation_points(
        process,
        gluon_count=0,
        model=model,
    )
    symbolic = ZeroGluonSymbolicEvaluator(model=model).metadata.to_json_dict()
    parameter_names = list(symbolic["parameter_names"])
    state_dir = output_path / "zero_gluon"
    state_dir.mkdir(parents=True, exist_ok=True)
    state_path = state_dir / "evaluator.bin"
    state_path.write_bytes(base64.b64decode(str(symbolic["evaluator_state_b64"])))
    pdgs = tuple(particle.pdg for particle in validation_points[0])
    z_left, z_right = model.z_fermion_coupling(abs(pdgs[0]))
    process_manifest = {
        "schema_version": 1,
        "kind": "pyamplicol-rusticol-process",
        "process": process,
        "family": "q-qbar-z-gluons-leading-color",
        "gluon_count": 0,
        "external_pdg_order": [int(particle.pdg) for particle in validation_points[0]],
        "external_particles": _rusticol_external_particles_manifest(
            validation_points[0]
        ),
        "momentum_conventions": _rusticol_momentum_conventions_manifest(
            validation_points[0]
        ),
        "model": _rusticol_model_manifest(model),
        "normalization": {
            "color_factor": model.leading_color_factor(pdgs),
            "average_factor": _initial_state_average_factor(pdgs[:2]),
            "identical_factor": _final_state_identical_factor(
                particle.pdg for particle in validation_points[0][2:]
            ),
            "coupling_factor": 2.0 * 4.0 * math.pi * model.alpha_ew,
        },
        "helicity_filter": {
            "strategy": "zero-gluon-symbolic-scalar",
            "retained_amplitude_count": 1,
            "amplitudes": [],
            "raw_sum_weights": [1.0],
        },
        "layout": {
            "parameter_count": len(parameter_names),
            "current_offsets": [],
            "momentum_offsets": {},
            "real_valued_inputs": list(range(len(parameter_names))),
            "momentum_offsets_and_labels": [],
        },
        "table": {
            "currents": [],
            "sources": [],
            "interactions": [],
            "interactions_by_result": [],
            "amplitudes": [],
        },
        "evaluator_manifest": None,
        "compiled": {
            "kind": "zero-gluon-symbolic-scalar",
            "stages": [],
            "amplitude_stage": {
                "kind": "amplitude-stage",
                "output_length": 1,
                "raw_sum_weights": [1.0],
                "amplitude_evaluator": None,
                "raw_sum_evaluator": None,
            },
            "zero_gluon": {
                "parameter_names": parameter_names,
                "evaluator_state_path": _artifact_path_for_manifest(
                    state_path,
                    output_path,
                ),
                "z_left": z_left,
                "z_right": z_right,
            },
        },
        "backend_settings": {
            "runtime": "rusticol",
            "compiled_kind": "zero-gluon-symbolic-scalar",
            "evaluator_manifest": None,
            "symbolica_evaluator_settings": {
                "backend": "symbolica-evaluator-state",
            },
        },
        "dependency_fingerprint": _rusticol_dependency_fingerprint(),
        "metadata": {
            "kernel": "symbolica-zero-gluon",
            "symbolica_evaluator_settings": {
                "backend": "symbolica-evaluator-state",
            },
            "symbolic_scalar_evaluator": {
                key: value
                for key, value in symbolic.items()
                if key != "evaluator_state_b64"
            },
        },
    }
    manifest_path = output_path / "process_manifest.json"
    manifest_path.write_text(
        json.dumps(process_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (output_path / "validation_momenta.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pyamplicol-rusticol-validation-momenta",
                "process": process,
                "points": [
                    [
                        {
                            "pdg": int(particle.pdg),
                            "momentum": [
                                _decimal_string(component)
                                for component in particle.momentum
                            ],
                        }
                        for particle in point
                    ]
                    for point in validation_points
                ],
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    (output_path / "check_standalone.py").write_text(
        _RUSTICOL_STANDALONE_CHECK_SCRIPT,
        encoding="utf-8",
    )
    return manifest_path


def _rusticol_validation_points(
    process: str,
    *,
    gluon_count: int,
    model: AmplicolSMLeadingColorModel,
) -> tuple[tuple[ExternalMomentum, ...], ...]:
    native = LeadingColorZJetsNativeEvaluator(model)
    if gluon_count == 0:
        return (
            native.canonical_zero_gluon_point(
                process,
                sqrt_s=1000.0,
            ),
        )
    return (
        native.canonical_z_gluon_point(
            process,
            gluon_count=gluon_count,
            sqrt_s=1000.0,
        ),
    )


def _rusticol_layout_manifest(
    table: SharedCurrentTable,
    layout: _SharedGlobalParameterLayout,
) -> dict[str, Any]:
    return {
        "parameter_count": layout.parameter_count,
        "current_offsets": list(layout.current_offsets),
        "momentum_offsets": {
            str(current_id): int(offset)
            for current_id, offset in sorted(layout.momentum_offsets.items())
        },
        "real_valued_inputs": list(layout.real_valued_inputs),
        "momentum_offsets_and_labels": [
            {
                "offset": int(offset),
                "labels": list(labels),
            }
            for offset, labels in _unique_momentum_offsets_and_labels(table, layout)
        ],
    }


def _rusticol_table_manifest(table: SharedCurrentTable) -> dict[str, Any]:
    return {
        "currents": [
            {
                "id": current.id,
                "pdg": current.pdg,
                "external_labels": list(current.external_labels),
                "chirality": current.chirality,
                "ext_source_bits": current.ext_source_bits,
                "source_ids": list(current.source_ids),
                "is_source": current.is_source,
                "needs_propagator": current.needs_propagator,
                "dimension": current.dimension,
            }
            for current in table.currents
        ],
        "sources": [
            {
                "current_id": source.current_id,
                "leg_label": source.leg_label,
                "helicity": source.helicity,
                "physical_helicity": source.physical_helicity,
                "chirality": source.chirality,
                "source_bit": source.source_bit,
            }
            for source in table.sources
        ],
        "interactions": [
            {
                "id": interaction.id,
                "vertex_kind": interaction.vertex_kind,
                "left_id": interaction.left_id,
                "right_id": interaction.right_id,
                "result_id": interaction.result_id,
                "coupling": list(interaction.coupling),
            }
            for interaction in table.interactions
        ],
        "interactions_by_result": [
            list(interaction_ids)
            for interaction_ids in table.interactions_by_result
        ],
        "amplitudes": [
            {
                "left_id": amplitude.left_id,
                "right_id": amplitude.right_id,
                "helicities": list(amplitude.helicities),
                "multiplicity": amplitude.multiplicity,
            }
            for amplitude in table.amplitudes
        ],
    }


def _shared_current_table_from_rusticol_manifest(
    manifest: dict[str, Any],
) -> SharedCurrentTable:
    table_payload = manifest.get("table")
    if not isinstance(table_payload, dict):
        raise NativeEvaluationError("process manifest is missing current-table metadata")

    current_payloads = table_payload.get("currents")
    if not isinstance(current_payloads, list):
        raise NativeEvaluationError("process manifest table is missing currents")
    currents: list[SharedCurrentNode] = []
    for expected_id, payload in enumerate(current_payloads):
        if not isinstance(payload, dict):
            raise NativeEvaluationError("current table entry is not a JSON object")
        current_id = int(payload["id"])
        if current_id != expected_id:
            raise NativeEvaluationError(
                "current table ids must be contiguous and ordered in process manifest"
            )
        key = CurrentKey(
            int(payload["pdg"]),
            tuple(int(label) for label in payload["external_labels"]),
            int(payload["chirality"]),
        )
        currents.append(
            SharedCurrentNode(
                id=current_id,
                key=key,
                ext_source_bits=int(payload["ext_source_bits"]),
                source_ids=tuple(int(source_id) for source_id in payload["source_ids"]),
                is_source=bool(payload["is_source"]),
                needs_propagator=bool(payload["needs_propagator"]),
                dimension=int(payload["dimension"]),
            )
        )

    source_payloads = table_payload.get("sources")
    if not isinstance(source_payloads, list):
        raise NativeEvaluationError("process manifest table is missing sources")
    sources = tuple(
        SharedSourceCurrent(
            current_id=int(payload["current_id"]),
            leg_label=int(payload["leg_label"]),
            helicity=int(payload["helicity"]),
            physical_helicity=int(payload["physical_helicity"]),
            chirality=int(payload["chirality"]),
            source_bit=int(payload["source_bit"]),
        )
        for payload in source_payloads
        if isinstance(payload, dict)
    )
    if len(sources) != len(source_payloads):
        raise NativeEvaluationError("source table entry is not a JSON object")

    interaction_payloads = table_payload.get("interactions")
    if not isinstance(interaction_payloads, list):
        raise NativeEvaluationError("process manifest table is missing interactions")
    interactions: list[SharedInteractionNode] = []
    for expected_id, payload in enumerate(interaction_payloads):
        if not isinstance(payload, dict):
            raise NativeEvaluationError("interaction table entry is not a JSON object")
        interaction_id = int(payload["id"])
        if interaction_id != expected_id:
            raise NativeEvaluationError(
                "interaction table ids must be contiguous and ordered in process manifest"
            )
        left_id = int(payload["left_id"])
        right_id = int(payload["right_id"])
        result_id = int(payload["result_id"])
        coupling = payload.get("coupling")
        if not isinstance(coupling, list) or len(coupling) != 2:
            raise NativeEvaluationError(
                "interaction table entry has invalid coupling metadata"
            )
        interactions.append(
            SharedInteractionNode(
                id=interaction_id,
                vertex_kind=int(payload["vertex_kind"]),
                left_id=left_id,
                right_id=right_id,
                result_id=result_id,
                left=currents[left_id].key,
                right=currents[right_id].key,
                result=currents[result_id].key,
                coupling=(float(coupling[0]), float(coupling[1])),
            )
        )

    interactions_by_result_payload = table_payload.get("interactions_by_result")
    if not isinstance(interactions_by_result_payload, list):
        raise NativeEvaluationError(
            "process manifest table is missing interactions_by_result"
        )
    interactions_by_result = tuple(
        tuple(int(interaction_id) for interaction_id in payload)
        for payload in interactions_by_result_payload
        if isinstance(payload, list)
    )
    if len(interactions_by_result) != len(interactions_by_result_payload):
        raise NativeEvaluationError("interactions_by_result entry is not a list")
    if len(interactions_by_result) != len(currents):
        raise NativeEvaluationError(
            "interactions_by_result length does not match current table length"
        )

    amplitude_payloads = table_payload.get("amplitudes")
    if not isinstance(amplitude_payloads, list):
        raise NativeEvaluationError("process manifest table is missing amplitudes")
    amplitudes = tuple(
        SharedAmplitudeRecord(
            left_id=int(payload["left_id"]),
            right_id=int(payload["right_id"]),
            helicities=tuple(int(helicity) for helicity in payload["helicities"]),
            multiplicity=int(payload["multiplicity"]),
        )
        for payload in amplitude_payloads
        if isinstance(payload, dict)
    )
    if len(amplitudes) != len(amplitude_payloads):
        raise NativeEvaluationError("amplitude table entry is not a JSON object")

    return SharedCurrentTable(
        currents=tuple(currents),
        sources=sources,
        interactions=tuple(interactions),
        interactions_by_result=interactions_by_result,
        amplitudes=amplitudes,
    )


def _decimal_string(value: float) -> str:
    return format(float(value), ".17g")


_RUSTICOL_STANDALONE_CHECK_SCRIPT = r'''#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import json
import math
import statistics
import sys
import time
from decimal import Decimal
from pathlib import Path

try:
    import numpy as np
except ImportError:  # pragma: no cover - standalone convenience path
    np = None


class Color:
    enabled = sys.stdout.isatty()
    cyan = "\033[36m" if enabled else ""
    green = "\033[32m" if enabled else ""
    yellow = "\033[33m" if enabled else ""
    red = "\033[31m" if enabled else ""
    bold = "\033[1m" if enabled else ""
    dim = "\033[2m" if enabled else ""
    reset = "\033[0m" if enabled else ""


def style(text, color: str):
    return f"{getattr(Color, color)}{text}{Color.reset}" if Color.enabled else str(text)


def format_float(value, digits: int = 4):
    return f"{float(value):.{digits}g}"


def format_measurement(value, error=None, unit: str = ""):
    formatted = format_float(value)
    if error is not None:
        formatted = f"{formatted} +/- {format_float(error, 2)}"
    return f"{formatted} {unit}".strip()


def render_table(title: str, rows):
    columns = ("Metric", "Value")
    string_rows = [(str(metric), str(value), row_style) for metric, value, row_style in rows]
    width_metric = max([len(columns[0])] + [len(metric) for metric, _, _ in string_rows])
    width_value = max([len(columns[1])] + [len(value) for _, value, _ in string_rows])
    border = f"+-{'-' * width_metric}-+-{'-' * width_value}-+"
    lines = [style(title, "bold"), border]
    lines.append(
        "| "
        + style(columns[0].ljust(width_metric), "bold")
        + " | "
        + style(columns[1].rjust(width_value), "bold")
        + " |"
    )
    lines.append(border)
    for metric, value, row_style in string_rows:
        metric_cell = metric.ljust(width_metric)
        value_cell = value.rjust(width_value)
        if row_style:
            metric_cell = style(metric_cell, row_style)
            value_cell = style(value_cell, row_style)
        lines.append(f"| {metric_cell} | {value_cell} |")
    lines.append(border)
    return "\n".join(lines)


def print_table(title: str, rows):
    print(render_table(title, rows))


def rusticol_candidate_paths(root: Path, explicit_folder: Path | None = None):
    seeds = []
    if explicit_folder is not None:
        seeds.append(explicit_folder.expanduser())
    seeds.extend([root, Path.cwd(), *root.parents])

    seen = set()
    for seed in seeds:
        for candidate in expand_rusticol_seed(seed):
            try:
                resolved = candidate.resolve()
            except OSError:
                resolved = candidate
            key = str(resolved)
            if key in seen:
                continue
            seen.add(key)
            yield resolved


def expand_rusticol_seed(seed: Path):
    yield seed
    if seed.name == "rusticol":
        yield seed.parent
    for pattern in (
        "lib/python*/site-packages",
        ".venv/lib/python*/site-packages",
        "dependencies/.venv/lib/python*/site-packages",
        "pyAmpliCol/dependencies/.venv/lib/python*/site-packages",
    ):
        yield from seed.glob(pattern)


def import_rusticol(root: Path, rusticol_folder: Path | None):
    errors = []
    searched = []
    for candidate in rusticol_candidate_paths(root, rusticol_folder):
        if not candidate.exists():
            continue
        sys.path.insert(0, str(candidate))
        searched.append(candidate)
        loaded = False
        try:
            import rusticol

            loaded = True
            return rusticol
        except ModuleNotFoundError as exc:
            if exc.name != "rusticol":
                raise
            errors.append(str(exc))
        finally:
            if not loaded:
                with contextlib.suppress(ValueError):
                    sys.path.remove(str(candidate))

    print(
        "Could not import rusticol. Re-run with --rusticol-folder pointing to "
        "the pyAmpliCol managed venv, a site-packages directory, or the "
        "rusticol package directory.",
        file=sys.stderr,
    )
    if searched:
        print("Searched import roots:", file=sys.stderr)
        for candidate in searched:
            print(f"  - {candidate}", file=sys.stderr)
    else:
        print("No existing candidate import roots were found.", file=sys.stderr)
    if errors:
        print(f"Last import error: {errors[-1]}", file=sys.stderr)
    return None


def load_points(root: Path, precision: int):
    payload = json.loads((root / "validation_momenta.json").read_text())
    points = [
        [[Decimal(component) for component in particle["momentum"]] for particle in point]
        for point in payload["points"]
    ]
    if precision == 16 and np is not None:
        return np.asarray(points, dtype=np.float64)
    if precision == 16:
        return [[[float(component) for component in particle] for particle in point] for point in points]
    return points


def repeat_points(points, count: int):
    if np is not None and hasattr(points, "shape"):
        reps = int(math.ceil(count / max(int(points.shape[0]), 1)))
        return np.tile(points, (reps, 1, 1))[:count]
    reps = int(math.ceil(count / max(len(points), 1)))
    return (points * reps)[:count]


def evaluate(runtime, points, precision: int):
    if precision == 16:
        return runtime.evaluate(points)
    return runtime.evaluate_with_prec(points, precision)


def profile(runtime, points, precision: int, target_s: float):
    estimate_points = repeat_points(points, 16)
    t0 = time.perf_counter()
    evaluate(runtime, estimate_points, precision)
    estimate_s = max(time.perf_counter() - t0, 1.0e-9)
    per_point = estimate_s / 16.0
    sample_count = max(16, int(target_s / per_point))
    block_count = 8
    block_size = max(1, sample_count // block_count)
    samples = []
    breakdown = None
    for _ in range(block_count):
        batch = repeat_points(points, block_size)
        t0 = time.perf_counter()
        profile_payload = runtime.profile(batch, precision=precision)
        elapsed = time.perf_counter() - t0
        samples.append(elapsed / block_size)
        breakdown = profile_payload
    mean = statistics.fmean(samples)
    stderr = statistics.stdev(samples) / math.sqrt(len(samples)) if len(samples) > 1 else 0.0
    return {
        "samples": block_count * block_size,
        "block_count": block_count,
        "block_size": block_size,
        "mean_us_per_point": mean * 1.0e6,
        "stderr_us_per_point": stderr * 1.0e6,
        "last_profile": breakdown,
    }


def profile_breakdown_rows(payload):
    last_profile = payload.get("last_profile")
    if not isinstance(last_profile, dict):
        return []
    points = max(int(last_profile.get("points", payload.get("block_size", 1))), 1)
    rows = []
    labels = {
        "source_fill_time_s": "Source fill",
        "parameter_pack_time_s": "Parameter pack",
        "stage_evaluator_time_s": "Stage evaluators",
        "amplitude_evaluator_time_s": "Amplitude evaluator",
        "me_reduction_time_s": "ME reduction",
        "output_transfer_time_s": "Output transfer",
        "wall_time_s": "Profile wall",
    }
    for key, label in labels.items():
        value = last_profile.get(key)
        if not isinstance(value, (int, float)):
            continue
        row_style = "cyan" if key in {"stage_evaluator_time_s", "amplitude_evaluator_time_s"} else None
        rows.append((label, format_measurement(float(value) * 1.0e6 / points, unit="us/point"), row_style))
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Standalone rusticol process check")
    parser.add_argument("--precision", type=int, default=16)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--target-runtime", type=float, default=10.0)
    parser.add_argument(
        "--rusticol-folder",
        type=Path,
        help=(
            "Optional rusticol location. Accepts the pyAmpliCol managed venv, "
            "a site-packages directory, or the rusticol package directory."
        ),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    rusticol = import_rusticol(root, args.rusticol_folder)
    if rusticol is None:
        return 1

    manifest = json.loads((root / "process_manifest.json").read_text())
    runtime = rusticol.Runtime.load(str(root))
    points = load_points(root, args.precision)
    values = evaluate(runtime, points, args.precision)

    backend = manifest.get("backend_settings", {})
    helicity_filter = manifest.get("helicity_filter", {})
    print_table(
        "Rusticol Standalone Check",
        [
            ("Process", manifest["process"], "bold"),
            ("family", manifest.get("family"), None),
            ("gluons", manifest.get("gluon_count"), None),
            (
                "backend",
                backend.get("compiled_kind") if isinstance(backend, dict) else None,
                None,
            ),
            (
                "retained amplitudes",
                helicity_filter.get("retained_amplitude_count")
                if isinstance(helicity_filter, dict)
                else None,
                None,
            ),
            ("precision", args.precision, None),
            ("points", len(values), None),
            ("matrix element(s)", values, "green"),
        ],
    )

    if args.profile:
        payload = profile(runtime, points, args.precision, args.target_runtime)
        print_table(
            "Rusticol Timing Summary",
            [
                (
                    "Wall runtime",
                    format_measurement(
                        payload["mean_us_per_point"],
                        payload["stderr_us_per_point"],
                        unit="us/point",
                    ),
                    "green",
                ),
                (
                    "Samples",
                    f"{payload['samples']} ({payload['block_count']} x {payload['block_size']})",
                    None,
                ),
                ("Target runtime", format_measurement(args.target_runtime, unit="s"), None),
            ],
        )
        rows = profile_breakdown_rows(payload)
        if rows:
            print_table("Rusticol Runtime Breakdown", rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def _artifact_path_for_manifest(path: Path, artifact_dir: Path) -> str:
    artifact_root = artifact_dir.resolve()
    source_path = path.expanduser()
    source_resolved = source_path.resolve()
    try:
        source_resolved.relative_to(artifact_root)
    except ValueError:
        if source_path.exists():
            target_dir = artifact_dir / "compiled" / "repackaged"
            target_dir.mkdir(parents=True, exist_ok=True)
            target = target_dir / source_path.name
            if source_resolved != target.resolve():
                shutil.copy2(source_path, target)
            source_path = target
        else:
            source_path = source_resolved
    try:
        return os.path.relpath(source_path, artifact_dir)
    except ValueError:
        return str(source_path)


def _artifact_subdirectory(artifact_dir: Path, name: str) -> Path:
    directory = artifact_dir / name
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def _artifact_path_from_manifest(path: str, artifact_dir: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return artifact_dir / candidate


def _contiguous_column_slice(columns: np.ndarray) -> slice | None:
    if columns.size == 0:
        return slice(0, 0)
    start = int(columns[0])
    stop = int(columns[-1]) + 1
    if stop - start != int(columns.size):
        return None
    if not np.array_equal(columns, np.arange(start, stop, dtype=columns.dtype)):
        return None
    return slice(start, stop)


def _build_shared_global_parameter_layout(
    table: SharedCurrentTable,
) -> _SharedGlobalParameterLayout:
    builder = ParamBuilder()
    current_offsets: list[int] = []
    for current in table.currents:
        head = ("shared_global", "current", str(current.id))
        builder.add_parameter_list(
            head,
            6,
            role="shared_global_current",
        )
        current_offsets.append(builder.positions[head][0])

    momentum_offsets = _register_shared_momentum_parameters(builder, table)

    return _SharedGlobalParameterLayout(
        parameter_symbols=tuple(builder.parameter_symbols()),
        parameter_count=len(builder.values),
        current_offsets=tuple(current_offsets),
        momentum_offsets=momentum_offsets,
        real_valued_inputs=tuple(builder.real_valued_inputs),
    )


def _build_shared_split_global_parameter_layout(
    table: SharedCurrentTable,
) -> _SharedGlobalParameterLayout:
    builder = ParamBuilder()
    current_offsets: list[int] = []
    for current in table.currents:
        head = ("shared_global", "current", str(current.id))
        builder.add_parameter_list(
            head,
            6,
            role="shared_global_current",
        )
        current_offsets.append(builder.positions[head][0])

    interaction_offsets: list[int] = []
    for interaction in table.interactions:
        head = ("shared_global", "interaction", str(interaction.id))
        builder.add_parameter_list(
            head,
            6,
            role="shared_global_interaction",
        )
        interaction_offsets.append(builder.positions[head][0])

    momentum_offsets = _register_shared_momentum_parameters(builder, table)

    return _SharedGlobalParameterLayout(
        parameter_symbols=tuple(builder.parameter_symbols()),
        parameter_count=len(builder.values),
        current_offsets=tuple(current_offsets),
        momentum_offsets=momentum_offsets,
        real_valued_inputs=tuple(builder.real_valued_inputs),
        interaction_offsets=tuple(interaction_offsets),
    )


def _register_shared_momentum_parameters(
    builder: ParamBuilder,
    table: SharedCurrentTable,
) -> dict[int, int]:
    momentum_offsets: dict[int, int] = {}
    offsets_by_labels: dict[tuple[int, ...], int] = {}
    for current_id in _shared_momentum_current_ids(table):
        labels = _momentum_labels_for_current(table.currents[current_id])
        offset = offsets_by_labels.get(labels)
        if offset is None:
            head = (
                "shared_global",
                "momentum",
                "_".join(str(label) for label in labels),
            )
            builder.add_parameter_list(
                head,
                4,
                role="shared_global_current_momentum",
                real_valued=True,
            )
            offset = builder.positions[head][0]
            offsets_by_labels[labels] = offset
        momentum_offsets[current_id] = offset
    return momentum_offsets


def _unique_momentum_offsets_and_labels(
    table: SharedCurrentTable,
    layout: _SharedGlobalParameterLayout,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    labels_by_offset: dict[int, tuple[int, ...]] = {}
    for current_id, offset in layout.momentum_offsets.items():
        labels_by_offset.setdefault(
            offset,
            _momentum_labels_for_current(table.currents[current_id]),
        )
    return tuple(sorted(labels_by_offset.items()))


def _momentum_labels_for_current(current: SharedCurrentNode) -> tuple[int, ...]:
    return tuple(sorted(int(label) for label in current.external_labels))


def _shared_momentum_current_ids(table: SharedCurrentTable) -> tuple[int, ...]:
    current_ids: set[int] = set()
    for current in table.currents:
        if current.needs_propagator:
            current_ids.add(current.id)
    for interaction in table.interactions:
        if interaction.vertex_kind == 0:
            current_ids.add(interaction.left_id)
            current_ids.add(interaction.right_id)
    return tuple(sorted(current_ids))


def _shared_stage_interaction_ids(
    table: SharedCurrentTable,
    current_ids: Sequence[int],
) -> tuple[int, ...]:
    return tuple(
        interaction_id
        for current_id in current_ids
        for interaction_id in table.interactions_by_result[current_id]
    )


def _interaction_dimension(
    table: SharedCurrentTable,
    interaction_id: int,
) -> int:
    return table.currents[table.interactions[interaction_id].result_id].dimension


class _CompiledCurrentStage:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        current_ids: tuple[int, ...],
        stage_index: int,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> None:
        self.layout = layout
        self.current_ids = current_ids
        self.input_current_ids: tuple[int, ...] = ()
        self.momentum_current_ids: tuple[int, ...] = ()
        self.output_slices: dict[int, tuple[int, int]] = {}
        self._register_parameters(table)
        outputs = self._build_outputs(table)
        self.output_length = len(outputs)
        self.output_slots = tuple(self.output_slices.items())
        (
            self.output_columns,
            self.output_current_ids,
            self.output_current_components,
        ) = _output_slot_index_arrays(self.output_slots)
        self.evaluator = _compile_symbolica_outputs(
            outputs,
            list(layout.parameter_symbols),
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            real_params=layout.real_valued_inputs,
            symbolica_settings=symbolica_settings,
            label=f"shared_stage_{stage_index}",
        )
        self.parameter_count = layout.parameter_count

    @classmethod
    def from_artifact(
        cls,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        manifest: dict[str, Any],
        artifact_dir: Path,
    ) -> "_CompiledCurrentStage":
        instance = cls.__new__(cls)
        instance.layout = layout
        instance.current_ids = tuple(int(value) for value in manifest["current_ids"])
        instance.input_current_ids = tuple(
            int(value) for value in manifest["input_current_ids"]
        )
        instance.momentum_current_ids = tuple(
            int(value) for value in manifest["momentum_current_ids"]
        )
        instance.output_slices = _output_slices_from_manifest(manifest)
        instance.output_length = int(manifest["output_length"])
        instance.output_slots = tuple(instance.output_slices.items())
        (
            instance.output_columns,
            instance.output_current_ids,
            instance.output_current_components,
        ) = _output_slot_index_arrays(instance.output_slots)
        instance.evaluator = _load_symbolica_evaluator_artifact(
            manifest["evaluator"],
            artifact_dir,
        )
        instance.parameter_count = layout.parameter_count
        return instance

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        return {
            "kind": "current-stage",
            "current_ids": list(self.current_ids),
            "input_current_ids": list(self.input_current_ids),
            "momentum_current_ids": list(self.momentum_current_ids),
            "output_length": self.output_length,
            "output_slots": _output_slots_to_manifest(self.output_slots),
            "evaluator": _symbolica_evaluator_artifact_manifest(
                self.evaluator,
                artifact_dir,
            ),
        }

    def evaluate_and_assign(
        self,
        state_rows: np.ndarray,
        current_rows: np.ndarray,
        table: SharedCurrentTable,
    ) -> DAGEvaluationTiming:
        pack_start = time.perf_counter()
        parameter_pack_time_s = time.perf_counter() - pack_start
        start = time.perf_counter()
        evaluated = _evaluate_complex_outputs(self.evaluator, state_rows)
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        if self.output_columns.size:
            _assign_complex_outputs(
                current_rows,
                self.output_current_ids,
                self.output_current_components,
                self.output_columns,
                evaluated,
            )
        output_transfer_time_s = time.perf_counter() - transfer_start
        return DAGEvaluationTiming(
            parameter_pack_time_s=parameter_pack_time_s,
            evaluator_time_s=evaluator_time_s,
            output_transfer_time_s=output_transfer_time_s,
        )

    def _register_parameters(
        self,
        table: SharedCurrentTable,
    ) -> None:
        input_current_ids = set()
        momentum_current_ids = set()
        for current_id in self.current_ids:
            current = table.currents[current_id]
            if current.needs_propagator:
                momentum_current_ids.add(current_id)
            for interaction_id in table.interactions_by_result[current_id]:
                interaction = table.interactions[interaction_id]
                input_current_ids.add(interaction.left_id)
                input_current_ids.add(interaction.right_id)
                if interaction.vertex_kind == 0:
                    momentum_current_ids.add(interaction.left_id)
                    momentum_current_ids.add(interaction.right_id)
        self.input_current_ids = tuple(sorted(input_current_ids))
        self.momentum_current_ids = tuple(sorted(momentum_current_ids))

    def _build_outputs(self, table: SharedCurrentTable) -> tuple[Any, ...]:
        parameter_symbols = self.layout.parameter_symbols
        input_expressions: dict[int, tuple[Any, ...]] = {}
        for current_id in self.input_current_ids:
            current = table.currents[current_id]
            start = self.layout.current_offsets[current_id]
            stop = start + current.dimension
            input_expressions[current_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )

        momentum_expressions: dict[int, tuple[Any, ...]] = {}
        for current_id in self.momentum_current_ids:
            start = self.layout.momentum_offsets[current_id]
            stop = start + 4
            momentum_expressions[current_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )

        outputs: list[Any] = []
        for current_id in self.current_ids:
            current = table.currents[current_id]
            total = tuple(0j for _ in range(current.dimension))
            for interaction_id in table.interactions_by_result[current_id]:
                interaction = table.interactions[interaction_id]
                contribution = _expr_vertex(
                    interaction,
                    input_expressions[interaction.left_id],
                    input_expressions[interaction.right_id],
                    momentum_expressions,
                )
                total = _sum_components(total, contribution)
            if current.needs_propagator:
                total = _expr_propagate(
                    current.key,
                    total,
                    momentum_expressions[current.id],
                )
            start = len(outputs)
            outputs.extend(total)
            self.output_slices[current_id] = (start, len(outputs))
        return tuple(outputs)


class _CompiledSplitCurrentStage:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        current_ids: tuple[int, ...],
        stage_index: int,
        stage_count: int,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
        progress_callback: ProgressCallback | None = None,
    ) -> None:
        interaction_ids = _shared_stage_interaction_ids(table, current_ids)
        _report_progress(
            progress_callback,
            stage="compile",
            item=f"vertices {stage_index}/{stage_count}",
        )
        self.vertex_stage = _CompiledVertexStage(
            table,
            layout=layout,
            interaction_ids=interaction_ids,
            stage_index=stage_index,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            symbolica_settings=symbolica_settings,
        )
        _report_progress(
            progress_callback,
            stage="compile",
            item=f"vertices {stage_index}/{stage_count}",
            increment=1,
        )
        _report_progress(
            progress_callback,
            stage="compile",
            item=f"currents {stage_index}/{stage_count}",
        )
        self.combine_stage = _CompiledCurrentCombineStage(
            table,
            layout=layout,
            current_ids=current_ids,
            stage_index=stage_index,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            symbolica_settings=symbolica_settings,
        )
        _report_progress(
            progress_callback,
            stage="compile",
            item=f"currents {stage_index}/{stage_count}",
            increment=1,
        )
        self.current_ids = current_ids
        self.output_length = self.vertex_stage.output_length + self.combine_stage.output_length

    def evaluate_and_assign(
        self,
        state_rows: np.ndarray,
        current_rows: np.ndarray,
        interaction_rows: np.ndarray,
        table: SharedCurrentTable,
    ) -> DAGEvaluationTiming:
        vertex_timing = self.vertex_stage.evaluate_and_assign(
            state_rows,
            interaction_rows,
            table,
        )
        combine_timing = self.combine_stage.evaluate_and_assign(
            state_rows,
            current_rows,
            table,
        )
        return vertex_timing + combine_timing


class _CompiledVertexStage:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        interaction_ids: tuple[int, ...],
        stage_index: int,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> None:
        self.layout = layout
        self.interaction_ids = interaction_ids
        self.input_current_ids: tuple[int, ...] = ()
        self.momentum_current_ids: tuple[int, ...] = ()
        self.output_slices: dict[int, tuple[int, int]] = {}
        self._register_parameters(table)
        outputs = self._build_outputs(table)
        self.output_length = len(outputs)
        self.output_slots = tuple(self.output_slices.items())
        (
            self.output_columns,
            self.output_interaction_ids,
            self.output_interaction_components,
        ) = _output_slot_index_arrays(self.output_slots)
        self.evaluator = _compile_symbolica_outputs(
            outputs,
            list(layout.parameter_symbols),
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            real_params=layout.real_valued_inputs,
            symbolica_settings=symbolica_settings,
            label=f"shared_stage_{stage_index}_vertices",
        )
        self.parameter_count = layout.parameter_count

    def evaluate_and_assign(
        self,
        state_rows: np.ndarray,
        interaction_rows: np.ndarray,
        table: SharedCurrentTable,
    ) -> DAGEvaluationTiming:
        pack_start = time.perf_counter()
        parameter_pack_time_s = time.perf_counter() - pack_start
        start = time.perf_counter()
        evaluated = _evaluate_complex_outputs(self.evaluator, state_rows)
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        if self.output_columns.size:
            _assign_complex_outputs(
                interaction_rows,
                self.output_interaction_ids,
                self.output_interaction_components,
                self.output_columns,
                evaluated,
            )
        output_transfer_time_s = time.perf_counter() - transfer_start
        return DAGEvaluationTiming(
            parameter_pack_time_s=parameter_pack_time_s,
            evaluator_time_s=evaluator_time_s,
            output_transfer_time_s=output_transfer_time_s,
        )

    def _register_parameters(
        self,
        table: SharedCurrentTable,
    ) -> None:
        input_current_ids = set()
        momentum_current_ids = set()
        for interaction_id in self.interaction_ids:
            interaction = table.interactions[interaction_id]
            input_current_ids.add(interaction.left_id)
            input_current_ids.add(interaction.right_id)
            if interaction.vertex_kind == 0:
                momentum_current_ids.add(interaction.left_id)
                momentum_current_ids.add(interaction.right_id)
        self.input_current_ids = tuple(sorted(input_current_ids))
        self.momentum_current_ids = tuple(sorted(momentum_current_ids))

    def _build_outputs(self, table: SharedCurrentTable) -> tuple[Any, ...]:
        parameter_symbols = self.layout.parameter_symbols
        input_expressions: dict[int, tuple[Any, ...]] = {}
        for current_id in self.input_current_ids:
            current = table.currents[current_id]
            start = self.layout.current_offsets[current_id]
            stop = start + current.dimension
            input_expressions[current_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )

        momentum_expressions: dict[int, tuple[Any, ...]] = {}
        for current_id in self.momentum_current_ids:
            start = self.layout.momentum_offsets[current_id]
            stop = start + 4
            momentum_expressions[current_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )

        outputs: list[Any] = []
        for interaction_id in self.interaction_ids:
            interaction = table.interactions[interaction_id]
            contribution = _expr_vertex(
                interaction,
                input_expressions[interaction.left_id],
                input_expressions[interaction.right_id],
                momentum_expressions,
            )
            start = len(outputs)
            outputs.extend(contribution)
            self.output_slices[interaction_id] = (start, len(outputs))
        return tuple(outputs)


class _CompiledCurrentCombineStage:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        current_ids: tuple[int, ...],
        stage_index: int,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> None:
        self.layout = layout
        self.current_ids = current_ids
        self.input_interaction_ids: tuple[int, ...] = ()
        self.momentum_current_ids: tuple[int, ...] = ()
        self.output_slices: dict[int, tuple[int, int]] = {}
        self._register_parameters(table)
        outputs = self._build_outputs(table)
        self.output_length = len(outputs)
        self.output_slots = tuple(self.output_slices.items())
        (
            self.output_columns,
            self.output_current_ids,
            self.output_current_components,
        ) = _output_slot_index_arrays(self.output_slots)
        self.evaluator = _compile_symbolica_outputs(
            outputs,
            list(layout.parameter_symbols),
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            real_params=layout.real_valued_inputs,
            symbolica_settings=symbolica_settings,
            label=f"shared_stage_{stage_index}_currents",
        )
        self.parameter_count = layout.parameter_count

    def evaluate_and_assign(
        self,
        state_rows: np.ndarray,
        current_rows: np.ndarray,
        table: SharedCurrentTable,
    ) -> DAGEvaluationTiming:
        pack_start = time.perf_counter()
        parameter_pack_time_s = time.perf_counter() - pack_start
        start = time.perf_counter()
        evaluated = _evaluate_complex_outputs(self.evaluator, state_rows)
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        if self.output_columns.size:
            _assign_complex_outputs(
                current_rows,
                self.output_current_ids,
                self.output_current_components,
                self.output_columns,
                evaluated,
            )
        output_transfer_time_s = time.perf_counter() - transfer_start
        return DAGEvaluationTiming(
            parameter_pack_time_s=parameter_pack_time_s,
            evaluator_time_s=evaluator_time_s,
            output_transfer_time_s=output_transfer_time_s,
        )

    def _register_parameters(
        self,
        table: SharedCurrentTable,
    ) -> None:
        interaction_ids = set()
        momentum_current_ids = set()
        for current_id in self.current_ids:
            current = table.currents[current_id]
            if current.needs_propagator:
                momentum_current_ids.add(current_id)
            for interaction_id in table.interactions_by_result[current_id]:
                interaction_ids.add(interaction_id)
        self.input_interaction_ids = tuple(sorted(interaction_ids))
        self.momentum_current_ids = tuple(sorted(momentum_current_ids))

    def _build_outputs(self, table: SharedCurrentTable) -> tuple[Any, ...]:
        parameter_symbols = self.layout.parameter_symbols
        interaction_expressions: dict[int, tuple[Any, ...]] = {}
        for interaction_id in self.input_interaction_ids:
            start = self.layout.interaction_offsets[interaction_id]
            stop = start + _interaction_dimension(table, interaction_id)
            interaction_expressions[interaction_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )

        momentum_expressions: dict[int, tuple[Any, ...]] = {}
        for current_id in self.momentum_current_ids:
            start = self.layout.momentum_offsets[current_id]
            stop = start + 4
            momentum_expressions[current_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )

        outputs: list[Any] = []
        for current_id in self.current_ids:
            current = table.currents[current_id]
            total = tuple(0j for _ in range(current.dimension))
            for interaction_id in table.interactions_by_result[current_id]:
                total = _sum_components(total, interaction_expressions[interaction_id])
            if current.needs_propagator:
                total = _expr_propagate(
                    current.key,
                    total,
                    momentum_expressions[current.id],
                )
            start = len(outputs)
            outputs.extend(total)
            self.output_slices[current_id] = (start, len(outputs))
        return tuple(outputs)


class _CompiledAmplitudeStage:
    def __init__(
        self,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        merge_evaluators_strategy: bool,
        verbose_evaluator_build: bool,
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> None:
        self.layout = layout
        self.input_current_ids: tuple[int, ...] = ()
        self._register_parameters(table)
        self.outputs = self._build_outputs(table)
        self.output_length = len(self.outputs)
        self.raw_sum_weights = np.asarray(
            [amplitude.multiplicity for amplitude in table.amplitudes],
            dtype=np.float64,
        )
        self.merge_evaluators_strategy = merge_evaluators_strategy
        self.verbose_evaluator_build = verbose_evaluator_build
        self.symbolica_settings = symbolica_settings
        self._amplitude_evaluator: Any | None = None
        self._raw_sum_evaluator = self._build_raw_sum_evaluator()
        self.parameter_count = layout.parameter_count

    @classmethod
    def from_artifact(
        cls,
        table: SharedCurrentTable,
        *,
        layout: _SharedGlobalParameterLayout,
        manifest: dict[str, Any],
        artifact_dir: Path,
        symbolica_settings: SymbolicaEvaluatorSettings,
    ) -> "_CompiledAmplitudeStage":
        instance = cls.__new__(cls)
        instance.layout = layout
        instance.input_current_ids = tuple(
            int(value) for value in manifest["input_current_ids"]
        )
        instance.outputs = ()
        instance.output_length = int(manifest["output_length"])
        instance.raw_sum_weights = np.asarray(
            manifest["raw_sum_weights"],
            dtype=np.float64,
        )
        instance.merge_evaluators_strategy = False
        instance.verbose_evaluator_build = False
        instance.symbolica_settings = symbolica_settings
        amplitude_manifest = manifest.get("amplitude_evaluator")
        instance._amplitude_evaluator = (
            None
            if amplitude_manifest is None
            else _load_symbolica_evaluator_artifact(amplitude_manifest, artifact_dir)
        )
        raw_sum_manifest = manifest.get("raw_sum_evaluator")
        instance._raw_sum_evaluator = (
            None
            if raw_sum_manifest is None
            else _load_symbolica_evaluator_artifact(raw_sum_manifest, artifact_dir)
        )
        if instance._amplitude_evaluator is None and instance._raw_sum_evaluator is None:
            raise NativeEvaluationError(
                "compiled evaluator artifact amplitude stage contains no evaluator"
            )
        instance.parameter_count = layout.parameter_count
        if len(table.amplitudes) != instance.output_length:
            raise NativeEvaluationError(
                "compiled evaluator artifact amplitude count does not match process"
            )
        return instance

    def materialize(self) -> None:
        self._ensure_amplitude_evaluator()

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        return {
            "kind": "amplitude-stage",
            "input_current_ids": list(self.input_current_ids),
            "output_length": self.output_length,
            "raw_sum_weights": self.raw_sum_weights.tolist(),
            "amplitude_evaluator": _symbolica_evaluator_artifact_manifest(
                self._ensure_amplitude_evaluator(),
                artifact_dir,
            ),
            "raw_sum_evaluator": (
                None
                if self._raw_sum_evaluator is None
                else _symbolica_evaluator_artifact_manifest(
                    self._raw_sum_evaluator,
                    artifact_dir,
                )
            ),
        }

    def evaluate_rows(
        self,
        state_rows: np.ndarray,
        current_rows: np.ndarray,
    ) -> tuple[tuple[tuple[complex, ...], ...], DAGEvaluationTiming]:
        pack_start = time.perf_counter()
        parameter_pack_time_s = time.perf_counter() - pack_start
        start = time.perf_counter()
        evaluated = _evaluate_complex_outputs(
            self._ensure_amplitude_evaluator(),
            state_rows,
        )
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        evaluated_array = _concatenate_complex_outputs(evaluated)
        rows = tuple(
            tuple(row[: self.output_length].tolist())
            for row in evaluated_array
        )
        output_transfer_time_s = time.perf_counter() - transfer_start
        return (
            rows,
            DAGEvaluationTiming(
                parameter_pack_time_s=parameter_pack_time_s,
                evaluator_time_s=evaluator_time_s,
                output_transfer_time_s=output_transfer_time_s,
            ),
        )

    def evaluate_raw_sums(
        self,
        state_rows: np.ndarray,
        current_rows: np.ndarray,
    ) -> tuple[tuple[float, ...], DAGEvaluationTiming]:
        pack_start = time.perf_counter()
        parameter_pack_time_s = time.perf_counter() - pack_start
        start = time.perf_counter()
        if self._raw_sum_evaluator is None:
            evaluated = _evaluate_complex_outputs(
                self._ensure_amplitude_evaluator(),
                state_rows,
            )
        else:
            evaluated = _evaluate_complex_outputs(
                self._raw_sum_evaluator,
                state_rows,
            )
        evaluator_time_s = time.perf_counter() - start
        transfer_start = time.perf_counter()
        if self._raw_sum_evaluator is None:
            raw_sums = _weighted_abs2_sums(
                evaluated,
                self.raw_sum_weights,
                self.output_length,
            )
        else:
            evaluated_array = _concatenate_complex_outputs(evaluated)
            raw_sums = tuple(float(value.real) for value in evaluated_array[:, 0].tolist())
        output_transfer_time_s = time.perf_counter() - transfer_start
        return (
            raw_sums,
            DAGEvaluationTiming(
                parameter_pack_time_s=parameter_pack_time_s,
                evaluator_time_s=evaluator_time_s,
                output_transfer_time_s=output_transfer_time_s,
            ),
        )

    def _ensure_amplitude_evaluator(self) -> Any:
        if self._amplitude_evaluator is None:
            if not self.outputs:
                raise NativeEvaluationError(
                    "cannot build amplitude evaluator from a loaded runtime-only artifact"
                )
            self._amplitude_evaluator = _compile_symbolica_outputs(
                self.outputs,
                list(self.layout.parameter_symbols),
                merge_evaluators_strategy=self.merge_evaluators_strategy,
                verbose_evaluator_build=self.verbose_evaluator_build,
                real_params=self.layout.real_valued_inputs,
                symbolica_settings=self.symbolica_settings,
                label="shared_amplitude",
            )
        return self._amplitude_evaluator

    def _build_raw_sum_evaluator(self) -> Any | None:
        if (
            not self.symbolica_settings.raw_sum_final_stage
            or self.symbolica_settings.backend == "jit"
        ):
            return None
        raw_sum = sum(
            weight * amplitude * amplitude.conj()
            for weight, amplitude in zip(
                self.raw_sum_weights.tolist(),
                self.outputs,
                strict=True,
            )
        )
        return _compile_symbolica_outputs(
            (raw_sum,),
            list(self.layout.parameter_symbols),
            merge_evaluators_strategy=self.merge_evaluators_strategy,
            verbose_evaluator_build=self.verbose_evaluator_build,
            real_params=self.layout.real_valued_inputs,
            symbolica_settings=self.symbolica_settings,
            jit_compile=False,
            label="shared_raw_sum",
        )

    def _register_parameters(self, table: SharedCurrentTable) -> None:
        self.input_current_ids = tuple(
            sorted(
                {
                    current_id
                    for amplitude in table.amplitudes
                    for current_id in (amplitude.left_id, amplitude.right_id)
                }
            )
        )

    def _build_outputs(self, table: SharedCurrentTable) -> tuple[Any, ...]:
        parameter_symbols = self.layout.parameter_symbols
        input_expressions: dict[int, tuple[Any, ...]] = {}
        for current_id in self.input_current_ids:
            current = table.currents[current_id]
            start = self.layout.current_offsets[current_id]
            stop = start + current.dimension
            input_expressions[current_id] = tuple(
                parameter_symbols[index] for index in range(start, stop)
            )
        return tuple(
            _expr_dot_weyl(
                input_expressions[amplitude.left_id],
                input_expressions[amplitude.right_id],
            )
            for amplitude in table.amplitudes
        )


def _shared_current_stage_ids(
    table: SharedCurrentTable,
) -> tuple[tuple[int, ...], ...]:
    grouped: dict[int, list[int]] = {}
    for current in table.currents:
        if current.is_source:
            continue
        grouped.setdefault(len(current.external_labels), []).append(current.id)
    return tuple(tuple(grouped[size]) for size in sorted(grouped))


def _compile_symbolica_outputs(
    outputs: tuple[Any, ...],
    params: list[Any],
    *,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    aliases: Sequence[tuple[Any, Any]] = (),
    real_params: Sequence[int] = (),
    symbolica_settings: SymbolicaEvaluatorSettings | None = None,
    jit_compile: bool = True,
    label: str = "symbolica",
) -> Any:
    if not outputs:
        raise NativeEvaluationError("cannot build evaluator with zero outputs")
    settings = symbolica_settings or SymbolicaEvaluatorSettings()
    outputs = tuple(_prepare_symbolica_output(output, settings) for output in outputs)
    chunk_size = settings.compiled_output_chunk_size
    if chunk_size is not None and len(outputs) > chunk_size:
        chunk_output_dir = (
            None
            if settings.compiled_output_dir is None
            else str(Path(settings.compiled_output_dir) / _safe_symbol_name(label))
        )
        unchunked_settings = replace(settings, compiled_output_chunk_size=None)
        if chunk_output_dir is not None:
            unchunked_settings = replace(
                unchunked_settings,
                compiled_output_dir=chunk_output_dir,
            )
        chunk_ranges = tuple(enumerate(range(0, len(outputs), chunk_size)))

        def compile_chunk(chunk_index: int, start: int) -> Any:
            stop = min(start + chunk_size, len(outputs))
            return _compile_symbolica_outputs(
                outputs[start:stop],
                params,
                merge_evaluators_strategy=merge_evaluators_strategy,
                verbose_evaluator_build=verbose_evaluator_build,
                aliases=aliases,
                real_params=real_params,
                symbolica_settings=unchunked_settings,
                jit_compile=jit_compile,
                label=f"{label}_chunk_{chunk_index}",
            )

        workers = min(settings.compiled_chunk_compile_workers, len(chunk_ranges))
        if workers <= 1:
            chunks = [
                compile_chunk(chunk_index, start)
                for chunk_index, start in chunk_ranges
            ]
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [
                    executor.submit(compile_chunk, chunk_index, start)
                    for chunk_index, start in chunk_ranges
                ]
                chunks = [future.result() for future in futures]
        return _ChunkedSymbolicaEvaluator(tuple(chunks))
    evaluator_kwargs = _symbolica_evaluator_kwargs(
        settings,
        verbose=verbose_evaluator_build,
        jit_compile=jit_compile,
    )
    alias_kwargs = {"aliases": list(aliases)} if aliases else {}
    if merge_evaluators_strategy:
        evaluator = outputs[0].evaluator(
            params,
            **alias_kwargs,
            **evaluator_kwargs,
        )
        for expression in _progress_outputs(
            outputs[1:],
            enabled=verbose_evaluator_build,
        ):
            other = expression.evaluator(
                params,
                **alias_kwargs,
                **evaluator_kwargs,
            )
            evaluator.merge(
                other,
                cpe_iterations=(
                    1
                    if settings.cpe_iterations is None
                    else settings.cpe_iterations
                ),
            )
        if real_params:
            evaluator.set_real_params(
                list(real_params),
                sqrt_real=settings.real_param_sqrt_real,
                log_real=settings.real_param_log_real,
                powf_real=settings.real_param_powf_real,
                real_if_args_real=settings.real_param_real_if_args_real,
                verbose=verbose_evaluator_build,
            )
        return _finalize_symbolica_evaluator(
            evaluator,
            settings,
            label,
            input_len=len(params),
            output_len=len(outputs),
        )

    from symbolica import Expression

    evaluator = Expression.evaluator_multiple(
        outputs,
        params,
        **alias_kwargs,
        **evaluator_kwargs,
    )
    if real_params:
        evaluator.set_real_params(
            list(real_params),
            sqrt_real=settings.real_param_sqrt_real,
            log_real=settings.real_param_log_real,
            powf_real=settings.real_param_powf_real,
            real_if_args_real=settings.real_param_real_if_args_real,
            verbose=verbose_evaluator_build,
        )
    return _finalize_symbolica_evaluator(
        evaluator,
        settings,
        label,
        input_len=len(params),
        output_len=len(outputs),
    )


def _prepare_symbolica_output(
    output: Any,
    settings: SymbolicaEvaluatorSettings,
) -> Any:
    if not settings.collect_factors:
        return output
    collect_factors = getattr(output, "collect_factors", None)
    if callable(collect_factors):
        return collect_factors()
    return output


def _symbolica_evaluator_kwargs(
    settings: SymbolicaEvaluatorSettings,
    *,
    verbose: bool,
    jit_compile: bool = True,
) -> dict[str, Any]:
    return {
        "iterations": settings.iterations,
        "cpe_iterations": settings.cpe_iterations,
        "n_cores": settings.n_cores,
        "verbose": verbose,
        "jit_compile": jit_compile,
        "direct_translation": settings.direct_translation,
        "jit_direct_translation": settings.jit_direct_translation,
        "jit_optimization_level": settings.jit_optimization_level,
        "max_horner_scheme_variables": settings.max_horner_scheme_variables,
        "max_common_pair_cache_entries": settings.max_common_pair_cache_entries,
        "max_common_pair_distance": settings.max_common_pair_distance,
    }


def _resolve_compiled_preset(
    preset: str,
    *,
    gluon_count: int | None = None,
    inline_asm: str,
    optimization_level: int,
    output_chunk_size: int | None = None,
) -> tuple[str, int, int | None]:
    if preset == "manual":
        return inline_asm, optimization_level, output_chunk_size
    if preset == "adaptive":
        if gluon_count is None:
            raise NativeEvaluationError(
                "adaptive symbolica compiled preset requires a gluon count"
            )
        chunk_size = output_chunk_size
        if chunk_size is None:
            if 5 <= gluon_count <= 6:
                chunk_size = 64
            elif gluon_count == 7:
                chunk_size = 96
        if gluon_count <= 7:
            return "none", 1, chunk_size
        return "default", 3, chunk_size
    if preset == "generation":
        return "default", 3, output_chunk_size
    if preset == "balanced":
        return "none", 1, output_chunk_size
    if preset == "runtime":
        if gluon_count is None:
            return "none", 3, output_chunk_size
        chunk_size = output_chunk_size
        if chunk_size is None:
            if 5 <= gluon_count <= 6:
                chunk_size = 64
            elif gluon_count == 7:
                chunk_size = 96
            elif gluon_count == 8:
                chunk_size = 96
        if 5 <= gluon_count <= 6:
            return "none", 2, chunk_size
        if 7 <= gluon_count <= 8:
            return "none", 1, chunk_size
        return "none", 3, chunk_size
    if preset == "runtime-o3":
        if gluon_count is None:
            return "none", 3, output_chunk_size
        chunk_size = output_chunk_size
        if chunk_size is None:
            if 5 <= gluon_count <= 6:
                chunk_size = 64
            elif gluon_count >= 7:
                chunk_size = 96
        return "none", 3, chunk_size
    raise NativeEvaluationError(
        "symbolica compiled preset must be 'manual', 'adaptive', "
        "'generation', 'balanced', 'runtime', or 'runtime-o3'"
    )


def _finalize_symbolica_evaluator(
    evaluator: Any,
    settings: SymbolicaEvaluatorSettings,
    label: str,
    *,
    input_len: int,
    output_len: int,
) -> Any:
    if settings.backend == "jit":
        return _JITSymbolicaEvaluatorAdapter(
            evaluator,
            settings,
            label,
            input_len=input_len,
            output_len=output_len,
        )
    if settings.backend in ("compiled-complex", "compiled-complex-4x"):
        return _CompiledComplexEvaluatorAdapter(
            evaluator,
            settings,
            label,
            input_len=input_len,
            output_len=output_len,
        )
    raise NativeEvaluationError(
        f"unsupported symbolica evaluator backend: {settings.backend}"
    )


class _JITSymbolicaEvaluatorAdapter:
    def __init__(
        self,
        evaluator: Any,
        settings: SymbolicaEvaluatorSettings,
        label: str,
        *,
        input_len: int,
        output_len: int,
    ) -> None:
        self.input_len = int(input_len)
        self.output_len = int(output_len)
        self.backend = settings.backend
        self.label = _safe_symbol_name(label)
        self._source_evaluator = evaluator
        self.evaluator_state_path: Path | None = None

    def evaluate_complex(self, parameter_rows: Any) -> Any:
        return self._source_evaluator.evaluate_complex(parameter_rows)

    def evaluate(self, parameter_rows: Any) -> Any:
        return self._source_evaluator.evaluate(parameter_rows)

    def _evaluate_complex_prepared(self, parameter_rows: np.ndarray) -> Any:
        return self._source_evaluator.evaluate_complex(parameter_rows)

    def materialize(self) -> None:
        self._ensure_jit_compiled()

    def _ensure_jit_compiled(self) -> None:
        # The Symbolica saved evaluator state includes a JIT payload only after
        # the complex evaluator has been called at least once.
        dummy = np.ones((1, self.input_len), dtype=np.complex128)
        self._source_evaluator.evaluate_complex(dummy)

    @classmethod
    def from_artifact(
        cls,
        manifest: dict[str, Any],
        artifact_dir: Path,
    ) -> "_JITSymbolicaEvaluatorAdapter":
        from symbolica import Evaluator

        instance = cls.__new__(cls)
        instance.input_len = int(manifest["input_len"])
        instance.output_len = int(manifest["output_len"])
        instance.backend = str(manifest["backend"])
        instance.label = str(manifest.get("label", "jit_symbolica_evaluator"))
        instance.evaluator_state_path = _artifact_path_from_manifest(
            str(manifest["evaluator_state_path"]),
            artifact_dir,
        )
        instance._source_evaluator = Evaluator.load(
            instance.evaluator_state_path.read_bytes()
        )
        return instance

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        self._ensure_jit_compiled()
        evaluator_dir = _artifact_subdirectory(artifact_dir, "evaluators")
        path = evaluator_dir / f"{self.label}_{uuid.uuid4().hex}.evaluator.bin"
        path.write_bytes(self._source_evaluator.save())
        self.evaluator_state_path = path
        return {
            "kind": "jit-symbolica-evaluator",
            "backend": self.backend,
            "label": self.label,
            "input_len": self.input_len,
            "output_len": self.output_len,
            "evaluator_state_path": _artifact_path_for_manifest(path, artifact_dir),
        }


class _CompiledComplexEvaluatorAdapter:
    def __init__(
        self,
        evaluator: Any,
        settings: SymbolicaEvaluatorSettings,
        label: str,
        *,
        input_len: int,
        output_len: int,
    ) -> None:
        if settings.backend == "compiled-complex-4x" and _aarch64_platform():
            raise NativeEvaluationError(
                "Symbolica compiled-complex-4x is currently unsafe on aarch64: "
                "the Rust side uses Complex<wide::f64x4>, while xsimd's "
                "std::complex<double> batch is two lanes on this architecture. "
                "Use compiled-complex with --symbolica-compiled-inline-asm "
                "default or none instead."
            )
        safe_label = _safe_symbol_name(label)
        unique = uuid.uuid4().hex
        function_name = f"pyamplicol_{safe_label}_{unique}"
        self.function_name = function_name
        self.input_len = int(input_len)
        self.output_len = int(output_len)
        self.backend = settings.backend
        self.number_type = (
            "complex_4x"
            if settings.backend == "compiled-complex-4x"
            else "complex"
        )
        self._source_evaluator = evaluator
        if settings.compiled_output_dir is None:
            self._tmpdir: tempfile.TemporaryDirectory[str] | None = (
                tempfile.TemporaryDirectory(prefix="pyamplicol-symbolica-")
            )
            path = Path(self._tmpdir.name)
        else:
            self._tmpdir = None
            path = Path(settings.compiled_output_dir).expanduser()
            path.mkdir(parents=True, exist_ok=True)
        self.source_path = path / f"{function_name}.cpp"
        self.library_path = path / f"lib{function_name}"
        self.evaluator_state_path: Path | None = path / f"{function_name}.evaluator.bin"
        save = getattr(evaluator, "save", None)
        if callable(save):
            self.evaluator_state_path.write_bytes(save())
        else:
            self.evaluator_state_path = None
        self._compiled = evaluator.compile(
            function_name,
            str(self.source_path),
            str(self.library_path),
            self.number_type,
            inline_asm=settings.compiled_inline_asm,
            optimization_level=settings.compiled_optimization_level,
            native=settings.compiled_native,
            compiler_path=settings.compiler_path,
            compiler_flags=_compiled_compiler_flags(settings),
        )

    def evaluate_complex(self, parameter_rows: Any) -> Any:
        return self._evaluate_complex_prepared(_complex128_parameter_rows(parameter_rows))

    def _evaluate_complex_prepared(self, parameter_rows: np.ndarray) -> Any:
        return self._compiled.evaluate(parameter_rows)

    @classmethod
    def from_artifact(
        cls,
        manifest: dict[str, Any],
        artifact_dir: Path,
    ) -> "_CompiledComplexEvaluatorAdapter":
        from symbolica import CompiledComplexEvaluator, CompiledSimdComplexEvaluator

        instance = cls.__new__(cls)
        instance._tmpdir = None
        instance.function_name = str(manifest["function_name"])
        instance.input_len = int(manifest["input_len"])
        instance.output_len = int(manifest["output_len"])
        instance.backend = str(manifest["backend"])
        instance.number_type = str(manifest["number_type"])
        instance._source_evaluator = None
        instance.source_path = _artifact_path_from_manifest(
            str(manifest["source_path"]),
            artifact_dir,
        )
        instance.library_path = _artifact_path_from_manifest(
            str(manifest["library_path"]),
            artifact_dir,
        )
        state_path = manifest.get("evaluator_state_path")
        instance.evaluator_state_path = (
            None
            if state_path is None
            else _artifact_path_from_manifest(str(state_path), artifact_dir)
        )
        if instance.number_type == "complex_4x":
            loader = CompiledSimdComplexEvaluator
        elif instance.number_type == "complex":
            loader = CompiledComplexEvaluator
        else:
            raise NativeEvaluationError(
                f"unsupported compiled evaluator number type: {instance.number_type!r}"
            )
        instance._compiled = loader.load(
            str(instance.library_path),
            instance.function_name,
            instance.input_len,
            instance.output_len,
        )
        return instance

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        return {
            "kind": "compiled-complex-evaluator",
            "backend": self.backend,
            "number_type": self.number_type,
            "function_name": self.function_name,
            "input_len": self.input_len,
            "output_len": self.output_len,
            "source_path": _artifact_path_for_manifest(
                self.source_path,
                artifact_dir,
            ),
            "library_path": _artifact_path_for_manifest(
                self.library_path,
                artifact_dir,
            ),
            "evaluator_state_path": (
                None
                if self.evaluator_state_path is None
                else _artifact_path_for_manifest(
                    self.evaluator_state_path,
                    artifact_dir,
                )
            ),
        }


def _compiled_compiler_flags(settings: SymbolicaEvaluatorSettings) -> tuple[str, ...]:
    flags = list(settings.compiler_flags)
    if settings.backend == "compiled-complex-4x":
        xsimd_include = (
            Path(__file__).resolve().parents[2] / "dependencies" / "xsimd" / "include"
        )
        include_flag = f"-I{xsimd_include}"
        if xsimd_include.exists() and include_flag not in flags:
            flags.append(include_flag)
    return tuple(flags)


def _aarch64_platform() -> bool:
    machine = platform.machine().lower()
    return machine in {"arm64", "aarch64"}


class _ChunkedSymbolicaEvaluator:
    def __init__(self, evaluators: tuple[Any, ...]) -> None:
        if not evaluators:
            raise NativeEvaluationError("chunked evaluator needs at least one chunk")
        self._evaluators = evaluators

    def evaluate_complex(self, parameter_rows: Any) -> Any:
        return np.concatenate(self.evaluate_complex_chunks(parameter_rows), axis=1)

    def evaluate_complex_chunks(self, parameter_rows: Any) -> tuple[Any, ...]:
        prepared_rows = _complex128_parameter_rows(parameter_rows)
        return tuple(
            _evaluate_prepared_complex(evaluator, prepared_rows)
            for evaluator in self._evaluators
        )

    @classmethod
    def from_artifact(
        cls,
        manifest: dict[str, Any],
        artifact_dir: Path,
    ) -> "_ChunkedSymbolicaEvaluator":
        chunks = manifest.get("chunks")
        if not isinstance(chunks, list):
            raise NativeEvaluationError("chunked evaluator artifact is missing chunks")
        return cls(
            tuple(
                _load_symbolica_evaluator_artifact(chunk, artifact_dir)
                for chunk in chunks
            )
        )

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        return {
            "kind": "chunked-symbolica-evaluator",
            "chunks": [
                _symbolica_evaluator_artifact_manifest(evaluator, artifact_dir)
                for evaluator in self._evaluators
            ],
        }


ComplexOutput = np.ndarray | tuple[np.ndarray, ...]


def _complex128_parameter_rows(parameter_rows: Any) -> np.ndarray:
    if (
        isinstance(parameter_rows, np.ndarray)
        and parameter_rows.dtype == np.complex128
        and parameter_rows.flags.c_contiguous
    ):
        return parameter_rows
    return np.asarray(parameter_rows, dtype=np.complex128)


def _evaluate_prepared_complex(evaluator: Any, parameter_rows: np.ndarray) -> Any:
    evaluate_prepared = getattr(evaluator, "_evaluate_complex_prepared", None)
    if callable(evaluate_prepared):
        return evaluate_prepared(parameter_rows)
    return evaluator.evaluate_complex(parameter_rows)


def _evaluate_complex_outputs(evaluator: Any, parameter_rows: Any) -> ComplexOutput:
    evaluate_chunks = getattr(evaluator, "evaluate_complex_chunks", None)
    if callable(evaluate_chunks):
        return tuple(
            np.asarray(chunk, dtype=np.complex128)
            for chunk in evaluate_chunks(parameter_rows)
        )
    return np.asarray(evaluator.evaluate_complex(parameter_rows), dtype=np.complex128)


def _symbolica_evaluator_artifact_manifest(
    evaluator: Any,
    artifact_dir: Path,
) -> dict[str, Any]:
    manifest = getattr(evaluator, "artifact_manifest", None)
    if not callable(manifest):
        raise NativeEvaluationError(
            "saving evaluator artifacts is currently supported only for "
            "compiled-complex shared-current DAG evaluators"
        )
    return manifest(artifact_dir)


def _load_symbolica_evaluator_artifact(
    manifest: Any,
    artifact_dir: Path,
) -> Any:
    if not isinstance(manifest, dict):
        raise NativeEvaluationError("compiled evaluator artifact entry is invalid")
    kind = manifest.get("kind")
    if kind == "jit-symbolica-evaluator":
        return _JITSymbolicaEvaluatorAdapter.from_artifact(manifest, artifact_dir)
    if kind == "compiled-complex-evaluator":
        return _CompiledComplexEvaluatorAdapter.from_artifact(manifest, artifact_dir)
    if kind == "chunked-symbolica-evaluator":
        return _ChunkedSymbolicaEvaluator.from_artifact(manifest, artifact_dir)
    raise NativeEvaluationError(f"unsupported evaluator artifact kind: {kind!r}")


def _concatenate_complex_outputs(evaluated: ComplexOutput) -> np.ndarray:
    if isinstance(evaluated, tuple):
        return np.concatenate(evaluated, axis=1)
    return evaluated


def _assign_complex_outputs(
    target_rows: np.ndarray,
    output_ids: np.ndarray,
    output_components: np.ndarray,
    output_columns: np.ndarray,
    evaluated: ComplexOutput,
) -> None:
    if isinstance(evaluated, tuple):
        evaluated_array = np.concatenate(evaluated, axis=1)
        target_rows[:, output_ids, output_components] = evaluated_array[:, output_columns]
        return
    target_rows[:, output_ids, output_components] = evaluated[:, output_columns]


def _selected_current_components(
    current_rows: np.ndarray,
    output_ids: np.ndarray,
    output_components: np.ndarray,
) -> np.ndarray:
    if output_ids.size == 0:
        return np.empty((current_rows.shape[0], 0), dtype=np.complex128)
    return current_rows[:, output_ids, output_components]


def _complex_array_checksum(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.complex128)
    if array.size == 0:
        return {
            "sum_re": 0.0,
            "sum_im": 0.0,
            "sum_abs2": 0.0,
            "max_abs": 0.0,
        }
    return {
        "sum_re": float(np.sum(array.real)),
        "sum_im": float(np.sum(array.imag)),
        "sum_abs2": float(np.sum(array.real * array.real + array.imag * array.imag)),
        "max_abs": float(np.max(np.abs(array))),
    }


def _weighted_abs2_sums(
    evaluated: ComplexOutput,
    weights: np.ndarray,
    output_length: int,
) -> tuple[float, ...]:
    return tuple(
        float(value)
        for value in _weighted_abs2_sums_array(
            evaluated,
            weights,
            output_length,
        ).tolist()
    )


def _weighted_abs2_sums_array(
    evaluated: ComplexOutput,
    weights: np.ndarray,
    output_length: int,
) -> np.ndarray:
    if isinstance(evaluated, tuple):
        if not evaluated:
            return np.empty(0, dtype=np.float64)
        raw_sums = np.zeros(evaluated[0].shape[0], dtype=np.float64)
        offset = 0
        for chunk in evaluated:
            width = min(int(chunk.shape[1]), output_length - offset)
            if width <= 0:
                break
            amplitudes = chunk[:, :width]
            squared = amplitudes.real * amplitudes.real + amplitudes.imag * amplitudes.imag
            raw_sums += squared @ weights[offset : offset + width]
            offset += width
        return raw_sums
    amplitudes = evaluated[:, :output_length]
    squared = amplitudes.real * amplitudes.real + amplitudes.imag * amplitudes.imag
    return np.dot(squared, weights)


def _safe_symbol_name(value: str) -> str:
    result = "".join(character if character.isalnum() else "_" for character in value)
    if not result:
        return "eval"
    if result[0].isdigit():
        return f"eval_{result}"
    return result


def _progress_outputs(
    outputs: tuple[Any, ...],
    *,
    enabled: bool,
) -> Iterable[Any]:
    if not enabled:
        return outputs
    try:
        from colorama import Fore, Style  # type: ignore[import-untyped]
        from progressbar import Bar, ETA, Percentage, ProgressBar  # type: ignore[import-not-found]
    except ImportError:
        return outputs
    widgets = [
        Fore.CYAN,
        " merging evaluators ",
        Percentage(),
        " ",
        Bar(),
        " ",
        ETA(),
        Style.RESET_ALL,
    ]
    return ProgressBar(max_value=len(outputs), widgets=widgets)(outputs)


def _expr_vertex(
    interaction: SharedInteractionNode,
    left: tuple[Any, ...],
    right: tuple[Any, ...],
    momentum_expressions: dict[int, tuple[Any, ...]],
) -> tuple[Any, ...]:
    kind = int(interaction.vertex_kind)
    if kind == 0:
        return _expr_three_gluon(
            left,
            momentum_expressions[interaction.left_id],
            right,
            momentum_expressions[interaction.right_id],
        )
    if kind == 1:
        return _expr_two_gluon_to_tensor(left, right)
    if kind == 2:
        return _expr_tensor_gluon_to_gluon(left, right)
    if kind == 3:
        return _expr_gluon_tensor_to_gluon(left, right)
    if kind == 6:
        return _expr_quark_vector_weyl(
            left,
            right,
            int(interaction.result.chirality),
        )
    if kind == 10:
        chirality = int(interaction.result.chirality)
        coupling = _weyl_coupling_for_chirality(chirality, interaction.coupling)
        return tuple(
            coupling * component
            for component in _expr_quark_vector_weyl(left, right, chirality)
        )
    raise NativeEvaluationError(f"unsupported shared DAG vertex kind: {kind}")


def _expr_three_gluon(
    left: tuple[Any, ...],
    left_momentum: tuple[Any, ...],
    right: tuple[Any, ...],
    right_momentum: tuple[Any, ...],
) -> tuple[Any, ...]:
    tmp1 = _expr_minkowski_dot(left, right)
    tmp2 = _expr_minkowski_dot_momentum(left, right_momentum)
    tmp3 = _expr_minkowski_dot_momentum(right, left_momentum)
    prefactor = 1j / math.sqrt(2.0)
    return tuple(
        prefactor
        * (
            tmp1 * (left_momentum[index] - right_momentum[index])
            + 2.0 * (tmp2 * right[index] - tmp3 * left[index])
        )
        for index in range(4)
    )


def _expr_two_gluon_to_tensor(
    left: tuple[Any, ...],
    right: tuple[Any, ...],
) -> tuple[Any, ...]:
    return (
        left[0] * right[1] - left[1] * right[0],
        left[0] * right[2] - left[2] * right[0],
        left[0] * right[3] - left[3] * right[0],
        left[1] * right[2] - left[2] * right[1],
        left[1] * right[3] - left[3] * right[1],
        left[2] * right[3] - left[3] * right[2],
    )


def _expr_tensor_gluon_to_gluon(
    tensor: tuple[Any, ...],
    gluon: tuple[Any, ...],
) -> tuple[Any, ...]:
    prefactor = 0.5j
    return (
        (tensor[0] * gluon[1] + tensor[1] * gluon[2] + tensor[2] * gluon[3])
        * prefactor,
        (tensor[0] * gluon[0] + tensor[3] * gluon[2] + tensor[4] * gluon[3])
        * prefactor,
        (tensor[1] * gluon[0] - tensor[3] * gluon[1] + tensor[5] * gluon[3])
        * prefactor,
        (tensor[2] * gluon[0] - tensor[4] * gluon[1] - tensor[5] * gluon[2])
        * prefactor,
    )


def _expr_gluon_tensor_to_gluon(
    gluon: tuple[Any, ...],
    tensor: tuple[Any, ...],
) -> tuple[Any, ...]:
    prefactor = 0.5j
    return (
        (-gluon[1] * tensor[0] - gluon[2] * tensor[1] - gluon[3] * tensor[2])
        * prefactor,
        (-gluon[0] * tensor[0] - gluon[2] * tensor[3] - gluon[3] * tensor[4])
        * prefactor,
        (-gluon[0] * tensor[1] + gluon[1] * tensor[3] - gluon[3] * tensor[5])
        * prefactor,
        (-gluon[0] * tensor[2] + gluon[1] * tensor[4] + gluon[2] * tensor[5])
        * prefactor,
    )


def _expr_quark_vector_weyl(
    quark: tuple[Any, ...],
    vector: tuple[Any, ...],
    chirality: int,
) -> tuple[Any, ...]:
    tmp1, tmp2, tmp3, tmp4 = _expr_vector_slash_terms(vector)
    prefactor = 1j / math.sqrt(2.0)
    q1, q2 = quark
    if chirality == 1:
        return (
            prefactor * (tmp2 * q1 - tmp3 * q2),
            prefactor * (tmp1 * q2 - tmp4 * q1),
        )
    if chirality == -1:
        return (
            prefactor * (tmp1 * q1 + tmp3 * q2),
            prefactor * (tmp2 * q2 + tmp4 * q1),
        )
    raise NativeEvaluationError("Weyl quark-vector expression needs nonzero chirality")


def _expr_propagate(
    current: CurrentKey,
    value: tuple[Any, ...],
    momentum: tuple[Any, ...],
) -> tuple[Any, ...]:
    pdg = int(current.pdg)
    if pdg == 21:
        denominator = _expr_minkowski_square(momentum)
        prefactor = -1j / denominator
        return tuple(component * prefactor for component in value)
    if 1 <= abs(pdg) <= 6 and int(current.chirality) != 0:
        return _expr_quark_propagator_weyl(
            value,
            momentum,
            int(current.chirality),
        )
    return value


def _expr_quark_propagator_weyl(
    quark: tuple[Any, ...],
    momentum: tuple[Any, ...],
    chirality: int,
) -> tuple[Any, ...]:
    energy, px, py, pz = momentum
    denominator = _expr_minkowski_square(momentum)
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


def _expr_vector_slash_terms(
    vector: tuple[Any, ...],
) -> tuple[Any, Any, Any, Any]:
    v0, v1, v2, v3 = vector
    return v0 + v3, v0 - v3, v1 + 1j * v2, v1 - 1j * v2


def _expr_dot_weyl(
    left: tuple[Any, ...],
    right: tuple[Any, ...],
) -> Any:
    return left[0] * right[0] + left[1] * right[1]


def _expr_minkowski_dot(
    left: tuple[Any, ...],
    right: tuple[Any, ...],
) -> Any:
    return left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3]


def _expr_minkowski_dot_momentum(
    vector: tuple[Any, ...],
    momentum: tuple[Any, ...],
) -> Any:
    return vector[0] * momentum[0] - vector[1] * momentum[1] - vector[2] * momentum[2] - vector[3] * momentum[3]


def _expr_minkowski_square(momentum: tuple[Any, ...]) -> Any:
    return (
        momentum[0] ** 2
        - momentum[1] ** 2
        - momentum[2] ** 2
        - momentum[3] ** 2
    )


def _build_shared_helicity_current_table(
    graph: RecursionGraph,
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
    progress_callback: ProgressCallback | None = None,
) -> SharedCurrentTable:
    """Build an AmpliCol-style imode=1 current table for q q~ -> Z+n gluons.

    External helicities become source-current ids.  Internal currents are then
    keyed by the ordinary current identity plus the bitset of source currents
    feeding them, matching the conceptual role of legacy AmpliCol's ``ext_cur``.
    """

    process = graph.process
    if len(process) != gluon_count + 3:
        raise NativeEvaluationError("shared-helicity table received a mismatched graph")
    incoming_quark_pdg = int(process[0])
    incoming_antiquark_pdg = int(process[1])
    quark_current_pdg = -incoming_antiquark_pdg
    anti_closure_pdg = -incoming_quark_pdg
    gluon_labels = tuple(range(3, 3 + gluon_count))
    z_label = gluon_count + 3

    currents: list[SharedCurrentNode] = []
    sources: list[SharedSourceCurrent] = []
    interactions: list[SharedInteractionNode] = []
    current_ids: dict[tuple[int, tuple[int, ...], int, int], int] = {}

    def add_current(
        key: CurrentKey,
        *,
        source_ids: tuple[int, ...],
        ext_source_bits: int,
        is_source: bool = False,
        needs_propagator: bool = False,
    ) -> int:
        identity = (
            int(key.pdg),
            tuple(int(label) for label in key.external_labels),
            int(key.chirality),
            int(ext_source_bits),
        )
        existing = current_ids.get(identity)
        if existing is not None:
            return existing
        current_id = len(currents)
        currents.append(
            SharedCurrentNode(
                id=current_id,
                key=key,
                ext_source_bits=ext_source_bits,
                source_ids=source_ids,
                is_source=is_source,
                needs_propagator=needs_propagator,
                dimension=_current_dimension(key),
            )
        )
        current_ids[identity] = current_id
        return current_id

    def add_source(
        key: CurrentKey,
        *,
        leg_label: int,
        helicity: int,
        physical_helicity: int,
        chirality: int,
    ) -> int:
        source_bit = 1 << len(sources)
        current_id = add_current(
            key,
            source_ids=(len(sources),),
            ext_source_bits=source_bit,
            is_source=True,
            needs_propagator=False,
        )
        sources.append(
            SharedSourceCurrent(
                current_id=current_id,
                leg_label=leg_label,
                helicity=helicity,
                physical_helicity=physical_helicity,
                chirality=chirality,
                source_bit=source_bit,
            )
        )
        _report_progress(
            progress_callback,
            stage="recursion",
            item=f"source leg {leg_label}",
            increment=1,
        )
        return current_id

    def add_interaction(
        vertex_kind: int,
        left_id: int,
        right_id: int,
        result_id: int,
        coupling: tuple[float, float] = (1.0, 0.0),
    ) -> None:
        left = currents[left_id].key
        right = currents[right_id].key
        result = currents[result_id].key
        interactions.append(
            SharedInteractionNode(
                id=len(interactions),
                vertex_kind=vertex_kind,
                left_id=left_id,
                right_id=right_id,
                result_id=result_id,
                left=left,
                right=right,
                result=result,
                coupling=coupling,
            )
        )

    anti_source_by_physical: dict[int, int] = {
        1: add_source(
            CurrentKey(anti_closure_pdg, (1,), 0),
            leg_label=1,
            helicity=1,
            physical_helicity=1,
            chirality=1,
        ),
        -1: add_source(
            CurrentKey(anti_closure_pdg, (1,), 0),
            leg_label=1,
            helicity=-1,
            physical_helicity=-1,
            chirality=-1,
        ),
    }
    quark_source_by_chirality: dict[int, int] = {
        1: add_source(
            CurrentKey(quark_current_pdg, (2,), 1),
            leg_label=2,
            helicity=-1,
            physical_helicity=-1,
            chirality=1,
        ),
        -1: add_source(
            CurrentKey(quark_current_pdg, (2,), -1),
            leg_label=2,
            helicity=1,
            physical_helicity=1,
            chirality=-1,
        ),
    }
    gluon_sources: dict[tuple[int, int], int] = {}
    for label in gluon_labels:
        for helicity in (-1, 1):
            gluon_sources[(label, helicity)] = add_source(
                CurrentKey(21, (label,), 0),
                leg_label=label,
                helicity=helicity,
                physical_helicity=helicity,
                chirality=0,
            )
    z_sources: dict[int, int] = {}
    for helicity in (-1, 0, 1):
        z_sources[helicity] = add_source(
            CurrentKey(23, (z_label,), 0),
            leg_label=z_label,
            helicity=helicity,
            physical_helicity=helicity,
            chirality=0,
        )

    def merge_source_ids(*current_ids_to_merge: int) -> tuple[int, ...]:
        merged: set[int] = set()
        for current_id in current_ids_to_merge:
            merged.update(currents[current_id].source_ids)
        return tuple(sorted(merged))

    def merge_source_bits(*current_ids_to_merge: int) -> int:
        bits = 0
        for current_id in current_ids_to_merge:
            bits |= currents[current_id].ext_source_bits
        return bits

    gluon_currents: dict[tuple[int, int, tuple[int, ...]], int] = {}
    tensor_currents: dict[tuple[int, int, tuple[int, ...]], int] = {}
    for index, label in enumerate(gluon_labels):
        for helicity in (-1, 1):
            gluon_currents[(index, index + 1, (helicity,))] = gluon_sources[
                (label, helicity)
            ]

    for length in range(2, gluon_count + 1):
        for start in range(0, gluon_count - length + 1):
            end = start + length
            labels = gluon_labels[start:end]
            for helicities in product((-1, 1), repeat=length):
                _report_progress(
                    progress_callback,
                    stage="recursion",
                    item=f"gluon L{length} {start + 1}-{end}",
                    increment=1,
                )
                if length < gluon_count:
                    first_left = gluon_currents[(start, start + 1, helicities[:1])]
                    first_right = gluon_currents[(start + 1, end, helicities[1:])]
                    tensor_id = add_current(
                        CurrentKey(-21, labels, 0),
                        source_ids=merge_source_ids(first_left, first_right),
                        ext_source_bits=merge_source_bits(first_left, first_right),
                        needs_propagator=False,
                    )
                    for split in range(start + 1, end):
                        left_helicities = helicities[: split - start]
                        right_helicities = helicities[split - start :]
                        left_id = gluon_currents[(start, split, left_helicities)]
                        right_id = gluon_currents[(split, end, right_helicities)]
                        add_interaction(1, left_id, right_id, tensor_id)
                    tensor_currents[(start, end, helicities)] = tensor_id

                first_left = gluon_currents[(start, start + 1, helicities[:1])]
                first_right = gluon_currents[(start + 1, end, helicities[1:])]
                gluon_id = add_current(
                    CurrentKey(21, labels, 0),
                    source_ids=merge_source_ids(first_left, first_right),
                    ext_source_bits=merge_source_bits(first_left, first_right),
                    needs_propagator=True,
                )
                for split in range(start + 1, end):
                    left_helicities = helicities[: split - start]
                    right_helicities = helicities[split - start :]
                    left_id = gluon_currents[(start, split, left_helicities)]
                    right_id = gluon_currents[(split, end, right_helicities)]
                    add_interaction(0, left_id, right_id, gluon_id)
                    if split - start >= 2:
                        add_interaction(
                            2,
                            tensor_currents[(start, split, left_helicities)],
                            right_id,
                            gluon_id,
                        )
                    if end - split >= 2:
                        add_interaction(
                            3,
                            left_id,
                            tensor_currents[(split, end, right_helicities)],
                            gluon_id,
                        )
                gluon_currents[(start, end, helicities)] = gluon_id

    z_coupling = model.z_fermion_coupling(abs(incoming_quark_pdg))
    quark_without_z: dict[tuple[int, int, tuple[int, ...]], int] = {}
    quark_with_z: dict[tuple[int, int, tuple[int, ...], int], int] = {}
    for chirality in (1, -1):
        source_id = quark_source_by_chirality[chirality]
        quark_without_z[(chirality, 0, ())] = source_id
        for end in range(1, gluon_count + 1):
            labels = (2, *gluon_labels[:end])
            for helicities in product((-1, 1), repeat=end):
                _report_progress(
                    progress_callback,
                    stage="recursion",
                    item=f"quark g{end}/{gluon_count}",
                    increment=1,
                )
                first_left = quark_without_z[(chirality, 0, ())]
                first_right = gluon_currents[(0, end, helicities)]
                result_id = add_current(
                    CurrentKey(quark_current_pdg, labels, chirality),
                    source_ids=merge_source_ids(first_left, first_right),
                    ext_source_bits=merge_source_bits(first_left, first_right),
                    needs_propagator=True,
                )
                for split in range(0, end):
                    left_helicities = helicities[:split]
                    right_helicities = helicities[split:]
                    left_id = quark_without_z[(chirality, split, left_helicities)]
                    right_id = gluon_currents[(split, end, right_helicities)]
                    add_interaction(6, left_id, right_id, result_id)
                quark_without_z[(chirality, end, helicities)] = result_id

        for end in range(0, gluon_count + 1):
            labels = (2, *gluon_labels[:end], z_label)
            for helicities in product((-1, 1), repeat=end):
                without_z_id = quark_without_z[(chirality, end, helicities)]
                for z_helicity, z_id in z_sources.items():
                    _report_progress(
                        progress_callback,
                        stage="recursion",
                        item=f"quark+Z g{end}/{gluon_count}",
                        increment=1,
                    )
                    result_id = add_current(
                        CurrentKey(quark_current_pdg, labels, chirality),
                        source_ids=merge_source_ids(without_z_id, z_id),
                        ext_source_bits=merge_source_bits(without_z_id, z_id),
                        needs_propagator=end != gluon_count,
                    )
                    add_interaction(10, without_z_id, z_id, result_id, z_coupling)
                    for split in range(0, end):
                        left_helicities = helicities[:split]
                        right_helicities = helicities[split:]
                        left_id = quark_with_z[
                            (chirality, split, left_helicities, z_helicity)
                        ]
                        right_id = gluon_currents[(split, end, right_helicities)]
                        add_interaction(6, left_id, right_id, result_id)
                    quark_with_z[(chirality, end, helicities, z_helicity)] = result_id

    amplitudes: list[SharedAmplitudeRecord] = []
    for physical_quark_helicity, physical_antiquark_helicity, chirality in (
        (1, -1, 1),
        (-1, 1, -1),
    ):
        anti_id = anti_source_by_physical[physical_quark_helicity]
        for z_helicity in (-1, 0, 1):
            for gluon_helicities in product((-1, 1), repeat=gluon_count):
                _report_progress(
                    progress_callback,
                    stage="recursion",
                    item="amplitude closure",
                    increment=1,
                )
                amplitudes.append(
                    SharedAmplitudeRecord(
                        left_id=quark_with_z[
                            (chirality, gluon_count, gluon_helicities, z_helicity)
                        ],
                        right_id=anti_id,
                        helicities=(
                            physical_quark_helicity,
                            physical_antiquark_helicity,
                            *gluon_helicities,
                            z_helicity,
                        ),
                    )
                )

    interactions_by_result_lists: list[list[int]] = [[] for _ in currents]
    for interaction in interactions:
        interactions_by_result_lists[interaction.result_id].append(interaction.id)
    return SharedCurrentTable(
        currents=tuple(currents),
        sources=tuple(sources),
        interactions=tuple(interactions),
        interactions_by_result=tuple(
            tuple(ids) for ids in interactions_by_result_lists
        ),
        amplitudes=tuple(amplitudes),
    )


def _fill_shared_source_currents(
    values: np.ndarray,
    table: SharedCurrentTable,
    particles: tuple[ExternalMomentum, ...],
    model: AmplicolSMLeadingColorModel,
    *,
    gluon_count: int,
    zero_current_values: np.ndarray | None = None,
    clear_values: bool = True,
) -> None:
    if clear_values:
        if zero_current_values is not None:
            values[:, :] = zero_current_values
        else:
            values[:, :] = 0.0
    incoming_quark = particles[0]
    incoming_antiquark = particles[1]
    anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
    quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
    gluons = particles[2 : 2 + gluon_count]
    z_boson = particles[-1]

    for source in table.sources:
        current = table.currents[source.current_id]
        if source.leg_label == 1:
            wavefunction: Sequence[complex]
            if source.physical_helicity == 1:
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
            values[current.id, : len(wavefunction)] = wavefunction
        elif source.leg_label == 2:
            if source.chirality == 1:
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
            values[current.id, : len(wavefunction)] = wavefunction
        elif 3 <= source.leg_label < gluon_count + 3:
            gluon = gluons[source.leg_label - 3]
            wavefunction = _ext_gluon_cmplx(gluon.momentum, source.helicity)
            values[current.id, : len(wavefunction)] = wavefunction
        elif source.leg_label == gluon_count + 3:
            wavefunction = _ext_massive_vector(
                z_boson.momentum,
                source.helicity,
                model.mass(23),
            )
            values[current.id, : len(wavefunction)] = wavefunction
        else:
            raise NativeEvaluationError(
                f"unexpected source leg label: {source.leg_label}"
            )


def _zero_current_tuple(dimension: int) -> tuple[complex, ...]:
    if dimension == 2:
        return (0j, 0j)
    if dimension == 4:
        return (0j, 0j, 0j, 0j)
    if dimension == 6:
        return (0j, 0j, 0j, 0j, 0j, 0j)
    return tuple(0j for _ in range(dimension))


class _DAGBlockEvaluators:
    def __init__(self, model: AmplicolSMLeadingColorModel) -> None:
        self.model = model
        self.three_gluon = _build_three_gluon_block(model)
        self.two_gluon_to_tensor = _build_two_gluon_to_tensor_block(model)
        self.tensor_gluon_to_gluon = _build_tensor_gluon_to_gluon_block(model)
        self.gluon_tensor_to_gluon = _build_gluon_tensor_to_gluon_block(model)
        self.quark_vector = {
            1: _build_quark_vector_block(model, chirality=1),
            -1: _build_quark_vector_block(model, chirality=-1),
        }
        self.gluon_propagator = _build_gluon_propagator_block(model)
        self.quark_propagator = {
            1: _build_quark_propagator_block(model, chirality=1),
            -1: _build_quark_propagator_block(model, chirality=-1),
        }

    def vertex(
        self,
        interaction: Any,
        left: tuple[complex, ...],
        right: tuple[complex, ...],
        momenta_by_label: dict[int, tuple[float, float, float, float]],
    ) -> tuple[complex, ...]:
        kind = int(interaction.vertex_kind)
        if kind == 0:
            return self.three_gluon.evaluate(
                left,
                right,
                _sum_momenta(
                    momenta_by_label[label]
                    for label in interaction.left.external_labels
                ),
                _sum_momenta(
                    momenta_by_label[label]
                    for label in interaction.right.external_labels
                ),
            )
        if kind == 1:
            return self.two_gluon_to_tensor.evaluate(left, right)
        if kind == 2:
            return self.tensor_gluon_to_gluon.evaluate(left, right)
        if kind == 3:
            return self.gluon_tensor_to_gluon.evaluate(left, right)
        if kind == 6:
            return self.quark_vector[int(interaction.result.chirality)].evaluate(
                left,
                right,
            )
        if kind == 10:
            chirality = int(interaction.result.chirality)
            coupling = _weyl_coupling_for_chirality(
                chirality,
                (float(interaction.coupling[0]), float(interaction.coupling[1])),
            )
            return tuple(
                coupling * value
                for value in self.quark_vector[chirality].evaluate(left, right)
            )
        raise NativeEvaluationError(f"unsupported DAG vertex kind: {kind}")

    def propagate(
        self,
        current: Any,
        value: tuple[complex, ...],
        momentum: tuple[float, float, float, float],
    ) -> tuple[complex, ...]:
        pdg = int(current.pdg)
        if pdg == 21:
            return self.gluon_propagator.evaluate(value, momentum)
        if 1 <= abs(pdg) <= 6 and int(current.chirality) != 0:
            return self.quark_propagator[int(current.chirality)].evaluate(
                value,
                momentum,
            )
        return value


class _TensorBlockEvaluator:
    def __init__(
        self,
        *,
        evaluator: Any,
        param_builder: ParamBuilder,
        heads: tuple[tuple[str, ...], ...],
        output_length: int,
    ) -> None:
        self.evaluator = evaluator
        self.param_builder = param_builder
        self.heads = heads
        self.output_length = output_length

    def evaluate(
        self,
        *values: Sequence[complex | float],
    ) -> tuple[complex, ...]:
        for head, head_values in zip(self.heads, values, strict=True):
            self.param_builder.set_parameter_values(
                head,
                head_values,
                check_phase_flags=False,
            )
        tensor = self.evaluator.evaluate_complex(
            [self.param_builder.get_complex_values()]
        )[0]
        return tuple(complex(tensor[index]) for index in range(self.output_length))


def _build_three_gluon_block(model: AmplicolSMLeadingColorModel) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import Representation, TensorName, TensorNetwork

    library = model.build_tensor_library()
    builder = ParamBuilder()
    mink = Representation.mink(4)
    left_head = ("block", "three_gluon", "left")
    right_head = ("block", "three_gluon", "right")
    left_momentum_head = ("block", "three_gluon", "left_momentum")
    right_momentum_head = ("block", "three_gluon", "right_momentum")
    builder.register_rank1_tensor(
        library,
        tensor_name="pyamplicol::block_three_gluon_left",
        representation=mink,
        head=left_head,
        length=4,
        role="block_current",
    )
    builder.register_rank1_tensor(
        library,
        tensor_name="pyamplicol::block_three_gluon_right",
        representation=mink,
        head=right_head,
        length=4,
        role="block_current",
    )
    builder.register_rank1_tensor(
        library,
        tensor_name="pyamplicol::block_three_gluon_left_momentum",
        representation=mink,
        head=left_momentum_head,
        length=4,
        role="block_momentum",
        real_valued=True,
    )
    builder.register_rank1_tensor(
        library,
        tensor_name="pyamplicol::block_three_gluon_right_momentum",
        representation=mink,
        head=right_momentum_head,
        length=4,
        role="block_momentum",
        real_valued=True,
    )
    expression = (
        model.three_gluon_current_expression(
            left_slot=mink("left"),
            right_slot=mink("right"),
            output_slot=mink("out"),
            left_momentum_tensor_name="pyamplicol::block_three_gluon_left_momentum",
            right_momentum_tensor_name="pyamplicol::block_three_gluon_right_momentum",
            dummy_prefix="block_three_gluon",
        )
        * TensorName("pyamplicol::block_three_gluon_left")(
            mink("left")
        ).to_expression()
        * TensorName("pyamplicol::block_three_gluon_right")(
            mink("right")
        ).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(left_head, right_head, left_momentum_head, right_momentum_head),
        output_length=4,
    )


def _build_two_gluon_to_tensor_block(
    model: AmplicolSMLeadingColorModel,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import Representation, TensorName

    library = model.build_tensor_library()
    builder = ParamBuilder()
    mink = Representation.mink(4)
    aux = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    left_head = ("block", "two_gluon_to_tensor", "left")
    right_head = ("block", "two_gluon_to_tensor", "right")
    _register_block_vector(library, builder, "block_two_left", mink, left_head, 4)
    _register_block_vector(library, builder, "block_two_right", mink, right_head, 4)
    expression = (
        TensorName(str(symbols.two_gluon_to_tensor))(
            mink("left"),
            mink("right"),
            aux("out"),
        ).to_expression()
        * TensorName("pyamplicol::block_two_left")(mink("left")).to_expression()
        * TensorName("pyamplicol::block_two_right")(mink("right")).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(left_head, right_head),
        output_length=6,
    )


def _build_tensor_gluon_to_gluon_block(
    model: AmplicolSMLeadingColorModel,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import Representation, TensorName

    library = model.build_tensor_library()
    builder = ParamBuilder()
    mink = Representation.mink(4)
    aux = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    tensor_head = ("block", "tensor_gluon_to_gluon", "tensor")
    gluon_head = ("block", "tensor_gluon_to_gluon", "gluon")
    _register_block_vector(library, builder, "block_tg_tensor", aux, tensor_head, 6)
    _register_block_vector(library, builder, "block_tg_gluon", mink, gluon_head, 4)
    expression = (
        TensorName(str(symbols.tensor_gluon_to_gluon))(
            aux("tensor"),
            mink("gluon"),
            mink("out"),
        ).to_expression()
        * TensorName("pyamplicol::block_tg_tensor")(aux("tensor")).to_expression()
        * TensorName("pyamplicol::block_tg_gluon")(mink("gluon")).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(tensor_head, gluon_head),
        output_length=4,
    )


def _build_gluon_tensor_to_gluon_block(
    model: AmplicolSMLeadingColorModel,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import Representation, TensorName

    library = model.build_tensor_library()
    builder = ParamBuilder()
    mink = Representation.mink(4)
    aux = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    gluon_head = ("block", "gluon_tensor_to_gluon", "gluon")
    tensor_head = ("block", "gluon_tensor_to_gluon", "tensor")
    _register_block_vector(library, builder, "block_gt_gluon", mink, gluon_head, 4)
    _register_block_vector(library, builder, "block_gt_tensor", aux, tensor_head, 6)
    expression = (
        TensorName(str(symbols.gluon_tensor_to_gluon))(
            mink("gluon"),
            aux("tensor"),
            mink("out"),
        ).to_expression()
        * TensorName("pyamplicol::block_gt_gluon")(mink("gluon")).to_expression()
        * TensorName("pyamplicol::block_gt_tensor")(aux("tensor")).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(gluon_head, tensor_head),
        output_length=4,
    )


def _build_quark_vector_block(
    model: AmplicolSMLeadingColorModel,
    *,
    chirality: int,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import Representation, TensorName

    library = model.build_tensor_library()
    builder = ParamBuilder()
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    mink = Representation.mink(4)
    quark_head = ("block", "quark_vector", str(chirality), "quark")
    vector_head = ("block", "quark_vector", str(chirality), "vector")
    _register_block_vector(library, builder, "block_qv_quark", weyl, quark_head, 2)
    _register_block_vector(library, builder, "block_qv_vector", mink, vector_head, 4)
    tensor_name = (
        str(symbols.quark_vector_weyl_plus)
        if chirality == 1
        else str(symbols.quark_vector_weyl_minus)
    )
    expression = (
        TensorName(tensor_name)(
            weyl("quark"),
            mink("vector"),
            weyl("out"),
        ).to_expression()
        * TensorName("pyamplicol::block_qv_quark")(weyl("quark")).to_expression()
        * TensorName("pyamplicol::block_qv_vector")(mink("vector")).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(quark_head, vector_head),
        output_length=2,
    )


def _build_gluon_propagator_block(
    model: AmplicolSMLeadingColorModel,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import LibraryTensor, Representation, TensorName

    library = model.build_tensor_library()
    builder = ParamBuilder()
    mink = Representation.mink(4)
    gluon_head = ("block", "gluon_propagator", "gluon")
    momentum_head = ("block", "gluon_propagator", "momentum")
    _register_block_vector(library, builder, "block_gp_gluon", mink, gluon_head, 4)
    momentum = builder.register_rank1_tensor(
        library,
        tensor_name="pyamplicol::block_gp_momentum",
        representation=mink,
        head=momentum_head,
        length=4,
        role="block_momentum",
        real_valued=True,
    )
    library.register(
        LibraryTensor.dense(
            TensorName("pyamplicol::block_gluon_propagator")(mink, mink),
            model.gluon_propagator_tensor_data(momentum),
        )
    )
    expression = (
        TensorName("pyamplicol::block_gluon_propagator")(
            mink("input"),
            mink("out"),
        ).to_expression()
        * TensorName("pyamplicol::block_gp_gluon")(mink("input")).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(gluon_head, momentum_head),
        output_length=4,
    )


def _build_quark_propagator_block(
    model: AmplicolSMLeadingColorModel,
    *,
    chirality: int,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import LibraryTensor, Representation, TensorName

    library = model.build_tensor_library()
    builder = ParamBuilder()
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    mink = Representation.mink(4)
    quark_head = ("block", "quark_propagator", str(chirality), "quark")
    momentum_head = ("block", "quark_propagator", str(chirality), "momentum")
    _register_block_vector(library, builder, "block_qp_quark", weyl, quark_head, 2)
    momentum = builder.register_rank1_tensor(
        library,
        tensor_name="pyamplicol::block_qp_momentum",
        representation=mink,
        head=momentum_head,
        length=4,
        role="block_momentum",
        real_valued=True,
    )
    library.register(
        LibraryTensor.dense(
            TensorName("pyamplicol::block_quark_propagator")(weyl, weyl),
            model.quark_weyl_propagator_tensor_data(
                momentum,
                chirality=chirality,
            ),
        )
    )
    expression = (
        TensorName("pyamplicol::block_quark_propagator")(
            weyl("input"),
            weyl("out"),
        ).to_expression()
        * TensorName("pyamplicol::block_qp_quark")(weyl("input")).to_expression()
    )
    return _tensor_block_from_expression(
        expression,
        library,
        builder,
        heads=(quark_head, momentum_head),
        output_length=2,
    )


def _register_block_vector(
    library: Any,
    builder: ParamBuilder,
    name: str,
    representation: Any,
    head: tuple[str, ...],
    length: int,
) -> None:
    builder.register_rank1_tensor(
        library,
        tensor_name=f"pyamplicol::{name}",
        representation=representation,
        head=head,
        length=length,
        role="block_current",
    )


def _tensor_block_from_expression(
    expression: Any,
    library: Any,
    builder: ParamBuilder,
    *,
    heads: tuple[tuple[str, ...], ...],
    output_length: int,
) -> _TensorBlockEvaluator:
    from symbolica.community.spenso import TensorNetwork

    network = TensorNetwork(expression, library)
    network.execute(library=library)
    tensor = network.result_tensor(library)
    evaluator = tensor.evaluator(
        constants={},
        funs={},
        params=builder.parameter_symbols(),
    )
    set_real_params = getattr(evaluator, "set_real_params", None)
    if builder.real_valued_inputs and callable(set_real_params):
        set_real_params(builder.real_valued_inputs)
    return _TensorBlockEvaluator(
        evaluator=evaluator,
        param_builder=builder,
        heads=heads,
        output_length=output_length,
    )


def _z_gluon_graph(
    process: str,
    model: AmplicolSMLeadingColorModel,
) -> RecursionGraph:
    result = NativeMatrixElementGenerator(model=model).generate(
        process,
        write_cache_metadata=False,
    )
    if result.graph is None:
        raise NativeEvaluationError(f"no native graph available for {process}")
    return result.graph


def _validate_z_gluon_graph(graph: RecursionGraph) -> int:
    gluon_count = sum(1 for pdg in graph.process[2:] if pdg == 21)
    if gluon_count < 1 or graph.process[-1] != 23:
        raise NativeEvaluationError(
            "DAG evaluator currently supports q q~ -> Z plus ordered gluons"
        )
    if any(pdg != 21 for pdg in graph.process[2:-1]):
        raise NativeEvaluationError("DAG graph is not an ordered Z+gluon graph")
    return gluon_count


def _validate_z_gluon_point(
    particles: Sequence[ExternalMomentum],
    *,
    gluon_count: int,
) -> tuple[ExternalMomentum, ...]:
    point = tuple(particles)
    expected = gluon_count + 3
    if len(point) != expected:
        raise NativeEvaluationError(
            f"q q~ -> Z + {gluon_count} gluons requires exactly {expected} external momenta"
        )
    if point[0].pdg + point[1].pdg != 0 or not 1 <= abs(point[0].pdg) <= 6:
        raise NativeEvaluationError("incoming particles must be a quark/antiquark pair")
    if any(particle.pdg != 21 for particle in point[2 : 2 + gluon_count]):
        raise NativeEvaluationError("final-state gluons must precede the Z boson")
    if point[-1].pdg != 23:
        raise NativeEvaluationError("final state must end with one Z boson")
    return point


def _current_momenta_by_label(
    particles: tuple[ExternalMomentum, ...],
) -> dict[int, tuple[float, float, float, float]]:
    momenta: dict[int, tuple[float, float, float, float]] = {}
    for index, particle in enumerate(particles, start=1):
        momenta[index] = (
            _negate_momentum(particle.momentum)
            if index <= 2
            else particle.momentum
        )
    return momenta


def _sum_components(
    left: tuple[complex, ...],
    right: tuple[complex, ...],
) -> tuple[complex, ...]:
    return tuple(
        left_value + right_value
        for left_value, right_value in zip(left, right, strict=True)
    )


def _dot_weyl(
    left: tuple[complex, ...],
    right: tuple[complex, ...],
) -> complex:
    return left[0] * right[0] + left[1] * right[1]


def _sum_momenta(
    momenta: Iterable[tuple[float, float, float, float]],
) -> tuple[float, float, float, float]:
    total = (0.0, 0.0, 0.0, 0.0)
    for momentum in momenta:
        total = (
            total[0] + momentum[0],
            total[1] + momentum[1],
            total[2] + momentum[2],
            total[3] + momentum[3],
        )
    return total


def _negate_momentum(
    momentum: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    return -momentum[0], -momentum[1], -momentum[2], -momentum[3]


__all__ = [
    "DAGEvaluatorMetadata",
    "ZGluonDAGEvaluator",
]
