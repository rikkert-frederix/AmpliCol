from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Literal, Sequence, cast

from .model import AmplicolSMLeadingColorModel
from .native import (
    ExternalMomentum,
    LeadingColorZJetsNativeEvaluator,
    MatrixElementEvaluation,
    NativeEvaluationError,
)

RuntimeBackend = Literal["auto", "python", "dag", "numeric-tensor-network"]


@dataclass(frozen=True)
class NativeRuntimeMetadata:
    process: str
    backend: str
    kernel: str
    setup_time_s: float
    fallback_reason: str | None = None
    evaluator_metadata: dict[str, object] | None = None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "backend": self.backend,
            "kernel": self.kernel,
            "setup_time_s": self.setup_time_s,
            "fallback_reason": self.fallback_reason,
            "evaluator_metadata": self.evaluator_metadata,
        }


class NativeRuntimeEvaluator:
    """Select the strongest available native evaluator for a process."""

    def __init__(
        self,
        process: str,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        runtime_backend: RuntimeBackend = "auto",
        allow_python_fallback: bool = True,
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
        symbolica_compiled_output_dir: str | None = None,
        symbolica_load_evaluator_dir: str | None = None,
        symbolica_raw_sum_final_stage: bool = False,
    ) -> None:
        setup_start = time.perf_counter()
        self.process = process
        self.model = model or AmplicolSMLeadingColorModel()
        self._python = LeadingColorZJetsNativeEvaluator(self.model)
        self._runtime: Any | None = None
        self.batch_size = batch_size
        self._gluon_count = self._python.supported_z_gluon_count(process)
        fallback_reason = None
        backend = "native-python-z-gluon"
        kernel = "staged-python-recursion"
        evaluator_metadata: dict[str, object] | None = None

        if self._python.supports_zero_gluon_z(process):
            if runtime_backend not in ("auto", "python"):
                raise NativeEvaluationError(
                    f"runtime backend {runtime_backend!r} is not available for q q~ -> Z"
                )
            backend = "native-python-zero-gluon"
            kernel = "python-zero-gluon"
        elif self._gluon_count is not None and self._gluon_count >= 1:
            if runtime_backend == "python":
                pass
            elif runtime_backend == "numeric-tensor-network":
                from .tensor_runtime import ZGluonNumericTensorNetworkEvaluator

                self._runtime = ZGluonNumericTensorNetworkEvaluator(
                    process,
                    model=self.model,
                )
                numeric_metadata = self._runtime.metadata
                backend = "native-spenso-numeric-tensor-network"
                kernel = numeric_metadata.kernel
                evaluator_metadata = numeric_metadata.to_json_dict()
            else:
                try:
                    from .dag_runtime import ZGluonDAGEvaluator

                    self._runtime = ZGluonDAGEvaluator(
                        process,
                        model=self.model,
                        batch_size=batch_size,
                        merge_evaluators_strategy=merge_evaluators_strategy,
                        split_vertex_current_stages=split_vertex_current_stages,
                        verbose_evaluator_build=verbose_evaluator_build,
                        symbolica_evaluator_backend=symbolica_evaluator_backend,
                        symbolica_iterations=symbolica_iterations,
                        symbolica_cpe_iterations=symbolica_cpe_iterations,
                        symbolica_n_cores=symbolica_n_cores,
                        symbolica_direct_translation=symbolica_direct_translation,
                        symbolica_jit_direct_translation=symbolica_jit_direct_translation,
                        symbolica_jit_optimization_level=symbolica_jit_optimization_level,
                        symbolica_max_horner_scheme_variables=(
                            symbolica_max_horner_scheme_variables
                        ),
                        symbolica_max_common_pair_cache_entries=(
                            symbolica_max_common_pair_cache_entries
                        ),
                        symbolica_max_common_pair_distance=(
                            symbolica_max_common_pair_distance
                        ),
                        symbolica_collect_factors=symbolica_collect_factors,
                        symbolica_compiled_preset=symbolica_compiled_preset,
                        symbolica_compiled_inline_asm=symbolica_compiled_inline_asm,
                        symbolica_compiled_optimization_level=(
                            symbolica_compiled_optimization_level
                        ),
                        symbolica_compiled_native=symbolica_compiled_native,
                        symbolica_compiler_path=symbolica_compiler_path,
                        symbolica_compiler_flags=symbolica_compiler_flags,
                        symbolica_compiled_output_chunk_size=(
                            symbolica_compiled_output_chunk_size
                        ),
                        symbolica_compiled_chunk_compile_workers=(
                            symbolica_compiled_chunk_compile_workers
                        ),
                        symbolica_compiled_output_dir=(
                            None
                            if symbolica_compiled_output_dir is None
                            else str(symbolica_compiled_output_dir)
                        ),
                        symbolica_load_evaluator_dir=(
                            None
                            if symbolica_load_evaluator_dir is None
                            else str(symbolica_load_evaluator_dir)
                        ),
                        symbolica_raw_sum_final_stage=symbolica_raw_sum_final_stage,
                    )
                    dag_metadata = self._runtime.metadata
                    backend = "native-spenso-symbolica-shared-helicity-current-dag"
                    kernel = dag_metadata.kernel
                    evaluator_metadata = dag_metadata.to_json_dict()
                except NativeEvaluationError as exc:
                    if runtime_backend == "dag" or not allow_python_fallback:
                        raise
                    fallback_reason = str(exc)
        else:
            backend = "not-implemented"
            kernel = "not-implemented"

        self.metadata = NativeRuntimeMetadata(
            process=process,
            backend=backend,
            kernel=kernel,
            setup_time_s=time.perf_counter() - setup_start,
            fallback_reason=fallback_reason,
            evaluator_metadata=evaluator_metadata,
        )

    def evaluate(
        self,
        *,
        particles: Sequence[ExternalMomentum] | None = None,
        sqrt_s: float | None = None,
    ) -> MatrixElementEvaluation:
        if self.metadata.backend == "not-implemented":
            raise NativeEvaluationError(
                "native numerical evaluation is currently implemented only for "
                "q q~ -> Z plus ordered gluons"
            )
        point = tuple(particles) if particles is not None else self._canonical_point(sqrt_s)
        if self._runtime is not None:
            return cast(MatrixElementEvaluation, self._runtime.evaluate(point))
        return self._python.evaluate(self.process, particles=point)

    def evaluate_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> tuple[MatrixElementEvaluation, ...]:
        if self.metadata.backend == "not-implemented":
            raise NativeEvaluationError(
                "native numerical evaluation is currently implemented only for "
                "q q~ -> Z plus ordered gluons"
            )
        if self._runtime is not None and hasattr(self._runtime, "evaluate_many"):
            return cast(
                tuple[MatrixElementEvaluation, ...],
                self._runtime.evaluate_many(points),
            )
        results: list[MatrixElementEvaluation] = []
        for start in range(0, len(points), self.batch_size):
            batch = points[start : start + self.batch_size]
            for point in batch:
                results.append(self.evaluate(particles=point))
        return tuple(results)

    def _canonical_point(
        self,
        sqrt_s: float | None,
    ) -> tuple[ExternalMomentum, ...]:
        if self._python.supports_zero_gluon_z(self.process):
            return self._python.canonical_zero_gluon_point(
                self.process,
                sqrt_s=sqrt_s,
            )
        if self._gluon_count is not None and self._gluon_count >= 1:
            return self._python.canonical_z_gluon_point(
                self.process,
                gluon_count=self._gluon_count,
                sqrt_s=sqrt_s,
            )
        raise NativeEvaluationError(
            "native numerical evaluation is currently implemented only for "
            "q q~ -> Z plus ordered gluons"
        )


__all__ = [
    "NativeRuntimeEvaluator",
    "NativeRuntimeMetadata",
    "RuntimeBackend",
]
