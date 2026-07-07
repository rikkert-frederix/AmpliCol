from __future__ import annotations

import math
import time
from dataclasses import dataclass, replace
from typing import Any, Literal, Sequence, cast

from .model import AmplicolSMLeadingColorModel
from .native import (
    ExternalMomentum,
    LeadingColorZJetsNativeEvaluator,
    MatrixElementEvaluation,
    NativeEvaluationError,
    _final_state_identical_factor,
    _initial_state_average_factor,
)

RuntimeBackend = Literal[
    "auto",
    "python",
    "dag",
    "numeric-tensor-network",
    "rusticol",
]


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
    """Reference-only selector for retired native/tensor/Z-specific evaluators.

    Production pyAmpliCol generation and evaluation are schema-v2 generic DAG
    process artifacts executed by Rusticol.  This class remains available only
    for migration diagnostics and historical reference tests, and callers must
    opt in explicitly with ``allow_reference_legacy=True``.
    """

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
        symbolica_jit_direct_translation: bool | None = None,
        symbolica_jit_optimization_level: int = 3,
        symbolica_max_horner_scheme_variables: int = 500,
        symbolica_max_common_pair_cache_entries: int = 1000000,
        symbolica_max_common_pair_distance: int | None = None,
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
        compiled_dag_lowering: str = "spenso",
        compiled_dag_cross_check_lowering: bool = False,
        compiled_dag_inline_external_wavefunctions: bool = True,
        compiled_dag_helicity_filter: bool = True,
        compiled_dag_helicity_filter_samples: int = 10,
        compiled_dag_helicity_filter_seed: int = 12345,
        compiled_dag_helicity_filter_relative_tolerance: float = 1.0e-12,
        compiled_dag_helicity_filter_zero_tolerance: float = 1.0e-300,
        compiled_dag_helicity_filter_phase_space: str = "rambo",
        allow_reference_legacy: bool = False,
    ) -> None:
        if not allow_reference_legacy:
            raise NativeEvaluationError(
                "NativeRuntimeEvaluator is a retired reference-only runtime. "
                "Production pyAmpliCol uses generic DAG process artifacts: run "
                "`pyamplicol generate-process PROCESS OUTPUT_DIR` and evaluate "
                "the resulting directory with Rusticol or `pyamplicol time-process`."
            )
        setup_start = time.perf_counter()
        self.process = process
        self.model = model or AmplicolSMLeadingColorModel()
        self._python = LeadingColorZJetsNativeEvaluator(self.model)
        self._runtime: Any | None = None
        self.batch_size = batch_size
        self._neutral_vector_target = (
            self._python.supported_electroweak_vector_gluon_process(process)
        )
        self._neutral_dilepton_target = (
            self._python.supported_neutral_dilepton_gluon_process(process)
        )
        self._charged_leptonic_w_target = (
            self._python.supported_charged_leptonic_w_gluon_process(process)
        )
        self._vector_pdg = (
            None if self._neutral_vector_target is None else self._neutral_vector_target[0]
        )
        self._gluon_count = (
            None if self._neutral_vector_target is None else self._neutral_vector_target[1]
        )
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
        elif self._neutral_vector_target is not None and self._gluon_count is not None:
            if self._vector_pdg != 23:
                if runtime_backend in (
                    "numeric-tensor-network",
                    "rusticol",
                ):
                    raise NativeEvaluationError(
                        f"runtime backend {runtime_backend!r} is not yet available for "
                        "one-quark-line non-Z vector plus ordered gluons"
                    )
                backend = (
                    "native-python-gamma-gluon"
                    if self._vector_pdg == 22
                    else "native-python-w-gluon"
                )
                kernel = "staged-python-recursion"
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
            elif runtime_backend == "rusticol":
                if symbolica_load_evaluator_dir is None:
                    raise NativeEvaluationError(
                        "runtime backend 'rusticol' requires --load-evaluator-dir "
                        "pointing at a generated process directory"
                    )
                import rusticol  # type: ignore[import-not-found]

                self._runtime = rusticol.Runtime.load(  # type: ignore[attr-defined]
                    str(symbolica_load_evaluator_dir)
                )
                backend = "rusticol-pyo3-shared-current-dag"
                kernel = "rusticol-pyo3-shared-current-dag"
                evaluator_metadata = dict(self._runtime.metadata())
            elif self._vector_pdg != 23 and runtime_backend == "auto":
                pass
            else:
                try:
                    from .dag_runtime import ZGluonDAGEvaluator

                    dag_common_pair_distance = (
                        100
                        if symbolica_max_common_pair_distance is None
                        else symbolica_max_common_pair_distance
                    )
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
                        symbolica_jit_direct_translation=(
                            False
                            if symbolica_jit_direct_translation is None
                            else symbolica_jit_direct_translation
                        ),
                        symbolica_jit_optimization_level=symbolica_jit_optimization_level,
                        symbolica_max_horner_scheme_variables=(
                            symbolica_max_horner_scheme_variables
                        ),
                        symbolica_max_common_pair_cache_entries=(
                            symbolica_max_common_pair_cache_entries
                        ),
                        symbolica_max_common_pair_distance=dag_common_pair_distance,
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
        elif self._neutral_dilepton_target is not None:
            if runtime_backend not in ("auto", "python"):
                raise NativeEvaluationError(
                    f"runtime backend {runtime_backend!r} is not yet available for "
                    "one-quark-line neutral dilepton plus ordered gluons"
                )
            self._gluon_count = self._neutral_dilepton_target[2]
            backend = "native-python-neutral-dilepton-gluon"
            kernel = "staged-python-recursion"
        elif self._charged_leptonic_w_target is not None:
            if runtime_backend not in ("auto", "python"):
                raise NativeEvaluationError(
                    f"runtime backend {runtime_backend!r} is not yet available for "
                    "one-quark-line charged-current leptonic W plus ordered gluons"
                )
            self._gluon_count = self._charged_leptonic_w_target[3]
            backend = "native-python-charged-leptonic-w-gluon"
            kernel = "staged-python-recursion"
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
                "one-quark-line electroweak vector, neutral dilepton, or "
                "charged-current leptonic W plus ordered gluons"
            )
        point = tuple(particles) if particles is not None else self._canonical_point(sqrt_s)
        if self.metadata.backend == "rusticol-pyo3-shared-current-dag":
            return self._evaluate_rusticol_many((point,))[0]
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
                "one-quark-line electroweak vector, neutral dilepton, or "
                "charged-current leptonic W plus ordered gluons"
            )
        if self.metadata.backend == "rusticol-pyo3-shared-current-dag":
            return self._evaluate_rusticol_many(points)
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

    def materialize_compiled_evaluators(self) -> None:
        if self._runtime is None:
            return
        materialize = getattr(self._runtime, "materialize_compiled_evaluators", None)
        if callable(materialize):
            materialize()

    def save_evaluator_artifact(self, output_dir: str) -> Any:
        if self._runtime is None:
            raise NativeEvaluationError(
                "selected native runtime does not expose a reusable evaluator artifact"
            )
        save = getattr(self._runtime, "save_evaluator_artifact", None)
        if not callable(save):
            raise NativeEvaluationError(
                "selected native runtime does not support evaluator artifact saving"
            )
        return save(output_dir)

    def refresh_metadata(self) -> NativeRuntimeMetadata:
        if self._runtime is None:
            return self.metadata
        runtime_metadata = getattr(self._runtime, "metadata", None)
        if runtime_metadata is None:
            return self.metadata
        to_json = getattr(runtime_metadata, "to_json_dict", None)
        if not callable(to_json):
            return self.metadata
        self.metadata = replace(
            self.metadata,
            evaluator_metadata=to_json(),
        )
        return self.metadata

    def stage_diagnostics_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> dict[str, object]:
        if self._runtime is None:
            raise NativeEvaluationError(
                "selected native runtime does not expose stage diagnostics"
            )
        if self.metadata.backend == "rusticol-pyo3-shared-current-dag":
            return dict(self._runtime.stage_diagnostics(_momenta_array(points)))
        diagnostics = getattr(self._runtime, "stage_diagnostics_many", None)
        if not callable(diagnostics):
            raise NativeEvaluationError(
                "selected native runtime does not expose stage diagnostics"
            )
        return cast(dict[str, object], diagnostics(points))

    def _canonical_point(
        self,
        sqrt_s: float | None,
    ) -> tuple[ExternalMomentum, ...]:
        if self._python.supports_zero_gluon_z(self.process):
            return self._python.canonical_zero_gluon_point(
                self.process,
                sqrt_s=sqrt_s,
            )
        if self._gluon_count is not None and self._gluon_count >= 0:
            if self._vector_pdg is None:
                if self._neutral_dilepton_target is not None:
                    lepton_pdg, antilepton_pdg, gluon_count = (
                        self._neutral_dilepton_target
                    )
                    return self._python.canonical_neutral_dilepton_gluon_point(
                        self.process,
                        lepton_pdg=lepton_pdg,
                        antilepton_pdg=antilepton_pdg,
                        gluon_count=gluon_count,
                        sqrt_s=sqrt_s,
                    )
                if self._charged_leptonic_w_target is not None:
                    (
                        vector_pdg,
                        first_lepton_pdg,
                        second_lepton_pdg,
                        gluon_count,
                    ) = self._charged_leptonic_w_target
                    return self._python.canonical_charged_leptonic_w_gluon_point(
                        self.process,
                        vector_pdg=vector_pdg,
                        first_lepton_pdg=first_lepton_pdg,
                        second_lepton_pdg=second_lepton_pdg,
                        gluon_count=gluon_count,
                        sqrt_s=sqrt_s,
                    )
                raise NativeEvaluationError(
                    "native numerical evaluation is currently implemented only for "
                    "one-quark-line electroweak vector, neutral dilepton, or "
                    "charged-current leptonic W plus ordered gluons"
                )
            return self._python.canonical_neutral_vector_gluon_point(
                self.process,
                vector_pdg=self._vector_pdg,
                gluon_count=self._gluon_count,
                sqrt_s=sqrt_s,
            )
        raise NativeEvaluationError(
            "native numerical evaluation is currently implemented only for "
            "one-quark-line electroweak vector, neutral dilepton, or charged-current "
            "leptonic W plus ordered gluons"
        )

    def _evaluate_rusticol_many(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> tuple[MatrixElementEvaluation, ...]:
        if self._runtime is None:
            raise NativeEvaluationError("rusticol runtime is not initialized")
        values = self._runtime.evaluate(_momenta_array(points))
        output: list[MatrixElementEvaluation] = []
        for point, value in zip(points, values, strict=True):
            particles = tuple(point)
            pdgs = tuple(particle.pdg for particle in particles)
            color_factor = self.model.leading_color_factor(pdgs)
            average_factor = _initial_state_average_factor(pdgs[:2])
            identical_factor = _final_state_identical_factor(
                particle.pdg for particle in particles[2:]
            )
            gluon_count = 0 if self._gluon_count is None else self._gluon_count
            coupling_factor = (
                (4.0 * math.pi * self.model.alpha_s_me_check) ** gluon_count
                * (2.0 * 4.0 * math.pi * self.model.alpha_ew)
            )
            output.append(
                MatrixElementEvaluation(
                    process=self.process,
                    particles=particles,
                    matrix_element=float(value),
                    raw_helicity_sum=float("nan"),
                    color_factor=color_factor,
                    average_factor=average_factor,
                    coupling_factor=coupling_factor,
                    helicity_contributions=(),
                    identical_factor=identical_factor,
                )
            )
        return tuple(output)


def _momenta_array(points: Sequence[Sequence[ExternalMomentum]]) -> Any:
    import numpy as np

    return np.asarray(
        [
            [
                [float(component) for component in particle.momentum]
                for particle in point
            ]
            for point in points
        ],
        dtype=np.float64,
    )


__all__ = [
    "NativeRuntimeEvaluator",
    "NativeRuntimeMetadata",
    "RuntimeBackend",
]
