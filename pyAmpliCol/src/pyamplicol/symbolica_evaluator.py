from __future__ import annotations

import os
import platform
import shutil
import tempfile
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np

from .core_types import NativeEvaluationError


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


@dataclass(frozen=True)
class SymbolicaEvaluatorSettings:
    backend: str = "jit"
    iterations: int = 10
    cpe_iterations: int | None = None
    n_cores: int = 4
    direct_translation: bool = True
    jit_direct_translation: bool = False
    jit_optimization_level: int = 3
    max_horner_scheme_variables: int = 1000
    max_common_pair_cache_entries: int = 5000000
    max_common_pair_distance: int = 1000
    collect_factors: bool = False
    compiled_preset: str = "adaptive"
    compiled_inline_asm: str = "default"
    compiled_optimization_level: int = 3
    compiled_native: bool = True
    compiler_path: str | None = None
    compiler_flags: tuple[str, ...] = ()
    compiled_output_chunk_size: int | None = None
    output_chunk_strategy: str = "uniform"
    output_chunk_autotune_batch_size: int = 128
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
        if self.output_chunk_strategy not in (
            "auto",
            "uniform",
            "tapered-stage",
            "measured-stage",
        ):
            raise NativeEvaluationError(
                "symbolica output chunk strategy must be 'auto', 'uniform', "
                "'tapered-stage', or 'measured-stage'"
            )
        if self.output_chunk_autotune_batch_size < 1:
            raise NativeEvaluationError(
                "symbolica output chunk autotune batch size must be positive"
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
            "output_chunk_strategy": self.output_chunk_strategy,
            "output_chunk_autotune_batch_size": (
                self.output_chunk_autotune_batch_size
            ),
            "compiled_chunk_compile_workers": self.compiled_chunk_compile_workers,
            "compiled_output_dir": self.compiled_output_dir,
            "raw_sum_final_stage": self.raw_sum_final_stage,
            "real_param_sqrt_real": self.real_param_sqrt_real,
            "real_param_log_real": self.real_param_log_real,
            "real_param_powf_real": self.real_param_powf_real,
            "real_param_real_if_args_real": self.real_param_real_if_args_real,
        }


def _compile_symbolica_outputs(
    outputs: tuple[Any, ...],
    params: list[Any],
    *,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    aliases: Sequence[tuple[Any, Any]] = (),
    functions: Mapping[tuple[Any, tuple[Any, ...]], Any] | None = None,
    real_params: Sequence[int] = (),
    symbolica_settings: SymbolicaEvaluatorSettings | None = None,
    jit_compile: bool = True,
    label: str = "symbolica",
    progress_callback: ProgressCallback | None = None,
) -> Any:
    if not outputs:
        raise NativeEvaluationError("cannot build evaluator with zero outputs")
    settings = symbolica_settings or SymbolicaEvaluatorSettings()
    progress_stage = _symbolica_progress_stage(settings, jit_compile=jit_compile)
    _report_progress(
        progress_callback,
        stage=progress_stage,
        item=f"{label} prepare {len(outputs)}",
    )
    total_started = time.perf_counter()
    prepare_started = time.perf_counter()
    outputs = tuple(_prepare_symbolica_output(output, settings) for output in outputs)
    output_prepare_s = time.perf_counter() - prepare_started
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
            _report_progress(
                progress_callback,
                stage=progress_stage,
                item=f"{label} chunk {chunk_index + 1}/{len(chunk_ranges)}",
            )
            return _compile_symbolica_outputs(
                outputs[start:stop],
                params,
                merge_evaluators_strategy=merge_evaluators_strategy,
                verbose_evaluator_build=verbose_evaluator_build,
                aliases=aliases,
                functions=functions,
                real_params=real_params,
                symbolica_settings=unchunked_settings,
                jit_compile=jit_compile,
                label=f"{label}_chunk_{chunk_index}",
                progress_callback=progress_callback,
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
        chunked = _ChunkedSymbolicaEvaluator(tuple(chunks))
        chunked.build_timing = {
            "output_prepare_s": output_prepare_s,
            "symbolica_evaluator_build_s": time.perf_counter() - total_started,
        }
        return chunked
    evaluator_kwargs = _symbolica_evaluator_kwargs(
        settings,
        verbose=verbose_evaluator_build,
        jit_compile=jit_compile,
    )
    alias_kwargs = {"aliases": list(aliases)} if aliases else {}
    function_kwargs = {"functions": dict(functions)} if functions else {}
    if merge_evaluators_strategy:
        _report_progress(
            progress_callback,
            stage=progress_stage,
            item=f"{label} evaluator 1/{len(outputs)}",
        )
        _report_jit_boundary(
            progress_callback,
            settings,
            jit_compile=jit_compile,
            phase="initialize",
            item=f"{label} eval 1/{len(outputs)} p={len(params)}",
        )
        with _JITBoundaryHeartbeat(
            progress_callback,
            settings,
            jit_compile=jit_compile,
            phase="initialize",
            item=f"{label} eval 1/{len(outputs)} p={len(params)}",
        ):
            evaluator_started = time.perf_counter()
            evaluator = outputs[0].evaluator(
                params,
                **alias_kwargs,
                **function_kwargs,
                **evaluator_kwargs,
            )
            evaluator_construct_s = time.perf_counter() - evaluator_started
        _report_jit_boundary(
            progress_callback,
            settings,
            jit_compile=jit_compile,
            phase="returned",
            item=f"{label} eval 1/{len(outputs)}",
        )
        for expression in _progress_outputs(
            outputs[1:],
            enabled=verbose_evaluator_build,
        ):
            _report_progress(
                progress_callback,
                stage=progress_stage,
                item=f"{label} merge",
            )
            _report_jit_boundary(
                progress_callback,
                settings,
                jit_compile=jit_compile,
                phase="initialize",
                item=f"{label} merge p={len(params)}",
            )
            with _JITBoundaryHeartbeat(
                progress_callback,
                settings,
                jit_compile=jit_compile,
                phase="initialize",
                item=f"{label} merge p={len(params)}",
            ):
                merge_construct_started = time.perf_counter()
                other = expression.evaluator(
                    params,
                    **alias_kwargs,
                    **function_kwargs,
                    **evaluator_kwargs,
                )
                evaluator_construct_s += time.perf_counter() - merge_construct_started
            _report_jit_boundary(
                progress_callback,
                settings,
                jit_compile=jit_compile,
                phase="returned",
                item=f"{label} merge",
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
            real_params_started = time.perf_counter()
            _report_jit_boundary(
                progress_callback,
                settings,
                jit_compile=jit_compile,
                phase="real params",
                item=f"{label} real={len(real_params)}",
            )
            evaluator.set_real_params(
                list(real_params),
                sqrt_real=settings.real_param_sqrt_real,
                log_real=settings.real_param_log_real,
                powf_real=settings.real_param_powf_real,
                real_if_args_real=settings.real_param_real_if_args_real,
                verbose=verbose_evaluator_build,
            )
            real_params_s = time.perf_counter() - real_params_started
        else:
            real_params_s = 0.0
        adapter = _finalize_symbolica_evaluator(
            evaluator,
            settings,
            label,
            input_len=len(params),
            output_len=len(outputs),
            progress_callback=progress_callback,
        )
        _set_evaluator_build_timing(
            adapter,
            {
                "output_prepare_s": output_prepare_s,
                "evaluator_construct_s": evaluator_construct_s,
                "real_params_s": real_params_s,
                "symbolica_evaluator_build_s": time.perf_counter() - total_started,
            },
        )
        return adapter

    from symbolica import Expression

    _report_progress(
        progress_callback,
        stage=progress_stage,
        item=f"{label} evaluator {len(outputs)}",
    )
    _report_jit_boundary(
        progress_callback,
        settings,
        jit_compile=jit_compile,
        phase="initialize",
        item=f"{label} out={len(outputs)} p={len(params)}",
    )
    with _JITBoundaryHeartbeat(
        progress_callback,
        settings,
        jit_compile=jit_compile,
        phase="initialize",
        item=f"{label} out={len(outputs)} p={len(params)}",
    ):
        evaluator_started = time.perf_counter()
        evaluator = Expression.evaluator_multiple(
            outputs,
            params,
            **alias_kwargs,
            **function_kwargs,
            **evaluator_kwargs,
        )
        evaluator_construct_s = time.perf_counter() - evaluator_started
    _report_jit_boundary(
        progress_callback,
        settings,
        jit_compile=jit_compile,
        phase="returned",
        item=f"{label} out={len(outputs)}",
    )
    if real_params:
        real_params_started = time.perf_counter()
        _report_jit_boundary(
            progress_callback,
            settings,
            jit_compile=jit_compile,
            phase="real params",
            item=f"{label} real={len(real_params)}",
        )
        evaluator.set_real_params(
            list(real_params),
            sqrt_real=settings.real_param_sqrt_real,
            log_real=settings.real_param_log_real,
            powf_real=settings.real_param_powf_real,
            real_if_args_real=settings.real_param_real_if_args_real,
            verbose=verbose_evaluator_build,
        )
        real_params_s = time.perf_counter() - real_params_started
    else:
        real_params_s = 0.0
    adapter = _finalize_symbolica_evaluator(
        evaluator,
        settings,
        label,
        input_len=len(params),
        output_len=len(outputs),
        progress_callback=progress_callback,
    )
    _set_evaluator_build_timing(
        adapter,
        {
            "output_prepare_s": output_prepare_s,
            "evaluator_construct_s": evaluator_construct_s,
            "real_params_s": real_params_s,
            "symbolica_evaluator_build_s": time.perf_counter() - total_started,
        },
    )
    return adapter


def _set_evaluator_build_timing(evaluator: Any, timing: Mapping[str, float]) -> None:
    try:
        current = getattr(evaluator, "build_timing", {})
        if not isinstance(current, dict):
            current = {}
        current.update({str(key): float(value) for key, value in timing.items()})
        setattr(evaluator, "build_timing", current)
    except Exception:
        return


def _report_jit_boundary(
    progress_callback: ProgressCallback | None,
    settings: SymbolicaEvaluatorSettings,
    *,
    jit_compile: bool,
    phase: str,
    item: str,
) -> None:
    if settings.backend != "jit" or not jit_compile:
        return
    _report_progress(
        progress_callback,
        stage=f"jit {phase}",
        item=item,
    )


class _JITBoundaryHeartbeat:
    """Emit progress while Symbolica is inside a blocking JIT build call."""

    def __init__(
        self,
        progress_callback: ProgressCallback | None,
        settings: SymbolicaEvaluatorSettings,
        *,
        jit_compile: bool,
        phase: str,
        item: str,
        interval_s: float = 5.0,
    ) -> None:
        self.progress_callback = progress_callback
        self.settings = settings
        self.jit_compile = jit_compile
        self.phase = phase
        self.item = item
        self.interval_s = max(float(interval_s), 0.01)
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._started_at = 0.0

    def __enter__(self) -> "_JITBoundaryHeartbeat":
        if (
            self.progress_callback is None
            or self.settings.backend != "jit"
            or not self.jit_compile
        ):
            return self
        self._started_at = time.perf_counter()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)

    def _run(self) -> None:
        while not self._stop.wait(self.interval_s):
            elapsed_s = time.perf_counter() - self._started_at
            _report_progress(
                self.progress_callback,
                stage=f"jit {self.phase}",
                item=f"{self.item} waiting {elapsed_s:.0f}s",
            )


def _symbolica_progress_stage(
    settings: SymbolicaEvaluatorSettings,
    *,
    jit_compile: bool,
) -> str:
    if settings.backend == "jit":
        return "jit compile" if jit_compile else "jit build"
    if settings.backend in ("compiled-complex", "compiled-complex-4x"):
        return "c++ build"
    return "eval build"


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
    progress_callback: ProgressCallback | None = None,
) -> Any:
    if settings.backend == "jit":
        _report_progress(
            progress_callback,
            stage="jit ready",
            item=label,
        )
        return _JITSymbolicaEvaluatorAdapter(
            evaluator,
            settings,
            label,
            input_len=input_len,
            output_len=output_len,
            progress_callback=progress_callback,
        )
    if settings.backend in ("compiled-complex", "compiled-complex-4x"):
        _report_progress(
            progress_callback,
            stage="c++ compile",
            item=label,
        )
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
        progress_callback: ProgressCallback | None = None,
    ) -> None:
        self.input_len = int(input_len)
        self.output_len = int(output_len)
        self.backend = settings.backend
        self.settings = settings.to_json_dict()
        self.label = _safe_symbol_name(label)
        self._source_evaluator = evaluator
        self._progress_callback = progress_callback
        self.evaluator_state_path: Path | None = None
        self.build_timing: dict[str, float] = {}

    def evaluate_complex(self, parameter_rows: Any) -> Any:
        return self._source_evaluator.evaluate_complex(parameter_rows)

    def evaluate_complex_profiled(
        self, parameter_rows: Any
    ) -> tuple[Any, tuple[float, float, float]]:
        return self._evaluate_complex_profiled_prepared(
            _complex128_parameter_rows(parameter_rows)
        )

    def supports_complex_profiled(self) -> bool:
        return callable(
            getattr(self._source_evaluator, "evaluate_complex_profiled", None)
        )

    def evaluate(self, parameter_rows: Any) -> Any:
        return self._source_evaluator.evaluate(parameter_rows)

    def _evaluate_complex_prepared(self, parameter_rows: np.ndarray) -> Any:
        return self._source_evaluator.evaluate_complex(parameter_rows)

    def _evaluate_complex_profiled_prepared(
        self, parameter_rows: np.ndarray
    ) -> tuple[Any, tuple[float, float, float]]:
        profile = getattr(self._source_evaluator, "evaluate_complex_profiled", None)
        if not callable(profile):
            raise NativeEvaluationError(
                "this Symbolica build does not expose native complex profiling"
            )
        output, timing = profile(parameter_rows)
        if not isinstance(timing, tuple) or len(timing) != 3:
            raise NativeEvaluationError(
                "Symbolica returned an invalid native complex profile"
            )
        return output, tuple(float(value) for value in timing)

    def materialize(self) -> None:
        self._ensure_jit_compiled()

    def _ensure_jit_compiled(self) -> None:
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
        instance.settings = dict(manifest.get("settings", {}))
        instance.label = str(manifest.get("label", "jit_symbolica_evaluator"))
        instance.evaluator_state_path = _artifact_path_from_manifest(
            str(manifest["evaluator_state_path"]),
            artifact_dir,
        )
        instance._source_evaluator = Evaluator.load(
            instance.evaluator_state_path.read_bytes()
        )
        instance._progress_callback = None
        instance.build_timing = {}
        return instance

    def artifact_manifest(self, artifact_dir: Path) -> dict[str, Any]:
        _report_progress(
            self._progress_callback,
            stage="jit materialize",
            item=self.label,
        )
        materialize_started = time.perf_counter()
        self._ensure_jit_compiled()
        jit_compile_s = time.perf_counter() - materialize_started
        _report_progress(
            self._progress_callback,
            stage="jit materialized",
            item=f"{self.label} {jit_compile_s:.3f}s",
            increment=1,
        )
        evaluator_dir = _artifact_subdirectory(artifact_dir, "evaluators")
        path = evaluator_dir / f"{self.label}_{uuid.uuid4().hex}.evaluator.bin"
        save_started = time.perf_counter()
        path.write_bytes(self._source_evaluator.save())
        evaluator_save_s = time.perf_counter() - save_started
        self.evaluator_state_path = path
        build_timing = dict(self.build_timing)
        build_timing["jit_materialize_s"] = jit_compile_s
        build_timing["evaluator_save_s"] = evaluator_save_s
        build_timing["artifact_manifest_s"] = jit_compile_s + evaluator_save_s
        return {
            "kind": "jit-symbolica-evaluator",
            "backend": self.backend,
            "label": self.label,
            "input_len": self.input_len,
            "output_len": self.output_len,
            "jit_payload_layout": _jit_payload_layout(),
            "settings": self.settings,
            "evaluator_state_path": _artifact_path_for_manifest(path, artifact_dir),
            "build_timing": build_timing,
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
        self.settings = settings.to_json_dict()
        self.number_type = (
            "complex_4x"
            if settings.backend == "compiled-complex-4x"
            else "complex"
        )
        self._source_evaluator = evaluator
        self.build_timing: dict[str, float] = {}
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
            save_started = time.perf_counter()
            self.evaluator_state_path.write_bytes(save())
            self.build_timing["evaluator_save_s"] = time.perf_counter() - save_started
        else:
            self.evaluator_state_path = None
        compile_started = time.perf_counter()
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
        self.build_timing["cxx_compile_s"] = time.perf_counter() - compile_started

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
        instance.settings = dict(manifest.get("settings", {}))
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
        instance.build_timing = {}
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
            "settings": self.settings,
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
            "build_timing": dict(self.build_timing),
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


def _jit_payload_layout() -> str:
    return "native-complex-f64x2" if _aarch64_platform() else "scalar-complex-f64"


class _ChunkedSymbolicaEvaluator:
    def __init__(self, evaluators: tuple[Any, ...]) -> None:
        if not evaluators:
            raise NativeEvaluationError("chunked evaluator needs at least one chunk")
        self._evaluators = evaluators
        self.build_timing: dict[str, float] = {}

    def evaluate_complex(self, parameter_rows: Any) -> Any:
        return np.concatenate(self.evaluate_complex_chunks(parameter_rows), axis=1)

    def evaluate_complex_chunks(self, parameter_rows: Any) -> tuple[Any, ...]:
        prepared_rows = _complex128_parameter_rows(parameter_rows)
        return tuple(
            _evaluate_prepared_complex(evaluator, prepared_rows)
            for evaluator in self._evaluators
        )

    def evaluate_complex_profiled(
        self, parameter_rows: Any
    ) -> tuple[tuple[Any, ...], tuple[float, float, float]]:
        prepared_rows = _complex128_parameter_rows(parameter_rows)
        outputs: list[Any] = []
        profiles: list[tuple[float, float, float]] = []
        for evaluator in self._evaluators:
            output, profile = _evaluate_prepared_complex_profiled(
                evaluator, prepared_rows
            )
            outputs.append(output)
            profiles.append(profile)
        pack_times = sorted(profile[0] for profile in profiles)
        shared_pack_s = pack_times[len(pack_times) // 2]
        return tuple(outputs), (
            shared_pack_s,
            sum(profile[1] for profile in profiles),
            sum(profile[2] for profile in profiles),
        )

    def supports_complex_profiled(self) -> bool:
        return all(
            bool(getattr(evaluator, "supports_complex_profiled", lambda: False)())
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
        worker_counts = []
        for evaluator in self._evaluators:
            settings = getattr(evaluator, "settings", None)
            if isinstance(settings, Mapping):
                workers = settings.get("compiled_chunk_compile_workers")
                if isinstance(workers, int):
                    worker_counts.append(workers)
        workers = min(max(worker_counts, default=1), len(self._evaluators))
        if workers <= 1:
            chunks = [
                _symbolica_evaluator_artifact_manifest(evaluator, artifact_dir)
                for evaluator in self._evaluators
            ]
        else:
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = [
                    executor.submit(
                        _symbolica_evaluator_artifact_manifest,
                        evaluator,
                        artifact_dir,
                    )
                    for evaluator in self._evaluators
                ]
                chunks = [future.result() for future in futures]
        timing_totals = dict(self.build_timing)
        timing_totals["chunk_count"] = float(len(chunks))
        for chunk in chunks:
            chunk_timing = chunk.get("build_timing") if isinstance(chunk, dict) else None
            if not isinstance(chunk_timing, Mapping):
                continue
            for key, value in chunk_timing.items():
                if isinstance(value, (float, int)):
                    timing_totals[str(key)] = timing_totals.get(str(key), 0.0) + float(
                        value
                    )
        return {
            "kind": "chunked-symbolica-evaluator",
            "chunks": chunks,
            "build_timing": timing_totals,
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


def _evaluate_prepared_complex_profiled(
    evaluator: Any,
    parameter_rows: np.ndarray,
) -> tuple[Any, tuple[float, float, float]]:
    evaluate_profiled = getattr(
        evaluator, "_evaluate_complex_profiled_prepared", None
    )
    if not callable(evaluate_profiled):
        raise NativeEvaluationError(
            "evaluator does not expose native complex profiling"
        )
    return evaluate_profiled(parameter_rows)


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
            "Symbolica evaluator adapters"
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


__all__ = [
    "ComplexOutput",
    "ProgressCallback",
    "SymbolicaEvaluatorSettings",
    "_ChunkedSymbolicaEvaluator",
    "_artifact_path_for_manifest",
    "_artifact_path_from_manifest",
    "_compile_symbolica_outputs",
    "_evaluate_complex_outputs",
    "_load_symbolica_evaluator_artifact",
    "_resolve_compiled_preset",
    "_symbolica_evaluator_artifact_manifest",
]
