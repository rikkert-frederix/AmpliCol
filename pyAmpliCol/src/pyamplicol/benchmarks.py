from __future__ import annotations

import json
import math
import multiprocessing as mp
import time
import traceback
from pathlib import Path
from typing import Any, Literal, Sequence, cast

from .native import ExternalMomentum
from .processes import ProcessOptions
from .reference import AmplicolAdapter, TimingRow


def benchmark_z_gluon_modes(
    *,
    min_gluons: int = 0,
    max_gluons: int = 6,
    points: int = 3,
    timing: int | None = None,
    amplicol_root: str | Path,
    jobs: int = 8,
    timeout: float | None = None,
    options: ProcessOptions | None = None,
    numeric_timeout: float = 120.0,
    parametric_timeout: float = 180.0,
    parametric_max_gluons: int = 4,
    tensor_strategy: str = "interleaved",
    output: str | Path | None = None,
    batch_size: int = 16,
    merge_evaluators_strategy: bool = False,
    verbose_evaluator_build: bool = False,
    evaluator_build_kwargs: dict[str, Any] | None = None,
    include_python: bool = True,
    include_numeric_tn: bool = True,
    include_parametric_tn: bool = True,
    include_shared_dag: bool = True,
) -> dict[str, Any]:
    """Benchmark legacy AmpliCol and pyamplicol runtime modes on Z+n gluons.

    The numerical and parametric full tensor-network modes intentionally run in
    child processes so large exploratory contractions can be bounded by a
    timeout. Use the repository memory watchdog outside this command for RAM
    limits, as documented in AGENTS.md.
    """

    if min_gluons < 0:
        raise ValueError("min_gluons must be non-negative")
    if max_gluons < min_gluons:
        raise ValueError("max_gluons must be greater than or equal to min_gluons")
    if points < 1:
        raise ValueError("points must be positive")
    shared_build_kwargs = (
        {
            "batch_size": batch_size,
            "merge_evaluators_strategy": merge_evaluators_strategy,
            "verbose_evaluator_build": verbose_evaluator_build,
        }
        if evaluator_build_kwargs is None
        else dict(evaluator_build_kwargs)
    )
    batch_size = int(shared_build_kwargs.get("batch_size", batch_size))
    merge_evaluators_strategy = bool(
        shared_build_kwargs.get(
            "merge_evaluators_strategy",
            merge_evaluators_strategy,
        )
    )
    verbose_evaluator_build = bool(
        shared_build_kwargs.get(
            "verbose_evaluator_build",
            verbose_evaluator_build,
        )
    )
    if include_shared_dag and batch_size < 16:
        raise ValueError("batch_size must be at least 16 for timing benchmarks")

    adapter = AmplicolAdapter(amplicol_root, jobs=jobs, timeout=timeout)
    timing_sample = 1 if timing is None else timing
    payload: dict[str, Any] = {
        "family": "d d~ -> Z + n g",
        "min_gluons": min_gluons,
        "max_gluons": max_gluons,
        "points": points,
        "timing": timing_sample,
        "tensor_strategy": tensor_strategy,
        "batch_size": batch_size,
        "merge_evaluators_strategy": merge_evaluators_strategy,
        "evaluator_build_kwargs": shared_build_kwargs,
        "definitions": {
            "generation_s": {
                "legacy": (
                    "make cleanlib + make amplicol_generate + library=create + "
                    "make amplicol_generate_library; n=0 fixed-probe setup uses "
                    "cleanlib + amplicol_generate"
                ),
                "python": "NativeMatrixElementGenerator recursion/lowering metadata time",
                "numeric_tn": "factorized skeleton build plus simplify_color",
                "parametric_tn": (
                    "parametric TensorNetwork reduction plus Symbolica evaluator build"
                ),
                "shared_dag": (
                    "AmpliCol-style shared helicity-current table plus staged "
                    "Symbolica evaluator build"
                ),
            },
            "runtime_s_per_point": {
                "legacy": (
                    "AmpliCol detailed timing row 'amplitude evaluation' divided "
                    "by the requested timing sample count; benchmark timing uses a "
                    "second quiet direct-probe run with --amplicol_probe_quiet. "
                    "The error field is the printed-timer quantization floor."
                ),
                "pyamplicol": "mean wall time over the same AmpliCol probe points",
                "pyamplicol_evaluator_only": (
                    "time spent inside Symbolica evaluator calls only, warmed up "
                    "and divided by evaluated phase-space points"
                ),
            },
        },
        "rows": [],
    }

    output_path = None if output is None else Path(output)
    for gluon_count in range(min_gluons, max_gluons + 1):
        process = z_gluon_family_process(gluon_count)
        row = _benchmark_one_gluon_count(
            adapter=adapter,
            process=process,
            gluon_count=gluon_count,
            points=points,
            timing_sample=timing_sample,
            options=options,
            numeric_timeout=numeric_timeout,
            parametric_timeout=parametric_timeout,
            parametric_max_gluons=parametric_max_gluons,
            tensor_strategy=tensor_strategy,
            batch_size=batch_size,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            evaluator_build_kwargs=shared_build_kwargs,
            include_python=include_python,
            include_numeric_tn=include_numeric_tn,
            include_parametric_tn=include_parametric_tn,
            include_shared_dag=include_shared_dag,
        )
        payload["rows"].append(row)
        if output_path is not None:
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))

    payload["summary"] = summarize_mode_benchmark(payload["rows"])
    if output_path is not None:
        output_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def z_gluon_family_process(gluon_count: int) -> str:
    return "d d~ > z" + ("" if gluon_count == 0 else " " + " ".join(["g"] * gluon_count))


def profile_z_gluon_tensor_evaluator(
    process: str,
    *,
    sqrt_s: float = 1000.0,
    repetitions: int = 100,
    evaluator_repetitions: int | None = None,
    tensor_strategy: str = "interleaved",
) -> dict[str, Any]:
    """Profile the parametric tensor-network evaluator hot path."""

    if repetitions < 1:
        raise ValueError("repetitions must be positive")
    if evaluator_repetitions is not None and evaluator_repetitions < 1:
        raise ValueError("evaluator_repetitions must be positive")
    if tensor_strategy not in ("interleaved", "monolithic"):
        raise ValueError("tensor_strategy must be 'interleaved' or 'monolithic'")
    tensor_strategy_value = cast(
        Literal["interleaved", "monolithic"],
        tensor_strategy,
    )

    from .native import (
        LeadingColorZJetsNativeEvaluator,
        _ext_antiquark_weyl,
        _ext_gluon_cmplx,
        _ext_massive_vector,
        _ext_quark_weyl,
    )
    from .tensor_runtime import ZGluonTensorNetworkEvaluator, _negate_momentum

    native = LeadingColorZJetsNativeEvaluator()
    gluon_count = native.supported_z_gluon_count(process)
    if gluon_count is None or gluon_count < 1:
        raise ValueError(
            "tensor evaluator profiling currently supports q q~ -> Z plus gluons"
        )
    point = native.canonical_z_gluon_point(
        process,
        gluon_count=gluon_count,
        sqrt_s=sqrt_s,
    )

    generation_start = time.perf_counter()
    evaluator = ZGluonTensorNetworkEvaluator(process, strategy=tensor_strategy_value)
    generation_s = time.perf_counter() - generation_start

    anti_momentum = _negate_momentum(point[0].momentum)
    quark_momentum = _negate_momentum(point[1].momentum)
    evaluator._populate_parameters(
        point,
        chirality=1,
        quark_wf=_ext_quark_weyl(quark_momentum, -1, 1),
        anti_wf=_ext_antiquark_weyl(anti_momentum, 1, -1),
        gluon_wfs=tuple(
            _ext_gluon_cmplx(particle.momentum, -1)
            for particle in point[2 : 2 + gluon_count]
        ),
        z_wf=_ext_massive_vector(point[-1].momentum, -1, evaluator.model.mass(23)),
    )
    param_values = evaluator.bundle.param_builder.values
    pure_repetitions = (
        evaluator_repetitions
        if evaluator_repetitions is not None
        else max(100, repetitions)
    )
    pure_evaluator_s = _time_repeated(
        lambda: evaluator.bundle.evaluator.evaluate_complex([param_values]),
        pure_repetitions,
    )
    full_evaluate_s = _time_repeated(lambda: evaluator.evaluate(point), repetitions)
    result = evaluator.evaluate(point)
    helicity_calls = 2 * 3 * (2**gluon_count)
    return {
        "process": process,
        "sqrt_s": sqrt_s,
        "tensor_strategy": tensor_strategy,
        "gluon_count": gluon_count,
        "generation_s": generation_s,
        "parameter_count": len(param_values),
        "helicity_calls": helicity_calls,
        "pure_evaluator_call_s": pure_evaluator_s,
        "pure_evaluator_call_us": pure_evaluator_s * 1.0e6,
        "full_evaluate_s": full_evaluate_s,
        "full_evaluate_ms": full_evaluate_s * 1.0e3,
        "full_evaluate_per_helicity_us": full_evaluate_s / helicity_calls * 1.0e6,
        "repetitions": repetitions,
        "evaluator_repetitions": pure_repetitions,
        "matrix_element": result.matrix_element,
        "raw_helicity_sum": result.raw_helicity_sum,
    }


def profile_z_gluon_dag_evaluator(
    process: str,
    *,
    sqrt_s: float = 1000.0,
    points: int = 16,
    repetitions: int = 100,
    evaluator_build_kwargs: dict[str, Any] | None = None,
    save_evaluator_dir: str | Path | None = None,
) -> dict[str, Any]:
    """Profile the shared-current D-mode evaluator without legacy rebuilds."""

    if points < 1:
        raise ValueError("points must be positive")
    if repetitions < 1:
        raise ValueError("repetitions must be positive")
    build_kwargs = (
        {"batch_size": 16}
        if evaluator_build_kwargs is None
        else dict(evaluator_build_kwargs)
    )
    if (
        save_evaluator_dir is not None
        and build_kwargs.get("symbolica_compiled_output_dir") is None
        and build_kwargs.get("symbolica_load_evaluator_dir") is None
    ):
        build_kwargs["symbolica_compiled_output_dir"] = str(
            Path(save_evaluator_dir).expanduser() / "compiled"
        )
    batch_size = int(build_kwargs.get("batch_size", 16))
    if batch_size < 16:
        raise ValueError("batch_size must be at least 16 for timing profiles")

    from .dag_runtime import ZGluonDAGEvaluator
    from .native import LeadingColorZJetsNativeEvaluator

    native = LeadingColorZJetsNativeEvaluator()
    gluon_count = native.supported_z_gluon_count(process)
    if gluon_count is None or gluon_count < 1:
        raise ValueError(
            "D-mode profiling currently supports q q~ -> Z plus at least one gluon"
        )
    point_list = tuple(
        native.canonical_z_gluon_point(
            process,
            gluon_count=gluon_count,
            sqrt_s=sqrt_s * (1.0 + 0.037 * index),
        )
        for index in range(points)
    )

    generation_start = time.perf_counter()
    evaluator = ZGluonDAGEvaluator(process, **build_kwargs)
    saved_evaluator_manifest = None
    if save_evaluator_dir is not None:
        saved_evaluator_manifest = str(
            evaluator.save_evaluator_artifact(save_evaluator_dir)
        )
    generation_s = time.perf_counter() - generation_start

    _time_shared_dag_evaluations(evaluator, point_list)
    total_runtime_s = 0.0
    evaluator_runtime_s = 0.0
    runtime_samples_s_per_point: list[float] = []
    evaluator_samples_s_per_point: list[float] = []
    timing_samples: dict[str, list[float]] = {}
    matrix_elements: list[float] = []
    for repetition in range(repetitions):
        values, runtime_s, evaluator_s, timing = _evaluate_shared_dag_points_once(
            evaluator,
            point_list,
        )
        total_runtime_s += runtime_s
        evaluator_runtime_s += evaluator_s
        runtime_samples_s_per_point.append(runtime_s / points)
        evaluator_samples_s_per_point.append(evaluator_s / points)
        for key, value in timing.to_json_dict().items():
            timing_samples.setdefault(key, []).append(value / points)
        if repetition == 0:
            matrix_elements = values

    denominator = max(points * repetitions, 1)
    runtime_s_per_point = total_runtime_s / denominator
    evaluator_s_per_point = evaluator_runtime_s / denominator
    runtime_error_s = _standard_error(runtime_samples_s_per_point)
    evaluator_error_s = _standard_error(evaluator_samples_s_per_point)
    timing_s_per_point = {
        key: sum(samples) / len(samples)
        for key, samples in timing_samples.items()
    }
    timing_error_s_per_point = {
        key: _standard_error(samples)
        for key, samples in timing_samples.items()
    }
    return {
        "process": process,
        "sqrt_s": sqrt_s,
        "gluon_count": gluon_count,
        "points": points,
        "repetitions": repetitions,
        "batch_size": batch_size,
        "generation_s": generation_s,
        "saved_evaluator_manifest": saved_evaluator_manifest,
        "loaded_evaluator_dir": (
            None
            if build_kwargs.get("symbolica_load_evaluator_dir") is None
            else str(build_kwargs["symbolica_load_evaluator_dir"])
        ),
        "metadata": evaluator.metadata.to_json_dict(),
        "runtime_total_s": total_runtime_s,
        "runtime_s_per_point": runtime_s_per_point,
        "runtime_s_per_point_error": runtime_error_s,
        "runtime_us_per_point": runtime_s_per_point * 1.0e6,
        "runtime_us_per_point_error": (
            None if runtime_error_s is None else runtime_error_s * 1.0e6
        ),
        "runtime_evaluator_only_total_s": evaluator_runtime_s,
        "runtime_evaluator_only_s_per_point": evaluator_s_per_point,
        "runtime_evaluator_only_s_per_point_error": evaluator_error_s,
        "runtime_evaluator_only_us_per_point": (
            evaluator_s_per_point * 1.0e6
        ),
        "runtime_evaluator_only_us_per_point_error": (
            None if evaluator_error_s is None else evaluator_error_s * 1.0e6
        ),
        "runtime_samples_s_per_point": runtime_samples_s_per_point,
        "runtime_evaluator_only_samples_s_per_point": evaluator_samples_s_per_point,
        "runtime_breakdown_s_per_point": timing_s_per_point,
        "runtime_breakdown_us_per_point": {
            key: value * 1.0e6
            for key, value in timing_s_per_point.items()
        },
        "runtime_breakdown_us_per_point_error": {
            key: (None if value is None else value * 1.0e6)
            for key, value in timing_error_s_per_point.items()
        },
        "runtime_breakdown_samples_s_per_point": timing_samples,
        "matrix_elements": matrix_elements,
    }


def summarize_mode_benchmark(rows: Sequence[dict[str, Any]]) -> dict[str, Any]:
    mode_keys = ("legacy", "python", "numeric_tn", "parametric_tn", "shared_dag")
    summary: dict[str, Any] = {
        "rows": len(rows),
        "mode_success_counts": {},
        "max_relative_difference_to_legacy": None,
        "all_four_modes_match_rows": [],
        "all_four_modes_match_for_all_rows": False,
    }
    max_rel = 0.0
    saw_rel = False
    all_four_rows: list[int] = []
    for mode in mode_keys:
        summary["mode_success_counts"][mode] = sum(
            1 for row in rows if _mode_status(row, mode) == "ok"
        )
    for row in rows:
        available_ok = True
        for mode in mode_keys:
            if _mode_status(row, mode) != "ok":
                available_ok = False
                break
        if available_ok:
            all_four_rows.append(int(row["gluon_count"]))
        for mode in ("python", "numeric_tn", "parametric_tn", "shared_dag"):
            mode_result = row.get(mode)
            if not isinstance(mode_result, dict):
                continue
            rel = mode_result.get("max_relative_difference_to_legacy")
            if isinstance(rel, (float, int)):
                saw_rel = True
                max_rel = max(max_rel, float(rel))
    summary["max_relative_difference_to_legacy"] = max_rel if saw_rel else None
    summary["all_four_modes_match_rows"] = all_four_rows
    summary["all_four_modes_match_for_all_rows"] = len(all_four_rows) == len(rows)
    return summary


def format_mode_benchmark_table(payload: dict[str, Any]) -> str:
    lines = [
        "| n | process | A legacy gen/run | B python gen/run | "
        "C numeric TN gen/run | D scalar TN gen/run | D-new shared gen/run/eval-only | max rel diff |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["rows"]:
        values = [
            str(row["gluon_count"]),
            str(row["process"]),
            _format_mode_cell(row.get("legacy")),
            _format_mode_cell(row.get("python")),
            _format_mode_cell(row.get("numeric_tn")),
            _format_mode_cell(row.get("parametric_tn")),
            _format_shared_dag_cell(row.get("shared_dag")),
            _format_optional_float(_row_max_rel_to_legacy(row), precision=3),
        ]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _benchmark_one_gluon_count(
    *,
    adapter: AmplicolAdapter,
    process: str,
    gluon_count: int,
    points: int,
    timing_sample: int,
    options: ProcessOptions | None,
    numeric_timeout: float,
    parametric_timeout: float,
    parametric_max_gluons: int,
    tensor_strategy: str,
    batch_size: int,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    evaluator_build_kwargs: dict[str, Any],
    include_python: bool,
    include_numeric_tn: bool,
    include_parametric_tn: bool,
    include_shared_dag: bool,
) -> dict[str, Any]:
    if gluon_count == 0:
        build = adapter.prepare_direct_probe(process, options=options)
        reference_run = adapter.run_amplicol_fixed_probe(
            process,
            points=points,
            options=options,
            timing_sample=timing_sample,
        )
        timing_run = adapter.run_amplicol_fixed_probe(
            process,
            points=timing_sample,
            options=options,
            timing_sample=1,
            quiet=True,
        )
    else:
        build = adapter.prepare_library(process, options=options)
        reference_run = adapter.run_amplicol_probe(
            process,
            points=points,
            options=options,
            timing_sample=timing_sample,
        )
        timing_run = adapter.run_amplicol_probe(
            process,
            points=timing_sample,
            options=options,
            timing_sample=1,
            quiet=True,
        )
    probe_points = tuple(reference_run.probe_points)
    if not probe_points:
        raise RuntimeError(f"AmpliCol probe produced no points for {process}")

    reference_values = [float(point.matrix_element) for point in probe_points]
    particles = [_particle_payload(point.particles) for point in probe_points]
    row: dict[str, Any] = {
        "gluon_count": gluon_count,
        "process": process,
        "probe_point_count": len(probe_points),
        "reference_matrix_elements": reference_values,
        "legacy": {
            "status": "ok",
            "generation_s": build.total_command_time_s,
            "runtime_s_per_point": _legacy_runtime_per_point(
                timing_run.timing_rows,
                timing_sample,
            ),
            "runtime_s_per_point_error": _legacy_runtime_per_point_error(
                timing_run.timing_rows,
                timing_sample,
            ),
            "runtime_wall_s_per_point": (
                timing_run.total_command_time_s / timing_sample
            ),
            "commands_s": [command.elapsed_s for command in build.commands],
            "run_command_s": timing_run.total_command_time_s,
            "reference_collection_s": reference_run.total_command_time_s,
            "timing_probe_quiet": True,
        },
    }

    if include_python:
        row["python"] = _run_mode_child(
            "python",
            gluon_count,
            process,
            particles,
            reference_values,
            tensor_strategy,
            batch_size,
            merge_evaluators_strategy,
            verbose_evaluator_build,
            evaluator_build_kwargs,
            timeout_s=120.0,
        )
    if include_shared_dag:
        row["shared_dag"] = _run_mode_child(
            "shared_dag",
            gluon_count,
            process,
            particles,
            reference_values,
            tensor_strategy,
            batch_size,
            merge_evaluators_strategy,
            verbose_evaluator_build,
            evaluator_build_kwargs,
            timeout_s=parametric_timeout,
        )
    if include_numeric_tn:
        row["numeric_tn"] = _run_mode_child(
            "numeric_tn",
            gluon_count,
            process,
            particles,
            reference_values,
            tensor_strategy,
            batch_size,
            merge_evaluators_strategy,
            verbose_evaluator_build,
            evaluator_build_kwargs,
            timeout_s=numeric_timeout,
        )
    if not include_parametric_tn:
        pass
    elif gluon_count > parametric_max_gluons:
        row["parametric_tn"] = {
            "status": "guarded",
            "reason": (
                "current pyamplicol artifact guard allows full tensor-network "
                f"evaluator artifacts only through {parametric_max_gluons} "
                "final-state gluons"
            ),
        }
    else:
        row["parametric_tn"] = _run_mode_child(
            "parametric_tn",
            gluon_count,
            process,
            particles,
            reference_values,
            tensor_strategy,
            batch_size,
            merge_evaluators_strategy,
            verbose_evaluator_build,
            evaluator_build_kwargs,
            timeout_s=parametric_timeout,
        )

    row["all_available_max_pairwise_rel_diff"] = _available_max_pairwise_rel_diff(row)
    return row


def _worker(
    mode: str,
    gluon_count: int,
    process: str,
    particles_payload: list[list[dict[str, Any]]],
    reference_values: list[float],
    tensor_strategy: str,
    batch_size: int,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    evaluator_build_kwargs: dict[str, Any],
    queue: Any,
) -> None:
    try:
        from .evaluation import NativeRuntimeEvaluator
        from .matrix import NativeMatrixElementGenerator
        from .tensor_runtime import (
            TensorNetworkStrategy,
            ZGluonNumericTensorNetworkEvaluator,
            ZGluonTensorNetworkEvaluator,
        )

        particles = [_particles_from_payload(point) for point in particles_payload]
        strategy: TensorNetworkStrategy
        if tensor_strategy == "interleaved":
            strategy = "interleaved"
        elif tensor_strategy == "monolithic":
            strategy = "monolithic"
        else:
            raise ValueError(
                "tensor_strategy must be 'interleaved' or 'monolithic'"
            )
        result: dict[str, Any] = {"status": "ok"}
        if mode == "python":
            generation_start = time.perf_counter()
            generation = NativeMatrixElementGenerator().generate(
                process,
                write_cache_metadata=False,
            )
            result["generation_s"] = time.perf_counter() - generation_start
            result["generation_report_s"] = generation.generation_time_s
            setup_start = time.perf_counter()
            python_evaluator = NativeRuntimeEvaluator(
                process,
                runtime_backend="python",
            )
            result["setup_s"] = time.perf_counter() - setup_start
            result["metadata"] = python_evaluator.metadata.to_json_dict()
            values, runtime_s = _time_evaluations(
                lambda point: python_evaluator.evaluate(particles=point).matrix_element,
                particles,
            )
        elif mode == "shared_dag":
            if gluon_count == 0:
                queue.put(
                    {
                        "status": "unsupported",
                        "reason": "Z-only path uses the scalar zero-gluon evaluator",
                    }
                )
                return
            from .dag_runtime import ZGluonDAGEvaluator

            setup_start = time.perf_counter()
            shared_evaluator = ZGluonDAGEvaluator(
                process,
                **evaluator_build_kwargs,
            )
            result["generation_s"] = time.perf_counter() - setup_start
            result["generation_report_s"] = (
                shared_evaluator.symbolica_evaluator_build_time_s
                + shared_evaluator.current_table_build_time_s
            )
            result["metadata"] = shared_evaluator.metadata.to_json_dict()
            values, runtime_s, evaluator_runtime_s, timing = _time_shared_dag_evaluations(
                shared_evaluator,
                particles,
            )
            result["runtime_evaluator_only_s_per_point"] = (
                evaluator_runtime_s / max(len(particles), 1)
            )
            result["runtime_evaluator_only_total_s"] = evaluator_runtime_s
            result["runtime_breakdown_total_s"] = timing.to_json_dict()
            result["runtime_breakdown_s_per_point"] = {
                key: value / max(len(particles), 1)
                for key, value in timing.to_json_dict().items()
            }
        elif mode == "numeric_tn":
            if gluon_count == 0:
                queue.put(
                    {
                        "status": "unsupported",
                        "reason": "no nontrivial tensor network for Z-only process",
                    }
                )
                return
            setup_start = time.perf_counter()
            numeric_evaluator = ZGluonNumericTensorNetworkEvaluator(
                process,
                strategy=strategy,
            )
            result["generation_s"] = time.perf_counter() - setup_start
            result["generation_report_s"] = (
                numeric_evaluator.metadata.expression_build_time_s
            )
            result["metadata"] = numeric_evaluator.metadata.to_json_dict()
            values, runtime_s = _time_evaluations(
                lambda point: numeric_evaluator.evaluate(point).matrix_element,
                particles,
            )
        elif mode == "parametric_tn":
            if gluon_count == 0:
                queue.put(
                    {
                        "status": "unsupported",
                        "reason": (
                            "Z-only path has a scalar Symbolica evaluator, not the "
                            "full tensor-network evaluator"
                        ),
                    }
                )
                return
            setup_start = time.perf_counter()
            parametric_evaluator = ZGluonTensorNetworkEvaluator(
                process,
                strategy=strategy,
                build_helicity_filter=True,
            )
            result["generation_s"] = time.perf_counter() - setup_start
            metadata = parametric_evaluator.bundle.metadata or {}
            result["generation_report_s"] = metadata.get("total_build_s")
            result["tensor_network_reduction_s"] = metadata.get(
                "tensor_network_reduction_s"
            )
            result["symbolica_evaluator_build_s"] = metadata.get(
                "symbolica_evaluator_build_s"
            )
            result["metadata"] = parametric_evaluator.metadata.to_json_dict()
            result["bundle_metadata"] = metadata
            if parametric_evaluator.helicity_filter is not None:
                result["helicity_filter"] = (
                    parametric_evaluator.helicity_filter.to_json_dict()
                )
            values, runtime_s = _time_evaluations(
                lambda point: parametric_evaluator.evaluate(point).matrix_element,
                particles,
            )
        else:
            raise ValueError(f"unknown benchmark mode: {mode}")

        result["matrix_elements"] = values
        result["runtime_s_per_point"] = runtime_s / max(len(particles), 1)
        result["runtime_total_s"] = runtime_s
        result["max_relative_difference_to_legacy"] = max(
            _relative_difference(value, reference)
            for value, reference in zip(values, reference_values, strict=True)
        )
        queue.put(result)
    except BaseException as exc:  # noqa: BLE001 - benchmark failures are serialized.
        queue.put(
            {
                "status": "error",
                "error": repr(exc),
                "traceback": traceback.format_exc(limit=12),
            }
        )


def _run_mode_child(
    mode: str,
    gluon_count: int,
    process: str,
    particles_payload: list[list[dict[str, Any]]],
    reference_values: list[float],
    tensor_strategy: str,
    batch_size: int,
    merge_evaluators_strategy: bool,
    verbose_evaluator_build: bool,
    evaluator_build_kwargs: dict[str, Any],
    *,
    timeout_s: float,
) -> dict[str, Any]:
    ctx = mp.get_context("fork")
    queue = ctx.Queue()
    process_handle = ctx.Process(
        target=_worker,
        args=(
            mode,
            gluon_count,
            process,
            particles_payload,
            reference_values,
            tensor_strategy,
            batch_size,
            merge_evaluators_strategy,
            verbose_evaluator_build,
            evaluator_build_kwargs,
            queue,
        ),
    )
    start = time.perf_counter()
    process_handle.start()
    process_handle.join(timeout_s)
    elapsed_s = time.perf_counter() - start
    if process_handle.is_alive():
        process_handle.terminate()
        process_handle.join(10.0)
        if process_handle.is_alive():
            process_handle.kill()
            process_handle.join(5.0)
        return {
            "status": "timeout",
            "timeout_s": timeout_s,
            "child_elapsed_s": elapsed_s,
        }
    result: dict[str, Any]
    if queue.empty():
        result = {
            "status": "error",
            "error": (
                f"child exited with code {process_handle.exitcode} and no result"
            ),
        }
    else:
        result = queue.get()
    result = dict(result)
    result["child_elapsed_s"] = elapsed_s
    result["exitcode"] = process_handle.exitcode
    return result


def _time_evaluations(
    evaluate: Any,
    particles: Sequence[tuple[ExternalMomentum, ...]],
) -> tuple[list[float], float]:
    values: list[float] = []
    start = time.perf_counter()
    for point in particles:
        values.append(float(evaluate(point)))
    return values, time.perf_counter() - start


def _time_shared_dag_evaluations(
    evaluator: Any,
    particles: Sequence[tuple[ExternalMomentum, ...]],
) -> tuple[list[float], float, float, Any]:
    if not particles:
        from .dag_runtime import DAGEvaluationTiming

        return [], 0.0, 0.0, DAGEvaluationTiming()
    if hasattr(evaluator, "evaluate_matrix_elements_many"):
        evaluator.evaluate_matrix_elements_many(particles[:1])
    else:
        evaluator.evaluate_many(particles[:1])
    return _evaluate_shared_dag_points_once(evaluator, particles)


def _evaluate_shared_dag_points_once(
    evaluator: Any,
    particles: Sequence[tuple[ExternalMomentum, ...]],
) -> tuple[list[float], float, float, Any]:
    from .dag_runtime import DAGEvaluationTiming

    values: list[float] = []
    evaluator_runtime_s = 0.0
    timing = DAGEvaluationTiming()
    start = time.perf_counter()
    for offset in range(0, len(particles), evaluator.batch_size):
        batch = particles[offset : offset + evaluator.batch_size]
        if hasattr(evaluator, "evaluate_matrix_elements_many"):
            values.extend(float(value) for value in evaluator.evaluate_matrix_elements_many(batch))
        else:
            results = evaluator.evaluate_many(batch)
            values.extend(float(result.matrix_element) for result in results)
        evaluator_runtime_s += evaluator.compiled.last_evaluator_time_s
        timing += evaluator.last_runtime_timing
    return values, time.perf_counter() - start, evaluator_runtime_s, timing


def _time_repeated(func: Any, repetitions: int) -> float:
    func()
    start = time.perf_counter()
    for _ in range(repetitions):
        func()
    return (time.perf_counter() - start) / repetitions


def _particle_payload(particles: Sequence[ExternalMomentum]) -> list[dict[str, Any]]:
    return [
        {"pdg": int(particle.pdg), "momentum": [float(x) for x in particle.momentum]}
        for particle in particles
    ]


def _particles_from_payload(
    payload: Sequence[dict[str, Any]],
) -> tuple[ExternalMomentum, ...]:
    particles: list[ExternalMomentum] = []
    for item in payload:
        components = tuple(float(component) for component in item["momentum"])
        if len(components) != 4:
            raise ValueError("external momentum payload must have four components")
        particles.append(
            ExternalMomentum(
                int(item["pdg"]),
                (components[0], components[1], components[2], components[3]),
            )
        )
    return tuple(particles)


def _legacy_runtime_per_point(
    rows: Sequence[TimingRow],
    timing_sample: int,
) -> float | None:
    for row in rows:
        if row.label == "amplitude evaluation":
            return row.seconds / timing_sample
    return None


def _legacy_runtime_per_point_error(
    rows: Sequence[TimingRow],
    timing_sample: int,
) -> float | None:
    for row in rows:
        if row.label == "amplitude evaluation":
            # Legacy timing summaries are printed with coarse decimal precision.
            # Treat half a millisecond on the printed total as the quantization
            # uncertainty floor for one timing run.
            return 0.5e-3 / timing_sample
    return None


def _standard_error(samples: Sequence[float]) -> float | None:
    if len(samples) < 2:
        return None
    mean = sum(samples) / len(samples)
    variance = sum((sample - mean) ** 2 for sample in samples) / (len(samples) - 1)
    return math.sqrt(variance / len(samples))


def _available_max_pairwise_rel_diff(row: dict[str, Any]) -> float | None:
    values: list[float] = []
    legacy_values = row.get("reference_matrix_elements")
    if isinstance(legacy_values, list) and legacy_values:
        values.append(float(legacy_values[0]))
    for mode in ("python", "numeric_tn", "parametric_tn", "shared_dag"):
        result = row.get(mode)
        if not isinstance(result, dict) or result.get("status") != "ok":
            continue
        mode_values = result.get("matrix_elements")
        if isinstance(mode_values, list) and mode_values:
            values.append(float(mode_values[0]))
    if len(values) < 2:
        return None
    return max(
        _relative_difference(left, right)
        for index, left in enumerate(values)
        for right in values[index + 1 :]
    )


def _row_max_rel_to_legacy(row: dict[str, Any]) -> float | None:
    max_rel: float | None = None
    for mode in ("python", "numeric_tn", "parametric_tn", "shared_dag"):
        result = row.get(mode)
        if not isinstance(result, dict):
            continue
        rel = result.get("max_relative_difference_to_legacy")
        if not isinstance(rel, (float, int)):
            continue
        value = float(rel)
        max_rel = value if max_rel is None else max(max_rel, value)
    return max_rel


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _mode_status(row: dict[str, Any], mode: str) -> str | None:
    result = row.get(mode)
    if isinstance(result, dict):
        status = result.get("status")
        if isinstance(status, str):
            return status
    return None


def _format_mode_cell(result: object) -> str:
    if not isinstance(result, dict):
        return "n/a"
    status = result.get("status")
    if status != "ok":
        reason = result.get("reason") or result.get("error")
        if isinstance(reason, str) and reason:
            return f"{status} ({reason})"
        return str(status)
    return (
        f"{_format_optional_float(result.get('generation_s'))} / "
        f"{_format_timing_with_error(result, 'runtime_s_per_point')}"
    )


def _format_shared_dag_cell(result: object) -> str:
    if not isinstance(result, dict):
        return "n/a"
    status = result.get("status")
    if status != "ok":
        reason = result.get("reason") or result.get("error")
        if isinstance(reason, str) and reason:
            return f"{status} ({reason})"
        return str(status)
    return (
        f"{_format_optional_float(result.get('generation_s'))} / "
        f"{_format_timing_with_error(result, 'runtime_s_per_point')} / "
        f"{_format_timing_with_error(result, 'runtime_evaluator_only_s_per_point')}"
    )


def _format_optional_float(value: object, *, precision: int = 4) -> str:
    if not isinstance(value, (float, int)):
        return "n/a"
    return f"{float(value):.{precision}g}"


def _format_timing_with_error(
    result: dict[str, Any],
    key: str,
    *,
    precision: int = 4,
) -> str:
    value = result.get(key)
    if not isinstance(value, (float, int)):
        return "n/a"
    error = result.get(f"{key}_error")
    formatted_value = _format_optional_float(value, precision=precision)
    if not isinstance(error, (float, int)):
        return formatted_value
    return f"{formatted_value} +/- {_format_optional_float(error, precision=2)}"


__all__ = [
    "benchmark_z_gluon_modes",
    "format_mode_benchmark_table",
    "profile_z_gluon_dag_evaluator",
    "profile_z_gluon_tensor_evaluator",
    "summarize_mode_benchmark",
    "z_gluon_family_process",
]
