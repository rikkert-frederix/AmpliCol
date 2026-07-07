#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import json
import math
import os
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path
from typing import Any, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
PYAMPLICOL_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PYAMPLICOL_ROOT / "src"
DEFAULT_JSON = PYAMPLICOL_ROOT / ".benchmarks" / "performance_comparison.json"
DOCS_DIR = PYAMPLICOL_ROOT / "docs"
DEFAULT_MARKDOWN = DOCS_DIR / "performance_comparison.md"
LEGACY_MARKDOWN = DOCS_DIR / "performance_compairson.md"


MODE_SPECS: dict[str, dict[str, Any]] = {
    "X": {
        "label": "pyamplicol C++ runtime preset",
        "cli": (
            "--runtime-backend dag "
            "--symbolica-evaluator-backend compiled-complex "
            "--symbolica-compiled-preset runtime"
        ),
        "kwargs": {
            "symbolica_evaluator_backend": "compiled-complex",
            "symbolica_compiled_preset": "runtime",
        },
    },
    "Y": {
        "label": "pyamplicol assembly generation preset",
        "cli": (
            "--runtime-backend dag "
            "--symbolica-evaluator-backend compiled-complex "
            "--symbolica-compiled-preset generation"
        ),
        "kwargs": {
            "symbolica_evaluator_backend": "compiled-complex",
            "symbolica_compiled_preset": "generation",
        },
    },
    "Z": {
        "label": "pyamplicol JIT",
        "cli": "--runtime-backend dag --symbolica-evaluator-backend jit",
        "kwargs": {
            "symbolica_evaluator_backend": "jit",
        },
    },
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Collect the A/X/Y/Z pyamplicol performance comparison table."
    )
    parser.add_argument("--min-gluons", type=int, default=1)
    parser.add_argument("--max-gluons", type=int, default=8)
    parser.add_argument(
        "--extra-z-gluons",
        type=int,
        default=9,
        help="Additional gluon multiplicity for A and Z only; set to -1 to disable.",
    )
    parser.add_argument("--points", type=int, default=16)
    parser.add_argument("--repetitions", type=int, default=20)
    parser.add_argument("--legacy-timing", type=int, default=10_000)
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--mode-timeout", type=float, default=900.0)
    parser.add_argument("--amplicol-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-md", type=Path, default=DEFAULT_MARKDOWN)
    parser.add_argument("--keep-going", action="store_true", default=True)
    parser.add_argument("--mode-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--worker-input", type=Path, help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    _ensure_import_path()
    if args.mode_worker:
        return _worker_main(args)
    return _collector_main(args)


def _ensure_import_path() -> None:
    src = str(SRC_DIR)
    if src not in sys.path:
        sys.path.insert(0, src)


def _collector_main(args: argparse.Namespace) -> int:
    if args.min_gluons < 1:
        raise ValueError("--min-gluons must be at least 1 for this comparison")
    if args.max_gluons < args.min_gluons:
        raise ValueError("--max-gluons must be >= --min-gluons")
    if args.points < 1:
        raise ValueError("--points must be positive")
    if args.repetitions < 1:
        raise ValueError("--repetitions must be positive")
    if args.legacy_timing < 1:
        raise ValueError("--legacy-timing must be positive")

    from pyamplicol.reference import AmplicolAdapter

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)

    gluon_counts = list(range(args.min_gluons, args.max_gluons + 1))
    if args.extra_z_gluons >= 0 and args.extra_z_gluons not in gluon_counts:
        gluon_counts.append(args.extra_z_gluons)

    payload: dict[str, Any] = {
        "generated_at_epoch_s": time.time(),
        "family": "d d~ -> Z + n g",
        "points": args.points,
        "repetitions": args.repetitions,
        "legacy_timing": args.legacy_timing,
        "batch_size": args.batch_size,
        "mode_specs": {
            key: {"label": spec["label"], "cli": spec["cli"]}
            for key, spec in MODE_SPECS.items()
        },
        "versions": _dependency_versions(),
        "rows": [],
    }

    adapter = AmplicolAdapter(
        args.amplicol_root,
        jobs=args.jobs,
        timeout=args.timeout,
    )

    for gluon_count in gluon_counts:
        process = _z_gluon_process(gluon_count)
        print(f"== n={gluon_count}: legacy AmpliCol ==", flush=True)
        row: dict[str, Any] = {
            "gluon_count": gluon_count,
            "process": process,
            "legacy": _collect_legacy(
                adapter,
                process,
                points=args.points,
                timing=args.legacy_timing,
            ),
            "modes": {},
        }
        particles = row["legacy"]["reference_particles"]
        reference_values = row["legacy"]["reference_matrix_elements"]
        mode_keys = ("Z",) if gluon_count > args.max_gluons else ("X", "Y", "Z")
        for mode_key in mode_keys:
            print(f"== n={gluon_count}: mode {mode_key} ==", flush=True)
            row["modes"][mode_key] = _run_mode_worker(
                mode_key,
                process,
                particles,
                reference_values,
                batch_size=args.batch_size,
                repetitions=args.repetitions,
                timeout=args.mode_timeout,
            )
            _write_outputs(payload, row, args.output_json, args.output_md)
        payload["rows"].append(row)
        _write_payload(payload, args.output_json, args.output_md)

    return 0


def _collect_legacy(
    adapter: Any,
    process: str,
    *,
    points: int,
    timing: int,
) -> dict[str, Any]:
    from pyamplicol.benchmarks import (
        _legacy_runtime_per_point,
        _legacy_runtime_per_point_error,
    )

    build = adapter.prepare_library(process)
    reference_run = adapter.run_amplicol_probe(
        process,
        points=points,
        process_file=build.process_file,
        timing_sample=1,
        use_library=True,
    )
    timing_run = adapter.run_library_use(
        process,
        nevents=timing,
        seed=101,
        process_file=build.process_file,
        timing_sample=1,
    )
    runtime_s = _legacy_runtime_per_point(timing_run.timing_rows, timing)
    runtime_error_s = _legacy_runtime_per_point_error(timing_run.timing_rows, timing)
    return {
        "status": "ok",
        "generation_s": build.total_command_time_s,
        "runtime_s_per_point": runtime_s,
        "runtime_s_per_point_error": runtime_error_s,
        "reference_collection_s": reference_run.total_command_time_s,
        "timing_probe_s": timing_run.total_command_time_s,
        "commands_s": [command.elapsed_s for command in build.commands],
        "reference_matrix_elements": [
            float(point.matrix_element) for point in reference_run.probe_points
        ],
        "reference_particles": [
            _particle_payload(point.particles) for point in reference_run.probe_points
        ],
        "timing_rows": [
            {"label": row.label, "seconds": row.seconds, "note": row.note}
            for row in timing_run.timing_rows
        ],
        "max_relative_difference": 0.0,
    }


def _run_mode_worker(
    mode_key: str,
    process: str,
    particles: list[list[dict[str, Any]]],
    reference_values: list[float],
    *,
    batch_size: int,
    repetitions: int,
    timeout: float,
) -> dict[str, Any]:
    worker_payload = {
        "mode_key": mode_key,
        "process": process,
        "particles": particles,
        "reference_values": reference_values,
        "batch_size": batch_size,
        "repetitions": repetitions,
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
        json.dump(worker_payload, handle)
        worker_input = Path(handle.name)
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC_DIR) + os.pathsep + env.get("PYTHONPATH", "")
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--mode-worker",
        "--worker-input",
        str(worker_input),
    ]
    start = time.perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=REPO_ROOT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = time.perf_counter() - start
        worker_input.unlink(missing_ok=True)
        return {
            "status": "timeout",
            "mode_key": mode_key,
            "label": MODE_SPECS[mode_key]["label"],
            "elapsed_s": elapsed,
            "timeout_s": timeout,
            "stdout": exc.stdout,
            "stderr": exc.stderr,
        }
    finally:
        worker_input.unlink(missing_ok=True)

    elapsed = time.perf_counter() - start
    if completed.returncode != 0:
        return {
            "status": "error",
            "mode_key": mode_key,
            "label": MODE_SPECS[mode_key]["label"],
            "elapsed_s": elapsed,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
        }
    try:
        result = json.loads(completed.stdout)
    except json.JSONDecodeError:
        return {
            "status": "error",
            "mode_key": mode_key,
            "label": MODE_SPECS[mode_key]["label"],
            "elapsed_s": elapsed,
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "error": "worker did not emit JSON",
        }
    result["child_elapsed_s"] = elapsed
    return result


def _worker_main(args: argparse.Namespace) -> int:
    if args.worker_input is None:
        raise ValueError("--worker-input is required in worker mode")
    data = json.loads(args.worker_input.read_text())
    mode_key = str(data["mode_key"])
    process = str(data["process"])
    particles = [_particles_from_payload(point) for point in data["particles"]]
    reference_values = [float(value) for value in data["reference_values"]]
    batch_size = int(data["batch_size"])
    repetitions = int(data["repetitions"])
    try:
        result = _collect_py_mode(
            mode_key,
            process,
            particles,
            reference_values,
            batch_size=batch_size,
            repetitions=repetitions,
        )
    except BaseException as exc:  # noqa: BLE001 - serialized for parent process.
        print(
            json.dumps(
                {
                    "status": "error",
                    "mode_key": mode_key,
                    "label": MODE_SPECS[mode_key]["label"],
                    "error": repr(exc),
                    "traceback": traceback.format_exc(limit=16),
                },
                sort_keys=True,
            )
        )
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


def _collect_py_mode(
    mode_key: str,
    process: str,
    particles: Sequence[tuple[Any, ...]],
    reference_values: Sequence[float],
    *,
    batch_size: int,
    repetitions: int,
) -> dict[str, Any]:
    from pyamplicol.benchmarks import _standard_error
    from pyamplicol.dag_runtime import ZGluonDAGEvaluator

    kwargs = {
        "batch_size": batch_size,
        **MODE_SPECS[mode_key]["kwargs"],
    }
    build_start = time.perf_counter()
    evaluator = ZGluonDAGEvaluator(process, **kwargs)
    generation_s = time.perf_counter() - build_start

    values = [float(value) for value in evaluator.evaluate_matrix_elements_many(particles)]
    max_rel = max(
        _relative_difference(value, reference)
        for value, reference in zip(values, reference_values, strict=True)
    )

    evaluator.evaluate_matrix_elements_many(particles[:batch_size])
    runtime_samples: list[float] = []
    evaluator_samples: list[float] = []
    breakdown_samples: dict[str, list[float]] = {}
    for _ in range(repetitions):
        start = time.perf_counter()
        evaluator.evaluate_matrix_elements_many(particles)
        runtime_s = time.perf_counter() - start
        evaluator_s = evaluator.compiled.last_evaluator_time_s
        denominator = max(len(particles), 1)
        runtime_samples.append(runtime_s / denominator)
        evaluator_samples.append(evaluator_s / denominator)
        for key, value in evaluator.last_runtime_timing.to_json_dict().items():
            breakdown_samples.setdefault(key, []).append(value / denominator)

    runtime_s_per_point = sum(runtime_samples) / len(runtime_samples)
    evaluator_s_per_point = sum(evaluator_samples) / len(evaluator_samples)
    result = {
        "status": "ok",
        "mode_key": mode_key,
        "label": MODE_SPECS[mode_key]["label"],
        "generation_s": generation_s,
        "runtime_s_per_point": runtime_s_per_point,
        "runtime_s_per_point_error": _standard_error(runtime_samples),
        "runtime_evaluator_only_s_per_point": evaluator_s_per_point,
        "runtime_evaluator_only_s_per_point_error": _standard_error(evaluator_samples),
        "runtime_breakdown_s_per_point": {
            key: sum(values_for_key) / len(values_for_key)
            for key, values_for_key in breakdown_samples.items()
        },
        "runtime_breakdown_s_per_point_error": {
            key: _standard_error(values_for_key)
            for key, values_for_key in breakdown_samples.items()
        },
        "max_relative_difference": max_rel,
        "matrix_elements": values,
        "metadata": evaluator.metadata.to_json_dict(),
    }
    del evaluator
    gc.collect()
    return result


def _write_outputs(
    payload: dict[str, Any],
    pending_row: dict[str, Any],
    output_json: Path,
    output_md: Path,
) -> None:
    rows = list(payload["rows"])
    if pending_row not in rows:
        rows.append(pending_row)
    tmp_payload = {**payload, "rows": rows}
    output_json.write_text(json.dumps(tmp_payload, indent=2, sort_keys=True) + "\n")
    _write_markdown(output_md, _render_markdown(tmp_payload))


def _write_payload(payload: dict[str, Any], output_json: Path, output_md: Path) -> None:
    output_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    _write_markdown(output_md, _render_markdown(payload))


def _write_markdown(output_md: Path, text: str) -> None:
    output_md.write_text(text, encoding="utf-8")
    resolved_output = output_md.resolve()
    default_output = DEFAULT_MARKDOWN.resolve()
    legacy_output = LEGACY_MARKDOWN.resolve()
    if resolved_output == default_output and default_output != legacy_output:
        LEGACY_MARKDOWN.write_text(text, encoding="utf-8")
    elif resolved_output == legacy_output and default_output != legacy_output:
        DEFAULT_MARKDOWN.write_text(text, encoding="utf-8")


def _render_markdown(payload: dict[str, Any]) -> str:
    lines = [
        "# pyamplicol Performance Comparison",
        "",
        "Process family: `d d~ -> Z + n g`.",
        "",
        "All heavy runs for this table were launched under `pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30`.",
        f"Validation points per multiplicity: `{payload['points']}`.",
        f"pyamplicol timing repetitions: `{payload['repetitions']}`.",
        f"Legacy timing sample: `{payload['legacy_timing']}`.",
        f"pyamplicol batch size: `{payload['batch_size']}`.",
        "",
        "## Modes",
        "",
        "- A: legacy AmpliCol generated Fortran library, timed with the `--library=use` integration path. Probe reference values are also collected with `--library=use`.",
        "- X: pyamplicol shared-helicity-current D-mode with Symbolica `compiled-complex` generic C++ output and the `runtime` preset.",
        "- Y: pyamplicol shared-helicity-current D-mode with Symbolica `compiled-complex` emitted assembly output and the `generation` preset.",
        "- Z: pyamplicol shared-helicity-current D-mode with Symbolica JIT compilation. This is the generation-time-oriented mode.",
        "",
        "For pyamplicol rows, `wall` is full `evaluate_matrix_elements_many` wall time per phase-space point and `eval` is time spent inside Symbolica evaluator calls only. For legacy AmpliCol, `eval` is the printed `amplitude evaluation` timing per point.",
        "",
        "## Table",
        "",
        "| n | mode | status | gen [s] | wall [us/pt] | eval [us/pt] | max rel diff | notes |",
        "|---:|---|---|---:|---:|---:|---:|---|",
    ]
    for row in payload["rows"]:
        n = int(row["gluon_count"])
        legacy = row["legacy"]
        lines.append(
            "| "
            + " | ".join(
                [
                    str(n),
                    "A",
                    str(legacy.get("status", "unknown")),
                    _format_seconds(legacy.get("generation_s")),
                    "n/a",
                    _format_timing_us(
                        legacy.get("runtime_s_per_point"),
                        legacy.get("runtime_s_per_point_error"),
                    ),
                    _format_float(legacy.get("max_relative_difference")),
                    "legacy AmpliCol",
                ]
            )
            + " |"
        )
        for mode_key in ("X", "Y", "Z"):
            mode = row.get("modes", {}).get(mode_key)
            if mode is None:
                continue
            lines.append(
                "| "
                + " | ".join(
                    [
                        str(n),
                        mode_key,
                        str(mode.get("status", "unknown")),
                        _format_seconds(mode.get("generation_s")),
                        _format_timing_us(
                            mode.get("runtime_s_per_point"),
                            mode.get("runtime_s_per_point_error"),
                        ),
                        _format_timing_us(
                            mode.get("runtime_evaluator_only_s_per_point"),
                            mode.get("runtime_evaluator_only_s_per_point_error"),
                        ),
                        _format_float(mode.get("max_relative_difference")),
                        _mode_note(mode),
                    ]
                )
                + " |"
            )
    lines.extend(
        [
            "",
            "## Reproduction Commands",
            "",
            "Set up the shell once:",
            "",
            "```bash",
            "cd /Users/vjhirsch/HEP_programs/AmpliCol",
            "export PYTHONPATH=pyAmpliCol/src",
            "PY=pyAmpliCol/dependencies/.venv/bin/python",
            "WATCH='pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30 --'",
            "```",
            "",
            "The full table can be regenerated with:",
            "",
            "```bash",
            "$WATCH $PY pyAmpliCol/scripts/collect_performance_comparison.py \\",
            "  --min-gluons 1 --max-gluons 8 --extra-z-gluons 9 \\",
            f"  --points {payload['points']} --repetitions {payload['repetitions']} \\",
            f"  --legacy-timing {payload['legacy_timing']} --batch-size {payload['batch_size']}",
            "```",
            "",
            "For an individual process, define for example:",
            "",
            "```bash",
            "PROCESS='d d~ > z g g g'",
            "```",
            "",
            "Legacy AmpliCol validation/timing:",
            "",
            "```bash",
            "$WATCH $PY -m pyamplicol compare-amplicol \"$PROCESS\" \\",
            f"  --amplicol-probe --points {payload['points']} --timing {payload['legacy_timing']} --json",
            "```",
            "",
            "Mode X profile:",
            "",
            "```bash",
            "$WATCH $PY -m pyamplicol profile-dag-evaluator \"$PROCESS\" \\",
            f"  --points {payload['points']} --repetitions {payload['repetitions']} --batch-size {payload['batch_size']} \\",
            "  --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset runtime --json",
            "```",
            "",
            "Mode Y profile:",
            "",
            "```bash",
            "$WATCH $PY -m pyamplicol profile-dag-evaluator \"$PROCESS\" \\",
            f"  --points {payload['points']} --repetitions {payload['repetitions']} --batch-size {payload['batch_size']} \\",
            "  --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset generation --json",
            "```",
            "",
            "Mode Z profile:",
            "",
            "```bash",
            "$WATCH $PY -m pyamplicol profile-dag-evaluator \"$PROCESS\" \\",
            f"  --points {payload['points']} --repetitions {payload['repetitions']} --batch-size {payload['batch_size']} \\",
            "  --symbolica-evaluator-backend jit --json",
            "```",
            "",
            "The mode-specific validation command is the same as the legacy command, with `--runtime-backend dag` and the mode-specific Symbolica options added.",
            "",
            "## Dependency Versions",
            "",
            "```json",
            json.dumps(payload.get("versions", {}), indent=2, sort_keys=True),
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def _mode_note(mode: dict[str, Any]) -> str:
    if mode.get("status") == "ok":
        settings = (
            mode.get("metadata", {})
            .get("symbolica_evaluator_settings", {})
        )
        backend = settings.get("backend")
        preset = settings.get("compiled_preset")
        inline_asm = settings.get("compiled_inline_asm")
        opt = settings.get("compiled_optimization_level")
        chunk = settings.get("compiled_output_chunk_size")
        if backend == "jit":
            return "jit"
        return f"{backend}, preset={preset}, asm={inline_asm}, O{opt}, chunk={chunk}"
    stderr = str(mode.get("stderr") or "").strip().splitlines()
    if stderr:
        return _escape_table(stderr[0][:160])
    error = mode.get("error")
    if error:
        return _escape_table(str(error)[:160])
    return _escape_table(str(mode.get("returncode", "")))


def _format_seconds(value: Any) -> str:
    if not isinstance(value, (float, int)) or not math.isfinite(float(value)):
        return "n/a"
    return f"{float(value):.3f}"


def _format_timing_us(value: Any, error: Any = None) -> str:
    if not isinstance(value, (float, int)) or not math.isfinite(float(value)):
        return "n/a"
    us = float(value) * 1.0e6
    if isinstance(error, (float, int)) and math.isfinite(float(error)):
        return f"{us:.3f} +/- {float(error) * 1.0e6:.3f}"
    return f"{us:.3f}"


def _format_float(value: Any) -> str:
    if not isinstance(value, (float, int)) or not math.isfinite(float(value)):
        return "n/a"
    return f"{float(value):.3e}"


def _escape_table(value: str) -> str:
    return value.replace("|", "\\|").replace("\n", " ")


def _z_gluon_process(gluon_count: int) -> str:
    return "d d~ > z" + ("" if gluon_count == 0 else " " + " ".join(["g"] * gluon_count))


def _particle_payload(particles: Sequence[Any]) -> list[dict[str, Any]]:
    return [
        {"pdg": int(particle.pdg), "momentum": [float(x) for x in particle.momentum]}
        for particle in particles
    ]


def _particles_from_payload(payload: Sequence[dict[str, Any]]) -> tuple[Any, ...]:
    from pyamplicol.native import ExternalMomentum

    particles = []
    for particle in payload:
        momentum = tuple(float(x) for x in particle["momentum"])
        particles.append(ExternalMomentum(int(particle["pdg"]), momentum))
    return tuple(particles)


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _dependency_versions() -> dict[str, Any]:
    versions: dict[str, Any] = {}
    try:
        import symbolica

        versions["symbolica_version"] = getattr(symbolica, "__version__", None)
        versions["symbolica_local_versions"] = getattr(symbolica, "LOCAL_VERSIONS", None)
    except Exception as exc:  # noqa: BLE001
        versions["symbolica_error"] = repr(exc)
    manifest = PYAMPLICOL_ROOT / "dependencies" / "install_manifest.json"
    if manifest.exists():
        try:
            versions["install_manifest"] = json.loads(manifest.read_text())
        except json.JSONDecodeError as exc:
            versions["install_manifest_error"] = str(exc)
    return versions


if __name__ == "__main__":
    raise SystemExit(main())
