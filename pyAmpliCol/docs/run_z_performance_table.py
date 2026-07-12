#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Sequence

DOCS_DIR = Path(__file__).resolve().parent
PROJECT_DIR = DOCS_DIR.parent
REPO_ROOT = PROJECT_DIR.parent
SRC_DIR = PROJECT_DIR / "src"
SCRIPTS_DIR = PROJECT_DIR / "scripts"
PYTHON = PROJECT_DIR / "dependencies" / ".venv" / "bin" / "python"
Z_DATA = DOCS_DIR / "z_performance_data.json"
LC_MATRIX_DATA = DOCS_DIR / "result_matrix_data.json"
OUTPUT_ROOT = DOCS_DIR / ".z_performance_outputs"
TABLE_SCRIPT = DOCS_DIR / "z_performance_table.py"

DEFAULT_BATCH_SIZE = 64
DEFAULT_CHUNK_SIZE = 128
DEFAULT_ALL_FLOW_BATCH_SIZE = 64
DEFAULT_ALL_FLOW_CHUNK_SIZE = 8192
DEFAULT_SYMBOLICA_ITERATIONS = 10
DEFAULT_TARGET_RUNTIME = 10.0
DEFAULT_N_CORES = 5
DEFAULT_CPP_TIME_LIMIT = 900.0
DEFAULT_MEMORY_LIMIT_GB = 30.0
HEARTBEAT_S = 30.0
MEMORY_POLL_S = 1.0

MODE_KEYS = ("amplicol", "jit_o1", "asm", "cpp_o3", "jit_o3")

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from memory_watch_support import ProcessTreeMemoryMonitor  # noqa: E402


class MemoryLimitExceeded(RuntimeError):
    def __init__(
        self,
        *,
        limit_gb: float,
        peak_memory_bytes: int,
        memory_metric: str,
        command: Sequence[str],
    ):
        self.limit_gb = float(limit_gb)
        self.peak_memory_bytes = int(peak_memory_bytes)
        self.peak_memory_gb = self.peak_memory_bytes / 1024**3
        self.memory_metric = str(memory_metric)
        self.command = tuple(str(item) for item in command)
        super().__init__(
            f"memory limit {self.limit_gb:g} GB exceeded "
            f"({self.memory_metric.replace('_', ' ')} "
            f"{self.peak_memory_gb:.3f} GiB): {' '.join(self.command)}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Populate docs/z_performance_data.json and z_performance_table.tex "
            "for the d d~ > Z + (n-1)*g dedicated comparison table."
        )
    )
    parser.add_argument(
        "--n",
        nargs="*",
        default=("1-9",),
        help="Final-state multiplicities to run, e.g. 1-6 or 7 8 9.",
    )
    parser.add_argument(
        "--modes",
        nargs="*",
        default=("all",),
        help=(
            "Rows to populate: amplicol, jit_o1, asm, cpp_o3, jit_o3, or all. "
            "Comma-separated lists are accepted."
        ),
    )
    parser.add_argument("--n-cores", type=int, default=DEFAULT_N_CORES)
    parser.add_argument("--target-runtime", type=float, default=DEFAULT_TARGET_RUNTIME)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--output-chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument(
        "--all-flow-batch-size",
        type=int,
        default=DEFAULT_ALL_FLOW_BATCH_SIZE,
        help="Batch size used by time-process for the all-flow fixed-helicity block.",
    )
    parser.add_argument(
        "--all-flow-output-chunk-size",
        type=int,
        default=DEFAULT_ALL_FLOW_CHUNK_SIZE,
        help="Symbolica output chunk size used when generating all-flow artifacts.",
    )
    parser.add_argument(
        "--cpp-time-limit",
        type=float,
        default=DEFAULT_CPP_TIME_LIMIT,
        help="Generation timeout in seconds for the C++ O3 row; use 0 to disable.",
    )
    parser.add_argument(
        "--only-missing",
        action="store_true",
        help="Skip rows already recorded with status ok in z_performance_data.json.",
    )
    parser.add_argument(
        "--no-seed-lc-cache",
        action="store_true",
        help=(
            "Do not seed AmpliCol and JIT O1 rows from result_matrix_data.json. "
            "Without a seed, amplicol rows are left untouched by this script."
        ),
    )
    parser.add_argument(
        "--force-pyamplicol-regeneration",
        action="store_true",
        help=(
            "Regenerate every requested pyAmpliCol selected/all-flow artifact "
            "with the installed code, while still seeding the AmpliCol "
            "reference from result_matrix_data.json."
        ),
    )
    parser.add_argument(
        "--retime-existing",
        action="store_true",
        help=(
            "Retain existing pyAmpliCol artifacts and generation times, but "
            "rerun both selected/all-flow time-process measurements."
        ),
    )
    args = parser.parse_args(argv)
    if args.retime_existing and args.force_pyamplicol_regeneration:
        parser.error(
            "--retime-existing and --force-pyamplicol-regeneration are mutually exclusive"
        )

    n_values = _parse_n_values(args.n)
    modes = _parse_modes(args.modes)
    for n_final in n_values:
        for mode in modes:
            if args.only_missing and _row_is_terminal(n_final, mode):
                print(
                    f"[z-table] skip n={n_final} mode={mode}: already recorded",
                    flush=True,
                )
                continue
            can_seed_mode = mode == "amplicol" or (
                mode == "jit_o1" and not args.force_pyamplicol_regeneration
            )
            if not args.no_seed_lc_cache and can_seed_mode:
                if _seed_from_lc_cache(n_final, mode, args=args):
                    continue
            if mode == "amplicol":
                print(
                    f"[z-table] n={n_final} mode=amplicol: no LC cache seed; "
                    "run docs/result_matrix.py for process id 1 first",
                    flush=True,
                )
                continue
            if args.retime_existing:
                _retime_pyamplicol_mode(n_final, mode, args)
                continue
            _run_pyamplicol_mode(n_final, mode, args)
    _render_table()
    return 0


def process_for_n(n_final: int) -> str:
    if n_final < 1:
        raise ValueError("n must be at least 1")
    return "d d~ > z" + (" " + " ".join("g" for _ in range(n_final - 1)) if n_final > 1 else "")


def _retime_pyamplicol_mode(
    n_final: int,
    mode: str,
    args: argparse.Namespace,
) -> None:
    existing = _existing_mode_record(n_final, mode)
    output_dir = _artifact_path(
        existing.get("output_dir"),
        fallback=OUTPUT_ROOT / f"n{n_final}" / mode,
    )
    all_flow_output_dir = _retime_all_flow_artifact_path(
        n_final,
        mode,
        existing=existing,
        args=args,
    )
    if not (output_dir / "process_manifest.json").is_file():
        _record_row(
            n_final,
            mode,
            status="error",
            error=f"existing selected artifact is missing: {output_dir}",
            output_dir=str(output_dir),
            all_flow_output_dir=str(all_flow_output_dir),
        )
        return

    print(f"[z-table] retime n={n_final} mode={mode}", flush=True)
    try:
        timed = _time_existing_artifact(
            output_dir,
            target_runtime=float(args.target_runtime),
            batch_size=int(args.batch_size),
        )
    except Exception as exc:  # noqa: BLE001 - benchmark script records failures.
        _record_row(
            n_final,
            mode,
            status="error",
            generation_s=_optional_float(existing.get("generation_s")),
            output_dir=str(output_dir),
            all_flow_output_dir=str(all_flow_output_dir),
            error=str(exc),
            notes=f"existing process kept at {_display_path(output_dir)}",
        )
        return

    profile = timed.get("profile", {})
    if not isinstance(profile, dict):
        profile = {}
    all_flow_record: dict[str, Any] = {
        "all_flow_status": str(existing.get("all_flow_status", "missing")),
        "all_flow_generation_s": _optional_float(
            existing.get("all_flow_generation_s")
        ),
        "all_flow_output_dir": str(all_flow_output_dir),
        "all_flow_time_batch_size": int(args.all_flow_batch_size),
        "all_flow_symbolica_output_chunk_size": int(
            args.all_flow_output_chunk_size
        ),
        "all_flow_notes": (
            f"retimed existing process kept at {_display_path(all_flow_output_dir)}; "
            f"time batch {int(args.all_flow_batch_size)}, output chunk "
            f"{int(args.all_flow_output_chunk_size)}"
        ),
    }
    if (all_flow_output_dir / "process_manifest.json").is_file():
        try:
            all_flow_timed = _time_existing_artifact(
                all_flow_output_dir,
                target_runtime=float(args.target_runtime),
                batch_size=int(args.all_flow_batch_size),
            )
            all_flow_profile = all_flow_timed.get("profile", {})
            if not isinstance(all_flow_profile, dict):
                all_flow_profile = {}
            all_flow_record.update(
                {
                    "all_flow_status": "ok",
                    "all_flow_wall_us_per_point": _optional_float(
                        all_flow_profile.get("wall_us_per_point")
                    ),
                    "all_flow_runtime_us_per_point": _optional_float(
                        all_flow_profile.get("core_evaluator_us_per_point")
                    ),
                }
            )
        except Exception as exc:  # noqa: BLE001 - preserve selected timing.
            all_flow_record.update(
                {
                    "all_flow_status": "error",
                    "all_flow_error": str(exc),
                }
            )

    _record_row(
        n_final,
        mode,
        status="ok",
        generation_s=_optional_float(existing.get("generation_s")),
        output_dir=str(output_dir),
        wall_us_per_point=_optional_float(profile.get("wall_us_per_point")),
        runtime_us_per_point=_optional_float(
            profile.get("core_evaluator_us_per_point")
        ),
        **all_flow_record,
        notes=(
            f"retimed existing process kept at {_display_path(output_dir)}"
        ),
    )


def _artifact_path(raw: object, *, fallback: Path) -> Path:
    if raw is None or not str(raw).strip():
        return fallback
    path = Path(str(raw))
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def _retime_all_flow_artifact_path(
    n_final: int,
    mode: str,
    *,
    existing: dict[str, Any],
    args: argparse.Namespace,
) -> Path:
    fallback = OUTPUT_ROOT / f"n{n_final}" / f"{mode}_all_flows"
    stored = existing.get("all_flow_output_dir")
    if stored is not None and str(stored).strip():
        return _artifact_path(stored, fallback=fallback)
    if mode != "jit_o1":
        return fallback

    lc_record = _jit_o1_all_flow_record_from_lc(n_final, args=args)
    if lc_record is None:
        return fallback
    lc_path = lc_record.get("all_flow_output_dir")
    if lc_path is None or not str(lc_path).strip():
        return fallback
    existing_generation = _optional_float(existing.get("all_flow_generation_s"))
    lc_generation = _optional_float(lc_record.get("all_flow_generation_s"))
    if (
        existing_generation is not None
        and lc_generation is not None
        and not _same_benchmark_value(existing_generation, lc_generation)
    ):
        return fallback
    recovered = _artifact_path(lc_path, fallback=fallback)
    print(
        f"[z-table] recovered n={n_final} mode=jit_o1 all-flow artifact "
        f"from LC cache: {_display_path(recovered)}",
        flush=True,
    )
    return recovered


def _same_benchmark_value(left: float, right: float) -> bool:
    return abs(left - right) <= max(1e-9, 1e-9 * max(abs(left), abs(right)))


def _display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT.resolve()))
    except ValueError:
        return str(path)


def _time_existing_artifact(
    output_dir: Path,
    *,
    target_runtime: float,
    batch_size: int,
) -> dict[str, Any]:
    return _run_json_command(
        [
            str(PYTHON),
            "-m",
            "pyamplicol",
            "time-process",
            "--target-runtime",
            str(float(target_runtime)),
            "--batch-size",
            str(int(batch_size)),
            "--json",
            str(output_dir),
        ],
        timeout=None,
        log_path=output_dir.with_name(f"{output_dir.name}.time.log"),
    )


def _existing_mode_record(n_final: int, mode: str) -> dict[str, Any]:
    try:
        data = json.loads(Z_DATA.read_text(encoding="utf-8"))
        row = data["entries"][str(n_final)]["modes"][mode]
    except (FileNotFoundError, KeyError, json.JSONDecodeError):
        return {}
    return dict(row) if isinstance(row, dict) else {}


def _run_pyamplicol_mode(n_final: int, mode: str, args: argparse.Namespace) -> None:
    process = process_for_n(n_final)
    output_dir = OUTPUT_ROOT / f"n{n_final}" / mode
    all_flow_output_dir = OUTPUT_ROOT / f"n{n_final}" / f"{mode}_all_flows"
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    _preserve_generated_output(output_dir)
    print(f"[z-table] start n={n_final} mode={mode}: {process}", flush=True)
    generate = _generate_command(
        process,
        output_dir,
        mode=mode,
        n_final=n_final,
        n_cores=max(1, int(args.n_cores)),
        batch_size=int(args.batch_size),
        output_chunk_size=int(args.output_chunk_size),
        all_flows=False,
    )
    all_flow_generate = _generate_command(
        process,
        all_flow_output_dir,
        mode=mode,
        n_final=n_final,
        n_cores=max(1, int(args.n_cores)),
        batch_size=int(args.all_flow_batch_size),
        output_chunk_size=int(args.all_flow_output_chunk_size),
        all_flows=True,
    )
    timeout = (
        None
        if mode != "cpp_o3" or float(args.cpp_time_limit) <= 0
        else float(args.cpp_time_limit)
    )
    try:
        gen = _run_json_command(
            generate,
            timeout=timeout,
            log_path=output_dir.with_name(f"{output_dir.name}.generate.log"),
        )
    except MemoryLimitExceeded as exc:
        _record_row(
            n_final,
            mode,
            status="ram_limit",
            notes=f"generated process kept at pyAmpliCol/docs/.z_performance_outputs/n{n_final}/{mode}",
            error=str(exc),
        )
        print(f"[z-table] RAM limit n={n_final} mode={mode}: {exc}", flush=True)
        return
    except subprocess.TimeoutExpired:
        _record_row(
            n_final,
            mode,
            status="timeout",
            notes=f"generated process kept at pyAmpliCol/docs/.z_performance_outputs/n{n_final}/{mode}",
        )
        print(f"[z-table] timeout n={n_final} mode={mode}", flush=True)
        return
    except Exception as exc:  # noqa: BLE001 - benchmark script records failures.
        _record_row(n_final, mode, status="error", error=str(exc))
        print(f"[z-table] error n={n_final} mode={mode}: {exc}", flush=True)
        return

    try:
        timed = _run_json_command(
            [
                str(PYTHON),
                "-m",
                "pyamplicol",
                "time-process",
                "--target-runtime",
                str(float(args.target_runtime)),
                "--batch-size",
                str(int(args.batch_size)),
                "--json",
                str(output_dir),
            ],
            timeout=None,
            log_path=output_dir.with_name(f"{output_dir.name}.time.log"),
        )
    except MemoryLimitExceeded as exc:
        _record_row(
            n_final,
            mode,
            status="ram_limit",
            notes=f"generated process kept at pyAmpliCol/docs/.z_performance_outputs/n{n_final}/{mode}",
            error=str(exc),
        )
        print(f"[z-table] timing RAM limit n={n_final} mode={mode}: {exc}", flush=True)
        return
    except Exception as exc:  # noqa: BLE001 - benchmark script records failures.
        _record_row(n_final, mode, status="error", error=str(exc))
        print(f"[z-table] timing error n={n_final} mode={mode}: {exc}", flush=True)
        return

    profile = timed.get("profile", {})
    if not isinstance(profile, dict):
        profile = {}
    generation_s = _optional_float(gen.get("_command_elapsed_s", gen.get("generation_s")))
    internal_generation_s = _optional_float(gen.get("generation_s"))
    wall = _optional_float(profile.get("wall_us_per_point"))
    runtime = _optional_float(profile.get("core_evaluator_us_per_point"))
    reused_all_flow_record = (
        _jit_o1_all_flow_record_from_lc(n_final, args=args)
        if mode == "jit_o1" and not args.force_pyamplicol_regeneration
        else None
    )
    all_flow_record: dict[str, Any] = reused_all_flow_record or {
        "all_flow_status": "missing",
        "all_flow_output_dir": str(all_flow_output_dir),
        "all_flow_notes": (
            f"generated process kept at "
            f"pyAmpliCol/docs/.z_performance_outputs/n{n_final}/{mode}_all_flows; "
            f"time batch {int(args.all_flow_batch_size)}, output chunk "
            f"{int(args.all_flow_output_chunk_size)}"
        ),
    }
    if reused_all_flow_record is not None:
        print(
            f"[z-table] reused matching LC all-flow O1 artifact for n={n_final}",
            flush=True,
        )
    else:
        try:
            _preserve_generated_output(all_flow_output_dir)
            all_flow_gen = _run_json_command(
                all_flow_generate,
                timeout=timeout,
                log_path=all_flow_output_dir.with_name(
                    f"{all_flow_output_dir.name}.generate.log"
                ),
            )
            all_flow_timed = _run_json_command(
                [
                    str(PYTHON),
                    "-m",
                    "pyamplicol",
                    "time-process",
                    "--target-runtime",
                    str(float(args.target_runtime)),
                    "--batch-size",
                    str(int(args.all_flow_batch_size)),
                    "--json",
                    str(all_flow_output_dir),
                ],
                timeout=None,
                log_path=all_flow_output_dir.with_name(
                    f"{all_flow_output_dir.name}.time.log"
                ),
            )
            all_flow_profile = all_flow_timed.get("profile", {})
            if not isinstance(all_flow_profile, dict):
                all_flow_profile = {}
            all_flow_record.update(
                {
                    "all_flow_status": "ok",
                    "all_flow_output_dir": str(all_flow_output_dir),
                    "all_flow_generation_s": _optional_float(
                        all_flow_gen.get(
                            "_command_elapsed_s",
                            all_flow_gen.get("generation_s"),
                        )
                    ),
                    "all_flow_wall_us_per_point": _optional_float(
                        all_flow_profile.get("wall_us_per_point")
                    ),
                    "all_flow_runtime_us_per_point": _optional_float(
                        all_flow_profile.get("core_evaluator_us_per_point")
                    ),
                    "all_flow_time_batch_size": int(args.all_flow_batch_size),
                    "all_flow_symbolica_output_chunk_size": int(
                        args.all_flow_output_chunk_size
                    ),
                }
            )
        except MemoryLimitExceeded as exc:
            all_flow_record.update(
                {
                    "all_flow_status": "ram_limit",
                    "all_flow_error": str(exc),
                }
            )
            print(
                f"[z-table] all-flow RAM limit n={n_final} mode={mode}: {exc}",
                flush=True,
            )
        except subprocess.TimeoutExpired:
            all_flow_record["all_flow_status"] = "timeout"
            print(f"[z-table] all-flow timeout n={n_final} mode={mode}", flush=True)
        except Exception as exc:  # noqa: BLE001 - benchmark script records failures.
            all_flow_record.update(
                {
                    "all_flow_status": "error",
                    "all_flow_error": str(exc),
                }
            )
            print(
                f"[z-table] all-flow error n={n_final} mode={mode}: {exc}",
                flush=True,
            )
    _record_row(
        n_final,
        mode,
        status="ok",
        generation_s=generation_s,
        output_dir=str(output_dir),
        wall_us_per_point=wall,
        runtime_us_per_point=runtime,
        **all_flow_record,
        notes=(
            f"generated process kept at "
            f"pyAmpliCol/docs/.z_performance_outputs/n{n_final}/{mode}"
        ),
    )
    print(
        f"[z-table] done n={n_final} mode={mode} "
        f"(selected gen {internal_generation_s or generation_s})",
        flush=True,
    )


def _generate_command(
    process: str,
    output_dir: Path,
    *,
    mode: str,
    n_final: int,
    n_cores: int,
    batch_size: int,
    output_chunk_size: int,
    all_flows: bool,
) -> list[str]:
    reference_order, sector_ids = _lc_cache_color_settings(n_final)
    if reference_order is None:
        reference_order = _default_reference_color_order(n_final)
    if sector_ids is None:
        sector_ids = [0]
    command = [
        str(PYTHON),
        "-m",
        "pyamplicol",
        "generate-process",
        process,
        str(output_dir),
        "--replace",
        "--n_cores",
        str(n_cores),
        "--color-accuracy",
        "lc",
        "--batch-size",
        str(batch_size),
        "--symbolica-n-cores",
        str(n_cores),
        "--json",
        "--monitor",
        "--symbolica-output-chunk-size",
        str(output_chunk_size),
        "--symbolica-output-chunk-strategy",
        "uniform",
        "--symbolica-stage-local-parameter-layout",
        "--symbolica-iterations",
        str(DEFAULT_SYMBOLICA_ITERATIONS),
        "--symbolica-max-horner-scheme-variables",
        "1000",
        "--symbolica-max-common-pair-cache-entries",
        "5000000",
        "--symbolica-max-common-pair-distance",
        "1000",
    ]
    if all_flows:
        fixed_helicity = _fixed_helicity_choice(process)
        command.extend(
            [
                "--lc-sector-strategy",
                "all",
                "--skip-generic-plan",
                "--no-runtime-lc-sector-selector",
                "--source-helicities",
                str(fixed_helicity["source_helicities_cli"]),
            ]
        )
    else:
        command.extend(
            [
                "--lc-sector-ids",
                ",".join(str(item) for item in sector_ids),
                "--reference-color-order",
                ",".join(str(item) for item in reference_order),
            ]
        )
    if mode == "jit_o1":
        command.extend(
            [
                "--symbolica-evaluator-backend",
                "jit",
                "--symbolica-jit-optimization-level",
                "1",
                "--symbolica-compiled-chunk-compile-workers",
                str(n_cores),
            ]
        )
    elif mode == "jit_o3":
        command.extend(
            [
                "--symbolica-evaluator-backend",
                "jit",
                "--symbolica-jit-optimization-level",
                "3",
                "--symbolica-compiled-chunk-compile-workers",
                str(n_cores),
            ]
        )
    elif mode == "asm":
        command.extend(
            [
                "--symbolica-evaluator-backend",
                "compiled-complex",
                "--symbolica-compiled-preset",
                "generation",
                "--symbolica-compiled-inline-asm",
                "default",
                "--symbolica-compiled-optimization-level",
                "3",
                "--symbolica-compiled-chunk-compile-workers",
                str(n_cores),
            ]
        )
    elif mode == "cpp_o3":
        command.extend(
            [
                "--symbolica-evaluator-backend",
                "compiled-complex",
                "--symbolica-compiled-preset",
                "runtime-o3",
                "--symbolica-compiled-inline-asm",
                "none",
                "--symbolica-compiled-optimization-level",
                "3",
                "--symbolica-compiled-chunk-compile-workers",
                str(n_cores),
            ]
        )
    else:
        raise ValueError(f"unsupported pyAmpliCol mode: {mode}")
    return command


def _seed_from_lc_cache(
    n_final: int,
    mode: str,
    *,
    args: argparse.Namespace,
) -> bool:
    if not LC_MATRIX_DATA.exists():
        return False
    try:
        data = json.loads(LC_MATRIX_DATA.read_text(encoding="utf-8"))
        case = data["entries"]["dd_z_jets"][str(n_final)]
    except (KeyError, json.JSONDecodeError):
        return False
    if mode == "amplicol":
        row = case.get("amplicol", {})
        if not isinstance(row, dict) or row.get("status") != "ok":
            return False
        _record_row(
            n_final,
            "amplicol",
            status="ok",
            generation_s=_optional_float(row.get("generation_s")),
            runtime_us_per_point=_optional_float(row.get("runtime_us_per_point")),
            all_flow_status=(
                "ok"
                if _optional_float(row.get("all_flow_runtime_us_per_point")) is not None
                else "missing"
            ),
            all_flow_generation_s=_optional_float(row.get("all_flow_generation_s")),
            all_flow_runtime_us_per_point=_optional_float(
                row.get("all_flow_runtime_us_per_point")
            ),
            all_flow_time_batch_size=_optional_int(row.get("all_flow_time_batch_size")),
            all_flow_symbolica_output_chunk_size=_optional_int(
                row.get("all_flow_symbolica_output_chunk_size")
            ),
            all_flow_notes="reused from LC result matrix cache",
            notes="reused from LC result matrix cache",
        )
        print(f"[z-table] seeded n={n_final} mode=amplicol from LC cache", flush=True)
        return True
    row = case.get("pyamplicol_jit", {})
    if not isinstance(row, dict) or row.get("status") != "ok":
        return False
    settings = row.get("matrix_settings")
    selected_batch_size = (
        _optional_int(settings.get("selected_runtime_batch_size"))
        if isinstance(settings, dict)
        else None
    )
    if selected_batch_size != int(args.batch_size):
        print(
            f"[z-table] LC JIT O1 seed n={n_final} uses selected batch "
            f"{selected_batch_size}; generating dedicated batch "
            f"{int(args.batch_size)} artifact",
            flush=True,
        )
        return False
    _record_row(
        n_final,
        "jit_o1",
        status="ok",
        generation_s=_selected_generation_from_lc_row(row),
        output_dir=str(row.get("selected_output_dir", "")),
        wall_us_per_point=_optional_float(row.get("wall_us_per_point")),
        runtime_us_per_point=_optional_float(row.get("runtime_us_per_point")),
        all_flow_status=(
            "ok"
            if _optional_float(row.get("all_flow_runtime_us_per_point")) is not None
            else "missing"
        ),
        all_flow_generation_s=(
            _optional_float(row.get("all_flow_generation_s"))
            or _optional_float(row.get("generation_s"))
        ),
        all_flow_output_dir=str(row.get("all_flow_output_dir", "")),
        all_flow_wall_us_per_point=_optional_float(
            row.get("all_flow_wall_us_per_point")
        ),
        all_flow_runtime_us_per_point=_optional_float(
            row.get("all_flow_runtime_us_per_point")
        ),
        all_flow_time_batch_size=_optional_int(row.get("all_flow_time_batch_size")),
        all_flow_symbolica_output_chunk_size=_optional_int(
            row.get("all_flow_symbolica_output_chunk_size")
        ),
        all_flow_notes="reused from LC result matrix cache",
        notes="reused from LC result matrix cache",
    )
    print(f"[z-table] seeded n={n_final} mode=jit_o1 from LC cache", flush=True)
    return True


def _jit_o1_all_flow_record_from_lc(
    n_final: int,
    *,
    args: argparse.Namespace,
) -> dict[str, Any] | None:
    if not LC_MATRIX_DATA.exists():
        return None
    try:
        data = json.loads(LC_MATRIX_DATA.read_text(encoding="utf-8"))
        row = data["entries"]["dd_z_jets"][str(n_final)]["pyamplicol_jit"]
    except (KeyError, json.JSONDecodeError):
        return None
    if not isinstance(row, dict) or row.get("status") != "ok":
        return None
    settings = row.get("matrix_settings")
    if not isinstance(settings, dict):
        return None
    if _optional_int(settings.get("all_flow_runtime_batch_size")) != int(
        args.all_flow_batch_size
    ):
        return None
    if _optional_int(settings.get("all_flow_symbolica_output_chunk_size")) != int(
        args.all_flow_output_chunk_size
    ):
        return None
    if _optional_int(settings.get("all_flow_symbolica_jit_optimization_level")) != 1:
        return None
    generation_s = _optional_float(row.get("all_flow_generation_s"))
    wall = _optional_float(row.get("all_flow_wall_us_per_point"))
    runtime = _optional_float(row.get("all_flow_runtime_us_per_point"))
    if generation_s is None or wall is None or runtime is None:
        return None
    return {
        "all_flow_status": "ok",
        "all_flow_generation_s": generation_s,
        "all_flow_output_dir": str(row.get("all_flow_output_dir", "")),
        "all_flow_wall_us_per_point": wall,
        "all_flow_runtime_us_per_point": runtime,
        "all_flow_time_batch_size": int(args.all_flow_batch_size),
        "all_flow_symbolica_output_chunk_size": int(
            args.all_flow_output_chunk_size
        ),
        "all_flow_notes": (
            "reused matching O1 all-flow artifact from LC result matrix cache; "
            f"time batch {int(args.all_flow_batch_size)}, output chunk "
            f"{int(args.all_flow_output_chunk_size)}"
        ),
    }


def _preserve_generated_output(output_dir: Path) -> None:
    stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
    if output_dir.exists():
        preserved = output_dir.with_name(f"{output_dir.name}_before_{stamp}")
        counter = 1
        while preserved.exists():
            preserved = output_dir.with_name(
                f"{output_dir.name}_before_{stamp}_{counter}"
            )
            counter += 1
        output_dir.rename(preserved)
        print(f"[z-table] preserved {output_dir} as {preserved}", flush=True)
    for suffix in ("generate.log", "time.log"):
        log_path = output_dir.with_name(f"{output_dir.name}.{suffix}")
        if not log_path.exists():
            continue
        preserved_log = log_path.with_name(
            f"{output_dir.name}.{suffix}.before_{stamp}"
        )
        counter = 1
        while preserved_log.exists():
            preserved_log = log_path.with_name(
                f"{output_dir.name}.{suffix}.before_{stamp}_{counter}"
            )
            counter += 1
        log_path.rename(preserved_log)


def _selected_generation_from_lc_row(row: dict[str, Any]) -> float | None:
    explicit = _optional_float(row.get("selected_generation_s"))
    if explicit is not None:
        return explicit
    selected = row.get("selected_generate_payload")
    if isinstance(selected, dict):
        generation = _optional_float(selected.get("_command_elapsed_s"))
        if generation is not None:
            return generation
        generation = _optional_float(selected.get("generation_s"))
        if generation is not None:
            return generation
    return _optional_float(row.get("generation_s"))


def _fixed_helicity_choice(process: str) -> dict[str, Any]:
    import result_matrix

    base = next(item for item in result_matrix.BASE_PROCESSES if item.key == "dd_z_jets")
    return result_matrix._fixed_helicity_choice(process, base)


def _lc_cache_color_settings(n_final: int) -> tuple[list[int] | None, list[int] | None]:
    if not LC_MATRIX_DATA.exists():
        return None, None
    try:
        data = json.loads(LC_MATRIX_DATA.read_text(encoding="utf-8"))
        case = data["entries"]["dd_z_jets"][str(n_final)]
    except (KeyError, json.JSONDecodeError):
        return None, None
    reference_order = None
    sector_ids = None
    for row_key in ("pyamplicol_jit", "amplicol"):
        row = case.get(row_key, {})
        if not isinstance(row, dict):
            continue
        raw_order = row.get("reference_color_order")
        if isinstance(raw_order, list) and raw_order:
            reference_order = [int(item) for item in raw_order]
        raw_sectors = row.get("selected_lc_sector_ids")
        if isinstance(raw_sectors, list) and raw_sectors:
            sector_ids = [int(item) for item in raw_sectors]
        if reference_order is not None and sector_ids is not None:
            break
    return reference_order, sector_ids


def _default_reference_color_order(n_final: int) -> list[int]:
    return [2, *range(4, n_final + 3), 1, 3]


def _record_row(
    n_final: int,
    mode: str,
    *,
    status: str,
    generation_s: float | None = None,
    output_dir: str = "",
    wall_us_per_point: float | None = None,
    runtime_us_per_point: float | None = None,
    all_flow_status: str | None = None,
    all_flow_generation_s: float | None = None,
    all_flow_output_dir: str = "",
    all_flow_wall_us_per_point: float | None = None,
    all_flow_runtime_us_per_point: float | None = None,
    all_flow_time_batch_size: int | None = None,
    all_flow_symbolica_output_chunk_size: int | None = None,
    all_flow_notes: str = "",
    all_flow_error: str = "",
    notes: str = "",
    error: str = "",
) -> None:
    command = [
        str(PYTHON),
        str(TABLE_SCRIPT),
        "record",
        "--n",
        str(n_final),
        "--mode",
        mode,
        "--status",
        status,
        "--notes",
        notes,
        "--error",
        error,
    ]
    if generation_s is not None:
        command.extend(["--generation-s", str(generation_s)])
    if output_dir:
        command.extend(["--output-dir", output_dir])
    if wall_us_per_point is not None:
        command.extend(["--wall-us-per-point", str(wall_us_per_point)])
    if runtime_us_per_point is not None:
        command.extend(["--runtime-us-per-point", str(runtime_us_per_point)])
    if all_flow_status is not None:
        command.extend(["--all-flow-status", str(all_flow_status)])
    if all_flow_generation_s is not None:
        command.extend(["--all-flow-generation-s", str(all_flow_generation_s)])
    if all_flow_output_dir:
        command.extend(["--all-flow-output-dir", all_flow_output_dir])
    if all_flow_wall_us_per_point is not None:
        command.extend(["--all-flow-wall-us-per-point", str(all_flow_wall_us_per_point)])
    if all_flow_runtime_us_per_point is not None:
        command.extend(
            ["--all-flow-runtime-us-per-point", str(all_flow_runtime_us_per_point)]
        )
    if all_flow_time_batch_size is not None:
        command.extend(["--all-flow-time-batch-size", str(all_flow_time_batch_size)])
    if all_flow_symbolica_output_chunk_size is not None:
        command.extend(
            [
                "--all-flow-symbolica-output-chunk-size",
                str(all_flow_symbolica_output_chunk_size),
            ]
        )
    if all_flow_notes:
        command.extend(["--all-flow-notes", all_flow_notes])
    if all_flow_error:
        command.extend(["--all-flow-error", all_flow_error])
    _run_checked(command)


def _render_table() -> None:
    _run_checked([str(PYTHON), str(TABLE_SCRIPT), "render"])


def _run_checked(command: Sequence[str]) -> None:
    env = _env()
    completed = subprocess.run(
        list(command),
        cwd=REPO_ROOT,
        env=env,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"command failed with exit code {completed.returncode}: {' '.join(command)}")


def _run_json_command(
    command: Sequence[str],
    *,
    timeout: float | None,
    log_path: Path,
) -> dict[str, Any]:
    env = _env()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[z-table] progress log: {log_path}", flush=True)
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log_file:
        log_file.write("# command: " + " ".join(str(arg) for arg in command) + "\n")
        process = subprocess.Popen(
            list(command),
            cwd=REPO_ROOT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=log_file,
            start_new_session=True,
        )
        try:
            stdout = _wait_with_heartbeat(process, command, started=started, timeout=timeout)
        except BaseException:
            _terminate_process_group(process.pid, signal.SIGTERM)
            try:
                process.communicate(timeout=10.0)
            except subprocess.TimeoutExpired:
                _terminate_process_group(process.pid, signal.SIGKILL)
                process.communicate()
            raise
    elapsed_s = time.perf_counter() - started
    if process.returncode != 0:
        raise RuntimeError(
            f"command failed with exit code {process.returncode}: {' '.join(command)}\n"
            f"stdout:\n{(stdout or '')[-4000:]}\n"
            f"stderr:\n{_read_tail(log_path)[-4000:]}"
        )
    payload = _parse_json_output(stdout or "")
    payload["_command_elapsed_s"] = elapsed_s
    payload["_progress_log"] = str(log_path)
    return payload


def _wait_with_heartbeat(
    process: subprocess.Popen[str],
    command: Sequence[str],
    *,
    started: float,
    timeout: float | None,
) -> str:
    deadline = None if timeout is None else started + timeout
    memory_limit_bytes = int(DEFAULT_MEMORY_LIMIT_GB * 1024**3)
    memory_monitor = ProcessTreeMemoryMonitor(
        physical_footprint_high_watermark_bytes=int(0.8 * memory_limit_bytes)
    )
    last_heartbeat = started
    while True:
        memory_sample = memory_monitor.sample(process.pid)
        memory_bytes = memory_sample.effective_bytes
        if memory_bytes is not None:
            if memory_bytes > memory_limit_bytes:
                raise MemoryLimitExceeded(
                    limit_gb=DEFAULT_MEMORY_LIMIT_GB,
                    peak_memory_bytes=memory_bytes,
                    memory_metric=memory_sample.effective_metric or "memory",
                    command=command,
                )
        wait_s = MEMORY_POLL_S
        if deadline is not None:
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                raise subprocess.TimeoutExpired(command, timeout)
            wait_s = min(wait_s, remaining)
        try:
            stdout, _ = process.communicate(timeout=wait_s)
            return stdout or ""
        except subprocess.TimeoutExpired:
            if deadline is not None and time.perf_counter() >= deadline:
                raise
            elapsed = time.perf_counter() - started
            if time.perf_counter() - last_heartbeat >= HEARTBEAT_S:
                memory_note = (
                    ""
                    if memory_monitor.peak_memory_bytes is None
                    else (
                        ", peak memory "
                        f"{memory_monitor.peak_memory_bytes / 1024**3:.2f} GiB"
                        f" ({(memory_monitor.peak_memory_metric or 'unknown').replace('_', ' ')})"
                    )
                )
                print(
                    f"[z-table] still running after {elapsed:.0f}s{memory_note}: "
                    f"{_command_label(command)}",
                    flush=True,
                )
                last_heartbeat = time.perf_counter()


def _command_label(command: Sequence[str]) -> str:
    if "generate-process" in command:
        idx = list(command).index("generate-process")
        return "generate-process " + " ".join(str(item) for item in command[idx + 1 : idx + 4])
    if "time-process" in command:
        return "time-process"
    return " ".join(str(item) for item in command[:5])


def _terminate_process_group(pid: int, sig: int) -> None:
    try:
        os.killpg(pid, sig)
    except ProcessLookupError:
        return


def _parse_json_output(text: str) -> dict[str, Any]:
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        start = text.find("{")
        end = text.rfind("}")
        if start < 0 or end < start:
            raise
        value = json.loads(text[start : end + 1])
    if not isinstance(value, dict):
        raise TypeError("expected JSON object output")
    return value


def _read_tail(path: Path, *, max_bytes: int = 16000) -> str:
    if not path.exists():
        return ""
    with path.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - max_bytes))
        return handle.read().decode("utf-8", errors="replace")


def _row_is_terminal(n_final: int, mode: str) -> bool:
    if not Z_DATA.exists():
        return False
    try:
        data = json.loads(Z_DATA.read_text(encoding="utf-8"))
        row = data["entries"][str(n_final)]["modes"][mode]
    except (KeyError, json.JSONDecodeError):
        return False
    if not isinstance(row, dict) or row.get("status") not in {"ok", "ram_limit"}:
        return False
    if row.get("status") == "ram_limit":
        return True
    return row.get("all_flow_status") in {"ok", "ram_limit"}


def _parse_n_values(values: Sequence[str]) -> list[int]:
    result: set[int] = set()
    for token in _split_tokens(values):
        if "-" in token:
            left, right = token.split("-", 1)
            start = int(left)
            end = int(right)
            if start > end:
                start, end = end, start
            result.update(range(start, end + 1))
        else:
            result.add(int(token))
    return [n for n in sorted(result) if 1 <= n <= 9]


def _parse_modes(values: Sequence[str]) -> list[str]:
    tokens = _split_tokens(values)
    if not tokens or "all" in tokens:
        return list(MODE_KEYS)
    unknown = sorted(set(tokens) - set(MODE_KEYS))
    if unknown:
        raise ValueError(f"unknown modes: {', '.join(unknown)}")
    return [mode for mode in MODE_KEYS if mode in tokens]


def _split_tokens(values: Sequence[str]) -> list[str]:
    tokens: list[str] = []
    for value in values:
        tokens.extend(part.strip() for part in str(value).split(",") if part.strip())
    return tokens


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _optional_int(value: object) -> int | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float) and value.is_integer():
        return int(value)
    return None


def _env() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = (
        f"{SRC_DIR}{os.pathsep}{env['PYTHONPATH']}"
        if env.get("PYTHONPATH")
        else str(SRC_DIR)
    )
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    return env


if __name__ == "__main__":
    raise SystemExit(main())
