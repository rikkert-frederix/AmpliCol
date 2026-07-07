#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

DOCS_DIR = Path(__file__).resolve().parent
PYAMPLICOL_DIR = DOCS_DIR.parent
REPO_ROOT = PYAMPLICOL_DIR.parent
SRC_DIR = PYAMPLICOL_DIR / "src"
DEFAULT_DATA = DOCS_DIR / "result_matrix_data.json"
DEFAULT_TABLE = DOCS_DIR / "result_matrix_table.tex"
DEFAULT_OUTPUT_ROOT = DOCS_DIR / ".result_matrix_outputs"
DEFAULT_TIME_LIMIT_S = 900.0
DEFAULT_TARGET_RUNTIME_S = 1.0
DEFAULT_AMPLICOL_POINTS = 10000
DEFAULT_JOBS = 4
DEFAULT_N_CORES = 4
DEFAULT_PROCESS_WORKERS = 1
DEFAULT_COLUMNS_PER_TABLE = 3
VALIDATION_REL_TOL = 1.0e-8
VALIDATION_ABS_TOL = 1.0e-16

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from pyamplicol.generic_validation import _library_runtime_per_point  # noqa: E402
from pyamplicol.core_types import ExternalMomentum  # noqa: E402
from pyamplicol.generic_artifact import (  # noqa: E402
    select_leading_color_sector_ids_from_plan,
)
from pyamplicol.phase_space import generic_validation_point  # noqa: E402
from pyamplicol.color_plan import (  # noqa: E402
    build_color_plan,
    lc_line_pairing_representative_ids,
)
from pyamplicol.process_ir import build_process_ir  # noqa: E402
from pyamplicol.processes import ProcessOptions  # noqa: E402
from pyamplicol.reference import (  # noqa: E402
    AmplicolAdapter,
    amplicol_process_file_entry,
    amplicol_process_file_integrals,
)


@dataclass(frozen=True)
class BaseProcess:
    key: str
    label: str
    initial: tuple[str, ...]
    base_final: tuple[str, ...]
    include_cc: bool = False
    include_3qqbar: bool = False
    include_resonance: bool = False
    max_quark_pairs: int | None = None
    max_currents: int | None = None
    max_color_sectors: int | None = None
    lc_sector_strategy: str | None = None

    @property
    def min_final_count(self) -> int:
        return len(self.base_final)

    def process_for_n(self, n_final: int) -> str | None:
        extra_gluons = n_final - len(self.base_final)
        if extra_gluons < 0:
            return None
        final = (*self.base_final, *(("g",) * extra_gluons))
        return f"{' '.join(self.initial)} > {' '.join(final)}"

    def process_options(self) -> ProcessOptions:
        return ProcessOptions(
            include_cc=self.include_cc,
            include_3qqbar=self.include_3qqbar,
            include_resonance=self.include_resonance,
        )

    def cli_flags(self) -> list[str]:
        flags: list[str] = []
        if self.include_cc:
            flags.append("--include-cc")
        if self.include_3qqbar:
            flags.append("--include-3qqbar")
        if self.include_resonance:
            flags.append("--include-resonance")
        if self.max_quark_pairs is not None:
            flags.extend(["--max-quark-pairs", str(self.max_quark_pairs)])
        if self.max_currents is not None:
            flags.extend(["--max-currents", str(self.max_currents)])
        if self.max_color_sectors is not None:
            flags.extend(["--max-color-sectors", str(self.max_color_sectors)])
        if self.lc_sector_strategy is not None:
            flags.extend(["--lc-sector-strategy", self.lc_sector_strategy])
        return flags


BASE_PROCESSES: tuple[BaseProcess, ...] = (
    BaseProcess(
        key="dd_z_jets",
        label=r"$d\bar d\to Z+(n-1)g$",
        initial=("d", "d~"),
        base_final=("z",),
    ),
    BaseProcess(
        key="ud_w_jets",
        label=r"$u\bar d\to W^++(n-1)g$",
        initial=("u", "d~"),
        base_final=("w+",),
        include_cc=True,
    ),
    BaseProcess(
        key="dd_epem_jets",
        label=r"$d\bar d\to e^+e^-+(n-2)g$",
        initial=("d", "d~"),
        base_final=("e+", "e-"),
    ),
    BaseProcess(
        key="ud_epve_jets",
        label=r"$u\bar d\to e^+\nu_e+(n-2)g$",
        initial=("u", "d~"),
        base_final=("e+", "ve"),
        include_cc=True,
    ),
    BaseProcess(
        key="dd_zz_jets",
        label=r"$d\bar d\to ZZ+(n-2)g$",
        initial=("d", "d~"),
        base_final=("z", "z"),
    ),
    BaseProcess(
        key="gg_tt_jets",
        label=r"$gg\to t\bar t+(n-2)g$",
        initial=("g", "g"),
        base_final=("t", "t~"),
    ),
    BaseProcess(
        key="dd_tt_jets",
        label=r"$d\bar d\to t\bar t+(n-2)g$",
        initial=("d", "d~"),
        base_final=("t", "t~"),
    ),
    BaseProcess(
        key="gg_gluons",
        label=r"$gg\to gg+(n-2)g$",
        initial=("g", "g"),
        base_final=("g", "g"),
    ),
    BaseProcess(
        key="dd_zzz_jets",
        label=r"$d\bar d\to ZZZ+(n-3)g$",
        initial=("d", "d~"),
        base_final=("z", "z", "z"),
    ),
    BaseProcess(
        key="dd_epemzh_jets",
        label=r"$d\bar d\to e^+e^-ZH+(n-4)g$",
        initial=("d", "d~"),
        base_final=("e+", "e-", "z", "h"),
    ),
    BaseProcess(
        key="dd_ttzh_jets",
        label=r"$d\bar d\to t\bar t ZH+(n-4)g$",
        initial=("d", "d~"),
        base_final=("t", "t~", "z", "h"),
    ),
    BaseProcess(
        key="dd_4l_jets",
        label=r"$d\bar d\to e^+e^-e^+e^-+(n-4)g$",
        initial=("d", "d~"),
        base_final=("e+", "e-", "e+", "e-"),
    ),
    BaseProcess(
        key="dd_3q_lines",
        label=r"$d\bar d\to u\bar{u}\,s\bar{s}+(n-4)g$",
        initial=("d", "d~"),
        base_final=("u", "u~", "s", "s~"),
        include_3qqbar=True,
        max_quark_pairs=3,
    ),
    BaseProcess(
        key="dd_4q_lines",
        label=r"$d\bar d\to u\bar{u}\,s\bar{s}\,c\bar{c}+(n-6)g$",
        initial=("d", "d~"),
        base_final=("u", "u~", "s", "s~", "c", "c~"),
        include_3qqbar=True,
        max_quark_pairs=4,
        max_currents=750000,
        max_color_sectors=50000,
        lc_sector_strategy="line-pairing-representatives",
    ),
)

BASE_PROCESS_IDS: dict[str, int] = {
    base.key: process_id for process_id, base in enumerate(BASE_PROCESSES, start=1)
}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Generate and render the pyAmpliCol LC performance matrix used by "
            "description.tex."
        )
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--generate-data",
        nargs="*",
        default=(),
        help="n values to regenerate, e.g. '1 2 3' or '1-3,6'.",
    )
    parser.add_argument(
        "--show-data",
        nargs="*",
        default=("1-9",),
        help="n values to show in the LaTeX table. Missing data is rendered N/A.",
    )
    parser.add_argument(
        "--time-limit",
        type=float,
        default=DEFAULT_TIME_LIMIT_S,
        help=(
            "Per-cell timeout in seconds for C++ O3 generation. "
            "AmpliCol and JIT rows are never capped by this option. "
            "Use 0 to disable the C++ O3 cap as well."
        ),
    )
    parser.add_argument(
        "--target-runtime",
        type=float,
        default=DEFAULT_TARGET_RUNTIME_S,
        help="Target runtime in seconds for pyAmpliCol time-process.",
    )
    parser.add_argument(
        "--amplicol-points",
        type=int,
        default=DEFAULT_AMPLICOL_POINTS,
        help="Repeated supplied-momenta probe points for AmpliCol runtime timing.",
    )
    parser.add_argument(
        "--columns-per-table",
        type=int,
        default=DEFAULT_COLUMNS_PER_TABLE,
        help="Maximum number of multiplicity columns per rendered landscape table.",
    )
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS)
    parser.add_argument("--n-cores", type=int, default=DEFAULT_N_CORES)
    parser.add_argument(
        "--process-workers",
        type=int,
        default=DEFAULT_PROCESS_WORKERS,
        help=(
            "Number of independent matrix cells to refresh concurrently. "
            "Fortran AmpliCol reference refreshes use shared top-level build "
            "files, so this is automatically forced to one unless "
            "--skip-amplicol is supplied."
        ),
    )
    parser.add_argument("--skip-amplicol", action="store_true")
    parser.add_argument("--skip-jit", action="store_true")
    parser.add_argument("--skip-cpp-o3", action="store_true")
    parser.add_argument(
        "--only-missing",
        action="store_true",
        help="Only run requested modes whose current table entry is missing or failed.",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help=(
            "Also run same-point validation against AmpliCol. The performance "
            "matrix refresh keeps this disabled by default so reference timings "
            "come from the direct generated-library benchmark only."
        ),
    )
    parser.add_argument(
        "--base-process",
        nargs="*",
        default=(),
        help="Optional base-process keys to regenerate, e.g. dd_z_jets gg_tt_jets.",
    )
    parser.add_argument(
        "--process-ids",
        nargs="*",
        default=(),
        help=(
            "Optional process-row ids to regenerate, in the order displayed in "
            "the matrix table. Accepts values like '1 3 7' or '1-4,8'."
        ),
    )
    parser.add_argument(
        "--no-recompile",
        action="store_true",
        help=(
            "Only write result_matrix_table.tex/result_matrix_data.json. By default "
            "the script also recompiles description.tex and refreshes pyAmpliCol.pdf."
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    data = _load_data(args.data)
    _refresh_data_metadata(data)
    generate_ns = _parse_n_values(args.generate_data)
    show_ns = _parse_n_values(args.show_data)
    selected_base_keys = {str(key) for key in args.base_process}
    try:
        selected_process_ids = set(_parse_process_ids(args.process_ids))
    except ValueError as exc:
        parser.error(str(exc))

    if generate_ns:
        work_items: list[tuple[BaseProcess, int, str]] = []
        for n_final in generate_ns:
            for process_id, base in enumerate(BASE_PROCESSES, start=1):
                if selected_base_keys and base.key not in selected_base_keys:
                    continue
                if selected_process_ids and process_id not in selected_process_ids:
                    continue
                process = base.process_for_n(n_final)
                if process is None:
                    _record_not_applicable(data, base, n_final)
                    _write_data(args.data, data)
                    continue
                if args.dry_run:
                    print(
                        f"[dry-run] id={process_id} n={n_final} "
                        f"{base.key}: {process}"
                    )
                    continue
                work_items.append((base, n_final, process))
        process_workers = max(1, int(args.process_workers))
        if process_workers > 1 and not args.skip_amplicol:
            print(
                "[matrix] --process-workers > 1 requested, but AmpliCol "
                "reference refreshes use shared top-level build files; "
                "forcing serial execution. Add --skip-amplicol for safe "
                "parallel pyAmpliCol-only refreshes.",
                flush=True,
            )
            process_workers = 1
        if process_workers > 1 and work_items:
            print(
                f"[matrix] running {len(work_items)} independent cells with "
                f"{process_workers} process workers",
                flush=True,
            )
            with ThreadPoolExecutor(max_workers=process_workers) as executor:
                futures = [
                    executor.submit(
                        _generate_case,
                        data,
                        base,
                        n_final,
                        process,
                        output_root=args.output_root,
                        time_limit=None if args.time_limit <= 0 else args.time_limit,
                        target_runtime=args.target_runtime,
                        amplicol_points=args.amplicol_points,
                        jobs=args.jobs,
                        n_cores=args.n_cores,
                        skip_amplicol=args.skip_amplicol,
                        skip_jit=args.skip_jit,
                        skip_cpp_o3=args.skip_cpp_o3,
                        only_missing=args.only_missing,
                        validate=args.validate,
                    )
                    for base, n_final, process in work_items
                ]
                for future in as_completed(futures):
                    future.result()
            _write_data(args.data, data)
        elif work_items:
            for base, n_final, process in work_items:
                _generate_case(
                    data,
                    base,
                    n_final,
                    process,
                    output_root=args.output_root,
                    time_limit=None if args.time_limit <= 0 else args.time_limit,
                    target_runtime=args.target_runtime,
                    amplicol_points=args.amplicol_points,
                    jobs=args.jobs,
                    n_cores=args.n_cores,
                    skip_amplicol=args.skip_amplicol,
                    skip_jit=args.skip_jit,
                    skip_cpp_o3=args.skip_cpp_o3,
                    only_missing=args.only_missing,
                    validate=args.validate,
                )
                _write_data(args.data, data)

    table = render_latex_table(data, show_ns, columns_per_table=args.columns_per_table)
    args.table.parent.mkdir(parents=True, exist_ok=True)
    args.table.write_text(table, encoding="utf-8")
    _write_data(args.data, data)
    print(f"wrote {args.table}")
    print(f"wrote {args.data}")
    if not args.no_recompile:
        _recompile_description_pdf(args.table)
    return 0


def _recompile_description_pdf(table_path: Path) -> None:
    if table_path.resolve() != DEFAULT_TABLE.resolve():
        print(
            "skipped PDF recompilation: --table does not point to "
            f"{DEFAULT_TABLE}",
            flush=True,
        )
        return
    description = DOCS_DIR / "description.tex"
    if not description.exists():
        print(f"skipped PDF recompilation: {description} does not exist", flush=True)
        return
    env = dict(os.environ)
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    print("recompiling pyAmpliCol/docs/description.tex", flush=True)
    completed = subprocess.run(
        ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "description.tex"],
        cwd=DOCS_DIR,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "description.tex recompilation failed with exit code "
            f"{completed.returncode}\nstdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    source = DOCS_DIR / "description.pdf"
    target = DOCS_DIR / "pyAmpliCol.pdf"
    target.write_bytes(source.read_bytes())
    print(f"wrote {source}")
    print(f"wrote {target}")


def _generate_case(
    data: dict[str, Any],
    base: BaseProcess,
    n_final: int,
    process: str,
    *,
    output_root: Path,
    time_limit: float | None,
    target_runtime: float,
    amplicol_points: int,
    jobs: int,
    n_cores: int,
    skip_amplicol: bool,
    skip_jit: bool,
    skip_cpp_o3: bool,
    only_missing: bool,
    validate: bool,
) -> None:
    case = _case_payload(data, base, n_final)
    case["process"] = process
    case["status"] = "running"
    case["updated_at"] = _now()
    started = time.monotonic()
    process_id = BASE_PROCESS_IDS.get(base.key, 0)
    print(
        f"[matrix] id={process_id} n={n_final} {base.key}: {process}",
        flush=True,
    )

    amplicol_settings = _amplicol_matrix_settings(
        amplicol_points=amplicol_points,
        jobs=jobs,
    )
    if not skip_amplicol and _should_run_mode(
        case,
        "amplicol",
        only_missing=only_missing,
        expected_settings=amplicol_settings,
    ):
        case["amplicol"] = _run_amplicol_case(
            process,
            base=base,
            deadline_started=started,
            time_limit=None,
            amplicol_points=amplicol_points,
            jobs=jobs,
            matrix_settings=amplicol_settings,
        )
    try:
        reference_color_order = _reference_color_order_for_case(case)
        selected_lc_sector_ids = _selected_lc_sector_ids_for_case(
            process,
            base,
            case,
        )
        sector_selection_error: Exception | None = None
    except Exception as exc:  # noqa: BLE001 - render as localized setup failure.
        reference_color_order = None
        selected_lc_sector_ids = None
        sector_selection_error = exc
    jit_settings = _pyamplicol_matrix_settings(
        backend_key="jit",
        target_runtime=target_runtime,
        n_cores=n_cores,
        selected_lc_sector_ids=selected_lc_sector_ids,
        reference_color_order=reference_color_order,
        base=base,
    )
    if sector_selection_error is not None and _ok(_mode(case, "amplicol")):
        sector_error_payload = _error_payload(
            RuntimeError(f"LC sector selection failed: {sector_selection_error}")
        )
        if not skip_jit and _should_run_mode(
            case,
            "pyamplicol_jit",
            only_missing=only_missing,
            expected_settings=jit_settings,
        ):
            case["pyamplicol_jit"] = dict(
                sector_error_payload,
                mode="pyAmpliCol - staged-DAG | JIT",
                matrix_settings=jit_settings,
            )
        if not skip_cpp_o3:
            cpp_o3_settings = _pyamplicol_matrix_settings(
                backend_key="cpp_o3",
                target_runtime=target_runtime,
                n_cores=n_cores,
                selected_lc_sector_ids=selected_lc_sector_ids,
                reference_color_order=reference_color_order,
                base=base,
            )
            if _should_run_mode(
                case,
                "pyamplicol_cpp_o3",
                only_missing=only_missing,
                expected_settings=cpp_o3_settings,
            ):
                case["pyamplicol_cpp_o3"] = dict(
                    sector_error_payload,
                    mode="pyAmpliCol - staged-DAG | C++ O3",
                    matrix_settings=cpp_o3_settings,
                )
        case["status"] = "done"
        case["updated_at"] = _now()
        return

    if not skip_jit and _should_run_mode(
        case,
        "pyamplicol_jit",
        only_missing=only_missing,
        expected_settings=jit_settings,
    ):
        case["pyamplicol_jit"] = _run_pyamplicol_case(
            process,
            base=base,
            n_final=n_final,
            backend_key="jit",
            output_dir=output_root / base.key / f"n{n_final}" / "jit",
            deadline_started=time.monotonic(),
            time_limit=None,
            target_runtime=target_runtime,
            n_cores=n_cores,
            selected_lc_sector_ids=selected_lc_sector_ids,
            reference_color_order=reference_color_order,
            matrix_settings=jit_settings,
        )
    cpp_o3_settings = _pyamplicol_matrix_settings(
        backend_key="cpp_o3",
        target_runtime=target_runtime,
        n_cores=n_cores,
        selected_lc_sector_ids=selected_lc_sector_ids,
        reference_color_order=reference_color_order,
        base=base,
    )
    if not skip_cpp_o3 and _should_run_mode(
        case,
        "pyamplicol_cpp_o3",
        only_missing=only_missing,
        expected_settings=cpp_o3_settings,
    ):
        case["pyamplicol_cpp_o3"] = _run_pyamplicol_case(
            process,
            base=base,
            n_final=n_final,
            backend_key="cpp_o3",
            output_dir=output_root / base.key / f"n{n_final}" / "cpp_o3",
            deadline_started=time.monotonic(),
            time_limit=time_limit,
            target_runtime=target_runtime,
            n_cores=n_cores,
            selected_lc_sector_ids=selected_lc_sector_ids,
            reference_color_order=reference_color_order,
            matrix_settings=cpp_o3_settings,
        )
    if validate and not skip_amplicol and _ok(_mode(case, "amplicol")):
        case["validation"] = _run_validation_case(
            process,
            base=base,
            case=case,
            deadline_started=started,
            time_limit=None,
            jobs=jobs,
        )
    elif not validate:
        case.pop("validation", None)
    case["status"] = "done"
    case["updated_at"] = _now()


def _run_amplicol_case(
    process: str,
    *,
    base: BaseProcess,
    deadline_started: float,
    time_limit: float | None,
    amplicol_points: int,
    jobs: int,
    matrix_settings: dict[str, Any],
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "mode": "AmpliCol - using-library",
        "status": "running",
        "started_at": _now(),
    }
    try:
        remaining = _remaining(deadline_started, time_limit)
        adapter = AmplicolAdapter(REPO_ROOT, jobs=jobs, timeout=remaining)
        try:
            point: tuple[ExternalMomentum, ...] | None = generic_validation_point(
                process
            )
        except Exception:  # noqa: BLE001 - benchmark harness records failures later.
            point = None
        process_list_backend = "python"
        build = adapter.prepare_library(
            process,
            options=base.process_options(),
            process_list_backend=process_list_backend,
            warmup_particles=point,
        )
        remaining = _remaining(deadline_started, time_limit)
        adapter = AmplicolAdapter(REPO_ROOT, jobs=jobs, timeout=remaining)
        run = adapter.run_library_benchmark(
            process,
            points=max(1, amplicol_points),
            process_file=build.process_file,
            options=base.process_options(),
            process_list_backend=process_list_backend,
        )
        reference_probe = "direct_generated_library_benchmark"
        runtime_s = _library_runtime_per_point(run, max(1, amplicol_points))
        raw_color_order = _first_process_file_color_order(build.process_file)
        mapped_color_order = _map_process_file_color_order_to_py_labels(
            process,
            build.process_file,
            raw_color_order,
            options=base.process_options(),
        )
        payload.update(
            {
                "status": "ok",
                "generation_s": build.total_command_time_s,
                "reference_probe": reference_probe,
                "process_file": str(build.process_file),
                "process_list_backend": (
                    process_list_backend
                ),
                "reference_color_order": mapped_color_order,
                "reference_color_order_process_file": raw_color_order,
                "runtime_us_per_point": (
                    None if runtime_s is None else 1.0e6 * runtime_s
                ),
                "commands": [
                    {
                        "args": command.args,
                        "elapsed_s": command.elapsed_s,
                        "returncode": command.returncode,
                    }
                    for command in (*build.commands, *run.commands)
                ],
                "commands_s": [
                    command.elapsed_s for command in (*build.commands, *run.commands)
                ],
                "timing_rows": [
                    {
                        "label": row.label,
                        "seconds": row.seconds,
                        "note": row.note,
                    }
                    for row in run.timing_rows
                ],
                "matrix_settings": matrix_settings,
                "finished_at": _now(),
            }
        )
    except Exception as exc:  # noqa: BLE001 - benchmark harness records failures.
        payload.update(_error_payload(exc))
        if _amplicol_error_is_unsupported(str(exc)):
            payload["status"] = "unsupported"
    return payload


def _first_process_file_color_order(process_file: str | Path) -> list[int] | None:
    first_entry = next(iter(amplicol_process_file_integrals(process_file)), None)
    if first_entry is None:
        return None
    group, integral = first_entry
    entry = amplicol_process_file_entry(process_file, group=group, integral=integral)
    if entry is None:
        return None
    return [int(label) for label in entry["color_order"]]


def _first_process_file_pdgs(process_file: str | Path) -> tuple[int, ...] | None:
    first_entry = next(iter(amplicol_process_file_integrals(process_file)), None)
    if first_entry is None:
        return None
    group, integral = first_entry
    entry = amplicol_process_file_entry(process_file, group=group, integral=integral)
    if entry is None:
        return None
    return tuple(int(pdg) for pdg in entry["process"])


def _map_process_file_color_order_to_py_labels(
    process: str,
    process_file: str | Path,
    color_order: list[int] | None,
    *,
    options: ProcessOptions,
) -> list[int] | None:
    if not color_order:
        return color_order
    process_file_pdgs = _first_process_file_pdgs(process_file)
    if process_file_pdgs is None:
        return color_order
    process_ir = build_process_ir(process, options=options)
    py_final_legs_by_pdg: dict[int, list[int]] = {}
    for leg in process_ir.final_legs:
        if leg.outgoing_pdg is None:
            continue
        py_final_legs_by_pdg.setdefault(int(leg.outgoing_pdg), []).append(int(leg.label))
    for labels in py_final_legs_by_pdg.values():
        labels.sort()
    used_by_pdg: dict[int, int] = {}
    process_file_to_py: dict[int, int] = {1: 1, 2: 2}
    for process_file_label, pdg in enumerate(process_file_pdgs, start=1):
        if process_file_label <= 2:
            continue
        candidates = py_final_legs_by_pdg.get(int(pdg), [])
        used = used_by_pdg.get(int(pdg), 0)
        if used >= len(candidates):
            return color_order
        process_file_to_py[process_file_label] = candidates[used]
        used_by_pdg[int(pdg)] = used + 1
    try:
        return [process_file_to_py[int(label)] for label in color_order]
    except KeyError:
        return color_order


def _selected_lc_sector_ids_for_case(
    process: str,
    base: BaseProcess,
    case: dict[str, Any],
) -> set[int] | None:
    reference_color_order = _reference_color_order_for_case(case)
    if reference_color_order is None:
        if base.lc_sector_strategy != "line-pairing-representatives":
            return None
        color_plan = build_color_plan(
            process,
            color_accuracy="lc",
            options=base.process_options(),
            max_sectors=(
                20000 if base.max_color_sectors is None else base.max_color_sectors
            ),
        )
        return set(lc_line_pairing_representative_ids(color_plan))
    color_plan = build_color_plan(
        process,
        color_accuracy="lc",
        options=base.process_options(),
        max_sectors=(
            20000 if base.max_color_sectors is None else base.max_color_sectors
        ),
        reference_color_order=reference_color_order,
    )
    return select_leading_color_sector_ids_from_plan(
        color_plan,
        reference_color_order=reference_color_order,
    )


def _reference_color_order_for_case(case: dict[str, Any]) -> tuple[int, ...] | None:
    amplicol = _mode(case, "amplicol")
    reference_color_order = amplicol.get("reference_color_order")
    if not isinstance(reference_color_order, list) or not reference_color_order:
        return None
    return tuple(int(label) for label in reference_color_order)


def _run_validation_case(
    process: str,
    *,
    base: BaseProcess,
    case: dict[str, Any],
    deadline_started: float,
    time_limit: float | None,
    jobs: int,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "status": "running",
        "started_at": _now(),
        "kind": "same-deterministic-point",
    }
    try:
        amplicol = _mode(case, "amplicol")
        process_file = amplicol.get("process_file")
        if not isinstance(process_file, str) or not process_file:
            raise RuntimeError("AmpliCol process_file missing from matrix payload")
        backend = str(amplicol.get("process_list_backend", "legacy"))
        if backend not in {"python", "legacy"}:
            backend = "legacy"
        point, point_source = _validation_particles_for_case(process, case)
        adapter = AmplicolAdapter(
            REPO_ROOT,
            jobs=jobs,
            timeout=_remaining(deadline_started, time_limit),
        )
        point, point_order_source = _reorder_validation_particles_for_fortran(
            adapter,
            process,
            particles=point,
            process_file=process_file,
            options=base.process_options(),
            backend=backend,
        )
        run = adapter.run_amplicol_momenta_probe(
            process,
            particles=point,
            points=1,
            process_file=process_file,
            options=base.process_options(),
            timing_sample=1,
            quiet=False,
            use_library=True,
            process_list_backend=backend,  # type: ignore[arg-type]
        )
        if not run.probe_points:
            raise RuntimeError("AmpliCol validation probe did not report a value")
        reference = float(run.probe_points[0].matrix_element)
        values: dict[str, float] = {}
        for key in ("pyamplicol_jit", "pyamplicol_cpp_o3"):
            value = _first_payload_value(_mode(case, key))
            if value is not None:
                values[key] = value
        if not values:
            payload.update(
                {
                    "status": "not_available",
                    "reference": reference,
                    "values": values,
                    "relative_differences": {},
                    "max_relative_difference": None,
                    "tolerance": VALIDATION_REL_TOL,
                    "absolute_tolerance": VALIDATION_ABS_TOL,
                    "reason": "No pyAmpliCol validation values available",
                    "amplicol_command": run.commands[0].args if run.commands else (),
                    "finished_at": _now(),
                }
            )
            return payload
        rel_diffs = {
            key: _relative_difference(reference, value)
            for key, value in values.items()
        }
        abs_diffs = {key: abs(reference - value) for key, value in values.items()}
        max_rel_diff = max(rel_diffs.values(), default=0.0)
        max_abs_diff = max(abs_diffs.values(), default=0.0)
        payload.update(
            {
                "status": (
                    "ok"
                    if max_rel_diff <= VALIDATION_REL_TOL
                    else "failed"
                ),
                "reference": reference,
                "values": values,
                "relative_differences": rel_diffs,
                "absolute_differences": abs_diffs,
                "max_relative_difference": max_rel_diff,
                "max_absolute_difference": max_abs_diff,
                "tolerance": VALIDATION_REL_TOL,
                "absolute_tolerance": VALIDATION_ABS_TOL,
                "point_source": point_source,
                "point_order_source": point_order_source,
                "amplicol_command": run.commands[0].args if run.commands else (),
                "finished_at": _now(),
            }
        )
    except Exception as exc:  # noqa: BLE001 - benchmark harness records failures.
        error_payload = _error_payload(exc)
        payload.update(error_payload)
        backend_limitation = _documented_backend_limitation_from_error(
            str(error_payload.get("error", "")),
            backend_key=backend_key,
        )
        if backend_limitation is not None:
            payload.update(
                {
                    "status": "backend_unsupported",
                    "backend_limitation": backend_limitation,
                    "note": (
                        "SymJIT 2.19.2 still aborts while materializing this "
                        "large JIT evaluator on AArch64; the entry is rendered "
                        "as a localized backend N/A rather than as a physics "
                        "validation or generated-library comparison failure."
                    ),
                }
            )
    return payload


def _validation_particles_for_case(
    process: str,
    case: dict[str, Any],
) -> tuple[tuple[ExternalMomentum, ...], str]:
    for key in ("pyamplicol_cpp_o3", "pyamplicol_jit"):
        mode = _mode(case, key)
        output_dir = mode.get("output_dir")
        if not isinstance(output_dir, str) or not output_dir:
            continue
        path = Path(output_dir) / "validation_momenta.json"
        if not path.exists():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("available") is False:
                continue
            points = payload.get("points")
            if not isinstance(points, list) or not points:
                continue
            first_point = points[0]
            if not isinstance(first_point, list):
                continue
            particles: list[ExternalMomentum] = []
            for item in first_point:
                if not isinstance(item, dict):
                    raise ValueError("validation particle entry is not an object")
                momentum = item.get("momentum")
                if not isinstance(momentum, list) or len(momentum) != 4:
                    raise ValueError("validation particle momentum must have 4 entries")
                particles.append(
                    ExternalMomentum(
                        pdg=int(item["pdg"]),
                        momentum=tuple(float(component) for component in momentum),
                    )
                )
            if particles:
                return tuple(particles), str(path)
        except Exception as exc:  # noqa: BLE001 - fall back to generated point.
            print(f"[matrix] warning: could not load {path}: {exc}", flush=True)
    return generic_validation_point(process), "generic_validation_point"


def _reorder_validation_particles_for_fortran(
    adapter: AmplicolAdapter,
    process: str,
    *,
    particles: tuple[ExternalMomentum, ...],
    process_file: str,
    options: ProcessOptions,
    backend: str,
) -> tuple[tuple[ExternalMomentum, ...], str]:
    first_entry = next(iter(amplicol_process_file_integrals(process_file)), None)
    if first_entry is not None:
        group, integral = first_entry
        entry = amplicol_process_file_entry(
            process_file,
            group=group,
            integral=integral,
        )
        if entry is not None:
            expected_pdgs = tuple(int(pdg) for pdg in entry["process"])
            return _reorder_particles_by_pdg(particles, expected_pdgs), (
                "process-file-pdg-order:" + ",".join(str(pdg) for pdg in expected_pdgs)
            )
    try:
        if process.strip().lower() == "d d~ > z":
            order_run = adapter.run_amplicol_fixed_probe(
                process,
                points=1,
                process_file=process_file,
                options=options,
                timing_sample=1,
                quiet=False,
                use_library=True,
                process_list_backend=backend,  # type: ignore[arg-type]
            )
        else:
            order_run = adapter.run_amplicol_probe(
                process,
                points=1,
                process_file=process_file,
                options=options,
                timing_sample=1,
                quiet=False,
                use_library=True,
                process_list_backend=backend,  # type: ignore[arg-type]
            )
        if not order_run.probe_points:
            return particles, "process-order"
        expected_pdgs = tuple(
            int(particle.pdg) for particle in order_run.probe_points[0].particles
        )
        return _reorder_particles_by_pdg(particles, expected_pdgs), (
            "fortran-probe-pdg-order:" + ",".join(str(pdg) for pdg in expected_pdgs)
        )
    except Exception as exc:  # noqa: BLE001 - validation can still try process order.
        print(
            f"[matrix] warning: could not determine Fortran external order for "
            f"{process}: {exc}",
            flush=True,
        )
        return particles, "process-order-fallback"


def _reorder_particles_by_pdg(
    particles: Sequence[ExternalMomentum],
    expected_pdgs: Sequence[int],
) -> tuple[ExternalMomentum, ...]:
    remaining = list(particles)
    ordered: list[ExternalMomentum] = []
    for expected in expected_pdgs:
        for index, particle in enumerate(remaining):
            if int(particle.pdg) == int(expected):
                ordered.append(particle)
                del remaining[index]
                break
        else:
            raise ValueError(
                f"validation point cannot be reordered to Fortran PDG order "
                f"{list(expected_pdgs)}"
            )
    if remaining:
        raise ValueError("validation point has extra particles after PDG reordering")
    return tuple(ordered)


def _run_pyamplicol_case(
    process: str,
    *,
    base: BaseProcess,
    n_final: int,
    backend_key: str,
    output_dir: Path,
    deadline_started: float,
    time_limit: float | None,
    target_runtime: float,
    n_cores: int,
    selected_lc_sector_ids: set[int] | None,
    reference_color_order: Sequence[int] | None,
    matrix_settings: dict[str, Any],
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "mode": (
            "pyAmpliCol - staged-DAG | JIT"
            if backend_key == "jit"
            else "pyAmpliCol - staged-DAG | C++ O3"
        ),
        "status": "running",
        "started_at": _now(),
        "output_dir": str(output_dir),
        "matrix_settings": matrix_settings,
    }
    generate = [
        sys.executable,
        "-m",
        "pyamplicol",
        "generate-process",
        *base.cli_flags(),
        "--replace",
        "--n_cores",
        str(n_cores),
        "--symbolica-n-cores",
        str(n_cores),
        "--json",
    ]
    if selected_lc_sector_ids:
        generate.extend(
            [
                "--lc-sector-ids",
                ",".join(str(sector_id) for sector_id in sorted(selected_lc_sector_ids)),
            ]
        )
    if reference_color_order:
        generate.extend(
            [
                "--reference-color-order",
                ",".join(str(label) for label in reference_color_order),
            ]
        )
    if backend_key == "jit":
        generate.extend(["--symbolica-evaluator-backend", "jit"])
    else:
        generate.extend(
            [
                "--symbolica-evaluator-backend",
                "compiled-complex",
                "--symbolica-compiled-preset",
                "runtime-o3",
                "--symbolica-compiled-chunk-compile-workers",
                str(max(1, n_cores)),
            ]
        )
    generate.extend([process, str(output_dir)])

    try:
        gen = _run_json_command(
            generate,
            timeout=_remaining(deadline_started, time_limit),
        )
        timed = _run_json_command(
            [
                sys.executable,
                "-m",
                "pyamplicol",
                "time-process",
                "--target-runtime",
                str(target_runtime),
                "--json",
                str(output_dir),
            ],
            timeout=_remaining(deadline_started, time_limit),
        )
        profile = timed.get("profile", {})
        payload.update(
            {
                "status": "ok",
                "generation_s": _optional_float(
                    gen.get("_command_elapsed_s", gen.get("generation_s"))
                ),
                "internal_generation_s": _optional_float(gen.get("generation_s")),
                "runtime_us_per_point": _optional_float(
                    profile.get("core_evaluator_us_per_point")
                ),
                "wall_us_per_point": _optional_float(
                    profile.get("wall_us_per_point")
                ),
                "samples": profile.get("samples"),
                "selected_lc_sector_ids": (
                    None
                    if selected_lc_sector_ids is None
                    else sorted(selected_lc_sector_ids)
                ),
                "reference_color_order": (
                    None
                    if reference_color_order is None
                    else [int(label) for label in reference_color_order]
                ),
                "generate_payload": _compact_payload(gen),
                "time_payload": _compact_payload(timed),
                "matrix_settings": matrix_settings,
                "finished_at": _now(),
            }
        )
    except Exception as exc:  # noqa: BLE001 - benchmark harness records failures.
        payload.update(_error_payload(exc))
        limitation = _documented_backend_limitation_from_error(
            str(payload.get("error", "")),
            backend_key=backend_key,
        )
        if limitation is not None:
            payload["status"] = "backend_unsupported"
            payload["backend_limitation"] = limitation
    return payload


def _run_json_command(args: Sequence[str], *, timeout: float | None) -> dict[str, Any]:
    env = dict(os.environ)
    env["PYTHONPATH"] = (
        f"{SRC_DIR}{os.pathsep}{env['PYTHONPATH']}"
        if env.get("PYTHONPATH")
        else str(SRC_DIR)
    )
    start = time.perf_counter()
    completed = subprocess.run(
        list(args),
        cwd=REPO_ROOT,
        env=env,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    elapsed_s = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{completed.returncode}: {' '.join(args)}\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    payload = _parse_json_output(completed.stdout)
    payload["_command_elapsed_s"] = elapsed_s
    payload["_command_args"] = list(args)
    return payload


def render_latex_table(
    data: dict[str, Any],
    n_values: Sequence[int],
    *,
    columns_per_table: int = 3,
) -> str:
    shown = tuple(sorted(dict.fromkeys(int(n) for n in n_values)))
    chunks = tuple(_chunks(shown, max(1, int(columns_per_table))))
    lines = [
        "% Generated by docs/result_matrix.py; edit result_matrix_data.json instead.",
        r"\providecommand{\matrixpunct}[1]{\textcolor{black}{\texttt{#1}}}",
        r"\providecommand{\matrixratio}[2]{\matrixpunct{(}\textcolor{#1}{\texttt{x#2}}\matrixpunct{)}}",
        r"\providecommand{\matrixratioinner}[2]{\textcolor{#1}{\texttt{#2}}}",
        r"\providecommand{\matrixna}[1]{\textcolor{#1}{\texttt{N/A}}}",
        r"\providecommand{\matrixnaratio}[1]{\matrixpunct{(}\matrixna{#1}\matrixpunct{)}}",
        (
            r"\providecommand{\matrixratiopair}[4]{"
            r"\matrixpunct{(}"
            r"\matrixratioinner{#1}{#2}"
            r"\matrixpunct{|}"
            r"\matrixratioinner{#3}{#4}"
            r"\matrixpunct{)}}"
        ),
        (
            r"\providecommand{\matrixslot}[2]{\makebox[#1][l]{#2}}"
        ),
        (
            r"\providecommand{\matrixcell}[6]{"
            r"\begin{tabular}[t]{@{}l@{\hspace{0.025in}}l@{\hspace{0.025in}}l@{}}"
            r"\matrixslot{0.62in}{#1}&\matrixslot{0.90in}{#2}&"
            r"\matrixslot{0.90in}{#3}\\"
            r"\matrixslot{0.62in}{#4}&\matrixslot{0.90in}{#5}&"
            r"\matrixslot{0.90in}{#6}"
            r"\end{tabular}}"
        ),
    ]
    entries = data.get("entries", {})
    if not isinstance(entries, dict):
        entries = {}
    validation_summary = _validation_summary(entries, shown)
    for chunk_index, chunk in enumerate(chunks):
        multiplicity_columns = (
            r"@{\hspace{0.055in}}".join("L{2.56in}" for _ in chunk)
        )
        colspec = (
            r"@{}r@{\hspace{0.055in}}L{1.42in}@{\hspace{0.075in}}"
            + multiplicity_columns
            + r"@{}"
        )
        lines.extend(
            [
                r"\begin{landscape}",
                (
                    r"\section{Generic LC Performance Matrix}"
                    if chunk_index == 0
                    else r"\subsection*{Generic LC Performance Matrix (continued)}"
                ),
                r"\begingroup",
                r"\scriptsize",
                r"\setlength{\tabcolsep}{2.2pt}",
                r"\renewcommand{\arraystretch}{1.23}",
            ]
        )
        if chunk_index == 0:
            lines.extend(
                [
                    r"\noindent\footnotesize Each cell compares generated-library \AC\ against "
                    r"\PAC\ staged-DAG JIT and staged-DAG C++ O3 for the same final-state "
                    r"multiplicity \(n\).  In each cell, the upper line is generation time; the lower "
                    r"line is wall runtime per phase-space point.  The left slot is the "
                    r"generated-library \AC\ reference value, the middle slot is \PAC\ JIT, "
                    r"and the right slot is \PAC\ C++ O3. Runtime multipliers are shown as "
                    r"\texttt{(wall|core)}. "
                    r"When no \AC\ reference exists but \PAC\ ran, the \PAC\ slots show "
                    r"absolute generation time and parenthesized \texttt{wall|core} "
                    r"runtime in microseconds. "
                    r"Colored "
                    r"ratios are \PAC/\AC: green below one, orange below two, red otherwise.  "
                    r"When rows are generated with \texttt{--validate}, \PAC\ JIT and "
                    r"C++ O3 values are compared to the same-point generated-library "
                    r"\AC\ value with a relative tolerance of \(10^{-8}\).",
                    r"\par\smallskip",
                    _validation_summary_latex(validation_summary),
                    r"\par\smallskip",
                ]
            )
        lines.extend(
            [
                r"\begin{longtable}{" + colspec + "}",
                r"\toprule",
                r"\textbf{ID} & \textbf{base process}"
                + "".join(f" & \\textbf{{n={n}}}" for n in chunk)
                + r" \\",
                r"\specialrule{0.85pt}{0pt}{0pt}",
                r"\endfirsthead",
                r"\toprule",
                r"\textbf{ID} & \textbf{base process}"
                + "".join(f" & \\textbf{{n={n}}}" for n in chunk)
                + r" \\",
                r"\specialrule{0.85pt}{0pt}{0pt}",
                r"\endhead",
            ]
        )
        for row_index, base in enumerate(BASE_PROCESSES):
            row_entries = entries.get(base.key, {})
            if not isinstance(row_entries, dict):
                row_entries = {}
            cells = [rf"\texttt{{{row_index + 1}}}", base.label]
            for n_final in chunk:
                if n_final < base.min_final_count:
                    cells.append(_structural_na())
                    continue
                case = row_entries.get(str(n_final), {})
                cells.append(_latex_cell(case if isinstance(case, dict) else {}))
            if row_index % 2 == 0:
                lines.append(r"\rowcolor{refblue}")
            lines.append(" & ".join(cells) + r" \\")
            lines.append(r"\addlinespace[0.22em]")
        lines.extend(
            [
                r"\bottomrule",
                r"\end{longtable}",
                r"\endgroup",
                r"\end{landscape}",
                "",
            ]
        )
    lines.extend(_matrix_run_settings_latex())
    lines.extend(_matrix_status_notes_latex(entries, shown))
    return "\n".join(lines)


def _matrix_run_settings_latex() -> list[str]:
    return [
        r"\paragraph{Matrix run settings.}",
        (
            r"The matrix table is generated with one uniform benchmark setup. "
            r"The full refresh command is:"
        ),
        r"\begin{lstlisting}[language=bash,breaklines=true,breakatwhitespace=true]",
        (
            "pyAmpliCol/dependencies/.venv/bin/python pyAmpliCol/docs/result_matrix.py "
            "--generate-data 1-9 --show-data 1-9 "
            f"--time-limit {DEFAULT_TIME_LIMIT_S:g} "
            f"--target-runtime {DEFAULT_TARGET_RUNTIME_S:g} "
            f"--amplicol-points {DEFAULT_AMPLICOL_POINTS} "
            f"--jobs {DEFAULT_JOBS} --n-cores {DEFAULT_N_CORES} "
            f"--process-workers {DEFAULT_PROCESS_WORKERS} "
            f"--columns-per-table {DEFAULT_COLUMNS_PER_TABLE}"
        ),
        r"\end{lstlisting}",
        (
            r"\noindent\footnotesize Partial refreshes use the same settings with "
            r"\texttt{--base-process}, \texttt{--skip-jit}, \texttt{--skip-cpp-o3}, "
            r"\texttt{--process-ids}, or \texttt{--skip-amplicol} added only to restrict which cells are "
            r"updated.  PyAmpliCol-only refreshes may additionally use "
            r"\texttt{--process-workers 5}; this is intentionally disabled when "
            r"the \AC\ reference is refreshed because that path uses shared "
            r"top-level generated-library files.  \AC\ reference timings use the "
            r"generated-library sequence "
            r"\texttt{amplicol\_generate --library=create}, "
            r"\texttt{make -j4 amplicol\_generate\_library}, and a direct "
            r"generated-library benchmark using \texttt{--library=use}; they do not "
            r"use integration-driver timings.  \PAC\ rows use the generic LC "
            r"staged-DAG artifact, Rusticol double precision, JIT or C++ O3 "
            r"Symbolica evaluators as labelled, and \texttt{time-process} with a "
            r"one-second target runtime.  The \texttt{--time-limit} setting applies "
            r"only to C++ O3 generation; \AC\ and JIT rows are left uncapped."
        ),
        r"\par\smallskip",
    ]


def _matrix_status_notes_latex(
    entries: dict[str, Any],
    n_values: Sequence[int],
) -> list[str]:
    unsupported_four_quark_line = False
    documented_backend_limits: list[str] = []
    non_structural: list[str] = []
    for base in BASE_PROCESSES:
        row_entries = entries.get(base.key, {})
        if not isinstance(row_entries, dict):
            continue
        for n_final in n_values:
            case = row_entries.get(str(n_final), {})
            if not isinstance(case, dict):
                continue
            ref = _mode(case, "amplicol")
            if ref.get("status") == "unsupported" and base.key == "dd_4q_lines":
                unsupported_four_quark_line = True
            for key, label in (
                ("amplicol", r"\AC"),
                ("pyamplicol_jit", r"\PAC\ JIT"),
                ("pyamplicol_cpp_o3", r"\PAC\ C++ O3"),
            ):
                mode = _mode(case, key)
                status = str(mode.get("status", ""))
                if _is_documented_backend_limitation(mode):
                    documented_backend_limits.append(
                        rf"{base.label}, \(n={n_final}\), {label}"
                    )
                    continue
                if status in {"", "ok", "missing", "not_applicable", "unsupported"}:
                    continue
                entry = (
                    rf"{base.label}, \(n={n_final}\), {label}: "
                    rf"\texttt{{{_latex_escape(status)}}}"
                )
                non_structural.append(
                    entry
                )
    lines: list[str] = []
    if unsupported_four_quark_line:
        lines.extend(
            [
                r"\paragraph{Structural reference limitations.}",
                (
                    r"For the four-quark-line row "
                    r"\(d\bar d\to u\bar u\,s\bar s\,c\bar c+(n-6)g\), "
                    r"the \AC\ reference column is a structural N/A: Fortran "
                    r"\AC\ stops with \texttt{Unknown number of quarks and "
                    r"anti-quarks}.  The \PAC\ cells are therefore shown as "
                    r"absolute generation/runtime timings when they are available."
                ),
                r"\par\smallskip",
            ]
        )
    if non_structural:
        shown = "; ".join(non_structural[:8])
        if len(non_structural) > 8:
            shown += rf"; plus {len(non_structural) - 8} more localized entries"
        lines.extend(
            [
                r"\paragraph{Current localized incomplete entries.}",
                (
                    r"The following entries are not structural N/A and should be "
                    r"resolved or explicitly documented before the table is treated "
                    r"as final: "
                    + shown
                    + "."
                ),
                r"\par\smallskip",
            ]
        )
    if documented_backend_limits:
        shown = "; ".join(documented_backend_limits[:8])
        if len(documented_backend_limits) > 8:
            shown += (
                rf"; plus {len(documented_backend_limits) - 8} more localized entries"
            )
        lines.extend(
            [
                r"\paragraph{Documented JIT backend limitations.}",
                (
                    r"The following \PAC\ JIT entries are rendered as N/A because "
                    r"they hit a localized Symbolica/SymJIT AArch64 complex-JIT "
                    r"materialization limitation, not a generated-library "
                    r"comparison or physics-validation failure: "
                    + shown
                    + r".  The raw assertion is retained in "
                    r"\texttt{result\_matrix\_data.json}. The "
                    r"current managed dependency set uses SymJIT \texttt{2.19.2}, "
                    r"which fixes the smaller historical \texttt{v220} AArch64 "
                    r"JIT segfault reproduced by \texttt{MRE\_symjit\_bug\_new.py}; "
                    r"the remaining listed AArch64 vector-offset assertions still "
                    r"require an upstream SymJIT fix."
                ),
                r"\par\smallskip",
            ]
        )
    return lines


def _is_symjit_backend_failure(mode: dict[str, Any]) -> bool:
    error = str(mode.get("error", ""))
    return (
        "exit code -11" in error
        or "exit code 139" in error
        or "SymJIT" in error
        or "symjit" in error
    )


def _is_documented_backend_limitation(mode: dict[str, Any]) -> bool:
    return str(mode.get("status", "")).lower() in {
        "backend_unsupported",
        "backend_limitation",
    }


def _documented_backend_limitation_from_error(
    error: str,
    *,
    backend_key: str,
) -> str | None:
    if backend_key != "jit":
        return None
    if (
        "symjit/rust/arm/vector.rs:795" in error
        and "assertion failed: x.abs() < 1048576" in error
    ):
        return "symjit_aarch64_vector_offset_assertion"
    return None


def _validation_summary(
    entries: dict[str, Any],
    n_values: Sequence[int],
) -> dict[str, Any]:
    validated = 0
    unvalidated = 0
    failed: list[str] = []
    max_rel = 0.0
    for base in BASE_PROCESSES:
        row_entries = entries.get(base.key, {})
        if not isinstance(row_entries, dict):
            row_entries = {}
        for n_final in n_values:
            if n_final < base.min_final_count:
                continue
            case = row_entries.get(str(n_final), {})
            if not isinstance(case, dict):
                continue
            if not _ok(_mode(case, "amplicol")):
                continue
            if not (_ok(_mode(case, "pyamplicol_jit")) or _ok(_mode(case, "pyamplicol_cpp_o3"))):
                continue
            validation = case.get("validation")
            if not isinstance(validation, dict):
                unvalidated += 1
                continue
            tolerance = _optional_float(validation.get("tolerance")) or 1.0e-8
            rel = _optional_float(validation.get("max_relative_difference"))
            if rel is not None:
                max_rel = max(max_rel, rel)
            if validation.get("status") == "ok" and (rel is None or rel <= tolerance):
                validated += 1
                continue
            failed.append(
                f"{base.key}, n={n_final}: status={validation.get('status')}, "
                f"max_rel={rel}, tol={tolerance}"
            )
    return {
        "validated": validated,
        "unvalidated": unvalidated,
        "failed": failed,
        "max_relative_difference": max_rel,
    }


def _validation_summary_latex(summary: dict[str, Any]) -> str:
    return (
        r"\noindent\footnotesize A red \textbf{VALIDATION FAILED} marker is shown "
        r"next to any validated entry whose relative difference exceeds \(10^{-8}\)."
    )


def _chunks(values: Sequence[int], size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(values[index : index + size]) for index in range(0, len(values), size))


def _latex_cell(case: dict[str, Any]) -> str:
    if not case:
        return _missing_na()
    if case.get("status") == "not_applicable":
        return _structural_na()
    ref = _mode(case, "amplicol")
    jit = _mode(case, "pyamplicol_jit")
    cpp = _mode(case, "pyamplicol_cpp_o3")
    if _ok(ref):
        gen = _metric_parts(
            ref,
            jit,
            cpp,
            metric="generation_s",
            formatter=_format_seconds,
        )
        run = _runtime_parts(ref, jit, cpp)
    else:
        gen = (
            _reference_metric_or_status(
                ref,
                metric="generation_s",
                formatter=_format_seconds,
            ),
            _absolute_generation(jit),
            _absolute_generation(cpp),
        )
        run = (
            _reference_metric_or_status(
                ref,
                metric="runtime_us_per_point",
                formatter=_format_us,
            ),
            _absolute_runtime_pair(jit),
            _absolute_runtime_pair(cpp),
        )
    cell = (
        rf"\matrixcell{{{gen[0]}}}{{{gen[1]}}}{{{gen[2]}}}"
        rf"{{{run[0]}}}{{{run[1]}}}{{{run[2]}}}"
    )
    marker = _validation_cell_marker(case)
    if marker:
        return rf"\begin{{tabular}}[t]{{@{{}}l@{{}}}}{cell}\\[-0.1em]{marker}\end{{tabular}}"
    return cell


def _latex_cell_without_reference(jit: dict[str, Any], cpp: dict[str, Any]) -> str:
    if not _ok(jit) and not _ok(cpp):
        if jit:
            return _latex_failure(jit)
        if cpp:
            return _latex_failure(cpp)
        return _missing_na()
    gen = (
        _missing_na(),
        _absolute_generation(jit),
        _absolute_generation(cpp),
    )
    run = (
        _missing_na(),
        _absolute_runtime_pair(jit),
        _absolute_runtime_pair(cpp),
    )
    return (
        rf"\matrixcell{{{gen[0]}}}{{{gen[1]}}}{{{gen[2]}}}"
        rf"{{{run[0]}}}{{{run[1]}}}{{{run[2]}}}"
    )


def _reference_metric_or_status(
    mode: dict[str, Any],
    *,
    metric: str,
    formatter,
) -> str:
    if _ok(mode):
        value = _optional_float(mode.get(metric))
        return _missing_na() if value is None else formatter(value)
    if mode:
        return _latex_failure(mode)
    return _missing_na()


def _absolute_generation(mode: dict[str, Any]) -> str:
    if _is_documented_backend_limitation(mode):
        return _structural_na()
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _latex_failure(mode)
    if not _ok(mode):
        return _missing_na()
    value = _optional_float(mode.get("generation_s"))
    return _missing_na() if value is None else _format_seconds(value)


def _absolute_runtime_pair(mode: dict[str, Any]) -> str:
    if _is_documented_backend_limitation(mode):
        return _structural_na()
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _missing_na()
    if not _ok(mode):
        return _missing_na()
    wall = _optional_float(mode.get("wall_us_per_point"))
    core = _optional_float(mode.get("runtime_us_per_point"))
    if wall is None or core is None:
        return _missing_na()
    return (
        rf"\matrixpunct{{(}}\texttt{{{_format_us_value(wall)}}}"
        rf"\matrixpunct{{|}}\texttt{{{_format_us_value(core)}}}"
        rf"\matrixpunct{{)}}"
    )


def _validation_cell_marker(case: dict[str, Any]) -> str:
    validation = case.get("validation")
    if not isinstance(validation, dict):
        return ""
    tolerance = _optional_float(validation.get("tolerance")) or 1.0e-8
    rel = _optional_float(validation.get("max_relative_difference"))
    if validation.get("status") == "ok" and (rel is None or rel <= tolerance):
        return ""
    rel_text = "N/A" if rel is None else _format_sig(rel, unit="")
    return rf"\textcolor{{speedred}}{{\scriptsize\textbf{{VALIDATION FAILED}} \texttt{{rel={rel_text}}}}}"


def _metric_parts(
    ref: dict[str, Any],
    jit: dict[str, Any],
    cpp: dict[str, Any],
    *,
    metric: str,
    formatter,
) -> tuple[str, str, str]:
    ref_value = _optional_float(ref.get(metric))
    if ref_value is None or ref_value <= 0.0:
        return (_missing_na(), _missing_ratio(), _missing_ratio())
    parts = [formatter(ref_value)]
    for mode in (jit, cpp):
        if _is_documented_backend_limitation(mode):
            parts.append(_missing_ratio())
            continue
        if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
            parts.append(_latex_failure(mode))
            continue
        if not _ok(mode):
            parts.append(_missing_ratio())
            continue
        value = _optional_float(mode.get(metric))
        if value is None:
            parts.append(_missing_ratio())
            continue
        ratio = value / ref_value
        parts.append(_ratio_latex(ratio))
    return (parts[0], parts[1], parts[2])


def _runtime_parts(
    ref: dict[str, Any],
    jit: dict[str, Any],
    cpp: dict[str, Any],
) -> tuple[str, str, str]:
    ref_value = _optional_float(ref.get("runtime_us_per_point"))
    if ref_value is None or ref_value <= 0.0:
        return (_missing_na(), _missing_runtime_pair(), _missing_runtime_pair())
    return (
        _format_us(ref_value),
        _runtime_ratio_pair(jit, ref_value),
        _runtime_ratio_pair(cpp, ref_value),
    )


def _format_us_value(value: float) -> str:
    return _format_compact_number(value)


def _format_compact_number(value: float) -> str:
    text = f"{value:.3g}"
    text = text.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")
    return text


def _runtime_ratio_pair(mode: dict[str, Any], ref_value: float) -> str:
    if _is_documented_backend_limitation(mode):
        return _missing_runtime_pair()
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _missing_runtime_pair()
    if not _ok(mode):
        return _missing_runtime_pair()
    wall = _optional_float(mode.get("wall_us_per_point"))
    core = _optional_float(mode.get("runtime_us_per_point"))
    return _ratio_pair_latex(
        None if wall is None else wall / ref_value,
        None if core is None else core / ref_value,
    )


def _ratio_latex(value: float) -> str:
    color = "speedgreen" if value < 1.0 else "speedorange" if value < 2.0 else "speedred"
    return rf"\matrixratio{{{color}}}{{{_format_ratio(value)}}}"


def _ratio_pair_latex(wall: float | None, core: float | None) -> str:
    wall_color, wall_text = _ratio_fragment(wall)
    core_color, core_text = _ratio_fragment(core)
    return rf"\matrixratiopair{{{wall_color}}}{{{wall_text}}}{{{core_color}}}{{{core_text}}}"


def _ratio_fragment(value: float | None) -> tuple[str, str]:
    if value is None:
        return ("speedred", "N/A")
    color = "speedgreen" if value < 1.0 else "speedorange" if value < 2.0 else "speedred"
    return (color, f"x{_format_ratio(value)}")


def _latex_failure(mode: dict[str, Any]) -> str:
    if not mode:
        return _missing_na()
    error = str(mode.get("error", ""))
    if "process_list.py did not produce processes.txt" in error:
        return _missing_na()
    if _is_documented_backend_limitation(mode):
        return _structural_na()
    if mode.get("status") == "unsupported" or _amplicol_error_is_unsupported(error):
        return _structural_na()
    status = str(mode.get("status", "missing"))
    if status == "missing":
        return _missing_na()
    return rf"\textcolor{{speedred}}{{\texttt{{{_latex_escape(status)}}}}}"


def _amplicol_error_is_unsupported(error: str) -> bool:
    unsupported_markers = (
        "Unknown number of quarks and anti-quarks",
    )
    return any(marker in error for marker in unsupported_markers)


def _structural_na() -> str:
    return r"\matrixna{black!45}"


def _missing_na() -> str:
    return r"\matrixna{speedred}"


def _missing_ratio() -> str:
    return r"\matrixnaratio{speedred}"


def _missing_runtime_pair() -> str:
    return _ratio_pair_latex(None, None)


def _mode(case: dict[str, Any], key: str) -> dict[str, Any]:
    value = case.get(key)
    return value if isinstance(value, dict) else {}


def _ok(mode: dict[str, Any]) -> bool:
    return mode.get("status") == "ok"


def _should_run_mode(
    case: dict[str, Any],
    key: str,
    *,
    only_missing: bool,
    expected_settings: dict[str, Any],
) -> bool:
    if not only_missing:
        return True
    mode = _mode(case, key)
    if _is_documented_backend_limitation(mode):
        return not _settings_match(mode.get("matrix_settings"), expected_settings)
    if not _ok(mode):
        return True
    return not _settings_match(mode.get("matrix_settings"), expected_settings)


def _settings_match(actual: object, expected: dict[str, Any]) -> bool:
    if not isinstance(actual, dict):
        return False
    return actual == expected


def _amplicol_matrix_settings(*, amplicol_points: int, jobs: int) -> dict[str, Any]:
    return {
        "amplicol_points": int(amplicol_points),
        "jobs": int(jobs),
        "process_list_backend": "python",
        "workflow": "library_create_make_library_use",
        "runtime_probe": "direct_generated_library_benchmark",
    }


def _pyamplicol_matrix_settings(
    *,
    backend_key: str,
    target_runtime: float,
    n_cores: int,
    selected_lc_sector_ids: set[int] | None,
    reference_color_order: Sequence[int] | None,
    base: BaseProcess,
) -> dict[str, Any]:
    settings: dict[str, Any] = {
        "runtime": "rusticol",
        "precision": 16,
        "target_runtime_s": float(target_runtime),
        "n_cores": int(n_cores),
        "symbolica_n_cores": int(n_cores),
        "selected_lc_sector_ids": (
            None
            if selected_lc_sector_ids is None
            else sorted(int(sector_id) for sector_id in selected_lc_sector_ids)
        ),
        "reference_color_order": (
            None
            if reference_color_order is None
            else [int(label) for label in reference_color_order]
        ),
        "max_currents": base.max_currents,
        "max_color_sectors": base.max_color_sectors,
        "max_quark_pairs": base.max_quark_pairs,
        "lc_sector_strategy": base.lc_sector_strategy,
    }
    if backend_key == "jit":
        settings.update(
            {
                "symbolica_evaluator_backend": "jit",
            }
        )
    elif backend_key == "cpp_o3":
        settings.update(
            {
                "symbolica_evaluator_backend": "compiled-complex",
                "symbolica_compiled_preset": "runtime-o3",
                "symbolica_compiled_chunk_compile_workers": int(n_cores),
            }
        )
    else:
        settings["symbolica_evaluator_backend"] = backend_key
    return settings


def _record_not_applicable(data: dict[str, Any], base: BaseProcess, n_final: int) -> None:
    case = _case_payload(data, base, n_final)
    case.clear()
    case.update(
        {
            "status": "not_applicable",
            "process": None,
            "updated_at": _now(),
        }
    )


def _case_payload(data: dict[str, Any], base: BaseProcess, n_final: int) -> dict[str, Any]:
    _refresh_data_metadata(data)
    entries = data.setdefault("entries", {})
    if not isinstance(entries, dict):
        raise TypeError("result matrix data 'entries' must be a JSON object")
    base_entries = entries.setdefault(base.key, {})
    if not isinstance(base_entries, dict):
        raise TypeError(f"entry for {base.key!r} must be a JSON object")
    case = base_entries.setdefault(str(n_final), {})
    if not isinstance(case, dict):
        raise TypeError(f"entry for {base.key!r}, n={n_final} must be an object")
    case.setdefault("base_key", base.key)
    case.setdefault("n_final", n_final)
    return case


def _refresh_data_metadata(data: dict[str, Any]) -> None:
    data.setdefault("schema_version", 1)
    data.setdefault("created_by", "pyAmpliCol/docs/result_matrix.py")
    data["updated_at"] = _now()
    data["base_processes"] = [
        base_payload(process_id, item)
        for process_id, item in enumerate(BASE_PROCESSES, start=1)
    ]


def base_payload(process_id: int, base: BaseProcess) -> dict[str, Any]:
    return {
        "process_id": int(process_id),
        "key": base.key,
        "label": base.label,
        "initial": list(base.initial),
        "base_final": list(base.base_final),
        "include_cc": base.include_cc,
        "include_3qqbar": base.include_3qqbar,
        "include_resonance": base.include_resonance,
        "max_quark_pairs": base.max_quark_pairs,
    }


def _load_data(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _write_data(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _parse_n_values(values: Iterable[str]) -> tuple[int, ...]:
    parsed: set[int] = set()
    for raw in values:
        for token in str(raw).replace(",", " ").split():
            if not token:
                continue
            if "-" in token:
                start_s, end_s = token.split("-", 1)
                start = int(start_s)
                end = int(end_s)
                step = 1 if end >= start else -1
                parsed.update(range(start, end + step, step))
            else:
                parsed.add(int(token))
    return tuple(sorted(value for value in parsed if value > 0))


def _parse_process_ids(values: Iterable[str]) -> tuple[int, ...]:
    parsed = _parse_n_values(values)
    max_id = len(BASE_PROCESSES)
    invalid = [value for value in parsed if value < 1 or value > max_id]
    if invalid:
        raise ValueError(
            f"--process-ids contains invalid ids {invalid}; valid range is 1-{max_id}"
        )
    return parsed


def _remaining(started: float, limit: float | None) -> float | None:
    if limit is None:
        return None
    remaining = limit - (time.monotonic() - started)
    if remaining <= 0:
        raise TimeoutError(f"time limit of {limit:g} seconds exceeded")
    return remaining


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


def _compact_payload(payload: dict[str, Any]) -> dict[str, Any]:
    keep = {
        "available",
        "artifact_class",
        "generation_s",
        "key",
        "kind",
        "manifest",
        "process",
        "process_dir",
        "runtime_available",
        "runtime_backend",
        "runtime_unavailable_message",
        "schema_version",
        "target_runtime_s",
        "values",
    }
    compact = {key: value for key, value in payload.items() if key in keep}
    profile = payload.get("profile")
    if isinstance(profile, dict):
        compact["profile"] = {
            key: profile.get(key)
            for key in (
                "samples",
                "wall_us_per_point",
                "wall_us_per_point_error",
                "core_evaluator_us_per_point",
                "core_evaluator_us_per_point_error",
            )
            if key in profile
        }
    return compact


def _error_payload(exc: Exception) -> dict[str, Any]:
    error = str(exc)[-4000:]
    timeout_markers = (
        "timed out",
        "time limit",
        "Terminated: 15",
        "Interrupted system call",
    )
    status = (
        "timeout"
        if isinstance(exc, (TimeoutError, subprocess.TimeoutExpired))
        or any(marker in error for marker in timeout_markers)
        else "failed"
    )
    return {
        "status": status,
        "error": error,
        "finished_at": _now(),
    }


def _optional_float(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return number


def _first_payload_value(mode: dict[str, Any]) -> float | None:
    payload = mode.get("time_payload")
    if not isinstance(payload, dict):
        return None
    values = payload.get("values")
    if not isinstance(values, list) or not values:
        return None
    return _optional_float(values[0])


def _relative_difference(reference: float, value: float) -> float:
    return abs(value - reference) / max(abs(reference), abs(value), 1.0e-300)


def _format_seconds(value: float) -> str:
    return rf"\texttt{{{_format_sig(value, unit='s')}}}"


def _format_us(value: float) -> str:
    return rf"\texttt{{{_format_sig(value, unit='us')}}}"


def _format_sig(value: float, *, unit: str) -> str:
    if value == 0.0:
        return f"0 {unit}"
    abs_value = abs(value)
    if abs_value >= 1000:
        text = f"{value:.3g}"
    elif abs_value >= 100:
        text = f"{value:.0f}"
    elif abs_value >= 10:
        text = f"{value:.2f}".rstrip("0").rstrip(".")
    elif abs_value >= 1:
        text = f"{value:.3f}".rstrip("0").rstrip(".")
    else:
        text = f"{value:.3g}"
    return f"{text} {unit}"


def _format_ratio(value: float) -> str:
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    if value >= 1:
        return f"{value:.2f}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _latex_escape(value: str) -> str:
    return (
        value.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
        .replace("$", r"\$")
        .replace("{", r"\{")
        .replace("}", r"\}")
    )


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


if __name__ == "__main__":
    raise SystemExit(main())
