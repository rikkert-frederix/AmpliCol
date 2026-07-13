#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import signal
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from itertools import islice, product
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

DOCS_DIR = Path(__file__).resolve().parent
PYAMPLICOL_DIR = DOCS_DIR.parent
REPO_ROOT = PYAMPLICOL_DIR.parent
SRC_DIR = PYAMPLICOL_DIR / "src"
SCRIPTS_DIR = PYAMPLICOL_DIR / "scripts"
DEFAULT_DATA = DOCS_DIR / "result_matrix_data.json"
DEFAULT_TABLE = DOCS_DIR / "result_matrix_table.tex"
DEFAULT_SHARED_LC_TABLE = DOCS_DIR / "shared_lc_recycling_table.tex"
DEFAULT_NLC_DATA = DOCS_DIR / "result_matrix_nlc_data.json"
DEFAULT_NLC_TABLE = DOCS_DIR / "result_matrix_nlc_table.tex"
DEFAULT_FULL_DATA = DOCS_DIR / "result_matrix_full_data.json"
DEFAULT_FULL_TABLE = DOCS_DIR / "result_matrix_full_table.tex"
DEFAULT_OUTPUT_ROOT = DOCS_DIR / ".result_matrix_outputs"
DEFAULT_UFO_SM_SOURCE = PYAMPLICOL_DIR / "assets" / "models" / "json" / "sm" / "sm.json"
DEFAULT_UFO_SM_DATA = DOCS_DIR / "result_matrix_ufo_sm_data.json"
DEFAULT_UFO_SM_TABLE = DOCS_DIR / "result_matrix_ufo_sm_table.tex"
DEFAULT_UFO_SM_SHARED_LC_TABLE = DOCS_DIR / "shared_lc_recycling_ufo_sm_table.tex"
DEFAULT_UFO_SM_NLC_DATA = DOCS_DIR / "result_matrix_ufo_sm_nlc_data.json"
DEFAULT_UFO_SM_NLC_TABLE = DOCS_DIR / "result_matrix_ufo_sm_nlc_table.tex"
DEFAULT_UFO_SM_FULL_DATA = DOCS_DIR / "result_matrix_ufo_sm_full_data.json"
DEFAULT_UFO_SM_FULL_TABLE = DOCS_DIR / "result_matrix_ufo_sm_full_table.tex"
DEFAULT_UFO_SM_OUTPUT_ROOT = DOCS_DIR / ".result_matrix_ufo_sm_outputs"
DEFAULT_TIME_LIMIT_S = 900.0
DEFAULT_TARGET_RUNTIME_S = 10.0
DEFAULT_AMPLICOL_POINTS = 10000
DEFAULT_JOBS = 5
DEFAULT_N_CORES = 5
DEFAULT_PROCESS_WORKERS = 1
DEFAULT_COLUMNS_PER_TABLE = 3
DEFAULT_BATCH_SIZE = 128
DEFAULT_NON_LC_BATCH_SIZE = 64
DEFAULT_ALL_FLOW_BATCH_SIZE = 64
DEFAULT_SYMBOLICA_ITERATIONS = 10
DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE = 128
DEFAULT_ALL_FLOW_SYMBOLICA_OUTPUT_CHUNK_SIZE = 8192
DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT = True
DEFAULT_MEMORY_LIMIT_GB = 30.0
VALIDATION_REL_TOL = 1.0e-8
COMMAND_HEARTBEAT_S = 30.0
MEMORY_POLL_S = 1.0
MEMORY_NEAR_LIMIT_POLL_S = 0.25
MEMORY_HIGH_WATERMARK_FRACTION = 0.8
VALIDATION_ABS_TOL = 1.0e-16
_FIXED_HELICITY_CACHE: dict[tuple[str, str], dict[str, Any]] = {}
QCD_ONLY_EXTERNAL_SM_FAMILIES = frozenset(
    ("gg_tt_jets", "dd_tt_jets", "gg_gluons", "dd_3q_lines", "dd_4q_lines")
)

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from memory_watch_support import ProcessTreeMemoryMonitor  # noqa: E402
from pyamplicol.generic_validation import _library_runtime_per_point  # noqa: E402
from pyamplicol.core_types import ExternalMomentum  # noqa: E402
from pyamplicol.generic_artifact import (  # noqa: E402
    _evaluate_dag_raw_sums_by_warmup,
    build_generic_process_manifest,
    select_leading_color_sector_ids_from_plan,
)
from pyamplicol.generic_dag import (  # noqa: E402
    _root_source_helicity_mapping,
    _source_helicity_signature_by_bit,
)
from pyamplicol.phase_space import generic_validation_point  # noqa: E402
from pyamplicol.color_plan import (  # noqa: E402
    build_color_plan,
    lc_line_pairing_representative_ids,
    lc_topology_replay_partitions,
)
from pyamplicol.model import AmplicolSMLeadingColorModel  # noqa: E402
from pyamplicol.process_ir import build_process_ir  # noqa: E402
from pyamplicol.processes import ProcessOptions  # noqa: E402
from pyamplicol.reference import (  # noqa: E402
    AmplicolAdapter,
    amplicol_process_file_entry,
    amplicol_process_file_integrals,
)


class MemoryLimitExceeded(RuntimeError):
    def __init__(
        self,
        *,
        limit_gb: float,
        peak_memory_bytes: int,
        memory_metric: str,
        command: Sequence[str],
        peak_rss_bytes: int | None = None,
        peak_physical_footprint_bytes: int | None = None,
    ):
        self.limit_gb = float(limit_gb)
        self.peak_memory_bytes = int(peak_memory_bytes)
        self.peak_memory_gb = self.peak_memory_bytes / 1024**3
        self.memory_metric = str(memory_metric)
        self.peak_rss_bytes = (
            None if peak_rss_bytes is None else int(peak_rss_bytes)
        )
        self.peak_rss_gb = (
            None
            if self.peak_rss_bytes is None
            else self.peak_rss_bytes / 1024**3
        )
        self.peak_physical_footprint_bytes = (
            None
            if peak_physical_footprint_bytes is None
            else int(peak_physical_footprint_bytes)
        )
        self.peak_physical_footprint_gb = (
            None
            if self.peak_physical_footprint_bytes is None
            else self.peak_physical_footprint_bytes / 1024**3
        )
        self.command = tuple(str(item) for item in command)
        metric_label = self.memory_metric.replace("_", " ")
        super().__init__(
            f"memory limit {self.limit_gb:g} GB exceeded "
            f"({metric_label} {self.peak_memory_gb:.3f} GiB): "
            f"{' '.join(self.command)}"
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
            flags.append(f"--max-currents={self.max_currents}")
        if self.max_color_sectors is not None:
            flags.append(f"--max-color-sectors={self.max_color_sectors}")
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
            "pyAmpliCol.tex."
        )
    )
    parser.add_argument(
        "--color-accuracy",
        choices=("lc", "nlc", "full"),
        default="lc",
        help="Colour accuracy represented by the generated matrix table.",
    )
    parser.add_argument(
        "--model",
        default="BUILTIN_SM",
        metavar="SOURCE",
        help=(
            "SM model source used by pyAmpliCol: BUILTIN_SM, a UFO directory, "
            "loader JSON, or compiled pyAmpliCol model. AmpliCol remains the "
            "common external reference."
        ),
    )
    parser.add_argument(
        "--model-label",
        default=None,
        help="Concise source label used in the generated table heading.",
    )
    parser.add_argument(
        "--reference-data",
        type=Path,
        default=None,
        help=(
            "Seed only AmpliCol records and process metadata from another matrix "
            "cache. This permits matched external-model runs without relabeling "
            "built-in pyAmpliCol timings."
        ),
    )
    parser.add_argument(
        "--validation-report",
        type=Path,
        default=None,
        help=(
            "Optional validate_ufo_sm_fixture.py report whose matched numerical "
            "results are attached to seeded cells."
        ),
    )
    parser.add_argument("--data", type=Path, default=None)
    parser.add_argument("--table", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument(
        "--shared-lc-table",
        type=Path,
        default=None,
        help="Override the generated shared-LC summary table path.",
    )
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
        "--selected-jit-opt-level",
        type=int,
        choices=(1, 3),
        default=1,
        help=(
            "SymJIT optimization level for the selected-flow JIT artifact. "
            "The LC all-flow artifact remains at O1; O3 selected-flow entries "
            "are marked with a star in the rendered table."
        ),
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
        "--reuse-lc-all-flow",
        action="store_true",
        help=(
            "When refreshing LC JIT cells, regenerate and retime only the "
            "selected-flow artifact and retain the existing validated all-flow "
            "artifact and cache fields."
        ),
    )
    parser.add_argument(
        "--retime-existing",
        action="store_true",
        help=(
            "Rerun time-process for existing pyAmpliCol artifacts without "
            "regenerating them. Generation metadata and validation records are "
            "preserved; AmpliCol reference rows are not rerun."
        ),
    )
    parser.add_argument(
        "--only-missing",
        action="store_true",
        help="Only run requested modes whose current table entry is missing or failed.",
    )
    parser.add_argument(
        "--reset-cache",
        action="store_true",
        help=(
            "Start from an empty schema-valid cache before rendering or "
            "regenerating requested cells. The current files should be archived "
            "before using this option."
        ),
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
            "the script also recompiles pyAmpliCol.tex and refreshes pyAmpliCol.pdf."
        ),
    )
    parser.add_argument(
        "--clean-output-artifacts",
        action="store_true",
        help=(
            "Remove per-cell .result_matrix_outputs artifacts after their timings "
            "have been recorded in the matrix JSON cache."
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.retime_existing and args.clean_output_artifacts:
        parser.error("--retime-existing cannot be combined with --clean-output-artifacts")
    if args.retime_existing and args.only_missing:
        parser.error("--retime-existing cannot be combined with --only-missing")
    if args.retime_existing and not args.skip_amplicol:
        print(
            "[matrix] --retime-existing leaves AmpliCol reference rows unchanged",
            flush=True,
        )
        args.skip_amplicol = True
    color_accuracy = str(args.color_accuracy).lower()
    model_source = _normalized_model_source(args.model)
    model_profile = _matrix_model_profile(model_source)
    model_label = str(args.model_label or _matrix_model_label(model_profile))
    data_path = args.data or _default_data_path(color_accuracy, model_profile=model_profile)
    table_path = args.table or _default_table_path(
        color_accuracy,
        model_profile=model_profile,
    )
    output_root = (
        args.output_root
        or (
            DEFAULT_OUTPUT_ROOT
            if model_profile == "built-in-sm"
            else DEFAULT_UFO_SM_OUTPUT_ROOT
        )
    ).expanduser().resolve()
    shared_lc_table_path = args.shared_lc_table or (
        DEFAULT_SHARED_LC_TABLE
        if model_profile == "built-in-sm"
        else DEFAULT_UFO_SM_SHARED_LC_TABLE
    )

    data = _load_data(data_path)
    if args.reset_cache:
        if args.dry_run:
            print(f"[dry-run] would reset {data_path} and {table_path}")
        else:
            data = {}
    if args.reference_data is not None and not args.dry_run:
        _seed_reference_cache(
            data,
            _load_data(args.reference_data),
            color_accuracy=color_accuracy,
            validation_report=args.validation_report,
        )
    generate_ns = _parse_n_values(args.generate_data)
    show_ns = _parse_n_values(args.show_data)
    benchmark_model_source = model_source
    model_preparation: dict[str, Any] | None = None
    if generate_ns and not args.dry_run and model_profile != "built-in-sm":
        benchmark_model_source, model_preparation = _prepare_benchmark_model_source(
            model_source,
            output_root=output_root,
        )
    selected_base_keys = {str(key) for key in args.base_process}
    try:
        selected_process_ids = set(_parse_process_ids(args.process_ids))
    except ValueError as exc:
        parser.error(str(exc))
    if args.dry_run and not generate_ns:
        print("[dry-run] no --generate-data values supplied; no writes performed")
        return 0
    if not args.dry_run:
        _refresh_data_metadata(data, color_accuracy=color_accuracy)
        data["model_profile"] = model_profile
        data["model_label"] = model_label
        data["model_source"] = model_source
        if model_preparation is not None:
            data["compiled_model_source"] = benchmark_model_source
            data["model_preparation"] = model_preparation

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
                    if args.dry_run:
                        print(
                            f"[dry-run] id={process_id} n={n_final} "
                            f"{base.key}: not applicable"
                        )
                        continue
                    _record_not_applicable(data, base, n_final)
                    continue
                if args.dry_run:
                    print(
                        f"[dry-run] id={process_id} n={n_final} "
                        f"{base.key}: {process}"
                    )
                    continue
                work_items.append((base, n_final, process))
        if args.dry_run:
            return 0
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
        if process_workers > 1 and not args.no_recompile:
            print(
                "[matrix] --process-workers > 1 requested with PDF refresh "
                "enabled; forcing serial execution so the JSON, table, and "
                "pyAmpliCol.pdf are refreshed after each completed cell.",
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
                        output_root=output_root,
                        time_limit=None if args.time_limit <= 0 else args.time_limit,
                        target_runtime=args.target_runtime,
                        amplicol_points=args.amplicol_points,
                        jobs=args.jobs,
                        n_cores=args.n_cores,
                        selected_jit_opt_level=args.selected_jit_opt_level,
                        skip_amplicol=args.skip_amplicol,
                        skip_jit=args.skip_jit,
                        skip_cpp_o3=args.skip_cpp_o3,
                        reuse_lc_all_flow=args.reuse_lc_all_flow,
                        retime_existing=args.retime_existing,
                        only_missing=args.only_missing,
                        validate=args.validate,
                        clean_output_artifacts=args.clean_output_artifacts,
                        color_accuracy=color_accuracy,
                        model_source=benchmark_model_source,
                        model_profile=model_profile,
                    )
                    for base, n_final, process in work_items
                ]
                for future in as_completed(futures):
                    future.result()
                    _write_table_data_and_maybe_pdf(
                        data_path,
                        table_path,
                        data,
                        show_ns,
                        columns_per_table=args.columns_per_table,
                        color_accuracy=color_accuracy,
                        shared_lc_table_path=shared_lc_table_path,
                        refresh_pdf=not args.no_recompile,
                    )
        elif work_items:
            for base, n_final, process in work_items:
                _generate_case(
                    data,
                    base,
                    n_final,
                    process,
                    output_root=output_root,
                    time_limit=None if args.time_limit <= 0 else args.time_limit,
                    target_runtime=args.target_runtime,
                    amplicol_points=args.amplicol_points,
                    jobs=args.jobs,
                    n_cores=args.n_cores,
                    selected_jit_opt_level=args.selected_jit_opt_level,
                    skip_amplicol=args.skip_amplicol,
                    skip_jit=args.skip_jit,
                    skip_cpp_o3=args.skip_cpp_o3,
                    reuse_lc_all_flow=args.reuse_lc_all_flow,
                    retime_existing=args.retime_existing,
                    only_missing=args.only_missing,
                    validate=args.validate,
                    clean_output_artifacts=args.clean_output_artifacts,
                    color_accuracy=color_accuracy,
                    model_source=benchmark_model_source,
                    model_profile=model_profile,
                )
                _write_table_data_and_maybe_pdf(
                    data_path,
                    table_path,
                    data,
                    show_ns,
                    columns_per_table=args.columns_per_table,
                    color_accuracy=color_accuracy,
                    shared_lc_table_path=shared_lc_table_path,
                    refresh_pdf=not args.no_recompile,
                )

    _write_table_data_and_maybe_pdf(
        data_path,
        table_path,
        data,
        show_ns,
        columns_per_table=args.columns_per_table,
        color_accuracy=color_accuracy,
        shared_lc_table_path=shared_lc_table_path,
        refresh_pdf=not args.no_recompile,
    )
    return 0


def _write_table_data_and_maybe_pdf(
    data_path: Path,
    table_path: Path,
    data: dict[str, Any],
    show_ns: Sequence[int],
    *,
    columns_per_table: int,
    color_accuracy: str,
    shared_lc_table_path: Path,
    refresh_pdf: bool,
) -> None:
    table = render_latex_table(
        data,
        show_ns,
        columns_per_table=columns_per_table,
        color_accuracy=color_accuracy,
    )
    table_path.parent.mkdir(parents=True, exist_ok=True)
    table_path.write_text(table, encoding="utf-8")
    if color_accuracy == "lc":
        shared_lc_table = render_shared_lc_recycling_table(data)
        shared_lc_table_path.write_text(shared_lc_table, encoding="utf-8")
        print(f"wrote {shared_lc_table_path}")
    _write_data(data_path, data)
    print(f"wrote {table_path}")
    print(f"wrote {data_path}")
    if refresh_pdf:
        _refresh_pyamplicol_pdf(table_path)


def _refresh_pyamplicol_pdf(table_path: Path) -> None:
    known_tables = {
        DEFAULT_TABLE.resolve(),
        DEFAULT_NLC_TABLE.resolve(),
        DEFAULT_FULL_TABLE.resolve(),
        DEFAULT_UFO_SM_TABLE.resolve(),
        DEFAULT_UFO_SM_NLC_TABLE.resolve(),
        DEFAULT_UFO_SM_FULL_TABLE.resolve(),
    }
    if table_path.resolve() not in known_tables:
        print(
            "skipped PDF recompilation: --table does not point to one of "
            f"{sorted(str(path) for path in known_tables)}",
            flush=True,
        )
        return
    document = DOCS_DIR / "pyAmpliCol.tex"
    if not document.exists():
        print(f"skipped PDF refresh: {document} does not exist", flush=True)
        return
    env = dict(os.environ)
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    print("refreshing pyAmpliCol/docs/pyAmpliCol.pdf", flush=True)
    completed = subprocess.run(
        [
            "latexmk",
            "-pdf",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "pyAmpliCol.tex",
        ],
        cwd=DOCS_DIR,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "pyAmpliCol.tex recompilation failed with exit code "
            f"{completed.returncode}\nstdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    target = DOCS_DIR / "pyAmpliCol.pdf"
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
    selected_jit_opt_level: int,
    skip_amplicol: bool,
    skip_jit: bool,
    skip_cpp_o3: bool,
    reuse_lc_all_flow: bool,
    retime_existing: bool,
    only_missing: bool,
    validate: bool,
    clean_output_artifacts: bool,
    color_accuracy: str,
    model_source: str,
    model_profile: str,
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
    structural_reason = _structural_color_accuracy_unsupported_reason(
        process,
        base=base,
        color_accuracy=color_accuracy,
    )
    if structural_reason is not None:
        unsupported = _unsupported_mode_payload(
            structural_reason,
            color_accuracy=color_accuracy,
        )
        case["amplicol"] = dict(
            unsupported,
            mode=(
                "AmpliCol - using-library"
                if color_accuracy == "lc"
                else f"AmpliCol - raw color library | {color_accuracy}"
            ),
        )
        case["pyamplicol_jit"] = dict(
            unsupported,
            mode="pyAmpliCol - staged-DAG | JIT",
        )
        case["pyamplicol_cpp_o3"] = dict(
            unsupported,
            mode="pyAmpliCol - staged-DAG | C++ O3",
        )
        case["validation"] = {
            "status": "unsupported",
            "kind": "same-deterministic-point",
            "reason": structural_reason,
            "tolerance": VALIDATION_REL_TOL,
            "relative_tolerance": VALIDATION_REL_TOL,
            "absolute_tolerance": VALIDATION_ABS_TOL,
            "finished_at": _now(),
        }
        case["status"] = "unsupported"
        case["updated_at"] = _now()
        return

    pyamplicol_base = _pyamplicol_base_for_case(
        base,
        n_final=n_final,
        color_accuracy=color_accuracy,
    )
    amplicol_refreshed = False
    amplicol_settings = _amplicol_matrix_settings(
        amplicol_points=amplicol_points,
        jobs=jobs,
        color_accuracy=color_accuracy,
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
            target_runtime=target_runtime,
            amplicol_points=amplicol_points,
            jobs=jobs,
            matrix_settings=amplicol_settings,
            color_accuracy=color_accuracy,
        )
        amplicol_refreshed = _ok(_mode(case, "amplicol"))
    try:
        reference_color_order = (
            _reference_color_order_for_case(case)
            if color_accuracy == "lc"
            else None
        )
        selected_lc_sector_ids = (
            _selected_lc_sector_ids_for_case(
                process,
                base,
                case,
            )
            if color_accuracy == "lc"
            else None
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
        selected_jit_opt_level=selected_jit_opt_level,
        selected_lc_sector_ids=selected_lc_sector_ids,
        reference_color_order=reference_color_order,
        base=pyamplicol_base,
        color_accuracy=color_accuracy,
        model_source=model_source,
        model_profile=model_profile,
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
                base=pyamplicol_base,
                color_accuracy=color_accuracy,
                model_source=model_source,
                model_profile=model_profile,
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

    existing_jit_mode = dict(_mode(case, "pyamplicol_jit"))
    if not skip_jit and _should_run_mode(
        case,
        "pyamplicol_jit",
        only_missing=only_missing,
        expected_settings=jit_settings,
    ):
        jit_output_dir = _pyamplicol_output_dir(
            output_root,
            color_accuracy=color_accuracy,
            base=base,
            n_final=n_final,
            backend_key="jit",
        )
        if retime_existing:
            refreshed_jit_mode = _retime_pyamplicol_case(
                existing_jit_mode,
                backend_key="jit",
                output_dir=jit_output_dir,
                target_runtime=target_runtime,
                color_accuracy=color_accuracy,
            )
        else:
            refreshed_jit_mode = _run_pyamplicol_case(
                process,
                base=pyamplicol_base,
                n_final=n_final,
                backend_key="jit",
                output_dir=jit_output_dir,
                deadline_started=time.monotonic(),
                time_limit=None,
                target_runtime=target_runtime,
                n_cores=n_cores,
                selected_jit_opt_level=selected_jit_opt_level,
                selected_lc_sector_ids=selected_lc_sector_ids,
                reference_color_order=reference_color_order,
                matrix_settings=jit_settings,
                color_accuracy=color_accuracy,
                model_source=model_source,
                model_profile=model_profile,
                reuse_lc_all_flow=reuse_lc_all_flow,
                existing_mode=existing_jit_mode,
            )
        if reuse_lc_all_flow and not retime_existing:
            refresh_validation = _selected_refresh_validation(
                existing_jit_mode,
                refreshed_jit_mode,
            )
            refreshed_jit_mode["selected_refresh_validation"] = refresh_validation
            if refresh_validation.get("status") == "failed":
                refreshed_jit_mode["status"] = "validation_failed"
                refreshed_jit_mode["error"] = (
                    "selected-flow refresh changed the validated value: "
                    f"relative difference "
                    f"{refresh_validation.get('relative_difference')}"
                )
        case["pyamplicol_jit"] = refreshed_jit_mode
        if clean_output_artifacts:
            shutil.rmtree(jit_output_dir, ignore_errors=True)
    cpp_o3_settings = _pyamplicol_matrix_settings(
        backend_key="cpp_o3",
        target_runtime=target_runtime,
        n_cores=n_cores,
        selected_lc_sector_ids=selected_lc_sector_ids,
        reference_color_order=reference_color_order,
        base=pyamplicol_base,
        color_accuracy=color_accuracy,
        model_source=model_source,
        model_profile=model_profile,
    )
    if not skip_cpp_o3 and _should_run_mode(
        case,
        "pyamplicol_cpp_o3",
        only_missing=only_missing,
        expected_settings=cpp_o3_settings,
    ):
        cpp_o3_output_dir = _pyamplicol_output_dir(
            output_root,
            color_accuracy=color_accuracy,
            base=base,
            n_final=n_final,
            backend_key="cpp_o3",
        )
        if retime_existing:
            case["pyamplicol_cpp_o3"] = _retime_pyamplicol_case(
                dict(_mode(case, "pyamplicol_cpp_o3")),
                backend_key="cpp_o3",
                output_dir=cpp_o3_output_dir,
                target_runtime=target_runtime,
                color_accuracy=color_accuracy,
            )
        else:
            case["pyamplicol_cpp_o3"] = _run_pyamplicol_case(
                process,
                base=pyamplicol_base,
                n_final=n_final,
                backend_key="cpp_o3",
                output_dir=cpp_o3_output_dir,
                deadline_started=time.monotonic(),
                time_limit=time_limit,
                target_runtime=target_runtime,
                n_cores=n_cores,
                selected_jit_opt_level=1,
                selected_lc_sector_ids=selected_lc_sector_ids,
                reference_color_order=reference_color_order,
                matrix_settings=cpp_o3_settings,
                color_accuracy=color_accuracy,
                model_source=model_source,
                model_profile=model_profile,
            )
        if clean_output_artifacts:
            shutil.rmtree(cpp_o3_output_dir, ignore_errors=True)
    should_validate = validate and not skip_amplicol and _ok(_mode(case, "amplicol"))
    if should_validate and only_missing and _validation_record_ok(case.get("validation")):
        should_validate = False
    if should_validate:
        case["validation"] = _run_validation_case(
            process,
            base=base,
            case=case,
            deadline_started=started,
            time_limit=None,
            jobs=jobs,
            library_current=amplicol_refreshed,
            color_accuracy=color_accuracy,
        )
    elif validate and color_accuracy != "lc":
        unavailable_reason = _fortran_color_reference_unavailable_reason(
            _mode(case, "amplicol"),
            _mode(case, "pyamplicol_jit"),
        )
        if unavailable_reason is not None:
            case["validation"] = {
                "status": "unsupported",
                "kind": "same-deterministic-point",
                "reason": unavailable_reason,
                "tolerance": VALIDATION_REL_TOL,
                "relative_tolerance": VALIDATION_REL_TOL,
                "absolute_tolerance": VALIDATION_ABS_TOL,
                "finished_at": _now(),
            }
    case["status"] = "done"
    case["updated_at"] = _now()


def _pyamplicol_base_for_case(
    base: BaseProcess,
    *,
    n_final: int,
    color_accuracy: str,
) -> BaseProcess:
    max_currents = -1
    max_color_sectors = -1
    lc_sector_strategy = "all" if color_accuracy == "lc" else base.lc_sector_strategy
    if color_accuracy != "lc" and n_final < 5:
        max_currents = base.max_currents
        max_color_sectors = base.max_color_sectors
    if (
        max_currents == base.max_currents
        and max_color_sectors == base.max_color_sectors
        and lc_sector_strategy == base.lc_sector_strategy
    ):
        return base
    return replace(
        base,
        max_currents=max_currents,
        max_color_sectors=max_color_sectors,
        lc_sector_strategy=lc_sector_strategy,
    )


def _run_amplicol_case(
    process: str,
    *,
    base: BaseProcess,
    deadline_started: float,
    time_limit: float | None,
    target_runtime: float,
    amplicol_points: int,
    jobs: int,
    matrix_settings: dict[str, Any],
    color_accuracy: str,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "mode": (
            "AmpliCol - using-library"
            if color_accuracy == "lc"
            else f"AmpliCol - raw color library | {color_accuracy}"
        ),
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
        if color_accuracy != "lc":
            if point is None:
                raise RuntimeError(
                    "NLC/full AmpliCol colour-library probe requires a deterministic "
                    "validation point"
                )
            build = adapter.prepare_library(
                process,
                options=base.process_options(),
                process_list_backend="python",
                warmup_particles=point,
                warmup_points=1,
                raw=True,
                color_complete=True,
            )
            remaining = _remaining(deadline_started, time_limit)
            adapter = AmplicolAdapter(REPO_ROOT, jobs=jobs, timeout=remaining)
            (
                run,
                probe_points,
                calibration,
                calibration_commands,
            ) = _run_color_library_probe_adaptive(
                adapter,
                process,
                color_accuracy=color_accuracy,
                particles=point,
                max_points=max(1, amplicol_points),
                target_runtime=target_runtime,
                process_file=build.process_file,
                options=base.process_options(),
                process_list_backend="python",
                color_complete=True,
            )
            runtime_s = _color_probe_runtime_per_point(
                run,
                probe_points,
            )
            payload.update(
                {
                    "status": "ok",
                    "generation_s": build.total_command_time_s,
                    "reference_probe": "generated_library_color_probe_raw_adaptive",
                    "process_file": str(run.process_file),
                    "process_list_backend": "python",
                    "reference_color_order": None,
                    "reference_color_order_process_file": None,
                    "runtime_probe_points": probe_points,
                    "runtime_probe_max_points": max(1, int(amplicol_points)),
                    "runtime_probe_target_runtime_s": float(target_runtime),
                    "runtime_probe_calibration": calibration,
                    "runtime_us_per_point": (
                        None if runtime_s is None else 1.0e6 * runtime_s
                    ),
                    "reference_value": run.first_point_matrix_element,
                    "color_probe_components": (
                        None
                        if run.color_probe_components is None
                        else list(run.color_probe_components)
                    ),
                    "color_probe_raw_components": (
                        None
                        if run.color_probe_raw_components is None
                        else list(run.color_probe_raw_components)
                    ),
                    "commands": [
                        {
                            "args": command.args,
                            "elapsed_s": command.elapsed_s,
                            "returncode": command.returncode,
                        }
                        for command in (
                            *build.commands,
                            *calibration_commands,
                            *run.commands,
                        )
                    ],
                    "commands_s": [
                        command.elapsed_s
                        for command in (
                            *build.commands,
                            *calibration_commands,
                            *run.commands,
                        )
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
            return payload
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
            reuse_process_file=True,
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
        helicity_started = time.perf_counter()
        print(
            f"[matrix] selecting replay-invariant fixed helicity: {process}",
            flush=True,
        )
        fixed_helicity = _fixed_helicity_choice(process, base)
        print(
            "[matrix] fixed helicity ready after "
            f"{time.perf_counter() - helicity_started:.3f}s: "
            f"{fixed_helicity['source_helicities_cli']}",
            flush=True,
        )
        all_flow_supported = _lc_fixed_helicity_all_flow_supported(process, base)
        all_flow_runtime = None
        all_flow_probe_points = None
        all_flow_runtime_s = None
        all_flow_generation_s = None
        if all_flow_supported:
            all_flow_runtime, all_flow_probe_points = _run_amplicol_all_flow_probe(
                adapter,
                process,
                base=base,
                process_file=build.process_file,
                point=point,
                helicities=fixed_helicity["amplicol_helicities"],
                target_runtime=target_runtime,
                amplicol_points=amplicol_points,
                process_list_backend=process_list_backend,
            )
            all_flow_runtime_s = _color_probe_runtime_per_point(
                all_flow_runtime,
                all_flow_probe_points,
            )
            all_flow_generation_s = _amplicol_color_probe_setup_s_from_run(
                all_flow_runtime
            )
        payload.update(
            {
                "status": "ok",
                "generation_s": build.total_command_time_s,
                "all_flow_generation_s": all_flow_generation_s,
                "all_flow_generation_source": (
                    "amplicol_color_probe_imode2_setup_build"
                    if all_flow_runtime is not None
                    else None
                ),
                "reference_probe": reference_probe,
                "process_file": str(build.process_file),
                "process_list_backend": (
                    process_list_backend
                ),
                "process_file_scope": "all_entries_generated_first_group_first_integral_timed",
                "reference_color_order": mapped_color_order,
                "reference_color_order_process_file": raw_color_order,
                "runtime_us_per_point": (
                    None if runtime_s is None else 1.0e6 * runtime_s
                ),
                "all_flow_reference_probe": "amplicol_color_probe_fixed_helicity_all_flows",
                "all_flow_helicity_mode": fixed_helicity["mode"],
                "all_flow_source_helicities": {
                    str(label): int(helicity)
                    for label, helicity in fixed_helicity[
                        "source_helicities"
                    ].items()
                },
                "all_flow_amplicol_helicities": list(
                    fixed_helicity["amplicol_helicities"]
                ),
                "all_flow_value_validation_enabled": bool(
                    fixed_helicity["value_validation_enabled"]
                    and all_flow_supported
                ),
                "all_flow_validation_note": (
                    fixed_helicity["validation_note"]
                    if all_flow_supported
                    else "AmpliCol fixed-helicity all-flow probe is unsupported for more than two quark lines"
                ),
                "all_flow_runtime_probe_points": all_flow_probe_points,
                "all_flow_runtime_probe_target_runtime_s": float(target_runtime),
                "all_flow_runtime_us_per_point": (
                    None if all_flow_runtime_s is None else 1.0e6 * all_flow_runtime_s
                ),
                "all_flow_reference_value": (
                    None
                    if all_flow_runtime is None
                    else all_flow_runtime.first_point_matrix_element
                ),
                "all_flow_current_count": (
                    None
                    if all_flow_runtime is None
                    else all_flow_runtime.color_probe_current_count
                ),
                "all_flow_vertex_count": (
                    None
                    if all_flow_runtime is None
                    else all_flow_runtime.color_probe_vertex_count
                ),
                "all_flow_amplitude_count": (
                    None
                    if all_flow_runtime is None
                    else all_flow_runtime.color_probe_amplitude_count
                ),
                "all_flow_color_order_count": (
                    None
                    if all_flow_runtime is None
                    else all_flow_runtime.color_probe_color_order_count
                ),
                "all_flow_color_probe_components": (
                    None
                    if all_flow_runtime is None
                    or all_flow_runtime.color_probe_components is None
                    else list(all_flow_runtime.color_probe_components)
                ),
                "all_flow_color_probe_raw_components": (
                    None
                    if all_flow_runtime is None
                    or all_flow_runtime.color_probe_raw_components is None
                    else list(all_flow_runtime.color_probe_raw_components)
                ),
                "commands": [
                    {
                        "args": command.args,
                        "elapsed_s": command.elapsed_s,
                        "returncode": command.returncode,
                    }
                        for command in (
                            *build.commands,
                            *run.commands,
                            *((all_flow_runtime.commands) if all_flow_runtime is not None else ()),
                        )
                    ],
                    "commands_s": [
                        command.elapsed_s
                        for command in (
                            *build.commands,
                            *run.commands,
                            *((all_flow_runtime.commands) if all_flow_runtime is not None else ()),
                        )
                    ],
                "timing_rows": [
                    {
                        "label": row.label,
                        "seconds": row.seconds,
                        "note": row.note,
                    }
                    for row in run.timing_rows
                ],
                "all_flow_timing_rows": [
                    {
                        "label": row.label,
                        "seconds": row.seconds,
                        "note": row.note,
                    }
                    for row in (
                        all_flow_runtime.timing_rows
                        if all_flow_runtime is not None
                        else ()
                    )
                ],
                "matrix_settings": matrix_settings,
                "finished_at": _now(),
            }
        )
    except Exception as exc:  # noqa: BLE001 - benchmark harness records failures.
        payload.update(_error_payload(exc))
        if color_accuracy != "lc":
            reference_gap = _fortran_color_reference_error_reason(str(exc))
            if reference_gap is not None:
                payload["status"] = "fortran_reference_unavailable"
                payload["reason"] = reference_gap
        if _amplicol_error_is_unsupported(str(exc)):
            payload["status"] = "unsupported"
    return payload


def _run_amplicol_all_flow_probe(
    adapter: AmplicolAdapter,
    process: str,
    *,
    base: BaseProcess,
    process_file: str | Path,
    point: tuple[ExternalMomentum, ...] | None,
    helicities: Sequence[int] | None,
    target_runtime: float,
    amplicol_points: int,
    process_list_backend: str,
) -> tuple[Any, int]:
    particles = point if point is not None else generic_validation_point(process)
    adaptive = target_runtime > 0.0
    max_points = max(1, int(amplicol_points)) if adaptive else 1
    run = adapter.run_color_probe(
        process,
        color_accuracy="lc",
        particles=particles,
        helicities=helicities,
        points=max_points,
        process_file=process_file,
        options=base.process_options(),
        process_list_backend=process_list_backend,  # type: ignore[arg-type]
        target_runtime=target_runtime if adaptive else None,
    )
    probe_points = getattr(run, "color_probe_points", None)
    if probe_points is None:
        probe_points = max_points
    return run, probe_points


def _fixed_helicity_choice(
    process: str,
    base: BaseProcess,
) -> dict[str, Any]:
    cache_key = (str(process), base.key)
    cached = _FIXED_HELICITY_CACHE.get(cache_key)
    if cached is not None:
        return cached

    model = AmplicolSMLeadingColorModel()
    process_ir = build_process_ir(
        process,
        color_accuracy="lc",
        options=base.process_options(),
    )
    fixed_root_choice = _fixed_helicity_from_available_lc_root(
        process,
        base,
        process_ir=process_ir,
        model=model,
    )
    fallback_used = False
    replay_helicity_invariant = False
    if fixed_root_choice is None:
        fallback_used = True
        source_helicities = {}
        for leg in sorted(process_ir.legs, key=lambda item: int(item.label)):
            if leg.outgoing_pdg is None:
                continue
            pdg = int(leg.outgoing_pdg)
            source_helicities[int(leg.label)] = _preferred_source_helicity_for_label(
                model,
                pdg,
                int(leg.label),
            )
    else:
        source_helicities = dict(fixed_root_choice[0])
        replay_helicity_invariant = bool(fixed_root_choice[1])

    amplicol_helicities: list[int] = []
    basis_comparable = True
    for leg in sorted(process_ir.legs, key=lambda item: int(item.label)):
        if leg.outgoing_pdg is None:
            continue
        pdg = int(leg.outgoing_pdg)
        helicity = source_helicities.get(int(leg.label))
        if helicity is None:
            helicity = _preferred_source_helicity_for_label(
                model,
                pdg,
                int(leg.label),
            )
            source_helicities[int(leg.label)] = helicity
        # amplitude_QCD crosses incoming momenta internally and expects the
        # same source-helicity convention that pyAmpliCol stores.
        amplicol_helicities.append(helicity)
        if model.mass(pdg) != 0.0:
            basis_comparable = False
    result = {
        "mode": "fixed-source-helicity",
        "source_helicities": source_helicities,
        "source_helicities_cli": ",".join(
            f"{label}={helicity}"
            for label, helicity in sorted(source_helicities.items())
        ),
        "amplicol_helicities": tuple(amplicol_helicities),
        "value_validation_enabled": basis_comparable,
        "validation_note": (
            (
                "fixed source-helicity basis is massless and value-comparable "
                "for the materialized all-sector LC artifact"
            )
            if basis_comparable
            else "fixed massive-spin basis is timing-only; selected-flow spin-summed validation remains authoritative"
        ),
        "replay_helicity_invariant": replay_helicity_invariant,
        "selection_source": (
            "fallback-preferred-source-helicity"
            if fallback_used
            else "selected-lc-root-signature"
        ),
    }
    _FIXED_HELICITY_CACHE[cache_key] = result
    return result


def _lc_fixed_helicity_all_flow_supported(
    process: str,
    base: BaseProcess,
) -> bool:
    process_ir = build_process_ir(
        process,
        color_accuracy="lc",
        options=base.process_options(),
    )
    return len(process_ir.quark_labels) <= 2


def _fixed_helicity_from_available_lc_root(
    process: str,
    base: BaseProcess,
    *,
    process_ir: Any | None = None,
    model: AmplicolSMLeadingColorModel | None = None,
) -> tuple[dict[int, int], bool] | None:
    model = model or AmplicolSMLeadingColorModel()
    process_ir = process_ir or build_process_ir(
        process,
        color_accuracy="lc",
        options=base.process_options(),
    )
    try:
        replay_equalities = _lc_replay_helicity_equalities(process, base)
    except Exception:
        replay_equalities = ()
    constrained_candidates = _replay_invariant_helicity_candidates(
        process_ir,
        model,
        replay_equalities,
    )
    if constrained_candidates:
        candidate_started = time.perf_counter()
        print(
            "[matrix] checking "
            f"{len(constrained_candidates)} constrained fixed-helicity candidates",
            flush=True,
        )
        point = generic_validation_point(process)
        best_candidate: dict[int, int] | None = None
        best_raw_sum = 0.0
        for candidate in constrained_candidates:
            try:
                candidate_manifest = build_generic_process_manifest(
                    process_ir,
                    model=model,
                    color_accuracy="lc",
                    options=base.process_options(),
                    max_color_sectors=-1,
                    selected_color_sector_ids={0},
                    selected_source_helicities=candidate,
                    numerical_filter_current=False,
                    numerical_current_merging=False,
                )
                if not candidate_manifest.dag.amplitude_roots:
                    continue
                raw_sums = _evaluate_dag_raw_sums_by_warmup(
                    candidate_manifest.dag,
                    model,
                    (point,),
                )
            except Exception:  # noqa: BLE001 - retain exact fallback below.
                continue
            raw_sum = abs(float(raw_sums[0])) if raw_sums else 0.0
            if raw_sum > best_raw_sum:
                best_candidate = candidate
                best_raw_sum = raw_sum
        if best_candidate is not None and best_raw_sum > 1.0e-280:
            invariant = _fixed_helicity_signature_is_replay_invariant(
                tuple(sorted(best_candidate.items())),
                replay_equalities,
            )
            print(
                "[matrix] constrained fixed-helicity search accepted "
                + ",".join(
                    f"{label}={helicity}"
                    for label, helicity in sorted(best_candidate.items())
                )
                + f" after {time.perf_counter() - candidate_started:.3f}s "
                + f"(raw sum {best_raw_sum:.6g})",
                flush=True,
            )
            return best_candidate, invariant
        print(
            "[matrix] constrained fixed-helicity search found no numerically "
            "nonzero root after "
            f"{time.perf_counter() - candidate_started:.3f}s; "
            "falling back to all-helicity root discovery",
            flush=True,
        )

    try:
        manifest = build_generic_process_manifest(
            process,
            color_accuracy="lc",
            options=base.process_options(),
            selected_color_sector_ids={0},
            max_color_sectors=-1,
            numerical_filter_current=False,
            numerical_current_merging=False,
        )
    except Exception:
        return None
    source_by_bit = _source_helicity_signature_by_bit(manifest.dag)
    candidates: set[tuple[tuple[int, int], ...]] = set()
    for root in manifest.dag.amplitude_roots:
        source_map = _root_source_helicity_mapping(
            manifest.dag,
            root,
            source_by_bit,
        )
        if source_map:
            candidates.add(tuple(sorted(source_map.items())))
    if not candidates:
        return None
    selected = min(
        candidates,
        key=lambda signature: (
            0
            if _fixed_helicity_signature_is_replay_invariant(
                signature,
                replay_equalities,
            )
            else 1,
            _fixed_helicity_signature_sort_key(signature),
        ),
    )
    invariant = _fixed_helicity_signature_is_replay_invariant(
        selected,
        replay_equalities,
    )
    return {int(label): int(helicity) for label, helicity in selected}, invariant


def _replay_invariant_helicity_candidates(
    process_ir: Any,
    model: AmplicolSMLeadingColorModel,
    replay_equalities: Sequence[tuple[int, int]],
) -> tuple[dict[int, int], ...]:
    allowed_by_label: dict[int, set[int]] = {}
    pdg_by_label: dict[int, int] = {}
    for leg in sorted(process_ir.legs, key=lambda item: int(item.label)):
        if leg.outgoing_pdg is None:
            continue
        label = int(leg.label)
        pdg = int(leg.outgoing_pdg)
        allowed = {
            int(state.helicity)
            for state in model.source_spin_states(pdg)
        }
        if not allowed:
            return ()
        allowed_by_label[label] = allowed
        pdg_by_label[label] = pdg
    if not allowed_by_label:
        return ()

    parents = {label: label for label in allowed_by_label}

    def find(label: int) -> int:
        root = label
        while parents[root] != root:
            root = parents[root]
        while parents[label] != label:
            parent = parents[label]
            parents[label] = root
            label = parent
        return root

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        if left_root < right_root:
            parents[right_root] = left_root
        else:
            parents[left_root] = right_root

    for left, right in replay_equalities:
        left = int(left)
        right = int(right)
        if left in parents and right in parents:
            union(left, right)

    components: dict[int, list[int]] = {}
    for label in sorted(allowed_by_label):
        components.setdefault(find(label), []).append(label)

    priority = {-1: 0, 1: 1, 0: 2}
    component_options: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for labels in components.values():
        common = set.intersection(*(allowed_by_label[label] for label in labels))
        if not common:
            return ()
        preferred = {
            label: _preferred_source_helicity_for_label(
                model,
                pdg_by_label[label],
                label,
            )
            for label in labels
        }
        ordered = tuple(
            sorted(
                common,
                key=lambda helicity: (
                    -sum(value == helicity for value in preferred.values()),
                    priority.get(int(helicity), 3),
                    abs(int(helicity)),
                    int(helicity),
                ),
            )
        )
        component_options.append((tuple(labels), ordered))

    combinations = product(
        *(helicities for _, helicities in component_options)
    )
    candidates: list[dict[int, int]] = []
    for selected in islice(combinations, 256):
        candidate: dict[int, int] = {}
        for (labels, _), helicity in zip(component_options, selected, strict=True):
            for label in labels:
                candidate[label] = int(helicity)
        candidates.append(candidate)
    candidates.sort(
        key=lambda candidate: _fixed_helicity_signature_sort_key(
            tuple(sorted(candidate.items()))
        )
    )
    return tuple(candidates)


def _lc_replay_helicity_equalities(
    process: str,
    base: BaseProcess,
) -> tuple[tuple[int, int], ...]:
    color_plan = build_color_plan(
        process,
        color_accuracy="lc",
        options=base.process_options(),
        max_sectors=-1,
    )
    equalities: set[tuple[int, int]] = set()
    for partition in lc_topology_replay_partitions(color_plan):
        for permutation in partition.label_permutations:
            for representative_label, sector_label in permutation:
                left = int(representative_label)
                right = int(sector_label)
                if left == right:
                    continue
                equalities.add(tuple(sorted((left, right))))
    return tuple(sorted(equalities))


def _fixed_helicity_signature_is_replay_invariant(
    signature: tuple[tuple[int, int], ...],
    equalities: Sequence[tuple[int, int]],
) -> bool:
    helicities = {int(label): int(helicity) for label, helicity in signature}
    return all(
        helicities.get(int(left)) == helicities.get(int(right))
        for left, right in equalities
        if int(left) in helicities and int(right) in helicities
    )


def _fixed_helicity_signature_sort_key(
    signature: tuple[tuple[int, int], ...],
) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    helicities = [int(helicity) for _, helicity in signature]
    return (
        abs(sum(helicities)),
        1 if len(set(helicities)) == 1 else 0,
        signature,
    )


def _preferred_source_helicity_for_label(
    model: AmplicolSMLeadingColorModel,
    pdg: int,
    label: int,
) -> int:
    helicities = {int(state.helicity) for state in model.source_spin_states(pdg)}
    if {-1, 1}.issubset(helicities):
        return -1 if int(label) % 2 else 1
    for preferred in (-1, 1, 0):
        if preferred in helicities:
            return preferred
    return next(iter(sorted(helicities)), 0)


def _first_process_file_color_order(process_file: str | Path) -> list[int] | None:
    first_entry = next(iter(amplicol_process_file_integrals(process_file)), None)
    if first_entry is None:
        return None
    group, integral = first_entry
    entry = amplicol_process_file_entry(process_file, group=group, integral=integral)
    if entry is None:
        return None
    return [int(label) for label in entry["color_order"]]


def _run_color_library_probe_adaptive(
    adapter: AmplicolAdapter,
    process: str,
    *,
    color_accuracy: str,
    particles: Sequence[ExternalMomentum],
    max_points: int,
    target_runtime: float,
    process_file: str | Path,
    options: ProcessOptions,
    process_list_backend: str,
    color_complete: bool,
) -> tuple[Any, int, dict[str, Any], tuple[Any, ...]]:
    calibration_points = 1
    calibration_run = adapter.run_color_library_probe(
        process,
        color_accuracy=color_accuracy,
        particles=particles,
        points=calibration_points,
        process_file=process_file,
        options=options,
        process_list_backend=process_list_backend,  # type: ignore[arg-type]
        color_complete=color_complete,
    )
    calibration_runtime_s = _color_probe_runtime_per_point(
        calibration_run,
        calibration_points,
    )
    if calibration_runtime_s is None or calibration_runtime_s <= 0.0:
        probe_points = max(1, int(max_points))
    else:
        requested_points = int(math.ceil(max(0.0, target_runtime) / calibration_runtime_s))
        probe_points = max(1, min(int(max_points), requested_points))
    if probe_points == calibration_points:
        run = calibration_run
        final_commands: tuple[Any, ...] = ()
    else:
        run = adapter.run_color_library_probe(
            process,
            color_accuracy=color_accuracy,
            particles=particles,
            points=probe_points,
            process_file=process_file,
            options=options,
            process_list_backend=process_list_backend,  # type: ignore[arg-type]
            color_complete=color_complete,
        )
        final_commands = run.commands
    calibration_payload = {
        "points": calibration_points,
        "runtime_us_per_point": (
            None if calibration_runtime_s is None else 1.0e6 * calibration_runtime_s
        ),
        "target_runtime_s": float(target_runtime),
        "max_points": int(max_points),
        "selected_points": int(probe_points),
    }
    return run, probe_points, calibration_payload, (
        calibration_run.commands if final_commands else ()
    )


def _color_probe_runtime_per_point(run: Any, points: int) -> float | None:
    for row in getattr(run, "timing_rows", ()):
        if str(getattr(row, "label", "")).strip().lower() == "total":
            return float(getattr(row, "seconds")) / max(1, int(points))
    commands = getattr(run, "commands", ())
    if commands:
        return float(commands[-1].elapsed_s) / max(1, int(points))
    return None


def _amplicol_color_probe_setup_s_from_run(run: Any) -> float | None:
    return _amplicol_color_probe_setup_s_from_records(
        getattr(run, "commands", ()),
        getattr(run, "timing_rows", ()),
    )


def _amplicol_color_probe_setup_s_from_records(
    commands: Sequence[Any],
    timing_rows: Sequence[Any],
) -> float | None:
    """Return the imode=2 colour-probe setup/build cost outside the timed loop."""

    probe_index: int | None = None
    probe_elapsed: float | None = None
    for index, command in enumerate(commands):
        args = _command_record_args(command)
        if args and str(args[0]) == "./amplicol_color_probe":
            elapsed = _command_record_elapsed_s(command)
            if elapsed is not None:
                probe_index = index
                probe_elapsed = elapsed
    if probe_index is None or probe_elapsed is None:
        return None
    build_elapsed = 0.0
    for command in commands[:probe_index]:
        args = _command_record_args(command)
        if (
            args
            and str(args[0]) == "make"
            and any(str(arg) == "amplicol_color_probe" for arg in args[1:])
        ):
            elapsed = _command_record_elapsed_s(command)
            if elapsed is not None:
                build_elapsed = elapsed
    timed_total = _timing_total_seconds(timing_rows)
    if timed_total is None:
        return build_elapsed if build_elapsed > 0.0 else None
    setup_elapsed = max(0.0, probe_elapsed - timed_total)
    total = build_elapsed + setup_elapsed
    return total if total > 0.0 else None


def _command_record_args(command: Any) -> list[str]:
    if isinstance(command, dict):
        args = command.get("args", ())
    else:
        args = getattr(command, "args", ())
    return [str(arg) for arg in (args or ())]


def _command_record_elapsed_s(command: Any) -> float | None:
    value = (
        command.get("elapsed_s")
        if isinstance(command, dict)
        else getattr(command, "elapsed_s", None)
    )
    return _optional_float(value)


def _timing_total_seconds(timing_rows: Sequence[Any]) -> float | None:
    for row in timing_rows:
        label = row.get("label") if isinstance(row, dict) else getattr(row, "label", "")
        if str(label).strip().lower() != "total":
            continue
        value = row.get("seconds") if isinstance(row, dict) else getattr(row, "seconds", None)
        return _optional_float(value)
    return None


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
            max_sectors=-1,
        )
        return set(lc_line_pairing_representative_ids(color_plan))
    color_plan = build_color_plan(
        process,
        color_accuracy="lc",
        options=base.process_options(),
        max_sectors=-1,
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
    library_current: bool,
    color_accuracy: str,
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
        setup_commands: tuple[Any, ...] = ()
        if color_accuracy != "lc":
            if not library_current:
                build = adapter.prepare_library(
                    process,
                    process_file=process_file,
                    options=base.process_options(),
                    process_list_backend=backend,  # type: ignore[arg-type]
                    warmup_particles=point,
                    warmup_points=1,
                    raw=True,
                    color_complete=True,
                )
                process_file = str(build.process_file)
                setup_commands = build.commands
                adapter = AmplicolAdapter(
                    REPO_ROOT,
                    jobs=jobs,
                    timeout=_remaining(deadline_started, time_limit),
                )
            run = adapter.run_color_library_probe(
                process,
                color_accuracy=color_accuracy,
                particles=point,
                points=1,
                process_file=process_file,
                options=base.process_options(),
                process_list_backend=backend,  # type: ignore[arg-type]
                color_complete=True,
            )
            if run.first_point_matrix_element is None:
                raise RuntimeError("AmpliCol colour validation probe did not report a value")
            reference = float(run.first_point_matrix_element)
            point_order_source = "amplicol-color-library-probe-pdg-order"
        else:
            if not library_current:
                build = adapter.prepare_library(
                    process,
                    process_file=process_file,
                    options=base.process_options(),
                    process_list_backend=backend,  # type: ignore[arg-type]
                    warmup_particles=point,
                )
                process_file = str(build.process_file)
                setup_commands = build.commands
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
                    "amplicol_command": run.commands[-1].args if run.commands else (),
                    "finished_at": _now(),
                }
            )
            return payload
        rel_diffs = {
            key: _relative_difference(reference, value)
            for key, value in values.items()
        }
        abs_diffs = {key: abs(reference - value) for key, value in values.items()}
        all_flow_reference = _optional_float(amplicol.get("all_flow_reference_value"))
        all_flow_value = _optional_float(_mode(case, "pyamplicol_jit").get("all_flow_value"))
        all_flow_validate_value = bool(
            amplicol.get("all_flow_value_validation_enabled")
        ) and bool(_mode(case, "pyamplicol_jit").get("all_flow_value_validation_enabled"))
        all_flow_rel_diff = (
            None
            if (
                not all_flow_validate_value
                or all_flow_reference is None
                or all_flow_value is None
            )
            else _relative_difference(all_flow_reference, all_flow_value)
        )
        all_flow_abs_diff = (
            None
            if (
                not all_flow_validate_value
                or all_flow_reference is None
                or all_flow_value is None
            )
            else abs(all_flow_reference - all_flow_value)
        )
        max_rel_diff = max(rel_diffs.values(), default=0.0)
        max_abs_diff = max(abs_diffs.values(), default=0.0)
        if all_flow_rel_diff is not None:
            max_rel_diff = max(max_rel_diff, all_flow_rel_diff)
        if all_flow_abs_diff is not None:
            max_abs_diff = max(max_abs_diff, all_flow_abs_diff)
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
                "all_flow_reference": all_flow_reference,
                "all_flow_value": all_flow_value,
                "all_flow_value_validation_enabled": all_flow_validate_value,
                "all_flow_validation_note": amplicol.get(
                    "all_flow_validation_note"
                ),
                "all_flow_relative_difference": all_flow_rel_diff,
                "all_flow_absolute_difference": all_flow_abs_diff,
                "max_relative_difference": max_rel_diff,
                "max_absolute_difference": max_abs_diff,
                "tolerance": VALIDATION_REL_TOL,
                "absolute_tolerance": VALIDATION_ABS_TOL,
                "point_source": point_source,
                "point_order_source": point_order_source,
                "amplicol_command": run.commands[-1].args if run.commands else (),
                "prepared_library_for_validation": not library_current,
                "commands": [
                    {
                        "args": command.args,
                        "elapsed_s": command.elapsed_s,
                        "returncode": command.returncode,
                    }
                    for command in (*setup_commands, *run.commands)
                ],
                "finished_at": _now(),
            }
        )
    except Exception as exc:  # noqa: BLE001 - benchmark harness records failures.
        payload.update(_error_payload(exc))
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
    selected_jit_opt_level: int,
    selected_lc_sector_ids: set[int] | None,
    reference_color_order: Sequence[int] | None,
    matrix_settings: dict[str, Any],
    color_accuracy: str,
    model_source: str,
    model_profile: str,
    reuse_lc_all_flow: bool = False,
    existing_mode: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    selected_batch_size = (
        DEFAULT_BATCH_SIZE if color_accuracy == "lc" else DEFAULT_NON_LC_BATCH_SIZE
    )
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
        "--color-accuracy",
        color_accuracy,
        "--batch-size",
        str(selected_batch_size),
        "--symbolica-n-cores",
        str(n_cores),
        "--json",
        "--monitor",
        "--model",
        model_source,
        "--symbolica-output-chunk-size",
        str(DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE),
        "--symbolica-output-chunk-strategy",
        "auto" if color_accuracy == "lc" else "uniform",
    ]
    if (
        model_profile != "built-in-sm"
        and base.key in QCD_ONLY_EXTERNAL_SM_FAMILIES
    ):
        generate.extend(["--max-coupling-order", "QED=0"])
    if model_profile != "built-in-sm":
        # AmpliCol's matrix workloads retain the hierarchy-preferred lowest
        # perturbative order. Derive the corresponding cap from model metadata
        # rather than encoding a process-family order table here.
        generate.extend(["--coupling-order-policy", "minimal"])
        model_cache_dir = (output_dir / ".model_cache").resolve()
        generate.extend(["--model-cache-dir", str(model_cache_dir)])
    if DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT:
        generate.append("--symbolica-stage-local-parameter-layout")
    else:
        generate.append("--no-symbolica-stage-local-parameter-layout")
    if not _matrix_uses_numerical_current_passes(
        base,
        color_accuracy=color_accuracy,
    ):
        generate.extend(
            [
                "--no-numerical-filter-current",
                "--no-numerical-current-merging",
            ]
        )
    fixed_helicity_all_flow_supported = (
        backend_key == "jit"
        and color_accuracy == "lc"
        and _lc_fixed_helicity_all_flow_supported(process, base)
    )
    fixed_helicity = (
        _fixed_helicity_choice(process, base)
        if fixed_helicity_all_flow_supported
        else None
    )
    selected_generate = list(generate)
    all_flow_generate = list(generate)
    if fixed_helicity_all_flow_supported:
        _replace_command_option(
            all_flow_generate,
            "--symbolica-output-chunk-size",
            str(DEFAULT_ALL_FLOW_SYMBOLICA_OUTPUT_CHUNK_SIZE),
        )
        _replace_command_option(
            all_flow_generate,
            "--batch-size",
            str(DEFAULT_ALL_FLOW_BATCH_SIZE),
        )
        _replace_command_option(
            all_flow_generate,
            "--symbolica-output-chunk-strategy",
            "uniform",
        )
    if reference_color_order:
        selected_generate.extend(
            [
                "--reference-color-order",
                ",".join(str(label) for label in reference_color_order),
            ]
        )
    if selected_lc_sector_ids:
        selected_generate.extend(
            [
                "--lc-sector-ids",
                ",".join(str(sector_id) for sector_id in sorted(selected_lc_sector_ids)),
            ]
        )
    run_lc_all_flow_generation = (
        backend_key == "jit"
        and color_accuracy == "lc"
        and fixed_helicity_all_flow_supported
    )
    if backend_key == "jit" and color_accuracy == "lc":
        selected_generate.append("--skip-generic-plan")
        selected_generate.append("--no-runtime-lc-sector-selector")
    if run_lc_all_flow_generation:
        all_flow_generate.append("--skip-generic-plan")
        all_flow_generate.append("--no-runtime-lc-sector-selector")
        if fixed_helicity is not None:
            fixed_helicity_flags = [
                "--source-helicities",
                str(fixed_helicity["source_helicities_cli"]),
            ]
            all_flow_generate.extend(fixed_helicity_flags)
    if backend_key == "jit":
        selected_evaluator_flags = [
            "--symbolica-evaluator-backend",
            "jit",
            "--symbolica-jit-optimization-level",
            str(selected_jit_opt_level),
            "--symbolica-iterations",
            str(DEFAULT_SYMBOLICA_ITERATIONS),
            "--symbolica-max-horner-scheme-variables",
            "1000",
            "--symbolica-max-common-pair-cache-entries",
            "5000000",
            "--symbolica-max-common-pair-distance",
            "1000",
            "--symbolica-compiled-chunk-compile-workers",
            str(max(1, n_cores)),
        ]
        all_flow_evaluator_flags = list(selected_evaluator_flags)
        _replace_command_option(
            all_flow_evaluator_flags,
            "--symbolica-jit-optimization-level",
            "1",
        )
    else:
        selected_evaluator_flags = [
            "--symbolica-evaluator-backend",
            "compiled-complex",
            "--symbolica-compiled-preset",
            "runtime-o3",
            "--symbolica-compiled-chunk-compile-workers",
            str(max(1, n_cores)),
        ]
        all_flow_evaluator_flags = list(selected_evaluator_flags)
    selected_generate.extend(selected_evaluator_flags)
    all_flow_generate.extend(all_flow_evaluator_flags)
    selected_output_dir = (
        output_dir / "selected_flow"
        if backend_key == "jit" and color_accuracy == "lc"
        else output_dir
    )
    all_flow_output_dir = output_dir / "all_flows"
    selected_generate.extend([process, str(selected_output_dir)])
    all_flow_generate.extend([process, str(all_flow_output_dir)])
    existing_mode = dict(existing_mode or {})
    existing_all_flow_dir = existing_mode.get("all_flow_output_dir")
    reuse_existing_all_flow = bool(
        reuse_lc_all_flow
        and run_lc_all_flow_generation
        and isinstance(existing_all_flow_dir, str)
        and existing_all_flow_dir
        and Path(existing_all_flow_dir).is_dir()
        and _optional_float(existing_mode.get("all_flow_generation_s")) is not None
        and _optional_float(existing_mode.get("all_flow_wall_us_per_point")) is not None
    )
    reuse_existing_all_flow_resource_status = bool(
        reuse_lc_all_flow
        and run_lc_all_flow_generation
        and str(existing_mode.get("all_flow_status", "")) == "ram_limit"
    )
    if reuse_lc_all_flow and run_lc_all_flow_generation:
        print(
            (
                f"[matrix] reusing LC all-flow artifact {existing_all_flow_dir}"
                if reuse_existing_all_flow
                else (
                    "[matrix] retaining known LC all-flow RAM-limit status"
                    if reuse_existing_all_flow_resource_status
                    else "[matrix] requested LC all-flow reuse is unavailable; "
                    "regenerating it"
                )
            ),
            flush=True,
        )

    try:
        _preserve_generated_output(selected_output_dir)
        selected_gen = _run_generation_json_command(
            selected_generate,
            timeout=_remaining(deadline_started, time_limit),
            log_path=_progress_log_path(selected_output_dir, "generate"),
            output_dir=selected_output_dir,
            median_samples=3,
        )
        all_flow_gen: dict[str, Any] | None = None
        all_flow_failure: dict[str, Any] | None = None
        all_flow_replay_fallback = False
        if (
            run_lc_all_flow_generation
            and not reuse_existing_all_flow
            and not reuse_existing_all_flow_resource_status
        ):
            _preserve_generated_output(all_flow_output_dir)
            try:
                all_flow_gen = _run_generation_json_command(
                    all_flow_generate,
                    timeout=_remaining(deadline_started, time_limit),
                    log_path=_progress_log_path(all_flow_output_dir, "generate"),
                    output_dir=all_flow_output_dir,
                    median_samples=3,
                )
            except MemoryLimitExceeded as exc:
                all_flow_failure = _error_payload(exc)
        timed = _run_json_command(
            [
                sys.executable,
                "-m",
                "pyamplicol",
                "time-process",
                "--target-runtime",
                str(target_runtime),
                "--batch-size",
                str(selected_batch_size),
                "--json",
                str(selected_output_dir),
            ],
            timeout=None,
            log_path=_progress_log_path(selected_output_dir, "time"),
        )
        all_flow_timed: dict[str, Any] | None = None
        if (
            fixed_helicity_all_flow_supported
            and not reuse_existing_all_flow
            and not reuse_existing_all_flow_resource_status
            and all_flow_failure is None
        ):
            try:
                all_flow_timed = _run_json_command(
                    [
                        sys.executable,
                        "-m",
                        "pyamplicol",
                        "time-process",
                        "--target-runtime",
                        str(target_runtime),
                        "--batch-size",
                        str(DEFAULT_ALL_FLOW_BATCH_SIZE),
                        "--json",
                        str(all_flow_output_dir),
                    ],
                    timeout=None,
                    log_path=_progress_log_path(all_flow_output_dir, "time_all_flows"),
                )
            except MemoryLimitExceeded as exc:
                all_flow_failure = _error_payload(exc)
        gen = all_flow_gen if all_flow_gen is not None else selected_gen
        profile = timed.get("profile", {})
        all_flow_profile = (
            all_flow_timed.get("profile", {}) if isinstance(all_flow_timed, dict) else {}
        )
        selected_generation_s = _optional_float(
            selected_gen.get("_command_elapsed_s", selected_gen.get("generation_s"))
        )
        all_flow_generation_s = (
            None
            if all_flow_gen is None
            else _optional_float(
                all_flow_gen.get("_command_elapsed_s", all_flow_gen.get("generation_s"))
            )
        )
        displayed_generation_s = (
            all_flow_generation_s
            if all_flow_generation_s is not None
            else selected_generation_s
        )
        lowering_status = gen.get("lowering_status")
        if not isinstance(lowering_status, dict):
            lowering_status = {}
        payload.update(
            {
                "status": "ok",
                "generation_s": displayed_generation_s,
                "selected_generation_s": selected_generation_s,
                "all_flow_generation_s": all_flow_generation_s,
                "internal_generation_s": _optional_float(gen.get("generation_s")),
                "selected_internal_generation_s": _optional_float(
                    selected_gen.get("generation_s")
                ),
                "all_flow_internal_generation_s": (
                    None
                    if all_flow_gen is None
                    else _optional_float(all_flow_gen.get("generation_s"))
                ),
                "jit_compile_s": _optional_float(gen.get("jit_compile_s")),
                "jit_fraction_of_generation": _optional_float(
                    gen.get("jit_fraction_of_generation")
                ),
                "current_count": lowering_status.get("current_count"),
                "interaction_count": lowering_status.get("interaction_count"),
                "amplitude_root_count": lowering_status.get("amplitude_root_count"),
                "amplitude_color_sector_count": lowering_status.get(
                    "amplitude_color_sector_count"
                ),
                "internal_current_color_sector_count": lowering_status.get(
                    "internal_current_color_sector_count"
                ),
                "runtime_us_per_point": _optional_float(
                    profile.get("core_evaluator_us_per_point")
                ),
                "wall_us_per_point": _optional_float(
                    profile.get("wall_us_per_point")
                ),
                "samples": profile.get("samples"),
                "all_flow_runtime_us_per_point": _optional_float(
                    all_flow_profile.get("core_evaluator_us_per_point")
                ),
                "all_flow_wall_us_per_point": _optional_float(
                    all_flow_profile.get("wall_us_per_point")
                ),
                "all_flow_samples": all_flow_profile.get("samples"),
                "all_flow_status": (
                    "ok"
                    if reuse_existing_all_flow
                    else (
                        str(existing_mode.get("all_flow_status"))
                        if reuse_existing_all_flow_resource_status
                        else (
                            str(all_flow_failure.get("status"))
                            if all_flow_failure is not None
                            else "ok"
                            if all_flow_timed is not None
                            else None
                        )
                    )
                ),
                "all_flow_error": (
                    existing_mode.get("all_flow_error")
                    if reuse_existing_all_flow_resource_status
                    else (
                        all_flow_failure.get("error")
                        if all_flow_failure is not None
                        else None
                    )
                ),
                "all_flow_memory_limit_gb": (
                    existing_mode.get("all_flow_memory_limit_gb")
                    if reuse_existing_all_flow_resource_status
                    else (
                        all_flow_failure.get("memory_limit_gb")
                        if all_flow_failure is not None
                        else None
                    )
                ),
                "all_flow_peak_memory_gb": (
                    existing_mode.get("all_flow_peak_memory_gb")
                    if reuse_existing_all_flow_resource_status
                    else (
                        all_flow_failure.get("peak_memory_gb")
                        if all_flow_failure is not None
                        else None
                    )
                ),
                "all_flow_memory_metric": (
                    existing_mode.get("all_flow_memory_metric")
                    if reuse_existing_all_flow_resource_status
                    else (
                        all_flow_failure.get("memory_metric")
                        if all_flow_failure is not None
                        else None
                    )
                ),
                "all_flow_time_batch_size": (
                    DEFAULT_ALL_FLOW_BATCH_SIZE
                    if all_flow_timed is not None
                    else None
                ),
                "all_flow_symbolica_output_chunk_size": (
                    DEFAULT_ALL_FLOW_SYMBOLICA_OUTPUT_CHUNK_SIZE
                    if all_flow_gen is not None
                    else None
                ),
                "all_flow_value": _payload_first_value(all_flow_timed or {}),
                "selected_jit_optimization_level": (
                    selected_jit_opt_level if backend_key == "jit" else None
                ),
                "all_flow_jit_optimization_level": (
                    1
                    if backend_key == "jit"
                    and color_accuracy == "lc"
                    and fixed_helicity_all_flow_supported
                    else None
                ),
                "all_flow_replay_fallback": all_flow_replay_fallback,
                "all_flow_helicity_mode": (
                    None if fixed_helicity is None else fixed_helicity["mode"]
                ),
                "all_flow_source_helicities": (
                    None
                    if fixed_helicity is None
                    else {
                        str(label): int(helicity)
                        for label, helicity in fixed_helicity[
                            "source_helicities"
                        ].items()
                    }
                ),
                "all_flow_amplicol_helicities": (
                    None
                    if fixed_helicity is None
                    else list(fixed_helicity["amplicol_helicities"])
                ),
                "all_flow_value_validation_enabled": (
                    False
                    if fixed_helicity is None
                    else bool(fixed_helicity["value_validation_enabled"])
                    and not all_flow_replay_fallback
                ),
                "all_flow_validation_note": (
                    (
                        fixed_helicity["validation_note"]
                        if not all_flow_replay_fallback
                        else (
                            f"{fixed_helicity['validation_note']}; "
                            "all-sector fallback scalar is timing-only"
                        )
                    )
                    if fixed_helicity is not None
                    else (
                        "AmpliCol fixed-helicity all-flow probe is unsupported "
                        "for more than two quark lines"
                    )
                    if backend_key == "jit" and color_accuracy == "lc"
                    else None
                ),
                "selected_lc_sector_ids": (
                    None
                    if selected_lc_sector_ids is None
                    else sorted(selected_lc_sector_ids)
                ),
                "runtime_lc_sector_ids": (
                    None
                    if selected_lc_sector_ids is None
                    else sorted(selected_lc_sector_ids)
                ),
                "reference_color_order": (
                    None
                    if reference_color_order is None
                    else [int(label) for label in reference_color_order]
                ),
                "selected_output_dir": str(selected_output_dir),
                "all_flow_output_dir": (
                    str(existing_all_flow_dir)
                    if reuse_existing_all_flow_resource_status
                    and isinstance(existing_all_flow_dir, str)
                    else str(all_flow_output_dir)
                    if all_flow_gen is not None or all_flow_failure is not None
                    else None
                ),
                "generate_payload": _compact_payload(gen),
                "selected_generate_payload": _compact_payload(selected_gen),
                "all_flow_generate_payload": (
                    None if all_flow_gen is None else _compact_payload(all_flow_gen)
                ),
                "time_payload": _compact_payload(timed),
                "all_flow_time_payload": (
                    None if all_flow_timed is None else _compact_payload(all_flow_timed)
                ),
                "generate_log": gen.get("_progress_log"),
                "selected_generate_log": selected_gen.get("_progress_log"),
                "all_flow_generate_log": (
                    None if all_flow_gen is None else all_flow_gen.get("_progress_log")
                ),
                "time_log": timed.get("_progress_log"),
                "all_flow_time_log": (
                    None if all_flow_timed is None else all_flow_timed.get("_progress_log")
                ),
                "matrix_settings": matrix_settings,
                "finished_at": _now(),
            }
        )
        if reuse_existing_all_flow:
            _restore_reused_lc_all_flow_fields(payload, existing_mode)
        elif reuse_existing_all_flow_resource_status:
            _restore_reused_lc_all_flow_fields(
                payload,
                existing_mode,
                include_primary_fields=False,
            )
    except Exception as exc:  # noqa: BLE001 - benchmark harness records failures.
        payload.update(_error_payload(exc))
        if reuse_existing_all_flow:
            _restore_reused_lc_all_flow_fields(payload, existing_mode)
            payload["all_flow_status"] = "ok"
        if _color_accuracy_error_is_structural_unsupported(str(payload.get("error", ""))):
            payload["status"] = "unsupported"
        limitation = _documented_backend_limitation_from_error(
            str(payload.get("error", "")),
            backend_key=backend_key,
        )
        if limitation is not None:
            payload["status"] = "backend_unsupported"
            payload["backend_limitation"] = limitation
    return payload


def _retime_pyamplicol_case(
    existing_mode: Mapping[str, Any],
    *,
    backend_key: str,
    output_dir: Path,
    target_runtime: float,
    color_accuracy: str,
) -> dict[str, Any]:
    payload = dict(existing_mode)
    if payload.get("status") != "ok":
        payload["retime_status"] = "skipped_non_ok"
        print(
            f"[matrix] skip retime {output_dir}: existing status "
            f"{payload.get('status', 'missing')}",
            flush=True,
        )
        return payload

    selected_default = (
        output_dir / "selected_flow"
        if backend_key == "jit" and color_accuracy == "lc"
        else output_dir
    )
    selected_output_dir = _cached_artifact_path(
        payload.get("selected_output_dir") or payload.get("output_dir"),
        fallback=selected_default,
    )
    if not (selected_output_dir / "process_manifest.json").is_file():
        payload.update(
            {
                "retime_status": "artifact_missing",
                "retime_error": f"missing artifact: {selected_output_dir}",
                "retime_attempted_at": _now(),
            }
        )
        print(f"[matrix] cannot retime: {payload['retime_error']}", flush=True)
        return payload

    settings = payload.get("matrix_settings")
    settings = dict(settings) if isinstance(settings, dict) else {}
    selected_batch_size = _optional_int(
        settings.get("selected_runtime_batch_size", settings.get("batch_size"))
    ) or (DEFAULT_BATCH_SIZE if color_accuracy == "lc" else DEFAULT_NON_LC_BATCH_SIZE)
    selected_log = _progress_log_path(selected_output_dir, "time")
    _preserve_timing_log(selected_log)
    print(
        f"[matrix] retime existing {selected_output_dir} "
        f"(batch {selected_batch_size}, target {target_runtime:g}s)",
        flush=True,
    )
    try:
        timed = _time_existing_pyamplicol_artifact(
            selected_output_dir,
            target_runtime=target_runtime,
            batch_size=selected_batch_size,
            log_path=selected_log,
        )
    except Exception as exc:  # noqa: BLE001 - keep the last valid table value.
        payload.update(
            {
                "retime_status": "error",
                "retime_error": str(exc),
                "retime_attempted_at": _now(),
            }
        )
        print(f"[matrix] retime failed for {selected_output_dir}: {exc}", flush=True)
        return payload

    profile = timed.get("profile", {})
    if not isinstance(profile, dict):
        profile = {}
    payload.update(
        {
            "selected_output_dir": str(selected_output_dir),
            "runtime_us_per_point": _optional_float(
                profile.get("core_evaluator_us_per_point")
            ),
            "wall_us_per_point": _optional_float(profile.get("wall_us_per_point")),
            "samples": profile.get("samples"),
            "time_payload": _compact_payload(timed),
            "time_log": timed.get("_progress_log"),
            "retime_status": "ok",
            "retimed_at": _now(),
        }
    )
    settings["target_runtime_s"] = float(target_runtime)
    settings["selected_runtime_batch_size"] = int(selected_batch_size)
    payload["matrix_settings"] = settings
    selected_validation = _selected_refresh_validation(existing_mode, payload)
    payload["retime_validation"] = selected_validation
    if selected_validation.get("status") == "failed":
        payload["status"] = "validation_failed"
        payload["retime_status"] = "validation_failed"
        payload["error"] = (
            "retiming existing artifact changed its selected-flow value: "
            f"relative difference {selected_validation.get('relative_difference')}"
        )
        return payload

    should_retime_all_flow = (
        backend_key == "jit"
        and color_accuracy == "lc"
        and str(payload.get("all_flow_status", "ok")) not in {"ram_limit", "timeout"}
        and _optional_float(payload.get("all_flow_wall_us_per_point")) is not None
    )
    if not should_retime_all_flow:
        return payload

    all_flow_output_dir = _cached_artifact_path(
        payload.get("all_flow_output_dir"),
        fallback=output_dir / "all_flows",
    )
    payload["all_flow_output_dir"] = str(all_flow_output_dir)
    if not (all_flow_output_dir / "process_manifest.json").is_file():
        payload["all_flow_retime_status"] = "artifact_missing"
        payload["all_flow_retime_error"] = f"missing artifact: {all_flow_output_dir}"
        print(
            f"[matrix] cannot retime all-flow: {payload['all_flow_retime_error']}",
            flush=True,
        )
        return payload

    all_flow_batch_size = _optional_int(
        settings.get("all_flow_runtime_batch_size")
    ) or DEFAULT_ALL_FLOW_BATCH_SIZE
    all_flow_log = _progress_log_path(all_flow_output_dir, "time_all_flows")
    _preserve_timing_log(all_flow_log)
    print(
        f"[matrix] retime existing all-flow {all_flow_output_dir} "
        f"(batch {all_flow_batch_size}, target {target_runtime:g}s)",
        flush=True,
    )
    try:
        all_flow_timed = _time_existing_pyamplicol_artifact(
            all_flow_output_dir,
            target_runtime=target_runtime,
            batch_size=all_flow_batch_size,
            log_path=all_flow_log,
        )
    except Exception as exc:  # noqa: BLE001 - selected timing remains valid.
        payload["all_flow_retime_status"] = "error"
        payload["all_flow_retime_error"] = str(exc)
        print(f"[matrix] all-flow retime failed for {all_flow_output_dir}: {exc}", flush=True)
        return payload

    all_flow_profile = all_flow_timed.get("profile", {})
    if not isinstance(all_flow_profile, dict):
        all_flow_profile = {}
    payload.update(
        {
            "all_flow_runtime_us_per_point": _optional_float(
                all_flow_profile.get("core_evaluator_us_per_point")
            ),
            "all_flow_wall_us_per_point": _optional_float(
                all_flow_profile.get("wall_us_per_point")
            ),
            "all_flow_samples": all_flow_profile.get("samples"),
            "all_flow_time_batch_size": int(all_flow_batch_size),
            "all_flow_value": _payload_first_value(all_flow_timed),
            "all_flow_time_payload": _compact_payload(all_flow_timed),
            "all_flow_time_log": all_flow_timed.get("_progress_log"),
            "all_flow_retime_status": "ok",
        }
    )
    settings["all_flow_runtime_batch_size"] = int(all_flow_batch_size)
    payload["matrix_settings"] = settings
    all_flow_validation = _all_flow_retime_validation(existing_mode, payload)
    payload["all_flow_retime_validation"] = all_flow_validation
    if all_flow_validation.get("status") == "failed":
        payload["status"] = "validation_failed"
        payload["all_flow_retime_status"] = "validation_failed"
        payload["error"] = (
            "retiming existing artifact changed its all-flow value: relative "
            f"difference {all_flow_validation.get('relative_difference')}"
        )
    return payload


def _time_existing_pyamplicol_artifact(
    output_dir: Path,
    *,
    target_runtime: float,
    batch_size: int,
    log_path: Path,
) -> dict[str, Any]:
    return _run_json_command(
        [
            sys.executable,
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
        log_path=log_path,
    )


def _cached_artifact_path(raw: object, *, fallback: Path) -> Path:
    if raw is None or not str(raw).strip():
        return fallback
    path = Path(str(raw))
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def _optional_int(value: object) -> int | None:
    if value is None or isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _all_flow_retime_validation(
    previous_mode: Mapping[str, Any],
    refreshed_mode: Mapping[str, Any],
) -> dict[str, Any]:
    previous_value = _optional_float(previous_mode.get("all_flow_value"))
    if previous_value is None:
        previous_payload = previous_mode.get("all_flow_time_payload")
        previous_value = _payload_first_value(
            previous_payload if isinstance(previous_payload, dict) else {}
        )
    refreshed_value = _optional_float(refreshed_mode.get("all_flow_value"))
    if previous_value is None or refreshed_value is None:
        return {
            "status": "unavailable",
            "previous_value": previous_value,
            "refreshed_value": refreshed_value,
            "relative_tolerance": 1.0e-10,
        }
    relative_difference = _relative_difference(previous_value, refreshed_value)
    return {
        "status": "ok" if relative_difference <= 1.0e-10 else "failed",
        "previous_value": previous_value,
        "refreshed_value": refreshed_value,
        "relative_difference": relative_difference,
        "relative_tolerance": 1.0e-10,
    }


def _preserve_timing_log(log_path: Path) -> None:
    if not log_path.exists():
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    preserved = log_path.with_name(f"{log_path.name}.before_{stamp}")
    counter = 1
    while preserved.exists():
        preserved = log_path.with_name(f"{log_path.name}.before_{stamp}_{counter}")
        counter += 1
    log_path.rename(preserved)


def _preserve_generated_output(output_dir: Path) -> Path | None:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    preserved_output: Path | None = None
    if output_dir.exists():
        preserved = output_dir.with_name(f"{output_dir.name}_before_{stamp}")
        counter = 1
        while preserved.exists():
            preserved = output_dir.with_name(
                f"{output_dir.name}_before_{stamp}_{counter}"
            )
            counter += 1
        output_dir.rename(preserved)
        preserved_output = preserved
        print(f"[matrix] preserved {output_dir} as {preserved}", flush=True)
    for suffix in ("generate.log", "time.log", "time_all_flows.log"):
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
    return preserved_output


def _restore_reused_lc_all_flow_fields(
    payload: dict[str, Any],
    existing_mode: Mapping[str, Any],
    *,
    include_primary_fields: bool = True,
) -> None:
    for key, value in existing_mode.items():
        if str(key).startswith("all_flow_"):
            payload[str(key)] = value
    if not include_primary_fields:
        return
    for key in (
        "generation_s",
        "internal_generation_s",
        "jit_compile_s",
        "jit_fraction_of_generation",
        "current_count",
        "interaction_count",
        "amplitude_root_count",
        "amplitude_color_sector_count",
        "internal_current_color_sector_count",
        "generate_payload",
        "generate_log",
    ):
        if key in existing_mode:
            payload[key] = existing_mode[key]


def _selected_refresh_validation(
    previous_mode: Mapping[str, Any],
    refreshed_mode: Mapping[str, Any],
) -> dict[str, Any]:
    previous_payload = previous_mode.get("time_payload")
    refreshed_payload = refreshed_mode.get("time_payload")
    previous_value = _payload_first_value(
        previous_payload if isinstance(previous_payload, dict) else {}
    )
    refreshed_value = _payload_first_value(
        refreshed_payload if isinstance(refreshed_payload, dict) else {}
    )
    if previous_value is None or refreshed_value is None:
        return {
            "status": "unavailable",
            "previous_value": previous_value,
            "refreshed_value": refreshed_value,
            "relative_tolerance": 1.0e-10,
        }
    absolute_difference = abs(refreshed_value - previous_value)
    scale = max(abs(previous_value), abs(refreshed_value), 1.0e-300)
    relative_difference = absolute_difference / scale
    return {
        "status": "ok" if relative_difference <= 1.0e-10 else "failed",
        "previous_value": previous_value,
        "refreshed_value": refreshed_value,
        "absolute_difference": absolute_difference,
        "relative_difference": relative_difference,
        "relative_tolerance": 1.0e-10,
    }


def _pyamplicol_output_dir(
    output_root: Path,
    *,
    color_accuracy: str,
    base: BaseProcess,
    n_final: int,
    backend_key: str,
) -> Path:
    color = color_accuracy.lower()
    root = output_root if color == "lc" else output_root / color
    return root / base.key / f"n{n_final}" / backend_key


def _progress_log_path(output_dir: Path, phase: str) -> Path:
    return output_dir.with_name(f"{output_dir.name}.{phase}.log")


def _replace_command_option(command: list[str], option: str, value: str) -> None:
    try:
        index = command.index(option)
    except ValueError:
        command.extend([option, value])
        return
    if index + 1 >= len(command):
        command.append(value)
    else:
        command[index + 1] = value


def _run_generation_json_command(
    args: Sequence[str],
    *,
    timeout: float | None,
    log_path: Path,
    output_dir: Path,
    median_samples: int,
    median_max_s: float = 90.0,
) -> dict[str, Any]:
    requested_samples = max(1, int(median_samples))
    payloads = [
        _run_json_command(
            args,
            timeout=timeout,
            log_path=log_path,
        )
    ]
    artifact_dirs: list[str | None] = []
    first_elapsed = _optional_float(payloads[0].get("_command_elapsed_s"))
    effective_samples = (
        requested_samples
        if first_elapsed is not None and first_elapsed < float(median_max_s)
        else 1
    )
    for _sample_index in range(1, effective_samples):
        preserved = _preserve_generated_output(output_dir)
        artifact_dirs.append(None if preserved is None else str(preserved))
        payloads.append(
            _run_json_command(
                args,
                timeout=timeout,
                log_path=log_path,
            )
        )
    artifact_dirs.append(str(output_dir))
    median_index = sorted(
        range(len(payloads)),
        key=lambda index: float(payloads[index]["_command_elapsed_s"]),
    )[len(payloads) // 2]
    selected = dict(payloads[median_index])
    artifact_payload = payloads[-1]
    selected["_generation_samples"] = [
        {
            "sample_index": index + 1,
            "artifact_dir": artifact_dirs[index],
            "command_elapsed_s": _optional_float(
                payload.get("_command_elapsed_s")
            ),
            "internal_generation_s": _optional_float(payload.get("generation_s")),
            "jit_compile_s": _optional_float(payload.get("jit_compile_s")),
        }
        for index, payload in enumerate(payloads)
    ]
    selected["_generation_sample_count"] = len(payloads)
    selected["_generation_median_sample_index"] = median_index + 1
    selected["_artifact_sample_index"] = len(payloads)
    selected["_artifact_command_elapsed_s"] = _optional_float(
        artifact_payload.get("_command_elapsed_s")
    )
    selected["_artifact_internal_generation_s"] = _optional_float(
        artifact_payload.get("generation_s")
    )
    return selected


def _run_json_command(
    args: Sequence[str],
    *,
    timeout: float | None,
    log_path: Path | None = None,
) -> dict[str, Any]:
    env = dict(os.environ)
    env["PYTHONPATH"] = (
        f"{SRC_DIR}{os.pathsep}{env['PYTHONPATH']}"
        if env.get("PYTHONPATH")
        else str(SRC_DIR)
    )
    start = time.perf_counter()
    stdout_text = ""
    stderr_text = ""
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"[matrix] progress log: {log_path}", flush=True)
        with log_path.open("w", encoding="utf-8") as log_file:
            log_file.write(f"# started: {_now()}\n")
            log_file.write("# command: " + " ".join(str(arg) for arg in args) + "\n")
            log_file.flush()
            completed, stdout_text, stderr_text = _run_process_group_command(
                list(args),
                timeout=timeout,
                env=env,
                stderr=log_file,
            )
        stderr_text = stderr_text or _read_tail(log_path)
    else:
        completed, stdout_text, stderr_text = _run_process_group_command(
            list(args),
            timeout=timeout,
            env=env,
            stderr=subprocess.PIPE,
        )
    elapsed_s = time.perf_counter() - start
    if completed.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{completed.returncode}: {' '.join(args)}\n"
            f"stdout:\n{stdout_text[-4000:]}\n"
            f"stderr:\n{stderr_text[-4000:]}"
        )
    payload = _parse_json_output(stdout_text)
    payload["_command_elapsed_s"] = elapsed_s
    payload["_command_args"] = list(args)
    if log_path is not None:
        payload["_progress_log"] = str(log_path)
    return payload


def _run_process_group_command(
    args: list[str],
    *,
    timeout: float | None,
    env: dict[str, str],
    stderr: Any,
) -> tuple[subprocess.CompletedProcess[str], str, str]:
    started = time.perf_counter()
    deadline = None if timeout is None else started + timeout
    memory_limit_bytes = int(DEFAULT_MEMORY_LIMIT_GB * 1024**3)
    memory_high_watermark_bytes = int(
        MEMORY_HIGH_WATERMARK_FRACTION * memory_limit_bytes
    )
    memory_monitor = ProcessTreeMemoryMonitor(
        physical_footprint_near_limit_poll_s=MEMORY_NEAR_LIMIT_POLL_S,
        physical_footprint_high_watermark_bytes=memory_high_watermark_bytes,
    )
    last_heartbeat = started
    process = subprocess.Popen(
        args,
        cwd=REPO_ROOT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=stderr,
        start_new_session=True,
    )
    try:
        while True:
            memory_sample = memory_monitor.sample(process.pid)
            memory_bytes = memory_sample.effective_bytes
            if memory_bytes is not None:
                if memory_bytes > memory_limit_bytes:
                    _terminate_process_group(process.pid, signal.SIGTERM)
                    try:
                        process.communicate(timeout=10.0)
                    except subprocess.TimeoutExpired:
                        _terminate_process_group(process.pid, signal.SIGKILL)
                        process.communicate()
                    raise MemoryLimitExceeded(
                        limit_gb=DEFAULT_MEMORY_LIMIT_GB,
                        peak_memory_bytes=memory_bytes,
                        memory_metric=memory_sample.effective_metric or "memory",
                        peak_rss_bytes=memory_monitor.peak_rss_bytes,
                        peak_physical_footprint_bytes=(
                            memory_monitor.peak_physical_footprint_bytes
                        ),
                        command=args,
                    )
            wait_s = _memory_poll_interval(
                memory_bytes,
                high_watermark_bytes=memory_high_watermark_bytes,
            )
            if deadline is not None:
                remaining_s = deadline - time.perf_counter()
                if remaining_s <= 0:
                    raise subprocess.TimeoutExpired(args, timeout)
                wait_s = min(wait_s, remaining_s)
            try:
                stdout_text, stderr_text = process.communicate(timeout=wait_s)
                break
            except subprocess.TimeoutExpired:
                if deadline is not None and time.perf_counter() >= deadline:
                    raise
                elapsed_s = time.perf_counter() - started
                if time.perf_counter() - last_heartbeat >= COMMAND_HEARTBEAT_S:
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
                        f"[matrix] still running after {elapsed_s:.0f}s{memory_note}: "
                        f"{_command_label(args)}",
                        flush=True,
                    )
                    last_heartbeat = time.perf_counter()
    except KeyboardInterrupt:
        _terminate_process_group(process.pid, signal.SIGTERM)
        try:
            stdout_text, stderr_text = process.communicate(timeout=10.0)
        except subprocess.TimeoutExpired:
            _terminate_process_group(process.pid, signal.SIGKILL)
            stdout_text, stderr_text = process.communicate()
        raise
    except subprocess.TimeoutExpired as exc:
        _terminate_process_group(process.pid, signal.SIGTERM)
        try:
            stdout_text, stderr_text = process.communicate(timeout=10.0)
        except subprocess.TimeoutExpired:
            _terminate_process_group(process.pid, signal.SIGKILL)
            stdout_text, stderr_text = process.communicate()
        raise subprocess.TimeoutExpired(
            cmd=args,
            timeout=timeout,
            output=stdout_text,
            stderr=stderr_text if isinstance(stderr_text, str) else None,
        ) from exc
    return (
        subprocess.CompletedProcess(
            args,
            process.returncode,
            stdout_text or "",
            stderr_text or "",
        ),
        stdout_text or "",
        stderr_text or "",
    )


def _memory_poll_interval(
    memory_bytes: int | None,
    *,
    high_watermark_bytes: int,
) -> float:
    if memory_bytes is not None and memory_bytes >= high_watermark_bytes:
        return MEMORY_NEAR_LIMIT_POLL_S
    return MEMORY_POLL_S


def _command_label(args: Sequence[str]) -> str:
    parts = [str(part) for part in args]
    if "generate-process" in parts and len(parts) >= 2:
        return f"generate-process {parts[-2]} -> {Path(parts[-1]).name}"
    if "time-process" in parts and parts:
        return f"time-process {Path(parts[-1]).name}"
    return " ".join(parts[:4])


def _terminate_process_group(
    pid: int,
    sig: signal.Signals = signal.SIGTERM,
) -> None:
    try:
        os.killpg(pid, sig)
    except ProcessLookupError:
        return
    if sig != signal.SIGTERM:
        return
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline:
        try:
            finished, _status = os.waitpid(pid, os.WNOHANG)
        except ChildProcessError:
            return
        if finished:
            return
        time.sleep(0.1)
    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        return


def _read_tail(path: Path, *, max_bytes: int = 16000) -> str:
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(size - max_bytes, 0), os.SEEK_SET)
            return handle.read().decode("utf-8", errors="replace")
    except OSError:
        return ""


def render_latex_table(
    data: dict[str, Any],
    n_values: Sequence[int],
    *,
    columns_per_table: int = 3,
    color_accuracy: str = "lc",
) -> str:
    shown = tuple(sorted(dict.fromkeys(int(n) for n in n_values)))
    chunks = tuple(_chunks(shown, max(1, int(columns_per_table))))
    lines = [
        "% Generated by docs/result_matrix.py; edit result_matrix_data.json instead.",
        r"\providecommand{\matrixpunct}[1]{\textcolor{black}{\texttt{#1}}}",
        r"\providecommand{\matrixratio}[2]{\matrixpunct{(}\textcolor{#1}{\texttt{x#2}}\matrixpunct{)}}",
        r"\providecommand{\matrixevalratio}[1]{\matrixpunct{(}\textcolor{black}{\texttt{x#1}}\matrixpunct{)}}",
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
            r"\begin{tabular}[t]{@{}l@{\hspace{0.012in}}l@{\hspace{0.012in}}l@{}}"
            r"\matrixslot{1.30in}{#1}&\matrixslot{0.53in}{#2}&"
            r"\matrixslot{0.53in}{#3}\\"
            r"\matrixslot{1.30in}{#4}&\matrixslot{0.53in}{#5}&"
            r"\matrixslot{0.53in}{#6}"
            r"\end{tabular}}"
        ),
        (
            r"\providecommand{\matrixcellnonlc}[6]{"
            r"\begin{tabular}[t]{@{}l@{\hspace{0.012in}}l@{\hspace{0.012in}}l@{}}"
            r"\matrixslot{0.52in}{#1}&\matrixslot{0.94in}{#2}&"
            r"\matrixslot{0.94in}{#3}\\"
            r"\matrixslot{0.52in}{#4}&\matrixslot{0.94in}{#5}&"
            r"\matrixslot{0.94in}{#6}"
            r"\end{tabular}}"
        ),
        r"\providecommand{\matrixrefslot}[1]{\makebox[0.59in][l]{#1}}",
        (
            r"\providecommand{\matrixrefpair}[2]{"
            r"\begin{tabular}[t]{@{}l@{\hspace{0.012in}\matrixpunct{/}\hspace{0.012in}}l@{}}"
            r"\matrixrefslot{#1}&\matrixrefslot{#2}"
            r"\end{tabular}}"
        ),
        (
            r"\providecommand{\matrixsummarycell}[2]{"
            r"\begin{tabular}[t]{@{}l@{}}#1\\#2\end{tabular}}"
        ),
        r"\providecommand{\matrixsummaryfield}[1]{\makebox[0.41in][r]{#1}}",
        (
            r"\providecommand{\matrixsummaryfour}[4]{"
            r"\matrixsummaryfield{#1}\matrixpunct{|}"
            r"\matrixsummaryfield{#2}\matrixpunct{|}"
            r"\matrixsummaryfield{#3}\matrixpunct{|}"
            r"\matrixsummaryfield{#4}}"
        ),
        (
            r"\providecommand{\matrixsummaryfive}[5]{"
            r"\matrixsummaryfield{#1}\matrixpunct{|}"
            r"\matrixsummaryfield{#2}\matrixpunct{|}"
            r"\matrixsummaryfield{#3}\matrixpunct{|}"
            r"\matrixsummaryfield{#4}\matrixpunct{|}"
            r"\matrixsummaryfield{#5}}"
        ),
    ]
    entries = data.get("entries", {})
    if not isinstance(entries, dict):
        entries = {}
    validation_summary = _validation_summary(entries, shown)
    model_profile = str(data.get("model_profile", "built-in-sm"))
    model_label = str(data.get("model_label", _matrix_model_label(model_profile)))
    title = _matrix_title(color_accuracy, model_profile=model_profile)
    for chunk_index, chunk in enumerate(chunks):
        multiplicity_columns = (
            r"@{\hspace{0.055in}}".join("L{2.51in}" for _ in chunk)
        )
        colspec = (
            r"@{}r@{\hspace{0.055in}}L{1.42in}@{\hspace{0.075in}}"
            + multiplicity_columns
            + r"@{}"
        )
        arraystretch = "1.13" if color_accuracy == "lc" else "1.23"
        row_space = "0.13em" if color_accuracy == "lc" else "0.22em"
        lines.extend(
            [
                r"\begin{landscape}",
                (
                    rf"\section{{{title}}}"
                    if chunk_index == 0
                    else rf"\subsection*{{{title} (continued)}}"
                ),
                r"\begingroup",
                r"\scriptsize",
                r"\setlength{\tabcolsep}{2.2pt}",
                rf"\renewcommand{{\arraystretch}}{{{arraystretch}}}",
            ]
        )
        if chunk_index == 0:
            lines.extend(
                _matrix_short_intro_latex(
                    color_accuracy,
                    model_label=model_label,
                    model_profile=model_profile,
                )
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
                cells.append(
                    _latex_cell(
                        case if isinstance(case, dict) else {},
                        color_accuracy=color_accuracy,
                    )
                )
            if row_index % 2 == 0:
                lines.append(r"\rowcolor{refblue}")
            lines.append(" & ".join(cells) + r" \\")
            lines.append(rf"\addlinespace[{row_space}]")
        lines.extend(_summary_rows_latex(entries, chunk, color_accuracy=color_accuracy))
        lines.extend(
            [
                r"\bottomrule",
                r"\end{longtable}",
                r"\endgroup",
            ]
        )
        if chunk_index == len(chunks) - 1:
            lines.extend(
                _matrix_long_intro_latex(
                    color_accuracy,
                    validation_summary,
                    model_label=model_label,
                    model_profile=model_profile,
                )
            )
            lines.extend(_matrix_run_settings_latex(color_accuracy))
            lines.extend(
                _matrix_status_notes_latex(
                    entries,
                    shown,
                    color_accuracy=color_accuracy,
                )
            )
        lines.extend([r"\end{landscape}", ""])
    return "\n".join(lines)


def render_shared_lc_recycling_table(data: dict[str, Any]) -> str:
    entries = data.get("entries", {})
    family = entries.get("gg_tt_jets", {}) if isinstance(entries, dict) else {}
    rows: list[str] = []
    if isinstance(family, dict):
        sortable: list[tuple[int, dict[str, Any]]] = []
        for n_text, case in family.items():
            try:
                n_final = int(n_text)
            except (TypeError, ValueError):
                continue
            if n_final >= 3 and isinstance(case, dict):
                sortable.append((n_final, case))
        for n_final, case in sorted(sortable):
            amplicol = case.get("amplicol", {})
            jit = case.get("pyamplicol_jit", {})
            validation = case.get("validation", {})
            if not isinstance(amplicol, dict) or not isinstance(jit, dict):
                continue
            amp_generation = _optional_positive_float(
                amplicol.get("all_flow_generation_s")
            )
            jit_generation = _optional_positive_float(
                jit.get("all_flow_generation_s")
            )
            amp_selected = _optional_positive_float(
                amplicol.get("runtime_us_per_point")
            )
            jit_selected = _optional_positive_float(jit.get("wall_us_per_point"))
            amp_all = _optional_positive_float(
                amplicol.get("all_flow_runtime_us_per_point")
            )
            jit_all = _optional_positive_float(
                jit.get("all_flow_wall_us_per_point")
            )
            jit_fraction = _optional_positive_float(
                jit.get("jit_fraction_of_generation")
            )
            current_count_value = _optional_positive_float(jit.get("current_count"))
            root_count_value = _optional_positive_float(
                jit.get("amplitude_root_count")
            )
            required = (
                amp_generation,
                jit_generation,
                amp_selected,
                jit_selected,
                amp_all,
                jit_all,
                jit_fraction,
                current_count_value,
                root_count_value,
            )
            if any(value is None for value in required):
                continue
            current_count = int(current_count_value)
            root_count = int(root_count_value)
            extra_gluons = n_final - 2
            gluon_term = "g" if extra_gluons == 1 else f"{extra_gluons}g"
            relative_difference = (
                _optional_float(validation.get("max_relative_difference"))
                if isinstance(validation, dict)
                else None
            )
            rows.append(
                " & ".join(
                    (
                        rf"\(gg\to t\bar t+{gluon_term}\)",
                        str(root_count),
                        str(current_count),
                        (
                            f"{_compact_latex_number(amp_generation)} / "
                            f"{_compact_latex_number(jit_generation)}"
                        ),
                        _compact_latex_number(jit_fraction),
                        rf"\({_format_ratio(jit_selected / amp_selected)}\times\)",
                        (
                            f"{_compact_latex_number(amp_all)} / "
                            f"{_compact_latex_number(jit_all)}"
                        ),
                        _scientific_latex(relative_difference),
                    )
                )
                + r" \\"
            )

    lines = [
        "% Generated by docs/result_matrix.py from result_matrix_data.json.",
        r"\Needspace{14\baselineskip}",
        r"\begin{center}",
        r"\begingroup",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{1.75pt}",
        r"\begin{tabular}{@{}L{1.2in}rrrrrrr@{}}",
        r"\toprule",
        (
            r"Process & Orders & Curr. & \AC/\PAC\ O1 gen & O1 JIT frac. & "
            r"Sel. ratio & \AC/\PAC\ all run & Rel. diff. \\"
        ),
        r"\midrule",
    ]
    if rows:
        lines.extend(rows)
    else:
        lines.append(r"\multicolumn{8}{c}{No completed shared-LC rows} \\")
    lines.extend(
        (
            r"\bottomrule",
            r"\end{tabular}",
            r"\endgroup",
            r"\end{center}",
        )
    )
    return "\n".join(lines) + "\n"


def _matrix_title(
    color_accuracy: str,
    *,
    model_profile: str = "built-in-sm",
) -> str:
    prefix = "External-SM " if model_profile != "built-in-sm" else "Generic "
    if color_accuracy == "nlc":
        return f"{prefix}NLC Performance Matrix"
    if color_accuracy == "full":
        return f"{prefix}Full-Colour Performance Matrix"
    return f"{prefix}LC Performance Matrix (SymJIT O1)"


def _matrix_short_intro_latex(
    color_accuracy: str,
    *,
    model_label: str | None = None,
    model_profile: str = "built-in-sm",
) -> list[str]:
    if model_label is None:
        model_label = _matrix_model_label(model_profile)
    reference = (
        r"\AC\ uses raw generated-library colour probes. "
        if color_accuracy in {"nlc", "full"}
        else r"\AC\ uses generated-library timings. "
    )
    slot_text = (
        r"slots are \AC, \PAC\ SymJIT O1 selected flow, "
        r"and \PAC\ SymJIT O1 all flows. "
        if color_accuracy == "lc"
        else r"slots are \AC, \PAC\ JIT O1, and \PAC\ C++ O3. "
    )
    return [
        (
            rf"\noindent\footnotesize \PAC\ model source: \texttt{{{_latex_escape(model_label)}}}. "
            if model_profile != "built-in-sm"
            else ""
        )
        +
        (
            r"Cell format: generation time above runtime "
            r"per phase-space point; "
            + slot_text
            + reference
            + r"Conventions and gaps are summarized after the table."
        ),
        r"\par\smallskip",
    ]


def _matrix_long_intro_latex(
    color_accuracy: str,
    validation_summary: dict[str, Any],
    *,
    model_label: str | None = None,
    model_profile: str = "built-in-sm",
) -> list[str]:
    if model_label is None:
        model_label = _matrix_model_label(model_profile)
    if color_accuracy == "lc":
        description = (
            r"\noindent\scriptsize Each \AC\ pair is selected-flow/helicity-summed "
            r"then all-flow/fixed-helicity; \PAC\ SymJIT O1 uses matching pruned "
            r"and shared artifacts.  Entries are generation seconds above runtime wall "
            r"microseconds per point; \PAC-only cells are absolute. "
        )
    else:
        description = (
            r"\noindent\scriptsize Each cell is raw-library \AC, \PAC\ JIT O1, "
            r"and \PAC\ C++ O3.  Entries are generation seconds above runtime "
            r"microseconds per point; \PAC\ runtime shows coloured wall and "
            r"black pure-evaluator ratios, or absolute time without an \AC\ reference. "
        )
    return [
        description
        + (
            rf"The \PAC\ entries use the compiled \texttt{{{_latex_escape(model_label)}}} "
            r"source; the \AC\ columns are the same matched external reference used by "
            r"the built-in-SM table. "
            if model_profile != "built-in-sm"
            else ""
        )
        + r"Green/orange/red means below one/below two/at least two.  Summaries "
        r"are \texttt{min|max|avg|med|sum}; \texttt{sum} compares paired totals. "
        + (
            r"\texttt{>30 GB RAM} marks only the affected selected/all-flow "
            r"workload. "
            if color_accuracy == "lc"
            else ""
        )
        + r"A \textbf{VALIDATION FAILED} marker means the same-point difference "
        r"exceeds the tolerance recorded for that validation result.",
        r"\par\vspace{0.15em}",
    ]


def _matrix_run_settings_latex(color_accuracy: str) -> list[str]:
    if color_accuracy == "lc":
        reference_text = (
            r"\AC\ uses its generated library for the selected workload and "
            r"the fixed-helicity \texttt{imode=2} colour probe for all flows. "
            r"\PAC\ writes separate \texttt{selected\_flow} and shared "
            r"\texttt{all\_flows} artifacts."
        )
    else:
        reference_text = (
            r"\AC\ uses the raw generated-library path "
            r"\texttt{amplicol\_generate --library=create-raw} followed by "
            r"\texttt{amplicol\_color\_library\_probe}, preserving NLC/full "
            r"colour-order interferences."
        )
    lines = [
        r"\paragraph{Settings.}",
        (
            rf"\noindent\scriptsize \path{{docs/result_matrix.py}} "
            rf"(\texttt{{--color-accuracy {color_accuracy}}}) writes this "
            r"table and its JSON cache. "
            + reference_text
            + (
                r" \PAC\ uses Rusticol double precision; selected/all-flow "
                r"batches are 128/64.  Selected JIT stages measure chunks "
                r"around base 128; all-flow outputs use uniform chunk 8192. "
                if color_accuracy == "lc"
                else r" \PAC\ uses Rusticol double precision, batch size 64, "
                r"output chunk size 128, "
            )
            + r" Stage-local inputs, ten Horner iterations, default CPE rounds, "
            r"expanded optimization limits, "
            + (
                r"SymJIT O1, "
                if color_accuracy == "lc"
                else r"SymJIT O1, "
            )
            + r"and "
            r"\texttt{--target-runtime 10} are used. "
            + (
                r"The two LC workloads are displayed side by side."
                if color_accuracy == "lc"
                else r"Only C++ O3 generation is subject to the 15-minute compile cap."
            )
        ),
        r"\par\vspace{0.15em}",
    ]
    if color_accuracy == "lc":
        lines.extend(
            [
                r"\paragraph{Runtime interpretation.}",
                (
                    r"\noindent\scriptsize On AArch64, Rusticol runs saved SymJIT "
                    r"stages as native two-lane complex evaluators and packs each "
                    r"stage input once across its chunks.  Wall time measures the "
                    r"public \texttt{evaluate} call; the profiler separates evaluator "
                    r"and packing work.  Low-multiplicity fixed overhead can dominate "
                    r"an \AC\ call measured in nanoseconds.  Ratios at or above two "
                    r"are reviewed throughout, with \(n\geq5\) as the scaling gate."
                ),
                r"\par\vspace{0.15em}",
            ]
        )
    return lines


def _matrix_status_notes_latex(
    entries: dict[str, Any],
    n_values: Sequence[int],
    *,
    color_accuracy: str,
) -> list[str]:
    unsupported_four_quark_line = False
    fortran_color_reference_gaps: list[str] = []
    documented_backend_limits: list[str] = []
    ram_limited: list[str] = []
    non_structural: list[str] = []
    cxx_not_run: list[str] = []
    cxx_over_budget: list[str] = []
    include_cxx_notes = color_accuracy != "lc"
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
            reference_gap = _fortran_color_reference_unavailable_reason(
                ref,
                _mode(case, "pyamplicol_jit"),
            )
            if reference_gap is not None:
                fortran_color_reference_gaps.append(
                    rf"{base.label}, \(n={n_final}\)"
                )
            for key, label in (
                ("amplicol", r"\AC"),
                ("pyamplicol_jit", r"\PAC\ JIT"),
                ("pyamplicol_cpp_o3", r"\PAC\ C++ O3"),
            ):
                mode = _mode(case, key)
                status = str(mode.get("status", ""))
                error = str(mode.get("error", ""))
                if (
                    include_cxx_notes
                    and
                    key == "pyamplicol_cpp_o3"
                    and not mode
                    and _ok(ref)
                    and _ok(_mode(case, "pyamplicol_jit"))
                ):
                    cxx_not_run.append(rf"{base.label}, \(n={n_final}\)")
                    continue
                if include_cxx_notes and key == "pyamplicol_cpp_o3" and status == "timeout":
                    cxx_over_budget.append(rf"{base.label}, \(n={n_final}\)")
                    continue
                if status == "ram_limit":
                    ram_limited.append(rf"{base.label}, \(n={n_final}\), {label}")
                    continue
                if (
                    key == "pyamplicol_jit"
                    and str(mode.get("all_flow_status", "")) == "ram_limit"
                ):
                    ram_limited.append(
                        rf"{base.label}, \(n={n_final}\), {label} all-flow"
                    )
                if _is_documented_backend_limitation(mode):
                    documented_backend_limits.append(
                        rf"{base.label}, \(n={n_final}\), {label}"
                    )
                    continue
                if (
                    _amplicol_error_is_unsupported(error)
                    or _color_accuracy_error_is_structural_unsupported(error)
                ):
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
    if unsupported_four_quark_line and color_accuracy != "lc":
        lines.extend(
            [
                r"\paragraph{Structural N/A.}",
                (
                    r"\noindent\scriptsize The four-quark-line row "
                    r"\(d\bar d\to u\bar u\,s\bar s\,c\bar c+(n-6)g\) "
                    r"exceeds the current \AC\ reference limit; completed "
                    r"\PAC\ cells are absolute timings."
                ),
                r"\par\vspace{0.15em}",
            ]
        )
    if fortran_color_reference_gaps:
        lines.extend(
            [
                r"\paragraph{Fortran reference gaps.}",
                (
                    rf"\noindent\scriptsize {len(fortran_color_reference_gaps)} "
                    r"completed \PAC\ cells have no direct Fortran colour reference "
                    r"because that path stops above two quark lines; they are shown "
                    r"as absolute timings."
                ),
                r"\par\vspace{0.15em}",
            ]
        )
    if ram_limited and color_accuracy != "lc":
        lines.extend(
            [
                r"\paragraph{RAM gaps.}",
                (
                    rf"\noindent\scriptsize {len(ram_limited)} entries exceeded "
                    r"the 30 GB watchdog and are shown as \texttt{>30 GB RAM}."
                ),
                r"\par\vspace{0.15em}",
            ]
        )
    if non_structural:
        shown = "; ".join(non_structural[:4])
        if len(non_structural) > 4:
            shown += rf"; plus {len(non_structural) - 4} more"
        lines.extend(
            [
                r"\paragraph{Incomplete entries.}",
                (
                    r"\noindent\scriptsize Non-structural incomplete cells: "
                    + shown
                    + "."
                ),
                r"\par\vspace{0.15em}",
            ]
        )
    cxx_compile_budget = [
        *(rf"{entry}: not run" for entry in cxx_not_run),
        *(rf"{entry}: \texttt{{t/o >15 min}}" for entry in cxx_over_budget),
    ]
    if cxx_compile_budget:
        lines.extend(
            [
                r"\paragraph{C++ O3 gaps.}",
                (
                    rf"\noindent\scriptsize {len(cxx_not_run)} optional cells were "
                    rf"not run and {len(cxx_over_budget)} exceeded the 15-minute "
                    r"compile budget; completed \AC/JIT data remain shown."
                ),
                r"\par\vspace{0.15em}",
            ]
        )
    if documented_backend_limits:
        lines.extend(
            [
                r"\paragraph{JIT backend limits.}",
                (
                    rf"\noindent\scriptsize {len(documented_backend_limits)} cells "
                    r"hit the documented AArch64 complex-JIT materialization "
                    r"limit; raw assertions remain in the JSON cache."
                ),
                r"\par\vspace{0.15em}",
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
            if validation.get("status") in {
                "not_available",
                "skipped",
                "unsupported",
            }:
                unvalidated += 1
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
        r"next to any validated entry whose relative difference exceeds its "
        r"recorded tolerance."
    )


def _chunks(values: Sequence[int], size: int) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(values[index : index + size]) for index in range(0, len(values), size))


def _latex_cell(case: dict[str, Any], *, color_accuracy: str) -> str:
    if not case:
        return _missing_na()
    if case.get("status") == "not_applicable":
        return _structural_na()
    ref = _mode(case, "amplicol")
    jit = _mode(case, "pyamplicol_jit")
    cpp = _mode(case, "pyamplicol_cpp_o3")
    if _case_is_structural_unsupported(ref, jit, cpp):
        return _structural_na()
    if color_accuracy == "lc":
        cell = _latex_lc_cell(case, ref, jit, cpp)
        marker = _validation_cell_marker(case)
        if marker:
            return rf"\begin{{tabular}}[t]{{@{{}}l@{{}}}}{cell}\\[-0.1em]{marker}\end{{tabular}}"
        return cell
    if _ok(ref):
        gen = _metric_parts(
            ref,
            jit,
            cpp,
            metric="generation_s",
            formatter=_format_seconds,
        )
        run = _runtime_parts(ref, jit, cpp, include_eval_ratio=True)
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
        rf"\matrixcellnonlc{{{gen[0]}}}{{{gen[1]}}}{{{gen[2]}}}"
        rf"{{{run[0]}}}{{{run[1]}}}{{{run[2]}}}"
    )
    marker = _validation_cell_marker(case)
    if marker:
        return rf"\begin{{tabular}}[t]{{@{{}}l@{{}}}}{cell}\\[-0.1em]{marker}\end{{tabular}}"
    return cell


def _latex_lc_cell(
    case: dict[str, Any],
    ref: dict[str, Any],
    jit: dict[str, Any],
    cpp: dict[str, Any],
) -> str:
    all_flow_failure = _lc_all_flow_failure_latex(jit)
    if _ok(ref):
        selected_ref_gen = _optional_float(ref.get("generation_s"))
        all_ref_gen = _lc_amplicol_all_flow_generation_s(ref)
        gen_ref = _reference_pair_from_values(
            selected_ref_gen,
            all_ref_gen,
            formatter=_format_seconds,
        )
        gen_jit = _generation_ratio_or_absolute_from_value(
            jit,
            _lc_pyamplicol_selected_generation_s(jit),
            selected_ref_gen,
            mark_selected_jit_level=True,
        )
        gen_all = (
            all_flow_failure
            or _generation_ratio_or_absolute_from_value(
                jit,
                _lc_pyamplicol_all_flow_generation_s(jit),
                all_ref_gen,
            )
        )
        run_ref = _reference_metric_pair_or_status(
            ref,
            metric="runtime_us_per_point",
            all_flow_metric="all_flow_runtime_us_per_point",
            formatter=_format_us,
        )
        selected_ref_runtime = _optional_float(ref.get("runtime_us_per_point"))
        all_ref_runtime = (
            _optional_float(ref.get("all_flow_runtime_us_per_point"))
            or selected_ref_runtime
        )
        run_jit = (
            _missing_runtime_pair()
            if selected_ref_runtime is None
            else _runtime_ratio_pair(jit, selected_ref_runtime)
        )
        run_all = (
            all_flow_failure
            or (
                _missing_runtime_pair()
                if all_ref_runtime is None
                else _runtime_ratio_pair(
                    jit,
                    all_ref_runtime,
                    wall_key="all_flow_wall_us_per_point",
                    core_key="all_flow_runtime_us_per_point",
                )
            )
        )
    else:
        gen_ref = _reference_metric_or_status(
            ref,
            metric="generation_s",
            formatter=_format_seconds,
        )
        gen_jit = _absolute_generation(
            {**jit, "generation_s": _lc_pyamplicol_selected_generation_s(jit)},
            mark_selected_jit_level=True,
        )
        gen_all = (
            all_flow_failure
            or _absolute_generation(
                {**jit, "generation_s": _lc_pyamplicol_all_flow_generation_s(jit)}
            )
        )
        run_ref = _reference_metric_or_status(
            ref,
            metric="runtime_us_per_point",
            formatter=_format_us,
        )
        run_jit = _absolute_runtime_pair(jit)
        run_all = (
            all_flow_failure
            or _absolute_runtime_pair(
                {
                    **jit,
                    "wall_us_per_point": jit.get("all_flow_wall_us_per_point"),
                    "runtime_us_per_point": jit.get("all_flow_runtime_us_per_point"),
                }
            )
        )
    return (
        rf"\matrixcell{{{gen_ref}}}{{{gen_jit}}}{{{gen_all}}}"
        rf"{{{run_ref}}}{{{run_jit}}}{{{run_all}}}"
    )


def _lc_all_flow_failure_latex(mode: dict[str, Any]) -> str:
    raw_status = mode.get("all_flow_status")
    status = "" if raw_status is None else str(raw_status).lower()
    if status in {"", "ok"}:
        return ""
    return _latex_failure(
        {
            "status": status,
            "error": mode.get("all_flow_error", ""),
        }
    )


def _case_is_structural_unsupported(
    ref: dict[str, Any],
    jit: dict[str, Any],
    cpp: dict[str, Any],
) -> bool:
    if any(mode.get("status") == "ram_limit" for mode in (ref, jit, cpp)):
        return False
    if _ok(ref) or _ok(jit) or _ok(cpp):
        return False
    return any(_mode_is_structural_unsupported(mode) for mode in (ref, jit, cpp))


def _mode_is_structural_unsupported(mode: dict[str, Any]) -> bool:
    if not mode:
        return False
    error = str(mode.get("error", ""))
    status = str(mode.get("status", "")).lower()
    return (
        status == "unsupported"
        or _amplicol_error_is_unsupported(error)
        or _color_accuracy_error_is_structural_unsupported(error)
    )


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
        rf"\matrixcellnonlc{{{gen[0]}}}{{{gen[1]}}}{{{gen[2]}}}"
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


def _reference_metric_pair_or_status(
    mode: dict[str, Any],
    *,
    metric: str,
    all_flow_metric: str,
    formatter,
) -> str:
    if not _ok(mode):
        return _latex_failure(mode) if mode else _missing_na()
    selected = _optional_float(mode.get(metric))
    all_flow = _optional_float(mode.get(all_flow_metric))
    if selected is None and all_flow is None:
        return _missing_na()
    if selected is None:
        return formatter(all_flow)  # type: ignore[arg-type]
    if all_flow is None:
        return formatter(selected)
    return rf"\matrixrefpair{{{formatter(selected)}}}{{{formatter(all_flow)}}}"


def _reference_pair_from_values(
    selected: float | None,
    all_flow: float | None,
    *,
    formatter,
) -> str:
    if selected is None and all_flow is None:
        return _missing_na()
    if selected is None:
        return formatter(all_flow)  # type: ignore[arg-type]
    if all_flow is None:
        return formatter(selected)
    return rf"\matrixrefpair{{{formatter(selected)}}}{{{formatter(all_flow)}}}"


def _generation_ratio_or_absolute_from_value(
    mode: dict[str, Any],
    value: float | None,
    ref_value: float | None,
    *,
    mark_selected_jit_level: bool = False,
) -> str:
    if _is_documented_backend_limitation(mode):
        return _structural_na()
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _latex_failure(mode)
    if not _ok(mode) or value is None:
        return _missing_na()
    if ref_value is None or ref_value <= 0.0:
        result = _format_seconds(value)
    else:
        result = _ratio_latex(value / ref_value)
    if mark_selected_jit_level and _selected_jit_is_o3(mode):
        result += r"\textsuperscript{*}"
    return result


def _lc_amplicol_all_flow_generation_s(mode: dict[str, Any]) -> float | None:
    source = str(mode.get("all_flow_generation_source", ""))
    value = _optional_float(mode.get("all_flow_generation_s"))
    if source.startswith("amplicol_color_probe_imode2"):
        return value
    probe = str(mode.get("all_flow_reference_probe", ""))
    if "amplicol_color_probe" in probe:
        derived = _amplicol_color_probe_setup_s_from_records(
            mode.get("commands", ()) if isinstance(mode.get("commands"), list) else (),
            mode.get("all_flow_timing_rows", ())
            if isinstance(mode.get("all_flow_timing_rows"), list)
            else (),
        )
        if derived is not None:
            return derived
        return None
    return value


def _lc_pyamplicol_selected_generation_s(mode: dict[str, Any]) -> float | None:
    explicit = _optional_float(mode.get("selected_generation_s"))
    if explicit is not None:
        return explicit
    payload = mode.get("selected_generate_payload")
    if isinstance(payload, dict):
        elapsed = _optional_float(payload.get("_command_elapsed_s"))
        if elapsed is not None:
            return elapsed
        internal = _optional_float(payload.get("generation_s"))
        if internal is not None:
            return internal
    return _optional_float(mode.get("generation_s"))


def _lc_pyamplicol_all_flow_generation_s(mode: dict[str, Any]) -> float | None:
    explicit = _optional_float(mode.get("all_flow_generation_s"))
    if explicit is not None:
        return explicit
    payload = mode.get("all_flow_generate_payload")
    if isinstance(payload, dict):
        elapsed = _optional_float(payload.get("_command_elapsed_s"))
        if elapsed is not None:
            return elapsed
        internal = _optional_float(payload.get("generation_s"))
        if internal is not None:
            return internal
    return _optional_float(mode.get("generation_s"))


def _absolute_generation(
    mode: dict[str, Any],
    *,
    mark_selected_jit_level: bool = False,
) -> str:
    if _is_documented_backend_limitation(mode):
        return _structural_na()
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _latex_failure(mode)
    if not _ok(mode):
        return _missing_na()
    value = _optional_float(mode.get("generation_s"))
    if value is None:
        return _missing_na()
    result = _format_seconds(value)
    if mark_selected_jit_level and _selected_jit_is_o3(mode):
        result += r"\textsuperscript{*}"
    return result


def _absolute_runtime_pair(mode: dict[str, Any]) -> str:
    if _is_documented_backend_limitation(mode):
        return _structural_na()
    if str(mode.get("status", "")).lower() == "ram_limit":
        return _latex_failure(mode)
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _missing_na()
    if not _ok(mode):
        return _missing_na()
    wall = _optional_float(mode.get("wall_us_per_point"))
    if wall is None:
        return _missing_na()
    return rf"\matrixpunct{{(}}\texttt{{{_format_us_value(wall)}}}\matrixpunct{{)}}"


def _validation_cell_marker(case: dict[str, Any]) -> str:
    validation = case.get("validation")
    if not isinstance(validation, dict):
        return ""
    status = str(validation.get("status", "")).lower()
    if status in {
        "",
        "missing",
        "not_applicable",
        "not_available",
        "unsupported",
        "skipped",
    }:
        return ""
    tolerance = _optional_float(validation.get("tolerance")) or 1.0e-8
    rel = _optional_float(validation.get("max_relative_difference"))
    if status == "ok" and (rel is None or rel <= tolerance):
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


def _metric_ratio_or_status(
    mode: dict[str, Any],
    metric: str,
    ref_value: float | None,
) -> str:
    if ref_value is None or ref_value <= 0.0:
        return _missing_ratio()
    if _is_documented_backend_limitation(mode):
        return _missing_ratio()
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _latex_failure(mode)
    if not _ok(mode):
        return _missing_ratio()
    value = _optional_float(mode.get(metric))
    if value is None:
        return _missing_ratio()
    return _ratio_latex(value / ref_value)


def _runtime_parts(
    ref: dict[str, Any],
    jit: dict[str, Any],
    cpp: dict[str, Any],
    *,
    include_eval_ratio: bool = False,
) -> tuple[str, str, str]:
    ref_value = _optional_float(ref.get("runtime_us_per_point"))
    if ref_value is None or ref_value <= 0.0:
        return (_missing_na(), _missing_runtime_pair(), _missing_runtime_pair())
    formatter = _runtime_ratio_wall_eval_pair if include_eval_ratio else _runtime_ratio_pair
    return (
        _format_us(ref_value),
        formatter(jit, ref_value),
        formatter(cpp, ref_value),
    )


def _summary_rows_latex(
    entries: dict[str, Any],
    n_values: Sequence[int],
    *,
    color_accuracy: str,
) -> list[str]:
    lc_summary = color_accuracy == "lc"
    generation_one_flow_cells = [
        (
            r"\multicolumn{2}{@{}L{1.74in}@{\hspace{0.075in}}}"
            + (
                r"{\textbf{gen O1 one-flow, hel. sum}}"
                if lc_summary
                else r"{\textbf{summary: gen}}"
            )
        )
    ]
    generation_all_flow_cells = [
        (
            r"\multicolumn{2}{@{}L{1.74in}@{\hspace{0.075in}}}"
            r"{\textbf{gen O1 all-flows, one-hel}}"
        )
    ]
    runtime_one_flow_cells = [
        (
            r"\multicolumn{2}{@{}L{1.74in}@{\hspace{0.075in}}}"
            + (
                r"{\textbf{run O1 one-flow, hel. sum}}"
                if lc_summary
                else r"{\textbf{summary: run}}"
            )
        )
    ]
    runtime_all_flow_cells = [
        (
            r"\multicolumn{2}{@{}L{1.74in}@{\hspace{0.075in}}}"
            r"{\textbf{run O1 all-flows, one-hel}}"
        )
    ]
    for n_final in n_values:
        summary = _column_summary(entries, n_final, color_accuracy=color_accuracy)
        generation_one_flow_cells.append(
            _summary_cell(
                _summary_numeric_stats_line(summary["amplicol_generation_one_flow"]),
                _summary_ratio_stats_line(
                    summary["jit_generation_one_flow_ratio"],
                    _summed_ratio(
                        summary["jit_generation_one_flow_paired"],
                        summary["jit_generation_one_flow_ref_paired"],
                    ),
                ),
            )
        )
        generation_all_flow_cells.append(
            _summary_cell(
                _summary_numeric_stats_line(summary["amplicol_generation_all_flow"]),
                _summary_ratio_stats_line(
                    summary["jit_generation_all_flow_ratio"],
                    _summed_ratio(
                        summary["jit_generation_all_flow_paired"],
                        summary["jit_generation_all_flow_ref_paired"],
                    ),
                ),
            )
        )
        runtime_one_flow_cells.append(
            _summary_cell(
                _summary_numeric_stats_line(summary["amplicol_runtime_one_flow"]),
                _summary_ratio_stats_line(
                    summary["jit_runtime_one_flow_ratio"],
                    _summed_ratio(
                        summary["jit_runtime_one_flow_paired"],
                        summary["jit_runtime_one_flow_ref_paired"],
                    ),
                ),
            )
        )
        runtime_all_flow_cells.append(
            _summary_cell(
                _summary_numeric_stats_line(summary["amplicol_runtime_all_flow"]),
                _summary_ratio_stats_line(
                    summary["jit_runtime_all_flow_ratio"],
                    _summed_ratio(
                        summary["jit_runtime_all_flow_paired"],
                        summary["jit_runtime_all_flow_ref_paired"],
                    ),
                ),
            )
        )
    rows = [
        r"\specialrule{1.05pt}{0.25em}{0.20em}",
        " & ".join(generation_one_flow_cells) + r" \\",
        r"\addlinespace[0.16em]",
    ]
    if lc_summary:
        rows.extend(
            [
                " & ".join(generation_all_flow_cells) + r" \\",
                r"\addlinespace[0.12em]",
                " & ".join(runtime_one_flow_cells) + r" \\",
                r"\addlinespace[0.12em]",
                " & ".join(runtime_all_flow_cells) + r" \\",
                r"\addlinespace[0.12em]",
            ]
        )
    else:
        rows.extend(
            [
                " & ".join(runtime_one_flow_cells) + r" \\",
                r"\addlinespace[0.12em]",
            ]
        )
    return rows


def _column_summary(
    entries: dict[str, Any],
    n_final: int,
    *,
    color_accuracy: str,
) -> dict[str, list[float]]:
    summary: dict[str, list[float]] = {
        "amplicol_generation_one_flow": [],
        "jit_generation_one_flow_ratio": [],
        "jit_generation_one_flow_paired": [],
        "jit_generation_one_flow_ref_paired": [],
        "amplicol_generation_all_flow": [],
        "jit_generation_all_flow_ratio": [],
        "jit_generation_all_flow_paired": [],
        "jit_generation_all_flow_ref_paired": [],
        "amplicol_runtime_one_flow": [],
        "jit_runtime_one_flow_ratio": [],
        "jit_runtime_one_flow_paired": [],
        "jit_runtime_one_flow_ref_paired": [],
        "amplicol_runtime_all_flow": [],
        "jit_runtime_all_flow_ratio": [],
        "jit_runtime_all_flow_paired": [],
        "jit_runtime_all_flow_ref_paired": [],
    }
    for base in BASE_PROCESSES:
        if n_final < base.min_final_count:
            continue
        row_entries = entries.get(base.key, {})
        if not isinstance(row_entries, dict):
            continue
        case = row_entries.get(str(n_final), {})
        if not isinstance(case, dict):
            continue
        ref = _mode(case, "amplicol")
        jit = _mode(case, "pyamplicol_jit")
        if not _ok(ref):
            continue
        ref_generation = _optional_positive_float(ref.get("generation_s"))
        if ref_generation is not None:
            summary["amplicol_generation_one_flow"].append(ref_generation)
            jit_generation = _optional_positive_float(
                _lc_pyamplicol_selected_generation_s(jit)
                if color_accuracy == "lc"
                else jit.get("generation_s")
            )
            if _ok(jit) and jit_generation is not None:
                summary["jit_generation_one_flow_ratio"].append(
                    jit_generation / ref_generation
                )
                summary["jit_generation_one_flow_paired"].append(jit_generation)
                summary["jit_generation_one_flow_ref_paired"].append(ref_generation)
        if color_accuracy == "lc":
            ref_all_flow_generation = _optional_positive_float(
                _lc_amplicol_all_flow_generation_s(ref)
            )
            if ref_all_flow_generation is not None:
                summary["amplicol_generation_all_flow"].append(ref_all_flow_generation)
                jit_all_flow_generation = _optional_positive_float(
                    _lc_pyamplicol_all_flow_generation_s(jit)
                )
                if _ok(jit) and jit_all_flow_generation is not None:
                    summary["jit_generation_all_flow_ratio"].append(
                        jit_all_flow_generation / ref_all_flow_generation
                    )
                    summary["jit_generation_all_flow_paired"].append(
                        jit_all_flow_generation
                    )
                    summary["jit_generation_all_flow_ref_paired"].append(
                        ref_all_flow_generation
                    )
        ref_runtime = _optional_positive_float(ref.get("runtime_us_per_point"))
        if ref_runtime is not None:
            summary["amplicol_runtime_one_flow"].append(ref_runtime)
            jit_runtime = _optional_positive_float(jit.get("wall_us_per_point"))
            if _ok(jit) and jit_runtime is not None:
                summary["jit_runtime_one_flow_ratio"].append(jit_runtime / ref_runtime)
                summary["jit_runtime_one_flow_paired"].append(jit_runtime)
                summary["jit_runtime_one_flow_ref_paired"].append(ref_runtime)
        if color_accuracy == "lc":
            ref_all_flow_runtime = _optional_positive_float(
                ref.get("all_flow_runtime_us_per_point")
            )
            if ref_all_flow_runtime is not None:
                summary["amplicol_runtime_all_flow"].append(ref_all_flow_runtime)
                jit_all_flow_runtime = _optional_positive_float(
                    jit.get("all_flow_wall_us_per_point")
                )
                if _ok(jit) and jit_all_flow_runtime is not None:
                    summary["jit_runtime_all_flow_ratio"].append(
                        jit_all_flow_runtime / ref_all_flow_runtime
                    )
                    summary["jit_runtime_all_flow_paired"].append(jit_all_flow_runtime)
                    summary["jit_runtime_all_flow_ref_paired"].append(
                        ref_all_flow_runtime
                    )
    return summary


def _summary_cell(first_line: str, second_line: str) -> str:
    return rf"\matrixsummarycell{{{first_line}}}{{{second_line}}}"


def _summary_numeric_stats_line(values: Sequence[float]) -> str:
    stats = _summary_stats(values)
    if stats is None:
        return _missing_na()
    fields = (rf"\texttt{{{_format_compact_number(value)}}}" for value in stats)
    return r"\matrixsummaryfour{" + "}{".join(fields) + "}"


def _summary_ratio_stats_line(
    values: Sequence[float],
    summed_ratio: float | None,
) -> str:
    stats = _summary_stats(values)
    if stats is None:
        return _missing_na()
    fields = [_summary_ratio_fragment(value) for value in stats]
    fields.append(_summary_ratio_fragment(summed_ratio) if summed_ratio is not None else _missing_na())
    return r"\matrixsummaryfive{" + "}{".join(fields) + "}"


def _summary_ratio_fragment(value: float) -> str:
    color = "speedgreen" if value < 1.0 else "speedorange" if value < 2.0 else "speedred"
    return rf"\textcolor{{{color}}}{{\texttt{{x{_format_ratio(value)}}}}}"


def _summed_ratio(numerators: Sequence[float], denominators: Sequence[float]) -> float | None:
    pairs = [
        (numerator, denominator)
        for numerator, denominator in zip(numerators, denominators, strict=False)
        if math.isfinite(numerator) and math.isfinite(denominator) and denominator > 0.0
    ]
    if not pairs:
        return None
    denominator_sum = math.fsum(denominator for _, denominator in pairs)
    if denominator_sum <= 0.0:
        return None
    return math.fsum(numerator for numerator, _ in pairs) / denominator_sum


def _summary_stats(values: Sequence[float]) -> tuple[float, float, float, float] | None:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return None
    count = len(finite)
    midpoint = count // 2
    median = (
        finite[midpoint]
        if count % 2
        else 0.5 * (finite[midpoint - 1] + finite[midpoint])
    )
    return (finite[0], finite[-1], sum(finite) / count, median)


def _optional_positive_float(value: object) -> float | None:
    number = _optional_float(value)
    if number is None or number <= 0.0:
        return None
    return number


def _format_us_value(value: float) -> str:
    return _format_compact_number(value)


def _format_compact_number(value: float) -> str:
    text = f"{value:.3g}"
    text = text.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")
    return text


def _runtime_ratio_pair(
    mode: dict[str, Any],
    ref_value: float,
    *,
    wall_key: str = "wall_us_per_point",
    core_key: str = "runtime_us_per_point",
) -> str:
    if _is_documented_backend_limitation(mode):
        return _missing_runtime_pair()
    if str(mode.get("status", "")).lower() == "ram_limit":
        return _latex_failure(mode)
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _missing_runtime_pair()
    if not _ok(mode):
        return _missing_runtime_pair()
    wall = _optional_float(mode.get(wall_key))
    if wall is None:
        return _missing_runtime_pair()
    ratio = _ratio_latex(wall / ref_value)
    if wall_key == "wall_us_per_point" and _selected_jit_is_o3(mode):
        return ratio + r"\textsuperscript{*}"
    return ratio


def _selected_jit_is_o3(mode: dict[str, Any]) -> bool:
    return _optional_float(mode.get("selected_jit_optimization_level")) == 3.0


def _runtime_ratio_wall_eval_pair(
    mode: dict[str, Any],
    ref_value: float,
    *,
    wall_key: str = "wall_us_per_point",
    core_key: str = "runtime_us_per_point",
) -> str:
    if _is_documented_backend_limitation(mode):
        return _missing_runtime_pair()
    if str(mode.get("status", "")).lower() == "ram_limit":
        return _latex_failure(mode)
    if mode and str(mode.get("status", "")).lower() not in {"", "ok"}:
        return _missing_runtime_pair()
    if not _ok(mode):
        return _missing_runtime_pair()
    wall = _optional_float(mode.get(wall_key))
    core = _optional_float(mode.get(core_key))
    if wall is None and core is None:
        return _missing_runtime_pair()
    wall_text = _missing_ratio() if wall is None else _ratio_latex(wall / ref_value)
    core_text = (
        _missing_ratio()
        if core is None
        else rf"\matrixevalratio{{{_format_ratio(core / ref_value)}}}"
    )
    return (
        r"\begin{tabular}[t]{@{}l@{\hspace{0.006in}}l@{}}"
        + wall_text
        + "&"
        + core_text
        + r"\end{tabular}"
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
    if str(mode.get("status", "")).lower() == "fortran_reference_unavailable":
        return _structural_na()
    if (
        mode.get("status") == "unsupported"
        or _amplicol_error_is_unsupported(error)
        or _color_accuracy_error_is_structural_unsupported(error)
    ):
        return _structural_na()
    status = str(mode.get("status", "missing"))
    if status == "missing":
        return _missing_na()
    if status == "ram_limit":
        return r"\textcolor{speedred}{\texttt{>30 GB RAM}}"
    if status == "timeout" and "C++ O3" in str(mode.get("mode", "")):
        return r"\textcolor{speedred}{\texttt{t/o >15 min}}"
    return rf"\textcolor{{speedred}}{{\texttt{{{_latex_escape(status)}}}}}"


def _amplicol_error_is_unsupported(error: str) -> bool:
    unsupported_markers = (
        "Unknown number of quarks and anti-quarks",
    )
    return any(marker in error for marker in unsupported_markers)


def _color_accuracy_error_is_structural_unsupported(error: str) -> bool:
    unsupported_markers = (
        "more than two quarks",
        "only for zero, one, or two quark pairs",
    )
    return any(marker in error for marker in unsupported_markers)


def _fortran_color_reference_unavailable_reason(
    amplicol: dict[str, Any],
    pyamplicol_jit: dict[str, Any],
) -> str | None:
    """Classify colour-reference gaps caused by Fortran AmpliCol limits.

    The generic pyAmpliCol colour contraction can handle multi-open-line
    sectors, but the current Fortran ``init_col`` reference matrix stops at
    more than two quark lines.  When the pyAmpliCol row itself ran, represent
    this as an unavailable reference rather than as a pyAmpliCol validation
    failure or unsupported process.
    """

    if not _ok(pyamplicol_jit):
        return None
    return _fortran_color_reference_error_reason(str(amplicol.get("error", "")))


def _fortran_color_reference_error_reason(error: str) -> str | None:
    if "more than two quarks" not in error:
        return None
    return (
        "Fortran AmpliCol init_col does not expose a direct NLC/full colour "
        "matrix for more than two quark lines; pyAmpliCol generated and "
        "evaluated the generic colour contraction, but this cell has no "
        "direct Fortran colour-contraction reference."
    )


def _structural_color_accuracy_unsupported_reason(
    process: str,
    *,
    base: BaseProcess,
    color_accuracy: str,
) -> str | None:
    if color_accuracy == "lc":
        return None
    return None


def _unsupported_mode_payload(
    reason: str,
    *,
    color_accuracy: str,
) -> dict[str, Any]:
    return {
        "status": "unsupported",
        "error": reason,
        "reason": reason,
        "color_accuracy": color_accuracy,
        "matrix_settings": {
            "color_accuracy": color_accuracy,
            "structural_unsupported": True,
        },
        "finished_at": _now(),
    }


def _structural_na() -> str:
    return r"\matrixna{black!45}"


def _missing_na() -> str:
    return r"\matrixna{speedred}"


def _missing_ratio() -> str:
    return r"\matrixnaratio{speedred}"


def _missing_runtime_pair() -> str:
    return _missing_ratio()


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
    if _is_documented_backend_limitation(mode) or mode.get("status") == "ram_limit":
        return not _settings_match(mode.get("matrix_settings"), expected_settings)
    if not _ok(mode):
        return True
    return not _settings_match(mode.get("matrix_settings"), expected_settings)


def _validation_record_ok(validation: object) -> bool:
    if not isinstance(validation, dict):
        return False
    tolerance = _optional_float(validation.get("tolerance")) or VALIDATION_REL_TOL
    rel = _optional_float(validation.get("max_relative_difference"))
    return validation.get("status") == "ok" and (rel is None or rel <= tolerance)


def _settings_match(actual: object, expected: dict[str, Any]) -> bool:
    if not isinstance(actual, dict):
        return False
    return actual == expected


def _amplicol_matrix_settings(
    *,
    amplicol_points: int,
    jobs: int,
    color_accuracy: str,
) -> dict[str, Any]:
    return {
        "amplicol_points": int(amplicol_points),
        "jobs": int(jobs),
        "color_accuracy": color_accuracy,
        "process_list_backend": "python",
        "workflow": (
            "full_library_create_make_library_first_entry_runtime_and_fixed_helicity_all_flow"
            if color_accuracy == "lc"
            else "library_create_raw_make_library_color_probe"
        ),
        "runtime_probe": (
            "direct_generated_library_benchmark_first_entry_plus_amplicol_color_probe_fixed_helicity_all_flow"
            if color_accuracy == "lc"
            else "amplicol_color_library_probe"
        ),
        "all_flow_runtime_probe": (
            "amplicol_color_probe_fixed_helicity" if color_accuracy == "lc" else None
        ),
    }


def _pyamplicol_matrix_settings(
    *,
    backend_key: str,
    target_runtime: float,
    n_cores: int,
    selected_lc_sector_ids: set[int] | None,
    reference_color_order: Sequence[int] | None,
    base: BaseProcess,
    color_accuracy: str,
    model_source: str = "built-in-sm",
    model_profile: str = "built-in-sm",
    selected_jit_opt_level: int = 1,
) -> dict[str, Any]:
    numerical_current_passes = _matrix_uses_numerical_current_passes(
        base,
        color_accuracy=color_accuracy,
    )
    settings: dict[str, Any] = {
        "runtime": "rusticol",
        "model_source": str(model_source),
        "model_profile": str(model_profile),
        "coupling_order_policy": (
            "minimal" if model_profile != "built-in-sm" else "all"
        ),
        "external_sm_qcd_only_cap": (
            "QED=0"
            if model_profile != "built-in-sm"
            and base.key in QCD_ONLY_EXTERNAL_SM_FAMILIES
            else None
        ),
        "precision": 16,
        "color_accuracy": color_accuracy,
        "target_runtime_s": float(target_runtime),
        "batch_size": (
            DEFAULT_BATCH_SIZE
            if color_accuracy == "lc"
            else DEFAULT_NON_LC_BATCH_SIZE
        ),
        "selected_runtime_batch_size": (
            DEFAULT_BATCH_SIZE
            if color_accuracy == "lc"
            else DEFAULT_NON_LC_BATCH_SIZE
        ),
        "all_flow_runtime_batch_size": (
            DEFAULT_ALL_FLOW_BATCH_SIZE if color_accuracy == "lc" else None
        ),
        "n_cores": int(n_cores),
        "symbolica_n_cores": int(n_cores),
        "symbolica_stage_local_parameter_layout": (
            DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT
        ),
        "symbolica_output_chunk_size": DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
        "selected_symbolica_output_chunk_size": DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
        "selected_symbolica_output_chunk_strategy": (
            "auto" if backend_key == "jit" and color_accuracy == "lc" else "uniform"
        ),
        "all_flow_symbolica_output_chunk_size": (
            DEFAULT_ALL_FLOW_SYMBOLICA_OUTPUT_CHUNK_SIZE
            if color_accuracy == "lc"
            else None
        ),
        "all_flow_symbolica_output_chunk_strategy": (
            "uniform" if color_accuracy == "lc" else None
        ),
        "generation_selected_lc_sector_ids": None,
        "runtime_lc_sector_ids": (
            None
            if selected_lc_sector_ids is None
            else sorted(int(sector_id) for sector_id in selected_lc_sector_ids)
        ),
        "runtime_lc_sector_selection_policy": (
            "matched-first-amplicol-group-integral"
            if color_accuracy == "lc"
            else None
        ),
        "runtime_lc_sector_selector": (
            "selected-flow-helicity-summed-artifact-plus-fixed-helicity-materialized-all-sector-artifact"
            if color_accuracy == "lc"
            else None
        ),
        "all_flow_runtime_probe": (
            "time-process-materialized-lc-all-sector-fixed-source-helicity"
            if color_accuracy == "lc"
            else None
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
        "zero_current_filter": {
            "enabled": numerical_current_passes,
            "sample_count": 10,
            "seed": 12345,
            "relative_tolerance": 1.0e-12,
            "zero_tolerance": 1.0e-300,
        },
        "current_merging": {
            "enabled": numerical_current_passes,
            "sample_count": 10,
            "seed": 12345,
            "relative_tolerance": 1.0e-12,
            "zero_tolerance": 1.0e-300,
        },
    }
    if backend_key == "jit":
        settings.update(
            {
                "symbolica_evaluator_backend": "jit",
                "symbolica_jit_optimization_level": selected_jit_opt_level,
                "selected_symbolica_jit_optimization_level": (
                    selected_jit_opt_level
                ),
                "all_flow_symbolica_jit_optimization_level": (
                    1 if color_accuracy == "lc" else None
                ),
                "symbolica_iterations": DEFAULT_SYMBOLICA_ITERATIONS,
                "symbolica_cpe_iterations": None,
                "symbolica_max_horner_scheme_variables": 1000,
                "symbolica_max_common_pair_cache_entries": 5000000,
                "symbolica_max_common_pair_distance": 1000,
                "symbolica_compiled_chunk_compile_workers": int(n_cores),
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


def _matrix_uses_numerical_current_passes(
    base: BaseProcess,
    *,
    color_accuracy: str,
) -> bool:
    return color_accuracy == "lc" and base.lc_sector_strategy != "all"


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


def _matrix_model_profile(model_source: str | Path) -> str:
    return (
        "built-in-sm"
        if str(model_source).casefold() in {"builtin_sm", "built-in-sm"}
        else "external-sm"
    )


def _normalized_model_source(model_source: str | Path) -> str:
    text = str(model_source)
    if text.casefold() in {"builtin_sm", "built-in-sm"}:
        return "BUILTIN_SM"
    return str(Path(text).expanduser().resolve())


def _prepare_benchmark_model_source(
    model_source: str,
    *,
    output_root: Path,
) -> tuple[str, dict[str, Any]]:
    from pyamplicol.model_source import compile_model_source

    started = time.perf_counter()
    compiled = compile_model_source(
        model_source,
        cache_dir=output_root / ".model_cache",
        use_cache=True,
        require_supported=True,
    )
    compiled_path = compiled.write(
        output_root / "external-sm.pyAmplicol-model.json"
    ).resolve()
    wall_s = time.perf_counter() - started
    print(
        f"[matrix] external model ready in {wall_s:.3f}s: {compiled_path}",
        flush=True,
    )
    return str(compiled_path), {
        "source": model_source,
        "compiled_model": str(compiled_path),
        "load_or_compile_wall_s": wall_s,
        "source_conversion_s": float(compiled.conversion_seconds),
        "phase_timings": dict(compiled.phase_timings),
        "compiled_model_name": compiled.name,
        "source_metadata": dict(compiled.source),
    }


def _matrix_model_label(model_profile: str) -> str:
    return "built-in-sm" if model_profile == "built-in-sm" else "external UFO/JSON SM"


def _default_data_path(color_accuracy: str, *, model_profile: str) -> Path:
    if model_profile != "built-in-sm":
        if color_accuracy == "nlc":
            return DEFAULT_UFO_SM_NLC_DATA
        if color_accuracy == "full":
            return DEFAULT_UFO_SM_FULL_DATA
        return DEFAULT_UFO_SM_DATA
    if color_accuracy == "nlc":
        return DEFAULT_NLC_DATA
    if color_accuracy == "full":
        return DEFAULT_FULL_DATA
    return DEFAULT_DATA


def _default_table_path(color_accuracy: str, *, model_profile: str) -> Path:
    if model_profile != "built-in-sm":
        if color_accuracy == "nlc":
            return DEFAULT_UFO_SM_NLC_TABLE
        if color_accuracy == "full":
            return DEFAULT_UFO_SM_FULL_TABLE
        return DEFAULT_UFO_SM_TABLE
    if color_accuracy == "nlc":
        return DEFAULT_NLC_TABLE
    if color_accuracy == "full":
        return DEFAULT_FULL_TABLE
    return DEFAULT_TABLE


def _seed_reference_cache(
    target: dict[str, Any],
    reference: dict[str, Any],
    *,
    color_accuracy: str,
    validation_report: Path | None,
) -> None:
    reference_entries = reference.get("entries", {})
    if not isinstance(reference_entries, dict):
        raise TypeError("reference matrix data 'entries' must be a JSON object")
    target_entries = target.setdefault("entries", {})
    if not isinstance(target_entries, dict):
        raise TypeError("target matrix data 'entries' must be a JSON object")
    fixture_results = _validation_report_results(validation_report)
    for base_key, raw_base_entries in reference_entries.items():
        if not isinstance(raw_base_entries, dict):
            continue
        seeded_base = target_entries.setdefault(str(base_key), {})
        if not isinstance(seeded_base, dict):
            raise TypeError(f"target entry for {base_key!r} must be a JSON object")
        for n_text, raw_case in raw_base_entries.items():
            if not isinstance(raw_case, dict):
                continue
            seeded_case = seeded_base.setdefault(str(n_text), {})
            if not isinstance(seeded_case, dict):
                raise TypeError(
                    f"target entry for {base_key!r}, n={n_text!r} must be an object"
                )
            for key in ("base_key", "n_final", "process"):
                if key in raw_case:
                    seeded_case[key] = copy.deepcopy(raw_case[key])
            amplicol = raw_case.get("amplicol")
            if isinstance(amplicol, dict):
                seeded_case["amplicol"] = copy.deepcopy(amplicol)
            fixture = fixture_results.get(
                f"{color_accuracy}:{base_key}:n{n_text}"
            )
            if fixture is not None:
                relative = _optional_float(fixture.get("relative_difference"))
                seeded_case["validation"] = {
                    "status": str(fixture.get("status", "unknown")),
                    "kind": "ufo-sm-vs-archived-built-in-same-point",
                    "max_relative_difference": relative,
                    "relative_differences": {"pyamplicol_jit": relative},
                    "tolerance": _optional_float(fixture.get("tolerance"))
                    or 1.0e-10,
                    "source_report": str(validation_report),
                }
            if "status" not in seeded_case:
                seeded_case["status"] = "reference_seeded"


def _validation_report_results(path: Path | None) -> dict[str, dict[str, Any]]:
    if path is None:
        return {}
    report = _load_data(path)
    raw_results = report.get("results", ())
    if not isinstance(raw_results, list):
        raise TypeError("validation report 'results' must be a JSON array")
    return {
        str(result["id"]): result
        for result in raw_results
        if isinstance(result, dict) and "id" in result
    }


def _case_payload(data: dict[str, Any], base: BaseProcess, n_final: int) -> dict[str, Any]:
    _refresh_data_metadata(
        data,
        color_accuracy=str(data.get("color_accuracy", "lc")),
    )
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


def _refresh_data_metadata(data: dict[str, Any], *, color_accuracy: str) -> None:
    data.setdefault("schema_version", 1)
    data.setdefault("created_by", "pyAmpliCol/docs/result_matrix.py")
    data["color_accuracy"] = color_accuracy
    data["updated_at"] = _now()
    data["base_processes"] = [
        base_payload(process_id, item)
        for process_id, item in enumerate(BASE_PROCESSES, start=1)
    ]
    data.setdefault("entries", {})


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
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        print(f"[matrix] warning: {path} is empty; starting a fresh cache", flush=True)
        return {}
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        backup = path.with_name(f"{path.name}.invalid")
        backup.write_text(text, encoding="utf-8")
        print(
            f"[matrix] warning: could not parse {path}: {exc}; "
            f"copied invalid cache to {backup} and starting fresh",
            flush=True,
        )
        return {}


def _write_data(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


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
        "jit_compile_s",
        "jit_fraction_of_generation",
        "key",
        "kind",
        "lowering_status",
        "manifest",
        "process",
        "process_dir",
        "requested_process_dir",
        "runtime_process_dir",
        "runtime_lc_sector_artifact",
        "_command_elapsed_s",
        "_progress_log",
        "_generation_samples",
        "_generation_sample_count",
        "_generation_median_sample_index",
        "_artifact_sample_index",
        "_artifact_command_elapsed_s",
        "_artifact_internal_generation_s",
        "runtime_available",
        "runtime_backend",
        "runtime_unavailable_message",
        "schema_version",
        "target_runtime_s",
        "lc_sector_ids",
        "effective_lc_sector_ids",
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
    if isinstance(exc, MemoryLimitExceeded):
        return {
            "status": "ram_limit",
            "error": error,
            "memory_limit_gb": exc.limit_gb,
            "peak_memory_gb": exc.peak_memory_gb,
            "peak_memory_bytes": exc.peak_memory_bytes,
            "memory_metric": exc.memory_metric,
            "peak_rss_gb": exc.peak_rss_gb,
            "peak_rss_bytes": exc.peak_rss_bytes,
            "peak_physical_footprint_gb": exc.peak_physical_footprint_gb,
            "peak_physical_footprint_bytes": (
                exc.peak_physical_footprint_bytes
            ),
            "finished_at": _now(),
        }
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
    return _payload_first_value(payload)


def _payload_first_value(payload: dict[str, Any]) -> float | None:
    values = payload.get("values")
    if not isinstance(values, list) or not values:
        return None
    return _optional_float(values[0])


def _relative_difference(reference: float, value: float) -> float:
    return abs(value - reference) / max(abs(reference), abs(value), 1.0e-300)


def _format_seconds(value: float) -> str:
    return rf"\texttt{{{_format_sig(value)}}}"


def _format_us(value: float) -> str:
    return rf"\texttt{{{_format_sig(value)}}}"


def _format_sig(value: float, *, unit: str = "") -> str:
    if value == 0.0:
        text = "0"
        return f"{text} {unit}" if unit else text
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
    return f"{text} {unit}" if unit else text


def _format_ratio(value: float) -> str:
    if value >= 100:
        return f"{value:.0f}"
    if value >= 10:
        return f"{value:.1f}"
    if value >= 1:
        return f"{value:.2f}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _scientific_latex(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return r"\texttt{N/A}"
    if value == 0.0:
        return r"\(0\)"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10.0**exponent)
    return rf"\({_format_compact_number(mantissa)}\times10^{{{exponent}}}\)"


def _compact_latex_number(value: float) -> str:
    text = _format_compact_number(value)
    if "e" not in text:
        return text
    mantissa, exponent = text.split("e", maxsplit=1)
    return rf"\({mantissa}\times10^{{{int(exponent)}}}\)"


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
