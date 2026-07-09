from __future__ import annotations

import argparse
import contextlib
import json
import math
import os
import queue
import signal
import shutil
import statistics
import subprocess
import sys
import threading
import time
from decimal import Decimal
from pathlib import Path
from pprint import pprint
from typing import Any, Literal, Mapping, Sequence, cast

from . import __version__
from .cli_display import (
    DisplayColumn,
    DisplayRow,
    default_display_for_args,
    format_measurement,
)
from .core_types import NativeEvaluationError
from .processes import (
    ANTI_PARTICLE,
    PDGS,
    ProcessEnumeration,
    ProcessEnumerator,
    ProcessOptions,
    ProcessSelectionReport,
    ProcessSetEntry,
    ProcessSetEnumeration,
    build_generic_process_selection_report,
    enumerate_process_set,
)
from .reference import (
    AmplicolAdapter,
    CommandResult,
    amplicol_process_file_entry,
    reference_color_order_for_run,
)

_RUNTIME_BACKENDS = (
    "auto",
    "python",
    "dag",
    "numeric-tensor-network",
    "rusticol",
)
_PRODUCTION_RUNTIME_BACKENDS = ("rusticol",)
RuntimeBackend = Literal[
    "auto",
    "python",
    "dag",
    "numeric-tensor-network",
    "rusticol",
]

_CHILD_PROGRESS_ENV = "PYAMPLICOL_CHILD_PROGRESS"
_CHILD_PROGRESS_PREFIX = "::pyamplicol-progress::"
_DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE = 128
_DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT = True


_PROCESS_SET_STANDALONE_CHECK_SCRIPT = r'''#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import runpy
import sys
from pathlib import Path


def resolve_process(root: Path, selected: str | None):
    manifest = json.loads((root / "process_set_manifest.json").read_text())
    entries = manifest.get("processes", [])
    if not isinstance(entries, list) or not entries:
        raise SystemExit(f"process-set artifact contains no subprocesses: {root}")
    selected_key = selected or str(manifest.get("default_process_key"))
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        key = str(entry.get("key"))
        process = str(entry.get("process"))
        if selected_key in {key, process}:
            path = Path(str(entry.get("path", "")))
            process_dir = path if path.is_absolute() else root / path
            return entry, process_dir
    available = ", ".join(
        str(entry.get("key"))
        for entry in entries
        if isinstance(entry, dict) and "key" in entry
    )
    raise SystemExit(
        f"process {selected_key!r} not found in {root}; available: {available}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Standalone rusticol process-set check")
    parser.add_argument(
        "--process",
        help=(
            "Canonical process key or process string to select from this "
            "process-set artifact. Defaults to the manifest default."
        ),
    )
    parser.add_argument("--precision", type=int, default=16)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--target-runtime", type=float, default=10.0)
    parser.add_argument(
        "--rusticol-folder",
        type=Path,
        help=(
            "Optional rusticol location forwarded to the selected subprocess "
            "checker."
        ),
    )
    args, passthrough = parser.parse_known_args()

    root = Path(__file__).resolve().parent
    _, process_dir = resolve_process(root, args.process)
    checker = process_dir / "check_standalone.py"
    if not checker.exists():
        raise SystemExit(f"selected subprocess has no standalone checker: {checker}")

    forwarded = [
        str(checker),
        "--precision",
        str(args.precision),
        "--target-runtime",
        str(args.target_runtime),
    ]
    if args.profile:
        forwarded.append("--profile")
    if args.rusticol_folder is not None:
        forwarded.extend(["--rusticol-folder", str(args.rusticol_folder)])
    forwarded.extend(passthrough)
    sys.argv = forwarded
    try:
        runpy.run_path(str(checker), run_name="__main__")
    except SystemExit as exc:
        code = exc.code
        if code is None:
            return 0
        if isinstance(code, int):
            return code
        print(code, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
'''


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="pyamplicol matrix-element generation and validation tooling."
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the pyamplicol package version and exit.",
    )

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    subparsers.add_parser(
        "inspect",
        help="Inspect the managed Symbolica/spenso/idenso environment.",
    )

    processes = subparsers.add_parser(
        "processes",
        help="Enumerate a process and optionally export legacy processes.txt.",
    )
    _add_process_options(processes)
    processes.add_argument("process", metavar="PROCESS")
    processes.add_argument("--legacy-output", type=Path)
    processes.add_argument("--json", action="store_true")

    process_plan = subparsers.add_parser(
        "process-plan",
        help=(
            "Write generic schema-v2 current-planning manifests without "
            "building runtime evaluators."
        ),
    )
    _add_process_options(process_plan)
    process_plan.add_argument("process", metavar="PROCESS")
    process_plan.add_argument("output_dir", type=Path, metavar="OUTPUT_DIR")
    process_plan.add_argument(
        "--color-accuracy",
        choices=("lc", "nlc", "full"),
        default="lc",
        help=(
            "Colour treatment to record in the planning manifest. LC is the "
            "default; NLC/full attach sparse colour-contraction metadata when "
            "supported."
        ),
    )
    process_plan.add_argument(
        "--max-currents",
        type=int,
        default=20000,
        help=(
            "Safety cap for generic current-plan construction; use a negative "
            "value to disable this internal cap."
        ),
    )
    process_plan.add_argument(
        "--max-color-sectors",
        type=int,
        default=20000,
        help=(
            "Safety cap for generic colour-flow sector enumeration; use a "
            "negative value to disable this internal cap."
        ),
    )
    _add_generic_dag_pruning_options(process_plan)
    process_plan.add_argument("--json", action="store_true")

    generate = subparsers.add_parser(
        "generate",
        help=argparse.SUPPRESS,
    )
    _add_process_options(generate)
    generate.add_argument("process", metavar="PROCESS")
    generate.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    generate.add_argument("--no-cache", action="store_true")
    _add_evaluator_build_options(generate)
    generate.add_argument("--json", action="store_true")

    evaluate = subparsers.add_parser(
        "evaluate",
        help=argparse.SUPPRESS,
    )
    _add_process_options(evaluate)
    evaluate.add_argument("process", metavar="PROCESS")
    evaluate.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    evaluate.add_argument("--sqrt-s", type=float)
    _add_runtime_backend_option(evaluate)
    _add_evaluator_build_options(evaluate)
    evaluate.add_argument("--json", action="store_true")

    compare = subparsers.add_parser(
        "compare-amplicol",
        help="Drive legacy AmpliCol reference generation/evaluation.",
    )
    _add_process_options(compare)
    compare.add_argument("process", metavar="PROCESS")
    compare.add_argument("--points", type=int, default=10)
    compare.add_argument("--timing", type=int)
    compare.add_argument("--amplicol-root", type=Path, default=_default_amplicol_root())
    compare.add_argument("--mg5-path", type=Path, default=_default_mg5_path())
    compare.add_argument("--jobs", type=int, default=8)
    compare.add_argument("--timeout", type=float)
    compare.add_argument("--process-file", type=Path)
    compare.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    _add_runtime_backend_option(compare, production=True)
    _add_evaluator_build_options(compare)
    _add_generic_dag_pruning_options(compare)
    compare.add_argument("--skip-build", action="store_true")
    compare.add_argument(
        "--amplicol-probe",
        action="store_true",
        help=(
            "Use AmpliCol's built-in generated-library random/fixed probe "
            "instead of the default deterministic supplied-momenta probe. "
            "When building is enabled, this creates and uses the generated "
            "Fortran amplitude library."
        ),
    )
    compare.add_argument(
        "--amplicol-momenta-probe",
        action="store_true",
        help=(
            "Use a deterministic pyAmpliCol validation point and evaluate it "
            "with the generated Fortran AmpliCol library. This is the default."
        ),
    )
    compare.add_argument(
        "--me-test",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    compare.add_argument("--dry-run", action="store_true")
    compare.add_argument("--json", action="store_true")

    family = subparsers.add_parser(
        "validate-z-gluon-family",
        help=argparse.SUPPRESS,
    )
    _add_process_options(family)
    family.add_argument("--min-gluons", type=int, default=0)
    family.add_argument("--max-gluons", type=int, default=6)
    family.add_argument("--points", type=int, default=10)
    family.add_argument("--rel-tol", type=float, default=1.0e-8)
    family.add_argument("--timing", type=int)
    family.add_argument("--amplicol-root", type=Path, default=_default_amplicol_root())
    family.add_argument("--jobs", type=int, default=8)
    family.add_argument("--timeout", type=float)
    family.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    _add_runtime_backend_option(family)
    _add_evaluator_build_options(family)
    family.add_argument("--dry-run", action="store_true")
    family.add_argument("--json", action="store_true")

    benchmark = subparsers.add_parser(
        "benchmark-z-gluon-modes",
        help=argparse.SUPPRESS,
    )
    _add_process_options(benchmark)
    benchmark.add_argument("--min-gluons", type=int, default=0)
    benchmark.add_argument("--max-gluons", type=int, default=6)
    benchmark.add_argument("--points", type=int, default=3)
    benchmark.add_argument("--timing", type=int)
    benchmark.add_argument("--amplicol-root", type=Path, default=_default_amplicol_root())
    benchmark.add_argument("--jobs", type=int, default=8)
    benchmark.add_argument("--timeout", type=float)
    benchmark.add_argument("--numeric-timeout", type=float, default=120.0)
    benchmark.add_argument("--parametric-timeout", type=float, default=180.0)
    benchmark.add_argument("--parametric-max-gluons", type=int, default=4)
    benchmark.add_argument(
        "--only-legacy-shared",
        action="store_true",
        help="Benchmark only legacy AmpliCol and the shared-current D-mode.",
    )
    _add_evaluator_build_options(benchmark)
    benchmark.add_argument(
        "--tensor-strategy",
        choices=("interleaved", "monolithic"),
        default="interleaved",
    )
    benchmark.add_argument("--output", type=Path)
    benchmark.add_argument("--json", action="store_true")

    tensor_profile = subparsers.add_parser(
        "profile-tensor-evaluator",
        help=argparse.SUPPRESS,
    )
    tensor_profile.add_argument("process", metavar="PROCESS")
    tensor_profile.add_argument("--sqrt-s", type=float, default=1000.0)
    tensor_profile.add_argument("--repetitions", type=int, default=100)
    tensor_profile.add_argument("--evaluator-repetitions", type=int)
    tensor_profile.add_argument(
        "--tensor-strategy",
        choices=("interleaved", "monolithic"),
        default="interleaved",
    )
    tensor_profile.add_argument("--json", action="store_true")

    dag_profile = subparsers.add_parser(
        "profile-dag-evaluator",
        help=argparse.SUPPRESS,
    )
    dag_profile.add_argument("process", metavar="PROCESS")
    dag_profile.add_argument("--sqrt-s", type=float, default=1000.0)
    dag_profile.add_argument("--points", type=int, default=16)
    dag_profile.add_argument("--repetitions", type=int, default=100)
    dag_profile.add_argument(
        "--generate-only",
        action="store_true",
        help=(
            "Generate and save a reusable DAG process artifact, then exit "
            "without timing runtime evaluation. This is primarily for the "
            "rusticol two-stage workflow."
        ),
    )
    _add_runtime_backend_option(dag_profile)
    _add_evaluator_build_options(dag_profile)
    _set_fast_rusticol_dag_defaults(dag_profile)
    dag_profile.add_argument("--json", action="store_true")

    generate_process = subparsers.add_parser(
        "generate-process",
        help=(
            "Generate the default fast Rusticol eager-DAG process directory. "
            "Defaults to JIT O3 stage evaluators."
        ),
    )
    _add_process_options(generate_process)
    generate_process.add_argument("process", metavar="PROCESS")
    generate_process.add_argument(
        "output_dir",
        type=Path,
        nargs="?",
        metavar="OUTPUT_DIR",
    )
    generate_process.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Enumerate and display selected concrete subprocesses without "
            "writing an artifact directory."
        ),
    )
    generate_process.add_argument(
        "--no-enumeration-prefilter",
        dest="enumeration_prefilter",
        action="store_false",
        default=True,
        help=(
            "Use the slower legacy-compatible inclusive enumeration path for "
            "diagnostics and before/after timing comparisons."
        ),
    )
    generate_process.add_argument(
        "--append",
        action="store_true",
        help="Append new process entries to an existing process-set output.",
    )
    generate_process.add_argument(
        "--replace",
        action="store_true",
        help="Replace existing process entries with matching canonical keys.",
    )
    generate_process.add_argument(
        "--color-accuracy",
        choices=("lc", "nlc", "full"),
        default="lc",
        help=(
            "Colour treatment requested for process generation. LC is the "
            "default; NLC/full use the same generic DAG with a sparse final "
            "colour contraction when supported."
        ),
    )
    generate_process.add_argument(
        "--n_cores",
        "--n-cores",
        "-n-cores",
        dest="n_cores",
        type=int,
        default=1,
        help=(
            "Number of subprocess artifacts to generate concurrently for a "
            "process-set request. Each subprocess still uses the Symbolica "
            "core settings specified by --symbolica-n-cores."
        ),
    )
    generate_process.add_argument(
        "--max-currents",
        type=int,
        default=50000,
        help=(
            "Safety cap for generic current-plan construction; use a negative "
            "value to disable this internal cap."
        ),
    )
    generate_process.add_argument(
        "--max-color-sectors",
        type=int,
        default=20000,
        help=(
            "Safety cap for generic colour-flow sector enumeration; use a "
            "negative value to disable this internal cap."
        ),
    )
    generate_process.add_argument(
        "--monitor",
        action="store_true",
        help=(
            "Emit fixed-width generation progress updates to stderr. This is "
            "safe to combine with --json because stdout remains machine JSON."
        ),
    )
    _add_evaluator_build_options(generate_process)
    _add_generic_dag_pruning_options(generate_process)
    _set_fast_rusticol_dag_defaults(generate_process)
    generate_process.add_argument("--json", action="store_true")

    time_process = subparsers.add_parser(
        "time-process",
        help=(
            "Load a generated process directory with Rusticol and time the "
            "validation momenta. Defaults to precision 16 and a 10 second target."
        ),
    )
    time_process.add_argument("process_dir", type=Path, metavar="PROCESS_DIR")
    time_process.add_argument(
        "--process",
        dest="process_key",
        help=(
            "Canonical process key or process string to select from a "
            "process-set artifact. Defaults to the first entry."
        ),
    )
    time_process.add_argument("--precision", type=int, default=16)
    time_process.add_argument("--target-runtime", type=float, default=10.0)
    time_process.add_argument(
        "--batch-size",
        type=int,
        help=(
            "Number of samples per repeated timing block. Defaults to the "
            "batch size recorded in the process artifact, which matches the "
            "performance-summary benchmark setup."
        ),
    )
    time_process.add_argument(
        "--model-parameters",
        type=Path,
        help=(
            "TOML file with runtime model-parameter overrides consumed by "
            "Rusticol. Keys must match the artifact model_parameter_names metadata."
        ),
    )
    time_process.add_argument("--json", action="store_true")

    profile = subparsers.add_parser(
        "profile",
        help=argparse.SUPPRESS,
    )
    _add_process_options(profile)
    profile.add_argument("process", metavar="PROCESS")
    profile.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    _add_runtime_backend_option(profile)
    _add_evaluator_build_options(profile)
    profile.add_argument("--json", action="store_true")

    for legacy_command in (
        "generate",
        "evaluate",
        "validate-z-gluon-family",
        "benchmark-z-gluon-modes",
        "profile-tensor-evaluator",
        "profile-dag-evaluator",
        "profile",
    ):
        _hide_subcommand_from_help(subparsers, legacy_command)

    return parser.parse_args(argv)


def _hide_subcommand_from_help(
    subparsers: Any,
    name: str,
) -> None:
    """Keep a compatibility subcommand parseable while removing it from help."""

    choices_actions = getattr(subparsers, "_choices_actions", None)
    if not isinstance(choices_actions, list):
        return
    choices_actions[:] = [
        action for action in choices_actions if getattr(action, "dest", None) != name
    ]


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.version:
        print(__version__)
        return 0
    args._display = default_display_for_args(args)

    if args.command in (None, "inspect"):
        return _cmd_inspect()
    if args.command == "processes":
        return _cmd_processes(args)
    if args.command == "process-plan":
        return _cmd_process_plan(args)
    if args.command == "generate":
        return _cmd_generate(args)
    if args.command == "evaluate":
        return _cmd_evaluate(args)
    if args.command == "compare-amplicol":
        return _cmd_compare_amplicol(args)
    if args.command == "validate-z-gluon-family":
        return _cmd_validate_z_gluon_family(args)
    if args.command == "benchmark-z-gluon-modes":
        return _cmd_benchmark_z_gluon_modes(args)
    if args.command == "profile-tensor-evaluator":
        return _cmd_profile_tensor_evaluator(args)
    if args.command == "profile-dag-evaluator":
        return _cmd_profile_dag_evaluator(args)
    if args.command == "generate-process":
        return _cmd_generate_process(args)
    if args.command == "time-process":
        return _cmd_time_process(args)
    if args.command == "profile":
        return _cmd_profile(args)
    raise ValueError(f"unknown command: {args.command}")


def _add_process_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-FS", "--flavour-scheme", type=int, default=5)
    parser.add_argument("-3", "--include-3qqbar", action="store_true")
    parser.add_argument("-cc", "--include-cc", action="store_true")
    parser.add_argument("-res", "--include-resonance", action="store_true")
    parser.add_argument(
        "--parallel-process-enumeration",
        action="store_true",
        help="Record that legacy-compatible enumeration may be parallelized later.",
    )


def _add_generic_dag_pruning_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--max-coupling-order",
        action="append",
        default=[],
        metavar="NAME=N",
        help=(
            "Generic model-level coupling-order cap, repeatable. Examples: "
            "--max-coupling-order QED=1 --max-coupling-order QCD=4. "
            "This is not a process-family tag; orders are accumulated from "
            "local model vertices."
        ),
    )
    parser.add_argument(
        "--max-qcd-order",
        type=int,
        help="Shortcut for --max-coupling-order QCD=N.",
    )
    parser.add_argument(
        "--max-qed-order",
        type=int,
        help="Shortcut for --max-coupling-order QED=N.",
    )
    parser.add_argument(
        "--coupling-order-policy",
        choices=("all", "minimal"),
        default="all",
        help=(
            "Generic coupling-order pruning policy. 'all' keeps every "
            "model-reachable order allowed by explicit caps. 'minimal' first "
            "infers the lowest contributing model coupling-order envelope and "
            "uses it as an additional cap; this is useful for fast leading-order "
            "generation without naming a process family."
        ),
    )
    parser.add_argument(
        "--max-lc-current-line-groups",
        type=int,
        help=(
            "Generic leading-colour pruning cap on how many colour-line groups "
            "one intermediate current may span. This can be used to restrict "
            "multi-quark-line warmups without naming a process family."
        ),
    )
    parser.add_argument(
        "--max-quark-lines",
        type=int,
        help=(
            "User-facing alias for --max-quark-pairs: generic cap on the "
            "number of external quark-antiquark colour lines allowed in a "
            "subprocess."
        ),
    )
    parser.add_argument(
        "--max-quark-pairs",
        type=int,
        help=(
            "Generic cap on the number of external quark-antiquark colour "
            "lines allowed in a subprocess. Useful for pruning inclusive "
            "p/j requests without naming a process family."
        ),
    )
    parser.add_argument(
        "--ignore-particles",
        default="",
        metavar="LIST",
        help=(
            "Comma-separated particle names or PDG ids to remove from generic "
            "DAG construction, e.g. 'h,a' or '25,22'."
        ),
    )
    parser.add_argument(
        "--ignore-vertex-kinds",
        default="",
        metavar="LIST",
        help="Comma-separated model vertex-kind ids to remove during DAG construction.",
    )
    parser.add_argument(
        "--no-closure-side-mask-pruning",
        dest="closure_side_mask_pruning",
        action="store_false",
        default=True,
        help=(
            "Disable generic closure-side subset pruning. This is intended for "
            "diagnostics; production generation keeps it enabled."
        ),
    )
    parser.add_argument(
        "--no-color-order-mask-pruning",
        dest="color_order_mask_pruning",
        action="store_false",
        default=True,
        help=(
            "Disable generic LC colour-order subset pruning. This is intended "
            "for diagnostics; production generation keeps it enabled."
        ),
    )
    parser.add_argument(
        "--no-species-reachability-pruning",
        dest="species_reachability_pruning",
        action="store_false",
        default=True,
        help=(
            "Disable generic particle-id reachability pruning. This is "
            "intended for diagnostics; production generation keeps it enabled."
        ),
    )
    parser.add_argument(
        "--numerical-filter-current",
        dest="numerical_filter_current",
        action="store_true",
        default=None,
        help=(
            "Force generation-time numerical warmup pruning of currents. "
            "By default this is enabled for LC and disabled for NLC/full "
            "all-sector colour contractions."
        ),
    )
    parser.add_argument(
        "--no-numerical-filter-current",
        dest="numerical_filter_current",
        action="store_false",
        default=None,
        help=(
            "Disable generation-time numerical warmup pruning of currents "
            "that are zero on deterministic random phase-space points."
        ),
    )
    parser.add_argument(
        "--numerical-current-merging",
        dest="numerical_current_merging",
        action="store_true",
        default=None,
        help=(
            "Force generation-time numerical detection and merging of "
            "identical currents. By default this is enabled for LC and "
            "disabled for NLC/full all-sector colour contractions."
        ),
    )
    parser.add_argument(
        "--no-numerical-current-merging",
        dest="numerical_current_merging",
        action="store_false",
        default=None,
        help=(
            "Disable generation-time numerical detection and merging of "
            "identical current values on deterministic random phase-space "
            "points."
        ),
    )
    parser.add_argument(
        "--numerical-current-samples",
        type=int,
        default=10,
        help=(
            "Number of deterministic random phase-space points used by the "
            "generation-time numerical current filter and merger."
        ),
    )
    parser.add_argument(
        "--numerical-current-seed",
        type=int,
        default=12345,
        help="Base seed for the generation-time numerical current warmup.",
    )
    parser.add_argument(
        "--numerical-current-relative-tolerance",
        type=float,
        default=1.0e-12,
        help=(
            "Relative threshold used by generation-time numerical current "
            "zero filtering and identical-current merging."
        ),
    )
    parser.add_argument(
        "--numerical-current-zero-tolerance",
        type=float,
        default=1.0e-300,
        help=(
            "Absolute floor used by generation-time numerical current zero "
            "filtering and identical-current merging."
        ),
    )
    parser.add_argument(
        "--lc-sector-strategy",
        choices=(
            "all",
            "topology-representatives",
            "line-pairing-representatives",
            "reference",
        ),
        default="topology-representatives",
        help=(
            "Leading-colour generation strategy. 'all' materializes every LC "
            "sector. 'topology-representatives' is the default: it compiles a "
            "single representative for replay-safe isomorphic LC topologies "
            "and falls back to all sectors when the topology is not safe. "
            "'line-pairing-representatives' compiles one representative per "
            "open quark-line pairing/allocation, dropping only permutations of "
            "whole open-line blocks. "
            "'reference' compiles only LC sector 0. Sibling-sector replay is "
            "opt-in via --lc-topology-replay."
        ),
    )
    parser.add_argument(
        "--lc-sector-ids",
        default="",
        metavar="LIST",
        help=(
            "Comma-separated explicit LC colour-sector ids to materialize. "
            "This is a generic colour-flow pruning hint for large multi-line "
            "processes and overrides --lc-sector-strategy when supplied."
        ),
    )
    parser.add_argument(
        "--reference-color-order",
        default="",
        metavar="LIST",
        help=(
            "Comma-separated external labels of a reference LC colour order. "
            "This is used for generic AmpliCol row comparisons so closure "
            "endpoints follow the supplied ordered amplitude."
        ),
    )
    parser.add_argument(
        "--lc-topology-replay",
        action="store_true",
        help=(
            "After compiling topology representatives, evaluate all replay-safe "
            "sibling LC sectors by external-label permutation and sum them. "
            "This is useful for full-sector diagnostics; the default matches "
            "the Fortran AmpliCol reference-order convention."
        ),
    )


def _add_runtime_backend_option(
    parser: argparse.ArgumentParser,
    *,
    production: bool = False,
) -> None:
    choices = _PRODUCTION_RUNTIME_BACKENDS if production else _RUNTIME_BACKENDS
    default = "rusticol" if production else "auto"
    help_text = (
        "Production runtime backend. The generic DAG workflow currently uses "
        "Rusticol process artifacts."
        if production
        else (
            "Legacy native runtime backend for compatibility commands. "
            "Production generation uses `generate-process` and Rusticol "
            "process artifacts."
        )
    )
    parser.add_argument(
        "--runtime-backend",
        choices=choices,
        default=default,
        help=help_text,
    )


def _add_evaluator_build_options(
    parser: argparse.ArgumentParser,
    *,
    include_legacy_compiled_dag_options: bool = False,
) -> None:
    parser.set_defaults(
        merge_evaluators_strategy=False,
        symbolica_direct_translation=True,
        symbolica_compiled_native=True,
        symbolica_jit_direct_translation=None,
    )
    parser.add_argument(
        "--merge-evaluators-strategy",
        dest="merge_evaluators_strategy",
        action="store_true",
        help=(
            "Build bounded evaluator pieces and merge them iteratively. This "
            "keeps peak memory lower but is currently not the default."
        ),
    )
    parser.add_argument(
        "--no-merge-evaluators-strategy",
        dest="merge_evaluators_strategy",
        action="store_false",
        help=(
            "Build bounded multi-output evaluator groups in a single call where "
            "supported. This can use more RAM but is currently the default."
        ),
    )
    parser.add_argument(
        "--symbolica-split-vertex-current-stages",
        action="store_true",
        help=(
            "Compile separate generic vertex and current-combine evaluators "
            "for each current-size stage. The default keeps the established "
            "generic Rusticol eager-DAG stage layout."
        ),
    )
    parser.add_argument(
        "--symbolica-stage-local-parameter-layout",
        action=argparse.BooleanOptionalAction,
        default=_DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT,
        help=(
            "Compile generic Rusticol stage evaluators with only the value, "
            "momentum, and model-parameter inputs required by that stage. This "
            "can reduce Symbolica evaluator input pressure at the cost of "
            "runtime stage-input packing. Enabled by default; use "
            "--no-symbolica-stage-local-parameter-layout to build global-input "
            "stage evaluators."
        ),
    )
    parser.set_defaults(compiled_dag_inline_external_wavefunctions=True)
    parser.set_defaults(compiled_dag_helicity_filter=True)
    parser.set_defaults(
        compiled_dag_lowering="spenso",
        compiled_dag_cross_check_lowering=False,
        compiled_dag_output_chunk_size=None,
        compiled_dag_helicity_filter_samples=10,
        compiled_dag_helicity_filter_seed=12345,
        compiled_dag_helicity_filter_relative_tolerance=1.0e-12,
        compiled_dag_helicity_filter_zero_tolerance=1.0e-300,
        compiled_dag_helicity_filter_phase_space="rambo",
    )
    if include_legacy_compiled_dag_options:
        parser.add_argument(
            "--compiled-dag-evaluator",
            action="store_true",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-lowering",
            choices=("spenso", "symbolic"),
            default="spenso",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-cross-check-lowering",
            action="store_true",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-output-chunk-size",
            type=int,
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-inline-external-wavefunctions",
            dest="compiled_dag_inline_external_wavefunctions",
            action="store_true",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--no-compiled-dag-inline-external-wavefunctions",
            dest="compiled_dag_inline_external_wavefunctions",
            action="store_false",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-helicity-filter",
            dest="compiled_dag_helicity_filter",
            action="store_true",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--no-compiled-dag-helicity-filter",
            dest="compiled_dag_helicity_filter",
            action="store_false",
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-helicity-filter-samples",
            type=int,
            default=10,
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-helicity-filter-seed",
            type=int,
            default=12345,
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-helicity-filter-relative-tolerance",
            type=float,
            default=1.0e-12,
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-helicity-filter-zero-tolerance",
            type=float,
            default=1.0e-300,
            help=argparse.SUPPRESS,
        )
        parser.add_argument(
            "--compiled-dag-helicity-filter-phase-space",
            choices=("rambo", "canonical"),
            default="rambo",
            help=argparse.SUPPRESS,
        )
    parser.add_argument(
        "--verbose-evaluator-build",
        action="store_true",
        help="Pass verbose evaluator-build progress through to supported backends.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=16,
        help="Number of phase-space samples to evaluate together at runtime.",
    )
    parser.add_argument(
        "--no-inlined-helicity-sum",
        action="store_true",
        help=(
            "Expose retained raw helicity amplitudes instead of only the final "
            "spin-summed matrix element."
        ),
    )
    parser.add_argument(
        "--symbolica-evaluator-backend",
        choices=("jit", "compiled-complex", "compiled-complex-4x"),
        default=argparse.SUPPRESS,
        help=(
            "Symbolica evaluator backend for generic DAG stage blocks. "
            "compiled-complex writes and compiles C++ complex evaluators; "
            "compiled-complex-4x uses Symbolica's SIMD complex backend."
        ),
    )
    parser.add_argument(
        "--symbolica-iterations",
        "--symbolica-n-iterations",
        dest="symbolica_iterations",
        type=int,
        default=10,
        help="Number of Horner-scheme optimization iterations.",
    )
    parser.add_argument(
        "--symbolica-cpe-iterations",
        "--symbolica-n-cpe-iterations",
        dest="symbolica_cpe_iterations",
        type=int,
        help=(
            "Number of common-pair elimination iterations for generic DAG "
            "stage evaluator optimization. Defaults to Symbolica's backend "
            "choice unless specified."
        ),
    )
    parser.add_argument(
        "--symbolica-n-cores",
        type=int,
        default=4,
        help="Number of cores used by Symbolica evaluator optimization.",
    )
    parser.add_argument(
        "--symbolica-direct-translation",
        dest="symbolica_direct_translation",
        action="store_true",
        help="Enable direct translation when Symbolica builds evaluator instructions.",
    )
    parser.add_argument(
        "--symbolica-no-direct-translation",
        dest="symbolica_direct_translation",
        action="store_false",
        help="Disable direct translation when Symbolica builds evaluator instructions.",
    )
    parser.add_argument(
        "--symbolica-jit-direct-translation",
        dest="symbolica_jit_direct_translation",
        action="store_true",
        help="Translate Symbolica instructions directly to SymJIT IR.",
    )
    parser.add_argument(
        "--symbolica-no-jit-direct-translation",
        dest="symbolica_jit_direct_translation",
        action="store_false",
        help="Disable direct translation to SymJIT IR.",
    )
    parser.add_argument(
        "--symbolica-jit-optimization-level",
        "--symbolica-jit-opt-level",
        dest="symbolica_jit_optimization_level",
        type=int,
        choices=(0, 1, 2, 3),
        default=3,
        help="SymJIT opt_level passed through Symbolica's JIT evaluator.",
    )
    parser.add_argument(
        "--symbolica-max-horner-scheme-variables",
        type=int,
        default=1000,
        help="Maximum number of variables considered in a Horner scheme.",
    )
    parser.add_argument(
        "--symbolica-max-common-pair-cache-entries",
        type=int,
        default=5_000_000,
        help="Maximum common-pair cache entries during Symbolica optimization.",
    )
    parser.add_argument(
        "--symbolica-max-common-pair-distance",
        type=int,
        help=(
            "Maximum distance between common pairs before cache eviction. "
            "Defaults to 1000 for generic DAG stage evaluators unless specified."
        ),
    )
    parser.add_argument(
        "--symbolica-collect-factors",
        action="store_true",
        help="Call collect_factors() on each block output before evaluator building.",
    )
    parser.add_argument(
        "--symbolica-compiled-preset",
        choices=(
            "manual",
            "adaptive",
            "generation",
            "balanced",
            "runtime",
            "runtime-o3",
        ),
        default="adaptive",
        help=(
            "Preset for compiled-complex code generation. adaptive uses the "
            "generic Rusticol stage-evaluator default, including output "
            "chunking when not set explicitly; manual respects the explicit "
            "inline-asm and optimization options; "
            "generation uses inline assembly; balanced uses generic C++ at "
            "-O1; runtime uses the measured generic C++ fast path; runtime-o3 "
            "uses generic C++ at -O3 with the runtime chunking policy."
        ),
    )
    parser.add_argument(
        "--symbolica-compiled-inline-asm",
        default="default",
        help="Inline assembly mode passed to Evaluator.compile for compiled-complex.",
    )
    parser.add_argument(
        "--symbolica-compiled-optimization-level",
        type=int,
        choices=(0, 1, 2, 3),
        default=3,
        help="C++ optimization level passed to Evaluator.compile.",
    )
    parser.add_argument(
        "--symbolica-compiled-native",
        dest="symbolica_compiled_native",
        action="store_true",
        help="Compile generated C++ for the native architecture.",
    )
    parser.add_argument(
        "--symbolica-no-compiled-native",
        dest="symbolica_compiled_native",
        action="store_false",
        help="Disable native-architecture flags in generated C++ compilation.",
    )
    parser.add_argument(
        "--symbolica-compiler-path",
        help="Optional compiler executable path for Evaluator.compile.",
    )
    parser.add_argument(
        "--symbolica-compiler-flag",
        dest="symbolica_compiler_flags",
        action="append",
        default=[],
        help=(
            "Additional compiler flag passed to Evaluator.compile. Repeat this "
            "option to pass multiple flags."
        ),
    )
    parser.add_argument(
        "--symbolica-output-chunk-size",
        "--symbolica-compiled-output-chunk-size",
        dest="symbolica_compiled_output_chunk_size",
        type=int,
        default=_DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
        help=(
            "Build Symbolica block outputs in chunks of this size. Applies to "
            "JIT and generated-code evaluators; the longer "
            "--symbolica-compiled-output-chunk-size spelling is kept as a "
            "compatibility alias. Defaults to 128."
        ),
    )
    parser.add_argument(
        "--symbolica-compiled-chunk-compile-workers",
        type=int,
        default=1,
        help=(
            "Number of worker threads used to compile chunked Symbolica "
            "generated-code evaluators."
        ),
    )
    parser.add_argument(
        "--symbolica-compiled-output-dir",
        type=Path,
        help=(
            "Directory where generated compiled Symbolica C++ sources and "
            "shared libraries are kept. By default they are written to a "
            "temporary directory."
        ),
    )
    parser.add_argument(
        "--save-evaluator-dir",
        type=Path,
        help=(
            "Save a reusable generic shared-current DAG process artifact to "
            "this directory after generation."
        ),
    )
    parser.add_argument(
        "--load-evaluator-dir",
        type=Path,
        help=(
            "Load a previously saved compiled shared-current DAG evaluator "
            "artifact from this directory instead of rebuilding it."
        ),
    )
    parser.add_argument(
        "--symbolica-raw-sum-final-stage",
        action="store_true",
        help=(
            "Compile the final helicity-weighted |amplitude|^2 sum as a "
            "single compiled Symbolica output. This is experimental and is "
            "off by default because it is not consistently faster."
        ),
    )


def _set_fast_rusticol_dag_defaults(parser: argparse.ArgumentParser) -> None:
    parser.set_defaults(
        runtime_backend="rusticol",
        symbolica_evaluator_backend="jit",
        symbolica_compiled_preset="runtime-o3",
        symbolica_n_cores=10,
        symbolica_compiled_chunk_compile_workers=10,
        batch_size=64,
        symbolica_stage_local_parameter_layout=(
            _DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT
        ),
        symbolica_compiled_output_chunk_size=_DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
    )


def _generation_build_kwargs(
    args: argparse.Namespace,
    runtime_backend: str,
    save_dir: Path,
) -> dict[str, Any]:
    build_kwargs = _runtime_evaluator_kwargs(args)
    if (
        build_kwargs.get("symbolica_load_evaluator_dir") is None
        and build_kwargs.get("symbolica_compiled_output_dir") is None
    ):
        build_kwargs["symbolica_compiled_output_dir"] = str(
            Path(save_dir).expanduser() / "compiled"
        )
    return build_kwargs


def _symbolica_settings_from_runtime_kwargs(
    values: dict[str, Any],
    *,
    process: str | None = None,
):
    from .symbolica_evaluator import (
        SymbolicaEvaluatorSettings,
        _resolve_compiled_preset,
    )

    (
        compiled_inline_asm,
        compiled_optimization_level,
        compiled_output_chunk_size,
    ) = _resolve_compiled_preset(
        str(values["symbolica_compiled_preset"]),
        gluon_count=(
            None if process is None else _external_gluon_count(process)
        ),
        inline_asm=str(values["symbolica_compiled_inline_asm"]),
        optimization_level=int(values["symbolica_compiled_optimization_level"]),
        output_chunk_size=values["symbolica_compiled_output_chunk_size"],
    )

    return SymbolicaEvaluatorSettings(
        backend=str(values["symbolica_evaluator_backend"]),
        iterations=int(values["symbolica_iterations"]),
        cpe_iterations=values["symbolica_cpe_iterations"],
        n_cores=int(values["symbolica_n_cores"]),
        direct_translation=bool(values["symbolica_direct_translation"]),
        jit_direct_translation=values["symbolica_jit_direct_translation"],
        jit_optimization_level=int(values["symbolica_jit_optimization_level"]),
        max_horner_scheme_variables=int(
            values["symbolica_max_horner_scheme_variables"]
        ),
        max_common_pair_cache_entries=int(
            values["symbolica_max_common_pair_cache_entries"]
        ),
        max_common_pair_distance=int(values["symbolica_max_common_pair_distance"]),
        collect_factors=bool(values["symbolica_collect_factors"]),
        compiled_preset=str(values["symbolica_compiled_preset"]),
        compiled_inline_asm=compiled_inline_asm,
        compiled_optimization_level=compiled_optimization_level,
        compiled_native=bool(values["symbolica_compiled_native"]),
        compiler_path=(
            None
            if values["symbolica_compiler_path"] is None
            else str(values["symbolica_compiler_path"])
        ),
        compiler_flags=tuple(str(flag) for flag in values["symbolica_compiler_flags"]),
        compiled_output_chunk_size=compiled_output_chunk_size,
        compiled_chunk_compile_workers=int(
            values["symbolica_compiled_chunk_compile_workers"]
        ),
        compiled_output_dir=(
            None
            if values["symbolica_compiled_output_dir"] is None
            else str(values["symbolica_compiled_output_dir"])
        ),
        raw_sum_final_stage=bool(values["symbolica_raw_sum_final_stage"]),
    )


def _external_gluon_count(process: str) -> int:
    from .process_ir import build_process_ir

    ir = build_process_ir(process)
    return sum(1 for pdg in ir.outgoing_pdgs if abs(int(pdg)) == 21)


def _runtime_backend(args: argparse.Namespace) -> RuntimeBackend:
    if bool(getattr(args, "compiled_dag_evaluator", False)):
        raise ValueError(
            "--compiled-dag-evaluator is deprecated; use generic DAG process "
            "artifacts with Rusticol instead."
        )
    value = getattr(args, "runtime_backend", "auto")
    if value not in _RUNTIME_BACKENDS:
        raise ValueError(f"unknown runtime backend: {value}")
    return cast(RuntimeBackend, value)


def _runtime_evaluator_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    runtime_backend = _runtime_backend(args)
    jit_direct_translation = getattr(args, "symbolica_jit_direct_translation", None)
    if jit_direct_translation is None:
        jit_direct_translation = False
    cpe_iterations = getattr(args, "symbolica_cpe_iterations", None)
    common_pair_distance = getattr(
        args,
        "symbolica_max_common_pair_distance",
        None,
    )
    if common_pair_distance is None:
        common_pair_distance = 1000
    return {
        "batch_size": int(getattr(args, "batch_size", 16)),
        "merge_evaluators_strategy": bool(
            getattr(args, "merge_evaluators_strategy", False)
        ),
        "split_vertex_current_stages": bool(
            getattr(args, "symbolica_split_vertex_current_stages", False)
        ),
        "stage_local_parameter_layout": bool(
            getattr(
                args,
                "symbolica_stage_local_parameter_layout",
                _DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT,
            )
        ),
        "verbose_evaluator_build": bool(
            getattr(args, "verbose_evaluator_build", False)
        ),
        "symbolica_evaluator_backend": str(
            getattr(args, "symbolica_evaluator_backend", "jit")
        ),
        "symbolica_iterations": int(getattr(args, "symbolica_iterations", 10)),
        "symbolica_cpe_iterations": cpe_iterations,
        "symbolica_n_cores": int(getattr(args, "symbolica_n_cores", 4)),
        "symbolica_direct_translation": bool(
            getattr(args, "symbolica_direct_translation", True)
        ),
        "symbolica_jit_direct_translation": jit_direct_translation,
        "symbolica_jit_optimization_level": int(
            getattr(args, "symbolica_jit_optimization_level", 3)
        ),
        "symbolica_max_horner_scheme_variables": int(
            getattr(args, "symbolica_max_horner_scheme_variables", 1000)
        ),
        "symbolica_max_common_pair_cache_entries": int(
            getattr(args, "symbolica_max_common_pair_cache_entries", 5_000_000)
        ),
        "symbolica_max_common_pair_distance": int(common_pair_distance),
        "symbolica_collect_factors": bool(
            getattr(args, "symbolica_collect_factors", False)
        ),
        "symbolica_compiled_preset": str(
            getattr(args, "symbolica_compiled_preset", "adaptive")
        ),
        "symbolica_compiled_inline_asm": str(
            getattr(args, "symbolica_compiled_inline_asm", "default")
        ),
        "symbolica_compiled_optimization_level": int(
            getattr(args, "symbolica_compiled_optimization_level", 3)
        ),
        "symbolica_compiled_native": bool(
            getattr(args, "symbolica_compiled_native", True)
        ),
        "symbolica_compiler_path": getattr(args, "symbolica_compiler_path", None),
        "symbolica_compiler_flags": tuple(
            getattr(args, "symbolica_compiler_flags", ())
        ),
        "symbolica_compiled_output_chunk_size": _compiled_output_chunk_size(args),
        "symbolica_compiled_chunk_compile_workers": int(
            getattr(args, "symbolica_compiled_chunk_compile_workers", 1)
        ),
        "symbolica_compiled_output_dir": getattr(
            args,
            "symbolica_compiled_output_dir",
            None,
        ),
        "symbolica_load_evaluator_dir": getattr(
            args,
            "load_evaluator_dir",
            None,
        ),
        "symbolica_raw_sum_final_stage": bool(
            getattr(args, "symbolica_raw_sum_final_stage", False)
        ),
        "compiled_dag_lowering": str(getattr(args, "compiled_dag_lowering", "spenso")),
        "compiled_dag_cross_check_lowering": bool(
            getattr(args, "compiled_dag_cross_check_lowering", False)
        ),
        "compiled_dag_inline_external_wavefunctions": bool(
            getattr(args, "compiled_dag_inline_external_wavefunctions", True)
        ),
        "compiled_dag_helicity_filter": bool(
            getattr(args, "compiled_dag_helicity_filter", True)
        ),
        "compiled_dag_helicity_filter_samples": int(
            getattr(args, "compiled_dag_helicity_filter_samples", 10)
        ),
        "compiled_dag_helicity_filter_seed": int(
            getattr(args, "compiled_dag_helicity_filter_seed", 12345)
        ),
        "compiled_dag_helicity_filter_relative_tolerance": float(
            getattr(args, "compiled_dag_helicity_filter_relative_tolerance", 1.0e-12)
        ),
        "compiled_dag_helicity_filter_zero_tolerance": float(
            getattr(args, "compiled_dag_helicity_filter_zero_tolerance", 1.0e-300)
        ),
        "compiled_dag_helicity_filter_phase_space": str(
            getattr(args, "compiled_dag_helicity_filter_phase_space", "rambo")
        ),
    }


def _compiled_output_chunk_size(args: argparse.Namespace) -> int | None:
    compiled_dag_chunk_size = getattr(args, "compiled_dag_output_chunk_size", None)
    if compiled_dag_chunk_size is not None:
        return int(compiled_dag_chunk_size)
    value = getattr(
        args,
        "symbolica_compiled_output_chunk_size",
        _DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
    )
    return None if value is None else int(value)


def _process_options(args: argparse.Namespace) -> ProcessOptions:
    return ProcessOptions(
        flavour_scheme=args.flavour_scheme,
        include_3qqbar=args.include_3qqbar,
        include_cc=args.include_cc,
        include_resonance=args.include_resonance,
        serial=not args.parallel_process_enumeration,
    )


def _max_quark_lines(args: argparse.Namespace) -> int | None:
    legacy = getattr(args, "max_quark_pairs", None)
    current = getattr(args, "max_quark_lines", None)
    if legacy is not None and current is not None and int(legacy) != int(current):
        raise ValueError("--max-quark-lines and --max-quark-pairs disagree")
    value = current if current is not None else legacy
    return None if value is None else int(value)


def _generic_dag_pruning_kwargs(
    args: argparse.Namespace,
    *,
    process: str | None = None,
) -> dict[str, Any]:
    color_accuracy = str(getattr(args, "color_accuracy", "lc")).lower()
    numerical_filter_current = getattr(args, "numerical_filter_current", None)
    if numerical_filter_current is None:
        numerical_filter_current = color_accuracy == "lc"
    numerical_current_merging = getattr(args, "numerical_current_merging", None)
    if numerical_current_merging is None:
        numerical_current_merging = color_accuracy == "lc"
    kwargs: dict[str, Any] = {
        "max_coupling_orders": _parse_max_coupling_orders(args),
        "max_lc_current_line_groups": getattr(
            args,
            "max_lc_current_line_groups",
            None,
        ),
        "max_quark_pairs": _max_quark_lines(args),
        "closure_side_mask_pruning": bool(
            getattr(args, "closure_side_mask_pruning", True)
        ),
        "color_order_mask_pruning": bool(
            getattr(args, "color_order_mask_pruning", True)
        ),
        "species_reachability_pruning": bool(
            getattr(args, "species_reachability_pruning", True)
        ),
        "ignored_particle_ids": _parse_ignored_particle_ids(
            str(getattr(args, "ignore_particles", ""))
        ),
        "ignored_vertex_kinds": _parse_int_list(
            str(getattr(args, "ignore_vertex_kinds", "")),
            option="--ignore-vertex-kinds",
        ),
        "numerical_filter_current": bool(numerical_filter_current),
        "numerical_current_merging": bool(numerical_current_merging),
        "numerical_current_samples": int(
            getattr(args, "numerical_current_samples", 10)
        ),
        "numerical_current_seed": int(
            getattr(args, "numerical_current_seed", 12345)
        ),
        "numerical_current_relative_tolerance": float(
            getattr(args, "numerical_current_relative_tolerance", 1.0e-12)
        ),
        "numerical_current_zero_tolerance": float(
            getattr(args, "numerical_current_zero_tolerance", 1.0e-300)
        ),
    }
    explicit_sector_ids = _parse_int_list(
        str(getattr(args, "lc_sector_ids", "")),
        option="--lc-sector-ids",
    )
    reference_color_order = _parse_int_list(
        str(getattr(args, "reference_color_order", "")),
        option="--reference-color-order",
    )
    if reference_color_order:
        kwargs["reference_color_order"] = tuple(reference_color_order)
    if explicit_sector_ids:
        kwargs["selected_color_sector_ids"] = set(explicit_sector_ids)
    elif (
        process is not None
        and str(getattr(args, "color_accuracy", "lc")) == "lc"
        and str(getattr(args, "lc_sector_strategy", "topology-representatives"))
        == "reference"
    ):
        kwargs["selected_color_sector_ids"] = {0}
    elif (
        process is not None
        and str(getattr(args, "color_accuracy", "lc")) == "lc"
        and str(getattr(args, "lc_sector_strategy", "topology-representatives"))
        == "line-pairing-representatives"
    ):
        representative_ids = _lc_line_pairing_representative_ids(
            process,
            args,
        )
        if representative_ids:
            kwargs["selected_color_sector_ids"] = representative_ids
    elif (
        process is not None
        and str(getattr(args, "color_accuracy", "lc")) == "lc"
        and str(getattr(args, "lc_sector_strategy", "topology-representatives"))
        == "topology-representatives"
    ):
        representative_ids = _lc_topology_representative_ids(
            process,
            args,
        )
        if representative_ids:
            kwargs["selected_color_sector_ids"] = representative_ids
    if (
        process is not None
        and str(getattr(args, "coupling_order_policy", "all")) == "minimal"
    ):
        from .generic_dag import infer_minimal_coupling_order_limits

        inferred_limits = infer_minimal_coupling_order_limits(
            process,
            color_accuracy=str(getattr(args, "color_accuracy", "lc")),
            options=_process_options(args),
            max_color_sectors=int(getattr(args, "max_color_sectors", 20000)),
            selected_color_sector_ids=kwargs.get("selected_color_sector_ids"),
            max_coupling_orders=kwargs["max_coupling_orders"],
            closure_side_mask_pruning=kwargs["closure_side_mask_pruning"],
            color_order_mask_pruning=kwargs["color_order_mask_pruning"],
            ignored_particle_ids=kwargs["ignored_particle_ids"],
            ignored_vertex_kinds=kwargs["ignored_vertex_kinds"],
        )
        kwargs["max_coupling_orders"] = _merge_coupling_order_limits(
            kwargs["max_coupling_orders"],
            inferred_limits,
        )
    return kwargs


def _compare_generic_dag_pruning_kwargs(
    args: argparse.Namespace,
    *,
    process: str,
) -> dict[str, Any]:
    """Generic DAG controls for AmpliCol comparison artifacts.

    Comparison defaults should stay tied to the LC sector selected by the
    Fortran process file.  Other generic pruning options still use the process
    string when needed, for example minimal coupling-order inference.
    """

    kwargs = _generic_dag_pruning_kwargs(args, process=process)
    explicit_sector_ids = _parse_int_list(
        str(getattr(args, "lc_sector_ids", "")),
        option="--lc-sector-ids",
    )
    strategy = str(getattr(args, "lc_sector_strategy", "topology-representatives"))
    if not explicit_sector_ids and strategy in {"topology-representatives", "reference"}:
        kwargs.pop("selected_color_sector_ids", None)
    return kwargs


def _normalize_generic_dag_pruning_kwargs(
    values: Mapping[str, Any] | None,
) -> dict[str, Any]:
    if values is None:
        return {}
    normalized = dict(values)
    selected_ids = normalized.get("selected_color_sector_ids")
    if selected_ids is not None:
        normalized["selected_color_sector_ids"] = {
            int(sector_id) for sector_id in selected_ids
        }
    coupling_orders = normalized.get("max_coupling_orders")
    if isinstance(coupling_orders, Mapping):
        normalized["max_coupling_orders"] = {
            str(name).upper(): int(value)
            for name, value in coupling_orders.items()
        }
    for key in ("ignored_particle_ids", "ignored_vertex_kinds"):
        if normalized.get(key) is not None:
            normalized[key] = [int(value) for value in normalized[key]]
    return normalized


def _merge_coupling_order_limits(
    explicit_limits: Mapping[str, int],
    inferred_limits: Mapping[str, int],
) -> dict[str, int]:
    merged = {str(name).upper(): int(value) for name, value in explicit_limits.items()}
    for name, value in inferred_limits.items():
        normalized = str(name).upper()
        if normalized in merged:
            merged[normalized] = min(merged[normalized], int(value))
        else:
            merged[normalized] = int(value)
    return merged


def _parse_max_coupling_orders(args: argparse.Namespace) -> dict[str, int]:
    limits: dict[str, int] = {}
    for item in getattr(args, "max_coupling_order", ()) or ():
        if "=" not in str(item):
            raise ValueError(
                "--max-coupling-order expects NAME=N, for example QED=1"
            )
        name, value = str(item).split("=", 1)
        limits[name.strip().upper()] = int(value)
    if getattr(args, "max_qcd_order", None) is not None:
        limits["QCD"] = int(args.max_qcd_order)
    if getattr(args, "max_qed_order", None) is not None:
        limits["QED"] = int(args.max_qed_order)
    return limits


def _parse_ignored_particle_ids(value: str) -> tuple[int, ...]:
    ids: list[int] = []
    for item in _split_cli_list(value):
        lowered = item.lower()
        if lowered in PDGS:
            ids.append(int(PDGS[lowered]))
        else:
            ids.append(int(item))
    return tuple(ids)


def _parse_int_list(value: str, *, option: str) -> tuple[int, ...]:
    try:
        return tuple(int(item) for item in _split_cli_list(value))
    except ValueError as exc:
        raise ValueError(f"{option} expects comma-separated integers") from exc


def _split_cli_list(value: str) -> tuple[str, ...]:
    return tuple(
        item.strip()
        for chunk in value.split(",")
        for item in chunk.split()
        if item.strip()
    )


def _lc_topology_representative_ids(
    process: str,
    args: argparse.Namespace,
) -> set[int]:
    from .color_plan import build_color_plan, lc_topology_replay_safe_groups

    plan = build_color_plan(
        process,
        color_accuracy=str(getattr(args, "color_accuracy", "lc")),
        options=_process_options(args),
        max_sectors=int(getattr(args, "max_color_sectors", 20000)),
    )
    return {
        int(group.representative_sector_id)
        for group in lc_topology_replay_safe_groups(plan)
    }


def _lc_line_pairing_representative_ids(
    process: str,
    args: argparse.Namespace,
) -> set[int]:
    from .color_plan import build_color_plan, lc_line_pairing_representative_ids

    plan = build_color_plan(
        process,
        color_accuracy=str(getattr(args, "color_accuracy", "lc")),
        options=_process_options(args),
        max_sectors=int(getattr(args, "max_color_sectors", 20000)),
    )
    return set(lc_line_pairing_representative_ids(plan))


def _cmd_inspect() -> int:
    import symbolica
    import symbolica.community.idenso  # noqa: F401
    import symbolica.community.spenso  # noqa: F401

    local_versions = getattr(symbolica, "LOCAL_VERSIONS", None)
    if not isinstance(local_versions, dict):
        raise RuntimeError("symbolica.LOCAL_VERSIONS is missing or is not a dict")

    expected = {"symbolica", "spenso", "idenso"}
    missing = expected.difference(local_versions)
    if missing:
        raise RuntimeError(f"symbolica.LOCAL_VERSIONS is missing keys: {sorted(missing)}")

    print("Loaded pyamplicol with Symbolica community modules: idenso, spenso")
    print("LOCAL_VERSIONS:")
    pprint(local_versions)
    return 0


def _cmd_processes(args: argparse.Namespace) -> int:
    try:
        process_set = enumerate_process_set(args.process, _process_options(args))
    except ValueError as exc:
        payload = {
            "available": False,
            "error": str(exc),
            "process": args.process,
            "runtime_backend": "rusticol",
        }
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            print(str(exc), file=sys.stderr)
        return 1
    if args.legacy_output is not None:
        if len(process_set.entries) != 1:
            raise ValueError(
                "--legacy-output is only supported for one concrete process entry"
            )
        enumerator = ProcessEnumerator(_process_options(args))
        enumerator.write_legacy_file(
            process_set.entries[0].enumeration,
            args.legacy_output,
        )

    payload = _process_set_enumeration_to_dict(process_set)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            f"{args.process}: {payload['n_entries']} process entries, "
            f"{payload['n_unique_processes']} unique processes, "
            f"{payload['n_groups']} phase-space groups, {payload['n_records']} records"
        )
        if args.legacy_output is not None:
            print(f"legacy export: {args.legacy_output}")
    return 0


def _cmd_process_plan(args: argparse.Namespace) -> int:
    from .generic_artifact import (
        build_generic_process_set_manifest,
        write_generic_process_manifest,
        write_generic_process_set_manifest,
    )

    manifest = build_generic_process_set_manifest(
        args.process,
        options=_process_options(args),
        color_accuracy=str(args.color_accuracy),
        max_currents=int(args.max_currents),
        max_color_sectors=int(args.max_color_sectors),
        **_generic_dag_pruning_kwargs(
            args,
            process=(
                args.process
                if not any(marker in args.process for marker in ("|", "[", "p", "j"))
                else None
            ),
        ),
    )
    output_dir = Path(args.output_dir).expanduser()
    if len(manifest.processes) == 1:
        subprocess_manifest = manifest.processes[0]
        subprocess_payload = subprocess_manifest.to_json_dict()
        manifest_path = write_generic_process_manifest(
            subprocess_manifest,
            output_dir,
        )
        payload = {
            "available": True,
            "kind": "pyamplicol-generic-process-plan",
            "manifest": str(manifest_path),
            "process": subprocess_manifest.process,
            "key": subprocess_manifest.key,
            "planning_status": subprocess_payload["planning_status"],
            "lowering_status": subprocess_payload["lowering_status"],
        }
    else:
        manifest_path = write_generic_process_set_manifest(manifest, output_dir)
        payload = {
            "available": True,
            "kind": "pyamplicol-generic-process-set-plan",
            "manifest": str(manifest_path),
            "request": manifest.request,
            "default_process_key": manifest.default_key,
            "generic_generation": manifest.generation_metadata,
            "n_processes": len(manifest.processes),
            "processes": [
                {
                    "key": subprocess_manifest.key,
                    "process": subprocess_manifest.process,
                    "manifest": str(
                        output_dir
                        / "subprocesses"
                        / subprocess_manifest.key
                        / "generic_process_manifest.json"
                    ),
                    "planning_status": subprocess_payload["planning_status"],
                    "lowering_status": subprocess_payload["lowering_status"],
                }
                for subprocess_manifest in manifest.processes
                for subprocess_payload in (subprocess_manifest.to_json_dict(),)
            ],
        }

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        rows = [
            DisplayRow({"metric": "Request", "value": args.process}, "bold"),
            {"metric": "Processes", "value": payload.get("n_processes", 1)},
            {"metric": "Colour", "value": args.color_accuracy},
            {"metric": "Manifest", "value": payload["manifest"]},
        ]
        if len(manifest.processes) == 1:
            status = payload["lowering_status"]
            if isinstance(status, dict):
                rows.extend(
                    [
                        {"metric": "Currents", "value": status["current_count"]},
                        {"metric": "Interactions", "value": status["interaction_count"]},
                        {"metric": "Closures", "value": status["closure_count"]},
                        {
                            "metric": "Tensor-ready",
                            "value": status["full_tensor_network_ready"],
                        },
                    ]
                )
        _display(args).print_table("Generic Process Plan", _kv_columns(), rows)
    return 0


def _cmd_generate(args: argparse.Namespace) -> int:
    return _legacy_native_command_unavailable(args, "generate")


def _cmd_generate_reference(args: argparse.Namespace) -> int:
    from .evaluation import NativeRuntimeEvaluator
    from .legacy_matrix import NativeMatrixElementGenerator

    generator = NativeMatrixElementGenerator(
        cache_dir=None if args.no_cache else args.cache_dir
    )
    result = generator.generate(
        args.process,
        options=_process_options(args),
        write_cache_metadata=not args.no_cache,
    )
    payload = result.to_json_dict()
    if bool(getattr(args, "compiled_dag_evaluator", False)) or getattr(
        args,
        "save_evaluator_dir",
        None,
    ) is not None:
        runtime = NativeRuntimeEvaluator(
            args.process,
            runtime_backend=_runtime_backend(args),
            allow_python_fallback=False,
            allow_reference_legacy=True,
            **_runtime_evaluator_kwargs(args),
        )
        payload["native_runtime_backend"] = runtime.metadata.to_json_dict()
        save_dir = getattr(args, "save_evaluator_dir", None)
        if save_dir is not None:
            manifest_path = runtime.save_evaluator_artifact(save_dir)
            payload["saved_evaluator_manifest"] = str(manifest_path)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        status = "supported" if result.supported_native_target else "unsupported"
        print(f"{args.process}: native target {status}; backend={result.backend}")
        runtime_payload = payload.get("native_runtime_backend")
        if isinstance(runtime_payload, dict):
            print(
                "runtime evaluator: "
                f"{runtime_payload.get('backend')} ({runtime_payload.get('kernel')})"
            )
        if "saved_evaluator_manifest" in payload:
            print(f"saved evaluator: {payload['saved_evaluator_manifest']}")
        if result.cache_file is not None:
            print(f"metadata cache: {result.cache_file}")
        for note in result.notes:
            print(note)
    return 0


def _cmd_evaluate(args: argparse.Namespace) -> int:
    return _legacy_native_command_unavailable(args, "evaluate")


def _cmd_evaluate_reference(args: argparse.Namespace) -> int:
    from .evaluation import NativeRuntimeEvaluator
    from .legacy_matrix import NativeMatrixElementGenerator

    artifact = _load_optional_evaluator_artifact(args.cache_dir, args.process)
    try:
        evaluator = NativeRuntimeEvaluator(
            args.process,
            runtime_backend=_runtime_backend(args),
            allow_reference_legacy=True,
            **_runtime_evaluator_kwargs(args),
        )
        evaluation_start = time.perf_counter()
        result = evaluator.evaluate(sqrt_s=args.sqrt_s)
        evaluation_time_s = time.perf_counter() - evaluation_start
    except NativeEvaluationError as exc:
        generator = NativeMatrixElementGenerator(cache_dir=args.cache_dir)
        generation = generator.generate(args.process, options=_process_options(args))
        payload = generation.to_json_dict()
        payload["evaluation_available"] = False
        payload["reason"] = str(exc)
        payload["evaluator_artifact"] = _artifact_status(
            args.cache_dir,
            args.process,
            artifact,
        )
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            print(payload["reason"])
        return 2

    runtime_metadata = evaluator.refresh_metadata().to_json_dict()
    payload = _evaluation_to_dict(result, runtime_metadata=runtime_metadata)
    payload["native_runtime_backend"] = runtime_metadata
    payload["native_evaluation_time_s"] = evaluation_time_s
    payload["batch_size"] = args.batch_size
    payload["merge_evaluators_strategy"] = args.merge_evaluators_strategy
    payload["verbose_evaluator_build"] = args.verbose_evaluator_build
    if args.no_inlined_helicity_sum:
        payload["raw_amplitude_vector"] = payload["helicity_contributions"]
    payload["evaluator_artifact"] = _artifact_status(
        args.cache_dir,
        args.process,
        artifact,
    )
    payload["symbolic_evaluation"] = _symbolic_artifact_evaluation(
        args.process,
        result,
        artifact,
    )
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"backend: {payload['backend']} ({payload['kernel']})")
        print(f"matrix element: {result.matrix_element:.16e}")
        print(f"raw helicity sum: {result.raw_helicity_sum:.16e}")
        artifact_payload = payload["evaluator_artifact"]
        if isinstance(artifact_payload, dict) and artifact_payload.get("available"):
            print(f"evaluator artifact: {artifact_payload['file']}")
        symbolic_payload = payload["symbolic_evaluation"]
        if isinstance(symbolic_payload, dict) and symbolic_payload.get("available"):
            print(
                "symbolica artifact matrix element: "
                f"{symbolic_payload['matrix_element']:.16e}"
            )
            print(
                "symbolica/native relative difference: "
                f"{symbolic_payload['relative_difference']:.6g}"
            )
    return 0


def _cmd_compare_amplicol(args: argparse.Namespace) -> int:
    options = _process_options(args)
    runtime_backend = _runtime_backend(args)
    runtime_evaluator_kwargs = _runtime_evaluator_kwargs(args)
    legacy_me_test = bool(getattr(args, "me_test", False))
    momenta_probe = bool(
        getattr(args, "amplicol_momenta_probe", False)
        or (not args.amplicol_probe and not legacy_me_test)
    )
    fixed_probe = (
        args.amplicol_probe
        and not momenta_probe
        and _should_use_fixed_amplicol_probe(args.process, options)
    )
    if args.dry_run:
        payload = _amplicol_dry_run_payload(
            args,
            options,
            fixed_probe=fixed_probe,
            momenta_probe=momenta_probe,
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    rusticol_process_dir: Path | None = None
    adapter = AmplicolAdapter(
        args.amplicol_root,
        jobs=args.jobs,
        timeout=args.timeout,
    )
    commands: list[CommandResult] = []
    if not args.skip_build:
        if args.amplicol_probe or momenta_probe:
            build = adapter.prepare_library(
                args.process,
                process_file=args.process_file,
                options=options,
            )
        elif legacy_me_test:
            build = adapter.prepare_direct_probe(
                args.process,
                process_file=args.process_file,
                options=options,
            )
        else:  # pragma: no cover - parser logic keeps one compare mode selected.
            raise RuntimeError("no AmpliCol comparison mode selected")
        commands.extend(build.commands)
        process_file = build.process_file
    else:
        process_file = args.process_file

    if momenta_probe:
        from .phase_space import generic_validation_point

        run = adapter.run_amplicol_momenta_probe(
            args.process,
            particles=generic_validation_point(args.process),
            points=args.points,
            process_file=process_file,
            options=options,
            timing_sample=args.timing,
            use_library=True,
        )
        mode = "amplicol_momenta_probe_library"
    elif args.amplicol_probe:
        if fixed_probe:
            run = adapter.run_amplicol_fixed_probe(
                args.process,
                points=args.points,
                process_file=process_file,
                options=options,
                timing_sample=args.timing,
                use_library=True,
            )
            mode = "amplicol_fixed_probe_library"
        else:
            run = adapter.run_amplicol_probe(
                args.process,
                points=args.points,
                process_file=process_file,
                options=options,
                timing_sample=args.timing,
                use_library=True,
            )
            mode = "amplicol_probe_library"
    elif legacy_me_test:
        run = adapter.run_me_test(
            args.process,
            points=args.points,
            process_file=process_file,
            options=options,
            mg5_path=args.mg5_path,
            timing_sample=args.timing,
        )
        mode = "me_test"
    else:  # pragma: no cover - parser logic keeps one compare mode selected.
        raise RuntimeError("no AmpliCol comparison mode selected")
    commands.extend(run.commands)
    if runtime_backend == "rusticol":
        native_generation, rusticol_process_dir = _rusticol_generation_profile(
            args.process,
            args.cache_dir,
            options,
            runtime_evaluator_kwargs=runtime_evaluator_kwargs,
            reference_color_order=_amplicol_reference_color_order_for_run(run),
            generic_dag_pruning_kwargs=_compare_generic_dag_pruning_kwargs(
                args,
                process=args.process,
            ),
        )
    else:
        native_generation = _native_generation_profile(
            args.process,
            args.cache_dir,
            options,
        )
    payload = {
        "mode": mode,
        "process": args.process,
        "process_file": str(run.process_file),
        "first_point_matrix_element": run.first_point_matrix_element,
        "first_phase_space_point": _first_point_to_dict(run.first_phase_space_point),
        "probe_points": [_probe_point_to_dict(point) for point in run.probe_points],
        "timing_rows": [
            {"label": row.label, "seconds": row.seconds, "note": row.note}
            for row in run.timing_rows
        ],
        "commands": [_command_summary(command) for command in commands],
        "total_command_time_s": sum(command.elapsed_s for command in commands),
        "pyamplicol_generation": native_generation,
    }
    if runtime_backend == "rusticol" and rusticol_process_dir is not None:
        _attach_rusticol_probe_comparison(
            run,
            payload,
            rusticol_process_dir,
            options=options,
            runtime_evaluator_kwargs=runtime_evaluator_kwargs,
            compare_args=args,
        )
    else:
        _attach_native_probe_comparison(
            args.process,
            run,
            payload,
            runtime_backend=runtime_backend,
            runtime_evaluator_kwargs=runtime_evaluator_kwargs,
        )
    payload["pyamplicol_performance"] = {
        "generation": native_generation,
        "runtime": payload.get("native_runtime"),
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        generation_time_value = native_generation.get("generation_time_s")
        generation_time_s = (
            generation_time_value
            if isinstance(generation_time_value, (float, int))
            else float("nan")
        )
        print(
            "pyamplicol generation: "
            f"{generation_time_s:.6g} s; "
            f"artifact cache hit={native_generation['artifact_cache_hit']}"
        )
        artifact_load_s = native_generation.get("artifact_load_s")
        if isinstance(artifact_load_s, (float, int)):
            print(f"pyamplicol artifact load: {artifact_load_s:.6g} s")
        print(f"process file: {payload['process_file']}")
        print(f"first-point AmpliCol ME: {payload['first_point_matrix_element']}")
        if payload.get("native_first_point_available"):
            print(f"first-point pyamplicol ME: {payload['native_first_point_matrix_element']}")
            if "native_first_point_relative_difference" in payload:
                print(
                    "first-point relative difference: "
                    f"{payload['native_first_point_relative_difference']:.6g}"
                )
        summary = payload.get("native_probe_summary")
        if isinstance(summary, dict):
            print(
                "native probe points: "
                f"{summary['available_points']}/{summary['reference_points']}; "
                f"max rel diff={summary.get('max_relative_difference')}"
            )
        runtime = payload.get("native_runtime")
        if isinstance(runtime, dict):
            backend = payload.get("native_runtime_backend")
            if isinstance(backend, dict):
                print(
                    "pyamplicol native backend: "
                    f"{backend.get('backend')} ({backend.get('kernel')})"
                )
            print(
                "pyamplicol native runtime: "
                f"total={runtime['total_s']:.6g} s; "
                f"mean={runtime['mean_per_point_s']:.6g} s/point"
            )
        print(f"total command time: {payload['total_command_time_s']:.6g} s")
        for row in run.timing_rows:
            print(f"{row.label}: {row.seconds:.6g} s {row.note}".rstrip())
    return 0


def _cmd_validate_z_gluon_family(args: argparse.Namespace) -> int:
    return _legacy_z_gluon_command_unavailable(args, "validate-z-gluon-family")


def _legacy_z_gluon_command_unavailable(
    args: argparse.Namespace,
    command: str,
) -> int:
    message = (
        f"{command} is a legacy Z+gluon-only command. Production generation, "
        "validation, and profiling now go through generic DAG process "
        "artifacts: use `process-plan`, `generate-process`, "
        "`time-process`, and `compare-amplicol --runtime-backend rusticol`."
    )
    if getattr(args, "json", False):
        print(
            json.dumps(
                {
                    "available": False,
                    "command": command,
                    "error": message,
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        print(message, file=sys.stderr)
    return 1


def _legacy_native_command_unavailable(
    args: argparse.Namespace,
    command: str,
) -> int:
    message = (
        f"{command} is a legacy native-kernel command. Production pyAmpliCol "
        "now uses generic DAG process artifacts: run `process-plan` to inspect "
        "support, `generate-process PROCESS OUTPUT_DIR` to build an artifact, "
        "`time-process OUTPUT_DIR` to evaluate/profile it through Rusticol, "
        "or `compare-amplicol --runtime-backend rusticol` for Fortran "
        "validation."
    )
    payload = {
        "available": False,
        "command": command,
        "error": message,
        "process": getattr(args, "process", None),
        "runtime_backend": _runtime_backend(args) if hasattr(args, "runtime_backend") else None,
    }
    if getattr(args, "json", False):
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(message, file=sys.stderr)
    return 1


def _cmd_validate_z_gluon_family_reference(args: argparse.Namespace) -> int:
    if args.min_gluons < 0:
        raise ValueError("--min-gluons must be non-negative")
    if args.max_gluons < args.min_gluons:
        raise ValueError("--max-gluons must be greater than or equal to --min-gluons")
    options = _process_options(args)
    if args.dry_run:
        payload = _z_gluon_family_dry_run_payload(args, options)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    adapter = AmplicolAdapter(
        args.amplicol_root,
        jobs=args.jobs,
        timeout=args.timeout,
    )
    rows: list[dict[str, object]] = []
    for gluon_count in range(args.min_gluons, args.max_gluons + 1):
        process = _z_gluon_family_process(gluon_count)
        native_generation = _native_generation_profile(
            process,
            args.cache_dir,
            options,
        )
        commands: list[CommandResult] = []
        build = adapter.prepare_library(process, options=options)
        if gluon_count == 0:
            run = adapter.run_amplicol_fixed_probe(
                process,
                points=args.points,
                process_file=build.process_file,
                options=options,
                timing_sample=args.timing,
                use_library=True,
            )
            mode = "amplicol_fixed_probe_library"
        else:
            run = adapter.run_amplicol_probe(
                process,
                points=args.points,
                process_file=build.process_file,
                options=options,
                timing_sample=args.timing,
                use_library=True,
            )
            mode = "amplicol_probe_library"
        timing_run = adapter.run_library_use(
            process,
            nevents=args.points if args.timing is None else args.timing,
            seed=101,
            options=options,
            timing_sample=1,
        )
        commands.extend(build.commands)
        commands.extend(run.commands)
        commands.extend(timing_run.commands)
        row: dict[str, object] = {
            "gluon_count": gluon_count,
            "process": process,
            "mode": mode,
            "fortran_timing_workflow": "generated_library_use",
            "runtime_backend": _z_gluon_family_runtime_backend(
                gluon_count,
                _runtime_backend(args),
            ),
            "process_file": str(run.process_file),
            "probe_points": [_probe_point_to_dict(point) for point in run.probe_points],
            "timing_rows": [
                {"label": timing.label, "seconds": timing.seconds, "note": timing.note}
                for timing in timing_run.timing_rows
            ],
            "commands": [_command_summary(command) for command in commands],
            "total_command_time_s": sum(command.elapsed_s for command in commands),
            "pyamplicol_generation": native_generation,
        }
        _attach_native_probe_comparison(
            process,
            run,
            row,
            runtime_backend=_z_gluon_family_runtime_backend(
                gluon_count,
                _runtime_backend(args),
            ),
            runtime_evaluator_kwargs=_runtime_evaluator_kwargs(args),
        )
        rows.append(row)

    summary = _z_gluon_family_summary(rows)
    passed = _z_gluon_family_passed(summary, rel_tol=args.rel_tol)
    payload = {
        "family": "d d~ -> Z + n g",
        "min_gluons": args.min_gluons,
        "max_gluons": args.max_gluons,
        "points": args.points,
        "rel_tol": args.rel_tol,
        "passed": passed,
        "runtime_backend": args.runtime_backend,
        "summary": summary,
        "rows": rows,
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "d d~ -> Z + n g validation: "
            f"n={args.min_gluons}..{args.max_gluons}, points={args.points}"
        )
        print(
            "validated rows: "
            f"{summary['validated_rows']}/{summary['requested_rows']}; "
            f"max rel diff={summary['max_relative_difference']}; "
            f"rel_tol={args.rel_tol}; passed={passed}"
        )
        for row in rows:
            probe_summary = row.get("native_probe_summary")
            max_rel = (
                probe_summary.get("max_relative_difference")
                if isinstance(probe_summary, dict)
                else None
            )
            runtime = row.get("native_runtime")
            mean_runtime = (
                runtime.get("mean_per_point_s")
                if isinstance(runtime, dict)
                else None
            )
            probe_points = row.get("probe_points")
            probe_point_count = len(probe_points) if isinstance(probe_points, list) else 0
            print(
                f"n={row['gluon_count']}: "
                f"points={probe_point_count}, "
                f"max_rel={max_rel}, "
                f"pyamplicol_mean_s={mean_runtime}"
            )
    return 0 if passed else 1


def _cmd_benchmark_z_gluon_modes(args: argparse.Namespace) -> int:
    return _legacy_z_gluon_command_unavailable(args, "benchmark-z-gluon-modes")


def _cmd_benchmark_z_gluon_modes_reference(args: argparse.Namespace) -> int:
    from .benchmarks import benchmark_z_gluon_modes, format_mode_benchmark_table

    evaluator_build_kwargs = _runtime_evaluator_kwargs(args)
    evaluator_build_kwargs["runtime_backend"] = "dag"
    payload = benchmark_z_gluon_modes(
        min_gluons=args.min_gluons,
        max_gluons=args.max_gluons,
        points=args.points,
        timing=args.timing,
        amplicol_root=args.amplicol_root,
        jobs=args.jobs,
        timeout=args.timeout,
        options=_process_options(args),
        numeric_timeout=args.numeric_timeout,
        parametric_timeout=args.parametric_timeout,
        parametric_max_gluons=args.parametric_max_gluons,
        tensor_strategy=args.tensor_strategy,
        output=args.output,
        batch_size=args.batch_size,
        merge_evaluators_strategy=args.merge_evaluators_strategy,
        verbose_evaluator_build=args.verbose_evaluator_build,
        evaluator_build_kwargs=evaluator_build_kwargs,
        include_python=not args.only_legacy_shared,
        include_numeric_tn=not args.only_legacy_shared,
        include_parametric_tn=not args.only_legacy_shared,
        include_shared_dag=True,
    )
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "d d~ -> Z + n g mode benchmark: "
            f"n={args.min_gluons}..{args.max_gluons}, points={args.points}, "
            f"tensor_strategy={args.tensor_strategy}, batch_size={args.batch_size}"
        )
        print(format_mode_benchmark_table(payload))
        summary = payload["summary"]
        print(
            "successful rows: "
            f"A={summary['mode_success_counts']['legacy']}, "
            f"B={summary['mode_success_counts']['python']}, "
            f"C={summary['mode_success_counts']['numeric_tn']}, "
            f"D-scalar={summary['mode_success_counts']['parametric_tn']}, "
            f"D-shared={summary['mode_success_counts']['shared_dag']}"
        )
        print(
            "max rel diff to legacy among available pyamplicol modes: "
            f"{summary['max_relative_difference_to_legacy']}"
        )
        print(
            "all four modes match for all rows: "
            f"{summary['all_four_modes_match_for_all_rows']}"
        )
    return 0


def _cmd_profile_tensor_evaluator(args: argparse.Namespace) -> int:
    return _legacy_z_gluon_command_unavailable(args, "profile-tensor-evaluator")


def _cmd_profile_tensor_evaluator_reference(args: argparse.Namespace) -> int:
    from .benchmarks import profile_z_gluon_tensor_evaluator

    payload = profile_z_gluon_tensor_evaluator(
        args.process,
        sqrt_s=args.sqrt_s,
        repetitions=args.repetitions,
        evaluator_repetitions=args.evaluator_repetitions,
        tensor_strategy=args.tensor_strategy,
    )
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"process: {payload['process']}")
        print(f"tensor strategy: {payload['tensor_strategy']}")
        print(f"gluon count: {payload['gluon_count']}")
        print(f"generation: {payload['generation_s']:.6g} s")
        print(
            "pure evaluator call: "
            f"{payload['pure_evaluator_call_us']:.6g} us"
        )
        print(f"full evaluate: {payload['full_evaluate_ms']:.6g} ms")
        print(
            "full evaluate per helicity: "
            f"{payload['full_evaluate_per_helicity_us']:.6g} us"
        )
    return 0


def _cmd_profile_dag_evaluator(args: argparse.Namespace) -> int:
    runtime_backend = _runtime_backend(args)
    if bool(getattr(args, "generate_only", False)):
        return _cmd_generate_dag_evaluator_artifact(args, runtime_backend)
    message = (
        "profile-dag-evaluator is a legacy Z+gluon profiler. Production "
        "profiling now uses generic DAG process artifacts: run "
        "`generate-process PROCESS OUTPUT_DIR` and then `time-process OUTPUT_DIR`."
    )
    if args.json:
        print(
            json.dumps(
                {
                    "available": False,
                    "error": message,
                    "process": args.process,
                    "runtime_backend": runtime_backend,
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        print(message, file=sys.stderr)
    return 1


def _cmd_generate_process(args: argparse.Namespace) -> int:
    try:
        max_quark_lines = _max_quark_lines(args)
        display = _display(args)
        with display.stage_progress(
            "Enumerating processes",
            total=4,
            metadata=str(args.process)[:48],
        ) as progress:
            progress.update(stage="parse", item=str(args.process)[:28])
            progress.update(stage="expand", item="variants", increment=1)
            report = build_generic_process_selection_report(
                args.process,
                _process_options(args),
                max_quark_pairs=max_quark_lines,
                use_prefilter=bool(getattr(args, "enumeration_prefilter", True)),
            )
            progress.update(
                stage="prefilter",
                item=f"cand={report.candidate_count} eval={report.evaluated_count}",
                increment=2,
            )
            progress.update(
                stage="canonicalize",
                item=f"sel={report.selected_count} rej={report.rejected_count}",
                increment=1,
            )
    except ValueError as exc:
        payload = {
            "available": False,
            "error": str(exc),
            "process": args.process,
            "runtime_backend": "rusticol",
        }
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            print(str(exc), file=sys.stderr)
        return 1
    if bool(getattr(args, "dry_run", False)):
        return _cmd_generate_process_dry_run(args, report)
    if getattr(args, "output_dir", None) is None:
        print("generate-process requires OUTPUT_DIR unless --dry-run is set", file=sys.stderr)
        return 2
    if not report.entries:
        message = f"no valid processes found for {args.process!r}"
        if args.json:
            print(json.dumps({"available": False, "error": message}, indent=2))
        else:
            print(message, file=sys.stderr)
        return 1
    process_set = ProcessSetEnumeration(
        request=args.process,
        options=_process_options(args),
        entries=report.entries,
        selection_report=report,
    )
    if len(process_set.entries) == 1 and not any(
        marker in args.process for marker in ("|", "[")
    ):
        args.process = process_set.entries[0].process
        return _cmd_generate_generic_dag_artifact(args, process_set.entries[0])

    return _cmd_generate_process_set(args, process_set)


def _cmd_generate_process_dry_run(
    args: argparse.Namespace,
    report: ProcessSelectionReport,
) -> int:
    payload = _process_selection_report_to_dict(report)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    display = _display(args)
    display.print_table(
        "Generic Process Enumeration Dry Run",
        [
            DisplayColumn("id", "ID", "right"),
            DisplayColumn("status", "Status"),
            DisplayColumn("process", "Process"),
            DisplayColumn("key", "Key"),
            DisplayColumn("quark_lines", "Quark lines", "right"),
            DisplayColumn("charge3", "Q3", "right"),
            DisplayColumn("source", "Source"),
        ],
        [
            DisplayRow(
                {
                    "id": index,
                    "status": record.status,
                    "process": record.process,
                    "key": record.key,
                    "quark_lines": record.quark_lines,
                    "charge3": record.charge3,
                    "source": record.source,
                },
                "green" if record.status == "selected" else "yellow",
            )
            for index, record in enumerate(
                (*report.selected_records, *report.duplicate_records),
                start=1,
            )
        ],
    )
    display.print_table(
        "Enumeration Summary",
        _kv_columns(),
        [
            DisplayRow({"metric": "Request", "value": report.request}, "bold"),
            {"metric": "Prefilter", "value": "enabled" if report.prefilter_enabled else "disabled"},
            {"metric": "Candidates", "value": report.candidate_count},
            {"metric": "Evaluated after cheap filters", "value": report.evaluated_count},
            {"metric": "Selected", "value": report.selected_count},
            {"metric": "Duplicates", "value": report.duplicate_count},
            {"metric": "Rejected", "value": report.rejected_count},
            {"metric": "Elapsed", "value": format_measurement(report.elapsed_s, unit="s")},
        ],
    )
    if report.rejection_counts:
        display.print_table(
            "Rejection Reasons",
            [
                DisplayColumn("reason", "Reason"),
                DisplayColumn("count", "Count", "right"),
            ],
            [
                {"reason": reason, "count": count}
                for reason, count in report.rejection_counts
                if count
            ],
        )
    return 0


def _process_selection_report_to_dict(
    report: ProcessSelectionReport,
) -> dict[str, object]:
    return {
        "available": bool(report.entries),
        "kind": "pyamplicol-generic-process-selection-report",
        "request": report.request,
        "prefilter_enabled": report.prefilter_enabled,
        "elapsed_s": report.elapsed_s,
        "candidate_count": report.candidate_count,
        "evaluated_count": report.evaluated_count,
        "selected_count": report.selected_count,
        "duplicate_count": report.duplicate_count,
        "rejected_count": report.rejected_count,
        "rejection_counts": {
            reason: count for reason, count in report.rejection_counts if count
        },
        "stage_timings": {
            stage: elapsed for stage, elapsed in report.stage_timings
        },
        "selected": [
            {
                "id": index,
                "key": record.key,
                "process": record.process,
                "source": record.source,
                "quark_lines": record.quark_lines,
                "charge3": record.charge3,
                "family": record.family,
            }
            for index, record in enumerate(report.selected_records, start=1)
        ],
        "duplicates": [
            {
                "key": record.key,
                "process": record.process,
                "source": record.source,
                "reason": record.reason,
            }
            for record in report.duplicate_records
        ],
    }


def _rusticol_artifact_unavailable_message(
    process: str,
    *,
    color_accuracy: str = "lc",
) -> str | None:
    """Return the current Rusticol artifact support diagnostic for a process."""

    report = _rusticol_artifact_support_report(
        process,
        color_accuracy=color_accuracy,
    )
    if bool(report.runtime_artifact_supported):
        return None
    return report.artifact_unavailable_message


def _rusticol_artifact_support_report(
    process: str,
    *,
    color_accuracy: str = "lc",
):
    """Return the current Rusticol artifact support report for a process."""

    from .process_support import classify_process_support

    return classify_process_support(process, color_accuracy=color_accuracy)


def _cmd_generate_generic_dag_artifact(
    args: argparse.Namespace,
    entry: ProcessSetEntry,
) -> int:
    from .generic_artifact import write_generic_dag_process_artifact

    output_dir = Path(args.output_dir).expanduser()
    if bool(getattr(args, "replace", False)) and output_dir.exists():
        shutil.rmtree(output_dir)
    generation_start = time.perf_counter()
    build_kwargs = _generation_build_kwargs(args, "rusticol", output_dir)
    symbolica_settings = _symbolica_settings_from_runtime_kwargs(
        build_kwargs,
        process=entry.process,
    )
    display = _display(args)
    with display.stage_progress(
        "Generating process",
        total=1,
        metadata=entry.key,
    ) as progress:
        progress_callback = _combined_progress_callback(
            progress.callback,
            _child_generation_progress_callback(entry.process),
        )
        manifest_path, manifest = write_generic_dag_process_artifact(
            entry.process,
            output_dir,
            options=_process_options(args),
            color_accuracy=str(getattr(args, "color_accuracy", "lc")),
            max_currents=int(getattr(args, "max_currents", 50000)),
            max_color_sectors=int(getattr(args, "max_color_sectors", 20000)),
            evaluator_backend=str(build_kwargs["symbolica_evaluator_backend"]),
            compiled_preset=str(build_kwargs["symbolica_compiled_preset"]),
            batch_size=int(build_kwargs["batch_size"]),
            emit_stage_evaluator_artifacts=True,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=bool(build_kwargs["merge_evaluators_strategy"]),
            stage_local_parameter_layout=bool(
                build_kwargs["stage_local_parameter_layout"]
            ),
            verbose_evaluator_build=bool(build_kwargs["verbose_evaluator_build"]),
            progress_callback=progress_callback,
            lc_topology_replay=bool(getattr(args, "lc_topology_replay", False)),
            **_generic_dag_pruning_kwargs(args, process=entry.process),
        )
        progress.update(stage="done", item=entry.key, increment=1)
    compiled = cast(dict[str, Any], manifest["compiled"])
    generation_s = time.perf_counter() - generation_start
    runtime_available = bool(compiled.get("runtime_available"))
    payload = {
        "available": True,
        "runtime_available": runtime_available,
        "runtime_backend": "rusticol",
        "kind": manifest["kind"],
        "artifact_class": manifest.get("artifact_class", "generic-dag-schema-v2"),
        "generation_s": generation_s,
        "process": entry.process,
        "key": entry.key,
        "saved_evaluator_manifest": str(manifest_path),
        "manifest": str(manifest_path),
        "planning_status": manifest["planning_status"],
        "lowering_status": manifest["lowering_status"],
        "runtime_unavailable_message": (
            None
            if runtime_available
            else compiled["runtime_unavailable_message"]
        ),
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        _display(args).print_table(
            "Generated Generic DAG Process Artifact",
            _kv_columns(),
            [
                DisplayRow({"metric": "Process", "value": entry.process}, "bold"),
                {
                    "metric": "Runtime",
                    "value": "rusticol schema-v2" if runtime_available else "pending schema-v2",
                },
                {"metric": "Manifest", "value": manifest_path},
                {
                    "metric": "Generation",
                    "value": format_measurement(generation_s, unit="s"),
                },
            ],
        )
    return 0


def _cmd_generate_process_set(args: argparse.Namespace, process_set: object) -> int:
    root = Path(args.output_dir).expanduser()
    manifest_path = root / "process_set_manifest.json"
    root.mkdir(parents=True, exist_ok=True)
    existing = _load_process_set_manifest(root)
    if existing is not None and not (args.append or args.replace):
        message = (
            f"process-set output already exists at {root}; use --append or "
            "--replace"
        )
        if args.json:
            print(json.dumps({"available": False, "error": message}, indent=2))
        else:
            print(message, file=sys.stderr)
        return 1

    existing_entries = {
        str(entry["key"]): dict(entry)
        for entry in _process_set_manifest_entries(existing or {})
        if isinstance(entry, dict) and "key" in entry
    }
    generated: list[dict[str, object]] = []
    skipped: list[str] = []
    subprocess_root = root / "subprocesses"
    subprocess_root.mkdir(parents=True, exist_ok=True)
    generation_metadata = _generic_process_set_generation_metadata(args)
    entries = cast(list[ProcessSetEntry], list(getattr(process_set, "entries")))
    work_items: list[tuple[ProcessSetEntry, Path]] = []
    representative_by_signature: dict[tuple[int, ...], ProcessSetEntry] = {}
    representative_for_key: dict[str, ProcessSetEntry] = {}
    crossing_signature_by_key: dict[str, tuple[int, ...]] = {}
    for existing_entry in existing_entries.values():
        if existing_entry.get("crossing_alias_of") is not None:
            continue
        existing_process = existing_entry.get("process")
        if isinstance(existing_process, str):
            try:
                signature = _process_crossing_reuse_signature(
                    existing_process,
                    color_accuracy=str(getattr(args, "color_accuracy", "lc")),
                    options=_process_options(args),
                )
            except ValueError:
                continue
            representative_by_signature.setdefault(
                signature,
                ProcessSetEntry(
                    key=str(existing_entry["key"]),
                    process=existing_process,
                    enumeration=entries[0].enumeration,
                ),
            )
    for entry in entries:
        signature = _process_crossing_reuse_signature(
            entry.process,
            color_accuracy=str(getattr(args, "color_accuracy", "lc")),
            options=_process_options(args),
        )
        representative = representative_by_signature.setdefault(signature, entry)
        representative_for_key[entry.key] = representative
        crossing_signature_by_key[entry.key] = signature
        if entry.key in existing_entries and args.append and not args.replace:
            skipped.append(entry.key)
            generated.append(existing_entries[entry.key])
            continue
        if representative.key != entry.key:
            continue
        subdir = subprocess_root / entry.key
        if subdir.exists() and not args.replace and entry.key in existing_entries:
            message = (
                f"process entry {entry.key!r} already exists; use --replace to rebuild it"
            )
            if args.json:
                print(json.dumps({"available": False, "error": message}, indent=2))
            else:
                print(message, file=sys.stderr)
            return 1
        if args.replace and subdir.exists():
            shutil.rmtree(subdir)
        work_items.append((entry, subdir))

    worker_count = max(int(getattr(args, "n_cores", 1)), 1)
    active: dict[int, tuple[subprocess.Popen[str], Any, Path, float]] = {}
    stderr_buffers: dict[int, list[str]] = {}
    stderr_threads: dict[int, threading.Thread] = {}
    child_progress_events: queue.SimpleQueue[tuple[int, dict[str, object]]] = (
        queue.SimpleQueue()
    )
    child_progress_seen: dict[int, tuple[str, str, float]] = {}
    failures: list[dict[str, object]] = []
    display = _display(args)
    try:
        with display.stage_progress(
            "Generating process set",
            total=len(entries),
            metadata=_process_set_progress_metadata(process_set, worker_count),
        ) as progress:
            if skipped:
                progress.update(
                    stage="skip",
                    item=f"{len(skipped)} existing",
                    increment=len(skipped),
                    ram=_rss_text_for_active_processes(active),
                )
            pending = list(work_items)
            while pending or active:
                while pending and len(active) < worker_count:
                    entry, subdir = pending.pop(0)
                    cmd = _generate_process_child_command(args, entry.process, subdir)
                    env = _child_generation_environment()
                    process = subprocess.Popen(
                        cmd,
                        cwd=Path.cwd(),
                        env=env,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        text=True,
                        start_new_session=True,
                    )
                    active[process.pid] = (process, entry, subdir, time.perf_counter())
                    stderr_buffers[process.pid] = []
                    stderr_threads[process.pid] = _start_child_stderr_reader(
                        process,
                        stderr_buffers[process.pid],
                        child_progress_events,
                    )
                    progress.update(
                        stage=f"run {len(active):02d}/{worker_count:02d}",
                        item=str(entry.key),
                        ram=_rss_text_for_active_processes(active),
                    )

                _drain_child_progress_events(
                    child_progress_events,
                    active,
                    progress,
                    child_progress_seen,
                )
                finished: list[int] = []
                for pid, (process, entry, subdir, started) in list(active.items()):
                    if process.poll() is None:
                        continue
                    stdout, stderr = _finish_generation_child_output(
                        process,
                        pid=pid,
                        stderr_buffers=stderr_buffers,
                        stderr_threads=stderr_threads,
                    )
                    finished.append(pid)
                    elapsed_s = time.perf_counter() - started
                    if process.returncode != 0:
                        child_payload = _parse_child_generation_payload(stdout)
                        child_error = str(
                            child_payload.get("error")
                            or stderr.strip()
                            or stdout.strip()
                            or f"subprocess exited with {process.returncode}"
                        )
                        failures.append(
                            {
                                "key": entry.key,
                                "process": entry.process,
                                "returncode": process.returncode,
                                "error": child_error,
                                "payload": child_payload,
                                "stdout": stdout,
                                "stderr": stderr,
                                "elapsed_s": elapsed_s,
                            }
                        )
                        progress.update(
                            stage="failed",
                            item=str(entry.key),
                            increment=1,
                            ram=_rss_text_for_active_processes(active),
                        )
                        continue
                    child_payload = _parse_child_generation_payload(stdout)
                    generated.append(
                        {
                            "key": entry.key,
                            "process": entry.process,
                            "path": str(Path("subprocesses") / entry.key),
                            "kind": str(
                                child_payload.get(
                                    "kind",
                                    "pyamplicol-generic-dag-process",
                                )
                            ),
                            "artifact_class": str(
                                child_payload.get(
                                    "artifact_class",
                                    "generic-dag-schema-v2",
                                )
                            ),
                            "generation_s": child_payload.get("generation_s", elapsed_s),
                            "runtime_available": bool(
                                child_payload.get("runtime_available", False)
                            ),
                            "generation_request": generation_metadata,
                        }
                    )
                    progress.update(
                        stage="done",
                        item=str(entry.key),
                        increment=1,
                        ram=_rss_text_for_active_processes(active),
                    )
                for pid in finished:
                    active.pop(pid, None)
                    stderr_buffers.pop(pid, None)
                    stderr_threads.pop(pid, None)
                _drain_child_progress_events(
                    child_progress_events,
                    active,
                    progress,
                    child_progress_seen,
                )
                if pending or active:
                    if not finished:
                        progress.update(
                            stage=f"run {len(active):02d}/{worker_count:02d}",
                            item=f"pending={len(pending):03d}",
                            ram=_rss_text_for_active_processes(active),
                        )
                        time.sleep(0.25)
            if failures:
                progress.update(
                    stage="failed",
                    item=f"{len(failures)} failures",
                    ram=_rss_text_for_active_processes(active),
                )
    except KeyboardInterrupt:
        _terminate_generation_children(active)
        _join_child_stderr_threads(stderr_threads)
        message = "process-set generation interrupted; active subprocesses were terminated"
        if args.json:
            print(json.dumps({"available": False, "error": message}, indent=2))
        else:
            print(message, file=sys.stderr)
        return 130
    except BaseException:
        _terminate_generation_children(active)
        _join_child_stderr_threads(stderr_threads)
        raise

    if failures:
        if args.json:
            print(
                json.dumps(
                    {"available": False, "failures": failures},
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            first = failures[0]
            print(
                f"failed to generate {first['key']}: {first.get('error')}",
                file=sys.stderr,
            )
        return 1

    representative_entries: dict[str, dict[str, object]] = dict(existing_entries)
    for item in generated:
        representative_entries[str(item["key"])] = item
    alias_generated: list[dict[str, object]] = []
    for entry in entries:
        if entry.key in existing_entries and args.append and not args.replace:
            continue
        alias_representative: ProcessSetEntry | None = representative_for_key.get(
            entry.key
        )
        if alias_representative is None or alias_representative.key == entry.key:
            continue
        representative_item = representative_entries.get(alias_representative.key)
        if representative_item is None:
            failures.append(
                {
                    "key": entry.key,
                    "process": entry.process,
                    "error": (
                        "crossing representative "
                        f"{alias_representative.key!r} was not generated"
                    ),
                }
            )
            continue
        alias_item = {
            "key": entry.key,
            "process": entry.process,
            "path": representative_item["path"],
            "kind": representative_item.get(
                "kind",
                "pyamplicol-generic-dag-process",
            ),
            "artifact_class": representative_item.get(
                "artifact_class",
                "generic-dag-schema-v2",
            ),
            "generation_s": 0.0,
            "runtime_available": bool(
                representative_item.get("runtime_available", False)
            ),
            "generation_request": generation_metadata,
            "crossing_alias_of": alias_representative.key,
            "crossing_signature": list(crossing_signature_by_key[entry.key]),
            "input_crossing_map": _process_input_crossing_map(
                alias_representative.process,
                entry.process,
                color_accuracy=str(getattr(args, "color_accuracy", "lc")),
                options=_process_options(args),
            ),
        }
        alias_generated.append(alias_item)
        generated.append(alias_item)
    if failures:
        if args.json:
            print(
                json.dumps(
                    {"available": False, "failures": failures},
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            first = failures[0]
            print(
                f"failed to generate {first['key']}: {first.get('error')}",
                file=sys.stderr,
            )
        return 1

    merged: dict[str, dict[str, object]] = dict(existing_entries)
    for item in generated:
        merged[str(item["key"])] = item
    if args.append and existing is not None and not args.replace:
        ordered_keys = [
            str(entry["key"])
            for entry in _process_set_manifest_entries(existing)
            if isinstance(entry, dict) and "key" in entry
        ]
        for entry in entries:
            if entry.key not in ordered_keys:
                ordered_keys.append(entry.key)
    else:
        ordered_keys = [entry.key for entry in entries]
        for key in existing_entries:
            if key not in ordered_keys:
                ordered_keys.append(key)
    process_entries = [merged[key] for key in ordered_keys if key in merged]
    existing_default_key = (
        str(existing.get("default_process_key"))
        if isinstance(existing, dict)
        and existing.get("default_process_key") in {entry["key"] for entry in process_entries}
        else None
    )
    default_process_key = (
        existing_default_key
        if args.append and existing is not None and not args.replace
        else (str(process_entries[0]["key"]) if process_entries else None)
    )
    runtime_available = all(
        bool(entry.get("runtime_available", False)) for entry in process_entries
    )
    process_set_manifest = {
        "schema_version": 2,
        "kind": "pyamplicol-generic-dag-process-set",
        "artifact_class": "generic-dag-schema-v2",
        "request": getattr(process_set, "request"),
        "default_process_key": default_process_key,
        "color_accuracy": getattr(args, "color_accuracy", "lc"),
        "generic_generation": generation_metadata,
        "runtime_available": runtime_available,
        "runtime_unavailable_message": (
            None
            if runtime_available
            else "one or more subprocesses do not contain serialized generic stage evaluators"
        ),
        "processes": process_entries,
    }
    manifest_path.write_text(
        json.dumps(process_set_manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (root / "check_standalone.py").write_text(
        _PROCESS_SET_STANDALONE_CHECK_SCRIPT,
        encoding="utf-8",
    )
    if args.json:
        print(
            json.dumps(
                {
                    "available": True,
                    "manifest": str(manifest_path),
                    "generated": generated,
                    "crossing_aliases": alias_generated,
                    "skipped": skipped,
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        _display(args).print_table(
            "Generated Process Set",
            _kv_columns(),
            [
                DisplayRow({"metric": "Request", "value": getattr(process_set, "request")}, "bold"),
                {"metric": "Entries", "value": len(process_entries)},
                {"metric": "Skipped", "value": len(skipped)},
                {"metric": "Manifest", "value": manifest_path},
            ],
        )
    return 0


def _generic_process_set_generation_metadata(
    args: argparse.Namespace,
) -> dict[str, object]:
    return {
        "process_options": {
            "flavour_scheme": int(getattr(args, "flavour_scheme", 5)),
            "include_3qqbar": bool(getattr(args, "include_3qqbar", False)),
            "include_cc": bool(getattr(args, "include_cc", False)),
            "include_resonance": bool(getattr(args, "include_resonance", False)),
        },
        "pruning": _json_ready_generic_dag_pruning_kwargs(
            _generic_dag_pruning_kwargs(args, process=None)
        ),
        "coupling_order_policy": str(getattr(args, "coupling_order_policy", "all")),
        "lc_sector_strategy": str(
            getattr(args, "lc_sector_strategy", "topology-representatives")
        ),
        "lc_topology_replay": bool(getattr(args, "lc_topology_replay", False)),
        "enumeration_prefilter": bool(
            getattr(args, "enumeration_prefilter", True)
        ),
        "n_cores": int(getattr(args, "n_cores", 1)),
    }


def _process_set_progress_metadata(process_set: object, worker_count: int) -> str:
    report = getattr(process_set, "selection_report", None)
    if report is None:
        return f"workers={worker_count}"
    return (
        f"workers={worker_count} "
        f"sel={int(getattr(report, 'selected_count', 0))} "
        f"rej={int(getattr(report, 'rejected_count', 0))}"
    )


def _generate_process_child_command(
    args: argparse.Namespace,
    process: str,
    output_dir: Path,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "pyamplicol",
        "generate-process",
        process,
        str(output_dir),
        "--json",
        "--color-accuracy",
        str(getattr(args, "color_accuracy", "lc")),
        "--flavour-scheme",
        str(int(getattr(args, "flavour_scheme", 5))),
        "--batch-size",
        str(int(getattr(args, "batch_size", 128))),
        "--symbolica-n-cores",
        str(int(getattr(args, "symbolica_n_cores", 4))),
        "--symbolica-evaluator-backend",
        str(getattr(args, "symbolica_evaluator_backend", "compiled-complex")),
        "--symbolica-iterations",
        str(int(getattr(args, "symbolica_iterations", 10))),
        "--symbolica-compiled-preset",
        str(getattr(args, "symbolica_compiled_preset", "runtime-o3")),
        "--symbolica-compiled-inline-asm",
        str(getattr(args, "symbolica_compiled_inline_asm", "default")),
        "--symbolica-compiled-optimization-level",
        str(int(getattr(args, "symbolica_compiled_optimization_level", 3))),
        "--symbolica-compiled-chunk-compile-workers",
        str(int(getattr(args, "symbolica_compiled_chunk_compile_workers", 1))),
        "--n_cores",
        "1",
        "--max-currents",
        str(int(getattr(args, "max_currents", 50000))),
        "--max-color-sectors",
        str(int(getattr(args, "max_color_sectors", 20000))),
    ]
    if bool(getattr(args, "include_3qqbar", False)):
        command.append("--include-3qqbar")
    if bool(getattr(args, "include_cc", False)):
        command.append("--include-cc")
    if bool(getattr(args, "include_resonance", False)):
        command.append("--include-resonance")
    for item in getattr(args, "max_coupling_order", ()) or ():
        command.extend(["--max-coupling-order", str(item)])
    if getattr(args, "max_qcd_order", None) is not None:
        command.extend(["--max-qcd-order", str(int(args.max_qcd_order))])
    if getattr(args, "max_qed_order", None) is not None:
        command.extend(["--max-qed-order", str(int(args.max_qed_order))])
    if str(getattr(args, "coupling_order_policy", "all")) != "all":
        command.extend(
            [
                "--coupling-order-policy",
                str(getattr(args, "coupling_order_policy", "all")),
            ]
        )
    if getattr(args, "max_lc_current_line_groups", None) is not None:
        command.extend(
            [
                "--max-lc-current-line-groups",
                str(int(args.max_lc_current_line_groups)),
            ]
        )
    max_quark_lines = _max_quark_lines(args)
    if max_quark_lines is not None:
        command.extend(["--max-quark-lines", str(max_quark_lines)])
    if not bool(getattr(args, "closure_side_mask_pruning", True)):
        command.append("--no-closure-side-mask-pruning")
    if not bool(getattr(args, "color_order_mask_pruning", True)):
        command.append("--no-color-order-mask-pruning")
    if not bool(getattr(args, "species_reachability_pruning", True)):
        command.append("--no-species-reachability-pruning")
    color_accuracy = str(getattr(args, "color_accuracy", "lc")).lower()
    numerical_filter_current = getattr(args, "numerical_filter_current", None)
    if numerical_filter_current is None:
        numerical_filter_current = color_accuracy == "lc"
    numerical_current_merging = getattr(args, "numerical_current_merging", None)
    if numerical_current_merging is None:
        numerical_current_merging = color_accuracy == "lc"
    if bool(numerical_filter_current) and color_accuracy != "lc":
        command.append("--numerical-filter-current")
    elif not bool(numerical_filter_current):
        command.append("--no-numerical-filter-current")
    if bool(numerical_current_merging) and color_accuracy != "lc":
        command.append("--numerical-current-merging")
    elif not bool(numerical_current_merging):
        command.append("--no-numerical-current-merging")
    if int(getattr(args, "numerical_current_samples", 10)) != 10:
        command.extend(
            [
                "--numerical-current-samples",
                str(int(getattr(args, "numerical_current_samples", 10))),
            ]
        )
    if int(getattr(args, "numerical_current_seed", 12345)) != 12345:
        command.extend(
            [
                "--numerical-current-seed",
                str(int(getattr(args, "numerical_current_seed", 12345))),
            ]
        )
    if (
        float(getattr(args, "numerical_current_relative_tolerance", 1.0e-12))
        != 1.0e-12
    ):
        command.extend(
            [
                "--numerical-current-relative-tolerance",
                str(
                    float(
                        getattr(
                            args,
                            "numerical_current_relative_tolerance",
                            1.0e-12,
                        )
                    )
                ),
            ]
        )
    if (
        float(getattr(args, "numerical_current_zero_tolerance", 1.0e-300))
        != 1.0e-300
    ):
        command.extend(
            [
                "--numerical-current-zero-tolerance",
                str(float(getattr(args, "numerical_current_zero_tolerance", 1.0e-300))),
            ]
        )
    ignore_particles = str(getattr(args, "ignore_particles", ""))
    if ignore_particles:
        command.extend(["--ignore-particles", ignore_particles])
    ignore_vertex_kinds = str(getattr(args, "ignore_vertex_kinds", ""))
    if ignore_vertex_kinds:
        command.extend(["--ignore-vertex-kinds", ignore_vertex_kinds])
    command.extend(
        [
            "--lc-sector-strategy",
            str(getattr(args, "lc_sector_strategy", "topology-representatives")),
        ]
    )
    lc_sector_ids = str(getattr(args, "lc_sector_ids", ""))
    if lc_sector_ids:
        command.extend(["--lc-sector-ids", lc_sector_ids])
    reference_color_order = str(getattr(args, "reference_color_order", ""))
    if reference_color_order:
        command.extend(["--reference-color-order", reference_color_order])
    if bool(getattr(args, "lc_topology_replay", False)):
        command.append("--lc-topology-replay")
    optional_integer_options = (
        ("symbolica_cpe_iterations", "--symbolica-cpe-iterations"),
        ("symbolica_jit_optimization_level", "--symbolica-jit-optimization-level"),
        ("symbolica_max_horner_scheme_variables", "--symbolica-max-horner-scheme-variables"),
        ("symbolica_max_common_pair_cache_entries", "--symbolica-max-common-pair-cache-entries"),
        ("symbolica_max_common_pair_distance", "--symbolica-max-common-pair-distance"),
        ("symbolica_compiled_output_chunk_size", "--symbolica-compiled-output-chunk-size"),
    )
    for attribute, option in optional_integer_options:
        value = getattr(args, attribute, None)
        if value is not None:
            command.extend([option, str(int(value))])
    compiler_path = getattr(args, "symbolica_compiler_path", None)
    if compiler_path:
        command.extend(["--symbolica-compiler-path", str(compiler_path)])
    for flag in getattr(args, "symbolica_compiler_flags", ()) or ():
        command.extend(["--symbolica-compiler-flag", str(flag)])
    if bool(getattr(args, "symbolica_collect_factors", False)):
        command.append("--symbolica-collect-factors")
    if bool(getattr(args, "symbolica_split_vertex_current_stages", False)):
        command.append("--symbolica-split-vertex-current-stages")
    if bool(
        getattr(
            args,
            "symbolica_stage_local_parameter_layout",
            _DEFAULT_SYMBOLICA_STAGE_LOCAL_PARAMETER_LAYOUT,
        )
    ):
        command.append("--symbolica-stage-local-parameter-layout")
    else:
        command.append("--no-symbolica-stage-local-parameter-layout")
    if bool(getattr(args, "merge_evaluators_strategy", False)):
        command.append("--merge-evaluators-strategy")
    else:
        command.append("--no-merge-evaluators-strategy")
    if bool(getattr(args, "symbolica_compiled_native", True)):
        command.append("--symbolica-compiled-native")
    else:
        command.append("--symbolica-no-compiled-native")
    jit_direct = getattr(args, "symbolica_jit_direct_translation", None)
    if jit_direct is True:
        command.append("--symbolica-jit-direct-translation")
    elif jit_direct is False:
        command.append("--symbolica-no-jit-direct-translation")
    if bool(getattr(args, "symbolica_direct_translation", True)):
        command.append("--symbolica-direct-translation")
    else:
        command.append("--symbolica-no-direct-translation")
    return command


def _process_crossing_reuse_signature(
    process: str,
    *,
    color_accuracy: str,
    options: ProcessOptions,
) -> tuple[int, ...]:
    del color_accuracy
    parsed = ProcessEnumerator(options).parse(process)
    initial_pdgs = tuple(int(PDGS[ANTI_PARTICLE[p]]) for p in parsed.initial_state)
    final_pdgs = tuple(sorted(int(PDGS[p]) for p in parsed.rest))
    return (
        len(initial_pdgs),
        *initial_pdgs,
        len(final_pdgs),
        *final_pdgs,
    )


def _process_input_crossing_map(
    representative_process: str,
    selected_process: str,
    *,
    color_accuracy: str,
    options: ProcessOptions,
) -> list[dict[str, int | float]]:
    from .process_ir import build_process_ir

    representative = build_process_ir(
        representative_process,
        color_accuracy=color_accuracy,
        options=options,
    )
    selected = build_process_ir(
        selected_process,
        color_accuracy=color_accuracy,
        options=options,
    )
    selected_legs = _external_all_outgoing_legs(selected)
    used: set[int] = set()
    mapping: list[dict[str, int | float]] = []
    for target_index, target in enumerate(_external_all_outgoing_legs(representative)):
        for source_index, source in enumerate(selected_legs):
            if source_index in used:
                continue
            if int(source["outgoing_pdg"]) != int(target["outgoing_pdg"]):
                continue
            used.add(source_index)
            mapping.append(
                {
                    "target_index": target_index,
                    "source_index": source_index,
                    "sign": float(target["sign"]) * float(source["sign"]),
                }
            )
            break
        else:
            raise ValueError(
                "could not build crossing map from "
                f"{selected_process!r} to {representative_process!r}"
            )
    return mapping


def _external_all_outgoing_legs(process_ir: Any) -> list[dict[str, int]]:
    legs: list[dict[str, int]] = []
    for index, leg in enumerate(process_ir.legs):
        outgoing_pdg = getattr(leg, "outgoing_pdg", None)
        if outgoing_pdg is None:
            continue
        legs.append(
            {
                "external_index": index,
                "outgoing_pdg": int(outgoing_pdg),
                "sign": -1 if bool(getattr(leg, "is_initial", False)) else 1,
            }
        )
    return legs


def _child_generation_environment() -> dict[str, str]:
    env = dict(os.environ)
    source_path = str(Path(__file__).resolve().parents[2])
    current = env.get("PYTHONPATH")
    if current:
        env["PYTHONPATH"] = f"{source_path}{os.pathsep}{current}"
    else:
        env["PYTHONPATH"] = source_path
    env.setdefault("PYAMPLICOL_NO_PROGRESS", "1")
    env[_CHILD_PROGRESS_ENV] = "1"
    return env


def _generation_progress_callback(args: argparse.Namespace, process: str):
    child_callback = _child_generation_progress_callback(process)
    monitor_callback = (
        _stderr_generation_monitor_callback(process)
        if bool(getattr(args, "monitor", False))
        or os.environ.get("PYAMPLICOL_MONITOR") == "1"
        else None
    )
    if monitor_callback is not None and child_callback is not None:
        return _combined_progress_callback(monitor_callback, child_callback)
    return monitor_callback or child_callback


def _child_generation_progress_callback(process: str):
    if os.environ.get(_CHILD_PROGRESS_ENV) != "1":
        return None

    def callback(event: dict[str, object]) -> None:
        payload: dict[str, object] = {
            "process": process,
            "stage": str(event.get("stage", "")),
            "item": str(event.get("item", "")),
        }
        for key in ("increment", "total", "ram"):
            if key in event:
                payload[key] = event[key]
        print(
            _CHILD_PROGRESS_PREFIX + json.dumps(payload, sort_keys=True),
            file=sys.stderr,
            flush=True,
        )

    return callback


def _stderr_generation_monitor_callback(process: str):
    started = time.perf_counter()
    last_seen: tuple[str, str] | None = None
    last_emitted = 0.0

    def callback(event: dict[str, object]) -> None:
        nonlocal last_seen, last_emitted
        stage = str(event.get("stage", ""))[:26]
        item = str(event.get("item", ""))[:66]
        now = time.perf_counter()
        signature = (stage, item)
        if signature == last_seen and now - last_emitted < 2.0:
            return
        last_seen = signature
        last_emitted = now
        wall = time.time()
        millis = int((wall - int(wall)) * 1000.0)
        stamp = time.strftime("%H:%M:%S", time.localtime(wall)) + f".{millis:03d}"
        elapsed = now - started
        extras: list[str] = []
        if "increment" in event:
            extras.append(f"inc={event['increment']}")
        if "total" in event:
            extras.append(f"total={event['total']}")
        if "ram" in event:
            extras.append(f"ram={event['ram']}")
        suffix = "" if not extras else " " + " ".join(str(item) for item in extras)
        print(
            f"[{stamp}] generate-process "
            f"stage={stage:<26} item={item:<66} elapsed={elapsed:8.1f}s"
            f"{suffix} process={process}",
            file=sys.stderr,
            flush=True,
        )

    return callback


def _combined_progress_callback(
    display_callback,
    child_callback,
):
    if child_callback is None:
        return display_callback

    def callback(event: dict[str, object]) -> None:
        display_callback(event)
        child_callback(event)

    return callback


def _parse_child_generation_payload(stdout: str) -> dict[str, object]:
    text = stdout.strip()
    if not text:
        return {}
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        start = text.find("{")
        end = text.rfind("}")
        if start < 0 or end <= start:
            return {"raw_stdout": stdout}
        payload = json.loads(text[start : end + 1])
    return payload if isinstance(payload, dict) else {}


def _start_child_stderr_reader(
    process: subprocess.Popen[str],
    buffer: list[str],
    events: queue.SimpleQueue[tuple[int, dict[str, object]]],
) -> threading.Thread:
    def run() -> None:
        stream = getattr(process, "stderr", None)
        if stream is None:
            return
        while True:
            try:
                line = stream.readline()
            except (OSError, ValueError):
                return
            if not line:
                return
            buffer.append(line)
            payload = _parse_child_progress_event(line)
            if payload is not None:
                events.put((int(process.pid), payload))

    thread = threading.Thread(
        target=run,
        name=f"pyamplicol-child-stderr-{process.pid}",
        daemon=True,
    )
    thread.start()
    return thread


def _finish_generation_child_output(
    process: subprocess.Popen[str],
    *,
    pid: int,
    stderr_buffers: dict[int, list[str]],
    stderr_threads: dict[int, threading.Thread],
) -> tuple[str, str]:
    stdout = ""
    stderr_tail = ""
    stdout_stream = getattr(process, "stdout", None)
    if stdout_stream is not None and hasattr(stdout_stream, "read"):
        stdout = str(stdout_stream.read())
    else:
        stdout, stderr_tail = process.communicate()
    thread = stderr_threads.get(pid)
    if thread is not None:
        thread.join(timeout=1.0)
    stderr = "".join(stderr_buffers.get(pid, ())) + stderr_tail
    return stdout, stderr


def _join_child_stderr_threads(
    stderr_threads: dict[int, threading.Thread],
) -> None:
    for thread in tuple(stderr_threads.values()):
        thread.join(timeout=1.0)
    stderr_threads.clear()


def _parse_child_progress_event(line: str) -> dict[str, object] | None:
    if not line.startswith(_CHILD_PROGRESS_PREFIX):
        return None
    text = line[len(_CHILD_PROGRESS_PREFIX) :].strip()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        return None
    return payload if isinstance(payload, dict) else None


def _drain_child_progress_events(
    events: queue.SimpleQueue[tuple[int, dict[str, object]]],
    active: dict[int, tuple[subprocess.Popen[str], object, Path, float]],
    progress: Any,
    seen: dict[int, tuple[str, str, float]] | None = None,
    *,
    min_interval_s: float = 0.5,
) -> None:
    while True:
        try:
            pid, event = events.get_nowait()
        except queue.Empty:
            return
        child = active.get(pid)
        if child is None:
            continue
        _, entry, _, _ = child
        stage = str(event.get("stage", "child"))
        item = str(event.get("item", ""))
        now = time.perf_counter()
        if seen is not None:
            previous_stage, previous_item, previous_at = seen.get(pid, ("", "", 0.0))
            changed_stage = stage != previous_stage
            changed_item = item != previous_item
            has_counter_update = "increment" in event or "total" in event
            if (
                not changed_stage
                and not has_counter_update
                and (not changed_item or now - previous_at < min_interval_s)
            ):
                continue
            if (
                not changed_stage
                and changed_item
                and not has_counter_update
                and now - previous_at < min_interval_s
            ):
                continue
            seen[pid] = (stage, item, now)
        progress.update(
            stage=f"child {stage}",
            item=f"{getattr(entry, 'key', pid)}:{item}",
            ram=_rss_text_for_active_processes(active),
        )


def _terminate_generation_children(
    active: dict[int, tuple[subprocess.Popen[str], object, Path, float]],
) -> None:
    root_pids = tuple(active)
    for process, _, _, _ in active.values():
        if process.poll() is None:
            with contextlib.suppress(ProcessLookupError, PermissionError):
                os.killpg(process.pid, signal.SIGTERM)
    for pid in sorted(_process_tree_pids(root_pids), reverse=True):
        if pid == os.getpid():
            continue
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.kill(pid, signal.SIGTERM)
    for process, _, _, _ in active.values():
        if process.poll() is None:
            with contextlib.suppress(ProcessLookupError, PermissionError):
                process.terminate()
    deadline = time.perf_counter() + 5.0
    for process, _, _, _ in active.values():
        remaining = max(deadline - time.perf_counter(), 0.0)
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            with contextlib.suppress(ProcessLookupError, PermissionError):
                os.killpg(process.pid, signal.SIGKILL)
            with contextlib.suppress(ProcessLookupError, PermissionError):
                process.kill()
    for process, _, _, _ in active.values():
        if process.poll() is None:
            with contextlib.suppress(ProcessLookupError, PermissionError):
                os.killpg(process.pid, signal.SIGKILL)
    for pid in sorted(_process_tree_pids(root_pids), reverse=True):
        if pid == os.getpid():
            continue
        with contextlib.suppress(ProcessLookupError, PermissionError):
            os.kill(pid, signal.SIGKILL)
    for process, _, _, _ in active.values():
        with contextlib.suppress(Exception):
            process.communicate(timeout=1.0)
    active.clear()


def _rss_text_for_active_processes(
    active: dict[int, tuple[subprocess.Popen[str], object, Path, float]],
) -> str:
    rss = _process_tree_rss_bytes((os.getpid(), *tuple(active)))
    return _format_bytes(rss) if rss is not None else "n/a"


def _process_tree_rss_bytes(root_pids: Sequence[int]) -> int | None:
    if not root_pids:
        return 0
    try:
        result = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return None
    if result.returncode != 0:
        return None
    rss_by_pid: dict[int, int] = {}
    children: dict[int, list[int]] = {}
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        try:
            pid = int(parts[0])
            ppid = int(parts[1])
            rss_kb = int(parts[2])
        except ValueError:
            continue
        rss_by_pid[pid] = max(rss_kb, 0) * 1024
        children.setdefault(ppid, []).append(pid)
    seen: set[int] = set()
    stack = list(root_pids)
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        stack.extend(children.get(pid, ()))
    return sum(rss_by_pid.get(pid, 0) for pid in seen)


def _process_tree_pids(root_pids: Sequence[int]) -> set[int]:
    if not root_pids:
        return set()
    try:
        result = subprocess.run(
            ["ps", "-axo", "pid=,ppid="],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return set(root_pids)
    if result.returncode != 0:
        return set(root_pids)
    children: dict[int, list[int]] = {}
    for line in result.stdout.splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        try:
            pid = int(parts[0])
            ppid = int(parts[1])
        except ValueError:
            continue
        children.setdefault(ppid, []).append(pid)
    seen: set[int] = set()
    stack = list(root_pids)
    while stack:
        pid = stack.pop()
        if pid in seen:
            continue
        seen.add(pid)
        stack.extend(children.get(pid, ()))
    return seen


def _format_bytes(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if amount < 1024.0 or unit == "TB":
            return f"{amount:.1f}{unit}" if unit != "B" else f"{int(amount)}B"
        amount /= 1024.0
    return f"{amount:.1f}TB"


def _load_process_set_manifest(root: Path) -> dict[str, object] | None:
    path = root / "process_set_manifest.json"
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"process-set manifest is not a JSON object: {path}")
    if payload.get("kind") not in {
        "pyamplicol-rusticol-process-set",
        "pyamplicol-generic-dag-process-set",
    }:
        raise ValueError(f"unsupported process-set artifact kind: {payload.get('kind')!r}")
    return payload


def _resolve_process_artifact_dir(root: Path, process_key: str | None = None) -> Path:
    if (root / "process_manifest.json").exists():
        if process_key is not None:
            raise ValueError("--process can only be used with a process-set artifact")
        return root
    manifest = _load_process_set_manifest(root)
    if manifest is None:
        raise ValueError(f"no process artifact manifest found in {root}")
    entries = _process_set_manifest_entries(manifest)
    if not entries:
        raise ValueError(f"process-set artifact contains no subprocesses: {root}")
    selected = process_key or str(manifest.get("default_process_key"))
    for entry in entries:
        key = str(entry.get("key"))
        process = str(entry.get("process"))
        if selected in {key, process}:
            path = Path(str(entry.get("path", "")))
            return path if path.is_absolute() else root / path
    available = ", ".join(str(entry.get("key")) for entry in entries)
    raise ValueError(f"process {selected!r} not found in {root}; available: {available}")


def _process_set_manifest_entries(payload: dict[str, object]) -> list[dict[str, object]]:
    entries = payload.get("processes", [])
    if not isinstance(entries, list):
        return []
    return [entry for entry in entries if isinstance(entry, dict)]


def _select_process_set_manifest_entry(
    payload: dict[str, object],
    process_key: str | None,
) -> dict[str, object]:
    entries = _process_set_manifest_entries(payload)
    if not entries:
        raise ValueError("process-set artifact contains no subprocesses")
    selected = process_key or str(payload.get("default_process_key"))
    for entry in entries:
        if selected in {str(entry.get("key")), str(entry.get("process"))}:
            return entry
    available = ", ".join(str(entry.get("key")) for entry in entries)
    raise ValueError(f"process {selected!r} not found; available: {available}")


def _cmd_time_process(args: argparse.Namespace) -> int:
    requested_root = Path(args.process_dir).expanduser()
    root = requested_root
    display = _display(args)
    try:
        import numpy as np
        import rusticol  # type: ignore[import-not-found]

        rusticol_build = _rusticol_build_metadata(rusticol)
        _require_release_rusticol(rusticol_build)
        with display.progress("Loading Rusticol process", metadata=str(root)):
            process_set_manifest = _load_process_set_manifest(requested_root)
            selected_process_entry = (
                _select_process_set_manifest_entry(
                    process_set_manifest,
                    args.process_key,
                )
                if process_set_manifest is not None
                else None
            )
            root = _resolve_process_artifact_dir(requested_root, args.process_key)
            manifest = json.loads((root / "process_manifest.json").read_text())
            unavailable = _generic_dag_runtime_unavailable_message(manifest)
            if unavailable is not None:
                raise RuntimeError(unavailable)
            runtime = _load_rusticol_runtime_compatible(
                rusticol,
                requested_root if process_set_manifest is not None else root,
                process_key=(
                    args.process_key if process_set_manifest is not None else None
                ),
                model_parameters=args.model_parameters,
            )
            runtime_metadata = dict(runtime.metadata())
            points = _load_rusticol_validation_momenta(root, int(args.precision), np)
            if selected_process_entry is not None:
                points = _validation_momenta_for_selected_crossing(
                    points,
                    selected_process_entry,
                    np,
                )
            values = _rusticol_evaluate(runtime, points, int(args.precision))
        with display.progress(
            "Profiling Rusticol runtime",
            metadata=f"precision={int(args.precision)}, target={float(args.target_runtime):.3g}s",
        ):
            profile_payload = _profile_rusticol_process(
                runtime,
                points,
                precision=int(args.precision),
                target_s=float(args.target_runtime),
                batch_size=_time_process_batch_size(args, manifest),
                np_module=np,
            )
    except Exception as exc:
        if args.json:
            print(
                json.dumps(
                    {
                        "available": False,
                        "error": str(exc),
                        "process_dir": str(requested_root),
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            print(str(exc), file=sys.stderr)
        return 1

    payload = {
        "available": True,
        "process_dir": str(root),
        "requested_process_dir": str(requested_root),
        "process": runtime_metadata.get("process") or manifest.get("process"),
        "artifact_class": (
            manifest.get("artifact_class")
            or runtime_metadata.get("artifact_class")
            or "generic-dag-schema-v2"
        ),
        "schema_version": runtime_metadata.get("schema_version"),
        "precision": int(args.precision),
        "target_runtime_s": float(args.target_runtime),
        "model_parameters": (
            None
            if args.model_parameters is None
            else str(Path(args.model_parameters).expanduser())
        ),
        "rusticol_build": rusticol_build,
        "values": [float(value) for value in values],
        "profile": profile_payload,
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        display.print_table(
            "Rusticol Timing Summary",
            _kv_columns(),
            [
                DisplayRow({"metric": "Process", "value": payload["process"]}, "bold"),
                {"metric": "Process dir", "value": root},
                {"metric": "Requested dir", "value": requested_root},
                {"metric": "Artifact class", "value": payload.get("artifact_class")},
                {
                    "metric": "Schema",
                    "value": payload.get("schema_version"),
                },
                {
                    "metric": "Rusticol build",
                    "value": (
                        f"{rusticol_build.get('profile', 'unknown')} "
                        f"[{rusticol_build.get('target', 'unknown')}]"
                    ),
                },
                {"metric": "Rusticol module", "value": rusticol_build.get("module_path")},
                {"metric": "Precision", "value": payload["precision"]},
                {
                    "metric": "Target runtime",
                    "value": format_measurement(
                        float(payload["target_runtime_s"]),
                        unit="s",
                    ),
                },
                {"metric": "Validation value(s)", "value": payload["values"]},
                DisplayRow(
                    {
                        "metric": "Wall runtime",
                        "value": format_measurement(
                            float(profile_payload["wall_us_per_point"]),
                            _float_or_none(
                                profile_payload.get("wall_us_per_point_error")
                            ),
                            unit="us/point",
                        ),
                    },
                    "green",
                ),
                DisplayRow(
                    {
                        "metric": "Evaluator + input pack",
                        "value": format_measurement(
                            float(profile_payload["core_evaluator_us_per_point"]),
                            _float_or_none(
                                profile_payload.get(
                                    "core_evaluator_us_per_point_error"
                                )
                            ),
                            unit="us/point",
                        ),
                    },
                    "cyan",
                ),
                DisplayRow(
                    {
                        "metric": "Pure evaluator calls",
                        "value": format_measurement(
                            float(profile_payload["pure_evaluator_us_per_point"]),
                            _float_or_none(
                                profile_payload.get(
                                    "pure_evaluator_us_per_point_error"
                                )
                            ),
                            unit="us/point",
                        ),
                    },
                    "cyan",
                ),
                {
                    "metric": "Input pack",
                    "value": format_measurement(
                        float(profile_payload["input_pack_us_per_point"]),
                        _float_or_none(
                            profile_payload.get("input_pack_us_per_point_error")
                        ),
                        unit="us/point",
                    ),
                },
                {
                    "metric": "Stage input pack",
                    "value": format_measurement(
                        float(profile_payload["stage_input_pack_us_per_point"]),
                        _float_or_none(
                            profile_payload.get(
                                "stage_input_pack_us_per_point_error"
                            )
                        ),
                        unit="us/point",
                    ),
                },
                {
                    "metric": "Amplitude input pack",
                    "value": format_measurement(
                        float(profile_payload["amplitude_input_pack_us_per_point"]),
                        _float_or_none(
                            profile_payload.get(
                                "amplitude_input_pack_us_per_point_error"
                            )
                        ),
                        unit="us/point",
                    ),
                },
                {
                    "metric": "Batch size",
                    "value": profile_payload["block_size"],
                },
                {
                    "metric": "Samples",
                    "value": (
                        f"{profile_payload['samples']} "
                        f"({profile_payload['block_count']} x "
                        f"{profile_payload['block_size']})"
                    ),
                },
            ],
        )
        breakdown_rows = _rusticol_breakdown_rows(profile_payload)
        if breakdown_rows:
            display.print_table(
                "Rusticol Runtime Breakdown",
                _kv_columns(),
                breakdown_rows,
            )
    if not bool(getattr(args, "_suppress_rusticol_exit", False)):
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(0)
    return 0


def _load_rusticol_runtime_compatible(
    rusticol_module: Any,
    process_dir: str | Path,
    *,
    process_key: str | None = None,
    model_parameters: str | Path | None = None,
) -> Any:
    """Load Rusticol while tolerating older local PyO3 signatures.

    Source builds accept keyword arguments for process-set selection and model
    parameters. Some existing editable installs expose an older positional-only
    loader. Keep `time-process` usable in both environments while still
    preferring the explicit new signature.
    """

    runtime_cls = getattr(rusticol_module, "Runtime")
    process_dir_text = str(Path(process_dir).expanduser())
    model_parameters_text = (
        None if model_parameters is None else str(Path(model_parameters).expanduser())
    )
    try:
        return runtime_cls.load(
            process_dir_text,
            process_key=process_key,
            model_parameters=model_parameters_text,
        )
    except TypeError as exc:
        if model_parameters_text is not None:
            raise
        if process_key is None:
            return runtime_cls.load(process_dir_text)
        try:
            return runtime_cls.load(process_dir_text, process_key)
        except TypeError:
            raise exc


def _rusticol_build_metadata(rusticol_module: Any) -> dict[str, str | None]:
    module = rusticol_module
    if not callable(getattr(module, "build_profile", None)):
        try:
            import importlib

            module = importlib.import_module("rusticol.rusticol")
        except Exception:  # noqa: BLE001 - metadata helper reports best effort.
            module = rusticol_module
    build_profile = getattr(module, "build_profile", None)
    build_target = getattr(module, "build_target", None)
    return {
        "profile": str(build_profile()) if callable(build_profile) else None,
        "target": str(build_target()) if callable(build_target) else None,
        "module_path": str(
            getattr(module, "__file__", getattr(rusticol_module, "__file__", ""))
        ),
    }


def _require_release_rusticol(build_metadata: Mapping[str, str | None]) -> None:
    profile = build_metadata.get("profile")
    if profile == "release":
        return
    if os.environ.get("PYAMPLICOL_ALLOW_DEBUG_RUSTICOL") == "1":
        return
    if profile is None:
        raise RuntimeError(
            "rusticol does not expose build-profile metadata. Reinstall it with "
            "`maturin develop --release` or run "
            "`python pyAmpliCol/dependencies/install_dependencies.py` before "
            "timing."
        )
    raise RuntimeError(
        "rusticol was built with Cargo profile "
        f"{profile!r}; timing requires the release PyO3 extension. Reinstall it "
        "with `maturin develop --release` or rerun the pyAmpliCol dependency "
        "installer."
    )


def _generic_dag_runtime_unavailable_message(
    manifest: dict[str, object],
) -> str | None:
    if manifest.get("kind") != "pyamplicol-generic-dag-process":
        return None
    compiled = manifest.get("compiled")
    if isinstance(compiled, dict):
        if bool(compiled.get("runtime_available", False)) and isinstance(
            compiled.get("stage_evaluators"),
            dict,
        ):
            return None
        message = compiled.get("runtime_unavailable_message")
        if isinstance(message, str) and message:
            return message
    return (
        "generic DAG evaluator stages are symbolically lowered, but serialized "
        "evaluator emission and Rusticol schema-v2 execution are not available"
    )


def _load_rusticol_validation_momenta(
    root: Path,
    precision: int,
    np_module: Any,
) -> Any:
    payload = json.loads((root / "validation_momenta.json").read_text())
    if payload.get("available") is False or not payload.get("points"):
        detail = payload.get("error") or "no validation momenta are bundled"
        raise RuntimeError(
            f"process artifact has no usable bundled validation momenta: {detail}"
        )
    points = [
        [
            [
                Decimal(_upcast_decimal_literal(str(component), precision))
                for component in particle["momentum"]
            ]
            for particle in point
        ]
        for point in payload["points"]
    ]
    if precision == 16:
        return np_module.asarray(points, dtype=np_module.float64)
    return points


def _upcast_decimal_literal(value: str, precision: int) -> str:
    if precision <= 16:
        return value
    text = value.strip()
    if not text:
        return text
    exponent = ""
    mantissa = text
    for marker in ("e", "E"):
        if marker in text:
            mantissa, exp = text.split(marker, 1)
            exponent = marker + exp
            break
    sign = ""
    if mantissa.startswith(("+", "-")):
        sign = mantissa[:1]
        mantissa = mantissa[1:]
    if "." in mantissa:
        integer, fraction = mantissa.split(".", 1)
    else:
        integer, fraction = mantissa, ""
    significant_digits = len((integer + fraction).lstrip("0"))
    if significant_digits == 0:
        significant_digits = 1
    padding = max(int(precision) - significant_digits, 0)
    return f"{sign}{integer}.{fraction}{'0' * padding}{exponent}"


def _validation_momenta_for_selected_crossing(
    points: Any,
    process_entry: dict[str, object],
    np_module: Any,
) -> Any:
    crossing_map = process_entry.get("input_crossing_map")
    if not isinstance(crossing_map, list) or not crossing_map:
        return points
    if hasattr(points, "shape"):
        mapped = np_module.empty_like(points)
        for item in crossing_map:
            if not isinstance(item, dict):
                raise RuntimeError("invalid input_crossing_map entry in process-set manifest")
            target_index = int(item["target_index"])
            source_index = int(item["source_index"])
            sign = float(item["sign"])
            mapped[:, source_index, :] = sign * points[:, target_index, :]
        return mapped

    mapped_points: list[list[list[Decimal]]] = []
    for point in points:
        mapped_point: list[list[Decimal] | None] = [None] * len(point)
        for item in crossing_map:
            if not isinstance(item, dict):
                raise RuntimeError("invalid input_crossing_map entry in process-set manifest")
            target_index = int(item["target_index"])
            source_index = int(item["source_index"])
            sign = -1 if float(item["sign"]) < 0.0 else 1
            mapped_point[source_index] = [
                (-component if sign < 0 else component)
                for component in point[target_index]
            ]
        if any(component is None for component in mapped_point):
            raise RuntimeError("input_crossing_map does not cover every selected leg")
        mapped_points.append(cast(list[list[Decimal]], mapped_point))
    return mapped_points


def _repeat_rusticol_points(points: Any, count: int, np_module: Any) -> Any:
    if hasattr(points, "shape"):
        reps = int(math.ceil(count / max(int(points.shape[0]), 1)))
        return np_module.tile(points, (reps, 1, 1))[:count]
    reps = int(math.ceil(count / max(len(points), 1)))
    return (points * reps)[:count]


def _rusticol_evaluate(runtime: Any, points: Any, precision: int) -> Sequence[Any]:
    if precision == 16:
        return runtime.evaluate(points)
    return runtime.evaluate_with_prec(points, precision)


def _time_process_batch_size(args: argparse.Namespace, manifest: dict[str, Any]) -> int:
    requested = getattr(args, "batch_size", None)
    if requested is not None:
        return max(int(requested), 1)
    metadata = manifest.get("metadata")
    if isinstance(metadata, dict):
        batch_size = metadata.get("batch_size")
        if isinstance(batch_size, int) and batch_size > 0:
            return batch_size
    return 128


def _mean_and_error(samples: Sequence[float]) -> tuple[float, float]:
    mean = statistics.fmean(samples) if samples else 0.0
    error = (
        statistics.stdev(samples) / math.sqrt(len(samples))
        if len(samples) > 1
        else 0.0
    )
    return mean, error


def _profile_rusticol_process(
    runtime: Any,
    points: Any,
    *,
    precision: int,
    target_s: float,
    batch_size: int,
    np_module: Any,
) -> dict[str, Any]:
    block_size = max(int(batch_size), 1)
    estimate_points = _repeat_rusticol_points(points, block_size, np_module)
    last_profile = dict(
        runtime.profile(estimate_points, precision=precision, include_values=False)
    )
    min_block_count = 8
    target_elapsed_s = max(float(target_s), 0.0)
    samples: list[float] = []
    core_samples: list[float] = []
    pure_evaluator_samples: list[float] = []
    stage_input_pack_samples: list[float] = []
    amplitude_input_pack_samples: list[float] = []
    input_pack_samples: list[float] = []
    elapsed_total_s = 0.0
    while len(samples) < min_block_count or elapsed_total_s < target_elapsed_s:
        batch = _repeat_rusticol_points(points, block_size, np_module)
        start = time.perf_counter()
        last_profile = dict(
            runtime.profile(batch, precision=precision, include_values=False)
        )
        elapsed_s = time.perf_counter() - start
        elapsed_total_s += elapsed_s
        samples.append(elapsed_s / block_size)
        core_samples.append(
            (
                float(last_profile.get("stage_evaluator_time_s", 0.0))
                + float(last_profile.get("amplitude_evaluator_time_s", 0.0))
            )
            / max(int(last_profile.get("points", block_size)), 1)
        )
        pure_evaluator_samples.append(
            (
                float(last_profile.get("stage_evaluator_call_time_s", 0.0))
                + float(
                    last_profile.get(
                        "amplitude_evaluator_call_time_s",
                        last_profile.get("amplitude_evaluator_time_s", 0.0),
                    )
                )
            )
            / max(int(last_profile.get("points", block_size)), 1)
        )
        points_in_profile = max(int(last_profile.get("points", block_size)), 1)
        stage_input_pack_s = float(last_profile.get("stage_input_pack_time_s", 0.0))
        amplitude_input_pack_s = float(
            last_profile.get("amplitude_input_pack_time_s", 0.0)
        )
        stage_input_pack_samples.append(stage_input_pack_s / points_in_profile)
        amplitude_input_pack_samples.append(
            amplitude_input_pack_s / points_in_profile
        )
        input_pack_samples.append(
            (stage_input_pack_s + amplitude_input_pack_s) / points_in_profile
        )
    wall_s, wall_error_s = _mean_and_error(samples)
    core_s, core_error_s = _mean_and_error(core_samples)
    pure_s, pure_error_s = _mean_and_error(pure_evaluator_samples)
    pack_s, pack_error_s = _mean_and_error(stage_input_pack_samples)
    amplitude_pack_s, amplitude_pack_error_s = _mean_and_error(
        amplitude_input_pack_samples
    )
    input_pack_s, input_pack_error_s = _mean_and_error(input_pack_samples)
    return {
        "samples": len(samples) * block_size,
        "block_count": len(samples),
        "block_size": block_size,
        "wall_us_per_point": wall_s * 1.0e6,
        "wall_us_per_point_error": wall_error_s * 1.0e6,
        "core_evaluator_us_per_point": core_s * 1.0e6,
        "core_evaluator_us_per_point_error": core_error_s * 1.0e6,
        "pure_evaluator_us_per_point": pure_s * 1.0e6,
        "pure_evaluator_us_per_point_error": pure_error_s * 1.0e6,
        "input_pack_us_per_point": input_pack_s * 1.0e6,
        "input_pack_us_per_point_error": input_pack_error_s * 1.0e6,
        "stage_input_pack_us_per_point": pack_s * 1.0e6,
        "stage_input_pack_us_per_point_error": pack_error_s * 1.0e6,
        "amplitude_input_pack_us_per_point": amplitude_pack_s * 1.0e6,
        "amplitude_input_pack_us_per_point_error": amplitude_pack_error_s * 1.0e6,
        "last_profile": _compact_rusticol_profile(last_profile),
    }


def _compact_rusticol_profile(profile: dict[str, Any]) -> dict[str, Any]:
    compact = dict(profile)
    values = compact.pop("values", None)
    if isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
        compact["value_count"] = len(values)
    return compact


def _cmd_generate_dag_evaluator_artifact(
    args: argparse.Namespace,
    runtime_backend: str,
) -> int:
    save_dir = getattr(args, "save_evaluator_dir", None)
    if save_dir is None:
        message = "--generate-only requires --save-evaluator-dir"
        if args.json:
            print(
                json.dumps(
                    {
                        "available": False,
                        "error": message,
                        "process": args.process,
                        "runtime_backend": runtime_backend,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            print(message, file=sys.stderr)
        return 1
    if runtime_backend != "rusticol":
        message = (
            "--generate-only only writes production generic DAG process "
            "artifacts through the Rusticol runtime. Use `generate-process` "
            "for new workflows."
        )
        if args.json:
            print(
                json.dumps(
                    {
                        "available": False,
                        "error": message,
                        "process": args.process,
                        "runtime_backend": runtime_backend,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
        else:
            print(message, file=sys.stderr)
        return 1

    delegate = argparse.Namespace(**vars(args))
    delegate.output_dir = Path(save_dir)
    delegate.color_accuracy = getattr(delegate, "color_accuracy", "lc")
    delegate.append = getattr(delegate, "append", False)
    delegate.replace = getattr(delegate, "replace", False)
    delegate.n_cores = getattr(delegate, "n_cores", 1)
    delegate.flavour_scheme = getattr(delegate, "flavour_scheme", 5)
    delegate.include_3qqbar = getattr(delegate, "include_3qqbar", False)
    delegate.include_cc = getattr(delegate, "include_cc", False)
    delegate.include_resonance = getattr(delegate, "include_resonance", False)
    delegate.parallel_process_enumeration = getattr(
        delegate,
        "parallel_process_enumeration",
        False,
    )
    return _cmd_generate_process(delegate)


def _format_error(value: object) -> str:
    if not isinstance(value, (float, int)):
        return "n/a"
    return f"{float(value):.3g}"


def _display(args: argparse.Namespace):
    return getattr(args, "_display", default_display_for_args(args))


def _kv_columns() -> tuple[DisplayColumn, DisplayColumn]:
    return (
        DisplayColumn("metric", "Metric"),
        DisplayColumn("value", "Value", align="right"),
    )


def _float_or_none(value: object) -> float | None:
    if isinstance(value, (float, int)):
        return float(value)
    return None


def _rusticol_breakdown_rows(
    profile_payload: dict[str, Any],
) -> list[DisplayRow | dict[str, object]]:
    last_profile = profile_payload.get("last_profile")
    if not isinstance(last_profile, dict):
        return []
    points = max(
        int(last_profile.get("points", profile_payload.get("block_size", 1))),
        1,
    )
    labels = {
        "source_fill_time_s": "Source fill",
        "parameter_pack_time_s": "Parameter pack",
        "stage_evaluator_time_s": "Stage evaluators",
        "amplitude_evaluator_time_s": "Amplitude evaluator",
        "me_reduction_time_s": "ME reduction",
        "output_transfer_time_s": "Output transfer",
        "wall_time_s": "Profile wall",
    }
    rows: list[DisplayRow | dict[str, object]] = []
    for key, label in labels.items():
        value = last_profile.get(key)
        if not isinstance(value, (float, int)):
            continue
        style = "cyan" if key in {"stage_evaluator_time_s", "amplitude_evaluator_time_s"} else None
        rows.append(
            DisplayRow(
                {
                    "metric": label,
                    "value": format_measurement(
                        float(value) * 1.0e6 / points,
                        unit="us/point",
                    ),
                },
                style,
            )
        )
    return rows


def _process_gluon_count_hint(process: str) -> int:
    if ">" not in process:
        return 0
    final_state = process.split(">", 1)[1]
    return sum(1 for token in final_state.lower().split() if token == "g")


def _charged_leptonic_w_effective_process(
    process: str,
    *,
    vector_pdg: int,
    gluon_count: int,
) -> str:
    initial, separator, _final = process.partition(">")
    if not separator:
        raise NativeEvaluationError(
            "invalid charged-current process string; expected 'initial > final'"
        )
    if vector_pdg == 24:
        vector_name = "w+"
    elif vector_pdg == -24:
        vector_name = "w-"
    else:
        raise NativeEvaluationError(
            f"charged-current effective process needs W+/W-, got {vector_pdg}"
        )
    final_tokens = [vector_name, *(["g"] * gluon_count)]
    return f"{initial.strip()} > {' '.join(final_tokens)}"


def _neutral_dilepton_effective_process(
    process: str,
    *,
    gluon_count: int,
) -> str:
    initial, separator, _final = process.partition(">")
    if not separator:
        raise NativeEvaluationError(
            "invalid neutral-dilepton process string; expected 'initial > final'"
        )
    final_tokens = ["z", *(["g"] * gluon_count)]
    return f"{initial.strip()} > {' '.join(final_tokens)}"


def _cmd_profile(args: argparse.Namespace) -> int:
    return _legacy_native_command_unavailable(args, "profile")


def _cmd_profile_reference(args: argparse.Namespace) -> int:
    from .evaluation import NativeRuntimeEvaluator
    from .legacy_matrix import NativeMatrixElementGenerator

    generator = NativeMatrixElementGenerator(cache_dir=args.cache_dir)
    result = generator.generate(args.process, options=_process_options(args))
    artifact = _load_optional_evaluator_artifact(args.cache_dir, args.process)
    per_point_runtime_s = None
    runtime_metadata = None
    try:
        evaluator = NativeRuntimeEvaluator(
            args.process,
            runtime_backend=_runtime_backend(args),
            allow_reference_legacy=True,
            **_runtime_evaluator_kwargs(args),
        )
        start = time.perf_counter()
        evaluator.evaluate()
        per_point_runtime_s = time.perf_counter() - start
        runtime_metadata = evaluator.metadata.to_json_dict()
    except NativeEvaluationError:
        per_point_runtime_s = None
    payload = result.to_json_dict()
    payload["profile"] = {
        "generation_time_s": result.generation_time_s,
        "tensor_network_reduction_s": _tensor_network_reduction_time(artifact),
        "symbolica_optimization_jit_s": _symbolic_evaluator_build_time(artifact),
        "artifact_cache_hit": result.artifact_cache_hit,
        "artifact_load_s": None if artifact is None else artifact.load_time_s,
        "per_point_runtime_s": per_point_runtime_s,
        "native_runtime_backend": runtime_metadata,
    }
    payload["profile_notes"] = (
        "Generation/cache/load/runtime timing is available for the staged native "
        "kernel; Symbolica evaluator build and tensor-network reduction timing are "
        "reported when an evaluator artifact is present."
    )
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(f"generation metadata time: {result.generation_time_s:.6g} s")
        if artifact is not None:
            print(f"artifact load time: {artifact.load_time_s:.6g} s")
        if runtime_metadata is not None:
            print(
                "native backend: "
                f"{runtime_metadata['backend']} ({runtime_metadata['kernel']})"
            )
        if per_point_runtime_s is not None:
            print(f"native per-point runtime: {per_point_runtime_s:.6g} s")
        print(payload["profile_notes"])
    return 0


def _process_enumeration_to_dict(enumeration: ProcessEnumeration) -> dict[str, object]:
    return {
        "request": {
            "initial_state": list(enumeration.request.initial_state),
            "jet_count": enumeration.request.jet_count,
            "rest": list(enumeration.request.rest),
            "leptons": list(enumeration.request.leptons),
        },
        "options": {
            "flavour_scheme": enumeration.options.flavour_scheme,
            "include_3qqbar": enumeration.options.include_3qqbar,
            "include_cc": enumeration.options.include_cc,
            "include_resonance": enumeration.options.include_resonance,
            "serial": enumeration.options.serial,
        },
        "n_external": enumeration.n_external,
        "n_unique_processes": len(enumeration.unique_processes),
        "n_groups": len(enumeration.groups),
        "n_records": enumeration.n_records,
        "unique_processes": [list(process) for process in enumeration.unique_processes],
        "groups": [
            {
                "group_id": group.group_id,
                "phase_space_order": [index + 1 for index in group.phase_space_order],
                "records": [
                    {
                        "process": list(record.process),
                        "color_order": [index + 1 for index in record.color_order],
                        "multichannel_partners": [
                            index + 1 for index in record.multichannel_partners
                        ],
                        "identical_factor": record.identical_factor,
                    }
                    for record in group.records
                ],
            }
            for group in enumeration.groups
        ],
    }


def _process_set_enumeration_to_dict(process_set: object) -> dict[str, object]:
    entries = getattr(process_set, "entries")
    entry_payloads = [
        {
            "key": entry.key,
            "process": entry.process,
            "enumeration": _process_enumeration_to_dict(entry.enumeration),
        }
        for entry in entries
    ]
    if len(entry_payloads) == 1:
        payload = dict(entry_payloads[0]["enumeration"])
        payload.update(
            {
                "process_set_request": getattr(process_set, "request"),
                "default_key": getattr(process_set, "default_key"),
                "n_entries": 1,
                "entries": entry_payloads,
            }
        )
        return payload
    return {
        "request": getattr(process_set, "request"),
        "default_key": getattr(process_set, "default_key"),
        "n_entries": len(entry_payloads),
        "n_unique_processes": sum(
            int(payload["enumeration"]["n_unique_processes"])
            for payload in entry_payloads
        ),
        "n_groups": sum(
            int(payload["enumeration"]["n_groups"]) for payload in entry_payloads
        ),
        "n_records": sum(
            int(payload["enumeration"]["n_records"]) for payload in entry_payloads
        ),
        "entries": entry_payloads,
    }


def _evaluation_to_dict(
    result: object,
    *,
    runtime_metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    particles = getattr(result, "particles")
    helicity_contributions = getattr(result, "helicity_contributions")
    inferred_backend = (
        "native-python-z-gluon"
        if any(particle.pdg == 21 for particle in particles[2:])
        else "native-python-zero-gluon"
    )
    backend = (
        runtime_metadata.get("backend", inferred_backend)
        if runtime_metadata is not None
        else inferred_backend
    )
    kernel = (
        runtime_metadata.get("kernel")
        if runtime_metadata is not None
        else None
    )
    return {
        "process": getattr(result, "process"),
        "evaluation_available": True,
        "backend": backend,
        "kernel": kernel,
        "matrix_element": getattr(result, "matrix_element"),
        "raw_helicity_sum": getattr(result, "raw_helicity_sum"),
        "color_factor": getattr(result, "color_factor"),
        "average_factor": getattr(result, "average_factor"),
        "identical_factor": getattr(result, "identical_factor", 1),
        "coupling_factor": getattr(result, "coupling_factor"),
        "particles": [
            {"pdg": particle.pdg, "momentum": list(particle.momentum)}
            for particle in particles
        ],
        "helicity_contributions": [
            {
                "helicities": list(contribution.helicities),
                "amplitude": {
                    "re": contribution.amplitude.real,
                    "im": contribution.amplitude.imag,
                },
                "squared": contribution.squared,
            }
            for contribution in helicity_contributions
        ],
    }


def _native_generation_profile(
    process: str,
    cache_dir: Path,
    options: ProcessOptions,
) -> dict[str, object]:
    from .legacy_matrix import NativeMatrixElementGenerator

    result = NativeMatrixElementGenerator(cache_dir=cache_dir).generate(
        process,
        options=options,
    )
    artifact = _load_optional_evaluator_artifact(cache_dir, process)
    artifact_payload = {} if artifact is None else artifact.payload
    full_me_ready = artifact_payload.get("full_me_tensor_network_ready")
    graph_counts = artifact_payload.get("graph_counts")
    blueprint = _tensor_network_blueprint_summary(artifact_payload)
    return {
        "cache_dir": str(cache_dir),
        "metadata_file": None if result.cache_file is None else str(result.cache_file),
        "artifact_file": (
            None if result.artifact_file is None else str(result.artifact_file)
        ),
        "backend": result.backend,
        "kernel": artifact_payload.get("kernel"),
        "artifact_fingerprint": result.artifact_fingerprint,
        "artifact_cache_hit": result.artifact_cache_hit,
        "artifact_load_s": None if artifact is None else artifact.load_time_s,
        "generation_time_s": result.generation_time_s,
        "supported_native_target": result.supported_native_target,
        "full_me_tensor_network_ready": full_me_ready,
        "graph_counts": graph_counts,
        "tensor_network_blueprint": blueprint,
    }


def _tensor_network_blueprint_summary(
    artifact_payload: dict[str, object],
) -> dict[str, object] | None:
    lowering = artifact_payload.get("symbolic_lowering")
    if not isinstance(lowering, dict):
        return None
    blueprint = lowering.get("tensor_network_blueprint")
    if not isinstance(blueprint, dict):
        return None
    keys = (
        "engine",
        "status",
        "expression_built",
        "expression_executed",
        "ready_interactions",
        "pending_interactions",
        "placeholder_vertex_kinds",
        "registered_tensor_names",
        "expression_length",
        "executed_expression_length",
        "execution_time_s",
        "full_me_tensor_network_ready",
        "propagator_lowering_ready",
    )
    return {key: blueprint.get(key) for key in keys}


def _load_optional_evaluator_artifact(
    cache_dir: Path,
    process: str,
) -> Any | None:
    from .legacy_matrix import load_evaluator_artifact

    try:
        return load_evaluator_artifact(cache_dir, process)
    except (OSError, ValueError):
        return None


def _artifact_status(
    cache_dir: Path,
    process: str,
    artifact: Any | None,
) -> dict[str, object]:
    from .legacy_matrix import evaluator_artifact_path

    if artifact is None:
        return {
            "available": False,
            "file": str(evaluator_artifact_path(cache_dir, process)),
        }
    payload = artifact.payload
    return {
        "available": True,
        "file": str(artifact.file),
        "load_time_s": artifact.load_time_s,
        "schema_version": payload.get("schema_version"),
        "artifact_fingerprint": payload.get("artifact_fingerprint"),
        "backend": payload.get("backend"),
        "kernel": payload.get("kernel"),
        "full_me_tensor_network_ready": payload.get("full_me_tensor_network_ready"),
        "symbolic_scalar_evaluator_ready": payload.get(
            "symbolic_scalar_evaluator_ready"
        ),
        "tensor_network_scalar_evaluator_ready": payload.get(
            "tensor_network_scalar_evaluator_ready"
        ),
    }


def _symbolic_artifact_evaluation(
    process: str,
    native_result: object,
    artifact: Any | None,
) -> dict[str, object]:
    if artifact is None:
        return {
            "available": False,
            "reason": "no evaluator artifact loaded",
        }
    kernel = artifact.payload.get("kernel")
    if kernel == "symbolica-zero-gluon":
        return _zero_gluon_artifact_evaluation(process, native_result, artifact)
    if kernel == "symbolica-one-gluon-tensor-network":
        return _z_gluon_tensor_network_artifact_evaluation(
            process,
            native_result,
            artifact,
            kernel=kernel,
        )
    if kernel == "symbolica-z-gluon-tensor-network":
        return _z_gluon_tensor_network_artifact_evaluation(
            process,
            native_result,
            artifact,
            kernel=kernel,
        )
    return {
        "available": False,
        "reason": f"artifact kernel is {kernel!r}",
    }


def _zero_gluon_artifact_evaluation(
    process: str,
    native_result: object,
    artifact: Any,
) -> dict[str, object]:
    from .symbolic import ZeroGluonSymbolicEvaluator

    evaluator_payload = artifact.payload.get("symbolic_scalar_evaluator")
    if not isinstance(evaluator_payload, dict):
        return {
            "available": False,
            "reason": "artifact has no symbolic scalar evaluator payload",
        }
    kernel = "symbolica-zero-gluon"
    try:
        start = time.perf_counter()
        symbolic = ZeroGluonSymbolicEvaluator.from_artifact_payload(evaluator_payload)
        symbolic_result = symbolic.evaluate(
            process,
            getattr(native_result, "particles"),
        )
        elapsed_s = time.perf_counter() - start
    except (KeyError, TypeError, ValueError, RuntimeError, NativeEvaluationError) as exc:
        return {
            "available": False,
            "kernel": kernel,
            "reason": str(exc),
        }

    return _symbolic_evaluation_success_payload(
        kernel=kernel,
        native_result=native_result,
        matrix_element=symbolic_result.matrix_element,
        raw_helicity_sum=symbolic_result.raw_helicity_sum,
        elapsed_s=elapsed_s,
    )


def _z_gluon_tensor_network_artifact_evaluation(
    process: str,
    native_result: object,
    artifact: Any,
    *,
    kernel: str,
) -> dict[str, object]:
    evaluator_payload = artifact.payload.get("tensor_network_scalar_evaluator")
    if not isinstance(evaluator_payload, dict):
        return {
            "available": False,
            "reason": "artifact has no tensor-network scalar evaluator payload",
        }
    try:
        from .tensor_runtime import ZGluonTensorNetworkEvaluator

        start = time.perf_counter()
        evaluator = ZGluonTensorNetworkEvaluator.from_artifact_payload(
            process,
            evaluator_payload,
        )
        tensor_result = evaluator.evaluate(
            getattr(native_result, "particles"),
        )
        elapsed_s = time.perf_counter() - start
    except (KeyError, TypeError, ValueError, RuntimeError, NativeEvaluationError) as exc:
        return {
            "available": False,
            "kernel": kernel,
            "reason": str(exc),
        }

    return _symbolic_evaluation_success_payload(
        kernel=kernel,
        native_result=native_result,
        matrix_element=tensor_result.matrix_element,
        raw_helicity_sum=tensor_result.raw_helicity_sum,
        elapsed_s=elapsed_s,
    )


def _symbolic_evaluation_success_payload(
    *,
    kernel: str,
    native_result: object,
    matrix_element: float,
    raw_helicity_sum: float,
    elapsed_s: float,
) -> dict[str, object]:
    native_me = float(getattr(native_result, "matrix_element"))
    denom = max(abs(native_me), abs(matrix_element), 1.0e-300)
    return {
        "available": True,
        "kernel": kernel,
        "matrix_element": matrix_element,
        "raw_helicity_sum": raw_helicity_sum,
        "relative_difference": abs(matrix_element - native_me) / denom,
        "evaluation_time_s": elapsed_s,
    }


def _symbolic_evaluator_build_time(
    artifact: Any | None,
) -> float | None:
    if artifact is None:
        return None
    evaluator_payload = artifact.payload.get("symbolic_scalar_evaluator")
    if not isinstance(evaluator_payload, dict):
        evaluator_payload = artifact.payload.get("tensor_network_scalar_evaluator")
    if not isinstance(evaluator_payload, dict):
        return None
    build_time = evaluator_payload.get("build_time_s")
    if not isinstance(build_time, (float, int)):
        metadata = evaluator_payload.get("metadata")
        if isinstance(metadata, dict):
            build_time = metadata.get("symbolica_evaluator_build_s")
    return float(build_time) if isinstance(build_time, (float, int)) else None


def _tensor_network_reduction_time(
    artifact: Any | None,
) -> float | None:
    if artifact is None:
        return None
    evaluator_payload = artifact.payload.get("tensor_network_scalar_evaluator")
    if not isinstance(evaluator_payload, dict):
        return None
    metadata = evaluator_payload.get("metadata")
    if not isinstance(metadata, dict):
        return None
    reduction_time = metadata.get("tensor_network_reduction_s")
    return float(reduction_time) if isinstance(reduction_time, (float, int)) else None


def _amplicol_dry_run_payload(
    args: argparse.Namespace,
    options: ProcessOptions,
    *,
    fixed_probe: bool = False,
    momenta_probe: bool = False,
) -> dict[str, object]:
    root = Path(args.amplicol_root)
    process_file = args.process_file or root / "processes.txt"
    legacy_me_test = bool(getattr(args, "me_test", False))
    commands = []
    if not args.skip_build:
        commands.extend(
            [["make", "cleanlib"], ["make", f"-j{args.jobs}", "amplicol_generate"]]
        )
        if args.amplicol_probe or momenta_probe:
            commands.extend(
                [
                    ["./amplicol_generate", "--library=create", f"--process={process_file}"],
                    ["make", f"-j{args.jobs}", "amplicol_generate_library"],
                ]
            )
    timing = args.points if args.timing is None else args.timing
    if momenta_probe:
        commands.append(
            [
                "./amplicol_generate",
                f"--amplicol_momenta_probe={args.points}",
                f"--timing={timing}",
                f"--process={process_file}",
                "--library=use",
            ]
        )
    elif args.amplicol_probe:
        probe_flag = "--amplicol_fixed_probe" if fixed_probe else "--amplicol_probe"
        commands.append(
            [
                "./amplicol_generate",
                f"{probe_flag}={args.points}",
                f"--timing={timing}",
                f"--process={process_file}",
                "--library=use",
            ]
        )
    elif legacy_me_test:
        commands.append(
            [
                "./amplicol_generate",
                f"--me_test={args.points}",
                f"--timing={timing}",
                f"--process={process_file}",
            ]
        )
    else:  # pragma: no cover - parser logic keeps one compare mode selected.
        raise RuntimeError("no AmpliCol comparison mode selected")
    mode = "me_test" if legacy_me_test else "amplicol_momenta_probe_library"
    if momenta_probe:
        mode = "amplicol_momenta_probe_library"
    elif args.amplicol_probe:
        mode = "amplicol_fixed_probe_library" if fixed_probe else "amplicol_probe_library"
    return {
        "mode": mode,
        "process": args.process,
        "options": {
            "flavour_scheme": options.flavour_scheme,
            "include_3qqbar": options.include_3qqbar,
            "include_cc": options.include_cc,
            "include_resonance": options.include_resonance,
            "serial": options.serial,
        },
        "amplicol_root": str(root),
        "process_file": str(process_file),
        "mg5_path": None if args.mg5_path is None else str(args.mg5_path),
        "pyamplicol_cache_dir": str(args.cache_dir),
        "pyamplicol_runtime_backend": args.runtime_backend,
        "commands": commands,
    }


def _z_gluon_family_process(gluon_count: int) -> str:
    if gluon_count == 0:
        return "d d~ > z"
    return "d d~ > z " + " ".join("g" for _ in range(gluon_count))


def _z_gluon_family_dry_run_payload(
    args: argparse.Namespace,
    options: ProcessOptions,
) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for gluon_count in range(args.min_gluons, args.max_gluons + 1):
        process = _z_gluon_family_process(gluon_count)
        timing = args.points if args.timing is None else args.timing
        probe_flag = "--amplicol_fixed_probe" if gluon_count == 0 else "--amplicol_probe"
        mode = "amplicol_fixed_probe_library" if gluon_count == 0 else "amplicol_probe_library"
        commands = [
            ["make", "cleanlib"],
            ["make", f"-j{args.jobs}", "amplicol_generate"],
            [
                "./amplicol_generate",
                "--library=create",
                "--process=processes.txt",
            ],
            ["make", f"-j{args.jobs}", "amplicol_generate_library"],
            [
                "./amplicol_generate",
                f"{probe_flag}={args.points}",
                f"--timing={timing}",
                "--process=processes.txt",
                "--library=use",
            ],
            [
                "./amplicol_generate",
                "--library=use",
                f"--nevents={timing}",
                "--seed=101",
                "--timing=1",
            ],
        ]
        rows.append(
            {
                "gluon_count": gluon_count,
                "process": process,
                "mode": mode,
                "fortran_timing_workflow": "generated_library_use",
                "runtime_backend": _z_gluon_family_runtime_backend(
                    gluon_count,
                    _runtime_backend(args),
                ),
                "commands": commands,
            }
        )
    return {
        "family": "d d~ -> Z + n g",
        "min_gluons": args.min_gluons,
        "max_gluons": args.max_gluons,
        "points": args.points,
        "rel_tol": args.rel_tol,
        "timing": args.points if args.timing is None else args.timing,
        "options": {
            "flavour_scheme": options.flavour_scheme,
            "include_3qqbar": options.include_3qqbar,
            "include_cc": options.include_cc,
            "include_resonance": options.include_resonance,
            "serial": options.serial,
        },
        "amplicol_root": str(args.amplicol_root),
        "pyamplicol_cache_dir": str(args.cache_dir),
        "pyamplicol_runtime_backend": args.runtime_backend,
        "rows": rows,
    }


def _z_gluon_family_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    max_relative_difference = None
    validated_rows = 0
    total_reference_points = 0
    total_native_points = 0
    total_command_time_s = 0.0
    pyamplicol_runtime_total_s = 0.0
    for row in rows:
        command_time = row.get("total_command_time_s")
        if isinstance(command_time, (float, int)):
            total_command_time_s += float(command_time)
        summary = row.get("native_probe_summary")
        if not isinstance(summary, dict):
            continue
        reference_points = int(summary.get("reference_points", 0))
        available_points = int(summary.get("available_points", 0))
        total_reference_points += reference_points
        total_native_points += available_points
        failures = summary.get("failures")
        if available_points == reference_points and not failures:
            validated_rows += 1
        row_max = summary.get("max_relative_difference")
        if isinstance(row_max, (float, int)):
            max_relative_difference = (
                float(row_max)
                if max_relative_difference is None
                else max(max_relative_difference, float(row_max))
            )
        runtime = row.get("native_runtime")
        if isinstance(runtime, dict):
            runtime_total = runtime.get("total_s")
            if isinstance(runtime_total, (float, int)):
                pyamplicol_runtime_total_s += float(runtime_total)
    return {
        "requested_rows": len(rows),
        "validated_rows": validated_rows,
        "reference_points": total_reference_points,
        "native_points": total_native_points,
        "max_relative_difference": max_relative_difference,
        "total_command_time_s": total_command_time_s,
        "pyamplicol_runtime_total_s": pyamplicol_runtime_total_s,
    }


def _z_gluon_family_passed(
    summary: dict[str, object],
    *,
    rel_tol: float,
) -> bool:
    requested_rows = summary.get("requested_rows")
    validated_rows = summary.get("validated_rows")
    max_relative_difference = summary.get("max_relative_difference")
    if not isinstance(requested_rows, int) or not isinstance(validated_rows, int):
        return False
    if requested_rows == 0 or validated_rows != requested_rows:
        return False
    if not isinstance(max_relative_difference, (float, int)):
        return False
    return float(max_relative_difference) <= rel_tol


def _z_gluon_family_runtime_backend(
    gluon_count: int,
    requested: str,
) -> str:
    if gluon_count == 0 and requested in ("dag", "numeric-tensor-network"):
        return "python"
    return requested


def _rusticol_generation_profile(
    process: str,
    cache_dir: Path,
    options: ProcessOptions,
    *,
    runtime_evaluator_kwargs: dict[str, Any],
    reference_color_order: Sequence[int] | None = None,
    generic_dag_pruning_kwargs: Mapping[str, Any] | None = None,
) -> tuple[dict[str, object], Path]:
    from .generic_artifact import (
        GenericProcessManifest,
        build_generic_process_manifest,
        select_leading_color_sector_ids_from_plan,
        write_generic_dag_process_artifact,
    )
    from .color_plan import build_color_plan
    from .process_ir import build_process_ir

    pruning_kwargs = _normalize_generic_dag_pruning_kwargs(generic_dag_pruning_kwargs)
    process_ir = build_process_ir(process, options=options)
    color_plan = build_color_plan(
        process_ir,
        color_accuracy=process_ir.color_accuracy,
        options=options,
    )
    selected_color_sector_ids = pruning_kwargs.pop("selected_color_sector_ids", None)
    if selected_color_sector_ids is None:
        selected_color_sector_ids = select_leading_color_sector_ids_from_plan(
            color_plan,
            reference_color_order=reference_color_order,
        )
    generic_manifest = build_generic_process_manifest(
        process_ir,
        options=options,
        reference_color_order=reference_color_order,
        selected_color_sector_ids=selected_color_sector_ids,
        **pruning_kwargs,
    )
    if selected_color_sector_ids is None:
        selected_color_sector_ids = _amplicol_reference_color_sector_ids(
            generic_manifest,
            reference_color_order=reference_color_order,
        )
    output_dir = Path(cache_dir).expanduser() / "generic-rusticol-compare" / process_ir.key
    if output_dir.exists():
        shutil.rmtree(output_dir)
    build_kwargs = dict(runtime_evaluator_kwargs)
    if (
        build_kwargs.get("symbolica_load_evaluator_dir") is None
        and build_kwargs.get("symbolica_compiled_output_dir") is None
    ):
        build_kwargs["symbolica_compiled_output_dir"] = str(output_dir / "compiled")
    settings = _symbolica_settings_from_runtime_kwargs(
        build_kwargs,
        process=process_ir.process,
    )
    start = time.perf_counter()
    manifest_path, manifest = write_generic_dag_process_artifact(
        generic_manifest,
        output_dir,
        options=options,
        evaluator_backend=str(build_kwargs["symbolica_evaluator_backend"]),
        compiled_preset=str(build_kwargs["symbolica_compiled_preset"]),
        batch_size=int(build_kwargs["batch_size"]),
        emit_stage_evaluator_artifacts=True,
        symbolica_settings=settings,
        merge_evaluators_strategy=bool(build_kwargs["merge_evaluators_strategy"]),
        stage_local_parameter_layout=bool(
            build_kwargs["stage_local_parameter_layout"]
        ),
        verbose_evaluator_build=bool(build_kwargs["verbose_evaluator_build"]),
        selected_color_sector_ids=selected_color_sector_ids,
        **pruning_kwargs,
    )
    generation_s = time.perf_counter() - start
    return (
        {
            "backend": "rusticol-generic-schema-v2",
            "generation_time_s": generation_s,
            "artifact_cache_hit": False,
            "artifact_load_s": None,
            "manifest": str(manifest_path),
            "process_dir": str(output_dir),
            "runtime_available": bool(
                cast(dict[str, Any], manifest["compiled"]).get("runtime_available")
            ),
            "evaluator_backend": str(build_kwargs["symbolica_evaluator_backend"]),
            "compiled_preset": str(build_kwargs["symbolica_compiled_preset"]),
            "validation_color_sector_ids": (
                None
                if selected_color_sector_ids is None
                else sorted(selected_color_sector_ids)
            ),
            "validation_color_sector_policy": (
                "first legacy LC ordering"
                if selected_color_sector_ids is not None
                else "full available LC sector set"
            ),
        },
        output_dir,
    )


def _amplicol_reference_color_sector_ids(
    manifest: Any,
    *,
    reference_color_order: Sequence[int] | None = None,
) -> set[int] | None:
    from .generic_artifact import select_leading_color_sector_ids

    return select_leading_color_sector_ids(
        manifest,
        reference_color_order=reference_color_order,
    )


def _amplicol_reference_color_order_for_run(run: object) -> tuple[int, ...] | None:
    if not hasattr(run, "process_file"):
        return None
    return reference_color_order_for_run(cast(Any, run))


def _amplicol_process_file_entry(
    process_file: Path,
    *,
    group: int,
    integral: int,
) -> dict[str, list[int]] | None:
    return amplicol_process_file_entry(process_file, group=group, integral=integral)


def _should_use_fixed_amplicol_probe(process: str, options: ProcessOptions) -> bool:
    parsed = ProcessEnumerator(options).parse(process)
    return parsed.jet_count == 0 and parsed.rest == ("z",)


def _first_point_to_dict(point: object | None) -> dict[str, object] | None:
    if point is None:
        return None
    return {
        "group": getattr(point, "group"),
        "integral": getattr(point, "integral"),
        "matrix_element": getattr(point, "matrix_element"),
        "particles": [
            {"pdg": particle.pdg, "momentum": list(particle.momentum)}
            for particle in getattr(point, "particles")
        ],
    }


def _probe_point_to_dict(point: object) -> dict[str, object]:
    return {
        "point": getattr(point, "point"),
        "group": getattr(point, "group"),
        "integral": getattr(point, "integral"),
        "matrix_element": getattr(point, "matrix_element"),
        "particles": [
            {"pdg": particle.pdg, "momentum": list(particle.momentum)}
            for particle in getattr(point, "particles")
        ],
    }


def _attach_native_probe_comparison(
    process: str,
    run: object,
    payload: dict[str, object],
    *,
    runtime_backend: str = "auto",
    runtime_evaluator_kwargs: dict[str, Any] | None = None,
) -> None:
    reference_points = list(getattr(run, "probe_points"))
    first_phase_space_point = getattr(run, "first_phase_space_point")
    if not reference_points and first_phase_space_point is not None:
        matrix_element = getattr(first_phase_space_point, "matrix_element")
        if matrix_element is not None:
            reference_points = [first_phase_space_point]
    if not reference_points:
        return

    from .evaluation import NativeRuntimeEvaluator

    evaluator = NativeRuntimeEvaluator(
        process,
        runtime_backend=cast(RuntimeBackend, runtime_backend),
        allow_reference_legacy=True,
        **(runtime_evaluator_kwargs or {}),
    )
    native_points: list[dict[str, int | float]] = []
    failures: list[str] = []
    timing_component_fields = (
        "source_fill_time_s",
        "momentum_setup_time_s",
        "parameter_pack_time_s",
        "evaluator_time_s",
        "output_transfer_time_s",
        "result_reduction_time_s",
    )
    timing_components = {field: 0.0 for field in timing_component_fields}
    timing_seen = False
    for point in reference_points:
        try:
            start = time.perf_counter()
            native = evaluator.evaluate(particles=getattr(point, "particles"))
            elapsed_s = time.perf_counter() - start
        except NativeEvaluationError as exc:
            failures.append(str(exc))
            continue
        reference_me = float(getattr(point, "matrix_element"))
        denom = max(abs(reference_me), abs(native.matrix_element), 1.0e-300)
        native_points.append(
            {
                "point": getattr(point, "point", 1),
                "reference_matrix_element": reference_me,
                "native_matrix_element": native.matrix_element,
                "relative_difference": abs(native.matrix_element - reference_me) / denom,
                "native_runtime_s": elapsed_s,
            }
        )
        runtime_impl = getattr(evaluator, "_runtime", None)
        timing = getattr(runtime_impl, "last_runtime_timing", None)
        timing_to_json = getattr(timing, "to_json_dict", None)
        if callable(timing_to_json):
            timing_payload = timing_to_json()
            timing_seen = True
            for field in timing_component_fields:
                timing_components[field] += float(timing_payload.get(field, 0.0))

    runtime_metadata = evaluator.refresh_metadata().to_json_dict()
    payload["native_runtime_backend"] = runtime_metadata

    if native_points:
        first = native_points[0]
        payload["native_first_point_available"] = True
        payload["native_first_point_matrix_element"] = first["native_matrix_element"]
        payload["native_first_point_relative_difference"] = first["relative_difference"]
    elif failures:
        payload["native_first_point_available"] = False
        payload["native_first_point_reason"] = failures[0]

    payload["native_probe_points"] = native_points
    payload["native_probe_summary"] = {
        "reference_points": len(reference_points),
        "available_points": len(native_points),
        "max_relative_difference": (
            max(float(point["relative_difference"]) for point in native_points)
            if native_points
            else None
        ),
        "failures": failures,
    }
    if native_points:
        runtimes = [float(point["native_runtime_s"]) for point in native_points]
        total_s = sum(runtimes)
        timing_breakdown: dict[str, object] | None = None
        if timing_seen:
            python_overhead = (
                timing_components["source_fill_time_s"]
                + timing_components["momentum_setup_time_s"]
                + timing_components["parameter_pack_time_s"]
                + timing_components["output_transfer_time_s"]
                + timing_components["result_reduction_time_s"]
            )
            timing_breakdown = {
                **timing_components,
                "python_overhead_time_s": python_overhead,
                "measured_time_s": python_overhead
                + timing_components["evaluator_time_s"],
            }
        evaluator_only_total_s = (
            _float_or_none(timing_breakdown.get("evaluator_time_s"))
            if isinstance(timing_breakdown, dict)
            else None
        )
        payload["native_runtime"] = {
            "evaluated_points": len(runtimes),
            "total_s": total_s,
            "mean_per_point_s": total_s / len(runtimes),
            "min_per_point_s": min(runtimes),
            "max_per_point_s": max(runtimes),
            "setup_time_s": evaluator.metadata.setup_time_s,
            "timing_breakdown_total_s": timing_breakdown,
            "evaluator_only_total_s": evaluator_only_total_s,
            "evaluator_only_mean_per_point_s": (
                evaluator_only_total_s / len(runtimes)
                if evaluator_only_total_s is not None
                else None
            ),
            "evaluations_per_second": (
                len(runtimes) / total_s if total_s > 0.0 else None
            ),
        }


def _attach_rusticol_probe_comparison(
    run: object,
    payload: dict[str, object],
    process_dir: Path,
    *,
    options: ProcessOptions,
    runtime_evaluator_kwargs: dict[str, Any],
    compare_args: argparse.Namespace | None = None,
) -> None:
    manifest_path = process_dir / "process_manifest.json"
    external_pdg_order: tuple[int, ...] | None = None
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        raw_order = manifest.get("external_pdg_order")
        if isinstance(raw_order, list):
            external_pdg_order = tuple(int(pdg) for pdg in raw_order)

    concrete_process = _concrete_probe_process_if_needed(
        run,
        external_pdg_order=external_pdg_order,
    )
    reference_color_order = _amplicol_reference_color_order_for_run(run)
    if concrete_process is not None:
        generation, process_dir = _rusticol_generation_profile_subprocess(
            concrete_process,
            process_dir.parent.parent,
            options,
            runtime_evaluator_kwargs=runtime_evaluator_kwargs,
            reference_color_order=reference_color_order,
            generic_dag_pruning_kwargs=(
                None
                if compare_args is None
                else _compare_generic_dag_pruning_kwargs(
                    compare_args,
                    process=concrete_process,
                )
            ),
        )
        payload["native_probe_generation"] = generation
        payload["native_probe_process"] = concrete_process
        manifest_path = process_dir / "process_manifest.json"
        external_pdg_order = None
        if manifest_path.exists():
            manifest = json.loads(manifest_path.read_text())
            raw_order = manifest.get("external_pdg_order")
            if isinstance(raw_order, list):
                external_pdg_order = tuple(int(pdg) for pdg in raw_order)

    import rusticol  # type: ignore[import-not-found]

    runtime = getattr(rusticol, "Runtime").load(str(process_dir))
    _attach_rusticol_probe_comparison_from_runtime(
        run,
        payload,
        runtime,
        external_pdg_order=external_pdg_order,
    )


def _rusticol_generation_profile_subprocess(
    process: str,
    cache_dir: Path,
    options: ProcessOptions,
    *,
    runtime_evaluator_kwargs: dict[str, Any],
    reference_color_order: Sequence[int] | None = None,
    generic_dag_pruning_kwargs: Mapping[str, Any] | None = None,
) -> tuple[dict[str, object], Path]:
    """Build a probe artifact in a fresh Python process.

    Symbolica's current Python bindings are not robust to building a second
    independent JIT artifact in the same interpreter after the first compare
    artifact has been built.  Probe-crossing fallbacks therefore keep the exact
    same generic generation function but run it in a short-lived interpreter.
    """

    code = r"""
import json
import os
from pathlib import Path

from pyamplicol.__main__ import _rusticol_generation_profile
from pyamplicol.processes import ProcessOptions

generation, process_dir = _rusticol_generation_profile(
    os.environ["PYAMPLICOL_PROBE_PROCESS"],
    Path(os.environ["PYAMPLICOL_PROBE_CACHE_DIR"]),
    ProcessOptions(**json.loads(os.environ["PYAMPLICOL_PROBE_OPTIONS"])),
    runtime_evaluator_kwargs=json.loads(os.environ["PYAMPLICOL_PROBE_RUNTIME"]),
    reference_color_order=json.loads(os.environ["PYAMPLICOL_REFERENCE_COLOR_ORDER"]),
    generic_dag_pruning_kwargs=json.loads(
        os.environ["PYAMPLICOL_PROBE_GENERIC_DAG_PRUNING"]
    ),
)
print(json.dumps({"generation": generation, "process_dir": str(process_dir)}, sort_keys=True))
"""
    env = _child_generation_environment()
    env.update(
        {
            "PYAMPLICOL_PROBE_PROCESS": process,
            "PYAMPLICOL_PROBE_CACHE_DIR": str(Path(cache_dir).expanduser()),
            "PYAMPLICOL_PROBE_OPTIONS": json.dumps(
                {
                    "flavour_scheme": options.flavour_scheme,
                    "include_3qqbar": options.include_3qqbar,
                    "include_cc": options.include_cc,
                    "include_resonance": options.include_resonance,
                    "serial": options.serial,
                },
                sort_keys=True,
            ),
            "PYAMPLICOL_PROBE_RUNTIME": json.dumps(
                _json_ready_runtime_kwargs(runtime_evaluator_kwargs),
                sort_keys=True,
            ),
            "PYAMPLICOL_REFERENCE_COLOR_ORDER": json.dumps(
                None
                if reference_color_order is None
                else [int(label) for label in reference_color_order],
                sort_keys=True,
            ),
            "PYAMPLICOL_PROBE_GENERIC_DAG_PRUNING": json.dumps(
                _json_ready_generic_dag_pruning_kwargs(
                    generic_dag_pruning_kwargs or {}
                ),
                sort_keys=True,
            ),
        }
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path.cwd(),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    payload = _parse_child_generation_payload(result.stdout)
    if result.returncode != 0:
        message = (
            result.stderr.strip()
            or result.stdout.strip()
            or f"probe subprocess generation exited with {result.returncode}"
        )
        raise RuntimeError(message)
    generation = payload.get("generation")
    process_dir = payload.get("process_dir")
    if not isinstance(generation, dict) or not isinstance(process_dir, str):
        raise RuntimeError(
            "probe subprocess generation did not return a valid artifact payload"
        )
    if result.stderr.strip():
        generation["subprocess_stderr"] = result.stderr.strip()
    return generation, Path(process_dir)


def _json_ready_runtime_kwargs(values: dict[str, Any]) -> dict[str, object]:
    converted: dict[str, object] = {}
    for key, value in values.items():
        if isinstance(value, Path):
            converted[key] = str(value)
        elif isinstance(value, tuple):
            converted[key] = list(value)
        elif isinstance(value, (str, int, float, bool)) or value is None:
            converted[key] = value
        elif isinstance(value, list):
            converted[key] = [
                str(item) if isinstance(item, Path) else item for item in value
            ]
        else:
            converted[key] = str(value)
    return converted


def _json_ready_generic_dag_pruning_kwargs(
    values: Mapping[str, Any],
) -> dict[str, object]:
    converted: dict[str, object] = {}
    for key, value in values.items():
        if isinstance(value, set):
            converted[key] = sorted(int(item) for item in value)
        elif isinstance(value, tuple):
            converted[key] = list(value)
        elif isinstance(value, list):
            converted[key] = value
        elif isinstance(value, Mapping):
            converted[key] = {
                str(item_key): item_value
                for item_key, item_value in value.items()
            }
        elif isinstance(value, (str, int, float, bool)) or value is None:
            converted[key] = value
        else:
            converted[key] = str(value)
    return converted


def _concrete_probe_process_if_needed(
    run: object,
    *,
    external_pdg_order: Sequence[int] | None,
) -> str | None:
    if external_pdg_order is None or len(external_pdg_order) < 2:
        return None
    reference_points = list(getattr(run, "probe_points"))
    first_phase_space_point = getattr(run, "first_phase_space_point", None)
    if not reference_points and first_phase_space_point is not None:
        matrix_element = getattr(first_phase_space_point, "matrix_element", None)
        if matrix_element is not None:
            reference_points = [first_phase_space_point]
    if not reference_points:
        return None
    particles = tuple(getattr(reference_points[0], "particles"))
    probe_order = tuple(int(getattr(particle, "pdg")) for particle in particles)
    if len(probe_order) != len(external_pdg_order):
        return None
    if probe_order[:2] == tuple(int(pdg) for pdg in external_pdg_order[:2]):
        return None
    try:
        initial = " ".join(_particle_name_from_pdg(pdg) for pdg in probe_order[:2])
        final = " ".join(_particle_name_from_pdg(pdg) for pdg in probe_order[2:])
    except KeyError:
        return None
    return f"{initial} > {final}"


def _particle_name_from_pdg(pdg: int) -> str:
    from .processes import PDGS

    names = {int(value): name for name, value in PDGS.items()}
    return names[int(pdg)]


def _attach_rusticol_probe_comparison_from_runtime(
    run: object,
    payload: dict[str, object],
    runtime: object,
    *,
    external_pdg_order: Sequence[int] | None = None,
) -> None:
    import numpy as np

    reference_points = list(getattr(run, "probe_points"))
    first_phase_space_point = getattr(run, "first_phase_space_point", None)
    if not reference_points and first_phase_space_point is not None:
        matrix_element = getattr(first_phase_space_point, "matrix_element", None)
        if matrix_element is not None:
            reference_points = [first_phase_space_point]
    if not reference_points:
        return

    native_points: list[dict[str, int | float]] = []
    failures: list[str] = []
    runtime_total_s = 0.0
    for point in reference_points:
        particles = tuple(getattr(point, "particles"))
        if external_pdg_order is not None:
            particles = _reorder_probe_particles_by_pdg(
                particles,
                external_pdg_order,
            )
        momenta = np.asarray(
            [
                [
                    [float(component) for component in particle.momentum]
                    for particle in particles
                ]
            ],
            dtype=np.float64,
        )
        try:
            start = time.perf_counter()
            values = runtime.evaluate(momenta)  # type: ignore[attr-defined]
            elapsed_s = time.perf_counter() - start
        except Exception as exc:  # pragma: no cover - surfaced in CLI payload
            failures.append(str(exc))
            continue
        pyamplicol_me = float(values[0])
        runtime_total_s += elapsed_s
        reference_me = float(getattr(point, "matrix_element"))
        denom = max(abs(reference_me), abs(pyamplicol_me), 1.0e-300)
        native_points.append(
            {
                "point": getattr(point, "point", 1),
                "reference_matrix_element": reference_me,
                "native_matrix_element": pyamplicol_me,
                "relative_difference": abs(pyamplicol_me - reference_me) / denom,
                "native_runtime_s": elapsed_s,
            }
        )

    metadata = dict(runtime.metadata())  # type: ignore[attr-defined]
    metadata["backend"] = "rusticol-generic-schema-v2"
    payload["native_runtime_backend"] = metadata

    if native_points:
        first = native_points[0]
        payload["native_first_point_available"] = True
        payload["native_first_point_matrix_element"] = first["native_matrix_element"]
        payload["native_first_point_relative_difference"] = first["relative_difference"]
    elif failures:
        payload["native_first_point_available"] = False
        payload["native_first_point_reason"] = failures[0]

    payload["native_probe_points"] = native_points
    payload["native_probe_summary"] = {
        "reference_points": len(reference_points),
        "available_points": len(native_points),
        "max_relative_difference": (
            max(float(point["relative_difference"]) for point in native_points)
            if native_points
            else None
        ),
        "failures": failures,
    }
    payload["native_runtime"] = {
        "total_s": runtime_total_s,
        "mean_per_point_s": runtime_total_s / max(len(native_points), 1),
    }


def _reorder_probe_particles_by_pdg(
    particles: Sequence[object],
    external_pdg_order: Sequence[int],
) -> tuple[object, ...]:
    if len(particles) != len(external_pdg_order):
        return tuple(particles)
    unused = list(range(len(particles)))
    reordered: list[object] = []
    for pdg in external_pdg_order:
        match_index = next(
            (
                index
                for index in unused
                if int(getattr(particles[index], "pdg")) == int(pdg)
            ),
            None,
        )
        if match_index is None:
            return tuple(particles)
        unused.remove(match_index)
        reordered.append(particles[match_index])
    return tuple(reordered)


def _command_summary(command: object) -> dict[str, object]:
    args = getattr(command, "args")
    return {
        "args": list(args),
        "cwd": str(getattr(command, "cwd")),
        "returncode": getattr(command, "returncode"),
        "elapsed_s": getattr(command, "elapsed_s"),
    }


def _default_amplicol_root() -> Path:
    cwd = Path.cwd()
    if (cwd / "amplicol_generate.f03").exists():
        return cwd
    if (cwd.parent / "amplicol_generate.f03").exists():
        return cwd.parent
    return Path(__file__).resolve().parents[3]


def _default_cache_dir() -> Path:
    return Path(".pyamplicol-cache")


def _default_mg5_path() -> Path | None:
    env_path = os.environ.get("MG5_PATH")
    if env_path:
        return Path(env_path)
    candidate = Path("/Users/vjhirsch/MG5/MG5_aMC_v3_6_0")
    return candidate if candidate.exists() else None


if __name__ == "__main__":
    raise SystemExit(main())
