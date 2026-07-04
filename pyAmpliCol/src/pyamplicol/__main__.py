from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from pprint import pprint
from typing import Any, Sequence, cast

from . import __version__
from .benchmarks import (
    benchmark_z_gluon_modes,
    format_mode_benchmark_table,
    profile_z_gluon_dag_evaluator,
    profile_z_gluon_tensor_evaluator,
)
from .evaluation import NativeRuntimeEvaluator, RuntimeBackend
from .matrix import (
    EvaluatorArtifact,
    NativeMatrixElementGenerator,
    evaluator_artifact_path,
    load_evaluator_artifact,
)
from .native import NativeEvaluationError
from .processes import ProcessEnumeration, ProcessEnumerator, ProcessOptions
from .reference import AmplicolAdapter, CommandResult
from .symbolic import ZeroGluonSymbolicEvaluator

_RUNTIME_BACKENDS = ("auto", "python", "dag", "numeric-tensor-network")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="pyamplicol matrix-element generation and validation tooling."
    )
    parser.add_argument(
        "--version",
        action="store_true",
        help="Print the pyamplicol package version and exit.",
    )

    subparsers = parser.add_subparsers(dest="command")

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

    generate = subparsers.add_parser(
        "generate",
        help="Generate native matrix-element metadata/evaluator artifacts.",
    )
    _add_process_options(generate)
    generate.add_argument("process", metavar="PROCESS")
    generate.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    generate.add_argument("--no-cache", action="store_true")
    _add_evaluator_build_options(generate)
    generate.add_argument("--json", action="store_true")

    evaluate = subparsers.add_parser(
        "evaluate",
        help="Evaluate a deterministic point with the native backend.",
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
    _add_runtime_backend_option(compare)
    _add_evaluator_build_options(compare)
    compare.add_argument("--skip-build", action="store_true")
    compare.add_argument(
        "--amplicol-probe",
        action="store_true",
        help="Use direct AmpliCol probing without MadGraph instead of --me_test.",
    )
    compare.add_argument("--dry-run", action="store_true")
    compare.add_argument("--json", action="store_true")

    family = subparsers.add_parser(
        "validate-z-gluon-family",
        help=(
            "Validate d d~ -> Z + n gluons against legacy AmpliCol direct probes."
        ),
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
        help=(
            "Benchmark legacy AmpliCol, native Python, numerical tensor-network, "
            "and parametric tensor-network modes on d d~ -> Z + n gluons."
        ),
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
        help="Profile the parametric Z+gluon tensor-network evaluator hot path.",
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
        help="Profile the shared-current D-mode evaluator hot path.",
    )
    dag_profile.add_argument("process", metavar="PROCESS")
    dag_profile.add_argument("--sqrt-s", type=float, default=1000.0)
    dag_profile.add_argument("--points", type=int, default=16)
    dag_profile.add_argument("--repetitions", type=int, default=100)
    _add_evaluator_build_options(dag_profile)
    dag_profile.add_argument("--json", action="store_true")

    profile = subparsers.add_parser(
        "profile",
        help="Profile native generation stages and cache behavior.",
    )
    _add_process_options(profile)
    profile.add_argument("process", metavar="PROCESS")
    profile.add_argument("--cache-dir", type=Path, default=_default_cache_dir())
    _add_runtime_backend_option(profile)
    _add_evaluator_build_options(profile)
    profile.add_argument("--json", action="store_true")

    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.version:
        print(__version__)
        return 0

    if args.command in (None, "inspect"):
        return _cmd_inspect()
    if args.command == "processes":
        return _cmd_processes(args)
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


def _add_runtime_backend_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--runtime-backend",
        choices=_RUNTIME_BACKENDS,
        default="auto",
        help=(
            "Native runtime backend: auto prefers the spenso/Symbolica current DAG; "
            "numeric-tensor-network executes the factorized network with fully "
            "numeric tensors in hep_lib."
        ),
    )


def _add_evaluator_build_options(parser: argparse.ArgumentParser) -> None:
    parser.set_defaults(
        merge_evaluators_strategy=False,
        symbolica_direct_translation=True,
        symbolica_compiled_native=True,
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
            "Experimental D-mode lowering that mirrors legacy AmpliCol more "
            "closely by compiling separate vertex and current-combine "
            "evaluators for each current-size stage."
        ),
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
        default="jit",
        help=(
            "Symbolica evaluator backend for D-mode blocks. compiled-complex "
            "writes and compiles C++ complex evaluators; compiled-complex-4x "
            "uses Symbolica's SIMD complex backend."
        ),
    )
    parser.add_argument(
        "--symbolica-iterations",
        "--symbolica-n-iterations",
        dest="symbolica_iterations",
        type=int,
        default=1,
        help="Number of Horner-scheme optimization iterations.",
    )
    parser.add_argument(
        "--symbolica-cpe-iterations",
        "--symbolica-n-cpe-iterations",
        dest="symbolica_cpe_iterations",
        type=int,
        help="Number of common-pair elimination iterations; default is unbounded.",
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
        action="store_true",
        help="Translate Symbolica instructions directly to SymJIT IR.",
    )
    parser.add_argument(
        "--symbolica-jit-optimization-level",
        type=int,
        choices=(0, 1, 2, 3),
        default=3,
        help="SymJIT optimization level.",
    )
    parser.add_argument(
        "--symbolica-max-horner-scheme-variables",
        type=int,
        default=500,
        help="Maximum number of variables considered in a Horner scheme.",
    )
    parser.add_argument(
        "--symbolica-max-common-pair-cache-entries",
        type=int,
        default=1_000_000,
        help="Maximum common-pair cache entries during Symbolica optimization.",
    )
    parser.add_argument(
        "--symbolica-max-common-pair-distance",
        type=int,
        default=100,
        help="Maximum distance between common pairs before cache eviction.",
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
            "best currently measured default by gluon multiplicity, including "
            "output chunking when not set explicitly; manual respects the "
            "explicit inline-asm and optimization options; "
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
        "--symbolica-compiled-output-chunk-size",
        type=int,
        help=(
            "Compile Symbolica block outputs in chunks of this size. "
            "Use 1 as a conservative workaround for the current AArch64 "
            "multi-output inline-assembly exporter issue."
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
            "Save a reusable compiled shared-current DAG evaluator artifact "
            "to this directory after generation. Supported for the "
            "compiled-complex D-mode used by profiling."
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


def _runtime_backend(args: argparse.Namespace) -> RuntimeBackend:
    value = getattr(args, "runtime_backend", "auto")
    if value not in _RUNTIME_BACKENDS:
        raise ValueError(f"unknown runtime backend: {value}")
    return cast(RuntimeBackend, value)


def _runtime_evaluator_kwargs(args: argparse.Namespace) -> dict[str, Any]:
    return {
        "batch_size": int(getattr(args, "batch_size", 16)),
        "merge_evaluators_strategy": bool(
            getattr(args, "merge_evaluators_strategy", False)
        ),
        "split_vertex_current_stages": bool(
            getattr(args, "symbolica_split_vertex_current_stages", False)
        ),
        "verbose_evaluator_build": bool(
            getattr(args, "verbose_evaluator_build", False)
        ),
        "symbolica_evaluator_backend": str(
            getattr(args, "symbolica_evaluator_backend", "jit")
        ),
        "symbolica_iterations": int(getattr(args, "symbolica_iterations", 1)),
        "symbolica_cpe_iterations": getattr(args, "symbolica_cpe_iterations", None),
        "symbolica_n_cores": int(getattr(args, "symbolica_n_cores", 4)),
        "symbolica_direct_translation": bool(
            getattr(args, "symbolica_direct_translation", True)
        ),
        "symbolica_jit_direct_translation": bool(
            getattr(args, "symbolica_jit_direct_translation", False)
        ),
        "symbolica_jit_optimization_level": int(
            getattr(args, "symbolica_jit_optimization_level", 3)
        ),
        "symbolica_max_horner_scheme_variables": int(
            getattr(args, "symbolica_max_horner_scheme_variables", 500)
        ),
        "symbolica_max_common_pair_cache_entries": int(
            getattr(args, "symbolica_max_common_pair_cache_entries", 1_000_000)
        ),
        "symbolica_max_common_pair_distance": int(
            getattr(args, "symbolica_max_common_pair_distance", 100)
        ),
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
        "symbolica_compiled_output_chunk_size": getattr(
            args,
            "symbolica_compiled_output_chunk_size",
            None,
        ),
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
    }


def _process_options(args: argparse.Namespace) -> ProcessOptions:
    return ProcessOptions(
        flavour_scheme=args.flavour_scheme,
        include_3qqbar=args.include_3qqbar,
        include_cc=args.include_cc,
        include_resonance=args.include_resonance,
        serial=not args.parallel_process_enumeration,
    )


def _cmd_inspect() -> int:
    import symbolica
    import symbolica.community.idenso  # noqa: F401
    import symbolica.community.spenso  # noqa: F401
    import symbolica.community.vakint  # noqa: F401

    local_versions = getattr(symbolica, "LOCAL_VERSIONS", None)
    if not isinstance(local_versions, dict):
        raise RuntimeError("symbolica.LOCAL_VERSIONS is missing or is not a dict")

    expected = {"symbolica", "spenso", "idenso", "vakint"}
    missing = expected.difference(local_versions)
    if missing:
        raise RuntimeError(f"symbolica.LOCAL_VERSIONS is missing keys: {sorted(missing)}")

    print("Loaded pyamplicol with Symbolica community modules: idenso, spenso, vakint")
    print("LOCAL_VERSIONS:")
    pprint(local_versions)
    return 0


def _cmd_processes(args: argparse.Namespace) -> int:
    enumerator = ProcessEnumerator(_process_options(args))
    enumeration = enumerator.enumerate(args.process)
    if args.legacy_output is not None:
        enumerator.write_legacy_file(enumeration, args.legacy_output)

    payload = _process_enumeration_to_dict(enumeration)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            f"{args.process}: {payload['n_unique_processes']} unique processes, "
            f"{payload['n_groups']} phase-space groups, {payload['n_records']} records"
        )
        if args.legacy_output is not None:
            print(f"legacy export: {args.legacy_output}")
    return 0


def _cmd_generate(args: argparse.Namespace) -> int:
    generator = NativeMatrixElementGenerator(
        cache_dir=None if args.no_cache else args.cache_dir
    )
    result = generator.generate(
        args.process,
        options=_process_options(args),
        write_cache_metadata=not args.no_cache,
    )
    payload = result.to_json_dict()
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        status = "supported" if result.supported_native_target else "unsupported"
        print(f"{args.process}: native target {status}; backend={result.backend}")
        if result.cache_file is not None:
            print(f"metadata cache: {result.cache_file}")
        for note in result.notes:
            print(note)
    return 0


def _cmd_evaluate(args: argparse.Namespace) -> int:
    artifact = _load_optional_evaluator_artifact(args.cache_dir, args.process)
    try:
        evaluator = NativeRuntimeEvaluator(
            args.process,
            runtime_backend=_runtime_backend(args),
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

    runtime_metadata = evaluator.metadata.to_json_dict()
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
    fixed_probe = args.amplicol_probe and _should_use_fixed_amplicol_probe(
        args.process, options
    )
    if args.dry_run:
        payload = _amplicol_dry_run_payload(args, options, fixed_probe=fixed_probe)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0

    native_generation = _native_generation_profile(
        args.process,
        args.cache_dir,
        options,
    )
    adapter = AmplicolAdapter(
        args.amplicol_root,
        jobs=args.jobs,
        timeout=args.timeout,
    )
    commands: list[CommandResult] = []
    if not args.skip_build:
        if fixed_probe:
            build = adapter.prepare_direct_probe(
                args.process,
                process_file=args.process_file,
                options=options,
            )
        else:
            build = adapter.prepare_library(
                args.process,
                process_file=args.process_file,
                options=options,
            )
        commands.extend(build.commands)
        process_file = build.process_file
    else:
        process_file = args.process_file

    if args.amplicol_probe:
        if fixed_probe:
            run = adapter.run_amplicol_fixed_probe(
                args.process,
                points=args.points,
                process_file=process_file,
                options=options,
                timing_sample=args.timing,
            )
            mode = "amplicol_fixed_probe"
        else:
            run = adapter.run_amplicol_probe(
                args.process,
                points=args.points,
                process_file=process_file,
                options=options,
                timing_sample=args.timing,
            )
            mode = "amplicol_probe"
    else:
        run = adapter.run_me_test(
            args.process,
            points=args.points,
            process_file=process_file,
            options=options,
            mg5_path=args.mg5_path,
            timing_sample=args.timing,
        )
        mode = "me_test"
    commands.extend(run.commands)
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
    _attach_native_probe_comparison(
        args.process,
        run,
        payload,
        runtime_backend=_runtime_backend(args),
        runtime_evaluator_kwargs=_runtime_evaluator_kwargs(args),
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
        if gluon_count == 0:
            build = adapter.prepare_direct_probe(process, options=options)
            run = adapter.run_amplicol_fixed_probe(
                process,
                points=args.points,
                options=options,
                timing_sample=args.timing,
            )
            mode = "amplicol_fixed_probe"
        else:
            build = adapter.prepare_library(process, options=options)
            run = adapter.run_amplicol_probe(
                process,
                points=args.points,
                options=options,
                timing_sample=args.timing,
            )
            mode = "amplicol_probe"
        commands.extend(build.commands)
        commands.extend(run.commands)
        row: dict[str, object] = {
            "gluon_count": gluon_count,
            "process": process,
            "mode": mode,
            "runtime_backend": _z_gluon_family_runtime_backend(
                gluon_count,
                _runtime_backend(args),
            ),
            "process_file": str(run.process_file),
            "probe_points": [_probe_point_to_dict(point) for point in run.probe_points],
            "timing_rows": [
                {"label": timing.label, "seconds": timing.seconds, "note": timing.note}
                for timing in run.timing_rows
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
        evaluator_build_kwargs=_runtime_evaluator_kwargs(args),
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
    payload = profile_z_gluon_dag_evaluator(
        args.process,
        sqrt_s=args.sqrt_s,
        points=args.points,
        repetitions=args.repetitions,
        evaluator_build_kwargs=_runtime_evaluator_kwargs(args),
        save_evaluator_dir=getattr(args, "save_evaluator_dir", None),
    )
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            f"{args.process}: shared D-mode profile, "
            f"points={payload['points']}, repetitions={payload['repetitions']}, "
            f"batch_size={payload['batch_size']}"
        )
        print(f"generation: {payload['generation_s']:.6g} s")
        print(
            "full runtime: "
            f"{payload['runtime_us_per_point']:.6g} "
            f"+/- {_format_error(payload.get('runtime_us_per_point_error'))} us/point"
        )
        print(
            "evaluator-only runtime: "
            f"{payload['runtime_evaluator_only_us_per_point']:.6g} "
            "+/- "
            f"{_format_error(payload.get('runtime_evaluator_only_us_per_point_error'))} "
            "us/point"
        )
        breakdown = payload.get("runtime_breakdown_us_per_point")
        breakdown_errors = payload.get("runtime_breakdown_us_per_point_error")
        if isinstance(breakdown, dict):
            print("runtime breakdown:")
            for key in (
                "source_fill_time_s",
                "momentum_setup_time_s",
                "parameter_pack_time_s",
                "evaluator_time_s",
                "output_transfer_time_s",
                "result_reduction_time_s",
                "python_overhead_time_s",
            ):
                value = breakdown.get(key)
                if not isinstance(value, (float, int)):
                    continue
                error = (
                    breakdown_errors.get(key)
                    if isinstance(breakdown_errors, dict)
                    else None
                )
                print(
                    f"  {key.removesuffix('_time_s')}: "
                    f"{float(value):.6g} +/- {_format_error(error)} us/point"
                )
    return 0


def _format_error(value: object) -> str:
    if not isinstance(value, (float, int)):
        return "n/a"
    return f"{float(value):.3g}"


def _cmd_profile(args: argparse.Namespace) -> int:
    generator = NativeMatrixElementGenerator(cache_dir=args.cache_dir)
    result = generator.generate(args.process, options=_process_options(args))
    artifact = _load_optional_evaluator_artifact(args.cache_dir, args.process)
    per_point_runtime_s = None
    runtime_metadata = None
    try:
        evaluator = NativeRuntimeEvaluator(
            args.process,
            runtime_backend=_runtime_backend(args),
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
) -> EvaluatorArtifact | None:
    try:
        return load_evaluator_artifact(cache_dir, process)
    except (OSError, ValueError):
        return None


def _artifact_status(
    cache_dir: Path,
    process: str,
    artifact: EvaluatorArtifact | None,
) -> dict[str, object]:
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
    artifact: EvaluatorArtifact | None,
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
    artifact: EvaluatorArtifact,
) -> dict[str, object]:
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
    artifact: EvaluatorArtifact,
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
    artifact: EvaluatorArtifact | None,
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
    artifact: EvaluatorArtifact | None,
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
) -> dict[str, object]:
    root = Path(args.amplicol_root)
    process_file = args.process_file or root / "processes.txt"
    commands = []
    if not args.skip_build:
        commands.extend(
            [["make", "cleanlib"], ["make", f"-j{args.jobs}", "amplicol_generate"]]
        )
        if not fixed_probe:
            commands.extend(
                [
                    ["./amplicol_generate", "--library=create", f"--process={process_file}"],
                    ["make", f"-j{args.jobs}", "amplicol_generate_library"],
                ]
            )
    timing = args.points if args.timing is None else args.timing
    if args.amplicol_probe:
        probe_flag = "--amplicol_fixed_probe" if fixed_probe else "--amplicol_probe"
        commands.append(
            [
                "./amplicol_generate",
                f"{probe_flag}={args.points}",
                f"--timing={timing}",
                f"--process={process_file}",
            ]
        )
    else:
        commands.append(
            [
                "./amplicol_generate",
                f"--me_test={args.points}",
                f"--timing={timing}",
                f"--process={process_file}",
            ]
        )
    mode = "me_test"
    if args.amplicol_probe:
        mode = "amplicol_fixed_probe" if fixed_probe else "amplicol_probe"
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
        commands = []
        if gluon_count == 0:
            commands.extend(
                [["make", "cleanlib"], ["make", f"-j{args.jobs}", "amplicol_generate"]]
            )
            commands.append(
                [
                    "./amplicol_generate",
                    f"--amplicol_fixed_probe={args.points}",
                    f"--timing={args.points if args.timing is None else args.timing}",
                    "--process=processes.txt",
                ]
            )
            mode = "amplicol_fixed_probe"
        else:
            commands.extend(
                [
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
                        f"--amplicol_probe={args.points}",
                        f"--timing={args.points if args.timing is None else args.timing}",
                        "--process=processes.txt",
                    ],
                ]
            )
            mode = "amplicol_probe"
        rows.append(
            {
                "gluon_count": gluon_count,
                "process": process,
                "mode": mode,
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
    requested: RuntimeBackend,
) -> RuntimeBackend:
    if gluon_count == 0 and requested in ("dag", "numeric-tensor-network"):
        return "python"
    return requested


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
    runtime_backend: RuntimeBackend = "auto",
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

    evaluator = NativeRuntimeEvaluator(
        process,
        runtime_backend=runtime_backend,
        **(runtime_evaluator_kwargs or {}),
    )
    payload["native_runtime_backend"] = evaluator.metadata.to_json_dict()
    native_points: list[dict[str, int | float]] = []
    failures: list[str] = []
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
        payload["native_runtime"] = {
            "evaluated_points": len(runtimes),
            "total_s": total_s,
            "mean_per_point_s": total_s / len(runtimes),
            "min_per_point_s": min(runtimes),
            "max_per_point_s": max(runtimes),
            "setup_time_s": evaluator.metadata.setup_time_s,
            "evaluations_per_second": (
                len(runtimes) / total_s if total_s > 0.0 else None
            ),
        }


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
