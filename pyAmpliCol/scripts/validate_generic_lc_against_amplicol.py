#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import cast


if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pyamplicol.generic_validation import (  # noqa: E402
    DEFAULT_LC_VALIDATION_PROCESSES,
    format_validation_table,
    summary_to_json,
    validate_generic_lc_processes,
)
from pyamplicol.processes import ProcessOptions  # noqa: E402


def _default_amplicol_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_output_dir() -> Path:
    return Path(__file__).resolve().parents[1] / ".generic-lc-validation"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate generic leading-colour schema-v2 pyAmpliCol artifacts "
            "against Fortran AmpliCol library-backed supplied-momenta probes."
        )
    )
    parser.add_argument(
        "processes",
        nargs="*",
        help="Concrete processes to validate. Defaults to a broad LC smoke matrix.",
    )
    parser.add_argument("--output-dir", type=Path, default=_default_output_dir())
    parser.add_argument("--amplicol-root", type=Path, default=_default_amplicol_root())
    parser.add_argument("--rel-tol", type=float, default=1.0e-8)
    parser.add_argument(
        "--symbolica-evaluator-backend",
        choices=("jit", "compiled-complex", "compiled-complex-4x"),
        default="jit",
    )
    parser.add_argument("--symbolica-compiled-preset", default="runtime-o3")
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--n-cores", type=int, default=4)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--timeout", type=float)
    parser.add_argument(
        "--library-timing-events",
        type=int,
        default=0,
        help=(
            "When positive, also collect Fortran runtime timing through "
            "--library=create, make amplicol_generate_library, and --library=use."
        ),
    )
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument(
        "--reference-process-list-backend",
        choices=("legacy", "python"),
        default="legacy",
        help=(
            "Process-list writer used to steer Fortran AmpliCol. The default "
            "uses the original process_list.py; 'python' exercises the ported "
            "pyAmpliCol writer."
        ),
    )
    parser.add_argument(
        "--flavour-scheme",
        type=int,
        default=5,
        choices=range(1, 7),
        metavar="[1-6]",
        help="Massless QCD flavour scheme used for p/j expansion.",
    )
    parser.add_argument(
        "--include-3qqbar",
        action="store_true",
        help="Allow legacy three-quark-line subprocess enumeration.",
    )
    parser.add_argument(
        "--include-cc",
        action="store_true",
        help="Allow flavour-changing charged-current subprocesses.",
    )
    parser.add_argument(
        "--include-resonance",
        action="store_true",
        help="Use legacy resonance grouping for explicit lepton pairs.",
    )
    parser.add_argument(
        "--no-isolate-processes",
        dest="isolate_processes",
        action="store_false",
        help=(
            "Validate all requested processes in one Python interpreter. "
            "The default isolates real multi-process sweeps because repeated "
            "Symbolica evaluator generation can abort the interpreter."
        ),
    )
    parser.set_defaults(isolate_processes=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--json", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    using_default_matrix = not args.processes
    processes = tuple(args.processes) or DEFAULT_LC_VALIDATION_PROCESSES
    if using_default_matrix:
        args.include_cc = True
        args.include_3qqbar = True
    if (
        not args.isolate_processes
        and not args.dry_run
        and len(processes) > 1
        and args.symbolica_evaluator_backend == "jit"
    ):
        parser.error(
            "--no-isolate-processes cannot be combined with multiple JIT "
            "validations because repeated Symbolica JIT construction can abort "
            "the interpreter; use the default isolated mode or validate one "
            "process at a time."
        )
    if args.isolate_processes and not args.dry_run and len(processes) > 1:
        payload = _run_isolated(args, processes)
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            from pyamplicol.generic_validation import GenericLCValidationSummary
            from pyamplicol.generic_validation import GenericLCValidationRow

            payload_rows = cast(list[dict[str, object]], payload["rows"])
            rows = tuple(
                GenericLCValidationRow(
                    process=str(row["process"]),
                    key=str(row["key"]),
                    supported=bool(row["supported"]),
                    status=str(row["status"]),
                    concrete_process=(
                        None
                        if row.get("concrete_process") is None
                        else str(row["concrete_process"])
                    ),
                    relative_difference=_optional_float(
                        row.get("relative_difference")
                    ),
                    reference_matrix_element=_optional_float(
                        row.get("reference_matrix_element")
                    ),
                    pyamplicol_matrix_element=_optional_float(
                        row.get("pyamplicol_matrix_element")
                    ),
                )
                for row in payload_rows
            )
            print(
                format_validation_table(
                    GenericLCValidationSummary(rows=rows, rel_tol=args.rel_tol)
                )
            )
        return 0 if bool(payload["passed"]) else 1
    summary = validate_generic_lc_processes(
        processes,
        output_dir=args.output_dir,
        amplicol_root=args.amplicol_root,
        rel_tol=args.rel_tol,
        options=_process_options(args),
        evaluator_backend=args.symbolica_evaluator_backend,
        compiled_preset=args.symbolica_compiled_preset,
        batch_size=args.batch_size,
        n_cores=args.n_cores,
        jobs=args.jobs,
        timeout=args.timeout,
        dry_run=args.dry_run,
        library_timing_events=args.library_timing_events,
        seed=args.seed,
        reference_process_list_backend=args.reference_process_list_backend,
    )
    if args.json:
        print(summary_to_json(summary))
    else:
        print(format_validation_table(summary))
    return 0 if summary.passed else 1


def _run_isolated(args: argparse.Namespace, processes: tuple[str, ...]) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    script = Path(__file__).resolve()
    for process in processes:
        command = [
            sys.executable,
            str(script),
            "--output-dir",
            str(args.output_dir),
            "--amplicol-root",
            str(args.amplicol_root),
            "--rel-tol",
            str(args.rel_tol),
            "--symbolica-evaluator-backend",
            str(args.symbolica_evaluator_backend),
            "--symbolica-compiled-preset",
            str(args.symbolica_compiled_preset),
            "--batch-size",
            str(args.batch_size),
            "--n-cores",
            str(args.n_cores),
            "--jobs",
            str(args.jobs),
            "--library-timing-events",
            str(args.library_timing_events),
            "--seed",
            str(args.seed),
            "--reference-process-list-backend",
            str(args.reference_process_list_backend),
            "--flavour-scheme",
            str(args.flavour_scheme),
            "--no-isolate-processes",
            "--json",
        ]
        if args.include_3qqbar:
            command.append("--include-3qqbar")
        if args.include_cc:
            command.append("--include-cc")
        if args.include_resonance:
            command.append("--include-resonance")
        if args.timeout is not None:
            command.extend(["--timeout", str(args.timeout)])
        command.append(process)
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode != 0:
            failures.append(
                {
                    "process": process,
                    "returncode": result.returncode,
                    "stdout_tail": result.stdout[-4000:],
                    "stderr_tail": result.stderr[-4000:],
                }
            )
            continue
        try:
            payload = json.loads(result.stdout)
        except json.JSONDecodeError as error:
            failures.append(
                {
                    "process": process,
                    "returncode": result.returncode,
                    "error": str(error),
                    "stdout_tail": result.stdout[-4000:],
                    "stderr_tail": result.stderr[-4000:],
                }
            )
            continue
        child_rows = payload.get("rows")
        if not isinstance(child_rows, list):
            failures.append(
                {
                    "process": process,
                    "returncode": result.returncode,
                    "error": "child validation payload did not contain rows",
                    "stdout_tail": result.stdout[-4000:],
                    "stderr_tail": result.stderr[-4000:],
                }
            )
            continue
        rows.extend(row for row in child_rows if isinstance(row, dict))
    max_relative_difference = None
    rel_values = [
        rel_value
        for row in rows
        if (rel_value := _optional_float(row.get("relative_difference")))
        is not None
    ]
    if rel_values:
        max_relative_difference = max(rel_values)
    passed = not failures and len(rows) == len(processes) and all(
        bool(row.get("supported"))
        and row.get("status") == "ok"
        and (
            (rel_value := _optional_float(row.get("relative_difference")))
            is None
            or rel_value <= float(args.rel_tol)
        )
        for row in rows
    )
    return {
        "kind": "pyamplicol-generic-lc-validation-summary",
        "isolated_processes": True,
        "process_count": len(rows),
        "requested_process_count": len(processes),
        "passed": passed,
        "rel_tol": float(args.rel_tol),
        "max_relative_difference": max_relative_difference,
        "rows": rows,
        "failures": failures,
    }


def _process_options(args: argparse.Namespace) -> ProcessOptions:
    return ProcessOptions(
        flavour_scheme=int(args.flavour_scheme),
        include_3qqbar=bool(args.include_3qqbar),
        include_cc=bool(args.include_cc),
        include_resonance=bool(args.include_resonance),
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, (float, int)):
        return float(value)
    return None


if __name__ == "__main__":
    raise SystemExit(main())
