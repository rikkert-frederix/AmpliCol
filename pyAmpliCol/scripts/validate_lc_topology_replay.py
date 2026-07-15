#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pyamplicol.core_types import ExternalMomentum  # noqa: E402
from pyamplicol.color_plan import lc_topology_replay_safe_groups  # noqa: E402
from pyamplicol.symbolica_evaluator import SymbolicaEvaluatorSettings  # noqa: E402
from pyamplicol.generic_artifact import (  # noqa: E402
    GenericProcessManifest,
    build_generic_process_manifest,
    write_generic_dag_process_artifact,
)
from pyamplicol.phase_space import generic_validation_point  # noqa: E402
from pyamplicol.processes import ProcessOptions  # noqa: E402


@dataclass(frozen=True)
class ReplayValidationRow:
    process: str
    key: str
    status: str
    full_value: float | None = None
    replay_value: float | None = None
    relative_difference: float | None = None
    full_generation_s: float | None = None
    replay_generation_s: float | None = None
    replay_sector_count: int | None = None
    representative_sector_ids: tuple[int, ...] = ()
    error: str | None = None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "key": self.key,
            "status": self.status,
            "full_value": self.full_value,
            "replay_value": self.replay_value,
            "relative_difference": self.relative_difference,
            "full_generation_s": self.full_generation_s,
            "replay_generation_s": self.replay_generation_s,
            "replay_sector_count": self.replay_sector_count,
            "representative_sector_ids": list(self.representative_sector_ids),
            "error": self.error,
        }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate LC topology-representative replay by comparing a fully "
            "materialized generic DAG artifact with a representative-sector "
            "artifact evaluated through Rusticol replay."
        )
    )
    parser.add_argument("processes", nargs="+")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_default_output_dir(),
    )
    parser.add_argument(
        "--symbolica-evaluator-backend",
        choices=("jit", "compiled-complex", "compiled-complex-4x"),
        default="jit",
    )
    parser.add_argument("--symbolica-compiled-preset", default="runtime-o3")
    parser.add_argument("--batch-size", type=int, default=16)
    parser.add_argument("--n-cores", type=int, default=2)
    parser.add_argument("--seed", type=int, default=101)
    parser.add_argument("--rel-tol", type=float, default=1.0e-10)
    parser.add_argument(
        "--flavour-scheme",
        type=int,
        default=5,
        choices=range(1, 7),
        metavar="[1-6]",
    )
    parser.add_argument("--include-3qqbar", action="store_true")
    parser.add_argument("--include-cc", action="store_true")
    parser.add_argument("--include-resonance", action="store_true")
    parser.add_argument(
        "--no-isolate-processes",
        dest="isolate_processes",
        action="store_false",
        help=(
            "Validate all requested processes in one Python interpreter. "
            "The default isolates multi-process sweeps because repeated "
            "Symbolica evaluator generation can abort the interpreter."
        ),
    )
    parser.set_defaults(isolate_processes=True)
    parser.add_argument("--json", action="store_true")
    return parser


def _default_output_dir() -> Path:
    return Path(__file__).resolve().parents[1] / ".lc-topology-replay-validation"


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.isolate_processes and len(args.processes) > 1:
        payload = _run_isolated(args)
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            rows = [
                ReplayValidationRow(
                    process=str(row["process"]),
                    key=str(row["key"]),
                    status=str(row["status"]),
                    full_value=_optional_float(row.get("full_value")),
                    replay_value=_optional_float(row.get("replay_value")),
                    relative_difference=_optional_float(
                        row.get("relative_difference")
                    ),
                    full_generation_s=_optional_float(row.get("full_generation_s")),
                    replay_generation_s=_optional_float(
                        row.get("replay_generation_s")
                    ),
                    replay_sector_count=(
                        None
                        if row.get("replay_sector_count") is None
                        else int(row["replay_sector_count"])
                    ),
                    representative_sector_ids=tuple(
                        int(value)
                        for value in row.get("representative_sector_ids", [])
                    ),
                    error=(
                        None
                        if row.get("error") is None
                        else str(row["error"])
                    ),
                )
                for row in payload["rows"]
                if isinstance(row, dict)
            ]
            print(_format_table(rows))
            print(
                "passed="
                + str(payload["passed"])
                + ", max_relative_difference="
                + _format_float(payload["max_relative_difference"])
            )
        return 0 if bool(payload["passed"]) else 1

    options = ProcessOptions(
        flavour_scheme=args.flavour_scheme,
        include_3qqbar=args.include_3qqbar,
        include_cc=args.include_cc,
        include_resonance=args.include_resonance,
    )
    rows = [
        validate_process(
            process,
            output_dir=args.output_dir,
            options=options,
            evaluator_backend=args.symbolica_evaluator_backend,
            compiled_preset=args.symbolica_compiled_preset,
            batch_size=args.batch_size,
            n_cores=args.n_cores,
            seed=args.seed,
        )
        for process in args.processes
    ]
    payload = {
        "kind": "pyamplicol-lc-topology-replay-validation",
        "passed": all(
            row.status == "not-applicable"
            or (
                row.status == "ok"
                and row.relative_difference is not None
                and row.relative_difference <= args.rel_tol
            )
            for row in rows
        ),
        "rel_tol": args.rel_tol,
        "max_relative_difference": max(
            (
                row.relative_difference
                for row in rows
                if row.relative_difference is not None
            ),
            default=None,
        ),
        "rows": [row.to_json_dict() for row in rows],
    }
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(_format_table(rows))
        print(
            "passed="
            + str(payload["passed"])
            + ", max_relative_difference="
            + _format_float(payload["max_relative_difference"])
        )
    return 0 if payload["passed"] else 1


def _run_isolated(args: argparse.Namespace) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    script = Path(__file__).resolve()
    for process in args.processes:
        command = [
            sys.executable,
            str(script),
            "--output-dir",
            str(args.output_dir),
            "--symbolica-evaluator-backend",
            str(args.symbolica_evaluator_backend),
            "--symbolica-compiled-preset",
            str(args.symbolica_compiled_preset),
            "--batch-size",
            str(args.batch_size),
            "--n-cores",
            str(args.n_cores),
            "--seed",
            str(args.seed),
            "--rel-tol",
            str(args.rel_tol),
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
    max_relative_difference = max(
        (
            float(row["relative_difference"])
            for row in rows
            if isinstance(row.get("relative_difference"), (float, int))
        ),
        default=None,
    )
    passed = not failures and len(rows) == len(args.processes) and all(
        row.get("status") == "not-applicable"
        or (
            row.get("status") == "ok"
            and row.get("relative_difference") is not None
            and float(row["relative_difference"]) <= float(args.rel_tol)
        )
        for row in rows
    )
    return {
        "kind": "pyamplicol-lc-topology-replay-validation",
        "isolated_processes": True,
        "requested_process_count": len(args.processes),
        "process_count": len(rows),
        "passed": passed,
        "rel_tol": args.rel_tol,
        "max_relative_difference": max_relative_difference,
        "rows": rows,
        "failures": failures,
    }


def validate_process(
    process: str,
    *,
    output_dir: Path,
    options: ProcessOptions,
    evaluator_backend: str,
    compiled_preset: str,
    batch_size: int,
    n_cores: int,
    seed: int,
) -> ReplayValidationRow:
    try:
        base_manifest = build_generic_process_manifest(process, options=options)
        representative_sector_ids = _representative_sector_ids(base_manifest)
        if representative_sector_ids is None:
            return ReplayValidationRow(
                process=process,
                key=base_manifest.key,
                status="not-applicable",
                error="process has no replay-safe reducible LC topology groups",
            )
        root = output_dir.expanduser() / base_manifest.key
        full_dir = root / "full"
        replay_dir = root / "topology-replay"
        for directory in (full_dir, replay_dir):
            if directory.exists():
                shutil.rmtree(directory)

        full_generation_s = _write_artifact(
            base_manifest,
            full_dir,
            options=options,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            n_cores=n_cores,
            selected_color_sector_ids=None,
            lc_topology_replay=False,
        )
        replay_generation_s = _write_artifact(
            base_manifest,
            replay_dir,
            options=options,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            n_cores=n_cores,
            selected_color_sector_ids=set(representative_sector_ids),
            lc_topology_replay=True,
        )
        point = _ordered_momenta(
            generic_validation_point(process, seed=seed),
            base_manifest.external_pdg_order,
        )
        full_value = _evaluate(full_dir, point)
        replay_value = _evaluate(replay_dir, point)
        return ReplayValidationRow(
            process=process,
            key=base_manifest.key,
            status="ok",
            full_value=full_value,
            replay_value=replay_value,
            relative_difference=_relative_difference(full_value, replay_value),
            full_generation_s=full_generation_s,
            replay_generation_s=replay_generation_s,
            replay_sector_count=_replay_sector_count(replay_dir),
            representative_sector_ids=representative_sector_ids,
        )
    except Exception as exc:  # pragma: no cover - script-level diagnostics
        return ReplayValidationRow(
            process=process,
            key="",
            status="failed",
            error=str(exc),
        )


def _write_artifact(
    manifest: GenericProcessManifest,
    output_dir: Path,
    *,
    options: ProcessOptions,
    evaluator_backend: str,
    compiled_preset: str,
    batch_size: int,
    n_cores: int,
    selected_color_sector_ids: set[int] | None,
    lc_topology_replay: bool,
) -> float:
    settings = SymbolicaEvaluatorSettings(
        backend=evaluator_backend,
        compiled_preset=compiled_preset,
        n_cores=n_cores,
        compiled_chunk_compile_workers=n_cores,
        compiled_output_dir=str(output_dir / "compiled"),
    )
    start = time.perf_counter()
    write_generic_dag_process_artifact(
        manifest,
        output_dir,
        options=options,
        evaluator_backend=evaluator_backend,
        compiled_preset=compiled_preset,
        batch_size=batch_size,
        emit_stage_evaluator_artifacts=True,
        symbolica_settings=settings,
        selected_color_sector_ids=selected_color_sector_ids,
        lc_topology_replay=lc_topology_replay,
    )
    return time.perf_counter() - start


def _representative_sector_ids(
    manifest: GenericProcessManifest,
) -> tuple[int, ...] | None:
    if manifest.color_plan.color_accuracy != "lc":
        return None
    groups = lc_topology_replay_safe_groups(manifest.color_plan)
    if not groups or manifest.color_plan.sector_count <= 1:
        return None
    representatives = tuple(
        sorted({int(group.representative_sector_id) for group in groups})
    )
    if len(representatives) >= manifest.color_plan.sector_count:
        return None
    return representatives


def _ordered_momenta(
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
            raise ValueError(f"could not order validation point as {list(expected_pdgs)}")
    return tuple(ordered)


def _evaluate(process_dir: Path, particles: Sequence[ExternalMomentum]) -> float:
    import rusticol  # type: ignore[import-not-found]

    momenta = np.asarray(
        [
            [
                [float(component) for component in particle.momentum]
                for particle in particles
            ]
        ],
        dtype=np.float64,
    )
    runtime = rusticol.Runtime.load(str(process_dir))
    values = runtime.evaluate(momenta)
    return float(values[0])


def _replay_sector_count(process_dir: Path) -> int:
    payload = json.loads((process_dir / "process_manifest.json").read_text())
    replay = payload.get("compiled", {}).get("lc_topology_replay", {})
    return int(replay.get("replayed_sector_count", 0))


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def _format_table(rows: Sequence[ReplayValidationRow]) -> str:
    table = [
        ("Process", "Status", "Rel. diff", "Full gen", "Replay gen", "Sectors"),
        ("-" * 24, "-" * 16, "-" * 12, "-" * 10, "-" * 10, "-" * 7),
    ]
    for row in rows:
        table.append(
            (
                row.process,
                row.status,
                _format_float(row.relative_difference),
                _format_float(row.full_generation_s),
                _format_float(row.replay_generation_s),
                "N/A" if row.replay_sector_count is None else str(row.replay_sector_count),
            )
        )
    widths = [max(len(str(row[index])) for row in table) for index in range(6)]
    return "\n".join(
        " | ".join(str(value).ljust(widths[index]) for index, value in enumerate(row))
        for row in table
    )


def _format_float(value: object) -> str:
    return "N/A" if value is None else f"{float(value):.6g}"


def _optional_float(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, (float, int)):
        return float(value)
    return float(str(value))


if __name__ == "__main__":
    raise SystemExit(main())
