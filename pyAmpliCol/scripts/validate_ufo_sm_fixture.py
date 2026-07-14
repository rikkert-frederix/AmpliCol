#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import deque
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from pyamplicol.model_assets import bundled_model_path


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixtures" / "result_matrix_references.json"
MATRIX_CACHES = {
    "lc": ROOT / "docs" / "result_matrix_data.json",
    "nlc": ROOT / "docs" / "result_matrix_nlc_data.json",
    "full": ROOT / "docs" / "result_matrix_full_data.json",
}
QCD_ONLY_FAMILIES = {
    "gg_tt_jets",
    "dd_tt_jets",
    "gg_gluons",
    # The archived multi-quark AmpliCol rows contain only strong vertices.
    # An unrestricted UFO SM also admits charged-current exchange, so matching
    # these benchmark workloads requires the same explicit coupling-order cap.
    "dd_3q_lines",
    "dd_4q_lines",
}
FAMILY_FLAGS = {
    "dd_3q_lines": ("--include-3qqbar", "--max-quark-pairs", "3"),
    "dd_4q_lines": ("--include-3qqbar", "--max-quark-pairs", "4"),
}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate UFO-SM artifacts for the archived result-matrix fixture "
            "and compare them at the identical stored phase-space points."
        )
    )
    parser.add_argument("--min-n", type=int, default=1)
    parser.add_argument("--max-n", type=int, default=4)
    parser.add_argument(
        "--color-accuracy",
        choices=("lc", "nlc", "full", "all"),
        default="all",
    )
    parser.add_argument("--case-id", action="append", default=[])
    parser.add_argument(
        "--output-root",
        type=Path,
        default=ROOT / ".ufo_support_outputs" / "ufo_sm_fixture_validation",
    )
    parser.add_argument(
        "--report-output",
        type=Path,
        help="Also write the cumulative portable report to this path.",
    )
    parser.add_argument("--n-cores", type=int, default=1)
    parser.add_argument("--tolerance", type=float, default=1.0e-10)
    parser.add_argument(
        "--regenerate",
        action="store_true",
        help="Archive an existing case artifact before regenerating it.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue after a generation or numerical failure.",
    )
    args = parser.parse_args(argv)

    if args.min_n < 1 or args.max_n < args.min_n:
        parser.error("require 1 <= --min-n <= --max-n")
    if args.n_cores < 1:
        parser.error("--n-cores must be positive")
    fixture = _read_json(FIXTURE)
    colors = (
        {"lc", "nlc", "full"}
        if args.color_accuracy == "all"
        else {args.color_accuracy}
    )
    selected_ids = set(args.case_id)
    cases = [
        case
        for case in fixture["cases"]
        if args.min_n <= int(case["n_final"]) <= args.max_n
        and str(case["color_accuracy"]) in colors
        and (not selected_ids or str(case["id"]) in selected_ids)
    ]
    cases.sort(
        key=lambda case: (
            int(case["n_final"]),
            str(case["color_accuracy"]),
            int(case["process_id"]),
        )
    )
    if selected_ids - {str(case["id"]) for case in cases}:
        missing = ", ".join(sorted(selected_ids - {str(case["id"]) for case in cases}))
        parser.error(f"unknown or filtered case ids: {missing}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    parameter_card = args.output_root / "matched-built-in-parameters.json"
    parameter_card.write_text(
        json.dumps(_matched_builtin_parameter_overrides(), indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    cache_payloads = {
        color: _read_json(path) for color, path in MATRIX_CACHES.items()
    }
    results: list[dict[str, Any]] = []
    failed = False
    for index, case in enumerate(cases, start=1):
        print(f"[{index}/{len(cases)}] {case['id']}: {case['process']}", flush=True)
        try:
            result = _validate_case(
                case,
                cache_payloads=cache_payloads,
                output_root=args.output_root,
                n_cores=args.n_cores,
                tolerance=args.tolerance,
                regenerate=args.regenerate,
                parameter_card=parameter_card,
            )
        except Exception as exc:
            result = {
                "id": str(case["id"]),
                "process": str(case["process"]),
                "color_accuracy": str(case["color_accuracy"]),
                "n_final": int(case["n_final"]),
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
            }
        results.append(result)
        _write_report(args.output_root / "report.json", results, args)
        if args.report_output is not None:
            _write_report(args.report_output, results, args)
        print(
            f"  {result['status']}"
            + (
                f" rel={float(result['relative_difference']):.3e}"
                if "relative_difference" in result
                else f" {result.get('error', '')}"
            ),
            flush=True,
        )
        if result["status"] != "ok":
            failed = True
            if not args.keep_going:
                break
    return 1 if failed else 0


def _validate_case(
    case: dict[str, Any],
    *,
    cache_payloads: dict[str, dict[str, Any]],
    output_root: Path,
    n_cores: int,
    tolerance: float,
    regenerate: bool,
    parameter_card: Path,
) -> dict[str, Any]:
    import rusticol

    case_id = str(case["id"])
    process_output = output_root / "processes" / case_id.replace(":", "_")
    generation_log = output_root / "logs" / f"{case_id.replace(':', '_')}.generate.log"
    if regenerate and process_output.exists():
        archive = _archive_process_output(output_root, process_output)
        print(f"    archived existing process output at {archive}", flush=True)
    elif process_output.exists() and not (
        process_output / "process_manifest.json"
    ).is_file():
        archive = _archive_process_output(output_root, process_output)
        print(f"    archived incomplete process output at {archive}", flush=True)
    if not (process_output / "process_manifest.json").is_file():
        command = _generation_command(
            case,
            cache_payloads=cache_payloads,
            process_output=process_output,
            model_cache=output_root / "model-cache",
            n_cores=n_cores,
        )
        environment = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
        _run_generation_command(
            command,
            cwd=ROOT,
            env=environment,
            log_path=generation_log,
        )
    point = np.asarray(
        [
            [
                [float(component) for component in particle["momentum"]]
                for particle in case["phase_space_point"]
            ]
        ],
        dtype=np.float64,
    )
    value = float(
        rusticol.Runtime.load(
            str(process_output),
            model_parameters=str(parameter_card),
        ).evaluate(point)[0]
    )
    reference = float(case["fortran"]["value"])
    relative_difference = abs(value - reference) / max(
        abs(value),
        abs(reference),
        1.0e-300,
    )
    return {
        "id": case_id,
        "process": str(case["process"]),
        "base_key": str(case["base_key"]),
        "color_accuracy": str(case["color_accuracy"]),
        "n_final": int(case["n_final"]),
        "status": "ok" if relative_difference <= tolerance else "mismatch",
        "value": value,
        "reference": reference,
        "relative_difference": relative_difference,
        "tolerance": tolerance,
        "process_output": os.path.relpath(process_output.resolve(), ROOT),
        "model_parameters": os.path.relpath(parameter_card.resolve(), ROOT),
        "generation_log": (
            os.path.relpath(generation_log.resolve(), ROOT)
            if generation_log.is_file()
            else None
        ),
    }


def _archive_process_output(output_root: Path, process_output: Path) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    archive = output_root / "archive" / stamp / process_output.name
    archive.parent.mkdir(parents=True, exist_ok=False)
    shutil.move(str(process_output), str(archive))
    return archive


def _run_generation_command(
    command: Sequence[str],
    *,
    cwd: Path,
    env: dict[str, str],
    log_path: Path,
) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    tail: deque[str] = deque(maxlen=80)
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"$ {shlex.join(command)}\n")
        log.flush()
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            log.write(line)
            log.flush()
            tail.append(line)
            if line.strip():
                print(f"    {line}", end="", flush=True)
        returncode = process.wait()
    if returncode:
        detail = "".join(tail).strip()
        raise RuntimeError(
            f"generation failed ({returncode}); log={log_path}"
            + (f": {detail}" if detail else "")
        )


def _generation_command(
    case: dict[str, Any],
    *,
    cache_payloads: dict[str, dict[str, Any]],
    process_output: Path,
    model_cache: Path,
    n_cores: int,
) -> list[str]:
    color = str(case["color_accuracy"])
    command = [
        sys.executable,
        "-m",
        "pyamplicol",
        "generate-process",
        "--model",
        str(bundled_model_path("sm", "json")),
        str(case["process"]),
        str(process_output),
        "--color-accuracy",
        color,
        "--symbolica-evaluator-backend",
        "jit",
        "--symbolica-jit-optimization-level",
        "1",
        "--symbolica-n-cores",
        str(n_cores),
        "--batch-size",
        "128" if color == "lc" else "64",
        "--symbolica-output-chunk-size",
        "128",
        "--symbolica-output-chunk-strategy",
        "auto" if color == "lc" else "uniform",
        "--model-cache-dir",
        str(model_cache),
        "--monitor",
    ]
    base_key = str(case["base_key"])
    command.extend(FAMILY_FLAGS.get(base_key, ()))
    if base_key in QCD_ONLY_FAMILIES:
        command.extend(("--max-coupling-order", "QED=0"))
    if color == "lc":
        command.extend(("--lc-sector-ids", "0"))
        entry = cache_payloads["lc"]["entries"][base_key][str(case["n_final"])]
        order = entry.get("amplicol", {}).get("reference_color_order")
        if isinstance(order, list) and order:
            command.extend(
                ("--reference-color-order", ",".join(str(value) for value in order))
            )
    return command


def _write_report(path: Path, results: list[dict[str, Any]], args: Any) -> None:
    payload = {
        "kind": "pyamplicol-ufo-sm-fixture-validation",
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "selection": {
            "min_n": int(args.min_n),
            "max_n": int(args.max_n),
            "color_accuracy": str(args.color_accuracy),
            "case_ids": list(args.case_id),
        },
        "summary": {
            "case_count": len(results),
            "ok": sum(result["status"] == "ok" for result in results),
            "failed": sum(result["status"] != "ok" for result in results),
            "max_relative_difference": max(
                (
                    float(result["relative_difference"])
                    for result in results
                    if "relative_difference" in result
                ),
                default=0.0,
            ),
        },
        "results": results,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _matched_builtin_parameter_overrides() -> dict[str, list[float]]:
    """Match the deliberately rounded constants in AmplicolSMLeadingColorModel."""

    default_alpha_ew = 1.0 / 132.507
    matched_alpha_ew = 0.007546771114
    # MW, and therefore sin(theta_W), depends on alpha_ew/Gf in the UFO SM.
    # Scale Gf with alpha_ew to retain the built-in model's unrounded weak angle.
    matched_gf = 1.16639e-5 * matched_alpha_ew / default_alpha_ew
    return {
        "aEWM1": [1.0 / matched_alpha_ew, 0.0],
        "Gf": [matched_gf, 0.0],
        "WH": [0.006382339, 0.0],
    }


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object in {path}")
    return payload


if __name__ == "__main__":
    raise SystemExit(main())
