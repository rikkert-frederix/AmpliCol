#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[2]
PYAMPLICOL_ROOT = REPO_ROOT / "pyAmpliCol"
DOCS_DIR = PYAMPLICOL_ROOT / "docs"
FIXTURE_PATH = PYAMPLICOL_ROOT / "tests" / "fixtures" / "result_matrix_references.json"

MATRIX_CACHES = {
    "lc": DOCS_DIR / "result_matrix_data.json",
    "nlc": DOCS_DIR / "result_matrix_nlc_data.json",
    "full": DOCS_DIR / "result_matrix_full_data.json",
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build the lightweight integration-test fixture from the documented "
            "LC/NLC/full result-matrix caches. With --refresh-caches, first "
            "rerun result_matrix.py for n<=4 so the fixture is regenerated from "
            "Fortran AmpliCol and pyAmpliCol validation runs."
        )
    )
    parser.add_argument("--output", type=Path, default=FIXTURE_PATH)
    parser.add_argument(
        "--max-n",
        type=int,
        default=4,
        help="Maximum final-state multiplicity archived in the fixture.",
    )
    parser.add_argument(
        "--color-accuracy",
        choices=("lc", "nlc", "full", "all"),
        default="all",
        help="Colour-accuracy cache(s) to extract.",
    )
    parser.add_argument(
        "--refresh-caches",
        action="store_true",
        help=(
            "Run docs/result_matrix.py before extracting fixtures. This requires "
            "a working Fortran AmpliCol build and pyAmpliCol runtime."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Fortran make jobs passed to result_matrix.py when refreshing.",
    )
    parser.add_argument(
        "--n-cores",
        type=int,
        default=4,
        help="Symbolica/Rusticol cores passed to result_matrix.py when refreshing.",
    )
    parser.add_argument(
        "--process-workers",
        type=int,
        default=1,
        help=(
            "Independent matrix-cell workers passed to result_matrix.py. Keep "
            "this at one when refreshing AmpliCol references."
        ),
    )
    parser.add_argument(
        "--skip-cpp-o3",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Skip C++ O3 rows when refreshing caches. The integration fixture "
            "only needs Fortran and JIT same-point values."
        ),
    )
    args = parser.parse_args(argv)

    colors = ("lc", "nlc", "full") if args.color_accuracy == "all" else (args.color_accuracy,)
    if args.refresh_caches:
        for color in colors:
            _refresh_cache(
                color,
                max_n=args.max_n,
                jobs=args.jobs,
                n_cores=args.n_cores,
                process_workers=args.process_workers,
                skip_cpp_o3=bool(args.skip_cpp_o3),
            )

    fixture = _build_fixture(colors, max_n=args.max_n)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(fixture, indent=2, sort_keys=True) + "\n")
    print(
        f"wrote {args.output} with {len(fixture['cases'])} validated cases, "
        f"{len(fixture['unsupported'])} structural unsupported cases, and "
        f"{len(fixture['unvalidated'])} unvalidated pyAmpliCol-supported cases"
    )
    return 0


def _refresh_cache(
    color: str,
    *,
    max_n: int,
    jobs: int,
    n_cores: int,
    process_workers: int,
    skip_cpp_o3: bool,
) -> None:
    cmd = [
        sys.executable,
        str(DOCS_DIR / "result_matrix.py"),
        "--color-accuracy",
        color,
        "--generate-data",
        f"1-{max_n}",
        "--show-data",
        f"1-{max_n}",
        "--validate",
        "--jobs",
        str(jobs),
        "--n-cores",
        str(n_cores),
        "--process-workers",
        str(process_workers),
    ]
    if skip_cpp_o3:
        cmd.append("--skip-cpp-o3")
    subprocess.run(cmd, cwd=REPO_ROOT, check=True)


def _build_fixture(colors: Iterable[str], *, max_n: int) -> dict[str, Any]:
    cases: list[dict[str, Any]] = []
    unsupported: list[dict[str, Any]] = []
    unvalidated: list[dict[str, Any]] = []
    sources: list[dict[str, Any]] = []
    for color in colors:
        path = MATRIX_CACHES[color]
        data = _read_json(path)
        process_ids = {
            str(base.get("key")): int(base.get("process_id", index))
            for index, base in enumerate(data.get("base_processes", []), start=1)
        }
        sources.append(
            {
                "color_accuracy": color,
                "path": _repo_relative(path),
                "updated_at": data.get("updated_at"),
                "schema_version": data.get("schema_version"),
            }
        )
        for base_key, by_n in sorted(data.get("entries", {}).items()):
            for n_text, entry in sorted(by_n.items(), key=lambda item: int(item[0])):
                n_final = int(n_text)
                if n_final > max_n:
                    continue
                validation = entry.get("validation") or {}
                amplicol = entry.get("amplicol") or {}
                jit = entry.get("pyamplicol_jit") or {}
                if validation.get("status") == "ok":
                    cases.append(
                        _case_from_entry(
                            color=color,
                            base_key=base_key,
                            process_id=process_ids.get(base_key),
                            n_final=n_final,
                            entry=entry,
                            validation=validation,
                            amplicol=amplicol,
                            jit=jit,
                            source_cache=path,
                        )
                    )
                elif _is_structural_unsupported(amplicol, jit):
                    payload = {
                        "id": f"{color}:{base_key}:n{n_final}",
                        "base_key": base_key,
                        "process_id": process_ids.get(base_key),
                        "n_final": n_final,
                        "process": entry.get("process"),
                        "color_accuracy": color,
                        "amplicol_status": amplicol.get("status"),
                        "pyamplicol_jit_status": jit.get("status"),
                        "reason": _structural_reason(amplicol, jit),
                        "source_cache": _repo_relative(path),
                    }
                    if _jit_supported_but_fortran_reference_missing(amplicol, jit):
                        payload["kind"] = "fortran-reference-unavailable"
                        payload["reason"] = _fortran_reference_missing_reason()
                        unvalidated.append(payload)
                    else:
                        payload["kind"] = "structural-unsupported"
                        unsupported.append(payload)

    cases.sort(key=lambda item: (item["color_accuracy"], item["process_id"], item["n_final"]))
    unsupported.sort(
        key=lambda item: (item["color_accuracy"], item.get("process_id") or 0, item["n_final"])
    )
    unvalidated.sort(
        key=lambda item: (item["color_accuracy"], item.get("process_id") or 0, item["n_final"])
    )
    summary = _summary(cases, unsupported, unvalidated)
    return {
        "schema_version": 1,
        "kind": "pyamplicol-result-matrix-reference-fixture",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_by": "pyAmpliCol/scripts/refresh_result_matrix_fixtures.py",
        "git": _git_metadata(),
        "max_n": max_n,
        "sources": sources,
        "summary": summary,
        "cases": cases,
        "unsupported": unsupported,
        "unvalidated": unvalidated,
    }


def _case_from_entry(
    *,
    color: str,
    base_key: str,
    process_id: int | None,
    n_final: int,
    entry: dict[str, Any],
    validation: dict[str, Any],
    amplicol: dict[str, Any],
    jit: dict[str, Any],
    source_cache: Path,
) -> dict[str, Any]:
    value = float(validation["values"]["pyamplicol_jit"])
    reference = float(validation["reference"])
    rel_diff = float(validation["relative_differences"]["pyamplicol_jit"])
    point = _read_validation_point(validation.get("point_source"))
    return {
        "id": f"{color}:{base_key}:n{n_final}",
        "base_key": base_key,
        "process_id": process_id,
        "n_final": n_final,
        "process": entry.get("process"),
        "color_accuracy": color,
        "phase_space_point": point,
        "fortran": {
            "value": reference,
            "mode": amplicol.get("mode"),
            "status": amplicol.get("status"),
            "reference_probe": amplicol.get("reference_probe"),
            "runtime_us_per_point": amplicol.get("runtime_us_per_point"),
        },
        "pyamplicol": {
            "mode": jit.get("mode"),
            "status": jit.get("status"),
            "value": value,
            "generation_s": jit.get("generation_s"),
            "wall_us_per_point": jit.get("wall_us_per_point"),
            "core_us_per_point": jit.get("runtime_us_per_point"),
        },
        "validation": {
            "status": validation.get("status"),
            "kind": validation.get("kind"),
            "tolerance": float(validation.get("tolerance", 1.0e-8)),
            "absolute_tolerance": float(validation.get("absolute_tolerance", 1.0e-16)),
            "relative_difference": rel_diff,
            "absolute_difference": abs(value - reference),
            "max_relative_difference": float(validation.get("max_relative_difference", rel_diff)),
            "point_order_source": validation.get("point_order_source"),
        },
        "metadata": {
            "source_cache": _repo_relative(source_cache),
            "point_source": _repo_relative(Path(str(validation.get("point_source", "")))),
            "amplicol_command": validation.get("amplicol_command"),
        },
    }


def _read_validation_point(path_value: object) -> list[dict[str, Any]]:
    if not isinstance(path_value, str) or not path_value:
        raise ValueError("validation record is missing point_source")
    path = Path(path_value)
    payload = _read_json(path)
    points = payload.get("points")
    if not isinstance(points, list) or not points:
        raise ValueError(f"validation point file has no points: {path}")
    point = points[0]
    if not isinstance(point, list) or not point:
        raise ValueError(f"validation point is malformed: {path}")
    return [
        {
            "pdg": int(entry["pdg"]),
            "momentum": [str(component) for component in entry["momentum"]],
        }
        for entry in point
    ]


def _summary(
    cases: list[dict[str, Any]],
    unsupported: list[dict[str, Any]],
    unvalidated: list[dict[str, Any]],
) -> dict[str, Any]:
    by_color: dict[str, dict[str, Any]] = {}
    for case in cases:
        bucket = by_color.setdefault(
            case["color_accuracy"],
            {
                "validated_cases": 0,
                "max_relative_difference": 0.0,
                "max_absolute_difference": 0.0,
            },
        )
        bucket["validated_cases"] += 1
        bucket["max_relative_difference"] = max(
            float(bucket["max_relative_difference"]),
            float(case["validation"]["relative_difference"]),
        )
        bucket["max_absolute_difference"] = max(
            float(bucket["max_absolute_difference"]),
            float(case["validation"]["absolute_difference"]),
        )
    for item in unsupported:
        bucket = by_color.setdefault(
            item["color_accuracy"],
            {
                "validated_cases": 0,
                "max_relative_difference": 0.0,
                "max_absolute_difference": 0.0,
            },
        )
        bucket["structural_unsupported_cases"] = (
            int(bucket.get("structural_unsupported_cases", 0)) + 1
        )
    for item in unvalidated:
        bucket = by_color.setdefault(
            item["color_accuracy"],
            {
                "validated_cases": 0,
                "max_relative_difference": 0.0,
                "max_absolute_difference": 0.0,
            },
        )
        bucket["unvalidated_pyamplicol_cases"] = (
            int(bucket.get("unvalidated_pyamplicol_cases", 0)) + 1
        )
    return {
        "validated_cases": len(cases),
        "structural_unsupported_cases": len(unsupported),
        "unvalidated_pyamplicol_cases": len(unvalidated),
        "by_color_accuracy": by_color,
    }


def _jit_supported_but_fortran_reference_missing(
    amplicol: dict[str, Any],
    jit: dict[str, Any],
) -> bool:
    if jit.get("status") != "ok":
        return False
    text = " ".join(
        str(value)
        for value in (
            amplicol.get("status"),
            amplicol.get("error"),
        )
    ).lower()
    return "more than two quarks" in text


def _fortran_reference_missing_reason() -> str:
    return (
        "Fortran AmpliCol init_col does not expose a direct NLC/full colour "
        "matrix for more than two quark lines; pyAmpliCol generated and "
        "evaluated the generic colour contraction, but this cell has no "
        "direct Fortran colour-contraction reference."
    )


def _is_structural_unsupported(amplicol: dict[str, Any], jit: dict[str, Any]) -> bool:
    text = " ".join(
        str(value)
        for value in (
            amplicol.get("status"),
            amplicol.get("error"),
            jit.get("status"),
            jit.get("error"),
        )
    ).lower()
    return (
        "unsupported" in text
        or "more than two quarks" in text
        or "one, or two quark pairs" in text
    )


def _structural_reason(amplicol: dict[str, Any], jit: dict[str, Any]) -> str:
    for payload in (jit, amplicol):
        error = str(payload.get("error", "")).strip()
        if error:
            for line in error.splitlines():
                if "quark" in line.lower() or "unsupported" in line.lower():
                    return line.strip()
            return error.splitlines()[0].strip()
    return "structurally unsupported"


def _git_metadata() -> dict[str, Any]:
    def run(*args: str) -> str | None:
        try:
            completed = subprocess.run(
                ["git", *args],
                cwd=REPO_ROOT,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
            )
            return completed.stdout.strip()
        except subprocess.CalledProcessError:
            return None

    status = run("status", "--short")
    return {
        "head": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(status),
    }


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _repo_relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


if __name__ == "__main__":
    raise SystemExit(main())
