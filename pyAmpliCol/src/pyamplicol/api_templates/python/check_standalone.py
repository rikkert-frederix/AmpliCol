#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import json
import math
import statistics
import sys
import time
from decimal import Decimal
from pathlib import Path
from typing import Any


def import_rusticol(root: Path, explicit_folder: str | None):
    candidates: list[Path] = []
    if explicit_folder:
        candidates.append(Path(explicit_folder).expanduser())
    candidates.extend([root, Path.cwd(), *root.parents])
    for seed in candidates:
        for candidate in (
            seed,
            *seed.glob("lib/python*/site-packages"),
            *seed.glob(".venv/lib/python*/site-packages"),
            *seed.glob("dependencies/.venv/lib/python*/site-packages"),
            *seed.glob("pyAmpliCol/dependencies/.venv/lib/python*/site-packages"),
        ):
            if candidate.exists():
                sys.path.insert(0, str(candidate))
    try:
        return importlib.import_module("rusticol")
    except ModuleNotFoundError as error:
        raise SystemExit(
            "could not import rusticol; activate the managed pyAmpliCol venv or "
            "pass --rusticol-folder"
        ) from error


def validation_point(path: Path, process_key: str, precision: int):
    if not path.exists():
        return None
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 2 or fields[0] != process_key:
            continue
        external_count = int(fields[1])
        components = fields[2:]
        if len(components) != external_count * 4:
            raise SystemExit(f"invalid validation row for {process_key!r}")
        convert = float if precision == 16 else Decimal
        return [
            [convert(value) for value in components[index : index + 4]]
            for index in range(0, len(components), 4)
        ]
    return None


def repeat_points(point, count: int):
    try:
        import numpy as np
    except ModuleNotFoundError:
        return [point for _ in range(count)]
    if point and isinstance(point[0][0], float):
        return np.repeat(np.asarray([point], dtype=np.float64), count, axis=0)
    return [point for _ in range(count)]


def flatten_values(values: Any) -> list[Any]:
    if hasattr(values, "reshape"):
        return values.reshape(-1).tolist()
    result: list[Any] = []

    def visit(value: Any) -> None:
        if isinstance(value, (list, tuple)):
            for item in value:
                visit(item)
        else:
            result.append(value)

    visit(values)
    return result


def scalar_for_json(value: Any, precision: int) -> float | str:
    return float(value) if precision == 16 else str(value)


def physics_payload(physics) -> dict[str, Any]:
    colors = [
        {"id": item.id, "kind": "lc-flow", "word": list(item.word)}
        for item in physics.color_flows
    ]
    colors.extend(
        {"id": item.id, "kind": "contracted", "word": []}
        for item in physics.contracted_color_components
    )
    return {
        "process": physics.process,
        "process_key": physics.process_key,
        "color_accuracy": physics.color_accuracy,
        "external_particles": [
            {
                "index": item.index,
                "pdg": item.pdg,
            }
            for item in physics.external_particles
        ],
        "helicities": [
            {"id": item.id, "helicities": list(item.helicities)}
            for item in physics.helicities
        ],
        "colors": colors,
    }


def profile(runtime, point, precision: int, target_runtime: float, batch_size: int):
    samples: list[float] = []
    elapsed = 0.0
    while len(samples) < 8 or elapsed < target_runtime:
        batch = repeat_points(point, batch_size)
        started = time.perf_counter()
        runtime.profile(batch, precision=precision)
        wall = time.perf_counter() - started
        elapsed += wall
        samples.append(wall / batch_size)
    return {
        "samples": len(samples) * batch_size,
        "wall_us_per_point": statistics.mean(samples) * 1.0e6,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Resolved Rusticol Python API example")
    parser.add_argument("--process")
    parser.add_argument("--model-parameters")
    parser.add_argument(
        "--set-parameter",
        nargs=3,
        action="append",
        default=[],
        metavar=("NAME", "REAL", "IMAG"),
    )
    parser.add_argument("--precision", type=int, default=16)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--target-runtime", type=float, default=10.0)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--rusticol-folder")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    rusticol = import_rusticol(root, args.rusticol_folder)
    runtime = rusticol.Runtime.load(
        str(root),
        process_key=args.process,
        model_parameters=args.model_parameters,
    )
    if args.set_parameter:
        runtime.set_model_parameters(
            {
                name: (float(real), float(imaginary))
                for name, real, imaginary in args.set_parameter
            }
        )
    physics = runtime.physics
    metadata = physics_payload(physics)
    point = validation_point(
        root / "API" / "validation_points.dat",
        physics.process_key,
        args.precision,
    )
    if point is None:
        payload = {
            "language": "python",
            "available": False,
            "diagnostic": "no bundled validation point is available",
            **metadata,
        }
        if args.json:
            print(json.dumps(payload, separators=(",", ":"), sort_keys=True))
        else:
            print(f"process: {physics.process}")
            print("no bundled validation point is available; metadata load succeeded")
        return 0

    batch = repeat_points(point, 1)
    if args.precision == 16:
        resolved = runtime.evaluate_resolved(batch)
        total = runtime.evaluate(batch)
    else:
        resolved = runtime.evaluate_resolved_with_prec(batch, args.precision)
        total = runtime.evaluate_with_prec(batch, args.precision)
    values = flatten_values(resolved.values)
    explicit_total = flatten_values(resolved.total())
    compatibility_total = flatten_values(total)
    payload = {
        "language": "python",
        "available": True,
        "precision": args.precision,
        **metadata,
        "shape": list(resolved.shape),
        "values": [scalar_for_json(value, args.precision) for value in values],
        "resolved_sum": [
            scalar_for_json(value, args.precision) for value in explicit_total
        ],
        "compatibility_total": [
            scalar_for_json(value, args.precision) for value in compatibility_total
        ],
    }
    if args.profile:
        payload["profile"] = profile(
            runtime,
            point,
            args.precision,
            max(args.target_runtime, 0.0),
            max(args.batch_size, 1),
        )
    if args.json:
        print(json.dumps(payload, separators=(",", ":"), sort_keys=True))
        return 0

    print(f"process: {physics.process} [{physics.process_key}]")
    print(f"resolved shape: {resolved.shape}")
    color_items = metadata["colors"]
    color_count = len(color_items)
    for helicity_index, helicity in enumerate(metadata["helicities"]):
        for color_index, color in enumerate(color_items):
            offset = helicity_index * color_count + color_index
            print(f"  {helicity['id']}  {color['id']}  {values[offset]}")
    print(f"explicit resolved sum: {explicit_total[0]}")
    print(f"compatibility total:   {compatibility_total[0]}")
    if not math.isclose(
        float(explicit_total[0]),
        float(compatibility_total[0]),
        rel_tol=1.0e-12,
        abs_tol=1.0e-15,
    ):
        raise SystemExit("resolved components do not reproduce the compatibility total")
    if args.profile:
        print(f"timing: {payload['profile']['wall_us_per_point']:.6g} us/point")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
