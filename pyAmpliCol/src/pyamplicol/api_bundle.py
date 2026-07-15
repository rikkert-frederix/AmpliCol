from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any, Mapping, Sequence


_BUNDLE_FILES = (
    ("python/check_standalone.py", "python/check_standalone.py"),
    ("cpp/check_standalone.cpp", "cpp/check_standalone.cpp"),
    ("cpp/Makefile", "cpp/Makefile"),
    ("fortran/check_standalone.f90", "fortran/check_standalone.f90"),
    ("fortran/Makefile", "fortran/Makefile"),
)


def write_api_bundle(output_dir: str | Path) -> Path:
    """Write the root-only, multi-language API examples for a schema-v2 artifact."""

    root = Path(output_dir).expanduser()
    bundle = root / "API"
    template_root = Path(__file__).with_name("api_templates")
    for source_name, target_name in _BUNDLE_FILES:
        source = template_root / source_name
        target = bundle / target_name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
    (bundle / "python" / "check_standalone.py").chmod(0o755)
    _write_validation_points(root, bundle / "validation_points.dat")
    legacy_checker = root / "check_standalone.py"
    if legacy_checker.exists():
        legacy_checker.unlink()
    return bundle


def remove_api_bundle(output_dir: str | Path) -> None:
    """Remove examples from an internal artifact that must not own a bundle."""

    root = Path(output_dir).expanduser()
    bundle = root / "API"
    if bundle.exists():
        shutil.rmtree(bundle)
    legacy_checker = root / "check_standalone.py"
    if legacy_checker.exists():
        legacy_checker.unlink()


def _write_validation_points(root: Path, target: Path) -> None:
    rows, unavailable = _validation_rows(root)
    lines = ["RUSTICOL_VALIDATION_POINTS_V1"]
    lines.extend(
        "\t".join((key, str(external_count), *components))
        for key, external_count, components in rows
    )
    lines.extend(f"# unavailable\t{key}\t{message}" for key, message in unavailable)
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _validation_rows(
    root: Path,
) -> tuple[list[tuple[str, int, list[str]]], list[tuple[str, str]]]:
    process_set_path = root / "process_set_manifest.json"
    if process_set_path.exists():
        manifest = _read_json(process_set_path)
        entries = manifest.get("processes", [])
        if not isinstance(entries, list):
            raise ValueError("process-set manifest processes must be a list")
        rows: list[tuple[str, int, list[str]]] = []
        unavailable: list[tuple[str, str]] = []
        for raw_entry in entries:
            if not isinstance(raw_entry, Mapping):
                continue
            key = str(raw_entry.get("key", ""))
            relative_path = Path(str(raw_entry.get("path", "")))
            process_root = relative_path if relative_path.is_absolute() else root / relative_path
            point, error = _load_first_validation_point(process_root)
            if point is None:
                unavailable.append((key, error or "no validation point"))
                continue
            crossing_map = raw_entry.get("input_crossing_map")
            if isinstance(crossing_map, list) and crossing_map:
                point = _invert_input_crossing_map(point, crossing_map)
            rows.append((key, len(point), _flatten_point(point)))
        return rows, unavailable

    manifest = _read_json(root / "process_manifest.json")
    key = str(manifest.get("key") or manifest.get("process") or "default")
    point, error = _load_first_validation_point(root)
    if point is None:
        return [], [(key, error or "no validation point")]
    return [(key, len(point), _flatten_point(point))], []


def _load_first_validation_point(
    root: Path,
) -> tuple[list[list[str]] | None, str | None]:
    path = root / "validation_momenta.json"
    if not path.exists():
        return None, f"missing {path.name}"
    payload = _read_json(path)
    points = payload.get("points")
    if payload.get("available") is False or not isinstance(points, list) or not points:
        return None, str(payload.get("error") or "no validation point is available")
    raw_point = points[0]
    if not isinstance(raw_point, list):
        return None, "validation point is not a particle list"
    point: list[list[str]] = []
    for particle in raw_point:
        if not isinstance(particle, Mapping):
            return None, "validation particle is not an object"
        momentum = particle.get("momentum")
        if not isinstance(momentum, list) or len(momentum) != 4:
            return None, "validation momentum does not have four components"
        point.append([str(component) for component in momentum])
    return point, None


def _invert_input_crossing_map(
    representative_point: Sequence[Sequence[str]],
    raw_map: Sequence[object],
) -> list[list[str]]:
    selected: list[list[str] | None] = [None] * len(representative_point)
    for raw_entry in raw_map:
        if not isinstance(raw_entry, Mapping):
            raise ValueError("input_crossing_map entry must be an object")
        target_index = int(raw_entry["target_index"])
        source_index = int(raw_entry["source_index"])
        sign = float(raw_entry["sign"])
        selected[source_index] = [
            _signed_decimal(component, sign) for component in representative_point[target_index]
        ]
    if any(momentum is None for momentum in selected):
        raise ValueError("input_crossing_map does not cover every selected external leg")
    return [momentum for momentum in selected if momentum is not None]


def _signed_decimal(value: str, sign: float) -> str:
    if sign >= 0.0 or _decimal_is_zero(value):
        return value
    return value[1:] if value.startswith("-") else f"-{value}"


def _decimal_is_zero(value: str) -> bool:
    try:
        return float(value) == 0.0
    except ValueError:
        return False


def _flatten_point(point: Sequence[Sequence[str]]) -> list[str]:
    return [component for momentum in point for component in momentum]


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON document is not an object: {path}")
    return payload


__all__ = ["remove_api_bundle", "write_api_bundle"]
