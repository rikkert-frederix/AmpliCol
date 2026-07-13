from __future__ import annotations

import hashlib
from importlib import resources
from pathlib import Path


BUNDLED_MODEL_NAMES = ("sm", "scalars", "scalar_gravity")
BUNDLED_MODEL_FORMATS = ("ufo", "json")


def packaged_models_root() -> Path:
    return Path(str(resources.files("pyamplicol").joinpath("assets", "models")))


def bundled_model_path(name: str, model_format: str = "ufo") -> Path:
    if name not in BUNDLED_MODEL_NAMES:
        raise ValueError(
            f"unknown bundled model {name!r}; expected one of "
            f"{', '.join(BUNDLED_MODEL_NAMES)}"
        )
    if model_format not in BUNDLED_MODEL_FORMATS:
        raise ValueError(
            f"unknown bundled model format {model_format!r}; expected ufo or json"
        )
    root = packaged_models_root() / model_format / name
    return root if model_format == "ufo" else root / f"{name}.json"


def verify_model_asset_manifest(root: Path | None = None) -> tuple[str, ...]:
    active_root = packaged_models_root() if root is None else Path(root)
    manifest_path = active_root / "MANIFEST.sha256"
    errors: list[str] = []
    seen: set[Path] = set()
    for line_number, line in enumerate(
        manifest_path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        if not line.strip():
            continue
        try:
            expected, relative_text = line.split(maxsplit=1)
        except ValueError:
            errors.append(f"line {line_number}: malformed manifest entry")
            continue
        relative = Path(relative_text)
        path = active_root / relative
        seen.add(relative)
        if not path.is_file():
            errors.append(f"missing: {relative.as_posix()}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            errors.append(
                f"digest mismatch: {relative.as_posix()} ({actual} != {expected})"
            )
    actual_files = {
        path.relative_to(active_root)
        for directory in BUNDLED_MODEL_FORMATS
        for path in (active_root / directory).rglob("*")
        if path.is_file() and path.suffix != ".pyc" and "__pycache__" not in path.parts
    }
    for relative in sorted(actual_files - seen):
        errors.append(f"unlisted: {relative.as_posix()}")
    return tuple(errors)
