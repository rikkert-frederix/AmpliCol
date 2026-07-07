from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from .generic_artifact import (
    GenericProcessManifest,
    build_generic_process_manifest,
    write_generic_process_manifest,
)
from .process_support import ProcessSupportReport, classify_process_support
from .processes import ProcessOptions


MatrixGenerationBackend = Literal[
    "generic-dag-schema-v2",
    "generic-dag-preflight-blocked",
]


@dataclass(frozen=True)
class GenerationResult:
    """Result of model-driven generic DAG matrix-element planning.

    This production-facing result is intentionally process-family agnostic.
    Retired staged native graph metadata remains in reference-only modules for
    migration diagnostics only.
    """

    process: str
    backend: MatrixGenerationBackend
    runtime_artifact_supported: bool
    generic_dag_runtime_supported: bool
    cache_file: Path | None
    generation_time_s: float
    manifest_file: Path | None
    manifest: GenericProcessManifest | None
    support_report: ProcessSupportReport
    notes: tuple[str, ...] = ()

    @property
    def artifact_file(self) -> Path | None:
        """Path to the generic process manifest emitted by this planning run."""

        return self.manifest_file

    @property
    def artifact_fingerprint(self) -> None:
        return None

    @property
    def artifact_cache_hit(self) -> bool:
        return False

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "backend": self.backend,
            "runtime_artifact_supported": self.runtime_artifact_supported,
            "generic_dag_runtime_supported": self.generic_dag_runtime_supported,
            "cache_file": None if self.cache_file is None else str(self.cache_file),
            "generation_time_s": self.generation_time_s,
            "manifest_file": (
                None if self.manifest_file is None else str(self.manifest_file)
            ),
            "artifact_file": (
                None if self.artifact_file is None else str(self.artifact_file)
            ),
            "manifest": (
                None if self.manifest is None else self.manifest.to_json_dict()
            ),
            "support_report": self.support_report.to_json_dict(),
            "notes": list(self.notes),
        }


class MatrixElementGenerator:
    """Process-generic, model-driven matrix-element planning facade."""

    def __init__(self, *, cache_dir: str | Path | None = None) -> None:
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir is not None else None

    def generate(
        self,
        process: str,
        *,
        options: ProcessOptions | None = None,
        color_accuracy: str = "lc",
        write_cache_metadata: bool = True,
        max_currents: int = 50000,
        max_color_sectors: int = 20000,
    ) -> GenerationResult:
        start = time.perf_counter()
        support_report = classify_process_support(
            process,
            color_accuracy=color_accuracy,
            options=options,
            color_plan_max_sectors=max_color_sectors,
            current_plan_max_currents=max_currents,
        )
        manifest: GenericProcessManifest | None = None
        manifest_file: Path | None = None
        cache_file: Path | None = None
        notes: list[str] = []
        if support_report.runtime_artifact_supported:
            manifest = build_generic_process_manifest(
                process,
                options=options,
                color_accuracy=color_accuracy,
                max_currents=max_currents,
                max_color_sectors=max_color_sectors,
            )
            backend: MatrixGenerationBackend = "generic-dag-schema-v2"
            notes.append(
                "Generic schema-v2 DAG planning succeeded. Use "
                "`generate-process PROCESS OUTPUT_DIR` to emit serialized "
                "Rusticol evaluator artifacts."
            )
        else:
            backend = "generic-dag-preflight-blocked"
            notes.append(
                support_report.artifact_unavailable_message
                or "generic DAG preflight did not produce an executable artifact"
            )
        result = GenerationResult(
            process=process,
            backend=backend,
            runtime_artifact_supported=support_report.runtime_artifact_supported,
            generic_dag_runtime_supported=support_report.generic_dag_runtime_supported,
            cache_file=None,
            generation_time_s=time.perf_counter() - start,
            manifest_file=None,
            manifest=manifest,
            support_report=support_report,
            notes=tuple(notes),
        )
        if write_cache_metadata and self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            cache_file = self.cache_dir / f"{_safe_process_name(process)}.generic.json"
            if manifest is not None:
                manifest_dir = self.cache_dir / _safe_process_name(process)
                manifest_file = write_generic_process_manifest(manifest, manifest_dir)
            result = GenerationResult(
                process=result.process,
                backend=result.backend,
                runtime_artifact_supported=result.runtime_artifact_supported,
                generic_dag_runtime_supported=result.generic_dag_runtime_supported,
                cache_file=cache_file,
                generation_time_s=result.generation_time_s,
                manifest_file=manifest_file,
                manifest=result.manifest,
                support_report=result.support_report,
                notes=result.notes,
            )
            cache_file.write_text(
                json.dumps(result.to_json_dict(), indent=2, sort_keys=True),
                encoding="utf-8",
            )
        return result


def _safe_process_name(process: str) -> str:
    return (
        process.replace(" ", "_")
        .replace(">", "to")
        .replace("~", "bar")
        .replace("+", "p")
        .replace("-", "m")
        .replace("/", "_")
    )


__all__ = [
    "GenerationResult",
    "MatrixElementGenerator",
    "MatrixGenerationBackend",
]
