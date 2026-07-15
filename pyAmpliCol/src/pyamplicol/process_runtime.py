from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Sequence, cast

from .core_types import NativeEvaluationError


ProcessRuntimeBackend = Literal["python", "rusticol"]


@dataclass(frozen=True)
class ProcessArtifactManifest:
    root: Path
    process: str
    external_pdg_order: tuple[int, ...]
    compiled_kind: str
    schema_version: int = 2
    kind: str = "pyamplicol-generic-dag-process"
    artifact_class: str = "generic-dag-schema-v2"
    key: str | None = None
    color_accuracy: str | None = None


class PythonProcessRuntime:
    """Schema-v2 Python facade backed by the shared Rusticol runtime core."""

    def __init__(
        self,
        process_dir: str | Path,
        *,
        process_key: str | None = None,
        model_parameters: str | Path | None = None,
    ) -> None:
        root = Path(process_dir).expanduser()
        self.manifest = load_process_manifest(
            _resolve_process_artifact_dir(root, process_key)
        )
        rusticol = cast(Any, __import__("rusticol"))
        self._runtime = rusticol.Runtime.load(
            str(root),
            process_key,
            str(Path(model_parameters).expanduser())
            if model_parameters is not None
            else None,
        )

    @property
    def process(self) -> str:
        return cast(str, self._runtime.process)

    @property
    def metadata(self) -> dict[str, object]:
        return cast(dict[str, object], self._runtime.metadata())

    @property
    def physics(self) -> Any:
        return self._runtime.physics

    def evaluate(self, momenta: Sequence[Sequence[Sequence[float]]] | Any) -> Any:
        return self._runtime.evaluate(momenta)

    def evaluate_with_prec(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | Any,
        precision: int,
    ) -> Any:
        return self._runtime.evaluate_with_prec(momenta, precision)

    def evaluate_resolved(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | Any,
        **selectors: object,
    ) -> Any:
        return self._runtime.evaluate_resolved(momenta, **selectors)

    def evaluate_resolved_with_prec(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | Any,
        precision: int,
        **selectors: object,
    ) -> Any:
        return self._runtime.evaluate_resolved_with_prec(
            momenta,
            precision,
            **selectors,
        )

    def profile(self, momenta: Sequence[Sequence[Sequence[float]]] | Any) -> Any:
        return self._runtime.profile(momenta)

    def set_model_parameters(self, parameters: dict[str, object]) -> None:
        self._runtime.set_model_parameters(parameters)

    def set_model_parameter(self, name: str, value: object) -> None:
        self._runtime.set_model_parameters({name: value})

    def mute_warnings(self) -> None:
        self._runtime.mute_warnings()

    def unmute_warnings(self) -> None:
        self._runtime.unmute_warnings()


def load_process(
    process_dir: str | Path,
    *,
    runtime: ProcessRuntimeBackend = "python",
    process_key: str | None = None,
    model_parameters: str | Path | None = None,
) -> Any:
    """Load a generated schema-v2 process through the Python or raw binding API."""

    if runtime == "python":
        return PythonProcessRuntime(
            process_dir,
            process_key=process_key,
            model_parameters=model_parameters,
        )
    if runtime == "rusticol":
        rusticol = cast(Any, __import__("rusticol"))
        return rusticol.Runtime.load(
            str(Path(process_dir).expanduser()),
            process_key,
            str(Path(model_parameters).expanduser())
            if model_parameters is not None
            else None,
        )
    raise ValueError(f"unknown process runtime: {runtime!r}")


def load_process_manifest(process_dir: str | Path) -> ProcessArtifactManifest:
    root = _resolve_process_artifact_dir(process_dir)
    payload = _load_process_manifest_payload(root)
    if (
        payload.get("kind") != "pyamplicol-generic-dag-process"
        or int(payload.get("schema_version", 0)) != 2
    ):
        raise NativeEvaluationError(
            "schema-v1 process artifacts are unsupported; regenerate this process "
            "as a schema-v2 generic DAG artifact"
        )
    compiled = payload.get("compiled")
    external_pdg_order = payload.get("external_pdg_order")
    if not isinstance(external_pdg_order, list):
        raise NativeEvaluationError(
            "generic DAG process artifact is missing external_pdg_order"
        )
    return ProcessArtifactManifest(
        root=root,
        process=str(payload["process"]),
        external_pdg_order=tuple(int(pdg) for pdg in external_pdg_order),
        compiled_kind=(
            str(compiled.get("kind")) if isinstance(compiled, dict) else "unknown"
        ),
        artifact_class=str(payload.get("artifact_class", "generic-dag-schema-v2")),
        key=str(payload.get("key")) if payload.get("key") is not None else None,
        color_accuracy=(
            str(payload.get("color_accuracy"))
            if payload.get("color_accuracy") is not None
            else None
        ),
    )


def _load_process_manifest_payload(process_dir: str | Path) -> dict[str, Any]:
    manifest_path = Path(process_dir).expanduser() / "process_manifest.json"
    if not manifest_path.exists():
        raise NativeEvaluationError(
            f"process artifact manifest does not exist: {manifest_path}"
        )
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise NativeEvaluationError("process artifact manifest is not a JSON object")
    return payload


def _resolve_process_artifact_dir(
    process_dir: str | Path,
    process_key: str | None = None,
) -> Path:
    root = Path(process_dir).expanduser()
    if (root / "process_manifest.json").exists():
        if process_key is not None:
            raise NativeEvaluationError(
                "process_key can only be used with a process-set artifact"
            )
        return root
    process_set_path = root / "process_set_manifest.json"
    if not process_set_path.exists():
        return root
    payload = json.loads(process_set_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise NativeEvaluationError("process-set manifest is not a JSON object")
    if (
        payload.get("kind") != "pyamplicol-generic-dag-process-set"
        or int(payload.get("schema_version", 0)) != 2
    ):
        raise NativeEvaluationError(
            "schema-v1 process sets are unsupported; regenerate this process set"
        )
    entries = [
        entry for entry in payload.get("processes", []) if isinstance(entry, dict)
    ]
    if not entries:
        raise NativeEvaluationError("process-set artifact contains no subprocesses")
    selected = process_key or str(payload.get("default_process_key"))
    for entry in entries:
        if selected in {str(entry.get("key")), str(entry.get("process"))}:
            path = Path(str(entry.get("path")))
            return path if path.is_absolute() else root / path
    available = ", ".join(str(entry.get("key")) for entry in entries)
    raise NativeEvaluationError(
        f"process {selected!r} not found in process set; available: {available}"
    )


__all__ = [
    "ProcessArtifactManifest",
    "ProcessRuntimeBackend",
    "PythonProcessRuntime",
    "load_process",
    "load_process_manifest",
]
