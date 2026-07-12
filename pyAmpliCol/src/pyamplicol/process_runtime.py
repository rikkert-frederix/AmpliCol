from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Sequence, cast

import numpy as np

from .core_types import ExternalMomentum, MatrixElementEvaluation, NativeEvaluationError
from .model import AmplicolSMLeadingColorModel

ProcessRuntimeBackend = Literal["python", "rusticol"]


@dataclass(frozen=True)
class ProcessArtifactManifest:
    root: Path
    process: str
    external_pdg_order: tuple[int, ...]
    compiled_kind: str
    schema_version: int = 1
    kind: str = "pyamplicol-rusticol-process"
    artifact_class: str = "legacy-schema-v1"
    family: str | None = None
    gluon_count: int | None = None
    key: str | None = None
    color_accuracy: str | None = None


class PythonProcessRuntime:
    """Load a self-contained pyAmpliCol process directory in Python."""

    def __init__(
        self,
        process_dir: str | Path,
        *,
        process_key: str | None = None,
        allow_legacy_schema_v1: bool = False,
    ) -> None:
        resolved = _resolve_process_artifact_dir(process_dir, process_key)
        self.manifest = load_process_manifest(resolved)
        self._payload = _load_process_manifest_payload(resolved)
        if not self._is_generic_schema_v2() and not allow_legacy_schema_v1:
            raise NativeEvaluationError(
                "schema-v1 process artifacts are retired from production; "
                "generate a schema-v2 generic DAG artifact or pass "
                "allow_legacy_schema_v1=True only for reference-only legacy "
                "artifact diagnostics"
            )
        self._compiled_sweep: Any | None = None
        self._model = AmplicolSMLeadingColorModel()
        self._zero_gluon_evaluator: Any | None = None

    @property
    def process(self) -> str:
        return self.manifest.process

    @property
    def metadata(self) -> dict[str, object]:
        if self._is_generic_schema_v2():
            return self._generic_metadata()
        if self._is_zero_gluon_symbolic_scalar():
            return {
                "process": self.manifest.process,
                "family": self.manifest.family,
                "gluon_count": self.manifest.gluon_count,
                "runtime": "python-zero-gluon-symbolic-scalar",
                "compiled_kind": self.manifest.compiled_kind,
                "evaluator": (
                    self._get_zero_gluon_evaluator().metadata.to_json_dict()
                ),
            }
        compiled = self._get_compiled_sweep()
        return {
            "process": self.manifest.process,
            "family": self.manifest.family,
            "gluon_count": self.manifest.gluon_count,
            "runtime": "python-eager-dag",
            "compiled_kind": self.manifest.compiled_kind,
            "parameter_count": compiled.parameter_count,
            "output_length": compiled.output_length,
            "current_count": len(compiled.table.currents),
            "stage_count": len(compiled.stages),
            "loaded_from_process_manifest": True,
        }

    def evaluate(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | np.ndarray,
    ) -> tuple[float, ...]:
        return tuple(
            evaluation.matrix_element
            for evaluation in self.evaluate_full(momenta)
        )

    def evaluate_full(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | np.ndarray,
    ) -> tuple[MatrixElementEvaluation, ...]:
        if self._is_generic_schema_v2():
            self._points_from_momenta(momenta)
            raise NativeEvaluationError(self._generic_execution_error())
        points = self._points_from_momenta(momenta)
        if self._is_zero_gluon_symbolic_scalar():
            return self._evaluate_zero_gluon_full(points)
        raw_sums = self._get_compiled_sweep().evaluate_raw_sum_rows(
            points,
            self._model,
            gluon_count=self._legacy_gluon_count(),
        )
        return tuple(
            self._evaluation_from_raw_sum(point, float(raw_sum))
            for point, raw_sum in zip(points, raw_sums, strict=True)
        )

    def profile(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | np.ndarray,
    ) -> dict[str, object]:
        points = self._points_from_momenta(momenta)
        start = time.perf_counter()
        values = self.evaluate(momenta)
        elapsed = time.perf_counter() - start
        timing = (
            None
            if self._is_zero_gluon_symbolic_scalar()
            else self._get_compiled_sweep().last_timing
        )
        return {
            "points": len(points),
            "values": list(values),
            "total_time_s": elapsed,
            "runtime_s_per_point": elapsed / max(len(points), 1),
            "timing": (
                None
                if timing is None or not hasattr(timing, "to_json_dict")
                else timing.to_json_dict()
            ),
        }

    def stage_diagnostics(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | np.ndarray,
    ) -> dict[str, object]:
        if self._is_generic_schema_v2():
            points = self._points_from_momenta(momenta)
            runtime_schema = self._generic_runtime_schema()
            return {
                "points": len(points),
                "schema_version": 2,
                "runtime": "python-generic-dag-schema-v2",
                "runtime_available": self._generic_runtime_available(),
                "process": self.manifest.process,
                "key": self.manifest.key,
                "stages": [
                    {
                        "stage": stage.get("stage_index"),
                        "stage_kind": stage.get("stage_kind"),
                        "subset_size": stage.get("subset_size"),
                        "interaction_count": stage.get(
                            "interaction_count",
                            len(stage.get("interactions", [])),
                        ),
                        "interaction_metadata_compacted": bool(
                            stage.get("interactions_compacted", False)
                        ),
                        "output_value_slot_count": len(
                            stage.get("output_value_slot_ids", [])
                        ),
                    }
                    for stage in runtime_schema.get("stages", [])
                    if isinstance(stage, dict)
                ],
                "amplitude_stage": runtime_schema.get("amplitude_stage"),
            }
        if self._is_zero_gluon_symbolic_scalar():
            points = self._points_from_momenta(momenta)
            return {
                "points": len(points),
                "stages": [
                    {
                        "stage": "zero_gluon_symbolic_scalar",
                        "output_len": len(points),
                    }
                ],
            }
        return self._get_compiled_sweep().stage_diagnostics(
            self._points_from_momenta(momenta),
            self._model,
            gluon_count=self._legacy_gluon_count(),
        )

    def _points_from_momenta(
        self,
        momenta: Sequence[Sequence[Sequence[float]]] | np.ndarray,
    ) -> tuple[tuple[ExternalMomentum, ...], ...]:
        array = np.asarray(momenta, dtype=np.float64)
        if array.ndim == 2:
            array = array.reshape((1, *array.shape))
        expected_shape_tail = (len(self.manifest.external_pdg_order), 4)
        if tuple(array.shape[1:]) != expected_shape_tail:
            raise NativeEvaluationError(
                "momenta must have shape "
                f"(batch, {expected_shape_tail[0]}, {expected_shape_tail[1]}), "
                f"got {tuple(array.shape)}"
            )
        return tuple(
            tuple(
                ExternalMomentum(
                    pdg,
                    cast(
                        tuple[float, float, float, float],
                        tuple(float(component) for component in row[index]),
                    ),
                )
                for index, pdg in enumerate(self.manifest.external_pdg_order)
            )
            for row in array
        )

    def _get_compiled_sweep(self) -> Any:
        if self._compiled_sweep is None:
            if self.manifest.compiled_kind != "shared-compiled-sweep":
                raise NativeEvaluationError(
                    "process artifact does not contain a shared-current sweep"
                )
            from .dag_runtime import (
                SymbolicaEvaluatorSettings,
                _SharedCompiledSweepEvaluator,
                _shared_current_table_from_rusticol_manifest,
            )

            table = _shared_current_table_from_rusticol_manifest(self._payload)
            self._compiled_sweep = _SharedCompiledSweepEvaluator.from_artifact(
                table,
                self.manifest.root,
                symbolica_settings=SymbolicaEvaluatorSettings(),
            )
        return self._compiled_sweep

    def _is_zero_gluon_symbolic_scalar(self) -> bool:
        return (
            self.manifest.gluon_count == 0
            and self.manifest.compiled_kind == "zero-gluon-symbolic-scalar"
        )

    def _is_generic_schema_v2(self) -> bool:
        return (
            self.manifest.schema_version == 2
            and self.manifest.kind == "pyamplicol-generic-dag-process"
        )

    def _legacy_gluon_count(self) -> int:
        if self.manifest.gluon_count is None:
            raise NativeEvaluationError(
                "generic schema-v2 artifacts do not expose a family-level gluon_count"
            )
        return self.manifest.gluon_count

    def _generic_runtime_schema(self) -> dict[str, Any]:
        runtime_schema = self._payload.get("runtime_schema")
        if not isinstance(runtime_schema, dict):
            raise NativeEvaluationError("generic schema-v2 artifact is missing runtime_schema")
        return runtime_schema

    def _generic_runtime_available(self) -> bool:
        compiled = self._payload.get("compiled")
        return bool(
            isinstance(compiled, dict)
            and compiled.get("runtime_available") is True
            and isinstance(compiled.get("stage_evaluators"), dict)
        )

    def _generic_execution_error(self) -> str:
        compiled = self._payload.get("compiled")
        if isinstance(compiled, dict):
            message = compiled.get("runtime_unavailable_message")
            if isinstance(message, str) and message:
                return message
        return (
            "generic schema-v2 Python execution is not available; use "
            "load_process(..., runtime='rusticol') for serialized evaluator execution"
        )

    def _generic_metadata(self) -> dict[str, object]:
        runtime_schema = self._generic_runtime_schema()
        dag_summary = self._payload.get("dag_summary")
        compiled = self._payload.get("compiled")
        parameter_layout = runtime_schema.get("parameter_layout")
        source_fill = runtime_schema.get("source_fill")
        stages = runtime_schema.get("stages")
        amplitude_stage = runtime_schema.get("amplitude_stage")
        return {
            "process": self.manifest.process,
            "key": self.manifest.key,
            "artifact_class": self.manifest.artifact_class,
            "schema_version": 2,
            "kind": self.manifest.kind,
            "color_accuracy": self.manifest.color_accuracy,
            "runtime": "python-generic-dag-schema-v2",
            "runtime_available": self._generic_runtime_available(),
            "compiled_kind": self.manifest.compiled_kind,
            "external_count": len(self.manifest.external_pdg_order),
            "current_count": (
                dag_summary.get("current_count")
                if isinstance(dag_summary, dict)
                else None
            ),
            "source_count": (
                dag_summary.get("source_count")
                if isinstance(dag_summary, dict)
                else None
            ),
            "interaction_count": (
                dag_summary.get("interaction_count")
                if isinstance(dag_summary, dict)
                else None
            ),
            "amplitude_root_count": (
                dag_summary.get("amplitude_root_count")
                if isinstance(dag_summary, dict)
                else None
            ),
            "stage_count": len(stages) if isinstance(stages, list) else None,
            "source_fill_count": (
                len(source_fill.get("sources", []))
                if isinstance(source_fill, dict)
                else None
            ),
            "parameter_count": (
                parameter_layout.get("parameter_count_if_flattened")
                if isinstance(parameter_layout, dict)
                else None
            ),
            "amplitude_output_count": (
                amplitude_stage.get("output_count")
                if isinstance(amplitude_stage, dict)
                else None
            ),
            "stage_evaluator_count": (
                len(compiled.get("stage_evaluators", {}).get("stages", [])) + 1
                if isinstance(compiled, dict)
                and isinstance(compiled.get("stage_evaluators"), dict)
                else 0
            ),
            "loaded_from_process_manifest": True,
        }

    def _get_zero_gluon_evaluator(self) -> Any:
        if self._zero_gluon_evaluator is None:
            from .symbolic import ZeroGluonSymbolicEvaluator

            zero_manifest = self._zero_gluon_manifest()
            state_path = _artifact_path_from_manifest(
                self.manifest.root,
                str(zero_manifest["evaluator_state_path"]),
            )
            parameter_names = zero_manifest.get("parameter_names")
            if not isinstance(parameter_names, list) or not all(
                isinstance(name, str) for name in parameter_names
            ):
                raise NativeEvaluationError(
                    "zero-gluon process artifact is missing parameter_names"
                )
            self._zero_gluon_evaluator = ZeroGluonSymbolicEvaluator(
                evaluator_state=state_path.read_bytes(),
                parameter_names=parameter_names,
            )
        return self._zero_gluon_evaluator

    def _zero_gluon_manifest(self) -> dict[str, object]:
        compiled = self._payload.get("compiled")
        if not isinstance(compiled, dict):
            raise NativeEvaluationError("process artifact is missing compiled metadata")
        zero = compiled.get("zero_gluon")
        if not isinstance(zero, dict):
            raise NativeEvaluationError("zero-gluon process artifact is missing zero_gluon")
        return zero

    def _evaluate_zero_gluon_full(
        self,
        points: Sequence[Sequence[ExternalMomentum]],
    ) -> tuple[MatrixElementEvaluation, ...]:
        evaluator = self._get_zero_gluon_evaluator()
        return tuple(
            evaluator.evaluate(self.manifest.process, tuple(point))
            for point in points
        )

    def _evaluation_from_raw_sum(
        self,
        point: Sequence[ExternalMomentum],
        raw_sum: float,
    ) -> MatrixElementEvaluation:
        norm = self._normalization_payload()
        color_factor = float(cast(Any, norm["color_factor"]))
        average_factor = float(cast(Any, norm["average_factor"]))
        identical_factor = float(cast(Any, norm["identical_factor"]))
        coupling_factor = float(cast(Any, norm["coupling_factor"]))
        matrix_element = (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
        )
        return MatrixElementEvaluation(
            process=self.manifest.process,
            particles=tuple(point),
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=int(color_factor),
            average_factor=int(average_factor),
            coupling_factor=coupling_factor,
            helicity_contributions=(),
            identical_factor=int(identical_factor),
        )

    def _normalization_payload(self) -> dict[str, object]:
        normalization = self._payload.get("normalization")
        if not isinstance(normalization, dict):
            raise NativeEvaluationError("process artifact is missing normalization")
        return normalization


def load_process(
    process_dir: str | Path,
    *,
    runtime: ProcessRuntimeBackend = "python",
    process_key: str | None = None,
    allow_legacy_schema_v1: bool = False,
) -> Any:
    """Load a generated process directory with either Python or rusticol."""

    if runtime == "python":
        return PythonProcessRuntime(
            process_dir,
            process_key=process_key,
            allow_legacy_schema_v1=allow_legacy_schema_v1,
        )
    if runtime == "rusticol":
        rusticol = cast(Any, __import__("rusticol"))

        return rusticol.Runtime.load(str(Path(process_dir).expanduser()), process_key)
    raise ValueError(f"unknown process runtime: {runtime!r}")


def load_process_manifest(process_dir: str | Path) -> ProcessArtifactManifest:
    root = _resolve_process_artifact_dir(process_dir)
    payload = _load_process_manifest_payload(root)
    if payload.get("kind") == "pyamplicol-generic-dag-process":
        return _generic_process_manifest_from_payload(root, payload)
    if payload.get("kind") != "pyamplicol-rusticol-process":
        raise NativeEvaluationError(
            f"unsupported process artifact kind: {payload.get('kind')!r}"
        )
    compiled = payload.get("compiled")
    compiled_kind = (
        str(compiled.get("kind"))
        if isinstance(compiled, dict)
        else "unknown"
    )
    return ProcessArtifactManifest(
        root=root,
        process=str(payload["process"]),
        external_pdg_order=tuple(int(pdg) for pdg in payload["external_pdg_order"]),
        compiled_kind=compiled_kind,
        artifact_class="legacy-schema-v1",
        family=str(payload["family"]),
        gluon_count=int(payload["gluon_count"]),
    )


def _generic_process_manifest_from_payload(
    root: Path,
    payload: dict[str, Any],
) -> ProcessArtifactManifest:
    if int(payload.get("schema_version", 0)) != 2:
        raise NativeEvaluationError(
            "generic DAG process artifact must use schema_version 2"
        )
    compiled = payload.get("compiled")
    compiled_kind = (
        str(compiled.get("kind"))
        if isinstance(compiled, dict)
        else "unknown"
    )
    external_pdg_order = payload.get("external_pdg_order")
    if not isinstance(external_pdg_order, list):
        raise NativeEvaluationError(
            "generic DAG process artifact is missing external_pdg_order"
        )
    return ProcessArtifactManifest(
        root=root,
        process=str(payload["process"]),
        external_pdg_order=tuple(int(pdg) for pdg in external_pdg_order),
        compiled_kind=compiled_kind,
        schema_version=2,
        kind="pyamplicol-generic-dag-process",
        artifact_class=str(payload.get("artifact_class", "generic-dag-schema-v2")),
        key=str(payload.get("key")) if payload.get("key") is not None else None,
        color_accuracy=(
            str(payload.get("color_accuracy"))
            if payload.get("color_accuracy") is not None
            else None
        ),
    )


def _load_process_manifest_payload(process_dir: str | Path) -> dict[str, Any]:
    root = Path(process_dir).expanduser()
    manifest_path = root / "process_manifest.json"
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
    if payload.get("kind") not in {
        "pyamplicol-rusticol-process-set",
        "pyamplicol-generic-dag-process-set",
    }:
        raise NativeEvaluationError(
            f"unsupported process-set artifact kind: {payload.get('kind')!r}"
        )
    entries = [entry for entry in payload.get("processes", []) if isinstance(entry, dict)]
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


def _artifact_path_from_manifest(root: Path, path: str) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return root / candidate


def _generic_dag_runtime_unavailable_message(
    payload: dict[str, Any],
) -> str | None:
    if payload.get("kind") != "pyamplicol-generic-dag-process":
        return None
    compiled = payload.get("compiled")
    if isinstance(compiled, dict):
        if bool(compiled.get("runtime_available", False)) and isinstance(
            compiled.get("stage_evaluators"),
            dict,
        ):
            return None
        message = compiled.get("runtime_unavailable_message")
        if isinstance(message, str) and message:
            return message
    return (
        "generic DAG evaluator stages are symbolically lowered, but serialized "
        "evaluator emission and Rusticol schema-v2 execution are not available"
    )


__all__ = [
    "ProcessArtifactManifest",
    "ProcessRuntimeBackend",
    "PythonProcessRuntime",
    "load_process",
    "load_process_manifest",
]
