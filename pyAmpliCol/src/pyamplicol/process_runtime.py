from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, Sequence

import numpy as np

from .model import AmplicolSMLeadingColorModel
from .native import ExternalMomentum, MatrixElementEvaluation, NativeEvaluationError

ProcessRuntimeBackend = Literal["python", "rusticol"]


@dataclass(frozen=True)
class ProcessArtifactManifest:
    root: Path
    process: str
    family: str
    gluon_count: int
    external_pdg_order: tuple[int, ...]
    compiled_kind: str


class PythonProcessRuntime:
    """Load a self-contained pyAmpliCol process directory in Python."""

    def __init__(self, process_dir: str | Path) -> None:
        self.manifest = load_process_manifest(process_dir)
        self._payload = _load_process_manifest_payload(process_dir)
        self._compiled_sweep: Any | None = None
        self._model = AmplicolSMLeadingColorModel()
        self._zero_gluon_evaluator: Any | None = None

    @property
    def process(self) -> str:
        return self.manifest.process

    @property
    def metadata(self) -> dict[str, object]:
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
        points = self._points_from_momenta(momenta)
        if self._is_zero_gluon_symbolic_scalar():
            return self._evaluate_zero_gluon_full(points)
        raw_sums = self._get_compiled_sweep().evaluate_raw_sum_rows(
            points,
            self._model,
            gluon_count=self.manifest.gluon_count,
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
            gluon_count=self.manifest.gluon_count,
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
                    tuple(float(component) for component in row[index]),
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
        color_factor = float(norm["color_factor"])
        average_factor = float(norm["average_factor"])
        identical_factor = float(norm["identical_factor"])
        coupling_factor = float(norm["coupling_factor"])
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
) -> Any:
    """Load a generated process directory with either Python or rusticol."""

    if runtime == "python":
        return PythonProcessRuntime(process_dir)
    if runtime == "rusticol":
        import rusticol  # type: ignore[import-not-found]

        return rusticol.Runtime.load(str(Path(process_dir).expanduser()))
    raise ValueError(f"unknown process runtime: {runtime!r}")


def load_process_manifest(process_dir: str | Path) -> ProcessArtifactManifest:
    root = Path(process_dir).expanduser()
    payload = _load_process_manifest_payload(root)
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
        family=str(payload["family"]),
        gluon_count=int(payload["gluon_count"]),
        external_pdg_order=tuple(int(pdg) for pdg in payload["external_pdg_order"]),
        compiled_kind=compiled_kind,
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


def _artifact_path_from_manifest(root: Path, path: str) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return root / candidate


__all__ = [
    "ProcessArtifactManifest",
    "ProcessRuntimeBackend",
    "PythonProcessRuntime",
    "load_process",
    "load_process_manifest",
]
