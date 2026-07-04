from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, Sequence

from .lowering import SymbolicLoweringReport, build_symbolic_lowering_report
from .model import AmplicolSMLeadingColorModel, Model
from .processes import ProcessEnumerator, ProcessOptions
from .symbolic import build_zero_gluon_symbolic_evaluator_payload

_AUTO_TENSOR_NETWORK_ARTIFACT_MAX_GLUONS = 4


@dataclass(frozen=True)
class CurrentKey:
    pdg: int
    external_labels: tuple[int, ...]
    chirality: int = 0


@dataclass(frozen=True)
class InteractionNode:
    vertex_kind: int
    left: CurrentKey
    right: CurrentKey
    result: CurrentKey
    coupling: tuple[float, float]


@dataclass(frozen=True)
class RecursionGraph:
    process: tuple[int, ...]
    color_order: tuple[int, ...]
    currents: tuple[CurrentKey, ...]
    interactions: tuple[InteractionNode, ...]
    amplitudes: tuple[tuple[CurrentKey, CurrentKey], ...]


@dataclass(frozen=True)
class GenerationResult:
    process: str
    backend: Literal[
        "native-python-recursion-staged",
        "native-symbolic-pending",
        "not-implemented",
    ]
    supported_native_target: bool
    cache_file: Path | None
    generation_time_s: float
    graph: RecursionGraph | None
    artifact_file: Path | None = None
    artifact_fingerprint: str | None = None
    artifact_cache_hit: bool = False
    symbolic_lowering: SymbolicLoweringReport | None = None
    notes: tuple[str, ...] = ()

    def to_json_dict(self) -> dict[str, object]:
        data = asdict(self)
        data["cache_file"] = None if self.cache_file is None else str(self.cache_file)
        data["artifact_file"] = (
            None if self.artifact_file is None else str(self.artifact_file)
        )
        return data


@dataclass(frozen=True)
class EvaluatorArtifact:
    file: Path
    payload: dict[str, Any]
    load_time_s: float


class NativeMatrixElementGenerator:
    """Generation facade for the staged native pyamplicol ME backend."""

    def __init__(
        self,
        *,
        model: Model | None = None,
        cache_dir: str | Path | None = None,
    ) -> None:
        self.model = model or AmplicolSMLeadingColorModel()
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None

    def generate(
        self,
        process: str,
        *,
        options: ProcessOptions | None = None,
        write_cache_metadata: bool = True,
    ) -> GenerationResult:
        start = time.perf_counter()
        enumeration = ProcessEnumerator(options).enumerate(process)
        supported = self._is_supported_native_milestone(process)
        graph = None
        symbolic_lowering = None
        notes: list[str] = []
        backend: Literal[
            "native-python-recursion-staged",
            "native-symbolic-pending",
            "not-implemented",
        ] = "not-implemented"
        if supported:
            if not isinstance(self.model, AmplicolSMLeadingColorModel):
                raise TypeError(
                    "staged native generation requires AmplicolSMLeadingColorModel"
                )
            first_record = enumeration.groups[0].records[0]
            graph = self._build_staged_graph(
                tuple(int(x) for x in _pdg_process(first_record.process)),
                first_record.color_order,
            )
            symbolic_lowering = build_symbolic_lowering_report(self.model, graph=graph)
            gluon_count = sum(1 for pdg in graph.process[2:] if pdg == 21)
            backend = "native-python-recursion-staged"
            if gluon_count <= _AUTO_TENSOR_NETWORK_ARTIFACT_MAX_GLUONS:
                notes.append(
                    "Native ordered Python recursion is available for the "
                    "q q~ -> Z + n g family; full spenso/Symbolica tensor-network "
                    "evaluator artifacts are currently generated through four "
                    "final-state gluons."
                )
            else:
                notes.append(
                    "Native ordered Python recursion is available for the "
                    "q q~ -> Z + n g family; full spenso/Symbolica tensor-network "
                    "blueprints are generated, while evaluator artifacts above "
                    "four final-state gluons remain guarded by expression-size "
                    "growth."
                )
        else:
            notes.append(
                "Full process architecture is available, but native matrix-element lowering is currently implemented only for the q q~ -> Z + n g milestone family."
            )

        result = GenerationResult(
            process=process,
            backend=backend,
            supported_native_target=supported,
            cache_file=None,
            generation_time_s=time.perf_counter() - start,
            graph=graph,
            artifact_file=None,
            artifact_fingerprint=None,
            artifact_cache_hit=False,
            symbolic_lowering=symbolic_lowering,
            notes=tuple(notes),
        )
        if write_cache_metadata and self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            cache_file, artifact_file = _cache_files(self.cache_dir, process)
            fingerprint = _result_fingerprint(result)
            cache_hit = _artifact_matches(artifact_file, fingerprint)
            if not cache_hit:
                artifact_payload = _evaluator_artifact_payload(result, fingerprint)
                artifact_file.write_text(
                    json.dumps(artifact_payload, indent=2, sort_keys=True)
                )
            cache_file.write_text(
                json.dumps(
                    _with_cache_files(
                        result,
                        cache_file=cache_file,
                        artifact_file=artifact_file,
                        artifact_fingerprint=fingerprint,
                        artifact_cache_hit=cache_hit,
                    ).to_json_dict(),
                    indent=2,
                    sort_keys=True,
                )
            )
            result = _with_cache_files(
                result,
                cache_file=cache_file,
                artifact_file=artifact_file,
                artifact_fingerprint=fingerprint,
                artifact_cache_hit=cache_hit,
            )
        return result

    def _is_supported_native_milestone(self, process: str) -> bool:
        try:
            parsed = ProcessEnumerator(ProcessOptions(flavour_scheme=5)).parse(process)
        except ValueError:
            return False
        incoming = parsed.initial_state
        rest = parsed.rest
        has_z = rest.count("z") == 1
        only_z_gluons = has_z and all(p in {"z", "g"} for p in rest)
        return (
            len(incoming) == 2
            and incoming[0].endswith("~")
            and incoming[1].replace("~", "") == incoming[0].replace("~", "")
            and only_z_gluons
        )

    def _build_staged_graph(
        self, process: tuple[int, ...], color_order: tuple[int, ...]
    ) -> RecursionGraph:
        if not isinstance(self.model, AmplicolSMLeadingColorModel):
            raise TypeError(
                "staged native graph construction requires AmplicolSMLeadingColorModel"
            )
        external = tuple(
            CurrentKey(
                pdg=_external_current_type(pdg) if i < 2 else pdg,
                external_labels=(i + 1,),
            )
            for i, pdg in enumerate(process)
        )
        interactions: list[InteractionNode] = []
        currents: list[CurrentKey] = list(external)
        gluon_labels = tuple(
            index + 1 for index, pdg in enumerate(process) if index >= 2 and pdg == 21
        )
        z_labels = tuple(
            index + 1 for index, pdg in enumerate(process) if index >= 2 and pdg == 23
        )
        if len(z_labels) != 1:
            raise ValueError("staged Z+gluon graph expects exactly one final-state Z")
        z_key = external[z_labels[0] - 1]
        anti_closure = external[0]
        quark_start = external[1]
        gluon_currents: dict[tuple[int, int], CurrentKey] = {
            (i, i + 1): external[label - 1] for i, label in enumerate(gluon_labels)
        }
        tensor_currents: dict[tuple[int, int], CurrentKey] = {}

        for length in range(2, len(gluon_labels) + 1):
            for start in range(0, len(gluon_labels) - length + 1):
                end = start + length
                if length < len(gluon_labels):
                    tensor_result = CurrentKey(
                        pdg=-21,
                        external_labels=gluon_labels[start:end],
                    )
                    currents.append(tensor_result)
                    for split in range(start + 1, end):
                        interactions.append(
                            InteractionNode(
                                1,
                                gluon_currents[(start, split)],
                                gluon_currents[(split, end)],
                                tensor_result,
                                (1.0, 0.0),
                            )
                        )
                    tensor_currents[(start, end)] = tensor_result

                result = CurrentKey(
                    pdg=21,
                    external_labels=gluon_labels[start:end],
                )
                currents.append(result)
                for split in range(start + 1, end):
                    interactions.append(
                        InteractionNode(
                            0,
                            gluon_currents[(start, split)],
                            gluon_currents[(split, end)],
                            result,
                            (1.0, 0.0),
                        )
                    )
                    if split - start >= 2:
                        interactions.append(
                            InteractionNode(
                                2,
                                tensor_currents[(start, split)],
                                gluon_currents[(split, end)],
                                result,
                                (1.0, 0.0),
                            )
                        )
                    if end - split >= 2:
                        interactions.append(
                            InteractionNode(
                                3,
                                gluon_currents[(start, split)],
                                tensor_currents[(split, end)],
                                result,
                                (1.0, 0.0),
                            )
                        )
                gluon_currents[(start, end)] = result

        amplitudes: list[tuple[CurrentKey, CurrentKey]] = []
        for chirality in (1, -1):
            quark_without_z: dict[int, CurrentKey] = {
                0: CurrentKey(
                    pdg=quark_start.pdg,
                    external_labels=quark_start.external_labels,
                    chirality=chirality,
                )
            }
            quark_with_z: dict[int, CurrentKey] = {}
            currents.append(quark_without_z[0])
            for end in range(1, len(gluon_labels) + 1):
                result = CurrentKey(
                    pdg=quark_start.pdg,
                    external_labels=quark_start.external_labels + gluon_labels[:end],
                    chirality=chirality,
                )
                currents.append(result)
                for split in range(0, end):
                    interactions.append(
                        InteractionNode(
                            6,
                            quark_without_z[split],
                            gluon_currents[(split, end)],
                            result,
                            (1.0, 0.0),
                        )
                    )
                quark_without_z[end] = result

            for end in range(0, len(gluon_labels) + 1):
                result = CurrentKey(
                    pdg=quark_start.pdg,
                    external_labels=quark_start.external_labels
                    + gluon_labels[:end]
                    + z_key.external_labels,
                    chirality=chirality,
                )
                currents.append(result)
                interactions.append(
                    InteractionNode(
                        10,
                        quark_without_z[end],
                        z_key,
                        result,
                        self.model.z_fermion_coupling(abs(process[0])),
                    )
                )
                for split in range(0, end):
                    interactions.append(
                        InteractionNode(
                            6,
                            quark_with_z[split],
                            gluon_currents[(split, end)],
                            result,
                            (1.0, 0.0),
                        )
                    )
                quark_with_z[end] = result
            amplitudes.append((quark_with_z[len(gluon_labels)], anti_closure))

        return RecursionGraph(
            process=process,
            color_order=tuple(index + 1 for index in color_order),
            currents=tuple(currents),
            interactions=tuple(interactions),
            amplitudes=tuple(amplitudes),
        )


def _safe_process_name(process: str) -> str:
    return (
        process.replace(" ", "_")
        .replace(">", "to")
        .replace("~", "bar")
        .replace("+", "p")
        .replace("-", "m")
    )


def evaluator_artifact_path(cache_dir: str | Path, process: str) -> Path:
    return _cache_files(Path(cache_dir), process)[1]


def load_evaluator_artifact(cache_dir: str | Path, process: str) -> EvaluatorArtifact:
    path = evaluator_artifact_path(cache_dir, process)
    start = time.perf_counter()
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"invalid evaluator artifact payload in {path}")
    return EvaluatorArtifact(
        file=path,
        payload=payload,
        load_time_s=time.perf_counter() - start,
    )


def _cache_files(cache_dir: Path, process: str) -> tuple[Path, Path]:
    stem = _safe_process_name(process)
    return (
        cache_dir / f"{stem}.metadata.json",
        cache_dir / f"{stem}.evaluator.json",
    )


def _with_cache_files(
    result: GenerationResult,
    *,
    cache_file: Path,
    artifact_file: Path,
    artifact_fingerprint: str,
    artifact_cache_hit: bool,
) -> GenerationResult:
    return GenerationResult(
        process=result.process,
        backend=result.backend,
        supported_native_target=result.supported_native_target,
        cache_file=cache_file,
        artifact_file=artifact_file,
        artifact_fingerprint=artifact_fingerprint,
        artifact_cache_hit=artifact_cache_hit,
        generation_time_s=result.generation_time_s,
        graph=result.graph,
        symbolic_lowering=result.symbolic_lowering,
        notes=result.notes,
    )


def _result_fingerprint(result: GenerationResult) -> str:
    payload = result.to_json_dict()
    payload["_artifact_kernel"] = _artifact_kernel(result)
    payload["_symbolic_evaluator_schema_version"] = (
        1 if _is_zero_gluon_z_graph(result.graph) else None
    )
    payload["_tensor_network_evaluator_schema_version"] = (
        1 if _has_z_gluon_tensor_network_evaluator(result) else None
    )
    _strip_unstable_fingerprint_fields(payload)
    for key in (
        "artifact_cache_hit",
        "artifact_file",
        "artifact_fingerprint",
        "cache_file",
        "generation_time_s",
    ):
        payload.pop(key, None)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _strip_unstable_fingerprint_fields(payload: dict[str, object]) -> None:
    symbolic_lowering = payload.get("symbolic_lowering")
    if not isinstance(symbolic_lowering, dict):
        return
    blueprint = symbolic_lowering.get("tensor_network_blueprint")
    if isinstance(blueprint, dict):
        blueprint["execution_time_s"] = None


def _artifact_matches(path: Path, fingerprint: str) -> bool:
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return False
    return (
        isinstance(payload, dict)
        and payload.get("artifact_kind") == "pyamplicol-evaluator"
        and payload.get("schema_version") == 1
        and payload.get("artifact_fingerprint") == fingerprint
    )


def _evaluator_artifact_payload(
    result: GenerationResult,
    fingerprint: str,
) -> dict[str, object]:
    payload = result.to_json_dict()
    graph = payload.get("graph")
    symbolic_lowering = payload.get("symbolic_lowering")
    symbolic_scalar_evaluator = (
        build_zero_gluon_symbolic_evaluator_payload()
        if _is_zero_gluon_z_graph(result.graph)
        else None
    )
    tensor_network_scalar_evaluator = _tensor_network_scalar_evaluator_payload(result)
    return {
        "artifact_kind": "pyamplicol-evaluator",
        "schema_version": 1,
        "process": result.process,
        "backend": result.backend,
        "kernel": _artifact_kernel(result),
        "supported_native_target": result.supported_native_target,
        "artifact_fingerprint": fingerprint,
        "full_me_tensor_network_ready": (
            tensor_network_scalar_evaluator is not None
            or _full_me_tensor_network_ready(symbolic_lowering)
        ),
        "symbolic_scalar_evaluator_ready": symbolic_scalar_evaluator is not None,
        "symbolic_scalar_evaluator": symbolic_scalar_evaluator,
        "tensor_network_scalar_evaluator_ready": (
            tensor_network_scalar_evaluator is not None
        ),
        "tensor_network_scalar_evaluator": tensor_network_scalar_evaluator,
        "graph_counts": _graph_counts(result.graph),
        "graph": graph,
        "symbolic_lowering": symbolic_lowering,
        "runtime_inputs": {
            "kinematics": "runtime",
            "external_wavefunctions": "runtime-boundary-cacheable",
        },
        "notes": list(result.notes),
    }


def _artifact_kernel(result: GenerationResult) -> str:
    if _is_zero_gluon_z_graph(result.graph):
        return "symbolica-zero-gluon"
    if _has_z_gluon_tensor_network_evaluator(result):
        gluon_count = _z_gluon_count(result.graph)
        if gluon_count == 1:
            return "symbolica-one-gluon-tensor-network"
        return "symbolica-z-gluon-tensor-network"
    if result.backend == "native-python-recursion-staged":
        return "staged-python-recursion"
    if result.backend == "native-symbolic-pending":
        return "symbolic-lowering-pending"
    return "not-implemented"


def _tensor_network_scalar_evaluator_payload(
    result: GenerationResult,
) -> dict[str, object] | None:
    if not _has_z_gluon_tensor_network_evaluator(result):
        return None
    from .tensor_runtime import build_z_gluon_tensor_network_evaluator_payload

    if result.graph is None:
        return None
    return build_z_gluon_tensor_network_evaluator_payload(
        result.process,
        result.graph,
    )


def _has_z_gluon_tensor_network_evaluator(result: GenerationResult) -> bool:
    gluon_count = _z_gluon_count(result.graph)
    return (
        gluon_count is not None
        and gluon_count <= _AUTO_TENSOR_NETWORK_ARTIFACT_MAX_GLUONS
        and result.symbolic_lowering is not None
        and _symbolic_lowering_supports_tensor_network_evaluator(result)
    )


def _symbolic_lowering_supports_tensor_network_evaluator(
    result: GenerationResult,
) -> bool:
    if result.symbolic_lowering is None:
        return False
    vertex_lowering = result.symbolic_lowering.vertex_lowering
    blueprint = result.symbolic_lowering.tensor_network_blueprint
    return (
        vertex_lowering is not None
        and vertex_lowering.full_tensor_network_ready is True
        and blueprint is not None
        and blueprint.propagator_lowering_ready is True
        and blueprint.pending_interactions == 0
    )


def _full_me_tensor_network_ready(symbolic_lowering: object) -> bool:
    return (
        isinstance(symbolic_lowering, dict)
        and symbolic_lowering.get("full_me_tensor_network_ready") is True
    )


def _graph_counts(graph: RecursionGraph | None) -> dict[str, int] | None:
    if graph is None:
        return None
    return {
        "currents": len(graph.currents),
        "interactions": len(graph.interactions),
        "amplitudes": len(graph.amplitudes),
    }


def _is_zero_gluon_z_graph(graph: RecursionGraph | None) -> bool:
    return (
        graph is not None
        and len(graph.process) == 3
        and graph.process[2] == 23
        and all(pdg != 21 for pdg in graph.process[2:])
    )


def _z_gluon_count(graph: RecursionGraph | None) -> int | None:
    if graph is None or len(graph.process) < 4 or graph.process[-1] != 23:
        return None
    final_without_z = graph.process[2:-1]
    if not final_without_z or any(pdg != 21 for pdg in final_without_z):
        return None
    return len(final_without_z)


def _pdg_process(process: Sequence[str]) -> tuple[str, ...]:
    from .processes import ANTI_PARTICLE, PDGS

    return tuple(PDGS[p] if i > 1 else PDGS[ANTI_PARTICLE[p]] for i, p in enumerate(process))


def _external_current_type(pdg: int) -> int:
    if 1 <= abs(pdg) <= 6:
        return -pdg
    return pdg


__all__ = [
    "CurrentKey",
    "EvaluatorArtifact",
    "GenerationResult",
    "InteractionNode",
    "NativeMatrixElementGenerator",
    "RecursionGraph",
    "evaluator_artifact_path",
    "load_evaluator_artifact",
]
