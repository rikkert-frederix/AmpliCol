from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal, Sequence

from .lowering import SymbolicLoweringReport, build_symbolic_lowering_report
from .model import AmplicolSMLeadingColorModel, Model
from .native import LeadingColorZJetsNativeEvaluator
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
        supported = self._is_supported_native_milestone(process)
        native = (
            LeadingColorZJetsNativeEvaluator(self.model)
            if isinstance(self.model, AmplicolSMLeadingColorModel)
            else None
        )
        dilepton_target = (
            None
            if native is None
            else native.supported_neutral_dilepton_gluon_process(process)
        )
        charged_leptonic_target = (
            None
            if native is None
            else native.supported_charged_leptonic_w_gluon_process(process)
        )
        supported_native_target = (
            supported
            or dilepton_target is not None
            or charged_leptonic_target is not None
        )
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
            enumeration = ProcessEnumerator(options).enumerate(process)
            record = _first_subprocess_record(enumeration)
            graph = self._build_staged_graph(
                tuple(int(pdg) for pdg in _pdg_process(record.process)),
                record.color_order,
            )
            symbolic_lowering = build_symbolic_lowering_report(self.model, graph=graph)
            gluon_count = sum(1 for pdg in graph.process[2:] if pdg == 21)
            vector_pdg = _neutral_vector_pdg(graph)
            backend = "native-python-recursion-staged"
            if vector_pdg == 23 and gluon_count <= _AUTO_TENSOR_NETWORK_ARTIFACT_MAX_GLUONS:
                notes.append(
                    "Native ordered Python recursion is available for the "
                    f"q q~ -> {_neutral_vector_name(vector_pdg)} + n g family; full spenso/Symbolica tensor-network "
                    "evaluator artifacts are currently generated through four "
                    "final-state gluons."
                )
            elif vector_pdg == 23:
                notes.append(
                    "Native ordered Python recursion is available for the "
                    f"q q~ -> {_neutral_vector_name(vector_pdg)} + n g family; full spenso/Symbolica tensor-network "
                    "blueprints are generated, while evaluator artifacts above "
                    "four final-state gluons remain guarded by expression-size "
                    "growth."
                )
            else:
                notes.append(
                    "Native ordered Python recursion and graph generation are "
                    f"available for q q~ -> {_neutral_vector_name(vector_pdg)} + n g; "
                    "compiled process artifacts still require the generic "
                    "neutral-vector source/runtime path."
                )
        elif dilepton_target is not None:
            backend = "native-python-recursion-staged"
            notes.append(
                "Native ordered Python recursion is available for the "
                "one-quark-line neutral dilepton + n g family, including "
                "gamma/Z interference and internal vector propagators. "
                "Shared-current graph lowering and reusable Rusticol artifacts "
                "for explicit lepton-source processes are still pending."
            )
        elif charged_leptonic_target is not None:
            backend = "native-python-recursion-staged"
            notes.append(
                "Native ordered Python recursion is available for the "
                "one-quark-line charged-current leptonic W + n g family, "
                "including internal W propagators. Shared-current graph lowering "
                "and reusable Rusticol artifacts for explicit lepton-source "
                "processes are still pending."
            )
        else:
            notes.append(
                "Full process architecture is available, but native matrix-element "
                "lowering is currently implemented only for one-quark-line "
                "electroweak vector + n g, neutral dilepton + n g, and "
                "charged-current leptonic W + n g families."
            )

        result = GenerationResult(
            process=process,
            backend=backend,
            supported_native_target=supported_native_target,
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
        zero_gluon_z = rest == ("z",)
        vector_names = tuple(name for name in ("a", "z", "w+", "w-") if rest.count(name) == 1)
        only_vector_gluons = (
            len(vector_names) == 1
            and all(p in {"a", "z", "w+", "w-", "g"} for p in rest)
            and rest.count("g") >= 1
        )
        if not (len(incoming) == 2 and (zero_gluon_z or only_vector_gluons)):
            return False

        from .processes import ANTI_PARTICLE, PDGS

        incoming_physical = tuple(
            int(PDGS[ANTI_PARTICLE[p]]) for p in incoming
        )
        if (
            any(not 1 <= abs(pdg) <= 6 for pdg in incoming_physical)
            or incoming_physical[0] * incoming_physical[1] >= 0
        ):
            return False
        if zero_gluon_z:
            return abs(incoming_physical[0]) == abs(incoming_physical[1])

        vector_pdg = int(PDGS[vector_names[0]])
        incoming_currents = tuple(-pdg for pdg in incoming_physical)
        quark_currents = [pdg for pdg in incoming_currents if pdg > 0]
        anti_closures = [pdg for pdg in incoming_currents if pdg < 0]
        if len(quark_currents) != 1 or len(anti_closures) != 1:
            return False
        try:
            vector_result = _electroweak_vector_result_pdg(
                quark_currents[0],
                vector_pdg,
            )
            return anti_closures[0] == -vector_result
        except ValueError:
            return False

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
        vector_labels = tuple(
            index + 1
            for index, pdg in enumerate(process)
            if index >= 2 and pdg in (22, 23, 24, -24)
        )
        if len(vector_labels) != 1:
            raise ValueError(
                "staged vector+gluon graph expects exactly one final-state photon, Z, W+, or W-"
            )
        vector_key = external[vector_labels[0] - 1]
        vector_pdg = int(process[vector_labels[0] - 1])
        quark_candidates = [key for key in external[:2] if key.pdg > 0]
        anti_candidates = [key for key in external[:2] if key.pdg < 0]
        if len(quark_candidates) != 1 or len(anti_candidates) != 1:
            raise ValueError(
                "staged vector+gluon graph expects one incoming quark and one antiquark"
            )
        quark_start = quark_candidates[0]
        anti_closure = anti_candidates[0]
        vector_result_pdg = _electroweak_vector_result_pdg(quark_start.pdg, vector_pdg)
        if anti_closure.pdg != -vector_result_pdg:
            raise ValueError(
                "incoming quark flavours are not compatible with the requested "
                f"vector emission: {quark_start.pdg} + {vector_pdg} -> "
                f"{vector_result_pdg}, but closure is {anti_closure.pdg}"
            )
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
                    pdg=vector_result_pdg,
                    external_labels=quark_start.external_labels
                    + gluon_labels[:end]
                    + vector_key.external_labels,
                    chirality=chirality,
                )
                currents.append(result)
                interactions.append(
                    InteractionNode(
                        10,
                        quark_without_z[end],
                        vector_key,
                        result,
                        _neutral_vector_coupling(
                            self.model,
                            vector_pdg=vector_pdg,
                            fermion_pdg=quark_start.pdg,
                        ),
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
        1
        if (
            _has_neutral_vector_gluon_tensor_network_evaluator(result)
            and result.graph is not None
            and _neutral_vector_pdg(result.graph) == 23
        )
        else None
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
    if (
        _has_neutral_vector_gluon_tensor_network_evaluator(result)
        and result.graph is not None
        and _neutral_vector_pdg(result.graph) == 23
    ):
        gluon_count = _neutral_vector_gluon_count(result.graph)
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
    if not _has_neutral_vector_gluon_tensor_network_evaluator(result):
        return None
    from .tensor_runtime import build_z_gluon_tensor_network_evaluator_payload

    if result.graph is None:
        return None
    if _neutral_vector_pdg(result.graph) != 23:
        return None
    return build_z_gluon_tensor_network_evaluator_payload(
        result.process,
        result.graph,
    )


def _has_neutral_vector_gluon_tensor_network_evaluator(result: GenerationResult) -> bool:
    gluon_count = _neutral_vector_gluon_count(result.graph)
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


def _neutral_vector_gluon_count(graph: RecursionGraph | None) -> int | None:
    if (
        graph is None
        or len(graph.process) < 4
        or graph.process[-1] not in (22, 23, 24, -24)
    ):
        return None
    final_without_vector = graph.process[2:-1]
    if not final_without_vector or any(pdg != 21 for pdg in final_without_vector):
        return None
    return len(final_without_vector)


def _neutral_vector_pdg(graph: RecursionGraph | None) -> int:
    if graph is None:
        raise ValueError("missing electroweak-vector graph")
    vectors = tuple(pdg for pdg in graph.process[2:] if pdg in (22, 23, 24, -24))
    if len(vectors) != 1:
        raise ValueError(f"expected one electroweak vector in graph, got {vectors}")
    return int(vectors[0])


def _neutral_vector_name(pdg: int) -> str:
    if pdg == 22:
        return "gamma"
    if pdg == 23:
        return "Z"
    if pdg == 24:
        return "W+"
    if pdg == -24:
        return "W-"
    return f"PDG {pdg}"


def _neutral_vector_coupling(
    model: AmplicolSMLeadingColorModel,
    *,
    vector_pdg: int,
    fermion_pdg: int,
) -> tuple[float, float]:
    if vector_pdg == 22:
        return model.photon_fermion_coupling(fermion_pdg)
    if vector_pdg == 23:
        return model.z_fermion_coupling(fermion_pdg)
    if abs(vector_pdg) == 24:
        _ = _electroweak_vector_result_pdg(fermion_pdg, vector_pdg)
        return (model.charged_current_coupling(), 0.0)
    raise ValueError(f"unsupported electroweak vector PDG: {vector_pdg}")


def _electroweak_vector_result_pdg(fermion_pdg: int, vector_pdg: int) -> int:
    if not 1 <= fermion_pdg <= 6:
        raise ValueError(
            f"charged-current quark-line support expects a quark, got {fermion_pdg}"
        )
    if vector_pdg in (22, 23):
        return fermion_pdg
    if vector_pdg == 24 and fermion_pdg in (1, 3, 5):
        return fermion_pdg + 1
    if vector_pdg == -24 and fermion_pdg in (2, 4, 6):
        return fermion_pdg - 1
    raise ValueError(f"unsupported charged-current transition {fermion_pdg} + {vector_pdg}")


def _pdg_process(process: Sequence[str]) -> tuple[str, ...]:
    from .processes import ANTI_PARTICLE, PDGS

    return tuple(PDGS[p] if i > 1 else PDGS[ANTI_PARTICLE[p]] for i, p in enumerate(process))


def _first_subprocess_record(enumeration: Any) -> Any:
    for group in enumeration.groups:
        if group.records:
            return group.records[0]
    raise ValueError("process enumeration did not produce subprocess records")


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
