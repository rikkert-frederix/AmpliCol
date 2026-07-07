from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Mapping, TypedDict

from .color_plan import GenericColorPlan, build_color_plan
from .generic_dag import GenericDAG, GenericDAGCompiler
from .model import AmplicolSMLeadingColorModel
from .process_ir import CanonicalProcessIR, build_process_ir
from .processes import (
    PDGS,
    SINGLETS,
    ParsedProcess,
    ProcessOptions,
    ProcessTuple,
)


_EW_VECTOR_NAMES = frozenset({"a", "z", "w+", "w-"})


@dataclass(frozen=True)
class ProcessContent:
    process: str
    key: str
    parsed: ParsedProcess
    ir: CanonicalProcessIR
    all_outgoing: ProcessTuple
    initial_pdgs: tuple[int, ...]
    final_pdgs: tuple[int, ...]
    quark_count: int
    antiquark_count: int
    quark_pair_count: int
    gluon_count: int
    final_singlets: ProcessTuple
    final_vectors: ProcessTuple
    final_higgs_count: int
    final_leptons: tuple[int, ...]

    @property
    def has_inclusive_initial_state(self) -> bool:
        return self.ir.has_inclusive_initial_state

    @property
    def has_multiple_nonleptonic_singlets(self) -> bool:
        return self.ir.has_multiple_nonleptonic_singlets


@dataclass(frozen=True)
class CurrentPlanSupportSummary:
    current_count: int
    source_count: int
    interaction_count: int
    closure_count: int
    color_sectors: tuple[int, ...]
    truncated: bool
    full_tensor_network_ready: bool
    required_vertex_kind_counts: tuple[tuple[int, int], ...]
    ready_vertex_kinds: tuple[int, ...]
    pending_vertex_kinds: tuple[int, ...]
    unimplemented_vertex_kinds: tuple[int, ...]
    required_propagator_kernel_counts: tuple[tuple[str, int], ...]
    ready_propagator_kernels: tuple[str, ...]
    pending_propagator_kernels: tuple[str, ...]
    unimplemented_propagator_kernels: tuple[str, ...]

    @property
    def has_closure(self) -> bool:
        return self.closure_count > 0

    def diagnostic_suffix(self) -> str:
        parts: list[str] = []
        if self.truncated:
            parts.append(
                f"generic current planning truncated at {self.current_count} currents"
            )
        if not self.has_closure:
            parts.append("generic current planning did not find amplitude closures")
        missing_parts: list[str] = []
        if self.pending_vertex_kinds:
            missing_parts.append(
                f"pending vertex kinds {list(self.pending_vertex_kinds)}"
            )
        if self.unimplemented_vertex_kinds:
            missing_parts.append(
                f"unimplemented vertex kinds {list(self.unimplemented_vertex_kinds)}"
            )
        if self.pending_propagator_kernels:
            missing_parts.append(
                "pending propagator kernels "
                f"{list(self.pending_propagator_kernels)}"
            )
        if self.unimplemented_propagator_kernels:
            missing_parts.append(
                "unimplemented propagator kernels "
                f"{list(self.unimplemented_propagator_kernels)}"
            )
        if missing_parts:
            parts.append(
                "generic current planning found closures but still needs "
                + " and ".join(missing_parts)
            )
        if not parts:
            parts.append(
                "generic current planning found a tensor-ready current graph; "
                "serialized schema-v2 evaluator emission and Rusticol runtime "
                "execution are available"
            )
        return "; ".join(parts)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "current_count": self.current_count,
            "source_count": self.source_count,
            "interaction_count": self.interaction_count,
            "closure_count": self.closure_count,
            "color_sectors": list(self.color_sectors),
            "color_sector_count": len(self.color_sectors),
            "truncated": self.truncated,
            "full_tensor_network_ready": self.full_tensor_network_ready,
            "required_vertex_kind_counts": [
                [kind, count] for kind, count in self.required_vertex_kind_counts
            ],
            "ready_vertex_kinds": list(self.ready_vertex_kinds),
            "pending_vertex_kinds": list(self.pending_vertex_kinds),
            "unimplemented_vertex_kinds": list(self.unimplemented_vertex_kinds),
            "required_propagator_kernel_counts": [
                [kernel, count]
                for kernel, count in self.required_propagator_kernel_counts
            ],
            "ready_propagator_kernels": list(self.ready_propagator_kernels),
            "pending_propagator_kernels": list(self.pending_propagator_kernels),
            "unimplemented_propagator_kernels": list(
                self.unimplemented_propagator_kernels
            ),
        }


class _PropagatorSupportSummary(TypedDict):
    required_propagator_kernel_counts: tuple[tuple[str, int], ...]
    ready_propagator_kernels: tuple[str, ...]
    pending_propagator_kernels: tuple[str, ...]
    unimplemented_propagator_kernels: tuple[str, ...]


@dataclass(frozen=True)
class ColorPlanSupportSummary:
    sector_count: int
    truncated: bool
    ready_for_leading_colour: bool
    idenso_required: bool
    diagnostics: tuple[str, ...]
    sector_kind_counts: tuple[tuple[str, int], ...]

    def blocking_diagnostic(self) -> str | None:
        if self.idenso_required:
            return "colour expansion requires Idenso basis/metric generation"
        if self.truncated:
            return f"colour planning truncated at {self.sector_count} sectors"
        if not self.sector_count:
            return "colour planning found no leading-colour sectors"
        if not self.ready_for_leading_colour:
            return "colour planning is not ready for leading colour"
        return None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "sector_count": self.sector_count,
            "truncated": self.truncated,
            "ready_for_leading_colour": self.ready_for_leading_colour,
            "idenso_required": self.idenso_required,
            "diagnostics": list(self.diagnostics),
            "sector_kind_counts": [
                [kind, count] for kind, count in self.sector_kind_counts
            ],
        }


@dataclass(frozen=True)
class ProcessSupportReport:
    process: str
    color_accuracy: str
    runtime_artifact_supported: bool
    generic_dag_runtime_supported: bool
    support_class: str
    missing_feature: str | None
    artifact_unavailable_message: str | None
    content: ProcessContent | None = None
    color_plan: ColorPlanSupportSummary | None = None
    current_plan: CurrentPlanSupportSummary | None = None

    def to_json_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {
            "runtime_artifact_supported": self.runtime_artifact_supported,
            "generic_dag_runtime_supported": self.generic_dag_runtime_supported,
            "process": self.process,
            "color_accuracy": self.color_accuracy,
            "support_class": self.support_class,
            "missing_feature": self.missing_feature,
            "artifact_unavailable_message": self.artifact_unavailable_message,
        }
        if self.content is not None:
            payload["content"] = {
                "key": self.content.key,
                "quark_pair_count": self.content.quark_pair_count,
                "gluon_count": self.content.gluon_count,
                "final_singlets": list(self.content.final_singlets),
                "final_vectors": list(self.content.final_vectors),
                "final_higgs_count": self.content.final_higgs_count,
                "final_leptons": list(self.content.final_leptons),
                "labels": {
                    "gluons": list(self.content.ir.gluon_labels),
                    "vectors": list(self.content.ir.vector_labels),
                    "leptons": list(self.content.ir.lepton_labels),
                    "higgs": list(self.content.ir.higgs_labels),
                },
            }
        if self.color_plan is not None:
            payload["color_plan"] = self.color_plan.to_json_dict()
        if self.current_plan is not None:
            payload["current_plan"] = self.current_plan.to_json_dict()
        return payload


def classify_process_support(
    process: str,
    *,
    color_accuracy: str = "lc",
    options: ProcessOptions | None = None,
    include_color_plan: bool = True,
    color_plan_max_sectors: int = 5000,
    include_current_plan: bool = True,
    current_plan_max_currents: int = 50000,
    selected_color_sector_ids: Iterable[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Iterable[int] | None = None,
    ignored_vertex_kinds: Iterable[int] | None = None,
) -> ProcessSupportReport:
    """Classify generic-DAG production readiness for one concrete process.

    This deliberately does not build a matrix element. It performs the
    production preflight required by the generic architecture: parse the
    process, build colour sectors, build the model-driven current DAG, then
    report the first missing layer. It must not classify support by matching
    whole process families.
    """

    try:
        content = describe_process_content(
            process,
            color_accuracy=color_accuracy,
            options=options,
        )
    except ValueError as exc:
        return ProcessSupportReport(
            process=process,
            color_accuracy=color_accuracy,
            runtime_artifact_supported=False,
            generic_dag_runtime_supported=False,
            support_class="invalid-process",
            missing_feature="process-parser",
            artifact_unavailable_message=str(exc),
        )

    color_plan = (
        _color_plan_summary(
            content,
            color_accuracy=color_accuracy,
            max_sectors=color_plan_max_sectors,
        )
        if include_color_plan
        else None
    )
    current_plan = (
        _current_plan_summary(
            content,
            max_currents=current_plan_max_currents,
            selected_color_sector_ids=selected_color_sector_ids,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
        )
        if include_current_plan
        else None
    )
    return _generic_preflight_report(
        process,
        content,
        color_plan=color_plan,
        current_plan=current_plan,
        selected_color_sector_ids=selected_color_sector_ids,
    )


def describe_process_content(
    process: str,
    *,
    color_accuracy: str = "lc",
    options: ProcessOptions | None = None,
) -> ProcessContent:
    ir = build_process_ir(process, color_accuracy=color_accuracy, options=options)
    parsed = ir.parsed
    all_outgoing = ir.outgoing_particles
    final_singlets = tuple(particle for particle in parsed.rest if particle in SINGLETS)
    final_vectors = tuple(
        particle for particle in parsed.rest if particle in _EW_VECTOR_NAMES
    )
    final_leptons = tuple(
        int(PDGS[particle])
        for particle in parsed.rest
        if _is_lepton_name(particle)
    )
    return ProcessContent(
        process=ir.process,
        key=ir.key,
        parsed=parsed,
        ir=ir,
        all_outgoing=tuple(all_outgoing),
        initial_pdgs=ir.initial_pdgs,
        final_pdgs=ir.final_pdgs,
        quark_count=ir.quark_lines.quark_count,
        antiquark_count=ir.quark_lines.antiquark_count,
        quark_pair_count=ir.quark_lines.quark_pair_count,
        gluon_count=len(ir.gluon_labels),
        final_singlets=final_singlets,
        final_vectors=final_vectors,
        final_higgs_count=len(ir.higgs_labels),
        final_leptons=final_leptons,
    )


def _supported_report(
    process: str,
    content: ProcessContent,
    *,
    support_class: str,
    color_plan: ColorPlanSupportSummary | None = None,
    current_plan: CurrentPlanSupportSummary | None = None,
) -> ProcessSupportReport:
    return ProcessSupportReport(
        process=process,
        color_accuracy=content.ir.color_accuracy,
        runtime_artifact_supported=True,
        generic_dag_runtime_supported=True,
        support_class=support_class,
        missing_feature=None,
        artifact_unavailable_message=None,
        content=content,
        color_plan=color_plan,
        current_plan=current_plan,
    )


def _unsupported_report(
    process: str,
    content: ProcessContent,
    *,
    support_class: str,
    missing_feature: str,
    message: str,
    generic_dag_runtime_supported: bool = True,
    color_plan: ColorPlanSupportSummary | None = None,
    current_plan: CurrentPlanSupportSummary | None = None,
) -> ProcessSupportReport:
    if color_plan is not None:
        color_diagnostic = color_plan.blocking_diagnostic()
        if color_diagnostic is not None:
            message = f"{message} ({color_diagnostic})"
    if current_plan is not None:
        current_diagnostic = current_plan.diagnostic_suffix()
        if current_diagnostic not in message:
            message = f"{message} ({current_diagnostic})"
    return ProcessSupportReport(
        process=process,
        color_accuracy=content.ir.color_accuracy,
        runtime_artifact_supported=False,
        generic_dag_runtime_supported=generic_dag_runtime_supported,
        support_class=support_class,
        missing_feature=missing_feature,
        artifact_unavailable_message=message,
        content=content,
        color_plan=color_plan,
        current_plan=current_plan,
    )


def _generic_preflight_report(
    process: str,
    content: ProcessContent,
    *,
    color_plan: ColorPlanSupportSummary | None,
    current_plan: CurrentPlanSupportSummary | None,
    selected_color_sector_ids: Iterable[int] | None,
) -> ProcessSupportReport:
    if content.has_inclusive_initial_state:
        return _unsupported_report(
            process,
            content,
            support_class="inclusive-process-request",
            missing_feature="process-set-expansion",
            generic_dag_runtime_supported=False,
            message=(
                "Rusticol process artifacts require concrete partonic subprocesses; "
                "expand inclusive labels through generate-process/process-set first."
            ),
            color_plan=color_plan,
            current_plan=None,
        )
    color_diagnostic = (
        color_plan.blocking_diagnostic() if color_plan is not None else None
    )
    if (
        color_diagnostic is not None
        and not _selected_lc_sector_truncation_is_safe(
            color_plan,
            current_plan,
            selected_color_sector_ids=selected_color_sector_ids,
        )
    ):
        return _unsupported_report(
            process,
            content,
            support_class="generic-dag-colour-preflight",
            missing_feature=(
                "colour-expansion"
                if color_plan is not None and color_plan.idenso_required
                else "colour-planning"
            ),
            generic_dag_runtime_supported=False,
            message=(
                f"{color_diagnostic} for "
                f"--color-accuracy={content.ir.color_accuracy}"
            ),
            color_plan=color_plan,
            current_plan=current_plan,
        )
    if current_plan is None:
        return _unsupported_report(
            process,
            content,
            support_class="generic-dag-current-preflight",
            missing_feature="generic-dag-current-planning",
            generic_dag_runtime_supported=False,
            message=(
                "generic model-driven current planning failed before an artifact "
                "could be built"
            ),
            color_plan=color_plan,
            current_plan=current_plan,
        )
    if current_plan.truncated:
        return _unsupported_report(
            process,
            content,
            support_class="generic-dag-current-preflight",
            missing_feature="generic-dag-current-cap",
            generic_dag_runtime_supported=False,
            message=current_plan.diagnostic_suffix(),
            color_plan=color_plan,
            current_plan=current_plan,
        )
    if not current_plan.has_closure:
        return _unsupported_report(
            process,
            content,
            support_class="generic-dag-current-preflight",
            missing_feature="generic-dag-amplitude-roots",
            generic_dag_runtime_supported=False,
            message=current_plan.diagnostic_suffix(),
            color_plan=color_plan,
            current_plan=current_plan,
        )
    if (
        current_plan.pending_vertex_kinds
        or current_plan.unimplemented_vertex_kinds
        or current_plan.pending_propagator_kernels
        or current_plan.unimplemented_propagator_kernels
    ):
        vertex_blocked = (
            bool(current_plan.pending_vertex_kinds)
            or bool(current_plan.unimplemented_vertex_kinds)
        )
        return _unsupported_report(
            process,
            content,
            support_class="generic-dag-lowering-preflight",
            missing_feature=(
                "generic-dag-vertex-lowering"
                if vertex_blocked
                else "generic-dag-propagator-lowering"
            ),
            generic_dag_runtime_supported=False,
            message=current_plan.diagnostic_suffix(),
            color_plan=color_plan,
            current_plan=current_plan,
        )
    return _supported_report(
        process,
        content,
        support_class="generic-dag-schema-v2",
        color_plan=color_plan,
        current_plan=current_plan,
    )


def _selected_lc_sector_truncation_is_safe(
    color_plan: ColorPlanSupportSummary | None,
    current_plan: CurrentPlanSupportSummary | None,
    *,
    selected_color_sector_ids: Iterable[int] | None,
) -> bool:
    if color_plan is None or current_plan is None:
        return False
    if not color_plan.truncated or color_plan.idenso_required:
        return False
    selected = {int(sector_id) for sector_id in (selected_color_sector_ids or ())}
    if not selected:
        return False
    if any(sector_id < 0 or sector_id >= color_plan.sector_count for sector_id in selected):
        return False
    if current_plan.truncated or not current_plan.has_closure:
        return False
    if not current_plan.full_tensor_network_ready:
        return False
    return set(current_plan.color_sectors).issubset(selected)


def _color_plan_summary(
    content: ProcessContent,
    *,
    color_accuracy: str,
    max_sectors: int,
) -> ColorPlanSupportSummary | None:
    try:
        plan = build_color_plan(
            content.ir,
            color_accuracy=color_accuracy,
            max_sectors=max_sectors,
        )
    except (KeyError, ValueError, RuntimeError):
        return None
    return _color_plan_summary_from_plan(plan)


def _color_plan_summary_from_plan(
    plan: GenericColorPlan,
) -> ColorPlanSupportSummary:
    sector_kind_counts = tuple(
        sorted(Counter(sector.kind for sector in plan.sectors).items())
    )
    return ColorPlanSupportSummary(
        sector_count=plan.sector_count,
        truncated=plan.truncated,
        ready_for_leading_colour=plan.ready_for_leading_colour,
        idenso_required=plan.idenso_required,
        diagnostics=plan.diagnostics,
        sector_kind_counts=sector_kind_counts,
    )


def _current_plan_summary(
    content: ProcessContent,
    *,
    max_currents: int,
    selected_color_sector_ids: Iterable[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Iterable[int] | None = None,
    ignored_vertex_kinds: Iterable[int] | None = None,
) -> CurrentPlanSupportSummary | None:
    try:
        dag = GenericDAGCompiler(
            max_currents=max_currents,
            selected_color_sector_ids=selected_color_sector_ids,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
        ).compile(content.ir)
    except (KeyError, ValueError, RuntimeError):
        return None
    return _current_plan_summary_from_dag(dag)


def _current_plan_summary_from_dag(
    dag: GenericDAG,
) -> CurrentPlanSupportSummary:
    vertex_counts = Counter(
        [interaction.vertex_kind for interaction in dag.interactions]
        + [
            int(root.vertex_kind)
            for root in dag.amplitude_roots
            if root.vertex_kind is not None
        ]
    )
    ready: set[int] = set()
    pending: set[int] = set()
    unimplemented: set[int] = set()
    model = AmplicolSMLeadingColorModel()
    for kind in vertex_counts:
        rule = model.vertex_lowering_rule(kind)
        if rule.backend == "unimplemented":
            unimplemented.add(kind)
        elif rule.full_tensor_network_ready:
            ready.add(kind)
        else:
            pending.add(kind)
    propagators = _current_plan_propagator_summary(dag, model)
    color_sectors = tuple(
        sorted({current.index.color_state.sector_id for current in dag.currents})
    )
    return CurrentPlanSupportSummary(
        current_count=len(dag.currents),
        source_count=len(dag.sources),
        interaction_count=len(dag.interactions),
        closure_count=len(dag.amplitude_roots),
        color_sectors=color_sectors,
        truncated=dag.truncated,
        full_tensor_network_ready=(
            dag.has_amplitudes
            and not pending
            and not unimplemented
            and not propagators["pending_propagator_kernels"]
            and not propagators["unimplemented_propagator_kernels"]
            and not dag.truncated
        ),
        required_vertex_kind_counts=tuple(sorted(vertex_counts.items())),
        ready_vertex_kinds=tuple(sorted(ready)),
        pending_vertex_kinds=tuple(sorted(pending)),
        unimplemented_vertex_kinds=tuple(sorted(unimplemented)),
        required_propagator_kernel_counts=propagators[
            "required_propagator_kernel_counts"
        ],
        ready_propagator_kernels=propagators["ready_propagator_kernels"],
        pending_propagator_kernels=propagators["pending_propagator_kernels"],
        unimplemented_propagator_kernels=propagators[
            "unimplemented_propagator_kernels"
        ],
    )


def _current_plan_propagator_summary(
    dag: GenericDAG,
    model: AmplicolSMLeadingColorModel,
) -> _PropagatorSupportSummary:
    used_as_input = {
        interaction.left_id for interaction in dag.interactions
    } | {
        interaction.right_id for interaction in dag.interactions
    }
    counts: Counter[str] = Counter()
    ready: set[str] = set()
    pending: set[str] = set()
    unimplemented: set[str] = set()
    for current_id in sorted(used_as_input):
        current = dag.currents[current_id]
        if current.is_source:
            continue
        rule = model.propagator_lowering_rule(
            current.index.particle_id,
            current.index.chirality,
        )
        counts[rule.kernel] += 1
        if rule.backend == "unimplemented":
            unimplemented.add(rule.kernel)
        elif rule.full_tensor_network_ready:
            ready.add(rule.kernel)
        else:
            pending.add(rule.kernel)
    return {
        "required_propagator_kernel_counts": tuple(sorted(counts.items())),
        "ready_propagator_kernels": tuple(sorted(ready)),
        "pending_propagator_kernels": tuple(sorted(pending)),
        "unimplemented_propagator_kernels": tuple(sorted(unimplemented)),
    }


def _is_lepton_name(particle: str) -> bool:
    try:
        pdg = abs(int(PDGS[particle]))
    except KeyError:
        return False
    return 11 <= pdg <= 16


__all__ = [
    "ColorPlanSupportSummary",
    "CurrentPlanSupportSummary",
    "ProcessContent",
    "ProcessSupportReport",
    "classify_process_support",
    "describe_process_content",
]
