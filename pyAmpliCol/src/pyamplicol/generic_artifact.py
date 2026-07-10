from __future__ import annotations

import json
import math
import shutil
import time
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping, Sequence, cast

from .color_plan import (
    GenericColorPlan,
    LCColorSectorReplayPartition,
    build_color_plan,
    lc_line_pairing_representative_ids,
    lc_topology_replay_partitions,
    lc_topology_replay_safe_groups,
)
from .color_contraction import (
    ColorGroupDescriptor,
    build_color_contraction_plan,
)
from .core_types import ExternalMomentum, NativeEvaluationError
from .generic_dag import (
    AmplitudeRoot,
    CurrentNode,
    GenericDAG,
    GenericDAGCompiler,
    InteractionNode,
    contributing_color_sector_ids,
    filter_dag_to_color_sectors,
    filter_dag_to_source_helicities,
    prune_global_helicity_flip_equivalent_roots,
    prune_dag_to_amplitude_roots,
)
from .model import AmplicolSMLeadingColorModel, Model
from .phase_space import massive_rambo_final_state
from .process_ir import CanonicalProcessIR
from .process_ir import ProcessSetIR, build_process_set_ir
from .processes import ProcessOptions

GENERIC_PROCESS_SCHEMA_VERSION = 2
GENERIC_PROCESS_MANIFEST_KIND = "pyamplicol-generic-process-plan"
GENERIC_PROCESS_SET_MANIFEST_KIND = "pyamplicol-generic-process-set-plan"
GENERIC_DAG_PROCESS_ARTIFACT_KIND = "pyamplicol-generic-dag-process"
GENERIC_DAG_PROCESS_SET_ARTIFACT_KIND = "pyamplicol-generic-dag-process-set"
GENERIC_LC_REPLAY_PARTITION_ARTIFACT_KIND = (
    "pyamplicol-lc-replay-partition-process"
)
LC_SECTOR_SELECTOR_PARAMETER = "runtime.lc_sector_id"

_FULL_COLOR_PLAN_SERIALIZATION_SECTOR_LIMIT = 1000
_NUMERICAL_REWRITE_VALIDATION_SAMPLE_LIMIT = 3
_NUMERICAL_REWRITE_VALIDATION_SEED_OFFSET = 1_000_003
_NUMERICAL_REWRITE_VALIDATION_RELATIVE_TOLERANCE = 1.0e-10
_NUMERICAL_REWRITE_VALIDATION_ZERO_TOLERANCE = 1.0e-24


@dataclass(frozen=True)
class GenericProcessManifest:
    """Process-generic planning manifest for the future schema-v2 runtime."""

    dag: GenericDAG
    model: Model
    color_plan: GenericColorPlan
    structural_current_aggregation: Mapping[str, object] | None = None
    zero_current_filter: Mapping[str, object] | None = None
    current_merging: Mapping[str, object] | None = None

    @property
    def process(self) -> str:
        return self.dag.process.process

    @property
    def key(self) -> str:
        return self.dag.process.key

    @property
    def external_pdg_order(self) -> tuple[int, ...]:
        return (*self.dag.process.initial_pdgs, *self.dag.process.final_pdgs)

    @property
    def outgoing_pdg_order(self) -> tuple[int, ...]:
        return self.dag.process.outgoing_pdgs

    def to_json_dict(
        self,
        *,
        selected_color_sector_ids: set[int] | None = None,
        enable_lc_sector_runtime_selector: bool = True,
    ) -> dict[str, object]:
        lowering = _dag_lowering_status(self.dag, self.model)
        current_ready = (
            self.dag.has_amplitudes
            and not lowering["pending_vertex_kinds"]
            and not lowering["unimplemented_vertex_kinds"]
            and not lowering["pending_propagator_kernels"]
            and not lowering["unimplemented_propagator_kernels"]
            and not self.dag.truncated
        )
        color_ready = self.color_plan.ready_for_requested_colour
        return {
            "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
            "kind": GENERIC_PROCESS_MANIFEST_KIND,
            "process": self.process,
            "key": self.key,
            "color_accuracy": self.dag.process.color_accuracy,
            "model": {
                "name": self.model.name,
                "vertex_lowering_coverage": (
                    self.model.vertex_lowering_coverage().to_json_dict()
                ),
            },
            "external_pdg_order": list(self.external_pdg_order),
            "outgoing_pdg_order": list(self.outgoing_pdg_order),
            "process_ir": self.dag.process.to_json_dict(),
            "color_plan": self.color_plan.to_json_dict(),
            "lc_topology_reuse": _lc_topology_reuse_payload(self),
            "generation_filters": {
                "structural_current_aggregation": dict(
                    self.structural_current_aggregation
                    or {
                        "enabled": False,
                        "reason": "structural current aggregation was not requested",
                    }
                ),
                "zero_current": dict(
                    self.zero_current_filter
                    or {
                        "enabled": False,
                        "reason": "zero-current warmup filter was not requested",
                    }
                ),
                "current_merging": dict(
                    self.current_merging
                    or {
                        "enabled": False,
                        "reason": "numerical current merging was not requested",
                    }
                ),
            },
            "currents": [
                current.to_json_dict() for current in self.dag.currents
            ],
            "sources": list(self.dag.sources),
            "interactions": [
                interaction.to_json_dict()
                for interaction in self.dag.interactions
            ],
            "amplitude_roots": [
                root.to_json_dict() for root in self.dag.amplitude_roots
            ],
            "stage_plan": _dag_stage_plan_payload(self.dag, self.model),
            "runtime_schema": _generic_runtime_schema_payload(
                self.dag,
                self.model,
                selected_color_sector_ids=selected_color_sector_ids,
                enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
            ),
            "planning_status": {
                "color_ready": color_ready,
                "color_sector_count": self.color_plan.sector_count,
                "color_truncated": self.color_plan.truncated,
                "idenso_required": self.color_plan.idenso_required,
                "current_ready": current_ready,
                "has_closure": self.dag.has_amplitudes,
                "has_amplitude_roots": self.dag.has_amplitudes,
                "generic_evaluator_ready": color_ready and current_ready,
            },
            "lowering_status": {
                **lowering,
                "closure_count": len(self.dag.amplitude_roots),
                "amplitude_root_count": len(self.dag.amplitude_roots),
                "truncated": self.dag.truncated,
                "has_closure": self.dag.has_amplitudes,
                "has_amplitude_roots": self.dag.has_amplitudes,
                "full_tensor_network_ready": current_ready,
            },
        }


@dataclass(frozen=True)
class GenericProcessSetManifest:
    """Process-set planning manifest with one generic plan per subprocess."""

    process_set: ProcessSetIR
    processes: tuple[GenericProcessManifest, ...]
    generation_metadata: dict[str, object] = field(default_factory=dict)

    @property
    def request(self) -> str:
        return self.process_set.request

    @property
    def color_accuracy(self) -> str:
        return self.process_set.color_accuracy

    @property
    def default_key(self) -> str:
        return self.process_set.default_key

    def to_json_dict(self) -> dict[str, object]:
        entries = []
        for manifest in self.processes:
            payload = manifest.to_json_dict()
            entries.append(
                {
                    "key": manifest.key,
                    "process": manifest.process,
                    "path": (
                        f"subprocesses/{manifest.key}/"
                        "generic_process_manifest.json"
                    ),
                    "generation_request": self.generation_metadata,
                    "planning_status": payload["planning_status"],
                    "lowering_status": payload["lowering_status"],
                }
            )
        return {
            "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
            "kind": GENERIC_PROCESS_SET_MANIFEST_KIND,
            "request": self.request,
            "default_process_key": self.default_key,
            "color_accuracy": self.color_accuracy,
            "generic_generation": self.generation_metadata,
            "processes": entries,
            "process_set_ir": self.process_set.to_json_dict(),
        }


def build_generic_process_manifest(
    process: str | CanonicalProcessIR | GenericDAG,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    color_accuracy: str = "lc",
    max_currents: int | None = None,
    max_color_sectors: int | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    selected_source_helicities: Mapping[int, int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
    progress_callback: Any | None = None,
) -> GenericProcessManifest:
    model = model or AmplicolSMLeadingColorModel()
    _emit_generation_progress(
        progress_callback,
        "dag compile",
        _process_progress_label(process),
    )
    phase_started = time.perf_counter()
    dag = (
        prune_dag_to_amplitude_roots(process)
        if isinstance(process, GenericDAG)
        else prune_dag_to_amplitude_roots(
            GenericDAGCompiler(
                model=model,
                color_accuracy=(
                    process.color_accuracy
                    if isinstance(process, CanonicalProcessIR)
                    else color_accuracy
                ),
                options=options,
                max_currents=max_currents,
                max_color_sectors=max_color_sectors,
                reference_color_order=(
                    None
                    if reference_color_order is None
                    else tuple(int(label) for label in reference_color_order)
                ),
                selected_color_sector_ids=selected_color_sector_ids,
                max_coupling_orders=max_coupling_orders,
                max_lc_current_line_groups=max_lc_current_line_groups,
                max_quark_pairs=max_quark_pairs,
                closure_side_mask_pruning=closure_side_mask_pruning,
                color_order_mask_pruning=color_order_mask_pruning,
                species_reachability_pruning=species_reachability_pruning,
                ignored_particle_ids=ignored_particle_ids,
                ignored_vertex_kinds=ignored_vertex_kinds,
            ).compile(process)
        )
    )
    _emit_generation_progress(
        progress_callback,
        "dag ready",
        _dag_progress_label(dag),
        increment=1,
        duration_s=time.perf_counter() - phase_started,
    )
    if selected_source_helicities:
        _emit_generation_progress(
            progress_callback,
            "helicity select",
            ",".join(
                f"{int(label)}:{int(helicity)}"
                for label, helicity in sorted(selected_source_helicities.items())
            ),
        )
        phase_started = time.perf_counter()
        dag = filter_dag_to_source_helicities(dag, selected_source_helicities)
        _emit_generation_progress(
            progress_callback,
            "helicity selected",
            _dag_progress_label(dag),
            increment=1,
            duration_s=time.perf_counter() - phase_started,
        )
    _emit_generation_progress(
        progress_callback,
        "helicity prune",
        "global flip representatives",
    )
    phase_started = time.perf_counter()
    dag = prune_global_helicity_flip_equivalent_roots(dag, model)
    _emit_generation_progress(
        progress_callback,
        "helicity ready",
        _dag_progress_label(dag),
        increment=1,
        duration_s=time.perf_counter() - phase_started,
    )
    _emit_generation_progress(
        progress_callback,
        "current aggregate",
        "lc gluon-flavour flow",
    )
    phase_started = time.perf_counter()
    dag, structural_current_aggregation = _aggregate_lc_gluon_flavour_flow_currents(
        dag,
        model,
        sample_count=numerical_current_samples,
        seed=numerical_current_seed,
    )
    _emit_generation_progress(
        progress_callback,
        "aggregate ready",
        _dag_progress_label(dag),
        increment=1,
        duration_s=time.perf_counter() - phase_started,
    )
    zero_current_filter: Mapping[str, object] | None = None
    if numerical_filter_current:
        _emit_generation_progress(
            progress_callback,
            "zero filter",
            f"{numerical_current_samples} samples",
        )
        phase_started = time.perf_counter()
        dag, zero_current_filter = _filter_zero_currents_by_warmup(
            dag,
            model,
            sample_count=numerical_current_samples,
            seed=numerical_current_seed,
            relative_tolerance=numerical_current_relative_tolerance,
            zero_tolerance=numerical_current_zero_tolerance,
        )
        _emit_generation_progress(
            progress_callback,
            "zero ready",
            _dag_progress_label(dag),
            increment=1,
            duration_s=time.perf_counter() - phase_started,
        )
    else:
        zero_current_filter = {
            "enabled": False,
            "reason": "disabled by --no-numerical-filter-current",
        }
        _emit_generation_progress(
            progress_callback,
            "zero filter",
            "disabled",
            increment=1,
            duration_s=0.0,
        )
    current_merging: Mapping[str, object] | None = None
    if numerical_current_merging:
        _emit_generation_progress(
            progress_callback,
            "current merge",
            f"{numerical_current_samples} samples",
        )
        phase_started = time.perf_counter()
        dag, current_merging = _merge_identical_currents_by_warmup(
            dag,
            model,
            sample_count=numerical_current_samples,
            seed=numerical_current_seed,
            relative_tolerance=numerical_current_relative_tolerance,
            zero_tolerance=numerical_current_zero_tolerance,
        )
        _emit_generation_progress(
            progress_callback,
            "merge ready",
            _dag_progress_label(dag),
            increment=1,
            duration_s=time.perf_counter() - phase_started,
        )
    else:
        current_merging = {
            "enabled": False,
            "reason": "disabled by --no-numerical-current-merging",
        }
        _emit_generation_progress(
            progress_callback,
            "current merge",
            "disabled",
            increment=1,
            duration_s=0.0,
        )
    _emit_generation_progress(
        progress_callback,
        "color plan",
        dag.process.color_accuracy,
    )
    phase_started = time.perf_counter()
    full_color_plan = build_color_plan(
        dag.process,
        color_accuracy=dag.process.color_accuracy,
        options=options,
        max_sectors=max_color_sectors,
        reference_color_order=reference_color_order,
    )
    _emit_generation_progress(
        progress_callback,
        "color ready",
        f"sectors={full_color_plan.sector_count}",
        increment=1,
        duration_s=time.perf_counter() - phase_started,
    )
    if (
        (selected_color_sector_ids is not None or reference_color_order is not None)
        and (
            full_color_plan.truncated
            or full_color_plan.sector_count
            > _FULL_COLOR_PLAN_SERIALIZATION_SECTOR_LIMIT
        )
    ):
        color_plan = dag.color_plan
    else:
        color_plan = full_color_plan
    return GenericProcessManifest(
        dag=dag,
        model=model,
        color_plan=color_plan,
        structural_current_aggregation=structural_current_aggregation,
        zero_current_filter=zero_current_filter,
        current_merging=current_merging,
    )


def build_generic_process_set_manifest(
    process_string: str,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    color_accuracy: str = "lc",
    max_currents: int | None = None,
    max_color_sectors: int | None = None,
    selected_color_sector_ids: set[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    selected_source_helicities: Mapping[int, int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
) -> GenericProcessSetManifest:
    model = model or AmplicolSMLeadingColorModel()
    process_set = build_process_set_ir(
        process_string,
        color_accuracy=color_accuracy,
        options=options,
        generic=True,
        max_quark_pairs=max_quark_pairs,
    )
    manifests = tuple(
        build_generic_process_manifest(
            entry.ir,
            model=model,
            options=options,
            max_currents=max_currents,
            max_color_sectors=max_color_sectors,
            selected_color_sector_ids=selected_color_sector_ids,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            selected_source_helicities=selected_source_helicities,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
        )
        for entry in process_set.entries
    )
    generation_metadata = _generic_process_set_generation_metadata(
        options=options,
        selected_color_sector_ids=selected_color_sector_ids,
        lc_topology_replay=False,
        max_coupling_orders=max_coupling_orders,
        max_lc_current_line_groups=max_lc_current_line_groups,
        max_quark_pairs=max_quark_pairs,
        closure_side_mask_pruning=closure_side_mask_pruning,
        color_order_mask_pruning=color_order_mask_pruning,
        species_reachability_pruning=species_reachability_pruning,
        ignored_particle_ids=ignored_particle_ids,
        ignored_vertex_kinds=ignored_vertex_kinds,
        selected_source_helicities=selected_source_helicities,
        numerical_filter_current=numerical_filter_current,
        numerical_current_merging=numerical_current_merging,
        numerical_current_samples=numerical_current_samples,
        numerical_current_seed=numerical_current_seed,
        numerical_current_relative_tolerance=numerical_current_relative_tolerance,
        numerical_current_zero_tolerance=numerical_current_zero_tolerance,
    )
    return GenericProcessSetManifest(
        process_set=process_set,
        processes=manifests,
        generation_metadata=generation_metadata,
    )


def _lc_topology_reuse_payload(
    manifest: GenericProcessManifest,
) -> dict[str, object]:
    if manifest.color_plan.color_accuracy != "lc":
        return {
            "available": False,
            "reason": "topology reuse is currently defined for leading-colour sectors",
        }
    groups = manifest.color_plan.topology_groups
    if not groups:
        return {
            "available": False,
            "reason": "no leading-colour topology groups were built",
        }
    dag_sector_ids = tuple(
        sorted(
            {
                current.index.color_state.sector_id
                for current in manifest.dag.currents
            }
        )
    )
    root_sector_ids = tuple(
        sorted(
            {
                _root_color_sector_id(manifest.dag, root)
                for root in manifest.dag.amplitude_roots
            }
        )
    )
    active_sector_ids = set(root_sector_ids or dag_sector_ids)
    line_pairing_representative_ids = lc_line_pairing_representative_ids(
        manifest.color_plan
    )
    safe_group_ids = {
        id(group)
        for group in lc_topology_replay_safe_groups(manifest.color_plan)
    }
    active_groups = []
    representative_ids: list[int] = []
    safe_representative_ids: list[int] = []
    for group in groups:
        active_members = tuple(
            sector_id
            for sector_id in group.sector_ids
            if sector_id in active_sector_ids
        )
        if not active_members:
            continue
        representative = (
            group.representative_sector_id
            if group.representative_sector_id in active_members
            else active_members[0]
        )
        representative_ids.append(representative)
        replay_safe = id(group) in safe_group_ids
        if replay_safe:
            safe_representative_ids.append(representative)
        active_groups.append(
            {
                **group.to_json_dict(),
                "active_sector_ids": list(active_members),
                "runtime_reuse_factor": len(active_members),
                "representative_sector_id": representative,
                "runtime_replay_safe": replay_safe,
            }
        )
    return {
        "available": bool(active_groups),
        "policy": (
            "compile LC topology representatives and evaluate sibling sectors "
            "by external-label/source-slot permutations"
        ),
        "full_color_sector_count": manifest.color_plan.sector_count,
        "dag_color_sector_ids": list(dag_sector_ids),
        "amplitude_color_sector_ids": list(root_sector_ids),
        "topology_group_count": len(groups),
        "active_topology_group_count": len(active_groups),
        "representative_sector_ids": representative_ids,
        "line_pairing_representative_sector_ids": list(
            line_pairing_representative_ids
        ),
        "line_pairing_representative_sector_count": len(
            line_pairing_representative_ids
        ),
        "replay_safe_topology_group_count": sum(
            1 for group in active_groups if group["runtime_replay_safe"]
        ),
        "replay_safe_representative_sector_ids": safe_representative_ids,
        "groups": active_groups,
    }


def _lc_topology_replay_payload(
    manifest: GenericProcessManifest,
    *,
    materialized_sector_ids: set[int] | None,
    enabled: bool,
) -> dict[str, object]:
    if not enabled:
        return {
            "enabled": False,
            "reason": "lc topology replay was not requested",
        }
    if manifest.color_plan.color_accuracy != "lc":
        return {
            "enabled": False,
            "reason": "lc topology replay is only defined for leading colour",
        }
    if materialized_sector_ids is None:
        return {
            "enabled": False,
            "reason": "all contributing colour sectors are already materialized",
        }
    groups = []
    replayed_sector_count = 0
    materialized = set(int(sector_id) for sector_id in materialized_sector_ids)
    safe_groups = lc_topology_replay_safe_groups(manifest.color_plan)
    if not safe_groups:
        return {
            "enabled": False,
            "reason": "no replay-safe topology groups are available",
            "materialized_sector_ids": sorted(materialized),
        }
    for group in safe_groups:
        if group.representative_sector_id not in materialized:
            continue
        sector_permutations = []
        for sector_id, permutation in zip(
            group.sector_ids,
            group.label_permutations,
            strict=True,
        ):
            sector_permutations.append(
                {
                    "sector_id": int(sector_id),
                    "label_permutation": [
                        {
                            "representative_label": int(representative_label),
                            "sector_label": int(sector_label),
                        }
                        for representative_label, sector_label in permutation
                    ],
                }
            )
        groups.append(
            {
                **group.to_json_dict(),
                "materialized_sector_id": int(group.representative_sector_id),
                "active_sector_ids": [int(sector_id) for sector_id in group.sector_ids],
                "sector_permutations": sector_permutations,
            }
        )
        replayed_sector_count += len(group.sector_ids)
    if not groups:
        return {
            "enabled": False,
            "reason": "no topology group representative is materialized",
            "materialized_sector_ids": sorted(materialized),
        }
    return {
        "enabled": True,
        "mode": "external-label-permutation",
        "normalization": "sum replayed sector raw sums before global normalization",
        "materialized_sector_ids": sorted(materialized),
        "replayed_sector_count": replayed_sector_count,
        "group_count": len(groups),
        "groups": groups,
    }


def _lc_topology_replay_payload_from_partition(
    partition: LCColorSectorReplayPartition,
) -> dict[str, object]:
    replay_weights = (
        partition.replay_weights
        if partition.replay_weights
        else tuple(1.0 for _ in partition.active_sector_ids)
    )
    sector_permutations = [
        {
            "sector_id": int(sector_id),
            "weight": float(weight),
            "label_permutation": [
                {
                    "representative_label": int(representative_label),
                    "sector_label": int(sector_label),
                }
                for representative_label, sector_label in permutation
            ],
        }
        for sector_id, permutation, weight in zip(
            partition.active_sector_ids,
            partition.label_permutations,
            replay_weights,
            strict=True,
        )
    ]
    group = {
        "representative_sector_id": int(partition.representative_sector_id),
        "materialized_sector_id": int(partition.representative_sector_id),
        "sector_ids": [int(sector_id) for sector_id in partition.active_sector_ids],
        "active_sector_ids": [
            int(sector_id) for sector_id in partition.active_sector_ids
        ],
        "label_permutations": [
            [[int(left), int(right)] for left, right in permutation]
            for permutation in partition.label_permutations
        ],
        "replay_weights": [float(weight) for weight in replay_weights],
        "sector_permutations": sector_permutations,
    }
    return {
        "enabled": True,
        "mode": "external-label-permutation",
        "normalization": "sum replayed sector raw sums before global normalization",
        "materialized_sector_ids": [int(partition.representative_sector_id)],
        "replayed_sector_count": len(partition.active_sector_ids),
        "group_count": 1,
        "groups": [group],
    }


def write_generic_process_manifest(
    manifest: GenericProcessManifest,
    output_dir: str | Path,
    *,
    filename: str = "generic_process_manifest.json",
) -> Path:
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    manifest_path = output_path / filename
    manifest_path.write_text(
        json.dumps(
            _json_safe_bigints(manifest.to_json_dict()),
            separators=(",", ":"),
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return manifest_path


def write_generic_process_set_manifest(
    manifest: GenericProcessSetManifest,
    output_dir: str | Path,
    *,
    filename: str = "generic_process_set_manifest.json",
) -> Path:
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    for subprocess_manifest in manifest.processes:
        write_generic_process_manifest(
            subprocess_manifest,
            output_path / "subprocesses" / subprocess_manifest.key,
        )
    manifest_path = output_path / filename
    manifest_path.write_text(
        json.dumps(_json_safe_bigints(manifest.to_json_dict()), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return manifest_path


def _emit_generation_progress(
    progress_callback: Any | None,
    stage: str,
    item: object = "",
    *,
    increment: int = 0,
    total: int | None = None,
    duration_s: float | None = None,
) -> None:
    if progress_callback is None:
        return
    event: dict[str, object] = {
        "stage": stage,
        "item": str(item),
    }
    if increment:
        event["increment"] = int(increment)
    if total is not None:
        event["total"] = int(total)
    if duration_s is not None:
        event["duration_s"] = float(duration_s)
    progress_callback(event)


def _process_progress_label(
    process: str | CanonicalProcessIR | GenericDAG,
) -> str:
    if isinstance(process, GenericDAG):
        return process.process.process
    if isinstance(process, CanonicalProcessIR):
        return process.process
    return str(process)


def _dag_progress_label(dag: GenericDAG) -> str:
    return (
        f"c={len(dag.currents)} i={len(dag.interactions)} "
        f"r={len(dag.amplitude_roots)}"
    )


def write_generic_dag_process_artifact(
    process: str | CanonicalProcessIR | GenericDAG | GenericProcessManifest,
    output_dir: str | Path,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    color_accuracy: str = "lc",
    max_currents: int | None = None,
    max_color_sectors: int | None = None,
    evaluator_backend: str = "compiled-complex",
    compiled_preset: str = "runtime-o3",
    batch_size: int = 64,
    emit_stage_evaluator_artifacts: bool = False,
    skip_main_stage_evaluator_artifacts: bool = False,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    stage_local_parameter_layout: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    write_generic_plan: bool = True,
    enable_lc_sector_runtime_selector: bool | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    runtime_lc_sector_ids: set[int] | None = None,
    lc_topology_replay: bool = False,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    selected_source_helicities: Mapping[int, int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
) -> tuple[Path, dict[str, object]]:
    """Write a schema-v2 generic DAG process artifact.

    This is the production-facing artifact for the generic refactor. It stores
    the full model-driven DAG, requested evaluator settings, symbolic
    stage-blueprint lowering, and optionally serialized generic stage
    evaluators loadable by Rusticol schema-v2 execution.
    """

    if enable_lc_sector_runtime_selector is None:
        enable_lc_sector_runtime_selector = selected_color_sector_ids is None
    generic_manifest = (
        process
        if isinstance(process, GenericProcessManifest)
        else build_generic_process_manifest(
            process,
            model=model,
            options=options,
            color_accuracy=color_accuracy,
            max_currents=max_currents,
            max_color_sectors=max_color_sectors,
            reference_color_order=reference_color_order,
            selected_color_sector_ids=selected_color_sector_ids,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            selected_source_helicities=selected_source_helicities,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
            progress_callback=progress_callback,
        )
    )
    if isinstance(process, GenericProcessManifest):
        _emit_generation_progress(
            progress_callback,
            "dag ready",
            _dag_progress_label(generic_manifest.dag),
            increment=1,
        )
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    if skip_main_stage_evaluator_artifacts:
        _emit_generation_progress(
            progress_callback,
            "main artifact",
            "compact manifest",
        )
        payload = _generic_dag_process_compact_artifact_payload(
            generic_manifest,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            stage_local_parameter_layout=stage_local_parameter_layout,
            reference_color_order=reference_color_order,
            selected_color_sector_ids=selected_color_sector_ids,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
        )
        _emit_generation_progress(
            progress_callback,
            "main artifact",
            "compact ready",
            increment=1,
        )
    else:
        _emit_generation_progress(
            progress_callback,
            "main plan",
            "serialize" if write_generic_plan else "skipped",
        )
        if write_generic_plan:
            phase_started = time.perf_counter()
            plan_path = write_generic_process_manifest(generic_manifest, output_path)
            plan_path_name: str | None = plan_path.name
            _emit_generation_progress(
                progress_callback,
                "main plan",
                plan_path.name,
                increment=1,
                duration_s=time.perf_counter() - phase_started,
            )
        else:
            plan_path_name = None
            _emit_generation_progress(
                progress_callback,
                "main plan",
                "skipped",
                increment=1,
                duration_s=0.0,
            )
        phase_started = time.perf_counter()
        payload = _generic_dag_process_artifact_payload(
            generic_manifest,
            plan_path=plan_path_name,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            stage_evaluator_artifact_dir=(
                output_path if emit_stage_evaluator_artifacts else None
            ),
            stage_evaluator_compiler=stage_evaluator_compiler,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            stage_local_parameter_layout=stage_local_parameter_layout,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
            reference_color_order=reference_color_order,
            selected_color_sector_ids=selected_color_sector_ids,
            enable_lc_sector_runtime_selector=bool(enable_lc_sector_runtime_selector),
            lc_topology_replay=lc_topology_replay,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
        )
        _emit_generation_progress(
            progress_callback,
            "main artifact",
            "runtime ready" if payload.get("runtime_schema") else "metadata ready",
            increment=1,
            duration_s=time.perf_counter() - phase_started,
        )
    sidecar_metadata: list[dict[str, object]] = []
    if runtime_lc_sector_ids and generic_manifest.dag.process.color_accuracy == "lc":
        _emit_generation_progress(
            progress_callback,
            "runtime sidecar",
            _runtime_lc_sector_artifact_key(runtime_lc_sector_ids),
        )
        sidecar_dir = (
            output_path
            / "runtime_lc_sectors"
            / _runtime_lc_sector_artifact_key(runtime_lc_sector_ids)
        )
        sidecar_dir.mkdir(parents=True, exist_ok=True)
        phase_started = time.perf_counter()
        sidecar_manifest = _runtime_manifest_for_color_sectors(
            generic_manifest,
            runtime_lc_sector_ids,
        )
        _emit_generation_progress(
            progress_callback,
            "sidecar dag",
            _dag_progress_label(sidecar_manifest.dag),
            increment=1,
            duration_s=time.perf_counter() - phase_started,
        )
        if write_generic_plan:
            phase_started = time.perf_counter()
            sidecar_plan_path = write_generic_process_manifest(
                sidecar_manifest,
                sidecar_dir,
            )
            sidecar_plan_path_name: str | None = sidecar_plan_path.name
            _emit_generation_progress(
                progress_callback,
                "sidecar plan",
                sidecar_plan_path.name,
                increment=1,
                duration_s=time.perf_counter() - phase_started,
            )
        else:
            sidecar_plan_path_name = None
            _emit_generation_progress(
                progress_callback,
                "sidecar plan",
                "skipped",
                increment=1,
                duration_s=0.0,
            )
        phase_started = time.perf_counter()
        sidecar_payload = _generic_dag_process_artifact_payload(
            sidecar_manifest,
            plan_path=sidecar_plan_path_name,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            stage_evaluator_artifact_dir=(
                sidecar_dir if emit_stage_evaluator_artifacts else None
            ),
            stage_evaluator_compiler=stage_evaluator_compiler,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            stage_local_parameter_layout=stage_local_parameter_layout,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
            reference_color_order=reference_color_order,
            selected_color_sector_ids=None,
            enable_lc_sector_runtime_selector=False,
            lc_topology_replay=lc_topology_replay,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
        )
        sidecar_manifest_path = sidecar_dir / "process_manifest.json"
        sidecar_manifest_path.write_text(
            json.dumps(sidecar_payload, separators=(",", ":"), sort_keys=True),
            encoding="utf-8",
        )
        _emit_generation_progress(
            progress_callback,
            "sidecar artifact",
            "process_manifest.json",
            increment=1,
            duration_s=time.perf_counter() - phase_started,
        )
        (sidecar_dir / "check_standalone.py").write_text(
            _GENERIC_DAG_CHECK_STANDALONE,
            encoding="utf-8",
        )
        _write_generic_validation_momenta(generic_manifest, sidecar_dir)
        sidecar_compiled = cast(dict[str, Any], sidecar_payload["compiled"])
        sidecar_metadata.append(
            {
                "color_sector_ids": sorted(
                    int(sector_id) for sector_id in runtime_lc_sector_ids
                ),
                "path": str(sidecar_dir.relative_to(output_path)),
                "kind": "selected-lc-runtime-sidecar",
                "runtime_available": bool(sidecar_compiled.get("runtime_available")),
                "runtime_selector": "none-specialized-runtime-artifact",
            }
        )
    if sidecar_metadata:
        compiled_payload = cast(dict[str, Any], payload["compiled"])
        compiled_payload["runtime_lc_sector_artifacts"] = sidecar_metadata
    manifest_path = output_path / "process_manifest.json"
    phase_started = time.perf_counter()
    manifest_path.write_text(
        json.dumps(payload, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )
    _emit_generation_progress(
        progress_callback,
        "manifest",
        "process_manifest.json",
        increment=1,
        duration_s=time.perf_counter() - phase_started,
    )
    (output_path / "check_standalone.py").write_text(
        _GENERIC_DAG_CHECK_STANDALONE,
        encoding="utf-8",
    )
    _write_generic_validation_momenta(generic_manifest, output_path)
    return manifest_path, payload


def write_lc_topology_replay_partition_artifact(
    process: str | CanonicalProcessIR,
    output_dir: str | Path,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    color_accuracy: str = "lc",
    max_currents: int | None = None,
    max_color_sectors: int | None = None,
    evaluator_backend: str = "compiled-complex",
    compiled_preset: str = "runtime-o3",
    batch_size: int = 64,
    emit_stage_evaluator_artifacts: bool = False,
    skip_main_stage_evaluator_artifacts: bool = False,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    stage_local_parameter_layout: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    write_generic_plan: bool = True,
    reference_color_order: Sequence[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    selected_source_helicities: Mapping[int, int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
) -> tuple[Path, dict[str, object]]:
    """Write an LC all-ordering artifact as exact replay partitions.

    Each representative sidecar is a normal selected-sector schema-v2 runtime.
    Its manifest is then annotated with the subset of LC sectors that can be
    replayed by external-label permutation while preserving the physical
    initial-state labels.  This avoids materializing a monolithic all-sector DAG
    and keeps every generated process output available for later re-timing.
    """

    if color_accuracy != "lc":
        raise ValueError("LC replay partition artifacts require color_accuracy='lc'")
    model = model or AmplicolSMLeadingColorModel()
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    _emit_generation_progress(
        progress_callback,
        "lc replay plan",
        _process_progress_label(process),
    )
    phase_started = time.perf_counter()
    color_plan = build_color_plan(
        process,
        color_accuracy="lc",
        options=options,
        max_sectors=max_color_sectors,
        reference_color_order=reference_color_order,
    )
    partitions = lc_topology_replay_partitions(color_plan)
    _emit_generation_progress(
        progress_callback,
        "lc replay plan",
        f"partitions={len(partitions)} sectors={color_plan.sector_count}",
        increment=1,
        duration_s=time.perf_counter() - phase_started,
    )
    if color_plan.truncated:
        raise ValueError(
            "LC replay partition generation requires the complete colour plan; "
            f"the current plan was truncated at max_color_sectors={max_color_sectors}"
        )
    if not partitions:
        raise ValueError(
            "no LC replay partitions are available for this process; use normal "
            "generic DAG generation or a selected LC sector instead"
        )

    representative_root = output_path / "lc_topology_replay_representatives"
    representative_root.mkdir(parents=True, exist_ok=True)
    generation_metadata = _generic_process_set_generation_metadata(
        options=options,
        selected_color_sector_ids=None,
        lc_topology_replay=True,
        max_coupling_orders=max_coupling_orders,
        max_lc_current_line_groups=max_lc_current_line_groups,
        max_quark_pairs=max_quark_pairs,
        closure_side_mask_pruning=closure_side_mask_pruning,
        color_order_mask_pruning=color_order_mask_pruning,
        species_reachability_pruning=species_reachability_pruning,
        ignored_particle_ids=ignored_particle_ids,
        ignored_vertex_kinds=ignored_vertex_kinds,
        selected_source_helicities=selected_source_helicities,
        numerical_filter_current=numerical_filter_current,
        numerical_current_merging=numerical_current_merging,
        numerical_current_samples=numerical_current_samples,
        numerical_current_seed=numerical_current_seed,
        numerical_current_relative_tolerance=numerical_current_relative_tolerance,
        numerical_current_zero_tolerance=numerical_current_zero_tolerance,
    )
    representative_artifacts: list[dict[str, object]] = []
    total_counts = Counter()
    replayed_root_count = 0
    runtime_available = True
    runtime_unavailable_messages: list[str] = []
    first_sidecar_dir: Path | None = None
    for index, partition in enumerate(partitions, start=1):
        representative_id = int(partition.representative_sector_id)
        _emit_generation_progress(
            progress_callback,
            "lc replay representative",
            f"{index}/{len(partitions)} sector={representative_id}",
        )
        sidecar_dir = representative_root / f"sector_{representative_id}"
        sidecar_manifest_path, sidecar_payload = write_generic_dag_process_artifact(
            process,
            sidecar_dir,
            model=model,
            options=options,
            color_accuracy="lc",
            max_currents=max_currents,
            max_color_sectors=max_color_sectors,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            emit_stage_evaluator_artifacts=emit_stage_evaluator_artifacts,
            skip_main_stage_evaluator_artifacts=skip_main_stage_evaluator_artifacts,
            stage_evaluator_compiler=stage_evaluator_compiler,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            stage_local_parameter_layout=stage_local_parameter_layout,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
            write_generic_plan=write_generic_plan,
            enable_lc_sector_runtime_selector=False,
            reference_color_order=reference_color_order,
            selected_color_sector_ids={representative_id},
            runtime_lc_sector_ids=None,
            lc_topology_replay=False,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            selected_source_helicities=selected_source_helicities,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
        )
        replay_payload = _lc_topology_replay_payload_from_partition(partition)
        compiled_payload = cast(dict[str, Any], sidecar_payload["compiled"])
        compiled_payload["lc_topology_replay"] = replay_payload
        sidecar_payload["runtime_lc_topology_replay"] = replay_payload
        sidecar_manifest_path.write_text(
            json.dumps(sidecar_payload, separators=(",", ":"), sort_keys=True),
            encoding="utf-8",
        )
        if first_sidecar_dir is None:
            first_sidecar_dir = sidecar_dir
        sidecar_runtime_available = bool(compiled_payload.get("runtime_available"))
        if not sidecar_runtime_available:
            runtime_available = False
            message = compiled_payload.get("runtime_unavailable_message")
            if isinstance(message, str) and message:
                runtime_unavailable_messages.append(message)
        dag_summary = sidecar_payload.get("dag_summary")
        if isinstance(dag_summary, Mapping):
            for key in (
                "current_count",
                "source_count",
                "interaction_count",
                "amplitude_root_count",
            ):
                value = dag_summary.get(key)
                if isinstance(value, int):
                    total_counts[key] += value
            root_count = dag_summary.get("amplitude_root_count")
            if isinstance(root_count, int):
                replayed_root_count += root_count * len(partition.active_sector_ids)
        representative_artifacts.append(
            {
                "representative_sector_id": representative_id,
                "active_sector_ids": [
                    int(sector_id) for sector_id in partition.active_sector_ids
                ],
                "replayed_sector_count": len(partition.active_sector_ids),
                "path": str(sidecar_dir.relative_to(output_path)),
                "manifest": str(sidecar_manifest_path.relative_to(output_path)),
                "kind": GENERIC_DAG_PROCESS_ARTIFACT_KIND,
                "artifact_class": "generic-dag-schema-v2",
                "runtime_available": sidecar_runtime_available,
                "dag_summary": dag_summary if isinstance(dag_summary, Mapping) else {},
                "lc_topology_replay": replay_payload,
            }
        )
        _emit_generation_progress(
            progress_callback,
            "lc replay representative",
            f"{index}/{len(partitions)} ready",
            increment=1,
        )
    if first_sidecar_dir is not None:
        validation_source = first_sidecar_dir / "validation_momenta.json"
        if validation_source.exists():
            shutil.copy2(validation_source, output_path / "validation_momenta.json")

    process_ir = color_plan.process
    payload: dict[str, object] = {
        "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
        "kind": GENERIC_LC_REPLAY_PARTITION_ARTIFACT_KIND,
        "artifact_class": "lc-replay-partition-schema-v2",
        "process": process_ir.process,
        "key": process_ir.key,
        "color_accuracy": "lc",
        "model": {
            "name": model.name,
            "vertex_lowering_coverage": (
                model.vertex_lowering_coverage().to_json_dict()
            ),
        },
        "external_pdg_order": [
            *process_ir.initial_pdgs,
            *process_ir.final_pdgs,
        ],
        "outgoing_pdg_order": list(process_ir.outgoing_pdgs),
        "process_ir": process_ir.to_json_dict(),
        "color_plan": color_plan.to_json_dict(),
        "generic_generation": generation_metadata,
        "planning_status": {
            "color_ready": color_plan.ready_for_requested_colour,
            "color_sector_count": color_plan.sector_count,
            "color_truncated": color_plan.truncated,
            "generic_evaluator_ready": runtime_available,
            "lc_replay_partition_count": len(partitions),
        },
        "lowering_status": {
            "current_color_sector_count": len(partitions),
            "amplitude_color_sector_count": color_plan.sector_count,
            "replay_partition_count": len(partitions),
            "replayed_color_sector_count": sum(
                len(partition.active_sector_ids) for partition in partitions
            ),
            "full_tensor_network_ready": runtime_available,
        },
        "compiled": {
            "kind": "lc-topology-replay-partition",
            "runtime_available": runtime_available,
            "runtime_unavailable_message": (
                None
                if runtime_available
                else "; ".join(sorted(set(runtime_unavailable_messages)))
                or "one or more LC replay representative artifacts are unavailable"
            ),
            "requested_evaluator_backend": evaluator_backend,
            "requested_compiled_preset": compiled_preset,
            "batch_size": batch_size,
            "stage_local_parameter_layout": stage_local_parameter_layout,
            "representative_count": len(representative_artifacts),
            "replayed_sector_count": sum(
                len(partition.active_sector_ids) for partition in partitions
            ),
            "representative_artifacts": representative_artifacts,
        },
        "dag_summary": {
            "current_count": int(total_counts["current_count"]),
            "source_count": int(total_counts["source_count"]),
            "interaction_count": int(total_counts["interaction_count"]),
            "amplitude_root_count": int(total_counts["amplitude_root_count"]),
            "materialized_color_sector_count": len(partitions),
            "replayed_color_sector_count": sum(
                len(partition.active_sector_ids) for partition in partitions
            ),
            "truncated": False,
        },
        "full_dag_summary": {
            "current_count": int(total_counts["current_count"]),
            "source_count": int(total_counts["source_count"]),
            "interaction_count": int(total_counts["interaction_count"]),
            "amplitude_root_count": replayed_root_count,
            "materialized_amplitude_root_count": int(
                total_counts["amplitude_root_count"]
            ),
            "materialized_color_sector_count": len(partitions),
            "replayed_color_sector_count": sum(
                len(partition.active_sector_ids) for partition in partitions
            ),
            "truncated": False,
        },
    }
    manifest_path = output_path / "process_manifest.json"
    manifest_path.write_text(
        json.dumps(payload, separators=(",", ":"), sort_keys=True),
        encoding="utf-8",
    )
    (output_path / "check_standalone.py").write_text(
        _GENERIC_DAG_CHECK_STANDALONE,
        encoding="utf-8",
    )
    return manifest_path, payload


def select_leading_color_sector_ids(
    manifest: GenericProcessManifest,
    *,
    reference_color_order: Sequence[int] | None = None,
) -> set[int] | None:
    """Select the LC sector matching a legacy AmpliCol colour ordering.

    For LC processes with multiple colour sectors, Fortran AmpliCol probe output
    reports one group/integral at a time. This helper selects the corresponding
    generic DAG sector so point-by-point validation compares the same ordered
    quantity. If no ordering is supplied, it falls back to the first
    contributing sector.
    """

    if (
        manifest.color_plan.color_accuracy != "lc"
        or manifest.color_plan.sector_count <= 1
    ):
        return None
    root_sectors = sorted(
        {
            _root_color_sector_id(manifest.dag, root)
            for root in manifest.dag.amplitude_roots
        }
    )
    if not root_sectors:
        return None
    if reference_color_order is not None:
        wanted = tuple(int(label) for label in reference_color_order)
        wanted_coloured = _reference_coloured_word(manifest.color_plan, wanted)
        if wanted_coloured:
            for sector in manifest.color_plan.sectors:
                if sector.id not in root_sectors:
                    continue
                if tuple(getattr(sector, "word_labels", ())) == wanted_coloured:
                    if wanted_coloured == wanted:
                        return {int(sector.id)}
                    return _lc_colored_word_sibling_sector_ids(
                        manifest.color_plan,
                        sector,
                        active_sector_ids=set(root_sectors),
                        reference_order=wanted,
                    )
        for sector in manifest.color_plan.sectors:
            if sector.id not in root_sectors:
                continue
            if wanted in sector.color_words:
                if wanted_coloured == wanted:
                    return {int(sector.id)}
                return _lc_colored_word_sibling_sector_ids(
                    manifest.color_plan,
                    sector,
                    active_sector_ids=set(root_sectors),
                    reference_order=wanted,
                )
            legacy_order_words = getattr(sector, "legacy_order_words", ())
            if wanted in legacy_order_words:
                if wanted_coloured == wanted:
                    return {int(sector.id)}
                return _lc_colored_word_sibling_sector_ids(
                    manifest.color_plan,
                    sector,
                    active_sector_ids=set(root_sectors),
                    reference_order=wanted,
                )
        for sector in manifest.color_plan.sectors:
            if sector.id not in root_sectors:
                continue
            if wanted in sector.compatibility_words:
                return _lc_colored_word_sibling_sector_ids(
                    manifest.color_plan,
                    sector,
                    active_sector_ids=set(root_sectors),
                )
        raise ValueError(
            "LC reference colour order does not match any generated colour sector: "
            f"{wanted}"
        )
    return {0 if 0 in root_sectors else root_sectors[0]}


def select_leading_color_sector_ids_from_plan(
    color_plan: GenericColorPlan,
    *,
    reference_color_order: Sequence[int] | None = None,
) -> set[int] | None:
    """Select an LC sector from colour-plan metadata before DAG construction."""

    if color_plan.color_accuracy != "lc" or color_plan.sector_count <= 1:
        return None
    if reference_color_order is not None:
        wanted = tuple(int(label) for label in reference_color_order)
        wanted_coloured = _reference_coloured_word(color_plan, wanted)
        if wanted_coloured:
            for sector in color_plan.sectors:
                if tuple(getattr(sector, "word_labels", ())) == wanted_coloured:
                    if wanted_coloured == wanted:
                        return {int(sector.id)}
                    return _lc_colored_word_sibling_sector_ids(
                        color_plan,
                        sector,
                        reference_order=wanted,
                    )
        for sector in color_plan.sectors:
            if wanted in sector.color_words:
                if wanted_coloured == wanted:
                    return {int(sector.id)}
                return _lc_colored_word_sibling_sector_ids(
                    color_plan,
                    sector,
                    reference_order=wanted,
                )
            legacy_order_words = getattr(sector, "legacy_order_words", ())
            if wanted in legacy_order_words:
                if wanted_coloured == wanted:
                    return {int(sector.id)}
                return _lc_colored_word_sibling_sector_ids(
                    color_plan,
                    sector,
                    reference_order=wanted,
                )
        for sector in color_plan.sectors:
            if wanted in sector.compatibility_words:
                return _lc_colored_word_sibling_sector_ids(color_plan, sector)
        raise ValueError(
            "LC reference colour order does not match any generated colour sector: "
            f"{wanted}"
        )
    return {0}


def _reference_coloured_word(
    color_plan,
    reference_order: Sequence[int],
) -> tuple[int, ...]:
    coloured = set(getattr(color_plan, "coloured_labels", ()) or ())
    return tuple(int(label) for label in reference_order if int(label) in coloured)


def _lc_colored_word_sibling_sector_ids(
    color_plan,
    sector,
    *,
    active_sector_ids: set[int] | None = None,
    reference_order: Sequence[int] | None = None,
) -> set[int]:
    """Return sectors with the same LC coloured word as ``sector``.

    Colour-singlet attachments can be assigned to different open fermion lines
    during current construction, but those assignments are diagrams of the same
    colour-ordered amplitude and must be generated together for coherent
    squaring.  The coloured word deliberately excludes those singlet labels.
    """

    word = tuple(getattr(sector, "word_labels", ()))
    if not word:
        return {int(sector.id)}
    siblings = set()
    for candidate in color_plan.sectors:
        if tuple(getattr(candidate, "word_labels", ())) != word:
            continue
        if active_sector_ids is not None and int(candidate.id) not in active_sector_ids:
            continue
        if reference_order is not None and not _lc_sector_preserves_reference_singlet_blocks(
            candidate,
            reference_order,
        ):
            continue
        siblings.add(int(candidate.id))
    return siblings or {int(sector.id)}


def _lc_sector_preserves_reference_singlet_blocks(
    sector,
    reference_order: Sequence[int],
) -> bool:
    """Return whether sector singlet allocations keep reference blocks intact.

    Fortran AmpliCol colour orders include colour-singlet particles in the
    ordered row.  The same ordered amplitude can emit an entire contiguous
    singlet block from either neighbouring open colour line, but one block is
    not split over several open lines.  This keeps sibling-sector replay
    generic while avoiding overcounting split allocations such as ``Z`` on one
    line and ``H`` on another for a reference order ending in ``Z H``.
    """

    if getattr(sector, "kind", "") != "open-lines":
        return True
    colored_labels = set(getattr(sector, "word_labels", ()) or ())
    if not colored_labels:
        return True
    blocks: list[set[int]] = []
    current: list[int] = []
    for raw_label in reference_order:
        label = int(raw_label)
        if label in colored_labels:
            if current:
                blocks.append(set(current))
                current = []
            continue
        current.append(label)
    if current:
        blocks.append(set(current))
    if not blocks:
        return True
    line_singlets = [set(line.singlet_labels) for line in getattr(sector, "quark_lines", ())]
    for block in blocks:
        if not block:
            continue
        if not any(block.issubset(singlets) for singlets in line_singlets):
            return False
    return True


def write_generic_dag_process_set_artifact(
    process_set: ProcessSetIR,
    output_dir: str | Path,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    max_currents: int | None = None,
    max_color_sectors: int | None = None,
    evaluator_backend: str = "compiled-complex",
    compiled_preset: str = "runtime-o3",
    batch_size: int = 64,
    emit_stage_evaluator_artifacts: bool = False,
    skip_main_stage_evaluator_artifacts: bool = False,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    stage_local_parameter_layout: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    enable_lc_sector_runtime_selector: bool | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    runtime_lc_sector_ids: set[int] | None = None,
    lc_topology_replay: bool = False,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    selected_source_helicities: Mapping[int, int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
) -> tuple[Path, dict[str, object]]:
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    generation_metadata = _generic_process_set_generation_metadata(
        options=options,
        selected_color_sector_ids=selected_color_sector_ids,
        lc_topology_replay=lc_topology_replay,
        max_coupling_orders=max_coupling_orders,
        max_lc_current_line_groups=max_lc_current_line_groups,
        max_quark_pairs=max_quark_pairs,
        closure_side_mask_pruning=closure_side_mask_pruning,
        color_order_mask_pruning=color_order_mask_pruning,
        species_reachability_pruning=species_reachability_pruning,
        ignored_particle_ids=ignored_particle_ids,
        ignored_vertex_kinds=ignored_vertex_kinds,
        selected_source_helicities=selected_source_helicities,
        numerical_filter_current=numerical_filter_current,
        numerical_current_merging=numerical_current_merging,
        numerical_current_samples=numerical_current_samples,
        numerical_current_seed=numerical_current_seed,
        numerical_current_relative_tolerance=numerical_current_relative_tolerance,
        numerical_current_zero_tolerance=numerical_current_zero_tolerance,
    )
    processes = []
    for entry in process_set.entries:
        subdir = output_path / "subprocesses" / entry.key
        manifest_path, payload = write_generic_dag_process_artifact(
            entry.ir,
            subdir,
            model=model,
            options=options,
            color_accuracy=process_set.color_accuracy,
            max_currents=max_currents,
            max_color_sectors=max_color_sectors,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            emit_stage_evaluator_artifacts=emit_stage_evaluator_artifacts,
            skip_main_stage_evaluator_artifacts=skip_main_stage_evaluator_artifacts,
            stage_evaluator_compiler=stage_evaluator_compiler,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            stage_local_parameter_layout=stage_local_parameter_layout,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
            enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
            selected_color_sector_ids=selected_color_sector_ids,
            runtime_lc_sector_ids=runtime_lc_sector_ids,
            lc_topology_replay=lc_topology_replay,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            selected_source_helicities=selected_source_helicities,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
        )
        compiled_payload = cast(dict[str, Any], payload["compiled"])
        processes.append(
            {
                "key": entry.key,
                "process": entry.process,
                "path": str(Path("subprocesses") / entry.key),
                "manifest": str(manifest_path.relative_to(output_path)),
                "kind": GENERIC_DAG_PROCESS_ARTIFACT_KIND,
                "artifact_class": "generic-dag-schema-v2",
                "planning_status": payload["planning_status"],
                "lowering_status": payload["lowering_status"],
                "runtime_available": bool(compiled_payload.get("runtime_available")),
                "generation_request": generation_metadata,
            }
        )
    runtime_available = all(
        bool(process["runtime_available"]) for process in processes
    )
    payload = {
        "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
        "kind": GENERIC_DAG_PROCESS_SET_ARTIFACT_KIND,
        "artifact_class": "generic-dag-schema-v2",
        "request": process_set.request,
        "default_process_key": process_set.default_key,
        "color_accuracy": process_set.color_accuracy,
        "generic_generation": generation_metadata,
        "runtime_available": runtime_available,
        "runtime_unavailable_message": (
            None
            if runtime_available
            else "one or more subprocesses do not contain serialized generic stage evaluators"
        ),
        "process_set_ir": process_set.to_json_dict(),
        "processes": processes,
    }
    manifest_path = output_path / "process_set_manifest.json"
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (output_path / "check_standalone.py").write_text(
        _GENERIC_DAG_CHECK_STANDALONE,
        encoding="utf-8",
    )
    return manifest_path, payload


def _generic_process_set_generation_metadata(
    *,
    options: ProcessOptions | None,
    selected_color_sector_ids: set[int] | None,
    lc_topology_replay: bool,
    max_coupling_orders: Mapping[str, int] | None,
    max_lc_current_line_groups: int | None,
    max_quark_pairs: int | None,
    closure_side_mask_pruning: bool,
    color_order_mask_pruning: bool,
    species_reachability_pruning: bool,
    ignored_particle_ids: Sequence[int] | None,
    ignored_vertex_kinds: Sequence[int] | None,
    selected_source_helicities: Mapping[int, int] | None,
    numerical_filter_current: bool,
    numerical_current_merging: bool,
    numerical_current_samples: int,
    numerical_current_seed: int,
    numerical_current_relative_tolerance: float,
    numerical_current_zero_tolerance: float,
) -> dict[str, object]:
    process_options = options or ProcessOptions()
    return {
        "process_options": {
            "flavour_scheme": int(process_options.flavour_scheme),
            "include_3qqbar": bool(process_options.include_3qqbar),
            "include_cc": bool(process_options.include_cc),
            "include_resonance": bool(process_options.include_resonance),
            "serial": bool(process_options.serial),
        },
        "pruning": {
            "max_coupling_orders": {
                str(name).upper(): int(value)
                for name, value in (max_coupling_orders or {}).items()
            },
            "max_lc_current_line_groups": max_lc_current_line_groups,
            "max_quark_pairs": max_quark_pairs,
            "closure_side_mask_pruning": bool(closure_side_mask_pruning),
            "color_order_mask_pruning": bool(color_order_mask_pruning),
            "species_reachability_pruning": bool(species_reachability_pruning),
            "ignored_particle_ids": [
                int(particle_id) for particle_id in (ignored_particle_ids or ())
            ],
            "ignored_vertex_kinds": [
                int(kind) for kind in (ignored_vertex_kinds or ())
            ],
            "selected_color_sector_ids": (
                None
                if selected_color_sector_ids is None
                else sorted(int(sector_id) for sector_id in selected_color_sector_ids)
            ),
            "selected_source_helicities": (
                None
                if selected_source_helicities is None
                else {
                    str(int(label)): int(helicity)
                    for label, helicity in sorted(
                        selected_source_helicities.items()
                    )
                }
            ),
            "zero_current_filter": {
                "enabled": bool(numerical_filter_current),
                "sample_count": int(numerical_current_samples),
                "seed": int(numerical_current_seed),
                "relative_tolerance": float(numerical_current_relative_tolerance),
                "zero_tolerance": float(numerical_current_zero_tolerance),
            },
            "current_merging": {
                "enabled": bool(numerical_current_merging),
                "sample_count": int(numerical_current_samples),
                "seed": int(numerical_current_seed),
                "relative_tolerance": float(numerical_current_relative_tolerance),
                "zero_tolerance": float(numerical_current_zero_tolerance),
            },
        },
        "lc_topology_replay": bool(lc_topology_replay),
    }


def _aggregate_lc_gluon_flavour_flow_currents(
    dag: GenericDAG,
    model: Model,
    *,
    sample_count: int,
    seed: int,
) -> tuple[GenericDAG, Mapping[str, object]]:
    before = _dag_count_payload(dag)
    report: dict[str, object] = {
        "enabled": dag.process.color_accuracy == "lc",
        "mode": "lc-gluon-flavour-flow-aggregation",
        "before": before,
    }
    if dag.process.color_accuracy != "lc":
        report.update(
            {
                "after": before,
                "skipped": True,
                "reason": "structural current aggregation is currently defined for LC DAGs",
            }
        )
        return dag, report

    def current_key(current: CurrentNode) -> tuple[object, ...]:
        index = current.index
        if current.is_source:
            return ("source", current.id)
        try:
            is_gluon = model.color_rep(int(index.particle_id)) == 8
        except KeyError:
            is_gluon = False
        if not is_gluon or index.auxiliary_kind is not None:
            return ("current", current.id)
        return (
            "lc-gluon",
            index.particle_id,
            index.external_mask,
            index.helicity_ancestry,
            index.chirality,
            index.spin_state,
            index.color_state,
            index.momentum_mask,
            index.coupling_orders,
            index.auxiliary_kind,
            index.ordered_external_labels,
        )

    groups: dict[tuple[object, ...], list[int]] = {}
    for current in dag.currents:
        groups.setdefault(current_key(current), []).append(current.id)
    merged_groups = tuple(tuple(ids) for ids in groups.values() if len(ids) > 1)
    merged_current_count = sum(len(ids) - 1 for ids in merged_groups)
    if merged_current_count == 0:
        report.update(
            {
                "after": before,
                "skipped": False,
                "merged_group_count": 0,
                "merged_current_count": 0,
                "removed_interaction_count": 0,
                "removed_amplitude_root_count": 0,
            }
        )
        return dag, report

    representative_by_old_id: dict[int, int] = {}
    for ids in groups.values():
        representative = min(ids)
        for current_id in ids:
            representative_by_old_id[current_id] = representative
    representatives = tuple(
        current.id
        for current in dag.currents
        if representative_by_old_id[current.id] == current.id
    )
    new_id_by_old_representative = {
        old_id: new_id for new_id, old_id in enumerate(representatives)
    }

    def mapped_current_id(old_id: int) -> int:
        return new_id_by_old_representative[representative_by_old_id[old_id]]

    currents = tuple(
        replace(dag.currents[old_id], id=new_id_by_old_representative[old_id])
        for old_id in representatives
    )
    interactions: list[InteractionNode] = []
    seen_interactions: set[tuple[object, ...]] = set()
    skipped_interaction_count = 0
    for interaction in dag.interactions:
        mapped = replace(
            interaction,
            id=0,
            left_id=mapped_current_id(interaction.left_id),
            right_id=mapped_current_id(interaction.right_id),
            result_id=mapped_current_id(interaction.result_id),
        )
        interaction_key = (
            mapped.vertex_kind,
            mapped.vertex_particles,
            mapped.left_id,
            mapped.right_id,
            mapped.result_id,
            mapped.coupling,
            mapped.color_weight,
            mapped.lowering_backend,
            mapped.full_tensor_network_ready,
        )
        if interaction_key in seen_interactions:
            skipped_interaction_count += 1
            continue
        seen_interactions.add(interaction_key)
        interactions.append(replace(mapped, id=len(interactions)))

    roots: list[AmplitudeRoot] = []
    seen_roots: set[tuple[object, ...]] = set()
    skipped_root_count = 0
    for root in dag.amplitude_roots:
        mapped = replace(
            root,
            id=0,
            left_id=mapped_current_id(root.left_id),
            right_id=mapped_current_id(root.right_id),
        )
        root_key = (
            mapped.kind,
            mapped.left_id,
            mapped.right_id,
            mapped.color_weight,
            mapped.color_sector_id,
            mapped.vertex_kind,
            mapped.vertex_particles,
            mapped.coupling,
            mapped.contraction,
            mapped.helicity_weight,
        )
        if root_key in seen_roots:
            skipped_root_count += 1
            continue
        seen_roots.add(root_key)
        roots.append(replace(mapped, id=len(roots)))

    sources: list[int] = []
    seen_sources: set[int] = set()
    for source_id in dag.sources:
        mapped_source = mapped_current_id(source_id)
        if mapped_source in seen_sources:
            continue
        seen_sources.add(mapped_source)
        sources.append(mapped_source)

    candidate = prune_dag_to_amplitude_roots(
        GenericDAG(
            process=dag.process,
            color_plan=dag.color_plan,
            currents=currents,
            sources=tuple(sources),
            interactions=tuple(interactions),
            amplitude_roots=tuple(roots),
            truncated=dag.truncated,
        )
    )
    validation = _validate_numerical_rewrite_preserves_raw_sums(
        dag,
        candidate,
        model,
        sample_count=sample_count,
        seed=seed,
    )
    accepted = bool(validation.get("accepted", False))
    effective_dag = candidate if accepted else dag
    after = _dag_count_payload(effective_dag)
    candidate_after = _dag_count_payload(candidate)
    report.update(
        {
            "after": after,
            "candidate_after": candidate_after,
            "merged_group_count": len(merged_groups),
            "merged_current_count": merged_current_count,
            "merged_current_ids": [list(ids[1:]) for ids in merged_groups],
            "removed_interaction_count": (
                before["interaction_count"] - after["interaction_count"]
            ),
            "removed_amplitude_root_count": (
                before["amplitude_root_count"] - after["amplitude_root_count"]
            ),
            "deduplicated_interaction_count": skipped_interaction_count,
            "deduplicated_amplitude_root_count": skipped_root_count,
            "validation": validation,
        }
    )
    if not accepted:
        report.update(
            {
                "skipped": True,
                "reason": "candidate structural aggregation failed raw-sum validation; leaving DAG unchanged",
            }
        )
        return dag, report
    report["skipped"] = False
    return candidate, report


def _filter_zero_currents_by_warmup(
    dag: GenericDAG,
    model: Model,
    *,
    sample_count: int,
    seed: int,
    relative_tolerance: float,
    zero_tolerance: float,
) -> tuple[GenericDAG, dict[str, object]]:
    before = _dag_count_payload(dag)
    report: dict[str, object] = {
        "enabled": True,
        "mode": "numeric-current-warmup",
        "sample_count": int(sample_count),
        "seed": int(seed),
        "relative_tolerance": float(relative_tolerance),
        "zero_tolerance": float(zero_tolerance),
        "before": before,
    }
    if sample_count <= 0:
        report.update(
            {
                "skipped": True,
                "reason": "sample_count must be positive",
                "after": before,
            }
        )
        return dag, report
    if not dag.interactions:
        report.update(
            {
                "skipped": False,
                "zero_current_ids": [],
                "removed_current_count": 0,
                "removed_interaction_count": 0,
                "removed_amplitude_root_count": 0,
                "after": before,
            }
        )
        return dag, report

    try:
        points = tuple(
            _generic_warmup_phase_space_point(dag, model, seed=seed + offset)
            for offset in range(sample_count)
        )
        maxima = _evaluate_current_component_maxima(dag, model, points)
    except Exception as error:  # noqa: BLE001 - filtering must be conservative.
        report.update(
            {
                "skipped": True,
                "reason": f"numeric warmup evaluation failed: {error}",
                "after": before,
            }
        )
        return dag, report

    source_ids = set(dag.sources)
    global_max = max(maxima.values(), default=0.0)
    threshold = max(float(zero_tolerance), abs(float(relative_tolerance)) * global_max)
    zero_ids = tuple(
        current.id
        for current in dag.currents
        if current.id not in source_ids and maxima.get(current.id, 0.0) <= threshold
    )
    if not zero_ids:
        report.update(
            {
                "skipped": False,
                "threshold": threshold,
                "global_max_abs": global_max,
                "zero_current_ids": [],
                "removed_current_count": 0,
                "removed_interaction_count": 0,
                "removed_amplitude_root_count": 0,
                "after": before,
            }
        )
        return dag, report

    filtered = _drop_currents_from_dag(dag, zero_ids)
    if not filtered.amplitude_roots:
        report.update(
            {
                "skipped": True,
                "reason": (
                    "zero-current candidates would remove every amplitude root; "
                    "leaving DAG unchanged"
                ),
                "threshold": threshold,
                "global_max_abs": global_max,
                "zero_current_ids": list(zero_ids),
                "after": before,
            }
        )
        return dag, report

    after = _dag_count_payload(filtered)
    validation = _validate_numerical_rewrite_preserves_raw_sums(
        dag,
        filtered,
        model,
        sample_count=sample_count,
        seed=seed,
    )
    if not bool(validation["accepted"]):
        report.update(
            {
                "skipped": True,
                "reason": (
                    "zero-current candidates failed raw-sum validation; "
                    "leaving DAG unchanged"
                ),
                "threshold": threshold,
                "global_max_abs": global_max,
                "zero_current_ids": list(zero_ids),
                "zero_current_count": len(zero_ids),
                "candidate_after": after,
                "validation": validation,
                "after": before,
            }
        )
        return dag, report

    report.update(
        {
            "skipped": False,
            "threshold": threshold,
            "global_max_abs": global_max,
            "zero_current_ids": list(zero_ids),
            "zero_current_count": len(zero_ids),
            "removed_current_count": before["current_count"] - after["current_count"],
            "removed_interaction_count": (
                before["interaction_count"] - after["interaction_count"]
            ),
            "removed_amplitude_root_count": (
                before["amplitude_root_count"] - after["amplitude_root_count"]
            ),
            "validation": validation,
            "after": after,
        }
    )
    return filtered, report


def _merge_identical_currents_by_warmup(
    dag: GenericDAG,
    model: Model,
    *,
    sample_count: int,
    seed: int,
    relative_tolerance: float,
    zero_tolerance: float,
) -> tuple[GenericDAG, dict[str, object]]:
    before = _dag_count_payload(dag)
    report: dict[str, object] = {
        "enabled": True,
        "mode": "numeric-current-signature-warmup",
        "sample_count": int(sample_count),
        "seed": int(seed),
        "relative_tolerance": float(relative_tolerance),
        "zero_tolerance": float(zero_tolerance),
        "before": before,
    }
    if sample_count <= 0:
        report.update(
            {
                "skipped": True,
                "reason": "sample_count must be positive",
                "after": before,
            }
        )
        return dag, report
    if not dag.interactions:
        report.update(
            {
                "skipped": False,
                "merged_current_ids": [],
                "merged_current_count": 0,
                "merged_group_count": 0,
                "removed_current_count": 0,
                "removed_interaction_count": 0,
                "removed_amplitude_root_count": 0,
                "after": before,
            }
        )
        return dag, report

    try:
        points = tuple(
            _generic_warmup_phase_space_point(dag, model, seed=seed + offset)
            for offset in range(sample_count)
        )
        maxima, signatures = _evaluate_current_warmup(dag, model, points)
    except Exception as error:  # noqa: BLE001 - merging must be conservative.
        report.update(
            {
                "skipped": True,
                "reason": f"numeric warmup evaluation failed: {error}",
                "after": before,
            }
        )
        return dag, report

    global_max = max(maxima.values(), default=0.0)
    threshold = max(float(zero_tolerance), abs(float(relative_tolerance)) * global_max)
    current_map, groups = _equivalent_current_representatives(
        dag,
        maxima=maxima,
        signatures=signatures,
        absolute_tolerance=threshold,
        relative_tolerance=abs(float(relative_tolerance)),
    )
    merged_ids = tuple(
        current_id
        for current_id, (representative_id, _sign) in sorted(current_map.items())
        if current_id != representative_id
    )
    if not merged_ids:
        report.update(
            {
                "skipped": False,
                "threshold": threshold,
                "global_max_abs": global_max,
                "merged_current_ids": [],
                "merged_current_count": 0,
                "merged_group_count": 0,
                "removed_current_count": 0,
                "removed_interaction_count": 0,
                "removed_amplitude_root_count": 0,
                "after": before,
            }
        )
        return dag, report

    merged = _merge_currents_in_dag(dag, current_map)
    if not merged.amplitude_roots:
        report.update(
            {
                "skipped": True,
                "reason": (
                    "current-merging candidates would remove every amplitude "
                    "root; leaving DAG unchanged"
                ),
                "threshold": threshold,
                "global_max_abs": global_max,
                "merged_current_ids": list(merged_ids),
                "after": before,
            }
        )
        return dag, report

    after = _dag_count_payload(merged)
    validation = _validate_numerical_rewrite_preserves_raw_sums(
        dag,
        merged,
        model,
        sample_count=sample_count,
        seed=seed,
    )
    if not bool(validation["accepted"]):
        report.update(
            {
                "skipped": True,
                "reason": (
                    "current-merging candidates failed raw-sum validation; "
                    "leaving DAG unchanged"
                ),
                "threshold": threshold,
                "global_max_abs": global_max,
                "merged_current_ids": list(merged_ids),
                "merged_current_count": len(merged_ids),
                "merged_group_count": len(groups),
                "candidate_after": after,
                "validation": validation,
                "after": before,
            }
        )
        return dag, report

    report.update(
        {
            "skipped": False,
            "threshold": threshold,
            "global_max_abs": global_max,
            "merged_current_ids": list(merged_ids),
            "merged_current_count": len(merged_ids),
            "merged_group_count": len(groups),
            "merged_groups": [
                {
                    "representative_current_id": representative_id,
                    "member_current_ids": [member_id for member_id, _sign in member_ids],
                    "member_signs": [
                        {
                            "current_id": member_id,
                            "sign": sign,
                        }
                        for member_id, sign in member_ids
                    ],
                }
                for representative_id, member_ids in groups
            ],
            "removed_current_count": before["current_count"] - after["current_count"],
            "removed_interaction_count": (
                before["interaction_count"] - after["interaction_count"]
            ),
            "removed_amplitude_root_count": (
                before["amplitude_root_count"] - after["amplitude_root_count"]
            ),
            "validation": validation,
            "after": after,
        }
    )
    return merged, report


def _dag_count_payload(dag: GenericDAG) -> dict[str, int]:
    return {
        "current_count": len(dag.currents),
        "source_count": len(dag.sources),
        "interaction_count": len(dag.interactions),
        "amplitude_root_count": len(dag.amplitude_roots),
    }


def _generic_warmup_phase_space_point(
    dag: GenericDAG,
    model: Model,
    *,
    seed: int,
) -> tuple[ExternalMomentum, ...]:
    initial_pdgs = tuple(int(pdg) for pdg in dag.process.initial_pdgs)
    final_pdgs = tuple(int(pdg) for pdg in dag.process.final_pdgs)
    if len(initial_pdgs) != 2:
        raise NativeEvaluationError(
            "zero-current filter currently requires a two-body initial state"
        )
    if not final_pdgs:
        raise NativeEvaluationError(
            "zero-current filter requires at least one final-state particle"
        )
    final_masses = tuple(float(_model_mass(model, pdg)) for pdg in final_pdgs)
    threshold = sum(final_masses)
    final_momenta: tuple[tuple[float, float, float, float], ...]
    if len(final_pdgs) == 1:
        if threshold <= 0.0:
            raise NativeEvaluationError(
                "no finite centre-of-mass warmup point exists for one massless final state"
            )
        sqrt_s = threshold
        final_momenta = ((float(sqrt_s), 0.0, 0.0, 0.0),)
    else:
        sqrt_s = max(1000.0, threshold + 100.0)
        final_momenta = massive_rambo_final_state(
            len(final_pdgs),
            sqrt_s=float(sqrt_s),
            masses=final_masses,
            seed=int(seed),
        )
    beam_energy = 0.5 * float(sqrt_s)
    return (
        ExternalMomentum(initial_pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
        ExternalMomentum(initial_pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
        *(
            ExternalMomentum(pdg, momentum)
            for pdg, momentum in zip(final_pdgs, final_momenta, strict=True)
        ),
    )


def _evaluate_current_component_maxima(
    dag: GenericDAG,
    model: Model,
    points: Sequence[Sequence[ExternalMomentum]],
) -> dict[int, float]:
    maxima, _ = _evaluate_current_warmup(
        dag,
        model,
        points,
        collect_signatures=False,
    )
    return maxima


def _evaluate_current_warmup(
    dag: GenericDAG,
    model: Model,
    points: Sequence[Sequence[ExternalMomentum]],
    *,
    collect_signatures: bool = True,
) -> tuple[
    dict[int, float],
    dict[int, tuple[tuple[complex, ...], ...]],
]:
    from .generic_stage_compiler import (
        _interaction_contribution,
        _momentum_components,
        _sum_components,
    )

    schema = cast(dict[str, Any], _generic_runtime_schema_payload(dag, model))
    value_slots = _schema_value_slots_by_id(schema)
    momentum_slots = _schema_momentum_slots_by_id(schema)
    current_slots = _schema_current_slots_by_id(schema)
    stages = cast(list[Any], schema["stages"])
    layout = cast(dict[str, Any], schema["parameter_layout"])
    value_count = int(layout["value_component_count"])
    momentum_count = int(layout["momentum_parameter_count"])
    stage_runtime: list[
        tuple[
            dict[str, Any],
            tuple[int, ...],
            dict[int, list[dict[str, Any]]],
            dict[int, list[dict[str, Any]]],
        ]
    ] = []
    for stage_obj in stages:
        stage = cast(dict[str, Any], stage_obj)
        slots_by_current: dict[int, list[dict[str, Any]]] = {}
        for slot_id in cast(list[Any], stage["output_value_slot_ids"]):
            slot = value_slots[int(slot_id)]
            slots_by_current.setdefault(int(slot["current_id"]), []).append(slot)
        interactions_by_result: dict[int, list[dict[str, Any]]] = {}
        for raw in cast(list[Any], stage["interactions"]):
            interaction = cast(dict[str, Any], raw)
            interactions_by_result.setdefault(
                int(interaction["result_current_id"]),
                [],
            ).append(interaction)
        stage_runtime.append(
            (
                stage,
                tuple(sorted(interactions_by_result)),
                interactions_by_result,
                slots_by_current,
            )
        )
    maxima = {current.id: 0.0 for current in dag.currents}
    signature_parts: dict[int, list[tuple[complex, ...]]] = (
        {current.id: [] for current in dag.currents} if collect_signatures else {}
    )
    for point in points:
        values = [0j for _ in range(value_count)]
        momenta = [0j for _ in range(momentum_count)]
        _fill_warmup_sources(schema, model, values, point)
        _fill_warmup_momenta(schema, momenta, point)
        for (
            _stage,
            current_ids,
            interactions_by_result,
            output_value_slots_by_current,
        ) in stage_runtime:
            for current_id in current_ids:
                current_slot = current_slots[current_id]
                dimension = int(current_slot["dimension"])
                total = tuple(0j for _ in range(dimension))
                for interaction in interactions_by_result[current_id]:
                    contribution = _interaction_contribution(
                        dag,
                        model,
                        interaction,
                        value_symbols=values,
                        momentum_symbols=momenta,
                        model_parameter_symbols={},
                        value_slots=value_slots,
                        momentum_slots=momentum_slots,
                    )
                    total = _sum_components(total, contribution)
                total_signature = tuple(complex(component) for component in total)
                if collect_signatures:
                    signature_parts[current_id].append(total_signature)
                for slot in output_value_slots_by_current.get(current_id, ()):
                    components = (
                        model.propagator_component_expression(
                            int(current_slot["particle_id"]),
                            total,
                            _momentum_components(
                                int(current_slot["momentum_mask"]),
                                momenta,
                                momentum_slots,
                            ),
                            chirality=int(current_slot["chirality"]),
                        )
                        if str(slot["variant"]) == "propagated"
                        else total
                    )
                    start = int(slot["component_start"])
                    for offset, component in enumerate(components):
                        value = complex(component)
                        values[start + offset] = value
                        maxima[current_id] = max(maxima[current_id], abs(value))
    return maxima, {
        current_id: tuple(parts)
        for current_id, parts in signature_parts.items()
    }


def _validate_numerical_rewrite_preserves_raw_sums(
    reference_dag: GenericDAG,
    candidate_dag: GenericDAG,
    model: Model,
    *,
    sample_count: int,
    seed: int,
) -> dict[str, object]:
    validation_count = max(
        1,
        min(
            int(sample_count),
            _NUMERICAL_REWRITE_VALIDATION_SAMPLE_LIMIT,
        ),
    )
    validation_seed = int(seed) + _NUMERICAL_REWRITE_VALIDATION_SEED_OFFSET
    report: dict[str, object] = {
        "mode": "independent-raw-sum-warmup",
        "accepted": False,
        "sample_count": validation_count,
        "seed": validation_seed,
        "relative_tolerance": _NUMERICAL_REWRITE_VALIDATION_RELATIVE_TOLERANCE,
        "zero_tolerance": _NUMERICAL_REWRITE_VALIDATION_ZERO_TOLERANCE,
    }
    try:
        points = tuple(
            _generic_warmup_phase_space_point(
                reference_dag,
                model,
                seed=validation_seed + offset,
            )
            for offset in range(validation_count)
        )
        reference_sums = _evaluate_dag_raw_sums_by_warmup(
            reference_dag,
            model,
            points,
        )
        candidate_sums = _evaluate_dag_raw_sums_by_warmup(
            candidate_dag,
            model,
            points,
        )
    except Exception as error:  # noqa: BLE001 - pruning must be conservative.
        report.update(
            {
                "reason": f"raw-sum validation failed: {error}",
                "reference_raw_sums": [],
                "candidate_raw_sums": [],
            }
        )
        return report

    differences = [
        abs(float(reference) - float(candidate))
        for reference, candidate in zip(reference_sums, candidate_sums, strict=True)
    ]
    scales = [
        max(abs(float(reference)), abs(float(candidate)), 1.0e-300)
        for reference, candidate in zip(reference_sums, candidate_sums, strict=True)
    ]
    relative_differences = [
        difference / scale
        for difference, scale in zip(differences, scales, strict=True)
    ]
    tolerances = [
        (
            _NUMERICAL_REWRITE_VALIDATION_ZERO_TOLERANCE
            if scale <= _NUMERICAL_REWRITE_VALIDATION_ZERO_TOLERANCE
            else _NUMERICAL_REWRITE_VALIDATION_RELATIVE_TOLERANCE * scale
        )
        for scale in scales
    ]
    accepted = all(
        difference <= tolerance
        for difference, tolerance in zip(differences, tolerances, strict=True)
    )
    report.update(
        {
            "accepted": accepted,
            "max_abs_difference": max(differences, default=0.0),
            "max_relative_difference": max(relative_differences, default=0.0),
            "reference_raw_sums": [float(value) for value in reference_sums],
            "candidate_raw_sums": [float(value) for value in candidate_sums],
        }
    )
    if not accepted:
        report["reason"] = "candidate raw sums differ from unmodified DAG"
    return report


def _evaluate_dag_raw_sums_by_warmup(
    dag: GenericDAG,
    model: Model,
    points: Sequence[Sequence[ExternalMomentum]],
) -> tuple[float, ...]:
    from .generic_stage_compiler import (
        _amplitude_root_expression,
        _interaction_contribution,
        _momentum_components,
        _sum_components,
    )

    schema = cast(dict[str, Any], _generic_runtime_schema_payload(dag, model))
    value_slots = _schema_value_slots_by_id(schema)
    momentum_slots = _schema_momentum_slots_by_id(schema)
    current_slots = _schema_current_slots_by_id(schema)
    stages = cast(list[Any], schema["stages"])
    amplitude_stage = cast(dict[str, Any], schema["amplitude_stage"])
    roots = [
        cast(dict[str, Any], raw)
        for raw in cast(list[Any], amplitude_stage["roots"])
    ]
    layout = cast(dict[str, Any], schema["parameter_layout"])
    value_count = int(layout["value_component_count"])
    momentum_count = int(layout["momentum_parameter_count"])
    raw_sums: list[float] = []
    for point in points:
        values = [0j for _ in range(value_count)]
        momenta = [0j for _ in range(momentum_count)]
        _fill_warmup_sources(schema, model, values, point)
        _fill_warmup_momenta(schema, momenta, point)
        for stage_obj in stages:
            stage = cast(dict[str, Any], stage_obj)
            interactions_by_result: dict[int, list[dict[str, Any]]] = {}
            for raw in cast(list[Any], stage["interactions"]):
                interaction = cast(dict[str, Any], raw)
                interactions_by_result.setdefault(
                    int(interaction["result_current_id"]),
                    [],
                ).append(interaction)
            for current_id in sorted(interactions_by_result):
                current_slot = current_slots[current_id]
                dimension = int(current_slot["dimension"])
                total = tuple(0j for _ in range(dimension))
                for interaction in interactions_by_result[current_id]:
                    contribution = _interaction_contribution(
                        dag,
                        model,
                        interaction,
                        value_symbols=values,
                        momentum_symbols=momenta,
                        model_parameter_symbols={},
                        value_slots=value_slots,
                        momentum_slots=momentum_slots,
                    )
                    total = _sum_components(total, contribution)
                for slot_id in cast(list[Any], stage["output_value_slot_ids"]):
                    slot = value_slots[int(slot_id)]
                    if int(slot["current_id"]) != current_id:
                        continue
                    components = (
                        model.propagator_component_expression(
                            int(current_slot["particle_id"]),
                            total,
                            _momentum_components(
                                int(current_slot["momentum_mask"]),
                                momenta,
                                momentum_slots,
                            ),
                            chirality=int(current_slot["chirality"]),
                        )
                        if str(slot["variant"]) == "propagated"
                        else total
                    )
                    start = int(slot["component_start"])
                    for offset, component in enumerate(components):
                        values[start + offset] = complex(component)
        amplitudes = tuple(
            complex(
                _amplitude_root_expression(
                    model,
                    root,
                    value_symbols=values,
                    model_parameter_symbols={},
                    value_slots=value_slots,
                )
            )
            for root in roots
        )
        weights = tuple(float(root["helicity_weight"]) for root in roots)
        group_ids = tuple(
            None
            if root.get("coherent_group_id") is None
            else int(root["coherent_group_id"])
            for root in roots
        )
        raw_sums.append(
            _coherent_weighted_abs2_sum(
                amplitudes,
                weights,
                group_ids,
            )
        )
    return tuple(raw_sums)


def _coherent_weighted_abs2_sum(
    amplitudes: Sequence[complex],
    weights: Sequence[float],
    group_ids: Sequence[int | None],
) -> float:
    if len(amplitudes) != len(weights) or len(amplitudes) != len(group_ids):
        raise NativeEvaluationError(
            "raw-sum validation metadata length does not match amplitude outputs"
        )
    raw_sum = 0.0
    grouped: dict[int, list[int]] = {}
    for index, group_id in enumerate(group_ids):
        if group_id is None:
            amplitude = amplitudes[index]
            raw_sum += float(weights[index]) * (
                amplitude.real * amplitude.real + amplitude.imag * amplitude.imag
            )
        else:
            grouped.setdefault(int(group_id), []).append(index)
    for group_id in sorted(grouped):
        indices = grouped[group_id]
        group_weights = {float(weights[index]) for index in indices}
        if len(group_weights) != 1:
            raise NativeEvaluationError(
                f"coherent amplitude group {group_id} has inconsistent raw-sum weights"
            )
        coherent = sum((amplitudes[index] for index in indices), start=0j)
        raw_sum += group_weights.pop() * (
            coherent.real * coherent.real + coherent.imag * coherent.imag
        )
    return float(raw_sum)


def _fill_warmup_sources(
    schema: Mapping[str, Any],
    model: Model,
    values: list[complex],
    point: Sequence[ExternalMomentum],
) -> None:
    source_fill = cast(dict[str, Any], schema["source_fill"])
    for raw_source in cast(list[Any], source_fill["sources"]):
        source = cast(dict[str, Any], raw_source)
        slot = cast(dict[str, Any], source["value_slot"])
        wave = _warmup_source_wavefunction(source, model, point)
        start = int(slot["component_start"])
        stop = int(slot["component_stop"])
        if stop - start != len(wave):
            raise NativeEvaluationError(
                f"source {source['source_id']} produced {len(wave)} components "
                f"for slot length {stop - start}"
            )
        values[start:stop] = [complex(component) for component in wave]


def _fill_warmup_momenta(
    schema: Mapping[str, Any],
    momenta: list[complex],
    point: Sequence[ExternalMomentum],
) -> None:
    conventions = cast(dict[str, Any], schema["momentum_conventions"])
    incoming = {int(label) for label in cast(list[Any], conventions["incoming_labels"])}
    for raw_slot in cast(list[Any], schema["momentum_slots"]):
        slot = cast(dict[str, Any], raw_slot)
        total = [0.0, 0.0, 0.0, 0.0]
        for raw_label in cast(list[Any], slot["external_labels"]):
            label = int(raw_label)
            try:
                physical = point[label - 1].momentum
            except IndexError as error:
                raise NativeEvaluationError(
                    f"momentum slot refers to unknown external label {label}"
                ) from error
            sign = -1.0 if label in incoming else 1.0
            for component in range(4):
                total[component] += sign * float(physical[component])
        start = int(slot["component_start"])
        for component in range(4):
            momenta[start + component] = complex(total[component], 0.0)


def _warmup_source_wavefunction(
    source: Mapping[str, Any],
    model: Model,
    point: Sequence[ExternalMomentum],
) -> tuple[complex, ...]:
    leg_label = int(source["leg_label"])
    try:
        momentum = tuple(float(component) for component in point[leg_label - 1].momentum)
    except IndexError as error:
        raise NativeEvaluationError(
            f"source refers to unknown external label {leg_label}"
        ) from error
    if str(source["crossing"]) == "negate-incoming-momentum":
        momentum = tuple(-component for component in momentum)
    dimension = int(source["dimension"])
    particle_id = int(source["particle_id"])
    helicity = int(source["source_helicity"])
    chirality = int(source["chirality"])
    if dimension == 1:
        return (1.0 + 0.0j,)
    if dimension == 2:
        return (
            _ext_antiquark_weyl_filter(momentum, helicity, chirality)
            if particle_id < 0
            else _ext_quark_weyl_filter(momentum, helicity, chirality)
        )
    if dimension == 4:
        if _is_fermion_pdg(particle_id):
            mass = _model_mass(model, particle_id)
            return (
                _ext_antiquark_dirac_massive(momentum, helicity, mass)
                if particle_id < 0
                else _ext_quark_dirac_massive(momentum, helicity, mass)
            )
        if abs(particle_id) == 21 or particle_id == 22:
            return _ext_gluon_filter(momentum, helicity)
        return _ext_massive_vector_filter(
            momentum,
            helicity,
            _model_mass(model, particle_id),
        )
    raise NativeEvaluationError(
        f"source dimension {dimension} is not implemented by zero-current filter"
    )


def _ext_quark_weyl_filter(
    momentum: Sequence[float],
    helicity: int,
    chirality: int,
) -> tuple[complex, complex]:
    energy, px, py, pz = (float(component) for component in momentum)
    if energy > 0.0:
        sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz < 0.0 else math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(helicity * px, -py) / sqp0p3
        )
        if helicity == 1 and chirality == 1:
            return chi1, chi2
        if helicity == -1 and chirality == -1:
            return chi2, chi1
        return 0j, 0j
    sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz > 0.0 else -math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(-helicity * (-px), py) / sqp0p3
    )
    if helicity == -1 and chirality == 1:
        return chi1, chi2
    if helicity == 1 and chirality == -1:
        return chi2, chi1
    return 0j, 0j


def _ext_antiquark_weyl_filter(
    momentum: Sequence[float],
    helicity: int,
    chirality: int,
) -> tuple[complex, complex]:
    energy, px, py, pz = (float(component) for component in momentum)
    if energy > 0.0:
        sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz < 0.0 else -math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(-helicity * px, py) / sqp0p3
        )
        if helicity == 1 and chirality == 1:
            return chi2, chi1
        if helicity == -1 and chirality == -1:
            return chi1, chi2
        return 0j, 0j
    sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz > 0.0 else math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(helicity * (-px), -py) / sqp0p3
    )
    if helicity == -1 and chirality == 1:
        return chi2, chi1
    if helicity == 1 and chirality == -1:
        return chi1, chi2
    return 0j, 0j


def _ext_quark_dirac_filter(
    momentum: Sequence[float],
    helicity: int,
) -> tuple[complex, complex, complex, complex]:
    energy, px, py, pz = (float(component) for component in momentum)
    if energy > 0.0:
        sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz < 0.0 else math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(helicity * px, -py) / sqp0p3
        )
        if helicity == 1:
            return chi1, chi2, 0j, 0j
        return 0j, 0j, chi2, chi1
    sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz > 0.0 else -math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(-helicity * (-px), py) / sqp0p3
    )
    if -helicity == 1:
        return chi1, chi2, 0j, 0j
    return 0j, 0j, chi2, chi1


def _ext_antiquark_dirac_filter(
    momentum: Sequence[float],
    helicity: int,
) -> tuple[complex, complex, complex, complex]:
    energy, px, py, pz = (float(component) for component in momentum)
    if energy > 0.0:
        sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz < 0.0 else -math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(-helicity * px, py) / sqp0p3
        )
        if -helicity == 1:
            return 0j, 0j, chi1, chi2
        return chi2, chi1, 0j, 0j
    sqp0p3 = 0.0 if px == 0.0 and py == 0.0 and pz > 0.0 else math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(helicity * (-px), -py) / sqp0p3
    )
    if helicity == 1:
        return 0j, 0j, chi1, chi2
    return chi2, chi1, 0j, 0j


def _ext_gluon_filter(
    momentum: Sequence[float],
    helicity: int,
) -> tuple[complex, complex, complex, complex]:
    energy, px, py, pz = (float(component) for component in momentum)
    if energy == 0.0:
        raise NativeEvaluationError("cannot generate external gluon with zero energy")
    sqh = math.sqrt(0.5)
    if energy > 0.0:
        hel = float(helicity)
        pp = energy
        pt = math.sqrt(px * px + py * py)
        wf3 = complex(hel * pt / pp * sqh)
        if pt != 0.0:
            pzpt = pz / (pp * pt) * sqh * hel
            wf1 = complex(-px * pzpt, -py / pt * sqh)
            wf2 = complex(-py * pzpt, px / pt * sqh)
        else:
            wf1 = complex(-hel * sqh)
            wf2 = complex(0.0, _fortran_sign(sqh, pz))
        return 0j, wf1, wf2, wf3
    hel = float(-helicity)
    pp = -energy
    pt = math.sqrt(px * px + py * py)
    wf3 = complex(hel * pt / pp * sqh)
    if pt != 0.0:
        pzpt = -pz / (pp * pt) * sqh * hel
        wf1 = complex(px * pzpt, py / pt * sqh)
        wf2 = complex(py * pzpt, -px / pt * sqh)
    else:
        wf1 = complex(-hel * sqh)
        wf2 = complex(0.0, -_fortran_sign(sqh, pz))
    return 0j, wf1, wf2, wf3


def _ext_massive_vector_filter(
    momentum: Sequence[float],
    helicity: int,
    mass: float,
) -> tuple[complex, complex, complex, complex]:
    energy, px, py, pz = (float(component) for component in momentum)
    if mass == 0.0:
        raise NativeEvaluationError("massless massive-vector source is invalid")
    sqh = math.sqrt(0.5)
    hel = float(helicity)
    nsvahl = abs(helicity)
    pt2 = px * px + py * py
    pp = min(energy, math.sqrt(pt2 + pz * pz))
    pt = min(pp, math.sqrt(pt2))
    hel0 = 1.0 - abs(hel)
    if pp == 0.0:
        return 0j, complex(-hel * sqh), complex(0.0, nsvahl * sqh), complex(hel0)
    emp = energy / (mass * pp)
    wf0 = complex(hel0 * pp / mass)
    wf3 = complex(hel0 * pz * emp + hel * pt / pp * sqh)
    if pt != 0.0:
        pzpt = pz / (pp * pt) * sqh * hel
        wf1 = complex(hel0 * px * emp - px * pzpt, -nsvahl * py / pt * sqh)
        wf2 = complex(hel0 * py * emp - py * pzpt, nsvahl * px / pt * sqh)
    else:
        wf1 = complex(-hel * sqh)
        wf2 = complex(0.0, nsvahl * _fortran_sign(sqh, pz))
    return wf0, wf1, wf2, wf3


def _ext_quark_dirac_massive(
    momentum: Sequence[float],
    helicity: int,
    mass: float,
) -> tuple[complex, complex, complex, complex]:
    if abs(mass) < 1.0e-8:
        return _ext_quark_dirac_filter(momentum, helicity)
    energy, px, py, pz = (float(component) for component in momentum)
    nsf = 1 if energy > 0.0 else -1
    nh = nsf * int(helicity)
    pp = math.sqrt(px * px + py * py + pz * pz)
    omega1 = math.sqrt(abs(energy) + pp)
    omega2 = mass / omega1
    omega = (omega1, omega2)
    sf1 = (1 + nsf + (1 - nsf) * nh) * 0.5
    sf2 = (1 + nsf - (1 - nsf) * nh) * 0.5
    ip = (3 + nh) // 2 - 1
    im = (3 - nh) // 2 - 1
    sfomeg = (sf1 * omega[ip], sf2 * omega[im])
    signed_px, signed_py, signed_pz = (
        (px, py, pz) if energy > 0.0 else (-px, -py, -pz)
    )
    pp3 = max(pp + signed_pz, 0.0)
    chi1 = complex(1.0 if pp == 0.0 else math.sqrt(pp3 * 0.5 / pp), 0.0)
    if pp3 == 0.0 or pp == 0.0:
        chi2 = complex(-nh, 0.0)
    else:
        denom = math.sqrt(2.0 * pp * pp3)
        chi2 = complex(nh * signed_px / denom, -signed_py / denom)
    chi = (chi1, chi2)
    return (
        chi[im] * sfomeg[1],
        chi[ip] * sfomeg[1],
        chi[im] * sfomeg[0],
        chi[ip] * sfomeg[0],
    )


def _ext_antiquark_dirac_massive(
    momentum: Sequence[float],
    helicity: int,
    mass: float,
) -> tuple[complex, complex, complex, complex]:
    if abs(mass) < 1.0e-8:
        return _ext_antiquark_dirac_filter(momentum, helicity)
    energy, px, py, pz = (float(component) for component in momentum)
    nsf = -1 if energy > 0.0 else 1
    nh = nsf * int(helicity)
    pp = math.sqrt(px * px + py * py + pz * pz)
    omega1 = math.sqrt(abs(energy) + pp)
    omega2 = mass / omega1
    omega = (omega1, omega2)
    sf1 = (1 + nsf + (1 - nsf) * nh) * 0.5
    sf2 = (1 + nsf - (1 - nsf) * nh) * 0.5
    ip = (3 + nh) // 2 - 1
    im = (3 - nh) // 2 - 1
    sfomeg = (sf1 * omega[ip], sf2 * omega[im])
    signed_px, signed_py, signed_pz = (
        (px, py, pz) if energy > 0.0 else (-px, -py, -pz)
    )
    pp3 = max(pp + signed_pz, 0.0)
    chi1 = complex(1.0 if pp == 0.0 else math.sqrt(pp3 * 0.5 / pp), 0.0)
    if pp3 == 0.0 or pp == 0.0:
        chi2 = complex(-nh, 0.0)
    else:
        denom = math.sqrt(2.0 * pp * pp3)
        chi2 = complex(nh * signed_px / denom, signed_py / denom)
    chi = (chi1, chi2)
    return (
        chi[im] * sfomeg[0],
        chi[ip] * sfomeg[0],
        chi[im] * sfomeg[1],
        chi[ip] * sfomeg[1],
    )


def _drop_currents_from_dag(
    dag: GenericDAG,
    current_ids: Sequence[int],
) -> GenericDAG:
    removed = set(int(current_id) for current_id in current_ids)
    source_ids = set(dag.sources)
    removed -= source_ids
    kept_ids = {current.id for current in dag.currents if current.id not in removed}
    current_id_map = {
        old_id: new_id
        for new_id, old_id in enumerate(sorted(kept_ids))
    }
    currents = tuple(
        CurrentNode(
            id=current_id_map[current.id],
            index=current.index,
            dimension=current.dimension,
            is_source=current.is_source,
            source_leg_label=current.source_leg_label,
            source_helicity=current.source_helicity,
        )
        for current in dag.currents
        if current.id in kept_ids
    )
    interactions: list[InteractionNode] = []
    for interaction in dag.interactions:
        if (
            interaction.left_id not in current_id_map
            or interaction.right_id not in current_id_map
            or interaction.result_id not in current_id_map
        ):
            continue
        interactions.append(
            InteractionNode(
                id=len(interactions),
                vertex_kind=interaction.vertex_kind,
                vertex_particles=interaction.vertex_particles,
                left_id=current_id_map[interaction.left_id],
                right_id=current_id_map[interaction.right_id],
                result_id=current_id_map[interaction.result_id],
                coupling=interaction.coupling,
                color_weight=interaction.color_weight,
                lowering_backend=interaction.lowering_backend,
                full_tensor_network_ready=interaction.full_tensor_network_ready,
            )
        )
    roots: list[AmplitudeRoot] = []
    for root in dag.amplitude_roots:
        if root.left_id not in current_id_map or root.right_id not in current_id_map:
            continue
        roots.append(
            AmplitudeRoot(
                id=len(roots),
                kind=root.kind,
                left_id=current_id_map[root.left_id],
                right_id=current_id_map[root.right_id],
                color_weight=root.color_weight,
                color_sector_id=root.color_sector_id,
                vertex_kind=root.vertex_kind,
                vertex_particles=root.vertex_particles,
                coupling=root.coupling,
                contraction=root.contraction,
                helicity_weight=root.helicity_weight,
            )
        )
    sources = tuple(
        current_id_map[source_id]
        for source_id in dag.sources
        if source_id in current_id_map
    )
    return prune_dag_to_amplitude_roots(
        GenericDAG(
            process=dag.process,
            color_plan=dag.color_plan,
            currents=currents,
            sources=sources,
            interactions=tuple(interactions),
            amplitude_roots=tuple(roots),
            truncated=dag.truncated,
        )
    )


def _equivalent_current_representatives(
    dag: GenericDAG,
    *,
    maxima: Mapping[int, float],
    signatures: Mapping[int, tuple[tuple[complex, ...], ...]],
    absolute_tolerance: float,
    relative_tolerance: float,
) -> tuple[
    dict[int, tuple[int, float]],
    tuple[tuple[int, tuple[tuple[int, float], ...]], ...],
]:
    source_ids = set(dag.sources)
    representative_by_current = {
        current.id: (current.id, 1.0) for current in dag.currents
    }
    representatives_by_key: dict[tuple[object, ...], list[int]] = {}
    group_members: dict[int, list[tuple[int, float]]] = {}
    for current in dag.currents:
        if current.id in source_ids:
            continue
        if maxima.get(current.id, 0.0) <= absolute_tolerance:
            continue
        signature = signatures.get(current.id, ())
        if not signature:
            continue
        key = _current_merge_key(current)
        matched: int | None = None
        relation_sign = 1.0
        for candidate_id in representatives_by_key.get(key, ()):
            relation = _current_signature_relation(
                signature,
                signatures.get(candidate_id, ()),
                absolute_tolerance=absolute_tolerance,
                relative_tolerance=relative_tolerance,
            )
            if relation is not None:
                matched = candidate_id
                relation_sign = relation
                break
        if matched is None:
            representatives_by_key.setdefault(key, []).append(current.id)
            continue
        representative_by_current[current.id] = (matched, relation_sign)
        group_members.setdefault(matched, []).append((current.id, relation_sign))
    groups = tuple(
        (representative_id, tuple(member_ids))
        for representative_id, member_ids in sorted(group_members.items())
    )
    return representative_by_current, groups


def _current_merge_key(current: CurrentNode) -> tuple[object, ...]:
    index = current.index
    return (
        index.particle_id,
        current.dimension,
        index.helicity_ancestry,
        index.chirality,
        index.spin_state,
        index.flavour_flow,
        index.charge_flow,
        index.color_state,
        index.momentum_mask,
        index.coupling_orders,
        index.auxiliary_kind,
    )


def _current_signature_relation(
    left: Sequence[Sequence[complex]],
    right: Sequence[Sequence[complex]],
    *,
    absolute_tolerance: float,
    relative_tolerance: float,
) -> float | None:
    if len(left) != len(right):
        return None
    relation: float | None = None
    for left_sample, right_sample in zip(left, right, strict=True):
        if len(left_sample) != len(right_sample):
            return None
        for left_value, right_value in zip(left_sample, right_sample, strict=True):
            scale = max(abs(left_value), abs(right_value))
            tolerance = max(absolute_tolerance, relative_tolerance * scale)
            if scale <= tolerance:
                continue
            if abs(left_value - right_value) <= tolerance:
                component_relation = 1.0
            elif abs(left_value + right_value) <= tolerance:
                component_relation = -1.0
            else:
                return None
            if relation is None:
                relation = component_relation
            elif relation != component_relation:
                return None
    return 1.0 if relation is None else relation


def _merge_currents_in_dag(
    dag: GenericDAG,
    representative_by_current: Mapping[int, tuple[int, float]],
) -> GenericDAG:
    representative = {
        current.id: int(representative_by_current.get(current.id, (current.id, 1.0))[0])
        for current in dag.currents
    }
    current_sign = {
        current.id: float(
            representative_by_current.get(current.id, (current.id, 1.0))[1]
        )
        for current in dag.currents
    }
    kept_old_ids = {
        current_id for current_id, representative_id in representative.items()
        if current_id == representative_id
    }
    current_id_map = {
        old_id: new_id
        for new_id, old_id in enumerate(sorted(kept_old_ids))
    }
    currents = tuple(
        CurrentNode(
            id=current_id_map[current.id],
            index=current.index,
            dimension=current.dimension,
            is_source=current.is_source,
            source_leg_label=current.source_leg_label,
            source_helicity=current.source_helicity,
        )
        for current in dag.currents
        if current.id in kept_old_ids
    )
    interactions: list[InteractionNode] = []
    for interaction in dag.interactions:
        result_representative = representative[interaction.result_id]
        if result_representative != interaction.result_id:
            continue
        left_id = representative[interaction.left_id]
        right_id = representative[interaction.right_id]
        result_id = result_representative
        if (
            left_id not in current_id_map
            or right_id not in current_id_map
            or result_id not in current_id_map
        ):
            continue
        input_sign = current_sign[interaction.left_id] * current_sign[interaction.right_id]
        interactions.append(
            InteractionNode(
                id=len(interactions),
                vertex_kind=interaction.vertex_kind,
                vertex_particles=interaction.vertex_particles,
                left_id=current_id_map[left_id],
                right_id=current_id_map[right_id],
                result_id=current_id_map[result_id],
                coupling=_scale_complex_tuple(interaction.coupling, input_sign),
                color_weight=interaction.color_weight,
                lowering_backend=interaction.lowering_backend,
                full_tensor_network_ready=interaction.full_tensor_network_ready,
            )
        )
    roots: list[AmplitudeRoot] = []
    for root in dag.amplitude_roots:
        left_id = representative[root.left_id]
        right_id = representative[root.right_id]
        if left_id not in current_id_map or right_id not in current_id_map:
            continue
        input_sign = current_sign[root.left_id] * current_sign[root.right_id]
        roots.append(
            AmplitudeRoot(
                id=len(roots),
                kind=root.kind,
                left_id=current_id_map[left_id],
                right_id=current_id_map[right_id],
                color_weight=_scale_complex_tuple(root.color_weight, input_sign),
                color_sector_id=root.color_sector_id,
                vertex_kind=root.vertex_kind,
                vertex_particles=root.vertex_particles,
                coupling=root.coupling,
                contraction=root.contraction,
                helicity_weight=root.helicity_weight,
            )
        )
    sources = tuple(
        current_id_map[source_id]
        for source_id in dag.sources
        if source_id in current_id_map
    )
    return prune_dag_to_amplitude_roots(
        GenericDAG(
            process=dag.process,
            color_plan=dag.color_plan,
            currents=currents,
            sources=sources,
            interactions=tuple(interactions),
            amplitude_roots=tuple(roots),
            truncated=dag.truncated,
        )
    )


def _scale_complex_tuple(value: tuple[float, float], factor: float) -> tuple[float, float]:
    return (float(value[0]) * float(factor), float(value[1]) * float(factor))


def _schema_value_slots_by_id(schema: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    storage = cast(dict[str, Any], schema["value_storage"])
    return {
        int(slot["value_slot_id"]): cast(dict[str, Any], slot)
        for slot in cast(list[Any], storage["value_slots"])
    }


def _schema_current_slots_by_id(schema: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    storage = cast(dict[str, Any], schema["current_storage"])
    return {
        int(slot["current_id"]): cast(dict[str, Any], slot)
        for slot in cast(list[Any], storage["current_slots"])
    }


def _schema_momentum_slots_by_id(schema: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    return {
        int(slot["momentum_slot_id"]): cast(dict[str, Any], slot)
        for slot in cast(list[Any], schema["momentum_slots"])
    }


def _is_fermion_pdg(particle_id: int) -> bool:
    abs_id = abs(int(particle_id))
    return 1 <= abs_id <= 6 or 11 <= abs_id <= 16


def _fortran_sign(value: float, sign_source: float) -> float:
    return abs(value) if sign_source >= 0.0 else -abs(value)


def _model_mass(model: Model, particle_id: int) -> float:
    try:
        return float(model.mass(particle_id))
    except KeyError:
        return float(model.mass(-particle_id))


def _write_generic_validation_momenta(
    manifest: GenericProcessManifest,
    output_path: Path,
) -> None:
    try:
        point = _generic_validation_point(manifest)
    except NativeEvaluationError as error:
        payload = {
            "schema_version": 1,
            "kind": "pyamplicol-rusticol-validation-momenta",
            "process": manifest.process,
            "available": False,
            "error": str(error),
            "points": [],
        }
    else:
        payload = {
            "schema_version": 1,
            "kind": "pyamplicol-rusticol-validation-momenta",
            "process": manifest.process,
            "available": True,
            "points": [
                [
                    {
                        "pdg": int(particle.pdg),
                        "momentum": [
                            _decimal_string(component)
                            for component in particle.momentum
                        ],
                    }
                    for particle in point
                ]
            ],
        }
    (output_path / "validation_momenta.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _generic_validation_point(
    manifest: GenericProcessManifest,
) -> tuple[ExternalMomentum, ...]:
    initial_pdgs = tuple(int(pdg) for pdg in manifest.dag.process.initial_pdgs)
    final_pdgs = tuple(int(pdg) for pdg in manifest.dag.process.final_pdgs)
    if len(initial_pdgs) != 2:
        raise NativeEvaluationError(
            "generic validation momenta currently require a two-body initial state"
        )
    if not final_pdgs:
        raise NativeEvaluationError(
            "generic validation momenta require at least one final-state particle"
        )
    final_masses = tuple(float(manifest.model.mass(pdg)) for pdg in final_pdgs)
    threshold = sum(final_masses)
    final_momenta: tuple[tuple[float, float, float, float], ...]
    if len(final_pdgs) == 1:
        if threshold <= 0.0:
            raise NativeEvaluationError(
                "no finite centre-of-mass validation point exists for one massless final state"
            )
        sqrt_s = threshold
        final_momenta = ((float(sqrt_s), 0.0, 0.0, 0.0),)
    else:
        sqrt_s = max(1000.0, threshold + 100.0)
        final_momenta = massive_rambo_final_state(
            len(final_pdgs),
            sqrt_s=sqrt_s,
            masses=final_masses,
            seed=101,
        )
    beam_energy = 0.5 * sqrt_s
    return (
        ExternalMomentum(initial_pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
        ExternalMomentum(initial_pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
        *(
            ExternalMomentum(pdg, momentum)
            for pdg, momentum in zip(final_pdgs, final_momenta, strict=True)
        ),
    )


def _decimal_string(value: float) -> str:
    return format(float(value), ".17g")


def _schema_int(value: object) -> int:
    return int(cast(Any, value))


def _schema_list(value: object) -> list[Any]:
    if not isinstance(value, list):
        raise TypeError("expected JSON array")
    return value


def _runtime_lc_sector_artifact_key(sector_ids: set[int]) -> str:
    return "lc_" + "_".join(str(sector_id) for sector_id in sorted(sector_ids))


def _generic_manifest_status_payloads(
    manifest: GenericProcessManifest,
    *,
    include_color_sector_summaries: bool = True,
) -> tuple[dict[str, object], dict[str, object]]:
    lowering = _dag_lowering_status(
        manifest.dag,
        manifest.model,
        include_color_sector_summaries=include_color_sector_summaries,
    )
    current_ready = (
        manifest.dag.has_amplitudes
        and not lowering["pending_vertex_kinds"]
        and not lowering["unimplemented_vertex_kinds"]
        and not lowering["pending_propagator_kernels"]
        and not lowering["unimplemented_propagator_kernels"]
        and not manifest.dag.truncated
    )
    color_ready = manifest.color_plan.ready_for_requested_colour
    planning_status = {
        "color_ready": color_ready,
        "color_sector_count": manifest.color_plan.sector_count,
        "color_truncated": manifest.color_plan.truncated,
        "idenso_required": manifest.color_plan.idenso_required,
        "current_ready": current_ready,
        "has_closure": manifest.dag.has_amplitudes,
        "has_amplitude_roots": manifest.dag.has_amplitudes,
        "generic_evaluator_ready": color_ready and current_ready,
    }
    lowering_status = {
        **lowering,
        "closure_count": len(manifest.dag.amplitude_roots),
        "amplitude_root_count": len(manifest.dag.amplitude_roots),
        "truncated": manifest.dag.truncated,
        "has_closure": manifest.dag.has_amplitudes,
        "has_amplitude_roots": manifest.dag.has_amplitudes,
        "full_tensor_network_ready": current_ready,
    }
    return planning_status, lowering_status


def _generation_filters_payload(
    manifest: GenericProcessManifest,
) -> dict[str, object]:
    return {
        "structural_current_aggregation": dict(
            manifest.structural_current_aggregation
            or {
                "enabled": False,
                "reason": "structural current aggregation was not requested",
            }
        ),
        "zero_current": dict(
            manifest.zero_current_filter
            or {
                "enabled": False,
                "reason": "zero-current warmup filter was not requested",
            }
        ),
        "current_merging": dict(
            manifest.current_merging
            or {
                "enabled": False,
                "reason": "numerical current merging was not requested",
            }
        ),
    }


def _generic_pruning_payload(
    manifest: GenericProcessManifest,
    *,
    max_coupling_orders: Mapping[str, int] | None,
    max_lc_current_line_groups: int | None,
    max_quark_pairs: int | None,
    closure_side_mask_pruning: bool,
    color_order_mask_pruning: bool,
    species_reachability_pruning: bool,
    ignored_particle_ids: Sequence[int] | None,
    ignored_vertex_kinds: Sequence[int] | None,
    numerical_filter_current: bool,
    numerical_current_merging: bool,
    numerical_current_samples: int,
    numerical_current_seed: int,
    numerical_current_relative_tolerance: float,
    numerical_current_zero_tolerance: float,
    reference_color_order: Sequence[int] | None,
) -> dict[str, object]:
    return {
        "max_coupling_orders": dict(max_coupling_orders or {}),
        "max_lc_current_line_groups": max_lc_current_line_groups,
        "max_quark_pairs": max_quark_pairs,
        "closure_side_mask_pruning": closure_side_mask_pruning,
        "color_order_mask_pruning": color_order_mask_pruning,
        "species_reachability_pruning": species_reachability_pruning,
        "ignored_particle_ids": list(ignored_particle_ids or ()),
        "ignored_vertex_kinds": list(ignored_vertex_kinds or ()),
        "structural_current_aggregation": dict(
            manifest.structural_current_aggregation
            or {
                "enabled": manifest.dag.process.color_accuracy == "lc",
                "mode": "lc-gluon-flavour-flow-aggregation",
                "reason": "structural current aggregation was not run",
            }
        ),
        "zero_current_filter": dict(
            manifest.zero_current_filter
            or {
                "enabled": bool(numerical_filter_current),
                "sample_count": int(numerical_current_samples),
                "seed": int(numerical_current_seed),
                "relative_tolerance": float(numerical_current_relative_tolerance),
                "zero_tolerance": float(numerical_current_zero_tolerance),
            }
        ),
        "current_merging": dict(
            manifest.current_merging
            or {
                "enabled": bool(numerical_current_merging),
                "sample_count": int(numerical_current_samples),
                "seed": int(numerical_current_seed),
                "relative_tolerance": float(numerical_current_relative_tolerance),
                "zero_tolerance": float(numerical_current_zero_tolerance),
            }
        ),
        "reference_color_order": (
            None
            if reference_color_order is None
            else [int(label) for label in reference_color_order]
        ),
    }


def _stage_evaluator_progress_callback(
    progress_callback: Any | None,
    *,
    completed_steps: int,
    trailing_steps: int,
):
    if progress_callback is None:
        return None

    def callback(event: dict[str, object]) -> None:
        forwarded = dict(event)
        total = forwarded.get("total")
        if isinstance(total, int):
            forwarded["total"] = max(
                int(total) + int(completed_steps) + int(trailing_steps),
                1,
            )
        progress_callback(forwarded)

    return callback


def _generic_dag_process_compact_artifact_payload(
    manifest: GenericProcessManifest,
    *,
    evaluator_backend: str,
    compiled_preset: str,
    batch_size: int,
    stage_local_parameter_layout: bool = False,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
) -> dict[str, object]:
    planning_status, lowering_status = _generic_manifest_status_payloads(
        manifest,
        include_color_sector_summaries=False,
    )
    topology_replay_payload = _lc_topology_replay_payload(
        manifest,
        materialized_sector_ids=selected_color_sector_ids,
        enabled=False,
    )
    compiled_payload: dict[str, object] = {
        "kind": "generic-dag-stage-blueprint",
        "runtime_available": False,
        "runtime_unavailable_message": (
            "main all-sector runtime evaluator was skipped; use a matching "
            "runtime_lc_sector_artifacts sidecar or regenerate without "
            "--skip-main-runtime-evaluator"
        ),
        "compact_main_artifact": True,
        "requested_evaluator_backend": evaluator_backend,
        "requested_compiled_preset": compiled_preset,
        "batch_size": batch_size,
        "stage_local_parameter_layout": stage_local_parameter_layout,
        "generic_pruning": _generic_pruning_payload(
            manifest,
            max_coupling_orders=max_coupling_orders,
            max_lc_current_line_groups=max_lc_current_line_groups,
            max_quark_pairs=max_quark_pairs,
            closure_side_mask_pruning=closure_side_mask_pruning,
            color_order_mask_pruning=color_order_mask_pruning,
            species_reachability_pruning=species_reachability_pruning,
            ignored_particle_ids=ignored_particle_ids,
            ignored_vertex_kinds=ignored_vertex_kinds,
            numerical_filter_current=numerical_filter_current,
            numerical_current_merging=numerical_current_merging,
            numerical_current_samples=numerical_current_samples,
            numerical_current_seed=numerical_current_seed,
            numerical_current_relative_tolerance=numerical_current_relative_tolerance,
            numerical_current_zero_tolerance=numerical_current_zero_tolerance,
            reference_color_order=reference_color_order,
        ),
        "lc_topology_replay": topology_replay_payload,
    }
    if selected_color_sector_ids is not None:
        compiled_payload["selected_color_sector_ids"] = sorted(
            int(sector_id) for sector_id in selected_color_sector_ids
        )
    return {
        "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
        "kind": GENERIC_DAG_PROCESS_ARTIFACT_KIND,
        "artifact_class": "generic-dag-schema-v2",
        "process": manifest.process,
        "key": manifest.key,
        "color_accuracy": manifest.dag.process.color_accuracy,
        "model": {
            "name": manifest.model.name,
            "vertex_lowering_coverage": (
                manifest.model.vertex_lowering_coverage().to_json_dict()
            ),
        },
        "external_pdg_order": list(manifest.external_pdg_order),
        "outgoing_pdg_order": list(manifest.outgoing_pdg_order),
        "process_ir": manifest.dag.process.to_json_dict(),
        "generic_plan_path": None,
        "generation_filters": _generation_filters_payload(manifest),
        "lc_topology_reuse": _lc_topology_reuse_payload(manifest),
        "runtime_lc_topology_reuse": _lc_topology_reuse_payload(manifest),
        "runtime_lc_topology_replay": topology_replay_payload,
        "planning_status": planning_status,
        "lowering_status": lowering_status,
        "compiled": compiled_payload,
        "dag_summary": _dag_count_payload(manifest.dag)
        | {"required_vertex_kinds": list(manifest.dag.required_vertex_kinds)},
        "full_dag_summary": _dag_count_payload(manifest.dag)
        | {"required_vertex_kinds": list(manifest.dag.required_vertex_kinds)},
    }


def _generic_dag_process_artifact_payload(
    manifest: GenericProcessManifest,
    *,
    plan_path: str | None,
    evaluator_backend: str,
    compiled_preset: str,
    batch_size: int,
    stage_evaluator_artifact_dir: Path | None = None,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    stage_local_parameter_layout: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    enable_lc_sector_runtime_selector: bool = True,
    lc_topology_replay: bool = False,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
    numerical_filter_current: bool = True,
    numerical_current_merging: bool = True,
    numerical_current_samples: int = 10,
    numerical_current_seed: int = 12345,
    numerical_current_relative_tolerance: float = 1.0e-12,
    numerical_current_zero_tolerance: float = 1.0e-300,
) -> dict[str, object]:
    from .generic_stage_compiler import (
        build_generic_stage_compiler_blueprint,
        write_generic_stage_evaluator_artifacts,
    )

    _emit_generation_progress(
        progress_callback,
        "main schema",
        "status",
    )
    full_planning_status, _full_lowering_status = _generic_manifest_status_payloads(
        manifest,
        include_color_sector_summaries=False,
    )
    _emit_generation_progress(
        progress_callback,
        "main schema",
        "status ready",
        increment=1,
    )
    runtime_sector_ids = _runtime_color_sector_ids(
        manifest.dag,
        selected_color_sector_ids=selected_color_sector_ids,
    )
    runtime_manifest = _runtime_manifest_for_color_sectors(
        manifest,
        runtime_sector_ids,
    )
    _emit_generation_progress(
        progress_callback,
        "runtime schema",
        (
            "all sectors"
            if runtime_sector_ids is None
            else f"{len(runtime_sector_ids)} sectors"
        ),
    )
    runtime_schema = _generic_runtime_schema_payload(
        runtime_manifest.dag,
        runtime_manifest.model,
        enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
    )
    _runtime_planning_status, runtime_lowering_status = _generic_manifest_status_payloads(
        runtime_manifest,
        include_color_sector_summaries=False,
    )
    _emit_generation_progress(
        progress_callback,
        "runtime schema",
        "ready",
        increment=1,
    )
    color_contraction_message = _runtime_color_contraction_unavailable_message(
        runtime_schema,
    )
    _emit_generation_progress(
        progress_callback,
        "stage blueprint",
        "build",
    )
    stage_blueprint = build_generic_stage_compiler_blueprint(
        runtime_manifest,
        stage_local_parameter_layout=stage_local_parameter_layout,
        enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
    )
    _emit_generation_progress(
        progress_callback,
        "stage blueprint",
        f"stages={stage_blueprint.stage_count}",
        increment=1,
    )
    compiled_payload: dict[str, object] = {
        "kind": "generic-dag-stage-blueprint",
        "runtime_available": False,
        "runtime_unavailable_message": (
            "generic DAG evaluator stages are symbolically lowered, but "
            "serialized evaluator emission and Rusticol schema-v2 execution "
            "are not implemented yet"
        ),
        "requested_evaluator_backend": evaluator_backend,
        "requested_compiled_preset": compiled_preset,
        "batch_size": batch_size,
        "stage_local_parameter_layout": stage_local_parameter_layout,
        "stage_compiler": stage_blueprint.to_json_dict(),
        "generic_pruning": {
            "max_coupling_orders": dict(max_coupling_orders or {}),
            "max_lc_current_line_groups": max_lc_current_line_groups,
            "max_quark_pairs": max_quark_pairs,
            "closure_side_mask_pruning": closure_side_mask_pruning,
            "color_order_mask_pruning": color_order_mask_pruning,
            "species_reachability_pruning": species_reachability_pruning,
            "ignored_particle_ids": list(ignored_particle_ids or ()),
            "ignored_vertex_kinds": list(ignored_vertex_kinds or ()),
            "structural_current_aggregation": dict(
                manifest.structural_current_aggregation
                or {
                    "enabled": manifest.dag.process.color_accuracy == "lc",
                    "mode": "lc-gluon-flavour-flow-aggregation",
                    "reason": "structural current aggregation was not run",
                }
            ),
            "zero_current_filter": dict(
                manifest.zero_current_filter
                or {
                    "enabled": bool(numerical_filter_current),
                    "sample_count": int(numerical_current_samples),
                    "seed": int(numerical_current_seed),
                    "relative_tolerance": float(numerical_current_relative_tolerance),
                    "zero_tolerance": float(numerical_current_zero_tolerance),
                }
            ),
            "current_merging": dict(
                manifest.current_merging
                or {
                    "enabled": bool(numerical_current_merging),
                    "sample_count": int(numerical_current_samples),
                    "seed": int(numerical_current_seed),
                    "relative_tolerance": float(numerical_current_relative_tolerance),
                    "zero_tolerance": float(numerical_current_zero_tolerance),
                }
            ),
            "reference_color_order": (
                None
                if reference_color_order is None
                else [int(label) for label in reference_color_order]
            ),
        },
    }
    topology_replay_payload = _lc_topology_replay_payload(
        manifest,
        materialized_sector_ids=runtime_sector_ids,
        enabled=lc_topology_replay,
    )
    compiled_payload["lc_topology_replay"] = topology_replay_payload
    if runtime_sector_ids is not None:
        compiled_payload["selected_color_sector_ids"] = sorted(runtime_sector_ids)
        compiled_payload["runtime_pruning"] = {
            "policy": "contributing-lc-colour-sectors",
            "full_current_count": len(manifest.dag.currents),
            "full_source_count": len(manifest.dag.sources),
            "full_interaction_count": len(manifest.dag.interactions),
            "full_amplitude_root_count": len(manifest.dag.amplitude_roots),
        }
    if stage_evaluator_artifact_dir is not None:
        stage_evaluator_payload = write_generic_stage_evaluator_artifacts(
            stage_blueprint,
            stage_evaluator_artifact_dir,
            compiler=stage_evaluator_compiler,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=_stage_evaluator_progress_callback(
                progress_callback,
                completed_steps=3,
                trailing_steps=2,
            ),
        )
        compiled_payload["stage_evaluators"] = stage_evaluator_payload
        if bool(stage_evaluator_payload.get("runtime_available", False)):
            compiled_payload["runtime_available"] = True
            compiled_payload["runtime_unavailable_message"] = None
    if color_contraction_message is not None:
        compiled_payload["runtime_available"] = False
        compiled_payload["runtime_unavailable_message"] = color_contraction_message
    return {
        "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
        "kind": GENERIC_DAG_PROCESS_ARTIFACT_KIND,
        "artifact_class": "generic-dag-schema-v2",
        "process": manifest.process,
        "key": manifest.key,
        "color_accuracy": manifest.dag.process.color_accuracy,
        "model": {
            "name": manifest.model.name,
            "vertex_lowering_coverage": (
                manifest.model.vertex_lowering_coverage().to_json_dict()
            ),
        },
        "external_pdg_order": list(manifest.external_pdg_order),
        "outgoing_pdg_order": list(manifest.outgoing_pdg_order),
        "process_ir": manifest.dag.process.to_json_dict(),
        "generic_plan_path": plan_path,
        "generation_filters": _generation_filters_payload(runtime_manifest),
        "lc_topology_reuse": _lc_topology_reuse_payload(manifest),
        "runtime_lc_topology_reuse": _lc_topology_reuse_payload(runtime_manifest),
        "runtime_lc_topology_replay": topology_replay_payload,
        "planning_status": full_planning_status,
        "lowering_status": runtime_lowering_status,
        "runtime_schema": runtime_schema,
        "normalization": runtime_schema["normalization"],
        "compiled": compiled_payload,
        "dag_summary": {
            "current_count": len(runtime_manifest.dag.currents),
            "source_count": len(runtime_manifest.dag.sources),
            "interaction_count": len(runtime_manifest.dag.interactions),
            "amplitude_root_count": len(
                runtime_manifest.dag.amplitude_roots
            ),
            "required_vertex_kinds": list(runtime_manifest.dag.required_vertex_kinds),
            "truncated": runtime_manifest.dag.truncated,
        },
        "full_dag_summary": {
            "current_count": len(manifest.dag.currents),
            "source_count": len(manifest.dag.sources),
            "interaction_count": len(manifest.dag.interactions),
            "amplitude_root_count": len(manifest.dag.amplitude_roots),
            "required_vertex_kinds": list(manifest.dag.required_vertex_kinds),
            "truncated": manifest.dag.truncated,
        },
    }


def _runtime_color_contraction_unavailable_message(
    runtime_schema: Mapping[str, Any],
) -> str | None:
    amplitude_stage = runtime_schema.get("amplitude_stage")
    if not isinstance(amplitude_stage, Mapping):
        return "generic runtime schema is missing amplitude_stage"
    color_contraction = amplitude_stage.get("color_contraction")
    if color_contraction is None:
        return None
    if not isinstance(color_contraction, Mapping):
        return "generic runtime schema has malformed color_contraction"
    if bool(color_contraction.get("supported", False)):
        return None
    reason = color_contraction.get("reason")
    return (
        str(reason)
        if reason
        else "requested colour contraction is not supported by this artifact"
    )


def _runtime_color_sector_ids(
    dag: GenericDAG,
    *,
    selected_color_sector_ids: set[int] | None,
) -> set[int] | None:
    if selected_color_sector_ids is not None:
        return set(selected_color_sector_ids)
    if dag.process.color_accuracy != "lc" or not dag.amplitude_roots:
        return None
    return set(contributing_color_sector_ids(dag))


def _runtime_manifest_for_color_sectors(
    manifest: GenericProcessManifest,
    sector_ids: set[int] | None,
) -> GenericProcessManifest:
    if sector_ids is None:
        return manifest
    filtered = filter_dag_to_color_sectors(manifest.dag, sector_ids)
    return GenericProcessManifest(
        dag=filtered,
        model=manifest.model,
        color_plan=manifest.color_plan,
        structural_current_aggregation=manifest.structural_current_aggregation,
        zero_current_filter=manifest.zero_current_filter,
        current_merging=manifest.current_merging,
    )


def _generic_runtime_schema_payload(
    dag: GenericDAG,
    model: Model,
    *,
    selected_color_sector_ids: set[int] | None = None,
    enable_lc_sector_runtime_selector: bool = True,
) -> dict[str, object]:
    current_slots = _runtime_current_slots(dag)
    slot_by_current_id = {
        _schema_int(slot["current_id"]): slot for slot in current_slots
    }
    current_usage = _runtime_current_usage(
        dag,
        selected_color_sector_ids=selected_color_sector_ids,
    )
    value_slots = _runtime_value_slots(
        dag,
        model=model,
        current_slots=current_slots,
        current_usage=current_usage,
    )
    value_slot_by_current_variant = {
        (_schema_int(slot["current_id"]), str(slot["variant"])): slot
        for slot in value_slots
    }
    momentum_slots = _runtime_momentum_slots(dag)
    momentum_slot_by_mask = {
        _schema_int(slot["momentum_mask"]): _schema_int(slot["momentum_slot_id"])
        for slot in momentum_slots
    }
    stage_payloads = _runtime_stage_payloads(
        dag,
        model,
        slot_by_current_id=slot_by_current_id,
        value_slot_by_current_variant=value_slot_by_current_variant,
        momentum_slot_by_mask=momentum_slot_by_mask,
    )
    amplitude_stage_payload = _runtime_amplitude_stage_payload(
        dag,
        slot_by_current_id=slot_by_current_id,
        value_slot_by_current_variant=value_slot_by_current_variant,
        selected_color_sector_ids=selected_color_sector_ids,
    )
    model_parameters = _runtime_model_parameter_records(
        dag,
        model,
        stage_payloads=stage_payloads,
        amplitude_stage_payload=amplitude_stage_payload,
        enable_lc_sector_runtime_selector=enable_lc_sector_runtime_selector,
    )
    source_component_count = sum(
        _schema_int(slot_by_current_id[source_id]["dimension"])
        for source_id in dag.sources
    )
    model_parameter_count = len(model_parameters)
    return {
        "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
        "kind": "pyamplicol-generic-dag-runtime-schema",
        "process_key": dag.process.key,
        "process": dag.process.process,
        "color_accuracy": dag.process.color_accuracy,
        "external_particles": _runtime_external_particles(dag),
        "momentum_conventions": _runtime_momentum_conventions(dag),
        "model": _runtime_model_payload(model),
        "normalization": _runtime_normalization_payload(dag, model),
        "parameter_layout": {
            "source_component_parameter_count": source_component_count,
            "momentum_parameter_count": 4 * len(momentum_slots),
            "model_parameter_count": model_parameter_count,
            "parameter_count_if_flattened": (
                source_component_count
                + 4 * len(momentum_slots)
                + model_parameter_count
            ),
            "value_component_count": (
                value_slots[-1]["component_stop"] if value_slots else 0
            ),
            "source_components_complex": True,
            "momentum_components_real": True,
            "real_valued_inputs": list(
                range(
                    source_component_count,
                    source_component_count
                    + 4 * len(momentum_slots)
                    + model_parameter_count,
                )
            ),
        },
        "model_parameters": model_parameters,
        "current_storage": {
            "component_count": (
                current_slots[-1]["component_stop"] if current_slots else 0
            ),
            "number_type": "complex",
            "current_slots": current_slots,
        },
        "value_storage": {
            "component_count": (
                value_slots[-1]["component_stop"] if value_slots else 0
            ),
            "number_type": "complex",
            "policy": (
                "current identity is independent of propagator state; runtime "
                "value slots store source, propagated, and unpropagated variants "
                "as required by interaction inputs and amplitude roots"
            ),
            "value_slots": value_slots,
        },
        "source_fill": {
            "source_count": len(dag.sources),
            "sources": [
                _runtime_source_record(
                    dag,
                    current_id=source_id,
                    source_index=source_index,
                    current_slot=slot_by_current_id[source_id],
                    value_slot=value_slot_by_current_variant[(source_id, "source")],
                )
                for source_index, source_id in enumerate(dag.sources)
            ],
        },
        "momentum_slots": momentum_slots,
        "stages": stage_payloads,
        "amplitude_stage": amplitude_stage_payload,
    }


def _runtime_external_particles(dag: GenericDAG) -> list[dict[str, object]]:
    return [
        {
            "label": leg.label,
            "index": leg.label - 1,
            "side": leg.side,
            "role": "initial" if leg.is_initial else "final",
            "particle": leg.particle,
            "outgoing_particle": leg.outgoing_particle,
            "pdg": leg.pdg,
            "outgoing_pdg": leg.outgoing_pdg,
            "particle_class": leg.particle_class,
            "momentum_slot": leg.label - 1,
            "momentum_components": ["E", "px", "py", "pz"],
        }
        for leg in dag.process.legs
    ]


def _runtime_momentum_conventions(dag: GenericDAG) -> dict[str, object]:
    incoming_labels = [leg.label for leg in dag.process.initial_legs]
    final_labels = [leg.label for leg in dag.process.final_legs]
    return {
        "input_shape": ["batch", len(dag.process.legs), 4],
        "component_order": ["E", "px", "py", "pz"],
        "input_momenta": "physical external four-momenta in process order",
        "incoming_labels": incoming_labels,
        "final_state_labels": final_labels,
        "all_outgoing_convention": {
            "crossed_incoming_labels": incoming_labels,
            "operation": "negate incoming four-vectors before current/source use",
        },
        "metric": "mostly-minus",
    }


def _runtime_model_payload(model: Model) -> dict[str, object]:
    return {
        "name": model.name,
        "parameters": {
            "alpha_s_me_check": getattr(model, "alpha_s_me_check", None),
            "alpha_ew": getattr(model, "alpha_ew", None),
            "sqrt_s_default": getattr(model, "sqrt_s", None),
        },
        "particles": [
            {
                "pdg": particle.pdg,
                "anti_pdg": particle.anti_pdg,
                "spin": particle.spin,
                "dimension": particle.dimension,
                "color_rep": particle.color_rep,
                "mass": particle.mass,
                "width": particle.width,
                "charge": particle.charge,
                "weak_isospin": list(particle.weak_isospin),
                "weak_hypercharge": list(particle.weak_hypercharge),
            }
            for particle in sorted(
                model.particles.values(),
                key=lambda item: item.pdg,
            )
        ],
        "vertices": [
            {
                "kind": vertex.kind,
                "particles": list(vertex.particles),
                "coupling": list(vertex.coupling),
                "lowering": (
                    model.vertex_lowering_rule(vertex.kind).to_json_dict()
                    if hasattr(model.vertex_lowering_rule(vertex.kind), "to_json_dict")
                    else {
                        "kind": vertex.kind,
                        "backend": model.vertex_lowering_rule(vertex.kind).backend,
                        "kernel": model.vertex_lowering_rule(vertex.kind).kernel,
                        "full_tensor_network_ready": (
                            model.vertex_lowering_rule(
                                vertex.kind
                            ).full_tensor_network_ready
                        ),
                    }
                ),
            }
            for vertex in model.vertices
        ],
    }


def _runtime_model_parameter_records(
    dag: GenericDAG,
    model: Model,
    *,
    stage_payloads: Sequence[Mapping[str, object]],
    amplitude_stage_payload: Mapping[str, object],
    enable_lc_sector_runtime_selector: bool = True,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    seen: set[str] = set()

    def add_record(
        name: str,
        kind: str,
        default: float,
        **metadata: object,
    ) -> None:
        if name in seen:
            return
        seen.add(name)
        records.append(
            {
                "name": name,
                "kind": kind,
                "parameter_index": len(records),
                "default": float(default),
                **metadata,
            }
        )

    alpha_s = getattr(model, "alpha_s_me_check", None)
    if alpha_s is not None:
        add_record(
            "normalization.alpha_s_me_check",
            "normalization",
            float(alpha_s),
        )
    alpha_ew = getattr(model, "alpha_ew", None)
    if alpha_ew is not None:
        add_record(
            "normalization.alpha_ew",
            "normalization",
            float(alpha_ew),
        )
    if (
        enable_lc_sector_runtime_selector
        and dag.process.color_accuracy == "lc"
        and dag.amplitude_roots
    ):
        add_record(
            LC_SECTOR_SELECTOR_PARAMETER,
            "runtime_control",
            -1.0,
            control="lc_sector_id",
            mode="minus_one_means_all",
        )

    for particle in sorted(model.particles.values(), key=lambda item: item.pdg):
        if float(particle.mass) != 0.0:
            add_record(
                _runtime_particle_parameter_name(particle.pdg, "mass"),
                "particle_mass",
                float(particle.mass),
                pdg=int(particle.pdg),
            )
        if float(particle.width) != 0.0:
            add_record(
                _runtime_particle_parameter_name(particle.pdg, "width"),
                "particle_width",
                float(particle.width),
                pdg=int(particle.pdg),
            )

    processed_coupling_name_sets: set[tuple[str | None, ...]] = set()
    for stage in stage_payloads:
        for interaction in _schema_list(stage["interactions"]):
            names = interaction.get("coupling_parameter_names")
            values = interaction.get("coupling")
            if not isinstance(names, list) or not isinstance(values, list):
                continue
            name_key = tuple(
                None if name is None else str(name) for name in names
            )
            if name_key in processed_coupling_name_sets:
                continue
            processed_coupling_name_sets.add(name_key)
            particles = tuple(int(pdg) for pdg in _schema_list(interaction["vertex_particles"]))
            for component, name in enumerate(names):
                if not isinstance(name, str):
                    continue
                add_record(
                    name,
                    "coupling_component",
                    float(values[component]),
                    vertex_kind=int(interaction["vertex_kind"]),
                    vertex_particles=list(particles),
                    component=component,
                )

    for root in _schema_list(amplitude_stage_payload["roots"]):
        names = root.get("coupling_parameter_names")
        values = root.get("coupling")
        particles = root.get("vertex_particles")
        if not isinstance(names, list) or not isinstance(values, list) or not isinstance(particles, list):
            continue
        name_key = tuple(None if name is None else str(name) for name in names)
        if name_key in processed_coupling_name_sets:
            continue
        processed_coupling_name_sets.add(name_key)
        for component, name in enumerate(names):
            if not isinstance(name, str):
                continue
            add_record(
                name,
                "coupling_component",
                float(values[component]),
                vertex_kind=int(root["vertex_kind"]),
                vertex_particles=[int(pdg) for pdg in particles],
                component=component,
            )

    return records


def _runtime_particle_parameter_name(pdg: int, field: str) -> str:
    return f"particle.{int(pdg)}.{field}"


def _runtime_coupling_parameter_names(
    vertex_kind: int,
    vertex_particles: Sequence[int],
    coupling: Sequence[object],
) -> list[str | None]:
    base = _runtime_coupling_parameter_base(vertex_kind, vertex_particles)
    names: list[str | None] = []
    for component, value in enumerate(coupling):
        numeric_value = float(value)
        if component == 1 and numeric_value == -10.0:
            # AmpliCol uses this sentinel as a structural phase flag in the
            # scalar kernel, not as a continuously variable model parameter.
            names.append(None)
        else:
            names.append(f"{base}.component_{component}")
    return names


def _runtime_coupling_parameter_base(
    vertex_kind: int,
    vertex_particles: Sequence[int],
) -> str:
    particles = "_".join(str(int(pdg)) for pdg in vertex_particles)
    return f"coupling.{int(vertex_kind)}.{particles}"


def _runtime_normalization_payload(
    dag: GenericDAG,
    model: Model,
) -> dict[str, object]:
    external_pdgs = [
        int(pdg)
        for pdg in (*dag.process.initial_pdgs, *dag.process.final_pdgs)
    ]
    try:
        leading_color_factor = getattr(model, "leading_color_factor")
        color_factor: int | None = int(leading_color_factor(external_pdgs))
    except Exception:
        color_factor = None
    initial_pdgs = [int(pdg) for pdg in dag.process.initial_pdgs]
    final_pdgs = [int(pdg) for pdg in dag.process.final_pdgs]
    final_state_identical_factor = _final_state_identical_factor(final_pdgs)
    quark_line_partner_factor = 1
    electroweak_power = _electroweak_coupling_power(dag)
    qcd_power = max(0, len(dag.process.legs) - 2 - electroweak_power)
    alpha_s_me_check = float(getattr(model, "alpha_s_me_check", 0.118))
    alpha_ew = float(getattr(model, "alpha_ew", 0.007546771114))
    global_coupling_factor = (
        (4.0 * math.pi * alpha_s_me_check) ** qcd_power
        * (2.0 * 4.0 * math.pi * alpha_ew) ** electroweak_power
    )
    return {
        "color_accuracy": dag.process.color_accuracy,
        "color_factor": color_factor,
        "average_factor": _initial_state_average_factor(initial_pdgs),
        "identical_factor": final_state_identical_factor,
        "final_state_identical_factor": final_state_identical_factor,
        "quark_line_partner_factor": quark_line_partner_factor,
        "global_coupling_factor": global_coupling_factor,
        "qcd_coupling_power": qcd_power,
        "electroweak_coupling_power": electroweak_power,
        "couplings_in_stage_evaluators": True,
        "coupling_policy": (
            "generic vertex couplings are stored on interaction/amplitude "
            "records and should be applied by the generated stage evaluators; "
            "dimensionful alpha_s/alpha_ew powers follow AmpliCol ME-check "
            "normalization"
        ),
    }


def _electroweak_coupling_power(dag: GenericDAG) -> int:
    if not dag.process.singlet_labels:
        return 0
    return max(1, len(dag.process.singlet_labels))


def _runtime_current_slots(dag: GenericDAG) -> list[dict[str, object]]:
    offset = 0
    slots: list[dict[str, object]] = []
    for current in dag.currents:
        start = offset
        stop = start + current.dimension
        offset = stop
        slots.append(
            {
                "current_id": current.id,
                "component_start": start,
                "component_stop": stop,
                "dimension": current.dimension,
                "is_source": current.is_source,
                "particle_id": current.index.particle_id,
                "external_mask": current.index.external_mask,
                "external_labels": list(current.index.external_labels),
                "momentum_mask": current.index.momentum_mask,
                "helicity_ancestry": _bigint_json(current.index.helicity_ancestry),
                "chirality": current.index.chirality,
                "spin_state": _spin_state_json(current.index.spin_state),
                "flavour_flow": list(current.index.flavour_flow),
                "charge_flow": current.index.charge_flow,
                "color_state": current.index.color_state.to_json_dict(),
                "auxiliary_kind": current.index.auxiliary_kind,
            }
        )
    return slots


def _selected_amplitude_roots(
    dag: GenericDAG,
    *,
    selected_color_sector_ids: set[int] | None = None,
):
    if selected_color_sector_ids is None:
        return tuple(dag.amplitude_roots)
    return tuple(
        root
        for root in dag.amplitude_roots
        if _root_color_sector_id(dag, root) in selected_color_sector_ids
    )


def _runtime_current_usage(
    dag: GenericDAG,
    *,
    selected_color_sector_ids: set[int] | None = None,
) -> dict[int, dict[str, bool]]:
    usage = {
        current.id: {
            "used_as_interaction_input": False,
            "used_as_amplitude_input": False,
        }
        for current in dag.currents
    }
    for interaction in dag.interactions:
        usage[interaction.left_id]["used_as_interaction_input"] = True
        usage[interaction.right_id]["used_as_interaction_input"] = True
    for root in _selected_amplitude_roots(
        dag,
        selected_color_sector_ids=selected_color_sector_ids,
    ):
        usage[root.left_id]["used_as_amplitude_input"] = True
        usage[root.right_id]["used_as_amplitude_input"] = True
    return usage


def _runtime_value_slots(
    dag: GenericDAG,
    *,
    model: Model,
    current_slots: list[dict[str, object]],
    current_usage: dict[int, dict[str, bool]],
) -> list[dict[str, object]]:
    offset = 0
    value_slot_id = 0
    value_slots: list[dict[str, object]] = []
    current_slot_by_id = {
        _schema_int(slot["current_id"]): slot for slot in current_slots
    }

    def add_slot(
        current_id: int,
        *,
        variant: str,
        applies_propagator: bool,
    ) -> None:
        nonlocal offset, value_slot_id
        current = dag.currents[current_id]
        current_slot = current_slot_by_id[current_id]
        usage = current_usage[current_id]
        propagator_rule = model.propagator_lowering_rule(
            current.index.particle_id,
            current.index.chirality,
        )
        start = offset
        stop = start + current.dimension
        value_slots.append(
            {
                "value_slot_id": value_slot_id,
                "current_id": current_id,
                "variant": variant,
                "component_start": start,
                "component_stop": stop,
                "dimension": current.dimension,
                "current_component_start": current_slot["component_start"],
                "current_component_stop": current_slot["component_stop"],
                "is_source": current.is_source,
                "applies_propagator": applies_propagator,
                "propagator": propagator_rule.to_json_dict(),
                "used_as_interaction_input": usage["used_as_interaction_input"],
                "used_as_amplitude_input": usage["used_as_amplitude_input"],
                "particle_id": current.index.particle_id,
                "external_mask": current.index.external_mask,
                "external_labels": list(current.index.external_labels),
                "momentum_mask": current.index.momentum_mask,
                "chirality": current.index.chirality,
            }
        )
        value_slot_id += 1
        offset = stop

    for current in dag.currents:
        usage = current_usage[current.id]
        if current.is_source:
            add_slot(current.id, variant="source", applies_propagator=False)
            continue
        propagator_rule = model.propagator_lowering_rule(
            current.index.particle_id,
            current.index.chirality,
        )
        needs_propagated = (
            usage["used_as_interaction_input"]
            and propagator_rule.applies_propagator
        )
        needs_unpropagated = usage["used_as_amplitude_input"]
        if (
            needs_unpropagated
            or not needs_propagated
            or (
                usage["used_as_interaction_input"]
                and not propagator_rule.applies_propagator
            )
        ):
            add_slot(
                current.id,
                variant="unpropagated",
                applies_propagator=False,
            )
        if needs_propagated:
            add_slot(
                current.id,
                variant="propagated",
                applies_propagator=True,
            )
    return value_slots


def _runtime_source_record(
    dag: GenericDAG,
    *,
    current_id: int,
    source_index: int,
    current_slot: dict[str, object],
    value_slot: dict[str, object],
) -> dict[str, object]:
    current = dag.currents[current_id]
    leg = _leg_by_label(dag, int(current.source_leg_label or 0))
    source_start = sum(
        dag.currents[source_id].dimension
        for source_id in dag.sources[:source_index]
    )
    return {
        "source_id": source_index,
        "current_id": current_id,
        "current_component_start": current_slot["component_start"],
        "current_component_stop": current_slot["component_stop"],
        "value_slot": _value_slot_ref(value_slot),
        "source_parameter_start": source_start,
        "source_parameter_stop": source_start + current.dimension,
        "leg_label": current.source_leg_label,
        "input_momentum_slot": None if leg is None else leg.label - 1,
        "side": None if leg is None else leg.side,
        "crossing": (
            "negate-incoming-momentum"
            if leg is not None and leg.is_initial
            else "identity"
        ),
        "physical_pdg": None if leg is None else leg.pdg,
        "outgoing_pdg": current.index.particle_id,
        "particle_id": current.index.particle_id,
        "source_kind": "external-wavefunction",
        "source_helicity": current.source_helicity,
        "chirality": current.index.chirality,
        "spin_state": _spin_state_json(current.index.spin_state),
        "dimension": current.dimension,
        "helicity_ancestry": _bigint_json(current.index.helicity_ancestry),
        "color_state": current.index.color_state.to_json_dict(),
    }


def _runtime_momentum_slots(dag: GenericDAG) -> list[dict[str, object]]:
    masks = sorted(
        {
            current.index.momentum_mask
            for current in dag.currents
        },
        key=lambda mask: (mask.bit_count(), mask),
    )
    initial_labels = {leg.label for leg in dag.process.initial_legs}
    slots = []
    for slot_id, mask in enumerate(masks):
        labels = _mask_labels(mask)
        offset = 4 * slot_id
        slots.append(
            {
                "momentum_slot_id": slot_id,
                "momentum_mask": mask,
                "external_labels": list(labels),
                "component_start": offset,
                "component_stop": offset + 4,
                "component_order": ["E", "px", "py", "pz"],
                "real_valued": True,
                "crossed_incoming_labels": [
                    label for label in labels if label in initial_labels
                ],
                "construction": (
                    "sum all-outgoing momenta for external_labels; incoming "
                    "physical momenta are negated first"
                ),
            }
        )
    return slots


def _runtime_stage_payloads(
    dag: GenericDAG,
    model: Model,
    *,
    slot_by_current_id: dict[int, dict[str, object]],
    value_slot_by_current_variant: dict[tuple[int, str], dict[str, object]],
    momentum_slot_by_mask: dict[int, int],
) -> list[dict[str, object]]:
    interactions_by_size: dict[int, list[int]] = {}
    for interaction in dag.interactions:
        result = dag.currents[interaction.result_id]
        size = len(result.index.external_labels)
        interactions_by_size.setdefault(size, []).append(interaction.id)
    stages = []
    coupling_name_cache: dict[tuple[int, tuple[int, ...], tuple[float, ...]], list[str | None]] = {}
    lowering_payload_cache: dict[int, dict[str, object]] = {}
    for stage_index, size in enumerate(sorted(interactions_by_size), start=1):
        interaction_ids = interactions_by_size[size]
        interactions: list[dict[str, object]] = []
        output_current_ids: set[int] = set()
        input_current_ids: set[int] = set()
        input_value_slot_ids: set[int] = set()
        output_value_slot_ids: set[int] = set()
        for interaction_id in interaction_ids:
            interaction = _runtime_interaction_record(
                dag,
                model,
                interaction_id=interaction_id,
                slot_by_current_id=slot_by_current_id,
                value_slot_by_current_variant=value_slot_by_current_variant,
                momentum_slot_by_mask=momentum_slot_by_mask,
                coupling_name_cache=coupling_name_cache,
                lowering_payload_cache=lowering_payload_cache,
            )
            interaction_any = cast(dict[str, Any], interaction)
            interactions.append(interaction)
            output_current_ids.add(_schema_int(interaction["result_current_id"]))
            input_current_ids.add(_schema_int(interaction["left_current_id"]))
            input_current_ids.add(_schema_int(interaction["right_current_id"]))
            input_value_slot_ids.add(
                _schema_int(interaction_any["left_value_slot"]["value_slot_id"])
            )
            input_value_slot_ids.add(
                _schema_int(interaction_any["right_value_slot"]["value_slot_id"])
            )
            for slot in interaction_any["result_value_slots"]:
                output_value_slot_ids.add(_schema_int(slot["value_slot_id"]))
        stages.append(
            {
                "stage_index": stage_index,
                "stage_kind": "current-combine",
                "subset_size": size,
                "input_current_ids": sorted(input_current_ids),
                "output_current_ids": sorted(output_current_ids),
                "input_value_slot_ids": sorted(input_value_slot_ids),
                "output_value_slot_ids": sorted(output_value_slot_ids),
                "interaction_count": len(interactions),
                "interactions": interactions,
            }
        )
    return stages


def _runtime_interaction_record(
    dag: GenericDAG,
    model: Model,
    *,
    interaction_id: int,
    slot_by_current_id: dict[int, dict[str, object]],
    value_slot_by_current_variant: dict[tuple[int, str], dict[str, object]],
    momentum_slot_by_mask: dict[int, int],
    coupling_name_cache: dict[tuple[int, tuple[int, ...], tuple[float, ...]], list[str | None]],
    lowering_payload_cache: dict[int, dict[str, object]],
) -> dict[str, object]:
    interaction = dag.interactions[interaction_id]
    left = dag.currents[interaction.left_id]
    right = dag.currents[interaction.right_id]
    result = dag.currents[interaction.result_id]
    rule = model.vertex_lowering_rule(interaction.vertex_kind)
    left_value_slot = _input_value_slot(
        left,
        model,
        value_slot_by_current_variant,
    )
    right_value_slot = _input_value_slot(
        right,
        model,
        value_slot_by_current_variant,
    )
    result_value_slots = _result_value_slots(
        result,
        value_slot_by_current_variant,
    )
    coupling_key = (
        int(interaction.vertex_kind),
        tuple(int(pdg) for pdg in interaction.vertex_particles),
        tuple(float(value) for value in interaction.coupling),
    )
    coupling_parameter_names = coupling_name_cache.get(coupling_key)
    if coupling_parameter_names is None:
        coupling_parameter_names = _runtime_coupling_parameter_names(
            interaction.vertex_kind,
            interaction.vertex_particles,
            interaction.coupling,
        )
        coupling_name_cache[coupling_key] = coupling_parameter_names
    lowering_payload = lowering_payload_cache.get(int(interaction.vertex_kind))
    if lowering_payload is None:
        lowering_payload = _lowering_rule_payload(rule)
        lowering_payload_cache[int(interaction.vertex_kind)] = lowering_payload
    return {
        "interaction_id": interaction.id,
        "vertex_kind": interaction.vertex_kind,
        "vertex_particles": list(interaction.vertex_particles),
        "left_current_id": interaction.left_id,
        "right_current_id": interaction.right_id,
        "result_current_id": interaction.result_id,
        "left_slot": _slot_ref(slot_by_current_id[interaction.left_id]),
        "right_slot": _slot_ref(slot_by_current_id[interaction.right_id]),
        "result_slot": _slot_ref(slot_by_current_id[interaction.result_id]),
        "left_value_slot": _value_slot_ref(left_value_slot),
        "right_value_slot": _value_slot_ref(right_value_slot),
        "result_value_slots": [
            _value_slot_ref(slot) for slot in result_value_slots
        ],
        "result_requires_propagated_value": any(
            str(slot["variant"]) == "propagated" for slot in result_value_slots
        ),
        "result_requires_unpropagated_value": any(
            str(slot["variant"]) == "unpropagated" for slot in result_value_slots
        ),
        "momentum_slots": {
            "left": momentum_slot_by_mask[left.index.momentum_mask],
            "right": momentum_slot_by_mask[right.index.momentum_mask],
            "result": momentum_slot_by_mask[result.index.momentum_mask],
        },
        "coupling": list(interaction.coupling),
        "coupling_parameter_names": list(coupling_parameter_names),
        "color_weight": list(interaction.color_weight),
        "accumulation": "sum-into-result-current",
        "lowering": lowering_payload,
        "full_tensor_network_ready": interaction.full_tensor_network_ready,
    }


def _runtime_amplitude_stage_payload(
    dag: GenericDAG,
    *,
    slot_by_current_id: dict[int, dict[str, object]],
    value_slot_by_current_variant: dict[tuple[int, str], dict[str, object]],
    selected_color_sector_ids: set[int] | None = None,
) -> dict[str, object]:
    roots = []
    selected_roots = _selected_amplitude_roots(
        dag,
        selected_color_sector_ids=selected_color_sector_ids,
    )
    coherent_group_ids, color_groups = _amplitude_group_metadata(dag, selected_roots)
    color_contraction = build_color_contraction_plan(dag.color_plan, color_groups)
    for output_index, root in enumerate(selected_roots):
        roots.append(
            {
                "output_index": output_index,
                "root_id": output_index,
                "dag_root_id": root.id,
                "kind": root.kind,
                "left_current_id": root.left_id,
                "right_current_id": root.right_id,
                "left_slot": _slot_ref(slot_by_current_id[root.left_id]),
                "right_slot": _slot_ref(slot_by_current_id[root.right_id]),
                "left_value_slot": _value_slot_ref(
                    _amplitude_value_slot(
                        dag.currents[root.left_id],
                        value_slot_by_current_variant,
                    )
                ),
                "right_value_slot": _value_slot_ref(
                    _amplitude_value_slot(
                        dag.currents[root.right_id],
                        value_slot_by_current_variant,
                    )
                ),
                "vertex_kind": root.vertex_kind,
                "vertex_particles": (
                    list(root.vertex_particles)
                    if root.vertex_particles is not None
                    else None
                ),
                "coupling": list(root.coupling),
                "coupling_parameter_names": (
                    None
                    if root.vertex_kind is None or root.vertex_particles is None
                    else _runtime_coupling_parameter_names(
                        int(root.vertex_kind),
                        tuple(int(pdg) for pdg in root.vertex_particles),
                        root.coupling,
                    )
                ),
                "color_weight": list(root.color_weight),
                "color_sector_id": _root_color_sector_id(dag, root),
                "contraction": root.contraction,
                "coherent_group_id": coherent_group_ids[root.id],
                "helicity_weight": root.helicity_weight,
            }
        )
    return {
        "stage_kind": "amplitude-roots",
        "output_count": len(roots),
        "selected_color_sector_ids": (
            None
            if selected_color_sector_ids is None
            else sorted(selected_color_sector_ids)
        ),
        "roots": roots,
        "final_reduction": {
            "status": (
                "sparse-color-contraction"
                if color_contraction is not None
                else "coherent-leading-color-diagonal"
            ),
            "operation": (
                "sum root outputs into coherent helicity/colour amplitudes, "
                "then apply the requested colour contraction"
            ),
        },
        "color_contraction": (
            None if color_contraction is None else color_contraction.to_json_dict()
        ),
    }


def _amplitude_group_metadata(
    dag: GenericDAG,
    roots,
) -> tuple[dict[int, int], tuple[ColorGroupDescriptor, ...]]:
    """Map root contributions to physical helicity/colour coherent sums.

    A schema-v2 amplitude root is one current-closure contribution.  Roots
    with the same source-current ancestry and LC colour sector are diagrams of
    the same helicity/colour amplitude and must be summed before squaring.
    """

    ids: dict[tuple[object, ...], int] = {}
    result: dict[int, int] = {}
    group_descriptors: dict[int, ColorGroupDescriptor] = {}
    source_by_ancestry: dict[int, tuple[object, ...]] = {}
    physical_sources_by_ancestry: dict[int, tuple[object, ...]] = {}
    for current in dag.currents:
        if not current.is_source:
            continue
        source_by_ancestry[int(current.index.helicity_ancestry)] = (
            int(current.source_leg_label or 0),
            int(current.index.particle_id),
            int(current.index.chirality),
            current.index.spin_state,
            current.source_helicity,
        )
    source_ancestry_is_single_bit = all(
        bit > 0 and bit & (bit - 1) == 0 for bit in source_by_ancestry
    )

    def physical_sources_for_ancestry(ancestry: int) -> tuple[object, ...]:
        cached = physical_sources_by_ancestry.get(ancestry)
        if cached is not None:
            return cached
        if source_ancestry_is_single_bit:
            sources: list[tuple[object, ...]] = []
            remaining = int(ancestry)
            while remaining:
                bit = remaining & -remaining
                source_key = source_by_ancestry.get(bit)
                if source_key is not None:
                    sources.append(source_key)
                remaining ^= bit
            physical_sources = tuple(sorted(sources))
        else:
            physical_sources = tuple(
                sorted(
                    source_key
                    for bit, source_key in source_by_ancestry.items()
                    if ancestry & bit
                )
            )
        physical_sources_by_ancestry[ancestry] = physical_sources
        return physical_sources

    for root in roots:
        left = dag.currents[root.left_id].index
        right = dag.currents[root.right_id].index
        ancestry = int(left.helicity_ancestry | right.helicity_ancestry)
        physical_sources = physical_sources_for_ancestry(ancestry)
        sector_id = _root_color_sector_id(dag, root)
        sector = dag.color_plan.sector(sector_id)
        sector_word = (
            tuple(sector.word_labels or sector.color_words[0])
            if sector is not None and sector.color_words
            else sector_id
        )
        color_key = (
            left.color_state.accuracy,
            sector_word,
            tuple(sorted(set(left.color_state.basis_key) | set(right.color_state.basis_key))),
            tuple(root.color_weight),
        )
        key = (
            physical_sources or ancestry,
            color_key,
        )
        group_id = ids.setdefault(key, len(ids))
        result[root.id] = group_id
        if group_id not in group_descriptors:
            word = ()
            if sector is not None:
                word = tuple(sector.word_labels or sector.color_words[0])
            group_descriptors[group_id] = ColorGroupDescriptor(
                group_id=group_id,
                helicity_key=physical_sources or (ancestry,),
                sector_id=int(sector_id),
                word=word,
                helicity_weight=float(root.helicity_weight),
            )
    return result, tuple(
        group_descriptors[group_id]
        for group_id in sorted(group_descriptors)
    )


def _root_color_sector_id(dag: GenericDAG, root) -> int:
    if getattr(root, "color_sector_id", None) is not None:
        return int(root.color_sector_id)
    return int(dag.currents[root.left_id].index.color_state.sector_id)


def _slot_ref(slot: dict[str, object]) -> dict[str, object]:
    return {
        "current_id": slot["current_id"],
        "component_start": slot["component_start"],
        "component_stop": slot["component_stop"],
        "dimension": slot["dimension"],
    }


def _value_slot_ref(slot: dict[str, object]) -> dict[str, object]:
    return {
        "value_slot_id": slot["value_slot_id"],
        "current_id": slot["current_id"],
        "variant": slot["variant"],
        "component_start": slot["component_start"],
        "component_stop": slot["component_stop"],
        "dimension": slot["dimension"],
    }


def _input_value_slot(
    current,
    model: Model,
    value_slot_by_current_variant: dict[tuple[int, str], dict[str, object]],
) -> dict[str, object]:
    if current.is_source:
        variant = "source"
    else:
        propagator = model.propagator_lowering_rule(
            current.index.particle_id,
            current.index.chirality,
        )
        variant = "propagated" if propagator.applies_propagator else "unpropagated"
    return value_slot_by_current_variant[(current.id, variant)]


def _result_value_slots(
    current,
    value_slot_by_current_variant: dict[tuple[int, str], dict[str, object]],
) -> tuple[dict[str, object], ...]:
    slots = []
    for variant in ("unpropagated", "propagated"):
        slot = value_slot_by_current_variant.get((current.id, variant))
        if slot is not None:
            slots.append(slot)
    if not slots and current.is_source:
        slots.append(value_slot_by_current_variant[(current.id, "source")])
    return tuple(slots)


def _amplitude_value_slot(
    current,
    value_slot_by_current_variant: dict[tuple[int, str], dict[str, object]],
) -> dict[str, object]:
    if current.is_source:
        return value_slot_by_current_variant[(current.id, "source")]
    return value_slot_by_current_variant[(current.id, "unpropagated")]


def _lowering_rule_payload(rule: object) -> dict[str, object]:
    return {
        "kind": getattr(rule, "kind"),
        "backend": getattr(rule, "backend"),
        "tensor_names": list(getattr(rule, "tensor_names", ())),
        "expression_head": getattr(rule, "expression_head", ""),
        "full_tensor_network_ready": getattr(rule, "full_tensor_network_ready"),
        "description": getattr(rule, "description", ""),
        "kernel": getattr(rule, "kernel", ""),
        "input_roles": list(getattr(rule, "input_roles", ())),
        "output_role": getattr(rule, "output_role", ""),
        "coupling_mode": getattr(rule, "coupling_mode", "none"),
    }


def _leg_by_label(dag: GenericDAG, label: int):
    for leg in dag.process.legs:
        if leg.label == label:
            return leg
    return None


def _spin_state_json(spin_state: object) -> object:
    if isinstance(spin_state, tuple):
        return list(spin_state)
    return spin_state


def _bigint_json(value: int) -> str:
    integer = int(value)
    if integer.bit_length() <= 63:
        return str(integer)
    return hex(integer)


def _json_safe_bigints(value: object) -> object:
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value if value.bit_length() <= 63 else hex(value)
    if isinstance(value, dict):
        return {key: _json_safe_bigints(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_json_safe_bigints(item) for item in value]
    if isinstance(value, tuple):
        return [_json_safe_bigints(item) for item in value]
    return value


def _mask_labels(mask: int) -> tuple[int, ...]:
    return tuple(index + 1 for index in range(mask.bit_length()) if mask & (1 << index))


def _initial_state_average_factor(initial_pdgs: list[int]) -> int:
    factor = 1
    for pdg in initial_pdgs:
        if pdg == 21:
            factor *= 2 * 8
        elif 1 <= abs(pdg) <= 6:
            factor *= 2 * 3
        else:
            factor *= 2
    return factor


def _final_state_identical_factor(final_pdgs: list[int]) -> int:
    factor = 1
    for multiplicity in Counter(final_pdgs).values():
        factor *= math.factorial(multiplicity)
    return factor


def _quark_line_partner_factor(process: CanonicalProcessIR) -> int:
    """Return the LC partner factor recorded for diagnostics.

    Generic schema-v2 matrix elements are normalized as partonic squared
    matrix elements.  Open-line partner or multichannel bookkeeping belongs to
    phase-space integration, not to the pointwise ME normalization.  Keep this
    as a stable manifest field for readers, but do not divide by it.
    """

    del process
    return 1


_GENERIC_DAG_CHECK_STANDALONE = """#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import json
import math
import os
import statistics
import sys
import time
from decimal import Decimal
from pathlib import Path


GREEN = "\\033[32m"
CYAN = "\\033[36m"
YELLOW = "\\033[33m"
RED = "\\033[31m"
RESET = "\\033[0m"


def color(text: object, code: str) -> str:
    return f"{code}{text}{RESET}" if sys.stdout.isatty() else str(text)


def import_rusticol(root: Path, explicit_folder: str | None):
    candidates = []
    if explicit_folder:
        candidates.append(Path(explicit_folder).expanduser())
    candidates.extend(
        [
            root / "rusticol",
            root.parent / "rusticol",
            root.parent.parent / "rusticol",
            root.parent.parent.parent / "rusticol",
        ]
    )
    for candidate in candidates:
        if candidate.exists():
            sys.path.insert(0, str(candidate))
    try:
        import rusticol  # type: ignore[import-not-found]
    except ModuleNotFoundError as error:
        searched = ", ".join(str(path) for path in candidates)
        raise SystemExit(
            "could not import rusticol; install pyAmpliCol dependencies or pass "
            f"--rusticol-folder. searched: {searched}"
        ) from error
    return rusticol


def rusticol_build_metadata(rusticol):
    module = rusticol
    if not callable(getattr(module, "build_profile", None)):
        try:
            module = importlib.import_module("rusticol.rusticol")
        except Exception:
            module = rusticol
    build_profile = getattr(module, "build_profile", None)
    build_target = getattr(module, "build_target", None)
    return {
        "profile": build_profile() if callable(build_profile) else None,
        "target": build_target() if callable(build_target) else None,
        "module_path": str(getattr(module, "__file__", getattr(rusticol, "__file__", ""))),
    }


def require_release_rusticol(metadata):
    if metadata.get("profile") == "release":
        return
    if os.environ.get("PYAMPLICOL_ALLOW_DEBUG_RUSTICOL") == "1":
        return
    if metadata.get("profile") is None:
        raise SystemExit(
            "rusticol does not expose build-profile metadata; reinstall with "
            "`maturin develop --release` or rerun the pyAmpliCol dependency installer"
        )
    raise SystemExit(
        "rusticol was built with Cargo profile "
        f"{metadata.get('profile')!r}; profiling requires the release PyO3 extension"
    )


def load_points(root: Path, precision: int):
    payload = json.loads((root / "validation_momenta.json").read_text())
    if payload.get("available") is False or not payload.get("points"):
        return None, payload.get("error", "no validation momenta are bundled")
    points = [
        [
            [Decimal(str(component)) for component in particle["momentum"]]
            for particle in point
        ]
        for point in payload["points"]
    ]
    if precision == 16:
        try:
            import numpy as np  # type: ignore[import-not-found]
        except ModuleNotFoundError:
            return points, None
        return np.asarray(points, dtype=np.float64), None
    return points, None


def repeat_points(points, count: int):
    if hasattr(points, "shape"):
        import numpy as np  # type: ignore[import-not-found]

        reps = int(math.ceil(count / max(int(points.shape[0]), 1)))
        return np.tile(points, (reps, 1, 1))[:count]
    reps = int(math.ceil(count / max(len(points), 1)))
    return (points * reps)[:count]


def evaluate(runtime, points, precision: int):
    if precision == 16:
        return runtime.evaluate(points)
    return runtime.evaluate_with_prec(points, precision)


def profile(runtime, points, precision: int, target_s: float, batch_size: int):
    samples = []
    core_samples = []
    elapsed = 0.0
    while len(samples) < 8 or elapsed < target_s:
        batch = repeat_points(points, batch_size)
        start = time.perf_counter()
        payload = dict(runtime.profile(batch, precision=precision))
        wall = time.perf_counter() - start
        elapsed += wall
        samples.append(wall / batch_size)
        core_samples.append(
            (
                float(payload.get("stage_evaluator_time_s", 0.0))
                + float(payload.get("amplitude_evaluator_time_s", 0.0))
            )
            / max(int(payload.get("points", batch_size)), 1)
        )
    wall_mean = statistics.mean(samples)
    core_mean = statistics.mean(core_samples)
    return {
        "samples": len(samples) * batch_size,
        "wall_us_per_point": wall_mean * 1.0e6,
        "core_us_per_point": core_mean * 1.0e6,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Standalone generic pyAmpliCol process check")
    parser.add_argument("--precision", type=int, default=16)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--target-runtime", type=float, default=10.0)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--rusticol-folder")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent
    rusticol = import_rusticol(root, args.rusticol_folder)
    build_metadata = rusticol_build_metadata(rusticol)
    require_release_rusticol(build_metadata)
    runtime = rusticol.Runtime.load(str(root))
    metadata = dict(runtime.metadata())
    print(color("pyAmpliCol generic process check", CYAN))
    print(f"process: {metadata.get('process')}")
    print(f"schema:  {metadata.get('schema_version')}")
    print(f"runtime: {color('available', GREEN)}")
    print(
        "rust:    "
        f"{build_metadata.get('profile') or 'unknown'} "
        f"[{build_metadata.get('target') or 'unknown'}]"
    )

    points, error = load_points(root, args.precision)
    if points is None:
        print(color(f"validation momenta unavailable: {error}", YELLOW))
        return 0
    values = evaluate(runtime, points, args.precision)
    print(f"values:  {[float(value) for value in values]}")
    if args.profile:
        payload = profile(
            runtime,
            points,
            args.precision,
            args.target_runtime,
            max(int(args.batch_size), 1),
        )
        print(
            "timing:  "
            f"wall={payload['wall_us_per_point']:.4g} us/point, "
            f"core={payload['core_us_per_point']:.4g} us/point, "
            f"samples={payload['samples']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
"""


def _dag_lowering_status(
    dag: GenericDAG,
    model: Model,
    *,
    include_color_sector_summaries: bool = True,
) -> dict[str, object]:
    vertex_kinds = [
        interaction.vertex_kind for interaction in dag.interactions
    ] + [
        int(root.vertex_kind)
        for root in dag.amplitude_roots
        if root.vertex_kind is not None
    ]
    counts = Counter(vertex_kinds)
    ready: set[int] = set()
    pending: set[int] = set()
    unimplemented: set[int] = set()
    for kind in counts:
        rule = model.vertex_lowering_rule(kind)
        if rule.backend == "unimplemented":
            unimplemented.add(kind)
        elif rule.full_tensor_network_ready:
            ready.add(kind)
        else:
            pending.add(kind)
    propagators = _dag_propagator_lowering_status(dag, model)
    internal_current_color_sectors = tuple(
        sorted({current.index.color_state.sector_id for current in dag.currents})
    )
    amplitude_color_sectors = contributing_color_sector_ids(dag)
    current_color_sectors = amplitude_color_sectors or internal_current_color_sectors
    payload: dict[str, object] = {
        "current_count": len(dag.currents),
        "source_count": len(dag.sources),
        "interaction_count": len(dag.interactions),
        "current_color_sectors": list(current_color_sectors),
        "current_color_sector_count": len(current_color_sectors),
        "internal_current_color_sectors": list(internal_current_color_sectors),
        "internal_current_color_sector_count": len(internal_current_color_sectors),
        "amplitude_color_sectors": list(amplitude_color_sectors),
        "amplitude_color_sector_count": len(amplitude_color_sectors),
        "required_vertex_kind_counts": [
            [kind, count] for kind, count in sorted(counts.items())
        ],
        "ready_vertex_kinds": sorted(ready),
        "pending_vertex_kinds": sorted(pending),
        "unimplemented_vertex_kinds": sorted(unimplemented),
        **propagators,
    }
    if include_color_sector_summaries:
        payload["color_sector_summaries"] = _dag_color_sector_summaries(dag, model)
    else:
        payload["color_sector_summaries"] = []
        payload["color_sector_summaries_omitted"] = True
        payload["color_sector_summaries_omitted_reason"] = (
            "compact main artifact; selected runtime sidecars keep detailed "
            "sector-local lowering manifests"
        )
    return payload


def _dag_propagator_lowering_status(
    dag: GenericDAG,
    model: Model,
) -> dict[str, object]:
    usage = _runtime_current_usage(dag)
    counts: Counter[str] = Counter()
    ready: set[str] = set()
    pending: set[str] = set()
    unimplemented: set[str] = set()
    for current in dag.currents:
        if current.is_source:
            continue
        if not usage[current.id]["used_as_interaction_input"]:
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
        "required_propagator_kernel_counts": [
            [kernel, count] for kernel, count in sorted(counts.items())
        ],
        "ready_propagator_kernels": sorted(ready),
        "pending_propagator_kernels": sorted(pending),
        "unimplemented_propagator_kernels": sorted(unimplemented),
    }


def _dag_color_sector_summaries(
    dag: GenericDAG,
    model: Model,
) -> list[dict[str, object]]:
    root_sectors = contributing_color_sector_ids(dag)
    sectors = (
        root_sectors
        if root_sectors
        else tuple(sorted({current.index.color_state.sector_id for current in dag.currents}))
    )
    summaries: list[dict[str, object]] = []
    for sector in sectors:
        sector_dag = (
            filter_dag_to_color_sectors(dag, {sector})
            if root_sectors
            else dag
        )
        vertex_kinds = [
            interaction.vertex_kind
            for interaction in sector_dag.interactions
        ] + [
            int(root.vertex_kind)
            for root in sector_dag.amplitude_roots
            if root.vertex_kind is not None
        ]
        counts = Counter(vertex_kinds)
        ready: set[int] = set()
        pending: set[int] = set()
        unimplemented: set[int] = set()
        for kind in counts:
            rule = model.vertex_lowering_rule(kind)
            if rule.backend == "unimplemented":
                unimplemented.add(kind)
            elif rule.full_tensor_network_ready:
                ready.add(kind)
            else:
                pending.add(kind)
        summaries.append(
            {
                "color_sector": sector,
                "current_count": len(sector_dag.currents),
                "interaction_count": len(sector_dag.interactions),
                "closure_count": len(sector_dag.amplitude_roots),
                "amplitude_root_count": len(sector_dag.amplitude_roots),
                "required_vertex_kind_counts": [
                    [kind, count] for kind, count in sorted(counts.items())
                ],
                "ready_vertex_kinds": sorted(ready),
                "pending_vertex_kinds": sorted(pending),
                "unimplemented_vertex_kinds": sorted(unimplemented),
            }
        )
    return summaries


def _dag_stage_plan_payload(dag: GenericDAG, model: Model) -> dict[str, object]:
    interactions_by_size: dict[int, list[int]] = {}
    currents_by_size: dict[int, set[int]] = {}
    for interaction in dag.interactions:
        result = dag.currents[interaction.result_id]
        size = len(result.index.external_labels)
        interactions_by_size.setdefault(size, []).append(interaction.id)
        currents_by_size.setdefault(size, set()).add(result.id)
    current_stages = []
    for stage_index, size in enumerate(sorted(interactions_by_size), start=1):
        counts = Counter(
            dag.interactions[interaction_id].vertex_kind
            for interaction_id in interactions_by_size[size]
        )
        ready: set[int] = set()
        pending: set[int] = set()
        unimplemented: set[int] = set()
        for kind in counts:
            rule = model.vertex_lowering_rule(kind)
            if rule.backend == "unimplemented":
                unimplemented.add(kind)
            elif rule.full_tensor_network_ready:
                ready.add(kind)
            else:
                pending.add(kind)
        current_stages.append(
            {
                "stage_index": stage_index,
                "subset_size": size,
                "current_ids": sorted(currents_by_size.get(size, ())),
                "interaction_ids": interactions_by_size[size],
                "interaction_count": len(interactions_by_size[size]),
                "required_vertex_kind_counts": [
                    [kind, count] for kind, count in sorted(counts.items())
                ],
                "ready_vertex_kinds": sorted(ready),
                "pending_vertex_kinds": sorted(pending),
                "unimplemented_vertex_kinds": sorted(unimplemented),
            }
        )
    return {
        "stage_count": len(current_stages) + 1,
        "current_stages": current_stages,
        "amplitude_stage": {
            "amplitude_root_ids": [root.id for root in dag.amplitude_roots],
            "closure_ids": [root.id for root in dag.amplitude_roots],
            "closure_count": len(dag.amplitude_roots),
            "amplitude_root_count": len(dag.amplitude_roots),
        },
    }


def load_generic_process_manifest(path: str | Path) -> dict[str, object]:
    manifest_path = Path(path).expanduser()
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("generic process manifest is not a JSON object")
    if payload.get("kind") != GENERIC_PROCESS_MANIFEST_KIND:
        raise ValueError(
            f"unsupported generic process manifest kind: {payload.get('kind')!r}"
        )
    if payload.get("schema_version") != GENERIC_PROCESS_SCHEMA_VERSION:
        raise ValueError(
            "unsupported generic process manifest schema version: "
            f"{payload.get('schema_version')!r}"
        )
    return payload


def load_generic_process_set_manifest(path: str | Path) -> dict[str, object]:
    manifest_path = Path(path).expanduser()
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("generic process-set manifest is not a JSON object")
    if payload.get("kind") != GENERIC_PROCESS_SET_MANIFEST_KIND:
        raise ValueError(
            "unsupported generic process-set manifest kind: "
            f"{payload.get('kind')!r}"
        )
    if payload.get("schema_version") != GENERIC_PROCESS_SCHEMA_VERSION:
        raise ValueError(
            "unsupported generic process-set manifest schema version: "
            f"{payload.get('schema_version')!r}"
        )
    return payload


__all__ = [
    "GENERIC_DAG_PROCESS_ARTIFACT_KIND",
    "GENERIC_DAG_PROCESS_SET_ARTIFACT_KIND",
    "GENERIC_LC_REPLAY_PARTITION_ARTIFACT_KIND",
    "GENERIC_PROCESS_MANIFEST_KIND",
    "GENERIC_PROCESS_SET_MANIFEST_KIND",
    "GENERIC_PROCESS_SCHEMA_VERSION",
    "GenericProcessManifest",
    "GenericProcessSetManifest",
    "build_generic_process_manifest",
    "build_generic_process_set_manifest",
    "load_generic_process_manifest",
    "load_generic_process_set_manifest",
    "select_leading_color_sector_ids",
    "select_leading_color_sector_ids_from_plan",
    "write_generic_process_manifest",
    "write_generic_process_set_manifest",
    "write_generic_dag_process_artifact",
    "write_generic_dag_process_set_artifact",
    "write_lc_topology_replay_partition_artifact",
]
