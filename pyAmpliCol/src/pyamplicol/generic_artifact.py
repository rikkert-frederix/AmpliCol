from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence, cast

from .color_plan import (
    GenericColorPlan,
    build_color_plan,
    lc_line_pairing_representative_ids,
    lc_topology_replay_safe_groups,
)
from .core_types import ExternalMomentum, NativeEvaluationError
from .generic_dag import (
    GenericDAG,
    GenericDAGCompiler,
    contributing_color_sector_ids,
    filter_dag_to_color_sectors,
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

_FULL_COLOR_PLAN_SERIALIZATION_SECTOR_LIMIT = 1000


@dataclass(frozen=True)
class GenericProcessManifest:
    """Process-generic planning manifest for the future schema-v2 runtime."""

    dag: GenericDAG
    model: Model
    color_plan: GenericColorPlan

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
        color_ready = self.color_plan.ready_for_leading_colour
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
    max_currents: int = 20000,
    max_color_sectors: int = 20000,
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
) -> GenericProcessManifest:
    model = model or AmplicolSMLeadingColorModel()
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
    full_color_plan = build_color_plan(
        dag.process,
        color_accuracy=dag.process.color_accuracy,
        options=options,
        max_sectors=max_color_sectors,
        reference_color_order=reference_color_order,
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
    return GenericProcessManifest(dag=dag, model=model, color_plan=color_plan)


def build_generic_process_set_manifest(
    process_string: str,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    color_accuracy: str = "lc",
    max_currents: int = 20000,
    max_color_sectors: int = 20000,
    selected_color_sector_ids: set[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
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
                manifest.dag.currents[root.left_id].index.color_state.sector_id
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
        json.dumps(_json_safe_bigints(manifest.to_json_dict()), indent=2, sort_keys=True),
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


def write_generic_dag_process_artifact(
    process: str | CanonicalProcessIR | GenericDAG | GenericProcessManifest,
    output_dir: str | Path,
    *,
    model: Model | None = None,
    options: ProcessOptions | None = None,
    color_accuracy: str = "lc",
    max_currents: int = 50000,
    max_color_sectors: int = 20000,
    evaluator_backend: str = "compiled-complex",
    compiled_preset: str = "runtime-o3",
    batch_size: int = 64,
    emit_stage_evaluator_artifacts: bool = False,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    lc_topology_replay: bool = False,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
) -> tuple[Path, dict[str, object]]:
    """Write a schema-v2 generic DAG process artifact.

    This is the production-facing artifact for the generic refactor. It stores
    the full model-driven DAG, requested evaluator settings, symbolic
    stage-blueprint lowering, and optionally serialized generic stage
    evaluators loadable by Rusticol schema-v2 execution.
    """

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
        )
    )
    output_path = Path(output_dir).expanduser()
    output_path.mkdir(parents=True, exist_ok=True)
    plan_path = write_generic_process_manifest(generic_manifest, output_path)
    payload = _generic_dag_process_artifact_payload(
        generic_manifest,
        plan_path=plan_path.name,
        evaluator_backend=evaluator_backend,
        compiled_preset=compiled_preset,
        batch_size=batch_size,
        stage_evaluator_artifact_dir=(
            output_path if emit_stage_evaluator_artifacts else None
        ),
        stage_evaluator_compiler=stage_evaluator_compiler,
        symbolica_settings=symbolica_settings,
        merge_evaluators_strategy=merge_evaluators_strategy,
        verbose_evaluator_build=verbose_evaluator_build,
        jit_compile=jit_compile,
        progress_callback=progress_callback,
        reference_color_order=reference_color_order,
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
    )
    manifest_path = output_path / "process_manifest.json"
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (output_path / "check_standalone.py").write_text(
        _GENERIC_DAG_CHECK_STANDALONE,
        encoding="utf-8",
    )
    _write_generic_validation_momenta(generic_manifest, output_path)
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
            manifest.dag.currents[root.left_id].index.color_state.sector_id
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
                return _lc_colored_word_sibling_sector_ids(
                    manifest.color_plan,
                    sector,
                    active_sector_ids=set(root_sectors),
                    reference_order=wanted,
                )
            legacy_order_words = getattr(sector, "legacy_order_words", ())
            if wanted in legacy_order_words:
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
                    return _lc_colored_word_sibling_sector_ids(
                        color_plan,
                        sector,
                        reference_order=wanted,
                    )
        for sector in color_plan.sectors:
            if wanted in sector.color_words:
                return _lc_colored_word_sibling_sector_ids(
                    color_plan,
                    sector,
                    reference_order=wanted,
                )
            legacy_order_words = getattr(sector, "legacy_order_words", ())
            if wanted in legacy_order_words:
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
    max_currents: int = 50000,
    max_color_sectors: int = 20000,
    evaluator_backend: str = "compiled-complex",
    compiled_preset: str = "runtime-o3",
    batch_size: int = 64,
    emit_stage_evaluator_artifacts: bool = False,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    lc_topology_replay: bool = False,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
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
            stage_evaluator_compiler=stage_evaluator_compiler,
            symbolica_settings=symbolica_settings,
            merge_evaluators_strategy=merge_evaluators_strategy,
            verbose_evaluator_build=verbose_evaluator_build,
            jit_compile=jit_compile,
            progress_callback=progress_callback,
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
        },
        "lc_topology_replay": bool(lc_topology_replay),
    }


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


def _generic_dag_process_artifact_payload(
    manifest: GenericProcessManifest,
    *,
    plan_path: str,
    evaluator_backend: str,
    compiled_preset: str,
    batch_size: int,
    stage_evaluator_artifact_dir: Path | None = None,
    stage_evaluator_compiler: Any | None = None,
    symbolica_settings: Any | None = None,
    merge_evaluators_strategy: bool = False,
    verbose_evaluator_build: bool = False,
    jit_compile: bool = True,
    progress_callback: Any | None = None,
    reference_color_order: Sequence[int] | None = None,
    selected_color_sector_ids: set[int] | None = None,
    lc_topology_replay: bool = False,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Sequence[int] | None = None,
    ignored_vertex_kinds: Sequence[int] | None = None,
) -> dict[str, object]:
    from .generic_stage_compiler import (
        build_generic_stage_compiler_blueprint,
        write_generic_stage_evaluator_artifacts,
    )

    full_plan_payload = manifest.to_json_dict()
    runtime_sector_ids = _runtime_color_sector_ids(
        manifest.dag,
        selected_color_sector_ids=selected_color_sector_ids,
    )
    runtime_manifest = _runtime_manifest_for_color_sectors(
        manifest,
        runtime_sector_ids,
    )
    plan_payload = runtime_manifest.to_json_dict()
    runtime_schema = cast(dict[str, Any], plan_payload["runtime_schema"])
    stage_blueprint = build_generic_stage_compiler_blueprint(
        runtime_manifest,
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
            progress_callback=progress_callback,
        )
        compiled_payload["stage_evaluators"] = stage_evaluator_payload
        if bool(stage_evaluator_payload.get("runtime_available", False)):
            compiled_payload["runtime_available"] = True
            compiled_payload["runtime_unavailable_message"] = None
    return {
        "schema_version": GENERIC_PROCESS_SCHEMA_VERSION,
        "kind": GENERIC_DAG_PROCESS_ARTIFACT_KIND,
        "artifact_class": "generic-dag-schema-v2",
        "process": manifest.process,
        "key": manifest.key,
        "color_accuracy": manifest.dag.process.color_accuracy,
        "model": plan_payload["model"],
        "external_pdg_order": list(manifest.external_pdg_order),
        "outgoing_pdg_order": list(manifest.outgoing_pdg_order),
        "process_ir": plan_payload["process_ir"],
        "generic_plan_path": plan_path,
        "lc_topology_reuse": full_plan_payload["lc_topology_reuse"],
        "runtime_lc_topology_reuse": plan_payload["lc_topology_reuse"],
        "runtime_lc_topology_replay": topology_replay_payload,
        "planning_status": full_plan_payload["planning_status"],
        "lowering_status": plan_payload["lowering_status"],
        "runtime_schema": plan_payload["runtime_schema"],
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
    )


def _generic_runtime_schema_payload(
    dag: GenericDAG,
    model: Model,
    *,
    selected_color_sector_ids: set[int] | None = None,
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
    source_component_count = sum(
        _schema_int(slot_by_current_id[source_id]["dimension"])
        for source_id in dag.sources
    )
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
            "parameter_count_if_flattened": (
                source_component_count + 4 * len(momentum_slots)
            ),
            "value_component_count": (
                value_slots[-1]["component_stop"] if value_slots else 0
            ),
            "source_components_complex": True,
            "momentum_components_real": True,
            "real_valued_inputs": list(
                range(
                    source_component_count,
                    source_component_count + 4 * len(momentum_slots),
                )
            ),
        },
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
        "stages": _runtime_stage_payloads(
            dag,
            model,
            slot_by_current_id=slot_by_current_id,
            value_slot_by_current_variant=value_slot_by_current_variant,
            momentum_slot_by_mask=momentum_slot_by_mask,
        ),
        "amplitude_stage": _runtime_amplitude_stage_payload(
            dag,
            slot_by_current_id=slot_by_current_id,
            value_slot_by_current_variant=value_slot_by_current_variant,
            selected_color_sector_ids=selected_color_sector_ids,
        ),
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
        if dag.currents[root.left_id].index.color_state.sector_id
        in selected_color_sector_ids
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
    for stage_index, size in enumerate(sorted(interactions_by_size), start=1):
        interaction_ids = interactions_by_size[size]
        interactions = [
            _runtime_interaction_record(
                dag,
                model,
                interaction_id=interaction_id,
                slot_by_current_id=slot_by_current_id,
                value_slot_by_current_variant=value_slot_by_current_variant,
                momentum_slot_by_mask=momentum_slot_by_mask,
            )
            for interaction_id in interaction_ids
        ]
        runtime_interactions = cast(list[dict[str, Any]], interactions)
        output_current_ids = sorted(
            {
                _schema_int(interaction["result_current_id"])
                for interaction in runtime_interactions
            }
        )
        input_current_ids = sorted(
            {
                _schema_int(interaction["left_current_id"])
                for interaction in runtime_interactions
            }
            | {
                _schema_int(interaction["right_current_id"])
                for interaction in runtime_interactions
            }
        )
        input_value_slot_ids = sorted(
            {
                _schema_int(interaction["left_value_slot"]["value_slot_id"])
                for interaction in runtime_interactions
            }
            | {
                _schema_int(interaction["right_value_slot"]["value_slot_id"])
                for interaction in runtime_interactions
            }
        )
        output_value_slot_ids = sorted(
            {
                _schema_int(slot["value_slot_id"])
                for interaction in runtime_interactions
                for slot in interaction["result_value_slots"]
            }
        )
        stages.append(
            {
                "stage_index": stage_index,
                "stage_kind": "current-combine",
                "subset_size": size,
                "input_current_ids": input_current_ids,
                "output_current_ids": output_current_ids,
                "input_value_slot_ids": input_value_slot_ids,
                "output_value_slot_ids": output_value_slot_ids,
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
        "color_weight": list(interaction.color_weight),
        "accumulation": "sum-into-result-current",
        "lowering": _lowering_rule_payload(rule),
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
    coherent_group_ids = _amplitude_coherent_group_ids(dag, selected_roots)
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
                "color_weight": list(root.color_weight),
                "contraction": root.contraction,
                "coherent_group_id": coherent_group_ids[root.id],
                "helicity_weight": 1.0,
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
            "status": "pending-rusticol-schema-v2",
            "operation": (
                "sum helicity/color weighted squared complex root outputs "
                "after evaluator execution"
            ),
        },
    }


def _amplitude_coherent_group_ids(
    dag: GenericDAG,
    roots,
) -> dict[int, int]:
    """Map root contributions to physical helicity/colour coherent sums.

    A schema-v2 amplitude root is one current-closure contribution.  Roots
    with the same source-current ancestry and LC colour sector are diagrams of
    the same helicity/colour amplitude and must be summed before squaring.
    """

    ids: dict[tuple[object, ...], int] = {}
    result: dict[int, int] = {}
    source_by_ancestry: dict[int, tuple[object, ...]] = {}
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
    for root in roots:
        left = dag.currents[root.left_id].index
        right = dag.currents[root.right_id].index
        ancestry = int(left.helicity_ancestry | right.helicity_ancestry)
        physical_sources = tuple(
            sorted(
                source_key
                for bit, source_key in source_by_ancestry.items()
                if ancestry & bit
            )
        )
        sector = (
            dag.color_plan.sector(left.color_state.sector_id)
            if left.color_state.accuracy == "lc"
            else None
        )
        color_key = (
            left.color_state.accuracy,
            (
                tuple(sector.word_labels)
                if sector is not None
                else left.color_state.sector_id
            ),
            tuple(sorted(set(left.color_state.basis_key) | set(right.color_state.basis_key))),
            tuple(root.color_weight),
        )
        key = (
            physical_sources or ancestry,
            color_key,
        )
        group_id = ids.setdefault(key, len(ids))
        result[root.id] = group_id
    return result


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


def _dag_lowering_status(dag: GenericDAG, model: Model) -> dict[str, object]:
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
    current_color_sectors = tuple(
        sorted({current.index.color_state.sector_id for current in dag.currents})
    )
    return {
        "current_count": len(dag.currents),
        "source_count": len(dag.sources),
        "interaction_count": len(dag.interactions),
        "current_color_sectors": list(current_color_sectors),
        "current_color_sector_count": len(current_color_sectors),
        "color_sector_summaries": _dag_color_sector_summaries(dag, model),
        "required_vertex_kind_counts": [
            [kind, count] for kind, count in sorted(counts.items())
        ],
        "ready_vertex_kinds": sorted(ready),
        "pending_vertex_kinds": sorted(pending),
        "unimplemented_vertex_kinds": sorted(unimplemented),
        **propagators,
    }


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
    current_counts = Counter(
        current.index.color_state.sector_id for current in dag.currents
    )
    interaction_counts = Counter(
        dag.currents[interaction.result_id].index.color_state.sector_id
        for interaction in dag.interactions
    )
    root_counts = Counter(
        dag.currents[root.left_id].index.color_state.sector_id
        for root in dag.amplitude_roots
    )
    sectors = sorted(set(current_counts) | set(interaction_counts) | set(root_counts))
    summaries: list[dict[str, object]] = []
    for sector in sectors:
        vertex_kinds = [
            interaction.vertex_kind
            for interaction in dag.interactions
            if dag.currents[interaction.result_id].index.color_state.sector_id
            == sector
        ] + [
            int(root.vertex_kind)
            for root in dag.amplitude_roots
            if root.vertex_kind is not None
            and dag.currents[root.left_id].index.color_state.sector_id == sector
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
                "current_count": current_counts[sector],
                "interaction_count": interaction_counts[sector],
                "closure_count": root_counts[sector],
                "amplitude_root_count": root_counts[sector],
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
]
