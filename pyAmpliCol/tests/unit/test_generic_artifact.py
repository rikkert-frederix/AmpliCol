from __future__ import annotations

import gzip
import json
import math
import pickle
import threading
from pathlib import Path

import pytest
import pyamplicol.generic_artifact as generic_artifact_module

from pyamplicol.generic_artifact import (
    GENERIC_DAG_PROCESS_ARTIFACT_KIND,
    GENERIC_DAG_PROCESS_SET_ARTIFACT_KIND,
    GENERIC_LC_REPLAY_PARTITION_ARTIFACT_KIND,
    GENERIC_PROCESS_MANIFEST_KIND,
    GENERIC_PROCESS_SET_MANIFEST_KIND,
    GENERIC_PROCESS_SCHEMA_VERSION,
    build_generic_process_set_manifest,
    build_generic_process_manifest,
    load_generic_process_manifest,
    load_generic_process_set_manifest,
    write_generic_dag_process_artifact,
    write_generic_dag_process_set_artifact,
    write_lc_topology_replay_partition_artifact,
    write_generic_process_set_manifest,
    write_generic_process_manifest,
    select_leading_color_sector_ids_from_plan,
    _current_signature_relation,
    _drop_currents_from_dag,
    _merge_currents_in_dag,
    _generic_runtime_schema_payload,
    _validate_numerical_rewrite_preserves_raw_sums,
    _json_safe_bigints,
)
from pyamplicol.generic_dag import (
    AmplitudeRoot,
    CurrentNode,
    GenericDAG,
    prune_dag_to_amplitude_roots,
)
from pyamplicol.color_plan import build_color_plan
from pyamplicol.generic_stage_compiler import (
    _parameter_builder,
    build_and_write_generic_stage_evaluator_artifacts,
    build_generic_stage_compiler_blueprint,
    write_generic_stage_evaluator_artifacts,
)
from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.process_ir import build_process_set_ir


def test_json_safe_bigints_serializes_large_bitsets_as_hex_strings() -> None:
    source = {
        "small": 42,
        "large": 1 << 20000,
        "nested": [True, 1 << 5000],
    }
    payload = _json_safe_bigints(source)

    assert payload is source
    assert payload["small"] == 42
    assert payload["large"].startswith("0x")
    assert payload["nested"][0] is True
    assert payload["nested"][1].startswith("0x")
    json.dumps(payload)


def test_generic_process_manifest_keeps_physical_and_outgoing_pdg_orders() -> None:
    manifest = build_generic_process_manifest("d d~ > z g")
    payload = manifest.to_json_dict()

    assert payload["schema_version"] == GENERIC_PROCESS_SCHEMA_VERSION
    assert payload["kind"] == GENERIC_PROCESS_MANIFEST_KIND
    assert payload["process"] == "d d~ > z g"
    assert payload["key"] == "d_dbar_to_z_g"
    assert payload["external_pdg_order"] == [1, -1, 23, 21]
    assert payload["outgoing_pdg_order"] == [-1, 1, 23, 21]
    assert payload["process_ir"]["labels"]["vectors"] == [3]
    assert payload["process_ir"]["labels"]["gluons"] == [4]
    assert payload["color_plan"]["sector_count"] == 1
    assert payload["color_plan"]["sectors"][0]["quark_lines"][0][
        "gluon_labels"
    ] == [4]
    assert payload["planning_status"] == {
        "color_ready": True,
        "color_sector_count": 1,
        "color_truncated": False,
        "current_ready": True,
        "generic_evaluator_ready": True,
        "has_amplitude_roots": True,
        "has_closure": True,
        "idenso_required": False,
    }
    assert payload["stage_plan"]["stage_count"] == (
        len(payload["stage_plan"]["current_stages"]) + 1
    )
    runtime_schema = payload["runtime_schema"]
    assert runtime_schema["kind"] == "pyamplicol-generic-dag-runtime-schema"
    assert runtime_schema["momentum_conventions"]["input_shape"] == [
        "batch",
        4,
        4,
    ]
    assert runtime_schema["momentum_conventions"]["incoming_labels"] == [1, 2]
    assert runtime_schema["normalization"]["color_factor"] == 9
    assert runtime_schema["normalization"]["average_factor"] == 36
    assert runtime_schema["normalization"]["couplings_in_stage_evaluators"] is True
    assert runtime_schema["current_storage"]["component_count"] > 0


def test_runtime_schema_reports_progress_for_long_interaction_passes() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    events: list[tuple[str, int, int]] = []

    _generic_runtime_schema_payload(
        manifest.dag,
        manifest.model,
        progress_callback=lambda phase, completed, total: events.append(
            (phase, completed, total)
        ),
    )

    interaction_total = len(manifest.dag.interactions)
    for phase in ("usage scan", "stage grouping", "interaction records"):
        phase_events = [event for event in events if event[0] == phase]
        assert phase_events
        assert phase_events[-1][1:] == (interaction_total, interaction_total)


def test_runtime_schema_precomputes_pure_gluon_all_sector_weight_context(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = build_generic_process_manifest(
        "g g > g g",
        selected_source_helicities={1: -1, 2: -1, 3: -1, 4: -1},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    original = generic_artifact_module._root_all_sector_weight
    contexts: list[bool | None] = []

    def recording_weight(dag, root, **kwargs):
        contexts.append(kwargs.get("has_multiple_lc_root_sectors"))
        return original(dag, root, **kwargs)

    monkeypatch.setattr(
        generic_artifact_module,
        "_root_all_sector_weight",
        recording_weight,
    )
    _generic_runtime_schema_payload(manifest.dag, manifest.model)

    assert len(contexts) == len(manifest.dag.amplitude_roots)
    assert set(contexts) == {True}


def test_zero_current_filter_prunes_massive_top_chirality_currents() -> None:
    unfiltered = build_generic_process_manifest(
        "g g > t t~ g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    filtered = build_generic_process_manifest(
        "g g > t t~ g",
        selected_color_sector_ids={0},
        numerical_current_merging=False,
    )
    report = filtered.to_json_dict()["generation_filters"]["zero_current"]

    assert report["enabled"] is True
    assert report["skipped"] is False
    assert report["validation"]["accepted"] is True
    assert report["zero_current_ids"] == [12, 13, 14, 15]
    assert report["removed_current_count"] == 4
    assert len(filtered.dag.currents) == len(unfiltered.dag.currents) - 4
    assert len(filtered.dag.interactions) == len(unfiltered.dag.interactions) - 16
    assert len(filtered.dag.amplitude_roots) == len(unfiltered.dag.amplitude_roots)


def test_raw_sum_validation_rejects_contributing_current_drop() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    candidate = _drop_currents_from_dag(
        manifest.dag,
        (manifest.dag.amplitude_roots[0].left_id,),
    )

    validation = _validate_numerical_rewrite_preserves_raw_sums(
        manifest.dag,
        candidate,
        manifest.model,
        sample_count=3,
        seed=12345,
    )

    assert validation["accepted"] is False
    assert validation["max_relative_difference"] > 1.0e-12


def test_raw_sum_validation_rejects_wrong_current_merge() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    current_map = {current.id: (current.id, 1.0) for current in manifest.dag.currents}
    current_map[manifest.dag.amplitude_roots[0].left_id] = (
        manifest.dag.amplitude_roots[1].left_id,
        1.0,
    )
    candidate = _merge_currents_in_dag(manifest.dag, current_map)

    validation = _validate_numerical_rewrite_preserves_raw_sums(
        manifest.dag,
        candidate,
        manifest.model,
        sample_count=3,
        seed=12345,
    )

    assert validation["accepted"] is False
    assert validation["max_relative_difference"] > 1.0e-12


def test_raw_sum_validation_rejects_weighted_single_initial_helicity() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > t t~ g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    dag = manifest.dag
    kept_bits: dict[int, int] = {}
    dropped_mask = 0
    for source_id in dag.sources:
        source = dag.currents[source_id]
        label = int(source.source_leg_label or 0)
        state = (
            int(source.source_helicity or 0),
            int(source.index.chirality),
        )
        bit = int(source.index.helicity_ancestry)
        if label == 1 and source.index.particle_id == -1:
            if state == (1, -1):
                kept_bits[label] = bit
            else:
                dropped_mask |= bit
        if label == 2 and source.index.particle_id == 1:
            if state == (-1, 1):
                kept_bits[label] = bit
            else:
                dropped_mask |= bit

    roots: list[AmplitudeRoot] = []
    for root in dag.amplitude_roots:
        left = dag.currents[root.left_id].index
        right = dag.currents[root.right_id].index
        ancestry = int(left.helicity_ancestry | right.helicity_ancestry)
        if ancestry & dropped_mask:
            continue
        if any(not (ancestry & bit) for bit in kept_bits.values()):
            continue
        roots.append(
            AmplitudeRoot(
                id=len(roots),
                kind=root.kind,
                left_id=root.left_id,
                right_id=root.right_id,
                color_weight=root.color_weight,
                color_sector_id=root.color_sector_id,
                vertex_kind=root.vertex_kind,
                vertex_particles=root.vertex_particles,
                coupling=root.coupling,
                contraction=root.contraction,
                helicity_weight=root.helicity_weight * 2.0,
            )
        )
    candidate = prune_dag_to_amplitude_roots(
        GenericDAG(
            process=dag.process,
            color_plan=dag.color_plan,
            currents=dag.currents,
            sources=dag.sources,
            interactions=dag.interactions,
            amplitude_roots=tuple(roots),
            truncated=dag.truncated,
        )
    )

    validation = _validate_numerical_rewrite_preserves_raw_sums(
        dag,
        candidate,
        manifest.model,
        sample_count=3,
        seed=12345,
    )

    assert len(candidate.currents) == 19
    assert validation["accepted"] is False
    assert validation["max_relative_difference"] > 1.0e-6


def test_lc_gluon_flavour_flow_aggregation_matches_amplicol_current_buckets() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > t t~ g g g g",
        selected_color_sector_ids={0},
        reference_color_order=(3, 5, 6, 7, 8, 1, 2, 4),
        numerical_filter_current=False,
        numerical_current_merging=False,
    )

    aggregation = manifest.structural_current_aggregation
    assert aggregation is not None
    assert aggregation["before"] == {
        "current_count": 438,
        "source_count": 16,
        "interaction_count": 2102,
        "amplitude_root_count": 128,
    }
    assert aggregation["after"] == {
        "current_count": 378,
        "source_count": 16,
        "interaction_count": 1590,
        "amplitude_root_count": 128,
    }
    assert aggregation["merged_current_count"] == 60
    assert aggregation["deduplicated_interaction_count"] == 512
    assert aggregation["validation"]["accepted"] is True
    assert aggregation["duration_s"] >= 0.0
    assert len(manifest.dag.currents) == 378
    assert len(manifest.dag.interactions) == 1590
    assert len(manifest.dag.amplitude_roots) == 128


def test_zero_current_filter_can_be_disabled_in_artifact(tmp_path: Path) -> None:
    manifest_path, payload = write_generic_dag_process_artifact(
        "d d~ > z g",
        tmp_path,
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    report = raw["generation_filters"]["zero_current"]
    merge_report = raw["generation_filters"]["current_merging"]

    assert payload == raw
    assert report == {
        "duration_s": 0.0,
        "enabled": False,
        "reason": "disabled by --no-numerical-filter-current",
    }
    assert merge_report == {
        "duration_s": 0.0,
        "enabled": False,
        "reason": "disabled by --no-numerical-current-merging",
    }
    assert raw["compiled"]["generic_pruning"]["zero_current_filter"] == report
    assert raw["compiled"]["generic_pruning"]["current_merging"] == merge_report
    plan = json.loads(
        (tmp_path / "generic_process_manifest.json").read_text(encoding="utf-8")
    )
    runtime_schema = plan["runtime_schema"]
    assert len(runtime_schema["current_storage"]["current_slots"]) == len(
        plan["currents"]
    )
    value_storage = runtime_schema["value_storage"]
    assert value_storage["component_count"] >= runtime_schema["current_storage"][
        "component_count"
    ]
    value_slots = value_storage["value_slots"]
    assert {slot["variant"] for slot in value_slots}.issuperset(
        {"source", "propagated", "unpropagated"}
    )
    assert len(runtime_schema["source_fill"]["sources"]) == len(plan["sources"])
    first_source = runtime_schema["source_fill"]["sources"][0]
    assert first_source["source_kind"] == "external-wavefunction"
    assert first_source["crossing"] == "negate-incoming-momentum"
    assert first_source["value_slot"]["variant"] == "source"
    first_stage = runtime_schema["stages"][0]
    assert first_stage["stage_kind"] == "current-combine"
    assert first_stage["interactions"][0]["left_slot"]["current_id"] in (
        first_stage["input_current_ids"]
    )
    first_interaction = first_stage["interactions"][0]
    assert first_interaction["left_value_slot"]["variant"] in {
        "source",
        "propagated",
    }
    assert first_interaction["right_value_slot"]["variant"] in {
        "source",
        "propagated",
    }
    assert first_interaction["result_value_slots"]
    assert set(first_stage["output_value_slot_ids"]).issuperset(
        {
            slot["value_slot_id"]
            for interaction in first_stage["interactions"]
            for slot in interaction["result_value_slots"]
        }
    )
    assert runtime_schema["amplitude_stage"]["output_count"] == len(
        plan["amplitude_roots"]
    )
    first_root = runtime_schema["amplitude_stage"]["roots"][0]
    assert first_root["left_value_slot"]["variant"] in {"source", "unpropagated"}
    assert first_root["right_value_slot"]["variant"] in {"source", "unpropagated"}
    assert plan["currents"][0]["index"]["particle_id"] == -1
    assert "color_state" in plan["currents"][0]["index"]
    assert plan["amplitude_roots"]


def test_numerical_current_merging_supports_opposite_sign_currents() -> None:
    assert _current_signature_relation(
        ((1.0 + 2.0j, -3.0 + 0.5j), (4.0 + 0.0j,)),
        ((-1.0 - 2.0j, 3.0 - 0.5j), (-4.0 + 0.0j,)),
        absolute_tolerance=1.0e-14,
        relative_tolerance=1.0e-14,
    ) == -1.0

    base = build_generic_process_manifest(
        "d d~ > z g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    ).dag
    representative = next(current for current in base.currents if not current.is_source)
    duplicate = CurrentNode(
        id=len(base.currents),
        index=representative.index,
        dimension=representative.dimension,
        is_source=False,
    )
    signed_root = AmplitudeRoot(
        id=len(base.amplitude_roots),
        kind="direct-contraction",
        left_id=duplicate.id,
        right_id=base.sources[0],
        color_weight=(1.0, 0.0),
        contraction="scalar",
    )
    dag = GenericDAG(
        process=base.process,
        color_plan=base.color_plan,
        currents=base.currents + (duplicate,),
        sources=base.sources,
        interactions=base.interactions,
        amplitude_roots=base.amplitude_roots + (signed_root,),
        truncated=base.truncated,
    )

    merged = _merge_currents_in_dag(
        dag,
        {
            current.id: (current.id, 1.0) for current in base.currents
        }
        | {duplicate.id: (representative.id, -1.0)},
    )

    assert len(merged.currents) == len(base.currents)
    assert any(root.color_weight == (-1.0, -0.0) for root in merged.amplitude_roots)


def test_generic_artifact_can_filter_amplitude_stage_to_one_lc_sector(
    tmp_path: Path,
) -> None:
    manifest = build_generic_process_manifest("g g > g g")
    full_payload = manifest.to_json_dict()
    filtered_payload = manifest.to_json_dict(selected_color_sector_ids={0})

    assert full_payload["color_plan"]["sector_count"] == 3
    assert full_payload["runtime_schema"]["amplitude_stage"]["output_count"] == 9
    assert filtered_payload["runtime_schema"]["amplitude_stage"]["output_count"] == 3
    assert filtered_payload["runtime_schema"]["amplitude_stage"][
        "selected_color_sector_ids"
    ] == [0]
    assert [
        root["root_id"]
        for root in filtered_payload["runtime_schema"]["amplitude_stage"]["roots"]
    ] == list(range(3))
    assert {
        root["helicity_weight"]
        for root in filtered_payload["runtime_schema"]["amplitude_stage"]["roots"]
    } == {2.0}
    assert all(
        root["dag_root_id"] >= 0
        for root in filtered_payload["runtime_schema"]["amplitude_stage"]["roots"]
    )

    _, artifact = write_generic_dag_process_artifact(
        manifest,
        tmp_path / "sector0",
        selected_color_sector_ids={0},
    )
    assert artifact["dag_summary"]["amplitude_root_count"] == 3
    assert artifact["compiled"]["selected_color_sector_ids"] == [0]


def test_low_n_pruning_policy_distinguishes_massless_qcd_and_massive_top() -> None:
    pure_gluon = build_generic_process_manifest(
        "g g > g g g g g",
        selected_color_sector_ids={0},
    ).to_json_dict()
    pure_gluon_unfiltered = build_generic_process_manifest(
        "g g > g g g g g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    ).to_json_dict()
    top_pair = build_generic_process_manifest(
        "g g > t t~ g g g",
        selected_color_sector_ids={0},
        numerical_current_merging=False,
    ).to_json_dict()
    top_pair_unfiltered = build_generic_process_manifest(
        "g g > t t~ g g g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    ).to_json_dict()

    assert len(pure_gluon_unfiltered["currents"]) == 310
    assert len(pure_gluon_unfiltered["interactions"]) == 1665
    assert len(pure_gluon["currents"]) == 308
    assert len(pure_gluon["interactions"]) == 1561
    assert len(pure_gluon["amplitude_roots"]) == 56
    assert pure_gluon["generation_filters"]["zero_current"][
        "removed_current_count"
    ] == 2
    assert {
        root["helicity_weight"]
        for root in pure_gluon["runtime_schema"]["amplitude_stage"]["roots"]
    } == {2.0}

    assert len(top_pair_unfiltered["currents"]) == 314
    assert len(top_pair_unfiltered["interactions"]) == 1332
    assert len(top_pair["currents"]) == 310
    assert len(top_pair["interactions"]) == 1232
    assert len(top_pair["amplitude_roots"]) == 128
    assert top_pair["generation_filters"]["zero_current"][
        "removed_current_count"
    ] == 4
    assert {
        root["helicity_weight"]
        for root in top_pair["runtime_schema"]["amplitude_stage"]["roots"]
    } == {1.0}


def test_generic_process_manifest_records_lc_topology_reuse_plan() -> None:
    payload = build_generic_process_manifest("d d~ > z g g").to_json_dict()
    reuse = payload["lc_topology_reuse"]

    assert reuse["available"] is True
    assert reuse["full_color_sector_count"] == 2
    assert reuse["topology_group_count"] == 1
    assert reuse["active_topology_group_count"] == 1
    assert reuse["representative_sector_ids"] == [0]
    assert reuse["groups"][0]["active_sector_ids"] == [0, 1]
    assert reuse["groups"][0]["runtime_reuse_factor"] == 2
    assert len(reuse["groups"][0]["label_permutations"]) == 2


def test_generic_process_manifest_records_line_pairing_representatives() -> None:
    payload = build_generic_process_manifest(
        "d d~ > u u~ s s~ c c~",
        selected_color_sector_ids={0},
        max_coupling_orders={"QCD": 6},
    ).to_json_dict()
    reuse = payload["lc_topology_reuse"]

    assert reuse["full_color_sector_count"] == 24
    assert reuse["line_pairing_representative_sector_count"] == 24
    assert reuse["line_pairing_representative_sector_ids"][:4] == [
        0,
        1,
        2,
        3,
    ]


def test_selected_sector_manifest_keeps_full_plan_and_filters_runtime_dag() -> None:
    payload = build_generic_process_manifest(
        "d d~ > u u~ s s~ c c~ b b~",
        selected_color_sector_ids={0},
        max_coupling_orders={"QCD": 8},
        max_quark_pairs=5,
        numerical_filter_current=False,
        numerical_current_merging=False,
    ).to_json_dict()

    assert payload["color_plan"]["sector_count"] == 120
    assert payload["color_plan"]["truncated"] is False
    assert payload["planning_status"]["color_ready"] is True
    assert payload["lowering_status"]["current_color_sector_count"] == 1
    assert payload["lowering_status"]["current_color_sectors"] == [0]
    assert payload["planning_status"]["generic_evaluator_ready"] is True
    assert payload["lowering_status"]["full_tensor_network_ready"] is True
    assert payload["runtime_schema"]["amplitude_stage"]["output_count"] == 16
    assert {
        root["helicity_weight"]
        for root in payload["runtime_schema"]["amplitude_stage"]["roots"]
    } == {2.0}


def test_reference_color_sector_can_be_selected_before_dag_construction() -> None:
    plan = build_color_plan("d d~ > u u~ s s~")

    assert select_leading_color_sector_ids_from_plan(
        plan,
        reference_color_order=(2, 4, 3, 6, 5, 1),
    ) == {3}
    assert select_leading_color_sector_ids_from_plan(plan) == {0}


def test_generic_process_manifest_can_build_representative_lc_sector() -> None:
    payload = build_generic_process_manifest(
        "d d~ > z g g",
        selected_color_sector_ids={0},
    ).to_json_dict()

    assert payload["color_plan"]["sector_count"] == 2
    assert {current["index"]["color_state"]["sector_id"] for current in payload["currents"]} == {
        0,
    }
    assert payload["lc_topology_reuse"]["representative_sector_ids"] == [0]
    assert payload["lc_topology_reuse"]["groups"][0]["active_sector_ids"] == [0]


def test_generic_dag_artifact_records_lc_topology_replay_for_representative_sector(
    tmp_path: Path,
) -> None:
    _, artifact = write_generic_dag_process_artifact(
        "d d~ > z g g",
        tmp_path,
        selected_color_sector_ids={0},
        lc_topology_replay=True,
    )

    replay = artifact["compiled"]["lc_topology_replay"]
    assert replay["enabled"] is True
    assert replay["materialized_sector_ids"] == [0]
    assert replay["replayed_sector_count"] == 2
    assert replay["groups"][0]["materialized_sector_id"] == 0
    assert replay["groups"][0]["active_sector_ids"] == [0, 1]
    assert len(replay["groups"][0]["sector_permutations"]) == 2
    assert artifact["runtime_lc_topology_replay"] == replay
    assert artifact["dag_summary"]["amplitude_root_count"] == 24


def test_lc_topology_replay_partition_artifact_writes_representative_sidecars(
    tmp_path: Path,
) -> None:
    manifest_path, artifact = write_lc_topology_replay_partition_artifact(
        "d d~ > z g g",
        tmp_path,
    )

    assert manifest_path == tmp_path / "process_manifest.json"
    assert artifact["kind"] == GENERIC_LC_REPLAY_PARTITION_ARTIFACT_KIND
    assert artifact["artifact_class"] == "lc-replay-partition-schema-v2"
    assert artifact["planning_status"]["color_sector_count"] == 2
    assert artifact["lowering_status"]["replayed_color_sector_count"] == 2
    assert artifact["compiled"]["representative_count"] == 1
    representatives = artifact["compiled"]["representative_artifacts"]
    assert len(representatives) == 1
    representative = representatives[0]
    assert representative["representative_sector_id"] == 0
    assert representative["active_sector_ids"] == [0, 1]
    sidecar_path = tmp_path / representative["path"]
    sidecar_manifest = json.loads(
        (sidecar_path / "process_manifest.json").read_text(encoding="utf-8")
    )
    replay = sidecar_manifest["compiled"]["lc_topology_replay"]
    assert replay["enabled"] is True
    assert replay["replayed_sector_count"] == 2
    assert replay["groups"][0]["active_sector_ids"] == [0, 1]
    assert (tmp_path / "validation_momenta.json").exists()


def test_generic_dag_artifact_rejects_unsafe_lc_topology_replay(
    tmp_path: Path,
) -> None:
    _, artifact = write_generic_dag_process_artifact(
        "g g > g g",
        tmp_path,
        selected_color_sector_ids={0},
        lc_topology_replay=True,
    )

    replay = artifact["compiled"]["lc_topology_replay"]
    assert replay["enabled"] is False
    assert replay["reason"] == "no replay-safe topology groups are available"
    assert artifact["lc_topology_reuse"]["replay_safe_representative_sector_ids"] == []
    assert artifact["lc_topology_reuse"]["groups"][0]["runtime_replay_safe"] is False


def test_generic_dag_artifact_keeps_single_sector_selection_without_replay(
    tmp_path: Path,
) -> None:
    _, artifact = write_generic_dag_process_artifact(
        "d d~ > z g g",
        tmp_path,
        selected_color_sector_ids={0},
    )

    assert artifact["compiled"]["lc_topology_replay"]["enabled"] is False
    assert artifact["compiled"]["selected_color_sector_ids"] == [0]
    assert artifact["dag_summary"]["amplitude_root_count"] == 24


def test_generic_dag_artifact_records_generic_pruning_options(
    tmp_path: Path,
) -> None:
    _, artifact = write_generic_dag_process_artifact(
        "d d~ > u u~",
        tmp_path,
        max_coupling_orders={"QED": 2},
        max_lc_current_line_groups=1,
        ignored_particle_ids=(25,),
        ignored_vertex_kinds=(16,),
    )

    pruning = artifact["compiled"]["generic_pruning"]
    assert pruning == {
        "max_coupling_orders": {"QED": 2},
        "max_lc_current_line_groups": 1,
        "max_quark_pairs": None,
        "closure_side_mask_pruning": True,
        "color_order_mask_pruning": True,
        "species_reachability_pruning": True,
        "ignored_particle_ids": [25],
        "ignored_vertex_kinds": [16],
        "structural_current_aggregation": artifact["generation_filters"][
            "structural_current_aggregation"
        ],
        "zero_current_filter": artifact["generation_filters"]["zero_current"],
        "current_merging": artifact["generation_filters"]["current_merging"],
        "reference_color_order": None,
    }


def test_generic_stage_blueprint_defaults_to_global_layout_but_supports_local_inputs() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g g",
        selected_color_sector_ids={0},
    )

    default_blueprint = build_generic_stage_compiler_blueprint(manifest)
    local_blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        stage_local_parameter_layout=True,
    )

    assert {
        stage.parameter_layout for stage in default_blueprint.stages
    } == {"global-value-momentum"}
    assert default_blueprint.amplitude_stage.parameter_layout == "global-value-momentum"
    assert all(not stage.input_components for stage in default_blueprint.stages)

    assert {
        stage.parameter_layout for stage in local_blueprint.stages
    } == {"stage-local-value-momentum"}
    assert local_blueprint.amplitude_stage.parameter_layout == (
        "stage-local-value-momentum"
    )
    assert all(stage.input_components for stage in local_blueprint.stages)
    assert max(stage.parameter_count for stage in local_blueprint.stages) < (
        default_blueprint.parameter_count
    )
    assert local_blueprint.amplitude_stage.parameter_count < (
        default_blueprint.parameter_count
    )
    assert local_blueprint.parameter_count == default_blueprint.parameter_count
    assert local_blueprint.parameter_symbols == ()


def test_generic_parameter_builder_accepts_empty_global_value_storage() -> None:
    builder = _parameter_builder(
        {
            "parameter_layout": {
                "value_component_count": 0,
                "momentum_parameter_count": 8,
                "model_parameter_count": 1,
            }
        }
    )

    assert len(builder.parameter_symbols()) == 9
    assert builder.real_valued_inputs == list(range(9))


def test_generic_stage_blueprint_keeps_four_quark_line_amplitude_outputs() -> None:
    manifest = build_generic_process_manifest(
        "d d~ > u u~ s s~ c c~ g g",
        selected_color_sector_ids={0},
        max_color_sectors=2000,
        max_currents=200000,
        max_quark_pairs=4,
    )

    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        selected_color_sector_ids={0},
        stage_local_parameter_layout=True,
    )

    assert blueprint.expression_ready is True
    assert blueprint.blockers == ()
    assert blueprint.amplitude_stage.output_length == 32
    assert blueprint.amplitude_stage.blockers == ()
    assert blueprint.amplitude_stage.output_expressions


def test_generic_artifact_defaults_to_contributing_lc_sectors(
    tmp_path: Path,
) -> None:
    _, artifact = write_generic_dag_process_artifact(
        "d d~ > u u~ g",
        tmp_path / "runtime",
    )

    assert artifact["compiled"]["selected_color_sector_ids"] == list(range(4))
    assert artifact["dag_summary"]["amplitude_root_count"] == 16
    assert artifact["dag_summary"]["current_count"] == 53
    assert artifact["dag_summary"]["interaction_count"] == 67
    assert artifact["dag_summary"]["source_count"] == 8
    assert artifact["full_dag_summary"]["current_count"] == 53
    assert artifact["full_dag_summary"]["interaction_count"] == 67
    assert artifact["full_dag_summary"]["source_count"] == 8
    aggregation = artifact["generation_filters"]["structural_current_aggregation"]
    assert aggregation["before"]["current_count"] == 61
    assert aggregation["before"]["interaction_count"] == 75
    assert aggregation["merged_current_count"] == 8
    assert aggregation["removed_interaction_count"] == 8
    assert aggregation["validation"]["accepted"] is True
    assert artifact["lowering_status"]["current_color_sectors"] == list(range(4))
    assert artifact["lowering_status"]["internal_current_color_sectors"] == [0]
    runtime_schema = artifact["runtime_schema"]
    assert runtime_schema["source_fill"]["source_count"] == 8
    assert runtime_schema["amplitude_stage"]["output_count"] == 16
    assert {
        root["helicity_weight"]
        for root in runtime_schema["amplitude_stage"]["roots"]
    } == {2.0}


def test_generic_runtime_schema_keeps_auxiliary_tensor_inputs_unpropagated() -> None:
    payload = build_generic_process_manifest("d d~ > z g g g").to_json_dict()
    runtime_schema = payload["runtime_schema"]
    value_slots = runtime_schema["value_storage"]["value_slots"]
    auxiliary_slots = [
        slot for slot in value_slots if slot["particle_id"] == -21
    ]

    assert auxiliary_slots
    assert {slot["variant"] for slot in auxiliary_slots} == {"unpropagated"}
    assert all(slot["applies_propagator"] is False for slot in auxiliary_slots)
    assert {
        slot["propagator"]["kernel"] for slot in auxiliary_slots
    } == {"auxiliary_tensor_embedded_propagator"}

    auxiliary_input_variants = [
        interaction["left_value_slot"]["variant"]
        for stage in runtime_schema["stages"]
        for interaction in stage["interactions"]
        if interaction["left_value_slot"]["current_id"]
        in {slot["current_id"] for slot in auxiliary_slots}
    ] + [
        interaction["right_value_slot"]["variant"]
        for stage in runtime_schema["stages"]
        for interaction in stage["interactions"]
        if interaction["right_value_slot"]["current_id"]
        in {slot["current_id"] for slot in auxiliary_slots}
    ]
    assert auxiliary_input_variants
    assert set(auxiliary_input_variants) == {"unpropagated"}


def test_generic_runtime_schema_keeps_higgsor_auxiliary_scalars_unpropagated() -> None:
    payload = build_generic_process_manifest("d d~ > h h z").to_json_dict()
    runtime_schema = payload["runtime_schema"]
    value_slots = runtime_schema["value_storage"]["value_slots"]
    higgsor_slots = [
        slot for slot in value_slots if slot["particle_id"] in {125, 126, 127}
    ]

    assert higgsor_slots
    assert {slot["variant"] for slot in higgsor_slots} == {"unpropagated"}
    assert all(slot["applies_propagator"] is False for slot in higgsor_slots)
    assert {
        slot["propagator"]["kernel"] for slot in higgsor_slots
    } == {"auxiliary_scalar_no_propagator"}


def test_generic_process_manifest_records_multi_boson_lowering_coverage() -> None:
    payload = build_generic_process_manifest("d d~ > z z g").to_json_dict()

    assert payload["lowering_status"]["has_closure"] is True
    assert payload["lowering_status"]["has_amplitude_roots"] is True
    assert payload["lowering_status"]["unimplemented_vertex_kinds"] == []
    assert payload["lowering_status"]["pending_vertex_kinds"] == []
    assert payload["lowering_status"]["full_tensor_network_ready"] is True
    assert payload["lowering_status"]["ready_vertex_kinds"] == [6, 10]
    assert payload["stage_plan"]["current_stages"][0][
        "unimplemented_vertex_kinds"
    ] == []
    assert payload["stage_plan"]["amplitude_stage"]["closure_count"] > 0


def test_generic_process_manifest_handles_multi_quark_line_plan() -> None:
    payload = build_generic_process_manifest("d d~ > u u~").to_json_dict()

    assert payload["process_ir"]["quark_lines"]["quark_pair_count"] == 2
    assert payload["planning_status"]["color_sector_count"] == 2
    assert payload["planning_status"]["color_ready"] is True
    assert payload["lowering_status"]["current_color_sector_count"] == 2
    assert payload["lowering_status"]["current_color_sectors"] == [0, 1]
    assert len(payload["lowering_status"]["color_sector_summaries"]) == 2
    assert payload["lowering_status"]["color_sector_summaries"][0][
        "pending_vertex_kinds"
    ] == []
    assert payload["lowering_status"]["unimplemented_vertex_kinds"] == []
    assert payload["lowering_status"]["pending_vertex_kinds"] == []
    assert payload["lowering_status"]["has_closure"] is True
    assert payload["lowering_status"]["full_tensor_network_ready"] is True
    runtime_schema = payload["runtime_schema"]
    assert runtime_schema["normalization"]["color_factor"] == 9
    assert runtime_schema["source_fill"]["source_count"] == len(payload["sources"])
    assert runtime_schema["value_storage"]["component_count"] > 0
    assert runtime_schema["amplitude_stage"]["output_count"] == len(
        payload["amplitude_roots"]
    )
    source_sectors = {
        source["color_state"]["sector_id"]
        for source in runtime_schema["source_fill"]["sources"]
    }
    assert source_sectors == {0}
    assert payload["stage_plan"]["current_stages"][-1]["subset_size"] == 3
    assert payload["stage_plan"]["amplitude_stage"]["closure_count"] > 0


def test_pure_gluon_lc_normalization_uses_model_leading_color_factor_once() -> None:
    payload = build_generic_process_manifest("g g > g g").to_json_dict()

    assert payload["runtime_schema"]["normalization"]["color_factor"] == 81


@pytest.mark.parametrize("color_accuracy", ["nlc", "full"])
def test_pure_gluon_subleading_colour_shares_currents_but_keeps_root_sectors(
    color_accuracy: str,
) -> None:
    payload = build_generic_process_manifest(
        "g g > g g g",
        color_accuracy=color_accuracy,
    ).to_json_dict()

    amplitude_stage = payload["runtime_schema"]["amplitude_stage"]
    color_contraction = amplitude_stage["color_contraction"]

    assert payload["color_plan"]["sector_count"] == 24
    assert payload["lowering_status"]["internal_current_color_sectors"] == [0]
    assert payload["lowering_status"]["current_color_sectors"] == list(range(24))
    assert {
        root["color_sector_id"]
        for root in amplitude_stage["roots"]
    } == set(range(24))
    assert color_contraction["supported"] is True
    assert color_contraction["group_count"] == amplitude_stage["output_count"]
    assert color_contraction["entry_count"] > color_contraction["group_count"]


def test_generic_process_manifest_does_not_divide_me_by_quark_line_partners() -> None:
    payload = build_generic_process_manifest("d d~ > u u~ s s~").to_json_dict()
    normalization = payload["runtime_schema"]["normalization"]

    assert payload["process_ir"]["quark_lines"]["quark_pair_count"] == 3
    assert normalization["final_state_identical_factor"] == 1
    assert normalization["quark_line_partner_factor"] == 1
    assert normalization["identical_factor"] == 1
    assert payload["lowering_status"]["has_closure"] is True
    assert payload["runtime_schema"]["amplitude_stage"]["output_count"] > 0


def test_generic_process_manifest_write_and_load_roundtrip(tmp_path: Path) -> None:
    manifest = build_generic_process_manifest("d d~ > e+ e- g")
    manifest_path = write_generic_process_manifest(manifest, tmp_path)

    loaded = load_generic_process_manifest(manifest_path)
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert loaded == raw
    assert loaded["process"] == "d d~ > e+ e- g"
    assert loaded["process_ir"]["labels"]["leptons"] == [3, 4]
    assert (tmp_path / "generic_process_manifest.json").exists()


def test_generic_process_set_manifest_writes_nested_subprocesses(
    tmp_path: Path,
) -> None:
    manifest = build_generic_process_set_manifest(
        "d d~ > z g | d d~ > z z g",
        selected_color_sector_ids={0},
        max_coupling_orders={"QED": 2},
        max_quark_pairs=1,
        species_reachability_pruning=False,
    )
    manifest_path = write_generic_process_set_manifest(manifest, tmp_path)

    loaded = load_generic_process_set_manifest(manifest_path)

    assert loaded["schema_version"] == GENERIC_PROCESS_SCHEMA_VERSION
    assert loaded["kind"] == GENERIC_PROCESS_SET_MANIFEST_KIND
    assert loaded["request"] == "d d~ > z g | d d~ > z z g"
    assert loaded["default_process_key"] == "d_dbar_to_z_g"
    assert loaded["generic_generation"]["pruning"]["selected_color_sector_ids"] == [0]
    assert loaded["generic_generation"]["pruning"]["max_coupling_orders"] == {
        "QED": 2
    }
    assert loaded["generic_generation"]["pruning"]["max_quark_pairs"] == 1
    assert loaded["generic_generation"]["pruning"]["species_reachability_pruning"] is False
    assert [entry["key"] for entry in loaded["processes"]] == [
        "d_dbar_to_z_g",
        "d_dbar_to_z_z_g",
    ]
    assert [entry["generation_request"] for entry in loaded["processes"]] == [
        loaded["generic_generation"],
        loaded["generic_generation"],
    ]
    assert (
        tmp_path
        / "subprocesses"
        / "d_dbar_to_z_g"
        / "generic_process_manifest.json"
    ).exists()
    assert loaded["processes"][1]["lowering_status"][
        "unimplemented_vertex_kinds"
    ] == []
    assert loaded["processes"][1]["planning_status"]["current_ready"] is True
    assert loaded["processes"][0]["planning_status"]["color_sector_count"] == 1
    assert loaded["processes"][0]["planning_status"]["color_ready"] is True
    nested = json.loads(
        (
            tmp_path
            / "subprocesses"
            / "d_dbar_to_z_g"
            / "generic_process_manifest.json"
        ).read_text(encoding="utf-8")
    )
    assert nested["color_plan"]["sector_count"] == 1


def test_generic_dag_process_artifact_writes_schema_v2_stage_blueprint_runtime(
    tmp_path: Path,
) -> None:
    manifest_path, payload = write_generic_dag_process_artifact(
        "d d~ > z g",
        tmp_path,
        evaluator_backend="jit",
        compiled_preset="generation",
        batch_size=16,
    )

    raw = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert payload == raw
    assert raw["schema_version"] == GENERIC_PROCESS_SCHEMA_VERSION
    assert raw["kind"] == GENERIC_DAG_PROCESS_ARTIFACT_KIND
    assert raw["artifact_class"] == "generic-dag-schema-v2"
    assert raw["compiled"]["kind"] == "generic-dag-stage-blueprint"
    assert raw["compiled"]["runtime_available"] is False
    assert raw["compiled"]["requested_evaluator_backend"] == "jit"
    assert raw["compiled"]["batch_size"] == 16
    assert raw["compiled"]["stage_compiler"]["runtime_available"] is False
    assert raw["compiled"]["stage_compiler"]["expression_ready"] is True
    assert raw["compiled"]["stage_compiler"]["stage_count"] == (
        len(raw["runtime_schema"]["stages"]) + 1
    )
    assert raw["compiled"]["stage_compiler"]["stages"][0]["output_slots"]
    assert raw["compiled"]["stage_compiler"]["amplitude_stage"]["output_length"] == len(
        raw["runtime_schema"]["amplitude_stage"]["roots"]
    )
    stage_compiler = raw["compiled"]["stage_compiler"]
    assert stage_compiler["real_valued_inputs"][0] == stage_compiler[
        "value_parameter_count"
    ]
    assert len(stage_compiler["real_valued_inputs"]) == (
        stage_compiler["momentum_parameter_count"]
        + stage_compiler["model_parameter_count"]
    )
    assert raw["planning_status"]["generic_evaluator_ready"] is True
    assert raw["dag_summary"]["current_count"] > 0
    assert raw["runtime_schema"]["kind"] == "pyamplicol-generic-dag-runtime-schema"
    assert raw["normalization"] == raw["runtime_schema"]["normalization"]
    assert raw["runtime_schema"]["parameter_layout"]["momentum_components_real"] is True
    assert raw["runtime_schema"]["amplitude_stage"]["roots"]
    assert (tmp_path / "generic_process_manifest.json").exists()
    assert (tmp_path / "check_standalone.py").exists()
    assert (tmp_path / "validation_momenta.json").exists()


def test_generic_dag_process_artifact_default_backend_is_cpp_o3(
    tmp_path: Path,
) -> None:
    def fake_compiler(stage, params, real_inputs):
        return {
            "kind": "compiled-complex-evaluator",
            "backend": "compiled-complex",
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "number_type": "f64",
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
            "source_path": f"compiled/{stage.evaluator_label}.cpp",
            "library_path": f"compiled/lib{stage.evaluator_label}",
            "function_name": stage.evaluator_label,
        }

    manifest_path, payload = write_generic_dag_process_artifact(
        "d d~ > z",
        tmp_path,
        emit_stage_evaluator_artifacts=True,
        stage_evaluator_compiler=fake_compiler,
    )
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert payload == raw
    assert raw["compiled"]["requested_evaluator_backend"] == "compiled-complex"
    assert raw["compiled"]["requested_compiled_preset"] == "runtime-o3"
    assert raw["compiled"]["runtime_available"] is True
    assert raw["compiled"]["stage_evaluators"]["runtime_available"] is True
    assert raw["compiled"]["stage_evaluators"]["stages"][0]["evaluator"][
        "kind"
    ] == "compiled-complex-evaluator"


def test_generic_dag_process_artifact_writes_massive_one_body_validation_point(
    tmp_path: Path,
) -> None:
    write_generic_dag_process_artifact("d d~ > z", tmp_path)

    payload = json.loads((tmp_path / "validation_momenta.json").read_text())
    point = payload["points"][0]

    assert payload["kind"] == "pyamplicol-rusticol-validation-momenta"
    assert payload["available"] is True
    assert [particle["pdg"] for particle in point] == [1, -1, 23]
    _assert_validation_point_is_conserved_and_on_shell(point)
    assert point[2]["momentum"][1:] == ["0", "0", "0"]


def test_generic_dag_process_artifact_writes_generic_rambo_validation_point(
    tmp_path: Path,
) -> None:
    write_generic_dag_process_artifact("d d~ > e+ e- g", tmp_path)

    payload = json.loads((tmp_path / "validation_momenta.json").read_text())
    point = payload["points"][0]

    assert payload["available"] is True
    assert [particle["pdg"] for particle in point] == [1, -1, -11, 11, 21]
    _assert_validation_point_is_conserved_and_on_shell(point)


def test_generic_dag_process_artifact_records_unavailable_validation_point(
    tmp_path: Path,
) -> None:
    manifest_path, raw = write_generic_dag_process_artifact("d d~ > a", tmp_path)
    validation = json.loads((tmp_path / "validation_momenta.json").read_text())

    assert manifest_path.exists()
    assert raw["process"] == "d d~ > a"
    assert validation["available"] is False
    assert validation["points"] == []
    assert "one massless final state" in validation["error"]


def test_generic_stage_compiler_keeps_expressions_out_of_json() -> None:
    manifest = build_generic_process_manifest("d d~ > z g")
    blueprint = build_generic_stage_compiler_blueprint(manifest)
    payload = blueprint.to_json_dict()

    assert blueprint.expression_ready is True
    assert len(blueprint.parameter_symbols) == blueprint.parameter_count
    assert blueprint.stages
    assert blueprint.stages[0].output_expressions
    assert len(blueprint.stages[0].output_expressions) == blueprint.stages[0].output_length
    assert blueprint.amplitude_stage.output_expressions
    assert len(blueprint.amplitude_stage.output_expressions) == (
        blueprint.amplitude_stage.output_length
    )
    assert blueprint.stages[0].evaluator_label.startswith("generic_stage_")
    assert blueprint.stages[0].parameter_layout == "global-value-momentum"
    assert blueprint.amplitude_stage.evaluator_label == "generic_amplitude_stage"
    assert "parameter_symbols" not in payload
    assert "output_expressions" not in payload["stages"][0]
    assert "output_expressions" not in payload["amplitude_stage"]
    assert payload["stages"][0]["evaluator_label"] == blueprint.stages[0].evaluator_label
    assert payload["stages"][0]["parameter_layout"] == "global-value-momentum"
    assert payload["amplitude_stage"]["evaluator_label"] == "generic_amplitude_stage"
    assert payload["stages"][0]["first_output_previews"]


def _assert_validation_point_is_conserved_and_on_shell(
    point: list[dict[str, object]],
) -> None:
    model = AmplicolSMLeadingColorModel()
    momenta = [
        tuple(float(component) for component in particle["momentum"])
        for particle in point
    ]
    initial = _sum_momenta(momenta[:2])
    final = _sum_momenta(momenta[2:])
    for observed, expected in zip(final, initial, strict=True):
        assert math.isclose(observed, expected, rel_tol=0.0, abs_tol=1.0e-9)
    for particle, momentum in zip(point, momenta, strict=True):
        mass = model.mass(int(particle["pdg"]))
        observed_m2 = _minkowski_square(momentum)
        expected_m2 = mass * mass
        assert math.isclose(
            observed_m2,
            expected_m2,
            rel_tol=1.0e-9,
            abs_tol=1.0e-7,
        )


def _sum_momenta(
    momenta: list[tuple[float, float, float, float]],
) -> tuple[float, float, float, float]:
    return (
        sum(momentum[0] for momentum in momenta),
        sum(momentum[1] for momentum in momenta),
        sum(momentum[2] for momentum in momenta),
        sum(momentum[3] for momentum in momenta),
    )


def _minkowski_square(momentum: tuple[float, float, float, float]) -> float:
    return (
        momentum[0] * momentum[0]
        - momentum[1] * momentum[1]
        - momentum[2] * momentum[2]
        - momentum[3] * momentum[3]
    )


def test_generic_stage_evaluator_artifact_writer_uses_stage_expressions(
    tmp_path: Path,
) -> None:
    manifest = build_generic_process_manifest("u d~ > e+ ve g")
    blueprint = build_generic_stage_compiler_blueprint(manifest)
    calls: list[dict[str, object]] = []

    def fake_compiler(stage, params, real_inputs):
        calls.append(
            {
                "label": stage.evaluator_label,
                "output_count": len(stage.output_expressions),
                "param_count": len(params),
                "real_count": len(real_inputs),
            }
        )
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    payload = write_generic_stage_evaluator_artifacts(
        blueprint,
        tmp_path / "generic-evaluators",
        compiler=fake_compiler,
    )

    assert payload["kind"] == "generic-dag-stage-evaluator-artifacts"
    assert payload["runtime_available"] is True
    assert payload["runtime_unavailable_message"] is None
    assert payload["parameter_layout"] == "global-value-momentum"
    assert payload["stage_count"] == len(blueprint.stages) + 1
    assert len(payload["stages"]) == len(blueprint.stages)
    assert payload["stages"][0]["evaluator"]["kind"] == "jit-symbolica-evaluator"
    assert payload["stages"][0]["evaluator"]["output_len"] == (
        blueprint.stages[0].output_length
    )
    assert payload["amplitude_stage"]["evaluator"]["label"] == (
        "generic_amplitude_stage"
    )
    assert "output_expressions" not in payload["stages"][0]
    assert "output_expressions" not in payload["amplitude_stage"]
    assert [call["label"] for call in calls] == [
        *(stage.evaluator_label for stage in blueprint.stages),
        blueprint.amplitude_stage.evaluator_label,
    ]
    assert all(call["param_count"] == blueprint.parameter_count for call in calls)
    assert all(
        call["real_count"]
        == blueprint.momentum_parameter_count + blueprint.model_parameter_count
        for call in calls
    )


def test_streaming_stage_evaluator_writer_releases_consumed_expressions(
    tmp_path: Path,
) -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )
    runtime_schema = _generic_runtime_schema_payload(
        manifest.dag,
        manifest.model,
    )
    calls: list[tuple[str, int]] = []

    def fake_compiler(stage, params, real_inputs):
        del params, real_inputs
        calls.append((stage.evaluator_label, len(stage.output_expressions)))
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": stage.parameter_count,
            "output_len": len(stage.output_expressions),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    blueprint, payload = build_and_write_generic_stage_evaluator_artifacts(
        manifest,
        runtime_schema,
        tmp_path / "streamed-evaluators",
        stage_local_parameter_layout=True,
        compiler=fake_compiler,
    )

    assert all(output_count > 0 for _label, output_count in calls)
    assert [label for label, _output_count in calls] == [
        *(stage.evaluator_label for stage in blueprint.stages),
        blueprint.amplitude_stage.evaluator_label,
    ]
    assert all(not stage.output_expressions for stage in blueprint.stages)
    assert not blueprint.amplitude_stage.output_expressions
    assert len(payload["stages"]) == len(blueprint.stages)
    assert payload["amplitude_stage"]["evaluator"]["label"] == (
        blueprint.amplitude_stage.evaluator_label
    )


def test_generic_stage_evaluator_writer_accepts_four_quark_line_amplitude_stage(
    tmp_path: Path,
) -> None:
    manifest = build_generic_process_manifest(
        "d d~ > u u~ s s~ c c~ g g",
        selected_color_sector_ids={0},
        max_color_sectors=2000,
        max_currents=200000,
        max_quark_pairs=4,
    )
    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        selected_color_sector_ids={0},
        stage_local_parameter_layout=True,
    )
    calls: list[dict[str, object]] = []

    def fake_compiler(stage, params, real_inputs):
        calls.append(
            {
                "label": stage.evaluator_label,
                "output_count": len(stage.output_expressions),
                "param_count": len(params),
                "real_count": len(real_inputs),
            }
        )
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    payload = write_generic_stage_evaluator_artifacts(
        blueprint,
        tmp_path / "four-quark-line-evaluators",
        compiler=fake_compiler,
    )

    assert payload["runtime_available"] is True
    assert payload["amplitude_stage"]["evaluator"]["label"] == (
        "generic_amplitude_stage"
    )
    assert payload["amplitude_stage"]["evaluator"]["output_len"] == 32
    assert [call["label"] for call in calls][-1] == "generic_amplitude_stage"
    assert [call["output_count"] for call in calls][-1] == 32
    assert all(call["output_count"] > 0 for call in calls)
    assert all(call["real_count"] <= call["param_count"] for call in calls)
    assert calls[-1]["real_count"] == 0


def test_generic_stage_evaluator_artifact_writer_default_symbolica_bridge(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pyamplicol import symbolica_evaluator

    manifest = build_generic_process_manifest("d d~ > z g")
    blueprint = build_generic_stage_compiler_blueprint(manifest)
    calls: list[dict[str, object]] = []

    def fake_compile_symbolica_outputs(
        outputs,
        params,
        *,
        merge_evaluators_strategy,
        verbose_evaluator_build,
        real_params,
        symbolica_settings,
        jit_compile,
        label,
        progress_callback,
        functions,
    ):
        del symbolica_settings, progress_callback, functions
        calls.append(
            {
                "label": label,
                "output_count": len(outputs),
                "param_count": len(params),
                "real_params": tuple(real_params),
                "merge": merge_evaluators_strategy,
                "verbose": verbose_evaluator_build,
                "jit": jit_compile,
            }
        )
        return {"label": label, "output_count": len(outputs)}

    def fake_artifact_manifest(evaluator, artifact_dir):
        return {
            "kind": "stubbed-symbolica-evaluator",
            "label": evaluator["label"],
            "input_len": blueprint.parameter_count,
            "output_len": evaluator["output_count"],
            "artifact_dir": str(artifact_dir),
        }

    monkeypatch.setattr(
        symbolica_evaluator,
        "_compile_symbolica_outputs",
        fake_compile_symbolica_outputs,
    )
    monkeypatch.setattr(
        symbolica_evaluator,
        "_symbolica_evaluator_artifact_manifest",
        fake_artifact_manifest,
    )

    payload = write_generic_stage_evaluator_artifacts(
        blueprint,
        tmp_path / "generic-default-writer",
        merge_evaluators_strategy=True,
        verbose_evaluator_build=True,
        jit_compile=False,
    )

    assert payload["runtime_available"] is True
    assert payload["stages"][0]["evaluator"]["kind"] == "stubbed-symbolica-evaluator"
    assert payload["amplitude_stage"]["evaluator"]["label"] == (
        "generic_amplitude_stage"
    )
    assert [call["label"] for call in calls] == [
        *(stage.evaluator_label for stage in blueprint.stages),
        blueprint.amplitude_stage.evaluator_label,
    ]
    assert all(call["param_count"] == blueprint.parameter_count for call in calls)
    assert all(
        call["real_params"] == blueprint.real_valued_inputs for call in calls
    )
    assert all(call["merge"] is True for call in calls)
    assert all(call["verbose"] is True for call in calls)
    assert all(call["jit"] is False for call in calls)


def test_generic_dag_process_artifact_can_embed_stage_evaluator_manifests(
    tmp_path: Path,
) -> None:
    def fake_compiler(stage, params, real_inputs):
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "real_input_count": len(real_inputs),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    manifest_path, payload = write_generic_dag_process_artifact(
        "u d~ > e+ ve g",
        tmp_path,
        emit_stage_evaluator_artifacts=True,
        stage_evaluator_compiler=fake_compiler,
    )
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    stage_evaluators = raw["compiled"]["stage_evaluators"]

    assert payload == raw
    assert raw["compiled"]["kind"] == "generic-dag-stage-blueprint"
    assert raw["compiled"]["runtime_available"] is True
    assert raw["compiled"]["runtime_unavailable_message"] is None
    assert stage_evaluators["kind"] == "generic-dag-stage-evaluator-artifacts"
    assert stage_evaluators["runtime_available"] is True
    assert stage_evaluators["runtime_unavailable_message"] is None
    assert stage_evaluators["stages"][0]["evaluator"]["kind"] == (
        "jit-symbolica-evaluator"
    )
    assert stage_evaluators["stages"][0]["evaluator"]["input_len"] == (
        stage_evaluators["parameter_count"]
    )
    assert stage_evaluators["amplitude_stage"]["evaluator"]["label"] == (
        "generic_amplitude_stage"
    )
    assert str(
        stage_evaluators["amplitude_stage"]["evaluator"]["evaluator_state_path"]
    ).startswith("evaluators/")
    runtime_stage = raw["runtime_schema"]["stages"][0]
    assert raw["runtime_schema"]["serialized_stage_metadata_compacted"] is True
    assert runtime_stage["interactions_compacted"] is True
    assert runtime_stage["interactions"] == []
    assert len(runtime_stage["interaction_ids"]) == runtime_stage["interaction_count"]
    assert sum(
        count for _vertex_kind, count in runtime_stage["interaction_vertex_kind_counts"]
    ) == runtime_stage["interaction_count"]
    stage_compiler = raw["compiled"]["stage_compiler"]
    assert stage_compiler["serialized_metadata_compacted"] is True
    assert "stages" not in stage_compiler
    assert "amplitude_stage" not in stage_compiler
    diagnostics = raw["compiled"]["diagnostic_runtime_schema"]
    diagnostics_path = tmp_path / diagnostics["path"]
    assert diagnostics["compression"] == "gzip"
    assert diagnostics["format"] == (
        "python-pickle-protocol-5-full-runtime-schema-v2"
    )
    assert diagnostics["interaction_count"] == raw["dag_summary"][
        "interaction_count"
    ]
    assert diagnostics["size_bytes"] == diagnostics_path.stat().st_size
    with gzip.open(diagnostics_path, mode="rb") as stream:
        detailed_runtime_schema = pickle.load(stream)
    assert detailed_runtime_schema["stages"][0]["interactions"]


def test_large_runtime_schema_uses_compact_storage_metadata(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_compiler(stage, params, real_inputs):
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "real_input_count": len(real_inputs),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    monkeypatch.setattr(
        generic_artifact_module,
        "_COMPACT_RUNTIME_STORAGE_MIN_CURRENTS",
        0,
    )
    manifest_path, _payload = write_generic_dag_process_artifact(
        "d d~ > z g g",
        tmp_path,
        emit_stage_evaluator_artifacts=True,
        stage_evaluator_compiler=fake_compiler,
    )

    runtime_schema = json.loads(manifest_path.read_text())["runtime_schema"]
    diagnostics = json.loads(manifest_path.read_text())["compiled"][
        "diagnostic_runtime_schema"
    ]
    current_storage = runtime_schema["current_storage"]
    value_storage = runtime_schema["value_storage"]
    runtime_stages = runtime_schema["stages"]
    assert current_storage["metadata_compacted"] is True
    assert value_storage["metadata_compacted"] is True
    assert "color_state" not in current_storage["current_slots"][0]
    assert "helicity_ancestry" not in current_storage["current_slots"][0]
    assert "propagator" not in value_storage["value_slots"][0]
    assert "external_labels" not in value_storage["value_slots"][0]
    assert all(stage["interactions_compacted"] is True for stage in runtime_stages)
    assert all(stage["interactions"] == [] for stage in runtime_stages)
    assert all(
        len(stage["interaction_ids"]) == stage["interaction_count"]
        for stage in runtime_stages
    )
    assert diagnostics["format"] == (
        "python-pickle-protocol-5-compact-runtime-schema-v2"
    )


def test_compact_interactions_preserve_stage_compiler_expressions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manifest = build_generic_process_manifest(
        "d d~ > z g g",
        selected_color_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )

    monkeypatch.setattr(
        generic_artifact_module,
        "_COMPACT_RUNTIME_STORAGE_MIN_CURRENTS",
        len(manifest.dag.currents) + 1,
    )
    verbose = build_generic_stage_compiler_blueprint(
        manifest,
        stage_local_parameter_layout=True,
    )
    monkeypatch.setattr(
        generic_artifact_module,
        "_COMPACT_RUNTIME_STORAGE_MIN_CURRENTS",
        0,
    )
    compact = build_generic_stage_compiler_blueprint(
        manifest,
        stage_local_parameter_layout=True,
    )

    assert compact.to_json_dict() == verbose.to_json_dict()
    for compact_stage, verbose_stage in zip(
        (*compact.stages, compact.amplitude_stage),
        (*verbose.stages, verbose.amplitude_stage),
        strict=True,
    ):
        assert tuple(map(str, compact_stage.output_expressions)) == tuple(
            map(str, verbose_stage.output_expressions)
        )


def test_large_runtime_schema_diagnostics_overlap_stage_compilation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    writer_threads: list[str] = []
    original_writer = generic_artifact_module._write_runtime_schema_diagnostics

    def recording_writer(runtime_schema, output_dir):
        writer_threads.append(threading.current_thread().name)
        return original_writer(runtime_schema, output_dir)

    def fake_compiler(stage, params, real_inputs):
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "real_input_count": len(real_inputs),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    monkeypatch.setattr(
        generic_artifact_module,
        "_ASYNC_RUNTIME_SCHEMA_DIAGNOSTICS_MIN_INTERACTIONS",
        0,
    )
    monkeypatch.setattr(
        generic_artifact_module,
        "_ASYNC_RUNTIME_SCHEMA_DIAGNOSTICS_TRIGGER_STAGE_S",
        0.0,
    )
    monkeypatch.setattr(
        generic_artifact_module,
        "_write_runtime_schema_diagnostics",
        recording_writer,
    )

    manifest_path, _payload = write_generic_dag_process_artifact(
        "u d~ > e+ ve g",
        tmp_path,
        emit_stage_evaluator_artifacts=True,
        stage_evaluator_compiler=fake_compiler,
    )

    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    diagnostics_path = tmp_path / raw["compiled"]["diagnostic_runtime_schema"]["path"]
    assert diagnostics_path.exists()
    assert writer_threads == ["pyamplicol-schema-diagnostics_0"]


def test_generic_dag_process_artifact_can_omit_lc_runtime_selector(
    tmp_path: Path,
) -> None:
    def fake_compiler(stage, params, real_inputs):
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "real_input_count": len(real_inputs),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    manifest_path, _ = write_generic_dag_process_artifact(
        "d d~ > z g g",
        tmp_path,
        emit_stage_evaluator_artifacts=True,
        enable_lc_sector_runtime_selector=False,
        stage_evaluator_compiler=fake_compiler,
        numerical_filter_current=False,
        numerical_current_merging=False,
    )

    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert raw["compiled"]["runtime_available"] is True
    assert not any(
        parameter["name"] == "runtime.lc_sector_id"
        for parameter in raw["runtime_schema"]["model_parameters"]
    )


def test_generic_dag_process_artifact_writes_compact_main_with_lc_runtime_sidecar(
    tmp_path: Path,
) -> None:
    def fake_compiler(stage, params, real_inputs):
        return {
            "kind": "jit-symbolica-evaluator",
            "label": stage.evaluator_label,
            "input_len": len(params),
            "output_len": len(stage.output_expressions),
            "real_input_count": len(real_inputs),
            "evaluator_state_path": f"evaluators/{stage.evaluator_label}.bin",
        }

    manifest_path, payload = write_generic_dag_process_artifact(
        "d d~ > z g g",
        tmp_path,
        emit_stage_evaluator_artifacts=True,
        skip_main_stage_evaluator_artifacts=True,
        stage_evaluator_compiler=fake_compiler,
        runtime_lc_sector_ids={0},
        numerical_filter_current=False,
        numerical_current_merging=False,
    )

    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    compiled = raw["compiled"]
    sidecars = compiled["runtime_lc_sector_artifacts"]
    sidecar_path = tmp_path / sidecars[0]["path"]
    sidecar = json.loads(
        (sidecar_path / "process_manifest.json").read_text(encoding="utf-8")
    )

    assert payload == raw
    assert raw["generic_plan_path"] is None
    assert "runtime_schema" not in raw
    assert not (tmp_path / "generic_process_manifest.json").exists()
    assert compiled["compact_main_artifact"] is True
    assert compiled["runtime_available"] is False
    assert raw["lowering_status"]["color_sector_summaries"] == []
    assert raw["lowering_status"]["color_sector_summaries_omitted"] is True
    assert sidecars == [
        {
            "color_sector_ids": [0],
            "path": "runtime_lc_sectors/lc_0",
            "kind": "selected-lc-runtime-sidecar",
            "runtime_available": True,
            "runtime_selector": "none-specialized-runtime-artifact",
        }
    ]
    assert sidecar["generic_plan_path"] == "generic_process_manifest.json"
    assert sidecar["compiled"]["runtime_available"] is True
    assert "runtime_schema" in sidecar
    assert not any(
        parameter["name"] == "runtime.lc_sector_id"
        for parameter in sidecar["runtime_schema"]["model_parameters"]
    )
    assert (sidecar_path / "generic_process_manifest.json").exists()


def test_generic_dag_process_set_artifact_writes_subprocess_manifests(
    tmp_path: Path,
) -> None:
    process_set = build_process_set_ir("d d~ > z g | u u~ > z g")
    manifest_path, payload = write_generic_dag_process_set_artifact(
        process_set,
        tmp_path,
        selected_color_sector_ids={0},
        max_coupling_orders={"QED": 1},
        max_quark_pairs=2,
        species_reachability_pruning=False,
    )

    raw = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert payload == raw
    assert raw["schema_version"] == GENERIC_PROCESS_SCHEMA_VERSION
    assert raw["kind"] == GENERIC_DAG_PROCESS_SET_ARTIFACT_KIND
    assert raw["artifact_class"] == "generic-dag-schema-v2"
    assert raw["runtime_available"] is False
    assert raw["default_process_key"] == "d_dbar_to_z_g"
    assert raw["generic_generation"]["pruning"]["selected_color_sector_ids"] == [0]
    assert raw["generic_generation"]["pruning"]["max_coupling_orders"] == {"QED": 1}
    assert raw["generic_generation"]["pruning"]["max_quark_pairs"] == 2
    assert raw["generic_generation"]["pruning"]["species_reachability_pruning"] is False
    assert [entry["kind"] for entry in raw["processes"]] == [
        GENERIC_DAG_PROCESS_ARTIFACT_KIND,
        GENERIC_DAG_PROCESS_ARTIFACT_KIND,
    ]
    assert [entry["generation_request"] for entry in raw["processes"]] == [
        raw["generic_generation"],
        raw["generic_generation"],
    ]
    assert [entry["artifact_class"] for entry in raw["processes"]] == [
        "generic-dag-schema-v2",
        "generic-dag-schema-v2",
    ]
    assert (
        tmp_path
        / "subprocesses"
        / "d_dbar_to_z_g"
        / "process_manifest.json"
    ).exists()


def test_generic_process_manifest_loader_rejects_wrong_kind(tmp_path: Path) -> None:
    path = tmp_path / "generic_process_manifest.json"
    path.write_text(
        json.dumps({"schema_version": GENERIC_PROCESS_SCHEMA_VERSION, "kind": "other"}),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="unsupported generic process manifest kind"):
        load_generic_process_manifest(path)


def test_generic_process_set_manifest_loader_rejects_wrong_kind(
    tmp_path: Path,
) -> None:
    path = tmp_path / "generic_process_set_manifest.json"
    path.write_text(
        json.dumps(
            {"schema_version": GENERIC_PROCESS_SCHEMA_VERSION, "kind": "other"}
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match="unsupported generic process-set manifest kind",
    ):
        load_generic_process_set_manifest(path)
