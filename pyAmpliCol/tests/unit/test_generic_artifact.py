from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from pyamplicol.generic_artifact import (
    GENERIC_DAG_PROCESS_ARTIFACT_KIND,
    GENERIC_DAG_PROCESS_SET_ARTIFACT_KIND,
    GENERIC_PROCESS_MANIFEST_KIND,
    GENERIC_PROCESS_SET_MANIFEST_KIND,
    GENERIC_PROCESS_SCHEMA_VERSION,
    build_generic_process_set_manifest,
    build_generic_process_manifest,
    load_generic_process_manifest,
    load_generic_process_set_manifest,
    write_generic_dag_process_artifact,
    write_generic_dag_process_set_artifact,
    write_generic_process_set_manifest,
    write_generic_process_manifest,
    select_leading_color_sector_ids_from_plan,
    _json_safe_bigints,
)
from pyamplicol.color_plan import build_color_plan
from pyamplicol.generic_stage_compiler import (
    build_generic_stage_compiler_blueprint,
    write_generic_stage_evaluator_artifacts,
)
from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.process_ir import build_process_set_ir


def test_json_safe_bigints_serializes_large_bitsets_as_hex_strings() -> None:
    payload = _json_safe_bigints(
        {
            "small": 42,
            "large": 1 << 20000,
            "nested": [True, 1 << 5000],
        }
    )

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
    assert len(runtime_schema["current_storage"]["current_slots"]) == len(
        payload["currents"]
    )
    value_storage = runtime_schema["value_storage"]
    assert value_storage["component_count"] >= runtime_schema["current_storage"][
        "component_count"
    ]
    value_slots = value_storage["value_slots"]
    assert {slot["variant"] for slot in value_slots}.issuperset(
        {"source", "propagated", "unpropagated"}
    )
    assert len(runtime_schema["source_fill"]["sources"]) == len(payload["sources"])
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
        payload["amplitude_roots"]
    )
    first_root = runtime_schema["amplitude_stage"]["roots"][0]
    assert first_root["left_value_slot"]["variant"] in {"source", "unpropagated"}
    assert first_root["right_value_slot"]["variant"] in {"source", "unpropagated"}
    assert payload["currents"][0]["index"]["particle_id"] == -1
    assert "color_state" in payload["currents"][0]["index"]
    assert payload["amplitude_roots"]


def test_generic_artifact_can_filter_amplitude_stage_to_one_lc_sector(
    tmp_path: Path,
) -> None:
    manifest = build_generic_process_manifest("g g > g g")
    full_payload = manifest.to_json_dict()
    filtered_payload = manifest.to_json_dict(selected_color_sector_ids={0})

    assert full_payload["color_plan"]["sector_count"] == 3
    assert full_payload["runtime_schema"]["amplitude_stage"]["output_count"] == 48
    assert filtered_payload["runtime_schema"]["amplitude_stage"]["output_count"] == 16
    assert filtered_payload["runtime_schema"]["amplitude_stage"][
        "selected_color_sector_ids"
    ] == [0]
    assert [
        root["root_id"]
        for root in filtered_payload["runtime_schema"]["amplitude_stage"]["roots"]
    ] == list(range(16))
    assert all(
        root["dag_root_id"] >= 0
        for root in filtered_payload["runtime_schema"]["amplitude_stage"]["roots"]
    )

    _, artifact = write_generic_dag_process_artifact(
        manifest,
        tmp_path / "sector0",
        selected_color_sector_ids={0},
    )
    assert artifact["dag_summary"]["amplitude_root_count"] == 16
    assert artifact["compiled"]["selected_color_sector_ids"] == [0]


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

    assert reuse["full_color_sector_count"] == 576
    assert reuse["line_pairing_representative_sector_count"] == 24
    assert reuse["line_pairing_representative_sector_ids"][:4] == [
        0,
        24,
        48,
        72,
    ]


def test_selected_sector_manifest_uses_filtered_color_plan_for_large_line_counts() -> None:
    payload = build_generic_process_manifest(
        "d d~ > u u~ s s~ c c~ b b~",
        selected_color_sector_ids={0},
        max_coupling_orders={"QCD": 8},
        max_quark_pairs=5,
    ).to_json_dict()

    assert payload["color_plan"]["sector_count"] == 1
    assert payload["color_plan"]["truncated"] is False
    assert payload["planning_status"]["color_ready"] is True
    assert payload["planning_status"]["generic_evaluator_ready"] is True
    assert payload["lowering_status"]["full_tensor_network_ready"] is True
    assert payload["runtime_schema"]["amplitude_stage"]["output_count"] == 32


def test_reference_color_sector_can_be_selected_before_dag_construction() -> None:
    plan = build_color_plan("d d~ > u u~ s s~")

    assert select_leading_color_sector_ids_from_plan(
        plan,
        reference_color_order=(2, 4, 3, 6, 5, 1),
    ) == {18}
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
    assert blueprint.amplitude_stage.output_length == 64
    assert blueprint.amplitude_stage.blockers == ()
    assert blueprint.amplitude_stage.output_expressions


def test_generic_artifact_defaults_to_contributing_lc_sectors(
    tmp_path: Path,
) -> None:
    _, artifact = write_generic_dag_process_artifact(
        "d d~ > u u~ g",
        tmp_path / "runtime",
    )

    assert artifact["compiled"]["selected_color_sector_ids"] == list(range(8))
    assert artifact["dag_summary"]["amplitude_root_count"] == 64
    assert artifact["dag_summary"]["current_count"] == 260
    assert artifact["dag_summary"]["interaction_count"] == 268
    assert artifact["dag_summary"]["source_count"] == 80
    assert artifact["full_dag_summary"]["current_count"] == 260
    assert artifact["full_dag_summary"]["interaction_count"] == 268
    assert artifact["full_dag_summary"]["source_count"] == 80
    assert artifact["lowering_status"]["current_color_sectors"] == list(range(8))
    runtime_schema = artifact["runtime_schema"]
    assert runtime_schema["source_fill"]["source_count"] == 80
    assert runtime_schema["amplitude_stage"]["output_count"] == 64


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
    assert payload["planning_status"]["color_sector_count"] == 4
    assert payload["planning_status"]["color_ready"] is True
    assert payload["lowering_status"]["current_color_sector_count"] == 4
    assert payload["lowering_status"]["current_color_sectors"] == [0, 1, 2, 3]
    assert len(payload["lowering_status"]["color_sector_summaries"]) == 4
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
    assert source_sectors == {0, 1, 2, 3}
    assert payload["stage_plan"]["current_stages"][-1]["subset_size"] == 3
    assert payload["stage_plan"]["amplitude_stage"]["closure_count"] > 0


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
    assert len(stage_compiler["real_valued_inputs"]) == stage_compiler[
        "momentum_parameter_count"
    ]
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
        call["real_count"] == blueprint.momentum_parameter_count
        for call in calls
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
    assert payload["amplitude_stage"]["evaluator"]["output_len"] == 64
    assert [call["label"] for call in calls][-1] == "generic_amplitude_stage"
    assert [call["output_count"] for call in calls][-1] == 64
    assert all(call["output_count"] > 0 for call in calls)
    assert all(call["real_count"] <= call["param_count"] for call in calls)
    assert calls[-1]["real_count"] == 0


def test_generic_stage_evaluator_artifact_writer_default_symbolica_bridge(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from pyamplicol import dag_runtime

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
    ):
        del symbolica_settings, progress_callback
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
        dag_runtime,
        "_compile_symbolica_outputs",
        fake_compile_symbolica_outputs,
    )
    monkeypatch.setattr(
        dag_runtime,
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
