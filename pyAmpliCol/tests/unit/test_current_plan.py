from __future__ import annotations

from pyamplicol.color_plan import build_color_plan
from pyamplicol.current_plan import (
    build_generic_current_plan,
    build_generic_stage_plan,
)
from pyamplicol.process_ir import build_process_ir


def test_generic_current_plan_builds_lc_vector_process_reachability() -> None:
    plan = build_generic_current_plan("d d~ > z g")

    assert plan.color_sectors == (0,)
    assert len(plan.sources) > 0
    assert plan.has_closure is True
    assert plan.truncated is False
    assert 8 not in {kind for kind, _ in plan.required_vertex_kind_counts}
    assert plan.unimplemented_vertex_kinds == ()
    assert set(plan.ready_vertex_kinds).issuperset({6, 10})
    assert plan.pending_vertex_kinds == ()
    assert plan.full_tensor_network_ready is True
    assert plan.process.outgoing_particles == ("d~", "d", "z", "g")
    assert all(hasattr(current, "index") for current in plan.currents)


def test_generic_current_plan_uses_model_owned_colour_filtering() -> None:
    lc_plan = build_generic_current_plan("d d~ > z g")
    full_plan = build_generic_current_plan(
        build_process_ir("d d~ > z g", color_accuracy="full")
    )

    lc_kinds = {kind for kind, _ in lc_plan.required_vertex_kind_counts}
    full_kinds = {kind for kind, _ in full_plan.required_vertex_kind_counts}
    assert 4 not in lc_kinds
    assert 4 in full_kinds
    assert full_plan.color_sectors == (0,)


def test_generic_current_plan_keys_currents_by_leading_colour_sector() -> None:
    plan = build_generic_current_plan("d d~ > z g g")
    payload = plan.to_json_dict()

    assert plan.color_sectors == (0, 1)
    assert {plan.currents[source].index.color_state.sector_id for source in plan.sources} == {0, 1}
    assert all(
        plan.currents[interaction.left_id].index.color_state.sector_id
        == plan.currents[interaction.right_id].index.color_state.sector_id
        == plan.currents[interaction.result_id].index.color_state.sector_id
        for interaction in plan.interactions
    )
    assert all(
        plan.currents[closure.left_id].index.color_state.sector_id
        == plan.currents[closure.right_id].index.color_state.sector_id
        for closure in plan.closures
    )
    summaries = payload["color_sector_summaries"]
    assert [summary["color_sector"] for summary in summaries] == [0, 1]
    assert summaries[0]["current_count"] == summaries[1]["current_count"]
    assert summaries[0]["interaction_count"] == summaries[1]["interaction_count"]
    assert summaries[0]["closure_count"] == summaries[1]["closure_count"]


def test_generic_current_plan_keeps_singlet_currents_in_colour_sector_support() -> None:
    color_plan = build_color_plan("d d~ > u u~ z")
    current_plan = build_generic_current_plan("d d~ > u u~ z")

    assert color_plan.sector_count == 2
    assert len(current_plan.currents) > 0
    assert len(current_plan.interactions) > 0
    assert len(current_plan.closures) > 0
    for current in current_plan.currents:
        if current.is_source:
            continue
        sector = color_plan.sector(current.index.color_state.sector_id)
        assert sector is not None
        labels = set(current.index.external_labels)
        sector_labels = set().union(*(set(group) for group in sector.line_label_groups))
        supported_labels = sector_labels.union(sector.singlet_labels)
        assert labels.issubset(supported_labels)


def test_generic_current_plan_finds_explicit_lepton_vector_currents() -> None:
    plan = build_generic_current_plan("d d~ > e+ e- g")

    assert plan.has_closure is True
    assert plan.pending_vertex_kinds == ()
    assert set(plan.ready_vertex_kinds).intersection({21, 22})
    assert plan.unimplemented_vertex_kinds == ()
    assert plan.full_tensor_network_ready is True
    assert any(
        current.index.particle_id in {22, 23}
        and set(current.index.external_labels) == {3, 4}
        for current in plan.currents
    )


def test_generic_current_plan_builds_multi_boson_reachability() -> None:
    plan = build_generic_current_plan("d d~ > z z g")

    assert plan.has_closure is True
    assert plan.unimplemented_vertex_kinds == ()
    assert plan.pending_vertex_kinds == ()
    assert plan.full_tensor_network_ready is True
    payload = plan.to_json_dict()
    assert payload["closure_count"] == len(plan.closures)
    assert payload["full_tensor_network_ready"] is True


def test_generic_current_plan_builds_higgs_vector_reachability() -> None:
    plan = build_generic_current_plan("d d~ > h z g")

    assert plan.has_closure is True
    assert plan.unimplemented_vertex_kinds == ()
    assert plan.pending_vertex_kinds == ()
    assert plan.full_tensor_network_ready is True
    assert 18 in {kind for kind, _ in plan.required_vertex_kind_counts}


def test_generic_current_plan_handles_more_than_three_quark_pair_candidates() -> None:
    plan = build_generic_current_plan("d d~ > u u~ s s~")
    payload = plan.to_json_dict()

    assert plan.process.quark_lines.quark_pair_count == 3
    assert plan.color_sectors == tuple(range(6))
    assert len(plan.sources) == 72
    assert plan.has_closure is True
    assert plan.unimplemented_vertex_kinds == ()
    assert plan.pending_vertex_kinds == ()
    assert all(
        len(current.index.external_labels) <= 5
        for current in plan.currents
        if not current.is_source
    )
    assert len(payload["color_sector_summaries"]) == 6
    assert {
        summary["color_sector"] for summary in payload["color_sector_summaries"]
    } == set(range(6))
    assert sum(summary["closure_count"] for summary in payload["color_sector_summaries"]) == len(plan.closures)


def test_generic_current_plan_honors_quark_pair_and_sector_caps() -> None:
    capped = build_generic_current_plan(
        "d d~ > u u~ s s~",
        max_quark_pairs=2,
    )
    selected = build_generic_current_plan(
        "d d~ > u u~ s s~ c c~",
        selected_color_sector_ids={0},
        max_color_sectors=2000,
        max_currents=5000,
    )

    assert capped.process.quark_lines.quark_pair_count == 3
    assert capped.currents == ()
    assert capped.interactions == ()
    assert capped.closures == ()
    assert capped.has_closure is False

    assert selected.process.quark_lines.quark_pair_count == 4
    assert selected.color_sectors == (0,)
    assert selected.has_closure is True
    assert selected.truncated is False
    assert selected.full_tensor_network_ready is True


def test_generic_current_plan_accepts_canonical_process_ir() -> None:
    ir = build_process_ir("u d~ > w+ g")

    plan = build_generic_current_plan(ir)

    assert plan.process is ir
    assert plan.has_closure is True
    assert plan.process.vector_labels == (3,)


def test_generic_stage_plan_groups_currents_by_external_subset_size() -> None:
    plan = build_generic_current_plan("d d~ > z g")
    stage_plan = build_generic_stage_plan(plan)

    assert stage_plan.stage_count == len(stage_plan.current_stages) + 1
    assert [stage.subset_size for stage in stage_plan.current_stages] == [2, 3]
    assert sum(len(stage.interaction_ids) for stage in stage_plan.current_stages) == len(
        plan.interactions
    )
    assert all(stage.current_ids for stage in stage_plan.current_stages)
    assert stage_plan.amplitude_stage.closure_ids == tuple(range(len(plan.closures)))
    assert stage_plan.full_tensor_network_ready is True


def test_generic_stage_plan_marks_higgs_vector_lowerings_ready() -> None:
    plan = build_generic_current_plan("d d~ > h z g")
    stage_plan = plan.build_stage_plan()

    assert all(
        stage.unimplemented_vertex_kinds == ()
        for stage in stage_plan.current_stages
    )
    assert any(
        18 in stage.ready_vertex_kinds
        for stage in stage_plan.current_stages
    )
    assert stage_plan.amplitude_stage.pending_vertex_kinds == ()
    assert stage_plan.amplitude_stage.full_tensor_network_ready is True
    payload = stage_plan.to_json_dict()
    assert payload["current_stages"][0]["subset_size"] == 2
    assert payload["amplitude_stage"]["closure_count"] == len(plan.closures)


def test_generic_stage_plan_handles_many_quark_line_stage_metadata() -> None:
    plan = build_generic_current_plan("d d~ > u u~ s s~")
    stage_plan = plan.build_stage_plan()
    payload = stage_plan.to_json_dict()

    assert [stage.subset_size for stage in stage_plan.current_stages] == [2, 3, 4, 5]
    assert stage_plan.amplitude_stage.unimplemented_vertex_kinds == ()
    assert stage_plan.amplitude_stage.pending_vertex_kinds == ()
    assert len(stage_plan.amplitude_stage.closure_ids) == len(plan.closures)
    assert stage_plan.amplitude_stage.full_tensor_network_ready is True
    stage_summaries = payload["current_stages"][0]["color_sector_summaries"]
    assert {summary["color_sector"] for summary in stage_summaries}.issubset(
        set(range(36))
    )
    assert len(stage_summaries) >= 1
    assert payload["amplitude_stage"]["color_sector_summaries"]
    assert all(summary["pending_vertex_kinds"] == [] for summary in stage_summaries)
