from __future__ import annotations

import pytest

from pyamplicol.color_plan import (
    build_color_plan,
    lc_line_pairing_representative_ids,
)
from pyamplicol.generic_dag import (
    ColorState,
    CurrentIndex,
    GenericDAGCompiler,
    compile_generic_dag,
    infer_minimal_coupling_order_limits,
    prune_dag_to_amplitude_roots,
)
from pyamplicol.model import AmplicolSMLeadingColorModel


def test_current_index_is_full_physics_state_identity() -> None:
    color = ColorState(accuracy="lc", sector_id=3, line_groups=(1,))
    kwargs = {
        "particle_id": 1,
        "external_mask": 0b101,
        "external_labels": (1, 3),
        "ordered_external_labels": (3, 1),
        "helicity_ancestry": 0b10010,
        "chirality": -1,
        "spin_state": -1,
        "flavour_flow": (1, 21, 1),
        "charge_flow": -1,
        "color_state": color,
        "momentum_mask": 0b101,
        "coupling_orders": (("qed", 1), ("QCD", 2)),
        "auxiliary_kind": "test-auxiliary",
    }
    index = CurrentIndex(**kwargs)

    assert index.overlaps(
        CurrentIndex(
            particle_id=21,
            external_mask=0b001,
            external_labels=(1,),
            helicity_ancestry=0b1,
            chirality=0,
            spin_state=1,
            flavour_flow=(21,),
            charge_flow=0,
            color_state=color,
            momentum_mask=0b001,
        )
    )
    payload = index.to_json_dict()
    assert set(payload).issuperset(
        {
            "particle_id",
            "external_mask",
            "external_labels",
            "ordered_external_labels",
            "helicity_ancestry",
            "chirality",
            "spin_state",
            "flavour_flow",
            "charge_flow",
            "color_state",
            "momentum_mask",
            "coupling_orders",
            "auxiliary_kind",
        }
    )
    assert payload["particle_id"] == kwargs["particle_id"]
    assert payload["external_labels"] == [1, 3]
    assert payload["ordered_external_labels"] == [3, 1]
    assert payload["helicity_ancestry"] == kwargs["helicity_ancestry"]
    assert payload["color_state"]["sector_id"] == 3
    assert payload["coupling_orders"] == [["QCD", 2], ["QED", 1]]
    assert payload["auxiliary_kind"] == "test-auxiliary"

    variants = [
        {**kwargs, "particle_id": 21},
        {
            **kwargs,
            "external_labels": (1, 2),
            "ordered_external_labels": (2, 1),
            "external_mask": 0b011,
            "momentum_mask": 0b011,
        },
        {**kwargs, "helicity_ancestry": 0b10011},
        {**kwargs, "chirality": 1},
        {**kwargs, "spin_state": (1, -1)},
        {**kwargs, "flavour_flow": (1, 22, 1)},
        {**kwargs, "charge_flow": 0},
        {
            **kwargs,
            "color_state": ColorState(
                accuracy="lc",
                sector_id=4,
                line_groups=(1,),
            ),
        },
        {**kwargs, "momentum_mask": 0b111},
        {**kwargs, "coupling_orders": (("QED", 2),)},
        {**kwargs, "auxiliary_kind": "other-auxiliary"},
        {**kwargs, "ordered_external_labels": (1, 3)},
    ]
    variant_indices = [CurrentIndex(**variant) for variant in variants]

    assert all(index != variant for variant in variant_indices)
    assert len({index, *variant_indices}) == 1 + len(variant_indices)
    assert "family" not in str(payload).lower()


def test_current_index_rejects_inconsistent_mask() -> None:
    with pytest.raises(ValueError, match="external_mask"):
        CurrentIndex(
            particle_id=21,
            external_mask=0b100,
            external_labels=(1,),
            helicity_ancestry=1,
            chirality=0,
            spin_state=1,
            flavour_flow=(21,),
            charge_flow=0,
            color_state=ColorState(accuracy="lc"),
            momentum_mask=0b1,
        )


def test_model_exposes_process_generic_vertex_lookup() -> None:
    model = AmplicolSMLeadingColorModel()

    assert model.vertices_accepting(1, 21) == model.vertices_for_inputs(1, 21)
    assert [vertex.kind for vertex in model.vertices_accepting(11, -11)] == [21, 21]


def test_generic_dag_compiler_builds_helicity_source_currents() -> None:
    dag = compile_generic_dag("d d~ > z g")

    assert dag.process.outgoing_particles == ("d~", "d", "z", "g")
    assert dag.color_plan.sector_count == 1
    assert len(dag.sources) == 9
    source_bits = [
        dag.currents[source_id].index.helicity_ancestry
        for source_id in dag.sources
    ]
    assert source_bits == [1 << index for index in range(len(dag.sources))]
    assert {dag.currents[source_id].source_helicity for source_id in dag.sources} == {
        -1,
        0,
        1,
    }
    assert dag.has_amplitudes is True
    assert 10 in dag.required_vertex_kinds


def test_generic_dag_can_compile_selected_lc_sector_without_family_branch() -> None:
    full = compile_generic_dag("d d~ > z g g g g", max_color_sectors=50)
    selected = compile_generic_dag(
        "d d~ > z g g g g",
        max_color_sectors=50,
        selected_color_sector_ids={0},
    )

    assert full.color_plan.sector_count == 24
    assert selected.color_plan.sector_count == 1
    assert {current.index.color_state.sector_id for current in selected.currents} == {0}
    assert len(selected.currents) * 24 == len(full.currents)
    assert len(selected.amplitude_roots) * 24 == len(full.amplitude_roots)


def test_generic_dag_prunes_dead_currents_after_amplitude_closure() -> None:
    raw = compile_generic_dag(
        "d d~ > z g g g g",
        max_color_sectors=50,
        selected_color_sector_ids={0},
        closure_side_mask_pruning=False,
        species_reachability_pruning=False,
    )
    pruned = prune_dag_to_amplitude_roots(raw)

    assert pruned.has_amplitudes is True
    assert len(pruned.amplitude_roots) == len(raw.amplitude_roots)
    assert len(pruned.currents) < len(raw.currents)
    assert len(pruned.interactions) < len(raw.interactions)
    assert {root.left_id for root in pruned.amplitude_roots}.issubset(
        {current.id for current in pruned.currents}
    )
    assert {root.right_id for root in pruned.amplitude_roots}.issubset(
        {current.id for current in pruned.currents}
    )


def test_generic_dag_uses_generic_closure_side_mask_pruning() -> None:
    unpruned = compile_generic_dag(
        "d d~ > z g g g",
        selected_color_sector_ids={0},
        closure_side_mask_pruning=False,
        species_reachability_pruning=False,
    )
    pruned = compile_generic_dag(
        "d d~ > z g g g",
        selected_color_sector_ids={0},
        species_reachability_pruning=False,
    )

    assert pruned.has_amplitudes is True
    assert len(pruned.currents) < len(unpruned.currents)
    assert len(pruned.interactions) < len(unpruned.interactions)
    assert len(pruned.amplitude_roots) == len(unpruned.amplitude_roots)
    sink_labels = {
        unpruned.currents[root.right_id].index.external_labels[0]
        for root in unpruned.amplitude_roots
        if len(unpruned.currents[root.right_id].index.external_labels) == 1
    }
    assert sink_labels
    assert all(
        not (
            len(current.index.external_labels) > 1
            and any(label in current.index.external_labels for label in sink_labels)
        )
        for current in pruned.currents
    )


def test_generic_dag_color_order_mask_pruning_preserves_lc_roots() -> None:
    unpruned = compile_generic_dag(
        "d d~ > z g g g g",
        selected_color_sector_ids={0},
        color_order_mask_pruning=False,
    )
    pruned = compile_generic_dag(
        "d d~ > z g g g g",
        selected_color_sector_ids={0},
    )
    word = pruned.color_plan.sectors[0].compatibility_words[0]
    word_positions = {label: index for index, label in enumerate(word)}

    assert pruned.has_amplitudes is True
    assert len(pruned.amplitude_roots) == len(unpruned.amplitude_roots)
    assert len(pruned.currents) == len(unpruned.currents)
    for current in pruned.currents:
        coloured_positions = tuple(
            sorted(
                word_positions[label]
                for label in current.index.external_labels
                if label in word_positions
            )
        )
        if not coloured_positions:
            continue
        assert coloured_positions == tuple(
            range(coloured_positions[0], coloured_positions[-1] + 1)
    )


def test_generic_dag_species_reachability_prunes_dead_particle_branches() -> None:
    without_reachability = compile_generic_dag(
        "d d~ > h h z",
        species_reachability_pruning=False,
    )
    with_reachability = compile_generic_dag("d d~ > h h z")
    backward_pruned = prune_dag_to_amplitude_roots(without_reachability)

    assert with_reachability.has_amplitudes is True
    assert len(with_reachability.amplitude_roots) == len(
        without_reachability.amplitude_roots
    )
    assert len(with_reachability.currents) < len(without_reachability.currents)
    assert len(with_reachability.interactions) < len(
        without_reachability.interactions
    )
    assert len(with_reachability.currents) == len(backward_pruned.currents)
    assert len(with_reachability.interactions) == len(backward_pruned.interactions)


def test_generic_dag_prunes_by_model_coupling_order_budget() -> None:
    no_ew = compile_generic_dag(
        "d d~ > z g g",
        max_coupling_orders={"QED": 0},
    )
    one_ew = compile_generic_dag(
        "d d~ > z g g",
        max_coupling_orders={"QED": 1},
    )

    assert no_ew.has_amplitudes is False
    assert no_ew.amplitude_roots == ()
    assert one_ew.has_amplitudes is True
    assert {tuple(order) for current in one_ew.currents for order in current.index.coupling_orders}
    assert max(
        dict(current.index.coupling_orders).get("QED", 0)
        for current in one_ew.currents
    ) <= 1


def test_generic_dag_reachability_tracks_coupling_order_states() -> None:
    broad = compile_generic_dag(
        "d d~ > z z z",
        max_coupling_orders={"QED": 3},
        species_reachability_pruning=False,
    )
    state_pruned = compile_generic_dag(
        "d d~ > z z z",
        max_coupling_orders={"QED": 3},
    )

    assert broad.has_amplitudes is True
    assert state_pruned.has_amplitudes is True
    assert len(state_pruned.amplitude_roots) == len(broad.amplitude_roots)
    assert len(state_pruned.currents) < len(broad.currents)
    assert len(state_pruned.interactions) < len(broad.interactions)
    assert max(
        dict(current.index.coupling_orders).get("QED", 0)
        for current in state_pruned.currents
    ) <= 3


def test_generic_dag_infers_minimal_coupling_order_envelope() -> None:
    limits = infer_minimal_coupling_order_limits(
        "d d~ > z g g",
        selected_color_sector_ids={0},
    )
    dag = compile_generic_dag(
        "d d~ > z g g",
        selected_color_sector_ids={0},
        max_coupling_orders=limits,
    )

    assert limits == {"QCD": 2, "QED": 1}
    assert dag.has_amplitudes is True
    assert all(
        dict(current.index.coupling_orders).get("QED", 0) <= 1
        and dict(current.index.coupling_orders).get("QCD", 0) <= 2
        for current in dag.currents
    )


def test_generic_dag_prunes_ignored_vertex_kinds_without_process_family_tag() -> None:
    dag = compile_generic_dag("d d~ > z g", ignored_vertex_kinds={10, 11})

    assert dag.has_amplitudes is False
    assert all(
        interaction.vertex_kind not in {10, 11}
        for interaction in dag.interactions
    )


def test_generic_dag_discovers_dilepton_vector_current_from_vertices() -> None:
    dag = compile_generic_dag("d d~ > e+ e- g")

    assert dag.has_amplitudes is True
    dilepton_currents = [
        current
        for current in dag.currents
        if current.index.particle_id in {22, 23}
        and set(current.index.external_labels) == {3, 4}
        and current.index.helicity_ancestry != 0
    ]
    assert dilepton_currents
    assert any(
        set(current.index.flavour_flow[:2]) == {-11, 11}
        and current.index.flavour_flow[-1] in {22, 23}
        for current in dilepton_currents
    )
    assert any(kind in dag.required_vertex_kinds for kind in (21, 22))


def test_generic_dag_keeps_one_mirrored_colour_singlet_current_orientation() -> None:
    dag = compile_generic_dag("d d~ > h z")

    assert dag.has_amplitudes is True
    scalar_vector_kinds_by_result: dict[int, set[int]] = {}
    for interaction in dag.interactions:
        result = dag.currents[interaction.result_id]
        if result.index.external_labels != (3, 4) or result.index.particle_id != 23:
            continue
        scalar_vector_kinds_by_result.setdefault(interaction.result_id, set()).add(
            interaction.vertex_kind
        )

    assert scalar_vector_kinds_by_result
    assert all(
        not {18, 19}.issubset(kinds)
        for kinds in scalar_vector_kinds_by_result.values()
    )


def test_generic_dag_uses_same_compiler_for_charged_current_w_process() -> None:
    dag = GenericDAGCompiler().compile("u d~ > w+ g")

    assert dag.has_amplitudes is True
    assert any(
        current.index.particle_id == 24
        and current.index.external_labels == (3,)
        for current in dag.currents
    )
    assert any(
        current.index.particle_id == 2
        and current.index.flavour_flow == (1, 2)
        for current in dag.currents
    )
    assert 10 in dag.required_vertex_kinds or 11 in dag.required_vertex_kinds


def test_generic_dag_carries_arbitrary_quark_line_colour_sectors() -> None:
    dag = compile_generic_dag(
        "d d~ > u u~ s s~",
        max_currents=10000,
    )

    assert dag.process.quark_lines.quark_pair_count == 3
    assert dag.color_plan.sector_count == 6
    assert {
        dag.currents[source_id].index.color_state.sector_id
        for source_id in dag.sources
    } == set(range(6))
    assert all(sector.compatibility_words for sector in dag.color_plan.sectors)
    assert any(
        current.index.ordered_external_labels != current.index.external_labels
        for current in dag.currents
        if not current.is_source
    )
    assert any(
        current.index.particle_id == 21
        and len(current.index.external_labels) == 2
        for current in dag.currents
        if not current.is_source
    )
    assert all(
        "q-qbar" not in str(current.index.to_json_dict())
        for current in dag.currents[:20]
    )


def test_generic_dag_builds_four_quark_line_selected_sector() -> None:
    process = "d d~ > u u~ s s~ c c~"
    color_plan = build_color_plan(process, max_sectors=2000)
    representative_ids = lc_line_pairing_representative_ids(color_plan)

    dag = compile_generic_dag(
        process,
        selected_color_sector_ids={representative_ids[0]},
        max_color_sectors=2000,
        max_currents=5000,
    )

    assert dag.process.quark_lines.quark_pair_count == 4
    assert color_plan.sector_count == 24
    assert len(representative_ids) == 24
    assert dag.color_plan.sector_count == 1
    assert dag.has_amplitudes is True
    assert dag.truncated is False
    assert dag.required_vertex_kinds
    assert all(
        "Z+gluons" not in str(current.index.to_json_dict())
        for current in dag.currents[:20]
    )


def test_generic_dag_prunes_by_lc_current_line_group_budget() -> None:
    full = compile_generic_dag("d d~ > u u~")
    pruned = compile_generic_dag(
        "d d~ > u u~",
        max_lc_current_line_groups=1,
    )

    assert full.has_amplitudes is True
    assert any(
        len(current.index.color_state.line_groups) > 1
        for current in full.currents
    )
    assert pruned.has_amplitudes is False
    assert all(
        len(current.index.color_state.line_groups) <= 1
        for current in pruned.currents
    )
    assert len(pruned.currents) < len(full.currents)


def test_lc_colour_engine_keeps_singlet_attachments_on_allocated_line() -> None:
    model = AmplicolSMLeadingColorModel()
    dag = compile_generic_dag(
        "d d~ > u u~ z",
        max_currents=50000,
        max_color_sectors=2000,
    )

    assert dag.has_amplitudes is True
    for interaction in dag.interactions:
        left = dag.currents[interaction.left_id]
        right = dag.currents[interaction.right_id]
        left_colored = (
            left.index.particle_id == 99
            or model.color_rep(left.index.particle_id) != 1
        )
        right_colored = (
            right.index.particle_id == 99
            or model.color_rep(right.index.particle_id) != 1
        )
        if left_colored == right_colored:
            continue
        colored = left if left_colored else right
        singlet = right if left_colored else left
        singlet_groups = set(singlet.index.color_state.line_groups)
        colored_groups = set(colored.index.color_state.line_groups)
        assert not singlet_groups or singlet_groups.issubset(colored_groups)
