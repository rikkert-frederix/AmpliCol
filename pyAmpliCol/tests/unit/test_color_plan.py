from __future__ import annotations

from pyamplicol.color_plan import (
    build_color_plan,
    lc_line_pairing_representative_ids,
    lc_topology_replay_safe_groups,
)


def test_color_plan_builds_open_line_gluon_orderings() -> None:
    plan = build_color_plan("d d~ > z g g")

    assert plan.color_accuracy == "lc"
    assert plan.ready_for_leading_colour is True
    assert plan.sector_count == 2
    assert [sector.quark_lines[0].gluon_labels for sector in plan.sectors] == [
        (4, 5),
        (5, 4),
    ]
    assert plan.sectors[0].quark_lines[0].quark_label == 2
    assert plan.sectors[0].quark_lines[0].antiquark_label == 1
    assert plan.sectors[0].quark_lines[0].singlet_labels == ()
    assert plan.sectors[0].singlet_labels == (3,)
    assert plan.sectors[0].coloured_label_groups == ((2, 4, 5, 1),)
    assert plan.sectors[0].line_label_groups == ((2, 4, 5, 1),)
    assert plan.coloured_labels == (1, 2, 4, 5)


def test_color_plan_generates_arbitrary_quark_line_pairings() -> None:
    plan = build_color_plan("d d~ > u u~ s s~")

    assert plan.sector_count == 6
    assert all(sector.kind == "open-lines" for sector in plan.sectors)
    assert all(len(sector.quark_lines) == 3 for sector in plan.sectors)
    assert all(len(sector.coloured_label_groups) == 3 for sector in plan.sectors)
    assert all(len(sector.color_words) == 1 for sector in plan.sectors)
    assert all(sector.word_labels for sector in plan.sectors)
    direct_plan = build_color_plan("s s~ > u u~ d~ d")
    assert (6, 1, 3, 4, 2, 5) in {
        word
        for sector in direct_plan.sectors
        for word in sector.compatibility_words
    }


def test_color_plan_line_pairing_representatives_drop_only_block_orderings() -> None:
    three_line = build_color_plan("d d~ > u u~ s s~")
    four_line = build_color_plan("d d~ > u u~ s s~ c c~")

    assert three_line.sector_count == 6
    assert lc_line_pairing_representative_ids(three_line) == (
        0,
        1,
        2,
        3,
        4,
        5,
    )
    assert four_line.sector_count == 24
    assert len(lc_line_pairing_representative_ids(four_line)) == 24
    assert lc_line_pairing_representative_ids(four_line)[:4] == (0, 1, 2, 3)


def test_open_line_compatibility_includes_complete_block_permutations() -> None:
    plan = build_color_plan("u u~ > d~ d g")
    sector = next(
        sector
        for sector in plan.sectors
        if sector.color_words == ((2, 3, 4, 5, 1),)
    )

    assert sector is not None
    assert sector.color_words == ((2, 3, 4, 5, 1),)
    assert (4, 5, 1, 2, 3) in sector.compatibility_words


def test_open_line_legacy_orders_include_antiquark_first_orientation() -> None:
    plan = build_color_plan("g g > u u~ g")
    sector = plan.sector(4)

    assert sector is not None
    assert sector.color_words == ((3, 5, 1, 2, 4),)
    assert (4, 5, 1, 2, 3) in sector.legacy_order_words
    assert (4, 5, 1, 2, 3) not in sector.compatibility_words


def test_color_plan_keeps_singlets_global_to_open_colour_sectors() -> None:
    plan = build_color_plan("d d~ > u u~ z")

    assert plan.sector_count == 2
    assert all(sector.kind == "open-lines" for sector in plan.sectors)
    assert {(5 in group) for sector in plan.sectors for group in sector.line_label_groups} == {
        False,
    }
    assert {sector.singlet_labels for sector in plan.sectors} == {(5,)}
    assert {
        tuple(line.singlet_labels for line in sector.quark_lines)
        for sector in plan.sectors
    } == {
        ((), ()),
    }


def test_color_plan_does_not_duplicate_singlet_orderings() -> None:
    plan = build_color_plan("d d~ > e+ e- g")

    assert plan.sector_count == 1
    line = plan.sectors[0].quark_lines[0]
    assert line.gluon_labels == (5,)
    assert line.singlet_labels == ()
    assert plan.sectors[0].singlet_labels == (3, 4)


def test_color_plan_allows_flavour_changing_open_line() -> None:
    plan = build_color_plan("u d~ > w+ g")

    assert plan.sector_count == 1
    line = plan.sectors[0].quark_lines[0]
    assert line.quark_label == 2
    assert line.antiquark_label == 1
    assert line.gluon_labels == (4,)


def test_color_plan_builds_pure_gluon_single_trace_sectors() -> None:
    plan = build_color_plan("g g > g g")

    assert plan.sector_count == 3
    assert all(sector.kind == "single-trace" for sector in plan.sectors)
    assert all(sector.trace_labels[0] == 1 for sector in plan.sectors)
    assert {
        sector.trace_labels
        for sector in plan.sectors
    } == {
        (1, 2, 3, 4),
        (1, 2, 4, 3),
        (1, 3, 2, 4),
    }


def test_color_plan_respects_sector_cap() -> None:
    plan = build_color_plan("d d~ > z g g g", max_sectors=3)

    assert plan.truncated is True
    assert plan.sector_count == 3
    assert "max_sectors=3" in plan.diagnostics[0]
    assert plan.ready_for_leading_colour is False


def test_color_plan_records_nlc_full_colour_scaffold() -> None:
    plan = build_color_plan("d d~ > z g", color_accuracy="full")

    assert plan.color_accuracy == "full"
    assert plan.sector_count == 0
    assert plan.idenso_required is True
    assert "requires Idenso" in plan.diagnostics[0]


def test_color_plan_groups_isomorphic_open_line_sectors() -> None:
    plan = build_color_plan("d d~ > z g g g g")

    assert plan.sector_count == 24
    assert len(plan.topology_groups) == 1
    group = plan.topology_groups[0]
    assert group.representative_sector_id == 0
    assert group.sector_ids == tuple(range(24))
    assert len(group.label_permutations) == 24
    assert group.label_permutations[0] == tuple(
        zip(
            plan.sectors[0].line_label_groups[0],
            plan.sectors[0].line_label_groups[0],
        )
    )


def test_color_plan_topology_groups_are_process_generic() -> None:
    plan = build_color_plan("d d~ > e+ e- g g")

    assert plan.sector_count == 2
    assert len(plan.topology_groups) == 1
    group = plan.topology_groups[0]
    assert group.sector_ids == (0, 1)
    assert len(group.label_permutations) == 2


def test_lc_topology_replay_safe_groups_preserve_initial_labels() -> None:
    safe = build_color_plan("d d~ > e+ e- g g")
    pure_gluon = build_color_plan("g g > g g")
    gluon_initiated = build_color_plan("g g > u u~ g")

    assert [group.representative_sector_id for group in lc_topology_replay_safe_groups(safe)] == [
        0,
    ]
    assert lc_topology_replay_safe_groups(pure_gluon) == ()
    assert lc_topology_replay_safe_groups(gluon_initiated) == ()
