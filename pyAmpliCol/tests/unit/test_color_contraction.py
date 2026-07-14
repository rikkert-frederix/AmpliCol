from __future__ import annotations

from pyamplicol.color_contraction import (
    ColorGroupDescriptor,
    amplicol_color_factors,
    build_color_contraction_plan,
)
from pyamplicol.color_plan import build_color_plan


def _groups_for_plan(process: str, *, color_accuracy: str) -> tuple[ColorGroupDescriptor, ...]:
    plan = build_color_plan(process, color_accuracy=color_accuracy)
    return tuple(
        ColorGroupDescriptor(
            group_id=index,
            helicity_key=("h",),
            sector_id=sector.id,
            word=tuple(sector.word_labels or sector.color_words[0]),
            helicity_weight=1.0,
        )
        for index, sector in enumerate(plan.sectors)
    )


def test_one_quark_line_nlc_and_full_factor_match_amplicol_convention() -> None:
    plan = build_color_plan("d d~ > z g", color_accuracy="full")
    sector = plan.sectors[0]

    assert amplicol_color_factors(plan, sector, sector) == (9.0, 8.0, 8.0)


def test_color_singlet_metric_is_unity_at_every_accuracy() -> None:
    for color_accuracy in ("lc", "nlc", "full"):
        plan = build_color_plan("z z > h h", color_accuracy=color_accuracy)
        sector = plan.sectors[0]

        assert amplicol_color_factors(plan, sector, sector) == (1.0, 1.0, 1.0)
        contraction = build_color_contraction_plan(
            plan,
            _groups_for_plan("z z > h h", color_accuracy=color_accuracy),
        )
        if color_accuracy == "lc":
            assert contraction is None
        else:
            assert contraction is not None
            assert tuple(entry.weight_re for entry in contraction.entries) == (1.0,)


def test_pure_gluon_full_colour_builds_sparse_off_diagonal_metric() -> None:
    plan = build_color_plan("g g > g g", color_accuracy="full")
    contraction = build_color_contraction_plan(
        plan,
        _groups_for_plan("g g > g g", color_accuracy="full"),
    )

    assert contraction is not None
    assert contraction.supported is True
    assert contraction.group_count == len(plan.sectors)
    assert any(entry.left_group_id != entry.right_group_id for entry in contraction.entries)
    assert {entry.weight_re for entry in contraction.entries}.issuperset(
        {50.666666666666664, -5.333333333333333}
    )


def test_one_quark_line_nlc_prunes_like_amplicol_for_two_gluons() -> None:
    plan = build_color_plan("g g > t t~ g g", color_accuracy="nlc")
    fortran_orders = (
        (3, 5, 2, 1, 6, 4),
        (3, 6, 1, 2, 5, 4),
        (3, 2, 1, 6, 5, 4),
        (3, 5, 6, 1, 2, 4),
        (3, 5, 1, 6, 2, 4),
        (3, 2, 6, 1, 5, 4),
        (3, 1, 6, 2, 5, 4),
        (3, 5, 2, 6, 1, 4),
        (3, 5, 1, 2, 6, 4),
        (3, 6, 2, 1, 5, 4),
        (3, 1, 2, 6, 5, 4),
        (3, 5, 6, 2, 1, 4),
        (3, 2, 5, 1, 6, 4),
        (3, 6, 1, 5, 2, 4),
        (3, 1, 6, 5, 2, 4),
        (3, 2, 5, 6, 1, 4),
        (3, 2, 1, 5, 6, 4),
        (3, 6, 5, 1, 2, 4),
        (3, 1, 5, 6, 2, 4),
        (3, 2, 6, 5, 1, 4),
        (3, 1, 5, 2, 6, 4),
        (3, 6, 2, 5, 1, 4),
        (3, 1, 2, 5, 6, 4),
        (3, 6, 5, 2, 1, 4),
    )
    fortran_first_row = {
        0: 151.7037037037037,
        3: 47.407407407407405,
        6: 36.74074074074074,
        7: -37.925925925925924,
        8: -37.925925925925924,
        9: 36.74074074074074,
        12: -37.925925925925924,
        13: 36.74074074074074,
        14: 42.074074074074076,
        22: 47.407407407407405,
    }
    sector_by_word = {tuple(sector.word_labels): sector for sector in plan.sectors}
    ordered_sectors = tuple(sector_by_word[word] for word in fortran_orders)

    left = ordered_sectors[0]
    for j, right in enumerate(ordered_sectors):
        py_value = amplicol_color_factors(plan, left, right)[1]
        if j != 0:
            py_value *= 2.0
        assert py_value == fortran_first_row.get(j, 0.0)


def test_two_quark_line_full_colour_keeps_open_line_block_interference() -> None:
    plan = build_color_plan("d d~ > u u~", color_accuracy="full")
    contraction = build_color_contraction_plan(
        plan,
        _groups_for_plan("d d~ > u u~", color_accuracy="full"),
    )

    assert plan.sector_count == 4
    assert contraction is not None
    assert contraction.supported is True
    assert any(
        entry.left_group_id != entry.right_group_id
        and entry.weight_re == -3.0
        and entry.symmetry_factor == 2.0
        for entry in contraction.entries
    )


def test_two_quark_line_nlc_prunes_nncl_entries_like_amplicol() -> None:
    plan = build_color_plan("d d~ > t t~ g g", color_accuracy="nlc")
    first_quark_label = plan.process.quark_labels[0]

    def rotate_to_reference(word: tuple[int, ...]) -> tuple[int, ...]:
        offset = word.index(first_quark_label)
        return word[offset:] + word[:offset]

    # This is the colour-order matrix emitted by amplicol_color_probe with
    # AMPICOL_COLOR_PROBE_MATRIX=1 for d d~ > t t~ g g at NLC.  Entries are
    # upper triangular and already include AmpliCol's off-diagonal factor 2.
    fortran_orders = (
        (2, 4, 3, 5, 6, 1),
        (2, 4, 3, 6, 5, 1),
        (2, 6, 4, 3, 5, 1),
        (2, 6, 5, 4, 3, 1),
        (2, 5, 4, 3, 6, 1),
        (2, 5, 6, 4, 3, 1),
        (2, 5, 6, 1, 3, 4),
        (2, 6, 5, 1, 3, 4),
        (2, 5, 1, 3, 6, 4),
        (2, 6, 1, 3, 5, 4),
        (2, 1, 3, 6, 5, 4),
        (2, 1, 3, 5, 6, 4),
    )
    fortran_entries = {
        (0, 0): 64.0,
        (0, 1): -16.0,
        (0, 3): 16.0,
        (0, 5): 16.0,
        (0, 6): -42.666666666666664,
        (0, 9): -42.666666666666664,
        (0, 11): -42.666666666666664,
    }
    sector_by_rotated_word = {
        rotate_to_reference(tuple(sector.word_labels)): sector
        for sector in plan.sectors
        if rotate_to_reference(tuple(sector.word_labels)) in fortran_orders
    }
    ordered_sectors = tuple(sector_by_rotated_word[word] for word in fortran_orders)

    left = ordered_sectors[0]
    for j, right in enumerate(ordered_sectors):
        py_value = amplicol_color_factors(plan, left, right)[1]
        if j != 0:
            py_value *= 2.0
        assert py_value == fortran_entries.get((0, j), 0.0)


def test_nlc_full_colour_supports_no_gluon_three_quark_lines() -> None:
    plan = build_color_plan("d d~ > u u~ s s~", color_accuracy="nlc")
    contraction = build_color_contraction_plan(
        plan,
        _groups_for_plan("d d~ > u u~ s s~", color_accuracy="nlc"),
    )

    assert contraction is not None
    assert contraction.supported is True
    assert contraction.reason is None
    assert plan.sector_count == 36
    assert {entry.weight_re for entry in contraction.entries} == {-9.0, 3.0, 27.0}

    full_plan = build_color_plan("d d~ > u u~ s s~", color_accuracy="full")
    full_contraction = build_color_contraction_plan(
        full_plan,
        _groups_for_plan("d d~ > u u~ s s~", color_accuracy="full"),
    )

    assert full_contraction is not None
    assert full_contraction.supported is True
    assert {entry.weight_re for entry in full_contraction.entries} == {
        -9.0,
        3.0,
        27.0,
    }


def test_nlc_full_colour_supports_multi_quark_lines_with_gluons() -> None:
    plan = build_color_plan("d d~ > u u~ s s~ g", color_accuracy="nlc")
    contraction = build_color_contraction_plan(
        plan,
        _groups_for_plan("d d~ > u u~ s s~ g", color_accuracy="nlc"),
    )

    assert contraction is not None
    assert contraction.supported is True
    assert contraction.reason is None
    assert {entry.weight_re for entry in contraction.entries} == {-27.0, 9.0, 72.0}

    full_plan = build_color_plan("d d~ > u u~ s s~ g", color_accuracy="full")
    full_contraction = build_color_contraction_plan(
        full_plan,
        _groups_for_plan("d d~ > u u~ s s~ g", color_accuracy="full"),
    )

    assert full_contraction is not None
    assert full_contraction.supported is True
    assert {entry.weight_re for entry in full_contraction.entries} == {
        -24.0,
        8.0,
        72.0,
    }
