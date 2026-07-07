from __future__ import annotations

from pyamplicol.process_ir import build_process_ir, build_process_set_ir
from pyamplicol.process_support import describe_process_content


def test_process_ir_keeps_physical_and_all_outgoing_conventions() -> None:
    ir = build_process_ir("d d~ > Z g")

    assert ir.process == "d d~ > z g"
    assert ir.key == "d_dbar_to_z_g"
    assert [leg.particle for leg in ir.legs] == ["d", "d~", "z", "g"]
    assert [leg.outgoing_particle for leg in ir.legs] == ["d~", "d", "z", "g"]
    assert ir.initial_pdgs == (1, -1)
    assert ir.final_pdgs == (23, 21)
    assert ir.outgoing_pdgs == (-1, 1, 23, 21)
    assert ir.vector_labels == (3,)
    assert ir.gluon_labels == (4,)
    assert ir.quark_lines.quark_count == 1
    assert ir.quark_lines.antiquark_count == 1
    assert ir.quark_lines.balanced is True


def test_process_ir_summarizes_arbitrary_quark_line_candidates() -> None:
    ir = build_process_ir("d d~ > u u~ s s~")

    assert ir.outgoing_particles == ("d~", "d", "u", "u~", "s", "s~")
    assert ir.quark_lines.quark_pair_count == 3
    assert ir.quark_lines.balanced is True
    assert ir.quark_labels == (2, 3, 5)
    assert ir.antiquark_labels == (1, 4, 6)
    assert ir.singlet_labels == ()


def test_process_ir_marks_multi_singlet_and_lepton_content() -> None:
    ir = build_process_ir("d d~ > e+ e- z g")

    assert ir.vector_labels == (5,)
    assert ir.lepton_labels == (3, 4)
    assert ir.gluon_labels == (6,)
    assert ir.has_multiple_nonleptonic_singlets is False

    multi = build_process_ir("d d~ > z z g")
    assert multi.vector_labels == (3, 4)
    assert multi.has_multiple_nonleptonic_singlets is True


def test_process_set_ir_uses_legacy_enumerator_entries() -> None:
    process_set = build_process_set_ir("p p > Z g")

    assert len(process_set.entries) == 10
    assert process_set.default_key == "d_dbar_to_g_z"
    first = process_set.entries[0]
    assert first.key == "d_dbar_to_g_z"
    assert first.process == "d d~ > g z"
    assert first.ir.outgoing_particles == ("d~", "d", "g", "z")
    assert first.n_groups == 1
    assert first.n_records == 1


def test_process_support_content_reuses_canonical_ir() -> None:
    content = describe_process_content("d d~ > z z g")

    assert content.key == "d_dbar_to_z_z_g"
    assert content.ir.vector_labels == (3, 4)
    assert content.ir.gluon_labels == (5,)
    assert content.has_multiple_nonleptonic_singlets is True
    assert content.all_outgoing == content.ir.outgoing_particles
