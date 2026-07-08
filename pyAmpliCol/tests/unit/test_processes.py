from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from pyamplicol.processes import (
    ANTIQUARKS,
    ProcessEnumerator,
    ProcessOptions,
    QUARKS,
    _process_dedupe_signature,
    _tokenize_side,
    build_generic_process_selection_report,
    enumerate_generic_process_set,
    enumerate_process_set,
    expand_process_variants,
    split_process_set,
)


def test_structured_process_enumeration_for_z_plus_two_gluons() -> None:
    enumeration = ProcessEnumerator().enumerate("d d~ > z g g")

    assert enumeration.request.initial_state == ("d~", "d")
    assert enumeration.request.rest == ("z", "g", "g")
    assert len(enumeration.unique_processes) == 1
    assert len(enumeration.groups) == 1
    assert enumeration.n_records == 1

    record = enumeration.groups[0].records[0]
    assert record.process == ("d~", "d", "g", "g", "z")
    assert record.identical_factor == 2.0


def test_structured_process_enumeration_prefilters_concrete_qcd_request() -> None:
    enumeration = ProcessEnumerator().enumerate("d d~ > z g g g g g g")

    assert len(enumeration.unique_processes) == 1
    assert len(enumeration.groups) == 1
    assert enumeration.n_records == 1
    assert enumeration.groups[0].records[0].process == (
        "d~",
        "d",
        "g",
        "g",
        "g",
        "g",
        "g",
        "g",
        "z",
    )


def test_process_variant_expansion_supports_repetition_and_anonymous_slots() -> None:
    assert expand_process_variants("d d~ > Z 4*g") == (
        "d d~ > z g g g g",
    )
    assert expand_process_variants("d d~ > Z 4g") == (
        "d d~ > z g g g g",
    )
    variants = expand_process_variants("d d~ > Z 2*[d g]")

    assert variants == (
        "d d~ > z d d",
        "d d~ > z d g",
        "d d~ > z g d",
        "d d~ > z g g",
    )
    assert expand_process_variants("d d~ > Z 3*[d g ]")[-1] == (
        "d d~ > z g g g"
    )


def test_process_set_supports_plan_syntax_with_anonymous_slots() -> None:
    request = "d d~ > e+ e- j j g g | p p > Z g [d  g]"

    assert split_process_set(request) == (
        "d d~ > e+ e- j j g g",
        "p p > Z g [d  g]",
    )
    assert expand_process_variants(request) == (
        "d d~ > e+ e- j j g g",
        "p p > z g d",
        "p p > z g g",
    )


def test_process_set_expands_builtin_p_to_concrete_partonic_entries() -> None:
    process_set = enumerate_process_set("p p > z g")

    assert len(process_set.entries) == 10
    assert process_set.default_key == "d_dbar_to_g_z"
    assert [entry.process for entry in process_set.entries[:5]] == [
        "d d~ > g z",
        "u u~ > g z",
        "s s~ > g z",
        "c c~ > g z",
        "b b~ > g z",
    ]
    assert [entry.process for entry in process_set.entries[5:]] == [
        "d~ d > g z",
        "u~ u > g z",
        "s~ s > g z",
        "c~ c > g z",
        "b~ b > g z",
    ]


def test_generic_process_set_uses_lightweight_concrete_entries() -> None:
    process_set = enumerate_generic_process_set("d d~ > z g g g g")

    assert len(process_set.entries) == 1
    entry = process_set.entries[0]
    assert entry.process == "d d~ > z g g g g"
    assert entry.enumeration.n_records == 1
    assert entry.enumeration.groups[0].records[0].process == (
        "d~",
        "d",
        "z",
        "g",
        "g",
        "g",
        "g",
    )


def test_generic_process_set_can_cap_inclusive_quark_pairs() -> None:
    unrestricted = enumerate_generic_process_set("p p > z 4j")
    capped = enumerate_generic_process_set("p p > z 4j", max_quark_pairs=1)

    assert len(capped.entries) < len(unrestricted.entries)
    assert capped.entries[0].process == "d d~ > g g g g z"
    assert all(
        min(
            sum(
                particle in QUARKS
                for particle in entry.enumeration.groups[0].records[0].process
            ),
            sum(
                particle in ANTIQUARKS
                for particle in entry.enumeration.groups[0].records[0].process
            ),
        )
        <= 1
        for entry in capped.entries
    )


def test_generic_process_set_supports_max_quark_line_report() -> None:
    unrestricted = build_generic_process_selection_report("p p > z 4j")
    capped = build_generic_process_selection_report("p p > z 4j", max_quark_pairs=1)

    assert capped.selected_count < unrestricted.selected_count
    assert dict(capped.rejection_counts)["max-quark-lines"] > 0
    assert all(record.quark_lines <= 1 for record in capped.selected_records)


def test_generic_process_set_deduplicates_final_permutations_only() -> None:
    process_set = enumerate_generic_process_set(
        "d d~ > z g | d d~ > g z | g d~ > d~ z"
    )

    assert [entry.process for entry in process_set.entries] == [
        "d d~ > z g",
        "g d~ > d~ z",
    ]
    assert process_set.selection_report is not None
    assert process_set.selection_report.duplicate_count == 1


def test_generic_process_set_canonicalizes_symmetric_pp_beam_order() -> None:
    report = build_generic_process_selection_report("p p > z g")

    assert report.selected_count == 5
    assert report.duplicate_count == 5
    assert [record.process for record in report.selected_records] == [
        "d d~ > g z",
        "u u~ > g z",
        "s s~ > g z",
        "c c~ > g z",
        "b b~ > g z",
    ]


def test_generic_process_selection_reports_higher_multiplicity_rejections() -> None:
    report = build_generic_process_selection_report("p p > e+ e- j j j")

    assert report.candidate_count > report.selected_count
    assert report.selected_count == 145
    reasons = dict(report.rejection_counts)
    assert reasons["charge"] > 0
    assert reasons["fermion-family"] > 0


def test_generic_process_selection_handles_mixed_inline_multiparticle_syntax() -> None:
    report = build_generic_process_selection_report(
        "d d~ > e+ e- j j g g | p p > Z g [d  g]"
    )

    assert report.selected_count == 12
    assert report.duplicate_count == 6
    assert report.selected_records[0].process == "d d~ > d d~ g g e+ e-"


@pytest.mark.parametrize(
    "process_request",
    [
        "p p > z j j",
        "p p > z 4j",
        "d d~ > e+ e- j j g g | p p > Z g [d  g]",
    ],
)
def test_generic_prefilter_preserves_legacy_selected_physical_signatures(
    process_request: str,
) -> None:
    legacy = build_generic_process_selection_report(
        process_request,
        use_prefilter=False,
    )
    prefiltered = build_generic_process_selection_report(
        process_request,
        use_prefilter=True,
    )

    assert {
        _selection_record_signature(record) for record in prefiltered.selected_records
    } == {_selection_record_signature(record) for record in legacy.selected_records}


def _selection_record_signature(record: object) -> str:
    initial, _, final = str(getattr(record, "process")).partition(">")
    source = str(getattr(record, "source"))
    source_request = ProcessEnumerator().parse(source)
    symmetric_initial = (
        len(source_request.initial_state) == 2
        and source_request.initial_state[0] == source_request.initial_state[1]
        and source_request.initial_state[0] in {"p", "j"}
    )
    return _process_dedupe_signature(
        tuple(_tokenize_side(initial.strip())),
        tuple(_tokenize_side(final.strip())),
        symmetric_initial=symmetric_initial,
    )


def test_process_set_expands_builtin_p_for_charged_current_entries() -> None:
    process_set = enumerate_process_set("p p > w+ g")

    assert process_set.default_key == "u_dbar_to_g_wplus"
    assert [entry.process for entry in process_set.entries] == [
        "u d~ > g w+",
        "c s~ > g w+",
        "d~ u > g w+",
        "s~ c > g w+",
    ]


def test_process_enumerator_allows_explicit_charged_leptonic_current() -> None:
    enumeration = ProcessEnumerator().enumerate("u d~ > e+ ve g")

    assert enumeration.request.initial_state == ("u~", "d")
    assert enumeration.request.rest == ("e+", "ve", "g")
    assert len(enumeration.unique_processes) == 1
    assert enumeration.n_records >= 1


def test_process_set_reports_empty_and_malformed_inputs() -> None:
    with pytest.raises(ValueError, match="empty process"):
        split_process_set("d d~ > z g | ")
    with pytest.raises(ValueError, match="anonymous multiparticle label cannot be empty"):
        expand_process_variants("d d~ > z []")
    with pytest.raises(ValueError, match="unmatched"):
        expand_process_variants("d d~ > z [d g")


def test_process_set_skips_invalid_anonymous_variants() -> None:
    process_set = enumerate_process_set("d d~ > Z 3*[d g]")

    assert [entry.process for entry in process_set.entries] == [
        "d d~ > z g g g",
    ]
    assert process_set.default_key == "d_dbar_to_z_g_g_g"


def test_process_enumerator_allows_more_than_three_quark_lines_in_formula() -> None:
    enumerator = ProcessEnumerator()

    assert enumerator.expected_number_of_dual_amplitudes(
        ("d", "d~", "u", "u~", "s", "s~", "c", "c~", "g")
    ) == 96.0


def test_process_enumerator_uses_generic_three_quark_line_orders() -> None:
    enumeration = ProcessEnumerator().enumerate("d d~ > u u~ s s~")

    assert len(enumeration.unique_processes) == 1
    assert len(enumeration.groups) == 3
    assert enumeration.n_records == 12


def test_color_complete_reference_enumeration_keeps_crossed_quark_pairings() -> None:
    enumeration = ProcessEnumerator().enumerate_color_complete("d d~ > t t~")

    assert len(enumeration.groups) == 1
    assert enumeration.n_records == 4
    color_orders = {
        tuple(label + 1 for label in record.color_order)
        for record in enumeration.groups[0].records
    }
    assert color_orders == {
        (2, 1, 3, 4),
        (2, 4, 3, 1),
        (3, 1, 2, 4),
        (3, 4, 2, 1),
    }


def test_color_complete_reference_enumeration_keeps_all_gluon_words() -> None:
    enumeration = ProcessEnumerator().enumerate_color_complete("g g > g g")

    assert len(enumeration.groups) == 1
    assert enumeration.n_records == 6
    assert all(record.identical_factor == 2.0 for record in enumeration.groups[0].records)


def test_process_export_matches_legacy_process_list_for_small_qcd_case(
    tmp_path: Path,
) -> None:
    assert_export_matches_legacy_process_list("d d~ > z g g", tmp_path)


@pytest.mark.parametrize(
    "process",
    [
        "d d~ > a g",
        "u d~ > w+ g",
        "d d~ > z z g",
        "d d~ > e+ e- g",
    ],
)
def test_process_export_matches_legacy_for_non_z_gluon_families(
    tmp_path: Path,
    process: str,
) -> None:
    assert_export_matches_legacy_process_list(process, tmp_path)


def assert_export_matches_legacy_process_list(process: str, tmp_path: Path) -> None:
    repo_root = Path(__file__).resolve().parents[3]
    subprocess.run(
        [
            "python3",
            str(repo_root / "process_list.py"),
            "--serial",
            process,
        ],
        cwd=tmp_path,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    legacy_text = (tmp_path / "processes.txt").read_text()

    enumerator = ProcessEnumerator(ProcessOptions(serial=True))
    enumeration = enumerator.enumerate(process)
    export = tmp_path / "pyamplicol-processes.txt"
    enumerator.write_legacy_file(enumeration, export)

    assert export.read_text() == legacy_text
