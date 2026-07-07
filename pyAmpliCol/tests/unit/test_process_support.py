from __future__ import annotations

import pytest

from pyamplicol.process_support import (
    classify_process_support,
    describe_process_content,
)


def test_process_content_describes_generic_quark_lines_and_singlets() -> None:
    content = describe_process_content("d d~ > z z g")

    assert content.quark_pair_count == 1
    assert content.gluon_count == 1
    assert content.final_singlets == ("z", "z")
    assert content.final_vectors == ("z", "z")
    assert content.final_higgs_count == 0


@pytest.mark.parametrize(
    "process",
    [
        "d d~ > z",
        "d d~ > z g g",
        "d d~ > e+ e-",
        "d d~ > e+ e- g",
        "d d~ > e+ e- g g",
        "d d~ > e+ e- a",
        "u d~ > e+ ve g",
        "u d~ > e+ ve",
        "d u~ > w- g",
        "d u~ > e- ve~ g",
        "u d~ > w+ g g",
        "d d~ > u u~",
        "d d~ > u u~ g",
        "d d~ > u u~ s s~",
        "g g > u u~ d d~",
        "g g > g g g",
        "g g > g g",
        "g g > u u~",
        "u d~ > w+ g",
        "u d~ > w+ z",
        "u d~ > w+ a",
        "d d~ > a z g",
        "d d~ > a a",
        "d d~ > a a g",
        "d d~ > z z g",
        "d d~ > z z g g",
        "d d~ > z z z",
        "d d~ > a z z",
        "d d~ > h h z",
        "d d~ > w+ w- g",
        "d d~ > w+ w-",
        "d d~ > h z",
        "u d~ > h w+",
    ],
)
def test_process_support_reports_generic_dag_ready_schema_v2_supported(
    process: str,
) -> None:
    report = classify_process_support(process)

    assert report.runtime_artifact_supported is True
    assert report.generic_dag_runtime_supported is True
    assert report.support_class == "generic-dag-schema-v2"
    assert report.missing_feature is None
    assert report.artifact_unavailable_message is None
    assert report.color_plan is not None
    assert report.color_plan.ready_for_leading_colour is True
    assert report.color_plan.sector_count > 0
    assert report.current_plan is not None
    assert report.current_plan.full_tensor_network_ready is True
    assert report.current_plan.has_closure is True
    assert report.current_plan.pending_vertex_kinds == ()
    assert report.current_plan.unimplemented_vertex_kinds == ()
    payload = report.to_json_dict()
    assert payload["support_class"] == "generic-dag-schema-v2"
    assert payload["generic_dag_runtime_supported"] is True
    assert "family" not in payload
    assert "artifact_supported" not in payload
    assert "native_supported" not in payload


@pytest.mark.parametrize(
    "process,missing,message",
    [
        (
            "d d~ > h g",
            "generic-dag-amplitude-roots",
            "did not find amplitude closures",
        ),
        (
            "d d~ > ve ve~ g",
            "generic-dag-amplitude-roots",
            "did not find amplitude closures",
        ),
        (
            "d d~ > vm vm~ g g",
            "generic-dag-amplitude-roots",
            "did not find amplitude closures",
        ),
    ],
)
def test_process_support_reports_generic_preflight_blockers(
    process: str,
    missing: str,
    message: str,
) -> None:
    report = classify_process_support(process)

    assert report.runtime_artifact_supported is False
    assert report.missing_feature == missing
    assert report.artifact_unavailable_message is not None
    assert message in report.artifact_unavailable_message


def test_process_support_reports_top_yukawa_as_generic_dag_supported() -> None:
    report = classify_process_support("t t~ > h")

    assert report.runtime_artifact_supported is True
    assert report.generic_dag_runtime_supported is True
    assert report.support_class == "generic-dag-schema-v2"
    assert report.missing_feature is None
    assert report.current_plan is not None
    assert report.current_plan.has_closure is True
    assert 16 in report.current_plan.ready_vertex_kinds
    assert report.current_plan.unimplemented_vertex_kinds == ()
    assert report.artifact_unavailable_message is None


def test_process_support_includes_generic_current_plan_for_multi_singlets() -> None:
    report = classify_process_support("d d~ > z z g")

    assert report.current_plan is not None
    assert report.current_plan.has_closure is True
    assert report.current_plan.unimplemented_vertex_kinds == ()
    assert report.current_plan.pending_vertex_kinds == ()
    assert report.current_plan.full_tensor_network_ready is True
    assert report.artifact_unavailable_message is None
    payload = report.to_json_dict()
    assert payload["current_plan"]["closure_count"] == report.current_plan.closure_count
    assert payload["current_plan"]["color_sector_count"] == 1


def test_process_support_includes_generic_current_plan_for_multi_quark_lines() -> None:
    report = classify_process_support("d d~ > u u~")

    assert report.color_plan is not None
    assert report.color_plan.sector_count == 4
    assert report.color_plan.sector_kind_counts == (("open-lines", 4),)
    assert report.current_plan is not None
    assert report.current_plan.has_closure is True
    assert report.current_plan.color_sectors == (0, 1, 2, 3)
    assert report.current_plan.unimplemented_vertex_kinds == ()
    assert report.current_plan.pending_vertex_kinds == ()
    assert report.artifact_unavailable_message is None


def test_process_support_honors_generic_reference_sector_pruning() -> None:
    report = classify_process_support(
        "d d~ > u u~ s s~ c c~",
        current_plan_max_currents=2000,
        selected_color_sector_ids={0},
        max_coupling_orders={"QCD": 6},
    )

    assert report.color_plan is not None
    assert report.color_plan.sector_count == 576
    assert report.current_plan is not None
    assert report.current_plan.truncated is False
    assert report.current_plan.color_sectors == (0,)
    assert report.current_plan.has_closure is True
    assert report.runtime_artifact_supported is True
    assert report.artifact_unavailable_message is None


def test_process_support_accepts_selected_sector_when_full_color_plan_truncates() -> None:
    report = classify_process_support(
        "d d~ > u u~ s s~ c c~ b b~ t t~",
        color_plan_max_sectors=20000,
        current_plan_max_currents=500000,
        selected_color_sector_ids={0},
        max_quark_pairs=6,
    )

    assert report.color_plan is not None
    assert report.color_plan.truncated is True
    assert report.current_plan is not None
    assert report.current_plan.color_sectors == (0,)
    assert report.current_plan.full_tensor_network_ready is True
    assert report.runtime_artifact_supported is True
    assert report.artifact_unavailable_message is None


def test_process_support_reports_colour_plan_truncation_without_hiding_current_plan() -> None:
    report = classify_process_support(
        "d d~ > u u~ s s~",
        color_plan_max_sectors=3,
    )

    assert report.color_plan is not None
    assert report.color_plan.truncated is True
    assert report.color_plan.sector_count == 3
    assert report.current_plan is not None
    assert report.artifact_unavailable_message is not None
    assert "colour planning truncated at 3 sectors" in report.artifact_unavailable_message
    assert "runtime execution are available" in report.artifact_unavailable_message


def test_process_support_can_skip_generic_current_plan_for_fast_diagnostics() -> None:
    report = classify_process_support(
        "d d~ > z z g",
        include_current_plan=False,
    )

    assert report.current_plan is None
    assert report.color_plan is not None
    assert report.artifact_unavailable_message is not None
    assert "vertex kinds" not in report.artifact_unavailable_message


def test_process_support_rejects_colour_expansion_without_lc_fallback() -> None:
    report = classify_process_support("d d~ > z g", color_accuracy="full")

    assert report.runtime_artifact_supported is False
    assert report.color_accuracy == "full"
    assert report.support_class == "generic-dag-colour-preflight"
    assert report.missing_feature == "colour-expansion"
    assert report.artifact_unavailable_message is not None
    assert "--color-accuracy=full" in report.artifact_unavailable_message
    assert "Idenso basis/metric" in report.artifact_unavailable_message
