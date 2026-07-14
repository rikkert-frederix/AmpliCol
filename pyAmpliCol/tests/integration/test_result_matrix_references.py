from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


FIXTURE_PATH = (
    Path(__file__).resolve().parents[1]
    / "fixtures"
    / "result_matrix_references.json"
)


def _load_fixture() -> dict[str, Any]:
    return json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def test_result_matrix_reference_fixture_is_complete_through_n6() -> None:
    fixture = _load_fixture()

    assert fixture["schema_version"] == 1
    assert fixture["kind"] == "pyamplicol-result-matrix-reference-fixture"
    assert fixture["max_n"] == 6
    assert fixture["summary"]["validated_cases"] == 144
    assert fixture["summary"]["structural_unsupported_cases"] == 1
    assert fixture["summary"]["unvalidated_pyamplicol_cases"] == 4

    counts = Counter(case["color_accuracy"] for case in fixture["cases"])
    assert counts == {"lc": 58, "nlc": 43, "full": 43}

    unsupported = Counter(item["color_accuracy"] for item in fixture["unsupported"])
    assert unsupported == {"lc": 1}
    unvalidated = Counter(item["color_accuracy"] for item in fixture["unvalidated"])
    assert unvalidated == {"nlc": 2, "full": 2}


def test_result_matrix_reference_fixture_values_match_fortran_references() -> None:
    fixture = _load_fixture()

    seen_ids: set[str] = set()
    for case in fixture["cases"]:
        assert case["id"] not in seen_ids
        seen_ids.add(case["id"])
        assert case["process"]
        assert case["color_accuracy"] in {"lc", "nlc", "full"}
        assert 1 <= int(case["n_final"]) <= 6
        assert case["fortran"]["status"] == "ok"
        assert case["pyamplicol"]["status"] == "ok"
        if case["color_accuracy"] == "lc":
            assert case["fortran"]["reference_probe"] == "direct_generated_library_benchmark"
        else:
            assert str(case["fortran"]["reference_probe"]).startswith(
                "generated_library_color_probe_raw"
            )

        reference = float(case["fortran"]["value"])
        value = float(case["pyamplicol"]["value"])
        rel_diff = abs(value - reference) / max(abs(value), abs(reference), 1.0e-300)
        tolerance = float(case["validation"]["tolerance"])

        assert math.isfinite(reference)
        assert math.isfinite(value)
        assert rel_diff <= tolerance
        assert rel_diff == case["validation"]["relative_difference"]
        assert case["validation"]["max_relative_difference"] <= tolerance
        assert case["phase_space_point"]
        for particle in case["phase_space_point"]:
            assert isinstance(particle["pdg"], int)
            assert len(particle["momentum"]) == 4


def test_result_matrix_reference_fixture_documents_unvalidated_colour_classes() -> None:
    fixture = _load_fixture()

    unsupported_ids = {item["id"] for item in fixture["unsupported"]}
    assert unsupported_ids == {"lc:dd_4q_lines:n6"}
    unsupported = fixture["unsupported"][0]
    assert unsupported["kind"] == "structural-unsupported"
    assert unsupported["pyamplicol_jit_status"] == "ok"
    assert "quarks" in unsupported["reason"].lower()
    unvalidated_ids = {item["id"] for item in fixture["unvalidated"]}
    assert unvalidated_ids == {
        "nlc:dd_3q_lines:n4",
        "full:dd_3q_lines:n4",
        "nlc:dd_3q_lines:n5",
        "full:dd_3q_lines:n5",
    }
    for item in fixture["unvalidated"]:
        assert str(item["process"]).startswith("d d~ > u u~ s s~")
        assert item["kind"] == "fortran-reference-unavailable"
        assert item["pyamplicol_jit_status"] == "ok"
        assert "quark" in item["reason"].lower()
