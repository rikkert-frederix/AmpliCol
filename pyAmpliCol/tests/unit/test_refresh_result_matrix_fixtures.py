from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _load_refresh_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "scripts"
        / "refresh_result_matrix_fixtures.py"
    )
    spec = importlib.util.spec_from_file_location(
        "pyamplicol_refresh_result_matrix_fixtures",
        path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _case(case_id: str, n_final: int, marker: str) -> dict[str, object]:
    return {
        "id": case_id,
        "color_accuracy": "lc",
        "process_id": 1,
        "n_final": n_final,
        "marker": marker,
        "validation": {
            "relative_difference": 1.0e-12,
            "absolute_difference": 2.0e-15,
        },
    }


def _structural(case_id: str, n_final: int, marker: str) -> dict[str, object]:
    return {
        "id": case_id,
        "color_accuracy": "nlc",
        "process_id": 2,
        "n_final": n_final,
        "marker": marker,
    }


def test_extend_fixture_preserves_lower_multiplicity_records_exactly() -> None:
    refresh = _load_refresh_module()
    old_case = _case("lc:old:n4", 4, "archived")
    old_unsupported = _structural("nlc:unsupported:n4", 4, "archived")
    old_unvalidated = _structural("nlc:unvalidated:n4", 4, "archived")
    existing = {
        "max_n": 4,
        "cases": [old_case],
        "unsupported": [old_unsupported],
        "unvalidated": [old_unvalidated],
    }
    new_case = _case("lc:new:n5", 5, "new")
    new_unsupported = _structural("nlc:unsupported:n5", 5, "new")
    new_unvalidated = _structural("nlc:unvalidated:n5", 5, "new")
    refreshed = {
        "max_n": 5,
        "created_at": "new-metadata",
        "cases": [_case("lc:old:n4", 4, "recomputed"), new_case],
        "unsupported": [
            _structural("nlc:unsupported:n4", 4, "recomputed"),
            new_unsupported,
        ],
        "unvalidated": [
            _structural("nlc:unvalidated:n4", 4, "recomputed"),
            new_unvalidated,
        ],
    }

    result = refresh._extend_fixture(existing, refreshed)

    assert result["cases"] == [old_case, new_case]
    assert result["unsupported"] == [old_unsupported, new_unsupported]
    assert result["unvalidated"] == [old_unvalidated, new_unvalidated]
    assert result["created_at"] == "new-metadata"
    assert result["summary"]["validated_cases"] == 2


def test_extend_fixture_requires_a_larger_multiplicity() -> None:
    refresh = _load_refresh_module()

    with pytest.raises(ValueError, match="requires --max-n to exceed"):
        refresh._extend_fixture(
            {"max_n": 4},
            {"max_n": 4},
        )
