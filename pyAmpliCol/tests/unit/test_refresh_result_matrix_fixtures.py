from __future__ import annotations

import importlib.util
import json
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
        "sources": [
            {"color_accuracy": "full", "marker": "archived"},
            {"color_accuracy": "lc", "marker": "archived"},
        ],
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
        "sources": [{"color_accuracy": "lc", "marker": "new"}],
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
    assert result["sources"] == [
        {"color_accuracy": "full", "marker": "archived"},
        {"color_accuracy": "lc", "marker": "new"},
    ]
    assert result["created_at"] == "new-metadata"
    assert result["summary"]["validated_cases"] == 2


def test_extend_fixture_requires_a_larger_multiplicity() -> None:
    refresh = _load_refresh_module()

    with pytest.raises(ValueError, match="requires --max-n to exceed"):
        refresh._extend_fixture(
            {"max_n": 4},
            {"max_n": 4},
        )


def test_build_fixture_min_n_does_not_read_retained_point_sources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    refresh = _load_refresh_module()
    selected_dir = tmp_path / "selected-flow"
    selected_dir.mkdir()
    point_path = selected_dir / "validation_momenta.json"
    point_path.write_text(
        json.dumps(
            {
                "points": [
                    [
                        {
                            "pdg": 1,
                            "momentum": ["1", "0", "0", "1"],
                        }
                    ]
                ]
            }
        ),
        encoding="utf-8",
    )

    def entry(
        n_final: int,
        point_source: str,
        selected_output_dir: Path | None = None,
    ) -> dict[str, object]:
        jit: dict[str, object] = {"status": "ok", "mode": "jit"}
        if selected_output_dir is not None:
            jit["selected_output_dir"] = str(selected_output_dir)
        return {
            "process": f"d d~ > z + {n_final - 1}g",
            "amplicol": {"status": "ok", "mode": "fortran"},
            "pyamplicol_jit": jit,
            "validation": {
                "status": "ok",
                "reference": 1.0,
                "values": {"pyamplicol_jit": 1.0},
                "relative_differences": {"pyamplicol_jit": 0.0},
                "point_source": point_source,
            },
        }

    cache_path = tmp_path / "result_matrix_data.json"
    cache_path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "base_processes": [{"key": "dd_z_jets", "process_id": 1}],
                "entries": {
                    "dd_z_jets": {
                        "4": entry(4, "stale-retained-point-source"),
                        "5": entry(
                            5,
                            "generic_validation_point",
                            selected_dir,
                        ),
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setitem(refresh.MATRIX_CACHES, "lc", cache_path)

    fixture = refresh._build_fixture(("lc",), min_n=5, max_n=5)

    assert [case["id"] for case in fixture["cases"]] == ["lc:dd_z_jets:n5"]
    assert fixture["cases"][0]["metadata"]["point_source"] == str(point_path)
    assert (
        fixture["cases"][0]["metadata"]["point_source_declared"]
        == "generic_validation_point"
    )
