from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pytest


def _load_model_results_module():
    path = Path(__file__).resolve().parents[2] / "docs" / "model_results.py"
    spec = importlib.util.spec_from_file_location("model_results_docs", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_scalar_table_renders_values_without_control_character_escapes() -> None:
    model_results = _load_model_results_module()
    data = {
        "models": {
            "scalars": {
                "entries": {
                    "2": {
                        "status": "ok",
                        "process": "scalar_0 scalar_0 > scalar_0 scalar_0",
                        "generation_s": 0.2,
                        "jit_compile_s": 0.01,
                        "dag": {
                            "currents": 8,
                            "interactions": 6,
                            "roots": 1,
                        },
                        "runtime": {
                            "pure_evaluator_us_per_point": 0.025,
                            "wall_us_per_point": 0.1,
                        },
                        "value": {
                            "double": 0.5,
                            "relative_difference": 0.0,
                        },
                    }
                },
                "source_measurements": {},
            }
        }
    }

    table = model_results.render_model_table(
        data,
        model_results.MODEL_PROFILES["scalars"],
    )

    assert "\t" not in table
    assert r"\texttt{0.5}" in table
    assert r"\phi_0\,\phi_0\to \phi_0\,\phi_0" in table
    assert r"eval [$\mu$s]" in table
    assert "0.025 & 0.1" in table
    assert "--max-coupling-order QCD=1" in table


def test_full_scalar_profile_keeps_all_tree_diagrams(tmp_path: Path) -> None:
    model_results = _load_model_results_module()
    contact_profile = model_results.MODEL_PROFILES["scalars"]
    full_profile = model_results.MODEL_PROFILES["scalars-full"]
    process = "scalar_0 scalar_0 > scalar_0 scalar_0 scalar_0"

    contact_command = model_results._generation_command(
        contact_profile,
        3,
        process=process,
        output_dir=tmp_path / "contact",
        model_artifact=tmp_path / "scalars.pyAmplicol-model.json",
    )
    full_command = model_results._generation_command(
        full_profile,
        3,
        process=process,
        output_dir=tmp_path / "full",
        model_artifact=tmp_path / "scalars.pyAmplicol-model.json",
    )

    assert full_profile.process(3) == process
    assert "--max-coupling-order" in contact_command
    assert "--max-coupling-order" not in full_command
    assert full_profile.expected_value(3) is None
    table = model_results.render_model_table(
        {},
        full_profile,
    )
    assert "Complete scalar trees" in table
    assert "seven-scalar entry is allowed one hour" in table


def test_generation_failure_records_generation_phase(monkeypatch) -> None:
    model_results = _load_model_results_module()
    data = {}

    def fail_generation(_command):
        raise RuntimeError("generation failed")

    monkeypatch.setattr(model_results, "_run_json", fail_generation)
    with pytest.raises(RuntimeError, match="generation failed"):
        model_results._run_case(
            data,
            model_results.MODEL_PROFILES["scalars-full"],
            7,
            output_root=model_results.PROJECT_ROOT / ".test_model_results_outputs",
            target_runtime_s=10.0,
            generation_timeout_s=300,
            timing_timeout_s=600,
            reuse_existing=False,
            generation_s=None,
            jit_compile_s=None,
        )

    entry = data["models"]["scalars_full"]["entries"]["7"]
    assert entry["status"] == "error"
    assert entry["failure_phase"] == "generation"


def test_scalar_gravity_timeout_label_uses_entry_limit() -> None:
    model_results = _load_model_results_module()
    table = model_results.render_model_table(
        {
            "models": {
                "scalar_gravity": {
                    "entries": {
                        "4": {
                            "status": "timeout",
                            "generation_timeout_s": 3600,
                        }
                    },
                    "source_measurements": {},
                }
            }
        },
        model_results.MODEL_PROFILES["scalar-gravity"],
    )

    assert r"t/o $>60$ min" in table


def test_scalar_gravity_timing_timeout_is_labeled_separately() -> None:
    model_results = _load_model_results_module()
    table = model_results.render_model_table(
        {
            "models": {
                "scalar_gravity": {
                    "entries": {
                        "4": {
                            "status": "timeout",
                            "failure_phase": "timing",
                            "generation_timeout_s": 3600,
                            "timing_timeout_s": 600,
                        }
                    },
                    "source_measurements": {},
                }
            }
        },
        model_results.MODEL_PROFILES["scalar-gravity"],
    )

    assert r"run t/o $>10$ min" in table
