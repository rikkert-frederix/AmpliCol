from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


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
                        "runtime": {"wall_us_per_point": 0.1},
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
