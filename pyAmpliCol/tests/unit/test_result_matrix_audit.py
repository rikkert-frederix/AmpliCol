from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_audit_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "docs"
        / "audit_result_matrix.py"
    )
    spec = importlib.util.spec_from_file_location("result_matrix_audit", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_result_matrix_gate_passes_with_documented_unsupported_cells() -> None:
    audit = _load_audit_module()
    data = {
        "entries": {
            "proc": {
                "2": {
                    "amplicol": {"status": "unsupported"},
                    "pyamplicol_jit": {"status": "backend_unsupported"},
                    "validation": {"status": "ok", "max_relative_difference": 0.0},
                }
            }
        }
    }

    report = audit.render_report(data, [])

    assert "Gate status: **PASS**." in report
    assert "`amplicol_unsupported`=1" in report
    assert "`jit_backend_unsupported`=1" in report
    assert "`missing_cpp_o3`=1" in report


def test_result_matrix_gate_fails_on_missing_required_jit_cell() -> None:
    audit = _load_audit_module()
    data = {
        "entries": {
            "proc": {
                "2": {
                    "amplicol": {"status": "ok"},
                    "validation": {"status": "ok", "max_relative_difference": 0.0},
                }
            }
        }
    }

    report = audit.render_report(data, [])

    assert "Gate status: **FAIL**." in report
    assert "`missing_jit`=1" in report
