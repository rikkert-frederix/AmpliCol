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

    assert "Gate status: **PASS_WITH_LIMITATIONS**." in report
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


def test_result_matrix_report_separates_low_n_mode_and_validation_coverage() -> None:
    audit = _load_audit_module()
    data = {
        "entries": {
            "proc": {
                "2": {
                    "amplicol": {"status": "ok"},
                    "pyamplicol_jit": {"status": "ok"},
                    "validation": {
                        "status": "ok",
                        "max_relative_difference": 1.0e-12,
                        "tolerance": 1.0e-8,
                    },
                },
                "3": {
                    "amplicol": {"status": "ok"},
                    "pyamplicol_jit": {"status": "ok"},
                },
                "6": {
                    "amplicol": {"status": "ok"},
                    "pyamplicol_jit": {"status": "failed"},
                    "validation": {"status": "failed"},
                },
            }
        }
    }

    report = audit.render_report(data, [])

    assert "Low-multiplicity coverage (`n <= 5`):" in report
    assert "`cases`=2" in report
    assert "`amplicol_ok`=2" in report
    assert "`jit_ok`=2" in report
    assert "`mode_failures`=0" in report
    assert "`validation_records`=1" in report
    assert "`validation_clean`=1" in report
    assert "`missing_validation_records`=1" in report
    assert "`validation_failures`=0" in report
