from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


def _load_audit_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "docs"
        / "low_n_pruning_audit.py"
    )
    spec = importlib.util.spec_from_file_location("low_n_pruning_audit", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_low_n_pruning_audit_parses_manifest_and_fortran_counts(tmp_path: Path) -> None:
    audit = _load_audit_module()
    manifest = tmp_path / "generic_process_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "process": "g g > t t~ g",
                "currents": [
                    {
                        "id": 0,
                        "index": {
                            "external_labels": [1],
                        },
                    },
                    {
                        "id": 1,
                        "index": {
                            "external_labels": [1, 2],
                        },
                    },
                ],
                "interactions": [
                    {
                        "vertex_kind": 0,
                        "result_id": 1,
                    },
                    {
                        "vertex_kind": 6,
                        "result_id": 1,
                    },
                ],
                "amplitude_roots": [{"id": 0}],
            }
        ),
        encoding="utf-8",
    )
    fortran = tmp_path / "amp1_1_lib.f03"
    fortran.write_text(
        """
module amp1_1_lib
contains
  subroutine evaluate_amp1_1(p,amps)
    complex(kind=8),dimension(1:6,2) :: val_c
    complex(kind=8),dimension(1:6,3) :: int_c
  end subroutine evaluate_amp1_1
  subroutine vertex_type2_0(pp,val_c,int_c)
    do i=1, 2
    enddo
  end subroutine vertex_type2_0
  subroutine vertex_type2_6(pp,val_c,int_c)
    do i=1, 1
    enddo
  end subroutine vertex_type2_6
  subroutine compute_amps(amps,val_c)
    amps(1)=val_c(1,1)
  end subroutine compute_amps
end module amp1_1_lib
""",
        encoding="utf-8",
    )

    py_counts = audit.pyamplicol_counts(manifest)
    ft_counts = audit.fortran_counts(fortran)

    assert py_counts.currents == 2
    assert py_counts.interactions == 2
    assert py_counts.amplitude_roots == 1
    assert py_counts.stage_kind_counts == {(2, 0): 1, (2, 6): 1}
    assert ft_counts.currents == 2
    assert ft_counts.interactions == 3
    assert ft_counts.amplitude_roots == 1
    assert ft_counts.stage_kind_counts == {(2, 0): 2, (2, 6): 1}
    report = audit.markdown_report(py_counts, ft_counts)
    assert "| delta | +0 | -1 | +0 |" in report
