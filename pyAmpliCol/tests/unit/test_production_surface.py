from __future__ import annotations

from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    "relative_path",
    [
        "src/pyamplicol/generic_dag.py",
        "src/pyamplicol/generic_artifact.py",
        "src/pyamplicol/generic_stage_compiler.py",
        "src/pyamplicol/process_support.py",
        "src/pyamplicol/current_plan.py",
    ],
)
def test_schema_v2_production_modules_do_not_reference_family_runtimes(
    relative_path: str,
) -> None:
    source = (PROJECT_ROOT / relative_path).read_text(encoding="utf-8")

    forbidden_tokens = (
        "ZGluon",
        "_z_gluon",
        "z_gluon_",
        "tensor_runtime",
        "compiled_dag_runtime",
        "legacy_matrix",
        "LeadingColorZJets",
    )

    assert {
        token for token in forbidden_tokens if token in source
    } == set()
