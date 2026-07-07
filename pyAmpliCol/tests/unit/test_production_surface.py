from __future__ import annotations

from pathlib import Path

import pytest

import pyamplicol


PROJECT_ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize(
    "relative_path",
    [
        "src/pyamplicol/generic_dag.py",
        "src/pyamplicol/generic_artifact.py",
        "src/pyamplicol/generic_stage_compiler.py",
        "src/pyamplicol/process_support.py",
        "src/pyamplicol/current_plan.py",
        "src/pyamplicol/matrix.py",
        "src/pyamplicol/symbolica_evaluator.py",
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
        "from .dag_runtime",
        "from .native",
        "from .evaluation",
        "from .tensor_runtime",
        "from .compiled_dag_runtime",
        "from .legacy_matrix",
        "LeadingColorZJets",
    )

    assert {
        token for token in forbidden_tokens if token in source
    } == set()


def test_package_root_does_not_export_retired_runtime_symbols() -> None:
    retired_symbols = {
        "LeadingColorZJetsNativeEvaluator",
        "NativeRuntimeEvaluator",
        "NativeMatrixElementGenerator",
        "ZGluonCompiledDAGEvaluator",
        "ZGluonDAGEvaluator",
        "ZGluonNumericTensorNetworkEvaluator",
        "ZGluonTensorNetworkEvaluator",
        "benchmark_z_gluon_modes",
        "profile_z_gluon_dag_evaluator",
        "profile_z_gluon_tensor_evaluator",
    }

    exported = set(pyamplicol.__all__)

    assert exported.isdisjoint(retired_symbols)
    assert {name for name in retired_symbols if hasattr(pyamplicol, name)} == set()
