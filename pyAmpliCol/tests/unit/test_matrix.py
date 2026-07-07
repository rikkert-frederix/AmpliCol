from __future__ import annotations

import json
from pathlib import Path

import pytest

from pyamplicol.core_types import NativeEvaluationError
from pyamplicol.evaluation import NativeRuntimeEvaluator
from pyamplicol.matrix import MatrixElementGenerator


def test_matrix_generator_uses_generic_dag_schema_without_family_graph(
    tmp_path: Path,
) -> None:
    result = MatrixElementGenerator(cache_dir=tmp_path).generate("d d~ > z z g")

    assert result.backend == "generic-dag-schema-v2"
    assert result.runtime_artifact_supported is True
    assert result.generic_dag_runtime_supported is True
    assert result.manifest is not None
    assert result.manifest.dag.has_amplitudes is True
    assert result.manifest_file is not None
    assert result.manifest_file.exists()
    assert result.cache_file is not None
    payload = json.loads(result.cache_file.read_text())
    assert payload["backend"] == "generic-dag-schema-v2"
    assert "supported_native_target" not in payload
    assert payload["manifest"]["kind"] == "pyamplicol-generic-process-plan"
    assert payload["support_report"]["support_class"] == "generic-dag-schema-v2"


def test_matrix_generator_handles_non_z_generic_process() -> None:
    result = MatrixElementGenerator().generate(
        "d d~ > e+ e- g",
        write_cache_metadata=False,
    )

    assert result.backend == "generic-dag-schema-v2"
    assert result.manifest is not None
    assert result.manifest.process == "d d~ > e+ e- g"


def test_native_runtime_evaluator_is_reference_only_by_default() -> None:
    with pytest.raises(NativeEvaluationError, match="retired reference-only runtime"):
        NativeRuntimeEvaluator("d d~ > z")


def test_native_runtime_evaluator_reference_opt_in_remains_available() -> None:
    evaluator = NativeRuntimeEvaluator(
        "d d~ > z",
        runtime_backend="python",
        allow_reference_legacy=True,
    )

    assert evaluator.metadata.backend == "native-python-zero-gluon"
