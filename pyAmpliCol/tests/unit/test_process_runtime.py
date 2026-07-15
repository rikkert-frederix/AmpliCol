from __future__ import annotations

import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from pyamplicol.core_types import NativeEvaluationError
from pyamplicol.process_runtime import (
    PythonProcessRuntime,
    load_process,
    load_process_manifest,
)


def _write_process_manifest(
    root: Path,
    *,
    schema_version: int = 2,
    kind: str = "pyamplicol-generic-dag-process",
) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "process_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": schema_version,
                "kind": kind,
                "artifact_class": "generic-dag-schema-v2",
                "key": "d_dbar_to_z",
                "process": "d d~ > z",
                "color_accuracy": "lc",
                "external_pdg_order": [1, -1, 23],
                "compiled": {"kind": "generic-dag-stage-evaluators"},
            }
        ),
        encoding="utf-8",
    )


def test_schema_v2_process_manifest_loads(tmp_path: Path) -> None:
    _write_process_manifest(tmp_path)

    manifest = load_process_manifest(tmp_path)

    assert manifest.schema_version == 2
    assert manifest.process == "d d~ > z"
    assert manifest.external_pdg_order == (1, -1, 23)
    assert manifest.compiled_kind == "generic-dag-stage-evaluators"


def test_schema_v1_process_manifest_is_rejected(tmp_path: Path) -> None:
    _write_process_manifest(
        tmp_path,
        schema_version=1,
        kind="pyamplicol-z-gluon-process",
    )

    with pytest.raises(NativeEvaluationError, match="schema-v1.*unsupported"):
        load_process_manifest(tmp_path)


def test_schema_v1_process_set_is_rejected(tmp_path: Path) -> None:
    (tmp_path / "process_set_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pyamplicol-process-set",
                "default_process_key": "d_dbar_to_z",
                "processes": [],
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(NativeEvaluationError, match="schema-v1.*unsupported"):
        load_process_manifest(tmp_path)


def test_python_facade_delegates_metadata_model_card_and_atomic_setters(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    _write_process_manifest(tmp_path)
    model_card = tmp_path / "parameters.json"
    model_card.write_text("{}\n", encoding="utf-8")
    loaded_with: list[object] = []

    class FakeHandle:
        process = "d d~ > z"
        physics = object()

        def __init__(self) -> None:
            self.parameter_updates: list[dict[str, object]] = []

        def metadata(self) -> dict[str, object]:
            return {"schema_version": 2, "process": self.process}

        def set_model_parameters(self, parameters: dict[str, object]) -> None:
            self.parameter_updates.append(parameters)

    handle = FakeHandle()

    class FakeRuntime:
        @staticmethod
        def load(*arguments: object) -> FakeHandle:
            loaded_with.extend(arguments)
            return handle

    monkeypatch.setitem(
        sys.modules,
        "rusticol",
        SimpleNamespace(Runtime=FakeRuntime),
    )

    runtime = load_process(tmp_path, model_parameters=model_card)

    assert isinstance(runtime, PythonProcessRuntime)
    assert runtime.metadata == {"schema_version": 2, "process": "d d~ > z"}
    assert loaded_with == [str(tmp_path), None, str(model_card)]
    runtime.set_model_parameter("normalization.alpha_s_me_check", 0.101)
    assert handle.parameter_updates == [
        {"normalization.alpha_s_me_check": 0.101}
    ]
