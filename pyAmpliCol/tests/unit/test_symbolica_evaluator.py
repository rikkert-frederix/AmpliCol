from __future__ import annotations

import platform
import threading

import numpy as np
import pytest
from symbolica import E, S

from pyamplicol.symbolica_evaluator import (
    SymbolicaEvaluatorSettings,
    _ChunkedSymbolicaEvaluator,
    _compile_symbolica_outputs,
    _jit_payload_layout,
    _symbolica_evaluator_artifact_manifest,
)


def test_jit_artifact_reports_materialization_progress(tmp_path) -> None:
    x = S("x")
    events: list[dict[str, object]] = []
    evaluator = _compile_symbolica_outputs(
        (x * x + 1,),
        [x],
        merge_evaluators_strategy=False,
        verbose_evaluator_build=False,
        real_params=(0,),
        symbolica_settings=SymbolicaEvaluatorSettings(
            iterations=1,
            n_cores=1,
            jit_optimization_level=0,
        ),
        jit_compile=True,
        label="progress_probe",
        progress_callback=events.append,
    )

    manifest = _symbolica_evaluator_artifact_manifest(evaluator, tmp_path)

    stages = [str(event["stage"]) for event in events]
    assert "jit materialize" in stages
    assert "jit materialized" in stages
    assert stages.index("jit materialize") < stages.index("jit materialized")
    assert manifest["kind"] == "jit-symbolica-evaluator"
    assert manifest["jit_payload_layout"] == _jit_payload_layout()
    assert (tmp_path / manifest["evaluator_state_path"]).is_file()


def test_symbolica_output_compiler_forwards_function_definitions() -> None:
    function, parameter, formal = S("compiled_function", "parameter", "formal")
    evaluator = _compile_symbolica_outputs(
        (function(parameter),),
        [parameter],
        merge_evaluators_strategy=False,
        verbose_evaluator_build=False,
        functions={(function, (formal,)): formal**2 + 1},
        real_params=(0,),
        symbolica_settings=SymbolicaEvaluatorSettings(
            iterations=1,
            n_cores=1,
            jit_optimization_level=0,
        ),
        jit_compile=True,
        label="function_definition_probe",
    )

    result = evaluator.evaluate(np.array([[2.0], [3.0]], dtype=np.float64))

    np.testing.assert_allclose(result, np.array([[5.0], [10.0]]))


def test_chunk_artifacts_materialize_concurrently_in_stable_order(tmp_path) -> None:
    barrier = threading.Barrier(3, timeout=2.0)

    class FakeEvaluator:
        settings = {"compiled_chunk_compile_workers": 3}

        def __init__(self, index: int) -> None:
            self.index = index

        def artifact_manifest(self, artifact_dir) -> dict[str, object]:
            assert artifact_dir == tmp_path
            barrier.wait()
            return {
                "kind": "fake-evaluator",
                "index": self.index,
                "build_timing": {"artifact_manifest_s": 1.0},
            }

    evaluator = _ChunkedSymbolicaEvaluator(
        tuple(FakeEvaluator(index) for index in range(3))
    )

    manifest = evaluator.artifact_manifest(tmp_path)

    assert [chunk["index"] for chunk in manifest["chunks"]] == [0, 1, 2]
    assert manifest["build_timing"]["artifact_manifest_s"] == 3.0


def test_real_jit_chunks_materialize_concurrently(tmp_path) -> None:
    x = S("x")
    evaluator = _compile_symbolica_outputs(
        (x + 1, x * x + 2),
        [x],
        merge_evaluators_strategy=False,
        verbose_evaluator_build=False,
        real_params=(0,),
        symbolica_settings=SymbolicaEvaluatorSettings(
            iterations=1,
            n_cores=1,
            jit_optimization_level=0,
            compiled_output_chunk_size=1,
            compiled_chunk_compile_workers=2,
        ),
        jit_compile=True,
        label="parallel_materialization_probe",
    )

    manifest = _symbolica_evaluator_artifact_manifest(evaluator, tmp_path)

    assert manifest["kind"] == "chunked-symbolica-evaluator"
    assert len(manifest["chunks"]) == 2
    assert all(
        (tmp_path / chunk["evaluator_state_path"]).is_file()
        for chunk in manifest["chunks"]
    )


@pytest.mark.skipif(
    platform.machine().lower() not in {"arm64", "aarch64"},
    reason="the external-call spill-offset regression is AArch64-specific",
)
@pytest.mark.parametrize("optimization_level", (1, 3))
def test_real_and_complex_jit_external_calls_on_aarch64(
    optimization_level: int,
) -> None:
    x = S("x")
    evaluator = E("0").evaluator_multiple(
        (E("sin(x) + cos(x)"), E("exp(x) / (1 + x^2)")),
        (x,),
        iterations=1,
        n_cores=1,
        jit_compile=True,
        jit_optimization_level=optimization_level,
    )

    complex_rows = np.array([[0.25 + 0.1j], [0.5 - 0.2j]], dtype=np.complex128)
    complex_expected = np.column_stack(
        (
            np.sin(complex_rows[:, 0]) + np.cos(complex_rows[:, 0]),
            np.exp(complex_rows[:, 0]) / (1 + complex_rows[:, 0] ** 2),
        )
    )
    np.testing.assert_allclose(
        evaluator.evaluate_complex(complex_rows),
        complex_expected,
        rtol=2e-13,
        atol=2e-13,
    )

    real_rows = np.array([[0.25], [0.5]], dtype=np.float64)
    real_expected = np.column_stack(
        (
            np.sin(real_rows[:, 0]) + np.cos(real_rows[:, 0]),
            np.exp(real_rows[:, 0]) / (1 + real_rows[:, 0] ** 2),
        )
    )
    np.testing.assert_allclose(
        evaluator.evaluate(real_rows),
        real_expected,
        rtol=2e-13,
        atol=2e-13,
    )


def test_chunk_profile_counts_shared_input_pack_once() -> None:
    class FakeEvaluator:
        def __init__(self, value: float, profile: tuple[float, float, float]) -> None:
            self.value = value
            self.profile = profile

        def supports_complex_profiled(self) -> bool:
            return True

        def _evaluate_complex_profiled_prepared(self, rows):
            return np.full((len(rows), 1), self.value), self.profile

    evaluator = _ChunkedSymbolicaEvaluator(
        (
            FakeEvaluator(1.0, (4.0, 5.0, 1.0)),
            FakeEvaluator(2.0, (2.0, 7.0, 2.0)),
            FakeEvaluator(3.0, (3.0, 11.0, 4.0)),
        )
    )

    outputs, profile = evaluator.evaluate_complex_profiled(
        np.ones((2, 1), dtype=np.complex128)
    )

    assert [output[0, 0] for output in outputs] == [1.0, 2.0, 3.0]
    assert profile == (3.0, 23.0, 7.0)


@pytest.mark.skipif(
    platform.machine().lower() not in {"arm64", "aarch64"},
    reason="native two-lane profiling is AArch64-specific",
)
def test_real_jit_chunk_profile_matches_complex_evaluation() -> None:
    x = S("x")
    evaluator = _compile_symbolica_outputs(
        (x + 1, x * x + 2),
        [x],
        merge_evaluators_strategy=False,
        verbose_evaluator_build=False,
        real_params=(0,),
        symbolica_settings=SymbolicaEvaluatorSettings(
            iterations=1,
            n_cores=1,
            jit_optimization_level=1,
            compiled_output_chunk_size=1,
        ),
        jit_compile=True,
        label="native_profile_probe",
    )
    rows = np.array([[0.5], [1.25]], dtype=np.complex128)

    expected = evaluator.evaluate_complex(rows)
    chunks, profile = evaluator.evaluate_complex_profiled(rows)

    assert evaluator.supports_complex_profiled()
    assert np.allclose(np.concatenate(chunks, axis=1), expected)
    assert all(value >= 0.0 for value in profile)
