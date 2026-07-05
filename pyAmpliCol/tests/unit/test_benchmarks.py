from __future__ import annotations

from pyamplicol.benchmarks import (
    _legacy_runtime_per_point,
    _legacy_runtime_per_point_error,
    _time_shared_dag_evaluations,
    format_mode_benchmark_table,
    summarize_mode_benchmark,
)
from pyamplicol.dag_runtime import DAGEvaluationTiming
from pyamplicol.reference import TimingRow


def test_mode_benchmark_summary_and_table_use_max_relative_difference() -> None:
    rows = [
        {
            "gluon_count": 1,
            "process": "d d~ > z g",
            "reference_matrix_elements": [1.0],
            "legacy": {"status": "ok", "generation_s": 2.0, "runtime_s_per_point": 1e-6},
            "python": {
                "status": "ok",
                "generation_s": 0.1,
                "runtime_s_per_point": 2e-4,
                "max_relative_difference_to_legacy": 1e-12,
            },
            "numeric_tn": {
                "status": "ok",
                "generation_s": 0.2,
                "runtime_s_per_point": 3e-2,
                "max_relative_difference_to_legacy": 3e-12,
            },
            "parametric_tn": {
                "status": "ok",
                "generation_s": 0.3,
                "runtime_s_per_point": 4e-3,
                "max_relative_difference_to_legacy": 2e-12,
            },
            "shared_dag": {
                "status": "ok",
                "generation_s": 0.4,
                "runtime_s_per_point": 5e-4,
                "runtime_evaluator_only_s_per_point": 5e-6,
                "max_relative_difference_to_legacy": 4e-12,
            },
        }
    ]

    summary = summarize_mode_benchmark(rows)
    table = format_mode_benchmark_table({"rows": rows})

    assert summary["all_four_modes_match_for_all_rows"] is True
    assert summary["max_relative_difference_to_legacy"] == 4e-12
    assert "4e-12" in table
    assert "5e-06" in table


def test_legacy_runtime_uses_probe_point_count() -> None:
    rows = (TimingRow("amplitude evaluation", 0.125, ""),)

    assert _legacy_runtime_per_point(rows, 1000) == 0.000125
    assert _legacy_runtime_per_point_error(rows, 1000) == 0.5e-6


def test_shared_dag_profile_warms_full_batch() -> None:
    class DummyEvaluator:
        batch_size = 4

        def __init__(self) -> None:
            self.call_sizes: list[int] = []
            self.last_runtime_timing = DAGEvaluationTiming(evaluator_time_s=0.0)

        def evaluate_matrix_elements_many(self, particles):
            self.call_sizes.append(len(particles))
            self.last_runtime_timing = DAGEvaluationTiming(evaluator_time_s=0.0)
            return [0.0 for _ in particles]

    evaluator = DummyEvaluator()
    particles = tuple(() for _ in range(10))

    _time_shared_dag_evaluations(evaluator, particles)

    assert evaluator.call_sizes[0] == evaluator.batch_size
    assert evaluator.call_sizes[1:] == [4, 4, 2]
