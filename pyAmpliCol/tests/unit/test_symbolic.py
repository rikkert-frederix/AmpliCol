from __future__ import annotations

import math

from pyamplicol.native import ExternalMomentum, LeadingColorZJetsNativeEvaluator
from pyamplicol.symbolic import ZeroGluonSymbolicEvaluator


def test_zero_gluon_symbolica_evaluator_matches_native_kernel() -> None:
    native_evaluator = LeadingColorZJetsNativeEvaluator()
    particles = native_evaluator.canonical_zero_gluon_point("d d~ > z")

    native = native_evaluator.evaluate("d d~ > z", particles=particles)
    symbolic = ZeroGluonSymbolicEvaluator().evaluate("d d~ > z", particles)

    assert math.isclose(symbolic.raw_helicity_sum, native.raw_helicity_sum, rel_tol=1e-12)
    assert math.isclose(symbolic.matrix_element, native.matrix_element, rel_tol=1e-12)


def test_zero_gluon_symbolica_evaluator_artifact_roundtrip() -> None:
    particles = (
        ExternalMomentum(1, (45.594, 0.0, 0.0, 45.594)),
        ExternalMomentum(-1, (45.594, 0.0, 0.0, -45.594)),
        ExternalMomentum(23, (91.188, 0.0, 0.0, 0.0)),
    )
    evaluator = ZeroGluonSymbolicEvaluator()
    original = evaluator.evaluate("d d~ > z", particles)

    loaded = ZeroGluonSymbolicEvaluator.from_artifact_payload(
        evaluator.metadata.to_json_dict(),
    )
    roundtrip = loaded.evaluate("d d~ > z", particles)

    assert math.isclose(roundtrip.raw_helicity_sum, original.raw_helicity_sum)
    assert math.isclose(roundtrip.matrix_element, original.matrix_element)
