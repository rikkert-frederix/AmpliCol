from __future__ import annotations

import math

from pyamplicol.native import (
    ExternalMomentum,
    LeadingColorZJetsNativeEvaluator,
)


def test_zero_gluon_z_evaluator_matches_ported_amplicol_conventions() -> None:
    result = LeadingColorZJetsNativeEvaluator().evaluate("d d~ > z")

    assert result.color_factor == 3
    assert result.average_factor == 36
    assert math.isclose(result.raw_helicity_sum, 8990.715478289452)
    assert math.isclose(result.matrix_element, 142.10653372872991)

    nonzero = [
        contribution
        for contribution in result.helicity_contributions
        if contribution.squared > 1.0e-12
    ]
    assert [contribution.helicities for contribution in nonzero] == [
        (-1, 1, 1),
        (1, -1, -1),
    ]


def test_zero_gluon_z_evaluator_accepts_explicit_kinematics() -> None:
    mass = 91.188
    particles = (
        ExternalMomentum(1, (mass / 2.0, 0.0, 0.0, mass / 2.0)),
        ExternalMomentum(-1, (mass / 2.0, 0.0, 0.0, -mass / 2.0)),
        ExternalMomentum(23, (mass, 0.0, 0.0, 0.0)),
    )

    result = LeadingColorZJetsNativeEvaluator().evaluate("d d~ > z", particles=particles)

    assert math.isclose(result.matrix_element, 142.10653372872991)


def test_one_gluon_z_evaluator_matches_amplicol_probe_point() -> None:
    particles = (
        ExternalMomentum(1, (1.96470320908030e01, 0.0, 0.0, 1.96470320908030e01)),
        ExternalMomentum(-1, (2.34889854851241e02, 0.0, 0.0, -2.34889854851241e02)),
        ExternalMomentum(
            21,
            (
                1.03330113525346e02,
                8.87977089800681e00,
                2.95057931637355e01,
                -9.86289521374870e01,
            ),
        ),
        ExternalMomentum(
            23,
            (
                1.51206773416698e02,
                -8.87977089800681e00,
                -2.95057931637355e01,
                -1.16613870622951e02,
            ),
        ),
    )

    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "d d~ > z g",
        particles=particles,
    )

    assert result.color_factor == 9
    assert result.average_factor == 36
    assert math.isclose(result.matrix_element, 1.62597906589508, rel_tol=1.0e-7)
    assert len(result.helicity_contributions) == 12


def test_one_gluon_z_evaluator_accepts_reversed_final_order() -> None:
    result = LeadingColorZJetsNativeEvaluator().evaluate("d d~ > g z", sqrt_s=1000.0)

    assert result.matrix_element > 0.0
    assert [particle.pdg for particle in result.particles] == [1, -1, 21, 23]


def test_two_gluon_z_evaluator_matches_amplicol_probe_point() -> None:
    particles = (
        ExternalMomentum(1, (1.2553060477027773e03, 0.0, 0.0, 1.2553060477027773e03)),
        ExternalMomentum(-1, (6.6059804963285214e03, 0.0, 0.0, -6.6059804963285214e03)),
        ExternalMomentum(
            21,
            (
                4.2543742815960177e03,
                -1.0638084437134454e02,
                -6.4917109296463082e01,
                -4.2525485785322971e03,
            ),
        ),
        ExternalMomentum(
            21,
            (
                1.2539358949931880e03,
                -1.0131709337979918e02,
                1.3090951115324427e02,
                1.2429612927324235e03,
            ),
        ),
        ExternalMomentum(
            23,
            (
                2.3529763674420933e03,
                2.0769793775114371e02,
                -6.5992401856781186e01,
                -2.3410871628258706e03,
            ),
        ),
    )

    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "d d~ > z g g",
        particles=particles,
    )

    assert result.color_factor == 27
    assert result.average_factor == 36
    assert result.identical_factor == 2
    assert math.isclose(result.matrix_element, 1.1320063844323158e-02, rel_tol=2.0e-8)
    assert len(result.helicity_contributions) == 24


def test_three_gluon_z_evaluator_matches_amplicol_probe_point() -> None:
    particles = (
        ExternalMomentum(1, (6.058305948333987e01, 0.0, 0.0, 6.058305948333987e01)),
        ExternalMomentum(
            -1,
            (1.9442845962702932e03, 0.0, 0.0, -1.9442845962702932e03),
        ),
        ExternalMomentum(
            21,
            (
                9.957785786477019e02,
                2.2161757176857336e01,
                2.4164979860623113e01,
                -9.95238608556163e02,
            ),
        ),
        ExternalMomentum(
            21,
            (
                6.831042131498102e02,
                5.5223129474791996e01,
                -6.592756778016911e01,
                -6.776690400192582e02,
            ),
        ),
        ExternalMomentum(
            21,
            (
                1.3566631028445107e02,
                5.16141623059848e01,
                -1.9288677346843286e01,
                -1.2397287171739545e02,
            ),
        ),
        ExternalMomentum(
            23,
            (
                1.9031855367167003e02,
                -1.2899904895763413e02,
                6.105126526638929e01,
                -8.682101649413711e01,
            ),
        ),
    )

    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "d d~ > z g g g",
        particles=particles,
    )

    assert result.color_factor == 81
    assert result.average_factor == 36
    assert result.identical_factor == 6
    assert math.isclose(result.matrix_element, 9.932526154443622e-06, rel_tol=6.0e-8)
    assert len(result.helicity_contributions) == 48


def test_evaluator_accepts_z_plus_seven_gluon_family() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()

    process = "d d~ > z g g g g g g g"
    point = evaluator.canonical_z_gluon_point(
        process,
        gluon_count=7,
        sqrt_s=1000.0,
    )

    assert evaluator.supported_z_gluon_count(process) == 7
    assert len(point) == 10
    assert [particle.pdg for particle in point] == [
        1,
        -1,
        21,
        21,
        21,
        21,
        21,
        21,
        21,
        23,
    ]


def test_evaluator_builds_z_plus_eight_gluon_canonical_point_fast_path() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()

    process = "d d~ > z g g g g g g g g"
    point = evaluator.canonical_z_gluon_point(
        process,
        gluon_count=8,
        sqrt_s=1000.0,
    )

    assert evaluator.supported_z_gluon_count(process) == 8
    assert len(point) == 11
    assert [particle.pdg for particle in point] == [
        1,
        -1,
        21,
        21,
        21,
        21,
        21,
        21,
        21,
        21,
        23,
    ]
