from __future__ import annotations

import math

import pytest

pytestmark = pytest.mark.skip(
    reason=(
        "legacy native Python kernels are retired from production; schema-v2 "
        "generic DAG/Rusticol tests cover the active path"
    )
)

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


def test_zero_gluon_z_evaluator_accepts_reversed_beam_order() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    forward = evaluator.evaluate("d d~ > z")
    reversed_beams = evaluator.evaluate("d~ d > z")

    assert evaluator.supports_zero_gluon_z("d~ d > z")
    assert [particle.pdg for particle in reversed_beams.particles] == [-1, 1, 23]
    assert math.isclose(
        reversed_beams.matrix_element,
        forward.matrix_element,
        rel_tol=1.0e-14,
    )


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


def test_one_gluon_z_evaluator_accepts_reversed_beam_order() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    forward = evaluator.evaluate("d d~ > z g", sqrt_s=1000.0)
    reversed_beams = evaluator.evaluate("d~ d > z g", sqrt_s=1000.0)

    assert evaluator.supported_z_gluon_count("d~ d > g z") == 1
    assert [particle.pdg for particle in reversed_beams.particles] == [
        -1,
        1,
        21,
        23,
    ]
    assert math.isclose(
        reversed_beams.matrix_element,
        forward.matrix_element,
        rel_tol=1.0e-14,
    )


def test_one_gluon_photon_evaluator_uses_neutral_vector_recursion() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    result = evaluator.evaluate("d d~ > a g", sqrt_s=1000.0)

    assert result.matrix_element > 0.0
    assert result.color_factor == 9
    assert result.average_factor == 36
    assert result.identical_factor == 1
    assert evaluator.supported_photon_gluon_count("d d~ > g a") == 1
    assert [particle.pdg for particle in result.particles] == [1, -1, 21, 22]
    assert len(result.helicity_contributions) == 8


def test_one_gluon_w_evaluator_uses_charged_current_recursion() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    result = evaluator.evaluate("u d~ > w+ g", sqrt_s=1000.0)

    assert result.matrix_element > 0.0
    assert result.color_factor == 9
    assert result.average_factor == 36
    assert result.identical_factor == 1
    assert evaluator.supported_w_gluon_count("u d~ > g w+") == 1
    assert [particle.pdg for particle in result.particles] == [2, -1, 21, 24]
    assert len(result.helicity_contributions) == 12
    assert sum(1 for contribution in result.helicity_contributions if contribution.squared > 0.0) == 6


def test_one_gluon_w_evaluator_accepts_reversed_beam_order() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    forward = evaluator.evaluate("u d~ > w+ g", sqrt_s=1000.0)
    reversed_beams = evaluator.evaluate("d~ u > w+ g", sqrt_s=1000.0)

    assert evaluator.supported_w_gluon_count("d~ u > g w+") == 1
    assert [particle.pdg for particle in reversed_beams.particles] == [
        -1,
        2,
        21,
        24,
    ]
    assert math.isclose(
        reversed_beams.matrix_element,
        forward.matrix_element,
        rel_tol=1.0e-14,
    )


def test_one_gluon_neutral_dilepton_evaluator_matches_amplicol_probe_point() -> None:
    particles = (
        ExternalMomentum(1, (4728.772310962675, 0.0, 0.0, 4728.772310962675)),
        ExternalMomentum(-1, (4945.245406507281, 0.0, 0.0, -4945.245406507281)),
        ExternalMomentum(
            21,
            (
                1734.133260051392,
                98.02973670128142,
                -202.56836593257958,
                1719.469217945512,
            ),
        ),
        ExternalMomentum(
            -11,
            (
                2999.0229133567823,
                108.56720098298311,
                14.681407346592067,
                2997.0211967865766,
            ),
        ),
        ExternalMomentum(
            11,
            (
                4940.861544061782,
                -206.59693768426453,
                187.88695858598751,
                -4932.963510276694,
            ),
        ),
    )

    evaluator = LeadingColorZJetsNativeEvaluator()
    result = evaluator.evaluate("d d~ > e+ e- g", particles=particles)

    assert result.color_factor == 9
    assert result.average_factor == 36
    assert result.identical_factor == 1
    assert math.isclose(result.matrix_element, 5.342225464923099e-08, rel_tol=1.0e-12)
    assert evaluator.supported_neutral_dilepton_gluon_process("d d~ > g e+ e-") == (
        11,
        -11,
        1,
    )
    assert [particle.pdg for particle in result.particles] == [1, -1, 21, 11, -11]
    assert len(result.helicity_contributions) == 16


def test_neutral_dilepton_evaluator_accepts_reversed_beam_order() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    reversed_beams = evaluator.evaluate("d~ d > e+ e- g", sqrt_s=1000.0)

    assert evaluator.supported_neutral_dilepton_gluon_process(
        "d~ d > g e+ e-"
    ) == (11, -11, 1)
    assert [particle.pdg for particle in reversed_beams.particles] == [
        -1,
        1,
        21,
        11,
        -11,
    ]
    assert reversed_beams.matrix_element > 0.0


def test_zero_gluon_neutral_dilepton_evaluator_uses_native_recursion() -> None:
    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "d d~ > e+ e-",
        sqrt_s=1000.0,
    )

    assert result.color_factor == 3
    assert result.average_factor == 36
    assert result.identical_factor == 1
    assert math.isclose(result.matrix_element, 0.0021350534019506165, rel_tol=1.0e-12)
    assert [particle.pdg for particle in result.particles] == [1, -1, 11, -11]
    assert len(result.helicity_contributions) == 8


def test_two_gluon_neutral_dilepton_evaluator_matches_amplicol_probe_point() -> None:
    particles = (
        ExternalMomentum(1, (109.12011508811767, 0.0, 0.0, 109.12011508811767)),
        ExternalMomentum(-1, (56.33553047687237, 0.0, 0.0, -56.33553047687237)),
        ExternalMomentum(
            21,
            (
                62.18121502859448,
                29.790727666707852,
                40.056737551665,
                -37.073896804050456,
            ),
        ),
        ExternalMomentum(
            21,
            (
                100.88787933531779,
                -29.56082674586922,
                -40.40294864455396,
                87.59065851849661,
            ),
        ),
        ExternalMomentum(
            -11,
            (
                1.5998233173131764,
                -4.008451474191155e-06,
                -3.7747315128622283e-06,
                1.5998233173037018,
            ),
        ),
        ExternalMomentum(
            11,
            (
                0.7867278837645935,
                -0.22989691238715626,
                0.34621486762047154,
                0.6679995794954566,
            ),
        ),
    )

    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "d d~ > e+ e- g g",
        particles=particles,
    )

    assert result.color_factor == 27
    assert result.average_factor == 36
    assert result.identical_factor == 2
    assert math.isclose(result.matrix_element, 0.11745153557252273, rel_tol=1.0e-12)
    assert [particle.pdg for particle in result.particles] == [
        1,
        -1,
        21,
        21,
        11,
        -11,
    ]
    assert len(result.helicity_contributions) == 32


def test_one_gluon_charged_leptonic_w_evaluator_matches_amplicol_probe_point() -> None:
    particles = (
        ExternalMomentum(2, (4728.772310962675, 0.0, 0.0, 4728.772310962675)),
        ExternalMomentum(-1, (4945.245406507281, 0.0, 0.0, -4945.245406507281)),
        ExternalMomentum(
            21,
            (
                1734.133260051392,
                98.02973670128142,
                -202.56836593257958,
                1719.469217945512,
            ),
        ),
        ExternalMomentum(
            -11,
            (
                4940.861544061782,
                -206.59693768426453,
                187.88695858598751,
                -4932.963510276694,
            ),
        ),
        ExternalMomentum(
            12,
            (
                2999.0229133567823,
                108.56720098298311,
                14.681407346592067,
                2997.0211967865766,
            ),
        ),
    )

    evaluator = LeadingColorZJetsNativeEvaluator()
    result = evaluator.evaluate("u d~ > e+ ve g", particles=particles)

    assert result.color_factor == 9
    assert result.average_factor == 36
    assert result.identical_factor == 1
    assert math.isclose(result.matrix_element, 2.9438706468720568e-06, rel_tol=1.0e-12)
    assert evaluator.supported_charged_leptonic_w_gluon_process(
        "u d~ > g e+ ve"
    ) == (24, -11, 12, 1)
    assert [particle.pdg for particle in result.particles] == [2, -1, 21, -11, 12]
    assert len(result.helicity_contributions) == 16


def test_charged_leptonic_w_evaluator_accepts_reversed_beam_order() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    reversed_beams = evaluator.evaluate("d~ u > e+ ve g", sqrt_s=1000.0)

    assert evaluator.supported_charged_leptonic_w_gluon_process(
        "d~ u > g e+ ve"
    ) == (24, -11, 12, 1)
    assert [particle.pdg for particle in reversed_beams.particles] == [
        -1,
        2,
        21,
        -11,
        12,
    ]
    assert reversed_beams.matrix_element > 0.0


def test_zero_gluon_charged_leptonic_w_evaluator_uses_native_recursion() -> None:
    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "u d~ > e+ ve",
        sqrt_s=1000.0,
    )

    assert result.color_factor == 3
    assert result.average_factor == 36
    assert result.identical_factor == 1
    assert math.isclose(result.matrix_element, 0.0009046372657017716, rel_tol=1.0e-12)
    assert [particle.pdg for particle in result.particles] == [2, -1, -11, 12]
    assert len(result.helicity_contributions) == 8


def test_one_gluon_charged_leptonic_w_minus_evaluator_uses_native_recursion() -> None:
    evaluator = LeadingColorZJetsNativeEvaluator()
    result = evaluator.evaluate("d u~ > e- ve~ g", sqrt_s=1000.0)

    assert result.matrix_element > 0.0
    assert result.color_factor == 9
    assert result.average_factor == 36
    assert evaluator.supported_charged_leptonic_w_gluon_process(
        "d u~ > g e- ve~"
    ) == (-24, 11, -12, 1)
    assert [particle.pdg for particle in result.particles] == [1, -2, 21, 11, -12]
    assert len(result.helicity_contributions) == 16


def test_two_gluon_charged_leptonic_w_evaluator_uses_ordered_gluon_recursion() -> None:
    result = LeadingColorZJetsNativeEvaluator().evaluate(
        "u d~ > e+ ve g g",
        sqrt_s=1000.0,
    )

    assert result.matrix_element > 0.0
    assert result.color_factor == 27
    assert result.average_factor == 36
    assert result.identical_factor == 2
    assert [particle.pdg for particle in result.particles] == [
        2,
        -1,
        21,
        21,
        -11,
        12,
    ]
    assert len(result.helicity_contributions) == 32


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
