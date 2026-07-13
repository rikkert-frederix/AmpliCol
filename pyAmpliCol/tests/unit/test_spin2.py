from __future__ import annotations

import math

import pytest

from pyamplicol.generic_artifact import _ext_spin2_filter
from pyamplicol.ufo_model import _expr_spin2_propagator


_METRIC = (1.0, -1.0, -1.0, -1.0)


def _tensor_inner(left: tuple[complex, ...], right: tuple[complex, ...]) -> complex:
    return sum(
        _METRIC[mu]
        * _METRIC[nu]
        * left[4 * mu + nu].conjugate()
        * right[4 * mu + nu]
        for mu in range(4)
        for nu in range(4)
    )


def _assert_physical_spin2_state(
    tensor: tuple[complex, ...],
    momentum: tuple[float, float, float, float],
) -> None:
    for mu in range(4):
        for nu in range(4):
            assert tensor[4 * mu + nu] == pytest.approx(tensor[4 * nu + mu])
    trace = sum(_METRIC[mu] * tensor[4 * mu + mu] for mu in range(4))
    assert trace == pytest.approx(0.0, abs=1.0e-12)
    for nu in range(4):
        contraction = sum(
            _METRIC[mu] * momentum[mu] * tensor[4 * mu + nu]
            for mu in range(4)
        )
        assert contraction == pytest.approx(0.0, abs=1.0e-12)


def test_massless_spin2_helicities_are_physical_and_orthonormal() -> None:
    momentum = (10.0, 3.0, 4.0, math.sqrt(75.0))
    states = {
        helicity: _ext_spin2_filter(momentum, helicity, 0.0)
        for helicity in (-2, 2)
    }

    for state in states.values():
        _assert_physical_spin2_state(state, momentum)
        assert _tensor_inner(state, state) == pytest.approx(1.0)
    assert _tensor_inner(states[-2], states[2]) == pytest.approx(0.0, abs=1.0e-12)


def test_massive_spin2_helicities_are_complete() -> None:
    mass = 4.0
    momentum = (5.0, 3.0, 0.0, 0.0)
    states = {
        helicity: _ext_spin2_filter(momentum, helicity, mass)
        for helicity in (-2, -1, 0, 1, 2)
    }

    for helicity, state in states.items():
        _assert_physical_spin2_state(state, momentum)
        for other_helicity, other_state in states.items():
            expected = 1.0 if helicity == other_helicity else 0.0
            assert _tensor_inner(state, other_state) == pytest.approx(
                expected,
                abs=1.0e-12,
            )

    theta = tuple(
        tuple(
            -(_METRIC[mu] if mu == nu else 0.0)
            + momentum[mu] * momentum[nu] / (mass * mass)
            for nu in range(4)
        )
        for mu in range(4)
    )
    for mu in range(4):
        for nu in range(4):
            for alpha in range(4):
                for beta in range(4):
                    actual = sum(
                        state[4 * mu + nu]
                        * state[4 * alpha + beta].conjugate()
                        for state in states.values()
                    )
                    expected = (
                        0.5
                        * (
                            theta[mu][alpha] * theta[nu][beta]
                            + theta[mu][beta] * theta[nu][alpha]
                        )
                        - theta[mu][nu] * theta[alpha][beta] / 3.0
                    )
                    assert actual == pytest.approx(expected, abs=1.0e-12)


def test_spin2_propagators_apply_de_donder_and_fierz_pauli_projectors() -> None:
    current = tuple(complex(index + 1, 0.5 * index) for index in range(16))
    massless_momentum = (6.0, 1.0, 2.0, 3.0)
    massless = _expr_spin2_propagator(
        current,
        massless_momentum,
        0.0,
        0.0,
        dimension=4.0,
        massive=False,
    )
    denominator = (
        massless_momentum[0] ** 2
        - massless_momentum[1] ** 2
        - massless_momentum[2] ** 2
        - massless_momentum[3] ** 2
    )
    numerator = tuple(component * denominator / 1j for component in massless)
    assert all(
        numerator[4 * mu + nu] == pytest.approx(numerator[4 * nu + mu])
        for mu in range(4)
        for nu in range(4)
    )

    massive_momentum = (5.0, 3.0, 0.0, 0.0)
    massive = _expr_spin2_propagator(
        current,
        massive_momentum,
        4.0,
        1.0,
        dimension=4.0,
        massive=True,
    )
    massive_denominator = 1j * 4.0
    massive_numerator = tuple(
        component * massive_denominator / 1j for component in massive
    )
    _assert_physical_spin2_state(massive_numerator, massive_momentum)
