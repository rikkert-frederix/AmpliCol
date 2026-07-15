from __future__ import annotations

import math

from .core_types import NativeEvaluationError


FourMomentum = tuple[float, float, float, float]
WaveFunction = tuple[complex, complex, complex, complex]
WeylWaveFunction = tuple[complex, complex]
TensorWaveFunction = tuple[complex, complex, complex, complex, complex, complex]


def _lepton_antilepton_to_vector_weyl(
    lepton: WeylWaveFunction,
    antilepton: WeylWaveFunction,
    coupling: tuple[float, float],
    lepton_chirality: int,
    antilepton_chirality: int,
) -> WaveFunction:
    prefactor = 1j / math.sqrt(2.0)
    left, right = coupling
    l1, l2 = lepton
    a1, a2 = antilepton
    if lepton_chirality == -1 and antilepton_chirality == 1:
        factor = prefactor * left
        return (
            factor * (l1 * a1 + l2 * a2),
            -factor * (l2 * a1 + l1 * a2),
            1j * factor * (-l2 * a1 + l1 * a2),
            factor * (-l1 * a1 + l2 * a2),
        )
    if lepton_chirality == 1 and antilepton_chirality == -1:
        factor = prefactor * right
        return (
            factor * (l1 * a1 + l2 * a2),
            factor * (l1 * a2 + l2 * a1),
            1j * factor * (-l1 * a2 + l2 * a1),
            factor * (l1 * a1 - l2 * a2),
        )
    return 0j, 0j, 0j, 0j


def _antilepton_lepton_to_vector_weyl(
    antilepton: WeylWaveFunction,
    lepton: WeylWaveFunction,
    coupling: tuple[float, float],
    antilepton_chirality: int,
    lepton_chirality: int,
) -> WaveFunction:
    prefactor = 1j / math.sqrt(2.0)
    left, right = coupling
    a1, a2 = antilepton
    l1, l2 = lepton
    if antilepton_chirality == 1 and lepton_chirality == -1:
        factor = prefactor * left
        return (
            factor * (l1 * a1 + l2 * a2),
            -factor * (l2 * a1 + l1 * a2),
            1j * factor * (-l2 * a1 + l1 * a2),
            factor * (-l1 * a1 + l2 * a2),
        )
    if antilepton_chirality == -1 and lepton_chirality == 1:
        factor = prefactor * right
        return (
            factor * (l1 * a1 + l2 * a2),
            factor * (l1 * a2 + l2 * a1),
            1j * factor * (-l1 * a2 + l2 * a1),
            factor * (l1 * a1 - l2 * a2),
        )
    return 0j, 0j, 0j, 0j


def _vector_slash_terms(
    vector: WaveFunction,
) -> tuple[complex, complex, complex, complex]:
    v0, v1, v2, v3 = vector
    return v0 + v3, v0 - v3, v1 + 1j * v2, v1 - 1j * v2


def _quark_gluon_to_quark_weyl(
    quark: WeylWaveFunction,
    gluon: WaveFunction,
    chirality: int,
) -> WeylWaveFunction:
    tmp1, tmp2, tmp3, tmp4 = _vector_slash_terms(gluon)
    prefactor = 1j / math.sqrt(2.0)
    q1, q2 = quark
    if chirality == 1:
        return (
            prefactor * (tmp2 * q1 - tmp3 * q2),
            prefactor * (tmp1 * q2 - tmp4 * q1),
        )
    if chirality == -1:
        return (
            prefactor * (tmp1 * q1 + tmp3 * q2),
            prefactor * (tmp2 * q2 + tmp4 * q1),
        )
    raise NativeEvaluationError("quark-vector kernel requires nonzero chirality")


def _quark_gluon_to_quark_coupl_weyl(
    quark: WeylWaveFunction,
    vector: WaveFunction,
    coupling: tuple[float, float],
    chirality: int,
) -> WeylWaveFunction:
    tmp1, tmp2, tmp3, tmp4 = _vector_slash_terms(vector)
    prefactor = 1j / math.sqrt(2.0)
    q1, q2 = quark
    left, right = coupling
    if chirality == 1:
        factor = prefactor * right
        return factor * (tmp2 * q1 - tmp3 * q2), factor * (tmp1 * q2 - tmp4 * q1)
    if chirality == -1:
        factor = prefactor * left
        return factor * (tmp1 * q1 + tmp3 * q2), factor * (tmp2 * q2 + tmp4 * q1)
    raise NativeEvaluationError("coupled quark-vector kernel requires nonzero chirality")


def _quark_propagator_weyl(
    quark: WeylWaveFunction,
    momentum: FourMomentum,
    chirality: int,
) -> WeylWaveFunction:
    energy, px, py, pz = momentum
    denominator = energy**2 - px**2 - py**2 - pz**2
    if denominator == 0.0:
        raise NativeEvaluationError("singular massless quark propagator")
    prefactor = 1j / denominator
    tmp1 = complex(energy + pz)
    tmp2 = complex(energy - pz)
    tmp3 = complex(px, py)
    tmp4 = complex(px, -py)
    q1, q2 = quark
    if chirality == 1:
        return (tmp1 * q1 + tmp3 * q2) * prefactor, (tmp2 * q2 + tmp4 * q1) * prefactor
    if chirality == -1:
        return (tmp2 * q1 - tmp3 * q2) * prefactor, (tmp1 * q2 - tmp4 * q1) * prefactor
    raise NativeEvaluationError("quark propagator requires nonzero chirality")


def _three_gluon(
    left: WaveFunction,
    left_momentum: FourMomentum,
    right: WaveFunction,
    right_momentum: FourMomentum,
) -> WaveFunction:
    tmp1 = _minkowski_dot(left, right)
    tmp2 = _minkowski_dot_momentum(left, right_momentum)
    tmp3 = _minkowski_dot_momentum(right, left_momentum)
    prefactor = 1j / math.sqrt(2.0)
    return tuple(
        prefactor
        * (
            tmp1 * (left_momentum[index] - right_momentum[index])
            + 2.0 * (tmp2 * right[index] - tmp3 * left[index])
        )
        for index in range(4)
    )  # type: ignore[return-value]


def _gluon_propagator(
    gluon: WaveFunction,
    momentum: FourMomentum,
) -> WaveFunction:
    denominator = _minkowski_square(momentum)
    if denominator == 0.0:
        raise NativeEvaluationError("singular massless gluon propagator")
    prefactor = -1j / denominator
    return tuple(value * prefactor for value in gluon)  # type: ignore[return-value]


def _massive_vector_propagator(
    vector: WaveFunction,
    momentum: FourMomentum,
    mass: float,
    width: float,
) -> WaveFunction:
    denominator = _minkowski_square(momentum) - mass**2 + 1j * mass * width
    if denominator == 0.0:
        raise NativeEvaluationError("singular massive vector propagator")
    prefactor = -1j / denominator
    longitudinal = _minkowski_dot_momentum(vector, momentum) / mass**2
    return tuple(
        (vector[index] - momentum[index] * longitudinal) * prefactor
        for index in range(4)
    )  # type: ignore[return-value]


def _two_gluon_to_tensor(
    left: WaveFunction,
    right: WaveFunction,
) -> TensorWaveFunction:
    return (
        left[0] * right[1] - left[1] * right[0],
        left[0] * right[2] - left[2] * right[0],
        left[0] * right[3] - left[3] * right[0],
        left[1] * right[2] - left[2] * right[1],
        left[1] * right[3] - left[3] * right[1],
        left[2] * right[3] - left[3] * right[2],
    )


def _tensor_gluon_to_gluon(
    tensor: TensorWaveFunction,
    gluon: WaveFunction,
) -> WaveFunction:
    prefactor = 0.5j
    return (
        (tensor[0] * gluon[1] + tensor[1] * gluon[2] + tensor[2] * gluon[3]) * prefactor,
        (tensor[0] * gluon[0] + tensor[3] * gluon[2] + tensor[4] * gluon[3]) * prefactor,
        (tensor[1] * gluon[0] - tensor[3] * gluon[1] + tensor[5] * gluon[3]) * prefactor,
        (tensor[2] * gluon[0] - tensor[4] * gluon[1] - tensor[5] * gluon[2]) * prefactor,
    )


def _gluon_tensor_to_gluon(
    gluon: WaveFunction,
    tensor: TensorWaveFunction,
) -> WaveFunction:
    prefactor = 0.5j
    return (
        (-gluon[1] * tensor[0] - gluon[2] * tensor[1] - gluon[3] * tensor[2]) * prefactor,
        (-gluon[0] * tensor[0] - gluon[2] * tensor[3] - gluon[3] * tensor[4]) * prefactor,
        (-gluon[0] * tensor[1] + gluon[1] * tensor[3] - gluon[3] * tensor[5]) * prefactor,
        (-gluon[0] * tensor[2] + gluon[1] * tensor[4] + gluon[2] * tensor[5]) * prefactor,
    )


def _minkowski_dot(left: WaveFunction, right: WaveFunction) -> complex:
    return left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3]


def _minkowski_dot_momentum(
    vector: WaveFunction,
    momentum: FourMomentum,
) -> complex:
    return (
        vector[0] * momentum[0]
        - vector[1] * momentum[1]
        - vector[2] * momentum[2]
        - vector[3] * momentum[3]
    )


def _minkowski_square(momentum: FourMomentum) -> float:
    return momentum[0] ** 2 - momentum[1] ** 2 - momentum[2] ** 2 - momentum[3] ** 2
