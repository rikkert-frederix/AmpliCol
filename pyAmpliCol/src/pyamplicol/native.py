from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import Iterable, Sequence

from .model import AmplicolSMLeadingColorModel, Model
from .processes import ANTI_PARTICLE, PDGS, ProcessEnumerator


FourMomentum = tuple[float, float, float, float]
WaveFunction = tuple[complex, complex, complex, complex]
WeylWaveFunction = tuple[complex, complex]
TensorWaveFunction = tuple[complex, complex, complex, complex, complex, complex]


@dataclass(frozen=True)
class ExternalMomentum:
    pdg: int
    momentum: FourMomentum


@dataclass(frozen=True)
class HelicityContribution:
    helicities: tuple[int, ...]
    amplitude: complex
    squared: float


@dataclass(frozen=True)
class MatrixElementEvaluation:
    process: str
    particles: tuple[ExternalMomentum, ...]
    matrix_element: float
    raw_helicity_sum: float
    color_factor: int
    average_factor: int
    coupling_factor: float
    helicity_contributions: tuple[HelicityContribution, ...]
    identical_factor: int = 1


class NativeEvaluationError(ValueError):
    pass


class LeadingColorZJetsNativeEvaluator:
    """Native evaluator entry point for staged q q~ -> Z + n g kernels."""

    def __init__(self, model: AmplicolSMLeadingColorModel | None = None) -> None:
        self.model: AmplicolSMLeadingColorModel = model or AmplicolSMLeadingColorModel()

    def evaluate(
        self,
        process: str,
        *,
        particles: Sequence[ExternalMomentum] | None = None,
        sqrt_s: float | None = None,
    ) -> MatrixElementEvaluation:
        if self.supports_zero_gluon_z(process):
            point = (
                tuple(particles)
                if particles is not None
                else self.canonical_zero_gluon_point(process, sqrt_s=sqrt_s)
            )
            return self.evaluate_zero_gluon_z(process, point)
        gluon_count = self.supported_z_gluon_count(process)
        if gluon_count is not None and gluon_count >= 1:
            point = (
                tuple(particles)
                if particles is not None
                else self.canonical_z_gluon_point(
                    process,
                    gluon_count=gluon_count,
                    sqrt_s=sqrt_s,
                )
            )
            return self.evaluate_z_gluons(process, point, gluon_count=gluon_count)
        raise NativeEvaluationError(
            "native numerical evaluation is currently implemented only for "
            "q q~ -> Z plus ordered gluons"
        )

    def supports_zero_gluon_z(self, process: str) -> bool:
        try:
            parsed = ProcessEnumerator().parse(process)
        except ValueError:
            return False
        if parsed.rest != ("z",):
            return False
        incoming = parsed.initial_state
        return (
            len(incoming) == 2
            and incoming[0].endswith("~")
            and incoming[0].replace("~", "") == incoming[1].replace("~", "")
        )

    def supports_one_gluon_z(self, process: str) -> bool:
        return self.supported_z_gluon_count(process) == 1

    def supported_z_gluon_count(self, process: str) -> int | None:
        try:
            parsed = ProcessEnumerator().parse(process)
        except ValueError:
            return None
        rest = Counter(parsed.rest)
        if rest.get("z", 0) != 1:
            return None
        gluon_count = rest.get("g", 0)
        if gluon_count == 0 or sum(rest.values()) != gluon_count + 1:
            return None
        incoming = parsed.initial_state
        if not (
            len(incoming) == 2
            and incoming[0].endswith("~")
            and incoming[0].replace("~", "") == incoming[1].replace("~", "")
        ):
            return None
        return gluon_count

    def canonical_zero_gluon_point(
        self,
        process: str,
        *,
        sqrt_s: float | None = None,
    ) -> tuple[ExternalMomentum, ...]:
        pdgs = _physical_pdgs(process)
        energy = (self.model.mass(23) if sqrt_s is None else sqrt_s) / 2.0
        return (
            ExternalMomentum(pdgs[0], (energy, 0.0, 0.0, energy)),
            ExternalMomentum(pdgs[1], (energy, 0.0, 0.0, -energy)),
            ExternalMomentum(23, (2.0 * energy, 0.0, 0.0, 0.0)),
        )

    def canonical_one_gluon_z_point(
        self,
        process: str,
        *,
        sqrt_s: float | None = None,
    ) -> tuple[ExternalMomentum, ...]:
        return self.canonical_z_gluon_point(process, gluon_count=1, sqrt_s=sqrt_s)

    def canonical_z_gluon_point(
        self,
        process: str,
        *,
        gluon_count: int,
        sqrt_s: float | None = None,
    ) -> tuple[ExternalMomentum, ...]:
        pdgs = _physical_pdgs(process)
        energy_cm = self.model.sqrt_s if sqrt_s is None else sqrt_s
        z_mass = self.model.mass(23)
        if energy_cm <= z_mass:
            raise NativeEvaluationError(
                "q q~ -> Z + gluons requires sqrt(s) above the Z mass"
            )
        beam_energy = 0.5 * energy_cm
        if gluon_count >= 2:
            gluons, z = _canonical_multi_gluon_recoil_point(
                energy_cm,
                z_mass,
                gluon_count,
            )
            return (
                ExternalMomentum(pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
                ExternalMomentum(pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
                *(ExternalMomentum(21, momentum) for momentum in gluons),
                ExternalMomentum(23, z),
            )
        if gluon_count != 1:
            raise NativeEvaluationError(
                "canonical native points are currently available for at least one gluon"
            )
        gluon_energy = (energy_cm**2 - z_mass**2) / (2.0 * energy_cm)
        z_energy = (energy_cm**2 + z_mass**2) / (2.0 * energy_cm)
        theta = 0.7
        phi = 0.3
        sin_theta = math.sin(theta)
        px = gluon_energy * sin_theta * math.cos(phi)
        py = gluon_energy * sin_theta * math.sin(phi)
        pz = gluon_energy * math.cos(theta)
        return (
            ExternalMomentum(pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
            ExternalMomentum(pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
            ExternalMomentum(21, (gluon_energy, px, py, pz)),
            ExternalMomentum(23, (z_energy, -px, -py, -pz)),
        )

    def evaluate_zero_gluon_z(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
    ) -> MatrixElementEvaluation:
        if len(particles) != 3:
            raise NativeEvaluationError("q q~ -> Z requires exactly three external momenta")
        if particles[2].pdg != 23:
            raise NativeEvaluationError("third external particle must be a Z boson")
        if particles[0].pdg + particles[1].pdg != 0:
            raise NativeEvaluationError("incoming particles must be a quark/antiquark pair")
        if not 1 <= abs(particles[0].pdg) <= 6:
            raise NativeEvaluationError("incoming pair must be quarks")

        current_inputs = [
            _external_current_input(self.model, particle, is_initial=i < 2)
            for i, particle in enumerate(particles)
        ]
        quark = next(item for item in current_inputs[:2] if item[0] > 0)
        antiquark = next(item for item in current_inputs[:2] if item[0] < 0)
        z_momentum = current_inputs[2][1]
        coupling = self.model.z_fermion_coupling(abs(particles[0].pdg))

        helicity_contributions: list[HelicityContribution] = []
        raw_sum = 0.0
        for h_quark in (-1, 1):
            for h_antiquark in (-1, 1):
                for h_z in (-1, 0, 1):
                    q_wf = _ext_quark(quark[1], h_quark, self.model.mass(quark[0]))
                    aq_wf = _ext_antiquark(
                        antiquark[1],
                        h_antiquark,
                        self.model.mass(antiquark[0]),
                    )
                    z_wf = _ext_massive_vector(z_momentum, h_z, self.model.mass(23))
                    vector_current = _fermion_antifermion_to_vector(q_wf, aq_wf, coupling)
                    amplitude = sum(
                        vector_current[index] * z_wf[index] for index in range(4)
                    )
                    squared = float((amplitude * amplitude.conjugate()).real)
                    raw_sum += squared
                    helicity_contributions.append(
                        HelicityContribution(
                            (h_quark, h_antiquark, h_z),
                            amplitude,
                            squared,
                        )
                    )

        pdgs = tuple(particle.pdg for particle in particles)
        color_factor = self.model.leading_color_factor(pdgs)
        average_factor = _initial_state_average_factor(pdgs[:2])
        coupling_factor = 2.0 * 4.0 * math.pi * self.model.alpha_ew
        matrix_element = raw_sum * color_factor * coupling_factor / average_factor
        return MatrixElementEvaluation(
            process=process,
            particles=tuple(particles),
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=color_factor,
            average_factor=average_factor,
            coupling_factor=coupling_factor,
            helicity_contributions=tuple(helicity_contributions),
        )

    def evaluate_one_gluon_z(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
    ) -> MatrixElementEvaluation:
        return self.evaluate_z_gluons(process, particles, gluon_count=1)

    def evaluate_z_gluons(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
        *,
        gluon_count: int,
    ) -> MatrixElementEvaluation:
        point = _normalise_z_gluon_particles(particles, gluon_count=gluon_count)
        if point[0].pdg + point[1].pdg != 0:
            raise NativeEvaluationError("incoming particles must be a quark/antiquark pair")
        if not 1 <= abs(point[0].pdg) <= 6:
            raise NativeEvaluationError("incoming pair must be quarks")
        if point[0].pdg < 0 or point[1].pdg > 0:
            raise NativeEvaluationError(
                "q q~ -> Z g native evaluator expects the quark before the antiquark"
            )

        amplitudes = _evaluate_ordered_z_gluon_amplicol_recursion(
            self.model,
            point,
            gluon_count=gluon_count,
        )
        raw_sum = sum(
            float((amplitude * amplitude.conjugate()).real)
            for _, amplitude in amplitudes
        )
        color_factor = self.model.leading_color_factor(particle.pdg for particle in point)
        average_factor = _initial_state_average_factor(particle.pdg for particle in point[:2])
        identical_factor = _final_state_identical_factor(particle.pdg for particle in point[2:])
        coupling_factor = (4.0 * math.pi * self.model.alpha_s_me_check) ** gluon_count * (
            2.0 * 4.0 * math.pi * self.model.alpha_ew
        )
        matrix_element = (
            raw_sum * color_factor * coupling_factor / (average_factor * identical_factor)
        )
        return MatrixElementEvaluation(
            process=process,
            particles=point,
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=color_factor,
            average_factor=average_factor,
            coupling_factor=coupling_factor,
            identical_factor=identical_factor,
            helicity_contributions=tuple(
                HelicityContribution(
                    helicities=helicities,
                    amplitude=amplitude,
                    squared=float((amplitude * amplitude.conjugate()).real),
                )
                for helicities, amplitude in amplitudes
            ),
        )


def _physical_pdgs(process: str) -> tuple[int, ...]:
    enumerator = ProcessEnumerator()
    try:
        parsed = enumerator.parse(process)
    except ValueError:
        parsed = None
    if parsed is not None and parsed.jet_count == 0:
        return tuple(
            int(PDGS[ANTI_PARTICLE[p]]) if i < 2 else int(PDGS[p])
            for i, p in enumerate((*parsed.initial_state, *parsed.rest))
        )

    enumeration = enumerator.enumerate(process)
    record = enumeration.groups[0].records[0]
    return tuple(
        int(PDGS[ANTI_PARTICLE[p]]) if i < 2 else int(PDGS[p])
        for i, p in enumerate(record.process)
    )


def _external_current_input(
    model: Model,
    particle: ExternalMomentum,
    *,
    is_initial: bool,
) -> tuple[int, FourMomentum]:
    current_type = model.anti_particle(particle.pdg)
    momentum = _negate_momentum(particle.momentum) if is_initial else particle.momentum
    return current_type, momentum


def _initial_state_average_factor(initial_pdgs: Iterable[int]) -> int:
    factor = 1
    for pdg in initial_pdgs:
        if pdg == 21:
            factor *= 2 * 8
        elif 1 <= abs(pdg) <= 6:
            factor *= 2 * 3
        else:
            factor *= 2
    return factor


def _final_state_identical_factor(final_pdgs: Iterable[int]) -> int:
    factor = 1
    for multiplicity in Counter(final_pdgs).values():
        factor *= math.factorial(multiplicity)
    return factor


def _normalise_one_gluon_z_particles(
    particles: Sequence[ExternalMomentum],
) -> tuple[ExternalMomentum, ExternalMomentum, ExternalMomentum, ExternalMomentum]:
    point = _normalise_z_gluon_particles(particles, gluon_count=1)
    return point[0], point[1], point[2], point[3]


def _normalise_z_gluon_particles(
    particles: Sequence[ExternalMomentum],
    *,
    gluon_count: int,
) -> tuple[ExternalMomentum, ...]:
    expected = gluon_count + 3
    if len(particles) != expected:
        raise NativeEvaluationError(
            f"q q~ -> Z + {gluon_count} gluons requires exactly {expected} external momenta"
        )
    if gluon_count < 1:
        raise NativeEvaluationError(
            "native ordered q q~ -> Z + gluon recursion needs at least one gluon"
        )
    initial = tuple(particles[:2])
    finals = tuple(particles[2:])
    gluons = [particle for particle in finals if particle.pdg == 21]
    z_bosons = [particle for particle in finals if particle.pdg == 23]
    if len(gluons) != gluon_count or len(z_bosons) != 1:
        raise NativeEvaluationError(
            f"final state must contain exactly {gluon_count} gluons and one Z"
        )
    return (*initial, *gluons, z_bosons[0])


def _canonical_multi_gluon_recoil_point(
    sqrt_s: float,
    z_mass: float,
    gluon_count: int,
) -> tuple[tuple[FourMomentum, ...], FourMomentum]:
    if gluon_count < 2:
        raise NativeEvaluationError("multi-gluon recoil point needs at least two gluons")
    fixed_gluons: list[FourMomentum] = []
    soft_scale = 0.035 * sqrt_s / gluon_count
    for i in range(gluon_count - 1):
        energy = soft_scale * (1.0 + 0.23 * i)
        direction = _unit_vector(theta=0.57 + 0.29 * i, phi=0.31 + 1.07 * i)
        fixed_gluons.append(
            (
                energy,
                energy * direction[0],
                energy * direction[1],
                energy * direction[2],
            )
        )

    fixed_total = _sum_momenta(fixed_gluons)
    recoil = (
        sqrt_s - fixed_total[0],
        -fixed_total[1],
        -fixed_total[2],
        -fixed_total[3],
    )
    recoil_mass2 = _minkowski_square(recoil)
    if recoil_mass2 <= z_mass**2:
        raise NativeEvaluationError("canonical multi-gluon recoil is below Z threshold")
    recoil_mass = math.sqrt(recoil_mass2)
    last_energy = (recoil_mass2 - z_mass**2) / (2.0 * recoil_mass)
    z_energy = (recoil_mass2 + z_mass**2) / (2.0 * recoil_mass)
    direction = _unit_vector(theta=1.13, phi=0.73)
    last_rest = (
        last_energy,
        last_energy * direction[0],
        last_energy * direction[1],
        last_energy * direction[2],
    )
    z_rest = (
        z_energy,
        -last_rest[1],
        -last_rest[2],
        -last_rest[3],
    )
    beta = (recoil[1] / recoil[0], recoil[2] / recoil[0], recoil[3] / recoil[0])
    return (
        (*fixed_gluons, _boost_from_rest(last_rest, beta)),
        _boost_from_rest(z_rest, beta),
    )


def _unit_vector(*, theta: float, phi: float) -> tuple[float, float, float]:
    sin_theta = math.sin(theta)
    return sin_theta * math.cos(phi), sin_theta * math.sin(phi), math.cos(theta)


def _boost_from_rest(momentum: FourMomentum, beta: tuple[float, float, float]) -> FourMomentum:
    beta2 = beta[0] ** 2 + beta[1] ** 2 + beta[2] ** 2
    if beta2 == 0.0:
        return momentum
    if beta2 >= 1.0:
        raise NativeEvaluationError("invalid canonical boost with beta >= 1")
    gamma = 1.0 / math.sqrt(1.0 - beta2)
    beta_dot_p = beta[0] * momentum[1] + beta[1] * momentum[2] + beta[2] * momentum[3]
    spatial_factor = ((gamma - 1.0) * beta_dot_p / beta2) + gamma * momentum[0]
    return (
        gamma * (momentum[0] + beta_dot_p),
        momentum[1] + spatial_factor * beta[0],
        momentum[2] + spatial_factor * beta[1],
        momentum[3] + spatial_factor * beta[2],
    )


def _negate_momentum(momentum: FourMomentum) -> FourMomentum:
    return (-momentum[0], -momentum[1], -momentum[2], -momentum[3])


def _ext_quark(momentum: FourMomentum, helicity: int, mass: float) -> WaveFunction:
    tiny = 1.0e-8
    if abs(mass) >= tiny:
        raise NativeEvaluationError("massive quark external wavefunctions are not implemented")
    energy, px, py, pz = momentum
    if energy > 0.0:
        if px == 0.0 and py == 0.0 and pz < 0.0:
            sqp0p3 = 0.0
        else:
            sqp0p3 = math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(helicity * px, -py) / sqp0p3
        )
        if helicity == 1:
            return chi1, chi2, 0j, 0j
        return 0j, 0j, chi2, chi1

    if px == 0.0 and py == 0.0 and pz > 0.0:
        sqp0p3 = 0.0
    else:
        sqp0p3 = -math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(-helicity * (-px), -(-py)) / sqp0p3
    )
    if -helicity == 1:
        return chi1, chi2, 0j, 0j
    return 0j, 0j, chi2, chi1


def _ext_antiquark(momentum: FourMomentum, helicity: int, mass: float) -> WaveFunction:
    tiny = 1.0e-8
    if abs(mass) >= tiny:
        raise NativeEvaluationError(
            "massive antiquark external wavefunctions are not implemented"
        )
    energy, px, py, pz = momentum
    if energy > 0.0:
        if px == 0.0 and py == 0.0 and pz < 0.0:
            sqp0p3 = 0.0
        else:
            sqp0p3 = -math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(-helicity * px, py) / sqp0p3
        )
        if -helicity == 1:
            return 0j, 0j, chi1, chi2
        return chi2, chi1, 0j, 0j

    if px == 0.0 and py == 0.0 and pz > 0.0:
        sqp0p3 = 0.0
    else:
        sqp0p3 = math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(helicity * (-px), -py) / sqp0p3
    )
    if helicity == 1:
        return 0j, 0j, chi1, chi2
    return chi2, chi1, 0j, 0j


def _ext_massive_vector(momentum: FourMomentum, helicity: int, mass: float) -> WaveFunction:
    energy, px, py, pz = momentum
    sqh = math.sqrt(0.5)
    hel = float(helicity)
    nsv = 1
    nsvahl = nsv * abs(helicity)
    pt2 = px**2 + py**2
    pp = min(energy, math.sqrt(pt2 + pz**2))
    pt = min(pp, math.sqrt(pt2))
    if mass == 0.0:
        raise NativeEvaluationError("massless vector wavefunction not expected for Z")
    hel0 = 1.0 - abs(hel)
    if pp == 0.0:
        return (
            0j,
            complex(-hel * sqh),
            complex(0.0, nsvahl * sqh),
            complex(hel0),
        )
    emp = energy / (mass * pp)
    wf0 = complex(hel0 * pp / mass)
    wf3 = complex(hel0 * pz * emp + hel * pt / pp * sqh)
    if pt != 0.0:
        pzpt = pz / (pp * pt) * sqh * hel
        wf1 = complex(hel0 * px * emp - px * pzpt, -nsvahl * py / pt * sqh)
        wf2 = complex(hel0 * py * emp - py * pzpt, nsvahl * px / pt * sqh)
    else:
        wf1 = complex(-hel * sqh)
        wf2 = complex(0.0, nsvahl * _fortran_sign(sqh, pz))
    return wf0, wf1, wf2, wf3


def _ext_gluon_cmplx(momentum: FourMomentum, helicity: int) -> WaveFunction:
    energy, px, py, pz = momentum
    if energy == 0.0:
        raise NativeEvaluationError("cannot generate external gluon with zero energy")
    sqh = math.sqrt(0.5)
    if energy > 0.0:
        hel = float(helicity)
        pp = energy
        pt = math.sqrt(px**2 + py**2)
        wf0 = 0j
        wf3 = complex(hel * pt / pp * sqh)
        if pt != 0.0:
            pzpt = pz / (pp * pt) * sqh * hel
            wf1 = complex(-px * pzpt, -py / pt * sqh)
            wf2 = complex(-py * pzpt, px / pt * sqh)
        else:
            wf1 = complex(-hel * sqh)
            wf2 = complex(0.0, _fortran_sign(sqh, pz))
        return wf0, wf1, wf2, wf3

    hel = float(-helicity)
    pp = -energy
    pt = math.sqrt(px**2 + py**2)
    wf0 = 0j
    wf3 = complex(hel * pt / pp * sqh)
    if pt != 0.0:
        pzpt = -pz / (pp * pt) * sqh * hel
        wf1 = complex(px * pzpt, py / pt * sqh)
        wf2 = complex(py * pzpt, -px / pt * sqh)
    else:
        wf1 = complex(-hel * sqh)
        wf2 = complex(0.0, -_fortran_sign(sqh, pz))
    return wf0, wf1, wf2, wf3


def _ext_quark_weyl(
    momentum: FourMomentum,
    helicity: int,
    chirality: int,
) -> WeylWaveFunction:
    energy, px, py, pz = momentum
    if energy > 0.0:
        if px == 0.0 and py == 0.0 and pz < 0.0:
            sqp0p3 = 0.0
        else:
            sqp0p3 = math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(helicity * px, -py) / sqp0p3
        )
        if helicity == 1 and chirality == 1:
            return chi1, chi2
        if helicity == -1 and chirality == -1:
            return chi2, chi1
        return 0j, 0j

    if px == 0.0 and py == 0.0 and pz > 0.0:
        sqp0p3 = 0.0
    else:
        sqp0p3 = -math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(-helicity * (-px), -(-py)) / sqp0p3
    )
    if helicity == -1 and chirality == 1:
        return chi1, chi2
    if helicity == 1 and chirality == -1:
        return chi2, chi1
    return 0j, 0j


def _ext_antiquark_weyl(
    momentum: FourMomentum,
    helicity: int,
    chirality: int,
) -> WeylWaveFunction:
    energy, px, py, pz = momentum
    if energy > 0.0:
        if px == 0.0 and py == 0.0 and pz < 0.0:
            sqp0p3 = 0.0
        else:
            sqp0p3 = -math.sqrt(max(energy + pz, 0.0))
        chi1 = complex(sqp0p3)
        chi2 = (
            complex(-helicity) * math.sqrt(2.0 * energy)
            if sqp0p3 == 0.0
            else complex(-helicity * px, py) / sqp0p3
        )
        if helicity == 1 and chirality == 1:
            return chi2, chi1
        if helicity == -1 and chirality == -1:
            return chi1, chi2
        return 0j, 0j

    if px == 0.0 and py == 0.0 and pz > 0.0:
        sqp0p3 = 0.0
    else:
        sqp0p3 = math.sqrt(max(-(energy + pz), 0.0))
    chi1 = complex(sqp0p3)
    chi2 = (
        complex(-helicity) * math.sqrt(2.0 * abs(energy))
        if sqp0p3 == 0.0
        else complex(helicity * (-px), -py) / sqp0p3
    )
    if helicity == -1 and chirality == 1:
        return chi2, chi1
    if helicity == 1 and chirality == -1:
        return chi1, chi2
    return 0j, 0j


def _fermion_antifermion_to_vector(
    quark: WaveFunction,
    antiquark: WaveFunction,
    coupling: tuple[float, float],
) -> WaveFunction:
    l1, l2, l3, l4 = quark
    a1, a2, a3, a4 = antiquark
    left, right = coupling
    prefactor = 1j / math.sqrt(2.0)
    return (
        prefactor * (left * (l3 * a1 + l4 * a2) + right * (l1 * a3 + l2 * a4)),
        prefactor * (left * (-l4 * a1 - l3 * a2) + right * (l1 * a4 + l2 * a3)),
        1j
        * prefactor
        * (left * (-l4 * a1 + l3 * a2) + right * (-l1 * a4 + l2 * a3)),
        prefactor * (left * (-l3 * a1 + l4 * a2) + right * (l1 * a3 - l2 * a4)),
    )


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
    raise NativeEvaluationError("QuarkGluontoQuark_weyl called with zero chirality")


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
        return (
            factor * (tmp2 * q1 - tmp3 * q2),
            factor * (tmp1 * q2 - tmp4 * q1),
        )
    if chirality == -1:
        factor = prefactor * left
        return (
            factor * (tmp1 * q1 + tmp3 * q2),
            factor * (tmp2 * q2 + tmp4 * q1),
        )
    raise NativeEvaluationError("QuarkGluontoQuark_coupl_weyl called with zero chirality")


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
        return (
            (tmp1 * q1 + tmp3 * q2) * prefactor,
            (tmp2 * q2 + tmp4 * q1) * prefactor,
        )
    if chirality == -1:
        return (
            (tmp2 * q1 - tmp3 * q2) * prefactor,
            (tmp1 * q2 - tmp4 * q1) * prefactor,
        )
    raise NativeEvaluationError("QuarkPropagator_weyl called with zero chirality")


def _vector_slash_terms(vector: WaveFunction) -> tuple[complex, complex, complex, complex]:
    v0, v1, v2, v3 = vector
    return v0 + v3, v0 - v3, v1 + 1j * v2, v1 - 1j * v2


def _weyl_dot(left: WeylWaveFunction, right: WeylWaveFunction) -> complex:
    return left[0] * right[0] + left[1] * right[1]


def _evaluate_ordered_z_gluon_amplicol_recursion(
    model: AmplicolSMLeadingColorModel,
    particles: tuple[ExternalMomentum, ...],
    *,
    gluon_count: int,
) -> tuple[tuple[tuple[int, ...], complex], ...]:
    if gluon_count < 1:
        raise NativeEvaluationError(
            "ordered q q~ -> Z + gluon recursion needs at least one gluon"
        )
    incoming_quark = particles[0]
    incoming_antiquark = particles[1]
    gluons = particles[2 : 2 + gluon_count]
    z_boson = particles[-1]

    anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
    quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
    gluon_momenta = tuple(gluon.momentum for gluon in gluons)
    z_momentum = z_boson.momentum

    anti_plus = _ext_antiquark_weyl(anti_closure_momentum, 1, -1)
    anti_minus = _ext_antiquark_weyl(anti_closure_momentum, -1, 1)
    quark_minus = _ext_quark_weyl(quark_start_momentum, -1, 1)
    quark_plus = _ext_quark_weyl(quark_start_momentum, 1, -1)
    z_vectors = {
        helicity: _ext_massive_vector(z_momentum, helicity, model.mass(23))
        for helicity in (-1, 0, 1)
    }
    coupling = model.z_fermion_coupling(abs(incoming_quark.pdg))

    amplitudes: list[tuple[tuple[int, ...], complex]] = []
    for physical_quark_helicity, physical_antiquark_helicity, chirality, quark_wf, anti_wf in (
        (1, -1, 1, quark_minus, anti_plus),
        (-1, 1, -1, quark_plus, anti_minus),
    ):
        for z_helicity in (-1, 0, 1):
            for gluon_helicities in product((-1, 1), repeat=gluon_count):
                gluon_vectors = tuple(
                    _ext_gluon_cmplx(momentum, helicity)
                    for momentum, helicity in zip(
                        gluon_momenta, gluon_helicities, strict=True
                    )
                )
                gluon_currents = _build_ordered_gluon_currents(
                    gluon_vectors,
                    gluon_momenta,
                )
                final_current = _build_ordered_quark_current(
                    quark_wf,
                    quark_start_momentum,
                    gluon_currents,
                    gluon_momenta,
                    z_vectors[z_helicity],
                    z_momentum,
                    coupling,
                    chirality,
                )
                amplitudes.append(
                    (
                        (
                            physical_quark_helicity,
                            physical_antiquark_helicity,
                            *gluon_helicities,
                            z_helicity,
                        ),
                        _weyl_dot(final_current, anti_wf),
                    )
                )
    return tuple(amplitudes)


def _build_ordered_gluon_currents(
    external_gluons: tuple[WaveFunction, ...],
    momenta: tuple[FourMomentum, ...],
) -> dict[tuple[int, int], WaveFunction]:
    currents: dict[tuple[int, int], WaveFunction] = {
        (i, i + 1): gluon for i, gluon in enumerate(external_gluons)
    }
    tensor_currents: dict[tuple[int, int], TensorWaveFunction] = {}
    gluon_count = len(external_gluons)
    for length in range(2, gluon_count + 1):
        for start in range(0, gluon_count - length + 1):
            end = start + length
            tensor_contributions = [
                _two_gluon_to_tensor(
                    currents[(start, split)],
                    currents[(split, end)],
                )
                for split in range(start + 1, end)
            ]
            tensor_currents[(start, end)] = _sum_tensor_wavefunctions(
                tensor_contributions
            )

            contributions: list[WaveFunction] = []
            for split in range(start + 1, end):
                left = currents[(start, split)]
                right = currents[(split, end)]
                contributions.append(
                    _three_gluon(
                        left,
                        _sum_momenta(momenta[start:split]),
                        right,
                        _sum_momenta(momenta[split:end]),
                    )
                )
                if split - start >= 2:
                    contributions.append(
                        _tensor_gluon_to_gluon(
                            tensor_currents[(start, split)],
                            right,
                        )
                    )
                if end - split >= 2:
                    contributions.append(
                        _gluon_tensor_to_gluon(
                            left,
                            tensor_currents[(split, end)],
                        )
                    )
            momentum = _sum_momenta(momenta[start:end])
            currents[(start, end)] = _gluon_propagator(
                _sum_wavefunctions(contributions),
                momentum,
            )
    return currents


def _build_ordered_quark_current(
    quark: WeylWaveFunction,
    quark_momentum: FourMomentum,
    gluon_currents: dict[tuple[int, int], WaveFunction],
    gluon_momenta: tuple[FourMomentum, ...],
    z_vector: WaveFunction,
    z_momentum: FourMomentum,
    coupling: tuple[float, float],
    chirality: int,
) -> WeylWaveFunction:
    gluon_count = len(gluon_momenta)
    quark_without_z: dict[int, WeylWaveFunction] = {0: quark}
    for end in range(1, gluon_count + 1):
        contributions = [
            _quark_gluon_to_quark_weyl(
                quark_without_z[split],
                gluon_currents[(split, end)],
                chirality,
            )
            for split in range(0, end)
        ]
        quark_without_z[end] = _quark_propagator_weyl(
            _sum_weyl_wavefunctions(contributions),
            _quark_sequence_momentum(quark_momentum, gluon_momenta, end, None),
            chirality,
        )

    quark_with_z: dict[int, WeylWaveFunction] = {}
    for end in range(0, gluon_count + 1):
        contributions = [
            _quark_gluon_to_quark_coupl_weyl(
                quark_without_z[end],
                z_vector,
                coupling,
                chirality,
            )
        ]
        contributions.extend(
            _quark_gluon_to_quark_weyl(
                quark_with_z[split],
                gluon_currents[(split, end)],
                chirality,
            )
            for split in range(0, end)
        )
        value = _sum_weyl_wavefunctions(contributions)
        if end != gluon_count:
            value = _quark_propagator_weyl(
                value,
                _quark_sequence_momentum(
                    quark_momentum,
                    gluon_momenta,
                    end,
                    z_momentum,
                ),
                chirality,
            )
        quark_with_z[end] = value
    return quark_with_z[gluon_count]


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
    component0 = tmp1 * (left_momentum[0] - right_momentum[0]) + 2.0 * (
        tmp2 * right[0] - tmp3 * left[0]
    )
    component1 = tmp1 * (left_momentum[1] - right_momentum[1]) + 2.0 * (
        tmp2 * right[1] - tmp3 * left[1]
    )
    component2 = tmp1 * (left_momentum[2] - right_momentum[2]) + 2.0 * (
        tmp2 * right[2] - tmp3 * left[2]
    )
    component3 = tmp1 * (left_momentum[3] - right_momentum[3]) + 2.0 * (
        tmp2 * right[3] - tmp3 * left[3]
    )
    return (
        prefactor * component0,
        prefactor * component1,
        prefactor * component2,
        prefactor * component3,
    )


def _gluon_propagator(gluon: WaveFunction, momentum: FourMomentum) -> WaveFunction:
    denominator = _minkowski_square(momentum)
    if denominator == 0.0:
        raise NativeEvaluationError("singular massless gluon propagator")
    prefactor = -1j / denominator
    return (
        gluon[0] * prefactor,
        gluon[1] * prefactor,
        gluon[2] * prefactor,
        gluon[3] * prefactor,
    )


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
        (tensor[0] * gluon[1] + tensor[1] * gluon[2] + tensor[2] * gluon[3])
        * prefactor,
        (tensor[0] * gluon[0] + tensor[3] * gluon[2] + tensor[4] * gluon[3])
        * prefactor,
        (tensor[1] * gluon[0] - tensor[3] * gluon[1] + tensor[5] * gluon[3])
        * prefactor,
        (tensor[2] * gluon[0] - tensor[4] * gluon[1] - tensor[5] * gluon[2])
        * prefactor,
    )


def _gluon_tensor_to_gluon(
    gluon: WaveFunction,
    tensor: TensorWaveFunction,
) -> WaveFunction:
    prefactor = 0.5j
    return (
        (-gluon[1] * tensor[0] - gluon[2] * tensor[1] - gluon[3] * tensor[2])
        * prefactor,
        (-gluon[0] * tensor[0] - gluon[2] * tensor[3] - gluon[3] * tensor[4])
        * prefactor,
        (-gluon[0] * tensor[1] + gluon[1] * tensor[3] - gluon[3] * tensor[5])
        * prefactor,
        (-gluon[0] * tensor[2] + gluon[1] * tensor[4] + gluon[2] * tensor[5])
        * prefactor,
    )


def _quark_sequence_momentum(
    quark_momentum: FourMomentum,
    gluon_momenta: tuple[FourMomentum, ...],
    gluon_count: int,
    z_momentum: FourMomentum | None,
) -> FourMomentum:
    momenta = (quark_momentum, *gluon_momenta[:gluon_count])
    if z_momentum is not None:
        momenta = (*momenta, z_momentum)
    return _sum_momenta(momenta)


def _sum_momenta(momenta: Sequence[FourMomentum]) -> FourMomentum:
    total = (0.0, 0.0, 0.0, 0.0)
    for momentum in momenta:
        total = _add_momenta(total, momentum)
    return total


def _sum_wavefunctions(wavefunctions: Sequence[WaveFunction]) -> WaveFunction:
    return (
        sum(wavefunction[0] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[1] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[2] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[3] for wavefunction in wavefunctions) + 0j,
    )


def _sum_weyl_wavefunctions(wavefunctions: Sequence[WeylWaveFunction]) -> WeylWaveFunction:
    return (
        sum(wavefunction[0] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[1] for wavefunction in wavefunctions) + 0j,
    )


def _sum_tensor_wavefunctions(
    wavefunctions: Sequence[TensorWaveFunction],
) -> TensorWaveFunction:
    return (
        sum(wavefunction[0] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[1] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[2] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[3] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[4] for wavefunction in wavefunctions) + 0j,
        sum(wavefunction[5] for wavefunction in wavefunctions) + 0j,
    )


def _minkowski_dot(left: WaveFunction, right: WaveFunction) -> complex:
    return left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3]


def _minkowski_dot_momentum(vector: WaveFunction, momentum: FourMomentum) -> complex:
    return vector[0] * momentum[0] - vector[1] * momentum[1] - vector[2] * momentum[2] - vector[3] * momentum[3]


def _minkowski_square(momentum: FourMomentum) -> float:
    return (
        momentum[0] ** 2
        - momentum[1] ** 2
        - momentum[2] ** 2
        - momentum[3] ** 2
    )


def _evaluate_one_gluon_z_amplicol_library(
    model: AmplicolSMLeadingColorModel,
    particles: tuple[ExternalMomentum, ExternalMomentum, ExternalMomentum, ExternalMomentum],
) -> tuple[tuple[tuple[int, int, int, int], complex], ...]:
    incoming_quark, incoming_antiquark, gluon, z_boson = particles
    pp1 = _negate_momentum(incoming_quark.momentum)
    pp2 = z_boson.momentum
    pp3 = gluon.momentum
    pp4 = _negate_momentum(incoming_antiquark.momentum)
    pp5 = _add_momenta(pp4, pp2)
    pp6 = _add_momenta(pp4, pp3)

    anti_plus = _ext_antiquark_weyl(pp1, 1, -1)
    anti_minus = _ext_antiquark_weyl(pp1, -1, 1)
    z_vectors = {
        -1: _ext_massive_vector(pp2, -1, model.mass(23)),
        0: _ext_massive_vector(pp2, 0, model.mass(23)),
        1: _ext_massive_vector(pp2, 1, model.mass(23)),
    }
    gluon_vectors = {
        -1: _ext_gluon_cmplx(pp3, -1),
        1: _ext_gluon_cmplx(pp3, 1),
    }
    quark_minus = _ext_quark_weyl(pp4, -1, 1)
    quark_plus = _ext_quark_weyl(pp4, 1, -1)
    coupling = model.z_fermion_coupling(abs(incoming_quark.pdg))

    amplitudes: list[tuple[tuple[int, int, int, int], complex]] = []
    for physical_quark_helicity, physical_antiquark_helicity, chirality, quark_wf, anti_wf in (
        (1, -1, 1, quark_minus, anti_plus),
        (-1, 1, -1, quark_plus, anti_minus),
    ):
        z_then_prop: dict[int, WeylWaveFunction] = {}
        g_then_prop: dict[int, WeylWaveFunction] = {}
        for z_helicity, z_wf in z_vectors.items():
            z_then = _quark_gluon_to_quark_coupl_weyl(
                quark_wf,
                z_wf,
                coupling,
                chirality,
            )
            z_then_prop[z_helicity] = _quark_propagator_weyl(
                z_then,
                pp5,
                chirality,
            )
        for gluon_helicity, gluon_wf in gluon_vectors.items():
            g_then = _quark_gluon_to_quark_weyl(quark_wf, gluon_wf, chirality)
            g_then_prop[gluon_helicity] = _quark_propagator_weyl(
                g_then,
                pp6,
                chirality,
            )

        for z_helicity in (-1, 0, 1):
            for gluon_helicity in (-1, 1):
                term_after_z = _quark_gluon_to_quark_weyl(
                    z_then_prop[z_helicity],
                    gluon_vectors[gluon_helicity],
                    chirality,
                )
                term_after_g = _quark_gluon_to_quark_coupl_weyl(
                    g_then_prop[gluon_helicity],
                    z_vectors[z_helicity],
                    coupling,
                    chirality,
                )
                final_current = (
                    term_after_z[0] + term_after_g[0],
                    term_after_z[1] + term_after_g[1],
                )
                amplitudes.append(
                    (
                        (
                            physical_quark_helicity,
                            physical_antiquark_helicity,
                            gluon_helicity,
                            z_helicity,
                        ),
                        _weyl_dot(final_current, anti_wf),
                    )
                )
    return tuple(amplitudes)


def _add_momenta(left: FourMomentum, right: FourMomentum) -> FourMomentum:
    return (
        left[0] + right[0],
        left[1] + right[1],
        left[2] + right[2],
        left[3] + right[3],
    )


def _fortran_sign(value: float, sign_source: float) -> float:
    return math.copysign(abs(value), sign_source)


__all__ = [
    "ExternalMomentum",
    "HelicityContribution",
    "LeadingColorZJetsNativeEvaluator",
    "MatrixElementEvaluation",
    "NativeEvaluationError",
]
