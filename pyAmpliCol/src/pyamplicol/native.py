from __future__ import annotations

import math
from collections import Counter
from itertools import product
from typing import Iterable, Sequence

from .core_types import (
    ExternalMomentum,
    FourMomentum,
    HelicityContribution,
    MatrixElementEvaluation,
    NativeEvaluationError,
    TensorWaveFunction,
    WaveFunction,
    WeylWaveFunction,
)
from .model import AmplicolSMLeadingColorModel, Model
from .processes import ANTI_PARTICLE, PDGS, ProcessEnumerator


class LeadingColorZJetsNativeEvaluator:
    """Native evaluator entry point for staged one-quark-line vector + n g kernels."""

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
        vector_target = self.supported_electroweak_vector_gluon_process(process)
        if vector_target is not None:
            vector_pdg, gluon_count = vector_target
            point = (
                tuple(particles)
                if particles is not None
                else self.canonical_neutral_vector_gluon_point(
                    process,
                    vector_pdg=vector_pdg,
                    gluon_count=gluon_count,
                    sqrt_s=sqrt_s,
                )
            )
            return self.evaluate_neutral_vector_gluons(
                process,
                point,
                vector_pdg=vector_pdg,
                gluon_count=gluon_count,
            )
        dilepton_target = self.supported_neutral_dilepton_gluon_process(process)
        if dilepton_target is not None:
            lepton_pdg, antilepton_pdg, gluon_count = dilepton_target
            point = (
                tuple(particles)
                if particles is not None
                else self.canonical_neutral_dilepton_gluon_point(
                    process,
                    lepton_pdg=lepton_pdg,
                    antilepton_pdg=antilepton_pdg,
                    gluon_count=gluon_count,
                    sqrt_s=sqrt_s,
                )
            )
            return self.evaluate_neutral_dilepton_gluons(
                process,
                point,
                lepton_pdg=lepton_pdg,
                antilepton_pdg=antilepton_pdg,
                gluon_count=gluon_count,
            )
        charged_leptonic_target = (
            self.supported_charged_leptonic_w_gluon_process(process)
        )
        if charged_leptonic_target is not None:
            vector_pdg, first_lepton_pdg, second_lepton_pdg, gluon_count = (
                charged_leptonic_target
            )
            point = (
                tuple(particles)
                if particles is not None
                else self.canonical_charged_leptonic_w_gluon_point(
                    process,
                    vector_pdg=vector_pdg,
                    first_lepton_pdg=first_lepton_pdg,
                    second_lepton_pdg=second_lepton_pdg,
                    gluon_count=gluon_count,
                    sqrt_s=sqrt_s,
                )
            )
            return self.evaluate_charged_leptonic_w_gluons(
                process,
                point,
                vector_pdg=vector_pdg,
                first_lepton_pdg=first_lepton_pdg,
                second_lepton_pdg=second_lepton_pdg,
                gluon_count=gluon_count,
            )
        raise NativeEvaluationError(
            "native numerical evaluation is currently implemented only for "
            "one-quark-line electroweak vector, neutral dilepton, or charged-current "
            "leptonic W plus ordered gluons"
        )

    def supports_zero_gluon_z(self, process: str) -> bool:
        try:
            parsed = ProcessEnumerator().parse(process)
        except ValueError:
            return False
        if parsed.rest != ("z",):
            return False
        physical_pdgs = _physical_pdgs(process)
        if len(physical_pdgs) != 3:
            return False
        initial_pdgs = physical_pdgs[:2]
        return (
            all(1 <= abs(pdg) <= 6 for pdg in initial_pdgs)
            and initial_pdgs[0] * initial_pdgs[1] < 0
            and abs(initial_pdgs[0]) == abs(initial_pdgs[1])
        )

    def supports_one_gluon_z(self, process: str) -> bool:
        return self.supported_z_gluon_count(process) == 1

    def supported_z_gluon_count(self, process: str) -> int | None:
        supported = self.supported_electroweak_vector_gluon_process(process)
        if supported is None or supported[0] != 23:
            return None
        return supported[1]

    def supported_photon_gluon_count(self, process: str) -> int | None:
        supported = self.supported_electroweak_vector_gluon_process(process)
        if supported is None or supported[0] != 22:
            return None
        return supported[1]

    def supported_neutral_vector_gluon_process(self, process: str) -> tuple[int, int] | None:
        supported = self.supported_electroweak_vector_gluon_process(process)
        if supported is None or supported[0] not in (22, 23):
            return None
        return supported

    def supported_w_gluon_count(self, process: str) -> int | None:
        supported = self.supported_electroweak_vector_gluon_process(process)
        if supported is None or abs(supported[0]) != 24:
            return None
        return supported[1]

    def supported_electroweak_vector_gluon_process(
        self,
        process: str,
    ) -> tuple[int, int] | None:
        try:
            parsed = ProcessEnumerator().parse(process)
        except ValueError:
            return None
        rest = Counter(parsed.rest)
        vector_names = [
            name for name in ("a", "z", "w+", "w-") if rest.get(name, 0) == 1
        ]
        if len(vector_names) != 1:
            return None
        vector_name = vector_names[0]
        gluon_count = rest.get("g", 0)
        if gluon_count == 0 or sum(rest.values()) != gluon_count + 1:
            return None
        physical_pdgs = _physical_pdgs(process)
        if len(physical_pdgs) != gluon_count + 3:
            return None
        initial_pdgs = physical_pdgs[:2]
        if (
            any(not 1 <= abs(pdg) <= 6 for pdg in initial_pdgs)
            or initial_pdgs[0] * initial_pdgs[1] >= 0
        ):
            return None
        vector_pdg = int(PDGS[vector_name])
        current_pdgs = tuple(-pdg for pdg in initial_pdgs)
        quark_currents = [pdg for pdg in current_pdgs if pdg > 0]
        anti_closures = [pdg for pdg in current_pdgs if pdg < 0]
        if len(quark_currents) != 1 or len(anti_closures) != 1:
            return None
        line_start = quark_currents[0]
        anti_closure = anti_closures[0]
        try:
            line_result = _electroweak_vector_result_pdg(line_start, vector_pdg)
            if anti_closure != -line_result:
                return None
        except NativeEvaluationError:
            return None
        return vector_pdg, gluon_count

    def supported_neutral_dilepton_gluon_process(
        self,
        process: str,
    ) -> tuple[int, int, int] | None:
        try:
            parsed = ProcessEnumerator().parse(process)
        except ValueError:
            return None
        rest = Counter(parsed.rest)
        gluon_count = rest.get("g", 0)
        if sum(rest.values()) != gluon_count + 2:
            return None
        physical_pdgs = _physical_pdgs(process)
        if len(physical_pdgs) != gluon_count + 4:
            return None
        initial_pdgs = physical_pdgs[:2]
        if (
            any(not 1 <= abs(pdg) <= 6 for pdg in initial_pdgs)
            or initial_pdgs[0] * initial_pdgs[1] >= 0
        ):
            return None
        if abs(initial_pdgs[0]) != abs(initial_pdgs[1]):
            return None
        final_leptons = [
            pdg
            for pdg in physical_pdgs[2:]
            if 11 <= abs(pdg) <= 16
        ]
        if len(final_leptons) != 2:
            return None
        positive = [pdg for pdg in final_leptons if pdg > 0]
        negative = [pdg for pdg in final_leptons if pdg < 0]
        if len(positive) != 1 or len(negative) != 1:
            return None
        lepton_pdg = positive[0]
        antilepton_pdg = negative[0]
        if lepton_pdg != -antilepton_pdg:
            return None
        return lepton_pdg, antilepton_pdg, gluon_count

    def supported_charged_leptonic_w_gluon_process(
        self,
        process: str,
    ) -> tuple[int, int, int, int] | None:
        try:
            parsed = ProcessEnumerator().parse(process)
        except ValueError:
            return None
        rest = Counter(parsed.rest)
        gluon_count = rest.get("g", 0)
        if sum(rest.values()) != gluon_count + 2:
            return None
        physical_pdgs = _physical_pdgs(process)
        if len(physical_pdgs) != gluon_count + 4:
            return None
        initial_pdgs = physical_pdgs[:2]
        if (
            any(not 1 <= abs(pdg) <= 6 for pdg in initial_pdgs)
            or initial_pdgs[0] * initial_pdgs[1] >= 0
        ):
            return None
        final_leptons = [
            pdg
            for pdg in physical_pdgs[2:]
            if 11 <= abs(pdg) <= 16
        ]
        if len(final_leptons) != 2:
            return None
        charge3 = sum(_pdg_charge3(pdg) for pdg in final_leptons)
        if charge3 == 3:
            vector_pdg = 24
            charged = [pdg for pdg in final_leptons if pdg in (-11, -13, -15)]
            neutrino = [pdg for pdg in final_leptons if pdg in (12, 14, 16)]
            if len(charged) != 1 or len(neutrino) != 1:
                return None
            first_lepton_pdg, second_lepton_pdg = charged[0], neutrino[0]
        elif charge3 == -3:
            vector_pdg = -24
            charged = [pdg for pdg in final_leptons if pdg in (11, 13, 15)]
            neutrino = [pdg for pdg in final_leptons if pdg in (-12, -14, -16)]
            if len(charged) != 1 or len(neutrino) != 1:
                return None
            first_lepton_pdg, second_lepton_pdg = charged[0], neutrino[0]
        else:
            return None
        if abs(abs(first_lepton_pdg) - abs(second_lepton_pdg)) != 1:
            return None
        current_pdgs = tuple(-pdg for pdg in initial_pdgs)
        quark_currents = [pdg for pdg in current_pdgs if pdg > 0]
        anti_closures = [pdg for pdg in current_pdgs if pdg < 0]
        if len(quark_currents) != 1 or len(anti_closures) != 1:
            return None
        try:
            expected_result = _electroweak_vector_result_pdg(
                quark_currents[0],
                vector_pdg,
            )
        except NativeEvaluationError:
            return None
        if anti_closures[0] != -expected_result:
            return None
        return vector_pdg, first_lepton_pdg, second_lepton_pdg, gluon_count

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
        return self.canonical_neutral_vector_gluon_point(
            process,
            vector_pdg=23,
            gluon_count=gluon_count,
            sqrt_s=sqrt_s,
        )

    def canonical_neutral_vector_gluon_point(
        self,
        process: str,
        *,
        vector_pdg: int,
        gluon_count: int,
        sqrt_s: float | None = None,
    ) -> tuple[ExternalMomentum, ...]:
        pdgs = _physical_pdgs(process)
        energy_cm = self.model.sqrt_s if sqrt_s is None else sqrt_s
        vector_mass = self.model.mass(vector_pdg)
        if energy_cm <= vector_mass:
            raise NativeEvaluationError(
                "q q~ -> electroweak vector + gluons requires sqrt(s) above the vector mass"
            )
        beam_energy = 0.5 * energy_cm
        if gluon_count >= 2:
            gluons, vector = _canonical_multi_gluon_recoil_point(
                energy_cm,
                vector_mass,
                gluon_count,
            )
            return (
                ExternalMomentum(pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
                ExternalMomentum(pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
                *(ExternalMomentum(21, momentum) for momentum in gluons),
                ExternalMomentum(vector_pdg, vector),
            )
        if gluon_count != 1:
            raise NativeEvaluationError(
                "canonical native points are currently available for at least one gluon"
            )
        gluon_energy = (energy_cm**2 - vector_mass**2) / (2.0 * energy_cm)
        vector_energy = (energy_cm**2 + vector_mass**2) / (2.0 * energy_cm)
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
            ExternalMomentum(vector_pdg, (vector_energy, -px, -py, -pz)),
        )

    def canonical_neutral_dilepton_gluon_point(
        self,
        process: str,
        *,
        lepton_pdg: int,
        antilepton_pdg: int,
        gluon_count: int,
        sqrt_s: float | None = None,
    ) -> tuple[ExternalMomentum, ...]:
        pdgs = _physical_pdgs(process)
        energy_cm = self.model.sqrt_s if sqrt_s is None else sqrt_s
        if energy_cm <= 0.0:
            raise NativeEvaluationError("sqrt(s) must be positive")
        beam_energy = 0.5 * energy_cm
        if gluon_count == 0:
            gluons: tuple[FourMomentum, ...] = ()
            dilepton_system = (energy_cm, 0.0, 0.0, 0.0)
        else:
            dilepton_mass = min(max(0.2 * energy_cm, 20.0), 120.0)
            if energy_cm <= dilepton_mass:
                raise NativeEvaluationError("sqrt(s) must be above the dilepton mass")
            if gluon_count >= 2:
                gluons, dilepton_system = _canonical_multi_gluon_recoil_point(
                    energy_cm,
                    dilepton_mass,
                    gluon_count,
                )
            elif gluon_count == 1:
                gluon_energy = (energy_cm**2 - dilepton_mass**2) / (2.0 * energy_cm)
                dilepton_energy = (energy_cm**2 + dilepton_mass**2) / (2.0 * energy_cm)
                direction = _unit_vector(theta=0.7, phi=0.3)
                gluon = (
                    gluon_energy,
                    gluon_energy * direction[0],
                    gluon_energy * direction[1],
                    gluon_energy * direction[2],
                )
                gluons = (gluon,)
                dilepton_system = (
                    dilepton_energy,
                    -gluon[1],
                    -gluon[2],
                    -gluon[3],
                )
            else:
                raise NativeEvaluationError(
                    "neutral dilepton plus gluons needs a non-negative gluon count"
                )
        lepton, antilepton = _decay_timelike_to_two_massless(
            dilepton_system,
            theta=1.03,
            phi=0.41,
        )
        return (
            ExternalMomentum(pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
            ExternalMomentum(pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
            *(ExternalMomentum(21, momentum) for momentum in gluons),
            ExternalMomentum(lepton_pdg, lepton),
            ExternalMomentum(antilepton_pdg, antilepton),
        )

    def canonical_charged_leptonic_w_gluon_point(
        self,
        process: str,
        *,
        vector_pdg: int,
        first_lepton_pdg: int,
        second_lepton_pdg: int,
        gluon_count: int,
        sqrt_s: float | None = None,
    ) -> tuple[ExternalMomentum, ...]:
        pdgs = _physical_pdgs(process)
        energy_cm = self.model.sqrt_s if sqrt_s is None else sqrt_s
        if energy_cm <= 0.0:
            raise NativeEvaluationError("sqrt(s) must be positive")
        beam_energy = 0.5 * energy_cm
        if gluon_count == 0:
            gluons = ()
            w_system = (energy_cm, 0.0, 0.0, 0.0)
        else:
            w_mass = min(max(0.2 * energy_cm, 20.0), 120.0)
            if energy_cm <= w_mass:
                raise NativeEvaluationError("sqrt(s) must be above the leptonic W mass")
            if gluon_count >= 2:
                gluons, w_system = _canonical_multi_gluon_recoil_point(
                    energy_cm,
                    w_mass,
                    gluon_count,
                )
            elif gluon_count == 1:
                gluon_energy = (energy_cm**2 - w_mass**2) / (2.0 * energy_cm)
                w_energy = (energy_cm**2 + w_mass**2) / (2.0 * energy_cm)
                direction = _unit_vector(theta=0.7, phi=0.3)
                gluon = (
                    gluon_energy,
                    gluon_energy * direction[0],
                    gluon_energy * direction[1],
                    gluon_energy * direction[2],
                )
                gluons = (gluon,)
                w_system = (
                    w_energy,
                    -gluon[1],
                    -gluon[2],
                    -gluon[3],
                )
            else:
                raise NativeEvaluationError(
                    "charged-current leptonic W plus gluons needs a non-negative gluon count"
                )
        first, second = _decay_timelike_to_two_massless(
            w_system,
            theta=1.03,
            phi=0.41,
        )
        return (
            ExternalMomentum(pdgs[0], (beam_energy, 0.0, 0.0, beam_energy)),
            ExternalMomentum(pdgs[1], (beam_energy, 0.0, 0.0, -beam_energy)),
            *(ExternalMomentum(21, momentum) for momentum in gluons),
            ExternalMomentum(first_lepton_pdg, first),
            ExternalMomentum(second_lepton_pdg, second),
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
        return self.evaluate_neutral_vector_gluons(
            process,
            particles,
            vector_pdg=23,
            gluon_count=gluon_count,
        )

    def evaluate_neutral_vector_gluons(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
        *,
        vector_pdg: int,
        gluon_count: int,
    ) -> MatrixElementEvaluation:
        point = _normalise_neutral_vector_gluon_particles(
            particles,
            gluon_count=gluon_count,
            vector_pdg=vector_pdg,
        )
        initial_currents = tuple(-particle.pdg for particle in point[:2])
        quark_currents = [pdg for pdg in initial_currents if pdg > 0]
        anti_closures = [pdg for pdg in initial_currents if pdg < 0]
        if len(quark_currents) != 1 or len(anti_closures) != 1:
            raise NativeEvaluationError("incoming pair must contain one quark and one antiquark")
        line_start = quark_currents[0]
        anti_closure = anti_closures[0]
        expected_result = _electroweak_vector_result_pdg(line_start, vector_pdg)
        if anti_closure != -expected_result:
            raise NativeEvaluationError(
                "incoming quark flavours are not compatible with the requested "
                f"vector emission: {line_start} + {vector_pdg} -> {expected_result}, "
                f"but closure requires {-anti_closure}"
            )

        amplitudes = _evaluate_ordered_neutral_vector_gluon_amplicol_recursion(
            self.model,
            point,
            vector_pdg=vector_pdg,
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

    def evaluate_neutral_dilepton_gluons(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
        *,
        lepton_pdg: int,
        antilepton_pdg: int,
        gluon_count: int,
    ) -> MatrixElementEvaluation:
        point = _normalise_neutral_dilepton_gluon_particles(
            particles,
            gluon_count=gluon_count,
            lepton_pdg=lepton_pdg,
            antilepton_pdg=antilepton_pdg,
        )
        incoming = tuple(particle.pdg for particle in point[:2])
        if (
            any(not 1 <= abs(pdg) <= 6 for pdg in incoming)
            or incoming[0] * incoming[1] >= 0
        ):
            raise NativeEvaluationError(
                "incoming pair must contain one quark and one antiquark"
            )
        if abs(incoming[0]) != abs(incoming[1]):
            raise NativeEvaluationError(
                "neutral dilepton native evaluator expects a flavour-preserving "
                "incoming quark/antiquark pair"
            )

        amplitudes = _evaluate_ordered_neutral_dilepton_gluon_amplicol_recursion(
            self.model,
            point,
            lepton_pdg=lepton_pdg,
            antilepton_pdg=antilepton_pdg,
            gluon_count=gluon_count,
        )
        raw_sum = sum(
            float((amplitude * amplitude.conjugate()).real)
            for _, amplitude in amplitudes
        )
        color_factor = self.model.leading_color_factor(
            particle.pdg for particle in point
        )
        average_factor = _initial_state_average_factor(
            particle.pdg for particle in point[:2]
        )
        identical_factor = _final_state_identical_factor(
            particle.pdg for particle in point[2:]
        )
        coupling_factor = (
            (4.0 * math.pi * self.model.alpha_s_me_check) ** gluon_count
            * (2.0 * 4.0 * math.pi * self.model.alpha_ew) ** 2
        )
        matrix_element = (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
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

    def evaluate_charged_leptonic_w_gluons(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
        *,
        vector_pdg: int,
        first_lepton_pdg: int,
        second_lepton_pdg: int,
        gluon_count: int,
    ) -> MatrixElementEvaluation:
        point = _normalise_charged_leptonic_w_gluon_particles(
            particles,
            gluon_count=gluon_count,
            vector_pdg=vector_pdg,
            first_lepton_pdg=first_lepton_pdg,
            second_lepton_pdg=second_lepton_pdg,
        )
        line_start, anti_closure = _incoming_current_line(point)
        expected_result = _electroweak_vector_result_pdg(line_start, vector_pdg)
        if anti_closure != -expected_result:
            raise NativeEvaluationError(
                "incoming quark flavours are not compatible with the requested "
                f"leptonic W current: {line_start} + {vector_pdg} -> "
                f"{expected_result}, but closure requires {-anti_closure}"
            )

        amplitudes = _evaluate_ordered_charged_leptonic_w_gluon_amplicol_recursion(
            self.model,
            point,
            vector_pdg=vector_pdg,
            first_lepton_pdg=first_lepton_pdg,
            second_lepton_pdg=second_lepton_pdg,
            gluon_count=gluon_count,
        )
        raw_sum = sum(
            float((amplitude * amplitude.conjugate()).real)
            for _, amplitude in amplitudes
        )
        color_factor = self.model.leading_color_factor(
            particle.pdg for particle in point
        )
        average_factor = _initial_state_average_factor(
            particle.pdg for particle in point[:2]
        )
        identical_factor = _final_state_identical_factor(
            particle.pdg for particle in point[2:]
        )
        coupling_factor = (
            (4.0 * math.pi * self.model.alpha_s_me_check) ** gluon_count
            * (2.0 * 4.0 * math.pi * self.model.alpha_ew) ** 2
        )
        matrix_element = (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
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


def _incoming_quark_antiquark(
    particles: Sequence[ExternalMomentum],
) -> tuple[ExternalMomentum, ExternalMomentum]:
    incoming = tuple(particles[:2])
    quarks = [particle for particle in incoming if particle.pdg > 0]
    antiquarks = [particle for particle in incoming if particle.pdg < 0]
    if (
        len(quarks) != 1
        or len(antiquarks) != 1
        or not 1 <= abs(quarks[0].pdg) <= 6
        or not 1 <= abs(antiquarks[0].pdg) <= 6
    ):
        raise NativeEvaluationError(
            "incoming pair must contain one quark and one antiquark"
        )
    return quarks[0], antiquarks[0]


def _incoming_current_line(
    particles: Sequence[ExternalMomentum],
) -> tuple[int, int]:
    incoming_currents = tuple(-particle.pdg for particle in particles[:2])
    quark_currents = [pdg for pdg in incoming_currents if pdg > 0]
    anti_closures = [pdg for pdg in incoming_currents if pdg < 0]
    if (
        len(quark_currents) != 1
        or len(anti_closures) != 1
        or not 1 <= quark_currents[0] <= 6
        or not 1 <= abs(anti_closures[0]) <= 6
    ):
        raise NativeEvaluationError(
            "incoming pair must contain one quark and one antiquark"
        )
    return quark_currents[0], anti_closures[0]


def _pdg_charge3(pdg: int) -> int:
    charges = {
        1: -1,
        2: 2,
        3: -1,
        4: 2,
        5: -1,
        6: 2,
        -1: 1,
        -2: -2,
        -3: 1,
        -4: -2,
        -5: 1,
        -6: -2,
        11: -3,
        -11: 3,
        13: -3,
        -13: 3,
        15: -3,
        -15: 3,
        12: 0,
        -12: 0,
        14: 0,
        -14: 0,
        16: 0,
        -16: 0,
        21: 0,
        22: 0,
        23: 0,
        24: 3,
        -24: -3,
        25: 0,
    }
    try:
        return charges[pdg]
    except KeyError as exc:
        raise NativeEvaluationError(f"unsupported PDG charge lookup: {pdg}") from exc


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
    return _normalise_neutral_vector_gluon_particles(
        particles,
        gluon_count=gluon_count,
        vector_pdg=23,
    )


def _normalise_neutral_vector_gluon_particles(
    particles: Sequence[ExternalMomentum],
    *,
    gluon_count: int,
    vector_pdg: int,
) -> tuple[ExternalMomentum, ...]:
    expected = gluon_count + 3
    if len(particles) != expected:
        raise NativeEvaluationError(
            f"q q~ -> electroweak vector + {gluon_count} gluons requires exactly {expected} external momenta"
        )
    if gluon_count < 1:
        raise NativeEvaluationError(
            "native ordered q q~ -> electroweak vector + gluon recursion needs at least one gluon"
        )
    initial = tuple(particles[:2])
    finals = tuple(particles[2:])
    gluons = [particle for particle in finals if particle.pdg == 21]
    vectors = [particle for particle in finals if particle.pdg == vector_pdg]
    if len(gluons) != gluon_count or len(vectors) != 1:
        raise NativeEvaluationError(
            f"final state must contain exactly {gluon_count} gluons and one requested electroweak vector"
        )
    return (*initial, *gluons, vectors[0])


def _normalise_neutral_dilepton_gluon_particles(
    particles: Sequence[ExternalMomentum],
    *,
    gluon_count: int,
    lepton_pdg: int,
    antilepton_pdg: int,
) -> tuple[ExternalMomentum, ...]:
    expected = gluon_count + 4
    if len(particles) != expected:
        raise NativeEvaluationError(
            "q q~ -> neutral dilepton + "
            f"{gluon_count} gluons requires exactly {expected} external momenta"
        )
    if gluon_count < 0:
        raise NativeEvaluationError(
            "native ordered neutral dilepton recursion needs a non-negative gluon count"
        )
    initial = tuple(particles[:2])
    finals = tuple(particles[2:])
    gluons = [particle for particle in finals if particle.pdg == 21]
    leptons = [particle for particle in finals if particle.pdg == lepton_pdg]
    antileptons = [particle for particle in finals if particle.pdg == antilepton_pdg]
    if len(gluons) != gluon_count or len(leptons) != 1 or len(antileptons) != 1:
        raise NativeEvaluationError(
            "final state must contain the requested gluons and exactly one "
            "neutral lepton/antilepton pair"
        )
    return (*initial, *gluons, leptons[0], antileptons[0])


def _normalise_charged_leptonic_w_gluon_particles(
    particles: Sequence[ExternalMomentum],
    *,
    gluon_count: int,
    vector_pdg: int,
    first_lepton_pdg: int,
    second_lepton_pdg: int,
) -> tuple[ExternalMomentum, ...]:
    expected = gluon_count + 4
    if len(particles) != expected:
        raise NativeEvaluationError(
            "q q~ -> charged leptonic W + "
            f"{gluon_count} gluons requires exactly {expected} external momenta"
        )
    if gluon_count < 0:
        raise NativeEvaluationError(
            "native ordered charged leptonic W recursion needs a non-negative gluon count"
        )
    if abs(vector_pdg) != 24:
        raise NativeEvaluationError("charged leptonic W normalization needs W+ or W-")
    initial = tuple(particles[:2])
    finals = tuple(particles[2:])
    gluons = [particle for particle in finals if particle.pdg == 21]
    first = [particle for particle in finals if particle.pdg == first_lepton_pdg]
    second = [particle for particle in finals if particle.pdg == second_lepton_pdg]
    if len(gluons) != gluon_count or len(first) != 1 or len(second) != 1:
        raise NativeEvaluationError(
            "final state must contain the requested gluons and one charged-current "
            "lepton/neutrino pair"
        )
    return (*initial, *gluons, first[0], second[0])


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


def _decay_timelike_to_two_massless(
    parent: FourMomentum,
    *,
    theta: float,
    phi: float,
) -> tuple[FourMomentum, FourMomentum]:
    parent_mass2 = _minkowski_square(parent)
    if parent_mass2 <= 0.0:
        raise NativeEvaluationError("dilepton system must be timelike")
    parent_mass = math.sqrt(parent_mass2)
    energy = 0.5 * parent_mass
    direction = _unit_vector(theta=theta, phi=phi)
    first_rest = (
        energy,
        energy * direction[0],
        energy * direction[1],
        energy * direction[2],
    )
    second_rest = (
        energy,
        -first_rest[1],
        -first_rest[2],
        -first_rest[3],
    )
    beta = (
        parent[1] / parent[0],
        parent[2] / parent[0],
        parent[3] / parent[0],
    )
    return _boost_from_rest(first_rest, beta), _boost_from_rest(second_rest, beta)


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
    return _evaluate_ordered_neutral_vector_gluon_amplicol_recursion(
        model,
        particles,
        vector_pdg=23,
        gluon_count=gluon_count,
    )


def _evaluate_ordered_neutral_vector_gluon_amplicol_recursion(
    model: AmplicolSMLeadingColorModel,
    particles: tuple[ExternalMomentum, ...],
    *,
    vector_pdg: int,
    gluon_count: int,
) -> tuple[tuple[tuple[int, ...], complex], ...]:
    if gluon_count < 1:
        raise NativeEvaluationError(
            "ordered q q~ -> neutral vector + gluon recursion needs at least one gluon"
        )
    incoming_quark, incoming_antiquark = _incoming_quark_antiquark(particles)
    gluons = particles[2 : 2 + gluon_count]
    vector_boson = particles[-1]

    anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
    quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
    gluon_momenta = tuple(gluon.momentum for gluon in gluons)
    vector_momentum = vector_boson.momentum

    anti_plus = _ext_antiquark_weyl(anti_closure_momentum, 1, -1)
    anti_minus = _ext_antiquark_weyl(anti_closure_momentum, -1, 1)
    quark_minus = _ext_quark_weyl(quark_start_momentum, -1, 1)
    quark_plus = _ext_quark_weyl(quark_start_momentum, 1, -1)
    vector_vectors = {
            helicity: _ext_neutral_vector(
                vector_momentum,
                helicity,
                vector_pdg=vector_pdg,
                model=model,
            )
        for helicity in _neutral_vector_helicities(vector_pdg)
    }
    coupling = _neutral_vector_fermion_coupling(
        model,
        vector_pdg=vector_pdg,
        fermion_pdg=-incoming_antiquark.pdg,
    )

    amplitudes: list[tuple[tuple[int, ...], complex]] = []
    for physical_quark_helicity, physical_antiquark_helicity, chirality, quark_wf, anti_wf in (
        (1, -1, 1, quark_minus, anti_plus),
        (-1, 1, -1, quark_plus, anti_minus),
    ):
        for vector_helicity in _neutral_vector_helicities(vector_pdg):
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
                    vector_vectors[vector_helicity],
                    vector_momentum,
                    coupling,
                    chirality,
                )
                amplitudes.append(
                    (
                        (
                            physical_quark_helicity,
                            physical_antiquark_helicity,
                            *gluon_helicities,
                            vector_helicity,
                        ),
                        _weyl_dot(final_current, anti_wf),
                    )
                )
    return tuple(amplitudes)


def _evaluate_ordered_neutral_dilepton_gluon_amplicol_recursion(
    model: AmplicolSMLeadingColorModel,
    particles: tuple[ExternalMomentum, ...],
    *,
    lepton_pdg: int,
    antilepton_pdg: int,
    gluon_count: int,
) -> tuple[tuple[tuple[int, ...], complex], ...]:
    if gluon_count < 0:
        raise NativeEvaluationError(
            "ordered q q~ -> neutral dilepton recursion needs a non-negative gluon count"
        )
    incoming_quark, incoming_antiquark = _incoming_quark_antiquark(particles)
    gluons = particles[2 : 2 + gluon_count]
    lepton = particles[-2]
    antilepton = particles[-1]

    anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
    quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
    gluon_momenta = tuple(gluon.momentum for gluon in gluons)
    dilepton_momentum = _sum_momenta((lepton.momentum, antilepton.momentum))

    anti_plus = _ext_antiquark_weyl(anti_closure_momentum, 1, -1)
    anti_minus = _ext_antiquark_weyl(anti_closure_momentum, -1, 1)
    quark_minus = _ext_quark_weyl(quark_start_momentum, -1, 1)
    quark_plus = _ext_quark_weyl(quark_start_momentum, 1, -1)

    neutral_vectors = (22, 23)
    quark_couplings = {
        vector_pdg: _neutral_vector_fermion_coupling(
            model,
            vector_pdg=vector_pdg,
            fermion_pdg=-incoming_antiquark.pdg,
        )
        for vector_pdg in neutral_vectors
    }
    lepton_couplings = {
        vector_pdg: _neutral_vector_fermion_coupling(
            model,
            vector_pdg=vector_pdg,
            fermion_pdg=lepton_pdg,
        )
        for vector_pdg in neutral_vectors
    }

    amplitudes: list[tuple[tuple[int, ...], complex]] = []
    for physical_quark_helicity, physical_antiquark_helicity, chirality, quark_wf, anti_wf in (
        (1, -1, 1, quark_minus, anti_plus),
        (-1, 1, -1, quark_plus, anti_minus),
    ):
        for lepton_helicity in (-1, 1):
            lepton_chirality = lepton_helicity
            lepton_wf = _ext_quark_weyl(
                lepton.momentum,
                lepton_helicity,
                lepton_chirality,
            )
            for antilepton_helicity in (-1, 1):
                antilepton_chirality = antilepton_helicity
                antilepton_wf = _ext_antiquark_weyl(
                    antilepton.momentum,
                    antilepton_helicity,
                    antilepton_chirality,
                )
                propagated_vectors = {
                    vector_pdg: _neutral_vector_propagator(
                        _lepton_antilepton_to_vector_weyl(
                            lepton_wf,
                            antilepton_wf,
                            lepton_couplings[vector_pdg],
                            lepton_chirality,
                            antilepton_chirality,
                        ),
                        dilepton_momentum,
                        vector_pdg=vector_pdg,
                        model=model,
                    )
                    for vector_pdg in neutral_vectors
                }
                for gluon_helicities in product((-1, 1), repeat=gluon_count):
                    gluon_vectors = tuple(
                        _ext_gluon_cmplx(momentum, helicity)
                        for momentum, helicity in zip(
                            gluon_momenta,
                            gluon_helicities,
                            strict=True,
                        )
                    )
                    gluon_currents = _build_ordered_gluon_currents(
                        gluon_vectors,
                        gluon_momenta,
                    )
                    amplitude = 0j
                    for vector_pdg in neutral_vectors:
                        final_current = _build_ordered_quark_current(
                            quark_wf,
                            quark_start_momentum,
                            gluon_currents,
                            gluon_momenta,
                            propagated_vectors[vector_pdg],
                            dilepton_momentum,
                            quark_couplings[vector_pdg],
                            chirality,
                        )
                        amplitude += _weyl_dot(final_current, anti_wf)
                    amplitudes.append(
                        (
                            (
                                physical_quark_helicity,
                                physical_antiquark_helicity,
                                *gluon_helicities,
                                lepton_helicity,
                                antilepton_helicity,
                            ),
                            amplitude,
                        )
                    )
    return tuple(amplitudes)


def _evaluate_ordered_charged_leptonic_w_gluon_amplicol_recursion(
    model: AmplicolSMLeadingColorModel,
    particles: tuple[ExternalMomentum, ...],
    *,
    vector_pdg: int,
    first_lepton_pdg: int,
    second_lepton_pdg: int,
    gluon_count: int,
) -> tuple[tuple[tuple[int, ...], complex], ...]:
    if gluon_count < 0:
        raise NativeEvaluationError(
            "ordered q q~ -> charged leptonic W recursion needs a non-negative gluon count"
        )
    if abs(vector_pdg) != 24:
        raise NativeEvaluationError("charged leptonic W recursion needs W+ or W-")

    incoming_quark, incoming_antiquark = _incoming_quark_antiquark(particles)
    gluons = particles[2 : 2 + gluon_count]
    first_lepton = particles[-2]
    second_lepton = particles[-1]

    anti_closure_momentum = _negate_momentum(incoming_quark.momentum)
    quark_start_momentum = _negate_momentum(incoming_antiquark.momentum)
    gluon_momenta = tuple(gluon.momentum for gluon in gluons)
    w_momentum = _sum_momenta((first_lepton.momentum, second_lepton.momentum))

    anti_plus = _ext_antiquark_weyl(anti_closure_momentum, 1, -1)
    anti_minus = _ext_antiquark_weyl(anti_closure_momentum, -1, 1)
    quark_minus = _ext_quark_weyl(quark_start_momentum, -1, 1)
    quark_plus = _ext_quark_weyl(quark_start_momentum, 1, -1)
    quark_coupling = _neutral_vector_fermion_coupling(
        model,
        vector_pdg=vector_pdg,
        fermion_pdg=-incoming_antiquark.pdg,
    )
    lepton_coupling = (model.charged_current_coupling(), 0.0)

    amplitudes: list[tuple[tuple[int, ...], complex]] = []
    for physical_quark_helicity, physical_antiquark_helicity, chirality, quark_wf, anti_wf in (
        (1, -1, 1, quark_minus, anti_plus),
        (-1, 1, -1, quark_plus, anti_minus),
    ):
        for first_helicity in (-1, 1):
            first_chirality = first_helicity
            for second_helicity in (-1, 1):
                second_chirality = second_helicity
                if vector_pdg == 24:
                    if first_lepton_pdg not in (-11, -13, -15) or second_lepton_pdg not in (
                        12,
                        14,
                        16,
                    ):
                        raise NativeEvaluationError(
                            "W+ leptonic recursion expects charged antilepton and neutrino"
                        )
                    first_wf = _ext_antiquark_weyl(
                        first_lepton.momentum,
                        first_helicity,
                        first_chirality,
                    )
                    second_wf = _ext_quark_weyl(
                        second_lepton.momentum,
                        second_helicity,
                        second_chirality,
                    )
                    w_current = _antilepton_lepton_to_vector_weyl(
                        first_wf,
                        second_wf,
                        lepton_coupling,
                        first_chirality,
                        second_chirality,
                    )
                else:
                    if first_lepton_pdg not in (11, 13, 15) or second_lepton_pdg not in (
                        -12,
                        -14,
                        -16,
                    ):
                        raise NativeEvaluationError(
                            "W- leptonic recursion expects charged lepton and antineutrino"
                        )
                    first_wf = _ext_quark_weyl(
                        first_lepton.momentum,
                        first_helicity,
                        first_chirality,
                    )
                    second_wf = _ext_antiquark_weyl(
                        second_lepton.momentum,
                        second_helicity,
                        second_chirality,
                    )
                    w_current = _lepton_antilepton_to_vector_weyl(
                        first_wf,
                        second_wf,
                        lepton_coupling,
                        first_chirality,
                        second_chirality,
                    )
                propagated_w = _neutral_vector_propagator(
                    w_current,
                    w_momentum,
                    vector_pdg=vector_pdg,
                    model=model,
                )
                for gluon_helicities in product((-1, 1), repeat=gluon_count):
                    gluon_vectors = tuple(
                        _ext_gluon_cmplx(momentum, helicity)
                        for momentum, helicity in zip(
                            gluon_momenta,
                            gluon_helicities,
                            strict=True,
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
                        propagated_w,
                        w_momentum,
                        quark_coupling,
                        chirality,
                    )
                    amplitudes.append(
                        (
                            (
                                physical_quark_helicity,
                                physical_antiquark_helicity,
                                *gluon_helicities,
                                first_helicity,
                                second_helicity,
                            ),
                            _weyl_dot(final_current, anti_wf),
                        )
                    )
    return tuple(amplitudes)


def _neutral_vector_helicities(vector_pdg: int) -> tuple[int, ...]:
    if vector_pdg == 22:
        return (-1, 1)
    if vector_pdg in (23, 24, -24):
        return (-1, 0, 1)
    raise NativeEvaluationError(f"unsupported electroweak vector PDG: {vector_pdg}")


def _neutral_vector_fermion_coupling(
    model: AmplicolSMLeadingColorModel,
    *,
    vector_pdg: int,
    fermion_pdg: int,
) -> tuple[float, float]:
    if vector_pdg == 22:
        return model.photon_fermion_coupling(fermion_pdg)
    if vector_pdg == 23:
        return model.z_fermion_coupling(fermion_pdg)
    if abs(vector_pdg) == 24:
        _ = _electroweak_vector_result_pdg(fermion_pdg, vector_pdg)
        return (model.charged_current_coupling(), 0.0)
    raise NativeEvaluationError(f"unsupported electroweak vector PDG: {vector_pdg}")


def _ext_neutral_vector(
    momentum: FourMomentum,
    helicity: int,
    *,
    vector_pdg: int,
    model: AmplicolSMLeadingColorModel,
) -> WaveFunction:
    if vector_pdg == 22:
        return _ext_gluon_cmplx(momentum, helicity)
    if vector_pdg in (23, 24, -24):
        return _ext_massive_vector(momentum, helicity, model.mass(vector_pdg))
    raise NativeEvaluationError(f"unsupported electroweak vector PDG: {vector_pdg}")


def _neutral_vector_propagator(
    vector: WaveFunction,
    momentum: FourMomentum,
    *,
    vector_pdg: int,
    model: AmplicolSMLeadingColorModel,
) -> WaveFunction:
    if vector_pdg == 22:
        return _gluon_propagator(vector, momentum)
    if vector_pdg == 23:
        return _massive_vector_propagator(
            vector,
            momentum,
            model.mass(23),
            model.width(23),
        )
    if abs(vector_pdg) == 24:
        return _massive_vector_propagator(
            vector,
            momentum,
            model.mass(vector_pdg),
            model.width(vector_pdg),
        )
    raise NativeEvaluationError(
        f"unsupported electroweak-vector propagator PDG: {vector_pdg}"
    )


def _electroweak_vector_result_pdg(fermion_pdg: int, vector_pdg: int) -> int:
    if not 1 <= fermion_pdg <= 6:
        raise NativeEvaluationError(
            f"charged-current quark-line support expects a quark, got {fermion_pdg}"
        )
    if vector_pdg in (22, 23):
        return fermion_pdg
    if vector_pdg == 24 and fermion_pdg in (1, 3, 5):
        return fermion_pdg + 1
    if vector_pdg == -24 and fermion_pdg in (2, 4, 6):
        return fermion_pdg - 1
    raise NativeEvaluationError(
        f"unsupported charged-current transition {fermion_pdg} + {vector_pdg}"
    )


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
    return (
        (vector[0] - momentum[0] * longitudinal) * prefactor,
        (vector[1] - momentum[1] * longitudinal) * prefactor,
        (vector[2] - momentum[2] * longitudinal) * prefactor,
        (vector[3] - momentum[3] * longitudinal) * prefactor,
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
