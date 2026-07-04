from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import cached_property
from typing import Any, Iterable, Sequence

from .symbols import symbols


@dataclass(frozen=True)
class Particle:
    pdg: int
    anti_pdg: int
    spin: int
    dimension: int
    color_rep: int
    mass: float = 0.0
    width: float = 0.0
    charge: float = 0.0
    weak_isospin: tuple[float, float] = (0.0, 0.0)
    weak_hypercharge: tuple[float, float] = (0.0, 0.0)


@dataclass(frozen=True)
class Vertex:
    kind: int
    particles: tuple[int, int, int]
    coupling: tuple[float, float] = (1.0, 0.0)


@dataclass(frozen=True)
class VertexLoweringRule:
    kind: int
    backend: str
    tensor_names: tuple[str, ...] = ()
    expression_head: str = ""
    full_tensor_network_ready: bool = False
    description: str = ""


@dataclass
class Model:
    name: str
    particles: dict[int, Particle] = field(default_factory=dict)
    vertices: tuple[Vertex, ...] = ()

    def particle(self, pdg: int) -> Particle:
        species = self._species_particle(pdg)
        if species is None:
            raise KeyError(f"particle not in model: {pdg}")
        return species

    def anti_particle(self, pdg: int) -> int:
        particle = self.particle(pdg)
        return particle.anti_pdg if particle.pdg == pdg else particle.pdg

    def mass(self, pdg: int) -> float:
        return self.particle(pdg).mass

    def width(self, pdg: int) -> float:
        return self.particle(pdg).width

    def spin(self, pdg: int) -> int:
        spin = self.particle(pdg).spin
        if spin < 0:
            raise ValueError(f"spin is ill-defined for particle {pdg}")
        return spin

    def dimension(self, pdg: int) -> int:
        return self.particle(pdg).dimension

    def charge(self, pdg: int) -> float:
        return self._property_sign(pdg) * self.particle(pdg).charge

    def weak_isospin_l(self, pdg: int) -> float:
        return self._property_sign(pdg) * self.particle(pdg).weak_isospin[0]

    def weak_isospin_r(self, pdg: int) -> float:
        return self._property_sign(pdg) * self.particle(pdg).weak_isospin[1]

    def color_rep(self, pdg: int) -> int:
        color = self.particle(pdg).color_rep
        if self._property_sign(pdg) < 0 and abs(color) == 3:
            return -color
        return color

    def color_dim(self, pdg: int) -> int:
        return abs(self.color_rep(pdg))

    def is_quark(self, pdg: int) -> bool:
        return 1 <= pdg <= 6

    def is_antiquark(self, pdg: int) -> bool:
        return -6 <= pdg <= -1

    def is_lepton(self, pdg: int) -> bool:
        return 11 <= pdg <= 16

    def is_antilepton(self, pdg: int) -> bool:
        return -16 <= pdg <= -11

    def is_fermion(self, pdg: int) -> bool:
        return (
            self.is_quark(pdg)
            or self.is_antiquark(pdg)
            or self.is_lepton(pdg)
            or self.is_antilepton(pdg)
        )

    def is_chiral_eligible(self, pdg: int) -> bool:
        return self.is_fermion(pdg) and self.mass(pdg) == 0.0

    def is_gluon(self, pdg: int) -> bool:
        return pdg in (21, 99)

    def is_singlet(self, pdg: int) -> bool:
        return not (abs(pdg) <= 6 or pdg == 21)

    def is_tensor(self, pdg: int) -> bool:
        return pdg in (-21, -23, 26, -26)

    def is_massive_boson(self, pdg: int) -> bool:
        return pdg == 23 or abs(pdg) == 24

    def is_photon(self, pdg: int) -> bool:
        return pdg == 22

    def is_higgs(self, pdg: int) -> bool:
        return pdg == 25

    def build_tensor_library(self) -> Any:
        raise NotImplementedError

    def vertex_lowering_rule(self, kind: int) -> VertexLoweringRule:
        raise NotImplementedError

    def three_gluon_current_expression(
        self,
        *,
        left_slot: Any,
        right_slot: Any,
        output_slot: Any,
        left_momentum_tensor_name: str,
        right_momentum_tensor_name: str,
        dummy_prefix: str,
    ) -> Any:
        raise NotImplementedError

    def gluon_propagator_tensor_data(
        self,
        momentum: Sequence[Any],
    ) -> list[Any]:
        raise NotImplementedError

    def quark_weyl_propagator_tensor_data(
        self,
        momentum: Sequence[Any],
        *,
        chirality: int,
    ) -> list[Any]:
        raise NotImplementedError

    def _species_particle(self, pdg: int) -> Particle | None:
        particle = self.particles.get(pdg)
        if particle is not None:
            return particle
        for candidate in self.particles.values():
            if candidate.anti_pdg == pdg:
                return candidate
        return None

    def _property_sign(self, pdg: int) -> int:
        particle = self.particles.get(pdg)
        if particle is not None and particle.pdg == pdg:
            return 1
        for candidate in self.particles.values():
            if candidate.anti_pdg == pdg:
                return -1
        raise KeyError(f"particle not in model: {pdg}")


@dataclass
class AmplicolSMLeadingColorModel(Model):
    """Model data and tensor registration matching legacy AmpliCol conventions."""

    name: str = "amplicol-sm-leading-color"
    alpha_s_mz: float = 0.119
    alpha_s_me_check: float = 0.118
    alpha_ew: float = 0.007546771114
    sin_weak: float = 0.47143025548407230
    sqrt_s: float = 14000.0

    def __post_init__(self) -> None:
        self.particles = {particle.pdg: particle for particle in self._build_particles()}
        self.vertices = tuple(self._build_vertices())

    @cached_property
    def cos_weak(self) -> float:
        return math.sqrt(1.0 - self.sin_weak**2)

    def weak_coupling(self) -> float:
        return 1.0 / self.sin_weak

    def neutral_gauge_coupling(self) -> float:
        return self.weak_coupling() * self.cos_weak

    def charged_current_coupling(self) -> float:
        return self.weak_coupling() / math.sqrt(2.0)

    def weak_coupling_over_cosine(self) -> float:
        return self.weak_coupling() / self.cos_weak

    def photon_fermion_coupling(self, pdg: int) -> tuple[float, float]:
        particle = self.particle(pdg)
        return particle.charge, particle.charge

    def z_fermion_coupling(self, pdg: int) -> tuple[float, float]:
        particle = self.particle(pdg)
        charge = particle.charge
        left = particle.weak_isospin[0]
        right = particle.weak_isospin[1]
        prefactor = self.weak_coupling_over_cosine()
        return (
            prefactor * (left - charge * self.sin_weak**2),
            prefactor * (right - charge * self.sin_weak**2),
        )

    def leading_color_factor(self, process: Iterable[int]) -> int:
        exponent_twice = 0
        for pdg in process:
            if pdg == 21:
                exponent_twice += 2
            elif 1 <= abs(pdg) <= 6:
                exponent_twice += 1
        if exponent_twice % 2:
            raise ValueError(f"non-integer leading-color exponent for {tuple(process)}")
        return 3 ** (exponent_twice // 2)

    def build_tensor_library(self) -> Any:
        from symbolica.community.spenso import (
            LibraryTensor,
            Representation,
            TensorLibrary,
            TensorName,
        )

        library = TensorLibrary.hep_lib_atom()
        mink = Representation.mink(4)
        antisym = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
        weyl = Representation("pyamplicol::weyl_spinor", 2)
        two_gluon_to_tensor = TensorName(str(symbols.two_gluon_to_tensor))
        tensor_gluon_to_gluon = TensorName(str(symbols.tensor_gluon_to_gluon))
        gluon_tensor_to_gluon = TensorName(str(symbols.gluon_tensor_to_gluon))
        quark_vector_weyl_plus = TensorName(str(symbols.quark_vector_weyl_plus))
        quark_vector_weyl_minus = TensorName(str(symbols.quark_vector_weyl_minus))

        library.register(
            LibraryTensor.dense(
                two_gluon_to_tensor(mink, mink, antisym),
                _two_gluon_to_tensor_data(),
            )
        )
        library.register(
            LibraryTensor.dense(
                tensor_gluon_to_gluon(antisym, mink, mink),
                _tensor_gluon_to_gluon_data(),
            )
        )
        library.register(
            LibraryTensor.dense(
                gluon_tensor_to_gluon(mink, antisym, mink),
                _gluon_tensor_to_gluon_data(),
            )
        )
        library.register(
            LibraryTensor.dense(
                quark_vector_weyl_plus(weyl, mink, weyl),
                _quark_vector_weyl_data(chirality=1),
            )
        )
        library.register(
            LibraryTensor.dense(
                quark_vector_weyl_minus(weyl, mink, weyl),
                _quark_vector_weyl_data(chirality=-1),
            )
        )
        return library

    def vertex_lowering_rule(self, kind: int) -> VertexLoweringRule:
        if kind == 0:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=("g", str(symbols.current_momentum)),
                expression_head="three_gluon_current",
                full_tensor_network_ready=True,
                description="color-ordered three-gluon current",
            )
        if kind == 1:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(str(symbols.two_gluon_to_tensor),),
                expression_head=str(symbols.two_gluon_to_tensor),
                full_tensor_network_ready=True,
                description="two gluons to auxiliary antisymmetric tensor",
            )
        if kind == 2:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(str(symbols.tensor_gluon_to_gluon),),
                expression_head=str(symbols.tensor_gluon_to_gluon),
                full_tensor_network_ready=True,
                description="auxiliary tensor and gluon to gluon current",
            )
        if kind == 3:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(str(symbols.gluon_tensor_to_gluon),),
                expression_head=str(symbols.gluon_tensor_to_gluon),
                full_tensor_network_ready=True,
                description="gluon and auxiliary tensor to gluon current",
            )
        if kind == 6:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(
                    str(symbols.quark_vector_weyl_minus),
                    str(symbols.quark_vector_weyl_plus),
                ),
                expression_head="quark_gluon_weyl_current",
                full_tensor_network_ready=True,
                description="Weyl quark-gluon current",
            )
        if kind in {4, 5, 7, 9}:
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica-pending",
                expression_head="quark_gluon_weyl_current",
                description="Weyl quark-gluon current",
            )
        if kind == 10:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(
                    str(symbols.quark_vector_weyl_minus),
                    str(symbols.quark_vector_weyl_plus),
                ),
                expression_head="fermion_gauge_weyl_current",
                full_tensor_network_ready=True,
                description="Weyl fermion electroweak current with graph coupling",
            )
        if kind in {11, 21, 22, 23, 24}:
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica-pending",
                expression_head="fermion_gauge_weyl_current",
                description="Weyl fermion electroweak current with runtime coupling",
            )
        return VertexLoweringRule(
            kind=kind,
            backend="unimplemented",
            description="no native pyamplicol lowering rule is registered yet",
        )

    def three_gluon_current_expression(
        self,
        *,
        left_slot: Any,
        right_slot: Any,
        output_slot: Any,
        left_momentum_tensor_name: str,
        right_momentum_tensor_name: str,
        dummy_prefix: str,
    ) -> Any:
        from symbolica import Expression
        from symbolica.community.spenso import Representation, TensorName

        mink = Representation.mink(4)
        metric = TensorName.g()
        left_momentum = TensorName(left_momentum_tensor_name)
        right_momentum = TensorName(right_momentum_tensor_name)
        left_dot_slot = mink(f"{dummy_prefix}_left_dot")
        right_dot_slot = mink(f"{dummy_prefix}_right_dot")
        prefactor = Expression.num(1j / math.sqrt(2.0))

        return prefactor * (
            metric(left_slot, right_slot).to_expression()
            * (
                left_momentum(output_slot).to_expression()
                - right_momentum(output_slot).to_expression()
            )
            + Expression.num(2.0)
            * (
                metric(left_slot, left_dot_slot).to_expression()
                * right_momentum(left_dot_slot).to_expression()
                * metric(right_slot, output_slot).to_expression()
                - metric(right_slot, right_dot_slot).to_expression()
                * left_momentum(right_dot_slot).to_expression()
                * metric(left_slot, output_slot).to_expression()
            )
        )

    def gluon_propagator_tensor_data(
        self,
        momentum: Sequence[Any],
    ) -> list[Any]:
        prefactor = _number(-1j) / _minkowski_square_expression(momentum)
        metric = (1.0, -1.0, -1.0, -1.0)
        data = [_number(0.0)] * (4 * 4)
        for index, sign in enumerate(metric):
            data[_flat_index((index, index), (4, 4))] = _number(sign) * prefactor
        return data

    def quark_weyl_propagator_tensor_data(
        self,
        momentum: Sequence[Any],
        *,
        chirality: int,
    ) -> list[Any]:
        prefactor = _number(1j) / _minkowski_square_expression(momentum)
        p0, p1, p2, p3 = (_as_expression(value) for value in momentum)
        data = [_number(0.0)] * (2 * 2)

        def set_entry(q_in: int, q_out: int, value: Any) -> None:
            data[_flat_index((q_in, q_out), (2, 2))] = prefactor * value

        if chirality == 1:
            set_entry(0, 0, p0 + p3)
            set_entry(1, 0, p1 + _number(1j) * p2)
            set_entry(0, 1, p1 - _number(1j) * p2)
            set_entry(1, 1, p0 - p3)
            return data
        if chirality == -1:
            set_entry(0, 0, p0 - p3)
            set_entry(1, 0, -(p1 + _number(1j) * p2))
            set_entry(0, 1, -(p1 - _number(1j) * p2))
            set_entry(1, 1, p0 + p3)
            return data
        raise ValueError(f"unsupported Weyl chirality: {chirality}")

    def _build_particles(self) -> list[Particle]:
        particles: list[Particle] = []
        for i in range(1, 6):
            if i % 2 == 0:
                particles.append(
                    Particle(i, -i, 2, 4, 3, charge=2.0 / 3.0, weak_isospin=(0.5, 0.0), weak_hypercharge=(1.0 / 3.0, 4.0 / 3.0))
                )
            else:
                particles.append(
                    Particle(i, -i, 2, 4, 3, charge=-1.0 / 3.0, weak_isospin=(-0.5, 0.0), weak_hypercharge=(1.0 / 3.0, -2.0 / 3.0))
                )
        particles.extend(
            [
                Particle(6, -6, 2, 4, 3, mass=173.0, width=1.491500, charge=2.0 / 3.0, weak_isospin=(0.5, 0.0), weak_hypercharge=(1.0 / 3.0, 4.0 / 3.0)),
                Particle(21, 21, 2, 4, 8),
                Particle(99, 99, -1, 4, 1),
                Particle(-21, -21, -1, 6, 8),
                Particle(22, 22, 2, 4, 1),
                Particle(23, 23, 3, 4, 1, mass=91.188, width=2.441404),
                Particle(-23, -23, -1, 6, 1),
                Particle(24, -24, 3, 4, 1, mass=80.419002445756163, width=2.0476, charge=1.0, weak_isospin=(1.0, 1.0)),
                Particle(25, 25, 1, 1, 1, mass=125.0, width=0.0063823389999999999, weak_isospin=(-0.5, -0.5), weak_hypercharge=(1.0, 1.0)),
                Particle(125, 125, -1, 1, 1),
                Particle(126, 126, -1, 1, 1),
                Particle(127, 127, -1, 1, 1),
                Particle(26, -26, -1, 6, 1, charge=1.0, weak_isospin=(1.0, 1.0)),
            ]
        )
        for i in range(1, 4):
            charged = 11 + (2 * i - 2)
            neutrino = 12 + (2 * i - 2)
            particles.append(
                Particle(charged, -charged, 2, 4, 1, charge=-1.0, weak_isospin=(-0.5, 0.0), weak_hypercharge=(-1.0, -2.0))
            )
            particles.append(
                Particle(neutrino, -neutrino, 2, 4, 1, weak_isospin=(0.5, 0.0), weak_hypercharge=(-1.0, 0.0))
            )
        return particles

    def _build_vertices(self) -> list[Vertex]:
        vertices: list[Vertex] = [
            Vertex(0, (21, 21, 21)),
            Vertex(1, (21, 21, -21)),
            Vertex(2, (-21, 21, 21)),
            Vertex(3, (21, -21, 21)),
        ]
        self._extend_quark_gluon_vertices(vertices)
        self._extend_electroweak_gauge_vertices(vertices)
        self._extend_higgs_vertices(vertices)
        self._extend_lepton_vertices(vertices)
        return vertices

    def _extend_quark_gluon_vertices(self, vertices: list[Vertex]) -> None:
        for i in range(1, 7):
            vertices.extend(
                [
                    Vertex(4, (21, i, i)),
                    Vertex(5, (21, -i, -i)),
                    Vertex(6, (i, 21, i)),
                    Vertex(7, (-i, 21, -i)),
                    Vertex(9, (-i, i, 21)),
                    Vertex(8, (i, -i, 99), (1.0 / 3.0, 0.0)),
                    Vertex(4, (99, i, i)),
                    Vertex(5, (99, -i, -i)),
                    Vertex(6, (i, 99, i)),
                    Vertex(7, (-i, 99, -i)),
                    Vertex(10, (i, 22, i), self.photon_fermion_coupling(i)),
                    Vertex(11, (-i, 22, -i), self.photon_fermion_coupling(-i)),
                    Vertex(10, (i, 23, i), self.z_fermion_coupling(i)),
                    Vertex(11, (-i, 23, -i), self.z_fermion_coupling(-i)),
                ]
            )
            if i % 2 == 0:
                vertices.append(Vertex(10, (i, -24, i - 1), (self.charged_current_coupling(), 0.0)))
                vertices.append(Vertex(11, (-i, 24, -i + 1), (self.charged_current_coupling(), 0.0)))
            else:
                vertices.append(Vertex(10, (i, 24, i + 1), (self.charged_current_coupling(), 0.0)))
                vertices.append(Vertex(11, (-i, -24, -i - 1), (self.charged_current_coupling(), 0.0)))

    def _extend_electroweak_gauge_vertices(self, vertices: list[Vertex]) -> None:
        ngc = self.neutral_gauge_coupling()
        wc = self.weak_coupling()
        vertices.extend(
            [
                Vertex(12, (24, -24, 23), (-ngc, 0.0)),
                Vertex(12, (-24, 24, 23), (ngc, 0.0)),
                Vertex(12, (24, -24, 22), (-self.charge(24), 0.0)),
                Vertex(12, (-24, 24, 22), (self.charge(24), 0.0)),
                Vertex(12, (24, 22, 24), (self.charge(24), 0.0)),
                Vertex(12, (-24, 22, -24), (self.charge(-24), 0.0)),
                Vertex(12, (22, 24, 24), (-self.charge(24), 0.0)),
                Vertex(12, (22, -24, -24), (-self.charge(-24), 0.0)),
                Vertex(12, (24, 23, 24), (ngc, 0.0)),
                Vertex(12, (-24, 23, -24), (-ngc, 0.0)),
                Vertex(12, (23, 24, 24), (-ngc, 0.0)),
                Vertex(12, (23, -24, -24), (ngc, 0.0)),
                Vertex(13, (24, -24, -23), (wc, 0.0)),
                Vertex(13, (-24, 24, -23), (-wc, 0.0)),
                Vertex(14, (-23, 24, 24), (wc, 0.0)),
                Vertex(14, (-23, -24, -24), (-wc, 0.0)),
                Vertex(15, (24, -23, 24), (-wc, 0.0)),
                Vertex(15, (-24, -23, -24), (wc, 0.0)),
                Vertex(13, (24, 22, 26), (self.charge(24), 0.0)),
                Vertex(13, (-24, 22, -26), (-self.charge(-24), 0.0)),
                Vertex(13, (22, 24, 26), (-self.charge(24), 0.0)),
                Vertex(13, (22, -24, -26), (self.charge(-24), 0.0)),
                Vertex(13, (24, 23, 26), (ngc, 0.0)),
                Vertex(13, (-24, 23, -26), (ngc, 0.0)),
                Vertex(13, (23, 24, 26), (-ngc, 0.0)),
                Vertex(13, (23, -24, -26), (-ngc, 0.0)),
                Vertex(14, (26, 22, 24), (self.charge(26), 0.0)),
                Vertex(14, (-26, 22, -24), (-self.charge(-26), 0.0)),
                Vertex(14, (26, -24, 22), (-self.charge(26), 0.0)),
                Vertex(14, (-26, 24, 22), (self.charge(-26), 0.0)),
                Vertex(14, (26, 23, 24), (ngc, 0.0)),
                Vertex(14, (-26, 23, -24), (ngc, 0.0)),
                Vertex(14, (26, -24, 23), (-ngc, 0.0)),
                Vertex(14, (-26, 24, 23), (-ngc, 0.0)),
                Vertex(15, (22, 26, 24), (-self.charge(26), 0.0)),
                Vertex(15, (22, -26, -24), (self.charge(-26), 0.0)),
                Vertex(15, (24, -26, 22), (self.charge(24), 0.0)),
                Vertex(15, (-24, 26, 22), (-self.charge(-24), 0.0)),
                Vertex(15, (23, 26, 24), (-ngc, 0.0)),
                Vertex(15, (23, -26, -24), (-ngc, 0.0)),
                Vertex(15, (24, -26, 23), (ngc, 0.0)),
                Vertex(15, (-24, 26, 23), (ngc, 0.0)),
            ]
        )

    def _extend_higgs_vertices(self, vertices: list[Vertex]) -> None:
        wc = self.weak_coupling()
        wc2 = wc**2
        cw2 = self.cos_weak**2
        for i in range(1, 7):
            if self.mass(i) == 0.0:
                continue
            coupling = self.mass(i) * wc / (self.mass(24) * 2.0)
            vertices.append(Vertex(16, (i, 25, i), (coupling, 0.0)))
            vertices.append(Vertex(16, (-i, 25, -i), (coupling, 0.0)))
        vertices.extend(
            [
                Vertex(17, (24, -24, 25), (self.mass(24) * wc, 0.0)),
                Vertex(17, (-24, 24, 25), (self.mass(24) * wc, 0.0)),
                Vertex(18, (25, 24, 24), (self.mass(24) * wc, 0.0)),
                Vertex(18, (25, -24, -24), (self.mass(24) * wc, 0.0)),
                Vertex(19, (24, 25, 24), (self.mass(24) * wc, 0.0)),
                Vertex(19, (-24, 25, -24), (self.mass(24) * wc, 0.0)),
                Vertex(17, (23, 23, 25), (self.mass(23) * self.weak_coupling_over_cosine(), 0.0)),
                Vertex(18, (25, 23, 23), (self.mass(23) * self.weak_coupling_over_cosine(), 0.0)),
                Vertex(19, (23, 25, 23), (self.mass(23) * self.weak_coupling_over_cosine(), 0.0)),
                Vertex(20, (25, 25, 25), ((-3.0 / 2.0) * wc * (self.mass(25) ** 2 / self.mass(24)), 0.0)),
                Vertex(17, (24, -24, 127), (0.5 * wc2, 0.0)),
                Vertex(17, (-24, 24, 127), (0.5 * wc2, 0.0)),
                Vertex(17, (23, 23, 127), (0.5 * wc2 / cw2, 0.0)),
                Vertex(20, (125, 25, 25), (1.0, 0.0)),
                Vertex(20, (25, 125, 25), (1.0, 0.0)),
                Vertex(20, (25, 25, 125), ((-3.0 / 4.0) * wc2 * self.mass(25) ** 2 / self.mass(24) ** 2, 0.0)),
                Vertex(20, (127, 25, 25), (1.0, 0.0)),
                Vertex(20, (25, 127, 25), (1.0, 0.0)),
                Vertex(20, (25, 25, 126), (1.0, -10.0)),
                Vertex(18, (126, 23, 23), (-0.5 * wc2 / cw2, 0.0)),
                Vertex(19, (23, 126, 23), (-0.5 * wc2 / cw2, 0.0)),
                Vertex(18, (126, 24, 24), (-0.5 * wc2, 0.0)),
                Vertex(18, (126, -24, -24), (-0.5 * wc2, 0.0)),
                Vertex(19, (24, 126, 24), (-0.5 * wc2, 0.0)),
                Vertex(19, (-24, 126, -24), (-0.5 * wc2, 0.0)),
            ]
        )

    def _extend_lepton_vertices(self, vertices: list[Vertex]) -> None:
        for i in range(1, 4):
            charged = 11 + (2 * i - 2)
            neutrino = 12 + (2 * i - 2)
            vertices.extend(
                [
                    Vertex(21, (charged, -charged, 22), self.photon_fermion_coupling(charged)),
                    Vertex(22, (-charged, charged, 22), self.photon_fermion_coupling(-charged)),
                    Vertex(21, (charged, -charged, 23), self.z_fermion_coupling(charged)),
                    Vertex(22, (-charged, charged, 23), self.z_fermion_coupling(-charged)),
                    Vertex(21, (charged, -neutrino, -24), (self.charged_current_coupling(), 0.0)),
                    Vertex(22, (-charged, neutrino, 24), (self.charged_current_coupling(), 0.0)),
                    Vertex(10, (charged, 22, charged), self.photon_fermion_coupling(charged)),
                    Vertex(11, (-charged, 22, -charged), self.photon_fermion_coupling(-charged)),
                    Vertex(23, (22, charged, charged), self.photon_fermion_coupling(charged)),
                    Vertex(24, (22, -charged, -charged), self.photon_fermion_coupling(-charged)),
                    Vertex(10, (charged, 23, charged), self.z_fermion_coupling(charged)),
                    Vertex(11, (-charged, 23, -charged), self.z_fermion_coupling(-charged)),
                    Vertex(23, (23, charged, charged), self.z_fermion_coupling(charged)),
                    Vertex(24, (23, -charged, -charged), self.z_fermion_coupling(-charged)),
                ]
            )


def _flat_index(indices: tuple[int, ...], dims: tuple[int, ...]) -> int:
    index = 0
    for value, dim in zip(indices, dims, strict=True):
        index = index * dim + value
    return index


def _number(value: complex | float) -> Any:
    from symbolica import Expression

    return Expression.num(value)


def _as_expression(value: Any) -> Any:
    if isinstance(value, int | float | complex):
        return _number(value)
    return value


def _minkowski_square_expression(momentum: Sequence[Any]) -> Any:
    if len(momentum) != 4:
        raise ValueError("Minkowski momentum needs four components")
    p0, p1, p2, p3 = (_as_expression(value) for value in momentum)
    return p0 * p0 - p1 * p1 - p2 * p2 - p3 * p3


def _two_gluon_to_tensor_data() -> list[complex]:
    data = [0j] * (6 * 4 * 4)
    metric = (1.0, -1.0, -1.0, -1.0)
    for tensor_index, (i, j) in enumerate(_ANTISYM_PAIRS):
        data[_flat_index((tensor_index, i, j), (6, 4, 4))] = (
            1.0 / (metric[i] * metric[j])
        ) + 0j
        data[_flat_index((tensor_index, j, i), (6, 4, 4))] = (
            -1.0 / (metric[j] * metric[i])
        ) + 0j
    return data


def _tensor_gluon_to_gluon_data() -> list[complex]:
    data = [0j] * (6 * 4 * 4)
    prefactor = 0.5j
    metric = (1.0, -1.0, -1.0, -1.0)
    rows = {
        0: ((0, 1, 1), (1, 2, 1), (2, 3, 1)),
        1: ((0, 0, 1), (3, 2, 1), (4, 3, 1)),
        2: ((1, 0, 1), (3, 1, -1), (5, 3, 1)),
        3: ((2, 0, 1), (4, 1, -1), (5, 2, -1)),
    }
    for out, entries in rows.items():
        for tensor_index, gluon_index, sign in entries:
            data[_flat_index((tensor_index, gluon_index, out), (6, 4, 4))] = (
                sign * prefactor / metric[gluon_index]
            )
    return data


def _gluon_tensor_to_gluon_data() -> list[complex]:
    data = [0j] * (6 * 4 * 4)
    prefactor = 0.5j
    metric = (1.0, -1.0, -1.0, -1.0)
    rows = {
        0: ((1, 0, -1), (2, 1, -1), (3, 2, -1)),
        1: ((0, 0, -1), (2, 3, -1), (3, 4, -1)),
        2: ((0, 1, -1), (1, 3, 1), (3, 5, -1)),
        3: ((0, 2, -1), (1, 4, 1), (2, 5, 1)),
    }
    for out, entries in rows.items():
        for gluon_index, tensor_index, sign in entries:
            data[_flat_index((tensor_index, gluon_index, out), (6, 4, 4))] = (
                sign * prefactor / metric[gluon_index]
            )
    return data


def _quark_vector_weyl_data(*, chirality: int) -> list[complex]:
    data = [0j] * (2 * 4 * 2)
    prefactor = 1j / math.sqrt(2.0)
    metric = (1.0, -1.0, -1.0, -1.0)

    def add(q_in: int, vector: int, q_out: int, coefficient: complex) -> None:
        # spenso canonicalizes T(weyl, mink, weyl) storage as
        # (weyl_in, weyl_out, mink), while expression calls keep the original
        # slot order.
        data[_flat_index((q_in, q_out, vector), (2, 2, 4))] = (
            coefficient / metric[vector]
        )

    if chirality == 1:
        add(0, 0, 0, prefactor)
        add(0, 3, 0, -prefactor)
        add(1, 1, 0, -prefactor)
        add(1, 2, 0, -1j * prefactor)
        add(1, 0, 1, prefactor)
        add(1, 3, 1, prefactor)
        add(0, 1, 1, -prefactor)
        add(0, 2, 1, 1j * prefactor)
        return data
    if chirality == -1:
        add(0, 0, 0, prefactor)
        add(0, 3, 0, prefactor)
        add(1, 1, 0, prefactor)
        add(1, 2, 0, 1j * prefactor)
        add(1, 0, 1, prefactor)
        add(1, 3, 1, -prefactor)
        add(0, 1, 1, prefactor)
        add(0, 2, 1, -1j * prefactor)
        return data
    raise ValueError(f"unsupported Weyl chirality: {chirality}")


_ANTISYM_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))


__all__ = [
    "AmplicolSMLeadingColorModel",
    "Model",
    "Particle",
    "Vertex",
    "VertexLoweringRule",
]
