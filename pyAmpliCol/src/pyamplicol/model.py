from __future__ import annotations

import math
from collections import Counter
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


CouplingOrders = tuple[tuple[str, int], ...]


@dataclass(frozen=True)
class SourceSpinState:
    helicity: int
    chirality: int
    spin_state: int | tuple[int, ...]


@dataclass(frozen=True)
class QuantumFlow:
    chirality: int
    spin_state: int | tuple[int, ...]
    flavour_flow: tuple[int, ...]
    charge_flow: int
    coupling: tuple[float, float]


@dataclass(frozen=True)
class VertexLoweringRule:
    kind: int
    backend: str
    tensor_names: tuple[str, ...] = ()
    expression_head: str = ""
    full_tensor_network_ready: bool = False
    description: str = ""
    kernel: str = ""
    input_roles: tuple[str, str] = ("", "")
    output_role: str = ""
    coupling_mode: str = "none"

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "backend": self.backend,
            "tensor_names": list(self.tensor_names),
            "expression_head": self.expression_head,
            "full_tensor_network_ready": self.full_tensor_network_ready,
            "description": self.description,
            "kernel": self.kernel,
            "input_roles": list(self.input_roles),
            "output_role": self.output_role,
            "coupling_mode": self.coupling_mode,
        }


@dataclass(frozen=True)
class PropagatorLoweringRule:
    particle_id: int
    chirality: int
    backend: str
    full_tensor_network_ready: bool
    applies_propagator: bool
    kernel: str
    description: str = ""

    def to_json_dict(self) -> dict[str, object]:
        return {
            "particle_id": self.particle_id,
            "chirality": self.chirality,
            "backend": self.backend,
            "full_tensor_network_ready": self.full_tensor_network_ready,
            "applies_propagator": self.applies_propagator,
            "kernel": self.kernel,
            "description": self.description,
        }


@dataclass(frozen=True)
class VertexLoweringCoverageEntry:
    kind: int
    vertex_count: int
    backend: str
    full_tensor_network_ready: bool
    tensor_names: tuple[str, ...]
    expression_head: str
    description: str
    kernel: str
    input_roles: tuple[str, str]
    output_role: str
    coupling_mode: str

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "vertex_count": self.vertex_count,
            "backend": self.backend,
            "full_tensor_network_ready": self.full_tensor_network_ready,
            "tensor_names": list(self.tensor_names),
            "expression_head": self.expression_head,
            "description": self.description,
            "kernel": self.kernel,
            "input_roles": list(self.input_roles),
            "output_role": self.output_role,
            "coupling_mode": self.coupling_mode,
        }


@dataclass(frozen=True)
class VertexLoweringCoverageReport:
    model: str
    entries: tuple[VertexLoweringCoverageEntry, ...]

    @property
    def ready_kinds(self) -> tuple[int, ...]:
        return tuple(
            entry.kind for entry in self.entries if entry.full_tensor_network_ready
        )

    @property
    def pending_kinds(self) -> tuple[int, ...]:
        return tuple(
            entry.kind
            for entry in self.entries
            if entry.backend != "unimplemented" and not entry.full_tensor_network_ready
        )

    @property
    def unimplemented_kinds(self) -> tuple[int, ...]:
        return tuple(
            entry.kind for entry in self.entries if entry.backend == "unimplemented"
        )

    def to_json_dict(self) -> dict[str, object]:
        return {
            "model": self.model,
            "ready_kinds": list(self.ready_kinds),
            "pending_kinds": list(self.pending_kinds),
            "unimplemented_kinds": list(self.unimplemented_kinds),
            "entries": [entry.to_json_dict() for entry in self.entries],
        }


@dataclass
class Model:
    name: str
    particles: dict[int, Particle] = field(default_factory=dict)
    vertices: tuple[Vertex, ...] = ()

    @cached_property
    def _species_by_pdg(self) -> dict[int, Particle]:
        species: dict[int, Particle] = {}
        for particle in self.particles.values():
            species.setdefault(particle.pdg, particle)
            species.setdefault(particle.anti_pdg, particle)
        return species

    @cached_property
    def _property_sign_by_pdg(self) -> dict[int, int]:
        signs: dict[int, int] = {}
        for particle in self.particles.values():
            signs.setdefault(particle.pdg, 1)
            if particle.anti_pdg != particle.pdg:
                signs.setdefault(particle.anti_pdg, -1)
            else:
                signs.setdefault(particle.anti_pdg, 1)
        return signs

    @cached_property
    def _color_rep_by_pdg(self) -> dict[int, int]:
        reps: dict[int, int] = {}
        for pdg, particle in self._species_by_pdg.items():
            color = particle.color_rep
            if self._property_sign_by_pdg[pdg] < 0 and abs(color) == 3:
                reps[pdg] = -color
            else:
                reps[pdg] = color
        return reps

    @cached_property
    def _vertices_by_input(self) -> dict[tuple[str, int, int], tuple[Vertex, ...]]:
        return {}

    def particle(self, pdg: int) -> Particle:
        species = self._species_by_pdg.get(pdg)
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
        try:
            return self._color_rep_by_pdg[pdg]
        except KeyError as exc:
            raise KeyError(f"particle not in model: {pdg}") from exc

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

    def propagator_lowering_rule(
        self,
        particle_id: int,
        chirality: int = 0,
    ) -> PropagatorLoweringRule:
        if self.is_tensor(particle_id):
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="identity",
                full_tensor_network_ready=True,
                applies_propagator=False,
                kernel="auxiliary_tensor_embedded_propagator",
                description=(
                    "auxiliary-tensor propagator factors are embedded in the "
                    "adjacent AmpliCol vertex kernels"
                ),
            )
        if particle_id in (125, 126, 127):
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="identity",
                full_tensor_network_ready=True,
                applies_propagator=False,
                kernel="auxiliary_scalar_no_propagator",
                description=(
                    "Higgsor auxiliary scalar currents are non-propagating in "
                    "legacy AmpliCol"
                ),
            )
        if particle_id in (21, 22, 99):
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="symbolica",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="massless_vector_feynman_gauge",
                description="massless vector propagator in mostly-minus metric",
            )
        if abs(particle_id) == 24 or particle_id == 23:
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="symbolica",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="massive_vector_unitary_gauge",
                description="massive vector propagator with width",
            )
        if self.is_chiral_eligible(particle_id) and chirality != 0:
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="spenso",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="weyl_fermion",
                description="massless Weyl fermion propagator",
            )
        if self.is_fermion(particle_id) and self.mass(particle_id) != 0.0:
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="symbolica",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="massive_dirac_fermion",
                description="massive Dirac fermion propagator",
            )
        if self.is_higgs(particle_id):
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="symbolica",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="scalar_with_width",
                description="scalar propagator with optional width",
            )
        return PropagatorLoweringRule(
            particle_id=particle_id,
            chirality=chirality,
            backend="unimplemented",
            full_tensor_network_ready=False,
            applies_propagator=True,
            kernel="unknown",
            description="no pyamplicol propagator lowering is registered",
        )

    def source_spin_states(self, particle_id: int) -> tuple[SourceSpinState, ...]:
        if self.is_chiral_eligible(particle_id):
            return (
                SourceSpinState(helicity=-1, chirality=-1, spin_state=-1),
                SourceSpinState(helicity=1, chirality=1, spin_state=1),
            )
        spin = self.spin(particle_id)
        if spin == 1:
            return (SourceSpinState(helicity=0, chirality=0, spin_state=0),)
        if spin == 2:
            return (
                SourceSpinState(helicity=-1, chirality=0, spin_state=-1),
                SourceSpinState(helicity=1, chirality=0, spin_state=1),
            )
        if spin == 3:
            return (
                SourceSpinState(helicity=-1, chirality=0, spin_state=-1),
                SourceSpinState(helicity=0, chirality=0, spin_state=0),
                SourceSpinState(helicity=1, chirality=0, spin_state=1),
            )
        return (SourceSpinState(helicity=0, chirality=0, spin_state=0),)

    def allowed_quantum_flows(
        self,
        vertex: Vertex,
        left_index: Any,
        right_index: Any,
    ) -> tuple[QuantumFlow, ...]:
        result_particle = vertex.particles[2]
        chiralities = _model_vertex_result_chiralities(
            self,
            vertex,
            left_index,
            right_index,
        )
        flows: list[QuantumFlow] = []
        for chirality in chiralities:
            flows.append(
                QuantumFlow(
                    chirality=chirality,
                    spin_state=self.result_spin_state(result_particle, chirality),
                    flavour_flow=self.combine_flavour_flow(
                        result_particle,
                        left_index,
                        right_index,
                    ),
                    charge_flow=self.charge_units(result_particle),
                    coupling=vertex.coupling,
                )
            )
        return tuple(flows)

    def combine_flavour_flow(
        self,
        result_particle: int,
        left_index: Any,
        right_index: Any,
    ) -> tuple[int, ...]:
        left_pdg = _index_particle_id(left_index)
        right_pdg = _index_particle_id(right_index)
        left_flow = _index_flavour_flow(left_index)
        right_flow = _index_flavour_flow(right_index)

        if self.is_fermion(result_particle):
            if self.is_fermion(left_pdg):
                return _append_flavour_transition(left_flow, result_particle)
            if self.is_fermion(right_pdg):
                return _append_flavour_transition(right_flow, result_particle)

        if self.is_fermion(left_pdg) and self.is_fermion(right_pdg):
            return (*left_flow, *right_flow, result_particle)

        return (result_particle,)

    def result_spin_state(self, particle_id: int, chirality: int) -> int:
        if self.is_fermion(particle_id):
            return chirality
        return 0

    def current_dimension(self, particle_id: int, chirality: int = 0) -> int:
        if chirality != 0 and self.is_chiral_eligible(particle_id):
            return 2
        try:
            return self.dimension(particle_id)
        except KeyError:
            return 0

    def charge_units(self, particle_id: int) -> int:
        return int(round(3.0 * self.charge(particle_id)))

    def auxiliary_kind(self, particle_id: int) -> str | None:
        if self.is_tensor(particle_id):
            return "antisymmetric-tensor"
        return None

    def vertex_component_expression(
        self,
        kind: int,
        left: Sequence[Any],
        right: Sequence[Any],
        *,
        result_particle_id: int,
        result_chirality: int,
        left_chirality: int = 0,
        right_chirality: int = 0,
        coupling: tuple[Any, Any] = (1.0, 0.0),
        left_momentum: Sequence[Any] | None = None,
        right_momentum: Sequence[Any] | None = None,
    ) -> tuple[Any, ...]:
        raise NotImplementedError

    def propagator_component_expression(
        self,
        particle_id: int,
        value: Sequence[Any],
        momentum: Sequence[Any],
        *,
        chirality: int = 0,
    ) -> tuple[Any, ...]:
        raise NotImplementedError

    def iter_vertices(self, *, color_accuracy: str = "lc") -> tuple[Vertex, ...]:
        if color_accuracy == "lc":
            return tuple(self.vertices)
        if color_accuracy in {"nlc", "full"}:
            return tuple(self.vertices)
        raise ValueError(f"unknown colour accuracy: {color_accuracy}")

    def vertices_for_inputs(
        self,
        left_pdg: int,
        right_pdg: int,
        *,
        color_accuracy: str = "lc",
    ) -> tuple[Vertex, ...]:
        key = (color_accuracy, int(left_pdg), int(right_pdg))
        if key not in self._vertices_by_input:
            self._vertices_by_input[key] = tuple(
                vertex
                for vertex in self.iter_vertices(color_accuracy=color_accuracy)
                if vertex.particles[0] == left_pdg
                and vertex.particles[1] == right_pdg
            )
        return self._vertices_by_input[key]

    def vertices_accepting(
        self,
        left_pdg: int,
        right_pdg: int,
        *,
        color_accuracy: str = "lc",
    ) -> tuple[Vertex, ...]:
        """Return model vertices for a local current-current combination.

        This is the process-generic name used by the DAG compiler.  It is a
        thin alias around the model table lookup, but keeping the name at the
        model boundary prevents production code from classifying whole process
        families before asking which local interactions are allowed.
        """

        return self.vertices_for_inputs(
            left_pdg,
            right_pdg,
            color_accuracy=color_accuracy,
        )

    def skip_duplicate_vertex_orientation(self, vertex: Vertex) -> bool:
        """Return whether a mirrored table entry should be skipped by DAG sweeps.

        Generic DAG generation asks the model because duplicated orientations are
        a model-table convention, not a process-family rule.
        """

        return False

    def vertex_coupling_orders(self, vertex: Vertex) -> CouplingOrders:
        """Return model-generic coupling-order increments for one vertex.

        The keys intentionally mirror UFO-style coupling-order names.  They are
        used only as local model metadata, so DAG pruning can cap e.g. QCD or
        QED order without recognizing whole process families.
        """

        return (("QED", 1),)

    def combine_coupling_orders(
        self,
        left_index: Any,
        right_index: Any,
        vertex: Vertex,
    ) -> CouplingOrders:
        totals: dict[str, int] = {}
        for orders in (
            _index_coupling_orders(left_index),
            _index_coupling_orders(right_index),
            self.vertex_coupling_orders(vertex),
        ):
            for name, value in orders:
                totals[str(name).upper()] = totals.get(str(name).upper(), 0) + int(value)
        return tuple(sorted((name, value) for name, value in totals.items() if value))

    def current_allowed(self, index: Any) -> bool:
        """Return whether a generated current index is valid in this model."""

        try:
            particle_id = int(getattr(index, "particle_id"))
            chirality = int(getattr(index, "chirality", 0))
            if chirality != 0 and self.is_chiral_eligible(particle_id):
                return True
            return self.dimension(particle_id) > 0
        except (KeyError, TypeError, ValueError):
            return False

    def vertex_lowering_coverage(self) -> VertexLoweringCoverageReport:
        counts = Counter(vertex.kind for vertex in self.vertices)
        entries: list[VertexLoweringCoverageEntry] = []
        for kind, count in sorted(counts.items()):
            rule = self.vertex_lowering_rule(kind)
            entries.append(
                VertexLoweringCoverageEntry(
                    kind=kind,
                    vertex_count=count,
                    backend=rule.backend,
                    full_tensor_network_ready=rule.full_tensor_network_ready,
                    tensor_names=rule.tensor_names,
                    expression_head=rule.expression_head,
                    description=rule.description,
                    kernel=rule.kernel,
                    input_roles=rule.input_roles,
                    output_role=rule.output_role,
                    coupling_mode=rule.coupling_mode,
                )
            )
        return VertexLoweringCoverageReport(
            model=self.name,
            entries=tuple(entries),
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
        return self._species_by_pdg.get(pdg)

    def _property_sign(self, pdg: int) -> int:
        try:
            return self._property_sign_by_pdg[pdg]
        except KeyError as exc:
            raise KeyError(f"particle not in model: {pdg}") from exc


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

    def skip_duplicate_vertex_orientation(self, vertex: Vertex) -> bool:
        """Skip mirrored model-table entries already covered by DAG sweeps."""

        return False

    def vertex_coupling_orders(self, vertex: Vertex) -> CouplingOrders:
        """Classify AmpliCol SM vertices by UFO-style coupling order."""

        if vertex.kind in {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}:
            return (("QCD", 1),)
        return (("QED", 1),)

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
                kernel="three_vector_current",
                input_roles=("vector", "vector"),
                output_role="vector",
                coupling_mode="fixed",
            )
        if kind == 1:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(str(symbols.two_gluon_to_tensor),),
                expression_head=str(symbols.two_gluon_to_tensor),
                full_tensor_network_ready=True,
                description="two gluons to auxiliary antisymmetric tensor",
                kernel="two_vector_to_tensor",
                input_roles=("vector", "vector"),
                output_role="antisymmetric_tensor",
                coupling_mode="fixed",
            )
        if kind == 2:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(str(symbols.tensor_gluon_to_gluon),),
                expression_head=str(symbols.tensor_gluon_to_gluon),
                full_tensor_network_ready=True,
                description="auxiliary tensor and gluon to gluon current",
                kernel="tensor_vector_to_vector",
                input_roles=("antisymmetric_tensor", "vector"),
                output_role="vector",
                coupling_mode="fixed",
            )
        if kind == 3:
            return VertexLoweringRule(
                kind=kind,
                backend="spenso",
                tensor_names=(str(symbols.gluon_tensor_to_gluon),),
                expression_head=str(symbols.gluon_tensor_to_gluon),
                full_tensor_network_ready=True,
                description="gluon and auxiliary tensor to gluon current",
                kernel="vector_tensor_to_vector",
                input_roles=("vector", "antisymmetric_tensor"),
                output_role="vector",
                coupling_mode="fixed",
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
                kernel="fermion_vector_to_fermion",
                input_roles=("fermion", "vector"),
                output_role="fermion",
                coupling_mode="fixed",
            )
        if kind in {4, 5, 7, 9}:
            qcd_role_map: dict[int, tuple[tuple[str, str], str, str]] = {
                4: (("vector", "fermion"), "fermion", "vector-fermion current"),
                5: (
                    ("vector", "antifermion"),
                    "antifermion",
                    "vector-antifermion current",
                ),
                7: (
                    ("antifermion", "vector"),
                    "antifermion",
                    "antifermion-vector current",
                ),
                9: (
                    ("antifermion", "fermion"),
                    "vector",
                    "antifermion-fermion vector current",
                ),
            }
            input_roles, output_role, description = qcd_role_map[kind]
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica",
                expression_head="quark_gluon_weyl_current",
                full_tensor_network_ready=True,
                description=f"Weyl QCD {description}",
                kernel=(
                    "fermion_pair_to_vector"
                    if kind == 9
                    else "fermion_vector_to_fermion"
                ),
                input_roles=input_roles,
                output_role=output_role,
                coupling_mode="fixed",
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
                kernel="fermion_vector_to_fermion",
                input_roles=("fermion", "vector"),
                output_role="fermion",
                coupling_mode="vertex",
            )
        if kind in {11, 21, 22, 23, 24}:
            fermion_gauge_role_map: dict[
                int,
                tuple[tuple[str, str], str, str, str],
            ] = {
                11: (
                    ("antifermion", "vector"),
                    "antifermion",
                    "antifermion electroweak current",
                    "fermion_vector_to_fermion",
                ),
                21: (
                    ("fermion", "antifermion"),
                    "vector",
                    "lepton-antilepton electroweak current",
                    "fermion_pair_to_vector",
                ),
                22: (
                    ("antifermion", "fermion"),
                    "vector",
                    "antilepton-lepton electroweak current",
                    "fermion_pair_to_vector",
                ),
                23: (
                    ("vector", "fermion"),
                    "fermion",
                    "vector-fermion electroweak current",
                    "fermion_vector_to_fermion",
                ),
                24: (
                    ("vector", "antifermion"),
                    "antifermion",
                    "vector-antifermion electroweak current",
                    "fermion_vector_to_fermion",
                ),
            }
            input_roles, output_role, description, kernel = fermion_gauge_role_map[kind]
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica",
                expression_head="fermion_gauge_weyl_current",
                full_tensor_network_ready=True,
                description=f"Weyl {description} with runtime coupling",
                kernel=kernel,
                input_roles=input_roles,
                output_role=output_role,
                coupling_mode="vertex",
            )
        if kind == 8:
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica",
                expression_head="qcd_u1_subtraction_current",
                full_tensor_network_ready=True,
                description="QCD U(1) subtraction current",
                kernel="fermion_pair_to_vector",
                input_roles=("fermion", "antifermion"),
                output_role="vector",
                coupling_mode="vertex",
            )
        if kind in {12, 13, 14, 15}:
            vector_role_map: dict[int, tuple[str, tuple[str, str], str, str]] = {
                12: (
                    "three_vector_current",
                    ("vector", "vector"),
                    "vector",
                    "electroweak three-vector current",
                ),
                13: (
                    "two_vector_to_tensor",
                    ("vector", "vector"),
                    "antisymmetric_tensor",
                    "two electroweak vectors to auxiliary tensor",
                ),
                14: (
                    "tensor_vector_to_vector",
                    ("antisymmetric_tensor", "vector"),
                    "vector",
                    "auxiliary tensor and vector to electroweak vector",
                ),
                15: (
                    "vector_tensor_to_vector",
                    ("vector", "antisymmetric_tensor"),
                    "vector",
                    "electroweak vector and auxiliary tensor to vector",
                ),
            }
            kernel, input_roles, output_role, description = vector_role_map[kind]
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica",
                expression_head=kernel,
                full_tensor_network_ready=True,
                description=f"{description} with runtime coupling",
                kernel=kernel,
                input_roles=input_roles,
                output_role=output_role,
                coupling_mode="vertex",
            )
        if kind == 16:
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica",
                expression_head="fermion_scalar_to_fermion",
                description=(
                    "massive Dirac fermion-scalar Yukawa current"
                ),
                full_tensor_network_ready=True,
                kernel="fermion_scalar_to_fermion",
                input_roles=("fermion", "scalar"),
                output_role="fermion",
                coupling_mode="vertex",
            )
        if kind in {17, 18, 19, 20}:
            scalar_role_map: dict[int, tuple[str, tuple[str, str], str, str]] = {
                17: (
                    "two_vector_to_scalar",
                    ("vector", "vector"),
                    "scalar",
                    "two vectors to scalar current",
                ),
                18: (
                    "scalar_vector_to_vector",
                    ("scalar", "vector"),
                    "vector",
                    "scalar-vector to vector current",
                ),
                19: (
                    "vector_scalar_to_vector",
                    ("vector", "scalar"),
                    "vector",
                    "vector-scalar to vector current",
                ),
                20: (
                    "two_scalar_to_scalar",
                    ("scalar", "scalar"),
                    "scalar",
                    "scalar-scalar to scalar current",
                ),
            }
            kernel, input_roles, output_role, description = scalar_role_map[kind]
            return VertexLoweringRule(
                kind=kind,
                backend="symbolica",
                expression_head=kernel,
                full_tensor_network_ready=True,
                description=f"{description} with runtime coupling",
                kernel=kernel,
                input_roles=input_roles,
                output_role=output_role,
                coupling_mode="vertex",
            )
        return VertexLoweringRule(
            kind=kind,
            backend="unimplemented",
            description="no native pyamplicol lowering rule is registered yet",
            kernel="unknown",
        )

    def vertex_component_expression(
        self,
        kind: int,
        left: Sequence[Any],
        right: Sequence[Any],
        *,
        result_particle_id: int,
        result_chirality: int,
        left_chirality: int = 0,
        right_chirality: int = 0,
        coupling: tuple[Any, Any] = (1.0, 0.0),
        left_momentum: Sequence[Any] | None = None,
        right_momentum: Sequence[Any] | None = None,
    ) -> tuple[Any, ...]:
        """Lower a local model vertex into component expressions.

        The inputs are already component tuples for the two parent currents.
        This method is intentionally process-blind: all decisions are local to
        the vertex kind, chirality labels, particle id, coupling, and optional
        current momenta.
        """

        if kind == 0:
            if left_momentum is None or right_momentum is None:
                raise ValueError("three-vector current requires parent momenta")
            return _expr_three_vector_current(
                tuple(left),
                tuple(left_momentum),
                tuple(right),
                tuple(right_momentum),
            )
        if kind == 1:
            return _expr_two_vector_to_tensor(tuple(left), tuple(right))
        if kind == 2:
            return _expr_tensor_vector_to_vector(tuple(left), tuple(right))
        if kind == 3:
            return _expr_vector_tensor_to_vector(tuple(left), tuple(right))
        if kind == 12:
            if left_momentum is None or right_momentum is None:
                raise ValueError("three-vector current requires parent momenta")
            return _expr_three_vector_current_coupled(
                tuple(left),
                tuple(left_momentum),
                tuple(right),
                tuple(right_momentum),
                coupling,
            )
        if kind == 13:
            return tuple(
                coupling[0] * component
                for component in _expr_two_vector_to_tensor(tuple(left), tuple(right))
            )
        if kind == 14:
            return tuple(
                coupling[0] * component
                for component in _expr_tensor_vector_to_vector(
                    tuple(left),
                    tuple(right),
                )
            )
        if kind == 15:
            return tuple(
                coupling[0] * component
                for component in _expr_vector_tensor_to_vector(
                    tuple(left),
                    tuple(right),
                )
            )
        if kind == 16:
            return _expr_fermion_scalar_to_fermion(
                tuple(left),
                tuple(right),
                coupling,
            )
        if kind == 17:
            return (
                (1j / math.sqrt(2.0))
                * coupling[0]
                * _expr_minkowski_dot(tuple(left), tuple(right)),
            )
        if kind == 18:
            return tuple(
                (1j / math.sqrt(2.0)) * coupling[0] * left[0] * component
                for component in tuple(right)
            )
        if kind == 19:
            return tuple(
                (1j / math.sqrt(2.0)) * coupling[0] * right[0] * component
                for component in tuple(left)
            )
        if kind == 20:
            phase = 1j if coupling[1] == -10.0 else 1.0
            return (
                (1j / math.sqrt(2.0))
                * phase
                * coupling[0]
                * left[0]
                * right[0],
            )
        if kind in {4, 6}:
            fermion, vector = (
                (tuple(right), tuple(left)) if kind == 4 else (tuple(left), tuple(right))
            )
            if len(fermion) == 4:
                return _expr_fermion_vector_dirac(
                    fermion,
                    vector,
                    antifermion=False,
                    coupling=None,
                )
            return _expr_fermion_vector_weyl(
                fermion,
                vector,
                result_chirality,
                antifermion=False,
                coupling=None,
            )
        if kind in {5, 7}:
            antifermion, vector = (
                (tuple(right), tuple(left)) if kind == 5 else (tuple(left), tuple(right))
            )
            if len(antifermion) == 4:
                return _expr_fermion_vector_dirac(
                    antifermion,
                    vector,
                    antifermion=True,
                    coupling=None,
                )
            return _expr_fermion_vector_weyl(
                antifermion,
                vector,
                result_chirality,
                antifermion=True,
                coupling=None,
            )
        if kind == 9:
            if len(tuple(right)) == 4 and len(tuple(left)) == 4:
                return _expr_fermion_antifermion_to_vector_dirac(
                    fermion=tuple(right),
                    antifermion=tuple(left),
                    coupling=(1.0, 1.0),
                )
            return _expr_fermion_antifermion_to_vector_weyl(
                fermion=tuple(right),
                antifermion=tuple(left),
                coupling=(1.0, 1.0),
                fermion_chirality=right_chirality,
                antifermion_chirality=left_chirality,
            )
        if kind == 8:
            qcd_coupling = (coupling[0], coupling[0])
            if len(tuple(left)) == 4 and len(tuple(right)) == 4:
                return _expr_fermion_antifermion_to_vector_dirac(
                    fermion=tuple(left),
                    antifermion=tuple(right),
                    coupling=qcd_coupling,
                )
            return _expr_fermion_antifermion_to_vector_weyl(
                fermion=tuple(left),
                antifermion=tuple(right),
                coupling=qcd_coupling,
                fermion_chirality=left_chirality,
                antifermion_chirality=right_chirality,
            )
        if kind in {10, 23}:
            fermion, vector = (
                (tuple(left), tuple(right))
                if kind == 10
                else (tuple(right), tuple(left))
            )
            if len(fermion) == 4:
                return _expr_fermion_vector_dirac(
                    fermion,
                    vector,
                    antifermion=False,
                    coupling=coupling,
                )
            return _expr_fermion_vector_weyl(
                fermion,
                vector,
                result_chirality,
                antifermion=False,
                coupling=coupling,
            )
        if kind in {11, 24}:
            antifermion, vector = (
                (tuple(left), tuple(right))
                if kind == 11
                else (tuple(right), tuple(left))
            )
            if len(antifermion) == 4:
                return _expr_fermion_vector_dirac(
                    antifermion,
                    vector,
                    antifermion=True,
                    coupling=coupling,
                )
            return _expr_fermion_vector_weyl(
                antifermion,
                vector,
                result_chirality,
                antifermion=True,
                coupling=coupling,
            )
        if kind == 21:
            return _expr_fermion_antifermion_to_vector_weyl(
                fermion=tuple(left),
                antifermion=tuple(right),
                coupling=coupling,
                fermion_chirality=left_chirality,
                antifermion_chirality=right_chirality,
            )
        if kind == 22:
            return _expr_fermion_antifermion_to_vector_weyl(
                fermion=tuple(right),
                antifermion=tuple(left),
                coupling=coupling,
                fermion_chirality=right_chirality,
                antifermion_chirality=left_chirality,
            )
        raise ValueError(f"vertex kind {kind} has no component expression lowering")

    def propagator_component_expression(
        self,
        particle_id: int,
        value: Sequence[Any],
        momentum: Sequence[Any],
        *,
        chirality: int = 0,
    ) -> tuple[Any, ...]:
        rule = self.propagator_lowering_rule(particle_id, chirality)
        components = tuple(value)
        current_momentum = tuple(momentum)
        if not rule.applies_propagator:
            return components
        if rule.kernel == "massless_vector_feynman_gauge":
            denominator = _minkowski_square_expression(current_momentum)
            prefactor = -1j / denominator
            return tuple(component * prefactor for component in components)
        if rule.kernel == "massive_vector_unitary_gauge":
            mass = self.mass(particle_id)
            width = self.width(particle_id)
            denominator = (
                _minkowski_square_expression(current_momentum)
                - mass * mass
                + 1j * mass * width
            )
            prefactor = -1j / denominator
            longitudinal = (
                _expr_minkowski_dot(components, current_momentum)
                / (mass * mass)
            )
            return tuple(
                (components[index] - current_momentum[index] * longitudinal)
                * prefactor
                for index in range(4)
            )
        if rule.kernel == "weyl_fermion":
            if particle_id < 0:
                return _expr_antiquark_propagator_weyl(
                    components,
                    current_momentum,
                    chirality,
                )
            return _expr_quark_propagator_weyl(
                components,
                current_momentum,
                chirality,
            )
        if rule.kernel == "massive_dirac_fermion":
            if particle_id < 0:
                return _expr_antiquark_propagator_dirac(
                    components,
                    current_momentum,
                    self.mass(particle_id),
                    self.width(particle_id),
                )
            return _expr_quark_propagator_dirac(
                components,
                current_momentum,
                self.mass(particle_id),
                self.width(particle_id),
            )
        if rule.kernel == "scalar_with_width":
            mass = self.mass(particle_id)
            width = self.width(particle_id)
            denominator = (
                _minkowski_square_expression(current_momentum)
                - mass * mass
                + 1j * mass * width
            )
            prefactor = 1j / denominator
            return tuple(component * prefactor for component in components)
        raise ValueError(
            f"propagator kernel {rule.kernel!r} is not lowered for particle "
            f"{particle_id}"
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


def _model_vertex_result_chiralities(
    model: Model,
    vertex: Vertex,
    left_index: Any,
    right_index: Any,
) -> tuple[int, ...]:
    result_pdg = vertex.particles[2]
    if model.is_chiral_eligible(result_pdg):
        input_chirality = _model_fermion_input_chirality(model, left_index, right_index)
        if input_chirality == 0:
            return (0,)
        if not _model_weyl_vertex_allowed(vertex, result_pdg, input_chirality):
            return ()
        return (input_chirality,)

    if _model_is_fermion_pair_to_vector_vertex(vertex.kind):
        left_pdg = _index_particle_id(left_index)
        right_pdg = _index_particle_id(right_index)
        left_chirality = (
            _index_chirality(left_index) if model.is_fermion(left_pdg) else 0
        )
        right_chirality = (
            _index_chirality(right_index) if model.is_fermion(right_pdg) else 0
        )
        if (
            left_chirality != 0
            and right_chirality != 0
            and left_chirality != -right_chirality
        ):
            return ()
        if not _model_fermion_pair_vector_coupling_allowed(
            vertex,
            left_chirality,
            right_chirality,
        ):
            return ()
    return (0,)


def _model_fermion_input_chirality(
    model: Model,
    left_index: Any,
    right_index: Any,
) -> int:
    left_pdg = _index_particle_id(left_index)
    right_pdg = _index_particle_id(right_index)
    left_chirality = _index_chirality(left_index)
    right_chirality = _index_chirality(right_index)
    if model.is_fermion(left_pdg) and left_chirality != 0:
        return left_chirality
    if model.is_fermion(right_pdg) and right_chirality != 0:
        return right_chirality
    return 0


def _model_weyl_vertex_allowed(
    vertex: Vertex,
    result_pdg: int,
    chirality: int,
) -> bool:
    if vertex.kind == 16:
        return False
    if vertex.kind in {10, 11, 23, 24}:
        index = _model_fermion_coupling_index(result_pdg, chirality)
        return vertex.coupling[index] != 0.0
    return True


def _model_fermion_coupling_index(pdg: int, chirality: int) -> int:
    if pdg > 0:
        return 0 if chirality == -1 else 1
    return 0 if chirality == 1 else 1


def _model_is_fermion_pair_to_vector_vertex(kind: int) -> bool:
    return kind in {8, 9, 21, 22}


def _model_fermion_pair_vector_coupling_allowed(
    vertex: Vertex,
    left_chirality: int,
    right_chirality: int,
) -> bool:
    if vertex.kind in {8, 9}:
        return True
    if left_chirality == 0 or right_chirality == 0:
        return any(component != 0.0 for component in vertex.coupling)
    if vertex.kind == 21:
        index = 0 if left_chirality == -1 and right_chirality == 1 else 1
        return vertex.coupling[index] != 0.0
    if vertex.kind == 22:
        index = 0 if left_chirality == 1 and right_chirality == -1 else 1
        return vertex.coupling[index] != 0.0
    return True


def _index_particle_id(index: Any) -> int:
    if hasattr(index, "particle_id"):
        return int(getattr(index, "particle_id"))
    return int(getattr(index, "pdg"))


def _index_chirality(index: Any) -> int:
    return int(getattr(index, "chirality", 0))


def _index_flavour_flow(index: Any) -> tuple[int, ...]:
    flow = getattr(index, "flavour_flow", None)
    if flow is None:
        return (_index_particle_id(index),)
    return tuple(int(value) for value in flow)


def _index_coupling_orders(index: Any) -> CouplingOrders:
    orders = getattr(index, "coupling_orders", None)
    if orders is None:
        return ()
    return tuple(
        sorted(
            (str(name).upper(), int(value))
            for name, value in orders
            if int(value) != 0
        )
    )


def _append_flavour_transition(
    flow: tuple[int, ...],
    result_particle: int,
) -> tuple[int, ...]:
    if flow and flow[-1] == result_particle:
        return flow
    return (*flow, result_particle)


def _expr_three_vector_current(
    left: tuple[Any, ...],
    left_momentum: tuple[Any, ...],
    right: tuple[Any, ...],
    right_momentum: tuple[Any, ...],
) -> tuple[Any, ...]:
    dot = _expr_minkowski_dot(left, right)
    left_dot_right_momentum = _expr_minkowski_dot(left, right_momentum)
    right_dot_left_momentum = _expr_minkowski_dot(right, left_momentum)
    prefactor = 1j / math.sqrt(2.0)
    return tuple(
        prefactor
        * (
            dot * (left_momentum[index] - right_momentum[index])
            + 2.0
            * (
                left_dot_right_momentum * right[index]
                - right_dot_left_momentum * left[index]
            )
        )
        for index in range(4)
    )


def _expr_three_vector_current_coupled(
    left: tuple[Any, ...],
    left_momentum: tuple[Any, ...],
    right: tuple[Any, ...],
    right_momentum: tuple[Any, ...],
    coupling: tuple[Any, Any],
) -> tuple[Any, ...]:
    dot = _expr_minkowski_dot(left, right)
    tmp2_momentum = tuple(
        2.0 * right_momentum[index] + left_momentum[index]
        for index in range(4)
    )
    tmp3_momentum = tuple(
        -2.0 * left_momentum[index] - right_momentum[index]
        for index in range(4)
    )
    tmp2 = _expr_minkowski_dot(left, tmp2_momentum)
    tmp3 = _expr_minkowski_dot(right, tmp3_momentum)
    prefactor = (1j / math.sqrt(2.0)) * coupling[0]
    return tuple(
        prefactor
        * (
            dot * (left_momentum[index] - right_momentum[index])
            + tmp2 * right[index]
            + tmp3 * left[index]
        )
        for index in range(4)
    )


def _expr_two_vector_to_tensor(
    left: tuple[Any, ...],
    right: tuple[Any, ...],
) -> tuple[Any, ...]:
    return (
        left[0] * right[1] - left[1] * right[0],
        left[0] * right[2] - left[2] * right[0],
        left[0] * right[3] - left[3] * right[0],
        left[1] * right[2] - left[2] * right[1],
        left[1] * right[3] - left[3] * right[1],
        left[2] * right[3] - left[3] * right[2],
    )


def _expr_tensor_vector_to_vector(
    tensor: tuple[Any, ...],
    vector: tuple[Any, ...],
) -> tuple[Any, ...]:
    prefactor = 0.5j
    return (
        (tensor[0] * vector[1] + tensor[1] * vector[2] + tensor[2] * vector[3])
        * prefactor,
        (tensor[0] * vector[0] + tensor[3] * vector[2] + tensor[4] * vector[3])
        * prefactor,
        (tensor[1] * vector[0] - tensor[3] * vector[1] + tensor[5] * vector[3])
        * prefactor,
        (tensor[2] * vector[0] - tensor[4] * vector[1] - tensor[5] * vector[2])
        * prefactor,
    )


def _expr_vector_tensor_to_vector(
    vector: tuple[Any, ...],
    tensor: tuple[Any, ...],
) -> tuple[Any, ...]:
    prefactor = 0.5j
    return (
        (-vector[1] * tensor[0] - vector[2] * tensor[1] - vector[3] * tensor[2])
        * prefactor,
        (-vector[0] * tensor[0] - vector[2] * tensor[3] - vector[3] * tensor[4])
        * prefactor,
        (-vector[0] * tensor[1] + vector[1] * tensor[3] - vector[3] * tensor[5])
        * prefactor,
        (-vector[0] * tensor[2] + vector[1] * tensor[4] + vector[2] * tensor[5])
        * prefactor,
    )


def _expr_fermion_vector_weyl(
    fermion: tuple[Any, ...],
    vector: tuple[Any, ...],
    chirality: int,
    *,
    antifermion: bool,
    coupling: tuple[Any, Any] | None,
) -> tuple[Any, ...]:
    tmp1, tmp2, tmp3, tmp4 = _expr_vector_slash_terms(vector)
    prefactor = 1j / math.sqrt(2.0)
    f1, f2 = fermion
    if antifermion:
        if chirality == 1:
            factor = prefactor if coupling is None else prefactor * coupling[0]
            return (
                factor * (tmp1 * f1 + tmp4 * f2),
                factor * (tmp2 * f2 + tmp3 * f1),
            )
        if chirality == -1:
            factor = prefactor if coupling is None else prefactor * coupling[1]
            return (
                factor * (tmp2 * f1 - tmp4 * f2),
                factor * (tmp1 * f2 - tmp3 * f1),
            )
    else:
        if chirality == 1:
            factor = prefactor if coupling is None else prefactor * coupling[1]
            return (
                factor * (tmp2 * f1 - tmp3 * f2),
                factor * (tmp1 * f2 - tmp4 * f1),
            )
        if chirality == -1:
            factor = prefactor if coupling is None else prefactor * coupling[0]
            return (
                factor * (tmp1 * f1 + tmp3 * f2),
                factor * (tmp2 * f2 + tmp4 * f1),
            )
    raise ValueError("fermion-vector Weyl kernel needs nonzero chirality")


def _expr_fermion_vector_dirac(
    fermion: tuple[Any, ...],
    vector: tuple[Any, ...],
    *,
    antifermion: bool,
    coupling: tuple[Any, Any] | None,
) -> tuple[Any, ...]:
    if len(fermion) != 4 or len(vector) != 4:
        raise ValueError("Dirac fermion-vector current expects dimensions 4 and 4")
    tmp1, tmp2, tmp3, tmp4 = _expr_vector_slash_terms(vector)
    prefactor = 1j / math.sqrt(2.0)
    f1, f2, f3, f4 = fermion
    if coupling is None:
        left_coupling = 1.0
        right_coupling = 1.0
    else:
        left_coupling, right_coupling = coupling
    if antifermion:
        upper = prefactor * right_coupling
        lower = prefactor * left_coupling
        return (
            upper * (tmp2 * f3 - tmp4 * f4),
            upper * (tmp1 * f4 - tmp3 * f3),
            lower * (tmp1 * f1 + tmp4 * f2),
            lower * (tmp2 * f2 + tmp3 * f1),
        )
    upper = prefactor * left_coupling
    lower = prefactor * right_coupling
    return (
        upper * (tmp1 * f3 + tmp3 * f4),
        upper * (tmp2 * f4 + tmp4 * f3),
        lower * (tmp2 * f1 - tmp3 * f2),
        lower * (tmp1 * f2 - tmp4 * f1),
    )


def _expr_fermion_antifermion_to_vector_weyl(
    *,
    fermion: tuple[Any, ...],
    antifermion: tuple[Any, ...],
    coupling: tuple[Any, Any],
    fermion_chirality: int,
    antifermion_chirality: int,
) -> tuple[Any, ...]:
    prefactor = 1j / math.sqrt(2.0)
    left, right = coupling
    f1, f2 = fermion
    a1, a2 = antifermion
    if fermion_chirality == -1 and antifermion_chirality == 1:
        factor = prefactor * left
        return (
            factor * (f1 * a1 + f2 * a2),
            -factor * (f2 * a1 + f1 * a2),
            1j * factor * (-f2 * a1 + f1 * a2),
            factor * (-f1 * a1 + f2 * a2),
        )
    if fermion_chirality == 1 and antifermion_chirality == -1:
        factor = prefactor * right
        return (
            factor * (f1 * a1 + f2 * a2),
            factor * (f1 * a2 + f2 * a1),
            1j * factor * (-f1 * a2 + f2 * a1),
            factor * (f1 * a1 - f2 * a2),
        )
    return (0j, 0j, 0j, 0j)


def _expr_fermion_antifermion_to_vector_dirac(
    *,
    fermion: tuple[Any, ...],
    antifermion: tuple[Any, ...],
    coupling: tuple[Any, Any],
) -> tuple[Any, ...]:
    if len(fermion) != 4 or len(antifermion) != 4:
        raise ValueError("Dirac fermion-antifermion vector current expects dimensions 4 and 4")
    prefactor = 1j / math.sqrt(2.0)
    left_coupling, right_coupling = coupling
    f1, f2, f3, f4 = fermion
    a1, a2, a3, a4 = antifermion
    left = (
        f3 * a1 + f4 * a2,
        -(f4 * a1 + f3 * a2),
        1j * (-f4 * a1 + f3 * a2),
        -f3 * a1 + f4 * a2,
    )
    right = (
        f1 * a3 + f2 * a4,
        f1 * a4 + f2 * a3,
        1j * (-f1 * a4 + f2 * a3),
        f1 * a3 - f2 * a4,
    )
    return tuple(
        prefactor * (left_coupling * left[index] + right_coupling * right[index])
        for index in range(4)
    )


def _expr_fermion_scalar_to_fermion(
    fermion: tuple[Any, ...],
    scalar: tuple[Any, ...],
    coupling: tuple[Any, Any],
) -> tuple[Any, ...]:
    if len(fermion) != 4 or len(scalar) != 1:
        raise ValueError("Dirac fermion-scalar current expects dimensions 4 and 1")
    prefactor = -1j / math.sqrt(2.0)
    return tuple(prefactor * coupling[0] * scalar[0] * component for component in fermion)


def _expr_vector_slash_terms(vector: tuple[Any, ...]) -> tuple[Any, Any, Any, Any]:
    v0, v1, v2, v3 = vector
    return v0 + v3, v0 - v3, v1 + 1j * v2, v1 - 1j * v2


def _expr_quark_propagator_weyl(
    quark: tuple[Any, ...],
    momentum: tuple[Any, ...],
    chirality: int,
) -> tuple[Any, ...]:
    energy, px, py, pz = momentum
    denominator = _minkowski_square_expression(momentum)
    prefactor = 1j / denominator
    tmp1 = energy + pz
    tmp2 = energy - pz
    tmp3 = px + 1j * py
    tmp4 = px - 1j * py
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
    raise ValueError("Weyl quark propagator expression needs nonzero chirality")


def _expr_antiquark_propagator_weyl(
    antiquark: tuple[Any, ...],
    momentum: tuple[Any, ...],
    chirality: int,
) -> tuple[Any, ...]:
    energy, px, py, pz = momentum
    denominator = _minkowski_square_expression(momentum)
    prefactor = 1j / denominator
    tmp1 = -(energy + pz)
    tmp2 = -(energy - pz)
    tmp3 = -(px + 1j * py)
    tmp4 = -(px - 1j * py)
    a1, a2 = antiquark
    if chirality == 1:
        return (
            (tmp2 * a1 - tmp4 * a2) * prefactor,
            (tmp1 * a2 - tmp3 * a1) * prefactor,
        )
    if chirality == -1:
        return (
            (tmp1 * a1 + tmp4 * a2) * prefactor,
            (tmp2 * a2 + tmp3 * a1) * prefactor,
        )
    raise ValueError("Weyl antiquark propagator expression needs nonzero chirality")


def _expr_quark_propagator_dirac(
    quark: tuple[Any, ...],
    momentum: tuple[Any, ...],
    mass: float,
    width: float,
) -> tuple[Any, ...]:
    if len(quark) != 4 or len(momentum) != 4:
        raise ValueError("Dirac quark propagator expects four components")
    energy, px, py, pz = momentum
    denominator = _minkowski_square_expression(momentum) - mass * mass + 1j * mass * width
    prefactor = 1j / denominator
    tmp1 = energy + pz
    tmp2 = energy - pz
    tmp3 = px + 1j * py
    tmp4 = px - 1j * py
    q1, q2, q3, q4 = quark
    return (
        (tmp1 * q3 + tmp3 * q4 + mass * q1) * prefactor,
        (tmp2 * q4 + tmp4 * q3 + mass * q2) * prefactor,
        (tmp2 * q1 - tmp3 * q2 + mass * q3) * prefactor,
        (tmp1 * q2 - tmp4 * q1 + mass * q4) * prefactor,
    )


def _expr_antiquark_propagator_dirac(
    antiquark: tuple[Any, ...],
    momentum: tuple[Any, ...],
    mass: float,
    width: float,
) -> tuple[Any, ...]:
    if len(antiquark) != 4 or len(momentum) != 4:
        raise ValueError("Dirac antiquark propagator expects four components")
    energy, px, py, pz = momentum
    denominator = _minkowski_square_expression(momentum) - mass * mass + 1j * mass * width
    prefactor = 1j / denominator
    tmp1 = -(energy + pz)
    tmp2 = -(energy - pz)
    tmp3 = -(px + 1j * py)
    tmp4 = -(px - 1j * py)
    a1, a2, a3, a4 = antiquark
    return (
        (tmp2 * a3 - tmp4 * a4 + mass * a1) * prefactor,
        (tmp1 * a4 - tmp3 * a3 + mass * a2) * prefactor,
        (tmp1 * a1 + tmp4 * a2 + mass * a3) * prefactor,
        (tmp2 * a2 + tmp3 * a1 + mass * a4) * prefactor,
    )


def _expr_minkowski_dot(
    left: tuple[Any, ...],
    right: tuple[Any, ...],
) -> Any:
    return left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3]


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
    "CouplingOrders",
    "Model",
    "Particle",
    "PropagatorLoweringRule",
    "QuantumFlow",
    "SourceSpinState",
    "Vertex",
    "VertexLoweringRule",
]
