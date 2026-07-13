from __future__ import annotations

import math
import re
from dataclasses import dataclass, replace
from typing import Mapping, Sequence

from symbolica import E, S, Expression, Replacement
from symbolica.community.spenso import (
    LibraryTensor,
    Representation,
    TensorLibrary,
    TensorName,
    TensorNetwork,
)

from .ufo_tensors import (
    classify_trilinear_color_expression,
    normalize_color_expression,
    normalize_lorentz_expression,
)


@dataclass(frozen=True)
class CompiledCouplingOrder:
    name: str
    expansion_order: int
    hierarchy: int

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "expansion_order": self.expansion_order,
            "hierarchy": self.hierarchy,
        }


@dataclass(frozen=True)
class CompiledParameterRecord:
    name: str
    nature: str
    parameter_type: str
    value: tuple[float, float] | None
    expression: str | None
    resolved_expression: str
    lhablock: str | None
    lhacode: tuple[int, ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "nature": self.nature,
            "parameter_type": self.parameter_type,
            "value": None if self.value is None else list(self.value),
            "expression": self.expression,
            "resolved_expression": self.resolved_expression,
            "lhablock": self.lhablock,
            "lhacode": list(self.lhacode),
        }


@dataclass(frozen=True)
class CompiledParticleRecord:
    name: str
    antiname: str
    pdg_code: int
    spin: int
    color: int
    mass: str
    width: str
    charge: float
    ghost_number: int
    propagating: bool
    goldstoneboson: bool
    propagator: str | None
    component_dimension: int | None = None
    auxiliary_kind: str | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "antiname": self.antiname,
            "pdg_code": self.pdg_code,
            "spin": self.spin,
            "color": self.color,
            "mass": self.mass,
            "width": self.width,
            "charge": self.charge,
            "ghost_number": self.ghost_number,
            "propagating": self.propagating,
            "goldstoneboson": self.goldstoneboson,
            "propagator": self.propagator,
            "component_dimension": self.component_dimension,
            "auxiliary_kind": self.auxiliary_kind,
        }


@dataclass(frozen=True)
class CompiledCouplingRecord:
    name: str
    expression: str
    resolved_expression: str
    value: tuple[float, float]
    orders: tuple[tuple[str, int], ...]

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "expression": self.expression,
            "resolved_expression": self.resolved_expression,
            "value": list(self.value),
            "orders": [[name, value] for name, value in self.orders],
        }


@dataclass(frozen=True)
class CompiledPropagatorRecord:
    name: str
    particle: str
    numerator: str
    denominator: str
    custom: bool

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "particle": self.particle,
            "numerator": self.numerator,
            "denominator": self.denominator,
            "custom": self.custom,
        }


@dataclass(frozen=True)
class CompiledVertexTerm:
    id: int
    vertex: str
    particles: tuple[str, ...]
    color_index: int
    lorentz_index: int
    color_source: str
    color_expression: str
    lorentz_name: str
    lorentz_source: str
    lorentz_expression: str
    coupling: str
    coupling_expression: str
    coupling_orders: tuple[tuple[str, int], ...]
    backend: str = "ufo"
    lc_color_normalization_power: int = 0

    @property
    def valence(self) -> int:
        return len(self.particles)

    def to_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "vertex": self.vertex,
            "particles": list(self.particles),
            "valence": self.valence,
            "color_index": self.color_index,
            "lorentz_index": self.lorentz_index,
            "color_source": self.color_source,
            "color_expression": self.color_expression,
            "lorentz_name": self.lorentz_name,
            "lorentz_source": self.lorentz_source,
            "lorentz_expression": self.lorentz_expression,
            "coupling": self.coupling,
            "coupling_expression": self.coupling_expression,
            "coupling_orders": [
                [name, value] for name, value in self.coupling_orders
            ],
            "backend": self.backend,
            "lc_color_normalization_power": self.lc_color_normalization_power,
        }


@dataclass(frozen=True)
class CompiledOrientedKernel:
    kind: int
    term_id: int
    vertex: str
    particles: tuple[str, str, str]
    source_particle_legs: tuple[int, int, int]
    component_expressions: tuple[str, ...]
    coupling_expression: str
    coupling_orders: tuple[tuple[str, int], ...]
    runtime_parameters: tuple[str, ...]
    color_source: str
    color_expression: str
    color_projection_structure: str | None = None
    color_projection_coefficient: tuple[float, float] | None = None
    lc_color_normalization_power: int = 0
    term_ids: tuple[int, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "term_id": self.term_id,
            "vertex": self.vertex,
            "particles": list(self.particles),
            "source_particle_legs": list(self.source_particle_legs),
            "component_expressions": list(self.component_expressions),
            "coupling_expression": self.coupling_expression,
            "coupling_orders": [
                [name, value] for name, value in self.coupling_orders
            ],
            "runtime_parameters": list(self.runtime_parameters),
            "color_source": self.color_source,
            "color_expression": self.color_expression,
            "color_projection_structure": self.color_projection_structure,
            "color_projection_coefficient": (
                None
                if self.color_projection_coefficient is None
                else list(self.color_projection_coefficient)
            ),
            "lc_color_normalization_power": self.lc_color_normalization_power,
            "term_ids": list(self.term_ids or (self.term_id,)),
        }


@dataclass(frozen=True)
class _ContactTreeNode:
    legs: tuple[int, ...]
    particle: CompiledParticleRecord
    physical_leg: int | None = None
    left: "_ContactTreeNode | None" = None
    right: "_ContactTreeNode | None" = None

    @property
    def is_leaf(self) -> bool:
        return self.physical_leg is not None


@dataclass(frozen=True)
class CompiledModelIR:
    name: str
    orders: tuple[CompiledCouplingOrder, ...]
    parameters: tuple[CompiledParameterRecord, ...]
    particles: tuple[CompiledParticleRecord, ...]
    couplings: tuple[CompiledCouplingRecord, ...]
    propagators: tuple[CompiledPropagatorRecord, ...]
    vertex_terms: tuple[CompiledVertexTerm, ...]
    oriented_kernels: tuple[CompiledOrientedKernel, ...]

    @property
    def max_vertex_valence(self) -> int:
        return max((term.valence for term in self.vertex_terms), default=0)

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "orders": [item.to_dict() for item in self.orders],
            "parameters": [item.to_dict() for item in self.parameters],
            "particles": [item.to_dict() for item in self.particles],
            "couplings": [item.to_dict() for item in self.couplings],
            "propagators": [item.to_dict() for item in self.propagators],
            "vertex_terms": [item.to_dict() for item in self.vertex_terms],
            "oriented_kernels": [item.to_dict() for item in self.oriented_kernels],
            "max_vertex_valence": self.max_vertex_valence,
        }

    @staticmethod
    def from_dict(payload: Mapping[str, object]) -> "CompiledModelIR":
        return CompiledModelIR(
            name=str(payload["name"]),
            orders=tuple(
                CompiledCouplingOrder(
                    name=str(item["name"]),
                    expansion_order=int(item["expansion_order"]),
                    hierarchy=int(item["hierarchy"]),
                )
                for item in _mappings(payload.get("orders"))
            ),
            parameters=tuple(
                CompiledParameterRecord(
                    name=str(item["name"]),
                    nature=str(item["nature"]),
                    parameter_type=str(item["parameter_type"]),
                    value=_optional_pair(item.get("value")),
                    expression=_optional_string(item.get("expression")),
                    resolved_expression=str(item["resolved_expression"]),
                    lhablock=_optional_string(item.get("lhablock")),
                    lhacode=tuple(int(value) for value in _sequence(item.get("lhacode"))),
                )
                for item in _mappings(payload.get("parameters"))
            ),
            particles=tuple(
                CompiledParticleRecord(
                    name=str(item["name"]),
                    antiname=str(item["antiname"]),
                    pdg_code=int(item["pdg_code"]),
                    spin=int(item["spin"]),
                    color=int(item["color"]),
                    mass=str(item["mass"]),
                    width=str(item["width"]),
                    charge=float(item["charge"]),
                    ghost_number=int(item["ghost_number"]),
                    propagating=bool(item["propagating"]),
                    goldstoneboson=bool(item["goldstoneboson"]),
                    propagator=_optional_string(item.get("propagator")),
                    component_dimension=(
                        None
                        if item.get("component_dimension") is None
                        else int(item["component_dimension"])
                    ),
                    auxiliary_kind=_optional_string(item.get("auxiliary_kind")),
                )
                for item in _mappings(payload.get("particles"))
            ),
            couplings=tuple(
                CompiledCouplingRecord(
                    name=str(item["name"]),
                    expression=str(item["expression"]),
                    resolved_expression=str(item["resolved_expression"]),
                    value=_pair(item.get("value")),
                    orders=_orders(item.get("orders")),
                )
                for item in _mappings(payload.get("couplings"))
            ),
            propagators=tuple(
                CompiledPropagatorRecord(
                    name=str(item["name"]),
                    particle=str(item["particle"]),
                    numerator=str(item["numerator"]),
                    denominator=str(item["denominator"]),
                    custom=bool(item["custom"]),
                )
                for item in _mappings(payload.get("propagators"))
            ),
            vertex_terms=tuple(
                CompiledVertexTerm(
                    id=int(item["id"]),
                    vertex=str(item["vertex"]),
                    particles=tuple(str(value) for value in _sequence(item["particles"])),
                    color_index=int(item["color_index"]),
                    lorentz_index=int(item["lorentz_index"]),
                    color_source=str(item["color_source"]),
                    color_expression=str(item["color_expression"]),
                    lorentz_name=str(item["lorentz_name"]),
                    lorentz_source=str(item["lorentz_source"]),
                    lorentz_expression=str(item["lorentz_expression"]),
                    coupling=str(item["coupling"]),
                    coupling_expression=str(item["coupling_expression"]),
                    coupling_orders=_orders(item.get("coupling_orders")),
                    backend=str(item.get("backend", "ufo")),
                    lc_color_normalization_power=int(
                        item.get("lc_color_normalization_power", 0)
                    ),
                )
                for item in _mappings(payload.get("vertex_terms"))
            ),
            oriented_kernels=tuple(
                CompiledOrientedKernel(
                    kind=int(item["kind"]),
                    term_id=int(item["term_id"]),
                    vertex=str(item["vertex"]),
                    particles=cast_tuple3(item["particles"]),
                    source_particle_legs=cast_int_tuple3(item["source_particle_legs"]),
                    component_expressions=tuple(
                        str(value) for value in _sequence(item["component_expressions"])
                    ),
                    coupling_expression=str(item["coupling_expression"]),
                    coupling_orders=_orders(item.get("coupling_orders")),
                    runtime_parameters=tuple(
                        str(value) for value in _sequence(item["runtime_parameters"])
                    ),
                    color_source=str(
                        item.get("color_source", item["color_expression"])
                    ),
                    color_expression=str(item["color_expression"]),
                    color_projection_structure=_optional_string(
                        item.get("color_projection_structure")
                    ),
                    color_projection_coefficient=(
                        None
                        if item.get("color_projection_coefficient") is None
                        else _pair(item.get("color_projection_coefficient"))
                    ),
                    lc_color_normalization_power=int(
                        item.get("lc_color_normalization_power", 0)
                    ),
                    term_ids=tuple(
                        int(value)
                        for value in _sequence(
                            item.get("term_ids", [int(item["term_id"])])
                        )
                    ),
                )
                for item in _mappings(payload.get("oriented_kernels"))
            ),
        )


def compile_ufo_model_ir(model: Mapping[str, object]) -> CompiledModelIR:
    particles = tuple(_particle(item) for item in _mappings(model.get("particles")))
    particle_by_name = {particle.name: particle for particle in particles}
    lorentz_by_name = {
        str(item["name"]): item
        for item in _mappings(model.get("lorentz_structures"))
    }
    parameter_records = tuple(
        _parameter(item) for item in _mappings(model.get("parameters"))
    )
    parameter_records = _resolve_parameter_records(parameter_records)
    coupling_records = tuple(
        _coupling(item) for item in _mappings(model.get("couplings"))
    )
    coupling_records = _resolve_coupling_records(
        coupling_records,
        parameter_records,
    )
    coupling_by_name = {coupling.name: coupling for coupling in coupling_records}
    terms: list[CompiledVertexTerm] = []
    for vertex in _mappings(model.get("vertex_rules")):
        particle_names = tuple(str(value) for value in _sequence(vertex["particles"]))
        try:
            vertex_particles = tuple(particle_by_name[name] for name in particle_names)
        except KeyError as exc:
            raise ValueError(
                f"vertex {vertex.get('name')} refers to unknown particle {exc.args[0]}"
            ) from exc
        colors = tuple(particle.color for particle in vertex_particles)
        color_sources = tuple(
            str(value) for value in _sequence(vertex["color_structures"])
        )
        lorentz_names = tuple(
            str(value) for value in _sequence(vertex["lorentz_structures"])
        )
        coupling_matrix = _sequence(vertex["couplings"])
        if len(coupling_matrix) != len(color_sources):
            raise ValueError(
                f"vertex {vertex.get('name')} coupling rows do not match color structures"
            )
        normalized_colors = tuple(
            normalize_color_expression(source, colors) for source in color_sources
        )
        normalized_lorentz = []
        for name in lorentz_names:
            try:
                lorentz = lorentz_by_name[name]
            except KeyError as exc:
                raise ValueError(
                    f"vertex {vertex.get('name')} refers to unknown Lorentz structure {name}"
                ) from exc
            normalized_lorentz.append(
                normalize_lorentz_expression(
                    str(lorentz["structure"]),
                    tuple(int(value) for value in _sequence(lorentz["spins"])),
                )
            )
        for color_index, row_value in enumerate(coupling_matrix):
            row = _sequence(row_value)
            if len(row) != len(lorentz_names):
                raise ValueError(
                    f"vertex {vertex.get('name')} coupling columns do not match "
                    "Lorentz structures"
                )
            for lorentz_index, coupling_value in enumerate(row):
                if coupling_value is None:
                    continue
                coupling_name = str(coupling_value)
                try:
                    coupling = coupling_by_name[coupling_name]
                except KeyError as exc:
                    raise ValueError(
                        f"vertex {vertex.get('name')} refers to unknown coupling "
                        f"{coupling_name}"
                    ) from exc
                source_lorentz = lorentz_by_name[lorentz_names[lorentz_index]]
                terms.append(
                    CompiledVertexTerm(
                        id=len(terms),
                        vertex=str(vertex["name"]),
                        particles=particle_names,
                        color_index=color_index,
                        lorentz_index=lorentz_index,
                        color_source=color_sources[color_index],
                        color_expression=normalized_colors[color_index].expression,
                        lorentz_name=lorentz_names[lorentz_index],
                        lorentz_source=str(source_lorentz["structure"]),
                        lorentz_expression=normalized_lorentz[lorentz_index].expression,
                        coupling=coupling.name,
                        coupling_expression=coupling.resolved_expression,
                        coupling_orders=coupling.orders,
                        lc_color_normalization_power=(
                            _lc_color_normalization_power(
                                color_sources[color_index]
                            )
                        ),
                    )
                )
    propagators = tuple(
        _propagator(item, particles) for item in _mappings(model.get("propagators"))
    )
    oriented_kernels = _compile_oriented_kernels(
        terms,
        particles,
        parameter_records,
    )
    contact_particles, contact_kernels = _compile_four_point_contact_kernels(
        terms,
        particles,
        start_kind=len(oriented_kernels),
    )
    contact_particles, contact_kernels = _deduplicate_contact_partials(
        contact_particles,
        contact_kernels,
        terms,
    )
    contact_kernels = _fuse_contact_finals(contact_kernels, terms)
    particles = (*particles, *contact_particles)
    oriented_kernels = (*oriented_kernels, *contact_kernels)
    tree_start_kind = (
        max((kernel.kind for kernel in oriented_kernels), default=-1) + 1
    )
    tree_particles, tree_kernels = _compile_color_singlet_contact_trees(
        terms,
        particles,
        start_kind=tree_start_kind,
    )
    particles = (*particles, *tree_particles)
    oriented_kernels = (*oriented_kernels, *tree_kernels)
    oriented_kernels = _annotate_oriented_kernel_color_projections(
        oriented_kernels,
        particles,
    )
    return CompiledModelIR(
        name=str(model.get("name", "unnamed-model")),
        orders=tuple(_order(item) for item in _mappings(model.get("orders"))),
        parameters=parameter_records,
        particles=particles,
        couplings=coupling_records,
        propagators=propagators,
        vertex_terms=tuple(terms),
        oriented_kernels=oriented_kernels,
    )


def _annotate_oriented_kernel_color_projections(
    kernels: Sequence[CompiledOrientedKernel],
    particles: Sequence[CompiledParticleRecord],
) -> tuple[CompiledOrientedKernel, ...]:
    particle_by_name = {particle.name: particle for particle in particles}
    annotated: list[CompiledOrientedKernel] = []
    for kernel in kernels:
        representations = tuple(
            particle_by_name[name].color for name in kernel.particles
        )
        structure, coefficient = classify_trilinear_color_expression(
            kernel.color_expression,
            kernel.color_source,
            representations,
        )
        annotated.append(
            replace(
                kernel,
                color_projection_structure=structure,
                color_projection_coefficient=(
                    float(coefficient.real),
                    float(coefficient.imag),
                ),
            )
        )
    return tuple(annotated)


def compile_builtin_model_ir(model: Mapping[str, object]) -> CompiledModelIR:
    particles = tuple(_particle(item) for item in _mappings(model.get("particles")))
    terms = tuple(
        CompiledVertexTerm(
            id=index,
            vertex=str(vertex["name"]),
            particles=tuple(str(value) for value in _sequence(vertex["particles"])),
            color_index=0,
            lorentz_index=0,
            color_source="built-in",
            color_expression="built-in",
            lorentz_name=str(vertex["builtin_kind"]),
            lorentz_source="built-in",
            lorentz_expression="built-in",
            coupling=f"builtin_coupling_{index}",
            coupling_expression=str(vertex.get("builtin_coupling", [1.0, 0.0])),
            coupling_orders=(),
            backend="built-in",
        )
        for index, vertex in enumerate(_mappings(model.get("vertex_rules")))
    )
    return CompiledModelIR(
        name=str(model.get("name", "built-in-sm")),
        orders=tuple(_order(item) for item in _mappings(model.get("orders"))),
        parameters=tuple(
            _parameter(item) for item in _mappings(model.get("parameters"))
        ),
        particles=particles,
        couplings=(),
        propagators=(),
        vertex_terms=terms,
        oriented_kernels=(),
    )


def _compile_oriented_kernels(
    terms: Sequence[CompiledVertexTerm],
    particles: Sequence[CompiledParticleRecord],
    parameters: Sequence[CompiledParameterRecord],
) -> tuple[CompiledOrientedKernel, ...]:
    particle_by_name = {particle.name: particle for particle in particles}
    parameter_by_name = {parameter.name: parameter for parameter in parameters}
    external_parameters = {
        parameter.name for parameter in parameters if parameter.nature == "external"
    }
    kernels: list[CompiledOrientedKernel] = []
    for term in terms:
        if term.valence != 3:
            continue
        oriented_result_particles: set[str] = set()
        for result_leg in range(3):
            result_source_name = term.particles[result_leg]
            if result_source_name in oriented_result_particles:
                continue
            oriented_result_particles.add(result_source_name)
            input_legs = tuple(leg for leg in range(3) if leg != result_leg)
            input_orders = (
                (input_legs,)
                if term.particles[input_legs[0]] == term.particles[input_legs[1]]
                else (input_legs, tuple(reversed(input_legs)))
            )
            for left_leg, right_leg in input_orders:
                result_source = particle_by_name[term.particles[result_leg]]
                try:
                    result_name = particle_by_name[result_source.antiname].name
                except KeyError as exc:
                    raise ValueError(
                        f"vertex {term.vertex} particle {result_source.name} refers to "
                        f"absent antiparticle {result_source.antiname}"
                    ) from exc
                components = _oriented_component_expressions(
                    term,
                    particle_by_name,
                    left_leg=left_leg,
                    right_leg=right_leg,
                    result_leg=result_leg,
                    kind=len(kernels),
                    use_transverse_massless_yang_mills=(
                        _is_single_structure_constant(term.color_expression)
                        and all(
                            _is_compile_time_zero_parameter(
                                particle_by_name[name].mass,
                                parameter_by_name,
                            )
                            for name in term.particles
                        )
                    ),
                )
                coupling_symbols = set(E(term.coupling_expression).get_all_symbols(False))
                runtime_parameters = tuple(
                    sorted(
                        name
                        for name in external_parameters
                        if S(f"UFO::{name}") in coupling_symbols
                    )
                )
                kernels.append(
                    CompiledOrientedKernel(
                        kind=len(kernels),
                        term_id=term.id,
                        vertex=term.vertex,
                        particles=(
                            term.particles[left_leg],
                            term.particles[right_leg],
                            result_name,
                        ),
                        source_particle_legs=(left_leg, right_leg, result_leg),
                        component_expressions=components,
                        coupling_expression=term.coupling_expression,
                        coupling_orders=term.coupling_orders,
                        runtime_parameters=runtime_parameters,
                        color_source=term.color_source,
                        color_expression=term.color_expression,
                        lc_color_normalization_power=(
                            term.lc_color_normalization_power
                        ),
                        term_ids=(term.id,),
                    )
                )
    return _fuse_oriented_kernels(kernels)


def _compile_four_point_contact_kernels(
    terms: Sequence[CompiledVertexTerm],
    particles: Sequence[CompiledParticleRecord],
    *,
    start_kind: int,
) -> tuple[tuple[CompiledParticleRecord, ...], tuple[CompiledOrientedKernel, ...]]:
    """Lower momentum-independent four-point tensors through dense auxiliaries."""

    particle_by_name = {particle.name: particle for particle in particles}
    used_pdgs = {abs(particle.pdg_code) for particle in particles}
    next_pdg = max(9_000_000, max(used_pdgs, default=0) + 1)
    auxiliary_particles: list[CompiledParticleRecord] = []
    kernels: list[CompiledOrientedKernel] = []

    def allocate_pdg() -> int:
        nonlocal next_pdg
        while next_pdg in used_pdgs:
            next_pdg += 1
        result = next_pdg
        used_pdgs.add(result)
        next_pdg += 1
        return result

    for term in terms:
        if term.valence != 4:
            continue
        if "ufo_momentum_" in term.lorentz_expression:
            continue
        source_particles = tuple(particle_by_name[name] for name in term.particles)
        oriented_result_particles: set[str] = set()
        for result_leg in range(4):
            source_result = source_particles[result_leg]
            if source_result.name in oriented_result_particles:
                continue
            oriented_result_particles.add(source_result.name)
            contact_split = _four_point_contact_color_split(term, result_leg)
            pair_legs = contact_split[0]
            remaining_leg = contact_split[1]
            outer_color_source = contact_split[2]
            final_color_source = contact_split[3]
            outer_color_power = contact_split[4]
            final_color_power = contact_split[5]
            outer_color_factor = contact_split[6]
            final_color_factor = contact_split[7]
            color_dummy = contact_split[8]
            assignment_multiplicity = sum(
                source_particles[leg].name
                == source_particles[remaining_leg].name
                for leg in range(4)
                if leg != result_leg
            )
            open_legs = tuple(sorted((remaining_leg, result_leg)))
            auxiliary_name = f"__pyamplicol_contact_{term.id}_r{result_leg}"
            canonical_kind = start_kind + len(kernels)
            canonical_components = _contact_partial_component_expressions(
                term,
                particle_by_name,
                left_leg=pair_legs[0],
                right_leg=pair_legs[1],
                open_legs=open_legs,
                kind=canonical_kind,
            )
            representative_indices, component_expansion = _compress_contact_components(
                canonical_components
            )
            auxiliary_dimension = len(representative_indices)
            auxiliary = CompiledParticleRecord(
                name=auxiliary_name,
                antiname=auxiliary_name,
                pdg_code=allocate_pdg(),
                spin=-1,
                color=_contact_auxiliary_color(
                    term,
                    source_particles,
                    remaining_leg=remaining_leg,
                    result_leg=result_leg,
                ),
                mass="ZERO",
                width="ZERO",
                charge=0.0,
                ghost_number=0,
                propagating=False,
                goldstoneboson=False,
                propagator=None,
                component_dimension=auxiliary_dimension,
                auxiliary_kind=(
                    f"ufo-contact:{term.id}:result-{result_leg}:"
                    + ",".join(str(leg) for leg in open_legs)
                ),
            )
            auxiliary_particles.append(auxiliary)

            pair_orders = (
                (pair_legs,)
                if source_particles[pair_legs[0]].name
                == source_particles[pair_legs[1]].name
                else (pair_legs, tuple(reversed(pair_legs)))
            )
            canonical_outer_sign = (
                _permutation_sign(
                    outer_color_factor,
                    (color_dummy, pair_legs[0] + 1, pair_legs[1] + 1),
                )
                if outer_color_factor
                else 1
            )
            for left_leg, right_leg in pair_orders:
                kind = start_kind + len(kernels)
                components = _contact_partial_component_expressions(
                    term,
                    particle_by_name,
                    left_leg=left_leg,
                    right_leg=right_leg,
                    open_legs=open_legs,
                    kind=kind,
                )
                outer_sign = (
                    _permutation_sign(
                        outer_color_factor,
                        (color_dummy, left_leg + 1, right_leg + 1),
                    )
                    if outer_color_factor
                    else 1
                )
                kernels.append(
                    CompiledOrientedKernel(
                        kind=kind,
                        term_id=term.id,
                        vertex=f"{term.vertex}::contact-partial",
                        particles=(
                            source_particles[left_leg].name,
                            source_particles[right_leg].name,
                            auxiliary.name,
                        ),
                        source_particle_legs=(left_leg, right_leg, -1),
                        component_expressions=tuple(
                            _canonicalize_oriented_kernel_component(
                                E(components[index])
                                * outer_sign
                                / canonical_outer_sign
                            ).to_canonical_string()
                            for index in representative_indices
                        ),
                        coupling_expression="1",
                        coupling_orders=(),
                        runtime_parameters=(),
                        color_source="1",
                        color_expression="1",
                        lc_color_normalization_power=0,
                        term_ids=(),
                    )
                )

            result_name = particle_by_name[source_result.antiname].name
            final_orders = (
                ((True, auxiliary.name, source_particles[remaining_leg].name),)
                if auxiliary.name == source_particles[remaining_leg].name
                else (
                    (True, auxiliary.name, source_particles[remaining_leg].name),
                    (False, source_particles[remaining_leg].name, auxiliary.name),
                )
            )
            for auxiliary_on_left, left_name, right_name in final_orders:
                kind = start_kind + len(kernels)
                final_input_tokens = (
                    (color_dummy, remaining_leg + 1)
                    if auxiliary_on_left
                    else (remaining_leg + 1, color_dummy)
                )
                final_sign = (
                    _permutation_sign(
                        final_color_factor,
                        (result_leg + 1, *final_input_tokens),
                    )
                    if final_color_factor
                    else 1
                )
                derived_coupling = S(f"UFO::derived_coupling_{term.id}")
                final_prefactor = (
                    canonical_outer_sign
                    * final_sign
                    * derived_coupling
                    / assignment_multiplicity
                )
                combined_color_source = (
                    "UFO::{}::f(1,2,3)*UFO::{}::f(1,2,3)"
                    if outer_color_factor
                    else term.color_source
                )
                kernels.append(
                    CompiledOrientedKernel(
                        kind=kind,
                        term_id=term.id,
                        vertex=f"{term.vertex}::contact-final",
                        particles=(left_name, right_name, result_name),
                        source_particle_legs=(
                            -1 if auxiliary_on_left else remaining_leg,
                            remaining_leg if auxiliary_on_left else -1,
                            result_leg,
                        ),
                        component_expressions=tuple(
                            (
                                E(component) * final_prefactor
                            ).to_canonical_string()
                            for component in _contact_final_component_expressions(
                                source_particles,
                                auxiliary,
                                open_legs=open_legs,
                                remaining_leg=remaining_leg,
                                result_leg=result_leg,
                                kind=kind,
                                auxiliary_on_left=auxiliary_on_left,
                                component_expansion=component_expansion,
                            )
                        ),
                        coupling_expression="1",
                        coupling_orders=term.coupling_orders,
                        runtime_parameters=(f"derived_coupling_{term.id}",),
                        color_source=combined_color_source,
                        color_expression=combined_color_source,
                        lc_color_normalization_power=(
                            outer_color_power + final_color_power
                        ),
                        term_ids=(term.id,),
                    )
                )
    return tuple(auxiliary_particles), tuple(kernels)


def _compile_color_singlet_contact_trees(
    terms: Sequence[CompiledVertexTerm],
    particles: Sequence[CompiledParticleRecord],
    *,
    start_kind: int,
) -> tuple[tuple[CompiledParticleRecord, ...], tuple[CompiledOrientedKernel, ...]]:
    """Lower arbitrary color-singlet contacts to balanced trivalent trees."""

    particle_by_name = {particle.name: particle for particle in particles}
    used_pdgs = {abs(particle.pdg_code) for particle in particles}
    next_pdg = max(9_000_000, max(used_pdgs, default=0) + 1)
    auxiliary_particles: list[CompiledParticleRecord] = []
    kernels: list[CompiledOrientedKernel] = []
    final_template_cache: dict[
        tuple[object, ...],
        tuple[int, tuple[Expression, ...]],
    ] = {}

    def allocate_pdg() -> int:
        nonlocal next_pdg
        while next_pdg in used_pdgs:
            next_pdg += 1
        result = next_pdg
        used_pdgs.add(result)
        next_pdg += 1
        return result

    for term in terms:
        if term.valence < 4 or not _contact_term_is_color_singlet(term):
            continue
        if term.valence == 4 and "ufo_momentum_" not in term.lorentz_expression:
            continue
        source_particles = tuple(particle_by_name[name] for name in term.particles)
        scalar_product_tree = all(particle.spin == 1 for particle in source_particles)
        oriented_result_particles: set[str] = set()
        for result_leg in range(term.valence):
            source_result = source_particles[result_leg]
            if source_result.name in oriented_result_particles:
                continue
            oriented_result_particles.add(source_result.name)
            input_legs = tuple(
                leg for leg in range(term.valence) if leg != result_leg
            )

            def build_node(legs: tuple[int, ...]) -> _ContactTreeNode:
                if len(legs) == 1:
                    leg = legs[0]
                    return _ContactTreeNode(
                        legs=legs,
                        particle=source_particles[leg],
                        physical_leg=leg,
                    )
                split = len(legs) // 2
                left_node = build_node(legs[:split])
                right_node = build_node(legs[split:])
                auxiliary_name = (
                    f"__pyamplicol_contact_tree_{term.id}_r{result_leg}_"
                    + "_".join(str(leg) for leg in legs)
                )
                auxiliary_dimension = (
                    1
                    if scalar_product_tree
                    else sum(
                        _spin_dimension(source_particles[leg].spin) + 4
                        for leg in legs
                    )
                )
                auxiliary = CompiledParticleRecord(
                    name=auxiliary_name,
                    antiname=auxiliary_name,
                    pdg_code=allocate_pdg(),
                    spin=-1,
                    color=1,
                    mass="ZERO",
                    width="ZERO",
                    charge=0.0,
                    ghost_number=0,
                    propagating=False,
                    goldstoneboson=False,
                    propagator=None,
                    component_dimension=auxiliary_dimension,
                    auxiliary_kind=(
                        f"ufo-contact-tree:{term.id}:result-{result_leg}:"
                        + ",".join(str(leg) for leg in legs)
                    ),
                )
                auxiliary_particles.append(auxiliary)
                node = _ContactTreeNode(
                    legs=legs,
                    particle=auxiliary,
                    left=left_node,
                    right=right_node,
                )
                _append_contact_tree_partial_kernels(
                    kernels,
                    term=term,
                    node=node,
                    source_particles=source_particles,
                    scalar_product_tree=scalar_product_tree,
                    start_kind=start_kind,
                )
                return node

            split = len(input_legs) // 2
            left_root = build_node(input_legs[:split])
            right_root = build_node(input_legs[split:])
            assignment_multiplicity = _contact_tree_assignment_multiplicity(
                input_legs,
                source_particles,
                left_root,
                right_root,
            )
            _append_contact_tree_final_kernels(
                kernels,
                term=term,
                source_particles=source_particles,
                left_node=left_root,
                right_node=right_root,
                result_leg=result_leg,
                result_name=particle_by_name[source_result.antiname].name,
                scalar_product_tree=scalar_product_tree,
                assignment_multiplicity=assignment_multiplicity,
                start_kind=start_kind,
                template_cache=final_template_cache,
            )
    return tuple(auxiliary_particles), tuple(kernels)


def _append_contact_tree_partial_kernels(
    kernels: list[CompiledOrientedKernel],
    *,
    term: CompiledVertexTerm,
    node: _ContactTreeNode,
    source_particles: Sequence[CompiledParticleRecord],
    scalar_product_tree: bool,
    start_kind: int,
) -> None:
    if node.left is None or node.right is None:
        raise ValueError("contact tree partial node is missing a child")
    orientations = (
        ((node.left, node.right),)
        if node.left.particle.name == node.right.particle.name
        else ((node.left, node.right), (node.right, node.left))
    )
    for actual_left, actual_right in orientations:
        kind = start_kind + len(kernels)
        canonical_left_side = "left" if actual_left is node.left else "right"
        canonical_right_side = "right" if actual_right is node.right else "left"
        left_payload = _contact_tree_node_payload(
            kind,
            canonical_left_side,
            node.left,
            source_particles,
            scalar_product_tree=scalar_product_tree,
        )
        right_payload = _contact_tree_node_payload(
            kind,
            canonical_right_side,
            node.right,
            source_particles,
            scalar_product_tree=scalar_product_tree,
        )
        components = (
            (left_payload[0] * right_payload[0],)
            if scalar_product_tree
            else (*left_payload, *right_payload)
        )
        kernels.append(
            CompiledOrientedKernel(
                kind=kind,
                term_id=term.id,
                vertex=f"{term.vertex}::contact-tree-partial",
                particles=(
                    actual_left.particle.name,
                    actual_right.particle.name,
                    node.particle.name,
                ),
                source_particle_legs=(
                    _contact_tree_source_leg(actual_left),
                    _contact_tree_source_leg(actual_right),
                    -1,
                ),
                component_expressions=tuple(
                    component.to_canonical_string() for component in components
                ),
                coupling_expression="1",
                coupling_orders=(),
                runtime_parameters=(),
                color_source="1",
                color_expression="1",
                lc_color_normalization_power=0,
                term_ids=(),
            )
        )


def _append_contact_tree_final_kernels(
    kernels: list[CompiledOrientedKernel],
    *,
    term: CompiledVertexTerm,
    source_particles: Sequence[CompiledParticleRecord],
    left_node: _ContactTreeNode,
    right_node: _ContactTreeNode,
    result_leg: int,
    result_name: str,
    scalar_product_tree: bool,
    assignment_multiplicity: int,
    start_kind: int,
    template_cache: dict[
        tuple[object, ...],
        tuple[int, tuple[Expression, ...]],
    ],
) -> None:
    orientations = (
        ((left_node, right_node),)
        if left_node.particle.name == right_node.particle.name
        else ((left_node, right_node), (right_node, left_node))
    )
    canonical_kind = start_kind + len(kernels)
    if scalar_product_tree:
        left_payload = _contact_tree_node_payload(
            canonical_kind,
            "left",
            left_node,
            source_particles,
            scalar_product_tree=True,
        )
        right_payload = _contact_tree_node_payload(
            canonical_kind,
            "right",
            right_node,
            source_particles,
            scalar_product_tree=True,
        )
        canonical_components = (
            left_payload[0]
            * right_payload[0]
            * E(term.lorentz_expression),
        )
    else:
        template_key = (
            term.lorentz_expression,
            tuple(particle.spin for particle in source_particles),
            result_leg,
            left_node.legs,
            right_node.legs,
        )
        cached = template_cache.get(template_key)
        if cached is None:
            canonical_components = _contact_tree_final_component_expressions(
                term,
                source_particles,
                left_node=left_node,
                right_node=right_node,
                result_leg=result_leg,
                kind=canonical_kind,
                canonical_left_side="left",
                canonical_right_side="right",
            )
            template_cache[template_key] = canonical_kind, canonical_components
        else:
            template_kind, template_components = cached
            canonical_components = tuple(
                _remap_kernel_symbols(
                    component,
                    old_kind=template_kind,
                    new_kind=canonical_kind,
                )
                for component in template_components
            )
    prefactor = S(f"UFO::derived_coupling_{term.id}") / assignment_multiplicity
    canonical_weighted_components = tuple(
        _canonicalize_oriented_kernel_component(component * prefactor)
        for component in canonical_components
    )
    for actual_left, actual_right in orientations:
        kind = start_kind + len(kernels)
        components = tuple(
            _remap_kernel_symbols(
                component,
                old_kind=canonical_kind,
                new_kind=kind,
                swap_sides=actual_left is right_node,
            )
            for component in canonical_weighted_components
        )
        kernels.append(
            CompiledOrientedKernel(
                kind=kind,
                term_id=term.id,
                vertex=f"{term.vertex}::contact-tree-final",
                particles=(
                    actual_left.particle.name,
                    actual_right.particle.name,
                    result_name,
                ),
                source_particle_legs=(
                    _contact_tree_source_leg(actual_left),
                    _contact_tree_source_leg(actual_right),
                    result_leg,
                ),
                component_expressions=tuple(
                    component.to_canonical_string() for component in components
                ),
                coupling_expression="1",
                coupling_orders=term.coupling_orders,
                runtime_parameters=(f"derived_coupling_{term.id}",),
                color_source=term.color_source,
                color_expression=term.color_expression,
                lc_color_normalization_power=term.lc_color_normalization_power,
                term_ids=(term.id,),
            )
        )


def _contact_tree_node_payload(
    kind: int,
    side: str,
    node: _ContactTreeNode,
    source_particles: Sequence[CompiledParticleRecord],
    *,
    scalar_product_tree: bool,
) -> tuple[Expression, ...]:
    if node.is_leaf:
        if node.physical_leg is None:
            raise ValueError("contact tree leaf has no source leg")
        dimension = _spin_dimension(source_particles[node.physical_leg].spin)
        components = tuple(
            S(f"pyamplicol::kernel_{kind}_{side}_{index}")
            for index in range(dimension)
        )
        if scalar_product_tree:
            return components
        momenta = tuple(
            S(f"pyamplicol::kernel_{kind}_{side}_momentum_{index}")
            for index in range(4)
        )
        return (*components, *momenta)
    dimension = node.particle.component_dimension
    if dimension is None:
        raise ValueError("contact tree auxiliary has no component dimension")
    return tuple(
        S(f"pyamplicol::kernel_{kind}_{side}_{index}")
        for index in range(dimension)
    )


def _contact_tree_final_component_expressions(
    term: CompiledVertexTerm,
    source_particles: Sequence[CompiledParticleRecord],
    *,
    left_node: _ContactTreeNode,
    right_node: _ContactTreeNode,
    result_leg: int,
    kind: int,
    canonical_left_side: str,
    canonical_right_side: str,
) -> tuple[Expression, ...]:
    payload_by_leg = {
        **_contact_tree_payload_by_leg(
            kind,
            canonical_left_side,
            left_node,
            source_particles,
        ),
        **_contact_tree_payload_by_leg(
            kind,
            canonical_right_side,
            right_node,
            source_particles,
        ),
    }
    momentum_by_leg = {
        leg: momentum for leg, (_components, momentum) in payload_by_leg.items()
    }
    momentum_by_leg[result_leg] = tuple(
        -sum(
            (momentum[component] for momentum in momentum_by_leg.values()),
            E("0"),
        )
        for component in range(4)
    )
    library = TensorLibrary.hep_lib_atom()
    expression = E(term.lorentz_expression)
    for leg, (components, _momentum) in sorted(payload_by_leg.items()):
        expression *= _contact_tree_physical_tensor_expression(
            library,
            kind=kind,
            leg=leg,
            spin=source_particles[leg].spin,
            components=components,
        )
    minkowski = Representation.mink(4)
    for leg, momentum in momentum_by_leg.items():
        library.register(
            LibraryTensor.dense(
                TensorName(f"pyamplicol::ufo_momentum_{leg + 1}")(minkowski),
                momentum,
            )
        )
    result = _execute_dense_tensor(expression, library)
    expected = _spin_dimension(source_particles[result_leg].spin)
    if len(result) != expected:
        raise ValueError(
            f"contact tree final {term.vertex}/{term.id} produced {len(result)} "
            f"components, expected {expected}"
        )
    return tuple(
        _replace_evaluator_constants(_as_expression(result[index]))
        for index in range(len(result))
    )


def _contact_tree_payload_by_leg(
    kind: int,
    side: str,
    node: _ContactTreeNode,
    source_particles: Sequence[CompiledParticleRecord],
) -> dict[int, tuple[tuple[Expression, ...], tuple[Expression, ...]]]:
    payload = _contact_tree_node_payload(
        kind,
        side,
        node,
        source_particles,
        scalar_product_tree=False,
    )
    result: dict[int, tuple[tuple[Expression, ...], tuple[Expression, ...]]] = {}
    cursor = 0
    for leg in node.legs:
        dimension = _spin_dimension(source_particles[leg].spin)
        components = tuple(payload[cursor : cursor + dimension])
        cursor += dimension
        momentum = tuple(payload[cursor : cursor + 4])
        cursor += 4
        result[leg] = components, momentum
    if cursor != len(payload):
        raise ValueError("contact tree payload layout mismatch")
    return result


def _contact_tree_physical_tensor_expression(
    library: TensorLibrary,
    *,
    kind: int,
    leg: int,
    spin: int,
    components: Sequence[Expression],
) -> Expression:
    representations = _spin_representations(spin)
    if not representations:
        if len(components) != 1:
            raise ValueError("scalar contact input must have one component")
        return components[0]
    name = TensorName(f"pyamplicol::contact_{kind}_leg_{leg}")
    library.register(LibraryTensor.dense(name(*representations), components))
    return name(*_spin_slots(spin, leg + 1)).to_expression()


def _contact_tree_assignment_multiplicity(
    input_legs: Sequence[int],
    source_particles: Sequence[CompiledParticleRecord],
    left_root: _ContactTreeNode,
    right_root: _ContactTreeNode,
) -> int:
    species_counts: dict[str, int] = {}
    for leg in input_legs:
        name = source_particles[leg].name
        species_counts[name] = species_counts.get(name, 0) + 1
    permutations = math.prod(math.factorial(count) for count in species_counts.values())
    symmetry_nodes = (
        _contact_tree_same_input_node_count(left_root)
        + _contact_tree_same_input_node_count(right_root)
        + int(left_root.particle.name == right_root.particle.name)
    )
    symmetry_divisor = 2**symmetry_nodes
    if permutations % symmetry_divisor:
        raise ValueError("contact tree assignment symmetry is not integral")
    return permutations // symmetry_divisor


def _contact_tree_same_input_node_count(node: _ContactTreeNode) -> int:
    if node.is_leaf:
        return 0
    if node.left is None or node.right is None:
        raise ValueError("contact tree internal node is missing a child")
    return (
        int(node.left.particle.name == node.right.particle.name)
        + _contact_tree_same_input_node_count(node.left)
        + _contact_tree_same_input_node_count(node.right)
    )


def _contact_tree_source_leg(node: _ContactTreeNode) -> int:
    return -1 if node.physical_leg is None else node.physical_leg


def _contact_term_is_color_singlet(term: CompiledVertexTerm) -> bool:
    return term.color_source in {"1", "UFO::{}::1"} or term.color_expression == "1"


def _deduplicate_contact_partials(
    auxiliary_particles: Sequence[CompiledParticleRecord],
    kernels: Sequence[CompiledOrientedKernel],
    terms: Sequence[CompiledVertexTerm],
) -> tuple[tuple[CompiledParticleRecord, ...], tuple[CompiledOrientedKernel, ...]]:
    term_by_id = {term.id: term for term in terms}
    auxiliary_by_name = {particle.name: particle for particle in auxiliary_particles}
    representative_by_signature: dict[tuple[object, ...], str] = {}
    replacement: dict[str, str] = {}
    retained: list[CompiledOrientedKernel] = []

    for kernel in kernels:
        if not kernel.vertex.endswith("::contact-partial"):
            retained.append(kernel)
            continue
        auxiliary = auxiliary_by_name[kernel.particles[2]]
        term = term_by_id[kernel.term_id]
        substitutions = {
            S(f"UFO::{name}"): E(term.coupling_expression)
            for name in kernel.runtime_parameters
            if name.startswith("derived_coupling_")
        }
        normalized_components = tuple(
            _replace_expression_symbols(
                _remap_kernel_symbols(
                    E(component),
                    old_kind=kernel.kind,
                    new_kind=0,
                ),
                substitutions,
            ).to_canonical_string()
            for component in kernel.component_expressions
        )
        signature = (
            kernel.particles[:2],
            auxiliary.component_dimension,
            auxiliary.color,
            kernel.coupling_orders,
            kernel.color_source,
            kernel.lc_color_normalization_power,
            normalized_components,
        )
        representative = representative_by_signature.get(signature)
        if representative is None:
            representative_by_signature[signature] = auxiliary.name
            retained.append(kernel)
        else:
            replacement[auxiliary.name] = representative

    if not replacement:
        return tuple(auxiliary_particles), tuple(kernels)

    rewritten = tuple(
        replace(
            kernel,
            particles=tuple(
                replacement.get(name, name) for name in kernel.particles
            ),
        )
        for kernel in retained
    )
    particles = tuple(
        particle
        for particle in auxiliary_particles
        if particle.name not in replacement
    )
    return particles, rewritten


def _fuse_contact_finals(
    kernels: Sequence[CompiledOrientedKernel],
    terms: Sequence[CompiledVertexTerm],
) -> tuple[CompiledOrientedKernel, ...]:
    coupling_by_term = {term.id: term.coupling_expression for term in terms}
    groups: dict[tuple[object, ...], list[CompiledOrientedKernel]] = {}
    passthrough: list[CompiledOrientedKernel] = []
    for kernel in kernels:
        if not kernel.vertex.endswith("::contact-final"):
            passthrough.append(kernel)
            continue
        key = (
            kernel.particles,
            kernel.coupling_orders,
            kernel.color_source,
            kernel.lc_color_normalization_power,
        )
        groups.setdefault(key, []).append(kernel)

    fused: list[CompiledOrientedKernel] = []
    for members in groups.values():
        first = members[0]
        runtime_aliases: dict[Expression, Expression] = {}
        canonical_runtime_name: str | None = None
        if all(
            coupling_by_term.get(member.term_id)
            == coupling_by_term.get(first.term_id)
            for member in members
        ):
            canonical_runtime_name = f"derived_coupling_{first.term_id}"
            canonical_runtime = S(f"UFO::{canonical_runtime_name}")
            runtime_aliases = {
                S(f"UFO::derived_coupling_{member.term_id}"): canonical_runtime
                for member in members
            }
        components = tuple(
            _canonicalize_oriented_kernel_component(
                sum(
                    (
                        _replace_expression_symbols(
                            _remap_kernel_symbols(
                                E(member.component_expressions[index]),
                                old_kind=member.kind,
                                new_kind=first.kind,
                            ),
                            runtime_aliases,
                        )
                        for member in members
                    ),
                    E("0"),
                )
            ).to_canonical_string()
            for index in range(len(first.component_expressions))
        )
        fused.append(
            replace(
                first,
                vertex=first.vertex.replace(
                    "::contact-final",
                    "::contact-final-fused",
                ),
                component_expressions=components,
                runtime_parameters=(
                    (canonical_runtime_name,)
                    if canonical_runtime_name is not None
                    else tuple(
                        sorted(
                            {
                                name
                                for member in members
                                for name in member.runtime_parameters
                            }
                        )
                    )
                ),
                term_ids=tuple(
                    term_id
                    for member in members
                    for term_id in (member.term_ids or (member.term_id,))
                ),
            )
        )
    return tuple(sorted((*passthrough, *fused), key=lambda kernel: kernel.kind))


def _contact_partial_component_expressions(
    term: CompiledVertexTerm,
    particle_by_name: Mapping[str, CompiledParticleRecord],
    *,
    left_leg: int,
    right_leg: int,
    open_legs: tuple[int, ...],
    kind: int,
) -> tuple[str, ...]:
    library = TensorLibrary.hep_lib_atom()
    expression = E(term.lorentz_expression)
    particles = tuple(particle_by_name[name] for name in term.particles)
    expression *= _input_tensor_expression(
        library,
        kind=kind,
        side="left",
        spin=particles[left_leg].spin,
        leg=left_leg + 1,
        components=_component_symbols(kind, "left", particles[left_leg].spin),
    )
    expression *= _input_tensor_expression(
        library,
        kind=kind,
        side="right",
        spin=particles[right_leg].spin,
        leg=right_leg + 1,
        components=_component_symbols(kind, "right", particles[right_leg].spin),
    )
    result = _execute_dense_tensor(expression, library)
    expected = math.prod(
        _spin_dimension(particles[leg].spin) for leg in open_legs
    )
    if len(result) != expected:
        raise ValueError(
            f"contact partial {term.vertex}/{term.id} produced {len(result)} "
            f"components, expected {expected}"
        )
    return tuple(
        _replace_evaluator_constants(_as_expression(result[index])).to_canonical_string()
        for index in range(len(result))
    )


def _contact_final_component_expressions(
    particles: Sequence[CompiledParticleRecord],
    auxiliary: CompiledParticleRecord,
    *,
    open_legs: tuple[int, ...],
    remaining_leg: int,
    result_leg: int,
    kind: int,
    auxiliary_on_left: bool,
    component_expansion: tuple[tuple[int, int] | None, ...],
) -> tuple[str, ...]:
    library = TensorLibrary.hep_lib_atom()
    auxiliary_side = "left" if auxiliary_on_left else "right"
    physical_side = "right" if auxiliary_on_left else "left"
    auxiliary_symbols = tuple(
        S(f"pyamplicol::kernel_{kind}_{auxiliary_side}_{component}")
        for component in range(auxiliary.component_dimension or 0)
    )
    expanded_auxiliary = tuple(
        E("0")
        if entry is None
        else entry[1] * auxiliary_symbols[entry[0]]
        for entry in component_expansion
    )
    expression = _contact_auxiliary_tensor_expression(
        library,
        kind=kind,
        side=auxiliary_side,
        particles=particles,
        open_legs=open_legs,
        components=expanded_auxiliary,
    )
    expression *= _input_tensor_expression(
        library,
        kind=kind,
        side=physical_side,
        spin=particles[remaining_leg].spin,
        leg=remaining_leg + 1,
        components=_component_symbols(
            kind,
            physical_side,
            particles[remaining_leg].spin,
        ),
    )
    result = _execute_dense_tensor(expression, library)
    expected = _spin_dimension(particles[result_leg].spin)
    if len(result) != expected:
        raise ValueError(
            f"contact final for source leg {result_leg} produced {len(result)} "
            f"components, expected {expected}"
        )
    return tuple(
        _replace_evaluator_constants(_as_expression(result[index])).to_canonical_string()
        for index in range(len(result))
    )


def _contact_auxiliary_tensor_expression(
    library: TensorLibrary,
    *,
    kind: int,
    side: str,
    particles: Sequence[CompiledParticleRecord],
    open_legs: tuple[int, ...],
    components: Sequence[Expression],
) -> Expression:
    representations = tuple(
        representation
        for leg in open_legs
        for representation in _spin_representations(particles[leg].spin)
    )
    if not representations:
        if len(components) != 1:
            raise ValueError("scalar contact auxiliary must have one component")
        return components[0]
    slots = tuple(
        slot
        for leg in open_legs
        for slot in _spin_slots(particles[leg].spin, leg + 1)
    )
    name = TensorName(f"pyamplicol::kernel_{kind}_{side}")
    library.register(LibraryTensor.dense(name(*representations), components))
    return name(*slots).to_expression()


def _execute_dense_tensor(expression: Expression, library: TensorLibrary):
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    result.to_dense()
    return result


def _contact_auxiliary_color(
    term: CompiledVertexTerm,
    particles: Sequence[CompiledParticleRecord],
    *,
    remaining_leg: int,
    result_leg: int,
) -> int:
    colors = tuple(particle.color for particle in particles)
    if all(color == 1 for color in colors):
        return 1
    if "f(" in term.color_source or "::f(" in term.color_source:
        return 8
    remaining = abs(colors[remaining_leg])
    result = abs(colors[result_leg])
    if remaining == 1:
        return colors[result_leg]
    if result == 1:
        return colors[remaining_leg]
    if remaining == result == 8:
        return 1
    return 1


def _compress_contact_components(
    components: Sequence[str],
) -> tuple[tuple[int, ...], tuple[tuple[int, int] | None, ...]]:
    representatives: list[Expression] = []
    representative_indices: list[int] = []
    expansion: list[tuple[int, int] | None] = []
    zero = E("0")
    for index, source in enumerate(components):
        expression = E(source)
        if str(expression) == "0":
            expansion.append(None)
            continue
        match: tuple[int, int] | None = None
        for basis_index, representative in enumerate(representatives):
            if str((expression - representative).expand()) == "0":
                match = (basis_index, 1)
                break
            if str((expression + representative).expand()) == "0":
                match = (basis_index, -1)
                break
        if match is None:
            match = (len(representatives), 1)
            representatives.append(expression)
            representative_indices.append(index)
        expansion.append(match)
    if not representatives:
        representatives.append(zero)
        representative_indices.append(0)
    return tuple(representative_indices), tuple(expansion)


def _four_point_contact_color_split(
    term: CompiledVertexTerm,
    result_leg: int,
) -> tuple[
    tuple[int, int],
    int,
    str,
    str,
    int,
    int,
    tuple[int, ...],
    tuple[int, ...],
    int,
]:
    factors = _normalized_structure_constant_factors(term.color_expression)
    if len(factors) == 2:
        shared_dummies = set(value for value in factors[0] if value < 0) & set(
            value for value in factors[1] if value < 0
        )
        if len(shared_dummies) == 1:
            dummy = next(iter(shared_dummies))
            result_index = result_leg + 1
            final_factor = next(
                (factor for factor in factors if result_index in factor),
                None,
            )
            if final_factor is not None:
                outer_factor = factors[1] if final_factor is factors[0] else factors[0]
                pair = tuple(value - 1 for value in outer_factor if value > 0)
                remaining = tuple(
                    value - 1
                    for value in final_factor
                    if value > 0 and value != result_index
                )
                if len(pair) == 2 and len(remaining) == 1:
                    canonical_f = "UFO::{}::f(1,2,3)"
                    return (
                        (pair[0], pair[1]),
                        remaining[0],
                        canonical_f,
                        canonical_f,
                        1,
                        1,
                        outer_factor,
                        final_factor,
                        dummy,
                    )

    input_legs = tuple(leg for leg in range(4) if leg != result_leg)
    return (
        (input_legs[0], input_legs[1]),
        input_legs[2],
        term.color_source,
        "1",
        term.lc_color_normalization_power,
        0,
        (),
        (),
        -1,
    )


def _normalized_structure_constant_factors(
    expression: str,
) -> tuple[tuple[int, ...], ...]:
    """Return typed f-tensor index words from a normalized color monomial."""

    if expression.count("::f(") != 2:
        return ()
    result: list[tuple[int, ...]] = []
    for arguments in _function_arguments(expression, "::f"):
        indices: list[int] = []
        for argument in arguments:
            dummy = re.search(r"ufo_c_dummy_([0-9]+)_adjoint", argument)
            if dummy is not None:
                indices.append(-int(dummy.group(1)))
                continue
            external = re.search(r"ufo_c_([0-9]+)", argument)
            if external is None:
                return ()
            indices.append(int(external.group(1)))
        if len(indices) != 3:
            return ()
        result.append(tuple(indices))
    return tuple(result)


def _function_arguments(source: str, head_suffix: str) -> tuple[tuple[str, ...], ...]:
    """Extract balanced canonical arguments for one fully qualified function."""

    marker = head_suffix + "("
    result: list[tuple[str, ...]] = []
    cursor = 0
    while True:
        marker_start = source.find(marker, cursor)
        if marker_start < 0:
            break
        open_index = marker_start + len(head_suffix)
        depth = 0
        brace_depth = 0
        argument_start = open_index + 1
        arguments: list[str] = []
        index = open_index
        while index < len(source):
            character = source[index]
            if character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    arguments.append(source[argument_start:index])
                    result.append(tuple(arguments))
                    cursor = index + 1
                    break
            elif character == "{":
                brace_depth += 1
            elif character == "}":
                brace_depth -= 1
            elif character == "," and depth == 1 and brace_depth == 0:
                arguments.append(source[argument_start:index])
                argument_start = index + 1
            index += 1
        else:
            raise ValueError(f"unbalanced canonical function {head_suffix}")
    return tuple(result)


def _permutation_sign(
    actual: tuple[int, ...],
    canonical: tuple[int, ...],
) -> int:
    if sorted(actual) != sorted(canonical):
        raise ValueError(
            f"color-factor indices {actual} do not match local orientation {canonical}"
        )
    positions = {value: index for index, value in enumerate(canonical)}
    permutation = tuple(positions[value] for value in actual)
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def _fuse_oriented_kernels(
    kernels: Sequence[CompiledOrientedKernel],
) -> tuple[CompiledOrientedKernel, ...]:
    groups: dict[tuple[object, ...], list[CompiledOrientedKernel]] = {}
    for kernel in kernels:
        key = (
            kernel.vertex,
            kernel.particles,
            kernel.source_particle_legs,
            kernel.coupling_orders,
            kernel.color_source,
            kernel.color_expression,
            kernel.lc_color_normalization_power,
        )
        groups.setdefault(key, []).append(kernel)

    fused: list[CompiledOrientedKernel] = []
    for members in groups.values():
        kind = len(fused)
        first = members[0]
        remapped_components = [
            tuple(
                _remap_kernel_symbols(
                    E(component),
                    old_kind=member.kind,
                    new_kind=kind,
                )
                for component in member.component_expressions
            )
            for member in members
        ]
        components = tuple(
            _canonicalize_oriented_kernel_component(
                sum(
                    (
                        remapped[index]
                        * S(f"UFO::derived_coupling_{member.term_id}")
                        for member, remapped in zip(
                            members,
                            remapped_components,
                            strict=True,
                        )
                    ),
                    E("0"),
                )
            )
            for index in range(len(remapped_components[0]))
        )
        coupling_expression = "1"
        fused.append(
            CompiledOrientedKernel(
                kind=kind,
                term_id=first.term_id,
                vertex=first.vertex,
                particles=first.particles,
                source_particle_legs=first.source_particle_legs,
                component_expressions=tuple(
                    component.to_canonical_string() for component in components
                ),
                coupling_expression=coupling_expression,
                coupling_orders=first.coupling_orders,
                runtime_parameters=tuple(
                    sorted(
                        {
                            f"derived_coupling_{member.term_id}"
                            for member in members
                        }
                    )
                ),
                color_source=first.color_source,
                color_expression=first.color_expression,
                lc_color_normalization_power=first.lc_color_normalization_power,
                term_ids=tuple(
                    term_id
                    for member in members
                    for term_id in (member.term_ids or (member.term_id,))
                ),
            )
        )
    return tuple(fused)


def _canonicalize_oriented_kernel_component(expression: Expression) -> Expression:
    """Cancel numeric sums, factor couplings, and group primitive inputs."""

    return expression.expand_num().collect_factors().collect_horner()


def _remap_kernel_symbols(
    expression: Expression,
    *,
    old_kind: int,
    new_kind: int,
    swap_sides: bool = False,
) -> Expression:
    source = expression.to_canonical_string()
    symbols = {
        (side, momentum_marker == "momentum_", int(index))
        for side, momentum_marker, index in re.findall(
            rf"kernel_{old_kind}_(left|right)_(momentum_)?([0-9]+)",
            source,
        )
    }
    replacements: list[Replacement] = []
    for side, is_momentum, index in symbols:
        target_side = (
            "right" if side == "left" else "left"
        ) if swap_sides else side
        if is_momentum:
            replacements.append(
                Replacement(
                    S(
                        f"pyamplicol::kernel_{old_kind}_{side}_momentum_{index}"
                    ),
                    S(
                        f"pyamplicol::kernel_{new_kind}_{target_side}_momentum_{index}"
                    ),
                )
            )
        else:
            replacements.append(
                Replacement(
                    S(f"pyamplicol::kernel_{old_kind}_{side}_{index}"),
                    S(f"pyamplicol::kernel_{new_kind}_{target_side}_{index}"),
                )
            )
    return expression.replace_multiple(replacements) if replacements else expression


def _replace_expression_symbols(
    expression: Expression,
    substitutions: Mapping[Expression, Expression],
) -> Expression:
    result = expression
    for source, target in substitutions.items():
        result = result.replace(source, target)
    return result


def _oriented_component_expressions(
    term: CompiledVertexTerm,
    particle_by_name: Mapping[str, CompiledParticleRecord],
    *,
    left_leg: int,
    right_leg: int,
    result_leg: int,
    kind: int,
    use_transverse_massless_yang_mills: bool = False,
) -> tuple[str, ...]:
    library = TensorLibrary.hep_lib_atom()
    expression = E(term.lorentz_expression)
    particles = tuple(particle_by_name[name] for name in term.particles)
    left_symbols = _component_symbols(kind, "left", particles[left_leg].spin)
    right_symbols = _component_symbols(kind, "right", particles[right_leg].spin)
    expression *= _input_tensor_expression(
        library,
        kind=kind,
        side="left",
        spin=particles[left_leg].spin,
        leg=left_leg + 1,
        components=left_symbols,
    )
    expression *= _input_tensor_expression(
        library,
        kind=kind,
        side="right",
        spin=particles[right_leg].spin,
        leg=right_leg + 1,
        components=right_symbols,
    )
    left_momentum = tuple(
        S(f"pyamplicol::kernel_{kind}_left_momentum_{component}")
        for component in range(4)
    )
    right_momentum = tuple(
        S(f"pyamplicol::kernel_{kind}_right_momentum_{component}")
        for component in range(4)
    )
    result_momentum = tuple(
        -(left_momentum[component] + right_momentum[component])
        for component in range(4)
    )
    momentum_by_leg = {
        left_leg: left_momentum,
        right_leg: right_momentum,
        result_leg: result_momentum,
    }
    minkowski = Representation.mink(4)
    for leg, momentum in momentum_by_leg.items():
        library.register(
            LibraryTensor.dense(
                TensorName(f"pyamplicol::ufo_momentum_{leg + 1}")(minkowski),
                momentum,
            )
        )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    result.to_dense()
    expected_dimension = _spin_dimension(particles[result_leg].spin)
    if len(result) != expected_dimension:
        raise ValueError(
            f"oriented kernel {term.vertex}/{term.id} produced {len(result)} "
            f"components for spin {particles[result_leg].spin}, expected "
            f"{expected_dimension}"
        )
    components = tuple(
        _replace_evaluator_constants(_as_expression(result[index]))
        for index in range(len(result))
    )
    if all(particle.spin == 3 for particle in particles):
        compact = _compact_yang_mills_three_vector_components(
            left_leg=left_leg,
            right_leg=right_leg,
            result_leg=result_leg,
            left_components=left_symbols,
            right_components=right_symbols,
            momentum_by_leg=momentum_by_leg,
        )
        scale = _equivalent_component_scale(components, compact)
        if scale is not None:
            if use_transverse_massless_yang_mills:
                compact = _transverse_yang_mills_three_vector_components(
                    left_components=left_symbols,
                    right_components=right_symbols,
                    left_momentum=left_momentum,
                    right_momentum=right_momentum,
                )
            components = tuple(scale * component for component in compact)
    return tuple(
        component.to_canonical_string() for component in components
    )


def _compact_yang_mills_three_vector_components(
    *,
    left_leg: int,
    right_leg: int,
    result_leg: int,
    left_components: Sequence[Expression],
    right_components: Sequence[Expression],
    momentum_by_leg: Mapping[int, Sequence[Expression]],
) -> tuple[Expression, ...]:
    """Return the compact canonical Yang-Mills three-vector contraction.

    The caller retains this form only after proving algebraic equivalence to the
    fully materialized spenso tensor. This makes the optimization independent
    of UFO Lorentz names and source-expression layout.
    """

    components_by_leg = {
        left_leg: tuple(left_components),
        right_leg: tuple(right_components),
    }
    terms = (
        ((0, 1), 0, 1, 2),
        ((0, 2), 2, 0, 1),
        ((1, 2), 1, 2, 0),
    )
    outputs: list[Expression] = []
    for component in range(4):
        output = E("0")
        for (
            metric_legs,
            positive_momentum_leg,
            negative_momentum_leg,
            vector_leg,
        ) in terms:
            momentum = tuple(
                momentum_by_leg[positive_momentum_leg][index]
                - momentum_by_leg[negative_momentum_leg][index]
                for index in range(4)
            )
            if result_leg == vector_leg:
                output += momentum[component] * _minkowski_dot(
                    components_by_leg[metric_legs[0]],
                    components_by_leg[metric_legs[1]],
                )
                continue
            other_metric_leg = (
                metric_legs[1]
                if metric_legs[0] == result_leg
                else metric_legs[0]
            )
            output += (
                components_by_leg[other_metric_leg][component]
                * _minkowski_dot(
                    momentum,
                    components_by_leg[vector_leg],
                )
            )
        outputs.append(output)
    return tuple(outputs)


def _transverse_yang_mills_three_vector_components(
    *,
    left_components: Sequence[Expression],
    right_components: Sequence[Expression],
    left_momentum: Sequence[Expression],
    right_momentum: Sequence[Expression],
) -> tuple[Expression, ...]:
    """Return the transverse Berends-Giele massless gauge current.

    The full Yang-Mills contraction differs only by terms proportional to each
    parent current's self-momentum contraction. AmpliCol removes those terms
    for massless adjoint gauge currents. Applying the same tensor-derived
    reduction preserves its fixed-width recursion convention without relying
    on model-specific particle or Lorentz names.
    """

    left = tuple(left_components)
    right = tuple(right_components)
    left_p = tuple(left_momentum)
    right_p = tuple(right_momentum)
    dot = _minkowski_dot(left, right)
    left_dot_right_p = _minkowski_dot(left, right_p)
    right_dot_left_p = _minkowski_dot(right, left_p)
    return tuple(
        dot * (left_p[index] - right_p[index])
        + 2 * (
            left_dot_right_p * right[index]
            - right_dot_left_p * left[index]
        )
        for index in range(4)
    )


def _equivalent_component_scale(
    materialized: Sequence[Expression],
    compact: Sequence[Expression],
) -> Expression | None:
    if len(materialized) != len(compact) or not materialized:
        return None
    for scale in (E("1"), E("-1")):
        if all(
            (dense - scale * candidate).expand() == E("0")
            for dense, candidate in zip(materialized, compact, strict=True)
        ):
            return scale

    first_dense, first_compact = next(
        (
            (dense, candidate)
            for dense, candidate in zip(materialized, compact, strict=True)
            if candidate != E("0")
        ),
        (E("0"), E("0")),
    )
    if first_compact == E("0"):
        return None
    scale = (first_dense / first_compact).cancel()
    if scale.get_all_symbols(False):
        return None
    if not all(
        (dense - scale * candidate).expand() == E("0")
        for dense, candidate in zip(materialized, compact, strict=True)
    ):
        return None
    return scale


def _is_compile_time_zero_parameter(
    name: str,
    parameters: Mapping[str, CompiledParameterRecord],
) -> bool:
    if name.upper() == "ZERO":
        return True
    parameter = parameters.get(name)
    if parameter is None or parameter.nature != "internal":
        return False
    expression = E(parameter.resolved_expression)
    if expression.get_all_symbols(False):
        return False
    try:
        return complex(expression.evaluate({})) == 0.0
    except (TypeError, ValueError):
        return False


def _is_single_structure_constant(expression: str) -> bool:
    return (
        expression.count("::f(") == 1
        and len(_function_arguments(expression, "::f")) == 1
        and "::t(" not in expression.lower()
        and "::d(" not in expression.lower()
    )


def _minkowski_dot(
    left: Sequence[Expression],
    right: Sequence[Expression],
) -> Expression:
    if len(left) != 4 or len(right) != 4:
        raise ValueError("Minkowski dot products require four components")
    return left[0] * right[0] - sum(
        (left[index] * right[index] for index in range(1, 4)),
        E("0"),
    )


def _input_tensor_expression(
    library: TensorLibrary,
    *,
    kind: int,
    side: str,
    spin: int,
    leg: int,
    components: Sequence[Expression],
) -> Expression:
    representations = _spin_representations(spin)
    if not representations:
        if len(components) != 1:
            raise ValueError("scalar current must have exactly one component")
        return components[0]
    name = TensorName(f"pyamplicol::kernel_{kind}_{side}")
    library.register(LibraryTensor.dense(name(*representations), components))
    slots = _spin_slots(spin, leg)
    return name(*slots).to_expression()


def _spin_representations(spin: int) -> tuple[Representation, ...]:
    minkowski = Representation.mink(4)
    if spin in {-1, 1}:
        return ()
    if spin == 2:
        return (Representation.bis(4),)
    if spin == 3:
        return (minkowski,)
    if spin == 5:
        return (minkowski, minkowski)
    raise ValueError(f"unsupported UFO spin code {spin}")


def _spin_slots(spin: int, leg: int):
    representations = _spin_representations(spin)
    if spin == 2:
        return (representations[0](f"ufo_s_1_{leg}"),)
    if spin == 3:
        return (representations[0](f"ufo_l_1_{leg}"),)
    if spin == 5:
        return (
            representations[0](f"ufo_l_1_{leg}"),
            representations[1](f"ufo_l_2_{leg}"),
        )
    return ()


def _spin_dimension(spin: int) -> int:
    return {-1: 1, 1: 1, 2: 4, 3: 4, 5: 16}[spin]


def _component_symbols(kind: int, side: str, spin: int) -> tuple[Expression, ...]:
    return tuple(
        S(f"pyamplicol::kernel_{kind}_{side}_{component}")
        for component in range(_spin_dimension(spin))
    )


def _lc_color_normalization_power(source: str) -> int:
    """Count normalized non-Abelian tensors in one UFO color monomial."""

    expression = E(source)
    return sum(
        len(list(expression.match(E(pattern))))
        for pattern in (
            "UFO::T(a_,b_,c_)",
            "UFO::f(a_,b_,c_)",
            "UFO::d(a_,b_,c_)",
        )
    )


def _as_expression(value: object) -> Expression:
    if isinstance(value, Expression):
        return value
    raise TypeError(f"spenso kernel component is not an Expression: {type(value).__name__}")


def cast_tuple3(value: object) -> tuple[str, str, str]:
    values = tuple(str(item) for item in _sequence(value))
    if len(values) != 3:
        raise ValueError("oriented kernel particles must have length three")
    return values[0], values[1], values[2]


def cast_int_tuple3(value: object) -> tuple[int, int, int]:
    values = tuple(int(item) for item in _sequence(value))
    if len(values) != 3:
        raise ValueError("oriented kernel source legs must have length three")
    return values[0], values[1], values[2]


def _order(item: Mapping[str, object]) -> CompiledCouplingOrder:
    return CompiledCouplingOrder(
        name=str(item["name"]),
        expansion_order=int(item["expansion_order"]),
        hierarchy=int(item["hierarchy"]),
    )


def _parameter(item: Mapping[str, object]) -> CompiledParameterRecord:
    expression = _optional_string(item.get("expression"))
    value = _optional_pair(item.get("value"))
    return CompiledParameterRecord(
        name=str(item["name"]),
        nature=str(item["nature"]),
        parameter_type=str(item["parameter_type"]),
        value=value,
        expression=expression,
        resolved_expression=(
            expression
            if expression is not None
            else _numeric_expression(value or (0.0, 0.0)).to_canonical_string()
        ),
        lhablock=_optional_string(item.get("lhablock")),
        lhacode=tuple(int(value) for value in _sequence(item.get("lhacode"))),
    )


def _particle(item: Mapping[str, object]) -> CompiledParticleRecord:
    return CompiledParticleRecord(
        name=str(item["name"]),
        antiname=str(item["antiname"]),
        pdg_code=int(item["pdg_code"]),
        spin=int(item["spin"]),
        color=int(item["color"]),
        mass=str(item["mass"]),
        width=str(item["width"]),
        charge=float(item.get("charge", 0.0)),
        ghost_number=int(item.get("ghost_number", 0)),
        propagating=bool(item.get("propagating", True)),
        goldstoneboson=bool(item.get("goldstoneboson", False)),
        propagator=_optional_string(item.get("propagator")),
    )


def _coupling(item: Mapping[str, object]) -> CompiledCouplingRecord:
    expression = str(item["expression"])
    return CompiledCouplingRecord(
        name=str(item["name"]),
        expression=expression,
        resolved_expression=expression,
        value=_pair(item.get("value")),
        orders=_orders(item.get("orders")),
    )


def _propagator(
    item: Mapping[str, object],
    particles: Sequence[CompiledParticleRecord],
) -> CompiledPropagatorRecord:
    particle_name = str(item["particle"])
    particle = next(
        (candidate for candidate in particles if candidate.name == particle_name),
        None,
    )
    linked_name = None if particle is None else particle.propagator
    name = str(item["name"])
    return CompiledPropagatorRecord(
        name=name,
        particle=particle_name,
        numerator=str(item["numerator"]),
        denominator=str(item["denominator"]),
        custom=linked_name == name and not name.endswith("_propFeynman"),
    )


def _orders(value: object) -> tuple[tuple[str, int], ...]:
    result = []
    for pair in _sequence(value):
        values = _sequence(pair)
        if len(values) != 2:
            raise ValueError("coupling order must be [name, value]")
        result.append((str(values[0]), int(values[1])))
    return tuple(result)


def _pair(value: object) -> tuple[float, float]:
    pair = _sequence(value)
    if len(pair) != 2:
        raise ValueError("complex value must be [real, imaginary]")
    return float(pair[0]), float(pair[1])


def _optional_pair(value: object) -> tuple[float, float] | None:
    return None if value is None else _pair(value)


def _optional_string(value: object) -> str | None:
    return None if value is None else str(value)


def _resolve_parameter_records(
    records: Sequence[CompiledParameterRecord],
) -> tuple[CompiledParameterRecord, ...]:
    by_name = {record.name: record for record in records}
    symbols = {name: S(f"UFO::{name}") for name in by_name}
    resolved: dict[str, Expression] = {}
    active: list[str] = []

    def resolve(name: str) -> Expression:
        if name in resolved:
            return resolved[name]
        if name in active:
            cycle = " -> ".join((*active, name))
            raise ValueError(f"cyclic UFO parameter definitions: {cycle}")
        record = by_name[name]
        if record.nature == "external":
            resolved[name] = symbols[name]
            return resolved[name]
        active.append(name)
        expression = (
            E(record.expression)
            if record.expression is not None
            else _numeric_expression(record.value or (0.0, 0.0))
        )
        expression_symbols = set(expression.get_all_symbols(False))
        for dependency, symbol in symbols.items():
            if dependency == name or symbol not in expression_symbols:
                continue
            expression = expression.replace(symbol, resolve(dependency))
        active.pop()
        resolved[name] = _replace_evaluator_constants(expression)
        return resolved[name]

    return tuple(
        replace(record, resolved_expression=resolve(record.name).to_canonical_string())
        for record in records
    )


def _resolve_coupling_records(
    records: Sequence[CompiledCouplingRecord],
    parameters: Sequence[CompiledParameterRecord],
) -> tuple[CompiledCouplingRecord, ...]:
    replacements = {
        S(f"UFO::{parameter.name}"): E(parameter.resolved_expression)
        for parameter in parameters
        if parameter.nature != "external"
    }
    result = []
    for record in records:
        expression = E(record.expression)
        expression_symbols = set(expression.get_all_symbols(False))
        for symbol, replacement_expression in replacements.items():
            if symbol in expression_symbols:
                expression = expression.replace(symbol, replacement_expression)
        expression = _replace_evaluator_constants(expression)
        result.append(
            replace(
                record,
                resolved_expression=expression.to_canonical_string(),
            )
        )
    return tuple(result)


def _replace_evaluator_constants(expression: Expression) -> Expression:
    # SymJIT 2.19.3 does not currently lower Symbolica's built-in Pi atom,
    # although ordinary numeric constants and powers are supported.
    return expression.replace(E("pi"), E(repr(math.pi)))


def _numeric_expression(value: tuple[float, float]) -> Expression:
    real, imaginary = value
    return E(repr(real)) + E(repr(imaginary)) * E("1𝑖")


def _sequence(value: object) -> list[object]:
    return list(value) if isinstance(value, (list, tuple)) else []


def _mappings(value: object) -> list[dict[str, object]]:
    result = []
    for item in _sequence(value):
        if not isinstance(item, dict):
            raise ValueError("compiled model list entry must be an object")
        result.append({str(key): element for key, element in item.items()})
    return result


__all__ = [
    "CompiledCouplingOrder",
    "CompiledCouplingRecord",
    "CompiledModelIR",
    "CompiledParameterRecord",
    "CompiledParticleRecord",
    "CompiledPropagatorRecord",
    "CompiledVertexTerm",
    "compile_builtin_model_ir",
    "compile_ufo_model_ir",
]
