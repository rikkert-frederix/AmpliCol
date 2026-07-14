from __future__ import annotations

import math
from collections import Counter
from copy import copy
from numbers import Number
from typing import Any, Mapping, Sequence, TYPE_CHECKING

from symbolica import E, S, Expression, Replacement
from symbolica.community.spenso import (
    LibraryTensor,
    Representation,
    TensorLibrary,
    TensorName,
    TensorNetwork,
)

from .model import (
    Model,
    Particle,
    PropagatorLoweringRule,
    QuantumFlow,
    SourceSpinState,
    Vertex,
    VertexEvaluationEquivalence,
    VertexLoweringRule,
    _expr_antiquark_propagator_dirac,
    _expr_antiquark_propagator_weyl,
    _expr_minkowski_dot,
    _expr_quark_propagator_dirac,
    _expr_quark_propagator_weyl,
    _minkowski_square_expression,
)
from .ufo_ir import (
    CompiledOrientedKernel,
    CompiledParameterRecord,
    _as_expression,
    _replace_evaluator_constants,
    _spin_representations,
    _spin_slots,
)
from .ufo_tensors import (
    classify_trilinear_color_expression,
    normalize_lorentz_expression,
)

if TYPE_CHECKING:
    from .model_source import CompiledModel


class CompiledUFOModel(Model):
    """Model boundary consumed by the existing generic DAG compiler."""

    def __init__(
        self,
        compiled: "CompiledModel",
        runtime_parameters: Mapping[str, Any] | None = None,
    ) -> None:
        self.compiled = compiled
        self.name = compiled.name
        self._particle_records_by_name = {
            particle.name: particle for particle in compiled.ir.particles
        }
        self._particle_records_by_pdg = {
            particle.pdg_code: particle for particle in compiled.ir.particles
        }
        self._parameter_records = {
            parameter.name: parameter for parameter in compiled.ir.parameters
        }
        self._propagator_records_by_particle_name = {
            propagator.particle: propagator
            for propagator in compiled.ir.propagators
        }
        self._vertex_terms = {
            term.id: term for term in compiled.ir.vertex_terms
        }
        self._coupling_records_by_name = {
            coupling.name: coupling for coupling in compiled.ir.couplings
        }
        defaults = {
            name: complex(value[0], value[1])
            for name, value in compiled.parameter_defaults.items()
        }
        self._runtime_parameters = {**defaults, **dict(runtime_parameters or {})}
        self.particles = {
            record.pdg_code: Particle(
                pdg=record.pdg_code,
                anti_pdg=self._particle_records_by_name[record.antiname].pdg_code,
                spin=record.spin,
                dimension=(
                    record.component_dimension
                    if record.component_dimension is not None
                    else _spin_dimension(record.spin)
                ),
                color_rep=record.color,
                mass=float(complex(self._parameter_default(record.mass)).real),
                width=float(complex(self._parameter_default(record.width)).real),
                charge=record.charge,
            )
            for record in compiled.ir.particles
        }
        self._kernels = {
            kernel.kind: kernel for kernel in compiled.ir.oriented_kernels
        }
        self._kernel_component_expression_cache: dict[
            int, tuple[Expression, ...]
        ] = {}
        self._kernel_coupling_expression_cache: dict[int, Expression] = {}
        self._expression_symbol_cache: dict[str, Expression] = {}
        self._kernel_function_specs: dict[
            tuple[int, int, int, int],
            tuple[tuple[Expression, tuple[int, ...]], ...],
        ] = {}
        self._symbolica_kernel_functions: dict[
            tuple[Expression, tuple[Expression, ...]], Expression
        ] = {}
        self._runtime_derived_definitions = {
            name: self._vertex_terms[int(name.rsplit("_", 1)[1])].coupling_expression
            for kernel in self._kernels.values()
            for name in kernel.runtime_parameters
            if name.startswith("derived_coupling_")
        }
        self._runtime_derived_expression_cache: dict[str, Expression] = {}
        self._runtime_derived_default_cache: dict[str, complex] = {}
        self._runtime_derived_domain_cache: dict[str, str] = {}
        self._runtime_parameter_domain_cache: dict[str, str] = {}
        self._weyl_projection_support: dict[
            tuple[int, int, int, int], bool
        ] = {}
        self._custom_propagator_expressions: dict[
            int,
            tuple[Expression, tuple[str, ...]],
        ] = {}
        self._custom_propagator_templates: dict[int, tuple[Expression, ...]] = {}
        self._color_projection_cache: dict[int, tuple[str, complex]] = {}
        self.inactive_goldstone_names = frozenset(
            record.name
            for record in compiled.ir.particles
            if record.goldstoneboson
            and self._goldstone_is_redundant_in_unitary_gauge(record)
        )
        self.vertices = tuple(
            Vertex(
                kind=kernel.kind,
                particles=tuple(
                    self._particle_records_by_name[name].pdg_code
                    for name in kernel.particles
                ),
                coupling=(1.0, 0.0),
            )
            for kernel in compiled.ir.oriented_kernels
            if not any(
                name in self.inactive_goldstone_names
                for name in kernel.particles
            )
        )
        vertices_by_input: dict[tuple[int, int], list[Vertex]] = {}
        for vertex in self.vertices:
            vertices_by_input.setdefault(
                (vertex.particles[0], vertex.particles[1]),
                [],
            ).append(vertex)
        self._compiled_vertices_by_input = {
            key: tuple(vertices) for key, vertices in vertices_by_input.items()
        }

    def with_runtime_parameters(
        self,
        parameters: Mapping[str, Any],
    ) -> "CompiledUFOModel":
        model = copy(self)
        model._runtime_parameters = {
            **self._runtime_parameters,
            **dict(parameters),
        }
        model._runtime_derived_default_cache = {}
        return model

    def particle(self, pdg: int) -> Particle:
        try:
            return self.particles[int(pdg)]
        except KeyError as exc:
            raise KeyError(f"particle not in model: {pdg}") from exc

    def vertices_for_inputs(
        self,
        left_pdg: int,
        right_pdg: int,
        *,
        color_accuracy: str = "lc",
    ) -> tuple[Vertex, ...]:
        if color_accuracy not in {"lc", "nlc", "full"}:
            raise ValueError(f"unknown colour accuracy: {color_accuracy}")
        return self._compiled_vertices_by_input.get(
            (int(left_pdg), int(right_pdg)),
            (),
        )

    def anti_particle(self, pdg: int) -> int:
        return self.particle(pdg).anti_pdg

    def mass(self, pdg: int) -> Any:
        return self._real_parameter_value(
            self._particle_records_by_pdg[int(pdg)].mass,
            field="mass",
        )

    def width(self, pdg: int) -> Any:
        return self._real_parameter_value(
            self._particle_records_by_pdg[int(pdg)].width,
            field="width",
        )

    def spin(self, pdg: int) -> int:
        return self._particle_records_by_pdg[int(pdg)].spin

    def dimension(self, pdg: int) -> int:
        record = self._particle_records_by_pdg[int(pdg)]
        if record.component_dimension is not None:
            return record.component_dimension
        return _spin_dimension(record.spin)

    def current_dimension(self, particle_id: int, chirality: int = 0) -> int:
        if chirality != 0 and self.is_chiral_eligible(particle_id):
            return 2
        return self.dimension(particle_id)

    def color_rep(self, pdg: int) -> int:
        return self._particle_records_by_pdg[int(pdg)].color

    def color_dim(self, pdg: int) -> int:
        return abs(self.color_rep(pdg))

    def charge(self, pdg: int) -> float:
        return self._particle_records_by_pdg[int(pdg)].charge

    def is_fermion(self, pdg: int) -> bool:
        return self.spin(pdg) == 2

    def is_quark(self, pdg: int) -> bool:
        return self.is_fermion(pdg) and self.color_rep(pdg) == 3

    def is_antiquark(self, pdg: int) -> bool:
        return self.is_fermion(pdg) and self.color_rep(pdg) == -3

    def is_lepton(self, pdg: int) -> bool:
        return self.is_fermion(pdg) and self.color_rep(pdg) == 1 and pdg > 0

    def is_antilepton(self, pdg: int) -> bool:
        return self.is_fermion(pdg) and self.color_rep(pdg) == 1 and pdg < 0

    def is_chiral_eligible(self, pdg: int) -> bool:
        if not self.is_fermion(pdg):
            return False
        propagator = self._propagator_record(pdg)
        if propagator is not None and propagator.custom:
            return False
        particle = self._particle_records_by_pdg[int(pdg)]
        if particle.mass.upper() == "ZERO":
            return True
        record = self._parameter_records.get(particle.mass)
        return (
            record is not None
            and record.nature != "external"
            and self._parameter_default(particle.mass) == 0.0
        )

    def is_gluon(self, pdg: int) -> bool:
        return self.spin(pdg) == 3 and self.color_rep(pdg) == 8

    def is_singlet(self, pdg: int) -> bool:
        return self.color_rep(pdg) == 1

    def is_tensor(self, pdg: int) -> bool:
        return self.spin(pdg) == 5

    def is_massive_boson(self, pdg: int) -> bool:
        return self.spin(pdg) in {3, 5} and self._parameter_default(
            self._particle_records_by_pdg[int(pdg)].mass
        ) != 0.0

    def is_photon(self, pdg: int) -> bool:
        return self.spin(pdg) == 3 and self.color_rep(pdg) == 1 and self.mass(pdg) == 0

    def is_higgs(self, pdg: int) -> bool:
        return self.spin(pdg) == 1

    def source_spin_states(self, particle_id: int) -> tuple[SourceSpinState, ...]:
        if self.is_chiral_eligible(particle_id):
            return super().source_spin_states(particle_id)
        spin = self.spin(particle_id)
        massive = complex(self._parameter_default(
            self._particle_records_by_pdg[int(particle_id)].mass
        )).real != 0.0
        if spin == 1:
            helicities = (0,)
        elif spin == 2:
            helicities = (-1, 1)
        elif spin == 3:
            helicities = (-1, 0, 1) if massive else (-1, 1)
        elif spin == 5:
            helicities = (-2, -1, 0, 1, 2) if massive else (-2, 2)
        else:
            raise ValueError(f"unsupported source spin {spin} for particle {particle_id}")
        return tuple(
            SourceSpinState(helicity=helicity, chirality=0, spin_state=helicity)
            for helicity in helicities
        )

    def allowed_quantum_flows(
        self,
        vertex: Vertex,
        left_index: Any,
        right_index: Any,
    ) -> tuple[QuantumFlow, ...]:
        result_particle = vertex.particles[2]
        left_chirality = int(getattr(left_index, "chirality", 0))
        right_chirality = int(getattr(right_index, "chirality", 0))
        result_chiralities = (
            (-1, 1) if self.is_chiral_eligible(result_particle) else (0,)
        )
        return tuple(
            QuantumFlow(
                chirality=result_chirality,
                spin_state=(
                    result_chirality
                    if self.is_chiral_eligible(result_particle)
                    else 0
                ),
                flavour_flow=self.combine_flavour_flow(
                    result_particle,
                    left_index,
                    right_index,
                ),
                charge_flow=self.charge_units(result_particle),
                coupling=(1.0, 0.0),
            )
            for result_chirality in result_chiralities
            if self._weyl_projection_is_nonzero(
                vertex.kind,
                left_chirality,
                right_chirality,
                result_chirality,
            )
        )

    def vertex_lowering_rule(self, kind: int) -> VertexLoweringRule:
        kernel = self._kernel(kind)
        return VertexLoweringRule(
            kind=kind,
            backend="spenso-ufo",
            tensor_names=(kernel.vertex,),
            expression_head="compiled_ufo_kernel",
            full_tensor_network_ready=True,
            description="typed and oriented UFO tensor kernel",
            kernel="compiled_ufo_kernel",
            input_roles=(kernel.particles[0], kernel.particles[1]),
            output_role=kernel.particles[2],
            coupling_mode="external-model-parameters",
        )

    def vertex_evaluation_equivalence(
        self,
        kind: int,
    ) -> VertexEvaluationEquivalence:
        kernel = self._kernel(kind)
        if not kernel.evaluation_equivalence_verified or not kernel.evaluation_class:
            return super().vertex_evaluation_equivalence(kind)
        input_order = tuple(int(value) for value in kernel.evaluation_input_order)
        if input_order not in {(0, 1), (1, 0)}:
            raise ValueError(
                f"compiled UFO kernel {kind} has invalid evaluation input order "
                f"{input_order}"
            )
        return VertexEvaluationEquivalence(
            class_id=kernel.evaluation_class,
            factor=kernel.evaluation_factor,
            input_order=input_order,
            verified=True,
        )

    def vertex_coupling_orders(self, vertex: Vertex):
        return self._kernel(vertex.kind).coupling_orders

    def coupling_order_hierarchies(self) -> dict[str, int]:
        return {
            str(order.name).upper(): max(1, int(order.hierarchy))
            for order in self.compiled.ir.orders
        }

    def vertex_color_weight(
        self,
        vertex: Vertex,
        *,
        color_accuracy: str,
    ) -> tuple[float, float]:
        if color_accuracy not in {"lc", "nlc", "full"}:
            raise ValueError(f"unknown colour accuracy: {color_accuracy}")
        power = self._kernel(vertex.kind).lc_color_normalization_power
        normalization = 2.0 ** (-0.5 * power)
        structure, coefficient = self._vertex_color_projection(vertex)
        if structure in {
            "adjoint-structure-constant",
            "adjoint-structure-constant-product",
        }:
            phase = (-1j) ** power
            weight = coefficient * normalization * phase
            return (weight.real, weight.imag)
        weight = coefficient * normalization
        return (weight.real, weight.imag)

    def vertex_color_structure(self, vertex: Vertex) -> str:
        return self._vertex_color_projection(vertex)[0]

    def _vertex_color_projection(self, vertex: Vertex) -> tuple[str, complex]:
        cached = self._color_projection_cache.get(vertex.kind)
        if cached is not None:
            return cached
        kernel = self._kernel(vertex.kind)
        if (
            kernel.color_projection_structure is not None
            and kernel.color_projection_coefficient is not None
        ):
            projected = (
                kernel.color_projection_structure,
                complex(*kernel.color_projection_coefficient),
            )
        else:
            projected = classify_trilinear_color_expression(
                kernel.color_expression,
                kernel.color_source,
                tuple(
                    self.color_rep(particle_id)
                    for particle_id in vertex.particles
                ),
            )
        self._color_projection_cache[vertex.kind] = projected
        return projected

    def vertex_is_internal_contact_fragment(self, vertex: Vertex) -> bool:
        return "::contact-" in self._kernel(vertex.kind).vertex

    def vertex_closure_allowed(self, vertex: Vertex) -> bool:
        del vertex
        return False

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
        del coupling
        kernel = self._kernel(kind)
        expected_result = self._particle_records_by_name[kernel.particles[2]].pdg_code
        if int(result_particle_id) != expected_result:
            raise ValueError(
                f"kernel {kind} returns particle {expected_result}, not {result_particle_id}"
            )
        left_momentum = tuple(left_momentum or (0.0, 0.0, 0.0, 0.0))
        right_momentum = tuple(right_momentum or (0.0, 0.0, 0.0, 0.0))
        if len(left_momentum) != 4 or len(right_momentum) != 4:
            raise ValueError("UFO kernels require four-component input momenta")
        runtime_parameter_values = {
            name: self._parameter_value(name)
            for name in kernel.runtime_parameters
        }
        if not all(
            _is_numeric(value)
            for value in (
                *left,
                *right,
                *left_momentum,
                *right_momentum,
                *runtime_parameter_values.values(),
            )
        ):
            return self._kernel_function_component_calls(
                kernel,
                left,
                right,
                left_chirality=left_chirality,
                right_chirality=right_chirality,
                result_chirality=result_chirality,
                left_momentum=left_momentum,
                right_momentum=right_momentum,
                runtime_parameter_values=runtime_parameter_values,
            )
        components = self._projected_kernel_components(
            kernel,
            left,
            right,
            left_chirality=left_chirality,
            right_chirality=right_chirality,
            result_chirality=result_chirality,
            left_momentum=left_momentum,
            right_momentum=right_momentum,
            runtime_parameter_values=runtime_parameter_values,
        )
        coupling_expression = self._resolved_kernel_coupling_expression(
            kernel,
            runtime_parameter_values,
        )
        return tuple(component * coupling_expression for component in components)

    def symbolica_function_definitions(
        self,
    ) -> Mapping[tuple[Expression, tuple[Expression, ...]], Expression]:
        """Return lazily materialized kernel functions used by stage expressions."""

        return self._symbolica_kernel_functions

    def propagator_lowering_rule(
        self,
        particle_id: int,
        chirality: int = 0,
    ) -> PropagatorLoweringRule:
        auxiliary_kind = self.auxiliary_kind(particle_id)
        if auxiliary_kind is not None:
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=0,
                backend="identity",
                full_tensor_network_ready=True,
                applies_propagator=False,
                kernel="ufo_contact_auxiliary_no_propagator",
                description="synthetic UFO contact current with no propagator",
            )
        if chirality != 0 and self.is_chiral_eligible(particle_id):
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=chirality,
                backend="spenso-ufo-weyl",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="ufo_weyl_fermion",
                description="UFO massless fermion projected to a Weyl current",
            )
        spin = self.spin(particle_id)
        custom = self._propagator_record(particle_id)
        if custom is not None and custom.custom:
            return PropagatorLoweringRule(
                particle_id=particle_id,
                chirality=0,
                backend="spenso-ufo-custom",
                full_tensor_network_ready=True,
                applies_propagator=True,
                kernel="ufo_custom_propagator",
                description="model-supplied UFO propagator lowered without gauge conversion",
            )
        kernels = {
            1: "ufo_scalar_propagator",
            2: "ufo_dirac_propagator",
            3: (
                "ufo_massive_vector_propagator"
                if self.is_massive_boson(particle_id)
                else "ufo_massless_vector_propagator"
            ),
            5: (
                "ufo_massive_spin2_fierz_pauli"
                if self.is_massive_boson(particle_id)
                else "ufo_massless_spin2_de_donder"
            ),
        }
        return PropagatorLoweringRule(
            particle_id=particle_id,
            chirality=0,
            backend="symbolica-ufo",
            full_tensor_network_ready=spin in {1, 2, 3, 5},
            applies_propagator=True,
            kernel=kernels[spin],
            description="UFO propagator with runtime mass and width",
        )

    def propagator_component_expression(
        self,
        particle_id: int,
        value: Sequence[Any],
        momentum: Sequence[Any],
        *,
        chirality: int = 0,
    ) -> tuple[Any, ...]:
        if self.auxiliary_kind(particle_id) is not None:
            return tuple(value)
        if chirality != 0 and self.is_chiral_eligible(particle_id):
            components = tuple(value)
            current_momentum = tuple(momentum)
            return (
                _expr_antiquark_propagator_weyl(
                    components,
                    current_momentum,
                    chirality,
                )
                if int(particle_id) < 0
                else _expr_quark_propagator_weyl(
                    components,
                    current_momentum,
                    chirality,
                )
            )
        custom = self._propagator_record(particle_id)
        if custom is not None and custom.custom:
            return self._evaluate_custom_propagator(
                particle_id,
                value,
                momentum,
            )
        spin = self.spin(particle_id)
        if len(momentum) != 4:
            raise ValueError("UFO propagators require four-momentum components")
        mass = self.mass(particle_id)
        width = self.width(particle_id)
        if spin == 1:
            if len(value) != 1:
                raise ValueError("scalar propagator expects one current component")
            denominator = (
                _minkowski_square_expression(momentum)
                - mass * mass
                + 1j * mass * width
            )
            return (1j * value[0] / denominator,)
        if spin == 2:
            components = tuple(value)
            current_momentum = tuple(momentum)
            return (
                _expr_antiquark_propagator_dirac(
                    components,
                    current_momentum,
                    mass,
                    width,
                )
                if int(particle_id) < 0
                else _expr_quark_propagator_dirac(
                    components,
                    current_momentum,
                    mass,
                    width,
                )
            )
        if spin == 3:
            if len(value) != 4:
                raise ValueError("vector propagator expects four current components")
            current = tuple(value)
            current_momentum = tuple(momentum)
            denominator = (
                _minkowski_square_expression(current_momentum)
                - mass * mass
                + 1j * mass * width
            )
            if mass == 0.0:
                return tuple(-1j * component / denominator for component in current)
            longitudinal = _expr_minkowski_dot(current, current_momentum) / (
                mass * mass
            )
            return tuple(
                -1j
                * (current[index] - current_momentum[index] * longitudinal)
                / denominator
                for index in range(4)
            )
        if spin == 5:
            dimension = self._runtime_parameters.get("dim", 4.0)
            return _expr_spin2_propagator(
                tuple(value),
                tuple(momentum),
                mass,
                width,
                dimension=dimension,
                massive=self.is_massive_boson(particle_id),
            )
        raise NotImplementedError(f"generic spin-{spin} propagator is not implemented")

    def runtime_parameter_names_for_vertex(self, kind: int) -> tuple[str, ...]:
        return self._kernel(kind).runtime_parameters

    def runtime_derived_parameter_definitions(self) -> dict[str, str]:
        return dict(self._runtime_derived_definitions)

    def runtime_derived_parameter_definitions_for(
        self,
        names: Sequence[str],
    ) -> dict[str, str]:
        return {
            name: self._runtime_derived_definitions[name]
            for name in names
            if name in self._runtime_derived_definitions
        }

    def runtime_derived_parameter_domains(self) -> dict[str, str]:
        """Return phases that follow exactly from UFO parameter declarations."""

        return self.runtime_derived_parameter_domains_for(
            tuple(self._runtime_derived_definitions)
        )

    def runtime_derived_parameter_domains_for(
        self,
        names: Sequence[str],
    ) -> dict[str, str]:
        """Return provable phases for only the requested derived couplings."""

        parameter_symbols = {
            S(f"UFO::{parameter.name}"): self._runtime_domain_symbol(
                parameter.name,
                self._runtime_parameter_domain(parameter.name),
            )
            for parameter in self._parameter_records.values()
        }
        domains: dict[str, str] = {}
        for name in names:
            cached = self._runtime_derived_domain_cache.get(name)
            if cached is not None:
                domains[name] = cached
                continue
            if name not in self._runtime_derived_definitions:
                continue
            term = self._vertex_terms[int(name.rsplit("_", 1)[1])]
            coupling = self._coupling_records_by_name[term.coupling]
            expression = E(coupling.expression)
            for source, replacement in parameter_symbols.items():
                expression = expression.replace(source, replacement)
            if expression.is_real():
                domain = "real"
            elif (-1j * expression).is_real():
                domain = "imaginary"
            else:
                domain = "complex"
            self._runtime_derived_domain_cache[name] = domain
            domains[name] = domain
        return domains

    def _runtime_parameter_domain(
        self,
        name: str,
        visiting: frozenset[str] = frozenset(),
    ) -> str:
        cached = self._runtime_parameter_domain_cache.get(name)
        if cached is not None:
            return cached
        parameter = self._parameter_records[name]
        if parameter.parameter_type.lower() == "real":
            domain = "real"
        elif parameter.nature.lower() != "internal" or name in visiting:
            domain = "complex"
        else:
            expression = E(parameter.resolved_expression)
            symbols = set(expression.get_all_symbols(False))
            nested_visiting = visiting | {name}
            for dependency in self._parameter_records.values():
                source = S(f"UFO::{dependency.name}")
                if source not in symbols:
                    continue
                replacement = self._runtime_domain_symbol(
                    dependency.name,
                    self._runtime_parameter_domain(
                        dependency.name,
                        nested_visiting,
                    ),
                )
                expression = expression.replace(source, replacement)
            if expression.is_real():
                domain = "real"
            elif (-1j * expression).is_real():
                domain = "imaginary"
            else:
                domain = "complex"
        self._runtime_parameter_domain_cache[name] = domain
        return domain

    @staticmethod
    def _runtime_domain_symbol(name: str, domain: str):
        if domain == "complex":
            return S(f"pyamplicol::runtime_domain::{name}")
        real_symbol = S(
            f"pyamplicol::runtime_domain::{name}",
            is_real=True,
        )
        return 1j * real_symbol if domain == "imaginary" else real_symbol

    def auxiliary_kind(self, particle_id: int) -> str | None:
        return self._particle_records_by_pdg[int(particle_id)].auxiliary_kind

    def runtime_derived_parameter_defaults(self) -> dict[str, complex]:
        return self.runtime_derived_parameter_defaults_for(
            tuple(self._runtime_derived_definitions)
        )

    def runtime_derived_parameter_defaults_for(
        self,
        names: Sequence[str],
    ) -> dict[str, complex]:
        substitutions: dict[Expression, Any] | None = None
        defaults: dict[str, complex] = {}
        for name in names:
            cached = self._runtime_derived_default_cache.get(name)
            if cached is not None:
                defaults[name] = cached
                continue
            expression = self._runtime_derived_definitions.get(name)
            if expression is None:
                continue
            template = self._runtime_derived_expression_cache.get(name)
            if template is None:
                template = E(expression)
                self._runtime_derived_expression_cache[name] = template
            if substitutions is None:
                substitutions = {
                    self._expression_symbol(f"UFO::{parameter_name}"): value
                    for parameter_name, value in self._runtime_parameters.items()
                }
            try:
                value = complex(template.evaluate(substitutions))
            except Exception:  # noqa: BLE001 - retain symbolic compatibility.
                value = complex(_replace_symbols(template, substitutions))
            self._runtime_derived_default_cache[name] = value
            defaults[name] = value
        return defaults

    def runtime_parameter_names_for_particle(self, pdg: int) -> tuple[str, ...]:
        particle = self._particle_records_by_pdg[int(pdg)]
        custom = self._propagator_record(pdg)
        if custom is not None and custom.custom:
            _expression, names = self._custom_propagator_expression(int(pdg))
            return names
        return tuple(
            name
            for name in (particle.mass, particle.width)
            if name.upper() != "ZERO" and name in self.compiled.parameter_defaults
        )

    def _custom_propagator_expression(
        self,
        particle_id: int,
    ) -> tuple[Expression, tuple[str, ...]]:
        cached = self._custom_propagator_expressions.get(int(particle_id))
        if cached is not None:
            return cached
        record = self._propagator_record(particle_id)
        if record is None or not record.custom:
            raise ValueError(f"particle {particle_id} has no custom UFO propagator")
        spin = self.spin(particle_id)
        normalized = normalize_lorentz_expression(
            f"({record.numerator})/({record.denominator})",
            (spin, spin),
        )
        expression = E(normalized.expression)
        for parameter in self._parameter_records.values():
            if parameter.nature == "external":
                continue
            expression = expression.replace(
                S(f"UFO::{parameter.name}"),
                E(parameter.resolved_expression),
            )
        expression = _replace_evaluator_constants(expression)
        symbols = set(expression.get_all_symbols(False))
        names = tuple(
            sorted(
                parameter.name
                for parameter in self._parameter_records.values()
                if parameter.nature == "external"
                and S(f"UFO::{parameter.name}") in symbols
            )
        )
        result = expression, names
        self._custom_propagator_expressions[int(particle_id)] = result
        return result

    def _custom_propagator_template(
        self,
        particle_id: int,
    ) -> tuple[Expression, ...]:
        cached = self._custom_propagator_templates.get(int(particle_id))
        if cached is not None:
            return cached
        expression, _names = self._custom_propagator_expression(particle_id)
        spin = self.spin(particle_id)
        components = tuple(
            S(f"pyamplicol::custom_propagator_input_{index}")
            for index in range(self.dimension(particle_id))
        )
        momenta = tuple(
            S(f"pyamplicol::custom_propagator_momentum_{index}")
            for index in range(4)
        )
        library = TensorLibrary.hep_lib_atom()
        representations = _spin_representations(spin)
        if representations:
            name = TensorName("pyamplicol::custom_propagator_input")
            library.register(LibraryTensor.dense(name(*representations), components))
            expression *= name(*_spin_slots(spin, 1)).to_expression()
        else:
            expression *= components[0]
        minkowski = Representation.mink(4)
        for leg in (1, 2):
            library.register(
                LibraryTensor.dense(
                    TensorName(f"pyamplicol::ufo_momentum_{leg}")(minkowski),
                    momenta,
                )
            )
        network = TensorNetwork(expression, library)
        network.execute(library=library)
        tensor = network.result_tensor(library)
        tensor.to_dense()
        expected = self.dimension(particle_id)
        if len(tensor) != expected:
            raise ValueError(
                f"custom propagator {particle_id} produced {len(tensor)} "
                f"components, expected {expected}"
            )
        template = tuple(
            _replace_evaluator_constants(_as_expression(tensor[index]))
            for index in range(len(tensor))
        )
        self._custom_propagator_templates[int(particle_id)] = template
        return template

    def _evaluate_custom_propagator(
        self,
        particle_id: int,
        value: Sequence[Any],
        momentum: Sequence[Any],
    ) -> tuple[Any, ...]:
        if len(value) != self.dimension(particle_id):
            raise ValueError(
                f"custom propagator for {particle_id} expects "
                f"{self.dimension(particle_id)} current components"
            )
        if len(momentum) != 4:
            raise ValueError("UFO custom propagators require four-momentum components")
        substitutions = {
            **{
                S(f"pyamplicol::custom_propagator_input_{index}"): component
                for index, component in enumerate(value)
            },
            **{
                S(f"pyamplicol::custom_propagator_momentum_{index}"): component
                for index, component in enumerate(momentum)
            },
            **{
                S(f"UFO::{name}"): self._parameter_value(name)
                for name in self._custom_propagator_expression(particle_id)[1]
            },
        }
        return tuple(
            _replace_symbols(component, substitutions)
            for component in self._custom_propagator_template(particle_id)
        )

    def runtime_mass_parameter_name(self, pdg: int) -> str | None:
        name = self._particle_records_by_pdg[int(pdg)].mass
        return name if name in self.compiled.parameter_defaults else None

    def runtime_width_parameter_name(self, pdg: int) -> str | None:
        name = self._particle_records_by_pdg[int(pdg)].width
        return name if name in self.compiled.parameter_defaults else None

    def runtime_parameter_defaults(self) -> dict[str, tuple[float, float]]:
        return dict(self.compiled.parameter_defaults)

    def runtime_parameter_type(self, name: str) -> str:
        try:
            return self._parameter_records[str(name)].parameter_type
        except KeyError as exc:
            raise KeyError(f"unknown runtime model parameter {name!r}") from exc

    def runtime_normalization_payload(self, dag: Any) -> dict[str, object]:
        initial_pdgs = tuple(int(pdg) for pdg in dag.process.initial_pdgs)
        final_pdgs = tuple(int(pdg) for pdg in dag.process.final_pdgs)
        average_factor = 1
        for pdg in initial_pdgs:
            average_factor *= len(self.source_spin_states(pdg))
            average_factor *= max(1, abs(self.color_rep(pdg)))
        identical_factor = math.prod(
            math.factorial(count) for count in Counter(final_pdgs).values()
        )
        color_factor = self.leading_color_factor((*initial_pdgs, *final_pdgs))
        return {
            "color_accuracy": dag.process.color_accuracy,
            "color_factor": float(color_factor),
            "average_factor": float(average_factor),
            "identical_factor": float(identical_factor),
            "final_state_identical_factor": float(identical_factor),
            "quark_line_partner_factor": 1,
            "global_coupling_factor": 1.0,
            "qcd_coupling_power": 0,
            "electroweak_coupling_power": 0,
            "couplings_in_stage_evaluators": True,
            "coupling_policy": (
                "UFO coupling expressions are fully included in generated stage "
                "evaluators; no built-in AmpliCol global coupling factor is applied"
            ),
        }

    def leading_color_factor(self, process: Sequence[int]) -> int:
        exponent_twice = 0
        for pdg in process:
            representation = abs(self.color_rep(int(pdg)))
            if representation == 8:
                exponent_twice += 2
            elif representation == 3:
                exponent_twice += 1
            elif representation != 1:
                raise ValueError(
                    f"unsupported leading-color representation {representation}"
                )
        if exponent_twice % 2:
            raise ValueError(
                f"non-integer leading-color exponent for {tuple(process)}"
            )
        return 3 ** (exponent_twice // 2)

    def _kernel(self, kind: int) -> CompiledOrientedKernel:
        try:
            return self._kernels[int(kind)]
        except KeyError as exc:
            raise KeyError(f"unknown compiled UFO kernel kind {kind}") from exc

    def _weyl_projection_is_nonzero(
        self,
        kind: int,
        left_chirality: int,
        right_chirality: int,
        result_chirality: int,
    ) -> bool:
        key = (kind, left_chirality, right_chirality, result_chirality)
        cached = self._weyl_projection_support.get(key)
        if cached is not None:
            return cached
        kernel = self._kernel(kind)
        left_pdg = self._particle_records_by_name[kernel.particles[0]].pdg_code
        right_pdg = self._particle_records_by_name[kernel.particles[1]].pdg_code
        left = tuple(
            S(f"pyamplicol::weyl_probe_{kind}_left_{index}")
            for index in range(self.current_dimension(left_pdg, left_chirality))
        )
        right = tuple(
            S(f"pyamplicol::weyl_probe_{kind}_right_{index}")
            for index in range(self.current_dimension(right_pdg, right_chirality))
        )
        components = self._projected_kernel_components(
            kernel,
            left,
            right,
            left_chirality=left_chirality,
            right_chirality=right_chirality,
            result_chirality=result_chirality,
            left_momentum=tuple(
                S(f"pyamplicol::weyl_probe_{kind}_left_momentum_{index}")
                for index in range(4)
            ),
            right_momentum=tuple(
                S(f"pyamplicol::weyl_probe_{kind}_right_momentum_{index}")
                for index in range(4)
            ),
        )
        supported = any(not _is_zero(component) for component in components)
        self._weyl_projection_support[key] = supported
        return supported

    def _projected_kernel_components(
        self,
        kernel: CompiledOrientedKernel,
        left: Sequence[Any],
        right: Sequence[Any],
        *,
        left_chirality: int,
        right_chirality: int,
        result_chirality: int,
        left_momentum: Sequence[Any],
        right_momentum: Sequence[Any],
        runtime_parameter_values: Mapping[str, Any] | None = None,
    ) -> tuple[Any, ...]:
        left_pdg = self._particle_records_by_name[kernel.particles[0]].pdg_code
        right_pdg = self._particle_records_by_name[kernel.particles[1]].pdg_code
        result_pdg = self._particle_records_by_name[kernel.particles[2]].pdg_code
        full_left = self._embed_weyl_current(left_pdg, left_chirality, left)
        full_right = self._embed_weyl_current(right_pdg, right_chirality, right)
        parameter_values = (
            {
                name: self._parameter_value(name)
                for name in kernel.runtime_parameters
            }
            if runtime_parameter_values is None
            else dict(runtime_parameter_values)
        )
        substitutions: dict[Expression, Any] = {}
        for index, value in enumerate(full_left):
            substitutions[
                self._expression_symbol(
                    f"pyamplicol::kernel_{kernel.kind}_left_{index}"
                )
            ] = value
        for index, value in enumerate(full_right):
            substitutions[
                self._expression_symbol(
                    f"pyamplicol::kernel_{kernel.kind}_right_{index}"
                )
            ] = value
        for index, value in enumerate(left_momentum):
            substitutions[
                self._expression_symbol(
                    f"pyamplicol::kernel_{kernel.kind}_left_momentum_{index}"
                )
            ] = value
        for index, value in enumerate(right_momentum):
            substitutions[
                self._expression_symbol(
                    f"pyamplicol::kernel_{kernel.kind}_right_momentum_{index}"
                )
            ] = value
        for name, value in parameter_values.items():
            substitutions[self._expression_symbol(f"UFO::{name}")] = value
        templates = self._kernel_component_expressions(kernel)
        if all(_is_numeric(value) for value in substitutions.values()):
            try:
                components = tuple(
                    complex(template.evaluate(substitutions))
                    for template in templates
                )
            except Exception:  # noqa: BLE001 - retain symbolic compatibility.
                components = tuple(
                    _replace_symbols(template, substitutions)
                    for template in templates
                )
        else:
            components = tuple(
                _replace_symbols(template, substitutions)
                for template in templates
            )
        if not self.is_chiral_eligible(result_pdg):
            return components
        if result_chirality == 1:
            return components[2:4]
        if result_chirality == -1:
            return components[0:2]
        raise ValueError("a projected Weyl result requires nonzero chirality")

    def _kernel_function_component_calls(
        self,
        kernel: CompiledOrientedKernel,
        left: Sequence[Any],
        right: Sequence[Any],
        *,
        left_chirality: int,
        right_chirality: int,
        result_chirality: int,
        left_momentum: Sequence[Any],
        right_momentum: Sequence[Any],
        runtime_parameter_values: Mapping[str, Any],
    ) -> tuple[Any, ...]:
        key = (
            int(kernel.kind),
            int(left_chirality),
            int(right_chirality),
            int(result_chirality),
        )
        function_specs = self._kernel_function_specs.get(key)
        if function_specs is None:
            function_specs = self._define_kernel_component_functions(
                kernel,
                left_dimension=len(left),
                right_dimension=len(right),
                left_chirality=left_chirality,
                right_chirality=right_chirality,
                result_chirality=result_chirality,
            )
            self._kernel_function_specs[key] = function_specs
        arguments = (
            *left,
            *right,
            *left_momentum,
            *right_momentum,
            *(runtime_parameter_values[name] for name in kernel.runtime_parameters),
        )
        return tuple(
            function(*(arguments[index] for index in argument_indices))
            for function, argument_indices in function_specs
        )

    def _define_kernel_component_functions(
        self,
        kernel: CompiledOrientedKernel,
        *,
        left_dimension: int,
        right_dimension: int,
        left_chirality: int,
        right_chirality: int,
        result_chirality: int,
    ) -> tuple[tuple[Expression, tuple[int, ...]], ...]:
        key_tag = "_".join(
            (
                str(int(kernel.kind)),
                _chirality_tag(left_chirality),
                _chirality_tag(right_chirality),
                _chirality_tag(result_chirality),
            )
        )
        formal_left = tuple(
            S(f"pyamplicol::kernel_function_{key_tag}_left_{index}")
            for index in range(left_dimension)
        )
        formal_right = tuple(
            S(f"pyamplicol::kernel_function_{key_tag}_right_{index}")
            for index in range(right_dimension)
        )
        formal_left_momentum = tuple(
            S(f"pyamplicol::kernel_function_{key_tag}_left_momentum_{index}")
            for index in range(4)
        )
        formal_right_momentum = tuple(
            S(f"pyamplicol::kernel_function_{key_tag}_right_momentum_{index}")
            for index in range(4)
        )
        formal_parameters = {
            name: S(f"pyamplicol::kernel_function_{key_tag}_parameter_{index}")
            for index, name in enumerate(kernel.runtime_parameters)
        }
        formal_arguments = (
            *formal_left,
            *formal_right,
            *formal_left_momentum,
            *formal_right_momentum,
            *(formal_parameters[name] for name in kernel.runtime_parameters),
        )
        components = self._projected_kernel_components(
            kernel,
            formal_left,
            formal_right,
            left_chirality=left_chirality,
            right_chirality=right_chirality,
            result_chirality=result_chirality,
            left_momentum=formal_left_momentum,
            right_momentum=formal_right_momentum,
            runtime_parameter_values=formal_parameters,
        )
        coupling_expression = self._resolved_kernel_coupling_expression(
            kernel,
            formal_parameters,
        )
        bodies = tuple(component * coupling_expression for component in components)
        functions = tuple(
            S(f"pyamplicol::kernel_function_{key_tag}_component_{index}")
            for index in range(len(bodies))
        )
        specs: list[tuple[Expression, tuple[int, ...]]] = []
        for function, body in zip(functions, bodies, strict=True):
            body_symbols = set(body.get_all_symbols(False))
            argument_indices = tuple(
                index
                for index, argument in enumerate(formal_arguments)
                if argument in body_symbols
            )
            used_arguments = tuple(
                formal_arguments[index] for index in argument_indices
            )
            self._symbolica_kernel_functions[(function, used_arguments)] = body
            specs.append((function, argument_indices))
        return tuple(specs)

    def _kernel_component_expressions(
        self,
        kernel: CompiledOrientedKernel,
    ) -> tuple[Expression, ...]:
        cached = self._kernel_component_expression_cache.get(kernel.kind)
        if cached is None:
            cached = tuple(E(component) for component in kernel.component_expressions)
            self._kernel_component_expression_cache[kernel.kind] = cached
        return cached

    def _resolved_kernel_coupling_expression(
        self,
        kernel: CompiledOrientedKernel,
        runtime_parameter_values: Mapping[str, Any],
    ) -> Any:
        template = self._kernel_coupling_expression_cache.get(kernel.kind)
        if template is None:
            template = E(kernel.coupling_expression)
            self._kernel_coupling_expression_cache[kernel.kind] = template
        substitutions = {
            self._expression_symbol(f"UFO::{name}"): value
            for name, value in runtime_parameter_values.items()
        }
        if all(_is_numeric(value) for value in substitutions.values()):
            try:
                return complex(template.evaluate(substitutions))
            except Exception:  # noqa: BLE001 - retain symbolic compatibility.
                pass
        return _replace_symbols(template, substitutions)

    def _expression_symbol(self, name: str) -> Expression:
        symbol = self._expression_symbol_cache.get(name)
        if symbol is None:
            symbol = S(name)
            self._expression_symbol_cache[name] = symbol
        return symbol

    def _embed_weyl_current(
        self,
        particle_id: int,
        chirality: int,
        values: Sequence[Any],
    ) -> tuple[Any, ...]:
        components = tuple(values)
        if not self.is_chiral_eligible(particle_id):
            return components
        if len(components) != 2:
            raise ValueError("a projected Weyl input must have two components")
        if chirality == 1:
            return (components[0], components[1], 0.0, 0.0)
        if chirality == -1:
            return (0.0, 0.0, components[0], components[1])
        raise ValueError("a projected Weyl input requires nonzero chirality")

    def _propagator_record(self, particle_id: int):
        name = self._particle_records_by_pdg[int(particle_id)].name
        return self._propagator_records_by_particle_name.get(name)

    def _goldstone_is_redundant_in_unitary_gauge(self, goldstone: Any) -> bool:
        """Return whether a unitary-gauge vector already carries this mode."""

        for vector in self.compiled.ir.particles:
            if (
                vector.spin != 3
                or vector.goldstoneboson
                or not vector.propagating
                or vector.mass != goldstone.mass
                or vector.color != goldstone.color
                or not math.isclose(
                    vector.charge,
                    goldstone.charge,
                    rel_tol=0.0,
                    abs_tol=1.0e-12,
                )
            ):
                continue
            propagator = self._propagator_records_by_particle_name.get(vector.name)
            # pyAmpliCol synthesizes the canonical unitary-gauge projector for
            # a massive vector unless the UFO supplied a genuinely custom
            # propagator. In that case a separate Goldstone would double count
            # the vector's longitudinal mode.
            if propagator is None or not propagator.custom:
                return True
        return False

    def _parameter_default(self, name: str) -> complex:
        if name.upper() == "ZERO":
            return 0.0 + 0.0j
        record = self._parameter_records.get(name)
        if record is None:
            raise ValueError(f"model parameter {name!r} is not defined")
        return _record_default(record)

    def _parameter_value(self, name: str) -> Any:
        if name in self._runtime_parameters:
            return self._runtime_parameters[name]
        definition = self._runtime_derived_definitions.get(name)
        if definition is not None:
            return _replace_symbols(
                E(definition),
                {
                    S(f"UFO::{parameter_name}"): value
                    for parameter_name, value in self._runtime_parameters.items()
                    if not parameter_name.startswith("derived_coupling_")
                },
            )
        return self._parameter_default(name)

    def _real_parameter_value(self, name: str, *, field: str) -> Any:
        value = self._parameter_value(name)
        if isinstance(value, int | float | complex):
            numeric = complex(value)
            if numeric.imag != 0.0:
                raise ValueError(f"particle {field} parameter {name!r} is not real")
            return numeric.real
        return value


def _record_default(record: CompiledParameterRecord) -> complex:
    if record.value is not None:
        return complex(record.value[0], record.value[1])
    return complex(E(record.resolved_expression).evaluate({}))


def _spin_dimension(spin: int) -> int:
    try:
        return {-1: 1, 1: 1, 2: 4, 3: 4, 5: 16}[int(spin)]
    except KeyError as exc:
        raise ValueError(f"unsupported UFO spin code {spin}") from exc


def _chirality_tag(chirality: int) -> str:
    value = int(chirality)
    if value < 0:
        return f"m{abs(value)}"
    if value > 0:
        return f"p{value}"
    return "z"


def _replace_symbols(expression: Expression, substitutions: Mapping[Expression, Any]) -> Any:
    symbols = set(expression.get_all_symbols(False))
    replacements = [
        Replacement(symbol, value)
        for symbol, value in substitutions.items()
        if symbol in symbols
    ]
    if not replacements:
        return expression
    return expression.replace_multiple(replacements)


def _is_numeric(value: Any) -> bool:
    return isinstance(value, Number)


def _is_zero(value: Any) -> bool:
    if isinstance(value, int | float | complex):
        return value == 0
    return isinstance(value, Expression) and value.to_canonical_string() == "0"


def _expr_spin2_propagator(
    value: Sequence[Any],
    momentum: Sequence[Any],
    mass: Any,
    width: Any,
    *,
    dimension: Any,
    massive: bool,
) -> tuple[Any, ...]:
    if len(value) != 16:
        raise ValueError("spin-2 propagator expects sixteen current components")
    if len(momentum) != 4:
        raise ValueError("spin-2 propagator expects four momentum components")
    tensor = tuple(tuple(value[4 * mu + nu] for nu in range(4)) for mu in range(4))
    metric = (1.0, -1.0, -1.0, -1.0)
    denominator = (
        _minkowski_square_expression(momentum)
        - mass * mass
        + 1j * mass * width
    )
    if not massive:
        trace = sum(
            (metric[index] * tensor[index][index] for index in range(4)),
            0.0,
        )
        trace_weight = 1.0 / (dimension - 2.0)
        projected = tuple(
            0.5 * (tensor[mu][nu] + tensor[nu][mu])
            - (
                metric[mu] * trace * trace_weight
                if mu == nu
                else 0.0
            )
            for mu in range(4)
            for nu in range(4)
        )
        return tuple(1j * component / denominator for component in projected)

    mass_squared = mass * mass
    first_projected = tuple(
        tuple(
            tensor[mu][nu]
            - momentum[mu]
            * sum(
                (
                    metric[alpha]
                    * momentum[alpha]
                    * tensor[alpha][nu]
                    for alpha in range(4)
                ),
                0.0,
            )
            / mass_squared
            for nu in range(4)
        )
        for mu in range(4)
    )
    transverse = tuple(
        tuple(
            first_projected[mu][nu]
            - momentum[nu]
            * sum(
                (
                    metric[beta]
                    * momentum[beta]
                    * first_projected[mu][beta]
                    for beta in range(4)
                ),
                0.0,
            )
            / mass_squared
            for nu in range(4)
        )
        for mu in range(4)
    )
    transverse_trace = sum(
        (metric[index] * transverse[index][index] for index in range(4)),
        0.0,
    )
    projected = tuple(
        0.5 * (transverse[mu][nu] + transverse[nu][mu])
        - (
            (metric[mu] if mu == nu else 0.0)
            - momentum[mu] * momentum[nu] / mass_squared
        )
        * transverse_trace
        / 3.0
        for mu in range(4)
        for nu in range(4)
    )
    return tuple(1j * component / denominator for component in projected)


__all__ = ["CompiledUFOModel"]
