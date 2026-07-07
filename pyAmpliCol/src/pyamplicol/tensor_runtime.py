from __future__ import annotations

import math
import time
from dataclasses import dataclass
from itertools import product
from typing import Any, Iterable, Literal, Sequence

from .lowering import (
    _GraphTensorExpressionBuilder,
    _current_dimension,
    _current_key_tuple,
    _current_momentum_currents,
    _current_momentum_parameter_head,
    _current_momentum_tensor_name,
    _current_parameter_head,
    _current_tensor_name,
    _is_weyl_fermion_current,
    _propagating_currents,
    _propagator_tensor_name,
    _source_currents,
    build_interleaved_tensor_network_scalar_bundle,
    build_tensor_network_scalar_bundle,
)
from .legacy_matrix import CurrentKey, NativeMatrixElementGenerator, RecursionGraph
from .model import AmplicolSMLeadingColorModel
from .native import (
    ExternalMomentum,
    HelicityContribution,
    LeadingColorZJetsNativeEvaluator,
    MatrixElementEvaluation,
    NativeEvaluationError,
    _ext_antiquark_weyl,
    _ext_gluon_cmplx,
    _ext_massive_vector,
    _ext_quark_weyl,
    _final_state_identical_factor,
    _initial_state_average_factor,
)
from .params import SymbolicaEvaluatorBundle

TensorNetworkStrategy = Literal["interleaved", "monolithic"]


@dataclass(frozen=True)
class TensorNetworkRuntimeMetadata:
    process: str
    kernel: str
    strategy: str
    gluon_count: int
    parameter_count: int
    graph_current_count: int
    graph_interaction_count: int
    graph_amplitude_count: int
    helicity_filter_original_count: int | None = None
    helicity_filter_kept_count: int | None = None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "kernel": self.kernel,
            "strategy": self.strategy,
            "gluon_count": self.gluon_count,
            "parameter_count": self.parameter_count,
            "graph_current_count": self.graph_current_count,
            "graph_interaction_count": self.graph_interaction_count,
            "graph_amplitude_count": self.graph_amplitude_count,
            "helicity_filter_original_count": self.helicity_filter_original_count,
            "helicity_filter_kept_count": self.helicity_filter_kept_count,
        }


@dataclass(frozen=True)
class HelicityFilterEntry:
    helicities: tuple[int, ...]
    multiplicity: int
    equivalents: tuple[tuple[int, ...], ...]

    def to_json_dict(self) -> dict[str, object]:
        return {
            "helicities": list(self.helicities),
            "multiplicity": self.multiplicity,
            "equivalents": [list(helicities) for helicities in self.equivalents],
        }

    @classmethod
    def from_json_dict(cls, payload: dict[str, object]) -> HelicityFilterEntry:
        helicities = _int_tuple(payload.get("helicities"))
        multiplicity = _int_payload_value(payload.get("multiplicity", 1))
        equivalents_payload = payload.get("equivalents", [])
        if not isinstance(equivalents_payload, list):
            raise ValueError("helicity filter entry equivalents must be a list")
        equivalents = tuple(_int_tuple(item) for item in equivalents_payload)
        return cls(
            helicities=helicities,
            multiplicity=multiplicity,
            equivalents=equivalents,
        )


@dataclass(frozen=True)
class HelicityFilter:
    entries: tuple[HelicityFilterEntry, ...]
    original_count: int
    sample_count: int
    relative_tolerance: float
    zero_tolerance: float
    sample_sqrt_s: tuple[float, ...]
    phase_space_mode: str = "canonical"
    seed: int | None = None

    @property
    def kept_count(self) -> int:
        return len(self.entries)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": "pyamplicol-z-gluon-helicity-filter",
            "entries": [entry.to_json_dict() for entry in self.entries],
            "original_count": self.original_count,
            "kept_count": self.kept_count,
            "sample_count": self.sample_count,
            "relative_tolerance": self.relative_tolerance,
            "zero_tolerance": self.zero_tolerance,
            "sample_sqrt_s": list(self.sample_sqrt_s),
            "phase_space_mode": self.phase_space_mode,
            "seed": self.seed,
        }

    @classmethod
    def from_json_dict(cls, payload: dict[str, object]) -> HelicityFilter:
        entries_payload = payload.get("entries")
        if not isinstance(entries_payload, list):
            raise ValueError("helicity filter payload is missing entries")
        sample_sqrt_s = payload.get("sample_sqrt_s", [])
        if not isinstance(sample_sqrt_s, list):
            raise ValueError("helicity filter sample_sqrt_s must be a list")
        return cls(
            entries=tuple(
                HelicityFilterEntry.from_json_dict(entry)
                for entry in entries_payload
                if isinstance(entry, dict)
            ),
            original_count=_int_payload_value(payload.get("original_count", 0)),
            sample_count=_int_payload_value(payload.get("sample_count", 0)),
            relative_tolerance=_float_payload_value(
                payload.get("relative_tolerance", 1.0e-12)
            ),
            zero_tolerance=_float_payload_value(
                payload.get("zero_tolerance", 1.0e-300)
            ),
            sample_sqrt_s=tuple(float(value) for value in sample_sqrt_s),
            phase_space_mode=str(payload.get("phase_space_mode", "canonical")),
            seed=(
                None
                if payload.get("seed") is None
                else _int_payload_value(payload.get("seed"))
            ),
        )


@dataclass(frozen=True)
class NumericTensorNetworkRuntimeMetadata:
    process: str
    kernel: str
    strategy: str
    gluon_count: int
    graph_current_count: int
    graph_interaction_count: int
    graph_amplitude_count: int
    expression_build_time_s: float

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "kernel": self.kernel,
            "strategy": self.strategy,
            "gluon_count": self.gluon_count,
            "graph_current_count": self.graph_current_count,
            "graph_interaction_count": self.graph_interaction_count,
            "graph_amplitude_count": self.graph_amplitude_count,
            "expression_build_time_s": self.expression_build_time_s,
        }


class ZGluonTensorNetworkEvaluator:
    """Propagated tensor-network scalar evaluator for q q~ -> Z + ordered gluons."""

    def __init__(
        self,
        process: str,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        graph: RecursionGraph | None = None,
        bundle: SymbolicaEvaluatorBundle | None = None,
        strategy: TensorNetworkStrategy = "interleaved",
        helicity_filter: HelicityFilter | None = None,
        build_helicity_filter: bool = False,
    ) -> None:
        self.process = process
        self.model = model or AmplicolSMLeadingColorModel()
        self.graph = graph or _z_gluon_graph(process, self.model)
        self.gluon_count = _validate_z_gluon_graph(self.graph)
        self.strategy: TensorNetworkStrategy = _validate_tensor_network_strategy(strategy)
        self.helicity_filter = helicity_filter
        if bundle is not None:
            self.bundle = bundle
        elif self.strategy == "interleaved":
            self.bundle = build_interleaved_tensor_network_scalar_bundle(
                self.model,
                self.graph,
                name=_tensor_network_bundle_name(process),
            )
        else:
            self.bundle = build_tensor_network_scalar_bundle(
                self.model,
                self.graph,
                name=_tensor_network_bundle_name(process),
            )
        if build_helicity_filter and self.helicity_filter is None:
            self.helicity_filter = build_z_gluon_helicity_filter(self)

    @property
    def metadata(self) -> TensorNetworkRuntimeMetadata:
        return TensorNetworkRuntimeMetadata(
            process=self.process,
            kernel=_tensor_network_kernel(self.gluon_count),
            strategy=self.strategy,
            gluon_count=self.gluon_count,
            parameter_count=len(self.bundle.param_builder.values),
            graph_current_count=len(self.graph.currents),
            graph_interaction_count=len(self.graph.interactions),
            graph_amplitude_count=len(self.graph.amplitudes),
            helicity_filter_original_count=(
                None if self.helicity_filter is None else self.helicity_filter.original_count
            ),
            helicity_filter_kept_count=(
                None if self.helicity_filter is None else self.helicity_filter.kept_count
            ),
        )

    def to_artifact_payload(self) -> dict[str, object]:
        payload = self.bundle.to_artifact_payload()
        payload["kind"] = "pyamplicol-z-gluon-tensor-network-evaluator"
        payload["runtime_metadata"] = self.metadata.to_json_dict()
        payload["helicity_filter"] = (
            None if self.helicity_filter is None else self.helicity_filter.to_json_dict()
        )
        return payload

    @classmethod
    def from_artifact_payload(
        cls,
        process: str,
        payload: dict[str, object],
        *,
        model: AmplicolSMLeadingColorModel | None = None,
    ) -> ZGluonTensorNetworkEvaluator:
        return cls(
            process,
            model=model,
            bundle=SymbolicaEvaluatorBundle.from_artifact_payload(payload),
            strategy=_strategy_from_bundle_payload(payload),
            helicity_filter=_helicity_filter_from_payload(payload),
        )

    def evaluate(
        self,
        particles: Sequence[ExternalMomentum],
    ) -> MatrixElementEvaluation:
        point = _validate_z_gluon_point(particles, gluon_count=self.gluon_count)
        raw_sum = 0.0
        helicity_contributions: list[HelicityContribution] = []
        for helicities, amplitude, multiplicity in self._weighted_helicity_amplitudes(point):
            squared = float((amplitude * amplitude.conjugate()).real)
            weighted_squared = squared * multiplicity
            raw_sum += weighted_squared
            helicity_contributions.append(
                HelicityContribution(
                    helicities=helicities,
                    amplitude=amplitude,
                    squared=weighted_squared,
                )
            )

        pdgs = tuple(particle.pdg for particle in point)
        color_factor = self.model.leading_color_factor(pdgs)
        average_factor = _initial_state_average_factor(pdgs[:2])
        identical_factor = _final_state_identical_factor(
            particle.pdg for particle in point[2:]
        )
        coupling_factor = (
            (4.0 * math.pi * self.model.alpha_s_me_check) ** self.gluon_count
            * (2.0 * 4.0 * math.pi * self.model.alpha_ew)
        )
        matrix_element = (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
        )
        return MatrixElementEvaluation(
            process=self.process,
            particles=point,
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=color_factor,
            average_factor=average_factor,
            coupling_factor=coupling_factor,
            helicity_contributions=tuple(helicity_contributions),
            identical_factor=identical_factor,
        )

    def _weighted_helicity_amplitudes(
        self,
        particles: tuple[ExternalMomentum, ...],
    ) -> tuple[tuple[tuple[int, ...], complex, int], ...]:
        if self.helicity_filter is None:
            entries = tuple(
                HelicityFilterEntry(
                    helicities=helicities,
                    multiplicity=1,
                    equivalents=(helicities,),
                )
                for helicities in _all_z_gluon_helicities(self.gluon_count)
            )
        else:
            entries = self.helicity_filter.entries
        return self._evaluate_helicity_entries(particles, entries)

    def _helicity_amplitudes(
        self,
        particles: tuple[ExternalMomentum, ...],
    ) -> tuple[tuple[tuple[int, ...], complex], ...]:
        return tuple(
            (helicities, amplitude)
            for helicities, amplitude, _ in self._evaluate_helicity_entries(
                particles,
                tuple(
                    HelicityFilterEntry(
                        helicities=helicities,
                        multiplicity=1,
                        equivalents=(helicities,),
                    )
                    for helicities in _all_z_gluon_helicities(self.gluon_count)
                ),
            )
        )

    def _evaluate_helicity_entries(
        self,
        particles: tuple[ExternalMomentum, ...],
        entries: tuple[HelicityFilterEntry, ...],
    ) -> tuple[tuple[tuple[int, ...], complex, int], ...]:
        incoming_quark = particles[0]
        incoming_antiquark = particles[1]
        gluons = particles[2 : 2 + self.gluon_count]
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
            helicity: _ext_massive_vector(z_momentum, helicity, self.model.mass(23))
            for helicity in (-1, 0, 1)
        }
        gluon_vectors = tuple(
            {
                helicity: _ext_gluon_cmplx(momentum, helicity)
                for helicity in (-1, 1)
            }
            for momentum in gluon_momenta
        )
        self._populate_momentum_parameters(particles)
        builder = self.bundle.param_builder
        anti_head = _current_parameter_head(
            CurrentKey(-incoming_quark.pdg, (1,), 0)
        )
        quark_heads = {
            chirality: _current_parameter_head(
                CurrentKey(-incoming_antiquark.pdg, (2,), chirality)
            )
            for chirality in (-1, 1)
        }
        gluon_heads = tuple(
            _current_parameter_head(CurrentKey(21, (offset,), 0))
            for offset in range(3, 3 + self.gluon_count)
        )
        z_head = _current_parameter_head(
            CurrentKey(23, (self.gluon_count + 3,), 0)
        )
        param_values = builder.values
        evaluator = self.bundle.evaluator
        fermion_states = {
            (1, -1): (1, quark_minus, anti_plus),
            (-1, 1): (-1, quark_plus, anti_minus),
        }
        amplitudes: list[tuple[tuple[int, ...], complex, int]] = []
        active_fermion_helicities: tuple[int, int] | None = None
        active_z_helicity: int | None = None
        for entry in entries:
            helicities = entry.helicities
            if len(helicities) != self.gluon_count + 3:
                raise NativeEvaluationError(
                    "helicity filter entry has incompatible helicity length"
                )
            fermion_helicities = (helicities[0], helicities[1])
            if fermion_helicities not in fermion_states:
                raise NativeEvaluationError(
                    f"unsupported q q~ helicity pair in filter: {fermion_helicities}"
                )
            chirality, quark_wf, anti_wf = fermion_states[fermion_helicities]
            if active_fermion_helicities != fermion_helicities:
                active_fermion_helicities = fermion_helicities
                active_z_helicity = None
                _set_parameter_values_unchecked(builder, anti_head, anti_wf)
                _set_parameter_values_unchecked(
                    builder,
                    quark_heads[chirality],
                    quark_wf,
                )
                _set_parameter_values_unchecked(
                    builder,
                    quark_heads[-chirality],
                    (0j, 0j),
                )
            z_helicity = helicities[-1]
            if z_helicity not in z_vectors:
                raise NativeEvaluationError(
                    f"unsupported Z helicity in filter: {z_helicity}"
                )
            if active_z_helicity != z_helicity:
                active_z_helicity = z_helicity
                _set_parameter_values_unchecked(builder, z_head, z_vectors[z_helicity])
            gluon_helicities = helicities[2:-1]
            for head, helicity, vectors in zip(
                gluon_heads,
                gluon_helicities,
                gluon_vectors,
                strict=True,
            ):
                _set_parameter_values_unchecked(
                    builder,
                    head,
                    vectors[helicity],
                )
            amplitude = complex(evaluator.evaluate_complex([param_values])[0][0])
            amplitudes.append((helicities, amplitude, entry.multiplicity))
        return tuple(amplitudes)

    def _populate_parameters(
        self,
        particles: tuple[ExternalMomentum, ...],
        *,
        chirality: int,
        quark_wf: tuple[complex, complex],
        anti_wf: tuple[complex, complex],
        gluon_wfs: tuple[tuple[complex, complex, complex, complex], ...],
        z_wf: tuple[complex, complex, complex, complex],
    ) -> None:
        self._populate_fermion_source_parameters(
            particles,
            chirality=chirality,
            quark_wf=quark_wf,
            anti_wf=anti_wf,
        )
        self._populate_gluon_source_values(gluon_wfs)
        self._populate_z_source_parameter(z_wf=z_wf)
        self._populate_momentum_parameters(particles)

    def _populate_fermion_source_parameters(
        self,
        particles: tuple[ExternalMomentum, ...],
        *,
        chirality: int,
        quark_wf: tuple[complex, complex],
        anti_wf: tuple[complex, complex],
    ) -> None:
        builder = self.bundle.param_builder
        incoming_quark = particles[0]
        incoming_antiquark = particles[1]
        anti_head = _current_parameter_head(
            CurrentKey(-incoming_quark.pdg, (1,), 0)
        )
        active_quark_head = _current_parameter_head(
            CurrentKey(-incoming_antiquark.pdg, (2,), chirality)
        )
        inactive_quark_head = _current_parameter_head(
            CurrentKey(-incoming_antiquark.pdg, (2,), -chirality)
        )
        _set_parameter_values_unchecked(builder, anti_head, anti_wf)
        _set_parameter_values_unchecked(builder, active_quark_head, quark_wf)
        _set_parameter_values_unchecked(builder, inactive_quark_head, (0j, 0j))

    def _populate_gluon_source_parameters(
        self,
        *,
        gluon_helicities: tuple[int, ...],
        gluon_vectors: tuple[dict[int, tuple[complex, complex, complex, complex]], ...],
    ) -> None:
        builder = self.bundle.param_builder
        for offset, (helicity, vectors) in enumerate(
            zip(gluon_helicities, gluon_vectors, strict=True),
            start=3,
        ):
            _set_parameter_values_unchecked(
                builder,
                _current_parameter_head(CurrentKey(21, (offset,), 0)),
                vectors[helicity],
            )

    def _populate_gluon_source_values(
        self,
        gluon_wfs: tuple[tuple[complex, complex, complex, complex], ...],
    ) -> None:
        builder = self.bundle.param_builder
        for offset, gluon_wf in enumerate(gluon_wfs, start=3):
            _set_parameter_values_unchecked(
                builder,
                _current_parameter_head(CurrentKey(21, (offset,), 0)),
                gluon_wf,
            )

    def _populate_z_source_parameter(
        self,
        *,
        z_wf: tuple[complex, complex, complex, complex],
    ) -> None:
        _set_parameter_values_unchecked(
            self.bundle.param_builder,
            _current_parameter_head(CurrentKey(23, (self.gluon_count + 3,), 0)),
            z_wf,
        )

    def _populate_momentum_parameters(
        self,
        particles: tuple[ExternalMomentum, ...],
    ) -> None:
        builder = self.bundle.param_builder
        momenta_by_label = _current_momenta_by_label(particles)
        for current in _current_momentum_currents(self.graph):
            momentum = _sum_momenta(
                momenta_by_label[label] for label in current.external_labels
            )
            _set_parameter_values_unchecked(
                builder,
                _current_momentum_parameter_head(current),
                tuple(complex(component) for component in momentum),
            )


class ZGluonNumericTensorNetworkEvaluator:
    """Validation evaluator using fully numeric tensors in a spenso library."""

    def __init__(
        self,
        process: str,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        graph: RecursionGraph | None = None,
        strategy: TensorNetworkStrategy = "interleaved",
    ) -> None:
        self.process = process
        self.model = model or AmplicolSMLeadingColorModel()
        self.graph = graph or _z_gluon_graph(process, self.model)
        self.gluon_count = _validate_z_gluon_graph(self.graph)
        self.strategy: TensorNetworkStrategy = _validate_tensor_network_strategy(strategy)
        start = time.perf_counter()
        self.expression = None
        if self.strategy == "monolithic":
            from symbolica.community.idenso import simplify_color

            raw_expression = _GraphTensorExpressionBuilder(
                self.model,
                self.graph,
            ).matrix_element_skeleton()
            self.expression = simplify_color(raw_expression)
        self.expression_build_time_s = time.perf_counter() - start
        self.last_interleaved_metadata: dict[str, float | int] | None = None

    @property
    def metadata(self) -> NumericTensorNetworkRuntimeMetadata:
        return NumericTensorNetworkRuntimeMetadata(
            process=self.process,
            kernel="spenso-numeric-tensor-network",
            strategy=self.strategy,
            gluon_count=self.gluon_count,
            graph_current_count=len(self.graph.currents),
            graph_interaction_count=len(self.graph.interactions),
            graph_amplitude_count=len(self.graph.amplitudes),
            expression_build_time_s=self.expression_build_time_s,
        )

    def evaluate(
        self,
        particles: Sequence[ExternalMomentum],
    ) -> MatrixElementEvaluation:
        point = _validate_z_gluon_point(particles, gluon_count=self.gluon_count)
        raw_sum = 0.0
        helicity_contributions: list[HelicityContribution] = []
        for helicities, amplitude in self._helicity_amplitudes(point):
            squared = float((amplitude * amplitude.conjugate()).real)
            raw_sum += squared
            helicity_contributions.append(
                HelicityContribution(
                    helicities=helicities,
                    amplitude=amplitude,
                    squared=squared,
                )
            )

        pdgs = tuple(particle.pdg for particle in point)
        color_factor = self.model.leading_color_factor(pdgs)
        average_factor = _initial_state_average_factor(pdgs[:2])
        identical_factor = _final_state_identical_factor(
            particle.pdg for particle in point[2:]
        )
        coupling_factor = (
            (4.0 * math.pi * self.model.alpha_s_me_check) ** self.gluon_count
            * (2.0 * 4.0 * math.pi * self.model.alpha_ew)
        )
        matrix_element = (
            raw_sum
            * color_factor
            * coupling_factor
            / (average_factor * identical_factor)
        )
        return MatrixElementEvaluation(
            process=self.process,
            particles=point,
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=color_factor,
            average_factor=average_factor,
            coupling_factor=coupling_factor,
            helicity_contributions=tuple(helicity_contributions),
            identical_factor=identical_factor,
        )

    def _helicity_amplitudes(
        self,
        particles: tuple[ExternalMomentum, ...],
    ) -> tuple[tuple[tuple[int, ...], complex], ...]:
        incoming_quark = particles[0]
        incoming_antiquark = particles[1]
        gluons = particles[2 : 2 + self.gluon_count]
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
            helicity: _ext_massive_vector(z_momentum, helicity, self.model.mass(23))
            for helicity in (-1, 0, 1)
        }
        amplitudes: list[tuple[tuple[int, ...], complex]] = []
        for (
            physical_quark_helicity,
            physical_antiquark_helicity,
            chirality,
            quark_wf,
            anti_wf,
        ) in (
            (1, -1, 1, quark_minus, anti_plus),
            (-1, 1, -1, quark_plus, anti_minus),
        ):
            for z_helicity in (-1, 0, 1):
                for gluon_helicities in product(
                    (-1, 1),
                    repeat=self.gluon_count,
                ):
                    gluon_wfs = tuple(
                        _ext_gluon_cmplx(momentum, helicity)
                        for momentum, helicity in zip(
                            gluon_momenta,
                            gluon_helicities,
                            strict=True,
                        )
                    )
                    source_values = _z_gluon_source_values(
                        self.graph,
                        particles,
                        gluon_count=self.gluon_count,
                        chirality=chirality,
                        quark_wf=quark_wf,
                        anti_wf=anti_wf,
                        gluon_wfs=gluon_wfs,
                        z_wf=z_vectors[z_helicity],
                    )
                    amplitudes.append(
                        (
                            (
                                physical_quark_helicity,
                                physical_antiquark_helicity,
                                *gluon_helicities,
                                z_helicity,
                            ),
                            self._evaluate_numeric_network(
                                particles,
                                source_values=source_values,
                            ),
                        )
                    )
        return tuple(amplitudes)

    def _evaluate_numeric_network(
        self,
        particles: tuple[ExternalMomentum, ...],
        *,
        source_values: dict[tuple[int, tuple[int, ...], int], tuple[complex, ...]],
    ) -> complex:
        from symbolica.community.spenso import TensorNetwork

        library = self.model.build_tensor_library()
        _register_numeric_source_currents(library, self.graph, source_values)
        _register_numeric_current_momenta_and_propagators(
            self.model,
            library,
            self.graph,
            particles,
        )
        if self.strategy == "interleaved":
            network, metadata = _GraphTensorExpressionBuilder(
                self.model,
                self.graph,
            ).matrix_element_interleaved_network(library, execute_between=True)
            self.last_interleaved_metadata = metadata
        else:
            if self.expression is None:
                raise NativeEvaluationError("monolithic numeric tensor expression is missing")
            network = TensorNetwork(self.expression, library)
        network.execute(library=library)
        scalar = network.result_scalar()
        evaluator = scalar.evaluator([])
        return complex(evaluator.evaluate_complex([[]])[0][0])


class OneGluonTensorNetworkEvaluator(ZGluonTensorNetworkEvaluator):
    """Backward-compatible one-gluon tensor-network evaluator wrapper."""

    def __init__(
        self,
        process: str,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        graph: RecursionGraph | None = None,
        bundle: SymbolicaEvaluatorBundle | None = None,
        strategy: TensorNetworkStrategy = "interleaved",
        helicity_filter: HelicityFilter | None = None,
        build_helicity_filter: bool = False,
    ) -> None:
        super().__init__(
            process,
            model=model,
            graph=graph,
            bundle=bundle,
            strategy=strategy,
            helicity_filter=helicity_filter,
            build_helicity_filter=build_helicity_filter,
        )
        if self.gluon_count != 1:
            raise NativeEvaluationError(
                "one-gluon tensor-network evaluator received a non-one-gluon graph"
            )


def build_z_gluon_tensor_network_evaluator_payload(
    process: str,
    graph: RecursionGraph,
) -> dict[str, object]:
    return ZGluonTensorNetworkEvaluator(
        process,
        graph=graph,
        build_helicity_filter=True,
    ).to_artifact_payload()


def build_one_gluon_tensor_network_evaluator_payload(
    process: str,
    graph: RecursionGraph,
) -> dict[str, object]:
    return OneGluonTensorNetworkEvaluator(
        process,
        graph=graph,
        build_helicity_filter=True,
    ).to_artifact_payload()


def build_z_gluon_helicity_filter(
    evaluator: ZGluonTensorNetworkEvaluator,
    *,
    sample_count: int = 10,
    relative_tolerance: float = 1.0e-12,
    zero_tolerance: float = 1.0e-300,
) -> HelicityFilter:
    """Build a legacy-style helicity filter for a generated tensor evaluator."""

    points = _helicity_filter_sample_points(
        evaluator.process,
        evaluator.model,
        evaluator.gluon_count,
        sample_count=sample_count,
    )
    all_helicities = _all_z_gluon_helicities(evaluator.gluon_count)
    squared_by_helicity: dict[tuple[int, ...], list[float]] = {
        helicities: [] for helicities in all_helicities
    }
    unfiltered_entries = tuple(
        HelicityFilterEntry(
            helicities=helicities,
            multiplicity=1,
            equivalents=(helicities,),
        )
        for helicities in all_helicities
    )
    for point in points:
        for helicities, amplitude, _ in evaluator._evaluate_helicity_entries(
            point,
            unfiltered_entries,
        ):
            squared_by_helicity[helicities].append(
                float((amplitude * amplitude.conjugate()).real)
            )

    active = [
        helicities
        for helicities in all_helicities
        if any(value > zero_tolerance for value in squared_by_helicity[helicities])
    ]

    assigned: set[tuple[int, ...]] = set()
    entries: list[HelicityFilterEntry] = []
    for helicities in active:
        if helicities in assigned:
            continue
        equivalents = [helicities]
        assigned.add(helicities)
        for candidate in active:
            if candidate in assigned:
                continue
            if _same_squared_helicity_signature(
                squared_by_helicity[helicities],
                squared_by_helicity[candidate],
                relative_tolerance=relative_tolerance,
                zero_tolerance=zero_tolerance,
            ):
                equivalents.append(candidate)
                assigned.add(candidate)
        entries.append(
            HelicityFilterEntry(
                helicities=helicities,
                multiplicity=len(equivalents),
                equivalents=tuple(equivalents),
            )
        )

    return HelicityFilter(
        entries=tuple(entries),
        original_count=len(all_helicities),
        sample_count=len(points),
        relative_tolerance=relative_tolerance,
        zero_tolerance=zero_tolerance,
        sample_sqrt_s=tuple(point[0].momentum[0] + point[1].momentum[0] for point in points),
        phase_space_mode="canonical",
        seed=None,
    )


def _z_gluon_graph(
    process: str,
    model: AmplicolSMLeadingColorModel,
) -> RecursionGraph:
    result = NativeMatrixElementGenerator(model=model).generate(
        process,
        write_cache_metadata=False,
    )
    if result.graph is None:
        raise NativeEvaluationError(f"no native graph available for {process}")
    return result.graph


def _validate_z_gluon_graph(graph: RecursionGraph) -> int:
    gluon_count = sum(1 for pdg in graph.process[2:] if pdg == 21)
    if gluon_count < 1 or gluon_count > 6 or graph.process[-1] != 23:
        raise NativeEvaluationError(
            "tensor-network evaluator currently supports q q~ -> Z plus one to six gluons"
        )
    if any(pdg != 21 for pdg in graph.process[2:-1]):
        raise NativeEvaluationError(
            "tensor-network graph is not an ordered Z+gluon graph"
        )
    return gluon_count


def _validate_z_gluon_point(
    particles: Sequence[ExternalMomentum],
    *,
    gluon_count: int,
) -> tuple[ExternalMomentum, ...]:
    point = tuple(particles)
    expected = gluon_count + 3
    if len(point) != expected:
        raise NativeEvaluationError(
            f"q q~ -> Z + {gluon_count} gluons requires exactly {expected} external momenta"
        )
    if point[0].pdg + point[1].pdg != 0 or not 1 <= abs(point[0].pdg) <= 6:
        raise NativeEvaluationError("incoming particles must be a quark/antiquark pair")
    if any(particle.pdg != 21 for particle in point[2 : 2 + gluon_count]):
        raise NativeEvaluationError("final-state gluons must precede the Z boson")
    if point[-1].pdg != 23:
        raise NativeEvaluationError("final state must end with one Z boson")
    return point


def _current_momenta_by_label(
    particles: tuple[ExternalMomentum, ...],
) -> dict[int, tuple[float, float, float, float]]:
    momenta: dict[int, tuple[float, float, float, float]] = {}
    for index, particle in enumerate(particles, start=1):
        momenta[index] = (
            _negate_momentum(particle.momentum)
            if index <= 2
            else particle.momentum
        )
    return momenta


def _z_gluon_source_values(
    graph: RecursionGraph,
    particles: tuple[ExternalMomentum, ...],
    *,
    gluon_count: int,
    chirality: int,
    quark_wf: tuple[complex, complex],
    anti_wf: tuple[complex, complex],
    gluon_wfs: tuple[tuple[complex, complex, complex, complex], ...],
    z_wf: tuple[complex, complex, complex, complex],
) -> dict[tuple[int, tuple[int, ...], int], tuple[complex, ...]]:
    source_values: dict[tuple[int, tuple[int, ...], int], tuple[complex, ...]] = {}
    for current in _source_currents(graph):
        source_values[_current_key_tuple(current)] = tuple(
            0j for _ in range(_current_dimension(current))
        )

    incoming_quark = particles[0]
    incoming_antiquark = particles[1]
    source_values[
        _current_key_tuple(CurrentKey(-incoming_quark.pdg, (1,), 0))
    ] = anti_wf
    source_values[
        _current_key_tuple(CurrentKey(-incoming_antiquark.pdg, (2,), chirality))
    ] = quark_wf
    for offset, gluon_wf in enumerate(gluon_wfs, start=3):
        source_values[_current_key_tuple(CurrentKey(21, (offset,), 0))] = gluon_wf
    source_values[
        _current_key_tuple(CurrentKey(23, (gluon_count + 3,), 0))
    ] = z_wf
    return source_values


def _register_numeric_source_currents(
    library: Any,
    graph: RecursionGraph,
    source_values: dict[tuple[int, tuple[int, ...], int], tuple[complex, ...]],
) -> None:
    from symbolica.community.spenso import LibraryTensor, TensorName

    for current in _source_currents(graph):
        library.register(
            LibraryTensor.dense(
                TensorName(_current_tensor_name(current))(
                    _current_representation(current)
                ),
                source_values[_current_key_tuple(current)],
            )
        )


def _register_numeric_current_momenta_and_propagators(
    model: AmplicolSMLeadingColorModel,
    library: Any,
    graph: RecursionGraph,
    particles: tuple[ExternalMomentum, ...],
) -> None:
    from symbolica.community.spenso import LibraryTensor, Representation, TensorName

    mink = Representation.mink(4)
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    momenta_by_label = _current_momenta_by_label(particles)
    momenta_by_current: dict[
        tuple[int, tuple[int, ...], int],
        tuple[float, float, float, float],
    ] = {}
    for current in _current_momentum_currents(graph):
        momentum = _sum_momenta(
            momenta_by_label[label] for label in current.external_labels
        )
        momenta_by_current[_current_key_tuple(current)] = momentum
        library.register(
            LibraryTensor.dense(
                TensorName(_current_momentum_tensor_name(current))(mink),
                momentum,
            )
        )

    for current in _propagating_currents(graph):
        momentum = momenta_by_current[_current_key_tuple(current)]
        pdg = int(current.pdg)
        if pdg == 21:
            library.register(
                LibraryTensor.dense(
                    TensorName(_propagator_tensor_name(current))(mink, mink),
                    model.gluon_propagator_tensor_data(momentum),
                )
            )
        elif _is_weyl_fermion_current(current):
            library.register(
                LibraryTensor.dense(
                    TensorName(_propagator_tensor_name(current))(weyl, weyl),
                    model.quark_weyl_propagator_tensor_data(
                        momentum,
                        chirality=int(current.chirality),
                    ),
                )
            )


def _current_representation(current: object) -> Any:
    from symbolica.community.spenso import Representation

    pdg = int(getattr(current, "pdg"))
    if pdg == -21:
        return Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    if pdg in (21, 22, 23):
        return Representation.mink(4)
    if 1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16:
        return Representation("pyamplicol::weyl_spinor", 2)
    raise NativeEvaluationError(f"unsupported source-current representation: {pdg}")


def _validate_tensor_network_strategy(strategy: str) -> TensorNetworkStrategy:
    if strategy == "interleaved":
        return "interleaved"
    if strategy == "monolithic":
        return "monolithic"
    raise NativeEvaluationError(
        "tensor-network strategy must be 'interleaved' or 'monolithic'"
    )


def _helicity_filter_from_payload(
    payload: dict[str, object],
) -> HelicityFilter | None:
    filter_payload = payload.get("helicity_filter")
    if filter_payload is None:
        return None
    if not isinstance(filter_payload, dict):
        raise NativeEvaluationError("helicity_filter payload must be a dictionary")
    return HelicityFilter.from_json_dict(filter_payload)


def _strategy_from_bundle_payload(payload: dict[str, object]) -> TensorNetworkStrategy:
    metadata = payload.get("metadata")
    if isinstance(metadata, dict):
        strategy = metadata.get("strategy")
        if strategy == "interleaved" or strategy == "monolithic":
            return strategy
    return "monolithic"


def _all_z_gluon_helicities(gluon_count: int) -> tuple[tuple[int, ...], ...]:
    helicities: list[tuple[int, ...]] = []
    for physical_quark_helicity, physical_antiquark_helicity in (
        (1, -1),
        (-1, 1),
    ):
        for z_helicity in (-1, 0, 1):
            for gluon_helicities in product((-1, 1), repeat=gluon_count):
                helicities.append(
                    (
                        physical_quark_helicity,
                        physical_antiquark_helicity,
                        *gluon_helicities,
                        z_helicity,
                    )
                )
    return tuple(helicities)


def _helicity_filter_sample_points(
    process: str,
    model: AmplicolSMLeadingColorModel,
    gluon_count: int,
    *,
    sample_count: int,
) -> tuple[tuple[ExternalMomentum, ...], ...]:
    native = LeadingColorZJetsNativeEvaluator(model)
    base_sqrt_s = max(model.sqrt_s, 1000.0)
    return tuple(
        native.canonical_z_gluon_point(
            process,
            gluon_count=gluon_count,
            sqrt_s=base_sqrt_s * (1.0 + 0.037 * index),
        )
        for index in range(sample_count)
    )


def _same_squared_helicity_signature(
    left: Sequence[float],
    right: Sequence[float],
    *,
    relative_tolerance: float,
    zero_tolerance: float,
) -> bool:
    if len(left) != len(right):
        return False
    for left_value, right_value in zip(left, right, strict=True):
        denominator = max(abs(left_value) + abs(right_value), zero_tolerance)
        if abs(left_value - right_value) / denominator > relative_tolerance:
            return False
    return True


def _int_tuple(value: object) -> tuple[int, ...]:
    if not isinstance(value, list):
        raise ValueError("expected a list of integers")
    return tuple(int(item) for item in value)


def _int_payload_value(value: object) -> int:
    if isinstance(value, bool):
        raise ValueError("expected integer payload value")
    if isinstance(value, int | str | bytes | bytearray):
        return int(value)
    raise ValueError("expected integer payload value")


def _float_payload_value(value: object) -> float:
    if isinstance(value, int | float | str | bytes | bytearray):
        return float(value)
    raise ValueError("expected floating-point payload value")


def _tensor_network_kernel(gluon_count: int) -> str:
    if gluon_count == 1:
        return "symbolica-one-gluon-tensor-network"
    return "symbolica-z-gluon-tensor-network"


def _sum_momenta(
    momenta: Iterable[tuple[float, float, float, float]],
) -> tuple[float, float, float, float]:
    total = (0.0, 0.0, 0.0, 0.0)
    for momentum in momenta:
        total = (
            total[0] + momentum[0],
            total[1] + momentum[1],
            total[2] + momentum[2],
            total[3] + momentum[3],
        )
    return total


def _set_parameter_values_unchecked(
    param_builder: Any,
    head: tuple[str, ...],
    values: Sequence[complex],
) -> None:
    """Set a known-valid evaluator parameter slice in the helicity hot path."""

    start, stop = param_builder.positions[head]
    param_builder.values[start:stop] = values


def _negate_momentum(
    momentum: tuple[float, float, float, float],
) -> tuple[float, float, float, float]:
    return -momentum[0], -momentum[1], -momentum[2], -momentum[3]


def _tensor_network_bundle_name(process: str) -> str:
    return (
        "tensor_network_"
        + process.replace(" ", "_")
        .replace(">", "to")
        .replace("~", "bar")
        .replace("+", "p")
        .replace("-", "m")
    )


__all__ = [
    "NumericTensorNetworkRuntimeMetadata",
    "OneGluonTensorNetworkEvaluator",
    "TensorNetworkRuntimeMetadata",
    "TensorNetworkStrategy",
    "HelicityFilter",
    "HelicityFilterEntry",
    "ZGluonNumericTensorNetworkEvaluator",
    "ZGluonTensorNetworkEvaluator",
    "build_z_gluon_helicity_filter",
    "build_one_gluon_tensor_network_evaluator_payload",
    "build_z_gluon_tensor_network_evaluator_payload",
]
