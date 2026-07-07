from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Literal, Mapping, Sequence, cast

from .color_plan import GenericColorPlan, LCColorSector, build_color_plan
from .model import (
    AmplicolSMLeadingColorModel,
    CouplingOrders,
    Model,
    QuantumFlow,
    Vertex,
)
from .process_ir import CanonicalProcessIR, ProcessLegIR, build_process_ir
from .processes import ProcessOptions

ColorAccuracy = Literal["lc", "nlc", "full"]


@dataclass(frozen=True)
class ColorState:
    """Colour identity carried by one current.

    The LC implementation stores the warmup sector and the open-line/trace
    groups touched by the current.  NLC/full colour will replace ``basis_key``
    with Idenso basis components without changing the current-index contract.
    """

    accuracy: str
    sector_id: int = 0
    line_groups: tuple[int, ...] = ()
    basis_key: tuple[str, ...] = ()
    _key: tuple[object, ...] = field(init=False, repr=False, compare=False)
    _hash: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        key = (
            self.accuracy,
            self.sector_id,
            self.line_groups,
            self.basis_key,
        )
        object.__setattr__(self, "_key", key)
        object.__setattr__(self, "_hash", hash(key))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ColorState):
            return NotImplemented
        return self._key == other._key

    def __hash__(self) -> int:
        return self._hash

    def to_json_dict(self) -> dict[str, object]:
        return {
            "accuracy": self.accuracy,
            "sector_id": self.sector_id,
            "line_groups": list(self.line_groups),
            "basis_key": list(self.basis_key),
        }


@dataclass(frozen=True)
class ColorFlow:
    state: ColorState
    weight: tuple[float, float] = (1.0, 0.0)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "state": self.state.to_json_dict(),
            "weight": list(self.weight),
        }


@dataclass(frozen=True)
class CurrentIndex:
    """Complete model-driven current identity.

    Equality of this dataclass is the only current deduplication rule.  It
    deliberately carries physics state, not a process-family label.
    """

    particle_id: int
    external_mask: int
    external_labels: tuple[int, ...]
    helicity_ancestry: int
    chirality: int
    spin_state: int | tuple[int, ...]
    flavour_flow: tuple[int, ...]
    charge_flow: int
    color_state: ColorState
    momentum_mask: int
    coupling_orders: tuple[tuple[str, int], ...] = ()
    auxiliary_kind: str | None = None
    ordered_external_labels: tuple[int, ...] = ()
    _key: tuple[object, ...] = field(init=False, repr=False, compare=False)
    _hash: int = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        labels = tuple(sorted(self.external_labels))
        if labels != self.external_labels:
            object.__setattr__(self, "external_labels", labels)
        if not self.ordered_external_labels:
            object.__setattr__(self, "ordered_external_labels", labels)
        elif tuple(sorted(self.ordered_external_labels)) != labels:
            raise ValueError(
                "ordered_external_labels must contain the same labels as "
                "external_labels"
            )
        expected_mask = _labels_mask(self.external_labels)
        if self.external_mask != expected_mask:
            raise ValueError(
                "external_mask does not match external_labels: "
                f"{self.external_mask} != {expected_mask}"
            )
        if self.momentum_mask == 0:
            raise ValueError("momentum_mask must be nonzero for a current")
        if self.coupling_orders:
            object.__setattr__(
                self,
                "coupling_orders",
                tuple(
                    sorted(
                        (str(name).upper(), int(value))
                        for name, value in self.coupling_orders
                        if int(value) != 0
                    )
                ),
            )
        key = (
            self.particle_id,
            self.external_mask,
            self.external_labels,
            self.helicity_ancestry,
            self.chirality,
            self.spin_state,
            self.flavour_flow,
            self.charge_flow,
            self.color_state,
            self.momentum_mask,
            self.coupling_orders,
            self.auxiliary_kind,
            self.ordered_external_labels,
        )
        object.__setattr__(self, "_key", key)
        object.__setattr__(self, "_hash", hash(key))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, CurrentIndex):
            return NotImplemented
        return self._key == other._key

    def __hash__(self) -> int:
        return self._hash

    def overlaps(self, other: CurrentIndex) -> bool:
        return bool(self.external_mask & other.external_mask)

    def to_json_dict(self) -> dict[str, object]:
        spin_state: object
        if isinstance(self.spin_state, tuple):
            spin_state = list(self.spin_state)
        else:
            spin_state = self.spin_state
        return {
            "particle_id": self.particle_id,
            "external_mask": self.external_mask,
            "external_labels": list(self.external_labels),
            "ordered_external_labels": list(self.ordered_external_labels),
            "helicity_ancestry": self.helicity_ancestry,
            "chirality": self.chirality,
            "spin_state": spin_state,
            "flavour_flow": list(self.flavour_flow),
            "charge_flow": self.charge_flow,
            "color_state": self.color_state.to_json_dict(),
            "momentum_mask": self.momentum_mask,
            "coupling_orders": [
                [name, value] for name, value in self.coupling_orders
            ],
            "auxiliary_kind": self.auxiliary_kind,
        }


@dataclass(frozen=True)
class CurrentNode:
    id: int
    index: CurrentIndex
    dimension: int
    is_source: bool
    source_leg_label: int | None = None
    source_helicity: int | None = None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "index": self.index.to_json_dict(),
            "dimension": self.dimension,
            "is_source": self.is_source,
            "source_leg_label": self.source_leg_label,
            "source_helicity": self.source_helicity,
        }


@dataclass(frozen=True)
class InteractionNode:
    id: int
    vertex_kind: int
    vertex_particles: tuple[int, int, int]
    left_id: int
    right_id: int
    result_id: int
    coupling: tuple[float, float]
    color_weight: tuple[float, float]
    lowering_backend: str
    full_tensor_network_ready: bool

    def to_json_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "vertex_kind": self.vertex_kind,
            "vertex_particles": list(self.vertex_particles),
            "left_id": self.left_id,
            "right_id": self.right_id,
            "result_id": self.result_id,
            "coupling": list(self.coupling),
            "color_weight": list(self.color_weight),
            "lowering_backend": self.lowering_backend,
            "full_tensor_network_ready": self.full_tensor_network_ready,
        }


@dataclass(frozen=True)
class AmplitudeRoot:
    id: int
    kind: str
    left_id: int
    right_id: int
    color_weight: tuple[float, float]
    vertex_kind: int | None = None
    vertex_particles: tuple[int, int, int] | None = None
    coupling: tuple[float, float] = (1.0, 0.0)
    contraction: str = ""
    helicity_weight: float = 1.0

    def to_json_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "kind": self.kind,
            "left_id": self.left_id,
            "right_id": self.right_id,
            "color_weight": list(self.color_weight),
            "vertex_kind": self.vertex_kind,
            "vertex_particles": (
                list(self.vertex_particles) if self.vertex_particles is not None else None
            ),
            "coupling": list(self.coupling),
            "contraction": self.contraction,
            "helicity_weight": self.helicity_weight,
        }


@dataclass(frozen=True)
class GenericDAG:
    process: CanonicalProcessIR
    color_plan: GenericColorPlan
    currents: tuple[CurrentNode, ...]
    sources: tuple[int, ...]
    interactions: tuple[InteractionNode, ...]
    amplitude_roots: tuple[AmplitudeRoot, ...]
    truncated: bool = False

    @property
    def has_amplitudes(self) -> bool:
        return bool(self.amplitude_roots)

    @property
    def required_vertex_kinds(self) -> tuple[int, ...]:
        return tuple(
            sorted(
                {
                    interaction.vertex_kind
                    for interaction in self.interactions
                }
                | {
                    root.vertex_kind
                    for root in self.amplitude_roots
                    if root.vertex_kind is not None
                }
            )
        )

    def currents_by_external_labels(
        self,
        labels: Iterable[int],
    ) -> tuple[CurrentNode, ...]:
        wanted = tuple(sorted(labels))
        return tuple(current for current in self.currents if current.index.external_labels == wanted)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process.to_json_dict(),
            "color_plan": self.color_plan.to_json_dict(),
            "current_count": len(self.currents),
            "source_count": len(self.sources),
            "interaction_count": len(self.interactions),
            "amplitude_root_count": len(self.amplitude_roots),
            "truncated": self.truncated,
            "required_vertex_kinds": list(self.required_vertex_kinds),
            "currents": [current.to_json_dict() for current in self.currents],
            "sources": list(self.sources),
            "interactions": [
                interaction.to_json_dict() for interaction in self.interactions
            ],
            "amplitude_roots": [root.to_json_dict() for root in self.amplitude_roots],
        }


class ColorEngine:
    """Local colour-flow engine used by the process-generic recursion."""

    def __init__(
        self,
        color_plan: GenericColorPlan,
        model: Model,
    ) -> None:
        self.color_plan = color_plan
        self.model = model
        self._sector_by_id = {sector.id: sector for sector in color_plan.sectors}
        self._leg_by_label = {leg.label: leg for leg in color_plan.process.legs}
        self._label_is_color_singlet: dict[int, bool] = {}
        for label, leg in self._leg_by_label.items():
            if leg.outgoing_pdg is None:
                self._label_is_color_singlet[label] = False
                continue
            try:
                self._label_is_color_singlet[label] = (
                    self.model.color_rep(leg.outgoing_pdg) == 1
                )
            except KeyError:
                self._label_is_color_singlet[label] = False
        self._vertex_has_colour_cache: dict[tuple[int, int, int, int], bool] = {}
        self._ordered_combination_labels_cache: dict[
            tuple[int, tuple[int, ...], tuple[int, ...], int, int, int, int],
            tuple[int, ...] | None,
        ] = {}

    def source_states_for_leg(self, leg: ProcessLegIR) -> tuple[ColorState, ...]:
        if self.color_plan.color_accuracy != "lc":
            return (
                ColorState(
                    accuracy=self.color_plan.color_accuracy,
                    basis_key=("idenso-basis-placeholder",),
                ),
            )
        if not self.color_plan.sectors:
            return (ColorState(accuracy="lc"),)
        states = []
        leg_is_singlet = self._label_is_color_singlet.get(leg.label, False)
        for sector in self.color_plan.sectors:
            groups = (
                ()
                if leg_is_singlet
                else _sector_group_indices_for_label(sector, leg.label)
            )
            if groups or leg_is_singlet or not sector.line_label_groups:
                states.append(
                    ColorState(
                        accuracy="lc",
                        sector_id=sector.id,
                        line_groups=groups,
                    )
                )
        return tuple(states)

    def combine(
        self,
        left: ColorState,
        right: ColorState,
        vertex: Vertex,
    ) -> tuple[ColorFlow, ...]:
        if left.accuracy != right.accuracy:
            return ()
        if left.accuracy != "lc":
            return (
                ColorFlow(
                    state=ColorState(
                        accuracy=left.accuracy,
                        basis_key=tuple(
                            sorted(set(left.basis_key) | set(right.basis_key))
                        ),
                    )
                ),
            )
        if left.sector_id != right.sector_id:
            return ()
        groups = self._lc_combined_line_groups(left, right, vertex)
        if groups is None:
            return ()
        return (
            ColorFlow(
                state=ColorState(
                    accuracy="lc",
                    sector_id=left.sector_id,
                    line_groups=groups,
                )
            ),
        )

    def _lc_combined_line_groups(
        self,
        left: ColorState,
        right: ColorState,
        vertex: Vertex,
    ) -> tuple[int, ...] | None:
        left_colored = self._particle_has_colour(vertex.particles[0])
        right_colored = self._particle_has_colour(vertex.particles[1])
        result_colored = self._particle_has_colour(vertex.particles[2])
        left_groups = set(left.line_groups)
        right_groups = set(right.line_groups)

        if not left_colored and not right_colored:
            if left_groups and right_groups and left_groups != right_groups:
                return None
            return tuple(sorted(left_groups | right_groups))

        if left_colored != right_colored:
            colored_groups = left_groups if left_colored else right_groups
            singlet_groups = right_groups if left_colored else left_groups
            if singlet_groups and colored_groups and not singlet_groups.issubset(
                colored_groups
            ):
                return None
            return tuple(sorted(colored_groups or singlet_groups))

        if not result_colored and left_groups and right_groups and left_groups != right_groups:
            return None
        return tuple(sorted(left_groups | right_groups))

    def _particle_has_colour(self, pdg: int) -> bool:
        if pdg == 99:
            return True
        try:
            return self.model.color_rep(pdg) != 1
        except KeyError:
            return True

    def closure_compatible(
        self,
        left: ColorState,
        right: ColorState,
        *,
        full_mask: int,
    ) -> tuple[ColorFlow, ...]:
        if left.accuracy != right.accuracy:
            return ()
        if left.accuracy != "lc":
            return (
                ColorFlow(
                    state=ColorState(
                        accuracy=left.accuracy,
                        basis_key=tuple(
                            sorted(set(left.basis_key) | set(right.basis_key))
                        ),
                    )
                ),
            )
        if left.sector_id != right.sector_id:
            return ()
        groups = tuple(sorted(set(left.line_groups) | set(right.line_groups)))
        sector = self._sector_by_id.get(left.sector_id)
        if sector is None:
            return (ColorFlow(state=left),)
        full_labels = set(_mask_labels(full_mask))
        sector_labels = set(sector.singlet_labels)
        for group in sector.line_label_groups:
            sector_labels.update(group)
        if full_labels and not full_labels.issubset(sector_labels):
            return ()
        return (
            ColorFlow(
                state=ColorState(
                    accuracy="lc",
                    sector_id=left.sector_id,
                    line_groups=groups,
                )
            ),
        )

    def _vertex_has_colour(self, vertex: Vertex) -> bool:
        key = (vertex.kind, *vertex.particles)
        cached = self._vertex_has_colour_cache.get(key)
        if cached is not None:
            return cached
        for pdg in vertex.particles:
            try:
                if self.model.color_rep(pdg) != 1:
                    self._vertex_has_colour_cache[key] = True
                    return True
            except KeyError:
                self._vertex_has_colour_cache[key] = True
                return True
        self._vertex_has_colour_cache[key] = False
        return False

    def vertex_allowed(self, vertex: Vertex) -> bool:
        if self.color_plan.color_accuracy != "lc":
            return True
        if 99 not in vertex.particles:
            return True
        return self.color_plan.process.quark_lines.quark_pair_count >= 2

    def ordered_combination_allowed(
        self,
        left_index: CurrentIndex,
        right_index: CurrentIndex,
        vertex: Vertex,
    ) -> bool:
        if self.color_plan.color_accuracy != "lc":
            return True
        if not self._vertex_has_colour(vertex):
            return True
        sector = self._sector_by_id.get(left_index.color_state.sector_id)
        if sector is None:
            return True
        return any(
            _ordered_combination_matches_word(
                left_index.ordered_external_labels,
                right_index.ordered_external_labels,
                word,
            )
            for word in _sector_intermediate_order_words(sector)
        )

    def ordered_combination_labels(
        self,
        left_index: CurrentIndex,
        right_index: CurrentIndex,
        vertex: Vertex,
    ) -> tuple[int, ...] | None:
        cache_key = (
            left_index.color_state.sector_id,
            left_index.ordered_external_labels,
            right_index.ordered_external_labels,
            vertex.kind,
            *vertex.particles,
        )
        if cache_key in self._ordered_combination_labels_cache:
            return self._ordered_combination_labels_cache[cache_key]
        proposed = (
            *left_index.ordered_external_labels,
            *right_index.ordered_external_labels,
        )
        singlet_order = self._singlet_order_allowed(
            left_index,
            right_index,
            vertex,
        )
        if not singlet_order:
            self._ordered_combination_labels_cache[cache_key] = None
            return None
        if self.color_plan.color_accuracy != "lc":
            self._ordered_combination_labels_cache[cache_key] = proposed
            return proposed
        sector = self._sector_by_id.get(left_index.color_state.sector_id)
        if sector is None:
            self._ordered_combination_labels_cache[cache_key] = proposed
            return proposed
        for word in _sector_intermediate_order_words(sector):
            segment = _ordered_combination_segment(
                left_index.ordered_external_labels,
                right_index.ordered_external_labels,
                word,
            )
            if segment is None:
                continue
            word_labels = set(word)
            extras = tuple(sorted(label for label in proposed if label not in word_labels))
            if not _line_local_singlet_extras_allowed(segment, extras, sector):
                continue
            labels = (*segment, *extras)
            self._ordered_combination_labels_cache[cache_key] = labels
            return labels
        self._ordered_combination_labels_cache[cache_key] = None
        return None

    def _singlet_order_allowed(
        self,
        left_index: CurrentIndex,
        right_index: CurrentIndex,
        vertex: Vertex,
    ) -> bool:
        left_singlet = self._all_external_labels_color_singlet(
            left_index.external_labels
        )
        right_singlet = self._all_external_labels_color_singlet(
            right_index.external_labels
        )
        if left_singlet and not right_singlet:
            return False
        if left_singlet and right_singlet:
            return max(left_index.external_labels) < max(right_index.external_labels)
        return True

    def _all_external_labels_color_singlet(self, labels: Iterable[int]) -> bool:
        for label in labels:
            is_singlet = self._label_is_color_singlet.get(label)
            if is_singlet is None:
                return False
            if not is_singlet:
                return False
        return True

    def ordered_closure_allowed(
        self,
        left_index: CurrentIndex,
        right_index: CurrentIndex,
    ) -> bool:
        if self.color_plan.color_accuracy != "lc":
            return True
        if left_index.color_state.sector_id != right_index.color_state.sector_id:
            return False
        sector = self._sector_by_id.get(left_index.color_state.sector_id)
        if sector is None:
            return True
        return any(
            _closure_combination_matches_word(
                _labels_projected_to_word(left_index.ordered_external_labels, word),
                _labels_projected_to_word(right_index.ordered_external_labels, word),
                word,
            )
            for word in sector.compatibility_words
        )


class GenericDAGCompiler:
    """Compile a concrete process into a model-driven current DAG.

    The compiler never classifies the whole process as a family.  It sweeps
    external subsets, asks the model which local vertices are valid for two
    current particle ids, asks the colour engine whether their colour states
    combine, and deduplicates solely by ``CurrentIndex`` equality.
    """

    def __init__(
        self,
        *,
        model: Model | None = None,
        color_accuracy: str = "lc",
        options: ProcessOptions | None = None,
        max_currents: int = 50000,
        max_color_sectors: int = 20000,
        reference_color_order: tuple[int, ...] | None = None,
        selected_color_sector_ids: Iterable[int] | None = None,
        max_coupling_orders: Mapping[str, int] | None = None,
        max_lc_current_line_groups: int | None = None,
        max_quark_pairs: int | None = None,
        closure_side_mask_pruning: bool = True,
        color_order_mask_pruning: bool = True,
        species_reachability_pruning: bool = True,
        ignored_particle_ids: Iterable[int] | None = None,
        ignored_vertex_kinds: Iterable[int] | None = None,
    ) -> None:
        self.model = model or AmplicolSMLeadingColorModel()
        self.color_accuracy = color_accuracy
        self.options = options
        self.max_currents = max_currents
        self.max_color_sectors = max_color_sectors
        self.reference_color_order = reference_color_order
        self.selected_color_sector_ids = (
            None
            if selected_color_sector_ids is None
            else frozenset(int(sector_id) for sector_id in selected_color_sector_ids)
        )
        self.max_coupling_orders = _normalize_coupling_order_limits(
            max_coupling_orders,
        )
        self.max_lc_current_line_groups = (
            None
            if max_lc_current_line_groups is None
            else max(0, int(max_lc_current_line_groups))
        )
        self.max_quark_pairs = (
            None
            if max_quark_pairs is None
            else max(0, int(max_quark_pairs))
        )
        self.closure_side_mask_pruning = bool(closure_side_mask_pruning)
        self.color_order_mask_pruning = bool(color_order_mask_pruning)
        self.species_reachability_pruning = bool(species_reachability_pruning)
        self.ignored_particle_ids = frozenset(
            int(particle_id) for particle_id in (ignored_particle_ids or ())
        )
        self.ignored_vertex_kinds = frozenset(
            int(kind) for kind in (ignored_vertex_kinds or ())
        )

    def compile(self, process: str | CanonicalProcessIR) -> GenericDAG:
        process_ir = (
            process
            if isinstance(process, CanonicalProcessIR)
            else build_process_ir(
                process,
                color_accuracy=self.color_accuracy,
                options=self.options,
            )
        )
        color_plan = build_color_plan(
            process_ir,
            color_accuracy=process_ir.color_accuracy,
            options=self.options,
            max_sectors=self.max_color_sectors,
            reference_color_order=self.reference_color_order,
        )
        if (
            self.max_quark_pairs is not None
            and process_ir.quark_lines.quark_pair_count > self.max_quark_pairs
        ):
            return GenericDAG(
                process=process_ir,
                color_plan=color_plan,
                currents=(),
                sources=(),
                interactions=(),
                amplitude_roots=(),
                truncated=False,
            )
        if self.selected_color_sector_ids is not None and color_plan.color_accuracy == "lc":
            selected_sectors = tuple(
                sector
                for sector in color_plan.sectors
                if sector.id in self.selected_color_sector_ids
            )
            missing_sector_ids = tuple(
                sorted(
                    int(sector_id)
                    for sector_id in self.selected_color_sector_ids
                    if all(sector.id != sector_id for sector in selected_sectors)
                )
            )
            diagnostics = color_plan.diagnostics
            if missing_sector_ids:
                diagnostics = (
                    *diagnostics,
                    "selected LC colour sector ids were not materialized: "
                    + ", ".join(str(sector_id) for sector_id in missing_sector_ids),
                )
            color_plan = GenericColorPlan(
                process=color_plan.process,
                color_accuracy=color_plan.color_accuracy,
                sectors=selected_sectors,
                diagnostics=diagnostics,
                truncated=bool(missing_sector_ids),
                idenso_required=color_plan.idenso_required,
            )
        color_engine = ColorEngine(color_plan, self.model)
        table = _CurrentTable(self.model)
        sources = self._build_sources(process_ir, color_engine, table)
        interactions: list[InteractionNode] = []
        interaction_keys: set[tuple[int, int, int, int]] = set()
        vertices_by_input: dict[tuple[int, int], tuple[Vertex, ...]] = {}
        right_particles_by_left = _right_particles_by_left(
            self.model,
            color_accuracy=process_ir.color_accuracy,
        )
        vertex_allowed_cache: dict[tuple[int, int, int, int], bool] = {}
        quantum_flow_cache: dict[tuple[object, ...], tuple[QuantumFlow, ...]] = {}
        full_mask = _labels_mask(leg.label for leg in process_ir.legs)
        closure_candidate_splits = _closure_candidate_splits(
            process_ir,
            self.model,
            color_engine,
            reference_color_order=self.reference_color_order,
        )
        closure_reachable_masks = (
            _closure_side_reachable_masks(
                full_mask,
                closure_candidate_splits,
            )
            if self.closure_side_mask_pruning
            else None
        )
        color_order_reachable_masks = (
            _lc_color_order_reachable_masks(
                process_ir,
                color_plan,
                self.model,
            )
            if self.color_order_mask_pruning
            else None
        )
        useful_states_by_mask = (
            _useful_states_by_mask(
                process_ir,
                self.model,
                color_engine,
                closure_candidate_splits,
                closure_reachable_masks,
                color_order_reachable_masks,
                max_coupling_orders=self.max_coupling_orders,
                ignored_particle_ids=self.ignored_particle_ids,
                ignored_vertex_kinds=self.ignored_vertex_kinds,
            )
            if self.species_reachability_pruning
            else None
        )
        truncated = False
        if any(
            int(leg.outgoing_pdg or 0) in self.ignored_particle_ids
            for leg in process_ir.legs
        ):
            return GenericDAG(
                process=process_ir,
                color_plan=color_plan,
                currents=tuple(table.currents),
                sources=tuple(sources),
                interactions=(),
                amplitude_roots=(),
                truncated=False,
            )

        for mask in _masks_by_size(full_mask):
            if mask & (mask - 1) == 0:
                continue
            if mask == full_mask:
                continue
            if not _mask_allowed_by_reachability(
                mask,
                closure_reachable_masks,
                color_order_reachable_masks,
            ):
                continue
            if useful_states_by_mask is not None and mask not in useful_states_by_mask:
                continue
            labels = _mask_labels(mask)
            for left_mask, right_mask in _ordered_splits(mask):
                if not (
                    _mask_allowed_by_reachability(
                        left_mask,
                        closure_reachable_masks,
                        color_order_reachable_masks,
                    )
                    and _mask_allowed_by_reachability(
                        right_mask,
                        closure_reachable_masks,
                        color_order_reachable_masks,
                    )
                ):
                    continue
                if useful_states_by_mask is not None and (
                    left_mask not in useful_states_by_mask
                    or right_mask not in useful_states_by_mask
                ):
                    continue
                left_ids = table.ids_by_mask(left_mask)
                if not left_ids or not table.has_mask(right_mask):
                    continue
                for left_id in left_ids:
                    left = table.current(left_id)
                    if (
                        useful_states_by_mask is not None
                        and not _state_allowed_by_reachability(
                            useful_states_by_mask,
                            left_mask,
                            left.index.particle_id,
                            left.index.coupling_orders,
                        )
                    ):
                        continue
                    possible_right_particles = right_particles_by_left.get(
                        left.index.particle_id,
                    )
                    if not possible_right_particles:
                        continue
                    candidate_right_ids = table.ids_by_mask_and_particles(
                        right_mask,
                        possible_right_particles,
                    )
                    if not candidate_right_ids:
                        continue
                    for right_id in candidate_right_ids:
                        right = table.current(right_id)
                        if left.index.overlaps(right.index):
                            continue
                        if (
                            useful_states_by_mask is not None
                            and not _state_allowed_by_reachability(
                                useful_states_by_mask,
                                right_mask,
                                right.index.particle_id,
                                right.index.coupling_orders,
                            )
                        ):
                            continue
                        vertex_lookup_key = (
                            left.index.particle_id,
                            right.index.particle_id,
                        )
                        if vertex_lookup_key in vertices_by_input:
                            vertices = vertices_by_input[vertex_lookup_key]
                        else:
                            vertices = self.model.vertices_accepting(
                                left.index.particle_id,
                                right.index.particle_id,
                                color_accuracy=process_ir.color_accuracy,
                            )
                            vertices_by_input[vertex_lookup_key] = vertices
                        for vertex in vertices:
                            if (
                                vertex.kind in self.ignored_vertex_kinds
                                or vertex.particles[2] in self.ignored_particle_ids
                            ):
                                continue
                            if (
                                useful_states_by_mask is not None
                                and vertex.particles[2]
                                not in useful_states_by_mask.get(mask, {})
                            ):
                                continue
                            vertex_key = (vertex.kind, *vertex.particles)
                            if vertex_key in vertex_allowed_cache:
                                vertex_allowed = vertex_allowed_cache[vertex_key]
                            else:
                                vertex_allowed = color_engine.vertex_allowed(vertex)
                                vertex_allowed_cache[vertex_key] = vertex_allowed
                            if not vertex_allowed:
                                continue
                            if self.model.skip_duplicate_vertex_orientation(vertex):
                                continue
                            ordered_external_labels = color_engine.ordered_combination_labels(
                                left.index,
                                right.index,
                                vertex,
                            )
                            if ordered_external_labels is None:
                                continue
                            coupling_orders = self.model.combine_coupling_orders(
                                left.index,
                                right.index,
                                vertex,
                            )
                            if not _coupling_orders_within_limits(
                                coupling_orders,
                                self.max_coupling_orders,
                            ):
                                continue
                            if (
                                useful_states_by_mask is not None
                                and not _state_allowed_by_reachability(
                                    useful_states_by_mask,
                                    mask,
                                    vertex.particles[2],
                                    coupling_orders,
                                )
                            ):
                                continue
                            quantum_flow_key = (
                                vertex.kind,
                                vertex.particles,
                                left.index.particle_id,
                                left.index.chirality,
                                left.index.flavour_flow,
                                right.index.particle_id,
                                right.index.chirality,
                                right.index.flavour_flow,
                            )
                            if quantum_flow_key in quantum_flow_cache:
                                quantum_flows = quantum_flow_cache[quantum_flow_key]
                            else:
                                quantum_flows = self.model.allowed_quantum_flows(
                                    vertex,
                                    left.index,
                                    right.index,
                                )
                                quantum_flow_cache[quantum_flow_key] = quantum_flows
                            for quantum_flow in quantum_flows:
                                for color_flow in color_engine.combine(
                                    left.index.color_state,
                                    right.index.color_state,
                                    vertex,
                                ):
                                    if not _lc_line_groups_within_limit(
                                        color_flow.state,
                                        self.max_lc_current_line_groups,
                                    ):
                                        continue
                                    out_index = CurrentIndex(
                                        particle_id=vertex.particles[2],
                                        external_mask=mask,
                                        external_labels=labels,
                                        ordered_external_labels=ordered_external_labels,
                                        helicity_ancestry=(
                                            left.index.helicity_ancestry
                                            | right.index.helicity_ancestry
                                        ),
                                        chirality=quantum_flow.chirality,
                                        spin_state=quantum_flow.spin_state,
                                        flavour_flow=quantum_flow.flavour_flow,
                                        charge_flow=quantum_flow.charge_flow,
                                        color_state=color_flow.state,
                                        momentum_mask=(
                                            left.index.momentum_mask
                                            | right.index.momentum_mask
                                        ),
                                        coupling_orders=coupling_orders,
                                        auxiliary_kind=self.model.auxiliary_kind(
                                            vertex.particles[2]
                                        ),
                                    )
                                    if not self.model.current_allowed(out_index):
                                        continue
                                    result = table.add_or_get(
                                        out_index,
                                        is_source=False,
                                    )
                                    key = (
                                        vertex.kind,
                                        left_id,
                                        right_id,
                                        result.id,
                                    )
                                    if key in interaction_keys:
                                        continue
                                    interaction_keys.add(key)
                                    rule = self.model.vertex_lowering_rule(vertex.kind)
                                    interactions.append(
                                        InteractionNode(
                                            id=len(interactions),
                                            vertex_kind=vertex.kind,
                                            vertex_particles=vertex.particles,
                                            left_id=left_id,
                                            right_id=right_id,
                                            result_id=result.id,
                                            coupling=quantum_flow.coupling,
                                            color_weight=color_flow.weight,
                                            lowering_backend=rule.backend,
                                            full_tensor_network_ready=(
                                                rule.full_tensor_network_ready
                                            ),
                                        )
                                    )
                                    if len(table.currents) > self.max_currents:
                                        truncated = True
                                        return GenericDAG(
                                            process=process_ir,
                                            color_plan=color_plan,
                                            currents=tuple(table.currents),
                                            sources=tuple(sources),
                                            interactions=tuple(interactions),
                                            amplitude_roots=tuple(
                                                self._build_amplitude_roots(
                                                    process_ir,
                                                    table,
                                                    color_engine,
                                                    candidate_splits=closure_candidate_splits,
                                                )
                                            ),
                                            truncated=truncated,
                                        )

        return GenericDAG(
            process=process_ir,
            color_plan=color_plan,
            currents=tuple(table.currents),
            sources=tuple(sources),
            interactions=tuple(interactions),
            amplitude_roots=tuple(
                self._build_amplitude_roots(
                    process_ir,
                    table,
                    color_engine,
                    candidate_splits=closure_candidate_splits,
                )
            ),
            truncated=truncated,
        )

    def _build_sources(
        self,
        process_ir: CanonicalProcessIR,
        color_engine: ColorEngine,
        table: "_CurrentTable",
    ) -> list[int]:
        sources: list[int] = []
        next_source_bit = 0
        for leg in process_ir.legs:
            if leg.outgoing_pdg is None:
                continue
            particle_id = int(leg.outgoing_pdg)
            for color_state in color_engine.source_states_for_leg(leg):
                if not _lc_line_groups_within_limit(
                    color_state,
                    self.max_lc_current_line_groups,
                ):
                    continue
                for source_state in self.model.source_spin_states(particle_id):
                    chirality = source_state.chirality
                    source_helicity = source_state.helicity
                    spin_state = source_state.spin_state
                    if leg.is_initial and self.model.is_chiral_eligible(particle_id):
                        chirality = -chirality
                        spin_state = chirality
                    elif leg.is_initial and self.model.is_gluon(particle_id):
                        source_helicity = -source_helicity
                        if not isinstance(spin_state, int):
                            raise TypeError("gluon source spin_state must be an integer")
                        spin_state = -spin_state
                    helicity_ancestry = 1 << next_source_bit
                    next_source_bit += 1
                    index = CurrentIndex(
                        particle_id=particle_id,
                        external_mask=1 << (leg.label - 1),
                        external_labels=(leg.label,),
                        ordered_external_labels=(leg.label,),
                        helicity_ancestry=helicity_ancestry,
                        chirality=chirality,
                        spin_state=spin_state,
                        flavour_flow=(particle_id,),
                        charge_flow=self.model.charge_units(particle_id),
                        color_state=color_state,
                        momentum_mask=1 << (leg.label - 1),
                        coupling_orders=(),
                        auxiliary_kind=self.model.auxiliary_kind(particle_id),
                    )
                    current = table.add_or_get(
                        index,
                        is_source=True,
                        source_leg_label=leg.label,
                        source_helicity=source_helicity,
                    )
                    sources.append(current.id)
        return sources

    def _build_amplitude_roots(
        self,
        process_ir: CanonicalProcessIR,
        table: "_CurrentTable",
        color_engine: ColorEngine,
        *,
        candidate_splits: tuple[tuple[int, int], ...] | None = None,
    ) -> list[AmplitudeRoot]:
        full_mask = _labels_mask(leg.label for leg in process_ir.legs)
        if candidate_splits is None:
            candidate_splits = _closure_candidate_splits(
                process_ir,
                self.model,
                color_engine,
            )
        roots: list[AmplitudeRoot] = []
        seen: set[tuple[object, ...]] = set()
        for left_mask, right_mask in candidate_splits:
            if left_mask == 0 or right_mask == 0:
                continue
            for left_id in table.ids_by_mask(left_mask):
                left = table.current(left_id)
                for right_id in table.ids_by_mask(right_mask):
                    right = table.current(right_id)
                    if left.index.overlaps(right.index):
                        continue
                    if not color_engine.ordered_closure_allowed(
                        left.index,
                        right.index,
                    ):
                        continue
                    if (
                        self.reference_color_order is not None
                        and process_ir.color_accuracy == "lc"
                    ):
                        sector = color_engine.color_plan.sector(
                            left.index.color_state.sector_id
                        )
                        if (
                            sector is not None
                            and self.reference_color_order in sector.compatibility_words
                            and not _closure_combination_matches_word(
                                _labels_projected_to_word(
                                    left.index.ordered_external_labels,
                                    self.reference_color_order,
                                ),
                                _labels_projected_to_word(
                                    right.index.ordered_external_labels,
                                    self.reference_color_order,
                                ),
                                self.reference_color_order,
                            )
                        ):
                            continue
                    color_flows = color_engine.closure_compatible(
                        left.index.color_state,
                        right.index.color_state,
                        full_mask=full_mask,
                    )
                    if not color_flows:
                        continue
                    direct = _direct_contraction_kind(
                        self.model,
                        left.index,
                        right.index,
                    )
                    for color_flow in color_flows:
                        if direct is not None:
                            direct_key: tuple[object, ...] = (
                                "direct",
                                direct,
                                left_id,
                                right_id,
                                color_flow.state,
                            )
                            if direct_key not in seen:
                                seen.add(direct_key)
                                roots.append(
                                    AmplitudeRoot(
                                        id=len(roots),
                                        kind="direct-contraction",
                                        left_id=left_id,
                                        right_id=right_id,
                                        color_weight=color_flow.weight,
                                        contraction=direct,
                                    )
                                )
                        for vertex in self.model.vertices_accepting(
                            left.index.particle_id,
                            right.index.particle_id,
                            color_accuracy=process_ir.color_accuracy,
                        ):
                            if (
                                vertex.kind in self.ignored_vertex_kinds
                                or vertex.particles[2] in self.ignored_particle_ids
                            ):
                                continue
                            if not color_engine.vertex_allowed(vertex):
                                continue
                            coupling_orders = self.model.combine_coupling_orders(
                                left.index,
                                right.index,
                                vertex,
                            )
                            if not _coupling_orders_within_limits(
                                coupling_orders,
                                self.max_coupling_orders,
                            ):
                                continue
                            closure_contraction = _closure_contraction_name(
                                self.model,
                                vertex.particles[2],
                            )
                            if closure_contraction != "scalar":
                                continue
                            if not self.model.allowed_quantum_flows(
                                vertex,
                                left.index,
                                right.index,
                            ):
                                continue
                            vertex_key: tuple[object, ...] = (
                                "vertex",
                                vertex.kind,
                                vertex.particles,
                                left_id,
                                right_id,
                                color_flow.state,
                            )
                            if vertex_key in seen:
                                continue
                            seen.add(vertex_key)
                            roots.append(
                                AmplitudeRoot(
                                    id=len(roots),
                                    kind="vertex-closure",
                                    left_id=left_id,
                                    right_id=right_id,
                                    color_weight=color_flow.weight,
                                    vertex_kind=vertex.kind,
                                    vertex_particles=vertex.particles,
                                    coupling=vertex.coupling,
                                    contraction=closure_contraction,
                                )
                            )
        return roots


class _CurrentTable:
    def __init__(self, model: Model) -> None:
        self.model = model
        self.currents: list[CurrentNode] = []
        self._ids: dict[CurrentIndex, int] = {}
        self._ids_by_mask: dict[int, list[int]] = {}
        self._ids_by_mask_particle: dict[tuple[int, int], list[int]] = {}

    def add_or_get(
        self,
        index: CurrentIndex,
        *,
        is_source: bool,
        source_leg_label: int | None = None,
        source_helicity: int | None = None,
    ) -> CurrentNode:
        current_id = self._ids.get(index)
        if current_id is not None:
            return self.currents[current_id]
        current_id = len(self.currents)
        node = CurrentNode(
            id=current_id,
            index=index,
            dimension=self.model.current_dimension(
                index.particle_id,
                index.chirality,
            ),
            is_source=is_source,
            source_leg_label=source_leg_label,
            source_helicity=source_helicity,
        )
        self._ids[index] = current_id
        self.currents.append(node)
        self._ids_by_mask.setdefault(index.external_mask, []).append(current_id)
        self._ids_by_mask_particle.setdefault(
            (index.external_mask, index.particle_id),
            [],
        ).append(current_id)
        return node

    def current(self, current_id: int) -> CurrentNode:
        return self.currents[current_id]

    def has_mask(self, mask: int) -> bool:
        return mask in self._ids_by_mask

    def ids_by_mask(self, mask: int) -> Sequence[int]:
        return self._ids_by_mask.get(mask, ())

    def ids_by_mask_and_particles(
        self,
        mask: int,
        particle_ids: Iterable[int],
    ) -> Sequence[int]:
        ids: list[int] = []
        for particle_id in particle_ids:
            ids.extend(self._ids_by_mask_particle.get((mask, particle_id), ()))
        return ids


def compile_generic_dag(
    process: str | CanonicalProcessIR,
    *,
    model: Model | None = None,
    color_accuracy: str = "lc",
    options: ProcessOptions | None = None,
    max_currents: int = 50000,
    max_color_sectors: int = 20000,
    reference_color_order: tuple[int, ...] | None = None,
    selected_color_sector_ids: Iterable[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    max_lc_current_line_groups: int | None = None,
    max_quark_pairs: int | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    species_reachability_pruning: bool = True,
    ignored_particle_ids: Iterable[int] | None = None,
    ignored_vertex_kinds: Iterable[int] | None = None,
) -> GenericDAG:
    return GenericDAGCompiler(
        model=model,
        color_accuracy=color_accuracy,
        options=options,
        max_currents=max_currents,
        max_color_sectors=max_color_sectors,
        reference_color_order=reference_color_order,
        selected_color_sector_ids=selected_color_sector_ids,
        max_coupling_orders=max_coupling_orders,
        max_lc_current_line_groups=max_lc_current_line_groups,
        max_quark_pairs=max_quark_pairs,
        closure_side_mask_pruning=closure_side_mask_pruning,
        color_order_mask_pruning=color_order_mask_pruning,
        species_reachability_pruning=species_reachability_pruning,
        ignored_particle_ids=ignored_particle_ids,
        ignored_vertex_kinds=ignored_vertex_kinds,
    ).compile(process)


def infer_minimal_coupling_order_limits(
    process: str | CanonicalProcessIR,
    *,
    model: Model | None = None,
    color_accuracy: str = "lc",
    options: ProcessOptions | None = None,
    max_color_sectors: int = 20000,
    selected_color_sector_ids: Iterable[int] | None = None,
    max_coupling_orders: Mapping[str, int] | None = None,
    closure_side_mask_pruning: bool = True,
    color_order_mask_pruning: bool = True,
    ignored_particle_ids: Iterable[int] | None = None,
    ignored_vertex_kinds: Iterable[int] | None = None,
) -> dict[str, int]:
    """Infer a generic lowest-order coupling envelope for a process.

    This is an opt-in generation accelerator.  It never recognizes a whole
    process family; it asks the model which local vertices can connect the
    external states, tracks UFO-style coupling orders, and returns the
    component-wise maximum over all minimal-total-order closure paths.  The
    returned dictionary can be used as ordinary ``max_coupling_orders``.
    """

    active_model = model or AmplicolSMLeadingColorModel()
    process_ir = (
        process
        if isinstance(process, CanonicalProcessIR)
        else build_process_ir(process, color_accuracy=color_accuracy, options=options)
    )
    color_plan = build_color_plan(
        process_ir,
        color_accuracy=process_ir.color_accuracy,
        options=options,
        max_sectors=max_color_sectors,
    )
    explicit_sector_ids = (
        None
        if selected_color_sector_ids is None
        else frozenset(int(sector_id) for sector_id in selected_color_sector_ids)
    )
    if explicit_sector_ids is not None and color_plan.color_accuracy == "lc":
        color_plan = GenericColorPlan(
            process=color_plan.process,
            color_accuracy=color_plan.color_accuracy,
            sectors=tuple(
                sector
                for sector in color_plan.sectors
                if sector.id in explicit_sector_ids
            ),
            diagnostics=color_plan.diagnostics,
            truncated=color_plan.truncated,
            idenso_required=color_plan.idenso_required,
        )
    color_engine = ColorEngine(color_plan, active_model)
    full_mask = _labels_mask(leg.label for leg in process_ir.legs)
    closure_candidate_splits = _closure_candidate_splits(
        process_ir,
        active_model,
        color_engine,
    )
    closure_reachable_masks = (
        _closure_side_reachable_masks(full_mask, closure_candidate_splits)
        if closure_side_mask_pruning
        else None
    )
    color_order_reachable_masks = (
        _lc_color_order_reachable_masks(process_ir, color_plan, active_model)
        if color_order_mask_pruning
        else None
    )
    totals = _closure_total_coupling_orders(
        process_ir,
        active_model,
        color_engine,
        closure_candidate_splits,
        closure_reachable_masks,
        color_order_reachable_masks,
        max_coupling_orders=_normalize_coupling_order_limits(max_coupling_orders),
        ignored_particle_ids=frozenset(
            int(particle_id) for particle_id in (ignored_particle_ids or ())
        ),
        ignored_vertex_kinds=frozenset(
            int(kind) for kind in (ignored_vertex_kinds or ())
        ),
    )
    if not totals:
        return {}
    minimum_total_order = min(_coupling_order_degree(total) for total in totals)
    minimal_totals = tuple(
        total
        for total in totals
        if _coupling_order_degree(total) == minimum_total_order
    )
    return _coupling_order_envelope(minimal_totals)


def prune_dag_to_amplitude_roots(dag: GenericDAG) -> GenericDAG:
    """Drop currents and interactions that cannot feed any amplitude root.

    The forward generic sweep intentionally over-generates valid local currents:
    it does not know which of them will survive closure until the full table is
    built.  Production artifacts should mirror AmpliCol's optimized library
    structure and keep only the backward-reachable sub-DAG feeding retained
    amplitude roots.
    """

    if not dag.amplitude_roots:
        return dag

    interactions_by_result: dict[int, list[InteractionNode]] = {}
    for interaction in dag.interactions:
        interactions_by_result.setdefault(interaction.result_id, []).append(interaction)

    required_current_ids: set[int] = set()
    stack: list[int] = []
    for root in dag.amplitude_roots:
        for current_id in (root.left_id, root.right_id):
            if current_id not in required_current_ids:
                required_current_ids.add(current_id)
                stack.append(current_id)

    required_interaction_ids: set[int] = set()
    while stack:
        current_id = stack.pop()
        for interaction in interactions_by_result.get(current_id, ()):
            if interaction.id in required_interaction_ids:
                continue
            required_interaction_ids.add(interaction.id)
            for parent_id in (interaction.left_id, interaction.right_id):
                if parent_id not in required_current_ids:
                    required_current_ids.add(parent_id)
                    stack.append(parent_id)

    if len(required_current_ids) == len(dag.currents) and len(required_interaction_ids) == len(dag.interactions):
        return dag

    current_id_map = {
        old_id: new_id
        for new_id, old_id in enumerate(sorted(required_current_ids))
    }
    pruned_currents = tuple(
        CurrentNode(
            id=current_id_map[current.id],
            index=current.index,
            dimension=current.dimension,
            is_source=current.is_source,
            source_leg_label=current.source_leg_label,
            source_helicity=current.source_helicity,
        )
        for current in dag.currents
        if current.id in required_current_ids
    )
    interaction_id_map = {
        old_id: new_id
        for new_id, old_id in enumerate(sorted(required_interaction_ids))
    }
    pruned_interactions = tuple(
        InteractionNode(
            id=interaction_id_map[interaction.id],
            vertex_kind=interaction.vertex_kind,
            vertex_particles=interaction.vertex_particles,
            left_id=current_id_map[interaction.left_id],
            right_id=current_id_map[interaction.right_id],
            result_id=current_id_map[interaction.result_id],
            coupling=interaction.coupling,
            color_weight=interaction.color_weight,
            lowering_backend=interaction.lowering_backend,
            full_tensor_network_ready=interaction.full_tensor_network_ready,
        )
        for interaction in dag.interactions
        if interaction.id in required_interaction_ids
    )
    pruned_roots = tuple(
        AmplitudeRoot(
            id=new_id,
            kind=root.kind,
            left_id=current_id_map[root.left_id],
            right_id=current_id_map[root.right_id],
            color_weight=root.color_weight,
            vertex_kind=root.vertex_kind,
            vertex_particles=root.vertex_particles,
            coupling=root.coupling,
            contraction=root.contraction,
            helicity_weight=root.helicity_weight,
        )
        for new_id, root in enumerate(dag.amplitude_roots)
    )
    pruned_sources = tuple(
        current_id_map[source_id]
        for source_id in dag.sources
        if source_id in required_current_ids
    )
    return GenericDAG(
        process=dag.process,
        color_plan=dag.color_plan,
        currents=pruned_currents,
        sources=pruned_sources,
        interactions=pruned_interactions,
        amplitude_roots=pruned_roots,
        truncated=dag.truncated,
    )


def prune_global_helicity_flip_equivalent_roots(
    dag: GenericDAG,
    model: Model,
) -> GenericDAG:
    """Group global-helicity-flip equivalent roots for parity-symmetric QCD.

    This is the safe, structural subset of AmpliCol's numerical helicity
    filtering.  For pure-QCD LC amplitudes, flipping all external helicities is
    a parity-equivalent contribution to the helicity-summed squared matrix
    element.  Keeping one representative with doubled helicity weight reduces
    amplitude roots and then dead-tree pruning removes currents that fed only
    the discarded partner roots.
    """

    if not _global_helicity_flip_equivalence_safe(dag, model):
        return dag
    if not dag.amplitude_roots:
        return dag

    source_by_bit = _source_helicity_signature_by_bit(dag)
    pure_gluon = _pure_gluon_tree_helicity_pruning_safe(dag, model)
    initial_leg_labels = {
        int(leg.label) for leg in dag.process.legs if leg.side == "initial"
    }
    roots_by_signature: dict[tuple[object, ...], list[AmplitudeRoot]] = {}
    zero_pruned = False
    for root in dag.amplitude_roots:
        signature = _root_physical_helicity_signature(dag, root, source_by_bit)
        if pure_gluon and _pure_gluon_tree_helicity_signature_is_zero(
            signature,
            initial_leg_labels,
        ):
            zero_pruned = True
            continue
        roots_by_signature.setdefault(
            signature,
            [],
        ).append(root)
    if not roots_by_signature:
        return dag

    handled: set[tuple[object, ...]] = set()
    retained: list[AmplitudeRoot] = []
    changed = False
    for signature in sorted(roots_by_signature):
        if signature in handled:
            continue
        flipped = _flip_root_physical_helicity_signature(signature)
        partner = roots_by_signature.get(flipped)
        handled.add(signature)
        if partner is not None:
            handled.add(flipped)
        weight = 1.0
        if partner is not None and flipped != signature:
            weight = 2.0
            changed = True
        for root in roots_by_signature[signature]:
            retained.append(
                AmplitudeRoot(
                    id=len(retained),
                    kind=root.kind,
                    left_id=root.left_id,
                    right_id=root.right_id,
                    color_weight=root.color_weight,
                    vertex_kind=root.vertex_kind,
                    vertex_particles=root.vertex_particles,
                    coupling=root.coupling,
                    contraction=root.contraction,
                    helicity_weight=root.helicity_weight * weight,
                )
            )

    if not changed and not zero_pruned:
        return dag
    return prune_dag_to_amplitude_roots(
        GenericDAG(
            process=dag.process,
            color_plan=dag.color_plan,
            currents=dag.currents,
            sources=dag.sources,
            interactions=dag.interactions,
            amplitude_roots=tuple(retained),
            truncated=dag.truncated,
        )
    )


def _global_helicity_flip_equivalence_safe(
    dag: GenericDAG,
    model: Model,
) -> bool:
    if dag.process.color_accuracy != "lc":
        return False
    for leg in dag.process.legs:
        if leg.outgoing_pdg is None:
            return False
        pdg = int(leg.outgoing_pdg)
        if not (model.is_gluon(pdg) or model.is_quark(abs(pdg))):
            return False
        if model.mass(pdg) != 0.0:
            return False
    for interaction in dag.interactions:
        orders = model.vertex_coupling_orders(
            Vertex(interaction.vertex_kind, interaction.vertex_particles)
        )
        if not orders or any(name != "QCD" for name, _value in orders):
            return False
    for root in dag.amplitude_roots:
        if root.vertex_kind is None or root.vertex_particles is None:
            continue
        orders = model.vertex_coupling_orders(
            Vertex(root.vertex_kind, root.vertex_particles)
        )
        if not orders or any(name != "QCD" for name, _value in orders):
            return False
    return True


def _pure_gluon_tree_helicity_pruning_safe(
    dag: GenericDAG,
    model: Model,
) -> bool:
    return all(
        leg.outgoing_pdg is not None and model.is_gluon(int(leg.outgoing_pdg))
        for leg in dag.process.legs
    )


def _pure_gluon_tree_helicity_signature_is_zero(
    signature: tuple[object, ...],
    initial_leg_labels: set[int],
) -> bool:
    _sector_id, source_helicities = signature
    helicities = []
    for label, helicity in cast(Sequence[tuple[int, int]], source_helicities):
        value = int(helicity)
        if int(label) in initial_leg_labels:
            value = -value
        helicities.append(value)
    return helicities.count(1) < 2 or helicities.count(-1) < 2


def _root_physical_helicity_signature(
    dag: GenericDAG,
    root: AmplitudeRoot,
    source_by_bit: Mapping[int, tuple[int, int]],
) -> tuple[object, ...]:
    left = dag.currents[root.left_id].index
    right = dag.currents[root.right_id].index
    ancestry = int(left.helicity_ancestry | right.helicity_ancestry)
    source_helicities = tuple(
        sorted(
            source
            for bit, source in source_by_bit.items()
            if ancestry & bit
        )
    )
    sector_id = (
        int(left.color_state.sector_id)
        if left.color_state.accuracy == "lc"
        else 0
    )
    return (sector_id, source_helicities)


def _source_helicity_signature_by_bit(
    dag: GenericDAG,
) -> dict[int, tuple[int, int]]:
    source_by_bit: dict[int, tuple[int, int]] = {}
    for current in dag.currents:
        if not current.is_source:
            continue
        source_by_bit[int(current.index.helicity_ancestry)] = (
            int(current.source_leg_label or 0),
            int(current.source_helicity or 0),
        )
    return source_by_bit


def _flip_root_physical_helicity_signature(
    signature: tuple[object, ...],
) -> tuple[object, ...]:
    sector_id, source_helicities = signature
    return (
        sector_id,
        tuple(
            sorted(
                (int(label), -int(helicity))
                for label, helicity in cast(
                    Sequence[tuple[int, int]],
                    source_helicities,
                )
            )
        ),
    )


def _normalize_coupling_order_limits(
    limits: Mapping[str, int] | None,
) -> dict[str, int]:
    if limits is None:
        return {}
    return {
        str(name).upper(): int(value)
        for name, value in limits.items()
        if int(value) >= 0
    }


def _coupling_orders_within_limits(
    orders: tuple[tuple[str, int], ...],
    limits: Mapping[str, int],
) -> bool:
    if not limits:
        return True
    order_map = {str(name).upper(): int(value) for name, value in orders}
    return all(order_map.get(name, 0) <= int(limit) for name, limit in limits.items())


def _combine_coupling_order_tuples(
    *orders: CouplingOrders,
) -> CouplingOrders:
    totals: dict[str, int] = {}
    for order_tuple in orders:
        for name, value in order_tuple:
            totals[str(name).upper()] = totals.get(str(name).upper(), 0) + int(value)
    return tuple(sorted((name, value) for name, value in totals.items() if value))


def _lc_line_groups_within_limit(
    color_state: ColorState,
    limit: int | None,
) -> bool:
    if limit is None or color_state.accuracy != "lc":
        return True
    return len(color_state.line_groups) <= limit


def _right_particles_by_left(
    model: Model,
    *,
    color_accuracy: str,
) -> dict[int, tuple[int, ...]]:
    rights: dict[int, set[int]] = {}
    for vertex in model.iter_vertices(color_accuracy=color_accuracy):
        rights.setdefault(vertex.particles[0], set()).add(vertex.particles[1])
    return {
        left: tuple(sorted(right_particles))
        for left, right_particles in rights.items()
    }


def contributing_color_sector_ids(dag: GenericDAG) -> tuple[int, ...]:
    """Return LC colour sectors that actually contribute amplitude roots."""

    return tuple(
        sorted(
            {
                dag.currents[root.left_id].index.color_state.sector_id
                for root in dag.amplitude_roots
            }
        )
    )


def filter_dag_to_color_sectors(
    dag: GenericDAG,
    sector_ids: Iterable[int],
) -> GenericDAG:
    """Return a dense-current DAG restricted to the requested colour sectors.

    Full DAG construction remains useful for diagnostics, but Rusticol schema-v2
    expects dense current ids.  This helper derives the runtime DAG by keeping
    only currents, interactions, sources, and roots whose LC colour sector is in
    ``sector_ids`` and remapping current/root/interaction ids densely.
    """

    selected = set(sector_ids)
    if not selected:
        return GenericDAG(
            process=dag.process,
            color_plan=dag.color_plan,
            currents=(),
            sources=(),
            interactions=(),
            amplitude_roots=(),
            truncated=dag.truncated,
        )
    current_id_map: dict[int, int] = {}
    currents: list[CurrentNode] = []
    for current in dag.currents:
        if current.index.color_state.sector_id not in selected:
            continue
        new_id = len(currents)
        current_id_map[current.id] = new_id
        currents.append(
            CurrentNode(
                id=new_id,
                index=current.index,
                dimension=current.dimension,
                is_source=current.is_source,
                source_leg_label=current.source_leg_label,
                source_helicity=current.source_helicity,
            )
        )

    sources = tuple(
        current_id_map[source_id]
        for source_id in dag.sources
        if source_id in current_id_map
    )

    interactions: list[InteractionNode] = []
    for interaction in dag.interactions:
        if (
            interaction.left_id not in current_id_map
            or interaction.right_id not in current_id_map
            or interaction.result_id not in current_id_map
        ):
            continue
        interactions.append(
            InteractionNode(
                id=len(interactions),
                vertex_kind=interaction.vertex_kind,
                vertex_particles=interaction.vertex_particles,
                left_id=current_id_map[interaction.left_id],
                right_id=current_id_map[interaction.right_id],
                result_id=current_id_map[interaction.result_id],
                coupling=interaction.coupling,
                color_weight=interaction.color_weight,
                lowering_backend=interaction.lowering_backend,
                full_tensor_network_ready=interaction.full_tensor_network_ready,
            )
        )

    amplitude_roots: list[AmplitudeRoot] = []
    for root in dag.amplitude_roots:
        if root.left_id not in current_id_map or root.right_id not in current_id_map:
            continue
        amplitude_roots.append(
            AmplitudeRoot(
                id=len(amplitude_roots),
                kind=root.kind,
                left_id=current_id_map[root.left_id],
                right_id=current_id_map[root.right_id],
                color_weight=root.color_weight,
                vertex_kind=root.vertex_kind,
                vertex_particles=root.vertex_particles,
                coupling=root.coupling,
                contraction=root.contraction,
                helicity_weight=root.helicity_weight,
            )
        )

    return GenericDAG(
        process=dag.process,
        color_plan=dag.color_plan,
        currents=tuple(currents),
        sources=sources,
        interactions=tuple(interactions),
        amplitude_roots=tuple(amplitude_roots),
        truncated=dag.truncated,
    )


def _direct_contraction_kind(
    model: Model,
    left: CurrentIndex,
    right: CurrentIndex,
) -> str | None:
    if model.anti_particle(left.particle_id) != right.particle_id:
        return None
    left_dimension = model.current_dimension(left.particle_id, left.chirality)
    right_dimension = model.current_dimension(right.particle_id, right.chirality)
    if left_dimension != right_dimension:
        return None
    if left_dimension == 1:
        return "scalar"
    if left_dimension == 2:
        if left.chirality != -right.chirality:
            return None
        return "weyl"
    if left_dimension == 4:
        if model.is_fermion(left.particle_id) and model.is_fermion(right.particle_id):
            return "dirac"
        return "lorentz"
    if left_dimension == 6:
        return "antisymmetric-tensor"
    return None


def _closure_contraction_name(model: Model, particle_id: int) -> str:
    dimension = model.current_dimension(particle_id, 0)
    if dimension == 1:
        return "scalar"
    if dimension == 2:
        return "weyl"
    if dimension == 4:
        return "lorentz"
    if dimension == 6:
        return "antisymmetric-tensor"
    return "model-vertex"


def _sector_group_indices_for_label(
    sector: LCColorSector,
    label: int,
) -> tuple[int, ...]:
    return tuple(
        index
        for index, group in enumerate(sector.line_label_groups)
        if label in set(group)
    )


def _line_positions(labels: Iterable[int], order: tuple[int, ...]) -> tuple[int, ...]:
    positions = {label: index for index, label in enumerate(order)}
    return tuple(positions[label] for label in labels if label in positions)


def _labels_projected_to_word(
    labels: Iterable[int],
    word: tuple[int, ...],
) -> tuple[int, ...]:
    word_labels = set(word)
    return tuple(label for label in labels if label in word_labels)


def _positions_contiguous(positions: tuple[int, ...]) -> bool:
    if not positions:
        return True
    return positions == tuple(range(positions[0], positions[-1] + 1))


def _ordered_combination_matches_word(
    left_labels: Iterable[int],
    right_labels: Iterable[int],
    word: tuple[int, ...],
) -> bool:
    return _ordered_combination_segment(left_labels, right_labels, word) is not None


def _closure_combination_matches_word(
    left_labels: Iterable[int],
    right_labels: Iterable[int],
    word: tuple[int, ...],
) -> bool:
    """Return whether two endpoint label sets close one colour word.

    Intermediate current construction must preserve the actual line order because
    it fixes where colour-singlet insertions may attach.  Final closure is
    different: the sink current can be oriented opposite to the sector's primary
    colour word, while still representing the same contiguous endpoint segment.
    Closure therefore checks contiguity of the projected coloured label sets, not
    their internal iteration order.
    """

    left_positions = tuple(sorted(_line_positions(left_labels, word)))
    right_positions = tuple(sorted(_line_positions(right_labels, word)))
    if not left_positions or not right_positions:
        return False
    if set(left_positions) & set(right_positions):
        return False
    union = tuple(sorted((*left_positions, *right_positions)))
    if (
        _positions_contiguous(left_positions)
        and _positions_contiguous(right_positions)
        and _positions_contiguous(union)
    ):
        return True
    return False


def _ordered_combination_segment(
    left_labels: Iterable[int],
    right_labels: Iterable[int],
    word: tuple[int, ...],
) -> tuple[int, ...] | None:
    left_positions = _line_positions(left_labels, word)
    right_positions = _line_positions(right_labels, word)
    if not left_positions and not right_positions:
        return ()
    if not left_positions or not right_positions:
        positions = left_positions or right_positions
        if not _positions_contiguous(positions):
            return None
        return tuple(word[positions[0] : positions[-1] + 1])
    union = tuple(sorted((*left_positions, *right_positions)))
    if (
        _positions_contiguous(left_positions)
        and _positions_contiguous(right_positions)
        and _positions_contiguous(union)
        and max(left_positions) < min(right_positions)
    ):
        return tuple(word[union[0] : union[-1] + 1])
    return None


def _canonical_lc_ordered_labels(
    labels: Iterable[int],
    sector: LCColorSector,
) -> tuple[int, ...]:
    label_tuple = tuple(labels)
    label_set = set(label_tuple)
    for word in _sector_intermediate_order_words(sector):
        positions = _line_positions(label_tuple, word)
        if not positions:
            continue
        if not _positions_contiguous(positions):
            continue
        word_labels = set(word)
        extras = tuple(sorted(label for label in label_set if label not in word_labels))
        segment = tuple(word[positions[0] : positions[-1] + 1])
        if not _line_local_singlet_extras_allowed(segment, extras, sector):
            continue
        return (*segment, *extras)
    return tuple(sorted(label_set))


def _sector_current_order_words(sector: LCColorSector) -> tuple[tuple[int, ...], ...]:
    """Return physical order words carried in diagnostics.

    Current construction uses ``sector.compatibility_words`` plus
    ``_line_local_singlet_extras_allowed``.  Keeping this helper as a diagnostic
    accessor makes it explicit that the full legacy words are not the pruning
    substrate: singlet insertions are line-local attachments to coloured
    segments, not fixed positions at the tail of a full word.
    """

    words = tuple(getattr(sector, "legacy_order_words", ()) or ())
    if words:
        return words
    return sector.compatibility_words or sector.color_words


def _sector_intermediate_order_words(
    sector: LCColorSector,
) -> tuple[tuple[int, ...], ...]:
    """Return LC words used for intermediate current construction.

    Multi-open-line sectors expose extra compatibility words so final closure
    can choose the opposite endpoint without duplicating physical sectors.
    Intermediate currents, however, must follow the sector's physical colour
    word; using all compatibility words here double-counts the same ordered
    current topology for multi-quark-line processes with singlet insertions.
    """

    return sector.color_words or sector.compatibility_words


def _line_local_singlet_extras_allowed(
    colored_segment: Iterable[int],
    extras: Iterable[int],
    sector: LCColorSector,
) -> bool:
    extra_set = set(extras)
    if not extra_set:
        return True
    if sector.kind != "open-lines":
        return extra_set.issubset(set(sector.singlet_labels))
    return extra_set.issubset(set(sector.singlet_labels))


def _labels_mask(labels: Iterable[int]) -> int:
    mask = 0
    for label in labels:
        mask |= 1 << (label - 1)
    return mask


def _mask_labels(mask: int) -> tuple[int, ...]:
    return tuple(index + 1 for index in range(mask.bit_length()) if mask & (1 << index))


def _canonical_sink_mask(process_ir: CanonicalProcessIR) -> int:
    fermion_classes = {"quark", "antiquark", "charged-lepton", "neutrino"}
    fermion_labels = [
        leg.label
        for leg in process_ir.legs
        if leg.outgoing_pdg is not None and leg.particle_class in fermion_classes
    ]
    if fermion_labels:
        return 1 << (min(fermion_labels) - 1)
    labels = [leg.label for leg in process_ir.legs if leg.outgoing_pdg is not None]
    if not labels:
        return 0
    return 1 << (min(labels) - 1)


def _closure_candidate_splits(
    process_ir: CanonicalProcessIR,
    model: Model,
    color_engine: ColorEngine,
    *,
    reference_color_order: tuple[int, ...] | None = None,
) -> tuple[tuple[int, int], ...]:
    full_mask = _labels_mask(leg.label for leg in process_ir.legs)
    splits: list[tuple[int, int]] = []
    split_seen: set[tuple[int, int]] = set()
    sink_labels: list[int] = []
    if reference_color_order:
        leg_by_label = {leg.label: leg for leg in process_ir.legs}
        colored_reference_labels: list[int] = []
        ordered_reference_labels: list[int] = []
        for raw_label in reference_color_order:
            label = int(raw_label)
            if not (full_mask & (1 << (label - 1))):
                continue
            ordered_reference_labels.append(label)
            leg = leg_by_label.get(label)
            if leg is None or leg.outgoing_pdg is None:
                continue
            try:
                is_colored = model.color_rep(int(leg.outgoing_pdg)) != 1
            except KeyError:
                is_colored = True
            if is_colored:
                colored_reference_labels.append(label)
        if colored_reference_labels:
            sink_labels.append(colored_reference_labels[-1])
        elif ordered_reference_labels:
            sink_labels.append(ordered_reference_labels[-1])
    if not sink_labels and color_engine.color_plan.color_accuracy == "lc":
        for sector in color_engine.color_plan.sectors:
            # Compatibility words are used while building intermediate currents:
            # complete open-line blocks may be traversed in several orders
            # without changing the physical LC sector.  Final amplitude
            # closure, however, must choose one physical sink per sector;
            # otherwise multi-open-line sectors are counted once per
            # compatible block ordering.
            for word in sector.color_words:
                if word:
                    sink_labels.append(int(word[-1]))
    if not sink_labels:
        sink_labels.append(_mask_labels(_canonical_sink_mask(process_ir))[0])
    seen: set[int] = set()
    for label in sink_labels:
        if label in seen:
            continue
        seen.add(label)
        sink_mask = 1 << (label - 1)
        if not (sink_mask & full_mask):
            continue
        split = (full_mask ^ sink_mask, sink_mask)
        if split in split_seen:
            continue
        split_seen.add(split)
        splits.append(split)
    return tuple(splits)


def _closure_side_reachable_masks(
    full_mask: int,
    candidate_splits: Iterable[tuple[int, int]],
) -> frozenset[int]:
    """Return current masks that can feed one configured amplitude closure.

    The generic forward sweep otherwise builds every locally valid current for
    every proper subset, then removes dead currents after closure.  AmpliCol's
    library generation is faster because it knows which endpoint side each
    current can ultimately feed.  This helper applies the same idea without
    naming a process family: once the generic closure splitter has selected the
    possible amplitude endpoint masks, any useful intermediate current must be
    a submask of one of those endpoint masks.
    """

    allowed: set[int] = set()
    for left_mask, right_mask in candidate_splits:
        for side_mask in (left_mask, right_mask):
            submask = side_mask & full_mask
            while submask:
                allowed.add(submask)
                submask = (submask - 1) & side_mask
    return frozenset(allowed)


def _lc_color_order_reachable_masks(
    process_ir: CanonicalProcessIR,
    color_plan: GenericColorPlan,
    model: Model,
) -> frozenset[int] | None:
    """Return subset masks compatible with at least one LC colour word.

    This is process-generic pruning.  It only uses model colour
    representations and the LC words produced by the colour planner.  Any
    useful coloured current in a colour-ordered recursion must cover a
    contiguous segment of one compatibility word.  Colour singlets are left as
    attachments to those segments because their allowed positions are governed
    by ordinary model vertices and the existing singlet-order rule.
    """

    if color_plan.color_accuracy != "lc" or not color_plan.sectors:
        return None

    full_mask = _labels_mask(leg.label for leg in process_ir.legs)
    colored_labels: set[int] = set()
    singlet_labels: set[int] = set()
    for leg in process_ir.legs:
        if leg.outgoing_pdg is None:
            continue
        try:
            is_colored = model.color_rep(int(leg.outgoing_pdg)) != 1
        except KeyError:
            is_colored = True
        if is_colored:
            colored_labels.add(leg.label)
        else:
            singlet_labels.add(leg.label)

    allowed: set[int] = set()

    if not colored_labels:
        return frozenset(_nonzero_submasks(full_mask))

    for sector in color_plan.sectors:
        if sector.kind == "open-lines":
            for singlet_submask in _submasks_for_labels(sector.singlet_labels):
                if singlet_submask:
                    allowed.add(singlet_submask)
        for raw_word in sector.compatibility_words:
            word = tuple(
                label
                for label in raw_word
                if label in colored_labels
            )
            for start in range(len(word)):
                segment_mask = 0
                for stop in range(start, len(word)):
                    segment_mask |= 1 << (word[stop] - 1)
                    allowed.add(segment_mask)
                    segment_labels = tuple(word[start : stop + 1])
                    line_singlets = _line_local_singlet_labels_for_segment(
                        sector,
                        segment_labels,
                    )
                    for singlet_submask in _submasks_for_labels(line_singlets):
                        allowed.add(segment_mask | singlet_submask)

    for label in colored_labels:
        allowed.add(1 << (label - 1))
    for label in singlet_labels:
        allowed.add(1 << (label - 1))
    return frozenset(mask for mask in allowed if mask & full_mask)


def _line_local_singlet_labels_for_segment(
    sector: LCColorSector,
    colored_segment: Iterable[int],
) -> tuple[int, ...]:
    if sector.kind != "open-lines":
        return sector.singlet_labels
    return sector.singlet_labels


def _submasks_for_labels(labels: Iterable[int]) -> tuple[int, ...]:
    label_tuple = tuple(labels)
    masks: list[int] = []
    count = len(label_tuple)
    for bits in range(1 << count):
        mask = 0
        for index, label in enumerate(label_tuple):
            if bits & (1 << index):
                mask |= 1 << (label - 1)
        masks.append(mask)
    return tuple(masks)


UsefulStateMap = dict[int, dict[int, frozenset[CouplingOrders]]]


def _useful_states_by_mask(
    process_ir: CanonicalProcessIR,
    model: Model,
    color_engine: ColorEngine,
    closure_candidate_splits: Iterable[tuple[int, int]],
    closure_reachable_masks: frozenset[int] | None,
    color_order_reachable_masks: frozenset[int] | None,
    *,
    max_coupling_orders: Mapping[str, int],
    ignored_particle_ids: frozenset[int],
    ignored_vertex_kinds: frozenset[int],
) -> UsefulStateMap:
    """Return current states that can feed at least one amplitude closure.

    The full current table is keyed by helicity, chirality, flavour flow and
    colour state.  This prepass intentionally ignores those expensive labels,
    but it keeps particle id and model coupling-order totals.  The order
    tracking lets user-supplied generic constraints such as ``QED=1`` prune
    impossible branches before the expensive helicity/current sweep, without
    recognizing any process family.  It is allowed to overestimate, but it must
    never depend on a process-family name.
    """

    full_mask = _labels_mask(leg.label for leg in process_ir.legs)
    possible: dict[int, dict[int, set[CouplingOrders]]] = {}
    for leg in process_ir.legs:
        if leg.outgoing_pdg is None:
            continue
        particle_id = int(leg.outgoing_pdg)
        if particle_id in ignored_particle_ids:
            continue
        possible.setdefault(1 << (leg.label - 1), {}).setdefault(
            particle_id,
            set(),
        ).add(())

    transitions: list[
        tuple[
            int,
            int,
            CouplingOrders,
            int,
            int,
            CouplingOrders,
            int,
            int,
            CouplingOrders,
        ]
    ] = []
    vertices_by_input: dict[tuple[int, int], tuple[Vertex, ...]] = {}
    for mask in _masks_by_size(full_mask):
        if mask & (mask - 1) == 0 or mask == full_mask:
            continue
        if not _mask_allowed_by_reachability(
            mask,
            closure_reachable_masks,
            color_order_reachable_masks,
        ):
            continue
        for left_mask, right_mask in _ordered_splits(mask):
            if not (
                _mask_allowed_by_reachability(
                    left_mask,
                    closure_reachable_masks,
                    color_order_reachable_masks,
                )
                and _mask_allowed_by_reachability(
                    right_mask,
                    closure_reachable_masks,
                    color_order_reachable_masks,
                )
            ):
                continue
            left_species = possible.get(left_mask)
            right_species = possible.get(right_mask)
            if not left_species or not right_species:
                continue
            for left_particle, left_orders_set in tuple(left_species.items()):
                for right_particle, right_orders_set in tuple(right_species.items()):
                    vertex_key = (left_particle, right_particle)
                    vertices = vertices_by_input.get(vertex_key)
                    if vertices is None:
                        vertices = model.vertices_accepting(
                            left_particle,
                            right_particle,
                            color_accuracy=process_ir.color_accuracy,
                        )
                        vertices_by_input[vertex_key] = vertices
                    for vertex in vertices:
                        if (
                            vertex.kind in ignored_vertex_kinds
                            or vertex.particles[2] in ignored_particle_ids
                            or not color_engine.vertex_allowed(vertex)
                            or model.skip_duplicate_vertex_orientation(vertex)
                        ):
                            continue
                        result_particle = vertex.particles[2]
                        for left_orders in tuple(left_orders_set):
                            for right_orders in tuple(right_orders_set):
                                coupling_orders = _combine_coupling_order_tuples(
                                    left_orders,
                                    right_orders,
                                    model.vertex_coupling_orders(vertex),
                                )
                                if not _coupling_orders_within_limits(
                                    coupling_orders,
                                    max_coupling_orders,
                                ):
                                    continue
                                possible.setdefault(mask, {}).setdefault(
                                    result_particle,
                                    set(),
                                ).add(coupling_orders)
                                transitions.append(
                                    (
                                        left_mask,
                                        left_particle,
                                        left_orders,
                                        right_mask,
                                        right_particle,
                                        right_orders,
                                        mask,
                                        result_particle,
                                        coupling_orders,
                                    )
                            )

    useful: dict[int, dict[int, set[CouplingOrders]]] = {}
    for left_mask, right_mask in closure_candidate_splits:
        left_species = possible.get(left_mask)
        right_species = possible.get(right_mask)
        if not left_species or not right_species:
            continue
        for left_particle, left_orders_set in left_species.items():
            for right_particle, right_orders_set in right_species.items():
                if _species_direct_contraction_kind(model, left_particle, right_particle):
                    for left_orders in left_orders_set:
                        for right_orders in right_orders_set:
                            total_orders = _combine_coupling_order_tuples(
                                left_orders,
                                right_orders,
                                (),
                            )
                            if not _coupling_orders_within_limits(
                                total_orders,
                                max_coupling_orders,
                            ):
                                continue
                            useful.setdefault(left_mask, {}).setdefault(
                                left_particle,
                                set(),
                            ).add(left_orders)
                            useful.setdefault(right_mask, {}).setdefault(
                                right_particle,
                                set(),
                            ).add(right_orders)
                for vertex in model.vertices_accepting(
                    left_particle,
                    right_particle,
                    color_accuracy=process_ir.color_accuracy,
                ):
                    if (
                        vertex.kind in ignored_vertex_kinds
                        or vertex.particles[2] in ignored_particle_ids
                        or not color_engine.vertex_allowed(vertex)
                        or model.skip_duplicate_vertex_orientation(vertex)
                    ):
                        continue
                    if _closure_contraction_name(model, vertex.particles[2]) != "scalar":
                        continue
                    vertex_orders = model.vertex_coupling_orders(vertex)
                    for left_orders in left_orders_set:
                        for right_orders in right_orders_set:
                            total_orders = _combine_coupling_order_tuples(
                                left_orders,
                                right_orders,
                                vertex_orders,
                            )
                            if not _coupling_orders_within_limits(
                                total_orders,
                                max_coupling_orders,
                            ):
                                continue
                            useful.setdefault(left_mask, {}).setdefault(
                                left_particle,
                                set(),
                            ).add(left_orders)
                            useful.setdefault(right_mask, {}).setdefault(
                                right_particle,
                                set(),
                            ).add(right_orders)

    changed = True
    while changed:
        changed = False
        for (
            left_mask,
            left_particle,
            left_orders,
            right_mask,
            right_particle,
            right_orders,
            mask,
            result_particle,
            result_orders,
        ) in reversed(transitions):
            if result_orders not in useful.get(mask, {}).get(result_particle, ()):
                continue
            left_useful = useful.setdefault(left_mask, {}).setdefault(
                left_particle,
                set(),
            )
            if left_orders not in left_useful:
                left_useful.add(left_orders)
                changed = True
            right_useful = useful.setdefault(right_mask, {}).setdefault(
                right_particle,
                set(),
            )
            if right_orders not in right_useful:
                right_useful.add(right_orders)
                changed = True

    return {
        mask: {
            particle: frozenset(orders)
            for particle, orders in species_orders.items()
            if orders
        }
        for mask, species_orders in useful.items()
        if species_orders
    }


def _closure_total_coupling_orders(
    process_ir: CanonicalProcessIR,
    model: Model,
    color_engine: ColorEngine,
    closure_candidate_splits: Iterable[tuple[int, int]],
    closure_reachable_masks: frozenset[int] | None,
    color_order_reachable_masks: frozenset[int] | None,
    *,
    max_coupling_orders: Mapping[str, int],
    ignored_particle_ids: frozenset[int],
    ignored_vertex_kinds: frozenset[int],
) -> frozenset[CouplingOrders]:
    """Return model-reachable total coupling orders for amplitude closures."""

    full_mask = _labels_mask(leg.label for leg in process_ir.legs)
    possible: dict[int, dict[int, set[CouplingOrders]]] = {}
    for leg in process_ir.legs:
        if leg.outgoing_pdg is None:
            continue
        particle_id = int(leg.outgoing_pdg)
        if particle_id in ignored_particle_ids:
            continue
        possible.setdefault(1 << (leg.label - 1), {}).setdefault(
            particle_id,
            set(),
        ).add(())

    vertices_by_input: dict[tuple[int, int], tuple[Vertex, ...]] = {}
    for mask in _masks_by_size(full_mask):
        if mask & (mask - 1) == 0 or mask == full_mask:
            continue
        if not _mask_allowed_by_reachability(
            mask,
            closure_reachable_masks,
            color_order_reachable_masks,
        ):
            continue
        for left_mask, right_mask in _ordered_splits(mask):
            if not (
                _mask_allowed_by_reachability(
                    left_mask,
                    closure_reachable_masks,
                    color_order_reachable_masks,
                )
                and _mask_allowed_by_reachability(
                    right_mask,
                    closure_reachable_masks,
                    color_order_reachable_masks,
                )
            ):
                continue
            left_species = possible.get(left_mask)
            right_species = possible.get(right_mask)
            if not left_species or not right_species:
                continue
            for left_particle, left_orders_set in tuple(left_species.items()):
                for right_particle, right_orders_set in tuple(right_species.items()):
                    vertex_key = (left_particle, right_particle)
                    vertices = vertices_by_input.get(vertex_key)
                    if vertices is None:
                        vertices = model.vertices_accepting(
                            left_particle,
                            right_particle,
                            color_accuracy=process_ir.color_accuracy,
                        )
                        vertices_by_input[vertex_key] = vertices
                    for vertex in vertices:
                        if (
                            vertex.kind in ignored_vertex_kinds
                            or vertex.particles[2] in ignored_particle_ids
                            or not color_engine.vertex_allowed(vertex)
                            or model.skip_duplicate_vertex_orientation(vertex)
                        ):
                            continue
                        result_particle = vertex.particles[2]
                        order_bucket = possible.setdefault(mask, {}).setdefault(
                            result_particle,
                            set(),
                        )
                        vertex_orders = model.vertex_coupling_orders(vertex)
                        for left_orders in tuple(left_orders_set):
                            for right_orders in tuple(right_orders_set):
                                coupling_orders = _combine_coupling_order_tuples(
                                    left_orders,
                                    right_orders,
                                    vertex_orders,
                                )
                                if not _coupling_orders_within_limits(
                                    coupling_orders,
                                    max_coupling_orders,
                                ):
                                    continue
                                order_bucket.add(coupling_orders)
                        possible[mask][result_particle] = set(
                            _pareto_minimal_coupling_orders(order_bucket)
                        )

    totals: set[CouplingOrders] = set()
    for left_mask, right_mask in closure_candidate_splits:
        left_species = possible.get(left_mask)
        right_species = possible.get(right_mask)
        if not left_species or not right_species:
            continue
        for left_particle, left_orders_set in left_species.items():
            for right_particle, right_orders_set in right_species.items():
                if _species_direct_contraction_kind(model, left_particle, right_particle):
                    for left_orders in left_orders_set:
                        for right_orders in right_orders_set:
                            total_orders = _combine_coupling_order_tuples(
                                left_orders,
                                right_orders,
                                (),
                            )
                            if _coupling_orders_within_limits(
                                total_orders,
                                max_coupling_orders,
                            ):
                                totals.add(total_orders)
                for vertex in model.vertices_accepting(
                    left_particle,
                    right_particle,
                    color_accuracy=process_ir.color_accuracy,
                ):
                    if (
                        vertex.kind in ignored_vertex_kinds
                        or vertex.particles[2] in ignored_particle_ids
                        or not color_engine.vertex_allowed(vertex)
                        or model.skip_duplicate_vertex_orientation(vertex)
                    ):
                        continue
                    if _closure_contraction_name(model, vertex.particles[2]) != "scalar":
                        continue
                    vertex_orders = model.vertex_coupling_orders(vertex)
                    for left_orders in left_orders_set:
                        for right_orders in right_orders_set:
                            total_orders = _combine_coupling_order_tuples(
                                left_orders,
                                right_orders,
                                vertex_orders,
                            )
                            if _coupling_orders_within_limits(
                                total_orders,
                                max_coupling_orders,
                            ):
                                totals.add(total_orders)
    return frozenset(_pareto_minimal_coupling_orders(totals))


def _coupling_order_degree(orders: CouplingOrders) -> int:
    return sum(int(value) for _, value in orders)


def _coupling_order_envelope(orders: Iterable[CouplingOrders]) -> dict[str, int]:
    envelope: dict[str, int] = {}
    for order_tuple in orders:
        for name, value in order_tuple:
            normalized = str(name).upper()
            envelope[normalized] = max(envelope.get(normalized, 0), int(value))
    return envelope


def _pareto_minimal_coupling_orders(
    orders: Iterable[CouplingOrders],
) -> tuple[CouplingOrders, ...]:
    normalized = tuple(sorted(set(orders)))
    minimal: list[CouplingOrders] = []
    for candidate in normalized:
        if any(
            _coupling_orders_dominate(other, candidate)
            for other in normalized
            if other != candidate
        ):
            continue
        minimal.append(candidate)
    return tuple(minimal)


def _coupling_orders_dominate(left: CouplingOrders, right: CouplingOrders) -> bool:
    left_map = dict(left)
    right_map = dict(right)
    keys = set(left_map) | set(right_map)
    return all(left_map.get(key, 0) <= right_map.get(key, 0) for key in keys) and any(
        left_map.get(key, 0) < right_map.get(key, 0) for key in keys
    )


def _species_direct_contraction_kind(
    model: Model,
    left_particle: int,
    right_particle: int,
) -> str | None:
    if model.anti_particle(left_particle) != right_particle:
        return None
    try:
        dimension = model.dimension(left_particle)
        right_dimension = model.dimension(right_particle)
    except KeyError:
        return None
    if dimension != right_dimension:
        return None
    if dimension == 1:
        return "scalar"
    if dimension == 2:
        return "weyl"
    if dimension == 4:
        return "lorentz"
    if dimension == 6:
        return "antisymmetric-tensor"
    return None


def _mask_allowed_by_reachability(
    mask: int,
    closure_reachable_masks: frozenset[int] | None,
    color_order_reachable_masks: frozenset[int] | None,
) -> bool:
    if closure_reachable_masks is not None and mask not in closure_reachable_masks:
        return False
    if color_order_reachable_masks is not None and mask not in color_order_reachable_masks:
        return False
    return True


def _state_allowed_by_reachability(
    useful_states_by_mask: UsefulStateMap,
    mask: int,
    particle_id: int,
    coupling_orders: CouplingOrders,
) -> bool:
    return coupling_orders in useful_states_by_mask.get(mask, {}).get(
        particle_id,
        (),
    )


def _nonzero_submasks(mask: int) -> tuple[int, ...]:
    submasks: list[int] = []
    submask = mask
    while submask:
        submasks.append(submask)
        submask = (submask - 1) & mask
    return tuple(submasks)


def _masks_by_size(full_mask: int) -> tuple[int, ...]:
    masks = []
    submask = full_mask
    while submask:
        masks.append(submask)
        submask = (submask - 1) & full_mask
    return tuple(sorted(masks, key=lambda value: (value.bit_count(), value)))


def _ordered_splits(mask: int) -> tuple[tuple[int, int], ...]:
    splits: list[tuple[int, int]] = []
    left = (mask - 1) & mask
    while left:
        right = mask ^ left
        if right:
            splits.append((left, right))
        left = (left - 1) & mask
    return tuple(splits)


__all__ = [
    "AmplitudeRoot",
    "ColorEngine",
    "ColorFlow",
    "ColorState",
    "CurrentIndex",
    "CurrentNode",
    "GenericDAG",
    "GenericDAGCompiler",
    "InteractionNode",
    "QuantumFlow",
    "compile_generic_dag",
    "contributing_color_sector_ids",
    "filter_dag_to_color_sectors",
    "infer_minimal_coupling_order_limits",
    "prune_global_helicity_flip_equivalent_roots",
]
