"""Physical Standard-Model tree discovery for phase-space maps.

The amplitude generator internally factorises several four-point rules with
auxiliary fields. Those fields are an implementation detail and are
deliberately absent here: diagram discovery uses an explicit catalogue of the
physical cubic and quartic interactions implemented by :mod:`particles`.

Tree states are projected immediately onto the information which affects the
numerical phase-space density. This keeps discovery finite for realistic
processes while retaining every distinct physical propagator/contact map.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import product
from typing import Callable, Iterable, Mapping


BREIT_WIGNER = "breit_wigner"
MASSLESS_POLE = "massless_pole"
MASSIVE_POWER = "massive_power"
FLAT_CONTACT = "flat_contact"

TRANSFORM_CODES = {
    BREIT_WIGNER: 1,
    MASSLESS_POLE: 2,
    MASSIVE_POWER: 3,
    FLAT_CONTACT: 4,
}

PHYSICAL_PDGS = frozenset(
    (*range(1, 7), *range(-6, 0), 21, 22, 23, 24, -24, 25,
     *range(11, 17), *range(-16, -10))
)
AUXILIARY_PDGS = frozenset((-21, -23, 26, -26, 99, 125, 126, 127))


def anti_pdg(pdg: int) -> int:
    """Return the antiparticle of a physical field."""

    if pdg not in PHYSICAL_PDGS:
        raise ValueError(f"auxiliary or unknown PDG {pdg} is not physical")
    return pdg if pdg in (21, 22, 23, 25) else -pdg


@dataclass(frozen=True)
class SMInteractionCatalogue:
    """Implemented physical SM vertices in an all-outgoing convention."""

    cubic: frozenset[tuple[int, int, int]]
    quartic: frozenset[tuple[int, int, int, int]]
    flavour_scheme: int

    def contains(self, particles: Iterable[int]) -> bool:
        vertex = tuple(sorted(particles))
        if len(vertex) == 3:
            return vertex in self.cubic
        if len(vertex) == 4:
            return vertex in self.quartic
        return False

    @property
    def vertices(self) -> frozenset[tuple[int, ...]]:
        return frozenset((*self.cubic, *self.quartic))


def build_sm_interaction_catalogue(flavour_scheme: int) -> SMInteractionCatalogue:
    """Build the physical tree-level interactions implemented by AmpliCol.

    Charged currents are diagonal in the three generations. A quark Yukawa is
    active precisely when that flavour is massive in the selected flavour
    scheme (top is always massive). No auxiliary amplitude field appears in
    either catalogue.
    """

    if not 1 <= flavour_scheme <= 5:
        raise ValueError("flavour scheme must be between 1 and 5")

    cubic: set[tuple[int, int, int]] = set()

    def add3(*particles: int) -> None:
        cubic.add(tuple(sorted(particles)))

    add3(21, 21, 21)
    add3(24, -24, 22)
    add3(24, -24, 23)
    add3(24, -24, 25)
    add3(23, 23, 25)
    add3(25, 25, 25)

    for flavour in range(1, 7):
        add3(21, flavour, -flavour)
        add3(22, flavour, -flavour)
        add3(23, flavour, -flavour)
    for flavour in range(flavour_scheme + 1, 7):
        add3(25, flavour, -flavour)
    for down in (1, 3, 5):
        up = down + 1
        add3(24, down, -up)
        add3(-24, -down, up)
    for charged, neutrino in ((11, 12), (13, 14), (15, 16)):
        add3(22, charged, -charged)
        add3(23, charged, -charged)
        add3(23, neutrino, -neutrino)
        add3(24, charged, -neutrino)
        add3(-24, -charged, neutrino)

    quartic: set[tuple[int, int, int, int]] = set()

    def add4(*particles: int) -> None:
        quartic.add(tuple(sorted(particles)))

    # Four-gluon, physical four-vector, vector-vector-Higgs-Higgs, and
    # four-Higgs contacts represented by auxiliary fields in particles.f03.
    add4(21, 21, 21, 21)
    add4(24, -24, 24, -24)
    add4(24, -24, 22, 22)
    add4(24, -24, 22, 23)
    add4(24, -24, 23, 23)
    add4(24, -24, 25, 25)
    add4(23, 23, 25, 25)
    add4(25, 25, 25, 25)

    catalogue = SMInteractionCatalogue(
        frozenset(cubic), frozenset(quartic), flavour_scheme
    )
    if any(AUXILIARY_PDGS & set(vertex) for vertex in catalogue.vertices):
        raise AssertionError("auxiliary field leaked into physical catalogue")
    return catalogue


def _is_quark(pdg: int) -> bool:
    return 1 <= abs(pdg) <= 6


def is_pure_qcd_vertex(vertex: Iterable[int]) -> bool:
    """Return whether a physical vertex belongs entirely to pure QCD."""

    particles = tuple(vertex)
    if len(particles) == 3:
        return (
            all(pdg == 21 for pdg in particles)
            or (
                particles.count(21) == 1
                and len([pdg for pdg in particles if _is_quark(pdg)]) == 2
                and sum(pdg for pdg in particles if _is_quark(pdg)) == 0
            )
        )
    return len(particles) == 4 and all(pdg == 21 for pdg in particles)


def transform_for_pdg(pdg: int, flavour_scheme: int) -> tuple[str, int]:
    """Return the numerical invariant-mass transform for a physical current.

    The integer parameter identifies numerical masses/widths where required.
    Charge conjugates deliberately share it. A zero parameter means that the
    density is independent of particle identity.
    """

    apdg = abs(pdg)
    if apdg in (6, 23, 24, 25):
        return BREIT_WIGNER, apdg
    if pdg in (21, 22) or 11 <= apdg <= 16 or (
        1 <= apdg <= flavour_scheme
    ):
        return MASSLESS_POLE, 0
    if 1 <= apdg <= 6:
        return MASSIVE_POWER, apdg
    raise ValueError(f"cannot choose a transform for non-physical PDG {pdg}")


@dataclass(frozen=True, order=True)
class TopologyNode:
    """One binary invariant-mass split in a projected phase-space tree."""

    pdg: int
    mask: int
    left: int
    right: int
    kind: str = ""
    parameter: int = 0

    def __post_init__(self) -> None:
        if self.mask <= 0 or self.left <= 0 or self.right <= 0:
            raise ValueError("topology masks must be non-zero")
        if self.left & self.right or self.left | self.right != self.mask:
            raise ValueError("topology children must be a disjoint partition")
        if self.pdg in AUXILIARY_PDGS:
            raise ValueError("auxiliary fields cannot appear in a phase map")
        if self.pdg == 0:
            if self.kind not in ("", FLAT_CONTACT):
                raise ValueError("PDG zero is reserved for flat contacts")
            object.__setattr__(self, "kind", FLAT_CONTACT)
            object.__setattr__(self, "parameter", 0)
        elif self.pdg not in PHYSICAL_PDGS:
            raise ValueError(f"non-physical topology PDG {self.pdg}")
        elif self.kind and self.kind not in TRANSFORM_CODES:
            raise ValueError(f"unknown phase-space transform {self.kind}")

    @property
    def density_signature(self) -> tuple[int, int, int, str, int]:
        return self.mask, self.left, self.right, self.kind, self.parameter


Topology = tuple[TopologyNode, ...]


@dataclass(frozen=True, order=True)
class TreeVertex:
    """One physical vertex used in a completed labelled tree."""

    parent: int
    mask: int
    children: tuple[int, ...]
    child_masks: tuple[int, ...]

    @property
    def particles(self) -> tuple[int, ...]:
        return tuple(sorted((anti_pdg(self.parent), *self.children)))


@dataclass(frozen=True, order=True)
class DiagramChannel:
    """One diagram-derived density in fixed external-leg labels."""

    topology: Topology
    order: tuple[int, ...]
    vertices: tuple[TreeVertex, ...] = ()
    multiplicity: int = 1

    @property
    def density_signature(self) -> tuple[tuple[int, ...], tuple[tuple, ...]]:
        return (
            self.order,
            tuple(node.density_signature for node in self.topology),
        )


def _canonical_node(
    pdg: int,
    mask: int,
    left: int,
    right: int,
    kind: str,
    parameter: int,
) -> TopologyNode:
    if (right.bit_count(), right) < (left.bit_count(), left):
        left, right = right, left
    return TopologyNode(pdg, mask, left, right, kind, parameter)


def _canonical_topology(nodes: Iterable[TopologyNode]) -> Topology:
    by_density: dict[tuple, TopologyNode] = {}
    for node in nodes:
        signature = node.density_signature
        representative = by_density.get(signature)
        if representative is None or node < representative:
            by_density[signature] = node
    return tuple(sorted(by_density.values(), key=lambda node: (
        node.mask.bit_count(), node.mask, node.left, node.right,
        node.kind, node.parameter, node.pdg,
    )))


def _mask_labels(mask: int) -> tuple[int, ...]:
    return tuple(label for label in range(mask.bit_length())
                 if mask & (1 << label))


def _three_partitions(mask: int) -> Iterable[tuple[int, int, int]]:
    """Yield each unordered partition of ``mask`` into three nonempty masks."""

    anchor = mask & -mask
    first = (mask - 1) & mask
    while first:
        remainder = mask ^ first
        if remainder and first & anchor:
            second_anchor = remainder & -remainder
            second = (remainder - 1) & remainder
            while second:
                third = remainder ^ second
                if third and second & second_anchor:
                    yield first, second, third
                second = (second - 1) & remainder
        first = (first - 1) & mask


def _build_combinations(catalogue: SMInteractionCatalogue, arity: int):
    combinations: dict[tuple[int, ...], set[int]] = {}
    vertices = catalogue.cubic if arity == 2 else catalogue.quartic
    for vertex in vertices:
        for current_index in range(arity + 1):
            parent = anti_pdg(vertex[current_index])
            children = tuple(sorted(
                vertex[index] for index in range(arity + 1)
                if index != current_index
            ))
            combinations.setdefault(children, set()).add(parent)
    return {key: tuple(sorted(value)) for key, value in combinations.items()}


@dataclass(frozen=True, order=True)
class _State:
    pdg: int
    topology: Topology
    spine: tuple[int, ...]
    vertices: tuple[TreeVertex, ...]


def _contact_nodes(
    parent: int,
    mask: int,
    child_masks: tuple[int, int, int],
    flavour_scheme: int,
) -> tuple[TopologyNode, TopologyNode]:
    """Return a physical parent split plus one deterministic flat scaffold."""

    ordered = tuple(sorted(child_masks, key=lambda item: (
        item.bit_count(), item
    )))
    first, second, third = ordered
    combined = second | third
    kind, parameter = transform_for_pdg(parent, flavour_scheme)
    physical = _canonical_node(
        parent, mask, first, combined, kind, parameter
    )
    scaffold = _canonical_node(
        0, combined, second, third, FLAT_CONTACT, 0
    )
    return physical, scaffold


def discover_diagram_channels(
    external_pdgs: tuple[int, ...],
    current_combinations: Mapping[tuple[int, int], Iterable[int]] | None = None,
    legacy_anti_pdg: Callable[[int], int] | None = None,
    *,
    catalogue: SMInteractionCatalogue | None = None,
    initial_labels: tuple[int, int] = (0, 1),
    mapped_final_mask: int | None = None,
    preserve_production_order: bool = True,
) -> tuple[DiagramChannel, ...]:
    """Return unique numerical maps obtained from complete physical trees.

    ``current_combinations`` and ``legacy_anti_pdg`` remain accepted for source
    compatibility, but quartic/validated discovery requires ``catalogue``.
    When only the legacy pair table is supplied, it is converted to a cubic
    catalogue containing physical fields only.
    """

    del legacy_anti_pdg
    nexternal = len(external_pdgs)
    if nexternal < 3:
        return ()
    if any(pdg not in PHYSICAL_PDGS for pdg in external_pdgs):
        raise ValueError("external states must all be physical particles")
    if catalogue is None:
        # Compatibility mode for callers of the former cubic API. Reconstruct
        # only vertices that contain physical fields; auxiliary pseudo-trees
        # can therefore never be introduced through this path.
        cubic: set[tuple[int, int, int]] = set()
        if current_combinations is not None:
            for children, parents in current_combinations.items():
                for parent in parents:
                    if parent not in PHYSICAL_PDGS:
                        continue
                    incident = tuple(sorted((*children, anti_pdg(parent))))
                    if all(pdg in PHYSICAL_PDGS for pdg in incident):
                        cubic.add(incident)
        catalogue = SMInteractionCatalogue(
            frozenset(cubic), frozenset(), 5
        )

    full_mask = (1 << nexternal) - 1
    first_initial, second_initial = initial_labels
    if first_initial == second_initial or any(
            label < 0 or label >= nexternal for label in initial_labels):
        raise ValueError("initial labels must be distinct external legs")
    initial_mask = sum(1 << label for label in initial_labels)
    all_final_mask = full_mask ^ initial_mask
    if mapped_final_mask is None:
        mapped_final_mask = all_final_mask
    if mapped_final_mask & ~all_final_mask:
        raise ValueError("mapped final-state mask contains an initial label")

    cubic_combinations = _build_combinations(catalogue, 2)
    quartic_combinations = _build_combinations(catalogue, 3)
    first_initial_bit = 1 << first_initial

    @lru_cache(maxsize=None)
    def currents(mask: int) -> tuple[_State, ...]:
        if mask & (mask - 1) == 0:
            label = mask.bit_length() - 1
            spine = (label,) if label == first_initial else ()
            return (_State(external_pdgs[label], (), spine, ()),)

        found: dict[tuple, _State] = {}

        def add_state(
            child_states: tuple[_State, ...],
            child_masks: tuple[int, ...],
            parent: int,
        ) -> None:
            incident = tuple(sorted((
                anti_pdg(parent), *(state.pdg for state in child_states)
            )))
            if not catalogue.contains(incident):
                return
            pure_qcd = is_pure_qcd_vertex(incident)
            nodes = tuple(
                node for state in child_states for node in state.topology
            )
            if mask & initial_mask == 0 and mask & mapped_final_mask == mask:
                if not pure_qcd:
                    if len(child_masks) == 2:
                        kind, parameter = transform_for_pdg(
                            parent, catalogue.flavour_scheme
                        )
                        nodes += (_canonical_node(
                            parent, mask, child_masks[0], child_masks[1],
                            kind, parameter,
                        ),)
                    else:
                        nodes += _contact_nodes(
                            parent, mask, child_masks,
                            catalogue.flavour_scheme,
                        )
            topology = _canonical_topology(nodes)

            if not preserve_production_order:
                spine: tuple[int, ...] = ()
            elif mask & first_initial_bit:
                initial_index = next(
                    index for index, child_mask in enumerate(child_masks)
                    if child_mask & first_initial_bit
                )
                spine = child_states[initial_index].spine
                branches = [
                    child_mask for index, child_mask in enumerate(child_masks)
                    if index != initial_index
                ]
                branches.sort(key=lambda item: (item.bit_count(), item))
                for branch in branches:
                    spine += _mask_labels(branch)
            else:
                spine = ()

            vertex = TreeVertex(
                parent,
                mask,
                tuple(state.pdg for state in child_states),
                child_masks,
            )
            vertices = tuple(sorted((
                *(vertex for state in child_states for vertex in state.vertices),
                vertex,
            ), key=lambda item: (
                item.mask.bit_count(), item.mask, item.parent,
                item.child_masks, item.children,
            )))
            state = _State(parent, topology, spine, vertices)
            # Physical trees which project to the same current, density and
            # production recipe are intentionally represented by one state.
            key = (parent, tuple(
                node.density_signature for node in topology
            ), spine)
            representative = found.get(key)
            if representative is None or state.vertices < representative.vertices:
                found[key] = state

        anchor = mask & -mask
        subset = (mask - 1) & mask
        while subset:
            other = mask ^ subset
            if other and subset & anchor:
                for left_state, right_state in product(
                    currents(subset), currents(other)
                ):
                    children = tuple(sorted((
                        left_state.pdg, right_state.pdg
                    )))
                    for parent in cubic_combinations.get(children, ()):
                        add_state(
                            (left_state, right_state),
                            (subset, other),
                            parent,
                        )
            subset = (subset - 1) & mask

        for first, second, third in _three_partitions(mask):
            for states in product(
                currents(first), currents(second), currents(third)
            ):
                children = tuple(sorted(state.pdg for state in states))
                for parent in quartic_combinations.get(children, ()):
                    add_state(states, (first, second, third), parent)

        return tuple(sorted(found.values()))

    rooted_mask = full_mask ^ (1 << second_initial)
    required_pdg = anti_pdg(external_pdgs[second_initial])
    canonical_order = (
        (first_initial,)
        + tuple(label for label in range(nexternal)
                if label not in initial_labels)
        + (second_initial,)
    )
    raw_channels = []
    for state in currents(rooted_mask):
        if state.pdg != required_pdg or not state.topology:
            continue
        order = state.spine + (second_initial,) \
            if preserve_production_order else canonical_order
        if len(order) != nexternal or set(order) != set(range(nexternal)):
            raise ValueError("diagram production path does not cover every leg")
        channel = DiagramChannel(
            state.topology, order, state.vertices, 1
        )
        if not validate_topology(channel.topology, all_final_mask):
            # A QCD subcurrent below an explicit EW decay boundary can be an
            # opaque ordinary-gen23 block. Such a tree is not yet a valid
            # exact binary recipe and is conservatively omitted.
            continue
        if not validate_tree_vertices(channel.vertices, catalogue):
            raise ValueError("completed tree failed physical vertex validation")
        raw_channels.append(channel)

    return canonicalize_channels_by_density(raw_channels)


def canonicalize_channels_by_density(
    channels: Iterable[DiagramChannel],
) -> tuple[DiagramChannel, ...]:
    """Collapse only channels with exactly identical numerical densities."""

    grouped: dict[tuple, list[DiagramChannel]] = {}
    for channel in channels:
        grouped.setdefault(channel.density_signature, []).append(channel)
    result = []
    for equivalent in grouped.values():
        representative = min(equivalent)
        result.append(replace(
            representative,
            multiplicity=sum(item.multiplicity for item in equivalent),
        ))
    return tuple(sorted(result))


def discover_diagram_topologies(
    external_pdgs: tuple[int, ...],
    current_combinations: Mapping[tuple[int, int], Iterable[int]] | None = None,
    legacy_anti_pdg: Callable[[int], int] | None = None,
    **kwargs,
) -> tuple[Topology, ...]:
    """Compatibility view returning unique projected topology trees."""

    return tuple(sorted({
        channel.topology
        for channel in discover_diagram_channels(
            external_pdgs,
            current_combinations,
            legacy_anti_pdg,
            **kwargs,
        )
    }))


def topology_roots(topology: Topology) -> tuple[int, ...]:
    """Return the maximal disjoint mapped masks in ``topology``."""

    masks = {node.mask for node in topology}
    roots = [
        mask for mask in masks
        if not any(mask != other and mask & other == mask for other in masks)
    ]
    return tuple(sorted(roots, key=lambda mask: (mask.bit_count(), mask)))


def validate_topology(topology: Topology, final_mask: int) -> bool:
    """Check physical fields, transforms, laminarity, and binary completion."""

    by_mask: dict[int, TopologyNode] = {}
    for node in topology:
        if node.pdg in AUXILIARY_PDGS:
            return False
        if node.pdg == 0 and node.kind != FLAT_CONTACT:
            return False
        if node.pdg != 0 and node.pdg not in PHYSICAL_PDGS:
            return False
        if node.kind not in TRANSFORM_CODES:
            return False
        if node.mask & final_mask != node.mask or node.mask in by_mask:
            return False
        by_mask[node.mask] = node
    masks = tuple(by_mask)
    for index, first in enumerate(masks):
        for second in masks[:index]:
            overlap = first & second
            if overlap and overlap != first and overlap != second:
                return False
    for root in topology_roots(topology):
        pending = [root]
        while pending:
            mask = pending.pop()
            if mask & (mask - 1) == 0:
                continue
            node = by_mask.get(mask)
            if node is None:
                return False
            pending.extend((node.left, node.right))
    return True


def validate_tree_vertices(
    vertices: Iterable[TreeVertex],
    catalogue: SMInteractionCatalogue,
) -> bool:
    """Validate every recorded completed-tree vertex against the catalogue."""

    for vertex in vertices:
        if any(pdg in AUXILIARY_PDGS for pdg in vertex.particles):
            return False
        if not catalogue.contains(vertex.particles):
            return False
        if len(vertex.children) not in (2, 3):
            return False
        if len(vertex.children) != len(vertex.child_masks):
            return False
        if any(mask <= 0 for mask in vertex.child_masks):
            return False
        combined = 0
        for mask in vertex.child_masks:
            if combined & mask:
                return False
            combined |= mask
        if combined != vertex.mask:
            return False
    return True


def has_physical_tree(
    external_pdgs: tuple[int, ...],
    catalogue: SMInteractionCatalogue,
) -> bool:
    """Return whether at least one complete cubic/quartic physical tree exists."""

    if any(pdg not in PHYSICAL_PDGS for pdg in external_pdgs):
        return False
    cubic = _build_combinations(catalogue, 2)
    quartic = _build_combinations(catalogue, 3)

    # Almost every collider subprocess in the implemented SM has a cubic
    # representation.  Prove that cheaper case first so quartic partitions
    # are only enumerated for genuinely contact-only candidates.
    @lru_cache(maxsize=None)
    def cubic_currents(mask: int) -> frozenset[int]:
        if mask & (mask - 1) == 0:
            return frozenset((external_pdgs[mask.bit_length() - 1],))
        found: set[int] = set()
        anchor = mask & -mask
        subset = (mask - 1) & mask
        while subset:
            other = mask ^ subset
            if other and subset & anchor:
                for first, second in product(
                    cubic_currents(subset), cubic_currents(other)
                ):
                    found.update(cubic.get(tuple(sorted((first, second))), ()))
            subset = (subset - 1) & mask
        return frozenset(found)

    full_mask = (1 << len(external_pdgs)) - 1
    anchor = full_mask & -full_mask
    subset = (full_mask - 1) & full_mask
    while subset:
        other = full_mask ^ subset
        if other and subset & anchor:
            if any(anti_pdg(current) in cubic_currents(other)
                   for current in cubic_currents(subset)):
                return True
        subset = (subset - 1) & full_mask

    @lru_cache(maxsize=None)
    def currents(mask: int) -> frozenset[int]:
        if mask & (mask - 1) == 0:
            return frozenset((external_pdgs[mask.bit_length() - 1],))
        found: set[int] = set()
        anchor = mask & -mask
        subset = (mask - 1) & mask
        while subset:
            other = mask ^ subset
            if other and subset & anchor:
                for first, second in product(currents(subset), currents(other)):
                    found.update(cubic.get(tuple(sorted((first, second))), ()))
            subset = (subset - 1) & mask
        for first_mask, second_mask, third_mask in _three_partitions(mask):
            for states in product(
                currents(first_mask), currents(second_mask), currents(third_mask)
            ):
                found.update(quartic.get(tuple(sorted(states)), ()))
        return frozenset(found)

    anchor = full_mask & -full_mask
    subset = (full_mask - 1) & full_mask
    while subset:
        other = full_mask ^ subset
        if other and subset & anchor:
            if any(anti_pdg(current) in currents(other)
                   for current in currents(subset)):
                return True
        subset = (subset - 1) & full_mask

    # A completed four-point root can connect four disjoint external currents.
    anchor = full_mask & -full_mask
    first = (full_mask - 1) & full_mask
    while first:
        rest = full_mask ^ first
        if rest and first & anchor:
            anchor2 = rest & -rest
            second = (rest - 1) & rest
            while second:
                rest2 = rest ^ second
                if rest2 and second & anchor2:
                    anchor3 = rest2 & -rest2
                    third = (rest2 - 1) & rest2
                    while third:
                        fourth = rest2 ^ third
                        if fourth and third & anchor3:
                            for states in product(
                                currents(first), currents(second),
                                currents(third), currents(fourth),
                            ):
                                if tuple(sorted(states)) in catalogue.quartic:
                                    return True
                        third = (third - 1) & rest2
                second = (second - 1) & rest
        first = (first - 1) & full_mask
    return False
