from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from itertools import permutations, product
from typing import Iterable, Literal, Sequence

from .process_ir import CanonicalProcessIR, ProcessLegIR, build_process_ir
from .processes import ProcessOptions

ColorAccuracy = Literal["lc", "nlc", "full"]
ColorSectorKind = Literal["singlet", "open-lines", "single-trace"]


@dataclass(frozen=True)
class LCQuarkLine:
    """One leading-colour open colour line in all-outgoing conventions."""

    quark_label: int
    antiquark_label: int
    gluon_labels: tuple[int, ...]
    singlet_labels: tuple[int, ...] = ()

    @property
    def coloured_labels(self) -> tuple[int, ...]:
        return (self.quark_label, *self.gluon_labels, self.antiquark_label)

    @property
    def line_labels(self) -> tuple[int, ...]:
        return (*self.coloured_labels, *self.singlet_labels)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "quark_label": self.quark_label,
            "antiquark_label": self.antiquark_label,
            "gluon_labels": list(self.gluon_labels),
            "singlet_labels": list(self.singlet_labels),
            "line_labels": list(self.line_labels),
        }


@dataclass(frozen=True)
class LCColorSector:
    """One colour-flow sector produced during generation warmup."""

    id: int
    kind: ColorSectorKind
    quark_lines: tuple[LCQuarkLine, ...] = ()
    trace_labels: tuple[int, ...] = ()
    singlet_labels: tuple[int, ...] = ()
    word_labels: tuple[int, ...] = ()

    @property
    def coloured_label_groups(self) -> tuple[tuple[int, ...], ...]:
        if self.kind == "open-lines":
            return tuple(line.coloured_labels for line in self.quark_lines)
        if self.kind == "single-trace":
            return (self.trace_labels,)
        return ()

    @property
    def line_label_groups(self) -> tuple[tuple[int, ...], ...]:
        if self.kind == "open-lines":
            return tuple(line.coloured_labels for line in self.quark_lines)
        if self.kind == "single-trace":
            return ((*self.trace_labels, *self.singlet_labels),)
        if self.kind == "singlet":
            return (self.singlet_labels,)
        return ()

    @cached_property
    def color_words(self) -> tuple[tuple[int, ...], ...]:
        if self.word_labels:
            return (self.word_labels,)
        if self.kind == "open-lines":
            return (
                tuple(label for line in self.quark_lines for label in line.coloured_labels),
            )
        if self.kind == "single-trace":
            return (self.trace_labels,)
        if self.kind == "singlet":
            return ((),)
        return ()

    @cached_property
    def compatibility_words(self) -> tuple[tuple[int, ...], ...]:
        """Colour words accepted while constructing currents for this sector.

        The sector itself has one physical colour word.  During current
        construction, however, complete open quark-line blocks can be traversed
        in different intermediate orders.  This lets the generic recursion
        reproduce AmpliCol's colour-ordered current closures without naming a
        process family.
        """

        if self.kind != "open-lines" or len(self.quark_lines) < 2:
            return self.color_words
        primary = self.color_words[0]
        blocks = _ordered_open_line_blocks(primary, self.quark_lines)
        if blocks is None or len(blocks) < 2:
            return self.color_words
        words: list[tuple[int, ...]] = [primary]
        seen = {primary}
        for block_permutation in permutations(blocks):
            word = tuple(
                label
                for block in block_permutation
                for label in block
            )
            if word in seen:
                continue
            seen.add(word)
            words.append(word)
        return tuple(words)

    @cached_property
    def legacy_order_words(self) -> tuple[tuple[int, ...], ...]:
        """Full phase-space orders, including colour-singlet attachments."""

        if self.kind == "open-lines":
            return _open_line_legacy_order_words(self)
        return self.color_words

    def to_json_dict(self) -> dict[str, object]:
        return {
            "id": self.id,
            "kind": self.kind,
            "quark_lines": [line.to_json_dict() for line in self.quark_lines],
            "trace_labels": list(self.trace_labels),
            "singlet_labels": list(self.singlet_labels),
            "word_labels": list(self.word_labels),
            "coloured_label_groups": [
                list(group) for group in self.coloured_label_groups
            ],
            "line_label_groups": [
                list(group) for group in self.line_label_groups
            ],
            "color_words": [list(word) for word in self.color_words],
            "compatibility_words": [
                list(word) for word in self.compatibility_words
            ],
            "legacy_order_words": [
                list(word) for word in self.legacy_order_words
            ],
        }


@dataclass(frozen=True)
class LCColorSectorTopologyGroup:
    """Colour sectors that can share one compiled current topology.

    The group key is built from model/process data carried by the external
    labels, not from a process-family name.  Runtime reuse still evaluates each
    sector with its own external-label permutation, so this is a generation-time
    sharing plan rather than a physics approximation.
    """

    signature: tuple[object, ...]
    representative_sector_id: int
    sector_ids: tuple[int, ...]
    label_permutations: tuple[tuple[tuple[int, int], ...], ...]

    def to_json_dict(self) -> dict[str, object]:
        return {
            "signature": _jsonable_signature(self.signature),
            "representative_sector_id": self.representative_sector_id,
            "sector_ids": list(self.sector_ids),
            "label_permutations": [
                [[left, right] for left, right in permutation]
                for permutation in self.label_permutations
            ],
        }


@dataclass(frozen=True)
class LCColorSectorReplayPartition:
    """Initial-label-safe replay block inside one LC topology group."""

    representative_sector_id: int
    active_sector_ids: tuple[int, ...]
    label_permutations: tuple[tuple[tuple[int, int], ...], ...]
    replay_weights: tuple[float, ...] = ()

    def to_json_dict(self) -> dict[str, object]:
        return {
            "representative_sector_id": self.representative_sector_id,
            "active_sector_ids": list(self.active_sector_ids),
            "label_permutations": [
                [[left, right] for left, right in permutation]
                for permutation in self.label_permutations
            ],
            "replay_weights": list(self.replay_weights),
        }


@dataclass(frozen=True)
class GenericColorPlan:
    """Colour-flow planning payload shared by future Python/Rust runtimes."""

    process: CanonicalProcessIR
    color_accuracy: str
    sectors: tuple[LCColorSector, ...]
    diagnostics: tuple[str, ...] = ()
    truncated: bool = False
    idenso_required: bool = False

    @property
    def sector_count(self) -> int:
        return len(self.sectors)

    @property
    def ready_for_leading_colour(self) -> bool:
        return self.color_accuracy == "lc" and bool(self.sectors) and not self.truncated

    @property
    def ready_for_requested_colour(self) -> bool:
        return bool(self.sectors) and not self.truncated and not self.idenso_required

    @property
    def coloured_labels(self) -> tuple[int, ...]:
        return tuple(
            sorted(
                {
                    *self.process.quark_labels,
                    *self.process.antiquark_labels,
                    *self.process.gluon_labels,
                }
            )
        )

    @cached_property
    def _sectors_by_id(self) -> dict[int, LCColorSector]:
        return {int(sector.id): sector for sector in self.sectors}

    def sector(self, color_sector: int) -> LCColorSector | None:
        return self._sectors_by_id.get(int(color_sector))

    @cached_property
    def topology_groups(self) -> tuple[LCColorSectorTopologyGroup, ...]:
        return _sector_topology_groups(self.process, self.sectors)

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process.to_json_dict(),
            "color_accuracy": self.color_accuracy,
            "sector_count": self.sector_count,
            "truncated": self.truncated,
            "idenso_required": self.idenso_required,
            "ready_for_leading_colour": self.ready_for_leading_colour,
            "ready_for_requested_colour": self.ready_for_requested_colour,
            "coloured_labels": list(self.coloured_labels),
            "diagnostics": list(self.diagnostics),
            "sectors": [sector.to_json_dict() for sector in self.sectors],
            "topology_groups": [
                group.to_json_dict() for group in self.topology_groups
            ],
        }


def build_color_plan(
    process: str | CanonicalProcessIR,
    *,
    color_accuracy: str = "lc",
    options: ProcessOptions | None = None,
    max_sectors: int | None = None,
    reference_color_order: Sequence[int] | None = None,
) -> GenericColorPlan:
    process_ir = (
        process
        if isinstance(process, CanonicalProcessIR)
        else build_process_ir(process, color_accuracy=color_accuracy, options=options)
    )
    max_sector_count = _normalize_sector_cap(max_sectors)
    if color_accuracy != process_ir.color_accuracy:
        process_ir = build_process_ir(
            process_ir.process,
            color_accuracy=color_accuracy,
            options=options,
        )
    quark_legs = _legs_by_labels(process_ir, process_ir.quark_labels)
    antiquark_legs = _legs_by_labels(process_ir, process_ir.antiquark_labels)
    gluon_labels = process_ir.gluon_labels
    singlet_labels = process_ir.singlet_labels
    if len(quark_legs) != len(antiquark_legs):
        return GenericColorPlan(
            process=process_ir,
            color_accuracy=color_accuracy,
            sectors=(),
            diagnostics=(
                "leading-colour open-line plan requires balanced outgoing "
                "quark and antiquark counts",
            ),
        )
    if not quark_legs:
        return _build_no_quark_color_plan(
            process_ir,
            gluon_labels=gluon_labels,
            singlet_labels=singlet_labels,
            max_sectors=max_sector_count,
            reference_color_order=reference_color_order,
        )

    sectors: list[LCColorSector] = list(
        _reference_lc_color_sectors(process_ir, reference_color_order)
    )
    seen_sector_keys = {_sector_dedup_key(sector) for sector in sectors}
    truncated = False
    for antiquark_permutation in permutations(antiquark_legs):
        for gluon_allocation in _iter_ordered_gluon_allocations(
            gluon_labels,
            len(quark_legs),
        ):
            lines = tuple(
                LCQuarkLine(
                    quark_label=quark.label,
                    antiquark_label=antiquark.label,
                    gluon_labels=tuple(gluon_allocation[index]),
                )
                for index, (quark, antiquark) in enumerate(
                    zip(quark_legs, antiquark_permutation, strict=True)
                )
            )
            for word_labels in _iter_open_line_color_words(
                lines,
                include_block_permutations=color_accuracy != "lc",
            ):
                candidate = LCColorSector(
                    id=len(sectors),
                    kind="open-lines",
                    quark_lines=lines,
                    singlet_labels=singlet_labels,
                    word_labels=word_labels,
                )
                key = _sector_dedup_key(candidate)
                if key in seen_sector_keys:
                    continue
                seen_sector_keys.add(key)
                sectors.append(candidate)
                if max_sector_count is not None and len(sectors) >= max_sector_count:
                    truncated = True
                    break
            if truncated:
                break
        if truncated:
            break

    diagnostics: tuple[str, ...] = ()
    if truncated:
        diagnostics = (
            f"leading-colour sector enumeration reached max_sectors={max_sector_count}",
        )
    return GenericColorPlan(
        process=process_ir,
        color_accuracy=color_accuracy,
        sectors=tuple(sectors),
        diagnostics=diagnostics,
        truncated=truncated,
    )


def _normalize_sector_cap(value: int | None) -> int | None:
    if value is None:
        return None
    normalized = int(value)
    return None if normalized < 0 else normalized


def _build_no_quark_color_plan(
    process: CanonicalProcessIR,
    *,
    gluon_labels: tuple[int, ...],
    singlet_labels: tuple[int, ...],
    max_sectors: int | None,
    reference_color_order: Sequence[int] | None = None,
) -> GenericColorPlan:
    if not gluon_labels:
        return GenericColorPlan(
            process=process,
            color_accuracy=process.color_accuracy,
            sectors=(
                LCColorSector(
                    id=0,
                    kind="singlet",
                    singlet_labels=singlet_labels,
                ),
            ),
        )

    first = min(gluon_labels)
    rest = tuple(label for label in gluon_labels if label != first)
    sectors: list[LCColorSector] = list(
        _reference_lc_color_sectors(process, reference_color_order)
    )
    seen_sector_keys = {_sector_dedup_key(sector) for sector in sectors}
    truncated = False
    fold_reflections = process.color_accuracy == "lc"
    seen_reversal_classes: set[tuple[int, ...]] = set()
    for ordered_rest in permutations(rest):
        if fold_reflections:
            canonical = min(ordered_rest, tuple(reversed(ordered_rest)))
            if canonical in seen_reversal_classes:
                continue
            seen_reversal_classes.add(canonical)
        candidate = LCColorSector(
            id=len(sectors),
            kind="single-trace",
            trace_labels=(first, *ordered_rest),
            singlet_labels=singlet_labels,
        )
        key = _sector_dedup_key(candidate)
        if key in seen_sector_keys:
            continue
        seen_sector_keys.add(key)
        sectors.append(candidate)
        if max_sectors is not None and len(sectors) >= max_sectors:
            truncated = True
            break
    diagnostics: tuple[str, ...] = ()
    if truncated:
        diagnostics = (
            f"leading-colour trace enumeration reached max_sectors={max_sectors}",
        )
    return GenericColorPlan(
        process=process,
        color_accuracy=process.color_accuracy,
        sectors=tuple(sectors),
        diagnostics=diagnostics,
        truncated=truncated,
    )


def _legs_by_labels(
    process: CanonicalProcessIR,
    labels: tuple[int, ...],
) -> tuple[ProcessLegIR, ...]:
    by_label = {leg.label: leg for leg in process.legs}
    return tuple(by_label[label] for label in labels)


def _reference_lc_color_sectors(
    process: CanonicalProcessIR,
    reference_color_order: Sequence[int] | None,
) -> tuple[LCColorSector, ...]:
    if reference_color_order is None:
        return ()
    reference = tuple(int(label) for label in reference_color_order)
    if not reference:
        return ()
    by_label = {leg.label: leg for leg in process.legs}
    coloured_labels = {
        *process.quark_labels,
        *process.antiquark_labels,
        *process.gluon_labels,
    }
    coloured_word = tuple(label for label in reference if label in coloured_labels)
    if sorted(coloured_word) != sorted(coloured_labels):
        return ()
    if not process.quark_labels:
        if tuple(label for label in coloured_word if label in process.gluon_labels) != coloured_word:
            return ()
        return (
            LCColorSector(
                id=0,
                kind="single-trace",
                trace_labels=coloured_word,
                singlet_labels=process.singlet_labels,
            ),
        )
    lines = _reference_open_lines(process, coloured_word, by_label)
    if lines is None:
        return ()
    word_labels = tuple(label for line in lines for label in line.coloured_labels)
    return (
        LCColorSector(
            id=0,
            kind="open-lines",
            quark_lines=lines,
            singlet_labels=process.singlet_labels,
            word_labels=word_labels,
        ),
    )


def _reference_open_lines(
    process: CanonicalProcessIR,
    coloured_word: tuple[int, ...],
    by_label: dict[int, ProcessLegIR],
) -> tuple[LCQuarkLine, ...] | None:
    quark_labels = set(process.quark_labels)
    antiquark_labels = set(process.antiquark_labels)
    gluon_labels = set(process.gluon_labels)
    lines: list[LCQuarkLine] = []
    offset = 0
    while offset < len(coloured_word):
        start = coloured_word[offset]
        if start in quark_labels:
            end_labels = antiquark_labels
            start_is_quark = True
        elif start in antiquark_labels:
            end_labels = quark_labels
            start_is_quark = False
        else:
            return None
        end = offset + 1
        while end < len(coloured_word) and coloured_word[end] in gluon_labels:
            end += 1
        if end >= len(coloured_word) or coloured_word[end] not in end_labels:
            return None
        middle = coloured_word[offset + 1 : end]
        if any(label not in gluon_labels for label in middle):
            return None
        stop = coloured_word[end]
        if start_is_quark:
            quark_label = start
            antiquark_label = stop
        else:
            quark_label = stop
            antiquark_label = start
        if by_label[quark_label].particle_class != "quark":
            return None
        if by_label[antiquark_label].particle_class != "antiquark":
            return None
        lines.append(
            LCQuarkLine(
                quark_label=quark_label,
                antiquark_label=antiquark_label,
                gluon_labels=middle,
            )
        )
        offset = end + 1
    if len(lines) != len(process.quark_labels):
        return None
    return tuple(lines)


def _sector_dedup_key(sector: LCColorSector) -> tuple[object, ...]:
    return (
        sector.kind,
        tuple(
            (
                line.quark_label,
                line.antiquark_label,
                line.gluon_labels,
                line.singlet_labels,
            )
            for line in sector.quark_lines
        ),
        sector.trace_labels,
        sector.singlet_labels,
        sector.word_labels,
    )


def _iter_open_line_color_words(
    lines: tuple[LCQuarkLine, ...],
    *,
    include_block_permutations: bool = False,
) -> tuple[tuple[int, ...], ...]:
    """Return explicit colour words for one open-line pairing/allocation."""

    if not include_block_permutations or len(lines) < 2:
        return (tuple(label for line in lines for label in line.coloured_labels),)
    words: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for line_permutation in permutations(lines):
        word = tuple(
            label
            for line in line_permutation
            for label in line.coloured_labels
        )
        if word in seen:
            continue
        seen.add(word)
        words.append(word)
    return tuple(words)


def _ordered_open_line_blocks(
    word: tuple[int, ...],
    lines: tuple[LCQuarkLine, ...],
) -> tuple[tuple[int, ...], ...] | None:
    remaining = list(lines)
    blocks: list[tuple[int, ...]] = []
    offset = 0
    while offset < len(word):
        matched_index = None
        for index, line in enumerate(remaining):
            block = line.coloured_labels
            if word[offset : offset + len(block)] == block:
                matched_index = index
                blocks.append(block)
                offset += len(block)
                break
        if matched_index is None:
            return None
        remaining.pop(matched_index)
    if remaining:
        return None
    return tuple(blocks)


def _open_line_legacy_order_words(
    sector: LCColorSector,
) -> tuple[tuple[int, ...], ...]:
    """Full open-line block orders accepted in legacy process rows."""

    if sector.word_labels:
        line_by_coloured = {line.coloured_labels: line for line in sector.quark_lines}
        ordered_blocks = _ordered_open_line_blocks(
            sector.word_labels,
            sector.quark_lines,
        )
        if ordered_blocks is None:
            return (sector.word_labels,)
        ordered_lines = tuple(line_by_coloured[block] for block in ordered_blocks)
    else:
        ordered_lines = sector.quark_lines

    line_orientations = tuple(
        _legacy_line_orientations(line) for line in ordered_lines
    )
    words: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for line_permutation in permutations(range(len(ordered_lines))):
        orientation_choices = (
            line_orientations[index] for index in line_permutation
        )
        for blocks in product(*orientation_choices):
            word = tuple(label for block in blocks for label in block)
            if word in seen:
                continue
            seen.add(word)
            words.append(word)
    return tuple(words)


def _legacy_line_orientations(line: LCQuarkLine) -> tuple[tuple[int, ...], ...]:
    canonical = line.line_labels
    antiquark_first = (
        line.antiquark_label,
        *line.gluon_labels,
        line.quark_label,
        *line.singlet_labels,
    )
    if antiquark_first == canonical:
        return (canonical,)
    return (canonical, antiquark_first)


def _iter_ordered_gluon_allocations(
    gluon_labels: tuple[int, ...],
    line_count: int,
) -> Iterable[tuple[tuple[int, ...], ...]]:
    if line_count <= 0:
        yield ()
        return
    for ordered_gluons in permutations(gluon_labels):
        yield from _iter_split_ordered_sequence(tuple(ordered_gluons), line_count)


def _iter_ordered_label_allocations(
    labels: tuple[int, ...],
    line_count: int,
) -> Iterable[tuple[tuple[int, ...], ...]]:
    if line_count <= 0:
        yield ()
        return
    if not labels:
        yield tuple(() for _ in range(line_count))
        return
    buckets: list[list[int]] = [[] for _ in range(line_count)]
    yield from _assign_unordered_labels_to_buckets(labels, buckets, 0)


def _assign_unordered_labels_to_buckets(
    labels: tuple[int, ...],
    buckets: list[list[int]],
    index: int,
) -> Iterable[tuple[tuple[int, ...], ...]]:
    if index == len(labels):
        yield tuple(tuple(bucket) for bucket in buckets)
        return
    label = labels[index]
    for bucket in buckets:
        bucket.append(label)
        yield from _assign_unordered_labels_to_buckets(labels, buckets, index + 1)
        bucket.pop()


def _iter_split_ordered_sequence(
    sequence: tuple[int, ...],
    bin_count: int,
) -> Iterable[tuple[tuple[int, ...], ...]]:
    if bin_count == 1:
        yield (sequence,)
        return
    for split_index in range(len(sequence) + 1):
        head = sequence[:split_index]
        for tail in _iter_split_ordered_sequence(
            sequence[split_index:],
            bin_count - 1,
        ):
            yield (head, *tail)


def _sector_topology_groups(
    process: CanonicalProcessIR,
    sectors: tuple[LCColorSector, ...],
) -> tuple[LCColorSectorTopologyGroup, ...]:
    by_signature: dict[tuple[object, ...], list[LCColorSector]] = {}
    for sector in sectors:
        by_signature.setdefault(
            _sector_topology_signature(process, sector),
            [],
        ).append(sector)

    groups: list[LCColorSectorTopologyGroup] = []
    for signature, sector_group in by_signature.items():
        representative = sector_group[0]
        representative_labels = _sector_topology_labels(representative)
        permutations: list[tuple[tuple[int, int], ...]] = []
        for sector in sector_group:
            sector_labels = _sector_topology_labels(sector)
            if len(sector_labels) != len(representative_labels):
                raise ValueError("isomorphic colour sectors have mismatched labels")
            permutations.append(tuple(zip(representative_labels, sector_labels)))
        groups.append(
            LCColorSectorTopologyGroup(
                signature=signature,
                representative_sector_id=representative.id,
                sector_ids=tuple(sector.id for sector in sector_group),
                label_permutations=tuple(permutations),
            )
        )
    return tuple(groups)


def _sector_topology_signature(
    process: CanonicalProcessIR,
    sector: LCColorSector,
) -> tuple[object, ...]:
    pdg_by_label = {
        leg.label: leg.outgoing_pdg
        for leg in process.legs
        if leg.outgoing_pdg is not None
    }
    if sector.kind == "open-lines":
        line_by_coloured = {line.coloured_labels: line for line in sector.quark_lines}
        ordered_blocks = (
            _ordered_open_line_blocks(sector.word_labels, sector.quark_lines)
            if sector.word_labels
            else None
        )
        ordered_lines = (
            tuple(line_by_coloured[block] for block in ordered_blocks)
            if ordered_blocks is not None
            else sector.quark_lines
        )
        return (
            sector.kind,
            tuple(
                (
                    pdg_by_label[line.quark_label],
                    tuple(pdg_by_label[label] for label in line.gluon_labels),
                    pdg_by_label[line.antiquark_label],
                    tuple(pdg_by_label[label] for label in line.singlet_labels),
                )
                for line in ordered_lines
            ),
        )
    if sector.kind == "single-trace":
        return (
            sector.kind,
            tuple(pdg_by_label[label] for label in sector.trace_labels),
            tuple(pdg_by_label[label] for label in sector.singlet_labels),
        )
    return (
        sector.kind,
        tuple(pdg_by_label[label] for label in sector.singlet_labels),
    )


def _sector_topology_labels(sector: LCColorSector) -> tuple[int, ...]:
    if sector.kind == "open-lines":
        line_by_coloured = {line.coloured_labels: line for line in sector.quark_lines}
        ordered_blocks = (
            _ordered_open_line_blocks(sector.word_labels, sector.quark_lines)
            if sector.word_labels
            else None
        )
        ordered_lines = (
            tuple(line_by_coloured[block] for block in ordered_blocks)
            if ordered_blocks is not None
            else sector.quark_lines
        )
        return tuple(label for line in ordered_lines for label in line.line_labels)
    if sector.kind == "single-trace":
        return (*sector.trace_labels, *sector.singlet_labels)
    return sector.singlet_labels


def _jsonable_signature(signature: tuple[object, ...]) -> list[object]:
    def convert(value: object) -> object:
        if isinstance(value, tuple):
            return [convert(item) for item in value]
        return value

    return [convert(item) for item in signature]


def lc_topology_replay_safe_groups(
    color_plan: GenericColorPlan,
) -> tuple[LCColorSectorTopologyGroup, ...]:
    """Return LC topology groups safe for physical-input replay.

    Runtime replay currently evaluates one representative sector with physical
    external momenta in the user's process ordering.  Therefore every replay
    permutation must preserve the initial-state label set.  Pure single-trace
    sectors are also kept out of this fast path for now because their trace
    symmetries need dedicated colour-convention validation.
    """

    if color_plan.color_accuracy != "lc":
        return ()
    return tuple(
        group
        for group in color_plan.topology_groups
        if _lc_topology_group_replay_safe(color_plan, group)
    )


def lc_topology_replay_partitions(
    color_plan: GenericColorPlan,
) -> tuple[LCColorSectorReplayPartition, ...]:
    """Partition LC topology groups into exact replay-safe representatives.

    The selected-sector ``lc_topology_replay_safe_groups`` helper only accepts
    groups whose full topology orbit preserves the initial-state label set.
    All-flow replay artifacts can be more general: they split one topology
    group into several initial-label-safe blocks and materialize one
    representative sidecar per block.  This also covers pure single-trace
    gluon sectors without falling back to an all-sector artifact.
    """

    if color_plan.color_accuracy != "lc":
        return ()
    initial_labels = {leg.label for leg in color_plan.process.initial_legs}
    if not initial_labels:
        return ()
    partitions: list[LCColorSectorReplayPartition] = []
    for group in color_plan.topology_groups:
        sectors = tuple(
            sector
            for sector_id in group.sector_ids
            if (sector := color_plan.sector(sector_id)) is not None
        )
        if len(sectors) != len(group.sector_ids):
            continue
        if any(sector.kind not in {"open-lines", "single-trace"} for sector in sectors):
            continue
        base_maps = {
            int(sector_id): {
                int(representative_label): int(sector_label)
                for representative_label, sector_label in permutation
            }
            for sector_id, permutation in zip(
                group.sector_ids,
                group.label_permutations,
                strict=True,
            )
        }
        sector_ids_by_initial_preimage: dict[tuple[int, ...], list[int]] = {}
        complete_initial_maps = True
        for sector_id in group.sector_ids:
            sector_map = base_maps[int(sector_id)]
            inverse_sector_map = {
                sector_label: representative_label
                for representative_label, sector_label in sector_map.items()
            }
            if not initial_labels.issubset(inverse_sector_map):
                complete_initial_maps = False
                break
            initial_preimage = tuple(
                sorted(inverse_sector_map[label] for label in initial_labels)
            )
            sector_ids_by_initial_preimage.setdefault(initial_preimage, []).append(
                int(sector_id)
            )
        if not complete_initial_maps:
            continue
        partition_sector_ids = sorted(
            sector_ids_by_initial_preimage.values(),
            key=lambda sector_ids: min(sector_ids),
        )
        for grouped_sector_ids in partition_sector_ids:
            representative = min(grouped_sector_ids)
            representative_map = base_maps[representative]
            inverse_representative_map = {
                sector_label: representative_label
                for representative_label, sector_label in representative_map.items()
            }
            active_sector_ids: list[int] = []
            relative_permutations: list[tuple[tuple[int, int], ...]] = []
            replay_weights: list[float] = []
            grouped_sector_id_set = set(grouped_sector_ids)
            for sector_id in group.sector_ids:
                if int(sector_id) not in grouped_sector_id_set:
                    continue
                sector = color_plan.sector(int(sector_id))
                if sector is None:
                    continue
                sector_map = base_maps[int(sector_id)]
                relative_map = {
                    representative_label: sector_map[
                        inverse_representative_map[representative_label]
                    ]
                    for representative_label in sorted(representative_map.values())
                }
                if {
                    relative_map[label]
                    for label in initial_labels
                } != initial_labels:
                    raise RuntimeError(
                        "internal LC replay partitioning error: sectors grouped "
                        "by initial-label preimage do not preserve the initial set"
                    )
                active_sector_ids.append(int(sector_id))
                relative_permutations.append(tuple(sorted(relative_map.items())))
                replay_weights.append(
                    _lc_topology_replay_sector_weight(color_plan, sector)
                )
            if not active_sector_ids:
                raise RuntimeError(
                    "internal LC replay partitioning error: representative sector "
                    f"{representative} did not produce a non-empty replay block"
                )
            partitions.append(
                LCColorSectorReplayPartition(
                    representative_sector_id=representative,
                    active_sector_ids=tuple(active_sector_ids),
                    label_permutations=tuple(relative_permutations),
                    replay_weights=tuple(replay_weights),
                )
            )
    return tuple(partitions)


def _lc_topology_replay_sector_weight(
    color_plan: GenericColorPlan,
    sector: LCColorSector,
) -> float:
    """Return the LC multiplicity represented by a materialized sector.

    Pure single-trace LC plans fold trace reflections during colour-sector
    enumeration.  AmpliCol's all-ordering ``imode=2`` basis keeps both trace
    orientations, while the reflected colour-ordered amplitude has the same
    squared contribution.  Replaying the folded sector with weight two preserves
    the full all-ordering sum without compiling or evaluating the reflected
    duplicate.
    """

    if (
        color_plan.color_accuracy == "lc"
        and sector.kind == "single-trace"
        and len(sector.trace_labels) > 2
    ):
        return 2.0
    return 1.0


def lc_line_pairing_representative_ids(
    color_plan: GenericColorPlan,
) -> tuple[int, ...]:
    """Return one sector per LC open-line pairing/allocation.

    This is a generic colour-flow pruning helper.  It keeps distinct
    quark-antiquark pairings and distinct gluon/singlet attachments, but drops
    duplicate sectors that only permute complete open-line blocks in the colour
    word.  It deliberately does not group sectors with different flavour
    pairings or different particles assigned to a line.
    """

    if color_plan.color_accuracy != "lc":
        return ()
    representatives: list[int] = []
    seen: set[tuple[object, ...]] = set()
    pdg_by_label = {
        leg.label: leg.outgoing_pdg
        for leg in color_plan.process.legs
        if leg.outgoing_pdg is not None
    }
    for sector in color_plan.sectors:
        signature = _sector_line_pairing_signature(pdg_by_label, sector)
        if signature in seen:
            continue
        seen.add(signature)
        representatives.append(int(sector.id))
    return tuple(representatives)


def _sector_line_pairing_signature(
    pdg_by_label: dict[int, int | None],
    sector: LCColorSector,
) -> tuple[object, ...]:
    if sector.kind == "open-lines":
        return (
            sector.kind,
            tuple(
                sorted(
                    (
                        line.quark_label,
                        pdg_by_label.get(line.quark_label),
                        line.antiquark_label,
                        pdg_by_label.get(line.antiquark_label),
                        tuple(
                            (label, pdg_by_label.get(label))
                            for label in line.gluon_labels
                        ),
                        tuple(
                            (label, pdg_by_label.get(label))
                            for label in line.singlet_labels
                        ),
                    )
                    for line in sector.quark_lines
                )
            ),
        )
    if sector.kind == "single-trace":
        canonical_trace = min(sector.trace_labels, tuple(reversed(sector.trace_labels)))
        return (
            sector.kind,
            tuple((label, pdg_by_label.get(label)) for label in canonical_trace),
            tuple((label, pdg_by_label.get(label)) for label in sector.singlet_labels),
        )
    return (
        sector.kind,
        tuple((label, pdg_by_label.get(label)) for label in sector.singlet_labels),
    )


def _lc_topology_group_replay_safe(
    color_plan: GenericColorPlan,
    group: LCColorSectorTopologyGroup,
) -> bool:
    sectors = tuple(
        sector
        for sector_id in group.sector_ids
        if (sector := color_plan.sector(sector_id)) is not None
    )
    if len(sectors) != len(group.sector_ids):
        return False
    if any(sector.kind != "open-lines" for sector in sectors):
        return False
    initial_labels = {leg.label for leg in color_plan.process.initial_legs}
    if not initial_labels:
        return False
    for permutation in group.label_permutations:
        mapping = {
            int(representative_label): int(sector_label)
            for representative_label, sector_label in permutation
        }
        if {
            mapping.get(label)
            for label in initial_labels
        } != initial_labels:
            return False
    return True


__all__ = [
    "ColorAccuracy",
    "ColorSectorKind",
    "GenericColorPlan",
    "LCColorSectorReplayPartition",
    "LCColorSectorTopologyGroup",
    "LCColorSector",
    "LCQuarkLine",
    "build_color_plan",
    "lc_line_pairing_representative_ids",
    "lc_topology_replay_partitions",
    "lc_topology_replay_safe_groups",
]
