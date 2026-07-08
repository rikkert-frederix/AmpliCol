from __future__ import annotations

import itertools
import math
import re
import time
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable, Sequence

ParticleName = str
ProcessTuple = tuple[ParticleName, ...]
OrderTuple = tuple[int, ...]


QUARKS = frozenset({"d", "u", "s", "c", "b", "t"})
ANTIQUARKS = frozenset({f"{q}~" for q in QUARKS})
SINGLETS = frozenset(
    {
        "a",
        "z",
        "w+",
        "w-",
        "e+",
        "e-",
        "mu+",
        "mu-",
        "ta+",
        "ta-",
        "ve",
        "ve~",
        "vm",
        "vm~",
        "vt",
        "vt~",
        "h",
    }
)
GLUONS = frozenset({"g"})
ALL_COLOURED = QUARKS | ANTIQUARKS | GLUONS
PDGS = {
    "g": "21",
    "d": "1",
    "u": "2",
    "s": "3",
    "c": "4",
    "b": "5",
    "t": "6",
    "d~": "-1",
    "u~": "-2",
    "s~": "-3",
    "c~": "-4",
    "b~": "-5",
    "t~": "-6",
    "a": "22",
    "z": "23",
    "w+": "24",
    "w-": "-24",
    "e+": "-11",
    "e-": "11",
    "mu+": "-13",
    "mu-": "13",
    "ta+": "-15",
    "ta-": "15",
    "ve": "12",
    "ve~": "-12",
    "vm": "14",
    "vm~": "-14",
    "vt": "16",
    "vt~": "-16",
    "h": "25",
}
ANTI_PARTICLE = {
    "g": "g",
    "d": "d~",
    "u": "u~",
    "s": "s~",
    "c": "c~",
    "b": "b~",
    "t": "t~",
    "d~": "d",
    "u~": "u",
    "s~": "s",
    "c~": "c",
    "b~": "b",
    "t~": "t",
    "a": "a",
    "z": "z",
    "w+": "w-",
    "w-": "w+",
    "e+": "e-",
    "e-": "e+",
    "mu+": "mu-",
    "mu-": "mu+",
    "ta+": "ta-",
    "ta-": "ta+",
    "ve": "ve~",
    "ve~": "ve",
    "vm": "vm~",
    "vm~": "vm",
    "vt": "vt~",
    "vt~": "vt",
    "h": "h",
}
SORT_PARTICLES = {
    "g": 13,
    "d": 1,
    "u": 2,
    "s": 3,
    "c": 4,
    "b": 5,
    "t": 6,
    "d~": 7,
    "u~": 8,
    "s~": 9,
    "c~": 10,
    "b~": 11,
    "t~": 12,
    "a": 80,
    "z": 81,
    "w+": 82,
    "w-": 83,
    "e+": 84,
    "e-": 85,
    "mu+": 86,
    "mu-": 87,
    "ta+": 88,
    "ta-": 89,
    "ve": 90,
    "ve~": 91,
    "vm": 92,
    "vm~": 93,
    "vt": 94,
    "vt~": 95,
    "h": 96,
}
CHARGES3 = {
    "g": 0,
    "d": -1,
    "u": 2,
    "s": -1,
    "c": 2,
    "b": -1,
    "t": 2,
    "d~": 1,
    "u~": -2,
    "s~": 1,
    "c~": -2,
    "b~": 1,
    "t~": -2,
    "a": 0,
    "z": 0,
    "w+": 3,
    "w-": -3,
    "e+": 3,
    "e-": -3,
    "mu+": 3,
    "mu-": -3,
    "ta+": 3,
    "ta-": -3,
    "ve": 0,
    "ve~": 0,
    "vm": 0,
    "vm~": 0,
    "vt": 0,
    "vt~": 0,
    "h": 0,
}
FAMILY = {
    "g": 0,
    "d": 1,
    "u": 1,
    "s": 11,
    "c": 11,
    "b": 21,
    "t": 21,
    "d~": -1,
    "u~": -1,
    "s~": -11,
    "c~": -11,
    "b~": -21,
    "t~": -21,
    "a": 0,
    "z": 0,
    "w+": 0,
    "w-": 0,
    "e+": -31,
    "e-": 31,
    "mu+": -41,
    "mu-": 41,
    "ta+": -51,
    "ta-": 51,
    "ve": 31,
    "ve~": -31,
    "vm": 41,
    "vm~": -41,
    "vt": 51,
    "vt~": -51,
    "h": 0,
}


@dataclass(frozen=True)
class ParticleSelectionMetadata:
    name: str
    pdg: int
    anti_name: str
    charge3: int
    family: int
    color_rep: int
    is_fermion: bool
    is_antiparticle: bool

    @property
    def is_colored(self) -> bool:
        return self.color_rep != 1


def _selection_metadata(name: str) -> ParticleSelectionMetadata:
    return _PARTICLE_SELECTION_METADATA[name]


def _metadata_color_rep(name: str) -> int:
    if name in QUARKS:
        return 3
    if name in ANTIQUARKS:
        return -3
    if name in GLUONS:
        return 8
    return 1


_PARTICLE_SELECTION_METADATA = {
    name: ParticleSelectionMetadata(
        name=name,
        pdg=int(PDGS[name]),
        anti_name=ANTI_PARTICLE[name],
        charge3=CHARGES3[name],
        family=FAMILY[name],
        color_rep=_metadata_color_rep(name),
        is_fermion=name in QUARKS or name in ANTIQUARKS or abs(int(PDGS[name])) in range(11, 17),
        is_antiparticle=name in ANTIQUARKS or int(PDGS[name]) < 0,
    )
    for name in PDGS
}


@dataclass(frozen=True)
class ProcessOptions:
    flavour_scheme: int = 5
    include_3qqbar: bool = False
    include_cc: bool = False
    include_resonance: bool = False
    serial: bool = True


@dataclass(frozen=True)
class ParsedProcess:
    initial_state: ProcessTuple
    jet_count: int
    rest: ProcessTuple
    leptons: ProcessTuple = ()


@dataclass(frozen=True)
class SubprocessRecord:
    process: ProcessTuple
    color_order: OrderTuple
    multichannel_partners: tuple[int, ...] = ()
    identical_factor: float = 1.0


@dataclass(frozen=True)
class PhaseSpaceGroup:
    group_id: int
    phase_space_order: OrderTuple
    records: tuple[SubprocessRecord, ...]


@dataclass(frozen=True)
class ProcessEnumeration:
    request: ParsedProcess
    options: ProcessOptions
    unique_processes: tuple[ProcessTuple, ...]
    groups: tuple[PhaseSpaceGroup, ...]

    @property
    def n_external(self) -> int:
        if self.unique_processes:
            return len(self.unique_processes[0])
        if self.groups and self.groups[0].records:
            return len(self.groups[0].records[0].process)
        return 0

    @property
    def n_records(self) -> int:
        return sum(len(group.records) for group in self.groups)


@dataclass(frozen=True)
class ProcessSetEntry:
    key: str
    process: str
    enumeration: ProcessEnumeration


@dataclass(frozen=True)
class ProcessSelectionRecord:
    source: str
    process: str
    key: str
    status: str
    reason: str
    quark_lines: int
    charge3: int
    family: int


@dataclass(frozen=True)
class ProcessSelectionReport:
    request: str
    options: ProcessOptions
    entries: tuple[ProcessSetEntry, ...]
    records: tuple[ProcessSelectionRecord, ...]
    candidate_count: int
    evaluated_count: int
    selected_count: int
    duplicate_count: int
    rejected_count: int
    rejection_counts: tuple[tuple[str, int], ...]
    stage_timings: tuple[tuple[str, float], ...]
    elapsed_s: float
    prefilter_enabled: bool

    @property
    def selected_records(self) -> tuple[ProcessSelectionRecord, ...]:
        return tuple(record for record in self.records if record.status == "selected")

    @property
    def duplicate_records(self) -> tuple[ProcessSelectionRecord, ...]:
        return tuple(record for record in self.records if record.status == "duplicate")


@dataclass(frozen=True)
class ProcessSetEnumeration:
    request: str
    options: ProcessOptions
    entries: tuple[ProcessSetEntry, ...]
    selection_report: ProcessSelectionReport | None = None

    @property
    def default_key(self) -> str:
        if not self.entries:
            raise ValueError("process set is empty")
        return self.entries[0].key


def split_process_set(process_string: str) -> tuple[str, ...]:
    """Split ``PROC | PROC`` input without treating bars inside brackets as separators."""

    parts: list[str] = []
    depth = 0
    start = 0
    for index, char in enumerate(process_string):
        if char == "[":
            depth += 1
        elif char == "]":
            depth -= 1
            if depth < 0:
                raise ValueError("unmatched ']' in process string")
        elif char == "|" and depth == 0:
            part = process_string[start:index].strip()
            if not part:
                raise ValueError("empty process in process set")
            parts.append(part)
            start = index + 1
    if depth != 0:
        raise ValueError("unmatched '[' in process string")
    tail = process_string[start:].strip()
    if not tail:
        raise ValueError("empty process in process set")
    parts.append(tail)
    return tuple(parts)


def expand_process_variants(process_string: str) -> tuple[str, ...]:
    """Expand anonymous multiparticle slots and repetition syntax.

    Built-in inclusive labels such as ``p`` and ``j`` are kept symbolic for the
    enumerator. Anonymous slots like ``[d g]`` are expanded by cartesian product,
    and each repeated slot in ``3*[d g]`` is treated independently.
    """

    variants: list[str] = []
    for process in split_process_set(process_string):
        parts = process.lower().replace("bar", "~").split(">")
        if len(parts) != 2:
            raise ValueError("invalid collision format; expected 'initial > final'")
        initial_options = _expand_side_tokens(_tokenize_side(parts[0].strip()))
        final_options = _expand_side_tokens(_tokenize_side(parts[1].strip()))
        for initial in itertools.product(*initial_options):
            for final in itertools.product(*final_options):
                variants.append(
                    f"{' '.join(initial)} > {' '.join(final)}"
                )
    return tuple(dict.fromkeys(variants))


def canonical_process_key(process: str) -> str:
    tokens = process.lower().replace("bar", "~").replace(">", " > ").split()
    safe = []
    for token in tokens:
        if token == ">":
            safe.append("to")
        else:
            safe.append(
                token.replace("~", "bar")
                .replace("+", "plus")
                .replace("-", "minus")
            )
    return "_".join(safe)


def _process_uses_inclusive_labels(process: str) -> bool:
    parts = process.lower().replace("bar", "~").split(">")
    if len(parts) != 2:
        raise ValueError("invalid collision format; expected 'initial > final'")
    for token in (*_tokenize_side(parts[0].strip()), *_tokenize_side(parts[1].strip())):
        _, item = _split_repeat_token(token)
        if item in {"p", "j"}:
            return True
        if item.startswith("[") and item.endswith("]"):
            if any(option in {"p", "j"} for option in _anonymous_options(item)):
                return True
    return False


def _request_allows_charged_current(request: ParsedProcess) -> bool:
    """Return whether request quantum numbers imply charged-current flow."""

    if "w+" in request.rest or "w-" in request.rest:
        return True
    non_qcd_rest = [
        particle
        for particle in request.rest
        if particle not in ALL_COLOURED and particle in CHARGES3
    ]
    return abs(sum(CHARGES3[particle] for particle in non_qcd_rest)) == 3


def _record_to_physical_process(record: SubprocessRecord) -> str:
    initial = tuple(ANTI_PARTICLE[p] for p in record.process[:2])
    final = record.process[2:]
    return f"{' '.join(initial)} > {' '.join(final)}"


def _concrete_processes_from_inclusive_enumeration(
    enumeration: ProcessEnumeration,
) -> tuple[str, ...]:
    processes: dict[str, None] = {}
    for group in enumeration.groups:
        for record in group.records:
            processes.setdefault(_record_to_physical_process(record), None)
    return tuple(sorted(processes, key=_concrete_process_sort_key))


def _concrete_process_sort_key(process: str) -> tuple[object, ...]:
    initial, _, final = process.partition(">")
    initial_tokens = tuple(_tokenize_side(initial.strip()))
    final_tokens = tuple(_tokenize_side(final.strip()))
    q_qbar_first = not (
        len(initial_tokens) == 2
        and not initial_tokens[0].endswith("~")
        and initial_tokens[1].endswith("~")
    )
    return (
        q_qbar_first,
        tuple(SORT_PARTICLES.get(token, 999) for token in initial_tokens),
        tuple(SORT_PARTICLES.get(token, 999) for token in final_tokens),
        initial_tokens,
        final_tokens,
    )


def _tokenize_side(side: str) -> list[str]:
    tokens: list[str] = []
    index = 0
    while index < len(side):
        if side[index].isspace():
            index += 1
            continue
        if side[index] == "[":
            end = side.find("]", index + 1)
            if end < 0:
                raise ValueError("unmatched '[' in process string")
            tokens.append(side[index : end + 1])
            index = end + 1
            continue
        end = index
        while end < len(side) and not side[end].isspace():
            if side[end] == "[":
                bracket = side.find("]", end + 1)
                if bracket < 0:
                    raise ValueError("unmatched '[' in process string")
                end = bracket + 1
            else:
                end += 1
        tokens.append(side[index:end])
        index = end
    return tokens


def _expand_side_tokens(tokens: Sequence[str]) -> tuple[tuple[str, ...], ...]:
    expanded: list[tuple[str, ...]] = []
    for token in tokens:
        repeat, item = _split_repeat_token(token)
        options = _anonymous_options(item)
        expanded.extend(options for _ in range(repeat))
    return tuple(expanded)


def _split_repeat_token(token: str) -> tuple[int, str]:
    match = re.fullmatch(r"(\d+)\*(.+)", token)
    if match:
        return int(match.group(1)), match.group(2)
    compact = re.fullmatch(r"(\d+)([A-Za-z][A-Za-z0-9+~\-]*)", token)
    if compact:
        return int(compact.group(1)), compact.group(2)
    return 1, token


def _anonymous_options(token: str) -> tuple[str, ...]:
    if token.startswith("[") and token.endswith("]"):
        options = tuple(_tokenize_side(token[1:-1].strip()))
        if not options:
            raise ValueError("anonymous multiparticle label cannot be empty")
        return options
    return (token,)


def _ordered_compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts <= 0:
        return ((),) if total == 0 else ()
    if parts == 1:
        return ((total,),)
    return tuple(
        (first, *rest)
        for first in range(total + 1)
        for rest in _ordered_compositions(total - first, parts - 1)
    )


def _chunk_sequence(
    values: Sequence[int],
    chunk_lengths: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    chunks: list[tuple[int, ...]] = []
    start = 0
    for length in chunk_lengths:
        chunks.append(tuple(values[start : start + length]))
        start += length
    return tuple(chunks)


class ProcessEnumerator:
    """Structured port of the legacy process_list.py enumeration logic."""

    def __init__(self, options: ProcessOptions | None = None) -> None:
        self.options = options or ProcessOptions()
        self.flavour_scheme = self._flavour_scheme(self.options.flavour_scheme)
        self.massless_qcd = self.flavour_scheme | frozenset(
            f"{q}~" for q in self.flavour_scheme
        ) | GLUONS
        self.proton = self.massless_qcd
        self.jet = self.massless_qcd
        self._phase_space_orders: dict[OrderTuple, list[SubprocessRecord]] = {}
        self._all_keys_sorted: list[OrderTuple] = []
        self._process_order_to_index: dict[tuple[ProcessTuple, OrderTuple], int] = {}

    def parse(self, process_string: str) -> ParsedProcess:
        variants = expand_process_variants(process_string)
        if len(variants) != 1:
            raise ValueError(
                "ProcessEnumerator.parse expects one concrete process; "
                "use enumerate_process_set for process sets or multiparticle expansion"
            )
        input_string = variants[0].lower().replace("bar", "~")
        parts = input_string.split(">")
        if len(parts) != 2:
            raise ValueError("invalid collision format; expected 'initial > final'")

        initial_state = _tokenize_side(parts[0].strip())
        if len(initial_state) != 2:
            raise ValueError("exactly two incoming particles are required")
        crossed_initial = tuple(
            ANTI_PARTICLE[p] if p not in {"p", "j"} else p for p in initial_state
        )

        final_state = _tokenize_side(parts[1].strip())
        jet_count = 0
        rest: list[str] = []
        for token in final_state:
            jet_match = re.fullmatch(r"(\d+)j", token)
            if token == "j":
                jet_count += 1
            elif jet_match:
                jet_count += int(jet_match.group(1))
            else:
                rest.append(token)
        leptons: list[str] = []
        if self.options.include_resonance:
            leptons = [p for p in rest if 11 <= abs(int(PDGS[p])) <= 16]
            for lepton in leptons:
                rest.remove(lepton)
            charge = sum(CHARGES3[p] for p in leptons)
            if charge == 0:
                rest.append("z")
            elif charge < 0:
                rest.append("w-")
            else:
                rest.append("w+")

        self._validate_particles([*crossed_initial, *rest])
        return ParsedProcess(crossed_initial, jet_count, tuple(rest), tuple(leptons))

    def enumerate(self, process_string: str) -> ProcessEnumeration:
        request = self.parse(process_string)
        unique_processes = self._generate_all_unique_processes(request)
        subprocesses = self._generate_all_processes(unique_processes, request)
        self._phase_space_orders = self._combine_results(
            self._process_subprocess(proc) for proc in subprocesses
        )
        self._all_keys_sorted = sorted(self._phase_space_orders)
        self._determine_multichannel_partners_and_symmetry_factor()
        self._check_consistency()

        groups = tuple(
            PhaseSpaceGroup(
                group_id=i + 1,
                phase_space_order=key,
                records=tuple(
                    sorted(self._phase_space_orders[key], key=self._sort_record)
                ),
            )
            for i, key in enumerate(self._all_keys_sorted)
        )
        return ProcessEnumeration(
            request=request,
            options=self.options,
            unique_processes=tuple(
                tuple(proc)
                for proc in sorted(
                    (tuple(sorted(p, key=lambda x: SORT_PARTICLES[x])) for p in unique_processes),
                    key=self._sort_process,
                )
            ),
            groups=groups,
        )

    def enumerate_color_complete(self, process_string: str) -> ProcessEnumeration:
        """Enumerate a reference-only colour-complete legacy process file.

        The production LC process list deliberately removes colour-order
        representatives related by legacy symmetries.  That is correct for
        ordinary integration and LC library timing, but NLC/full-colour
        generated-library validation needs raw amplitudes for every colour
        basis row used by AmpliCol's colour matrix.  This reference path keeps
        all model-compatible candidate colour words in a single phase-space
        group and avoids process-family assumptions.
        """

        request = self.parse(process_string)
        unique_processes = self._generate_all_unique_processes(request)
        subprocesses = self._generate_all_processes(unique_processes, request)
        records: list[SubprocessRecord] = []
        seen: set[tuple[ProcessTuple, OrderTuple]] = set()
        for proc in sorted(subprocesses, key=self._sort_process):
            for perm in self._candidate_color_orders(proc):
                order = tuple(perm)
                key = (proc, order)
                if key in seen:
                    continue
                seen.add(key)
                records.append(
                    SubprocessRecord(
                        process=proc,
                        color_order=order,
                        multichannel_partners=(0,),
                        identical_factor=self._identical_particle_symmetry_factor(proc),
                    )
                )
        if not records:
            return ProcessEnumeration(
                request=request,
                options=self.options,
                unique_processes=(),
                groups=(),
            )
        return ProcessEnumeration(
            request=request,
            options=self.options,
            unique_processes=tuple(
                tuple(proc)
                for proc in sorted(
                    (
                        tuple(sorted(p, key=lambda x: SORT_PARTICLES[x]))
                        for p in unique_processes
                    ),
                    key=self._sort_process,
                )
            ),
            groups=(
                PhaseSpaceGroup(
                    group_id=1,
                    phase_space_order=tuple(range(len(records[0].process))),
                    records=tuple(records),
                ),
            ),
        )

    def to_legacy_lines(self, enumeration: ProcessEnumeration) -> list[str]:
        unique_lines = self._unique_process_lines(enumeration)
        group_lines = self._phase_space_group_lines(enumeration)
        return [*unique_lines, *group_lines]

    def write_legacy_file(self, enumeration: ProcessEnumeration, path: str | Path) -> None:
        Path(path).write_text("\n".join(self.to_legacy_lines(enumeration)) + "\n")

    def _flavour_scheme(self, value: int) -> frozenset[str]:
        flavours = ("d", "u", "s", "c", "b", "t")
        if not 1 <= value <= 6:
            raise ValueError(f"unknown flavour scheme: {value}")
        return frozenset(flavours[:value])

    def _validate_particles(self, particles: Iterable[str]) -> None:
        unknown = sorted(set(particles).difference(PDGS).difference({"p", "j"}))
        if unknown:
            raise ValueError(f"unknown particle name(s): {unknown}")

    def _valid_color_order(self, proc: ProcessTuple, perm: OrderTuple) -> bool:
        found_quark = found_antiquark = found_singlet = found_gluon = found_first = False
        quark_idx = -1
        for idx in perm:
            if idx == 0:
                found_first = True
            particle = proc[idx]
            if particle in QUARKS:
                if found_quark or found_gluon:
                    return False
                found_quark = True
                found_antiquark = found_singlet = found_gluon = False
                quark_idx = idx
            elif particle in ANTIQUARKS:
                if found_antiquark or found_singlet or not found_quark:
                    return False
                if not found_first:
                    return False
                found_antiquark = True
                found_quark = found_singlet = found_gluon = False
                _ = quark_idx
            elif particle in GLUONS:
                if found_antiquark or found_singlet:
                    return False
                found_gluon = True
            else:
                if found_quark or found_gluon or not found_antiquark:
                    return False
                found_singlet = True
        return not (found_gluon and perm[0] != 0)

    def _unique_color_order(self, proc: ProcessTuple, perm: OrderTuple) -> bool:
        zero = perm.index(0)
        perm_mapped = perm[zero:] + perm[:zero]
        for particle in ALL_COLOURED:
            particle_positions = [
                i + 2 for i, part in enumerate(proc[2:]) if part == particle
            ]
            previous_position = 0
            for position in particle_positions:
                mapped_position = perm_mapped.index(position)
                if mapped_position < previous_position:
                    return False
                previous_position = mapped_position
        return True

    def _order_proc_perm(
        self, proc: ProcessTuple, perm: OrderTuple
    ) -> tuple[ProcessTuple, OrderTuple]:
        zero = perm.index(0)
        perm_mapped = list(perm[zero:] + perm[:zero])
        elements_to_order = [
            i for i in perm_mapped if proc[i] in self.massless_qcd and i > 1
        ]
        indices = [perm_mapped.index(x) for x in elements_to_order]
        for idx, value in zip(indices, sorted(elements_to_order), strict=True):
            perm_mapped[idx] = value
        perm_ordered = tuple(
            perm_mapped[len(perm) - zero :] + perm_mapped[: len(perm) - zero]
        )
        proc_ordered: list[str | None] = [None] * len(proc)
        for i, perm_index in enumerate(perm_ordered):
            proc_ordered[perm_index] = proc[perm[i]]
        return tuple(x for x in proc_ordered if x is not None), perm_ordered

    def _process_subprocess(
        self, proc: ProcessTuple
    ) -> dict[OrderTuple, list[SubprocessRecord]]:
        local: dict[OrderTuple, list[SubprocessRecord]] = {}
        for perm in self._candidate_color_orders(proc):
            order = tuple(perm)
            if not self._valid_color_order(proc, order):
                continue
            if not self._unique_color_order(proc, order):
                continue
            ordered_proc, ordered_perm = self._order_proc_perm(proc, order)
            zero = ordered_perm.index(0)
            perm_mapped = tuple(ordered_perm[zero:] + ordered_perm[:zero])
            local.setdefault(perm_mapped, []).append(
                SubprocessRecord(ordered_proc, ordered_perm)
            )
        return local

    def _candidate_color_orders(self, proc: ProcessTuple) -> Iterable[OrderTuple]:
        quark_indices = tuple(i for i, particle in enumerate(proc) if particle in QUARKS)
        anti_indices = tuple(i for i, particle in enumerate(proc) if particle in ANTIQUARKS)
        gluon_indices = tuple(i for i, particle in enumerate(proc) if particle in GLUONS)
        singlet_indices = tuple(i for i, particle in enumerate(proc) if particle in SINGLETS)
        if not quark_indices:
            for tail in itertools.permutations(index for index in gluon_indices if index != 0):
                yield (0, *tail)
            return

        if singlet_indices:
            singlet_perms = tuple(itertools.permutations(singlet_indices))
        else:
            singlet_perms = ((),)
        gluon_compositions = _ordered_compositions(
            len(gluon_indices),
            len(quark_indices),
        )
        singlet_compositions = _ordered_compositions(
            len(singlet_indices),
            len(quark_indices),
        )
        for quark_order in itertools.permutations(quark_indices):
            for anti_order in itertools.permutations(anti_indices):
                for gluon_perm in itertools.permutations(gluon_indices):
                    for gluon_chunks in gluon_compositions:
                        gluon_by_line = _chunk_sequence(gluon_perm, gluon_chunks)
                        for singlet_perm in singlet_perms:
                            for singlet_chunks in singlet_compositions:
                                singlet_by_line = _chunk_sequence(
                                    singlet_perm,
                                    singlet_chunks,
                                )
                                order: list[int] = []
                                for line in range(len(quark_indices)):
                                    order.append(quark_order[line])
                                    order.extend(gluon_by_line[line])
                                    order.append(anti_order[line])
                                    order.extend(singlet_by_line[line])
                                yield tuple(order)

    def _valid_process(
        self,
        proc: Sequence[str],
        *,
        allow_charged_current: bool | None = None,
    ) -> bool:
        allow_cc = (
            self.options.include_cc
            if allow_charged_current is None
            else allow_charged_current
        )
        nq = self._count_matching(proc, QUARKS)
        naq = self._count_matching(proc, ANTIQUARKS)
        if nq != naq:
            return False
        if sum(CHARGES3[x] for x in proc) != 0:
            return False
        if not self.options.include_cc and sum(FAMILY[x] for x in proc) != 0:
            return False
        if not allow_cc:
            if "w+" not in proc and "w-" not in proc:
                for q in QUARKS:
                    if self._count_matching(proc, [q]) != self._count_matching(
                        proc, [f"{q}~"]
                    ):
                        return False
        return not (nq == 0 and self._count_matching(proc, SINGLETS) > 0)

    def _compatible_unique_process(
        self, request: ParsedProcess, proc: Sequence[str]
    ) -> bool:
        mandatory = [p for p in request.initial_state if p not in {"p", "j"}]
        mandatory.extend(request.rest)
        proc_local = list(proc)
        try:
            for particle in mandatory:
                proc_local.remove(particle)
        except ValueError:
            return False
        return True

    def _compatible_process(self, request: ParsedProcess, proc: ProcessTuple) -> bool:
        proc_local = list(proc[2:])
        for i in (0, 1):
            if request.initial_state[i] not in {"p", "j"} and proc[i] != request.initial_state[i]:
                return False
        try:
            for particle in request.rest:
                proc_local.remove(particle)
        except ValueError:
            return False
        return True

    def _generate_all_unique_processes(self, request: ParsedProcess) -> set[ProcessTuple]:
        processes: list[list[str]] = [[]]
        for part in request.initial_state:
            if part not in {"p", "j"} and part not in self.massless_qcd:
                raise ValueError("initial state should be a proton or massless QCD parton")

        if request.jet_count == 0 and all(
            part not in {"p", "j"} for part in request.initial_state
        ):
            proc = tuple(
                sorted(
                    (*request.initial_state, *request.rest),
                    key=lambda item: SORT_PARTICLES[item],
                )
            )
            if self._valid_process(
                proc,
                allow_charged_current=_request_allows_charged_current(request),
            ) and self._compatible_unique_process(request, proc):
                return {proc}
            return set()

        qcd_rest = sum(1 for part in request.rest if part in self.massless_qcd)
        for part_index in range(request.jet_count + 2 + qcd_rest):
            new_processes: list[list[str]] = []
            for proc in processes:
                for particle in self.massless_qcd:
                    new_processes.append(sorted([*proc, particle]))
            processes = new_processes

        for part in request.rest:
            if part in self.massless_qcd:
                continue
            for proc in processes:
                proc.append(part)

        return {
            tuple(proc)
            for proc in processes
            if self._valid_process(
                proc,
                allow_charged_current=_request_allows_charged_current(request),
            )
            and self._compatible_unique_process(request, proc)
        }

    def _generate_all_processes(
        self, unique_processes: set[ProcessTuple], request: ParsedProcess
    ) -> set[ProcessTuple]:
        processes: set[ProcessTuple] = set()
        for proc in unique_processes:
            for i, j in itertools.combinations(range(len(proc)), 2):
                if proc[i] in self.jet and proc[j] in self.jet:
                    remaining = [entry for k, entry in enumerate(proc) if k not in (i, j)]
                    proc1 = (proc[i], proc[j], *remaining)
                    proc2 = (proc[j], proc[i], *remaining)
                    if self._compatible_process(request, proc1):
                        processes.add(proc1)
                    if self._compatible_process(request, proc2):
                        processes.add(proc2)
        return processes

    def _combine_results(
        self, results: Iterable[dict[OrderTuple, list[SubprocessRecord]]]
    ) -> dict[OrderTuple, list[SubprocessRecord]]:
        combined: dict[OrderTuple, list[SubprocessRecord]] = {}
        for result in results:
            for key, value in result.items():
                combined.setdefault(key, []).extend(value)
        return combined

    def _determine_multichannel_partners_and_symmetry_factor(self) -> None:
        self._build_process_index()
        for key in self._all_keys_sorted:
            for index, record in enumerate(tuple(self._phase_space_orders[key])):
                self._phase_space_orders[key][index] = self._record_with_partners(record)

        for key in self._all_keys_sorted:
            for index, record in enumerate(tuple(self._phase_space_orders[key])):
                self._phase_space_orders[key][index] = replace(
                    record,
                    identical_factor=record.identical_factor
                    * self._identical_particle_symmetry_factor(record.process),
                )

    def _build_process_index(self) -> None:
        self._process_order_to_index = {}
        for i, key in enumerate(self._all_keys_sorted):
            for record in self._phase_space_orders[key]:
                self._process_order_to_index[(record.process, record.color_order)] = i

    def _record_with_partners(self, record: SubprocessRecord) -> SubprocessRecord:
        proc = record.process
        perm = record.color_order
        singlet_indices = [perm.index(i) for i, p in enumerate(proc) if p in SINGLETS]
        anti_quark_indices = tuple(
            perm.index(i) for i, p in enumerate(proc) if p in ANTIQUARKS
        )
        nsinglets = len(singlet_indices)
        if len(singlet_indices) > 1:
            singlet_perms = tuple(itertools.permutations(singlet_indices))
        elif len(singlet_indices) == 1:
            singlet_perms = (tuple(singlet_indices),)
        else:
            singlet_perms = ()

        iden = float(math.factorial(max(len(anti_quark_indices) - 1, 0)))
        possible: list[tuple[OrderTuple, ProcessTuple]] = []
        if not singlet_perms:
            possible.append((perm, proc))
        elif anti_quark_indices:
            for singlets in singlet_perms:
                for chunks in _ordered_compositions(
                    nsinglets,
                    len(anti_quark_indices),
                ):
                    starts = [0]
                    for chunk in chunks[:-1]:
                        starts.append(starts[-1] + chunk)
                    singlet_chunks = tuple(
                        singlets[start : start + chunk]
                        for start, chunk in zip(starts, chunks, strict=True)
                    )
                    order_list: list[int] = []
                    for i in range(len(perm)):
                        if i in anti_quark_indices:
                            anti_index = anti_quark_indices.index(i)
                            order_list.extend(
                                perm[p]
                                for p in (
                                    (anti_quark_indices[anti_index],)
                                    + singlet_chunks[anti_index]
                                )
                            )
                        elif i not in singlet_indices:
                            order_list.append(perm[i])
                    possible.append((tuple(order_list), proc))
        else:
            raise RuntimeError("colour-singlet partners require at least one quark line")

        partners: list[int] = []
        for possible_order, process in possible:
            partner = self._process_order_to_index.get((process, possible_order))
            if partner is None:
                raise RuntimeError(
                    f"expected multichannel partner not found for {process} {possible_order}"
                )
            if partner in partners:
                raise RuntimeError(
                    f"duplicate multichannel partner for {process} {possible_order}"
                )
            partners.append(partner)
        return replace(
            record,
            multichannel_partners=tuple(sorted(partners)),
            identical_factor=1.0 / iden,
        )

    def _identical_particle_symmetry_factor(self, proc: ProcessTuple) -> float:
        factor = 1.0
        for particle in ALL_COLOURED:
            factor *= max(1, math.factorial(proc[2:].count(particle)))
        return factor

    def _check_consistency(self) -> None:
        all_processes: dict[ProcessTuple, float] = {}
        for key in self._all_keys_sorted:
            for record in self._phase_space_orders[key]:
                proc = (record.process[0], record.process[1], *sorted(record.process[2:], key=lambda x: SORT_PARTICLES[x]))
                all_processes[proc] = all_processes.get(proc, 0.0) + (
                    record.identical_factor / len(record.multichannel_partners)
                )
        for proc, count in all_processes.items():
            expected = self.expected_number_of_dual_amplitudes(proc)
            if abs(count - expected) > 1e-5:
                raise RuntimeError(
                    f"inconsistent number of dual amplitudes for {proc}: {count}, expected {expected}"
                )

    def expected_number_of_dual_amplitudes(self, proc: Sequence[str]) -> float:
        nq = self._count_matching(proc, QUARKS)
        ng = self._count_matching(proc, GLUONS)
        if nq == 0:
            return math.factorial(ng - 1)
        return (
            math.factorial(ng)
            * math.comb(ng + nq - 1, nq - 1)
            * math.factorial(nq)
        )

    def _unique_process_lines(self, enumeration: ProcessEnumeration) -> list[str]:
        sorted_processes = [list(proc) for proc in enumeration.unique_processes]
        if enumeration.options.include_resonance:
            for proc in sorted_processes:
                if "z" in proc:
                    index = proc.index("z")
                    proc[index : index + 1] = enumeration.request.leptons
        sorted_processes = self._add_different_flavour_quark_line_processes(sorted_processes)
        if not sorted_processes:
            raise ValueError("no processes found")

        lines = [f"{len(sorted_processes[0])} {len(sorted_processes)}"]
        lines.extend(" ".join(PDGS[p] for p in proc) for proc in sorted_processes)
        lines.extend(["", ""])
        return lines

    def _phase_space_group_lines(self, enumeration: ProcessEnumeration) -> list[str]:
        lines = [str(len(enumeration.groups)), ""]
        groups = list(enumeration.groups)
        if enumeration.options.include_resonance:
            groups = self._resonance_adjusted_groups(enumeration)

        for group in groups:
            max_channels = max(len(record.multichannel_partners) for record in group.records)
            order = " ".join(str(k + 1) for k in group.phase_space_order)
            lines.append(
                f"{group.group_id}   {len(group.records)}   {max_channels}   {order}"
            )
            for record in sorted(group.records, key=self._sort_record):
                lines.append(self._convert_record_to_legacy_line(record))
            lines.extend(["", "", ""])
        return lines

    def _resonance_adjusted_groups(
        self, enumeration: ProcessEnumeration
    ) -> list[PhaseSpaceGroup]:
        adjusted = []
        for group in enumeration.groups:
            first = group.records[0]
            phase_order = group.phase_space_order
            records = []
            if "z" in first.process:
                z_index = first.process.index("z")
                j = phase_order.index(z_index)
                phase_order = phase_order[:j] + (z_index, z_index + 1) + phase_order[j + 1 :]
            for record in group.records:
                if "z" not in record.process:
                    records.append(record)
                    continue
                z_index = record.process.index("z")
                process = (
                    record.process[:z_index]
                    + enumeration.request.leptons
                    + record.process[z_index + 1 :]
                )
                j = record.color_order.index(z_index)
                color_order = (
                    record.color_order[:j]
                    + (z_index, z_index + 1)
                    + record.color_order[j + 1 :]
                )
                records.append(
                    replace(record, process=process, color_order=color_order)
                )
            adjusted.append(replace(group, phase_space_order=phase_order, records=tuple(records)))
        return adjusted

    def _convert_record_to_legacy_line(self, record: SubprocessRecord) -> str:
        crossed = [
            PDGS[p] if i > 1 else PDGS[ANTI_PARTICLE[p]]
            for i, p in enumerate(record.process)
        ]
        parts = [
            str(len(record.multichannel_partners)),
            " ".join(str(m + 1) for m in record.multichannel_partners),
            " ".join(crossed),
            " ".join(str(o + 1) for o in record.color_order),
            str(record.identical_factor),
        ]
        return "   ".join(parts)

    def _add_different_flavour_quark_line_processes(
        self, sorted_processes: list[list[str]]
    ) -> list[list[str]]:
        index = 0
        while index < len(sorted_processes):
            proc = sorted_processes[index]
            nq = self._count_matching(proc, QUARKS)
            if nq >= 2:
                quarks = proc[0:nq]
                antiquarks = proc[nq : 2 * nq]
                for q_perm in itertools.permutations(quarks):
                    for aq_perm in itertools.permutations(antiquarks):
                        swapped = [*q_perm, *aq_perm, *proc[2 * nq :]]
                        if swapped not in sorted_processes:
                            sorted_processes.insert(index + 1, swapped)
                            index += 1
            index += 1
        return sorted_processes

    def _sort_process(self, process: Sequence[str]) -> tuple[int, int, list[int]]:
        nq = self._count_matching(process, QUARKS)
        same_flavour = max(self._count_matching(process, [q]) for q in QUARKS)
        return nq, same_flavour, [SORT_PARTICLES[p] for p in process]

    def _sort_record(self, record: SubprocessRecord) -> tuple[int, int, list[int]]:
        return self._sort_process(record.process)

    def _count_matching(self, main: Sequence[str], check: Iterable[str]) -> int:
        counts = Counter(main)
        return sum(counts[item] for item in check if item in counts)


def enumerate_processes(
    process_string: str, options: ProcessOptions | None = None
) -> ProcessEnumeration:
    return ProcessEnumerator(options).enumerate(process_string)


def enumerate_process_set(
    process_string: str,
    options: ProcessOptions | None = None,
) -> ProcessSetEnumeration:
    active_options = options or ProcessOptions()
    entries: list[ProcessSetEntry] = []
    seen: set[str] = set()
    for process in expand_process_variants(process_string):
        inclusive_enumeration: ProcessEnumeration | None = None
        if _process_uses_inclusive_labels(process):
            inclusive_enumeration = ProcessEnumerator(active_options).enumerate(process)
            concrete_processes = _concrete_processes_from_inclusive_enumeration(
                inclusive_enumeration
            )
        else:
            concrete_processes = (process,)
        for concrete_process in concrete_processes:
            key = canonical_process_key(concrete_process)
            if key in seen:
                continue
            seen.add(key)
            enumeration = ProcessEnumerator(active_options).enumerate(concrete_process)
            if enumeration.n_records == 0:
                continue
            entries.append(
                ProcessSetEntry(
                    key=key,
                    process=concrete_process,
                    enumeration=enumeration,
                )
            )
    if not entries:
        raise ValueError(f"no valid processes found for {process_string!r}")
    return ProcessSetEnumeration(
        request=process_string,
        options=active_options,
        entries=tuple(entries),
    )


def enumerate_generic_process_set(
    process_string: str,
    options: ProcessOptions | None = None,
    *,
    max_quark_pairs: int | None = None,
    use_prefilter: bool = True,
) -> ProcessSetEnumeration:
    """Expand process-set syntax for generic DAG generation.

    The generic compiler only needs concrete subprocess keys and parser
    metadata.  It should not pay the legacy phase-space/color-order enumeration
    cost for an already concrete partonic process; colour sectors and closures
    are discovered later by the model-driven DAG.  Inclusive ``p``/``j``
    requests still use the legacy enumerator once to obtain concrete children.
    """

    active_options = options or ProcessOptions()
    report = build_generic_process_selection_report(
        process_string,
        active_options,
        max_quark_pairs=max_quark_pairs,
        use_prefilter=use_prefilter,
    )
    if not report.entries:
        raise ValueError(f"no valid processes found for {process_string!r}")
    return ProcessSetEnumeration(
        request=process_string,
        options=active_options,
        entries=report.entries,
        selection_report=report,
    )


def build_generic_process_selection_report(
    process_string: str,
    options: ProcessOptions | None = None,
    *,
    max_quark_pairs: int | None = None,
    use_prefilter: bool = True,
) -> ProcessSelectionReport:
    """Return selected generic subprocesses plus enumeration diagnostics."""

    started = time.perf_counter()
    active_options = options or ProcessOptions()
    entries: list[ProcessSetEntry] = []
    seen: set[str] = set()
    records: list[ProcessSelectionRecord] = []
    rejection_counts: Counter[str] = Counter()
    candidate_count = 0
    evaluated_count = 0
    stage_timings: list[tuple[str, float]] = []

    stage_started = time.perf_counter()
    variants = expand_process_variants(process_string)
    stage_timings.append(("expand", time.perf_counter() - stage_started))

    stage_started = time.perf_counter()
    for process in variants:
        if use_prefilter:
            result = _generic_selection_records_from_request(
                process,
                active_options,
                max_quark_pairs=max_quark_pairs,
                seen=seen,
            )
        else:
            result = _legacy_selection_records_from_request(
                process,
                active_options,
                max_quark_pairs=max_quark_pairs,
                seen=seen,
            )
        candidate_count += result[0]
        evaluated_count += result[1]
        rejection_counts.update(result[2])
        records.extend(result[3])
    stage_timings.append(("prefilter", time.perf_counter() - stage_started))

    stage_started = time.perf_counter()
    for record in records:
        if record.status != "selected":
            continue
        enumeration = _lightweight_concrete_process_enumeration(
            record.process,
            active_options,
        )
        if enumeration.n_records == 0:
            rejection_counts["model-reachability"] += 1
            continue
        entries.append(
            ProcessSetEntry(
                key=record.key,
                process=record.process,
                enumeration=enumeration,
            )
        )
    stage_timings.append(("canonicalize", time.perf_counter() - stage_started))
    duplicate_count = sum(1 for record in records if record.status == "duplicate")
    selected_count = len(entries)
    rejected_count = candidate_count - selected_count - duplicate_count
    return ProcessSelectionReport(
        request=process_string,
        options=active_options,
        entries=tuple(entries),
        records=tuple(records),
        candidate_count=candidate_count,
        evaluated_count=evaluated_count,
        selected_count=selected_count,
        duplicate_count=duplicate_count,
        rejected_count=max(rejected_count, 0),
        rejection_counts=tuple(sorted(rejection_counts.items())),
        stage_timings=tuple(stage_timings),
        elapsed_s=time.perf_counter() - started,
        prefilter_enabled=use_prefilter,
    )


def _generic_selection_records_from_request(
    process: str,
    options: ProcessOptions,
    *,
    max_quark_pairs: int | None,
    seen: set[str],
) -> tuple[int, int, Counter[str], list[ProcessSelectionRecord]]:
    if _process_uses_inclusive_labels(process):
        return _prefilter_inclusive_request(
            process,
            options,
            max_quark_pairs=max_quark_pairs,
            seen=seen,
        )
    return _classify_concrete_request(
        process,
        options,
        max_quark_pairs=max_quark_pairs,
        seen=seen,
        source=process,
        symmetric_initial=False,
    )


def _legacy_selection_records_from_request(
    process: str,
    options: ProcessOptions,
    *,
    max_quark_pairs: int | None,
    seen: set[str],
) -> tuple[int, int, Counter[str], list[ProcessSelectionRecord]]:
    if _process_uses_inclusive_labels(process):
        enumeration = ProcessEnumerator(options).enumerate(process)
        concrete_processes = _concrete_processes_from_inclusive_enumeration(enumeration)
    else:
        concrete_processes = (process,)
    records: list[ProcessSelectionRecord] = []
    rejection_counts: Counter[str] = Counter()
    candidate_count = 0
    evaluated_count = 0
    for concrete_process in concrete_processes:
        result = _classify_concrete_request(
            concrete_process,
            options,
            max_quark_pairs=max_quark_pairs,
            seen=seen,
            source=process,
            symmetric_initial=False,
        )
        candidate_count += result[0]
        evaluated_count += result[1]
        rejection_counts.update(result[2])
        records.extend(result[3])
    return candidate_count, evaluated_count, rejection_counts, records


def _classify_concrete_request(
    process: str,
    options: ProcessOptions,
    *,
    max_quark_pairs: int | None,
    seen: set[str],
    source: str,
    symmetric_initial: bool,
) -> tuple[int, int, Counter[str], list[ProcessSelectionRecord]]:
    enumerator = ProcessEnumerator(options)
    request = enumerator.parse(process)
    physical_initial = tuple(ANTI_PARTICLE[p] for p in request.initial_state)
    final_state = tuple(request.rest)
    return _classify_physical_candidate(
        enumerator,
        physical_initial,
        final_state,
        source=source,
        allow_charged_current=_request_allows_charged_current(request),
        max_quark_pairs=max_quark_pairs,
        seen=seen,
        symmetric_initial=symmetric_initial,
    )


def _prefilter_inclusive_request(
    process: str,
    options: ProcessOptions,
    *,
    max_quark_pairs: int | None,
    seen: set[str],
) -> tuple[int, int, Counter[str], list[ProcessSelectionRecord]]:
    enumerator = ProcessEnumerator(options)
    request = enumerator.parse(process)
    initial_options: list[tuple[str, ...]] = []
    parton_options = tuple(
        sorted(enumerator.massless_qcd, key=lambda p: SORT_PARTICLES[p])
    )
    for crossed_particle in request.initial_state:
        if crossed_particle in {"p", "j"}:
            initial_options.append(parton_options)
        else:
            initial_options.append((ANTI_PARTICLE[crossed_particle],))

    symmetric_initial = (
        len(request.initial_state) == 2
        and request.initial_state[0] == request.initial_state[1]
        and request.initial_state[0] in {"p", "j"}
    )
    final_candidates = tuple(
        _generic_final_candidate(request.rest, jets)
        for jets in itertools.combinations_with_replacement(
            parton_options,
            request.jet_count,
        )
    )
    final_by_charge: Counter[int] = Counter(
        final_charge for _, final_charge, _, _, _ in final_candidates
    )
    final_by_charge_family: Counter[tuple[int, int]] = Counter(
        (final_charge, final_family)
        for _, final_charge, final_family, _, _ in final_candidates
    )
    bucketed_finals: dict[tuple[int, int | None], list[tuple[ProcessTuple, int, int, int, int]]] = {}
    for candidate in final_candidates:
        final_state, final_charge, final_family, _, _ = candidate
        del final_state
        bucketed_finals.setdefault((final_charge, None), []).append(candidate)
        bucketed_finals.setdefault((final_charge, final_family), []).append(candidate)

    allow_charged_current = _request_allows_charged_current(request)
    records: list[ProcessSelectionRecord] = []
    rejection_counts: Counter[str] = Counter()
    candidate_count = 0
    evaluated_count = 0
    for initial_state in itertools.product(*initial_options):
        crossed_initial = tuple(ANTI_PARTICLE[p] for p in initial_state)
        initial_charge = _charge3_sum(crossed_initial)
        initial_family = _family_sum(crossed_initial)
        candidate_count += len(final_candidates)
        charge_key = -initial_charge
        charge_matches = final_by_charge[charge_key]
        rejection_counts["charge"] += len(final_candidates) - charge_matches
        if options.include_cc:
            candidate_bucket = bucketed_finals.get((charge_key, None), ())
        else:
            family_key = -initial_family
            family_matches = final_by_charge_family[(charge_key, family_key)]
            rejection_counts["fermion-family"] += charge_matches - family_matches
            candidate_bucket = bucketed_finals.get((charge_key, family_key), ())
        for final_state, _, _, _, _ in candidate_bucket:
            result = _classify_physical_candidate(
                enumerator,
                initial_state,
                final_state,
                source=process,
                allow_charged_current=allow_charged_current,
                max_quark_pairs=max_quark_pairs,
                seen=seen,
                symmetric_initial=symmetric_initial,
                prechecked_charge_family=True,
            )
            candidate_count += result[0] - 1
            evaluated_count += result[1]
            rejection_counts.update(result[2])
            records.extend(result[3])
    return candidate_count, evaluated_count, rejection_counts, records


def _classify_physical_candidate(
    enumerator: ProcessEnumerator,
    physical_initial: Sequence[str],
    final_state: Sequence[str],
    *,
    source: str,
    allow_charged_current: bool,
    max_quark_pairs: int | None,
    seen: set[str],
    symmetric_initial: bool,
    prechecked_charge_family: bool = False,
) -> tuple[int, int, Counter[str], list[ProcessSelectionRecord]]:
    crossed_initial = tuple(ANTI_PARTICLE[p] for p in physical_initial)
    display_process = _canonical_physical_process(
        physical_initial,
        final_state,
        symmetric_initial=symmetric_initial,
    )
    canonical_initial, _, canonical_final = display_process.partition(">")
    canonical_crossed_initial = tuple(
        ANTI_PARTICLE[p] for p in _tokenize_side(canonical_initial.strip())
    )
    canonical_final_state = tuple(_tokenize_side(canonical_final.strip()))
    all_outgoing = (*canonical_crossed_initial, *canonical_final_state)
    charge3 = _charge3_sum(all_outgoing)
    family = _family_sum(all_outgoing)
    quark_lines = _quark_pair_count(all_outgoing)
    records: list[ProcessSelectionRecord] = []
    rejection_counts: Counter[str] = Counter()
    evaluated_count = 1
    reason = ""
    if not prechecked_charge_family and charge3 != 0:
        reason = "charge"
    elif not prechecked_charge_family and not enumerator.options.include_cc and family != 0:
        reason = "fermion-family"
    elif max_quark_pairs is not None and quark_lines > int(max_quark_pairs):
        reason = "max-quark-lines"
    elif _single_impossible_colored_state(all_outgoing):
        reason = "single-coloured-state"
    elif not enumerator._valid_process(
        all_outgoing,
        allow_charged_current=allow_charged_current,
    ):
        reason = "model-reachability"

    key = canonical_process_key(display_process)
    dedupe_key = _process_dedupe_signature(
        physical_initial,
        final_state,
        symmetric_initial=symmetric_initial,
    )
    if reason:
        rejection_counts[reason] += 1
        return 1, evaluated_count, rejection_counts, []
    if dedupe_key in seen:
        records.append(
            ProcessSelectionRecord(
                source=source,
                process=display_process,
                key=key,
                status="duplicate",
                reason="canonical-duplicate",
                quark_lines=quark_lines,
                charge3=charge3,
                family=family,
            )
        )
        return 1, evaluated_count, rejection_counts, records
    seen.add(dedupe_key)
    records.append(
        ProcessSelectionRecord(
            source=source,
            process=display_process,
            key=key,
            status="selected",
            reason="selected",
            quark_lines=quark_lines,
            charge3=charge3,
            family=family,
        )
    )
    return 1, evaluated_count, rejection_counts, records


def _canonical_physical_process(
    initial_state: Sequence[str],
    final_state: Sequence[str],
    *,
    symmetric_initial: bool,
) -> str:
    initial = tuple(initial_state)
    if symmetric_initial:
        initial = tuple(sorted(initial, key=lambda particle: SORT_PARTICLES[particle]))
    final = tuple(final_state)
    return f"{' '.join(initial)} > {' '.join(final)}"


def _process_dedupe_signature(
    initial_state: Sequence[str],
    final_state: Sequence[str],
    *,
    symmetric_initial: bool,
) -> str:
    initial = tuple(initial_state)
    if symmetric_initial:
        initial = tuple(sorted(initial, key=lambda particle: SORT_PARTICLES[particle]))
    final = tuple(sorted(final_state, key=lambda particle: SORT_PARTICLES[particle]))
    return f"{' '.join(initial)} > {' '.join(final)}"


def _single_impossible_colored_state(process: Sequence[str]) -> bool:
    return sum(1 for particle in process if _selection_metadata(particle).is_colored) == 1


def _quark_pair_count(process: Sequence[str]) -> int:
    quarks, antiquarks = _quark_counts(process)
    return min(quarks, antiquarks)


def _charge3_sum(process: Sequence[str]) -> int:
    return sum(_selection_metadata(particle).charge3 for particle in process)


def _family_sum(process: Sequence[str]) -> int:
    return sum(_selection_metadata(particle).family for particle in process)


def _quark_counts(process: Sequence[str]) -> tuple[int, int]:
    counts = Counter(process)
    quarks = sum(counts[item] for item in QUARKS if item in counts)
    antiquarks = sum(counts[item] for item in ANTIQUARKS if item in counts)
    return quarks, antiquarks


def _generic_final_candidate(
    rest: Sequence[str],
    jets: Sequence[str],
) -> tuple[ProcessTuple, int, int, int, int]:
    final_state = tuple(
        sorted(
            (*rest, *jets),
            key=lambda particle: SORT_PARTICLES[particle],
        )
    )
    quarks, antiquarks = _quark_counts(final_state)
    return (
        final_state,
        _charge3_sum(final_state),
        _family_sum(final_state),
        quarks,
        antiquarks,
    )


def _lightweight_concrete_process_enumeration(
    process: str,
    options: ProcessOptions,
) -> ProcessEnumeration:
    enumerator = ProcessEnumerator(options)
    request = enumerator.parse(process)
    all_outgoing = (*request.initial_state, *request.rest)
    if not enumerator._valid_process(
        all_outgoing,
        allow_charged_current=_request_allows_charged_current(request),
    ):
        return ProcessEnumeration(
            request=request,
            options=options,
            unique_processes=(),
            groups=(),
        )
    order = tuple(range(len(all_outgoing)))
    return ProcessEnumeration(
        request=request,
        options=options,
        unique_processes=(
            tuple(sorted(all_outgoing, key=lambda item: SORT_PARTICLES[item])),
        ),
        groups=(
            PhaseSpaceGroup(
                group_id=1,
                phase_space_order=order,
                records=(
                    SubprocessRecord(
                        process=tuple(all_outgoing),
                        color_order=order,
                        multichannel_partners=(0,),
                        identical_factor=1.0,
                    ),
                ),
            ),
        ),
    )


def write_legacy_process_file(
    process_string: str, path: str | Path, options: ProcessOptions | None = None
) -> ProcessEnumeration:
    enumerator = ProcessEnumerator(options)
    enumeration = enumerator.enumerate(process_string)
    enumerator.write_legacy_file(enumeration, path)
    return enumeration


def write_color_complete_legacy_process_file(
    process_string: str, path: str | Path, options: ProcessOptions | None = None
) -> ProcessEnumeration:
    enumerator = ProcessEnumerator(options)
    enumeration = enumerator.enumerate_color_complete(process_string)
    enumerator.write_legacy_file(enumeration, path)
    return enumeration


__all__ = [
    "ANTI_PARTICLE",
    "PDGS",
    "ParsedProcess",
    "PhaseSpaceGroup",
    "ProcessEnumeration",
    "ProcessEnumerator",
    "ProcessOptions",
    "ProcessSelectionRecord",
    "ProcessSelectionReport",
    "ProcessSetEntry",
    "ProcessSetEnumeration",
    "SubprocessRecord",
    "build_generic_process_selection_report",
    "canonical_process_key",
    "enumerate_generic_process_set",
    "enumerate_processes",
    "enumerate_process_set",
    "expand_process_variants",
    "split_process_set",
    "write_color_complete_legacy_process_file",
    "write_legacy_process_file",
]
