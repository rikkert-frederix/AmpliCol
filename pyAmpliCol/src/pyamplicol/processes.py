from __future__ import annotations

import itertools
import math
import re
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
        input_string = process_string.lower().replace("bar", "~").replace(" j", " 1j")
        parts = input_string.split(">")
        if len(parts) != 2:
            raise ValueError("invalid collision format; expected 'initial > final'")

        initial_state = parts[0].strip().split()
        if len(initial_state) != 2:
            raise ValueError("exactly two incoming particles are required")
        crossed_initial = tuple(
            ANTI_PARTICLE[p] if p != "p" else p for p in initial_state
        )

        final_state = parts[1].strip().split()
        jet_match = re.match(r"(\d+)j", final_state[-1]) if final_state else None
        jet_count = int(jet_match.group(1)) if jet_match else 0
        rest = list(final_state[:-1] if jet_match else final_state)
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
        unknown = sorted(set(particles).difference(PDGS).difference({"p"}))
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
        for perm in itertools.permutations(range(len(proc))):
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

    def _valid_process(self, proc: Sequence[str]) -> bool:
        nq = self._count_matching(proc, QUARKS)
        naq = self._count_matching(proc, ANTIQUARKS)
        if self.options.include_3qqbar:
            if nq > 3 or naq > 3:
                return False
        elif nq > 2 or naq > 2:
            return False
        if nq != naq:
            return False
        if sum(CHARGES3[x] for x in proc) != 0:
            return False
        if not self.options.include_cc:
            if sum(FAMILY[x] for x in proc) != 0:
                return False
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
        mandatory = [p for p in request.initial_state if p != "p"]
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
            if request.initial_state[i] != "p" and proc[i] != request.initial_state[i]:
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
            if part != "p" and part not in self.massless_qcd:
                raise ValueError("initial state should be a proton or massless QCD parton")

        qcd_rest = sum(1 for part in request.rest if part in self.massless_qcd)
        for part_index in range(request.jet_count + 2 + qcd_rest):
            new_processes: list[list[str]] = []
            for proc in processes:
                for particle in self.massless_qcd:
                    if (
                        not self.options.include_3qqbar
                        and part_index > 3
                        and particle not in GLUONS
                    ):
                        continue
                    if part_index > 5 and particle not in GLUONS:
                        continue
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
            if self._valid_process(proc) and self._compatible_unique_process(request, proc)
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

        iden = 2.0 if len(anti_quark_indices) == 3 else 1.0
        possible: list[tuple[OrderTuple, ProcessTuple]] = []
        if not singlet_perms:
            possible.append((perm, proc))
        elif len(anti_quark_indices) == 1:
            for singlets in singlet_perms:
                order: list[int] = []
                for i in range(len(perm)):
                    if i == anti_quark_indices[0]:
                        order.extend(perm[p] for p in (anti_quark_indices + singlets))
                    elif i not in singlet_indices:
                        order.append(perm[i])
                possible.append((tuple(order), proc))
        elif len(anti_quark_indices) == 2:
            for j in range(nsinglets + 1):
                for singlets in singlet_perms:
                    order_list: list[int] = []
                    for i in range(len(perm)):
                        if i == anti_quark_indices[0]:
                            order_list.extend(
                                perm[p] for p in ((anti_quark_indices[0],) + singlets[:j])
                            )
                        elif i == anti_quark_indices[1]:
                            order_list.extend(
                                perm[p] for p in ((anti_quark_indices[1],) + singlets[j:])
                            )
                        elif i not in singlet_indices:
                            order_list.append(perm[i])
                    possible.append((tuple(order_list), proc))
        elif len(anti_quark_indices) == 3:
            for j1 in range(nsinglets + 1):
                for j2 in range(nsinglets - j1 + 1):
                    for singlets in singlet_perms:
                        order_list = []
                        for i in range(len(perm)):
                            if i == anti_quark_indices[0]:
                                order_list.extend(
                                    perm[p]
                                    for p in ((anti_quark_indices[0],) + singlets[:j1])
                                )
                            elif i == anti_quark_indices[1]:
                                order_list.extend(
                                    perm[p]
                                    for p in (
                                        (anti_quark_indices[1],)
                                        + singlets[j1 : j1 + j2]
                                    )
                                )
                            elif i == anti_quark_indices[2]:
                                order_list.extend(
                                    perm[p]
                                    for p in ((anti_quark_indices[2],) + singlets[j1 + j2 :])
                                )
                            elif i not in singlet_indices:
                                order_list.append(perm[i])
                        possible.append((tuple(order_list), proc))

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
        if nq == 1:
            return math.factorial(ng)
        if nq == 2:
            return math.factorial(ng) * (ng + 1) * 2
        if nq == 3:
            return math.factorial(ng) * ((ng + 2) * (ng + 1) / 2) * 6
        raise ValueError(f"unknown number of quark lines: {nq}")

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


def write_legacy_process_file(
    process_string: str, path: str | Path, options: ProcessOptions | None = None
) -> ProcessEnumeration:
    enumerator = ProcessEnumerator(options)
    enumeration = enumerator.enumerate(process_string)
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
    "SubprocessRecord",
    "enumerate_processes",
    "write_legacy_process_file",
]
