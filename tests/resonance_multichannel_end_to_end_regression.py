#!/usr/bin/env python3
"""Fixed-seed 4FS acceptance test for compact integration families."""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path


PROCESS = "p p > e+ e- ve ve~ b b~"
SEED = 24681357
MAX_CONSTRUCTION_SECONDS = 10.0
MAX_STARTUP_SECONDS = 15.0
MAX_PROCESS_FILE_BYTES = 5_000_000
BOTTOM_MASS = 4.7
AUXILIARY_PDGS = {-21, -23, 26, -26, 99, 125, 126, 127}
RESULT_RE = re.compile(
    r"Integral\s+\(production/accumulated\):\s*"
    r"([-+0-9.Ee]+)\s+\+/-\s*([-+0-9.Ee]+)"
)


def run_checked(command: list[str], cwd: Path, timeout: float = 900) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout


def section_header(lines: list[str], name: str) -> tuple[int, int]:
    for index, line in enumerate(lines):
        fields = line.split()
        if fields and fields[0] == name:
            if len(fields) != 2:
                raise AssertionError(f"malformed {name} header: {line}")
            return index, int(fields[1])
    raise AssertionError(f"missing {name} catalogue")


def parse_v6(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    if len(raw) >= MAX_PROCESS_FILE_BYTES:
        raise AssertionError(
            f"processes.txt is {len(raw)} bytes, expected less than "
            f"{MAX_PROCESS_FILE_BYTES}"
        )
    lines = raw.decode("ascii").splitlines()
    header = [int(value) for value in lines[0].split()]
    if len(header) != 3 or header[2] != 6:
        raise AssertionError(f"expected process-file version 6, found {header}")
    if "flavour_scheme=4" not in lines[2]:
        raise AssertionError("target process file is not in the four-flavour scheme")

    p_index, permutation_count = section_header(lines, "PERMUTATIONS")
    m_index, map_count = section_header(lines, "PHASE_MAPS")
    s_index, partner_set_count = section_header(lines, "PARTNER_SETS")
    f_index, family_count = section_header(lines, "INTEGRATION_FAMILIES")
    if family_count > 20:
        raise AssertionError(f"expected at most 20 families, found {family_count}")
    if partner_set_count != family_count:
        raise AssertionError("partner sets and integration families are not one-to-one")

    permutations = {}
    for line in lines[p_index + 1 : m_index]:
        fields = line.split()
        if fields and fields[0] == "P":
            permutations[int(fields[1])] = tuple(map(int, fields[2:]))
    if len(permutations) != permutation_count:
        raise AssertionError("permutation catalogue count mismatch")

    maps: dict[int, tuple[tuple[int, ...], list[tuple[int, ...]]]] = {}
    current_map = None
    for line in lines[m_index + 1 : s_index]:
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "M":
            current_map = int(fields[1])
            maps[current_map] = (tuple(map(int, fields[3:])), [])
            if int(fields[2]) < 0:
                raise AssertionError("negative topology-node count")
        elif fields[0] == "N":
            if current_map is None or len(fields) != 7:
                raise AssertionError(f"orphan or malformed node: {line}")
            node = tuple(map(int, fields[1:]))
            pdg, kind = node[:2]
            if pdg in AUXILIARY_PDGS:
                raise AssertionError(f"auxiliary PDG {pdg} leaked into a map")
            if (pdg == 0) != (kind == 4):
                raise AssertionError(f"invalid flat-contact encoding: {line}")
            maps[current_map][1].append(node)
    if len(maps) != map_count:
        raise AssertionError("phase-map catalogue count mismatch")

    partner_sets = {}
    for line in lines[s_index + 1 : f_index]:
        fields = line.split()
        if not fields:
            continue
        if fields[0] != "S":
            raise AssertionError(f"unexpected partner-set record: {line}")
        set_id, pair_count = map(int, fields[1:3])
        values = list(map(int, fields[3:]))
        if len(values) != 2 * pair_count:
            raise AssertionError(f"malformed partner set {set_id}")
        pairs = tuple(zip(values[0::2], values[1::2]))
        if pairs != tuple(sorted(set(pairs))):
            raise AssertionError(f"partner set {set_id} is not canonical")
        if any(map_id not in maps or permutation_id not in permutations
               for map_id, permutation_id in pairs):
            raise AssertionError(f"partner set {set_id} has a dangling reference")
        partner_sets[set_id] = pairs
    if len(partner_sets) != partner_set_count:
        raise AssertionError("partner-set catalogue count mismatch")

    families = {}
    cursor = f_index + 1
    while cursor < len(lines):
        fields = lines[cursor].split()
        cursor += 1
        if not fields:
            continue
        if fields[0] == "END_PROCESSES":
            break
        if fields[0] != "F" or len(fields) != 4:
            raise AssertionError(f"unexpected family record: {' '.join(fields)}")
        family_id, set_id, row_count = map(int, fields[1:])
        rows = []
        while len(rows) < row_count:
            fields = lines[cursor].split()
            cursor += 1
            if fields:
                if fields[0] != "C":
                    raise AssertionError("family coefficient row is not tagged C")
                rows.append(tuple(fields[1:]))
        if set_id not in partner_sets or not rows:
            raise AssertionError(f"invalid integration family {family_id}")
        if len(set(rows)) != len(rows):
            raise AssertionError(f"family {family_id} repeats a coefficient row")
        families[family_id] = (set_id, tuple(rows))
    if len(families) != family_count:
        raise AssertionError("integration-family catalogue count mismatch")

    coefficient_rows = [
        row for _, rows in families.values() for row in rows
    ]
    for row in coefficient_rows:
        if len(row) != 2 * header[0] + 1:
            raise AssertionError("malformed integration-family coefficient row")
        incoming = tuple(map(int, row[:2]))
        if any(abs(pdg) == 5 for pdg in incoming):
            raise AssertionError(
                f"four-flavour catalogue contains initial bottom {incoming}"
            )

    topologies = [nodes for _, nodes in maps.values()]

    def has_nodes(*requirements: tuple[int, int]) -> bool:
        return any(
            all(any(node[0] == pdg and node[3] == mask for node in nodes)
                for pdg, mask in requirements)
            for nodes in topologies
        )

    # External labels are fixed by PROCESS: e+, e-, ve, ve~, b, b~ occupy
    # zero-based labels 2--7.  These signatures distinguish the required
    # physical histories rather than merely checking that a PDG occurs.
    expected_histories = {
        "top pair": has_nodes((6, 84), (-6, 168)),
        "double W": has_nodes((24, 20), (-24, 40)),
        "gamma to bottoms": has_nodes((22, 192)),
        "Z to bottoms": has_nodes((23, 192)),
        "Higgs to bottoms": has_nodes((25, 192)),
        "charged-lepton neutral current": has_nodes((22, 12))
        or has_nodes((23, 12)),
        "neutrino neutral current": has_nodes((23, 48)),
        "bottom radiation": any(
            abs(node[0]) == 5 and node[3].bit_count() > 2
            for nodes in topologies for node in nodes
        ),
        "lepton radiation": any(
            11 <= abs(node[0]) <= 12 and node[3].bit_count() > 2
            for nodes in topologies for node in nodes
        ),
        "quartic contact": any(
            node[0] == 0 and node[1] == 4
            for nodes in topologies for node in nodes
        ),
    }
    missing = [name for name, present in expected_histories.items() if not present]
    if missing:
        raise AssertionError(
            "target catalogue is missing physical histories: "
            + ", ".join(missing)
        )

    return {
        "bytes": len(raw),
        "prefix": lines[:p_index],
        "permutations": permutations,
        "maps": maps,
        "partner_sets": partner_sets,
        "families": families,
    }


def write_single_map_control(catalogues: dict[str, object], cwd: Path) -> None:
    """Write a cheap full-run control while retaining the fixed target state."""

    permutations = catalogues["permutations"]
    maps = catalogues["maps"]
    families = catalogues["families"]
    permutation = permutations[min(permutations)]
    order, nodes = maps[min(maps)]
    _, rows = families[min(families)]
    row = rows[0]
    lines = [
        *catalogues["prefix"],
        "PERMUTATIONS 1",
        "P 1 " + " ".join(map(str, permutation)),
        "",
        "PHASE_MAPS 1",
        f"M 1 {len(nodes)} " + " ".join(map(str, order)),
        *("N " + " ".join(map(str, node)) for node in nodes),
        "",
        "PARTNER_SETS 1",
        "S 1 1 1 1",
        "",
        "INTEGRATION_FAMILIES 1",
        "F 1 1 1",
        "C " + " ".join(row),
        "",
        "END_PROCESSES",
    ]
    (cwd / "processes.txt").write_text(
        "\n".join(lines) + "\n", encoding="ascii"
    )


def generate_process_file(process_list: Path, cwd: Path) -> tuple[Path, float]:
    started = time.perf_counter()
    run_checked([
        sys.executable,
        str(process_list),
        "--serial",
        "--flavour_scheme=4",
        PROCESS,
    ], cwd)
    elapsed = time.perf_counter() - started
    if elapsed >= MAX_CONSTRUCTION_SECONDS:
        raise AssertionError(
            f"process construction took {elapsed:.3f}s; expected <"
            f"{MAX_CONSTRUCTION_SECONDS:.0f}s"
        )
    path = cwd / "processes.txt"
    if not path.exists():
        raise AssertionError("process-list generator did not create processes.txt")
    return path, elapsed


def generator_command(
    generator: Path,
    input_card: Path,
    *,
    tag: str,
    itmax: int,
    nevents: int,
    combine_subprocesses: bool,
) -> list[str]:
    command = [
        str(generator),
        "--process=processes.txt",
        f"--input={input_card}",
        f"--nevents={nevents}",
        f"--itmax={itmax}",
        f"--seed={SEED}",
        "--phasespace=1",
        f"--tag={tag}",
        "--timing=none",
    ]
    if combine_subprocesses:
        command.append("--combine_subprocesses")
    return command


def check_startup(
    generator: Path,
    input_card: Path,
    cwd: Path,
    family_count: int,
) -> float:
    (cwd / "Outputs").mkdir(exist_ok=True)
    log_path = cwd / "Outputs" / "startup_log_file.txt"
    environment = os.environ.copy()
    environment["GFORTRAN_UNBUFFERED_ALL"] = "1"
    started = time.perf_counter()
    process = subprocess.Popen(
        generator_command(
            generator, input_card, tag="startup", itmax=2, nevents=1,
            combine_subprocesses=True,
        ),
        cwd=cwd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=environment,
    )
    try:
        while time.perf_counter() - started < MAX_STARTUP_SECONDS:
            if process.poll() is not None:
                raise AssertionError(
                    f"generator exited during startup with code {process.returncode}"
                )
            if log_path.exists():
                log = log_path.read_text(encoding="ascii", errors="replace")
                if "Start phase-space integration" in log:
                    elapsed = time.perf_counter() - started
                    expected_initialisations = family_count + 1
                    if log.count("Initialising amplitude for:") != \
                            expected_initialisations:
                        raise AssertionError(
                            "expected one amplitude initialisation per family "
                            "plus the unique-process target"
                        )
                    mass_rows = [
                        [float(value) for value in line.split()[1:]]
                        for line in log.splitlines()
                        if line.strip().startswith("masses:")
                    ]
                    if len(mass_rows) != family_count or any(
                        sum(abs(value - BOTTOM_MASS) < 1.0e-12
                            for value in masses) != 2
                        for masses in mass_rows
                    ):
                        raise AssertionError(
                            "integration families lost their two massive-bottom "
                            "thresholds"
                        )
                    return elapsed
            time.sleep(0.05)
        raise AssertionError(
            f"generator startup exceeded {MAX_STARTUP_SECONDS:.0f}s"
        )
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)


def validate_lhe_events(event_file: Path, require_bottoms: bool) -> list[float]:
    lines = event_file.read_text(encoding="ascii").splitlines()
    event_weights = []
    cursor = 0
    while cursor < len(lines):
        if lines[cursor].strip() != "<event>":
            cursor += 1
            continue
        cursor += 1
        header = lines[cursor].split()
        if len(header) < 3:
            raise AssertionError("malformed LHE event header")
        n_particles = int(header[0])
        event_weight = float(header[2])
        if not math.isfinite(event_weight) or event_weight <= 0.0:
            raise AssertionError(f"non-finite LHE event weight {event_weight}")
        event_weights.append(event_weight)
        cursor += 1
        bottom_masses = []
        for _ in range(n_particles):
            fields = lines[cursor].split()
            cursor += 1
            if len(fields) < 13:
                raise AssertionError("malformed LHE particle row")
            pdg, status = int(fields[0]), int(fields[1])
            momentum = [float(value) for value in fields[6:13]]
            if any(not math.isfinite(value) for value in momentum):
                raise AssertionError("non-finite momentum or weight in LHE event")
            if abs(pdg) == 5 and status == 1:
                px, py, pz, energy, stored_mass = momentum[:5]
                shell_mass = math.sqrt(max(
                    0.0, energy * energy - px * px - py * py - pz * pz
                ))
                if abs(stored_mass - BOTTOM_MASS) > 1.0e-10 or \
                        abs(shell_mass - BOTTOM_MASS) > 2.0e-8:
                    raise AssertionError(
                        "massive bottom left its finite phase-space threshold: "
                        f"stored={stored_mass}, shell={shell_mass}"
                    )
                bottom_masses.append(stored_mass)
        if require_bottoms and len(bottom_masses) != 2:
            raise AssertionError(
                f"expected two massive final bottoms, found {len(bottom_masses)}"
            )
    if not event_weights:
        raise AssertionError("fixed-seed run did not produce an event")
    return event_weights


def check_fixed_seed_run(
    generator: Path,
    input_card: Path,
    cwd: Path,
    *,
    tag: str,
    combine_subprocesses: bool,
    require_bottoms: bool = True,
    timeout: float = 900,
) -> tuple[tuple[float, float], Path, list[float]]:
    run_checked(
        generator_command(
            generator, input_card, tag=tag, itmax=1, nevents=1,
            combine_subprocesses=combine_subprocesses,
        ),
        cwd,
        timeout=timeout,
    )
    log_path = cwd / "Outputs" / f"{tag}_log_file.txt"
    log = log_path.read_text(encoding="ascii")
    matches = RESULT_RE.findall(log)
    if not matches:
        raise AssertionError("fixed-seed run has no production result")
    result, uncertainty = map(float, matches[-1])
    if not math.isfinite(result) or not math.isfinite(uncertainty) or result <= 0:
        raise AssertionError(
            f"fixed-seed run returned {result} +/- {uncertainty}"
        )
    spam = (
        "LUP decomposition failure",
        "Warning: gram4",
        "Warning: variable not between varmin and varmax",
    )
    if any(marker in log for marker in spam):
        raise AssertionError("invalid partner inversion leaked warning spam")
    if log.count("Unexpected failure of selected family phase map") > 1:
        raise AssertionError("selected-map failure diagnostic was not rate-limited")
    event_file = cwd / "Outputs" / f"{tag}_events.lhe"
    event_weights = validate_lhe_events(event_file, require_bottoms)
    return (result, uncertainty), event_file, event_weights


def assert_three_sigma(
    first: tuple[float, float],
    second: tuple[float, float],
) -> None:
    difference = abs(first[0] - second[0])
    uncertainty = math.hypot(first[1], second[1])
    if uncertainty <= 0.0 or difference > 3.0 * uncertainty:
        raise AssertionError(
            "separate/combined subprocess estimates disagree: "
            f"{first[0]:.8g} +/- {first[1]:.3g} versus "
            f"{second[0]:.8g} +/- {second[1]:.3g}"
        )


def check_colour_expansion(
    reweighter: Path,
    event_file: Path,
    input_card: Path,
    cwd: Path,
    event_weights: list[float],
) -> tuple[tuple[float, float, float], ...]:
    run_checked(
        [str(reweighter), str(event_file), f"--input={input_card}"],
        cwd,
        timeout=300,
    )
    reweighted = Path(f"{event_file}.rwgt")
    rows = [
        tuple(float(value) for value in line.split()[1:4])
        for line in reweighted.read_text(encoding="ascii").splitlines()
        if line.startswith("#color_expansion")
    ]
    if len(rows) != len(event_weights):
        raise AssertionError(
            "colour reweighting did not preserve the fixed event sample"
        )
    for event_weight, values in zip(event_weights, rows):
        if any(not math.isfinite(value) or value <= 0.0 for value in values):
            raise AssertionError(f"non-finite LC/NLC/full-colour values {values}")
        scale = max(1.0, abs(event_weight), *(abs(value) for value in values))
        if abs(values[0] - event_weight) > 5.0e-8 * scale:
            raise AssertionError(
                "stored event value disagrees with its LC reweighting: "
                f"{event_weight}, {values[0]}"
            )
        if abs(values[1] - values[2]) > 5.0e-12 * scale:
            raise AssertionError(
                f"NLC and full-colour event values disagree: {values}"
            )
    return tuple(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--reweighter", type=Path, required=True)
    parser.add_argument("--process-list", type=Path, required=True)
    parser.add_argument("--input-card", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generator = args.generator.resolve()
    reweighter = args.reweighter.resolve()
    process_list = args.process_list.resolve()
    input_card = args.input_card.resolve()
    source_root = process_list.parent

    multichannel_source = (source_root / "multichannel.f03").read_text(
        encoding="ascii"
    )
    family_routine = multichannel_source.split(
        "subroutine compute_family_multichannel_weight", 1
    )[1].split("end subroutine compute_family_multichannel_weight", 1)[0]
    if "compute_x_from_momenta" in family_routine:
        raise AssertionError("family density loop uses the legacy inverse path")
    if "compute_x_from_cache" not in family_routine:
        raise AssertionError("family density loop does not use the batch cache")

    with tempfile.TemporaryDirectory(prefix="amplicol_v6_family_") as tmp:
        root = Path(tmp)
        first = root / "first"
        second = root / "second"
        first.mkdir()
        second.mkdir()
        process_file, construction_seconds = generate_process_file(
            process_list, first
        )
        catalogues = parse_v6(process_file)
        second_file, _ = generate_process_file(process_list, second)
        if process_file.read_bytes() != second_file.read_bytes():
            raise AssertionError("version-6 process construction is not deterministic")
        (second / "Outputs").mkdir()

        startup_seconds = check_startup(
            generator,
            input_card,
            first,
            len(catalogues["families"]),
        )
        combined_result, combined_events, combined_weights = \
            check_fixed_seed_run(
                generator,
                input_card,
                first,
                tag="combined",
                combine_subprocesses=True,
            )
        separate_result, separate_events, separate_weights = \
            check_fixed_seed_run(
                generator,
                input_card,
                second,
                tag="separate",
                combine_subprocesses=False,
            )
        assert_three_sigma(combined_result, separate_result)
        check_colour_expansion(
            reweighter, combined_events, input_card, first, combined_weights
        )
        check_colour_expansion(
            reweighter, separate_events, input_card, second, separate_weights
        )

        control = root / "control"
        control.mkdir()
        (control / "Outputs").mkdir()
        write_single_map_control(catalogues, control)
        control_result, _, _ = check_fixed_seed_run(
            generator,
            input_card,
            control,
            tag="control",
            combine_subprocesses=True,
        )

    print(
        "v6 family end-to-end regression passed: "
        f"{len(catalogues['families'])} families, "
        f"{len(catalogues['maps'])} maps, "
        f"{catalogues['bytes']} bytes, "
        f"construction={construction_seconds:.3f}s, "
        f"startup={startup_seconds:.3f}s, "
        f"combined={combined_result[0]:.8g} +/- {combined_result[1]:.3g}, "
        f"separate={separate_result[0]:.8g} +/- {separate_result[1]:.3g}, "
        f"control={control_result[0]:.8g} +/- {control_result[1]:.3g}"
    )


if __name__ == "__main__":
    main()
