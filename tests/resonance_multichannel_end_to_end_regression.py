#!/usr/bin/env python3
"""Fixed-seed regression for adaptive resonance multichannel integration."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path


PROCESS = "u u~ > e+ e- ve ve~"
MADGRAPH_REFERENCE_PB = 0.05525988
SEED = 24681357
NEVENTS = 100
ITMAX = 5
RESULT_RE = re.compile(
    r"Integral\s+\(production/accumulated\):\s*"
    r"([-+0-9.Ee]+)\s+\+/-\s*([-+0-9.Ee]+)"
)


def run_checked(command: list[str], cwd: Path) -> None:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        timeout=900,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


def parse_groups(process_file: Path) -> tuple[list[str], int, list[dict[str, object]]]:
    lines = process_file.read_text(encoding="ascii").splitlines()
    ngroups_index = next(
        index
        for index in range(4, len(lines))
        if lines[index].strip().isdigit()
    )
    ngroups = int(lines[ngroups_index].strip())
    groups: list[dict[str, object]] = []
    cursor = ngroups_index + 1
    for _ in range(ngroups):
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        header = lines[cursor]
        header_fields = header.split()
        group_id = int(header_fields[0])
        nprocesses = int(header_fields[1])
        nresonances = int(header_fields[-1])
        cursor += 1
        resonances = lines[cursor : cursor + nresonances]
        cursor += nresonances
        rows = lines[cursor : cursor + nprocesses]
        cursor += nprocesses
        groups.append(
            {
                "id": group_id,
                "header": header,
                "resonances": resonances,
                "rows": rows,
            }
        )
    return lines[:ngroups_index], ngroups, groups


def write_one_ww_control(full_process_file: Path, control_file: Path) -> None:
    prefix, ngroups, groups = parse_groups(full_process_file)
    if ngroups != 37:
        raise AssertionError(f"expected 37 resonance maps, found {ngroups}")

    wanted_resonances = {(24, (3, 5)), (-24, (4, 6))}
    matches = []
    for group in groups:
        parsed = set()
        for line in group["resonances"]:
            fields = [int(value) for value in str(line).split()]
            parsed.add((fields[0], tuple(fields[2:])))
        if parsed == wanted_resonances:
            matches.append(group)
    if len(matches) != 1:
        raise AssertionError(
            f"expected one double-W resonance map, found {len(matches)}"
        )
    group = matches[0]
    rows = list(group["rows"])
    if len(rows) != 1:
        raise AssertionError("one-WW control expects one full-amplitude row")

    header_fields = str(group["header"]).split()
    header_fields[0] = "1"
    header_fields[2] = "1"

    row_fields = rows[0].split()
    number_of_channels = int(row_fields[0])
    channel_ids = [int(value) for value in row_fields[1 : 1 + number_of_channels]]
    offset = 1 + number_of_channels
    process = row_fields[offset : offset + 6]
    offset += 6
    colour_order = row_fields[offset : offset + 6]
    offset += 6
    identical_factor = row_fields[offset]
    offset += 1
    phase_permutation = row_fields[offset : offset + 6]
    offset += 6
    permutations = row_fields[offset:]
    if len(permutations) != 6 * number_of_channels:
        raise AssertionError("malformed multichannel permutation list")
    self_position = channel_ids.index(int(group["id"]))
    self_permutation = permutations[6 * self_position : 6 * (self_position + 1)]
    control_row = [
        "1",
        "1",
        *process,
        *colour_order,
        identical_factor,
        *phase_permutation,
        *self_permutation,
    ]

    output = [
        *prefix,
        "1",
        "",
        " ".join(header_fields),
        *[str(line) for line in group["resonances"]],
        " ".join(control_row),
        "",
        "",
        "# end",
    ]
    control_file.write_text("\n".join(output) + "\n", encoding="ascii")


def run_generator(
    generator: Path,
    process_file: Path,
    input_card: Path,
    workdir: Path,
    tag: str,
    combine_subprocesses: bool,
) -> tuple[Path, Path]:
    (workdir / "Outputs").mkdir(parents=True)
    command = [
        str(generator),
        # The generator's legacy process-filename buffer is only 80 characters.
        # The process file lives in workdir, so a relative name also exercises
        # the normal command-line interface without depending on /tmp length.
        f"--process={process_file.name}",
        f"--input={input_card}",
        f"--nevents={NEVENTS}",
        f"--itmax={ITMAX}",
        f"--seed={SEED}",
        "--phasespace=1",
        f"--tag={tag}",
        "--timing=none",
    ]
    if combine_subprocesses:
        command.append("--combine_subprocesses")
    run_checked(command, workdir)
    return (
        workdir / "Outputs" / f"{tag}_log_file.txt",
        workdir / "Outputs" / f"{tag}_events.lhe",
    )


def production_results(log_file: Path) -> tuple[str, list[tuple[float, float]]]:
    log = log_file.read_text(encoding="ascii")
    results = [
        (float(match.group(1)), float(match.group(2)))
        for match in RESULT_RE.finditer(log)
    ]
    if not results:
        raise AssertionError(f"no production result found in {log_file}")
    if "warm-up/current" not in log or "production/accumulated" not in log:
        raise AssertionError("warm-up/production log labels are missing")
    spam_markers = (
        "LUP decomposition failure",
        "Warning: gram4",
        "Warning: variable not between varmin and varmax",
    )
    for marker in spam_markers:
        if marker in log:
            raise AssertionError(f"low-level inverse-map warning leaked into log: {marker}")
    if log.count("unexpected failure of a forward phase-space map") > 1:
        raise AssertionError("forward-map diagnostic was not rate-limited")
    return log, results


def assert_three_sigma(
    first: tuple[float, float],
    second: tuple[float, float],
    label: str,
) -> None:
    difference = abs(first[0] - second[0])
    combined_uncertainty = math.hypot(first[1], second[1])
    if difference > 3.0 * combined_uncertainty:
        raise AssertionError(
            f"{label}: {first[0]:.8g} +/- {first[1]:.3g} versus "
            f"{second[0]:.8g} +/- {second[1]:.3g}"
        )


def check_colour_expansion(
    reweighter: Path,
    event_file: Path,
    input_card: Path,
    workdir: Path,
) -> None:
    run_checked(
        [str(reweighter), str(event_file), f"--input={input_card}"],
        workdir,
    )
    reweighted = Path(f"{event_file}.rwgt")
    colour_rows = []
    for line in reweighted.read_text(encoding="ascii").splitlines():
        if line.startswith("#color_expansion"):
            colour_rows.append([float(value) for value in line.split()[1:]])
    if len(colour_rows) != NEVENTS:
        raise AssertionError(
            f"expected {NEVENTS} colour-expansion rows, found {len(colour_rows)}"
        )
    for event_number, values in enumerate(colour_rows, start=1):
        if values[0] != values[1] or values[1] != values[2]:
            raise AssertionError(
                f"LC/NLC/full-colour mismatch in event {event_number}: {values}"
            )


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

    with tempfile.TemporaryDirectory(prefix="amplicol_resonance_multichannel_") as tmp:
        root = Path(tmp)
        full_dir = root / "full"
        control_dir = root / "control"
        full_dir.mkdir()
        control_dir.mkdir()
        full_process_file = full_dir / "processes.txt"
        run_checked(
            [sys.executable, str(process_list), "--serial", PROCESS],
            full_dir,
        )
        if not full_process_file.exists():
            raise AssertionError("process-list generator did not create processes.txt")
        _, ngroups, _ = parse_groups(full_process_file)
        if ngroups != 37:
            raise AssertionError(f"expected 37 maps, found {ngroups}")

        control_process_file = control_dir / "processes.txt"
        write_one_ww_control(full_process_file, control_process_file)

        full_log, _ = run_generator(
            generator,
            full_process_file,
            input_card,
            full_dir,
            "full",
            combine_subprocesses=False,
        )
        control_log, control_events = run_generator(
            generator,
            control_process_file,
            input_card,
            control_dir,
            "control",
            combine_subprocesses=True,
        )

        _, full_history = production_results(full_log)
        _, control_history = production_results(control_log)
        if len(full_history) < 2:
            raise AssertionError("37-map run did not accumulate two production iterations")
        full_result = full_history[-1]
        control_result = control_history[-1]
        reference_result = (MADGRAPH_REFERENCE_PB, 0.0)
        assert_three_sigma(full_result, control_result, "37-map/control mismatch")
        assert_three_sigma(full_result, reference_result, "37-map/MadGraph mismatch")
        assert_three_sigma(control_result, reference_result, "control/MadGraph mismatch")

        first_production = full_history[0]
        relative_shift = abs(full_result[0] - first_production[0]) / abs(
            first_production[0]
        )
        if relative_shift >= 0.08:
            raise AssertionError(
                f"production accumulation drifted by {100.0 * relative_shift:.2f}%"
            )
        assert_three_sigma(
            full_result,
            first_production,
            "production iterations are statistically incompatible",
        )

        check_colour_expansion(
            reweighter,
            control_events,
            input_card,
            control_dir,
        )

    print(
        "resonance multichannel end-to-end regression passed: "
        f"37-map={full_result[0]:.8f} +/- {full_result[1]:.8f} pb, "
        f"control={control_result[0]:.8f} +/- {control_result[1]:.8f} pb"
    )


if __name__ == "__main__":
    main()
