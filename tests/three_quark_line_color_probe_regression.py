#!/usr/bin/env python3
"""Exercise the direct all-flow oracle for three open quark lines."""

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

COMPONENTS_RE = re.compile(
    r"^AMPICOL_COLOR_PROBE_COMPONENTS\s+(\S+)\s+(\S+)\s+(\S+)$",
    re.MULTILINE,
)


@dataclass(frozen=True)
class Case:
    process: str
    pdgs: tuple[int, ...]
    helicities: tuple[int, ...]
    momenta: tuple[tuple[float, float, float, float], ...]
    color_orders: int
    matrix_entries: dict[str, Counter[float]]
    components: tuple[float, float, float]


CASES = (
    Case(
        process="d d~ > u u~ s s~",
        pdgs=(1, -1, 2, -2, 3, -3),
        helicities=(-1, 1, -1, 1, -1, 1),
        momenta=(
            (4671.200996478833, 0.0, 0.0, 4671.200996478833),
            (5452.62449645975, 0.0, 0.0, -5452.62449645975),
            (
                3848.769685279069,
                -1105.1428951212934,
                -1737.2006769281086,
                -3251.7412381317486,
            ),
            (
                3660.1488063565753,
                2228.296688374595,
                1630.5712325370819,
                2402.6278548445193,
            ),
            (
                1275.9167404758375,
                -470.2433927087442,
                -436.62760397349416,
                1102.8105076070963,
            ),
            (
                1338.9902608271016,
                -652.9104005445573,
                543.2570483645209,
                -1035.1206243007837,
            ),
        ),
        color_orders=6,
        matrix_entries={
            "lc": Counter({27.0: 6}),
            "nlc": Counter({-18.0: 9, 6.0: 6, 27.0: 6}),
            "full": Counter({-18.0: 9, 6.0: 6, 27.0: 6}),
        },
        components=(
            1.1063946294179874e-14,
            6.482595953811198e-15,
            6.482595953811198e-15,
        ),
    ),
    Case(
        process="d d~ > u u~ s s~ g",
        pdgs=(1, -1, 2, -2, 3, -3, 21),
        helicities=(-1, 1, -1, 1, -1, 1, -1),
        momenta=(
            (500.0, 0.0, 0.0, 500.0),
            (500.0, 0.0, 0.0, -500.0),
            (
                80.89693031749577,
                -35.363962752382264,
                25.38917257558109,
                -68.18426056773842,
            ),
            (
                251.36895499042657,
                -179.6909561218757,
                131.2696963946096,
                -116.90072125291724,
            ),
            (
                344.0807935954407,
                219.66991720940467,
                -259.57359242327067,
                52.519235628094506,
            ),
            (
                131.41282860492473,
                125.97195147403914,
                -1.2352831219401352,
                37.401511191104284,
            ),
            (
                192.24049249171213,
                -130.58694980918583,
                104.15000657501999,
                95.16423500145686,
            ),
        ),
        color_orders=18,
        matrix_entries={
            "lc": Counter({81.0: 18}),
            "nlc": Counter({16.0: 54, -48.0: 45, 72.0: 18}),
            "full": Counter({16.0: 54, -48.0: 45, 72.0: 18}),
        },
        components=(
            1.933598375456697e-12,
            1.5674722235444195e-12,
            1.5674722235444195e-12,
        ),
    ),
)


def run(command: list[str], *, cwd: Path, env: dict[str, str]) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        capture_output=True,
    )
    if completed.returncode:
        output = "\n".join(
            part for part in (completed.stdout, completed.stderr) if part
        )
        raise RuntimeError(f"{' '.join(command)} failed:\n{output}")
    return completed.stdout


def selected_entry(process_file: Path, pdgs: tuple[int, ...]) -> tuple[int, int]:
    lines = process_file.read_text(encoding="ascii").splitlines()
    n_external, n_unique = (int(value) for value in lines[0].split())
    cursor = 1 + n_unique

    def next_nonempty() -> list[str]:
        nonlocal cursor
        while cursor < len(lines) and not lines[cursor].strip():
            cursor += 1
        if cursor >= len(lines):
            raise RuntimeError("truncated generated process file")
        fields = lines[cursor].split()
        cursor += 1
        return fields

    n_groups = int(next_nonempty()[0])
    for _ in range(n_groups):
        header = next_nonempty()
        group, n_integrals = int(header[0]), int(header[1])
        for integral in range(1, n_integrals + 1):
            fields = next_nonempty()
            process_start = 1 + int(fields[0])
            row_pdgs = tuple(
                int(value)
                for value in fields[process_start : process_start + n_external]
            )
            if row_pdgs == pdgs:
                return group, integral
    raise RuntimeError(f"generated process file has no exact row for {pdgs}")


def parse_matrix(output: str) -> dict[str, Counter[float]]:
    result = {accuracy: Counter() for accuracy in ("lc", "nlc", "full")}
    for line in output.splitlines():
        fields = line.split()
        if len(fields) != 5 or fields[0] != "AMPICOL_COLOR_MATRIX_ENTRY":
            continue
        result[fields[1]][float(fields[4])] += 1
    return result


def parse_components(output: str) -> tuple[float, float, float]:
    matches = COMPONENTS_RE.findall(output)
    if len(matches) != 1:
        raise RuntimeError("probe did not emit exactly one component record")
    return tuple(float(value) for value in matches[0])  # type: ignore[return-value]


def check_case(
    repository: Path, temporary: Path, case: Case, env: dict[str, str]
) -> None:
    run(
        [
            sys.executable,
            "process_list.py",
            "--serial",
            "--include_3qqbar",
            case.process,
        ],
        cwd=repository,
        env=env,
    )
    group, integral = selected_entry(repository / "processes.txt", case.pdgs)
    momenta_file = temporary / f"momenta-{case.color_orders}.txt"
    momenta_file.write_text(
        "\n".join(" ".join(f"{value:.17g}" for value in row) for row in case.momenta)
        + "\n",
        encoding="ascii",
    )
    output = run(
        [
            "./amplicol_color_probe",
            "1",
            str(group),
            str(integral),
            "full",
            "processes.txt",
            str(momenta_file),
            *(str(value) for value in case.helicities),
        ],
        cwd=repository,
        env={**env, "AMPICOL_COLOR_PROBE_MATRIX": "1"},
    )
    expected_count = f"AMPICOL_COLOR_PROBE_COLOR_ORDERS {case.color_orders}"
    if expected_count not in output:
        raise AssertionError(f"{case.process}: missing {expected_count!r}")
    expected_amplitudes = f"AMPICOL_COLOR_PROBE_AMPLITUDES {case.color_orders}"
    if expected_amplitudes not in output:
        raise AssertionError(f"{case.process}: missing {expected_amplitudes!r}")
    if parse_matrix(output) != case.matrix_entries:
        raise AssertionError(f"{case.process}: unexpected LC/NLC/full color matrices")
    for expected, actual in zip(case.components, parse_components(output), strict=True):
        if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-24):
            raise AssertionError(
                f"{case.process}: component differs: "
                f"expected={expected}, actual={actual}"
            )
    print(f"PASS {case.process}: {case.color_orders} all-flow color orders")


def main() -> int:
    repository = Path(__file__).resolve().parents[1]
    env = dict(os.environ)
    process_file = repository / "processes.txt"
    saved_processes = process_file.read_bytes() if process_file.exists() else None
    try:
        run(
            ["make", "-B", "-j4", "amplicol_color_probe"],
            cwd=repository,
            env=env,
        )
        with tempfile.TemporaryDirectory(
            dir="/tmp", prefix="amplicol-three-line-color-"
        ) as raw:
            temporary = Path(raw)
            for case in CASES:
                check_case(repository, temporary, case, env)
    finally:
        if saved_processes is None:
            process_file.unlink(missing_ok=True)
        else:
            process_file.write_bytes(saved_processes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
