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
    # Component values are independently cross-checked against pyAmpliCol's
    # compiled and recurrence evaluators, including 32-decimal-digit precision.
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
            1.3039052156690463e-14,
            5.849282219821445e-15,
            5.849282219821445e-15,
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
            1.7878308543514494e-12,
            8.295533767448123e-13,
            8.295533767448123e-13,
        ),
    ),
    Case(
        process="d d~ > u u~ u u~",
        pdgs=(1, -1, 2, -2, 2, -2),
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
            2.6038594934460394e-14,
            1.3395024066340542e-14,
            1.3395024066340542e-14,
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


def parse_indexed_matrix(
    output: str,
) -> tuple[dict[int, tuple[int, ...]], dict[tuple[str, int, int], float]]:
    permutations: dict[int, tuple[int, ...]] = {}
    entries: dict[tuple[str, int, int], float] = {}
    for line in output.splitlines():
        fields = line.split()
        if fields and fields[0] == "AMPICOL_COLOR_MATRIX_PERM":
            permutations[int(fields[1])] = tuple(int(value) for value in fields[2:])
        elif fields and fields[0] == "AMPICOL_COLOR_MATRIX_ENTRY":
            entries[(fields[1], int(fields[2]), int(fields[3]))] = float(fields[4])
    return permutations, entries


def _is_quark_label(label: int, pdgs: tuple[int, ...]) -> bool:
    pdg = pdgs[label - 1]
    return (label <= 2 and -6 <= pdg <= -1) or (label > 2 and 1 <= pdg <= 6)


def _is_antiquark_label(label: int, pdgs: tuple[int, ...]) -> bool:
    pdg = pdgs[label - 1]
    return (label <= 2 and 1 <= pdg <= 6) or (label > 2 and -6 <= pdg <= -1)


def endpoint_permutation(
    word: tuple[int, ...], pdgs: tuple[int, ...]
) -> tuple[int, ...]:
    endpoints: dict[int, int] = {}
    quark = 0
    for label in word:
        if _is_quark_label(label, pdgs):
            quark = label
        elif _is_antiquark_label(label, pdgs):
            if quark == 0:
                raise AssertionError(f"antiquark {label} has no open string in {word}")
            endpoints[quark] = label
            quark = 0
    if quark or len(endpoints) != 3:
        raise AssertionError(f"incomplete three-line flow word {word}")
    return tuple(endpoints[label] for label in sorted(endpoints))


def check_no_gluon_indexed_matrix(output: str, case: Case) -> None:
    permutations, entries = parse_indexed_matrix(output)
    endpoint_orders = {
        row: endpoint_permutation(word, case.pdgs)
        for row, word in permutations.items()
    }
    expected: dict[tuple[str, int, int], float] = {}
    for row, row_order in endpoint_orders.items():
        for col, col_order in endpoint_orders.items():
            if col < row:
                continue
            if row == col:
                expected[("lc", row, col)] = 27.0
                coefficient = 27.0
            else:
                relative = tuple(col_order.index(endpoint) for endpoint in row_order)
                fixed_points = sum(index == value for index, value in enumerate(relative))
                coefficient = -18.0 if fixed_points == 1 else 6.0
            expected[("nlc", row, col)] = coefficient
            expected[("full", row, col)] = coefficient
    if entries != expected:
        raise AssertionError(f"{case.process}: indexed color matrix is incorrect")


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
    if case.color_orders == 6:
        check_no_gluon_indexed_matrix(output, case)
    full_components = parse_components(output)
    for expected, actual in zip(case.components, full_components, strict=True):
        if not math.isclose(actual, expected, rel_tol=1.0e-12, abs_tol=1.0e-24):
            raise AssertionError(
                f"{case.process}: component differs: "
                f"expected={expected}, actual={actual}"
            )

    lc_output = run(
        [
            "./amplicol_color_probe",
            "1",
            str(group),
            str(integral),
            "lc",
            "processes.txt",
            str(momenta_file),
            *(str(value) for value in case.helicities),
        ],
        cwd=repository,
        env={**env, "AMPICOL_COLOR_PROBE_MATRIX": "1"},
    )
    expected_lc_matrix = {
        "lc": case.matrix_entries["lc"],
        "nlc": Counter(),
        "full": Counter(),
    }
    if parse_matrix(lc_output) != expected_lc_matrix:
        raise AssertionError(
            f"{case.process}: sparse LC matrix differs from the generic result"
        )
    lc_components = parse_components(lc_output)
    if not math.isclose(
        lc_components[0],
        full_components[0],
        rel_tol=1.0e-12,
        abs_tol=1.0e-24,
    ) or lc_components[1:] != (0.0, 0.0):
        raise AssertionError(
            f"{case.process}: sparse LC contraction differs: {lc_components}"
        )
    print(f"PASS {case.process}: {case.color_orders} all-flow color orders")


def check_singlet_rejected(
    repository: Path, temporary: Path, env: dict[str, str]
) -> None:
    process = "d d~ > u u~ s s~ a"
    pdgs = (1, -1, 2, -2, 3, -3, 22)
    run(
        [
            sys.executable,
            "process_list.py",
            "--serial",
            "--include_3qqbar",
            process,
        ],
        cwd=repository,
        env=env,
    )
    group, integral = selected_entry(repository / "processes.txt", pdgs)
    momenta_file = temporary / "momenta-singlet.txt"
    momenta_file.write_text(
        "\n".join(
            " ".join(f"{value:.17g}" for value in row)
            for row in CASES[1].momenta
        )
        + "\n",
        encoding="ascii",
    )
    completed = subprocess.run(
        [
            "./amplicol_color_probe",
            "1",
            str(group),
            str(integral),
            "full",
            "processes.txt",
            str(momenta_file),
            *(str(value) for value in CASES[1].helicities),
        ],
        cwd=repository,
        env=env,
        text=True,
        capture_output=True,
    )
    if completed.returncode == 0:
        raise AssertionError(f"{process}: unsupported direct probe unexpectedly ran")
    if "does not support colour singlets" not in completed.stdout:
        raise AssertionError(
            f"{process}: direct probe did not fail with its supported-scope message"
        )
    print(f"PASS {process}: unsupported singlet scope rejected before allocation")


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
            temporary = Path(raw) / ("long-momenta-path-" + "x" * 80)
            temporary.mkdir()
            for case in CASES:
                check_case(repository, temporary, case, env)
            check_singlet_rejected(repository, temporary, env)
    finally:
        if saved_processes is None:
            process_file.unlink(missing_ok=True)
        else:
            process_file.write_bytes(saved_processes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
