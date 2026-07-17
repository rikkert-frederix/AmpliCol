#!/usr/bin/env python3
"""Compare generated-library color rows with the direct legacy probe."""

from __future__ import annotations

import array
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


CASES = (
    ("d d~ > z g g", ("nlc", "full")),
    ("u d~ > w+ g g", ("nlc", "full")),
    ("g g > g g g", ("nlc",)),
)
COMPONENTS_RE = re.compile(
    r"^AMPICOL_COLOR_PROBE_COMPONENTS\s+(\S+)\s+(\S+)\s+(\S+)$",
    re.MULTILINE,
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
        output = "\n".join(part for part in (completed.stdout, completed.stderr) if part)
        raise RuntimeError(f"{' '.join(command)} failed:\n{output}")
    return completed.stdout


def probe_components(output: str) -> tuple[float, float, float]:
    matches = COMPONENTS_RE.findall(output)
    if len(matches) != 1:
        raise RuntimeError("probe did not emit exactly one component record")
    return tuple(float(value) for value in matches[0])  # type: ignore[return-value]


def write_library_momenta(repository: Path, destination: Path) -> None:
    external_count = int((repository / "processes.txt").read_text().split()[0])
    values = array.array("d")
    byte_count = 4 * external_count * values.itemsize
    payload = (repository / "Library" / "amp1_1_lib.data").read_bytes()
    values.frombytes(payload[:byte_count])
    if values.itemsize != 8 or len(values) != 4 * external_count:
        raise RuntimeError("unexpected generated-library momentum encoding")
    rows = (
        " ".join(f"{component:.17g}" for component in values[4 * index : 4 * index + 4])
        for index in range(external_count)
    )
    destination.write_text("\n".join(rows) + "\n", encoding="ascii")


def assert_components_close(
    direct: tuple[float, float, float],
    library: tuple[float, float, float],
    *,
    context: str,
) -> None:
    for index, (expected, actual) in enumerate(zip(direct, library, strict=True), 1):
        if not math.isclose(actual, expected, rel_tol=1.0e-8, abs_tol=1.0e-15):
            raise AssertionError(
                f"{context}: component {index} differs: direct={expected}, library={actual}"
            )


def main() -> int:
    repository = Path(__file__).resolve().parents[1]
    env = dict(os.environ)
    process_file = repository / "processes.txt"
    saved_processes = process_file.read_bytes() if process_file.exists() else None

    try:
        run(
            ["make", "-B", "-j4", "amplicol_generate", "amplicol_color_probe"],
            cwd=repository,
            env=env,
        )
        with tempfile.TemporaryDirectory(dir="/tmp", prefix="amplicol-color-map-") as temporary:
            momenta_file = Path(temporary) / "momenta.txt"
            for process, accuracies in CASES:
                run(["make", "cleanlib"], cwd=repository, env=env)
                run(
                    [sys.executable, "process_list.py", "--serial", process],
                    cwd=repository,
                    env=env,
                )
                run(
                    [
                        "./amplicol_generate",
                        "--library=create-raw",
                        "--process=processes.txt",
                    ],
                    cwd=repository,
                    env=env,
                )
                run(
                    ["make", "-j4", "amplicol_color_library_probe"],
                    cwd=repository,
                    env=env,
                )
                write_library_momenta(repository, momenta_file)
                for accuracy in accuracies:
                    direct = probe_components(
                        run(
                            [
                                "./amplicol_color_probe",
                                "1",
                                "1",
                                "1",
                                accuracy,
                                "processes.txt",
                                str(momenta_file),
                            ],
                            cwd=repository,
                            env=env,
                        )
                    )
                    library = probe_components(
                        run(
                            [
                                "./amplicol_color_library_probe",
                                "1",
                                "1",
                                "1",
                                accuracy,
                                str(momenta_file),
                            ],
                            cwd=repository,
                            env=env,
                        )
                    )
                    assert_components_close(
                        direct,
                        library,
                        context=f"{process} ({accuracy})",
                    )
                    print(f"PASS {process} {accuracy}")
    finally:
        if saved_processes is None:
            process_file.unlink(missing_ok=True)
        else:
            process_file.write_bytes(saved_processes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
