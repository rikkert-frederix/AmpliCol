#!/usr/bin/env python3
"""Check that the generated-library probe returns one physical helicity."""

from __future__ import annotations

import math
import os
import re
import subprocess
import sys
from pathlib import Path

PROCESS = "d d~ > z"
HELICITIES = (-1, 1, -1)
MOMENTA = (
    (45.594, 0.0, 0.0, 45.594),
    (45.594, 0.0, 0.0, -45.594),
    (91.188, 0.0, 0.0, 0.0),
)


def run(command: list[str], *, cwd: Path) -> str:
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=os.environ,
        text=True,
        capture_output=True,
    )
    if completed.returncode:
        output = "\n".join(
            part for part in (completed.stdout, completed.stderr) if part
        )
        raise RuntimeError(f"{' '.join(command)} failed:\n{output}")
    return completed.stdout


def scalar(output: str, label: str) -> float:
    matches = re.findall(
        rf"^{re.escape(label)}(?:\s+\d+\s+\d+)?\s+([+\-0-9.Ee]+)$",
        output,
        re.MULTILINE,
    )
    if len(matches) != 1:
        raise AssertionError(f"expected exactly one {label} record")
    return float(matches[0])


def main() -> int:
    repository = Path(__file__).resolve().parents[1]
    process_file = repository / "processes.txt"
    momentum_file = repository / "Utilities" / "ME_checks" / "momenta_1_1.txt"
    saved_process = process_file.read_bytes() if process_file.exists() else None
    saved_momentum = momentum_file.read_bytes() if momentum_file.exists() else None
    momentum_argument = os.fspath(momentum_file.relative_to(repository))
    try:
        run([sys.executable, "process_list.py", "--serial", PROCESS], cwd=repository)
        momentum_file.write_text(
            "\n".join(
                " ".join(f"{component:.17e}" for component in vector)
                for vector in MOMENTA
            )
            + "\n",
            encoding="ascii",
        )
        run(["make", "-j4", "amplicol_generate"], cwd=repository)
        run(
            [
                "./amplicol_generate",
                "--library=create",
                "--process=processes.txt",
                "--amplicol_momenta_probe=10",
                "--amplicol_probe_quiet",
                "--timing=none",
            ],
            cwd=repository,
        )
        run(["make", "-j4", "amplicol_generate_library"], cwd=repository)
        run(
            ["make", "-j4", "amplicol_library_benchmark", "amplicol_color_probe"],
            cwd=repository,
        )

        generated = run(
            [
                "./amplicol_library_benchmark",
                "1",
                "1",
                "1",
                momentum_argument,
                *(str(value) for value in HELICITIES),
            ],
            cwd=repository,
        )
        direct = run(
            [
                "./amplicol_color_probe",
                "1",
                "1",
                "1",
                "lc",
                "processes.txt",
                momentum_argument,
                *(str(value) for value in HELICITIES),
            ],
            cwd=repository,
        )
        fixed = scalar(
            generated,
            "AMPICOL_SELECTED_FLOW_PROBE_FIXED_HELICITY_VALUE",
        )
        aggregate = scalar(generated, "AMPICOL_SELECTED_FLOW_PROBE_VALUE")
        direct_value = scalar(direct, "AMPICOL_COLOR_PROBE_VALUE lc")
        expected_helicity = "AMPICOL_SELECTED_FLOW_PROBE_HELICITIES 3 -1 1 -1"
        if expected_helicity not in generated:
            raise AssertionError("generated probe did not authenticate its helicity")
        if not math.isclose(fixed, direct_value, rel_tol=1.0e-13, abs_tol=1.0e-15):
            raise AssertionError(
                f"fixed library value {fixed} differs from direct value {direct_value}"
            )
        if math.isclose(fixed, aggregate, rel_tol=1.0e-13, abs_tol=1.0e-15):
            raise AssertionError(
                "fixed component still contains the summed-helicity multiplicity"
            )
        print("PASS selected-flow fixed-helicity generated-library probe")
    finally:
        if saved_process is None:
            process_file.unlink(missing_ok=True)
        else:
            process_file.write_bytes(saved_process)
        if saved_momentum is None:
            momentum_file.unlink(missing_ok=True)
        else:
            momentum_file.write_bytes(saved_momentum)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
