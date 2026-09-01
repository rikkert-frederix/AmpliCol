#!/usr/bin/env python3
"""Exercise valid and malformed command-line parser inputs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run(executable: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(executable), *arguments],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def require_success(executable: Path, arguments: list[str], marker: str) -> None:
    result = run(executable, *arguments)
    if result.returncode != 0 or marker not in result.stdout:
        raise AssertionError(
            f"valid parser invocation failed ({result.returncode}): {arguments}\n{result.stdout}"
        )


def require_failure(executable: Path, *arguments: str) -> None:
    result = run(executable, *arguments)
    if result.returncode == 0:
        raise AssertionError(f"malformed parser invocation succeeded: {arguments}\n{result.stdout}")
    if "Fortran runtime error" in result.stdout:
        raise AssertionError(f"parser leaked a runtime conversion error: {arguments}\n{result.stdout}")


def main() -> int:
    executable = Path(sys.argv[1] if len(sys.argv) > 1 else "tests/command_line_parser_driver.exe").resolve()
    require_success(executable, [], "PARSER_OK 10000 128 1 0 100 F F")
    require_success(
        executable,
        [
            "--process=born.dat",
            "--real-process=real.dat",
            "--input=card.dat",
            "--nevents=7",
            "--seed=904866561",
            "--phasespace=4",
            "--itmax=3",
            "--library=none",
            "--tag=safe-tag",
            "--me_test=2",
            "--timing=detailed",
            "--timing-sample=1",
            "--accuracy=.2d0",
            "--dim-reg=fdh",
            "--tail-replay=tail.dat",
            "--migration-tail-fraction=0",
        ],
        "PARSER_OK 7 3 4 2 1 T T",
    )

    malformed = [
        ("--helpful",),
        ("--nevents=",),
        ("--nevents=abc",),
        ("--nevents=1 2",),
        ("--nevents=999999999999999999999",),
        ("--nevents=0",),
        ("--seed=-1",),
        ("--seed=904866562",),
        ("--phasespace=0",),
        ("--itmax=0",),
        ("--me_test=0",),
        ("--timing-sample=0",),
        ("--accuracy=NaN",),
        ("--accuracy=1",),
        ("--accuracy=.5junk",),
        ("--migration-tail-fraction=NaN",),
        ("--migration-tail-fraction=1.1",),
        ("--process=",),
        ("--tag=../escape",),
        ("--tag=" + "x" * 80,),
        ("--process=" + "x" * 81,),
        ("--input=" + "x" * 260,),
    ]
    for arguments in malformed:
        require_failure(executable, *arguments)

    print("Command-line parser regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
