#!/usr/bin/env python3
"""Check that malformed helicity-filter metadata is rejected deliberately."""

from pathlib import Path
import subprocess
import sys


def main() -> None:
    executable = Path(
        sys.argv[1] if len(sys.argv) > 1 else "./amplitude_optimisation_regression"
    ).resolve()
    failures = {
        "missing-filter-offsets": "Helicity-filter subprocess offsets are missing",
        "nonmonotone-filter-offsets": "Invalid helicity-filter subprocess offsets",
    }
    for mode, expected in failures.items():
        result = subprocess.run(
            [str(executable), mode], text=True, capture_output=True, check=False
        )
        output = result.stdout + result.stderr
        if result.returncode == 0:
            raise AssertionError(f"{mode} was accepted\n{output}")
        if expected not in output:
            raise AssertionError(
                f"{mode} failed without the intended diagnostic\n{output}"
            )
    print("Amplitude optimisation hostile-state regression passed")


if __name__ == "__main__":
    main()
