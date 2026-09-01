#!/usr/bin/env python3
"""Exercise successful and rejected amplitude-serialization states."""

from __future__ import annotations

import subprocess
import sys


def run(executable: str, mode: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [executable, mode],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: amplitude_serialization_regression.py EXECUTABLE")
        return 2
    executable = sys.argv[1]
    successful = run(executable, "roundtrip")
    if successful.returncode != 0 or "amplitude serialization regression: PASS" not in successful.stdout:
        print(successful.stdout, end="")
        print("valid amplitude serialization failed")
        return 1

    rejected_modes = (
        "missing_current_ranges",
        "short_current_vertices",
        "malformed_singlet_map",
        "invalid_permutation",
        "invalid_process",
        "bad_offsets",
        "bad_process_mask",
        "orphan_same_flavour_operation",
        "extra_spin_multiplicity",
    )
    for mode in rejected_modes:
        result = run(executable, mode)
        if result.returncode == 0 or "Invalid amplitude serialization state" not in result.stdout:
            print(result.stdout, end="")
            print(f"malformed state was not rejected safely: {mode}")
            return 1

    print("amplitude serialization malformed-state regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
