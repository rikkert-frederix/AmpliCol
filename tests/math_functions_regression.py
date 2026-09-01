#!/usr/bin/env python3
"""Exercise checked integer arithmetic and permutation boundaries."""

from __future__ import annotations

import pathlib
import subprocess
import sys


def run(executable: pathlib.Path, mode: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(executable), mode],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    success = run(executable, "success")
    if success.returncode != 0 or "Math-functions regression: PASS" not in success.stdout:
        raise RuntimeError(f"checked arithmetic success path failed\n{success.stdout}")
    failures = {
        "multiply-overflow": "integer overflow",
        "multiply8-overflow": "64-bit integer overflow",
        "add-overflow": "integer overflow",
        "add8-overflow": "64-bit integer overflow",
        "power-overflow": "integer overflow",
        "factorial-overflow": "factorial does not fit",
        "negative": "negative integer factor",
        "final-permutation": "successor of the final",
    }
    for mode, expected in failures.items():
        result = run(executable, mode)
        if result.returncode == 0 or expected not in result.stdout:
            raise RuntimeError(
                f"arithmetic failure mode {mode!r} was not rejected\n{result.stdout}"
            )
    print("Checked-arithmetic failure regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
