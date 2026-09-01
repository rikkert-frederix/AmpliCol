#!/usr/bin/env python3
"""Run bitset success and fail-closed regression modes."""

from __future__ import annotations

import pathlib
import subprocess
import sys


def invoke(executable: pathlib.Path, mode: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(executable), mode],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    success = invoke(executable, "roundtrip")
    if success.returncode != 0 or "Bitset regression: PASS" not in success.stdout:
        raise RuntimeError(f"bitset round trip failed\n{success.stdout}")
    failures = {
        "mismatch": "different lengths",
        "out-of-range": "outside the bitset",
        "malformed": "malformed bitset",
        "too-large-integer": "too large to fit",
    }
    for mode, expected in failures.items():
        result = invoke(executable, mode)
        if result.returncode == 0 or expected not in result.stdout:
            raise RuntimeError(
                f"bitset failure mode {mode!r} was not rejected\n{result.stdout}"
            )
    print("Bitset malformed-state regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
