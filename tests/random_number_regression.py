#!/usr/bin/env python3
"""Exercise RANMAR state, maximum-seed, and failure paths."""

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile


def run(executable: pathlib.Path, mode: str, cwd: pathlib.Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(executable), mode],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="amplicol_rng_") as directory:
        cwd = pathlib.Path(directory) / "isolated"
        cwd.mkdir()
        first = run(executable, "max", cwd)
        second = run(executable, "max", cwd)
        if first.returncode != 0 or second.returncode != 0:
            raise RuntimeError(f"maximum seed failed\n{first.stdout}\n{second.stdout}")
        if first.stdout != second.stdout or "RNG_MAX" not in first.stdout:
            raise RuntimeError("maximum seed is not deterministic without offset files")

        repeat = run(executable, "repeat", cwd)
        if repeat.returncode != 0 or "RNG_REPEAT PASS" not in repeat.stdout:
            raise RuntimeError(f"RANMAR reinitialization failed\n{repeat.stdout}")

        failures = {
            "uninitialized": "before initialization",
            "invalid-seed": "Bad initialization value",
            "wide-interval": "interval is too wide",
        }
        for mode, expected in failures.items():
            result = run(executable, mode, cwd)
            if result.returncode == 0 or expected not in result.stdout:
                raise RuntimeError(
                    f"RNG failure mode {mode!r} was not rejected as expected\n{result.stdout}"
                )
    print("Random-number regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
