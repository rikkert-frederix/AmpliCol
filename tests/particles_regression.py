#!/usr/bin/env python3
"""Exercise particle-model integer boundaries and invalid lifecycle states."""

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
    if success.returncode != 0 or "Particles regression: PASS" not in success.stdout:
        raise RuntimeError(f"particle-model success path failed\n{success.stdout}")

    failures = {
        "nan-width": "invalid particle width",
        "huge-width": "invalid particle width",
        "tiny-width": "invalid particle width",
        "vertices-before-particles": "before particles",
        "stale-vertices": "has not been initialised",
    }
    for mode, expected in failures.items():
        result = run(executable, mode)
        if result.returncode == 0 or expected not in result.stdout:
            raise RuntimeError(
                f"particle-model failure mode {mode!r} was not rejected\n{result.stdout}"
            )
    print("Particle-model failure regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
