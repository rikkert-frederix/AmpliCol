#!/usr/bin/env python3
"""Exercise successful and controlled-failure event-handling paths."""

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
    if len(sys.argv) != 2:
        raise SystemExit("usage: handling_events_regression.py EXECUTABLE")
    executable = pathlib.Path(sys.argv[1]).resolve()

    success = run(executable, "success")
    if success.returncode != 0 or "Handling-events regression: PASS" not in success.stdout:
        raise AssertionError(
            f"successful event regression failed ({success.returncode}):\n{success.stdout}"
        )

    failure_modes = (
        "nan-weight",
        "nan-results",
        "header-missing",
        "header-extra",
        "header-nan",
        "header-huge",
        "header-process",
        "init-nproc",
        "init-nan-beam",
        "init-extra",
        "init-bad-process",
        "bad-multiplicity",
        "shifted-process-bounds",
        "stale-process-group",
    )
    for mode in failure_modes:
        result = run(executable, mode)
        if result.returncode == 0:
            raise AssertionError(f"{mode} was unexpectedly accepted:\n{result.stdout}")
        # A negative return code denotes a signal on POSIX.  Invalid data must
        # reach an intentional STOP, not SIGFPE/SIGSEGV/SIGABRT.
        if result.returncode < 0:
            raise AssertionError(
                f"{mode} terminated by signal {-result.returncode}:\n{result.stdout}"
            )
        if "unexpectedly accepted" in result.stdout:
            raise AssertionError(f"{mode} reached the regression fallthrough:\n{result.stdout}")

    print("Handling-events failure regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
