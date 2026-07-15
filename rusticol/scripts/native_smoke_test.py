#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import os
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str | Path]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(part) for part in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Build and run Rusticol C++/Fortran smoke tests.")
    parser.add_argument("process_output", type=Path)
    parser.add_argument("--config", default="rusticol-config")
    args = parser.parse_args()

    config = shutil.which(args.config) or args.config
    cflags = shlex.split(run([config, "--cflags"]).stdout)
    libs = shlex.split(run([config, "--libs"]).stdout)
    fortran_module = Path(run([config, "--fortran-module"]).stdout.strip())
    cxx = os.environ.get("CXX", "c++")
    fc = os.environ.get("FC", "gfortran")

    with tempfile.TemporaryDirectory(prefix="rusticol-native-smoke-") as temporary:
        build = Path(temporary)
        cpp_exe = build / "cpp-smoke"
        fortran_exe = build / "fortran-smoke"
        run(
            [
                cxx,
                "-std=c++17",
                "-O2",
                *cflags,
                ROOT / "tests" / "native_smoke.cpp",
                *libs,
                "-o",
                cpp_exe,
            ]
        )
        run(
            [
                fc,
                "-std=f2008",
                "-O2",
                fortran_module,
                ROOT / "tests" / "native_smoke.f90",
                *libs,
                "-o",
                fortran_exe,
            ]
        )
        cpp_values = parse_values(run([cpp_exe, args.process_output]).stdout)
        fortran_values = parse_values(run([fortran_exe, args.process_output]).stdout)
    for left, right in zip(cpp_values, fortran_values, strict=True):
        if not math.isclose(left, right, rel_tol=1.0e-12, abs_tol=1.0e-15):
            raise RuntimeError(f"native smoke mismatch: C++={left:.17e}, Fortran={right:.17e}")
    print(f"native smoke passed: {cpp_values[0]:.17e}")
    return 0


def parse_values(output: str) -> tuple[float, float]:
    fields = output.strip().split()
    if len(fields) != 2:
        raise RuntimeError(f"expected two smoke-test values, got {output!r}")
    return float(fields[0]), float(fields[1])


if __name__ == "__main__":
    raise SystemExit(main())
