#!/usr/bin/env python3
"""Validate MadGraph reference parsing and relative comparisons synthetically."""

from __future__ import annotations

import math
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


def write_reference(path: pathlib.Path, count: int, values: list[str], momentum: str = "0 0 0 0") -> None:
    rows = [momentum] * 4 + [str(count)] + values
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    normalization = (4.0 * math.pi * 0.118) ** 2
    with tempfile.TemporaryDirectory(prefix="amplicol_mg_") as directory:
        cwd = pathlib.Path(directory)
        reference_dir = cwd / "Utilities" / "ME_checks"
        reference_dir.mkdir(parents=True)
        reference = reference_dir / "momenta_1_1.txt"

        write_reference(reference, 2, ["0.0", f"{normalization * (1.0 + 5.0e-5):.17e}"])
        passed = run(executable, "pass", cwd)
        if passed.returncode != 0 or "comparison regression: PASS" not in passed.stdout:
            raise RuntimeError(f"valid MadGraph comparison failed\n{passed.stdout}")

        write_reference(reference, 1, ["0.0"])
        zero = run(executable, "zero", cwd)
        if zero.returncode != 0 or "zero comparison regression: PASS" not in zero.stdout:
            raise RuntimeError(f"zero MadGraph comparison failed\n{zero.stdout}")

        cases: list[tuple[str, str]] = []
        write_reference(reference, 1, [f"{normalization * 1.01:.17e}"])
        cases.append(("mismatch", "disagreement with MadGraph"))
        for mode, expected in cases:
            result = run(executable, mode, cwd)
            if result.returncode == 0 or expected not in result.stdout:
                raise RuntimeError(f"{mode} was not rejected\n{result.stdout}")

        write_reference(reference, 1, ["-1.0"])
        result = run(executable, "bad-reference", cwd)
        if result.returncode == 0 or "invalid matrix elements" not in result.stdout:
            raise RuntimeError(f"negative reference was not rejected\n{result.stdout}")

        write_reference(reference, 1, ["NaN"])
        result = run(executable, "bad-reference", cwd)
        if result.returncode == 0 or "invalid matrix elements" not in result.stdout:
            raise RuntimeError(f"non-finite reference was not rejected\n{result.stdout}")

        write_reference(reference, 1, ["1.0"], momentum="NaN 0 0 0")
        result = run(executable, "bad-momentum", cwd)
        if result.returncode == 0 or "non-finite momentum" not in result.stdout:
            raise RuntimeError(f"non-finite momentum was not rejected\n{result.stdout}")

        write_reference(reference, 1_000_001, [])
        result = run(executable, "bad-count", cwd)
        if result.returncode == 0 or "Invalid MadGraph reference count" not in result.stdout:
            raise RuntimeError(f"oversized reference count was not rejected\n{result.stdout}")

        write_reference(reference, 1, [])
        result = run(executable, "truncated", cwd)
        if result.returncode == 0 or "Could not read MadGraph reference value" not in result.stdout:
            raise RuntimeError(f"truncated reference was not rejected\n{result.stdout}")

        write_reference(reference, 1, ["1.0", "2.0"])
        result = run(executable, "trailing", cwd)
        if result.returncode == 0 or "contains trailing data" not in result.stdout:
            raise RuntimeError(f"trailing reference data was not rejected\n{result.stdout}")

        for mode, expected in {
            "overflow": "overflows MadGraph comparison normalization",
            "coupling-overflow": "normalization cannot be represented",
            "negative-code": "Invalid normalization metadata",
            "bad-next": "Invalid process dimensions",
        }.items():
            result = run(executable, mode, cwd)
            if result.returncode == 0 or expected not in result.stdout:
                raise RuntimeError(f"{mode} state was not rejected\n{result.stdout}")
    print("MadGraph malformed-reference regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
