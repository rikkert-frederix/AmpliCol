#!/usr/bin/env python3
"""Exercise failure paths in the fixed-size legacy internal-PDF reader."""

from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile


def run_failure(executable: pathlib.Path, grid: pathlib.Path, expected: str, member: int = 0) -> None:
    result = subprocess.run(
        [str(executable), str(grid), str(member)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode == 0:
        raise RuntimeError(f"malformed grid unexpectedly succeeded: {grid}\n{result.stdout}")
    if expected not in result.stdout:
        raise RuntimeError(
            f"failure for {grid} did not report {expected!r}\n{result.stdout}"
        )


def header(nx: int) -> str:
    return "\n".join(
        [
            "'Parameterlist:'",
            "'list', 0, 1, 0.119",
            "'NNPDF20int'",
            "0, 1",
            str(nx),
        ]
    ) + "\n"


def complete_grid(x_values: list[float], alternating_huge: bool = False) -> str:
    rows = [
        "'Parameterlist:'",
        "'list', 0, 1, 0.119",
        "'NNPDF20int'",
        "0, 1",
        str(len(x_values)),
        *(f"{value:.17e}" for value in x_values),
        "2",
        "'Q2grid'",
        "1.0",
        "1.0e8",
        "'Pdfdata'",
    ]
    for _ in x_values:
        for iq in range(2):
            if alternating_huge:
                value = "1.7e308" if iq == 0 else "-1.7e308"
            else:
                value = "1.0"
            rows.append(" ".join([value] * 13))
    return "\n".join(rows) + "\n"


def main() -> int:
    executable = pathlib.Path(sys.argv[1]).resolve()
    real_grid = pathlib.Path(sys.argv[2]).resolve()
    subprocess.run([str(executable), str(real_grid), "0"], check=True)
    with tempfile.TemporaryDirectory(prefix="amplicol_pdf_") as directory:
        root = pathlib.Path(directory)
        oversized = root / "oversized.grid"
        oversized.write_text(header(101), encoding="utf-8")
        run_failure(executable, oversized, "x-grid size exceeds")

        truncated = root / "truncated.grid"
        truncated.write_text(header(4) + "0.1\n", encoding="utf-8")
        run_failure(executable, truncated, "reading the PDF x grid")

        unsorted = root / "unsorted.grid"
        unsorted.write_text(header(4) + "0.1\n0.2\n0.15\n0.3\n", encoding="utf-8")
        run_failure(executable, unsorted, "x grid is not strictly increasing")

        unphysical_x = root / "unphysical_x.grid"
        unphysical_x.write_text(
            complete_grid([0.01, 0.1, 0.5, 1.1]), encoding="utf-8"
        )
        run_failure(executable, unphysical_x, "outside (0,1]")

        nan_alpha = root / "nan_alpha.grid"
        nan_alpha.write_text(
            "'Parameterlist:'\n'list', 0, 1, NaN\n", encoding="utf-8"
        )
        run_failure(executable, nan_alpha, "invalid alpha_s")

        nan_x = root / "nan_x.grid"
        nan_x.write_text(header(4) + "NaN\n", encoding="utf-8")
        run_failure(executable, nan_x, "outside (0,1]")

        nan_q2 = root / "nan_q2.grid"
        nan_q2.write_text(
            header(4) + "0.01\n0.1\n0.5\n1.0\n2\n'Q2grid'\nNaN\n",
            encoding="utf-8",
        )
        run_failure(executable, nan_q2, "nonpositive or non-finite")

        overflowing = root / "overflowing.grid"
        overflowing.write_text(
            complete_grid([0.01, 0.1, 0.5, 1.0], alternating_huge=True),
            encoding="utf-8",
        )
        run_failure(executable, overflowing, "unsafe interpolation arithmetic")

        run_failure(executable, real_grid, "replica out of range", member=1)
    print("Internal-PDF malformed-input regression: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
