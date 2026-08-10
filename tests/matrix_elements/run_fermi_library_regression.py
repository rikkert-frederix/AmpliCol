#!/usr/bin/env python3
"""Build and run the generated identical-lepton amplitude-library check."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import tempfile
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", required=True)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--fflags", default="")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    generator = Path(args.generator).resolve()
    verifier = root / "tests/matrix_elements/fermi_statistics_library_regression.f03"
    with tempfile.TemporaryDirectory(prefix="amplicol-fermi-library-") as tmp_name:
        tmp = Path(tmp_name)
        (tmp / "Library").mkdir()
        subprocess.run([str(generator), "--emit-library"], cwd=tmp, check=True)

        generated = tmp / "Library/amp1_1_lib.f03"
        checker = tmp / "fermi_statistics_library_regression"
        command = [
            *shlex.split(args.compiler),
            *shlex.split(args.fflags),
            "-I",
            str(root),
            "-o",
            str(checker),
            str(root / "feynmanrules.o"),
            str(generated),
            str(verifier),
        ]
        subprocess.run(command, cwd=tmp, check=True)
        subprocess.run([str(checker)], cwd=tmp, check=True)


if __name__ == "__main__":
    main()
