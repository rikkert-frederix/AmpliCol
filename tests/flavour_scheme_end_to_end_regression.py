#!/usr/bin/env python3
"""Exercise flavour-scheme metadata through process and library generation."""

from __future__ import annotations

import argparse
import struct
import subprocess
import sys
import tempfile
from pathlib import Path


def generate_scheme(
    root: Path,
    generator: Path,
    process_list: Path,
    input_card: Path,
    scheme: int,
    flavour: str,
    mass: float,
) -> None:
    workdir = root / f"fs{scheme}"
    workdir.mkdir()
    (workdir / "Outputs").mkdir()
    (workdir / "Library").mkdir()

    subprocess.run(
        [
            sys.executable,
            str(process_list),
            "--serial",
            "-FS",
            str(scheme),
            f"p p > {flavour} {flavour}~ h",
        ],
        cwd=workdir,
        check=True,
    )

    process_lines = (workdir / "processes.txt").read_text(
        encoding="ascii"
    ).splitlines()
    if process_lines[0].split()[2] != "4":
        raise AssertionError("flavour-scheme process file is not version 4")
    if f"flavour_scheme={scheme}" not in process_lines[2]:
        raise AssertionError("process file lost flavour-scheme metadata")

    tag = f"fs{scheme}"
    subprocess.run(
        [
            str(generator),
            "--process=processes.txt",
            f"--input={input_card}",
            "--nevents=1",
            "--phasespace=1",
            "--library=create",
            f"--tag={tag}",
        ],
        cwd=workdir,
        check=True,
        capture_output=True,
        text=True,
    )

    log_lines = (workdir / "Outputs" / f"{tag}_log_file.txt").read_text(
        encoding="ascii"
    ).splitlines()
    if f"Active flavour scheme: {scheme}" not in log_lines:
        raise AssertionError("event generator initialized the wrong flavour scheme")

    phase_space_masses = []
    for line in log_lines:
        if line.strip().startswith("masses:"):
            phase_space_masses.append(
                [float(value) for value in line.split()[1:]]
            )
    if not any(
        sum(abs(value - mass) < 1.0e-12 for value in masses) >= 2
        for masses in phase_space_masses
    ):
        raise AssertionError(
            f"FS{scheme} did not use massive {flavour} external legs"
        )

    metadata = (workdir / "Library" / "amplitudes.bin").read_bytes()[:8]
    library_version, library_scheme = struct.unpack("=ii", metadata)
    if (library_version, library_scheme) != (6, scheme):
        raise AssertionError(
            "amplitude library lost format/flavour-scheme metadata: "
            f"{library_version}, {library_scheme}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--process-list", type=Path, required=True)
    parser.add_argument("--input-card", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with tempfile.TemporaryDirectory(prefix="amplicol-flavour-scheme-") as tmp:
        root = Path(tmp)
        generate_scheme(
            root,
            args.generator.resolve(),
            args.process_list.resolve(),
            args.input_card.resolve(),
            4,
            "b",
            4.7,
        )
        generate_scheme(
            root,
            args.generator.resolve(),
            args.process_list.resolve(),
            args.input_card.resolve(),
            3,
            "c",
            1.42,
        )
    print("Flavour-scheme process/library end-to-end regression passed")


if __name__ == "__main__":
    main()
