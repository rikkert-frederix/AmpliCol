#!/usr/bin/env python3
"""Exercise 2->1 process construction, generation, and colour reweighting."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Particle:
    pdg: int
    status: int
    colours: tuple[int, int]
    px: float
    py: float
    pz: float
    energy: float
    mass: float


@dataclass(frozen=True)
class Event:
    weight: float
    scale: float
    particles: tuple[Particle, ...]
    colour_expansion: tuple[float, float, float] | None


def as_float(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def run(command: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stdout}"
        )
    return completed


def internal_pdf_card(source: Path, destination: Path, grid: Path) -> None:
    text = source.read_text(encoding="utf-8")
    text, count = re.subn(
        r"(?m)^(\s*use_lhapdf\s*=\s*)\.true\.\s*$",
        r"\1.false.",
        text,
        count=1,
    )
    if count != 1:
        raise AssertionError("could not select the internal PDF in the run card")
    text, count = re.subn(
        r"(?m)^(\s*internal_pdf_grid\s*=\s*)'[^']*'\s*$",
        lambda match: f"{match.group(1)}'{grid.resolve()}'",
        text,
        count=1,
    )
    if count != 1:
        raise AssertionError("could not set the internal PDF grid in the run card")
    destination.write_text(text, encoding="utf-8")


def process_summary(path: Path) -> tuple[tuple[int, ...], int]:
    lines = [
        line.strip()
        for line in path.read_text(encoding="ascii").splitlines()
        if line.strip()
    ]
    header = lines[0].split()
    n_external, n_unique, version = (int(value) for value in header[:3])
    if n_external != 3 or n_unique != 1 or version != 7:
        raise AssertionError(f"unexpected one-body process header: {header}")
    if not lines[1].startswith("# process:") or not lines[2].startswith("# options:"):
        raise AssertionError("one-body process metadata is missing")
    process = tuple(int(value) for value in lines[3].split())
    family_line = next(
        (line for line in lines if line.startswith("INTEGRATION_FAMILIES ")),
        None,
    )
    if family_line is None:
        raise AssertionError("one-body integration-family catalogue is missing")
    families = int(family_line.split()[1])
    return process, families


def parse_lhe(path: Path) -> tuple[float, float, int, list[Event]]:
    text = path.read_text(encoding="ascii")
    declared_match = re.search(r"<nevents>\s+(\d+)\s+</nevents>", text)
    if declared_match is None:
        raise AssertionError(f"missing event count in {path}")
    init_match = re.search(r"<init>\s*[^\n]+\n\s*([^\n]+)", text)
    if init_match is None:
        raise AssertionError(f"missing cross section in {path}")
    rate_fields = init_match.group(1).split()
    cross_section, error = map(as_float, rate_fields[:2])

    events: list[Event] = []
    for match in re.finditer(r"<event>\s*(.*?)\s*</event>", text, re.DOTALL):
        lines = [line.strip() for line in match.group(1).splitlines() if line.strip()]
        header = lines[0].split()
        nup = int(header[0])
        particles: list[Particle] = []
        for line in lines[1 : 1 + nup]:
            fields = line.split()
            particles.append(
                Particle(
                    int(fields[0]),
                    int(fields[1]),
                    (int(fields[4]), int(fields[5])),
                    *(as_float(value) for value in fields[6:11]),
                )
            )
        expansion = None
        for line in lines[1 + nup :]:
            if line.startswith("#color_expansion"):
                expansion = tuple(as_float(value) for value in line.split()[1:4])
        events.append(
            Event(
                as_float(header[2]),
                as_float(header[3]),
                tuple(particles),
                expansion,
            )
        )
    return cross_section, error, int(declared_match.group(1)), events


def assert_close(label: str, actual: float, expected: float, tolerance: float) -> None:
    if not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise AssertionError(f"{label}: {actual} != {expected}")


def check_kinematics(events: list[Event], pdgs: tuple[int, int, int], mass: float) -> None:
    if not events:
        raise AssertionError("generator wrote no one-body events")
    tau = (mass / 14000.0) ** 2
    for event in events:
        if tuple(particle.pdg for particle in event.particles) != pdgs:
            raise AssertionError("unexpected particles in one-body event")
        incoming1, incoming2, final = event.particles
        if (incoming1.status, incoming2.status, final.status) != (-1, -1, 1):
            raise AssertionError("unexpected one-body LHE particle statuses")
        for component in ("px", "py", "pz", "energy"):
            assert_close(
                f"four-momentum conservation ({component})",
                getattr(final, component),
                getattr(incoming1, component) + getattr(incoming2, component),
                2.0e-10,
            )
        invariant = final.energy**2 - final.px**2 - final.py**2 - final.pz**2
        assert_close("final-state invariant mass", invariant, mass**2, 2.0e-10)
        assert_close("reported final-state mass", final.mass, mass, 2.0e-12)
        assert_close("one-body transverse momentum", math.hypot(final.px, final.py), 0.0, 1.0e-12)
        x1 = incoming1.energy / 7000.0
        x2 = incoming2.energy / 7000.0
        if not (0.0 < x1 < 1.0 and 0.0 < x2 < 1.0):
            raise AssertionError("one-body Bjorken x is outside (0,1)")
        assert_close("fixed Bjorken product", x1 * x2, tau, 2.0e-11)
        assert_close("one-body scale", event.scale, mass / 2.0, 2.0e-10)


def generate_and_reweight(
    directory: Path,
    process_list: Path,
    generator: Path,
    reweighter: Path,
    input_card: Path,
    process: str,
    expected_pdgs: tuple[int, int, int],
    mass: float,
    expected_nlc: float,
    expected_full: float,
    heft: bool,
) -> None:
    directory.mkdir()
    (directory / "Outputs").mkdir()
    process_command = [sys.executable, str(process_list)]
    if heft:
        process_command.append("--heft")
    process_command.extend(["--serial", process])
    run(process_command, directory)
    serialized_process, groups = process_summary(directory / "processes.txt")
    if serialized_process != expected_pdgs or groups != 1:
        raise AssertionError(
            f"expected one channel for {process}, got {groups}: {serialized_process}"
        )

    tag = "higgs" if heft else "zboson"
    event_file = directory / "Outputs" / f"{tag}_events.lhe"
    run(
        [
            str(generator),
            "--process=processes.txt",
            f"--input={input_card}",
            "--nevents=32",
            "--itmax=8",
            "--seed=260822",
            f"--tag={tag}",
            "--timing=none",
        ],
        directory,
    )
    leading_rate, leading_error, declared, leading_events = parse_lhe(event_file)
    if declared != 32 or len(leading_events) != 32:
        raise AssertionError("one-body LHE event count does not match the request")
    if not (math.isfinite(leading_rate) and leading_rate > 0.0):
        raise AssertionError("one-body integration did not produce a positive rate")
    if not (math.isfinite(leading_error) and leading_error >= 0.0):
        raise AssertionError("one-body integration produced an invalid error")
    check_kinematics(leading_events, expected_pdgs, mass)

    run([str(reweighter), str(event_file), f"--input={input_card}"], directory)
    full_rate, _, full_declared, full_events = parse_lhe(
        Path(f"{event_file}.rwgt")
    )
    if full_declared != 32 or len(full_events) != 32:
        raise AssertionError("reweighted one-body LHE lost events")
    assert_close(
        "full/leading cross section", full_rate / leading_rate, expected_full, 2.0e-7
    )
    for event in full_events:
        if event.colour_expansion is None:
            raise AssertionError("reweighted event has no colour expansion")
        lc, nlc, full = event.colour_expansion
        assert_close("event NLC/LC ratio", nlc / lc, expected_nlc, 2.0e-7)
        assert_close("event full/LC ratio", full / lc, expected_full, 2.0e-7)
    check_kinematics(full_events, expected_pdgs, mass)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--generator", type=Path, required=True)
    parser.add_argument("--reweighter", type=Path, required=True)
    parser.add_argument("--process-list", type=Path, required=True)
    parser.add_argument("--input-card", type=Path, required=True)
    parser.add_argument("--pdf-grid", type=Path, required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="amplicol-onebody-") as raw:
        root = Path(raw)
        input_card = root / "run_card.dat"
        internal_pdf_card(args.input_card, input_card, args.pdf_grid)
        generate_and_reweight(
            root / "higgs",
            args.process_list.resolve(),
            args.generator.resolve(),
            args.reweighter.resolve(),
            input_card,
            "g g > h",
            (21, 21, 25),
            125.0,
            7.0 / 9.0,
            8.0 / 9.0,
            True,
        )
        generate_and_reweight(
            root / "zboson",
            args.process_list.resolve(),
            args.generator.resolve(),
            args.reweighter.resolve(),
            input_card,
            "u u~ > z",
            (2, -2, 23),
            91.188,
            1.0,
            1.0,
            False,
        )

    print("One-body generation and reweighting regression passed")


if __name__ == "__main__":
    main()
