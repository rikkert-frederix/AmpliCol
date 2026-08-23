#!/usr/bin/env python3
"""Summarize the matched AmpliCol/MadGraph gg > h g g g rates."""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class LHERate:
    cross_section: float
    integration_error: float
    declared_events: int
    events: int
    weight_mean: float
    weight_sem: float


def parse_lhe(path: Path) -> LHERate:
    cross_section = None
    integration_error = None
    declared_events = None
    events = 0
    mean = 0.0
    moment2 = 0.0
    need_init_beams = False
    need_init_rate = False
    need_event_header = False

    with path.open(encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if line.startswith("<nevents>"):
                declared_events = int(line.split()[1])
            if line == "<init>":
                need_init_beams = True
                continue
            if need_init_beams and line:
                need_init_beams = False
                need_init_rate = True
                continue
            if need_init_rate and line:
                fields = line.split()
                cross_section = float(fields[0])
                integration_error = float(fields[1])
                need_init_rate = False
                continue
            if line == "<event>":
                need_event_header = True
                continue
            if need_event_header and line:
                weight = float(line.split()[2])
                events += 1
                delta = weight - mean
                mean += delta / events
                moment2 += delta * (weight - mean)
                need_event_header = False

    if cross_section is None or integration_error is None:
        raise ValueError(f"no LHE init rate found in {path}")
    if declared_events is None:
        raise ValueError(f"no <nevents> metadata found in {path}")
    if events < 2:
        raise ValueError(f"fewer than two events found in {path}")
    sem = math.sqrt(moment2 / (events - 1) / events)
    return LHERate(
        cross_section,
        integration_error,
        declared_events,
        events,
        mean,
        sem,
    )


def parse_madgraph_results(path: Path) -> tuple[float, float]:
    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    fields = first_line.split()
    if len(fields) < 2:
        raise ValueError(f"invalid MadGraph results file {path}")
    # MadEvent writes the final combined central value in field 10, while
    # fields 1 and 2 contain the integration estimate and its uncertainty.
    # Survey-only files can leave field 10 at zero.
    combined = float(fields[9]) if len(fields) >= 10 else 0.0
    cross_section = combined if combined != 0.0 else float(fields[0])
    return cross_section, float(fields[1])


def parse_madgraph_banner(path: Path) -> tuple[int, float]:
    text = path.read_text(encoding="utf-8")
    event_match = re.search(r"Number of Events\s*:\s*(\d+)", text)
    rate_match = re.search(
        r"Integrated weight \(pb\)\s*:\s*([+\-0-9.eE]+)", text
    )
    if event_match is None or rate_match is None:
        raise ValueError(f"no MadGraph generation summary found in {path}")
    return int(event_match.group(1)), float(rate_match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("leading_colour_lhe", type=Path)
    parser.add_argument("full_colour_lhe", type=Path)
    parser.add_argument("madgraph_results", type=Path)
    parser.add_argument("--madgraph-banner", type=Path)
    parser.add_argument("--expected-events", type=int, default=100_000)
    args = parser.parse_args()

    leading = parse_lhe(args.leading_colour_lhe)
    full = parse_lhe(args.full_colour_lhe)
    madgraph_xsec, madgraph_error = parse_madgraph_results(
        args.madgraph_results
    )
    madgraph_events = None
    if args.madgraph_banner is not None:
        madgraph_events, banner_xsec = parse_madgraph_banner(
            args.madgraph_banner
        )
        if not math.isclose(
            madgraph_xsec, banner_xsec, rel_tol=1e-4, abs_tol=1e-8
        ):
            raise SystemExit(
                "MadGraph results and generated-event banner disagree: "
                f"{madgraph_xsec} versus {banner_xsec} pb"
            )
        madgraph_xsec = banner_xsec

    if leading.events != args.expected_events:
        raise SystemExit(
            f"expected {args.expected_events} leading-colour events, "
            f"found {leading.events}"
        )
    if full.events != args.expected_events:
        raise SystemExit(
            f"expected {args.expected_events} full-colour events, "
            f"found {full.events}"
        )
    if leading.declared_events != leading.events:
        raise SystemExit("leading-colour LHE event metadata does not match")
    if full.declared_events != full.events:
        raise SystemExit("full-colour LHE event metadata does not match")

    combined_error = math.hypot(full.integration_error, madgraph_error)
    difference = full.cross_section - madgraph_xsec
    pull = difference / combined_error
    weight_ratio = full.weight_mean / leading.weight_mean

    print(f"events:                    {full.events}")
    print(
        "AmpliCol leading colour:  "
        f"{leading.cross_section:.10g} +/- "
        f"{leading.integration_error:.4g} pb"
    )
    print(
        "AmpliCol full colour:     "
        f"{full.cross_section:.10g} +/- {full.integration_error:.4g} pb"
    )
    print(f"event-weight SEM:          {full.weight_sem:.4g} pb")
    print(f"full/leading weight ratio: {weight_ratio:.12g}")
    print(
        "MadGraph full colour:     "
        f"{madgraph_xsec:.10g} +/- {madgraph_error:.4g} pb"
    )
    if madgraph_events is not None:
        print(f"MadGraph events:           {madgraph_events}")
    print(f"AmpliCol - MadGraph:       {difference:+.6g} pb")
    print(f"combined-error pull:       {pull:+.4f} sigma")


if __name__ == "__main__":
    main()
