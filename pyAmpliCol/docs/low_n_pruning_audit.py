#!/usr/bin/env python3
"""Compare pyAmpliCol generic-DAG counts with a generated AmpliCol library.

This helper is intentionally read-only.  It does not drive AmpliCol generation;
prepare the Fortran library separately, behind the usual RAM watchdog for heavy
cases, then point this script at a `generic_process_manifest.json` and one
`Library/amp*_lib.f03` module.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


@dataclass(frozen=True)
class CountSummary:
    process: str
    currents: int
    interactions: int
    amplitude_roots: int
    stage_kind_counts: dict[tuple[int, int], int]


def pyamplicol_counts(manifest_path: Path) -> CountSummary:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    currents = payload.get("currents", ())
    interactions = payload.get("interactions", ())
    roots = payload.get("amplitude_roots", ())
    current_by_id = {int(current["id"]): current for current in currents}
    stage_kind_counts: Counter[tuple[int, int]] = Counter()
    for interaction in interactions:
        result = current_by_id[int(interaction["result_id"])]
        stage_size = len(result["index"]["external_labels"])
        stage_kind_counts[(stage_size, int(interaction["vertex_kind"]))] += 1
    return CountSummary(
        process=str(payload.get("process", manifest_path.stem)),
        currents=len(currents),
        interactions=len(interactions),
        amplitude_roots=len(roots),
        stage_kind_counts=dict(sorted(stage_kind_counts.items())),
    )


def fortran_counts(module_path: Path) -> CountSummary:
    source = module_path.read_text(encoding="utf-8", errors="ignore")
    dimensions = [
        int(value)
        for value in re.findall(
            r"complex\(kind=8\),dimension\(1:6,\s*(\d+)\)",
            source,
        )
    ]
    if len(dimensions) < 2:
        raise ValueError(f"could not read val_c/int_c dimensions in {module_path}")
    stage_kind_counts: Counter[tuple[int, int]] = Counter()
    for match in re.finditer(
        r"subroutine vertex_type(\d+)_(\d+)(?:_v\d+)?\(.*?\n(.*?)end subroutine",
        source,
        flags=re.S,
    ):
        stage_size = int(match.group(1))
        vertex_kind = int(match.group(2))
        body = match.group(3)
        loop_counts = [int(value) for value in re.findall(r"do i=1,\s*(\d+)", body)]
        stage_kind_counts[(stage_size, vertex_kind)] += sum(loop_counts)
    amplitude_roots = len(re.findall(r"^\s*amps\(\d+\)\s*=", source, flags=re.M))
    return CountSummary(
        process=module_path.stem,
        currents=dimensions[0],
        interactions=dimensions[1],
        amplitude_roots=amplitude_roots,
        stage_kind_counts=dict(sorted(stage_kind_counts.items())),
    )


def markdown_report(py_counts: CountSummary, ft_counts: CountSummary) -> str:
    lines = [
        f"# Low-n Count Delta: `{py_counts.process}`",
        "",
        "| source | currents | interactions | amplitude roots |",
        "| --- | ---: | ---: | ---: |",
        (
            f"| pyAmpliCol | {py_counts.currents} | {py_counts.interactions} "
            f"| {py_counts.amplitude_roots} |"
        ),
        (
            f"| Fortran AmpliCol | {ft_counts.currents} | {ft_counts.interactions} "
            f"| {ft_counts.amplitude_roots} |"
        ),
        (
            f"| delta | {py_counts.currents - ft_counts.currents:+d} "
            f"| {py_counts.interactions - ft_counts.interactions:+d} "
            f"| {py_counts.amplitude_roots - ft_counts.amplitude_roots:+d} |"
        ),
        "",
        "| stage size | vertex kind | pyAmpliCol | Fortran | delta |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    keys = sorted(set(py_counts.stage_kind_counts) | set(ft_counts.stage_kind_counts))
    for stage_size, vertex_kind in keys:
        py_value = py_counts.stage_kind_counts.get((stage_size, vertex_kind), 0)
        ft_value = ft_counts.stage_kind_counts.get((stage_size, vertex_kind), 0)
        lines.append(
            f"| {stage_size} | {vertex_kind} | {py_value} | {ft_value} "
            f"| {py_value - ft_value:+d} |"
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--fortran-module", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)

    report = markdown_report(
        pyamplicol_counts(args.manifest),
        fortran_counts(args.fortran_module),
    )
    if args.output is None:
        print(report)
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report, encoding="utf-8")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
