#!/usr/bin/env python3
"""Apply the matched HEFT comparison settings to a MadGraph output."""

from __future__ import annotations

import argparse
import shlex
from pathlib import Path


def parse_settings(path: Path) -> tuple[dict[str, str], dict[tuple[str, str], str]]:
    """Parse `set run_card key value` and `set param_card block id value`."""
    run_updates: dict[str, str] = {}
    param_updates: dict[tuple[str, str], str] = {}
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        fields = shlex.split(raw_line, comments=True)
        if not fields:
            continue
        if fields[:2] == ["set", "run_card"] and len(fields) == 4:
            run_updates[fields[2].lower()] = fields[3]
        elif fields[:2] == ["set", "param_card"] and len(fields) == 5:
            param_updates[(fields[2].lower(), fields[3])] = fields[4]
        else:
            raise ValueError(f"invalid setting at {path}:{line_number}")
    return run_updates, param_updates


def update_run_card(path: Path, updates: dict[str, str]) -> None:
    found: set[str] = set()
    output: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        data, separator, comment = raw_line.partition("!")
        if "=" not in data:
            output.append(raw_line)
            continue
        _, right = data.split("=", 1)
        fields = right.split()
        key = fields[0].lower() if fields else ""
        if key not in updates:
            output.append(raw_line)
            continue
        ending = "\n" if raw_line.endswith("\n") else ""
        suffix = f" !{comment.rstrip()}" if separator else ""
        output.append(f"  {updates[key]}\t= {key}{suffix}{ending}")
        found.add(key)
    missing = set(updates) - found
    if missing:
        raise ValueError(f"run-card keys not found in {path}: {sorted(missing)}")
    path.write_text("".join(output), encoding="utf-8")


def replace_lha_value(line: str, value: str, value_field: int = 1) -> str:
    ending = "\n" if line.endswith("\n") else ""
    body = line[:-1] if ending else line
    data, separator, comment = body.partition("#")
    fields = data.split()
    fields[value_field] = value
    suffix = f" #{comment}" if separator else ""
    return f"    {' '.join(fields)}{suffix}{ending}"


def update_param_card(path: Path, updates: dict[tuple[str, str], str]) -> None:
    found: set[tuple[str, str]] = set()
    current_block = ""
    output: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        fields = raw_line.split()
        if len(fields) >= 2 and fields[0].lower() == "block":
            current_block = fields[1].lower()
            output.append(raw_line)
            continue
        if len(fields) >= 3 and fields[0].lower() == "decay":
            target = ("decay", fields[1])
            if target in updates:
                output.append(
                    replace_lha_value(raw_line, updates[target], value_field=2)
                )
                found.add(target)
            else:
                output.append(raw_line)
            current_block = ""
            continue
        if current_block and len(fields) >= 2 and not raw_line.lstrip().startswith("#"):
            target = (current_block, fields[0])
            if target in updates:
                output.append(replace_lha_value(raw_line, updates[target]))
                found.add(target)
                continue
        output.append(raw_line)
    missing = set(updates) - found
    if missing:
        raise ValueError(f"parameter-card entries not found in {path}: {sorted(missing)}")
    path.write_text("".join(output), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("madgraph_output", type=Path)
    parser.add_argument(
        "--settings",
        type=Path,
        default=Path(__file__).with_name("madgraph_settings.txt"),
    )
    args = parser.parse_args()

    run_updates, param_updates = parse_settings(args.settings)
    cards = args.madgraph_output / "Cards"
    update_run_card(cards / "run_card.dat", run_updates)
    update_param_card(cards / "param_card.dat", param_updates)
    print(f"configured {args.madgraph_output}")


if __name__ == "__main__":
    main()
