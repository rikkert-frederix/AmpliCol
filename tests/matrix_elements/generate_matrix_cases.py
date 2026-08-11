#!/usr/bin/env python
"""Generate deterministic matrix-element regression cases.

The process expansion stays delegated to process_list.py. This script only
collects the generated process/order rows into a compact fixture format that
the Fortran regression harness can read.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path


FAMILIES = (
    ("pp_2j", None, ("-FS", "4", "p p > 2j")),
    ("pp_3j", None, ("-FS", "3", "p p > 3j")),
    ("pp_4j", None, ("-FS", "2", "p p > 4j")),
    ("pp_a_2j", None, ("-FS", "4", "p p > a 2j")),
    ("pp_z_2j", None, ("-FS", "4", "p p > z 2j")),
    ("pp_wp_2j", None, ("-FS", "4", "p p > w+ 2j")),
    ("pp_az_2j", None, ("-FS", "4", "p p > a z 2j")),
    ("pp_zz_2j", None, ("-FS", "4", "p p > z z 2j")),
    ("pp_wpz_2j", None, ("-FS", "4", "p p > w+ z 2j")),
    ("pp_awp_2j", None, ("-FS", "4", "p p > a w+ 2j")),
    ("pp_zwp_2j", None, ("-FS", "4", "p p > z w+ 2j")),
    ("pp_wpwp_2j", None, ("-FS", "4", "p p > w+ w+ 2j")),
    ("pp_wpwm_2j", None, ("-FS", "4", "p p > w+ w- 2j")),
    ("pp_az_3j", None, ("-FS", "4", "p p > a z 3j")),
    ("pp_zz_3j", None, ("-FS", "4", "p p > z z 3j")),
    ("pp_wpz_3j", None, ("-FS", "4", "p p > w+ z 3j")),
    ("pp_awp_3j", None, ("-FS", "4", "p p > a w+ 3j")),
    ("pp_zwp_3j", None, ("-FS", "4", "p p > z w+ 3j")),
    ("pp_wpwp_3j", None, ("-FS", "4", "p p > w+ w+ 3j")),
    ("pp_wpwm_3j", None, ("-FS", "4", "p p > w+ w- 3j")),
    ("pp_ttbar_0j", None, ("-FS", "3", "p p > t t~")),
    ("pp_ttbar_1j", None, ("-FS", "3", "p p > t t~ 1j")),
    ("pp_ttbar_2j", None, ("-FS", "2", "p p > t t~ 2j")),
    ("pp_ttbar_3j", None, ("-FS", "1", "p p > t t~ 3j")),
    ("pp_ttbar_4j", None, ("-FS", "1", "p p > t t~ 4j")),
    ("pp_aa_0j", None, ("-FS", "3", "p p > a a")),
    ("pp_aa_1j", None, ("-FS", "3", "p p > a a 1j")),
    ("pp_aa_2j", None, ("-FS", "2", "p p > a a 2j")),
    ("pp_aa_3j", None, ("-FS", "2", "p p > a a 3j")),
    ("pp_aa_4j", None, ("-FS", "2", "p p > a a 4j")),
    ("pp_aaa_0j", None, ("-FS", "2", "p p > a a a")),
    ("pp_aaa_1j", None, ("-FS", "2", "p p > a a a 1j")),
    ("pp_aaa_2j", None, ("-FS", "2", "p p > a a a 2j")),
    ("pp_aaa_3j", 80, ("-FS", "2", "p p > a a a 3j")),
    ("pp_aaa_4j", 40, ("-FS", "2", "p p > a a a 4j")),
    ("pp_ttbar_a_0j", None, ("-FS", "2", "p p > t t~ a")),
    ("pp_ttbar_a_1j", None, ("-FS", "2", "p p > t t~ a 1j")),
    ("pp_ttbar_a_2j", None, ("-FS", "2", "p p > t t~ a 2j")),
    ("pp_ttbar_a_3j", 120, ("-FS", "2", "p p > t t~ a 3j")),
    ("pp_ttbar_a_4j", 60, ("-FS", "2", "p p > t t~ a 4j")),
    ("pp_ttbar_aa_0j", None, ("-FS", "1", "p p > t t~ a a")),
    ("pp_ttbar_aa_1j", None, ("-FS", "1", "p p > t t~ a a 1j")),
    ("pp_ttbar_aa_2j", 80, ("-FS", "1", "p p > t t~ a a 2j")),
    ("pp_ttbar_aa_3j", 40, ("-FS", "1", "p p > t t~ a a 3j")),
    ("pp_ttbar_aa_4j", 20, ("-FS", "1", "p p > t t~ a a 4j")),
    # Electroweak stress cases.  In particular, the three- and four-boson
    # families exercise every auxiliary-field decomposition of a four-point
    # vector/Higgs interaction before the FD-gauge implementation is enabled.
    ("pp_zzz_0j", None, ("-FS", "2", "p p > z z z")),
    ("pp_wpwmz_0j", None, ("-FS", "2", "p p > w+ w- z")),
    ("pp_wpwma_0j", None, ("-FS", "2", "p p > w+ w- a")),
    ("pp_wpwmh_0j", None, ("-FS", "2", "p p > w+ w- h")),
    ("pp_zhh_0j", None, ("-FS", "2", "p p > z h h")),
    ("pp_zzzz_0j", None, ("-FS", "2", "p p > z z z z")),
    ("pp_wpwmwpwm_0j", None, ("-FS", "2", "p p > w+ w- w+ w-")),
    ("pp_ttbar_z_0j", None, ("-FS", "2", "p p > t t~ z")),
    ("pp_ttbar_h_0j", None, ("-FS", "2", "p p > t t~ h")),
    ("pp_tbbar_wm_0j", None, ("-FS", "2", "p p > t b~ w-")),
    ("pp_tbarb_wp_0j", None, ("-FS", "2", "p p > t~ b w+")),
)


def parse_process_file(path: Path, family: str) -> list[tuple[str, int, int, int, str, list[int], list[int]]]:
    lines = path.read_text().splitlines()
    n_external, n_unique = map(int, lines[0].split()[:2])

    idx = 1 + n_unique
    while idx < len(lines) and not lines[idx].strip():
        idx += 1
    n_groups = int(lines[idx].split()[0])
    idx += 1

    cases = []
    for _ in range(n_groups):
        while idx < len(lines) and not lines[idx].strip():
            idx += 1
        header = lines[idx].split()
        group_id = int(header[0])
        n_rows = int(header[1])
        idx += 1

        for row_id in range(1, n_rows + 1):
            parts = lines[idx].split()
            idx += 1
            n_channels = int(parts[0])
            proc_start = 1 + n_channels
            process = [int(x) for x in parts[proc_start : proc_start + n_external]]
            order = [int(x) for x in parts[proc_start + n_external : proc_start + 2 * n_external]]
            order = direct_amplitude_order(process, order)
            point = point_name(n_external, process)
            cases.append((family, n_external, group_id, row_id, point, process, order))

    return cases


def direct_amplitude_order(process: list[int], order: list[int]) -> list[int]:
    """Move colour singlets before the closing coloured leg.

    process_list.py writes singlets after the anti-quark in the colour string.
    amplitude_QCD%init uses the last entry to close the current and requires it
    to be coloured, so direct matrix-element tests need the equivalent local
    convention used by handling_processes.setup_color_order.
    """
    singlets = [idx for idx in order if is_singlet(process[idx - 1])]
    if not singlets:
        return order
    coloured = [idx for idx in order if not is_singlet(process[idx - 1])]
    return coloured[:-1] + singlets + coloured[-1:]


def is_singlet(pdg: int) -> bool:
    return not (abs(pdg) <= 6 or pdg == 21)


def point_name(n_external: int, process: list[int]) -> str:
    return "generic"


def select_cases(
    cases: list[tuple[str, int, int, int, str, list[int], list[int]]],
    limit: int | None,
) -> list[tuple[str, int, int, int, str, list[int], list[int]]]:
    if limit is None or len(cases) <= limit:
        return cases
    if limit <= 1:
        return cases[:limit]
    indices = sorted({round(i * (len(cases) - 1) / (limit - 1)) for i in range(limit)})
    return [cases[i] for i in indices]


def generate_cases(repo_root: Path) -> list[tuple[str, int, int, int, str, list[int], list[int]]]:
    process_list = repo_root / "process_list.py"
    all_cases = []
    with tempfile.TemporaryDirectory(prefix="weyl-me-cases-") as tmp:
        tmp_path = Path(tmp)
        for family, limit, args in FAMILIES:
            workdir = tmp_path / family
            workdir.mkdir()
            command = [sys.executable, str(process_list), "--serial", *args]
            subprocess.run(command, cwd=workdir, check=True)
            family_cases = parse_process_file(workdir / "processes.txt", family)
            family_cases = select_cases(family_cases, limit)
            all_cases.extend(family_cases)
    return all_cases


def write_cases(path: Path, cases: list[tuple[str, int, int, int, str, list[int], list[int]]]) -> None:
    with path.open("w") as handle:
        handle.write(f"1 {len(cases)}\n")
        for case_id, (family, n_external, group_id, row_id, point, process, order) in enumerate(cases, 1):
            handle.write(f"{case_id} {family} {n_external} {group_id} {row_id} {point}\n")
            handle.write(" ".join(str(x) for x in process) + "\n")
            handle.write(" ".join(str(x) for x in order) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/matrix_elements/cases.dat"),
        help="Path for the generated case fixture.",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Repository root containing process_list.py.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output = args.output if args.output.is_absolute() else repo_root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    cases = generate_cases(repo_root)
    write_cases(output, cases)
    print(f"wrote {len(cases)} matrix-element cases to {output}")


if __name__ == "__main__":
    main()
