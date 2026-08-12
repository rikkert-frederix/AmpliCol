#!/usr/bin/env python3
"""Exercise coupling-order metadata through two consecutive LHE reads."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile


EVENT = re.compile(r"^\s*<event>\s*$", re.MULTILINE)
DECLARED_EVENTS = re.compile(r"<nevents>\s*(\d+)\s*</nevents>")

FIXTURE_HEADER = """<LesHouchesEvents version="3.0">
<header>
           6           4
          -1   1.0000000000000000               -1          -2          -1          -2          21          21
          -1   1.0000000000000000               -1          -2          -2          -1          21          21
           2   1.0000000000000000               -2          -1          -1          -2          21          21
           1   1.0000000000000000               -2          -1          -2          -1          21          21
<nevents>           30 </nevents>
<seed>          970300 </seed>
</header>
<init>
   2212   2212 0.70000000E+04 0.70000000E+04 -1 -1   244800   244800 -3   1
 0.16540831E+03 0.81790887E+01 0.16540831E+03      1
<generator name='AmpliCol' version='1.0'>please cite arXiv:2601.19483</generator>
</init>
"""

FIXTURE_EVENT = """<event>
  6      1 0.16540831E+03 0.91188000E+02 0.75467711E-02 0.11900000E+00
        2 -1    0    0  503    0  0.00000000000000000E+00  0.00000000000000000E+00  0.68228618054761409E+03  0.68228618054761409E+03  0.00000000000000000E+00 0.0000E+00 -.1000E+01
        1 -1    0    0  504    0  0.00000000000000000E+00  0.00000000000000000E+00 -0.16112480598564357E+04  0.16112480598564357E+04  0.00000000000000000E+00 0.0000E+00 -.1000E+01
        1  1    1    2  504    0 -0.25551012650678423E+02 -0.43930923015685686E+02 -0.15986667139721549E+04  0.15994743019526800E+04  0.00000000000000000E+00 0.0000E+00 -.1000E+01
        2  1    1    2  501    0  0.52475825344633265E+02  0.10744569147905978E+02  0.27891142395702150E+03  0.28400834569693166E+03  0.00000000000000000E+00 0.0000E+00 -.1000E+01
       21  1    1    2  502  501 -0.36716031804542560E+02 -0.30435832176319082E+02  0.84149505606467727E+02  0.96724072316823907E+02  0.00000000000000000E+00 0.0000E+00 -.1000E+01
       21  1    1    2  503  502  0.97912191105877184E+01  0.63622186044098790E+02  0.30664390509984418E+03  0.31332752043761457E+03  0.00000000000000000E+00 0.0000E+00 -.1000E+01
#color  4  5  6  1  3  2
#overwgt 0.16540831E+03 0.16540831E+03 0.00000000E+00
</event>
"""

PURE_EW_LEGACY = """<LesHouchesEvents version="3.0">
<header>
4 1
-1 1.0 -2 2 24 -24
<nevents> 1 </nevents>
</header>
<init>
2212 2212 6500.0 6500.0 -1 -1 0 0 -3 1
1.0 0.0 1.0 1
<generator name='audit'>legacy pure-EW smoke</generator>
</init>
<event>
4 1 1.0 200.0 0.0075467711 0.118
2 -1 0 0 501 0 0.0 0.0 100.0 100.0 0.0 0.0 -1.0
-2 -1 0 0 0 501 0.0 0.0 -100.0 100.0 0.0 0.0 1.0
24 1 1 2 0 0 0.0 59.496638655 0.0 100.0 80.379 0.0 -1.0
-24 1 1 2 0 0 0.0 -59.496638655 0.0 100.0 80.379 0.0 1.0
#color 2 3 4 1
</event>
</LesHouchesEvents>
"""


def run_reweighter(command: list[str], cwd: Path) -> None:
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            f"command failed with exit status {completed.returncode}:\n"
            f"{' '.join(command)}\n{completed.stdout}"
        )


def check_lhe(path: Path, resolved_as2: int) -> int:
    if not path.is_file():
        raise RuntimeError(f"reweighter did not create {path}")
    contents = path.read_text(encoding="utf-8")
    coupling_tag = re.compile(
        rf"<coupling_orders>\s*0\s+{resolved_as2}\s+"
        r"-1\s+-1\s+-1\s+-1\s*</coupling_orders>"
    )
    if not coupling_tag.search(contents):
        raise RuntimeError(f"{path} is missing the expected coupling-order tag")
    match = DECLARED_EVENTS.search(contents)
    if match is None:
        raise RuntimeError(f"{path} is missing its declared event count")
    declared = int(match.group(1))
    actual = len(EVENT.findall(contents))
    if declared <= 0 or declared != actual:
        raise RuntimeError(
            f"{path} declares {declared} events but contains {actual}"
        )
    return actual


def check_weight_strategy(path: Path) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        init_line = lines.index("<init>")
    except ValueError as error:
        raise RuntimeError(f"{path} is missing its init block") from error
    beam = lines[init_line + 1].split()
    process = lines[init_line + 2].split()
    if beam[-2:] != ["-4", "1"] or process[-1] != "1":
        raise RuntimeError(
            f"{path} did not retain IDWTUP=-4, NPRUP=1, LPRUP=1"
        )


def write_fixture(path: Path) -> None:
    path.write_text(
        FIXTURE_HEADER + FIXTURE_EVENT * 30 + "</LesHouchesEvents>\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reweighter", type=Path, required=True)
    parser.add_argument("--input-card", type=Path, required=True)
    args = parser.parse_args()

    reweighter = args.reweighter.resolve()
    input_card = args.input_card.resolve()
    with tempfile.TemporaryDirectory(prefix="amplicol-coupling-lhe-") as temp:
        workdir = Path(temp)
        first_input = workdir / "events.lhe"
        write_fixture(first_input)

        run_reweighter(
            [
                str(reweighter),
                str(first_input),
                "--unwgt",
                f"--input={input_card}",
            ],
            workdir,
        )
        first_output = Path(f"{first_input}.rwgt")
        first_count = check_lhe(first_output, 8)

        run_reweighter(
            [str(reweighter), str(first_output), f"--input={input_card}"],
            workdir,
        )
        second_output = Path(f"{first_output}.rwgt")
        second_count = check_lhe(second_output, 8)
        if second_count != first_count:
            raise RuntimeError(
                "the second coupling-order read changed the event count: "
                f"{first_count} -> {second_count}"
            )

        legacy_input = workdir / "pure_ew_legacy.lhe"
        legacy_input.write_text(PURE_EW_LEGACY, encoding="utf-8")
        run_reweighter(
            [str(reweighter), str(legacy_input), f"--input={input_card}"],
            workdir,
        )
        legacy_output = Path(f"{legacy_input}.rwgt")
        if check_lhe(legacy_output, 0) != 1:
            raise RuntimeError("legacy pure-EW inference changed its event count")
        check_weight_strategy(legacy_output)

    print("Coupling-order LHE roundtrip regression passed")


if __name__ == "__main__":
    main()
