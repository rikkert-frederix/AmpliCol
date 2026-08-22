#!/usr/bin/env python3
"""Exercise HEFT LHE metadata and full-colour reweighting."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import tempfile
from pathlib import Path


W_MASS = 80.419002445756163
Z_MASS = 91.188
ALPHA_EW = 0.007546771114
SW = math.sqrt(1.0 - (W_MASS / Z_MASS) ** 2)
HEFT_VEV = 2.0 * W_MASS * SW / math.sqrt(4.0 * math.pi * ALPHA_EW)


EVENT = f"""<LesHouchesEvents version="3.0">
<header>
4 0 2 1 1.0 {HEFT_VEV:.17g}
<nevents> 1 </nevents>
<seed>   31 </seed>
</header>
<init>
2212 2212 7000 7000 -1 -1 0 0 -4 1
1.0 0.0 1.0 1
<generator name='AmpliCol' version='HEFT regression'></generator>
</init>
<event>
4 1 1.0 500.0 0.007546771114 0.118
21 -1 0 0 501 502 0 0 250 250 0 0 -1
21 -1 0 0 503 501 0 0 -250 250 0 0 -1
21 1 1 2 502 503 234.375 0 0 234.375 0 0 -1
25 1 1 2 0 0 -234.375 0 0 265.625 125 0 0
#color 1 2 3 4
#overwgt 1 1 1
</event>
</LesHouchesEvents>
"""

EXPECTED_NLC = 7.0 / 9.0
EXPECTED_FULL = 8.0 / 9.0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reweighter", required=True)
    parser.add_argument("--input-card", required=True)
    args = parser.parse_args()
    reweighter = Path(args.reweighter).resolve()
    input_card = Path(args.input_card).resolve()

    with tempfile.TemporaryDirectory(prefix="amplicol-heft-reweight-") as raw:
        event_file = Path(raw) / "event.lhe"
        event_file.write_text(EVENT, encoding="ascii")
        completed = subprocess.run(
            [str(reweighter), str(event_file), f"--input={input_card}"],
            check=True,
            text=True,
            capture_output=True,
        )
        match = re.search(r"Total FC cross section:\s+(\S+)", completed.stdout)
        if match is None:
            raise RuntimeError("reweighter did not report its HEFT result")
        if not math.isclose(
            float(match.group(1)), EXPECTED_FULL, rel_tol=1.0e-12
        ):
            raise AssertionError("unexpected HEFT full-colour cross section")

        output = Path(f"{event_file}.rwgt").read_text(encoding="ascii")
        header = re.search(
            r"<header>\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)",
            output,
        )
        if header is None or header.groups()[:4] != ("4", "0", "2", "1"):
            raise AssertionError("HEFT metadata was not preserved in reweighted LHE")
        kappa, vev = (float(value) for value in header.groups()[4:])
        if not math.isclose(kappa, 1.0) or not math.isclose(
            vev, HEFT_VEV, rel_tol=1.0e-14
        ):
            raise AssertionError("HEFT Wilson-coefficient metadata changed")

        expansion = re.search(
            r"#color_expansion\s+(\S+)\s+(\S+)\s+(\S+)", output
        )
        if expansion is None:
            raise RuntimeError("reweighted HEFT LHE has no colour expansion")
        lc, nlc, full = (float(value) for value in expansion.groups())
        if not math.isclose(lc, 1.0, rel_tol=1.0e-12):
            raise AssertionError("HEFT leading-colour input weight changed")
        if not math.isclose(nlc, EXPECTED_NLC, rel_tol=1.0e-7):
            raise AssertionError(f"unexpected HEFT NLC weight: {nlc}")
        if not math.isclose(full, EXPECTED_FULL, rel_tol=1.0e-7):
            raise AssertionError(f"unexpected HEFT full-colour weight: {full}")

    print("HEFT LHE full-colour reweight regression passed")


if __name__ == "__main__":
    main()
