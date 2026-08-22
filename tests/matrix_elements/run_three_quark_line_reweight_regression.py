#!/usr/bin/env python3
"""Exercise the LC-LHE to full-colour three-quark-line reweight path."""

from __future__ import annotations

import argparse
import math
import re
import subprocess
import tempfile
from pathlib import Path


EVENT = """<LesHouchesEvents version="3.0">
<header>
<flavour_scheme> 4 </flavour_scheme>
6 0
<nevents> 1 </nevents>
</header>
<init>
2212 2212 7000 7000 -1 -1 0 0 -4 1
1.0 0.0 1.0 1
<generator name='AmpliCol' version='regression'></generator>
</init>
<event>
6 1 1.0 100.0 0.0075467711 0.118
1 -1 0 0 501 0 0 0 4671.200996478833 4671.200996478833 0 0 -1
-1 -1 0 0 0 501 0 0 -5452.624496459750 5452.624496459750 0 0 1
2 1 1 2 502 0 -1105.1428951212934 -1737.2006769281086 -3251.7412381317486 3848.769685279069 0 0 -1
-2 1 1 2 0 502 2228.296688374595 1630.5712325370819 2402.6278548445193 3660.1488063565753 0 0 1
3 1 1 2 503 0 -470.2433927087442 -436.62760397349416 1102.8105076070963 1275.9167404758375 0 0 -1
-3 1 1 2 0 503 -652.9104005445573 543.2570483645209 -1035.1206243007837 1338.9902608271016 0 0 1
#color 1 2 3 4 5 6
#overwgt 1 1 1
</event>
</LesHouchesEvents>
"""

EXPECTED_REWEIGHT = 0.4485971947600595


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reweighter", required=True)
    parser.add_argument("--input-card", required=True)
    args = parser.parse_args()
    reweighter = Path(args.reweighter).resolve()
    input_card = Path(args.input_card).resolve()

    with tempfile.TemporaryDirectory(prefix="amplicol-three-line-reweight-") as raw:
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
            raise RuntimeError("reweighter did not report its full-colour result")
        reported = float(match.group(1))
        if not math.isclose(reported, EXPECTED_REWEIGHT, rel_tol=1.0e-11):
            raise AssertionError(
                f"full-colour reweight mismatch: {reported} != {EXPECTED_REWEIGHT}"
            )

        output = Path(f"{event_file}.rwgt").read_text(encoding="ascii")
        if "<flavour_scheme>  4 </flavour_scheme>" not in output and \
                "<flavour_scheme> 4 </flavour_scheme>" not in output:
            raise AssertionError("reweighted LHE did not preserve flavour-scheme metadata")
        match = re.search(
            r"#color_expansion\s+(\S+)\s+(\S+)\s+(\S+)", output
        )
        if match is None:
            raise RuntimeError("reweighted LHE has no colour-expansion record")
        lc, nlc, full = (float(value) for value in match.groups())
        if not math.isclose(lc, 1.0, rel_tol=1.0e-12):
            raise AssertionError(f"input LC weight changed unexpectedly: {lc}")
        if not math.isclose(nlc, EXPECTED_REWEIGHT, rel_tol=1.0e-7) or not math.isclose(
            full, EXPECTED_REWEIGHT, rel_tol=1.0e-7
        ):
            raise AssertionError(f"unexpected LHE colour expansion: {(lc, nlc, full)}")

    print("Three-quark-line LHE reweight regression passed")


if __name__ == "__main__":
    main()
