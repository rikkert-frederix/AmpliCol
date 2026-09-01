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
6 1
-1 1.0 -1 -2 3 -1 -2 -3
<nevents> 1 </nevents>
<seed> 1 </seed>
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


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise AssertionError(f"regression fixture does not contain exactly one {old!r}")
    return source.replace(old, new, 1)


def expect_rejected(
    reweighter: Path, input_card: Path, directory: Path, name: str, contents: str
) -> None:
    event_file = directory / f"invalid-{name}.lhe"
    event_file.write_text(contents, encoding="ascii")
    completed = subprocess.run(
        [str(reweighter), str(event_file), f"--input={input_card}"],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if completed.returncode == 0:
        raise AssertionError(f"invalid reweight input {name} was accepted")
    if completed.returncode < 0:
        raise AssertionError(
            f"invalid reweight input {name} terminated by signal "
            f"{-completed.returncode}:\n{completed.stdout}"
        )


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
        init_match = re.search(
            r"<init>\s*\n([^\n]+)\n([^\n]+)", output, flags=re.MULTILINE
        )
        if init_match is None:
            raise RuntimeError("reweighted LHE has no parseable init block")
        beam_fields = init_match.group(1).split()
        process_fields = init_match.group(2).split()
        if int(beam_fields[-2]) != -4:
            raise AssertionError("weighted output did not set IDWTUP=-4")
        if int(process_fields[-1]) != 1:
            raise AssertionError("reweighted output corrupted LPRUP")
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

        invalid_inputs = {
            "opening-tag": replace_once(
                EVENT,
                '<LesHouchesEvents version="3.0">',
                '<LesHouchesEvents version="3.0"> trailing',
            ),
            "overlong-record": replace_once(
                EVENT,
                '<LesHouchesEvents version="3.0">',
                '<LesHouchesEvents version="3.0">' + " " * 1100 + "trailing",
            ),
            "missing-seed": replace_once(EVENT, "<seed> 1 </seed>\n", ""),
            "zero-process-count": replace_once(EVENT, "6 1\n", "6 0\n"),
            "unique-row-extra": replace_once(
                EVENT,
                "-1 1.0 -1 -2 3 -1 -2 -3",
                "-1 1.0 -1 -2 3 -1 -2 -3 extra",
            ),
            "cyclic-map": replace_once(
                EVENT,
                "-1 1.0 -1 -2 3 -1 -2 -3",
                "1 1.0 -1 -2 3 -1 -2 -3",
            ),
            "duplicate-process": replace_once(
                replace_once(EVENT, "6 1\n", "6 2\n"),
                "-1 1.0 -1 -2 3 -1 -2 -3\n",
                "-1 1.0 -1 -2 3 -1 -2 -3\n"
                "-1 1.0 -1 -2 3 -1 -2 -3\n",
            ),
            "event-workspace": replace_once(
                EVENT, "<nevents> 1 </nevents>", "<nevents> 999999999 </nevents>"
            ),
            "init-process-count": replace_once(
                EVENT,
                "2212 2212 7000 7000 -1 -1 0 0 -4 1",
                "2212 2212 7000 7000 -1 -1 0 0 -4 2",
            ),
            "nan-beam": replace_once(
                EVENT,
                "2212 2212 7000 7000 -1 -1 0 0 -4 1",
                "2212 2212 NaN 7000 -1 -1 0 0 -4 1",
            ),
            "generator-trailing": replace_once(
                EVENT,
                "<generator name='AmpliCol' version='regression'></generator>",
                "<generator name='AmpliCol' version='regression'></generator> trailing",
            ),
            "generator-unclosed-opening": replace_once(
                EVENT,
                "<generator name='AmpliCol' version='regression'></generator>",
                "<generator name='AmpliCol' version='regression'</generator>",
            ),
            "nan-event-header": replace_once(
                EVENT,
                "6 1 1.0 100.0 0.0075467711 0.118",
                "6 1 NaN 100.0 0.0075467711 0.118",
            ),
            "particle-row-extra": replace_once(
                EVENT,
                "1 -1 0 0 501 0 0 0 4671.200996478833 4671.200996478833 0 0 -1",
                "1 -1 0 0 501 0 0 0 4671.200996478833 4671.200996478833 0 0 -1 extra",
            ),
            "quark-colour-slot": replace_once(
                EVENT,
                "1 -1 0 0 501 0 0 0 4671.200996478833 4671.200996478833 0 0 -1",
                "1 -1 0 0 0 501 0 0 4671.200996478833 4671.200996478833 0 0 -1",
            ),
            "unpaired-colour": replace_once(
                EVENT,
                "2 1 1 2 502 0 -1105.1428951212934",
                "2 1 1 2 504 0 -1105.1428951212934",
            ),
            "off-shell": replace_once(
                EVENT,
                "4671.200996478833 4671.200996478833 0 0 -1",
                "4671.200996478833 4672.200996478833 0 0 -1",
            ),
            "weight-mismatch": replace_once(
                EVENT, "#overwgt 1 1 1", "#overwgt 2 1 1"
            ),
            "missing-document-close": replace_once(
                EVENT, "</LesHouchesEvents>\n", ""
            ),
            "extra-event": replace_once(
                EVENT, "</LesHouchesEvents>", "<event>\n</LesHouchesEvents>"
            ),
        }
        for name, contents in invalid_inputs.items():
            expect_rejected(reweighter, input_card, Path(raw), name, contents)

        unweighted_input = Path(raw) / "event-unweighted.lhe"
        unweighted_input.write_text(EVENT, encoding="ascii")
        subprocess.run(
            [
                str(reweighter),
                str(unweighted_input),
                "--unwgt",
                f"--input={input_card}",
            ],
            check=True,
            text=True,
            capture_output=True,
        )
        unweighted_output = Path(f"{unweighted_input}.rwgt").read_text(
            encoding="ascii"
        )
        init_match = re.search(
            r"<init>\s*\n([^\n]+)\n([^\n]+)",
            unweighted_output,
            flags=re.MULTILINE,
        )
        event_match = re.search(r"<event>\s*\n([^\n]+)", unweighted_output)
        if init_match is None or event_match is None:
            raise RuntimeError("unweighted reweight output is incomplete")
        if int(init_match.group(1).split()[-2]) != -3:
            raise AssertionError("unweighted output did not set IDWTUP=-3")
        if int(init_match.group(2).split()[-1]) != 1:
            raise AssertionError("unweighted output corrupted LPRUP")
        unweighted_event_weight = float(event_match.group(1).split()[2])
        if not math.isclose(
            unweighted_event_weight, EXPECTED_REWEIGHT, rel_tol=1.0e-11
        ):
            raise AssertionError(
                f"unexpected unweighted event weight: {unweighted_event_weight}"
            )

    print("Three-quark-line LHE reweight/parser regression passed")


if __name__ == "__main__":
    main()
