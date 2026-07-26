#!/usr/bin/env python3
"""Check inclusive HwU closure for a combined Born plus real-emission run."""

from __future__ import annotations

import pathlib
import re
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
TAG = "integration_histogram_regression"
HWU = ROOT / "Outputs" / f"{TAG}_histograms.HwU"
LOG = ROOT / "Outputs" / f"{TAG}_log_file.txt"


def close(actual: float, expected: float, tolerance: float, label: str) -> None:
    if abs(actual - expected) > tolerance * max(1.0, abs(expected)):
        raise AssertionError(f"{label}: {actual} != {expected}")


def result_line(text: str, label: str) -> tuple[float, float]:
    match = re.search(
        rf"^{re.escape(label)}:\s+([+-]?\d+\.\d+E[+-]\d+)\s+\+/-\s+([+-]?\d+\.\d+E[+-]\d+)",
        text,
        re.MULTILINE,
    )
    if not match:
        raise AssertionError(f"missing result line: {label}")
    return float(match.group(1)), float(match.group(2))


def hwu_curve(text: str, curve: str) -> tuple[float, float]:
    match = re.search(
        rf'<histogram>\s+\d+\s+".*?\|T@{curve}"\s*\n'
        rf"\s*[+-]?\d+\.\d+E[+-]\d+\s+[+-]?\d+\.\d+E[+-]\d+\s+"
        rf"([+-]?\d+\.\d+E[+-]\d+)\s+([+-]?\d+\.\d+E[+-]\d+)",
        text,
    )
    if not match:
        raise AssertionError(f"missing HwU {curve} curve")
    return float(match.group(1)), float(match.group(2))


def main() -> None:
    command = [
        str(ROOT / "amplicol_generate"),
        "--process=processes_zz.txt",
        "--real-process=processes_zzj.txt",
        "--accuracy=0.9",
        "--itmax=1",
        "--seed=12345",
        f"--tag={TAG}",
    ]
    try:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        output = completed.stdout
        hwu = HWU.read_text(encoding="utf-8")

        signed = result_line(output, "    Integral     (accum)")
        born = result_line(output, "Born")
        nlo_hwu = hwu_curve(hwu, "NLO")
        lo_hwu = hwu_curve(hwu, "LO")

        close(nlo_hwu[0], signed[0], 5e-6, "inclusive NLO central value")
        close(nlo_hwu[1], signed[1], 5e-5, "inclusive NLO uncertainty")
        close(lo_hwu[0], born[0], 5e-6, "inclusive Born central value")
        close(lo_hwu[1], born[1], 5e-5, "inclusive Born uncertainty")

        required_histograms = [
            "selected jet multiplicity",
            "jet 1 pT [GeV]",
            "jet 6 eta",
            "jet HT [GeV]",
            "leading dijet invariant mass [GeV]",
            "charged-lepton multiplicity",
            "leading dilepton invariant mass [GeV]",
            "photon multiplicity",
            "missing transverse momentum [GeV]",
            "top/W/Z/H multiplicity",
        ]
        for title in required_histograms:
            if f'" {title} |T@NLO"' not in hwu or f'" {title} |T@LO"' not in hwu:
                raise AssertionError(f"missing default histogram: {title}")
        if hwu.count("<histogram>") != 66:
            raise AssertionError("unexpected number of default LO/NLO histogram blocks")

        with tempfile.TemporaryDirectory() as tmpdir:
            subprocess.run(
                [
                    "python3",
                    str(ROOT / "Utilities/plot_events/internal/histograms.py"),
                    str(HWU),
                    "--gnuplot",
                    f"--out={pathlib.Path(tmpdir) / 'default_analysis'}",
                    "--no_open",
                ],
                cwd=ROOT,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.STDOUT,
            )
    finally:
        HWU.unlink(missing_ok=True)
        LOG.unlink(missing_ok=True)

    print("integration histogram physics closure: PASS")


if __name__ == "__main__":
    main()
