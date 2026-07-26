#!/usr/bin/env python3
"""Regression checks for the standalone HwU histogram converter."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPT = REPOSITORY / "Utilities" / "plot_events" / "internal" / "histograms.py"

HWU_DATA = r"""##& xmin & xmax & central value & dy & muR=0.5 muF=0.5 & muR=2.0 muF=2.0

<histogram> 4 "Transverse momentum |TYPE@NLO |X_AXIS@LIN |Y_AXIS@LOG"
 +0.0000000e+00 +1.0000000e+00 +1.0000000e+01 +1.0000000e+00 +9.0000000e+00 +1.1000000e+01
 +1.0000000e+00 +2.0000000e+00 +2.0000000e+01 +2.0000000e+00 +1.8000000e+01 +2.2000000e+01
 +2.0000000e+00 +3.0000000e+00 +3.0000000e+01 +3.0000000e+00 +2.7000000e+01 +3.3000000e+01
 +3.0000000e+00 +4.0000000e+00 +4.0000000e+01 +4.0000000e+00 +3.6000000e+01 +4.4000000e+01
<\histogram>
"""


class HistogramCliTest(unittest.TestCase):
    def run_converter(self, *arguments):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), *map(str, arguments)],
            cwd=self.workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if result.returncode:
            self.fail(
                f"converter exited with {result.returncode}\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        return result

    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.workdir = Path(self.temporary_directory.name)
        self.input_a = self.workdir / "input_a.HwU"
        self.input_b = self.workdir / "input_b.HwU"
        self.input_a.write_text(HWU_DATA)
        self.input_b.write_text(HWU_DATA)

    def test_help_has_no_package_dependency(self):
        result = subprocess.run(
            [str(SCRIPT), "--help"],
            cwd=self.workdir,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        self.assertIn("python histograms.py", result.stdout)
        self.assertNotIn("ModuleNotFoundError", result.stderr)

    def test_gnuplot_and_hwu_outputs(self):
        output = self.workdir / "converted"
        self.run_converter(self.input_a, "--gnuplot", "--no_open", f"--out={output}")
        self.assertTrue(output.with_suffix(".HwU").is_file())
        self.assertTrue(output.with_suffix(".gnuplot").is_file())
        self.assertIn("Transverse momentum", output.with_suffix(".HwU").read_text())
        self.assertIn('set output "converted.ps"', output.with_suffix(".gnuplot").read_text())

    def test_raw_hwu_round_trip(self):
        output = self.workdir / "roundtrip"
        self.run_converter(self.input_a, "--HwU", f"--out={output}")
        content = output.with_suffix(".HwU").read_text()
        self.assertIn("central value", content)
        self.assertIn("<histogram> 4", content)

    def test_multi_file_sum_multiply_and_rebin(self):
        output = self.workdir / "combined"
        self.run_converter(
            self.input_a,
            self.input_b,
            "--sum",
            "--multiply=1,2",
            "--rebin=2",
            "--no_stat",
            "--no_open",
            f"--out={output}",
        )
        content = output.with_suffix(".HwU").read_text()
        self.assertIn("<histogram> 2", content)
        self.assertIn("+9.0000000e+01", content)


if __name__ == "__main__":
    unittest.main()
