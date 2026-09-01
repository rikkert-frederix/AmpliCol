#!/usr/bin/env python3
"""Check inclusive HwU closure for a combined Born plus real-emission run."""

from __future__ import annotations

import pathlib
import re
import subprocess
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[1]
TAG = "integration_histogram_regression"
REPLAY_TAG = "integration_histogram_replay"
RESIDUAL_REPLAY_TAG = "integration_histogram_residual_replay"
HWU = ROOT / "Outputs" / f"{TAG}_histograms.HwU"
LOG = ROOT / "Outputs" / f"{TAG}_log_file.txt"
TAIL_LOG = ROOT / "Outputs" / f"{TAG}_tail_diagnostics.log"
TAIL_REPLAY = ROOT / "Outputs" / f"{TAG}_tail_replay.dat"
TAIL_RESIDUAL_REPLAY = ROOT / "Outputs" / f"{TAG}_tail_residual_replay.dat"
REPLAY_LOG = ROOT / "Outputs" / f"{REPLAY_TAG}_log_file.txt"
RESIDUAL_REPLAY_LOG = ROOT / "Outputs" / f"{RESIDUAL_REPLAY_TAG}_log_file.txt"
REPLAY_TAIL_LOG = ROOT / "Outputs" / f"{REPLAY_TAG}_tail_diagnostics.log"
REPLAY_TAIL_REPLAY = ROOT / "Outputs" / f"{REPLAY_TAG}_tail_replay.dat"
INVALID_CUT_TAGS = ("integration_invalid_nlo_pt", "integration_invalid_nlo_radius")
TEST_BORN_PROCESS = ROOT / "tests_nf_born_processes.txt"
TEST_REAL_PROCESS = ROOT / "tests_nf_real_processes.txt"
WJ_MIGRATION_REPLAY = ROOT / "tests" / "wj_migration_tail_replay.dat"
WJ_AUXILIARY_VECTOR_REPLAY = ROOT / "tests" / "wj_auxiliary_vector_soft_tail_replay.dat"
ALPHA01_CARD = ROOT / "tests" / "input" / "alpha01_run_card.dat"


def generate_process(process: str, destination: pathlib.Path) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(
            ["python3", str(ROOT / "process_list.py"), "--serial", process],
            cwd=tmpdir,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
        )
        destination.write_bytes((pathlib.Path(tmpdir) / "processes.txt").read_bytes())


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
    generate_process("p p > z z", TEST_BORN_PROCESS)
    generate_process("p p > z z 1j", TEST_REAL_PROCESS)
    command = [
        str(ROOT / "amplicol_generate"),
        f"--process={TEST_BORN_PROCESS.name}",
        f"--real-process={TEST_REAL_PROCESS.name}",
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
        tail_log = TAIL_LOG.read_text(encoding="utf-8")

        signed = result_line(output, "    Integral     (accum)")
        born = result_line(output, "Born")
        real_total = result_line(output, "Real - local dipoles")
        real_regular = result_line(output, "Real - local dipoles, regular")
        real_migration = result_line(output, "Real - local dipoles, migration")
        nlo_hwu = hwu_curve(hwu, "NLO")
        lo_hwu = hwu_curve(hwu, "LO")

        if real_total == (0.0, 0.0):
            raise AssertionError("all sampled real-subtraction points were lost")

        close(nlo_hwu[0], signed[0], 5e-6, "inclusive NLO central value")
        close(nlo_hwu[1], signed[1], 5e-5, "inclusive NLO uncertainty")
        close(lo_hwu[0], born[0], 5e-6, "inclusive Born central value")
        close(lo_hwu[1], born[1], 5e-5, "inclusive Born uncertainty")
        close(
            real_total[0],
            real_regular[0] + real_migration[0],
            5e-6,
            "regular plus migration real-subtraction closure",
        )
        close(
            real_total[1],
            (real_regular[1] ** 2 + real_migration[1] ** 2) ** 0.5,
            5e-5,
            "regular plus migration uncertainty closure",
        )

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
        if "# AmpliCol subtracted-real tail diagnostics v2" not in tail_log:
            raise AssertionError("missing subtracted-real tail diagnostics")
        if (
            "scores total residual component" not in tail_log
            or " residual_records" not in tail_log
            or " component_records" not in tail_log
            or " dipole" not in tail_log
            or " stratum " not in tail_log
            or "migration variance_proxy max_point_proxy fraction" not in tail_log
        ):
            raise AssertionError("tail diagnostics omit component or dipole details")
        if "Migration largest-point variance fraction:" not in output:
            raise AssertionError("combined run did not report the migration-tail convergence gate")
        if not TAIL_REPLAY.is_file():
            raise AssertionError("missing deterministic tail replay file")
        if not TAIL_RESIDUAL_REPLAY.is_file():
            raise AssertionError("missing deterministic residual-tail replay file")
        if not TAIL_REPLAY.read_text(encoding="utf-8").startswith("# AmpliCol tail replay v3"):
            raise AssertionError("tail replay fixture does not encode its integration stratum")

        for replay_file, replay_tag in (
            (TAIL_REPLAY, REPLAY_TAG),
            (TAIL_RESIDUAL_REPLAY, RESIDUAL_REPLAY_TAG),
        ):
            replay = subprocess.run(
                [
                    str(ROOT / "amplicol_generate"),
                    f"--process={TEST_BORN_PROCESS.name}",
                    f"--real-process={TEST_REAL_PROCESS.name}",
                    "--accuracy=0.9",
                    "--itmax=1",
                    "--seed=12345",
                    f"--tag={replay_tag}",
                    f"--tail-replay={replay_file}",
                ],
                cwd=ROOT,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            if "Tail replay: PASS" not in replay.stdout:
                raise AssertionError(f"saved tail point did not replay: {replay_file}")

        generate_process("p p > w+ 1j", TEST_BORN_PROCESS)
        generate_process("p p > w+ 2j", TEST_REAL_PROCESS)
        with tempfile.TemporaryDirectory() as tmpdir:
            for tag, setting in zip(INVALID_CUT_TAGS, ("pTj_min=0d0", "DRjj_min=0d0")):
                invalid_card = pathlib.Path(tmpdir) / f"{tag}.dat"
                invalid_card.write_text(f"&amplicol\n  {setting}\n/\n", encoding="utf-8")
                invalid = subprocess.run(
                    [
                        str(ROOT / "amplicol_generate"),
                        f"--process={TEST_BORN_PROCESS.name}",
                        f"--real-process={TEST_REAL_PROCESS.name}",
                        "--accuracy=0.9",
                        "--itmax=1",
                        f"--input={invalid_card}",
                        f"--tag={tag}",
                    ],
                    cwd=ROOT,
                    check=False,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                if invalid.returncode == 0 or "requires positive pTj_min and DRjj_min" not in invalid.stdout:
                    raise AssertionError(f"unsafe NLO jet cuts were accepted: {setting}")
        for replay_file in (WJ_MIGRATION_REPLAY, WJ_AUXILIARY_VECTOR_REPLAY):
            migration_replay = subprocess.run(
                [
                    str(ROOT / "amplicol_generate"),
                    f"--process={TEST_BORN_PROCESS.name}",
                    f"--real-process={TEST_REAL_PROCESS.name}",
                    "--accuracy=0.9",
                    "--itmax=1",
                    "--seed=13579",
                    f"--input={ALPHA01_CARD}",
                    f"--tag={RESIDUAL_REPLAY_TAG}",
                    f"--tail-replay={replay_file}",
                ],
                cwd=ROOT,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            if "Tail replay: PASS" not in migration_replay.stdout:
                raise AssertionError(f"saved W+j cut-migration point did not replay: {replay_file}")
            if "stratum migration" not in migration_replay.stdout:
                raise AssertionError(f"saved W+j point was not assigned to migration: {replay_file}")

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
        TAIL_LOG.unlink(missing_ok=True)
        TAIL_REPLAY.unlink(missing_ok=True)
        TAIL_RESIDUAL_REPLAY.unlink(missing_ok=True)
        REPLAY_LOG.unlink(missing_ok=True)
        RESIDUAL_REPLAY_LOG.unlink(missing_ok=True)
        REPLAY_TAIL_LOG.unlink(missing_ok=True)
        REPLAY_TAIL_REPLAY.unlink(missing_ok=True)
        for tag in INVALID_CUT_TAGS:
            (ROOT / "Outputs" / f"{tag}_log_file.txt").unlink(missing_ok=True)
        TEST_BORN_PROCESS.unlink(missing_ok=True)
        TEST_REAL_PROCESS.unlink(missing_ok=True)

    print("integration histogram physics closure: PASS")


if __name__ == "__main__":
    main()
