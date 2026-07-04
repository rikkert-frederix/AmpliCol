from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

from pyamplicol.reference import (
    AmplicolAdapter,
    CommandResult,
    parse_amplicol_probe_points,
    parse_first_phase_space_point,
    parse_first_matrix_element,
    parse_timing_rows,
)


def test_reference_adapter_prepares_library_with_expected_command_order(
    tmp_path: Path,
) -> None:
    calls: list[tuple[str, ...]] = []

    def runner(
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult:
        calls.append(tuple(args))
        return CommandResult(tuple(args), cwd, 0, "", "", 0.01)

    adapter = AmplicolAdapter(tmp_path, runner=runner, jobs=3)
    result = adapter.prepare_library("d d~ > z g g")

    assert result.process_file == tmp_path / "processes.txt"
    assert calls == [
        ("make", "cleanlib"),
        ("make", "-j3", "amplicol_generate"),
        (
            "./amplicol_generate",
            "--library=create",
            f"--process={tmp_path / 'processes.txt'}",
        ),
        ("make", "-j3", "amplicol_generate_library"),
    ]
    assert (tmp_path / "processes.txt").exists()


def test_reference_adapter_me_test_automatically_enables_matching_timing_sample(
    tmp_path: Path,
) -> None:
    calls: list[tuple[tuple[str, ...], Mapping[str, str] | None]] = []

    def runner(
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult:
        calls.append((tuple(args), env))
        stdout = "  AmpliCol matrix element:   1.82022213266464E-05\n"
        return CommandResult(tuple(args), cwd, 0, stdout, "", 0.02)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    result = adapter.run_me_test(
        "d d~ > z g g",
        points=7,
        mg5_path=tmp_path / "MG5",
    )

    assert calls == [
        (
            (
                "./amplicol_generate",
                "--me_test=7",
                "--timing=7",
                f"--process={tmp_path / 'processes.txt'}",
            ),
            {"MG5_PATH": str(tmp_path / "MG5")},
        )
    ]
    assert result.first_point_matrix_element == 1.82022213266464e-05


def test_reference_adapter_amplicol_probe_uses_direct_probe_flag(tmp_path: Path) -> None:
    calls: list[tuple[str, ...]] = []

    def runner(
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult:
        calls.append(tuple(args))
        stdout = """
AMPICOL_PROBE_VALUE 1 1 1   1.42000000000000E+02
AMPICOL_PROBE_MOM 1 1 1   5.00000000000000E+02   0.00000000000000E+00   0.00000000000000E+00   5.00000000000000E+02
AMPICOL_PROBE_MOM 1 2 -1   5.00000000000000E+02   0.00000000000000E+00   0.00000000000000E+00  -5.00000000000000E+02
AMPICOL_PROBE_MOM 1 3 23   1.00000000000000E+03   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00
 AmpliCol probe first phase-space point
   group, integral:           1           1
     i      pdg                      E                     px                     py                     pz
     1        1   5.00000000000000E+02   0.00000000000000E+00   0.00000000000000E+00   5.00000000000000E+02
     2       -1   5.00000000000000E+02   0.00000000000000E+00   0.00000000000000E+00  -5.00000000000000E+02
     3       23   1.00000000000000E+03   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00
  AmpliCol matrix element:   1.42000000000000E+02
"""
        return CommandResult(tuple(args), cwd, 0, stdout, "", 0.02)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    result = adapter.run_amplicol_probe("d d~ > z", points=3)

    assert calls == [
        (
            "./amplicol_generate",
            "--amplicol_probe=3",
            "--timing=3",
            f"--process={tmp_path / 'processes.txt'}",
        )
    ]
    assert result.first_phase_space_point is not None
    assert result.first_phase_space_point.matrix_element == 142.0
    assert len(result.probe_points) == 1
    assert result.probe_points[0].matrix_element == 142.0


def test_reference_adapter_fixed_amplicol_probe_uses_fixed_probe_flag(
    tmp_path: Path,
) -> None:
    calls: list[tuple[str, ...]] = []

    def runner(
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult:
        calls.append(tuple(args))
        stdout = """
AMPICOL_PROBE_VALUE 1 1 1   1.4210653372872991E+02
AMPICOL_PROBE_MOM 1 1 1   4.55940000000000E+01   0.0   0.0   4.55940000000000E+01
AMPICOL_PROBE_MOM 1 2 -1  4.55940000000000E+01   0.0   0.0  -4.55940000000000E+01
AMPICOL_PROBE_MOM 1 3 23  9.11880000000000E+01   0.0   0.0   0.0
"""
        return CommandResult(tuple(args), cwd, 0, stdout, "", 0.02)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    result = adapter.run_amplicol_fixed_probe("d d~ > z", points=4, timing_sample=2)

    assert calls == [
        (
            "./amplicol_generate",
            "--amplicol_fixed_probe=4",
            "--timing=2",
            f"--process={tmp_path / 'processes.txt'}",
        )
    ]
    assert len(result.probe_points) == 1
    assert result.probe_points[0].matrix_element == 142.10653372872991


def test_reference_adapter_quiet_probe_adds_legacy_quiet_flag(tmp_path: Path) -> None:
    calls: list[tuple[str, ...]] = []

    def runner(
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult:
        calls.append(tuple(args))
        stdout = "Timing summary                           seconds    percent  note\n"
        return CommandResult(tuple(args), cwd, 0, stdout, "", 0.02)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    adapter.run_amplicol_probe("d d~ > z g", points=5, quiet=True)

    assert calls == [
        (
            "./amplicol_generate",
            "--amplicol_probe=5",
            "--timing=5",
            f"--process={tmp_path / 'processes.txt'}",
            "--amplicol_probe_quiet",
        )
    ]


def test_parse_first_matrix_element_returns_none_when_missing() -> None:
    assert parse_first_matrix_element("no matrix element here") is None


def test_parse_timing_rows_uses_only_timing_summary_table() -> None:
    output = """
  AmpliCol matrix element:   1.62597906589508E+00
------------------------------------------------------------------------------
Timing summary                           seconds    percent  note
------------------------------------------------------------------------------
phase-space initialisation              0.000015      0.20%
amplitude evaluation                    0.000006      0.08%  sampled
total                                   0.007684    100.00%
------------------------------------------------------------------------------
"""

    rows = parse_timing_rows(output)

    assert [row.label for row in rows] == [
        "phase-space initialisation",
        "amplitude evaluation",
        "total",
    ]
    assert rows[0].seconds == 0.000015
    assert rows[1].note == "sampled"


def test_parse_amplicol_probe_points_extracts_all_compact_records() -> None:
    output = """
AMPICOL_PROBE_VALUE 2 1 1   2.50000000000000E+00
AMPICOL_PROBE_MOM 2 2 -1   2.00000000000000E+00 0.0 0.0 -2.0
AMPICOL_PROBE_MOM 2 1 1    2.00000000000000E+00 0.0 0.0  2.0
AMPICOL_PROBE_VALUE 1 1 1   1.50000000000000E+00
AMPICOL_PROBE_MOM 1 1 1    1.00000000000000E+00 0.0 0.0  1.0
AMPICOL_PROBE_MOM 1 2 -1   1.00000000000000E+00 0.0 0.0 -1.0
"""

    points = parse_amplicol_probe_points(output)

    assert [point.point for point in points] == [1, 2]
    assert points[0].matrix_element == 1.5
    assert [particle.pdg for particle in points[1].particles] == [1, -1]


def test_parse_first_phase_space_point_extracts_kinematics_and_me() -> None:
    output = """
 ME-check first phase-space point
   group, integral:           1           1
     i      pdg                      E                     px                     py                     pz
     1        1   5.00000000000000E+02   0.00000000000000E+00   0.00000000000000E+00   5.00000000000000E+02
     2       -1   5.00000000000000E+02   0.00000000000000E+00   0.00000000000000E+00  -5.00000000000000E+02
     3       23   1.00000000000000E+03   0.00000000000000E+00   0.00000000000000E+00   0.00000000000000E+00
  AmpliCol matrix element:   1.42000000000000E+02
"""

    point = parse_first_phase_space_point(output)

    assert point is not None
    assert point.group == 1
    assert point.integral == 1
    assert point.matrix_element == 142.0
    assert [particle.pdg for particle in point.particles] == [1, -1, 23]
    assert point.particles[0].momentum == (500.0, 0.0, 0.0, 500.0)
