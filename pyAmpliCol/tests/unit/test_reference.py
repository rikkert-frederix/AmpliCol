from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import pytest

import pyamplicol.reference as reference_module
from pyamplicol.core_types import ExternalMomentum
from pyamplicol.reference import (
    AmplicolAdapter,
    AmplicolProbePoint,
    AmplicolWorkflowResult,
    CommandResult,
    amplicol_process_file_integrals,
    parse_amplicol_probe_points,
    parse_first_phase_space_point,
    parse_first_matrix_element,
    parse_timing_rows,
    reference_color_order_for_run,
    reorder_external_momenta_by_pdg,
)


def test_popen_runner_terminates_process_group_on_keyboard_interrupt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    terminated: list[tuple[int, object]] = []

    class FakeProcess:
        pid = 1234
        returncode = None
        calls = 0

        def communicate(self, timeout: float | None = None):  # noqa: ANN001
            self.calls += 1
            if self.calls == 1:
                raise KeyboardInterrupt
            return "", ""

    monkeypatch.setattr(reference_module.subprocess, "Popen", lambda *args, **kwargs: FakeProcess())
    monkeypatch.setattr(
        reference_module,
        "_terminate_process_group",
        lambda pid, sig: terminated.append((pid, sig)),
    )

    with pytest.raises(KeyboardInterrupt):
        reference_module.popen_runner(["make"], cwd=tmp_path)

    assert terminated == [(1234, reference_module.signal.SIGTERM)]


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


def test_reference_adapter_can_warm_library_creation_with_supplied_momenta(
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

    particles = (
        ExternalMomentum(1, (500.0, 0.0, 0.0, 500.0)),
        ExternalMomentum(-1, (500.0, 0.0, 0.0, -500.0)),
        ExternalMomentum(23, (1000.0, 0.0, 0.0, 0.0)),
    )
    adapter = AmplicolAdapter(tmp_path, runner=runner, jobs=3)

    adapter.prepare_library(
        "d d~ > z",
        warmup_particles=particles,
        warmup_points=10,
    )

    assert calls[2] == (
        "./amplicol_generate",
        "--library=create",
        f"--process={tmp_path / 'processes.txt'}",
        "--amplicol_momenta_probe=10",
        "--amplicol_probe_quiet",
    )
    assert (
        tmp_path / "Utilities" / "ME_checks" / "momenta_1_1.txt"
    ).exists()


def test_process_file_integral_parser_lists_group_rows(tmp_path: Path) -> None:
    path = tmp_path / "processes.txt"
    path.write_text(
        "\n".join(
            [
                "4 1",
                "1 -1 23 21",
                "2",
                "1 2",
                "1 1 1 -1 23 21 1 2 3 4",
                "1 2 1 -1 21 23 1 2 4 3",
                "3 1",
                "1 3 2 -2 23 21 1 2 3 4",
            ]
        ),
        encoding="utf-8",
    )

    assert amplicol_process_file_integrals(path) == ((1, 1), (1, 2), (3, 1))


def test_reorder_external_momenta_by_fortran_pdg_order() -> None:
    particles = (
        ExternalMomentum(1, (500.0, 0.0, 0.0, 500.0)),
        ExternalMomentum(-1, (500.0, 0.0, 0.0, -500.0)),
        ExternalMomentum(23, (504.0, -304.0, 209.0, 331.0)),
        ExternalMomentum(21, (496.0, 304.0, -209.0, -331.0)),
    )

    ordered = reorder_external_momenta_by_pdg(particles, (1, -1, 21, 23))

    assert tuple(particle.pdg for particle in ordered) == (1, -1, 21, 23)
    assert ordered[2].momentum == particles[3].momentum
    assert ordered[3].momentum == particles[2].momentum


def test_reference_adapter_reorders_supplied_momenta_probe_to_process_file_order(
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

    particles = (
        ExternalMomentum(1, (500.0, 0.0, 0.0, 500.0)),
        ExternalMomentum(-1, (500.0, 0.0, 0.0, -500.0)),
        ExternalMomentum(23, (504.0, -304.0, 209.0, 331.0)),
        ExternalMomentum(21, (496.0, 304.0, -209.0, -331.0)),
    )
    adapter = AmplicolAdapter(tmp_path, runner=runner)

    adapter.run_amplicol_momenta_probe("d d~ > z g", particles=particles)

    written = (
        tmp_path / "Utilities" / "ME_checks" / "momenta_1_1.txt"
    ).read_text(encoding="utf-8")
    lines = written.splitlines()
    assert calls == [
        (
            "./amplicol_generate",
            "--amplicol_momenta_probe=1",
            "--nevents=1",
            "--timing=1",
            f"--process={tmp_path / 'processes.txt'}",
        )
    ]
    assert lines[2].startswith("4.96000000000000000e+02")
    assert lines[3].startswith("5.04000000000000000e+02")


def test_reference_adapter_stops_library_preparation_after_failed_command(
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
        returncode = 2 if tuple(args) == ("make", "-j3", "amplicol_generate") else 0
        return CommandResult(tuple(args), cwd, returncode, "make stdout", "make stderr", 0.01)

    adapter = AmplicolAdapter(tmp_path, runner=runner, jobs=3)

    with pytest.raises(RuntimeError, match="amplicol_generate"):
        adapter.prepare_library("d d~ > z g g")

    assert calls == [
        ("make", "cleanlib"),
        ("make", "-j3", "amplicol_generate"),
    ]


def test_reference_adapter_can_use_legacy_process_list_backend(
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
        (cwd / "processes.txt").write_text("4 1\n1 -1 21 23\n", encoding="utf-8")
        return CommandResult(tuple(args), cwd, 0, "", "", 0.01)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    process_file = adapter.write_process_file(
        "d d~ > z g",
        process_list_backend="legacy",
    )

    assert process_file == tmp_path / "processes.txt"
    assert process_file.read_text(encoding="utf-8").startswith("4 1")
    assert calls == [("python3", "process_list.py", "--serial", "d d~ > z g")]


def test_reference_adapter_rejects_missing_legacy_process_file(
    tmp_path: Path,
) -> None:
    (tmp_path / "processes.txt").write_text("stale\n", encoding="utf-8")

    def runner(
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult:
        return CommandResult(tuple(args), cwd, 0, "", "", 0.01)

    adapter = AmplicolAdapter(tmp_path, runner=runner)

    with pytest.raises(RuntimeError, match="did not produce processes.txt"):
        adapter.write_process_file(
            "u d~ > e+ ve",
            process_list_backend="legacy",
        )

    assert not (tmp_path / "processes.txt").exists()


def test_reference_adapter_runs_generated_library_use_path(tmp_path: Path) -> None:
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
------------------------------------------------------------------------------
Timing summary                           seconds    percent  note
------------------------------------------------------------------------------
amplitude evaluation                    0.000042      1.00%
------------------------------------------------------------------------------
"""
        return CommandResult(tuple(args), cwd, 0, stdout, "", 0.01)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    result = adapter.run_library_use(
        "d d~ > z g",
        nevents=1234,
        seed=101,
        timing_sample=17,
    )

    assert calls == [
        (
            "./amplicol_generate",
            "--library=use",
            "--nevents=1234",
            "--seed=101",
            "--timing=17",
        )
    ]
    assert result.process_file == tmp_path / "processes.txt"
    assert [row.label for row in result.timing_rows] == ["amplitude evaluation"]
    assert result.timing_rows[0].seconds == 0.000042


def test_fortran_probe_normalization_uses_phase_space_group_external_count() -> None:
    """Guard against stale module-level ``next`` under ``--library=use``.

    In library mode the AmpliCol probe can evaluate a loaded phase-space group
    while the module variable named ``next`` is zero.  The ME-check and direct
    probe normalization must therefore use ``pgl(ichan)%next``.
    """

    repo_root = Path(__file__).resolve().parents[3]
    source = (repo_root / "mg_checks.f03").read_text(encoding="utf-8")

    assert "**(pgl(ichan)%next-2-pgl(ichan)%amps(iint)%n_sing(1))" in source
    assert "**(next-2-pgl(ichan)%amps(iint)%n_sing(1))" not in source


def test_fortran_supplied_momenta_probe_can_create_libraries() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    source = (repo_root / "amplicol_generate.f03").read_text(encoding="utf-8")

    assert "call optimise_the_amplitudes(iint,igroup,done)" in source
    assert "call create_amplitude_lib()" in source
    assert "Supplied-momenta library creation did not create all amplitude libraries" in source


def test_fortran_scalar_current_library_emitter_handles_single_vertex_sum() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    source = (repo_root / "amplitude_QCD.f03").read_text(encoding="utf-8")

    assert "val_c(1,int1(0,i))=int_c(1,int1(1,i))" in source
    assert "val_c(1,int1(0,i))=sum(int_c(1,int1(1:" in source
    assert "val_c(1,int1(0,i))=sum(int_c(1,int1(1:" in source
    assert "val_c(1,int1(0,i))=sum(int_c(1,int1(1:'//trim(adjustl(tmp))//',i)),dim=2)" not in source


def test_reference_color_order_for_run_reads_process_file_entry(tmp_path: Path) -> None:
    adapter = AmplicolAdapter(tmp_path)
    process_file = adapter.write_process_file("d d~ > z g")
    run = AmplicolWorkflowResult(
        commands=(),
        process_file=process_file,
        probe_points=(
            AmplicolProbePoint(
                point=1,
                group=1,
                integral=1,
                particles=(),
                matrix_element=1.0,
            ),
        ),
    )

    color_order = reference_color_order_for_run(run)

    assert color_order is not None
    assert len(color_order) == 4
    assert sorted(color_order) == [1, 2, 3, 4]


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
            "--nevents=3",
            "--timing=3",
            f"--process={tmp_path / 'processes.txt'}",
        )
    ]
    assert result.first_phase_space_point is not None
    assert result.first_phase_space_point.matrix_element == 142.0
    assert len(result.probe_points) == 1
    assert result.probe_points[0].matrix_element == 142.0


def test_reference_adapter_amplicol_probe_can_use_generated_library(
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
        return CommandResult(tuple(args), cwd, 0, "", "", 0.02)

    adapter = AmplicolAdapter(tmp_path, runner=runner)
    adapter.run_amplicol_probe("d d~ > z", points=3, use_library=True)

    assert calls == [
        (
            "./amplicol_generate",
            "--amplicol_probe=3",
            "--nevents=3",
            "--timing=3",
            f"--process={tmp_path / 'processes.txt'}",
            "--library=use",
        )
    ]


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
            "--nevents=4",
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
            "--nevents=5",
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
