from __future__ import annotations

import os
import re
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Protocol, Sequence

from .native import ExternalMomentum
from .processes import ProcessOptions, ProcessEnumerator


@dataclass(frozen=True)
class CommandResult:
    args: tuple[str, ...]
    cwd: Path
    returncode: int
    stdout: str
    stderr: str
    elapsed_s: float

    def check_returncode(self) -> None:
        if self.returncode != 0:
            raise RuntimeError(
                f"command failed with exit code {self.returncode}: {' '.join(self.args)}\n"
                f"stdout:\n{self.stdout}\n\nstderr:\n{self.stderr}"
            )


class CommandRunner(Protocol):
    def __call__(
        self,
        args: Sequence[str],
        *,
        cwd: Path,
        env: Mapping[str, str] | None = None,
        timeout: float | None = None,
    ) -> CommandResult: ...


def popen_runner(
    args: Sequence[str],
    *,
    cwd: Path,
    env: Mapping[str, str] | None = None,
    timeout: float | None = None,
) -> CommandResult:
    start = time.perf_counter()
    process = subprocess.Popen(
        list(args),
        cwd=cwd,
        env=None if env is None else {**os.environ, **dict(env)},
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
        text=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        _terminate_process_group(process.pid, signal.SIGTERM)
        try:
            stdout, stderr = process.communicate(timeout=5.0)
        except subprocess.TimeoutExpired:
            _terminate_process_group(process.pid, signal.SIGKILL)
            stdout, stderr = process.communicate()
    return CommandResult(
        args=tuple(args),
        cwd=cwd,
        returncode=process.returncode,
        stdout=stdout,
        stderr=stderr,
        elapsed_s=time.perf_counter() - start,
    )


def _terminate_process_group(pid: int, sig: signal.Signals) -> None:
    try:
        os.killpg(pid, sig)
    except ProcessLookupError:
        return


@dataclass(frozen=True)
class TimingRow:
    label: str
    seconds: float
    note: str = ""


@dataclass(frozen=True)
class AmplicolWorkflowResult:
    commands: tuple[CommandResult, ...]
    process_file: Path
    timing_rows: tuple[TimingRow, ...] = ()
    first_point_matrix_element: float | None = None
    first_phase_space_point: "AmplicolFirstPoint | None" = None
    probe_points: tuple["AmplicolProbePoint", ...] = ()

    @property
    def total_command_time_s(self) -> float:
        return sum(command.elapsed_s for command in self.commands)


@dataclass(frozen=True)
class AmplicolFirstPoint:
    group: int
    integral: int
    particles: tuple[ExternalMomentum, ...]
    matrix_element: float | None = None


@dataclass(frozen=True)
class AmplicolProbePoint:
    point: int
    group: int
    integral: int
    particles: tuple[ExternalMomentum, ...]
    matrix_element: float


class AmplicolAdapter:
    """Python steering layer for the legacy Fortran AmpliCol executable."""

    def __init__(
        self,
        repo_root: str | Path,
        *,
        runner: CommandRunner = popen_runner,
        jobs: int = 8,
        timeout: float | None = None,
    ) -> None:
        self.repo_root = Path(repo_root).resolve()
        self.runner = runner
        self.jobs = jobs
        self.timeout = timeout

    def write_process_file(
        self,
        process: str,
        path: str | Path | None = None,
        *,
        options: ProcessOptions | None = None,
    ) -> Path:
        output = self.repo_root / "processes.txt" if path is None else Path(path)
        if not output.is_absolute():
            output = self.repo_root / output
        enumerator = ProcessEnumerator(options)
        enumeration = enumerator.enumerate(process)
        enumerator.write_legacy_file(enumeration, output)
        return output

    def prepare_library(
        self,
        process: str,
        *,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(process, process_file, options=options)
        commands = [
            self._run(["make", "cleanlib"]),
            self._run(["make", f"-j{self.jobs}", "amplicol_generate"]),
            self._run(["./amplicol_generate", "--library=create", f"--process={path}"]),
            self._run(["make", f"-j{self.jobs}", "amplicol_generate_library"]),
        ]
        for command in commands:
            command.check_returncode()
        return AmplicolWorkflowResult(commands=tuple(commands), process_file=path)

    def prepare_direct_probe(
        self,
        process: str,
        *,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(process, process_file, options=options)
        commands = [
            self._run(["make", "cleanlib"]),
            self._run(["make", f"-j{self.jobs}", "amplicol_generate"]),
        ]
        for command in commands:
            command.check_returncode()
        return AmplicolWorkflowResult(commands=tuple(commands), process_file=path)

    def run_me_test(
        self,
        process: str,
        *,
        points: int = 10,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        mg5_path: str | Path | None = None,
        timing_sample: int | None = None,
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(process, process_file, options=options)
        env = {}
        if mg5_path is not None:
            env["MG5_PATH"] = str(mg5_path)
        timing = points if timing_sample is None else timing_sample
        command = self._run(
            [
                "./amplicol_generate",
                f"--me_test={points}",
                f"--timing={timing}",
                f"--process={path}",
            ],
            env=env,
        )
        command.check_returncode()
        output = "\n".join([command.stdout, command.stderr])
        probe_points = parse_amplicol_probe_points(output)
        return AmplicolWorkflowResult(
            commands=(command,),
            process_file=path,
            timing_rows=parse_timing_rows(output),
            first_point_matrix_element=parse_first_matrix_element(output),
            first_phase_space_point=parse_first_phase_space_point(output),
            probe_points=probe_points,
        )

    def run_amplicol_probe(
        self,
        process: str,
        *,
        points: int = 10,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        timing_sample: int | None = None,
        quiet: bool = False,
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(process, process_file, options=options)
        timing = points if timing_sample is None else timing_sample
        args = [
            "./amplicol_generate",
            f"--amplicol_probe={points}",
            f"--timing={timing}",
            f"--process={path}",
        ]
        if quiet:
            args.append("--amplicol_probe_quiet")
        command = self._run(args)
        command.check_returncode()
        output = "\n".join([command.stdout, command.stderr])
        probe_points = parse_amplicol_probe_points(output)
        return AmplicolWorkflowResult(
            commands=(command,),
            process_file=path,
            timing_rows=parse_timing_rows(output),
            first_point_matrix_element=parse_first_matrix_element(output),
            first_phase_space_point=parse_first_phase_space_point(output),
            probe_points=probe_points,
        )

    def run_amplicol_fixed_probe(
        self,
        process: str,
        *,
        points: int = 10,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        timing_sample: int | None = None,
        quiet: bool = False,
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(process, process_file, options=options)
        timing = points if timing_sample is None else timing_sample
        args = [
            "./amplicol_generate",
            f"--amplicol_fixed_probe={points}",
            f"--timing={timing}",
            f"--process={path}",
        ]
        if quiet:
            args.append("--amplicol_probe_quiet")
        command = self._run(args)
        command.check_returncode()
        output = "\n".join([command.stdout, command.stderr])
        probe_points = parse_amplicol_probe_points(output)
        return AmplicolWorkflowResult(
            commands=(command,),
            process_file=path,
            timing_rows=parse_timing_rows(output),
            first_point_matrix_element=parse_first_matrix_element(output),
            first_phase_space_point=parse_first_phase_space_point(output),
            probe_points=probe_points,
        )

    def run_legacy_process_list(
        self,
        process: str,
        *,
        options: ProcessOptions | None = None,
        output: str | Path = "processes.txt",
    ) -> CommandResult:
        opts = options or ProcessOptions()
        args = ["python3", "process_list.py", "--serial"]
        if opts.flavour_scheme != 5:
            args.extend(["-FS", str(opts.flavour_scheme)])
        if opts.include_3qqbar:
            args.append("-3")
        if opts.include_cc:
            args.append("-cc")
        if opts.include_resonance:
            args.append("-res")
        args.append(process)
        command = self._run(args)
        command.check_returncode()
        target = self.repo_root / output
        if target.name != "processes.txt":
            (self.repo_root / "processes.txt").replace(target)
        return command

    def _run(
        self,
        args: Sequence[str],
        *,
        env: Mapping[str, str] | None = None,
    ) -> CommandResult:
        return self.runner(args, cwd=self.repo_root, env=env, timeout=self.timeout)


def parse_timing_rows(output: str) -> tuple[TimingRow, ...]:
    rows: list[TimingRow] = []
    in_timing_summary = False
    pattern = re.compile(
        r"^\s*(?P<label>[A-Za-z][A-Za-z0-9 /_-]+?)\s+"
        r"(?P<seconds>[0-9.Ee+-]+)"
        r"(?:\s+(?P<percent>[0-9.]+%))?"
        r"(?:\s+(?P<note>.*?))?\s*$"
    )
    for line in output.splitlines():
        if "Timing summary" in line:
            in_timing_summary = True
            continue
        if not in_timing_summary:
            continue
        stripped = line.strip()
        if not stripped or set(stripped) == {"-"}:
            continue
        match = pattern.match(line)
        if match is None:
            continue
        try:
            seconds = float(match.group("seconds"))
        except ValueError:
            continue
        rows.append(TimingRow(match.group("label").strip(), seconds, match.group("note") or ""))
    return tuple(rows)


def parse_first_matrix_element(output: str) -> float | None:
    match = re.search(r"AmpliCol matrix element:\s*([0-9.Ee+-]+)", output)
    return None if match is None else float(match.group(1))


def parse_first_phase_space_point(output: str) -> AmplicolFirstPoint | None:
    header = re.search(
        r"(?:ME-check|AmpliCol probe) first phase-space point\s+group,\s+integral:\s*(\d+)\s+(\d+)",
        output,
        re.MULTILINE,
    )
    if header is None:
        return None

    particles: list[ExternalMomentum] = []
    lines = output[header.end() :].splitlines()
    row_pattern = re.compile(
        r"^\s*(?P<i>\d+)\s+(?P<pdg>[+-]?\d+)\s+"
        r"(?P<e>[0-9.Ee+-]+)\s+"
        r"(?P<px>[0-9.Ee+-]+)\s+"
        r"(?P<py>[0-9.Ee+-]+)\s+"
        r"(?P<pz>[0-9.Ee+-]+)\s*$"
    )
    for line in lines:
        row = row_pattern.match(line)
        if row is not None:
            particles.append(
                ExternalMomentum(
                    int(row.group("pdg")),
                    (
                        float(row.group("e")),
                        float(row.group("px")),
                        float(row.group("py")),
                        float(row.group("pz")),
                    ),
                )
            )
            continue
        if particles and "AmpliCol matrix element:" in line:
            break

    if not particles:
        return None
    return AmplicolFirstPoint(
        group=int(header.group(1)),
        integral=int(header.group(2)),
        particles=tuple(particles),
        matrix_element=parse_first_matrix_element(output),
    )


def parse_amplicol_probe_points(output: str) -> tuple[AmplicolProbePoint, ...]:
    values: dict[int, tuple[int, int, float]] = {}
    particles: dict[int, list[tuple[int, ExternalMomentum]]] = {}
    value_pattern = re.compile(
        r"^\s*AMPICOL_PROBE_VALUE\s+"
        r"(?P<point>\d+)\s+"
        r"(?P<group>\d+)\s+"
        r"(?P<integral>\d+)\s+"
        r"(?P<me>[0-9.Ee+-]+)\s*$"
    )
    momentum_pattern = re.compile(
        r"^\s*AMPICOL_PROBE_MOM\s+"
        r"(?P<point>\d+)\s+"
        r"(?P<i>\d+)\s+"
        r"(?P<pdg>[+-]?\d+)\s+"
        r"(?P<e>[0-9.Ee+-]+)\s+"
        r"(?P<px>[0-9.Ee+-]+)\s+"
        r"(?P<py>[0-9.Ee+-]+)\s+"
        r"(?P<pz>[0-9.Ee+-]+)\s*$"
    )
    for line in output.splitlines():
        value = value_pattern.match(line)
        if value is not None:
            point = int(value.group("point"))
            values[point] = (
                int(value.group("group")),
                int(value.group("integral")),
                float(value.group("me")),
            )
            continue
        momentum = momentum_pattern.match(line)
        if momentum is None:
            continue
        point = int(momentum.group("point"))
        particles.setdefault(point, []).append(
            (
                int(momentum.group("i")),
                ExternalMomentum(
                    int(momentum.group("pdg")),
                    (
                        float(momentum.group("e")),
                        float(momentum.group("px")),
                        float(momentum.group("py")),
                        float(momentum.group("pz")),
                    ),
                ),
            )
        )

    parsed: list[AmplicolProbePoint] = []
    for point in sorted(values):
        group, integral, matrix_element = values[point]
        point_particles = tuple(
            particle for _, particle in sorted(particles.get(point, ()), key=lambda item: item[0])
        )
        if not point_particles:
            continue
        parsed.append(
            AmplicolProbePoint(
                point=point,
                group=group,
                integral=integral,
                particles=point_particles,
                matrix_element=matrix_element,
            )
        )
    return tuple(parsed)


__all__ = [
    "AmplicolAdapter",
    "AmplicolFirstPoint",
    "AmplicolProbePoint",
    "AmplicolWorkflowResult",
    "CommandResult",
    "TimingRow",
    "parse_amplicol_probe_points",
    "parse_first_phase_space_point",
    "parse_first_matrix_element",
    "parse_timing_rows",
    "popen_runner",
]
