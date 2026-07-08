from __future__ import annotations

import os
import re
import signal
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Mapping, Protocol, Sequence

from .core_types import ExternalMomentum
from .processes import ProcessOptions, ProcessEnumerator

ProcessListBackend = Literal["python", "legacy"]


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
    except KeyboardInterrupt:
        _terminate_process_group(process.pid, signal.SIGTERM)
        try:
            process.communicate(timeout=5.0)
        except subprocess.TimeoutExpired:
            _terminate_process_group(process.pid, signal.SIGKILL)
            process.communicate()
        raise
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
    color_probe_components: tuple[float, float, float] | None = None
    color_probe_raw_components: tuple[float, float, float] | None = None

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


def amplicol_process_file_entry(
    process_file: str | Path,
    *,
    group: int,
    integral: int,
) -> dict[str, list[int]] | None:
    """Return the legacy process-file row for a group/integral pair."""

    path = Path(process_file)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None
    cursor = 0

    def next_nonempty() -> str | None:
        nonlocal cursor
        while cursor < len(lines):
            line = lines[cursor].strip()
            cursor += 1
            if line:
                return line
        return None

    header = next_nonempty()
    if header is None:
        return None
    try:
        external_count, unique_count = (int(value) for value in header.split()[:2])
    except (ValueError, IndexError):
        return None
    for _ in range(unique_count):
        if next_nonempty() is None:
            return None
    group_count_line = next_nonempty()
    if group_count_line is None:
        return None
    try:
        group_count = int(group_count_line.split()[0])
    except (ValueError, IndexError):
        return None
    for _ in range(group_count):
        group_header = next_nonempty()
        if group_header is None:
            return None
        try:
            group_tokens = [int(token) for token in group_header.split()]
            group_id = group_tokens[0]
            process_count = group_tokens[1]
        except (ValueError, IndexError):
            return None
        for row_index in range(1, process_count + 1):
            row = next_nonempty()
            if row is None:
                return None
            tokens = row.split()
            try:
                channel_count = int(tokens[0])
            except (ValueError, IndexError):
                return None
            start = 1 + channel_count
            stop_process = start + external_count
            stop_order = stop_process + external_count
            if len(tokens) < stop_order:
                return None
            try:
                process = [int(token) for token in tokens[start:stop_process]]
                color_order = [int(token) for token in tokens[stop_process:stop_order]]
            except ValueError:
                return None
            if group_id == group and row_index == integral:
                return {
                    "process": process,
                    "color_order": color_order,
                }
    return None


def amplicol_process_file_integrals(process_file: str | Path) -> tuple[tuple[int, int], ...]:
    """Return ``(group, integral)`` entries present in a legacy process file."""

    path = Path(process_file)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return ()
    cursor = 0

    def next_nonempty() -> str | None:
        nonlocal cursor
        while cursor < len(lines):
            line = lines[cursor].strip()
            cursor += 1
            if line:
                return line
        return None

    header = next_nonempty()
    if header is None:
        return ()
    try:
        _, unique_count = (int(value) for value in header.split()[:2])
    except (ValueError, IndexError):
        return ()
    for _ in range(unique_count):
        if next_nonempty() is None:
            return ()
    group_count_line = next_nonempty()
    if group_count_line is None:
        return ()
    try:
        group_count = int(group_count_line.split()[0])
    except (ValueError, IndexError):
        return ()
    entries: list[tuple[int, int]] = []
    for _ in range(group_count):
        group_header = next_nonempty()
        if group_header is None:
            return ()
        try:
            group_tokens = [int(token) for token in group_header.split()]
            group_id = group_tokens[0]
            process_count = group_tokens[1]
        except (ValueError, IndexError):
            return ()
        for row_index in range(1, process_count + 1):
            if next_nonempty() is None:
                return ()
            entries.append((group_id, row_index))
    return tuple(entries)


def reorder_external_momenta_by_pdg(
    particles: Sequence[ExternalMomentum],
    expected_pdgs: Sequence[int],
) -> tuple[ExternalMomentum, ...]:
    """Return ``particles`` in ``expected_pdgs`` order using stable PDG matching.

    The generated Fortran process file is allowed to use an external ordering
    different from pyAmpliCol's canonical process order, for example
    ``d d~ > g z`` for a user request ``d d~ > z g``.  Supplied momenta used
    for AmpliCol library warmup/probing must therefore be reordered before
    being written to ``Utilities/ME_checks/momenta_*.txt``.  Equal-PDG legs are
    matched stably, preserving the user's order among identical particles.
    """

    remaining = list(particles)
    ordered: list[ExternalMomentum] = []
    for expected in expected_pdgs:
        for index, particle in enumerate(remaining):
            if int(particle.pdg) == int(expected):
                ordered.append(particle)
                del remaining[index]
                break
        else:
            raise ValueError(
                "external momenta cannot be reordered to Fortran process "
                f"PDG order {list(expected_pdgs)}"
            )
    if remaining:
        raise ValueError(
            "external momenta contain extra particles after Fortran PDG reordering"
        )
    return tuple(ordered)


def reference_color_order_for_run(
    run: AmplicolWorkflowResult,
) -> tuple[int, ...] | None:
    """Extract the leading-colour ordering used by a probe/ME-check run."""

    reference_points: list[AmplicolProbePoint | AmplicolFirstPoint] = list(
        run.probe_points
    )
    if not reference_points and run.first_phase_space_point is not None:
        reference_points = [run.first_phase_space_point]
    if not reference_points:
        return None
    point = reference_points[0]
    entry = amplicol_process_file_entry(
        run.process_file,
        group=point.group,
        integral=point.integral,
    )
    if entry is None:
        return None
    return tuple(int(label) for label in entry["color_order"])


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
        process_list_backend: ProcessListBackend = "python",
    ) -> Path:
        output = self.repo_root / "processes.txt" if path is None else Path(path)
        if not output.is_absolute():
            output = self.repo_root / output
        if process_list_backend == "legacy":
            self.run_legacy_process_list(
                process,
                options=options,
                output=output,
            )
            return output
        if process_list_backend != "python":
            raise ValueError(
                "process_list_backend must be either 'python' or 'legacy', "
                f"got {process_list_backend!r}"
            )
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
        process_list_backend: ProcessListBackend = "python",
        warmup_particles: Sequence[ExternalMomentum] | None = None,
        warmup_points: int = 10,
    ) -> AmplicolWorkflowResult:
        """Create and compile the generated Fortran amplitude library.

        When ``warmup_particles`` are supplied, ``--library=create`` is driven
        by deterministic supplied momenta instead of the legacy phase-space
        integrator.  This still exercises the intended create/compile/use
        library chain, but avoids reference-side phase-space failures in
        point-by-point validation jobs.
        """

        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        create_args = ["./amplicol_generate", "--library=create", f"--process={path}"]
        if warmup_particles is not None:
            entries = amplicol_process_file_integrals(path) or ((1, 1),)
            for group, integral in entries:
                entry = amplicol_process_file_entry(
                    path,
                    group=group,
                    integral=integral,
                )
                ordered_particles: Sequence[ExternalMomentum] = warmup_particles
                if entry is not None:
                    ordered_particles = reorder_external_momenta_by_pdg(
                        warmup_particles,
                        entry["process"],
                    )
                self.write_momenta_probe_file(
                    ordered_particles,
                    group=group,
                    integral=integral,
                )
            create_args.extend(
                [
                    f"--amplicol_momenta_probe={warmup_points}",
                    "--amplicol_probe_quiet",
                ]
            )
        commands: list[CommandResult] = []
        for args in (
            ["make", "cleanlib"],
            ["make", f"-j{self.jobs}", "amplicol_generate"],
            create_args,
            ["make", f"-j{self.jobs}", "amplicol_generate_library"],
        ):
            command = self._run(args)
            commands.append(command)
            command.check_returncode()
        return AmplicolWorkflowResult(commands=tuple(commands), process_file=path)

    def prepare_direct_probe(
        self,
        process: str,
        *,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        commands: list[CommandResult] = []
        for args in (
            ["make", "cleanlib"],
            ["make", f"-j{self.jobs}", "amplicol_generate"],
        ):
            command = self._run(args)
            commands.append(command)
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
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
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

    def run_library_use(
        self,
        process: str,
        *,
        nevents: int = 10000,
        seed: int | None = None,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        timing_sample: int | None = None,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        """Run the generated-library integration path.

        This is the benchmark path recommended by the Fortran AmpliCol author:
        first call :meth:`prepare_library`, then run ``./amplicol_generate
        --library=use``.  Direct probes remain useful for deterministic
        point-by-point debugging, but they do not exercise the generated
        library.
        """

        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        args = [
            "./amplicol_generate",
            "--library=use",
            f"--nevents={nevents}",
        ]
        if seed is not None:
            args.append(f"--seed={seed}")
        if timing_sample is not None:
            args.append(f"--timing={timing_sample}")
        if process_file is not None:
            args.append(f"--process={path}")
        command = self._run(args)
        command.check_returncode()
        output = "\n".join([command.stdout, command.stderr])
        return AmplicolWorkflowResult(
            commands=(command,),
            process_file=path,
            timing_rows=parse_timing_rows(output),
            first_point_matrix_element=parse_first_matrix_element(output),
            first_phase_space_point=parse_first_phase_space_point(output),
            probe_points=parse_amplicol_probe_points(output),
        )

    def run_library_benchmark(
        self,
        process: str,
        *,
        points: int = 100000,
        group: int = 1,
        integral: int = 1,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        """Run the direct generated-library timing driver.

        The driver calls the generated ``amp_lib:evaluate_amp`` dispatcher
        directly on the saved library test momentum.  It intentionally bypasses
        the AmpliCol integration/probe control flow, phase-space generation,
        PDFs, cuts, and event-weight bookkeeping.  This is the cleanest runtime
        benchmark for the generated Fortran amplitude library itself.
        """

        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        build = self._run(["make", f"-j{self.jobs}", "amplicol_library_benchmark"])
        build.check_returncode()
        command = self._run(
            [
                "./amplicol_library_benchmark",
                str(max(1, points)),
                str(group),
                str(integral),
            ]
        )
        command.check_returncode()
        output = "\n".join([command.stdout, command.stderr])
        return AmplicolWorkflowResult(
            commands=(build, command),
            process_file=path,
            timing_rows=parse_timing_rows(output),
        )

    def run_color_probe(
        self,
        process: str,
        *,
        color_accuracy: str,
        particles: Sequence[ExternalMomentum],
        points: int = 1,
        group: int = 1,
        integral: int = 1,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        """Run the direct AmpliCol LC/NLC/full colour reference probe.

        The probe uses AmpliCol's ``imode=2`` all-colour-order amplitude and
        ``init_col`` sparse colour matrix, then applies the same coupling and
        averaging normalization as ``AMPICOL_PROBE_VALUE``.  It is the intended
        point-reference path for NLC/full validation.
        """

        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        entry = amplicol_process_file_entry(path, group=group, integral=integral)
        ordered_particles: Sequence[ExternalMomentum] = particles
        if entry is not None:
            ordered_particles = reorder_external_momenta_by_pdg(
                particles,
                entry["process"],
            )
        momenta_path = self.write_momenta_probe_file(
            ordered_particles,
            group=group,
            integral=integral,
        )
        build = self._run(["make", f"-j{self.jobs}", "amplicol_color_probe"])
        build.check_returncode()
        command = self._run(
            [
                "./amplicol_color_probe",
                str(max(1, points)),
                str(group),
                str(integral),
                color_accuracy,
                str(path),
                str(momenta_path),
            ]
        )
        command.check_returncode()
        output = "\n".join([command.stdout, command.stderr])
        return AmplicolWorkflowResult(
            commands=(build, command),
            process_file=path,
            timing_rows=parse_timing_rows(output),
            first_point_matrix_element=parse_color_probe_value(output),
            color_probe_components=parse_color_probe_components(output),
            color_probe_raw_components=parse_color_probe_raw_components(output),
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
        use_library: bool = False,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        timing = points if timing_sample is None else timing_sample
        args = [
            "./amplicol_generate",
            f"--amplicol_probe={points}",
            f"--nevents={max(1, points)}",
            f"--timing={timing}",
            f"--process={path}",
        ]
        if use_library:
            args.append("--library=use")
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
        use_library: bool = False,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        timing = points if timing_sample is None else timing_sample
        args = [
            "./amplicol_generate",
            f"--amplicol_fixed_probe={points}",
            f"--nevents={max(1, points)}",
            f"--timing={timing}",
            f"--process={path}",
        ]
        if use_library:
            args.append("--library=use")
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

    def write_momenta_probe_file(
        self,
        particles: Sequence[ExternalMomentum],
        *,
        group: int = 1,
        integral: int = 1,
    ) -> Path:
        directory = self.repo_root / "Utilities" / "ME_checks"
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"momenta_{group}_{integral}.txt"
        with path.open("w", encoding="utf-8") as handle:
            for particle in particles:
                handle.write(
                    " ".join(f"{component:.17e}" for component in particle.momentum)
                    + "\n"
                )
        return path

    def run_amplicol_momenta_probe(
        self,
        process: str,
        *,
        particles: Sequence[ExternalMomentum],
        points: int = 1,
        process_file: str | Path | None = None,
        options: ProcessOptions | None = None,
        timing_sample: int | None = None,
        quiet: bool = False,
        use_library: bool = False,
        process_list_backend: ProcessListBackend = "python",
    ) -> AmplicolWorkflowResult:
        path = self.write_process_file(
            process,
            process_file,
            options=options,
            process_list_backend=process_list_backend,
        )
        ordered_particles: Sequence[ExternalMomentum] = particles
        entry = amplicol_process_file_entry(path, group=1, integral=1)
        if entry is not None:
            ordered_particles = reorder_external_momenta_by_pdg(
                particles,
                entry["process"],
            )
        self.write_momenta_probe_file(ordered_particles, group=1, integral=1)
        timing = points if timing_sample is None else timing_sample
        args = [
            "./amplicol_generate",
            f"--amplicol_momenta_probe={points}",
            f"--nevents={max(1, points)}",
            f"--timing={timing}",
            f"--process={path}",
        ]
        if use_library:
            args.append("--library=use")
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
        legacy_output = self.repo_root / "processes.txt"
        try:
            legacy_output.unlink()
        except FileNotFoundError:
            pass
        command = self._run(args)
        command.check_returncode()
        target = self.repo_root / output
        if not legacy_output.exists():
            raise RuntimeError(
                "legacy process_list.py did not produce processes.txt for "
                f"process {process!r}; required process-list options may be missing"
            )
        if target.name != "processes.txt":
            legacy_output.replace(target)
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


def parse_color_probe_value(output: str) -> float | None:
    match = re.search(
        r"^\s*AMPICOL_COLOR_PROBE_VALUE\s+"
        r"\S+\s+\d+\s+\d+\s+([0-9.Ee+-]+)\s*$",
        output,
        re.MULTILINE,
    )
    return None if match is None else float(match.group(1))


def parse_color_probe_components(output: str) -> tuple[float, float, float] | None:
    return _parse_color_probe_component_line(output, "AMPICOL_COLOR_PROBE_COMPONENTS")


def parse_color_probe_raw_components(output: str) -> tuple[float, float, float] | None:
    return _parse_color_probe_component_line(
        output,
        "AMPICOL_COLOR_PROBE_RAW_COMPONENTS",
    )


def _parse_color_probe_component_line(
    output: str,
    label: str,
) -> tuple[float, float, float] | None:
    escaped = re.escape(label)
    match = re.search(
        rf"^\s*{escaped}\s+"
        r"([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s+([0-9.Ee+-]+)\s*$",
        output,
        re.MULTILINE,
    )
    if match is None:
        return None
    return (float(match.group(1)), float(match.group(2)), float(match.group(3)))


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
    "amplicol_process_file_entry",
    "amplicol_process_file_integrals",
    "parse_color_probe_components",
    "parse_color_probe_raw_components",
    "parse_color_probe_value",
    "parse_amplicol_probe_points",
    "parse_first_phase_space_point",
    "parse_first_matrix_element",
    "parse_timing_rows",
    "popen_runner",
    "reference_color_order_for_run",
    "reorder_external_momenta_by_pdg",
]
