from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence, cast

import numpy as np

from .core_types import ExternalMomentum
from .symbolica_evaluator import SymbolicaEvaluatorSettings
from .generic_artifact import (
    build_generic_process_manifest,
    select_leading_color_sector_ids,
    write_generic_dag_process_artifact,
)
from .phase_space import generic_validation_point
from .process_ir import build_process_ir
from .process_support import classify_process_support
from .processes import ProcessOptions
from .processes import PDGS
from .reference import (
    AmplicolAdapter,
    AmplicolWorkflowResult,
    ProcessListBackend,
    reference_color_order_for_run,
)


DEFAULT_LC_VALIDATION_PROCESSES: tuple[str, ...] = (
    "d d~ > z",
    "d d~ > z g",
    "d d~ > z g g",
    "d d~ > e+ e-",
    "d d~ > e+ e- g",
    "d d~ > e+ e- g g",
    "d d~ > e+ e- a",
    "d d~ > mu+ mu- g",
    "d d~ > mu+ mu- a",
    "u d~ > w+",
    "u d~ > w+ g",
    "u d~ > w+ g g",
    "d u~ > w-",
    "d u~ > w- g",
    "u d~ > e+ ve",
    "u d~ > e+ ve g",
    "u d~ > e+ ve z",
    "u d~ > e+ ve a",
    "u d~ > mu+ vm g",
    "d u~ > e- ve~",
    "d u~ > e- ve~ g",
    "d u~ > e- ve~ z",
    "d u~ > e- ve~ a",
    "d d~ > z z",
    "d d~ > z z g",
    "d d~ > z z g g",
    "d d~ > a z",
    "d d~ > a z g",
    "d d~ > a a",
    "d d~ > a a g",
    "d d~ > h z",
    "d d~ > h z g",
    "d d~ > h h z",
    "u d~ > h w+",
    "d u~ > h w-",
    "u d~ > w+ z",
    "u d~ > w+ a",
    "d u~ > w- z",
    "d u~ > w- a",
    "d d~ > w+ w- g",
    "d d~ > z z z",
    "d d~ > a z z",
    "g g > g g",
    "g g > g g g",
    "g g > u u~",
    "g g > u u~ g",
    "g g > d d~ g",
    "g g > t t~",
    "g g > t t~ g",
    "g g > t t~ h",
    "d d~ > t t~",
    "d d~ > t t~ g",
    "d d~ > u u~",
    "d d~ > u u~ g",
    "u d~ > u d~ g",
    "u u~ > d d~ g",
    "g g > u u~ d d~",
    "d d~ > u u~ s s~",
)


@dataclass(frozen=True)
class GenericLCValidationRow:
    process: str
    key: str
    supported: bool
    status: str
    concrete_process: str | None = None
    relative_difference: float | None = None
    reference_matrix_element: float | None = None
    pyamplicol_matrix_element: float | None = None
    reference_color_order: tuple[int, ...] | None = None
    selected_color_sector_ids: tuple[int, ...] | None = None
    artifact_dir: Path | None = None
    generation_s: float | None = None
    rusticol_runtime_s: float | None = None
    fortran_validation_workflow: str | None = None
    fortran_timing_workflow: str | None = None
    fortran_library_runtime_s_per_point: float | None = None
    error: str | None = None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "process": self.process,
            "key": self.key,
            "supported": self.supported,
            "status": self.status,
            "concrete_process": self.concrete_process,
            "relative_difference": self.relative_difference,
            "reference_matrix_element": self.reference_matrix_element,
            "pyamplicol_matrix_element": self.pyamplicol_matrix_element,
            "reference_color_order": (
                None
                if self.reference_color_order is None
                else list(self.reference_color_order)
            ),
            "selected_color_sector_ids": (
                None
                if self.selected_color_sector_ids is None
                else list(self.selected_color_sector_ids)
            ),
            "artifact_dir": None if self.artifact_dir is None else str(self.artifact_dir),
            "generation_s": self.generation_s,
            "rusticol_runtime_s": self.rusticol_runtime_s,
            "fortran_validation_workflow": self.fortran_validation_workflow,
            "fortran_timing_workflow": self.fortran_timing_workflow,
            "fortran_library_runtime_s_per_point": (
                self.fortran_library_runtime_s_per_point
            ),
            "error": self.error,
        }


@dataclass(frozen=True)
class GenericLCValidationSummary:
    rows: tuple[GenericLCValidationRow, ...]
    rel_tol: float

    @property
    def passed(self) -> bool:
        if not self.rows:
            return False
        return all(
            row.status in {"dry-run", "ok"}
            and row.supported
            and (
                row.relative_difference is None
                or row.relative_difference <= self.rel_tol
            )
            for row in self.rows
        )

    @property
    def max_relative_difference(self) -> float | None:
        values = [
            row.relative_difference
            for row in self.rows
            if row.relative_difference is not None
        ]
        return max(values) if values else None

    def to_json_dict(self) -> dict[str, object]:
        return {
            "kind": "pyamplicol-generic-lc-validation-summary",
            "process_count": len(self.rows),
            "passed": self.passed,
            "rel_tol": self.rel_tol,
            "max_relative_difference": self.max_relative_difference,
            "rows": [row.to_json_dict() for row in self.rows],
        }


def validate_generic_lc_processes(
    processes: Sequence[str] = DEFAULT_LC_VALIDATION_PROCESSES,
    *,
    output_dir: str | Path,
    amplicol_root: str | Path,
    rel_tol: float = 1.0e-8,
    options: ProcessOptions | None = None,
    evaluator_backend: str = "jit",
    compiled_preset: str = "runtime-o3",
    batch_size: int = 64,
    n_cores: int = 4,
    jobs: int = 8,
    timeout: float | None = None,
    dry_run: bool = False,
    library_timing_events: int = 0,
    seed: int = 101,
    reference_process_list_backend: str = "legacy",
) -> GenericLCValidationSummary:
    """Validate generic LC schema-v2 artifacts against Fortran AmpliCol.

    Deterministic supplied-momenta comparisons are run through the generated
    Fortran library workflow: ``--library=create``, ``make
    amplicol_generate_library``, then ``--library=use``.  If
    ``library_timing_events`` is positive, the same generated library is also
    used for a warm Fortran runtime timing run.
    """

    if reference_process_list_backend not in {"legacy", "python"}:
        raise ValueError(
            "reference_process_list_backend must be either 'legacy' or 'python'"
        )
    typed_reference_process_list_backend = cast(
        ProcessListBackend,
        reference_process_list_backend,
    )
    adapter = AmplicolAdapter(amplicol_root, jobs=jobs, timeout=timeout)
    rows = [
        _validate_one_process(
            process,
            output_dir=Path(output_dir).expanduser(),
            adapter=adapter,
            rel_tol=rel_tol,
            options=options,
            evaluator_backend=evaluator_backend,
            compiled_preset=compiled_preset,
            batch_size=batch_size,
            n_cores=n_cores,
            dry_run=dry_run,
            library_timing_events=library_timing_events,
            seed=seed,
            reference_process_list_backend=typed_reference_process_list_backend,
        )
        for process in processes
    ]
    return GenericLCValidationSummary(rows=tuple(rows), rel_tol=rel_tol)


def _validate_one_process(
    process: str,
    *,
    output_dir: Path,
    adapter: AmplicolAdapter,
    rel_tol: float,
    options: ProcessOptions | None,
    evaluator_backend: str,
    compiled_preset: str,
    batch_size: int,
    n_cores: int,
    dry_run: bool,
    library_timing_events: int,
    seed: int,
    reference_process_list_backend: ProcessListBackend,
) -> GenericLCValidationRow:
    try:
        ir = build_process_ir(process, options=options)
    except ValueError as exc:
        return GenericLCValidationRow(
            process=process,
            key="",
            supported=False,
            status="unsupported",
            error=str(exc),
        )
    support = classify_process_support(process, options=options)
    if not support.runtime_artifact_supported:
        return GenericLCValidationRow(
            process=process,
            key=ir.key,
            supported=False,
            status="unsupported",
            error=support.artifact_unavailable_message,
        )
    if dry_run:
        return GenericLCValidationRow(
            process=process,
            key=ir.key,
            supported=True,
            status="dry-run",
            artifact_dir=output_dir / ir.key,
        )

    point = generic_validation_point(process, seed=seed)
    build = adapter.prepare_library(
        process,
        options=options,
        process_list_backend=reference_process_list_backend,
        warmup_particles=point,
    )
    run = adapter.run_amplicol_momenta_probe(
        process,
        particles=point,
        points=1,
        process_file=build.process_file,
        options=options,
        timing_sample=1,
        use_library=True,
        process_list_backend=reference_process_list_backend,
    )
    if not run.probe_points:
        return GenericLCValidationRow(
            process=process,
            key=ir.key,
            supported=True,
            status="failed",
            fortran_validation_workflow="generated_library_momenta_probe",
            error="Fortran AmpliCol did not return a probe point",
        )
    reference_point = run.probe_points[0]
    reference_color_order = reference_color_order_for_run(run)

    base_manifest = build_generic_process_manifest(process, options=options)
    concrete_process = _concrete_process_from_probe_if_needed(
        reference_point.particles,
        base_manifest.external_pdg_order,
    )
    generic_manifest = build_generic_process_manifest(
        process if concrete_process is None else concrete_process,
        options=options,
        reference_color_order=reference_color_order,
    )
    selected_color_sector_ids = select_leading_color_sector_ids(
        generic_manifest,
        reference_color_order=reference_color_order,
    )
    process_output_dir = output_dir / generic_manifest.key
    settings = SymbolicaEvaluatorSettings(
        backend=evaluator_backend,
        compiled_preset=compiled_preset,
        n_cores=n_cores,
        compiled_chunk_compile_workers=n_cores,
        compiled_output_dir=str(process_output_dir / "compiled"),
    )
    start = time.perf_counter()
    _, payload = write_generic_dag_process_artifact(
        generic_manifest,
        process_output_dir,
        options=options,
        evaluator_backend=evaluator_backend,
        compiled_preset=compiled_preset,
        batch_size=batch_size,
        emit_stage_evaluator_artifacts=True,
        symbolica_settings=settings,
        selected_color_sector_ids=selected_color_sector_ids,
    )
    generation_s = time.perf_counter() - start
    compiled = payload.get("compiled")
    if not isinstance(compiled, dict) or not compiled.get("runtime_available"):
        return GenericLCValidationRow(
            process=process,
            key=generic_manifest.key,
            supported=True,
            status="failed",
            concrete_process=concrete_process,
            reference_color_order=reference_color_order,
            selected_color_sector_ids=_sorted_sector_tuple(selected_color_sector_ids),
            artifact_dir=process_output_dir,
            generation_s=generation_s,
            fortran_validation_workflow="generated_library_momenta_probe",
            error=str(
                compiled.get("runtime_unavailable_message")
                if isinstance(compiled, dict)
                else "generic artifact did not contain compiled stage evaluators"
            ),
        )

    pyamplicol_value, rusticol_runtime_s = _evaluate_rusticol(
        process_output_dir,
        _reorder_particles_by_pdg(
            reference_point.particles,
            generic_manifest.external_pdg_order,
        ),
    )
    reference_value = float(reference_point.matrix_element)
    relative_difference = _relative_difference(reference_value, pyamplicol_value)
    timing_s_per_point = None
    timing_workflow = None
    if library_timing_events > 0:
        timing_run = adapter.run_library_use(
            process,
            nevents=library_timing_events,
            seed=seed,
            process_file=build.process_file,
            options=options,
            timing_sample=1,
            process_list_backend=reference_process_list_backend,
        )
        timing_s_per_point = _library_runtime_per_point(
            timing_run,
            library_timing_events,
        )
        timing_workflow = "generated_library_use"
    return GenericLCValidationRow(
        process=process,
        key=generic_manifest.key,
        supported=True,
        status="ok" if relative_difference <= rel_tol else "failed",
        concrete_process=concrete_process,
        relative_difference=relative_difference,
        reference_matrix_element=reference_value,
        pyamplicol_matrix_element=pyamplicol_value,
        reference_color_order=reference_color_order,
        selected_color_sector_ids=_sorted_sector_tuple(selected_color_sector_ids),
        artifact_dir=process_output_dir,
        generation_s=generation_s,
        rusticol_runtime_s=rusticol_runtime_s,
        fortran_validation_workflow="generated_library_momenta_probe",
        fortran_timing_workflow=timing_workflow,
        fortran_library_runtime_s_per_point=timing_s_per_point,
    )


def _evaluate_rusticol(
    process_dir: Path,
    particles: Sequence[ExternalMomentum],
) -> tuple[float, float]:
    import rusticol as rusticol_module  # type: ignore[import-not-found]

    momenta = np.asarray(
        [
            [
                [float(component) for component in particle.momentum]
                for particle in particles
            ]
        ],
        dtype=np.float64,
    )
    runtime_factory: Any = getattr(rusticol_module, "Runtime")
    runtime = runtime_factory.load(str(process_dir))
    start = time.perf_counter()
    values = runtime.evaluate(momenta)
    return float(values[0]), time.perf_counter() - start


def _concrete_process_from_probe_if_needed(
    particles: Sequence[ExternalMomentum],
    expected_pdgs: Sequence[int],
) -> str | None:
    probe_order = tuple(int(particle.pdg) for particle in particles)
    expected = tuple(int(pdg) for pdg in expected_pdgs)
    if len(probe_order) != len(expected) or probe_order == expected:
        return None
    initial = " ".join(_particle_name_from_pdg(pdg) for pdg in probe_order[:2])
    final = " ".join(_particle_name_from_pdg(pdg) for pdg in probe_order[2:])
    return f"{initial} > {final}"


def _particle_name_from_pdg(pdg: int) -> str:
    names = {int(value): name for name, value in PDGS.items()}
    return names[int(pdg)]


def _reorder_particles_by_pdg(
    particles: Sequence[ExternalMomentum],
    expected_pdgs: Sequence[int],
) -> tuple[ExternalMomentum, ...]:
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
                "Fortran probe point cannot be reordered to the pyAmpliCol "
                f"external PDG order {list(expected_pdgs)}"
            )
    if remaining:
        raise ValueError(
            "Fortran probe point contains extra particles after reordering"
        )
    return tuple(ordered)


def _library_runtime_per_point(
    run: AmplicolWorkflowResult,
    nevents: int,
) -> float | None:
    for row in run.timing_rows:
        if row.label.lower().strip() == "amplitude evaluation":
            return row.seconds / max(1, nevents)
    for row in run.timing_rows:
        if "amplitude" in row.label.lower():
            return row.seconds / max(1, nevents)
    return None


def _relative_difference(reference: float, actual: float) -> float:
    return abs(actual - reference) / max(abs(actual), abs(reference), 1.0e-300)


def _sorted_sector_tuple(values: set[int] | None) -> tuple[int, ...] | None:
    return None if values is None else tuple(sorted(int(value) for value in values))


def format_validation_table(summary: GenericLCValidationSummary) -> str:
    rows = [
        ("Process", "Status", "Rel. diff", "Reference", "pyAmpliCol"),
        ("-" * 28, "-" * 10, "-" * 12, "-" * 12, "-" * 12),
    ]
    for row in summary.rows:
        rows.append(
            (
                row.process,
                row.status,
                _format_float(row.relative_difference),
                _format_float(row.reference_matrix_element),
                _format_float(row.pyamplicol_matrix_element),
            )
        )
    widths = [max(len(str(row[i])) for row in rows) for i in range(5)]
    lines = [
        " | ".join(str(value).ljust(widths[index]) for index, value in enumerate(row))
        for row in rows
    ]
    lines.append(
        "passed="
        + str(summary.passed)
        + ", max_relative_difference="
        + _format_float(summary.max_relative_difference)
    )
    return "\n".join(lines)


def _format_float(value: float | None) -> str:
    return "N/A" if value is None else f"{value:.6g}"


def summary_to_json(summary: GenericLCValidationSummary) -> str:
    return json.dumps(summary.to_json_dict(), indent=2, sort_keys=True)


__all__ = [
    "DEFAULT_LC_VALIDATION_PROCESSES",
    "GenericLCValidationRow",
    "GenericLCValidationSummary",
    "format_validation_table",
    "summary_to_json",
    "validate_generic_lc_processes",
]
