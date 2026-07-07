from __future__ import annotations

import inspect

from pyamplicol.generic_validation import (
    DEFAULT_LC_VALIDATION_PROCESSES,
    _concrete_process_from_probe_if_needed,
    _reorder_particles_by_pdg,
    format_validation_table,
    summary_to_json,
    validate_generic_lc_processes,
)
import pyamplicol.generic_validation as generic_validation
from pyamplicol.core_types import ExternalMomentum
from pyamplicol.process_support import classify_process_support


def test_default_generic_lc_validation_processes_cover_broad_classes() -> None:
    processes = set(DEFAULT_LC_VALIDATION_PROCESSES)

    assert "d d~ > z g g" in processes
    assert "d d~ > e+ e- g g" in processes
    assert "d d~ > mu+ mu- g" in processes
    assert "u d~ > e+ ve g" in processes
    assert "u d~ > e+ ve z" in processes
    assert "u d~ > e+ ve a" in processes
    assert "u d~ > mu+ vm g" in processes
    assert "u d~ > w+ g g" in processes
    assert "d u~ > w-" in processes
    assert "d u~ > w- g" in processes
    assert "d u~ > e- ve~" in processes
    assert "d u~ > e- ve~ g" in processes
    assert "d u~ > e- ve~ z" in processes
    assert "d u~ > e- ve~ a" in processes
    assert "d d~ > z z g" in processes
    assert "d d~ > z z g g" in processes
    assert "d d~ > w+ w- g" in processes
    assert "d d~ > h z" in processes
    assert "d d~ > h z g" in processes
    assert "d d~ > h h z" in processes
    assert "d u~ > h w-" in processes
    assert "d u~ > w- z" in processes
    assert "d u~ > w- a" in processes
    assert "g g > g g" in processes
    assert "g g > g g g" in processes
    assert "g g > u u~ g" in processes
    assert "g g > t t~ g" in processes
    assert "d d~ > t t~" in processes
    assert "d d~ > t t~ g" in processes
    assert "u d~ > u d~ g" in processes
    assert "d d~ > u u~ s s~" in processes


def test_default_generic_lc_validation_processes_are_preflight_supported() -> None:
    reports = {
        process: classify_process_support(
            process,
            color_plan_max_sectors=2000,
            current_plan_max_currents=200000,
            max_quark_pairs=4,
        )
        for process in DEFAULT_LC_VALIDATION_PROCESSES
    }
    unsupported = {
        process: report.artifact_unavailable_message
        for process, report in reports.items()
        if not report.runtime_artifact_supported
    }

    assert unsupported == {}
    assert all(report.current_plan is not None for report in reports.values())
    assert all(
        report.current_plan is not None
        and report.current_plan.full_tensor_network_ready
        for report in reports.values()
    )
    assert all("family" not in report.to_json_dict() for report in reports.values())


def test_generic_lc_validation_dry_run_uses_support_preflight(tmp_path) -> None:
    summary = validate_generic_lc_processes(
        ("d d~ > z g", "u d~ > w+ z"),
        output_dir=tmp_path,
        amplicol_root=tmp_path,
        dry_run=True,
    )

    assert summary.passed is True
    assert summary.max_relative_difference is None
    assert [row.status for row in summary.rows] == ["dry-run", "dry-run"]
    assert [row.supported for row in summary.rows] == [True, True]
    assert summary.rows[0].artifact_dir == tmp_path / "d_dbar_to_z_g"


def test_generic_lc_validation_summary_renderers(tmp_path) -> None:
    summary = validate_generic_lc_processes(
        ("d d~ > z g",),
        output_dir=tmp_path,
        amplicol_root=tmp_path,
        dry_run=True,
    )

    table = format_validation_table(summary)
    payload = summary_to_json(summary)

    assert "d d~ > z g" in table
    assert "passed=True" in table
    assert '"passed": true' in payload
    assert '"fortran_validation_workflow": null' in payload


def test_generic_lc_validation_uses_library_backed_momenta_probe() -> None:
    source = inspect.getsource(generic_validation._validate_one_process)

    assert "adapter.prepare_library(" in source
    assert "warmup_particles=point" in source
    assert "adapter.prepare_direct_probe(" not in source
    assert "use_library=True" in source
    assert "generated_library_momenta_probe" in source


def test_generic_lc_validation_reorders_fortran_probe_particles_by_pdg() -> None:
    particles = (
        ExternalMomentum(1, (1.0, 0.0, 0.0, 1.0)),
        ExternalMomentum(-1, (1.0, 0.0, 0.0, -1.0)),
        ExternalMomentum(21, (1.0, 0.0, 1.0, 0.0)),
        ExternalMomentum(23, (1.0, 0.0, -1.0, 0.0)),
    )

    reordered = _reorder_particles_by_pdg(particles, (1, -1, 23, 21))

    assert [particle.pdg for particle in reordered] == [1, -1, 23, 21]


def test_generic_lc_validation_reconstructs_crossed_probe_subprocess() -> None:
    particles = (
        ExternalMomentum(2, (1.0, 0.0, 0.0, 1.0)),
        ExternalMomentum(-2, (1.0, 0.0, 0.0, -1.0)),
        ExternalMomentum(-1, (1.0, 0.0, 1.0, 0.0)),
        ExternalMomentum(1, (1.0, 0.0, -1.0, 0.0)),
        ExternalMomentum(21, (1.0, 1.0, 0.0, 0.0)),
    )

    concrete = _concrete_process_from_probe_if_needed(
        particles,
        (1, -1, 2, -2, 21),
    )

    assert concrete == "u u~ > d~ d g"


def test_generic_lc_validation_reconstructs_final_state_probe_reordering() -> None:
    particles = (
        ExternalMomentum(2, (1.0, 0.0, 0.0, 1.0)),
        ExternalMomentum(-2, (1.0, 0.0, 0.0, -1.0)),
        ExternalMomentum(-1, (1.0, 0.0, 1.0, 0.0)),
        ExternalMomentum(1, (1.0, 0.0, -1.0, 0.0)),
        ExternalMomentum(21, (1.0, 1.0, 0.0, 0.0)),
    )

    concrete = _concrete_process_from_probe_if_needed(
        particles,
        (2, -2, 1, -1, 21),
    )

    assert concrete == "u u~ > d~ d g"
