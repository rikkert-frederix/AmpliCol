from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_result_matrix_module():
    path = Path(__file__).resolve().parents[2] / "docs" / "result_matrix.py"
    spec = importlib.util.spec_from_file_location("result_matrix_docs", path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_fixed_helicity_search_rejects_numerically_zero_pure_gluon_root() -> None:
    result_matrix = _load_result_matrix_module()
    base = next(
        item
        for item in result_matrix.BASE_PROCESSES
        if item.key == "gg_gluons"
    )

    result = result_matrix._fixed_helicity_choice(
        "g g > g g g g g g",
        base,
    )

    assert result["source_helicities"] == {
        label: -1 for label in range(1, 9)
    }
    assert result["replay_helicity_invariant"] is True


def test_memory_monitor_polling_accelerates_near_the_watchdog_limit() -> None:
    result_matrix = _load_result_matrix_module()

    assert (
        result_matrix._memory_poll_interval(None, high_watermark_bytes=100)
        == result_matrix.MEMORY_POLL_S
    )
    assert (
        result_matrix._memory_poll_interval(99, high_watermark_bytes=100)
        == result_matrix.MEMORY_POLL_S
    )
    assert (
        result_matrix._memory_poll_interval(100, high_watermark_bytes=100)
        == result_matrix.MEMORY_NEAR_LIMIT_POLL_S
    )


def test_shared_lc_recycling_table_is_rendered_from_cache() -> None:
    result_matrix = _load_result_matrix_module()
    data = {
        "entries": {
            "gg_tt_jets": {
                "3": {
                    "amplicol": {
                        "all_flow_generation_s": 0.1,
                        "runtime_us_per_point": 2.0,
                        "all_flow_runtime_us_per_point": 4.0,
                    },
                    "pyamplicol_jit": {
                        "all_flow_generation_s": 0.2,
                        "wall_us_per_point": 3.0,
                        "all_flow_wall_us_per_point": 5.0,
                        "jit_fraction_of_generation": 0.75,
                        "current_count": 29,
                        "amplitude_root_count": 6,
                    },
                    "validation": {"max_relative_difference": 1.0e-12},
                }
            }
        }
    }

    table = result_matrix.render_shared_lc_recycling_table(data)

    assert r"\(gg\to t\bar t+g\)" in table
    assert "6 & 29 & 0.1 / 0.2 & 0.75" in table
    assert r"\(1.50\times\)" in table
    assert r"\(1\times10^{-12}\)" in table


def test_reused_lc_all_flow_fields_do_not_overwrite_selected_measurements() -> None:
    result_matrix = _load_result_matrix_module()
    payload = {
        "selected_generation_s": 12.0,
        "wall_us_per_point": 3.0,
        "generation_s": 12.0,
        "current_count": 10,
    }
    existing = {
        "generation_s": 100.0,
        "current_count": 200,
        "all_flow_generation_s": 100.0,
        "all_flow_wall_us_per_point": 7.0,
        "all_flow_output_dir": "/tmp/preserved-all-flow",
        "wall_us_per_point": 9.0,
        "selected_generation_s": 90.0,
    }

    result_matrix._restore_reused_lc_all_flow_fields(payload, existing)

    assert payload["selected_generation_s"] == 12.0
    assert payload["wall_us_per_point"] == 3.0
    assert payload["generation_s"] == 100.0
    assert payload["current_count"] == 200
    assert payload["all_flow_wall_us_per_point"] == 7.0
    assert payload["all_flow_output_dir"] == "/tmp/preserved-all-flow"


def test_reused_lc_ram_limit_does_not_restore_primary_measurements() -> None:
    result_matrix = _load_result_matrix_module()
    payload = {
        "generation_s": 12.0,
        "current_count": 10,
        "all_flow_status": "ram_limit",
    }
    existing = {
        "generation_s": 100.0,
        "current_count": 200,
        "all_flow_status": "ram_limit",
        "all_flow_output_dir": "/tmp/preserved-all-flow",
        "all_flow_peak_memory_gb": 38.2,
    }

    result_matrix._restore_reused_lc_all_flow_fields(
        payload,
        existing,
        include_primary_fields=False,
    )

    assert payload["generation_s"] == 12.0
    assert payload["current_count"] == 10
    assert payload["all_flow_status"] == "ram_limit"
    assert payload["all_flow_output_dir"] == "/tmp/preserved-all-flow"
    assert payload["all_flow_peak_memory_gb"] == 38.2


def test_lc_cell_renders_all_flow_ram_limit_without_hiding_selected_result() -> None:
    result_matrix = _load_result_matrix_module()
    reference = {
        "status": "ok",
        "generation_s": 10.0,
        "all_flow_generation_s": 20.0,
        "all_flow_generation_source": "amplicol_color_probe_imode2_setup_build",
        "runtime_us_per_point": 100.0,
        "all_flow_runtime_us_per_point": 200.0,
    }
    jit = {
        "status": "ok",
        "selected_generation_s": 5.0,
        "wall_us_per_point": 150.0,
        "all_flow_status": "ram_limit",
    }

    cell = result_matrix._latex_lc_cell({}, reference, jit, {})

    assert r"\matrixratio{speedorange}{1.50}" in cell
    assert cell.count(r">30 GB RAM") == 2


def test_lc_cell_treats_legacy_null_all_flow_status_as_success() -> None:
    result_matrix = _load_result_matrix_module()
    reference = {
        "status": "ok",
        "generation_s": 10.0,
        "all_flow_generation_s": 20.0,
        "all_flow_generation_source": "amplicol_color_probe_imode2_setup_build",
        "runtime_us_per_point": 100.0,
        "all_flow_runtime_us_per_point": 200.0,
    }
    jit = {
        "status": "ok",
        "selected_generation_s": 5.0,
        "all_flow_generation_s": 10.0,
        "wall_us_per_point": 150.0,
        "all_flow_wall_us_per_point": 180.0,
        "all_flow_status": None,
    }

    cell = result_matrix._latex_lc_cell({}, reference, jit, {})

    assert r"\matrixratio{speedgreen}{0.9}" in cell
    assert r"\texttt{none}" not in cell


def test_lc_ram_limit_is_not_hidden_by_unsupported_reference() -> None:
    result_matrix = _load_result_matrix_module()
    case = {
        "amplicol": {"status": "unsupported"},
        "pyamplicol_jit": {
            "status": "ram_limit",
            "error": "memory limit 30 GB exceeded",
        },
    }

    cell = result_matrix._latex_cell(case, color_accuracy="lc")

    assert cell.count(r">30 GB RAM") >= 2
    assert cell != result_matrix._structural_na()


def test_matrix_preserves_generated_output_and_logs(tmp_path: Path) -> None:
    result_matrix = _load_result_matrix_module()
    output_dir = tmp_path / "selected_flow"
    output_dir.mkdir()
    (output_dir / "process_manifest.json").write_text("{}\n")
    generate_log = tmp_path / "selected_flow.generate.log"
    generate_log.write_text("old\n")

    result_matrix._preserve_generated_output(output_dir)

    assert not output_dir.exists()
    preserved_dirs = list(tmp_path.glob("selected_flow_before_*"))
    preserved_logs = list(tmp_path.glob("selected_flow.generate.log.before_*"))
    assert len(preserved_dirs) == 1
    assert len(preserved_logs) == 1
    assert (preserved_dirs[0] / "process_manifest.json").is_file()
    assert preserved_logs[0].read_text() == "old\n"


def test_matrix_retime_preserves_generation_and_validates_both_workloads(
    monkeypatch,
    tmp_path: Path,
) -> None:
    result_matrix = _load_result_matrix_module()
    selected = tmp_path / "jit" / "selected_flow"
    all_flow = tmp_path / "jit" / "all_flows"
    selected.mkdir(parents=True)
    all_flow.mkdir(parents=True)
    (selected / "process_manifest.json").write_text("{}\n")
    (all_flow / "process_manifest.json").write_text("{}\n")
    existing = {
        "status": "ok",
        "generation_s": 12.0,
        "selected_generation_s": 2.0,
        "all_flow_generation_s": 12.0,
        "selected_output_dir": str(selected),
        "all_flow_output_dir": str(all_flow),
        "wall_us_per_point": 8.0,
        "runtime_us_per_point": 7.0,
        "all_flow_status": "ok",
        "all_flow_wall_us_per_point": 10.0,
        "all_flow_runtime_us_per_point": 9.0,
        "time_payload": {"values": [2.0]},
        "all_flow_time_payload": {"values": [3.0]},
        "all_flow_value": 3.0,
        "matrix_settings": {
            "selected_runtime_batch_size": 128,
            "all_flow_runtime_batch_size": 64,
        },
    }
    payloads = iter(
        (
            {
                "values": [2.0],
                "profile": {
                    "wall_us_per_point": 5.0,
                    "core_evaluator_us_per_point": 4.0,
                    "samples": 100,
                },
            },
            {
                "values": [3.0],
                "profile": {
                    "wall_us_per_point": 7.0,
                    "core_evaluator_us_per_point": 6.0,
                    "samples": 200,
                },
            },
        )
    )
    calls: list[tuple[Path, int]] = []

    def fake_time(output_dir, *, target_runtime, batch_size, log_path):
        assert target_runtime == 10.0
        assert log_path.name.endswith(".log")
        calls.append((output_dir, batch_size))
        return next(payloads)

    monkeypatch.setattr(
        result_matrix,
        "_time_existing_pyamplicol_artifact",
        fake_time,
    )

    refreshed = result_matrix._retime_pyamplicol_case(
        existing,
        backend_key="jit",
        output_dir=tmp_path / "jit",
        target_runtime=10.0,
        color_accuracy="lc",
    )

    assert refreshed["generation_s"] == 12.0
    assert refreshed["selected_generation_s"] == 2.0
    assert refreshed["all_flow_generation_s"] == 12.0
    assert refreshed["wall_us_per_point"] == 5.0
    assert refreshed["runtime_us_per_point"] == 4.0
    assert refreshed["all_flow_wall_us_per_point"] == 7.0
    assert refreshed["all_flow_runtime_us_per_point"] == 6.0
    assert refreshed["retime_validation"]["status"] == "ok"
    assert refreshed["all_flow_retime_validation"]["status"] == "ok"
    assert calls == [(selected, 128), (all_flow, 64)]


def test_matrix_retime_missing_artifact_keeps_last_valid_measurement(
    tmp_path: Path,
) -> None:
    result_matrix = _load_result_matrix_module()
    existing = {
        "status": "ok",
        "generation_s": 12.0,
        "wall_us_per_point": 8.0,
        "selected_output_dir": str(tmp_path / "missing"),
    }

    refreshed = result_matrix._retime_pyamplicol_case(
        existing,
        backend_key="jit",
        output_dir=tmp_path / "jit",
        target_runtime=10.0,
        color_accuracy="lc",
    )

    assert refreshed["status"] == "ok"
    assert refreshed["generation_s"] == 12.0
    assert refreshed["wall_us_per_point"] == 8.0
    assert refreshed["retime_status"] == "artifact_missing"


def test_selected_refresh_validation_rejects_changed_value() -> None:
    result_matrix = _load_result_matrix_module()
    previous = {"time_payload": {"values": [2.0]}}

    matching = result_matrix._selected_refresh_validation(
        previous,
        {"time_payload": {"values": [2.0 * (1.0 + 1.0e-12)]}},
    )
    changed = result_matrix._selected_refresh_validation(
        previous,
        {"time_payload": {"values": [2.1]}},
    )

    assert matching["status"] == "ok"
    assert changed["status"] == "failed"
    assert changed["relative_difference"] > changed["relative_tolerance"]


def test_lc_matrix_settings_record_selected_o3_and_all_flow_o1() -> None:
    result_matrix = _load_result_matrix_module()
    base = next(
        item
        for item in result_matrix.BASE_PROCESSES
        if item.key == "gg_tt_jets"
    )

    settings = result_matrix._pyamplicol_matrix_settings(
        backend_key="jit",
        target_runtime=10.0,
        n_cores=5,
        selected_lc_sector_ids={0},
        reference_color_order=(3, 5, 1, 2, 4),
        base=base,
        color_accuracy="lc",
        selected_jit_opt_level=3,
    )

    assert settings["symbolica_jit_optimization_level"] == 3
    assert settings["selected_symbolica_jit_optimization_level"] == 3
    assert settings["all_flow_symbolica_jit_optimization_level"] == 1


def test_selected_o3_runtime_ratio_is_starred() -> None:
    result_matrix = _load_result_matrix_module()
    mode = {
        "status": "ok",
        "wall_us_per_point": 1.5,
        "all_flow_wall_us_per_point": 2.0,
        "selected_jit_optimization_level": 3,
    }

    selected = result_matrix._runtime_ratio_pair(mode, 1.0)
    all_flow = result_matrix._runtime_ratio_pair(
        mode,
        1.0,
        wall_key="all_flow_wall_us_per_point",
    )

    assert selected.endswith(r"\textsuperscript{*}")
    assert r"\textsuperscript{*}" not in all_flow


def test_selected_o3_generation_ratio_is_starred() -> None:
    result_matrix = _load_result_matrix_module()
    mode = {
        "status": "ok",
        "selected_jit_optimization_level": 3,
    }

    selected = result_matrix._generation_ratio_or_absolute_from_value(
        mode,
        1.5,
        1.0,
        mark_selected_jit_level=True,
    )
    all_flow = result_matrix._generation_ratio_or_absolute_from_value(
        mode,
        2.0,
        1.0,
    )

    assert selected.endswith(r"\textsuperscript{*}")
    assert r"\textsuperscript{*}" not in all_flow


def test_lc_matrix_notes_match_runtime_batches_and_native_path() -> None:
    result_matrix = _load_result_matrix_module()

    short_intro = " ".join(result_matrix._matrix_short_intro_latex("lc"))
    text = " ".join(result_matrix._matrix_run_settings_latex("lc"))

    assert result_matrix._matrix_title("lc").endswith("(SymJIT O1)")
    assert "O1 selected flow" in short_intro
    assert "starred" not in short_intro
    assert "O1 all flows" in short_intro
    assert "batches are 128/64" in text
    assert "native two-lane complex evaluators" in text
    assert r"\(n\geq5\)" in text


def test_lc_structural_na_is_kept_in_compact_table_notes() -> None:
    result_matrix = _load_result_matrix_module()
    entries = {
        "dd_4q_lines": {
            "7": {
                "amplicol": {"status": "unsupported"},
                "pyamplicol_jit": {"status": "ok"},
            }
        }
    }

    intro = " ".join(result_matrix._matrix_long_intro_latex("lc", {}))
    status = " ".join(
        result_matrix._matrix_status_notes_latex(
            entries,
            [7],
            color_accuracy="lc",
        )
    )

    assert r"\PAC-only cells are absolute" in intro
    assert "Structural N/A" not in status


def test_lc_ram_limit_is_kept_in_compact_table_notes() -> None:
    result_matrix = _load_result_matrix_module()
    entries = {
        "gg_gluons": {
            "8": {
                "amplicol": {"status": "ok"},
                "pyamplicol_jit": {
                    "status": "ok",
                    "all_flow_status": "ram_limit",
                },
            }
        }
    }

    intro = " ".join(result_matrix._matrix_long_intro_latex("lc", {}))
    status = " ".join(
        result_matrix._matrix_status_notes_latex(
            entries,
            [8],
            color_accuracy="lc",
        )
    )

    assert r"\texttt{>30 GB RAM}" in intro
    assert "RAM gaps" not in status


def test_not_applicable_rows_only_trigger_final_writer(
    monkeypatch,
    tmp_path: Path,
) -> None:
    result_matrix = _load_result_matrix_module()
    writes: list[tuple[object, ...]] = []

    monkeypatch.setattr(
        result_matrix,
        "_write_table_data_and_maybe_pdf",
        lambda *args, **kwargs: writes.append((*args, kwargs)),
    )

    assert (
        result_matrix.main(
            [
                "--data",
                str(tmp_path / "data.json"),
                "--table",
                str(tmp_path / "table.tex"),
                "--generate-data",
                "1",
                "--show-data",
                "1",
                "--base-process",
                "dd_4q_lines",
                "--skip-amplicol",
                "--skip-jit",
                "--skip-cpp-o3",
                "--no-recompile",
            ]
        )
        == 0
    )

    assert len(writes) == 1
