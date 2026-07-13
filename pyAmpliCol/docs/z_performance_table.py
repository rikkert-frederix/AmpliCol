#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

DOCS_DIR = Path(__file__).resolve().parent
DEFAULT_DATA = DOCS_DIR / "z_performance_data.json"
DEFAULT_TABLE = DOCS_DIR / "z_performance_table.tex"
DEFAULT_TARGET_RUNTIME_S = 10.0
DEFAULT_BATCH_SIZE = 64
DEFAULT_ALL_FLOW_BATCH_SIZE = 64
DEFAULT_SYMBOLICA_ITERATIONS = 10
DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE = 128
DEFAULT_ALL_FLOW_SYMBOLICA_OUTPUT_CHUNK_SIZE = 8192

MODES: tuple[dict[str, Any], ...] = (
    {
        "key": "amplicol",
        "route": "reference",
        "setup": r"\AC",
        "row_color": "refblue",
    },
    {
        "key": "jit_o1",
        "route": "Rusticol",
        "setup": r"\PAC\ JIT \(\mathrm{O}1\)",
        "row_color": "bestgreen",
    },
    {
        "key": "asm",
        "route": "Rusticol",
        "setup": r"\PAC\ ASM \(\mathrm{O}3\)",
    },
    {
        "key": "cpp_o3",
        "route": "Rusticol",
        "setup": r"\PAC\ C++ \(\mathrm{O}3\)",
    },
    {
        "key": "jit_o3",
        "route": "Rusticol",
        "setup": r"\PAC\ JIT \(\mathrm{O}3\)",
    },
)
MODE_KEYS = tuple(str(mode["key"]) for mode in MODES)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Render or update the manually-steered d d~ > Z + gluons "
            "performance table used by pyAmpliCol.tex."
        )
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument(
        "--built-in-data",
        type=Path,
        default=DEFAULT_DATA,
        help=(
            "Built-in-SM Z-table cache used for external-model generation- and "
            "wall-time ratios."
        ),
    )
    parser.add_argument("--model-source", default="BUILTIN_SM")
    parser.add_argument(
        "--model-profile",
        choices=("built-in-sm", "external-sm"),
        default="built-in-sm",
    )
    parser.add_argument("--model-label", default="built-in-sm")
    parser.add_argument("--no-recompile", action="store_true")
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("reset", help="Reset the cache and render an empty table.")
    subparsers.add_parser("render", help="Render the table from the current cache.")

    record = subparsers.add_parser(
        "record",
        help="Record one manually measured table row and refresh the table/PDF.",
    )
    record.add_argument("--n", type=int, required=True, choices=range(1, 10))
    record.add_argument("--mode", choices=MODE_KEYS, required=True)
    record.add_argument(
        "--status",
        choices=("ok", "timeout", "ram_limit", "error", "not_run"),
        default="ok",
    )
    record.add_argument("--generation-s", type=float)
    record.add_argument("--output-dir", default="")
    record.add_argument("--wall-us-per-point", type=float)
    record.add_argument("--runtime-us-per-point", type=float)
    record.add_argument("--generation-measurement-json", default="")
    record.add_argument(
        "--all-flow-status",
        choices=("ok", "missing", "timeout", "ram_limit", "error", "not_run"),
        default=None,
    )
    record.add_argument("--all-flow-generation-s", type=float)
    record.add_argument("--all-flow-output-dir", default="")
    record.add_argument("--all-flow-wall-us-per-point", type=float)
    record.add_argument("--all-flow-runtime-us-per-point", type=float)
    record.add_argument("--all-flow-generation-measurement-json", default="")
    record.add_argument("--all-flow-time-batch-size", type=int)
    record.add_argument("--all-flow-symbolica-output-chunk-size", type=int)
    record.add_argument("--all-flow-notes", default="")
    record.add_argument("--all-flow-error", default="")
    record.add_argument("--notes", default="")
    record.add_argument("--error", default="")

    args = parser.parse_args(argv)
    command = args.command or "render"
    if command == "reset":
        data: dict[str, Any] = {}
    else:
        data = _load_data(args.data)
    _refresh_metadata(
        data,
        model_source=str(args.model_source),
        model_profile=str(args.model_profile),
        model_label=str(args.model_label),
    )

    if command == "record":
        _record_row(
            data,
            int(args.n),
            str(args.mode),
            status=str(args.status),
            generation_s=args.generation_s,
            output_dir=str(args.output_dir),
            wall_us_per_point=args.wall_us_per_point,
            runtime_us_per_point=args.runtime_us_per_point,
            generation_measurement=_json_object_argument(
                args.generation_measurement_json,
                option="--generation-measurement-json",
            ),
            all_flow_status=args.all_flow_status,
            all_flow_generation_s=args.all_flow_generation_s,
            all_flow_output_dir=str(args.all_flow_output_dir),
            all_flow_wall_us_per_point=args.all_flow_wall_us_per_point,
            all_flow_runtime_us_per_point=args.all_flow_runtime_us_per_point,
            all_flow_generation_measurement=_json_object_argument(
                args.all_flow_generation_measurement_json,
                option="--all-flow-generation-measurement-json",
            ),
            all_flow_time_batch_size=args.all_flow_time_batch_size,
            all_flow_symbolica_output_chunk_size=args.all_flow_symbolica_output_chunk_size,
            all_flow_notes=str(args.all_flow_notes),
            all_flow_error=str(args.all_flow_error),
            notes=str(args.notes),
            error=str(args.error),
        )

    _write_data(args.data, data)
    built_in_data = (
        _load_data(args.built_in_data)
        if str(args.model_profile) != "built-in-sm"
        else None
    )
    table = render_table(data, built_in_data=built_in_data)
    args.table.write_text(table, encoding="utf-8")
    print(f"wrote {args.table}")
    print(f"wrote {args.data}")
    if not args.no_recompile:
        _refresh_pdf()
    return 0


def _load_data(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    text = path.read_text(encoding="utf-8")
    if not text.strip():
        return {}
    value = json.loads(text)
    if not isinstance(value, dict):
        raise TypeError(f"{path} must contain a JSON object")
    return value


def _json_object_argument(text: str, *, option: str) -> dict[str, Any] | None:
    if not text:
        return None
    value = json.loads(text)
    if not isinstance(value, dict):
        raise TypeError(f"{option} must encode a JSON object")
    return value


def _write_data(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _refresh_metadata(
    data: dict[str, Any],
    *,
    model_source: str,
    model_profile: str,
    model_label: str,
) -> None:
    data["schema_version"] = 1
    data["created_by"] = "pyAmpliCol/docs/z_performance_table.py"
    data["updated_at"] = datetime.now(timezone.utc).isoformat()
    data["multiplicity_convention"] = "final_state"
    data["process_family"] = "d d~ > z + (n-1)g"
    data["model_source"] = model_source
    data["model_profile"] = model_profile
    data["model_label"] = model_label
    data["settings"] = {
        "target_runtime_s": DEFAULT_TARGET_RUNTIME_S,
        "batch_size": DEFAULT_BATCH_SIZE,
        "selected_runtime_batch_size": DEFAULT_BATCH_SIZE,
        "all_flow_runtime_batch_size": DEFAULT_ALL_FLOW_BATCH_SIZE,
        "symbolica_output_chunk_size": DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
        "selected_symbolica_output_chunk_size": DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
        "all_flow_symbolica_output_chunk_size": (
            DEFAULT_ALL_FLOW_SYMBOLICA_OUTPUT_CHUNK_SIZE
        ),
        "symbolica_stage_local_parameter_layout": True,
        "symbolica_iterations": DEFAULT_SYMBOLICA_ITERATIONS,
        "symbolica_cpe_iterations": None,
        "symbolica_max_horner_scheme_variables": 1000,
        "symbolica_max_common_pair_cache_entries": 5_000_000,
        "symbolica_max_common_pair_distance": 1000,
    }
    data.setdefault("entries", {})


def _record_row(
    data: dict[str, Any],
    n_final: int,
    mode_key: str,
    *,
    status: str,
    generation_s: float | None,
    output_dir: str,
    wall_us_per_point: float | None,
    runtime_us_per_point: float | None,
    generation_measurement: dict[str, Any] | None,
    all_flow_status: str | None,
    all_flow_generation_s: float | None,
    all_flow_output_dir: str,
    all_flow_wall_us_per_point: float | None,
    all_flow_runtime_us_per_point: float | None,
    all_flow_generation_measurement: dict[str, Any] | None,
    all_flow_time_batch_size: int | None,
    all_flow_symbolica_output_chunk_size: int | None,
    all_flow_notes: str,
    all_flow_error: str,
    notes: str,
    error: str,
) -> None:
    entries = data.setdefault("entries", {})
    if not isinstance(entries, dict):
        raise TypeError("'entries' must be a JSON object")
    case = entries.setdefault(str(n_final), {})
    if not isinstance(case, dict):
        raise TypeError(f"entry for n={n_final} must be a JSON object")
    case["process"] = process_for_n(n_final)
    modes = case.setdefault("modes", {})
    if not isinstance(modes, dict):
        raise TypeError(f"modes for n={n_final} must be a JSON object")
    row: dict[str, Any] = {
        "status": status,
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }
    if generation_s is not None:
        row["generation_s"] = float(generation_s)
    if output_dir:
        row["output_dir"] = output_dir
    if wall_us_per_point is not None:
        row["wall_us_per_point"] = float(wall_us_per_point)
    if runtime_us_per_point is not None:
        row["runtime_us_per_point"] = float(runtime_us_per_point)
    if generation_measurement is not None:
        row["generation_measurement"] = generation_measurement
    if all_flow_status is not None:
        row["all_flow_status"] = str(all_flow_status)
    if all_flow_generation_s is not None:
        row["all_flow_generation_s"] = float(all_flow_generation_s)
    if all_flow_output_dir:
        row["all_flow_output_dir"] = all_flow_output_dir
    if all_flow_wall_us_per_point is not None:
        row["all_flow_wall_us_per_point"] = float(all_flow_wall_us_per_point)
    if all_flow_runtime_us_per_point is not None:
        row["all_flow_runtime_us_per_point"] = float(all_flow_runtime_us_per_point)
    if all_flow_generation_measurement is not None:
        row["all_flow_generation_measurement"] = all_flow_generation_measurement
    if all_flow_time_batch_size is not None:
        row["all_flow_time_batch_size"] = int(all_flow_time_batch_size)
    if all_flow_symbolica_output_chunk_size is not None:
        row["all_flow_symbolica_output_chunk_size"] = int(
            all_flow_symbolica_output_chunk_size
        )
    if all_flow_notes:
        row["all_flow_notes"] = all_flow_notes
    if all_flow_error:
        row["all_flow_error"] = all_flow_error
    if notes:
        row["notes"] = notes
    if error:
        row["error"] = error
    modes[mode_key] = row


def process_for_n(n_final: int) -> str:
    return f"d d~ > Z + ({n_final}-1)*g"


def _reproduction_n(entries: dict[str, Any]) -> int | None:
    preferred = (7, 9, 8, 6, 5, 4, 3, 2, 1)
    for n_final in preferred:
        case = entries.get(str(n_final), {})
        if not isinstance(case, dict):
            continue
        modes = case.get("modes", {})
        if not isinstance(modes, dict):
            continue
        jit_o1 = modes.get("jit_o1", {})
        if isinstance(jit_o1, dict) and jit_o1.get("status") == "ok":
            return n_final
    return None


def _explicit_process_for_n(n_final: int) -> str:
    final_state = ["z", *(["g"] * (n_final - 1))]
    return "d d~ > " + " ".join(final_state)


def _reference_color_order_for_n(n_final: int) -> str:
    order = [2, *range(4, n_final + 3), 1, 3]
    return ",".join(str(label) for label in order)


def _display_n_values(
    entries: dict[str, Any],
    *,
    external_model: bool,
) -> tuple[int, ...]:
    if not external_model:
        return tuple(range(1, 10))
    populated = [
        n_final
        for n_final in range(1, 10)
        if isinstance(entries.get(str(n_final)), dict)
        and isinstance(entries[str(n_final)].get("modes"), dict)
        and any(
            isinstance(row, dict)
            and (row.get("status") == "ok" or row.get("all_flow_status") == "ok")
            for row in entries[str(n_final)]["modes"].values()
        )
    ]
    return tuple(range(1, max(populated, default=1) + 1))


def render_table(
    data: dict[str, Any],
    *,
    built_in_data: dict[str, Any] | None = None,
) -> str:
    entries = data.get("entries", {})
    if not isinstance(entries, dict):
        entries = {}
    model_profile = str(data.get("model_profile", "built-in-sm"))
    model_label = str(data.get("model_label", "built-in-sm"))
    external_prefix = "External-SM " if model_profile != "built-in-sm" else ""
    compare_to_built_in = model_profile != "built-in-sm"
    built_in_entries = (
        built_in_data.get("entries", {}) if isinstance(built_in_data, dict) else {}
    )
    if not isinstance(built_in_entries, dict):
        built_in_entries = {}
    display_n_values = _display_n_values(
        entries,
        external_model=compare_to_built_in,
    )
    output_root = (
        "docs/.z_performance_ufo_sm_outputs"
        if model_profile != "built-in-sm"
        else "docs/.z_performance_outputs"
    )
    model_option = (
        "--model assets/models/json/sm/sm.json "
        if model_profile != "built-in-sm"
        else ""
    )
    lines = [
        "% Generated by docs/z_performance_table.py; edit z_performance_data.json instead.",
        r"\begin{landscape}",
        (
            rf"\section{{\texorpdfstring{{{external_prefix}Dedicated \(d\bar d\to Z+\) Gluon "
            rf"Performance}}{{{external_prefix}Dedicated d dbar to Z plus Gluon Performance}}}}"
        ),
        r"\begingroup",
        r"\tiny",
        r"\setlength{\tabcolsep}{1.55pt}",
        r"\renewcommand{\arraystretch}{1.03}",
        (
            r"\begin{longtable}{@{}r "
            + (
                r"L{0.78in} L{0.40in} L{0.94in} "
                r"r L{0.42in} r L{0.42in} r "
                r"@{}p{0.10in}@{} "
                r"r L{0.42in} r L{0.42in} r@{}}"
                if compare_to_built_in
                else (
                    r"L{0.92in} L{0.48in} L{0.86in} r r r "
                    r"@{}p{0.10in}@{} r r r@{}}"
                )
            )
        ),
        r"\toprule",
        (
            r"\textbf{n} & \textbf{process} & \textbf{route} & \textbf{setup} "
            + rf"& \multicolumn{{{5 if compare_to_built_in else 3}}}{{c}}{{\textbf{{selected flow, helicity sum}}}} "
            + r"& "
            + rf"& \multicolumn{{{5 if compare_to_built_in else 3}}}{{c}}{{\textbf{{all flows, fixed helicity}}}} "
            r"\\"
        ),
        (
            r"& & & & \textbf{gen [s]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{wall [us/pt]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{eval [us/pt]} & & \textbf{gen [s]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{wall [us/pt]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{eval [us/pt]} \\"
        ),
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        (
            r"\textbf{n} & \textbf{process} & \textbf{route} & \textbf{setup} "
            + rf"& \multicolumn{{{5 if compare_to_built_in else 3}}}{{c}}{{\textbf{{selected flow, helicity sum}}}} "
            + r"& "
            + rf"& \multicolumn{{{5 if compare_to_built_in else 3}}}{{c}}{{\textbf{{all flows, fixed helicity}}}} "
            r"\\"
        ),
        (
            r"& & & & \textbf{gen [s]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{wall [us/pt]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{eval [us/pt]} & & \textbf{gen [s]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{wall [us/pt]} "
            + (r"& \textbf{vs blt-in} " if compare_to_built_in else "")
            + r"& \textbf{eval [us/pt]} \\"
        ),
        r"\midrule",
        r"\endhead",
    ]
    for n_final in display_n_values:
        case = entries.get(str(n_final), {})
        if not isinstance(case, dict):
            case = {}
        process = process_for_n(n_final)
        modes = case.get("modes", {})
        if not isinstance(modes, dict):
            modes = {}
        reference = modes.get("amplicol", {})
        if not isinstance(reference, dict):
            reference = {}
        ref_generation = _optional_float(reference.get("generation_s"))
        ref_runtime = _optional_float(reference.get("runtime_us_per_point"))
        ref_all_flow_generation = _optional_float(reference.get("all_flow_generation_s"))
        ref_all_flow_runtime = _optional_float(reference.get("all_flow_runtime_us_per_point"))
        built_in_case = built_in_entries.get(str(n_final), {})
        if not isinstance(built_in_case, dict):
            built_in_case = {}
        built_in_modes = built_in_case.get("modes", {})
        if not isinstance(built_in_modes, dict):
            built_in_modes = {}
        for mode in MODES:
            mode_key = str(mode["key"])
            row = modes.get(mode_key, {})
            if not isinstance(row, dict):
                row = {}
            row_color = mode.get("row_color")
            if row_color:
                lines.append(rf"\rowcolor{{{row_color}}}")
            cells = _render_mode_cells(
                row,
                mode_key=mode_key,
                ref_generation=ref_generation,
                ref_runtime=ref_runtime,
                ref_all_flow_generation=ref_all_flow_generation,
                ref_all_flow_runtime=ref_all_flow_runtime,
            )
            if compare_to_built_in:
                built_in_row = built_in_modes.get(mode_key, {})
                if not isinstance(built_in_row, dict):
                    built_in_row = {}
                cells = [
                    cells[0],
                    _ratio_against_built_in(
                        row,
                        built_in_row,
                        mode_key=mode_key,
                        status_key="status",
                        value_key="generation_s",
                    ),
                    cells[1],
                    _ratio_against_built_in(
                        row,
                        built_in_row,
                        mode_key=mode_key,
                        status_key="status",
                        value_key="wall_us_per_point",
                    ),
                    cells[2],
                    cells[3],
                    _ratio_against_built_in(
                        row,
                        built_in_row,
                        mode_key=mode_key,
                        status_key="all_flow_status",
                        value_key="all_flow_generation_s",
                    ),
                    cells[4],
                    _ratio_against_built_in(
                        row,
                        built_in_row,
                        mode_key=mode_key,
                        status_key="all_flow_status",
                        value_key="all_flow_wall_us_per_point",
                    ),
                    cells[5],
                ]
            cells = (
                [*cells[:5], "", *cells[5:10]]
                if compare_to_built_in
                else [*cells[:3], "", *cells[3:6]]
            )
            lines.append(
                " & ".join(
                    [
                        rf"\textbf{{{n_final}}}" if mode_key == "amplicol" else str(n_final),
                        rf"\texttt{{{_latex_escape(process)}}}",
                        str(mode["route"]),
                        str(mode["setup"]),
                        *cells,
                    ]
                )
                + r" \\"
            )
        if n_final != display_n_values[-1]:
            lines.append(r"\midrule[0.45pt]")
    lines.extend([r"\bottomrule", r"\end{longtable}", r"\endgroup"])
    reproduction_n = _reproduction_n(entries)
    if reproduction_n is not None:
        process = _explicit_process_for_n(reproduction_n)
        reference_order = _reference_color_order_for_n(reproduction_n)
        output_dir = f"{output_root}/manual/n{reproduction_n}/jit_o1"
        lines.extend(
            [
                r"\par\smallskip",
                (
                    rf"\noindent\footnotesize Reproduce the \(n={reproduction_n}\) "
                    r"JIT O1 selected-flow entry with the following generation and "
                    r"timing commands."
                ),
                r"\begin{lstlisting}[language=bash,basicstyle=\ttfamily\tiny,breaklines=true,breakatwhitespace=true]",
                (
                    "env PYTHONPATH=src dependencies/.venv/bin/python -m pyamplicol generate-process "
                    + model_option
                    + f"'{process}' {output_dir} "
                    "--replace --n_cores 5 --color-accuracy lc --batch-size 64 "
                    "--symbolica-n-cores 5 --symbolica-output-chunk-size 128 "
                    "--symbolica-output-chunk-strategy uniform "
                    "--symbolica-stage-local-parameter-layout --symbolica-iterations 10 "
                    "--symbolica-max-horner-scheme-variables 1000 "
                    "--symbolica-max-common-pair-cache-entries 5000000 "
                    "--symbolica-max-common-pair-distance 1000 --lc-sector-ids 0 "
                    f"--reference-color-order {reference_order} "
                    "--symbolica-evaluator-backend jit --symbolica-jit-optimization-level 1"
                ),
                (
                    "env PYTHONPATH=src dependencies/.venv/bin/python -m pyamplicol time-process "
                    "--target-runtime 10 --batch-size 64 --json "
                    + output_dir
                ),
                r"\end{lstlisting}",
            ]
        )
    lines.extend(
        [
            r"\par\smallskip",
            (
                rf"\noindent\footnotesize \PAC\ model source: \texttt{{{_latex_escape(model_label)}}}. "
                r"Here \(n\) is final-state multiplicity. "
                r"Each block reports generation seconds and wall/evaluator "
                r"microseconds per point: one selected flow with the helicity sum "
                r"on the left, all flows at one fixed helicity on the right. "
                r"Wall measures Rusticol's public \texttt{evaluate} call; eval "
                r"separately profiles evaluator plus input packing.  Both use batch "
                r"64, stage-local inputs, ten Horner iterations, and a ten-second "
                r"timing target; output chunks are 128 selected and 8192 all-flow. "
                r"The green O1 row is the measured default: repeated O3 timings "
                r"showed no significant Z-family gain.  ASM remains a scalar "
                r"AArch64 diagnostic.  The driver keeps generated artifacts under "
                + rf"\texttt{{{_latex_escape(output_root)}}}."
                + (
                    r"  \AC's generated-library time is reported under wall; its "
                    r"separate eval metric is N/A."
                )
                + (
                    r"  Its \texttt{vs blt-in} cells are left blank."
                    if compare_to_built_in
                    else ""
                )
                + (
                    r"  In this external-model table, each \texttt{vs blt-in} entry "
                    r"is the adjacent generation- or wall-time ratio to the matching "
                    r"built-in-SM row at the same multiplicity, backend, and "
                    r"flow/helicity workload; values below one mean that the external "
                    r"model is faster."
                    if compare_to_built_in
                    else ""
                )
            ),
            r"\end{landscape}",
            "",
        ]
    )
    return "\n".join(lines)


def _render_mode_cells(
    row: dict[str, Any],
    *,
    mode_key: str,
    ref_generation: float | None,
    ref_runtime: float | None,
    ref_all_flow_generation: float | None,
    ref_all_flow_runtime: float | None,
) -> list[str]:
    status = str(row.get("status", "missing"))
    all_flow_status = str(row.get("all_flow_status", "missing"))
    if status == "timeout" and mode_key == "cpp_o3":
        return [
            r"\textcolor{speedred}{\texttt{t/o >15 min}}",
            _missing(),
            _missing(),
            *_render_timing_triplet(
                row,
                mode_key=mode_key,
                status=all_flow_status,
                generation_key="all_flow_generation_s",
                wall_key="all_flow_wall_us_per_point",
                runtime_key="all_flow_runtime_us_per_point",
                ref_generation=ref_all_flow_generation,
                ref_runtime=ref_all_flow_runtime,
                notes_key="all_flow_notes",
                error_key="all_flow_error",
            ),
            _notes(row),
        ]
    if status == "ram_limit":
        return [
            r"\textcolor{speedred}{\texttt{>30 GB RAM}}",
            _missing(color="black!45") if mode_key == "amplicol" else _missing(),
            _missing(),
            *_render_timing_triplet(
                row,
                mode_key=mode_key,
                status=all_flow_status,
                generation_key="all_flow_generation_s",
                wall_key="all_flow_wall_us_per_point",
                runtime_key="all_flow_runtime_us_per_point",
                ref_generation=ref_all_flow_generation,
                ref_runtime=ref_all_flow_runtime,
                notes_key="all_flow_notes",
                error_key="all_flow_error",
            ),
            _notes(row),
        ]
    if status not in {"ok", "missing"}:
        return [
            rf"\textcolor{{speedred}}{{\texttt{{{_latex_escape(status)}}}}}",
            _missing(color="black!45") if mode_key == "amplicol" else _missing(),
            _missing(),
            *_render_timing_triplet(
                row,
                mode_key=mode_key,
                status=all_flow_status,
                generation_key="all_flow_generation_s",
                wall_key="all_flow_wall_us_per_point",
                runtime_key="all_flow_runtime_us_per_point",
                ref_generation=ref_all_flow_generation,
                ref_runtime=ref_all_flow_runtime,
                notes_key="all_flow_notes",
                error_key="all_flow_error",
            ),
            _notes(row),
        ]
    if status == "missing":
        if mode_key == "amplicol":
            selected = [_missing(), _missing(color="black!45"), _missing()]
        else:
            selected = [_missing(), _missing(), _missing()]
        return [
            *selected,
            *_render_timing_triplet(
                row,
                mode_key=mode_key,
                status=all_flow_status,
                generation_key="all_flow_generation_s",
                wall_key="all_flow_wall_us_per_point",
                runtime_key="all_flow_runtime_us_per_point",
                ref_generation=ref_all_flow_generation,
                ref_runtime=ref_all_flow_runtime,
                notes_key="all_flow_notes",
                error_key="all_flow_error",
            ),
            "",
        ]
    generation = _optional_float(row.get("generation_s"))
    wall = _optional_float(row.get("wall_us_per_point"))
    runtime = _optional_float(row.get("runtime_us_per_point"))
    if mode_key == "amplicol":
        selected = [
            _format_plain(generation),
            _format_plain(runtime),
            _missing(color="black!45"),
        ]
    else:
        selected = [
            _format_with_ratio(generation, ref_generation),
            _format_with_ratio(wall, ref_runtime),
            _format_with_ratio(runtime, ref_runtime),
        ]
    return [
        *selected,
        *_render_timing_triplet(
            row,
            mode_key=mode_key,
            status=all_flow_status,
            generation_key="all_flow_generation_s",
            wall_key="all_flow_wall_us_per_point",
            runtime_key="all_flow_runtime_us_per_point",
            ref_generation=ref_all_flow_generation,
            ref_runtime=ref_all_flow_runtime,
            notes_key="all_flow_notes",
            error_key="all_flow_error",
        ),
        _joined_notes(row),
    ]


def _render_timing_triplet(
    row: dict[str, Any],
    *,
    mode_key: str,
    status: str,
    generation_key: str,
    wall_key: str,
    runtime_key: str,
    ref_generation: float | None,
    ref_runtime: float | None,
    notes_key: str,
    error_key: str,
) -> list[str]:
    if status == "missing":
        if mode_key == "amplicol":
            return [_missing(), _missing(color="black!45"), _missing()]
        return [_missing(), _missing(), _missing()]
    if status == "timeout" and mode_key == "cpp_o3":
        return [r"\textcolor{speedred}{\texttt{t/o >15 min}}", _missing(), _missing()]
    if status == "ram_limit":
        return [
            r"\textcolor{speedred}{\texttt{>30 GB RAM}}",
            _missing(color="black!45") if mode_key == "amplicol" else _missing(),
            _missing(),
        ]
    if status != "ok":
        return [
            rf"\textcolor{{speedred}}{{\texttt{{{_latex_escape(status)}}}}}",
            _missing(color="black!45") if mode_key == "amplicol" else _missing(),
            _missing(),
        ]
    generation = _optional_float(row.get(generation_key))
    wall = _optional_float(row.get(wall_key))
    runtime = _optional_float(row.get(runtime_key))
    if mode_key == "amplicol":
        return [
            _format_plain(generation),
            _format_plain(runtime),
            _missing(color="black!45"),
        ]
    return [
        _format_with_ratio(generation, ref_generation),
        _format_with_ratio(wall, ref_runtime),
        _format_with_ratio(runtime, ref_runtime),
    ]


def _ratio_against_built_in(
    row: dict[str, Any],
    built_in_row: dict[str, Any],
    *,
    mode_key: str,
    status_key: str,
    value_key: str,
) -> str:
    if mode_key == "amplicol":
        return ""
    if row.get(status_key) != "ok" or built_in_row.get(status_key) != "ok":
        return _missing(color="black!45")
    value = _optional_float(row.get(value_key))
    reference = _optional_float(built_in_row.get(value_key))
    if value is None or reference is None or reference <= 0.0:
        return _missing(color="black!45")
    ratio = value / reference
    color = "speedgreen" if ratio < 1.0 else "speedorange" if ratio < 2.0 else "speedred"
    return rf"\textcolor{{{color}}}{{\texttt{{x{_format_ratio(ratio)}}}}}"


def _joined_notes(row: dict[str, Any]) -> str:
    notes = [note for note in (_notes(row), _notes(row, notes_key="all_flow_notes", error_key="all_flow_error")) if note]
    return "; ".join(notes)


def _notes(
    row: dict[str, Any],
    *,
    notes_key: str = "notes",
    error_key: str = "error",
) -> str:
    notes = str(row.get(notes_key, ""))
    error = str(row.get(error_key, ""))
    if notes:
        if (
            notes.startswith("generated process kept at ")
            or notes.startswith("retimed existing process kept at ")
            or notes.startswith("reused matching O1 all-flow artifact from ")
            or notes == "reused from LC result matrix cache"
        ):
            return ""
        notes = notes.replace(
            "generated process kept at pyAmpliCol/docs/.z_performance_outputs/",
            "kept: .z_performance_outputs/",
        )
        notes = notes.replace(
            "reused from LC result matrix cache",
            "reused from LC matrix cache",
        )
        return _latex_escape(notes)
    if error:
        return _latex_escape(error[-220:])
    return ""


def _format_plain(value: float | None) -> str:
    if value is None:
        return _missing()
    return rf"\texttt{{{_format_number(value)}}}"


def _format_with_ratio(value: float | None, reference: float | None) -> str:
    if value is None:
        return _missing()
    text = _format_plain(value)
    if reference is None or reference <= 0.0:
        return text
    ratio = value / reference
    color = "speedgreen" if ratio < 1.0 else "speedorange" if ratio < 2.0 else "speedred"
    return text + rf" \textcolor{{{color}}}{{\texttt{{(x{_format_ratio(ratio)})}}}}"


def _format_number(value: float) -> str:
    text = f"{value:.3g}"
    return text.replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")


def _format_ratio(value: float) -> str:
    if value < 0.095:
        return f"{value:.2f}"
    if value < 9.95:
        return f"{value:.2f}"
    return _format_number(value)


def _missing(*, color: str = "speedred") -> str:
    return rf"\textcolor{{{color}}}{{\texttt{{N/A}}}}"


def _optional_float(value: object) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _latex_escape(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(char, char) for char in value)


def _refresh_pdf() -> None:
    document = DOCS_DIR / "pyAmpliCol.tex"
    if not document.exists():
        print(f"skipped PDF refresh: {document} does not exist", flush=True)
        return
    env = dict(os.environ)
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    print("refreshing pyAmpliCol/docs/pyAmpliCol.pdf", flush=True)
    completed = subprocess.run(
        [
            "latexmk",
            "-pdf",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "pyAmpliCol.tex",
        ],
        cwd=DOCS_DIR,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "pyAmpliCol.tex recompilation failed with exit code "
            f"{completed.returncode}\nstdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    print(f"wrote {DOCS_DIR / 'pyAmpliCol.pdf'}")


if __name__ == "__main__":
    raise SystemExit(main())
