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
DEFAULT_SYMBOLICA_ITERATIONS = 10
DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE = 128

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
        "row_color": "bestgreen",
    },
)
MODE_KEYS = tuple(str(mode["key"]) for mode in MODES)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Render or update the manually-steered d d~ > Z + gluons "
            "performance table used by pyAmpliCol.tex."
        )
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
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
    record.add_argument("--wall-us-per-point", type=float)
    record.add_argument("--runtime-us-per-point", type=float)
    record.add_argument("--notes", default="")
    record.add_argument("--error", default="")

    args = parser.parse_args(argv)
    command = args.command or "render"
    if command == "reset":
        data: dict[str, Any] = {}
    else:
        data = _load_data(args.data)
    _refresh_metadata(data)

    if command == "record":
        _record_row(
            data,
            int(args.n),
            str(args.mode),
            status=str(args.status),
            generation_s=args.generation_s,
            wall_us_per_point=args.wall_us_per_point,
            runtime_us_per_point=args.runtime_us_per_point,
            notes=str(args.notes),
            error=str(args.error),
        )

    _write_data(args.data, data)
    table = render_table(data)
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


def _write_data(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _refresh_metadata(data: dict[str, Any]) -> None:
    data["schema_version"] = 1
    data["created_by"] = "pyAmpliCol/docs/z_performance_table.py"
    data["updated_at"] = datetime.now(timezone.utc).isoformat()
    data["multiplicity_convention"] = "final_state"
    data["process_family"] = "d d~ > z + (n-1)g"
    data["settings"] = {
        "target_runtime_s": DEFAULT_TARGET_RUNTIME_S,
        "batch_size": DEFAULT_BATCH_SIZE,
        "symbolica_output_chunk_size": DEFAULT_SYMBOLICA_OUTPUT_CHUNK_SIZE,
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
    wall_us_per_point: float | None,
    runtime_us_per_point: float | None,
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
    if wall_us_per_point is not None:
        row["wall_us_per_point"] = float(wall_us_per_point)
    if runtime_us_per_point is not None:
        row["runtime_us_per_point"] = float(runtime_us_per_point)
    if notes:
        row["notes"] = notes
    if error:
        row["error"] = error
    modes[mode_key] = row


def process_for_n(n_final: int) -> str:
    return f"d d~ > Z + ({n_final}-1)*g"


def render_table(data: dict[str, Any]) -> str:
    entries = data.get("entries", {})
    if not isinstance(entries, dict):
        entries = {}
    lines = [
        "% Generated by docs/z_performance_table.py; edit z_performance_data.json instead.",
        r"\begin{landscape}",
        (
            r"\section{\texorpdfstring{Dedicated \(d\bar d\to Z+\) Gluon "
            r"Performance}{Dedicated d dbar to Z plus Gluon Performance}}"
        ),
        r"\begingroup",
        r"\tiny",
        r"\setlength{\tabcolsep}{1.55pt}",
        r"\renewcommand{\arraystretch}{1.03}",
        (
            r"\begin{longtable}{@{}r L{1.02in} L{0.60in} L{1.00in} "
            r"r r r L{1.30in}@{}}"
        ),
        r"\toprule",
        (
            r"\textbf{n} & \textbf{process} & \textbf{route} & \textbf{setup} "
            r"& \textbf{gen [s]} & \textbf{wall [us/pt]} "
            r"& \textbf{eval [us/pt]} & \textbf{notes} \\"
        ),
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        (
            r"\textbf{n} & \textbf{process} & \textbf{route} & \textbf{setup} "
            r"& \textbf{gen [s]} & \textbf{wall [us/pt]} "
            r"& \textbf{eval [us/pt]} & \textbf{notes} \\"
        ),
        r"\midrule",
        r"\endhead",
    ]
    for n_final in range(1, 10):
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
        if n_final != 9:
            lines.append(r"\midrule[0.45pt]")
    lines.extend(
        [
            r"\bottomrule",
            r"\end{longtable}",
            r"\endgroup",
            r"\par\smallskip",
            (
                r"\noindent\footnotesize Final-state multiplicity \(n\) means "
                r"\(d\bar d\to Z+(n-1)g\).  \PAC\ rows use Rusticol, batch size "
                r"64, output chunk size 128, stage-local evaluator inputs, ten "
                r"Horner iterations, and \texttt{time-process --target-runtime 10}. "
                r"The green row is the default \PAC\ route, SymJIT O3.  Generated "
                r"processes are kept under \texttt{docs/.z\_performance\_outputs}; "
                r"the driver is \texttt{docs/run\_z\_performance\_table.py}."
            ),
            r"\par\smallskip",
            (
                r"\noindent\footnotesize To reproduce one \PAC\ row, run "
                r"generation followed by timing.  This example is the \(n=7\) "
                r"JIT O3 row; use optimization level 1 for the JIT O1 row, or "
                r"the compiled/ASM flags from \texttt{run\_z\_performance\_table.py} "
                r"for the ASM row."
            ),
            r"\begin{lstlisting}[language=bash,basicstyle=\ttfamily\tiny,breaklines=true,breakatwhitespace=true]",
            (
                "env PYTHONPATH=src dependencies/.venv/bin/python -m pyamplicol generate-process "
                "'d d~ > z g g g g g g' docs/.z_performance_outputs/manual/n7/jit_o3 "
                "--replace --n_cores 5 --color-accuracy lc --batch-size 64 "
                "--symbolica-n-cores 5 --symbolica-output-chunk-size 128 "
                "--symbolica-stage-local-parameter-layout --symbolica-iterations 10 "
                "--symbolica-max-horner-scheme-variables 1000 "
                "--symbolica-max-common-pair-cache-entries 5000000 "
                "--symbolica-max-common-pair-distance 1000 --lc-sector-ids 0 "
                "--reference-color-order 2,4,5,6,7,8,9,1,3 "
                "--symbolica-evaluator-backend jit --symbolica-jit-optimization-level 3"
            ),
            (
                "env PYTHONPATH=src dependencies/.venv/bin/python -m pyamplicol time-process "
                "--target-runtime 10 --batch-size 64 --json "
                "docs/.z_performance_outputs/manual/n7/jit_o3"
            ),
            r"\end{lstlisting}",
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
) -> list[str]:
    status = str(row.get("status", "missing"))
    if status == "timeout" and mode_key == "cpp_o3":
        return [
            r"\textcolor{speedred}{\texttt{t/o >15 min}}",
            _missing(),
            _missing(),
            _notes(row),
        ]
    if status == "ram_limit":
        return [
            r"\textcolor{speedred}{\texttt{>30 GB RAM}}",
            _missing(color="black!45") if mode_key == "amplicol" else _missing(),
            _missing(),
            _notes(row),
        ]
    if status not in {"ok", "missing"}:
        return [
            rf"\textcolor{{speedred}}{{\texttt{{{_latex_escape(status)}}}}}",
            _missing(color="black!45") if mode_key == "amplicol" else _missing(),
            _missing(),
            _notes(row),
        ]
    if status == "missing":
        if mode_key == "amplicol":
            return [_missing(), _missing(color="black!45"), _missing(), ""]
        return [_missing(), _missing(), _missing(), ""]
    generation = _optional_float(row.get("generation_s"))
    wall = _optional_float(row.get("wall_us_per_point"))
    runtime = _optional_float(row.get("runtime_us_per_point"))
    if mode_key == "amplicol":
        return [
            _format_plain(generation),
            _missing(color="black!45"),
            _format_plain(runtime),
            _notes(row),
        ]
    return [
        _format_with_ratio(generation, ref_generation),
        _format_with_ratio(wall, ref_runtime),
        _format_with_ratio(runtime, ref_runtime),
        _notes(row),
    ]


def _notes(row: dict[str, Any]) -> str:
    notes = str(row.get("notes", ""))
    error = str(row.get("error", ""))
    if notes:
        if (
            notes.startswith("generated process kept at ")
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
