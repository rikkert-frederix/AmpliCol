#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Sequence


DOCS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = DOCS_DIR.parent
DEFAULT_DATA = DOCS_DIR / "model_results_data.json"
DEFAULT_SCALARS_TABLE = DOCS_DIR / "scalars_model_table.tex"
DEFAULT_SCALARS_FULL_TABLE = DOCS_DIR / "scalars_full_model_table.tex"
DEFAULT_GRAVITY_TABLE = DOCS_DIR / "scalar_gravity_model_table.tex"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / ".ufo_support_outputs" / "model_matrices"
DEFAULT_TARGET_RUNTIME_S = 10.0
DEFAULT_MEMORY_LIMIT_GB = 30.0
DEFAULT_GENERATION_TIMEOUT_S = 300
DEFAULT_TIMING_TIMEOUT_S = 600
DEFAULT_BATCH_SIZE = 128


@dataclass(frozen=True)
class ModelProfile:
    key: str
    label: str
    source_name: str
    minimum_n: int
    maximum_n: int

    def process(self, n_final: int) -> str:
        if self.source_name == "scalars":
            final_state = " ".join("scalar_0" for _ in range(n_final))
            return f"scalar_0 scalar_0 > {final_state}"
        final_state = " ".join("graviton" for _ in range(n_final))
        return f"scalar_0 scalar_0 > {final_state}"

    def expected_value(self, n_final: int) -> float | None:
        if self.key == "scalars":
            return 1.0 / math.factorial(n_final)
        return None


MODEL_PROFILES = {
    "scalars": ModelProfile(
        key="scalars",
        label="massless scalar contact model",
        source_name="scalars",
        minimum_n=2,
        maximum_n=8,
    ),
    "scalars-full": ModelProfile(
        key="scalars_full",
        label="massless scalar all-tree model",
        source_name="scalars",
        minimum_n=2,
        maximum_n=7,
    ),
    "scalar-gravity": ModelProfile(
        key="scalar_gravity",
        label="scalar-gravity model",
        source_name="scalar_gravity",
        minimum_n=2,
        maximum_n=4,
    ),
}


class CommandFailure(RuntimeError):
    def __init__(self, command: Sequence[str], completed: subprocess.CompletedProcess[str]):
        self.command = tuple(command)
        self.returncode = int(completed.returncode)
        self.stdout = completed.stdout
        self.stderr = completed.stderr
        super().__init__(
            f"command exited with {self.returncode}: {' '.join(self.command)}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Generate and render the contact-scalar, full-scalar, and "
            "scalar-gravity UFO validation matrices."
        ),
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--scalars-table", type=Path, default=DEFAULT_SCALARS_TABLE)
    parser.add_argument(
        "--scalars-full-table",
        type=Path,
        default=DEFAULT_SCALARS_FULL_TABLE,
    )
    parser.add_argument(
        "--scalar-gravity-table",
        type=Path,
        default=DEFAULT_GRAVITY_TABLE,
    )
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--no-recompile", action="store_true")
    subparsers = parser.add_subparsers(dest="command")

    subparsers.add_parser("render", help="Render all tables from the JSON cache.")

    run = subparsers.add_parser(
        "run",
        help="Generate, validate, time, record, and render one model result.",
    )
    run.add_argument("--model", choices=tuple(MODEL_PROFILES), required=True)
    run.add_argument("--n", type=int, required=True)
    run.add_argument(
        "--target-runtime",
        type=float,
        default=DEFAULT_TARGET_RUNTIME_S,
    )
    run.add_argument(
        "--generation-timeout",
        type=int,
        default=DEFAULT_GENERATION_TIMEOUT_S,
        help=(
            "Per-process generation timeout in seconds (default: 300). "
            "All generation remains subject to the 30 GB memory watchdog."
        ),
    )
    run.add_argument(
        "--timing-timeout",
        type=int,
        default=DEFAULT_TIMING_TIMEOUT_S,
        help=(
            "Timeout in seconds for timing and 50-digit validation of an "
            "already generated artifact (default: 600)."
        ),
    )
    run.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Retain an existing process directory and only rerun validation/timing.",
    )
    run.add_argument(
        "--generation-s",
        type=float,
        help="Required with --reuse-existing unless already present in the cache.",
    )
    run.add_argument(
        "--jit-compile-s",
        type=float,
        help="JIT time associated with a reused process directory.",
    )

    sources = subparsers.add_parser(
        "measure-sources",
        help="Measure UFO/JSON conversion and exact compiled-model loading separately.",
    )
    sources.add_argument("--model", choices=tuple(MODEL_PROFILES), required=True)
    sources.add_argument(
        "--source",
        choices=("all", "ufo", "json", "compiled"),
        default="all",
    )

    args = parser.parse_args(argv)
    data = _load_data(args.data)
    _refresh_metadata(data)
    command = args.command or "render"

    failure: Exception | None = None
    try:
        if command == "run":
            profile = MODEL_PROFILES[str(args.model)]
            _validate_n(profile, int(args.n))
            _run_case(
                data,
                profile,
                int(args.n),
                output_root=args.output_root,
                target_runtime_s=float(args.target_runtime),
                generation_timeout_s=int(args.generation_timeout),
                timing_timeout_s=int(args.timing_timeout),
                reuse_existing=bool(args.reuse_existing),
                generation_s=args.generation_s,
                jit_compile_s=args.jit_compile_s,
            )
        elif command == "measure-sources":
            _measure_sources(
                data,
                MODEL_PROFILES[str(args.model)],
                output_root=args.output_root,
                source_choice=str(args.source),
            )
    except Exception as error:
        failure = error
    finally:
        _write_results(
            data,
            data_path=args.data,
            scalars_table=args.scalars_table,
            scalars_full_table=args.scalars_full_table,
            gravity_table=args.scalar_gravity_table,
            refresh_pdf=not args.no_recompile,
        )
    if failure is not None:
        raise failure
    return 0


def _validate_n(profile: ModelProfile, n_final: int) -> None:
    if not profile.minimum_n <= n_final <= profile.maximum_n:
        raise ValueError(
            f"{profile.key} n must be in "
            f"[{profile.minimum_n}, {profile.maximum_n}], got {n_final}"
        )


def _load_data(path: Path) -> dict[str, Any]:
    if not path.exists() or not path.read_text(encoding="utf-8").strip():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise TypeError(f"{path} must contain a JSON object")
    return payload


def _refresh_metadata(data: dict[str, Any]) -> None:
    data["schema_version"] = 1
    data["created_by"] = "pyAmpliCol/docs/model_results.py"
    data["updated_at"] = datetime.now(timezone.utc).isoformat()
    data["generation_timeout_s"] = DEFAULT_GENERATION_TIMEOUT_S
    data["memory_limit_gb"] = DEFAULT_MEMORY_LIMIT_GB
    data["target_runtime_s"] = DEFAULT_TARGET_RUNTIME_S
    data.setdefault("models", {})


def _model_data(data: dict[str, Any], profile: ModelProfile) -> dict[str, Any]:
    models = data.setdefault("models", {})
    if not isinstance(models, dict):
        raise TypeError("model-results cache field 'models' must be an object")
    model = models.setdefault(
        profile.key,
        {
            "label": profile.label,
            "entries": {},
            "source_measurements": {},
        },
    )
    if not isinstance(model, dict):
        raise TypeError(f"model-results cache for {profile.key} must be an object")
    model["label"] = profile.label
    model.setdefault("entries", {})
    model.setdefault("source_measurements", {})
    return model


def _run_case(
    data: dict[str, Any],
    profile: ModelProfile,
    n_final: int,
    *,
    output_root: Path,
    target_runtime_s: float,
    generation_timeout_s: int,
    timing_timeout_s: int,
    reuse_existing: bool,
    generation_s: float | None,
    jit_compile_s: float | None,
) -> None:
    model_data = _model_data(data, profile)
    entries = model_data["entries"]
    assert isinstance(entries, dict)
    existing_entry = entries.get(str(n_final), {})
    if not isinstance(existing_entry, dict):
        existing_entry = {}
    output_dir = output_root / profile.key / f"n{n_final}_lc"
    model_artifact = (
        output_root / "models" / f"{profile.source_name}.pyAmplicol-model.json"
    )
    process = profile.process(n_final)
    entry: dict[str, Any] = {
        "status": "running",
        "process": process,
        "n_final": n_final,
        "color_accuracy": "lc",
        "output_dir": str(output_dir.relative_to(PROJECT_ROOT)),
        "model_artifact": str(model_artifact.relative_to(PROJECT_ROOT)),
        "target_runtime_s": target_runtime_s,
        "generation_timeout_s": generation_timeout_s,
        "timing_timeout_s": timing_timeout_s,
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }
    entries[str(n_final)] = entry
    active_phase = "generation"

    try:
        if reuse_existing:
            if not (output_dir / "process_manifest.json").is_file():
                raise FileNotFoundError(
                    f"cannot reuse missing process artifact {output_dir}"
                )
            recorded_generation = (
                generation_s
                if generation_s is not None
                else existing_entry.get("generation_s")
            )
            if recorded_generation is None:
                raise ValueError(
                    "--reuse-existing requires --generation-s when the cache "
                    "does not already contain it"
                )
            entry["generation_s"] = float(recorded_generation)
            recorded_jit = (
                jit_compile_s
                if jit_compile_s is not None
                else existing_entry.get("jit_compile_s")
            )
            if recorded_jit is not None:
                entry["jit_compile_s"] = float(recorded_jit)
            entry["generation_measurement"] = "reused preserved artifact"
        else:
            _archive_existing_output(output_dir)
            generation = _run_json(
                _watched_command(
                    _generation_command(
                        profile,
                        n_final,
                        process=process,
                        output_dir=output_dir,
                        model_artifact=model_artifact,
                    ),
                    timeout_s=generation_timeout_s,
                )
            )
            entry["generation_s"] = float(generation["generation_s"])
            entry["jit_compile_s"] = float(generation.get("jit_compile_s", 0.0))
            entry["generation_measurement"] = "fresh watched generation"

        manifest = json.loads(
            (output_dir / "process_manifest.json").read_text(encoding="utf-8")
        )
        dag = manifest["dag_summary"]
        entry["dag"] = {
            "currents": int(dag["current_count"]),
            "interactions": int(dag["interaction_count"]),
            "roots": int(dag["amplitude_root_count"]),
        }
        active_phase = "timing"
        timing = _run_json(
            _watched_command(
                _timing_command(output_dir, target_runtime_s=target_runtime_s),
                timeout_s=timing_timeout_s,
            )
        )
        timing_profile = timing["profile"]
        entry["runtime"] = {
            "wall_us_per_point": float(timing_profile["wall_us_per_point"]),
            "wall_error_us_per_point": float(
                timing_profile["wall_us_per_point_error"]
            ),
            "pure_evaluator_us_per_point": float(
                timing_profile["pure_evaluator_us_per_point"]
            ),
            "pure_evaluator_error_us_per_point": float(
                timing_profile["pure_evaluator_us_per_point_error"]
            ),
            "samples": int(timing_profile["samples"]),
            "batch_size": int(timing_profile["block_size"]),
        }
        double_value = float(timing["values"][0])
        active_phase = "50-digit validation"
        high_precision_value = _high_precision_value(
            output_dir,
            timeout_s=timing_timeout_s,
        )
        denominator = max(abs(double_value), abs(high_precision_value), 1.0e-300)
        relative_difference = abs(double_value - high_precision_value) / denominator
        entry["value"] = {
            "double": double_value,
            "precision_50": high_precision_value,
            "relative_difference": relative_difference,
        }
        expected = profile.expected_value(n_final)
        if expected is not None:
            expected_relative = abs(double_value - expected) / max(
                abs(double_value), abs(expected), 1.0e-300
            )
            entry["value"]["expected"] = expected
            entry["value"]["expected_relative_difference"] = expected_relative
            if expected_relative > 1.0e-12:
                raise ValueError(
                    f"{process} returned {double_value}, expected {expected} "
                    f"(relative difference {expected_relative:.3g})"
                )
        if relative_difference > 1.0e-12:
            raise ValueError(
                f"{process} double/50-digit relative difference is "
                f"{relative_difference:.3g}"
            )
        entry["status"] = "ok"
    except Exception as exc:
        entry["status"] = _failure_status(exc)
        entry["failure_phase"] = active_phase
        entry["error"] = _failure_message(exc)
        raise
    finally:
        entry["updated_at"] = datetime.now(timezone.utc).isoformat()


def _generation_command(
    profile: ModelProfile,
    n_final: int,
    *,
    process: str,
    output_dir: Path,
    model_artifact: Path,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "pyamplicol",
        "generate-process",
        process,
        str(output_dir),
        "--model",
        str(model_artifact),
        "--color-accuracy",
        "lc",
        "--symbolica-evaluator-backend",
        "jit",
        "--symbolica-jit-optimization-level",
        "1",
        "--symbolica-n-cores",
        "5",
        "--n-cores",
        "1",
        "--batch-size",
        str(DEFAULT_BATCH_SIZE),
        "--symbolica-output-chunk-size",
        "128",
        "--monitor",
        "--json",
    ]
    if profile.key == "scalars":
        command.extend(("--max-coupling-order", "QCD=1"))
    return command


def _timing_command(output_dir: Path, *, target_runtime_s: float) -> list[str]:
    return [
        sys.executable,
        "-m",
        "pyamplicol",
        "time-process",
        str(output_dir),
        "--precision",
        "16",
        "--target-runtime",
        str(target_runtime_s),
        "--batch-size",
        str(DEFAULT_BATCH_SIZE),
        "--json",
    ]


def _high_precision_value(output_dir: Path, *, timeout_s: int) -> float:
    command = _watched_command(
        [
            sys.executable,
            str(output_dir / "check_standalone.py"),
            "--precision",
            "50",
        ],
        timeout_s=timeout_s,
    )
    completed = _run_command(command)
    match = re.search(r"^values:\s*\[([^,\]]+)", completed.stdout, re.MULTILINE)
    if match is None:
        raise ValueError("could not parse 50-digit value from check_standalone.py")
    return float(match.group(1))


def _measure_sources(
    data: dict[str, Any],
    profile: ModelProfile,
    *,
    output_root: Path,
    source_choice: str,
) -> None:
    model_data = _model_data(data, profile)
    measurements = model_data["source_measurements"]
    assert isinstance(measurements, dict)
    source_root = output_root / "models" / "source_measurements"
    source_root.mkdir(parents=True, exist_ok=True)
    source_kinds = (
        ("ufo", "json")
        if source_choice == "all"
        else tuple(
            source
            for source in ("ufo", "json")
            if source == source_choice
        )
    )
    for source_kind in source_kinds:
        if source_kind == "ufo":
            source = PROJECT_ROOT / "assets" / "models" / "ufo" / profile.source_name
        else:
            source = (
                PROJECT_ROOT
                / "assets"
                / "models"
                / "json"
                / profile.source_name
                / f"{profile.source_name}.json"
            )
        output = source_root / (
            f"{profile.source_name}-{source_kind}.pyAmplicol-model.json"
        )
        _archive_existing_file(output)
        payload = _run_json(
            _watched_command(
                [
                    sys.executable,
                    "-m",
                    "pyamplicol",
                    "compile-model",
                    "--model",
                    str(source),
                    "--output",
                    str(output),
                    "--no-model-cache",
                    "--json",
                ],
                timeout_s=DEFAULT_GENERATION_TIMEOUT_S,
            )
        )
        measurements[source_kind] = {
            "status": "ok",
            "seconds": float(payload["conversion_seconds"]),
            "output": str(output.relative_to(PROJECT_ROOT)),
            "digest": payload["source"]["digest"],
        }

    if source_choice in {"all", "compiled"}:
        compiled_model = (
            output_root / "models" / f"{profile.source_name}.pyAmplicol-model.json"
        )
        code = (
            "import json,time;"
            "from pyamplicol.model_source import compile_model_source;"
            f"path={str(compiled_model)!r};"
            "started=time.perf_counter();"
            "compile_model_source(path);"
            "print(json.dumps({'seconds':time.perf_counter()-started}))"
        )
        payload = _run_json(
            _watched_command(
                [sys.executable, "-c", code],
                timeout_s=60,
            )
        )
        measurements["compiled"] = {
            "status": "ok",
            "seconds": float(payload["seconds"]),
            "output": str(compiled_model.relative_to(PROJECT_ROOT)),
            "measurement": "compile_model_source after Python imports",
        }


def _watched_command(command: Sequence[str], *, timeout_s: int) -> list[str]:
    return [
        sys.executable,
        str(PROJECT_ROOT / "scripts" / "run_with_memory_watch.py"),
        "--limit-gb",
        str(DEFAULT_MEMORY_LIMIT_GB),
        "--",
        "gtimeout",
        str(timeout_s),
        *command,
    ]


def _run_json(command: Sequence[str]) -> dict[str, Any]:
    completed = _run_command(command)
    payload = json.loads(completed.stdout)
    if not isinstance(payload, dict):
        raise TypeError("command JSON output must be an object")
    if payload.get("available") is False:
        raise RuntimeError(str(payload.get("error", "command reported unavailable")))
    return payload


def _run_command(command: Sequence[str]) -> subprocess.CompletedProcess[str]:
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    completed = subprocess.run(
        list(command),
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=None,
        check=False,
    )
    if completed.returncode != 0:
        raise CommandFailure(command, completed)
    return completed


def _archive_existing_output(path: Path) -> None:
    if not path.exists():
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    destination = path.with_name(f"{path.name}_before_{stamp}")
    counter = 1
    while destination.exists():
        destination = path.with_name(f"{path.name}_before_{stamp}_{counter}")
        counter += 1
    path.rename(destination)


def _archive_existing_file(path: Path) -> None:
    if not path.exists():
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive = path.parent / "archive"
    archive.mkdir(parents=True, exist_ok=True)
    destination = archive / f"{path.stem}_before_{stamp}{path.suffix}"
    path.rename(destination)


def _failure_status(error: Exception) -> str:
    if isinstance(error, CommandFailure):
        if error.returncode == 124:
            return "timeout"
        if error.returncode in {137, -9}:
            return "ram_limit_or_killed"
    return "error"


def _failure_message(error: Exception) -> str:
    if not isinstance(error, CommandFailure):
        return str(error)
    if "generate-process" in error.command:
        phase = "generate-process"
    elif "time-process" in error.command:
        phase = "time-process"
    elif "check_standalone.py" in error.command:
        phase = "50-digit validation"
    else:
        phase = "subprocess"
    if error.returncode == 124:
        return f"{phase} exceeded its configured timeout"
    if error.returncode in {137, -9}:
        return f"{phase} was killed at the configured memory limit"
    return f"{phase} exited with status {error.returncode}"


def _write_results(
    data: dict[str, Any],
    *,
    data_path: Path,
    scalars_table: Path,
    scalars_full_table: Path,
    gravity_table: Path,
    refresh_pdf: bool,
) -> None:
    _atomic_write_json(data_path, data)
    scalars_table.write_text(
        render_model_table(data, MODEL_PROFILES["scalars"]),
        encoding="utf-8",
    )
    scalars_full_table.write_text(
        render_model_table(data, MODEL_PROFILES["scalars-full"]),
        encoding="utf-8",
    )
    gravity_table.write_text(
        render_model_table(data, MODEL_PROFILES["scalar-gravity"]),
        encoding="utf-8",
    )
    print(f"wrote {data_path}")
    print(f"wrote {scalars_table}")
    print(f"wrote {scalars_full_table}")
    print(f"wrote {gravity_table}")
    if refresh_pdf:
        _refresh_pdf()


def _atomic_write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def render_model_table(data: dict[str, Any], profile: ModelProfile) -> str:
    models = data.get("models", {})
    model = models.get(profile.key, {}) if isinstance(models, dict) else {}
    entries = model.get("entries", {}) if isinstance(model, dict) else {}
    measurements = (
        model.get("source_measurements", {}) if isinstance(model, dict) else {}
    )
    if profile.key == "scalars_full" and not measurements:
        contact_model = models.get("scalars", {}) if isinstance(models, dict) else {}
        measurements = (
            contact_model.get("source_measurements", {})
            if isinstance(contact_model, dict)
            else {}
        )
    lines = [
        "% Generated by docs/model_results.py; edit model_results_data.json instead.",
        rf"\subsection*{{{_tex(profile.label.title())} Validation Ladder}}",
        r"\begingroup",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\renewcommand{\arraystretch}{1.12}",
        r"\begin{longtable}{@{}rL{2.35in}rrrrL{0.78in}rrr@{}}",
        r"\toprule",
        r"$n$ & process & gen [s] & JIT [s] & currents & interactions "
        r"& value & eval [$\mu$s] & wall [$\mu$s] & rel. 50d \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"$n$ & process & gen [s] & JIT [s] & currents & interactions "
        r"& value & eval [$\mu$s] & wall [$\mu$s] & rel. 50d \\",
        r"\midrule",
        r"\endhead",
    ]
    for n_final in range(profile.minimum_n, profile.maximum_n + 1):
        entry = entries.get(str(n_final), {}) if isinstance(entries, dict) else {}
        lines.append(_render_entry(profile, n_final, entry))
    lines.extend(
        [
            r"\bottomrule",
            r"\end{longtable}",
            _render_source_measurements(measurements),
            _render_entry_notes(entries),
            r"\par\smallskip",
            _render_model_method(profile),
            (
                r" Values use the retained RAMBO point generated with seed 101 at "
                r"$\sqrt{s}=1\,\mathrm{TeV}$. "
                r"Generation starts from the compiled model and uses SymJIT O1, "
                r"batch/chunk size 128 and a 30 GB watchdog."
            ),
            _render_generation_limit(profile),
            (
                r" The eval column isolates generated evaluators; wall includes "
                r"Rusticol orchestration. Runtime uses double precision and "
                r"\texttt{--target-runtime 10}; the final column compares the "
                r"same point with a 50-digit reevaluation. Superseded process "
                r"directories are archived beside the active output."
            ),
            r"\endgroup",
            "",
        ]
    )
    return "\n".join(lines)


def _render_entry(
    profile: ModelProfile,
    n_final: int,
    entry: object,
) -> str:
    if not isinstance(entry, dict) or not entry:
        return rf"{n_final} & ${_tex_process(profile.process(n_final))}$ & \multicolumn{{8}}{{c}}{{\textcolor{{black!45}}{{not run}}}} \\"
    status = str(entry.get("status", "not_run"))
    if status != "ok":
        label = {
            "timeout": _timeout_label(entry),
            "ram_limit_or_killed": r"$>30$ GB RAM",
        }.get(status, _tex(status))
        return rf"{n_final} & ${_tex_process(profile.process(n_final))}$ & \multicolumn{{8}}{{c}}{{\textcolor{{speedred}}{{{label}}}}} \\"
    dag = entry["dag"]
    runtime = entry["runtime"]
    value = entry["value"]
    return " & ".join(
        (
            str(n_final),
            f"${_tex_process(str(entry['process']))}$",
            _format_number(float(entry["generation_s"])),
            _format_number(float(entry.get("jit_compile_s", 0.0))),
            f"{int(dag['currents']):,}",
            f"{int(dag['interactions']):,}",
            _format_value(float(value["double"])),
            _format_number(float(runtime["pure_evaluator_us_per_point"])),
            _format_number(float(runtime["wall_us_per_point"])),
            _format_scientific(float(value["relative_difference"])),
        )
    ) + r" \\"


def _timeout_label(entry: dict[str, Any]) -> str:
    phase = str(entry.get("failure_phase", "generation"))
    if phase == "generation":
        timeout_s = int(
            entry.get("generation_timeout_s", DEFAULT_GENERATION_TIMEOUT_S)
        )
        prefix = ""
    else:
        timeout_s = int(entry.get("timing_timeout_s", DEFAULT_TIMING_TIMEOUT_S))
        prefix = "run " if phase == "timing" else "50d "
    if timeout_s % 60 == 0:
        return rf"{prefix}t/o $>{timeout_s // 60}$ min"
    return rf"{prefix}t/o $>{timeout_s}$ s"


def _render_model_method(profile: ModelProfile) -> str:
    if profile.key == "scalars":
        return (
            r"\noindent\textit{Contact oracle.} The scalar ladder explicitly "
            r"sets \texttt{--max-coupling-order QCD=1}. With $\lambda=1$ this "
            r"retains only the direct $(n+2)$-point vertex, so no propagator "
            r"denominator enters and the identical-final-state normalization "
            r"gives $1/n!$. It tests generic high-rank contact decomposition, "
            r"not the sum of all scalar-model tree diagrams; the uncapped "
            r"all-tree ladder below provides that complementary test."
        )
    if profile.key == "scalars_full":
        return (
            r"\noindent\textit{Complete scalar trees.} These entries omit the "
            r"coupling-order cap and therefore sum the direct contact term with "
            r"every allowed multi-vertex tree and internal scalar propagator. "
            r"They use the same model parameters and phase-space convention as "
            r"the contact oracle above."
        )
    return (
        r"\noindent\textit{Full-model ladder.} The scalar-gravity entries retain "
        r"the complete tree-level model, including internal scalar and graviton "
        r"propagators."
    )


def _render_generation_limit(profile: ModelProfile) -> str:
    if profile.key == "scalar_gravity":
        return (
            r" The ordinary generation cap is five minutes; the four-graviton "
            r"entry is allowed one hour."
        )
    if profile.key == "scalars_full":
        return (
            r" The ordinary generation cap is five minutes; the seven-scalar "
            r"entry is allowed one hour."
        )
    return r" The scalar rows use the ordinary five-minute generation cap."


def _render_source_measurements(measurements: object) -> str:
    if not isinstance(measurements, dict) or not measurements:
        return (
            r"\noindent Source conversion/loading timings have not yet been recorded."
        )
    fields = []
    for key, label in (("ufo", "UFO conversion"), ("json", "JSON conversion"), ("compiled", "compiled load")):
        payload = measurements.get(key)
        if isinstance(payload, dict) and payload.get("status") == "ok":
            fields.append(
                rf"{label}: \texttt{{{_format_number(float(payload['seconds']))} s}}"
            )
    if not fields:
        return r"\noindent Source conversion/loading timings are unavailable."
    return r"\noindent " + r"; \quad ".join(fields) + "."


def _render_entry_notes(entries: object) -> str:
    if not isinstance(entries, dict):
        return ""
    notes = [
        str(entry["notes"])
        for entry in entries.values()
        if isinstance(entry, dict) and entry.get("notes")
    ]
    if not notes:
        return ""
    return r"\par\smallskip\noindent\textit{Limit note.} " + " ".join(
        _tex(note) for note in notes
    )


def _format_number(value: float) -> str:
    if value == 0.0:
        return "0"
    magnitude = abs(value)
    if magnitude < 1.0e-3 or magnitude >= 1.0e4:
        return f"{value:.3g}"
    return f"{value:.4g}"


def _format_value(value: float) -> str:
    return rf"\texttt{{{value:.8g}}}"


def _format_scientific(value: float) -> str:
    if value == 0.0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10.0**exponent)
    return rf"${mantissa:.2f}\times10^{{{exponent}}}$"


def _tex_process(process: str) -> str:
    initial, final = process.split(">", 1)
    replacements = {
        "scalar_0": r"\phi_0",
        "graviton": r"G",
    }
    left = r"\,".join(
        replacements.get(token, _tex(token)) for token in initial.split()
    )
    right = r"\,".join(
        replacements.get(token, _tex(token)) for token in final.split()
    )
    return rf"{left}\to {right}"


def _tex(value: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
    }
    return "".join(replacements.get(character, character) for character in value)


def _refresh_pdf() -> None:
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    environment["LANG"] = "C"
    completed = subprocess.run(
        [
            "latexmk",
            "-pdf",
            "-interaction=nonstopmode",
            "-halt-on-error",
            "pyAmpliCol.tex",
        ],
        cwd=DOCS_DIR,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "pyAmpliCol.tex recompilation failed\n"
            f"stdout:\n{completed.stdout[-4000:]}\n"
            f"stderr:\n{completed.stderr[-4000:]}"
        )
    print(f"wrote {DOCS_DIR / 'pyAmpliCol.pdf'}")


if __name__ == "__main__":
    raise SystemExit(main())
