#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


DOCS_DIR = Path(__file__).resolve().parent
DEFAULT_DATA = DOCS_DIR / "result_matrix_data.json"
DEFAULT_REPORT = DOCS_DIR / "result_matrix_audit.md"


@dataclass(frozen=True)
class Finding:
    severity: str
    process_id: int
    process_key: str
    n_final: int
    mode: str
    message: str

    def markdown_row(self) -> str:
        return (
            f"| `{self.severity}` | {self.process_id} | `{self.process_key}` | "
            f"{self.n_final} | `{self.mode}` | {self.message} |"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit pyAmpliCol result-matrix timings for stale or odd cells."
    )
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--high-ratio", type=float, default=10.0)
    parser.add_argument("--nonmonotonic-factor", type=float, default=3.0)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()

    data = json.loads(args.data.read_text(encoding="utf-8"))
    findings = audit_matrix(
        data,
        high_ratio=float(args.high_ratio),
        nonmonotonic_factor=float(args.nonmonotonic_factor),
    )
    report = render_report(data, findings)
    if args.write:
        args.report.write_text(report, encoding="utf-8")
        print(f"wrote {args.report}")
    else:
        print(report)
    return 0


def audit_matrix(
    data: dict[str, Any],
    *,
    high_ratio: float,
    nonmonotonic_factor: float,
) -> list[Finding]:
    findings: list[Finding] = []
    base_processes = data.get("base_processes", [])
    entries = data.get("entries", {})
    if not isinstance(base_processes, list) or not isinstance(entries, dict):
        return [
            Finding("error", 0, "matrix", 0, "schema", "Malformed matrix data.")
        ]
    for base in base_processes:
        if not isinstance(base, dict):
            continue
        process_id = int(base.get("process_id", 0))
        process_key = str(base.get("key", "unknown"))
        row = entries.get(process_key, {})
        if not isinstance(row, dict):
            continue
        findings.extend(_audit_row_staleness(process_id, process_key, row))
        findings.extend(
            _audit_row_monotonicity(
                process_id,
                process_key,
                row,
                factor=nonmonotonic_factor,
            )
        )
        findings.extend(
            _audit_row_ratios(
                process_id,
                process_key,
                row,
                threshold=high_ratio,
            )
        )
        findings.extend(_audit_row_validation(process_id, process_key, row))
        findings.extend(_audit_row_statuses(process_id, process_key, row))
    return sorted(
        findings,
        key=lambda item: (
            {"error": 0, "warning": 1, "info": 2}.get(item.severity, 3),
            item.process_id,
            item.n_final,
            item.mode,
        ),
    )


def _audit_row_staleness(
    process_id: int,
    process_key: str,
    row: dict[str, Any],
) -> Iterable[Finding]:
    for n_text, case in _iter_cases(row):
        for mode_key, mode in _iter_modes(case):
            if mode.get("status") == "ok" and mode.get("matrix_settings") is None:
                yield Finding(
                    "info",
                    process_id,
                    process_key,
                    n_text,
                    mode_key,
                    "OK cell predates explicit `matrix_settings`; rerun if this "
                    "cell is used for a precise comparison.",
                )


def _audit_row_monotonicity(
    process_id: int,
    process_key: str,
    row: dict[str, Any],
    *,
    factor: float,
) -> Iterable[Finding]:
    previous: dict[str, tuple[int, float]] = {}
    for n_final, case in _iter_cases(row):
        for mode_key, mode in _iter_modes(case):
            if mode.get("status") != "ok":
                continue
            for metric_key, label in (
                ("generation_s", "generation"),
                ("runtime_us_per_point", "runtime"),
                ("wall_us_per_point", "wall runtime"),
            ):
                value = _float(mode.get(metric_key))
                if value is None:
                    continue
                previous_key = f"{mode_key}:{metric_key}"
                if previous_key in previous:
                    previous_n, previous_value = previous[previous_key]
                    if previous_value > 0.0 and value < previous_value / factor:
                        yield Finding(
                            "warning",
                            process_id,
                            process_key,
                            n_final,
                            mode_key,
                            (
                                f"{label} drops by more than x{factor:g}: "
                                f"n={previous_n} has {previous_value:.4g}, "
                                f"n={n_final} has {value:.4g}."
                            ),
                        )
                previous[previous_key] = (n_final, value)


def _audit_row_ratios(
    process_id: int,
    process_key: str,
    row: dict[str, Any],
    *,
    threshold: float,
) -> Iterable[Finding]:
    for n_final, case in _iter_cases(row):
        reference = _mode(case, "amplicol")
        ref_runtime = _float(reference.get("runtime_us_per_point"))
        if reference.get("status") != "ok" or not ref_runtime:
            continue
        for mode_key in ("pyamplicol_jit", "pyamplicol_cpp_o3"):
            mode = _mode(case, mode_key)
            if mode.get("status") != "ok":
                continue
            runtime = _float(mode.get("wall_us_per_point")) or _float(
                mode.get("runtime_us_per_point")
            )
            if runtime is None:
                continue
            ratio = runtime / ref_runtime
            if ratio >= threshold:
                yield Finding(
                    "warning",
                    process_id,
                    process_key,
                    n_final,
                    mode_key,
                    (
                        f"wall/runtime ratio to AmpliCol is x{ratio:.2g}; "
                        f"inspect current recycling and backend codegen."
                    ),
                )


def _audit_row_validation(
    process_id: int,
    process_key: str,
    row: dict[str, Any],
) -> Iterable[Finding]:
    for n_final, case in _iter_cases(row):
        validation = _mode(case, "validation")
        if not validation:
            continue
        tolerance = _float(validation.get("tolerance")) or 1.0e-8
        rel = _float(validation.get("max_relative_difference"))
        if rel is None or rel <= tolerance:
            continue
        yield Finding(
            "error",
            process_id,
            process_key,
            n_final,
            "validation",
            (
                f"relative validation difference {rel:.3g} exceeds "
                f"tolerance {tolerance:.3g}."
            ),
        )


def _audit_row_statuses(
    process_id: int,
    process_key: str,
    row: dict[str, Any],
) -> Iterable[Finding]:
    for n_final, case in _iter_cases(row):
        validation = _mode(case, "validation")
        status = str(validation.get("status", ""))
        if status and status not in {"ok", "not_available"}:
            yield Finding(
                "error",
                process_id,
                process_key,
                n_final,
                "validation",
                f"validation status is `{status}`.",
            )
        for mode_key, mode in _iter_modes(case):
            status = str(mode.get("status", ""))
            if status in {"", "ok", "unsupported", "backend_unsupported"}:
                continue
            yield Finding(
                "error",
                process_id,
                process_key,
                n_final,
                mode_key,
                f"unexpected status `{status}`.",
            )


def render_report(data: dict[str, Any], findings: list[Finding]) -> str:
    gate = _gate_summary(data, findings)
    lines = [
        "# Result Matrix Audit",
        "",
        "This report is generated by `pyAmpliCol/docs/audit_result_matrix.py`.",
        "It flags stale metadata, nonmonotonic timings, large runtime ratios, and",
        "unresolved statuses in `result_matrix_data.json`.",
        "",
        f"Matrix updated at: `{data.get('updated_at', 'unknown')}`",
        "",
        f"Gate status: **{gate['status']}**.",
        "",
        (
            "Gate counts: "
            f"`errors`={gate['error_findings']}, "
            f"`validation_failures`={gate['validation_failures']}, "
            f"`missing_amplicol`={gate['missing_amplicol']}, "
            f"`missing_jit`={gate['missing_jit']}, "
            f"`amplicol_unsupported`={gate['amplicol_unsupported']}, "
            f"`jit_backend_unsupported`={gate['jit_backend_unsupported']}, "
            f"`missing_cpp_o3`={gate['missing_cpp_o3']}."
        ),
        "",
        (
            "`missing_cpp_o3` is informational: C++ O3 cells are intentionally "
            "filled only where generation was feasible within the matrix-run "
            "time budget. `amplicol_unsupported` records processes outside "
            "Fortran AmpliCol's supported quark-line range; pyAmpliCol entries "
            "for those cells are still shown as absolute timings."
        ),
        "",
    ]
    if not findings:
        lines.append("No findings.")
        return "\n".join(lines) + "\n"
    counts: dict[str, int] = {}
    for finding in findings:
        counts[finding.severity] = counts.get(finding.severity, 0) + 1
    lines.append(
        "Summary: "
        + ", ".join(f"`{key}`={counts[key]}" for key in sorted(counts))
        + "."
    )
    lines.extend(
        [
            "",
            "| Severity | ID | Process | n | Mode | Finding |",
            "| --- | ---: | --- | ---: | --- | --- |",
        ]
    )
    lines.extend(finding.markdown_row() for finding in findings)
    return "\n".join(lines) + "\n"


def _gate_summary(data: dict[str, Any], findings: list[Finding]) -> dict[str, int | str]:
    entries = data.get("entries", {})
    error_findings = sum(1 for finding in findings if finding.severity == "error")
    validation_failures = 0
    missing_amplicol = 0
    missing_jit = 0
    amplicol_unsupported = 0
    jit_backend_unsupported = 0
    missing_cpp_o3 = 0
    if isinstance(entries, dict):
        for row in entries.values():
            if not isinstance(row, dict):
                continue
            for _, case in _iter_cases(row):
                if case.get("status") == "not_applicable":
                    continue
                validation = _mode(case, "validation")
                if validation:
                    tolerance = _float(validation.get("tolerance")) or 1.0e-8
                    rel = _float(validation.get("max_relative_difference"))
                    if validation.get("status") != "ok" or (
                        rel is not None and rel > tolerance
                    ):
                        validation_failures += 1
                amplicol = _mode(case, "amplicol")
                if not amplicol:
                    missing_amplicol += 1
                elif amplicol.get("status") == "unsupported":
                    amplicol_unsupported += 1
                jit = _mode(case, "pyamplicol_jit")
                if not jit:
                    missing_jit += 1
                elif jit.get("status") == "backend_unsupported":
                    jit_backend_unsupported += 1
                cpp = _mode(case, "pyamplicol_cpp_o3")
                if not cpp or cpp.get("status") == "missing":
                    missing_cpp_o3 += 1
    hard_failures = (
        error_findings
        + validation_failures
        + missing_amplicol
        + missing_jit
    )
    documented_limitations = (
        amplicol_unsupported
        + jit_backend_unsupported
        + missing_cpp_o3
    )
    if hard_failures:
        status = "FAIL"
    elif documented_limitations:
        status = "PASS_WITH_LIMITATIONS"
    else:
        status = "PASS"
    return {
        "status": status,
        "error_findings": error_findings,
        "validation_failures": validation_failures,
        "missing_amplicol": missing_amplicol,
        "missing_jit": missing_jit,
        "amplicol_unsupported": amplicol_unsupported,
        "jit_backend_unsupported": jit_backend_unsupported,
        "missing_cpp_o3": missing_cpp_o3,
    }


def _iter_cases(row: dict[str, Any]) -> Iterable[tuple[int, dict[str, Any]]]:
    for n_text, case in sorted(row.items(), key=lambda item: int(item[0])):
        if isinstance(case, dict):
            yield int(n_text), case


def _iter_modes(case: dict[str, Any]) -> Iterable[tuple[str, dict[str, Any]]]:
    for mode_key in ("amplicol", "pyamplicol_jit", "pyamplicol_cpp_o3"):
        mode = _mode(case, mode_key)
        if mode:
            yield mode_key, mode


def _mode(case: dict[str, Any], key: str) -> dict[str, Any]:
    value = case.get(key)
    return value if isinstance(value, dict) else {}


def _float(value: object) -> float | None:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


if __name__ == "__main__":
    raise SystemExit(main())
