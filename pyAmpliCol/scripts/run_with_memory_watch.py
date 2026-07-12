#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence

from memory_watch_support import ProcessTreeMemoryMonitor, process_tree_pids


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a command while limiting total process-tree memory. On macOS "
            "this includes compressed physical footprint as well as RSS."
        )
    )
    parser.add_argument("--limit-gb", type=float, required=True)
    parser.add_argument("--poll-s", type=float, default=1.0)
    parser.add_argument(
        "--physical-footprint-poll-s",
        type=float,
        default=5.0,
        help="macOS physical-footprint polling interval (default: 5 seconds).",
    )
    parser.add_argument("--stop-file", type=Path, default=Path("stop.order"))
    parser.add_argument(
        "--report-json",
        type=Path,
        default=None,
        help="Optional path where peak memory and exit status are written as JSON.",
    )
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if not args.command or args.command[0] != "--":
        parser.error("command must be passed after --")
    args.command = args.command[1:]
    if not args.command:
        parser.error("missing command after --")
    if args.poll_s <= 0:
        parser.error("--poll-s must be positive")
    if args.physical_footprint_poll_s <= 0:
        parser.error("--physical-footprint-poll-s must be positive")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    limit_bytes = int(args.limit_gb * 1024**3)
    child = subprocess.Popen(
        args.command,
        start_new_session=True,
        preexec_fn=lambda: _apply_process_memory_limit(limit_bytes),
    )
    monitor = ProcessTreeMemoryMonitor(
        physical_footprint_poll_s=args.physical_footprint_poll_s,
        physical_footprint_high_watermark_bytes=int(0.8 * limit_bytes),
    )
    warned_no_memory = False
    try:
        while True:
            returncode = child.poll()
            if returncode is not None:
                _write_report(
                    args.report_json,
                    limit_gb=args.limit_gb,
                    limit_bytes=limit_bytes,
                    monitor=monitor,
                    returncode=returncode,
                )
                return returncode

            if args.stop_file.exists():
                tracked_pids = _terminate_tree(child.pid, signal.SIGTERM)
                returncode = _wait_or_kill(child, tracked_pids)
                _write_report(
                    args.report_json,
                    limit_gb=args.limit_gb,
                    limit_bytes=limit_bytes,
                    monitor=monitor,
                    returncode=returncode,
                )
                return returncode

            sample = monitor.sample(child.pid)
            memory_bytes = sample.effective_bytes
            if memory_bytes is None:
                if not warned_no_memory:
                    print(
                        "memory watchdog: process-tree memory polling is "
                        "unavailable; relying on per-process OS memory limits",
                        file=sys.stderr,
                        flush=True,
                    )
                    warned_no_memory = True
            elif memory_bytes > limit_bytes:
                metric = (
                    "physical footprint"
                    if sample.effective_metric == "physical_footprint"
                    else "RSS"
                )
                print(
                    f"memory watchdog: {metric} "
                    f"{memory_bytes / 1024**3:.3f} GiB exceeded "
                    f"limit {args.limit_gb:.3f} GiB",
                    file=sys.stderr,
                    flush=True,
                )
                tracked_pids = _terminate_tree(child.pid, signal.SIGTERM)
                _wait_or_kill(child, tracked_pids)
                _write_report(
                    args.report_json,
                    limit_gb=args.limit_gb,
                    limit_bytes=limit_bytes,
                    monitor=monitor,
                    returncode=137,
                )
                return 137

            time.sleep(args.poll_s)
    except KeyboardInterrupt:
        tracked_pids = _terminate_tree(child.pid, signal.SIGINT)
        returncode = _wait_or_kill(child, tracked_pids)
        _write_report(
            args.report_json,
            limit_gb=args.limit_gb,
            limit_bytes=limit_bytes,
            monitor=monitor,
            returncode=returncode,
        )
        return returncode


def _write_report(
    path: Path | None,
    *,
    limit_gb: float,
    limit_bytes: int,
    monitor: ProcessTreeMemoryMonitor,
    returncode: int,
) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "limit_gb": limit_gb,
        "limit_bytes": limit_bytes,
        "peak_rss_bytes": monitor.peak_rss_bytes,
        "peak_rss_gb": (
            None
            if monitor.peak_rss_bytes is None
            else monitor.peak_rss_bytes / 1024**3
        ),
        "peak_physical_footprint_bytes": monitor.peak_physical_footprint_bytes,
        "peak_physical_footprint_gb": (
            None
            if monitor.peak_physical_footprint_bytes is None
            else monitor.peak_physical_footprint_bytes / 1024**3
        ),
        "peak_memory_bytes": monitor.peak_memory_bytes,
        "peak_memory_gb": (
            None
            if monitor.peak_memory_bytes is None
            else monitor.peak_memory_bytes / 1024**3
        ),
        "peak_memory_metric": monitor.peak_memory_metric,
        "returncode": returncode,
        "rss_polling_available": monitor.rss_polling_available,
        "physical_footprint_supported": monitor.physical_footprint_supported,
        "physical_footprint_polling_available": (
            monitor.physical_footprint_polling_available
        ),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _wait_or_kill(
    child: subprocess.Popen[bytes],
    tracked_pids: Sequence[int],
) -> int:
    try:
        returncode = child.wait(timeout=10.0)
    except subprocess.TimeoutExpired:
        _signal_pids(tracked_pids, signal.SIGKILL, exclude=child.pid)
        _terminate_group(child.pid, signal.SIGKILL)
        returncode = child.wait()
    else:
        # A child may have called setsid(), so waiting for the root does not
        # prove every process sampled by the watchdog has exited.
        _signal_pids(tracked_pids, signal.SIGKILL, exclude=child.pid)
    return returncode


def _apply_process_memory_limit(limit_bytes: int) -> None:
    for resource_name in ("RLIMIT_AS", "RLIMIT_DATA", "RLIMIT_RSS"):
        resource_id = getattr(resource, resource_name, None)
        if resource_id is None:
            continue
        try:
            soft, hard = resource.getrlimit(resource_id)
            new_soft = _bounded_limit(limit_bytes, hard)
            resource.setrlimit(resource_id, (new_soft, hard))
        except (OSError, ValueError):
            continue


def _bounded_limit(limit_bytes: int, hard: int) -> int:
    if hard == resource.RLIM_INFINITY:
        return limit_bytes
    return min(limit_bytes, hard)


def _terminate_group(pid: int, sig: signal.Signals) -> None:
    try:
        os.killpg(pid, sig)
    except ProcessLookupError:
        return


def _terminate_tree(pid: int, sig: signal.Signals) -> tuple[int, ...]:
    tracked_pids = process_tree_pids(pid)
    _signal_pids(tracked_pids, sig, exclude=pid)
    _terminate_group(pid, sig)
    return tracked_pids


def _signal_pids(
    pids: Sequence[int],
    sig: signal.Signals,
    *,
    exclude: int | None = None,
) -> None:
    for process_pid in reversed(tuple(pids)):
        if process_pid == exclude:
            continue
        try:
            os.kill(process_pid, sig)
        except ProcessLookupError:
            continue


if __name__ == "__main__":
    raise SystemExit(main())
