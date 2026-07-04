#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib
import os
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a command while limiting total process-tree RSS."
    )
    parser.add_argument("--limit-gb", type=float, required=True)
    parser.add_argument("--poll-s", type=float, default=1.0)
    parser.add_argument("--stop-file", type=Path, default=Path("stop.order"))
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if not args.command or args.command[0] != "--":
        parser.error("command must be passed after --")
    args.command = args.command[1:]
    if not args.command:
        parser.error("missing command after --")
    return args


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    limit_bytes = int(args.limit_gb * 1024**3)
    child = subprocess.Popen(
        args.command,
        start_new_session=True,
        preexec_fn=lambda: _apply_process_memory_limit(limit_bytes),
    )
    warned_no_rss = False
    try:
        while True:
            returncode = child.poll()
            if returncode is not None:
                return returncode

            if args.stop_file.exists():
                _terminate_group(child.pid, signal.SIGTERM)
                return _wait_or_kill(child)

            rss = _process_tree_rss_bytes(child.pid)
            if rss is None:
                if not warned_no_rss:
                    print(
                        "memory watchdog: RSS polling is unavailable; relying on "
                        "per-process OS memory limits",
                        file=sys.stderr,
                        flush=True,
                    )
                    warned_no_rss = True
            elif rss > limit_bytes:
                print(
                    f"memory watchdog: RSS {rss / 1024**3:.3f} GiB exceeded "
                    f"limit {args.limit_gb:.3f} GiB",
                    file=sys.stderr,
                    flush=True,
                )
                _terminate_group(child.pid, signal.SIGTERM)
                _wait_or_kill(child)
                return 137

            time.sleep(args.poll_s)
    except KeyboardInterrupt:
        _terminate_group(child.pid, signal.SIGINT)
        return _wait_or_kill(child)


def _wait_or_kill(child: subprocess.Popen[bytes]) -> int:
    try:
        return child.wait(timeout=10.0)
    except subprocess.TimeoutExpired:
        _terminate_group(child.pid, signal.SIGKILL)
        return child.wait()


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


def _process_tree_rss_bytes(pid: int) -> int | None:
    psutil_rss = _process_tree_rss_bytes_psutil(pid)
    if psutil_rss is not None:
        return psutil_rss
    if sys.platform == "darwin":
        return _process_tree_rss_bytes_ps(pid)
    return _process_tree_rss_bytes_proc(pid)


def _process_tree_rss_bytes_psutil(pid: int) -> int | None:
    try:
        psutil = importlib.import_module("psutil")
    except ImportError:
        return None
    try:
        root = psutil.Process(pid)
        processes = [root, *root.children(recursive=True)]
        return sum(process.memory_info().rss for process in processes)
    except Exception:
        return None


def _process_tree_rss_bytes_ps(pid: int) -> int | None:
    try:
        output = subprocess.run(
            ["ps", "-axo", "pid=,ppid=,rss="],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    children: dict[int, list[int]] = {}
    rss_kb: dict[int, int] = {}
    for line in output.splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        proc_pid, parent_pid, rss = (int(part) for part in parts)
        children.setdefault(parent_pid, []).append(proc_pid)
        rss_kb[proc_pid] = rss

    total = 0
    stack = [pid]
    while stack:
        current = stack.pop()
        total += rss_kb.get(current, 0)
        stack.extend(children.get(current, ()))
    return total * 1024


def _process_tree_rss_bytes_proc(pid: int) -> int | None:
    proc = Path("/proc")
    if not proc.exists():
        return _process_tree_rss_bytes_ps(pid)

    children: dict[int, list[int]] = {}
    rss_bytes: dict[int, int] = {}
    page_size = os.sysconf("SC_PAGE_SIZE")
    for stat_path in proc.glob("[0-9]*/stat"):
        try:
            stat = stat_path.read_text().split()
            proc_pid = int(stat[0])
            parent_pid = int(stat[3])
            rss_pages = int(stat[23])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(parent_pid, []).append(proc_pid)
        rss_bytes[proc_pid] = rss_pages * page_size

    total = 0
    stack = [pid]
    while stack:
        current = stack.pop()
        total += rss_bytes.get(current, 0)
        stack.extend(children.get(current, ()))
    return total


if __name__ == "__main__":
    raise SystemExit(main())
