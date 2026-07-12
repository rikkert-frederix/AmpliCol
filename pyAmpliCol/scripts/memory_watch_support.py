from __future__ import annotations

import importlib
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


_SUMMARY_FOOTPRINT_RE = re.compile(r"Summary Footprint:\s*(\d+) B")
_PHYSICAL_FOOTPRINT_RE = re.compile(r"phys_footprint:\s*(\d+) B")


@dataclass(frozen=True)
class ProcessTreeMemorySample:
    rss_bytes: int | None
    physical_footprint_bytes: int | None

    @property
    def effective_bytes(self) -> int | None:
        values = tuple(
            value
            for value in (self.rss_bytes, self.physical_footprint_bytes)
            if value is not None
        )
        return max(values) if values else None

    @property
    def effective_metric(self) -> str | None:
        effective = self.effective_bytes
        if effective is None:
            return None
        if self.physical_footprint_bytes == effective:
            return "physical_footprint"
        return "rss"


class ProcessTreeMemoryMonitor:
    """Sample process-tree RSS and macOS compressed physical footprint."""

    def __init__(
        self,
        *,
        physical_footprint_poll_s: float = 5.0,
        physical_footprint_near_limit_poll_s: float = 1.0,
        physical_footprint_high_watermark_bytes: int | None = None,
    ) -> None:
        if physical_footprint_poll_s <= 0:
            raise ValueError("physical_footprint_poll_s must be positive")
        if physical_footprint_near_limit_poll_s <= 0:
            raise ValueError(
                "physical_footprint_near_limit_poll_s must be positive"
            )
        self.physical_footprint_poll_s = float(physical_footprint_poll_s)
        self.physical_footprint_near_limit_poll_s = float(
            physical_footprint_near_limit_poll_s
        )
        self.physical_footprint_high_watermark_bytes = (
            None
            if physical_footprint_high_watermark_bytes is None
            else int(physical_footprint_high_watermark_bytes)
        )
        self.latest_physical_footprint_bytes: int | None = None
        self.peak_rss_bytes: int | None = None
        self.peak_physical_footprint_bytes: int | None = None
        self.peak_memory_bytes: int | None = None
        self.peak_memory_metric: str | None = None
        self.rss_polling_available = True
        self.physical_footprint_supported = (
            sys.platform == "darwin" and Path("/usr/bin/footprint").exists()
        )
        self.physical_footprint_polling_available = self.physical_footprint_supported
        self._last_physical_footprint_poll = float("-inf")

    def sample(
        self,
        pid: int,
        *,
        now: float | None = None,
    ) -> ProcessTreeMemorySample:
        timestamp = time.monotonic() if now is None else float(now)
        rss_bytes = process_tree_rss_bytes(pid)
        if rss_bytes is None:
            self.rss_polling_available = False
        else:
            self.peak_rss_bytes = _maximum(self.peak_rss_bytes, rss_bytes)

        physical_poll_s = self.physical_footprint_poll_s
        if (
            self.physical_footprint_high_watermark_bytes is not None
            and self.latest_physical_footprint_bytes is not None
            and self.latest_physical_footprint_bytes
            >= self.physical_footprint_high_watermark_bytes
        ):
            physical_poll_s = self.physical_footprint_near_limit_poll_s
        if (
            sys.platform == "darwin"
            and timestamp - self._last_physical_footprint_poll
            >= physical_poll_s
        ):
            self._last_physical_footprint_poll = timestamp
            footprint_bytes = process_tree_physical_footprint_bytes(pid)
            if footprint_bytes is not None:
                self.latest_physical_footprint_bytes = footprint_bytes
                self.peak_physical_footprint_bytes = _maximum(
                    self.peak_physical_footprint_bytes,
                    footprint_bytes,
                )
                self.physical_footprint_polling_available = True
            elif self.peak_physical_footprint_bytes is None:
                self.physical_footprint_polling_available = False

        sample = ProcessTreeMemorySample(
            rss_bytes=rss_bytes,
            physical_footprint_bytes=self.latest_physical_footprint_bytes,
        )
        effective_bytes = sample.effective_bytes
        if effective_bytes is not None and (
            self.peak_memory_bytes is None
            or effective_bytes > self.peak_memory_bytes
        ):
            self.peak_memory_bytes = effective_bytes
            self.peak_memory_metric = sample.effective_metric
        return sample


def _maximum(current: int | None, value: int) -> int:
    return value if current is None else max(current, value)


def process_tree_rss_bytes(pid: int) -> int | None:
    psutil_rss = _process_tree_rss_bytes_psutil(pid)
    if psutil_rss is not None:
        return psutil_rss
    if sys.platform == "darwin":
        return _process_tree_rss_bytes_ps(pid)
    return _process_tree_rss_bytes_proc(pid)


def process_tree_physical_footprint_bytes(pid: int) -> int | None:
    if sys.platform != "darwin" or not Path("/usr/bin/footprint").exists():
        return None
    pids = _process_tree_pids(pid)
    if not pids:
        return None
    try:
        completed = subprocess.run(
            [
                "/usr/bin/footprint",
                "--swapped",
                "-f",
                "bytes",
                *(str(process_pid) for process_pid in pids),
            ],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=30.0,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    return parse_physical_footprint_bytes(completed.stdout)


def process_tree_pids(pid: int) -> tuple[int, ...]:
    """Return the root and all descendants, including detached process groups."""

    return _process_tree_pids(pid)


def parse_physical_footprint_bytes(output: str) -> int | None:
    summary_match = _SUMMARY_FOOTPRINT_RE.search(output)
    if summary_match is not None:
        return int(summary_match.group(1))
    values = tuple(
        int(match.group(1)) for match in _PHYSICAL_FOOTPRINT_RE.finditer(output)
    )
    return sum(values) if values else None


def _process_tree_pids(pid: int) -> tuple[int, ...]:
    psutil_pids = _process_tree_pids_psutil(pid)
    if psutil_pids is not None:
        return psutil_pids
    ps_pids = _process_tree_pids_ps(pid)
    if ps_pids is not None:
        return ps_pids
    return (pid,)


def _process_tree_pids_psutil(pid: int) -> tuple[int, ...] | None:
    try:
        psutil = importlib.import_module("psutil")
    except ImportError:
        return None
    try:
        root = psutil.Process(pid)
        return tuple(
            sorted({root.pid, *(child.pid for child in root.children(recursive=True))})
        )
    except Exception:
        return None


def _process_tree_pids_ps(pid: int) -> tuple[int, ...] | None:
    try:
        output = subprocess.run(
            ["ps", "-axo", "pid=,ppid="],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    children: dict[int, list[int]] = {}
    seen: set[int] = set()
    for line in output.splitlines():
        parts = line.split()
        if len(parts) != 2:
            continue
        process_pid, parent_pid = (int(part) for part in parts)
        seen.add(process_pid)
        children.setdefault(parent_pid, []).append(process_pid)
    if pid not in seen:
        return None
    return _walk_process_tree(pid, children)


def _walk_process_tree(
    pid: int,
    children: dict[int, list[int]],
) -> tuple[int, ...]:
    result: list[int] = []
    stack = [pid]
    while stack:
        current = stack.pop()
        result.append(current)
        stack.extend(children.get(current, ()))
    return tuple(sorted(result))


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
        process_pid, parent_pid, rss = (int(part) for part in parts)
        children.setdefault(parent_pid, []).append(process_pid)
        rss_kb[process_pid] = rss
    if pid not in rss_kb:
        return None
    return sum(rss_kb.get(item, 0) for item in _walk_process_tree(pid, children)) * 1024


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
            process_pid = int(stat[0])
            parent_pid = int(stat[3])
            rss_pages = int(stat[23])
        except (OSError, ValueError, IndexError):
            continue
        children.setdefault(parent_pid, []).append(process_pid)
        rss_bytes[process_pid] = rss_pages * page_size
    if pid not in rss_bytes:
        return None
    return sum(
        rss_bytes.get(item, 0) for item in _walk_process_tree(pid, children)
    )
