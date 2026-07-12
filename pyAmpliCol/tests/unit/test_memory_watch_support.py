from __future__ import annotations

import sys
import signal
from pathlib import Path

import pytest


SCRIPTS_DIR = Path(__file__).resolve().parents[2] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from memory_watch_support import (  # noqa: E402
    ProcessTreeMemoryMonitor,
    ProcessTreeMemorySample,
    parse_physical_footprint_bytes,
)
import memory_watch_support  # noqa: E402
import run_with_memory_watch  # noqa: E402


def test_parse_physical_footprint_prefers_deduplicated_summary() -> None:
    output = """
    phys_footprint: 100 B
    phys_footprint: 200 B
    Summary Footprint: 250 B
    """

    assert parse_physical_footprint_bytes(output) == 250


def test_parse_physical_footprint_sums_single_process_sections() -> None:
    output = """
    Auxiliary data:
        phys_footprint: 100 B
    Auxiliary data:
        phys_footprint: 200 B
    """

    assert parse_physical_footprint_bytes(output) == 300


def test_process_tree_memory_sample_uses_larger_metric() -> None:
    sample = ProcessTreeMemorySample(
        rss_bytes=10,
        physical_footprint_bytes=20,
    )

    assert sample.effective_bytes == 20
    assert sample.effective_metric == "physical_footprint"


def test_monitor_polls_faster_above_physical_high_watermark(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    footprint_calls: list[int] = []
    monkeypatch.setattr(memory_watch_support.sys, "platform", "darwin")
    monkeypatch.setattr(
        memory_watch_support,
        "process_tree_rss_bytes",
        lambda _pid: 10,
    )
    monkeypatch.setattr(
        memory_watch_support,
        "process_tree_physical_footprint_bytes",
        lambda pid: footprint_calls.append(pid) or 80,
    )
    monitor = ProcessTreeMemoryMonitor(
        physical_footprint_poll_s=5.0,
        physical_footprint_near_limit_poll_s=1.0,
        physical_footprint_high_watermark_bytes=50,
    )

    monitor.sample(123, now=0.0)
    monitor.sample(123, now=2.0)

    assert footprint_calls == [123, 123]


def test_watchdog_termination_signals_detached_descendants(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, int, signal.Signals]] = []
    monkeypatch.setattr(
        run_with_memory_watch,
        "process_tree_pids",
        lambda _pid: (100, 101, 102),
    )
    monkeypatch.setattr(
        run_with_memory_watch.os,
        "kill",
        lambda pid, sig: calls.append(("pid", pid, sig)),
    )
    monkeypatch.setattr(
        run_with_memory_watch.os,
        "killpg",
        lambda pid, sig: calls.append(("group", pid, sig)),
    )

    tracked = run_with_memory_watch._terminate_tree(100, signal.SIGTERM)

    assert tracked == (100, 101, 102)
    assert calls == [
        ("pid", 102, signal.SIGTERM),
        ("pid", 101, signal.SIGTERM),
        ("group", 100, signal.SIGTERM),
    ]


def test_watchdog_kills_tracked_descendants_after_root_exits(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[tuple[int, ...], signal.Signals, int | None]] = []

    class FinishedChild:
        pid = 100

        @staticmethod
        def wait(*, timeout: float | None = None) -> int:
            assert timeout == 10.0
            return 0

    monkeypatch.setattr(
        run_with_memory_watch,
        "_signal_pids",
        lambda pids, sig, *, exclude=None: calls.append(
            (tuple(pids), sig, exclude)
        ),
    )

    returncode = run_with_memory_watch._wait_or_kill(
        FinishedChild(),
        (100, 101, 102),
    )

    assert returncode == 0
    assert calls == [((100, 101, 102), signal.SIGKILL, 100)]
