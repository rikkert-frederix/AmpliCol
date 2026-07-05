from __future__ import annotations

import contextlib
import logging
import os
import sys
import threading
import time
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from typing import Any, TextIO

from .logging import DEFAULT_LOG_DATE_FORMAT, configure_logging, get_logger


@dataclass(frozen=True)
class DisplayColumn:
    key: str
    title: str
    align: str = "left"


@dataclass(frozen=True)
class DisplayRow:
    values: dict[str, object]
    style: str | None = None


class CliDisplay:
    """Small color/progress/table facade for command-line paths."""

    def __init__(
        self,
        *,
        enabled: bool = True,
        color: bool | None = None,
        progress: bool | None = None,
        stream: TextIO | None = None,
        log_level: int | str = logging.INFO,
    ) -> None:
        self.enabled = enabled
        self.stream = sys.stdout if stream is None else stream
        self.error_stream = sys.stderr
        self.color = _color_enabled(self.stream) if color is None else bool(color)
        self.progress_enabled = (
            _progress_enabled(self.error_stream) if progress is None else bool(progress)
        )
        self._colorama = _load_colorama()
        self._progressbar = _load_progressbar()
        self._tabled: Any | None = None
        self._tabled_checked = False
        if self._colorama is not None:
            self._colorama.just_fix_windows_console()
        configure_logging(
            log_level,
            fmt="%(asctime)s.%(msecs)03d %(message)s",
            datefmt=DEFAULT_LOG_DATE_FORMAT,
            stream=self.error_stream,
            force=True,
        )
        self.logger = get_logger("cli")

    @property
    def has_tabled(self) -> bool:
        self._ensure_tabled_checked()
        return self._tabled is not None

    def style(self, text: object, style: str | None = None) -> str:
        value = str(text)
        if not self.color or style is None:
            return value
        codes = _STYLE_CODES.get(style, "")
        reset = _STYLE_CODES.get("reset", "")
        return f"{codes}{value}{reset}" if codes else value

    def info(self, message: str) -> None:
        if self.enabled:
            self.logger.info(self.style(message, "cyan"))

    def success(self, message: str) -> None:
        if self.enabled:
            self.logger.info(self.style(message, "green"))

    def warning(self, message: str) -> None:
        if self.enabled:
            self.logger.warning(self.style(message, "yellow"))

    def error(self, message: str) -> None:
        if self.enabled:
            self.logger.error(self.style(message, "red"))

    @contextlib.contextmanager
    def progress(self, label: str, *, metadata: str | None = None) -> Iterator[None]:
        if not self.enabled:
            yield
            return

        title = label if metadata is None else f"{label} [{metadata}]"
        start = time.perf_counter()
        spinner = _ProgressSpinner(
            title,
            display=self,
            progressbar_module=self._progressbar,
        )
        try:
            spinner.start()
            yield
        except Exception:
            elapsed = time.perf_counter() - start
            spinner.finish(success=False, elapsed_s=elapsed)
            raise
        else:
            elapsed = time.perf_counter() - start
            spinner.finish(success=True, elapsed_s=elapsed)

    @contextlib.contextmanager
    def stage_progress(
        self,
        label: str,
        *,
        total: int,
        metadata: str | None = None,
    ) -> Iterator["StageProgress"]:
        if not self.enabled:
            yield StageProgress.disabled()
            return

        title = label if metadata is None else f"{label} [{metadata}]"
        progress = StageProgress(
            title,
            total=max(int(total), 1),
            display=self,
            progressbar_module=self._progressbar,
        )
        try:
            progress.start()
            yield progress
        except Exception:
            progress.finish(success=False)
            raise
        else:
            progress.finish(success=True)

    def table(
        self,
        title: str,
        columns: Sequence[DisplayColumn],
        rows: Sequence[DisplayRow | dict[str, object]],
    ) -> str:
        self._ensure_tabled_checked()
        normalized = [
            row if isinstance(row, DisplayRow) else DisplayRow(dict(row))
            for row in rows
        ]
        return _render_table(title, columns, normalized, self)

    def print_table(
        self,
        title: str,
        columns: Sequence[DisplayColumn],
        rows: Sequence[DisplayRow | dict[str, object]],
    ) -> None:
        print(self.table(title, columns, rows), file=self.stream)

    def _ensure_tabled_checked(self) -> None:
        if self._tabled_checked:
            return
        self._tabled = _load_tabled()
        self._tabled_checked = True


class _ProgressSpinner:
    def __init__(
        self,
        label: str,
        *,
        display: CliDisplay,
        progressbar_module: Any,
    ) -> None:
        self.label = label
        self.display = display
        self.progressbar_module = progressbar_module
        self._bar: Any = None
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._started = False

    def start(self) -> None:
        if not self.display.progress_enabled or self.progressbar_module is None:
            self.display.info(f"{self.label} ...")
            return

        widgets = [
            self.display.style(self.label, "cyan"),
            " ",
            self.progressbar_module.AnimatedMarker(),
            " ",
            self.progressbar_module.Timer(format="elapsed %(elapsed)s"),
        ]
        self._bar = self.progressbar_module.ProgressBar(
            max_value=self.progressbar_module.UnknownLength,
            widgets=widgets,
            fd=self.display.error_stream,
        )
        self._bar.start()
        self._started = True
        self._thread = threading.Thread(target=self._drive, daemon=True)
        self._thread.start()

    def finish(self, *, success: bool, elapsed_s: float) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=1.0)
        if self._started and self._bar is not None:
            self._bar.finish()
        status = "done" if success else "failed"
        color = "green" if success else "red"
        self.display.logger.info(
            "%s",
            self.display.style(f"{self.label} {status} in {elapsed_s:.2f}s", color),
        )

    def _drive(self) -> None:
        assert self._bar is not None
        count = 0
        while not self._stop.wait(0.2):
            count += 1
            with contextlib.suppress(Exception):
                self._bar.update(count)


class StageProgress:
    """Progressbar2-backed stage/counter progress with fixed-width metadata."""

    def __init__(
        self,
        label: str = "",
        *,
        total: int = 1,
        display: CliDisplay | None = None,
        progressbar_module: Any = None,
        enabled: bool = True,
    ) -> None:
        self.label = label
        self.total = max(int(total), 1)
        self.display = display
        self.progressbar_module = progressbar_module
        self.enabled = enabled
        self.current = 0
        self.stage = ""
        self.item = ""
        self.start_time = 0.0
        self._bar: Any = None
        self._started = False
        self._last_logged_stage: str | None = None

    @classmethod
    def disabled(cls) -> "StageProgress":
        return cls(enabled=False)

    def start(self) -> None:
        self.start_time = time.perf_counter()
        if not self.enabled or self.display is None:
            return
        if not self.display.progress_enabled or self.progressbar_module is None:
            self.display.info(self._line("starting"))
            return

        widgets = [
            self.display.style(self.label, "cyan"),
            " ",
            self.progressbar_module.Percentage(),
            " ",
            self.progressbar_module.Bar(marker="#", left="|", right="|"),
            " ",
            self.progressbar_module.SimpleProgress(),
            " ",
            self.progressbar_module.Variable(
                "stage",
                format="stage:{formatted_value}",
                width=18,
            ),
            " ",
            self.progressbar_module.Variable(
                "item",
                format="item:{formatted_value}",
                width=28,
            ),
            " ",
            self.progressbar_module.Timer(format="elapsed %(elapsed)s"),
            " ",
            self.progressbar_module.ETA(),
        ]
        self._bar = self.progressbar_module.ProgressBar(
            max_value=self.total,
            widgets=widgets,
            fd=self.display.error_stream,
        )
        self._bar.start(
            stage=_fixed_field("starting", 18),
            item=_fixed_field("initializing", 28),
        )
        with contextlib.suppress(Exception):
            self._bar.update(
                0,
                force=True,
                stage=_fixed_field("starting", 18),
                item=_fixed_field("initializing", 28),
            )
        self._started = True

    def callback(self, event: dict[str, object]) -> None:
        increment = event.get("increment", 0)
        total = event.get("total")
        self.update(
            stage=_optional_str(event.get("stage")),
            item=_optional_str(event.get("item")),
            increment=increment if isinstance(increment, int) else 0,
            total=total if isinstance(total, int) else None,
        )

    def update(
        self,
        *,
        stage: str | None = None,
        item: str | None = None,
        increment: int = 0,
        total: int | None = None,
    ) -> None:
        if total is not None and total > self.total:
            self.total = int(total)
            if self._bar is not None:
                with contextlib.suppress(Exception):
                    self._bar.max_value = self.total
        if stage is not None:
            self.stage = stage
        if item is not None:
            self.item = item
        if increment:
            self.current = min(self.total, self.current + int(increment))

        if not self.enabled or self.display is None:
            return
        fixed_stage = _fixed_field(self.stage, 18)
        fixed_item = _fixed_field(self.item, 28)
        if self._started and self._bar is not None:
            with contextlib.suppress(Exception):
                self._bar.update(
                    self.current,
                    force=True,
                    stage=fixed_stage,
                    item=fixed_item,
                )
            return
        if self.stage != self._last_logged_stage:
            self._last_logged_stage = self.stage
            self.display.info(self._line("running"))

    def finish(self, *, success: bool) -> None:
        if self.current < self.total:
            self.current = self.total
        elapsed = time.perf_counter() - self.start_time if self.start_time else 0.0
        if self._started and self._bar is not None:
            with contextlib.suppress(Exception):
                self._bar.update(
                    self.total,
                    force=True,
                    stage=_fixed_field(self.stage, 18),
                    item=_fixed_field(self.item, 28),
                )
                self._bar.finish()
        if self.enabled and self.display is not None:
            status = "done" if success else "failed"
            color = "green" if success else "red"
            self.display.logger.info(
                "%s",
                self.display.style(
                    f"{self.label} {status} in {elapsed:.2f}s",
                    color,
                ),
            )

    def _line(self, status: str) -> str:
        return (
            f"{self.label} {status}: "
            f"stage={_fixed_field(self.stage, 18)} "
            f"item={_fixed_field(self.item, 28)} "
            f"{self.current:>6d}/{self.total:<6d}"
        )


def default_display_for_args(args: object) -> CliDisplay:
    json_output = bool(getattr(args, "json", False))
    no_color = _truthy_env("PYAMPLICOL_NO_COLOR")
    force_color = _truthy_env("PYAMPLICOL_FORCE_COLOR")
    no_progress = _truthy_env("PYAMPLICOL_NO_PROGRESS")
    force_progress = _truthy_env("PYAMPLICOL_FORCE_PROGRESS")
    color = True if force_color else False if no_color else None
    progress = True if force_progress else False if no_progress else None
    log_level = os.environ.get("PYAMPLICOL_LOG_LEVEL", "INFO")
    return CliDisplay(
        enabled=not json_output,
        color=color,
        progress=progress,
        log_level=log_level,
    )


def format_measurement(value: float, error: float | None = None, *, unit: str = "") -> str:
    formatted = _format_float(value)
    if error is not None and math_is_finite(error):
        formatted = f"{formatted} +/- {_format_float(error, digits=2)}"
    return f"{formatted} {unit}".strip()


def _render_table(
    title: str,
    columns: Sequence[DisplayColumn],
    rows: Sequence[DisplayRow],
    display: CliDisplay,
) -> str:
    rendered_rows = [
        {
            column.key: str(row.values.get(column.key, ""))
            for column in columns
        }
        for row in rows
    ]
    widths = {
        column.key: max(
            len(column.title),
            *(len(row[column.key]) for row in rendered_rows),
            0,
        )
        for column in columns
    }

    def border(left: str, sep: str, right: str) -> str:
        return (
            left
            + sep.join("-" * (widths[column.key] + 2) for column in columns)
            + right
        )

    def render_line(values: dict[str, str], *, style: str | None = None) -> str:
        cells: list[str] = []
        for column in columns:
            text = values[column.key]
            width = widths[column.key]
            if column.align == "right":
                cell = f" {text.rjust(width)} "
            else:
                cell = f" {text.ljust(width)} "
            cells.append(display.style(cell, style))
        return "|" + "|".join(cells) + "|"

    lines = [display.style(title, "bold")]
    lines.append(border("+", "+", "+"))
    lines.append(
        render_line(
            {column.key: column.title for column in columns},
            style="bold",
        )
    )
    lines.append(border("+", "+", "+"))
    for row, rendered in zip(rows, rendered_rows, strict=True):
        lines.append(render_line(rendered, style=row.style))
    lines.append(border("+", "+", "+"))
    return "\n".join(lines)


def _format_float(value: float, *, digits: int = 4) -> str:
    return f"{float(value):.{digits}g}"


def _optional_str(value: object) -> str | None:
    if value is None:
        return None
    return str(value)


def _fixed_field(value: str, width: int) -> str:
    if len(value) > width:
        return value[: max(width - 1, 0)] + ">"
    return value.ljust(width)


def math_is_finite(value: float) -> bool:
    try:
        return value == value and value not in (float("inf"), float("-inf"))
    except TypeError:
        return False


def _color_enabled(stream: TextIO) -> bool:
    if _truthy_env("PYAMPLICOL_NO_COLOR"):
        return False
    if _truthy_env("PYAMPLICOL_FORCE_COLOR"):
        return True
    return bool(getattr(stream, "isatty", lambda: False)())


def _progress_enabled(stream: TextIO) -> bool:
    if _truthy_env("PYAMPLICOL_NO_PROGRESS"):
        return False
    if _truthy_env("PYAMPLICOL_FORCE_PROGRESS"):
        return True
    return bool(getattr(stream, "isatty", lambda: False)())


def _truthy_env(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in {"1", "true", "yes", "on"}


def _load_colorama() -> Any | None:
    with contextlib.suppress(Exception):
        import colorama  # type: ignore[import-untyped]

        return colorama
    return None


def _load_progressbar() -> Any | None:
    with contextlib.suppress(Exception):
        import progressbar

        return progressbar
    return None


def _load_tabled() -> Any | None:
    with contextlib.suppress(Exception):
        import tabled  # type: ignore[import-untyped]

        return tabled
    return None


_STYLE_CODES = {
    "bold": "\033[1m",
    "cyan": "\033[36m",
    "green": "\033[32m",
    "yellow": "\033[33m",
    "red": "\033[31m",
    "dim": "\033[2m",
    "reset": "\033[0m",
}


__all__ = [
    "CliDisplay",
    "DisplayColumn",
    "DisplayRow",
    "StageProgress",
    "default_display_for_args",
    "format_measurement",
]
