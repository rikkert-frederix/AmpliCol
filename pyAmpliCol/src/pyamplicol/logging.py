from __future__ import annotations

import logging
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from typing import TextIO

LOGGER_NAME = "pyamplicol"
DEFAULT_LOG_FORMAT = "%(asctime)s.%(msecs)03d %(levelname)s:%(name)s:%(message)s"
DEFAULT_LOG_DATE_FORMAT = "%H:%M:%S"

_HANDLER_MARKER = "_pyamplicol_configured_handler"


def get_logger(name: str | None = None) -> logging.Logger:
    """Return the package logger or one of its named children."""

    if name is None or name == "":
        return logging.getLogger(LOGGER_NAME)
    return logging.getLogger(f"{LOGGER_NAME}.{name}")


def configure_logging(
    level: int | str = logging.INFO,
    *,
    stream: TextIO | None = None,
    fmt: str = DEFAULT_LOG_FORMAT,
    datefmt: str | None = DEFAULT_LOG_DATE_FORMAT,
    force: bool = False,
) -> logging.Logger:
    """Configure pyamplicol progress logging."""

    selected_level = _normalize_level(level)
    logger = get_logger()
    if force:
        _remove_configured_handlers(logger)

    handler = _configured_handler(logger)
    if handler is None:
        handler = logging.StreamHandler(sys.stderr if stream is None else stream)
        setattr(handler, _HANDLER_MARKER, True)
        logger.addHandler(handler)
    elif stream is not None and isinstance(handler, logging.StreamHandler):
        handler.setStream(stream)

    handler.setFormatter(logging.Formatter(fmt=fmt, datefmt=datefmt))
    logger.setLevel(selected_level)
    logger.propagate = False
    return logger


def disable_logging() -> None:
    """Remove pyamplicol-installed logging handlers and restore quiet defaults."""

    logger = get_logger()
    _remove_configured_handlers(logger)
    _ensure_null_handler(logger)
    logger.propagate = False


@contextmanager
def progress(
    message: str,
    *,
    logger: logging.Logger | None = None,
    level: int | str = logging.INFO,
    done_message: str | None = None,
) -> Iterator[None]:
    """Log start, finish, and failure messages around an operation."""

    selected_logger = get_logger() if logger is None else logger
    log_level = _normalize_level(level)
    start = time.perf_counter()
    if selected_logger.isEnabledFor(log_level):
        selected_logger.log(log_level, "%s ...", message)
    try:
        yield
    except Exception:
        elapsed = time.perf_counter() - start
        selected_logger.exception("%s failed after %.2fs", message, elapsed)
        raise
    else:
        if selected_logger.isEnabledFor(log_level):
            selected_logger.log(
                log_level,
                "%s done in %.2fs",
                done_message or message,
                time.perf_counter() - start,
            )


def _normalize_level(level: int | str) -> int:
    if isinstance(level, int):
        return level
    value = logging.getLevelName(level.upper())
    if not isinstance(value, int):
        raise ValueError(f"unknown logging level {level!r}")
    return value


def _configured_handler(logger: logging.Logger) -> logging.Handler | None:
    for handler in logger.handlers:
        if getattr(handler, _HANDLER_MARKER, False):
            return handler
    return None


def _remove_configured_handlers(logger: logging.Logger) -> None:
    for handler in tuple(logger.handlers):
        if getattr(handler, _HANDLER_MARKER, False):
            logger.removeHandler(handler)
            try:
                handler.close()
            except ValueError:
                pass


def _ensure_null_handler(logger: logging.Logger) -> None:
    if not any(isinstance(handler, logging.NullHandler) for handler in logger.handlers):
        logger.addHandler(logging.NullHandler())


_ensure_null_handler(get_logger())


__all__ = [
    "DEFAULT_LOG_FORMAT",
    "DEFAULT_LOG_DATE_FORMAT",
    "LOGGER_NAME",
    "configure_logging",
    "disable_logging",
    "get_logger",
    "progress",
]
