from __future__ import annotations

import io
import logging

import pyamplicol
from pyamplicol.__main__ import main


def test_package_exports_version_and_logger_name() -> None:
    assert pyamplicol.__version__
    assert pyamplicol.LOGGER_NAME == "pyamplicol"


def test_configure_logging_uses_package_logger() -> None:
    stream = io.StringIO()
    logger = pyamplicol.configure_logging(stream=stream, force=True)
    logger.info("ready")

    assert logger.name == "pyamplicol"
    assert "INFO:pyamplicol:ready" in stream.getvalue()

    pyamplicol.disable_logging()
    assert any(isinstance(handler, logging.NullHandler) for handler in logger.handlers)


def test_cli_version_mode_does_not_import_native_dependencies(capsys) -> None:
    assert main(["--version"]) == 0

    captured = capsys.readouterr()
    assert captured.out.strip() == pyamplicol.__version__
