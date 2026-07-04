from __future__ import annotations

from pathlib import Path

import pytest


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    for item in items:
        path = Path(str(item.path))
        path_parts = path.parts
        if _path_has_parts(path_parts, ("tests", "unit")):
            item.add_marker(pytest.mark.unit)
        if _path_has_parts(path_parts, ("tests", "integration")):
            item.add_marker(pytest.mark.integration)
        if path.name == "test_static_typing.py":
            item.add_marker(pytest.mark.typing)
        if path.name == "test_dependencies.py":
            item.add_marker(pytest.mark.dependencies)


def _path_has_parts(path_parts: tuple[str, ...], expected: tuple[str, ...]) -> bool:
    for index in range(0, len(path_parts) - len(expected) + 1):
        if path_parts[index : index + len(expected)] == expected:
            return True
    return False
