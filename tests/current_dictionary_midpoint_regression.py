#!/usr/bin/env python3
"""Guard the large-current dictionary against binary-search overflow."""

from __future__ import annotations

import ctypes
import re
from pathlib import Path
from typing import Callable


SOURCE = Path(__file__).resolve().parents[1] / "amplitude_QCD.f03"
Z_N9_MAX_KEY = 1_646_119_211


def int32(value: int) -> int:
    return ctypes.c_int32(value).value


def unsafe_midpoint(left: int, right: int) -> int:
    return int32(left + right) // 2


def safe_midpoint(left: int, right: int) -> int:
    return left + (right - left) // 2


def find_last_key(midpoint: Callable[[int, int], int]) -> tuple[int | None, int]:
    left = 1
    right = Z_N9_MAX_KEY
    steps = 0
    while left <= right:
        middle = midpoint(left, right)
        steps += 1
        if middle < 1 or middle > Z_N9_MAX_KEY:
            return None, middle
        if middle == Z_N9_MAX_KEY:
            return middle, steps
        left = middle + 1
    return None, steps


def main() -> int:
    source = re.sub(r"\s+", "", SOURCE.read_text(encoding="ascii").lower())
    if "middle=left+(right-left)/2" not in source:
        raise AssertionError("solve_dict must use an overflow-safe midpoint")
    if "middle=(right+left)/2" in source or "middle=(left+right)/2" in source:
        raise AssertionError("solve_dict still contains an overflowing midpoint")

    unsafe_result, unsafe_failure = find_last_key(unsafe_midpoint)
    if unsafe_result is not None or unsafe_failure >= 1:
        raise AssertionError("the regression setup did not reproduce int32 overflow")

    safe_result, safe_steps = find_last_key(safe_midpoint)
    if safe_result != Z_N9_MAX_KEY:
        raise AssertionError("overflow-safe binary search did not find the final key")

    print(
        "PASS current dictionary midpoint "
        f"max_key={Z_N9_MAX_KEY} steps={safe_steps} "
        f"old_midpoint={unsafe_failure}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
