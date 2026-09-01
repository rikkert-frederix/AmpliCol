#!/usr/bin/env python3
"""Check that hostile amplitude-library dimensions fail before allocation."""

from __future__ import annotations

import struct
import subprocess
import sys
import tempfile
from pathlib import Path


def write_prefix(path: Path, n_external: int, n_processes: int) -> None:
    with path.open("wb") as stream:
        stream.write(struct.pack("=i", 4))
        stream.write(struct.pack("=9d", *([1.0] * 9)))
        stream.write(struct.pack("=2i", n_external, n_processes))


def require_rejection(executable: Path, fixture: Path, expected: str) -> None:
    result = subprocess.run(
        [str(executable), str(fixture)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode == 0 or expected not in result.stdout:
        raise AssertionError(
            f"fixture {fixture.name} was not rejected as expected\n"
            f"return code: {result.returncode}\noutput:\n{result.stdout}"
        )


def main() -> None:
    executable = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="amplicol-library-metadata-") as tmp:
        root = Path(tmp)

        oversized = root / "oversized.bin"
        write_prefix(oversized, 20, 200_000_000)
        require_rejection(executable, oversized, "unique-process array shape")

        truncated = root / "truncated.bin"
        write_prefix(truncated, 4, 1)
        require_rejection(executable, truncated, "payload is truncated before unique-process map")

        unsupported_multiplicity = root / "two_to_one.bin"
        write_prefix(unsupported_multiplicity, 3, 1)
        require_rejection(
            executable,
            unsupported_multiplicity,
            "Invalid unique-process dimensions",
        )

    print("amplitude-library metadata regression: PASS")


if __name__ == "__main__":
    main()
