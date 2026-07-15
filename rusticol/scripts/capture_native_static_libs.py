#!/usr/bin/env python3
from __future__ import annotations

import platform
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "target" / "release" / "rusticol-native-static-libs.txt"


def main() -> int:
    completed = subprocess.run(
        [
            "cargo",
            "rustc",
            "--package",
            "rusticol-capi",
            "--release",
            "--",
            "--print",
            "native-static-libs",
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    )
    matches = re.findall(r"native-static-libs:\s*(.+)", completed.stdout)
    if not matches:
        raise RuntimeError("rustc did not report native-static-libs")
    flags = matches[-1].split()
    if platform.system() == "Darwin" and "-lgcc_s" in flags:
        replacement = macos_libgcc_path()
        if replacement is not None:
            flags = [str(replacement) if flag == "-lgcc_s" else flag for flag in flags]
            flags.insert(0, f"-Wl,-rpath,{replacement.parent}")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(" ".join(flags) + "\n")
    print(f"wrote {OUTPUT}")
    return 0


def macos_libgcc_path() -> Path | None:
    compiler = shutil.which("gfortran") or shutil.which("gcc")
    if compiler is None:
        return None
    completed = subprocess.run(
        [compiler, "-print-file-name=libgcc_s.1.1.dylib"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    candidate = Path(completed.stdout.strip())
    if completed.returncode != 0 or candidate.name == str(candidate) or not candidate.exists():
        return None
    return candidate.resolve()


if __name__ == "__main__":
    raise SystemExit(main())
