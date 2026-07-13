#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
AMPLICOL_ROOT = REPO_ROOT.parent
DEPS_DIR = Path(__file__).resolve().parent
VENV_DIR = DEPS_DIR / ".venv"
WHEEL_DIR = DEPS_DIR / "wheels"
DEPENDENCY_MANIFEST = DEPS_DIR / "install_manifest.json"
PATCHES_DIR = DEPS_DIR / "patches"
SYMBOLICA_WHEEL_DIR = WHEEL_DIR / "symbolica"
GAMMALOOP_WHEEL_DIR = WHEEL_DIR / "gammaloop"
RUSTICOL_WHEEL_DIR = WHEEL_DIR / "rusticol"
UFO_MODEL_LOADER_WHEEL_DIR = WHEEL_DIR / "ufo-model-loader"

SYMBOLICA_COMMUNITY_DIR = DEPS_DIR / "symbolica-community"
SYMBOLICA_DIR = DEPS_DIR / "symbolica"
SYMJIT_DIR = DEPS_DIR / "symjit"
GAMMALOOP_DIR = DEPS_DIR / "gammaloop"
XSIMD_DIR = DEPS_DIR / "xsimd"
UFO_MODEL_LOADER_DIR = DEPS_DIR / "ufo-model-loader"
RUSTICOL_DIR = AMPLICOL_ROOT / "rusticol"

SYMBOLICA_COMMUNITY_URL = "https://github.com/symbolica-dev/symbolica-community.git"
SYMBOLICA_URL = "https://github.com/symbolica-dev/symbolica.git"
SYMBOLICA_REF = "dev"
SYMJIT_URL = "https://github.com/siravan/symjit.git"
SYMJIT_REF = "7fb09d1cb2a943c25a6fd71a208af44fcc6d813d"
SYMJIT_REV = "7fb09d1cb2a943c25a6fd71a208af44fcc6d813d"
SYMJIT_VERSION = "2.19.3"
GAMMALOOP_URL = "https://github.com/alphal00p/gammaloop.git"
XSIMD_URL = "https://github.com/xtensor-stack/xsimd.git"
UFO_MODEL_LOADER_URL = "https://github.com/alphal00p/ufo_model_loader.git"
UFO_MODEL_LOADER_REV = "9cb4deeae40ddd64184049af07ac1d03ce5f6162"
UFO_MODEL_LOADER_VERSION = "0.1.7"
SYMBOLICA_COMMUNITY_REF = "main"
XSIMD_REF = "master"
DEFAULT_GAMMALOOP_REF = "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08"

MIN_RUST_VERSION = (1, 89, 0)
MIN_PYTHON_VERSION = (3, 11)
REEXEC_SENTINEL = "PYAMPLICOL_REEXECED_OUTSIDE_VENV"
BOOTSTRAP_REQUIREMENTS = (
    "pip",
    "maturin",
    "pytest",
    "mypy",
    "numpy<2.5",
    "ipykernel",
    "colorama",
    "progressbar2",
    "tabled",
    "setuptools>=70",
    "wheel",
)
BOOTSTRAP_IMPORTS = (
    "pip",
    "maturin",
    "pytest",
    "mypy",
    "numpy",
    "ipykernel",
    "colorama",
    "progressbar",
    "tabled",
    "setuptools",
    "wheel",
)

MANAGED_PATHS = (
    VENV_DIR,
    WHEEL_DIR,
    DEPENDENCY_MANIFEST,
    SYMBOLICA_COMMUNITY_DIR,
    SYMBOLICA_DIR,
    SYMJIT_DIR,
    GAMMALOOP_DIR,
    XSIMD_DIR,
    UFO_MODEL_LOADER_DIR,
)


class DependencySetupError(RuntimeError):
    pass


def quote_command(command: list[str | Path]) -> str:
    return shlex.join(str(part) for part in command)


def run(
    command: list[str | Path],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    capture: bool = False,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    if not capture:
        prefix = f" ({cwd})" if cwd else ""
        print(f"$ {quote_command(command)}{prefix}")

    completed = subprocess.run(
        [str(part) for part in command],
        cwd=str(cwd) if cwd else None,
        env=env,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
    )

    if check and completed.returncode != 0:
        if capture:
            if completed.stdout:
                print(completed.stdout, end="")
            if completed.stderr:
                print(completed.stderr, end="", file=sys.stderr)
        raise DependencySetupError(
            f"Command failed with exit code {completed.returncode}: {quote_command(command)}"
        )

    return completed


def is_inside_virtualenv() -> bool:
    return bool(os.environ.get("VIRTUAL_ENV")) or sys.prefix != sys.base_prefix


def base_python_executable() -> Path:
    candidate = getattr(sys, "_base_executable", None)
    if candidate and Path(candidate).exists():
        return Path(candidate)
    return Path(sys.executable)


def reexec_outside_virtualenv_if_needed(argv: list[str]) -> None:
    if not is_inside_virtualenv() or os.environ.get(REEXEC_SENTINEL):
        return

    base_python = base_python_executable()
    if base_python.resolve() == Path(sys.executable).resolve():
        return

    env = os.environ.copy()
    env.pop("VIRTUAL_ENV", None)
    env.pop("PYTHONHOME", None)
    env[REEXEC_SENTINEL] = "1"
    os.execvpe(str(base_python), [str(base_python), str(Path(__file__).resolve()), *argv], env)


def venv_bin_dir() -> Path:
    return VENV_DIR / ("Scripts" if os.name == "nt" else "bin")


def venv_python() -> Path:
    exe = "python.exe" if os.name == "nt" else "python"
    return venv_bin_dir() / exe


def venv_environment() -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONHOME", None)
    env["VIRTUAL_ENV"] = str(VENV_DIR)
    env["PATH"] = str(venv_bin_dir()) + os.pathsep + env.get("PATH", "")
    env["PYTHON_BIN_PATH"] = str(venv_python())
    return env


def require_tool(name: str) -> None:
    if shutil.which(name) is None:
        raise DependencySetupError(
            f"`{name}` is required but was not found on PATH. Install it, then rerun this script."
        )


def parse_rust_version(output: str) -> tuple[int, int, int] | None:
    match = re.search(r"rustc\s+(\d+)\.(\d+)\.(\d+)", output)
    if not match:
        return None
    return tuple(int(part) for part in match.groups())


def format_version(version: tuple[int, ...]) -> str:
    return ".".join(str(part) for part in version)


def resolved_executable(command: str) -> str | None:
    if os.sep in command:
        return command if Path(command).exists() else None
    return shutil.which(command)


def command_version(command: str) -> str | None:
    completed = run([command, "--version"], capture=True, check=False)
    if completed.returncode != 0:
        return None
    output = (completed.stdout or completed.stderr).strip().splitlines()
    return output[0] if output else None


def system_prerequisite_issues() -> list[str]:
    issues: list[str] = []
    current_python = sys.version_info[:3]
    if current_python < MIN_PYTHON_VERSION:
        issues.append(
            "Python "
            f"{format_version(MIN_PYTHON_VERSION)} or newer is required "
            f"for pyAmpliCol/rusticol, but this interpreter is "
            f"{format_version(current_python)} at {sys.executable}."
        )

    for module in ("venv", "ensurepip"):
        if importlib.util.find_spec(module) is None:
            issues.append(
                f"Python module `{module}` is required to create/bootstrap "
                f"{display_path(VENV_DIR)}. Install a Python distribution that "
                f"ships `{module}` or install the corresponding OS package."
            )

    required_tools = (
        (
            "git",
            "clone and update Symbolica, symbolica-community, GammaLoop, and xsimd",
            "Install Git and make sure `git` is on PATH.",
        ),
        (
            "cargo",
            "build the Symbolica community and rusticol PyO3 wheels",
            "Install Rust with rustup or your system package manager.",
        ),
        (
            "rustc",
            "compile Symbolica and rusticol",
            "Install Rust with rustup or your system package manager.",
        ),
        (
            "make",
            "build and validate the Fortran AmpliCol reference code",
            "Install GNU Make or a compatible make implementation.",
        ),
        (
            "gfortran",
            "compile the Fortran AmpliCol reference and generated libraries",
            "Install GCC/GFortran, for example with your system package manager.",
        ),
    )
    for name, purpose, hint in required_tools:
        if shutil.which(name) is None:
            issues.append(f"`{name}` is required to {purpose}. {hint}")

    cxx = os.environ.get("CXX", "g++")
    if resolved_executable(cxx) is None:
        if "CXX" in os.environ:
            issues.append(
                f"`CXX={cxx}` does not resolve to an executable C++ compiler. "
                "Set CXX to a valid compiler, or unset it and install `g++`."
            )
        else:
            issues.append(
                "`g++` is required to compile legacy C++ shims and Symbolica "
                "compiled-complex evaluator code. Install g++, or set CXX to "
                "a valid C++ compiler before running the installer."
            )

    if shutil.which("rustc") is not None:
        rustc = run(["rustc", "--version"], capture=True, check=False)
        if rustc.returncode != 0:
            issues.append(
                "`rustc --version` failed; make sure the Rust compiler works "
                "before running the installer."
            )
        else:
            rust_version = parse_rust_version(rustc.stdout)
            if rust_version is None:
                issues.append(
                    f"Could not parse rustc version from {rustc.stdout.strip()!r}."
                )
            elif rust_version < MIN_RUST_VERSION:
                issues.append(
                    "Rust "
                    f"{format_version(MIN_RUST_VERSION)} or newer is required, "
                    f"but rustc reports {format_version(rust_version)}. "
                    "Update Rust with `rustup update stable`, then rerun this script."
                )
    return issues


def ensure_system_prerequisites() -> None:
    issues = system_prerequisite_issues()
    if issues:
        formatted = "\n".join(f"  - {issue}" for issue in issues)
        raise DependencySetupError(
            "Missing or incompatible system prerequisites:\n"
            f"{formatted}\n\n"
            "After installing the missing tools, rerun this script."
        )

    print("System prerequisite check passed:")
    for command in ("git", "cargo", "rustc", "make", "gfortran"):
        version = command_version(command)
        print(f"  {command}: {version or shutil.which(command)}")
    cxx = os.environ.get("CXX", "g++")
    cxx_source = "CXX" if "CXX" in os.environ else "default"
    print(f"  C++ compiler ({cxx_source}): {command_version(cxx) or cxx}")
    print(f"  maturin: {maturin_status()}")
    print(f"  python: {format_version(sys.version_info[:3])} ({sys.executable})")


def create_venv() -> None:
    if venv_python().exists():
        return

    DEPS_DIR.mkdir(parents=True, exist_ok=True)
    run([base_python_executable(), "-m", "venv", str(VENV_DIR)])


def ensure_maturin() -> None:
    create_venv()
    env = venv_environment()
    run([venv_python(), "-m", "pip", "install", "--upgrade", *BOOTSTRAP_REQUIREMENTS], env=env)
    run([venv_python(), "-m", "maturin", "--version"], env=env)


def remove_path(path: Path) -> None:
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def reset_managed_state() -> None:
    print("Removing managed dependency checkouts, wheel output, and venv.")
    for path in MANAGED_PATHS:
        remove_path(path)


def clone_branch(url: str, branch: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    run(["git", "clone", "--depth", "1", "--branch", branch, url, destination])


def clone_ref(url: str, ref: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    run(["git", "clone", "--filter=blob:none", "--no-checkout", url, destination])
    run(["git", "checkout", ref], cwd=destination)


def ensure_branch_checkout(
    url: str,
    branch: str,
    destination: Path,
    *,
    update_existing: bool,
) -> None:
    if not destination.exists():
        clone_branch(url, branch, destination)
        return

    if update_existing:
        run(["git", "fetch", "--depth", "1", "origin", branch], cwd=destination)
        run(["git", "checkout", branch], cwd=destination)
        run(["git", "reset", "--hard", f"origin/{branch}"], cwd=destination)


def ensure_ref_checkout(
    url: str,
    ref: str,
    destination: Path,
    *,
    update_existing: bool,
) -> None:
    if not destination.exists():
        clone_ref(url, ref, destination)
        return

    if update_existing:
        run(["git", "fetch", "--depth", "1", "origin", ref], cwd=destination)
        run(["git", "checkout", ref], cwd=destination)


def discover_gammaloop_ref() -> str:
    cargo_toml = SYMBOLICA_COMMUNITY_DIR / "Cargo.toml"
    if not cargo_toml.exists():
        return DEFAULT_GAMMALOOP_REF

    text = cargo_toml.read_text(encoding="utf-8")
    for line in text.splitlines():
        if "github.com/alphal00p/gammaloop" in line and "rev" in line:
            match = re.search(r'rev\s*=\s*"([^"]+)"', line)
            if match:
                return match.group(1)

    return DEFAULT_GAMMALOOP_REF


def ensure_sources(*, update_existing: bool) -> None:
    ensure_branch_checkout(
        SYMBOLICA_COMMUNITY_URL,
        SYMBOLICA_COMMUNITY_REF,
        SYMBOLICA_COMMUNITY_DIR,
        update_existing=update_existing,
    )
    ensure_branch_checkout(
        SYMBOLICA_URL,
        SYMBOLICA_REF,
        SYMBOLICA_DIR,
        update_existing=update_existing,
    )
    ensure_ref_checkout(
        SYMJIT_URL,
        SYMJIT_REV,
        SYMJIT_DIR,
        update_existing=update_existing,
    )
    ensure_ref_checkout(
        GAMMALOOP_URL,
        discover_gammaloop_ref(),
        GAMMALOOP_DIR,
        update_existing=update_existing,
    )
    ensure_branch_checkout(
        XSIMD_URL,
        XSIMD_REF,
        XSIMD_DIR,
        update_existing=update_existing,
    )
    ensure_ref_checkout(
        UFO_MODEL_LOADER_URL,
        UFO_MODEL_LOADER_REV,
        UFO_MODEL_LOADER_DIR,
        update_existing=update_existing,
    )


def xsimd_headers_are_ready() -> bool:
    return (XSIMD_DIR / "include" / "xsimd" / "xsimd.hpp").exists()


def replace_toml_section(text: str, section_name: str, replacement_body: str) -> str:
    header = f"[{section_name}]"
    pattern = re.compile(
        rf"(?ms)^{re.escape(header)}\n.*?(?=^\[[^\n]+\]\n|\Z)"
    )
    replacement = f"{header}\n{replacement_body.strip()}\n\n"
    if pattern.search(text):
        return pattern.sub(replacement, text, count=1)
    return text.rstrip() + "\n\n" + replacement


def patch_symbolica_community_cargo() -> None:
    cargo_toml = SYMBOLICA_COMMUNITY_DIR / "Cargo.toml"
    text = cargo_toml.read_text(encoding="utf-8")

    dependencies = """
example_extension = { path = "example_extension" }
idenso = { path = "../gammaloop/crates/idenso", features = ["bincode", "python"] }
spenso = { path = "../gammaloop/crates/spenso", features = ["shadowing", "python"] }
spynso3 = { path = "../gammaloop/crates/spynso3" }
symbolica = { path = "../symbolica", features = ["python_export"] }
pyo3 = { version = "0.28", features = ["abi3"] }
pyo3-stub-gen = { version = "0.17", optional = true, default-features = false, features = [
    "numpy",
] }
mimalloc = { version = "0.1", features = [
    "local_dynamic_tls",
] } # prevent TLS allocation errors in conjunction with numpy
"""

    patches = """
graphica = { path = "../symbolica/lib/graphica" }
idenso = { path = "../gammaloop/crates/idenso" }
linnet = { path = "../gammaloop/crates/linnet" }
linnest = { path = "../gammaloop/crates/linnest" }
numerica = { path = "../symbolica/lib/numerica" }
spenso = { path = "../gammaloop/crates/spenso" }
spenso-hep-lib = { path = "../gammaloop/crates/spenso-hep-lib" }
spenso-macros = { path = "../gammaloop/crates/spenso-macros" }
spynso3 = { path = "../gammaloop/crates/spynso3" }
symbolica = { path = "../symbolica" }
"""

    text = replace_toml_section(text, "dependencies", dependencies)
    text = replace_toml_section(text, "patch.crates-io", patches)
    text = text.replace('    "vakint/python_stubgen",\n', "")
    text = re.sub(
        r'(?m)^numerica\s*=\s*\{[^\n]*\}\s*$',
        'numerica = { path = "../symbolica/lib/numerica" }',
        text,
        count=1,
    )
    cargo_toml.write_text(text.rstrip() + "\n", encoding="utf-8")

    example_cargo_toml = SYMBOLICA_COMMUNITY_DIR / "example_extension" / "Cargo.toml"
    example_text = example_cargo_toml.read_text(encoding="utf-8")
    example_text = re.sub(
        r'(?m)^symbolica\s*=\s*\{[^\n]*\}\s*$',
        'symbolica = { path = "../../symbolica", features = ["python_export"] }',
        example_text,
        count=1,
    )
    example_cargo_toml.write_text(example_text.rstrip() + "\n", encoding="utf-8")


def patch_gammaloop_cargo() -> None:
    cargo_toml = GAMMALOOP_DIR / "Cargo.toml"
    text = cargo_toml.read_text(encoding="utf-8")

    text = re.sub(
        r'(?m)^symbolica\s*=\s*\{[^\n]*\}\s*$',
        'symbolica = { path = "../symbolica", default-features = false, features = ["gmp"] }',
        text,
        count=1,
    )

    patches = """
graphica = { path = "../symbolica/lib/graphica" }
numerica = { path = "../symbolica/lib/numerica" }
symbolica = { path = "../symbolica" }
"""

    text = replace_toml_section(text, "patch.crates-io", patches)
    cargo_toml.write_text(text.rstrip() + "\n", encoding="utf-8")


def dependency_patch_specs() -> tuple[tuple[str, Path, Path], ...]:
    return (
        ("gammaloop", GAMMALOOP_DIR, PATCHES_DIR / "gammaloop"),
        ("symjit", SYMJIT_DIR, PATCHES_DIR / "symjit"),
        ("symbolica", SYMBOLICA_DIR, PATCHES_DIR / "symbolica"),
        ("symbolica_community", SYMBOLICA_COMMUNITY_DIR, PATCHES_DIR / "symbolica-community"),
    )


def dependency_patch_manifest_entries() -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    for dependency, _target_dir, patch_dir in dependency_patch_specs():
        if not patch_dir.exists():
            continue
        for patch_path in sorted(patch_dir.glob("*.patch")):
            entries.append(
                {
                    "dependency": dependency,
                    "path": display_path(patch_path),
                    "sha256": hashlib.sha256(patch_path.read_bytes()).hexdigest(),
                }
            )
    return entries


def dependency_patch_is_recorded(patch_path: Path) -> bool:
    if not DEPENDENCY_MANIFEST.exists():
        return False
    try:
        manifest = json.loads(DEPENDENCY_MANIFEST.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    entries = manifest.get("dependency_patches")
    if not isinstance(entries, list):
        return False
    expected_path = display_path(patch_path)
    expected_sha256 = hashlib.sha256(patch_path.read_bytes()).hexdigest()
    for raw_entry in entries:
        if not isinstance(raw_entry, dict) or raw_entry.get("path") != expected_path:
            continue
        recorded_sha256 = raw_entry.get("sha256")
        return recorded_sha256 is None or recorded_sha256 == expected_sha256
    return False


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        try:
            return os.path.relpath(path, REPO_ROOT)
        except ValueError:
            return str(path)


def apply_dependency_patch(target_dir: Path, patch_path: Path) -> bool:
    check = run(
        ["git", "apply", "--check", patch_path],
        cwd=target_dir,
        capture=True,
        check=False,
    )
    if check.returncode == 0:
        run(["git", "apply", patch_path], cwd=target_dir)
        print(f"Applied dependency patch: {display_path(patch_path)}")
        return True

    reverse = run(
        ["git", "apply", "--reverse", "--check", patch_path],
        cwd=target_dir,
        capture=True,
        check=False,
    )
    if reverse.returncode == 0:
        print(f"Dependency patch already applied: {display_path(patch_path)}")
        return False

    if dependency_patch_is_recorded(patch_path):
        print(
            "Dependency patch recorded as applied despite overlapping later "
            f"edits: {display_path(patch_path)}"
        )
        return False

    diagnostics = "".join(
        part
        for part in (check.stdout, check.stderr, reverse.stdout, reverse.stderr)
        if part
    ).strip()
    raise DependencySetupError(
        f"Could not apply dependency patch {display_path(patch_path)}"
        + (f":\n{diagnostics}" if diagnostics else ".")
    )


def apply_dependency_patches() -> None:
    for _dependency, target_dir, patch_dir in dependency_patch_specs():
        if not patch_dir.exists():
            continue
        if not target_dir.exists():
            raise DependencySetupError(
                f"Cannot apply patches from {patch_dir.relative_to(REPO_ROOT)} because "
                f"{target_dir.relative_to(REPO_ROOT)} is missing."
            )
        for patch_path in sorted(patch_dir.glob("*.patch")):
            apply_dependency_patch(target_dir, patch_path)


def git_head(path: Path) -> str:
    completed = run(["git", "rev-parse", "HEAD"], cwd=path, capture=True)
    return completed.stdout.strip()


def optional_git_head(path: Path) -> str | None:
    try:
        return git_head(path)
    except DependencySetupError:
        return None


def write_dependency_manifest(
    *,
    gammaloop_requested: bool,
    symbolica_installed: bool,
    gammaloop_installed: bool,
) -> None:
    manifest = {
        "schema_version": 1,
        "symbolica": {
            "requested": True,
            "installed": symbolica_installed,
            "source_ref": SYMBOLICA_REF,
            "source_path": str(SYMBOLICA_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(SYMBOLICA_DIR),
            "source_url": SYMBOLICA_URL,
            "symjit_version": SYMJIT_VERSION,
            "symjit_source_ref": SYMJIT_REF,
            "symjit_source_path": str(SYMJIT_DIR.relative_to(REPO_ROOT)),
            "symjit_source_rev": optional_git_head(SYMJIT_DIR),
            "symjit_source_url": SYMJIT_URL,
        },
        "symbolica_community": {
            "requested": True,
            "installed": symbolica_installed,
            "source_ref": SYMBOLICA_COMMUNITY_REF,
            "source_path": str(SYMBOLICA_COMMUNITY_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(SYMBOLICA_COMMUNITY_DIR),
            "source_url": SYMBOLICA_COMMUNITY_URL,
        },
        "gammaloop": {
            "requested": gammaloop_requested,
            "installed": gammaloop_installed,
            "source_path": str(GAMMALOOP_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(GAMMALOOP_DIR),
            "features": ["ufo_support", "python_abi"] if gammaloop_requested else [],
            "symbolica_features": ["gmp"] if gammaloop_requested else [],
        },
        "xsimd": {
            "requested": True,
            "installed": xsimd_headers_are_ready(),
            "source_ref": XSIMD_REF,
            "source_path": str(XSIMD_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(XSIMD_DIR),
            "source_url": XSIMD_URL,
            "usage": "header-only include path for Symbolica complex_4x compiled evaluators",
        },
        "ufo_model_loader": {
            "requested": True,
            "installed": python_package_is_installed("ufo-model-loader"),
            "version": UFO_MODEL_LOADER_VERSION,
            "source_path": str(UFO_MODEL_LOADER_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(UFO_MODEL_LOADER_DIR),
            "source_url": UFO_MODEL_LOADER_URL,
        },
        "rusticol": {
            "requested": True,
            "installed": python_package_is_installed("rusticol"),
            "source_path": display_path(RUSTICOL_DIR),
            "source_rev": optional_git_head(RUSTICOL_DIR),
            "usage": "PyO3 runtime for pyAmpliCol eager-DAG process artifacts",
        },
        "dependency_patches": dependency_patch_manifest_entries(),
    }
    DEPENDENCY_MANIFEST.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def local_versions() -> dict[str, str]:
    gammaloop_head = git_head(GAMMALOOP_DIR)
    return {
        "symbolica": git_head(SYMBOLICA_DIR),
        "spenso": gammaloop_head,
        "idenso": gammaloop_head,
    }


def render_local_versions_block(versions: dict[str, str]) -> str:
    lines = [
        "    // pyamplicol: expose managed source revisions.",
        "    let local_versions = PyDict::new(m.py());",
    ]
    for name, rev in versions.items():
        lines.append(f'    local_versions.set_item("{name}", "{rev}")?;')
    lines.extend(
        [
            '    m.add("LOCAL_VERSIONS", local_versions)?;',
            "    // pyamplicol: end managed source revisions.",
        ]
    )
    return "\n".join(lines)


def patch_symbolica_community_module() -> None:
    lib_rs = SYMBOLICA_COMMUNITY_DIR / "src" / "lib.rs"
    text = lib_rs.read_text(encoding="utf-8")

    text = text.replace(
        "types::{PyAnyMethods, PyModule, PyModuleMethods},",
        "types::{PyAnyMethods, PyDict, PyDictMethods, PyModule, PyModuleMethods},",
    )

    block = render_local_versions_block(local_versions())
    marker_pattern = re.compile(
        r"(?ms)^    // py(?:amplicol|chete): expose managed source revisions\.\n"
        r".*?"
        r"^    // py(?:amplicol|chete): end managed source revisions\.\n"
    )
    if marker_pattern.search(text):
        text = marker_pattern.sub(block + "\n", text, count=1)
    else:
        anchor = "    create_symbolica_module(m)?;\n"
        if anchor not in text:
            raise DependencySetupError(
                f"Could not find `create_symbolica_module(m)?;` in {lib_rs}"
            )
        text = text.replace(anchor, anchor + "\n" + block + "\n", 1)

    text = text.replace(
        "    register_module!(m, vakint::symbolica_community_module::VakintWrapper);\n",
        "",
    )

    lib_rs.write_text(text, encoding="utf-8")


def patch_symbolica_symjit_constraint() -> None:
    """Point Symbolica at the managed SymJIT checkout with AArch64 fixes."""

    cargo_toml = SYMBOLICA_DIR / "Cargo.toml"
    text = cargo_toml.read_text(encoding="utf-8")
    patched, replacements = re.subn(
        r'(?m)^symjit\s*=.*$',
        'symjit = { path = "../symjit" }',
        text,
        count=1,
    )
    if replacements != 1:
        raise DependencySetupError(f"Could not find the symjit dependency in {cargo_toml}")
    cargo_toml.write_text(patched, encoding="utf-8")


def patch_symjit_crate_type() -> None:
    """Make the managed SymJIT checkout usable as a Rust library dependency."""

    cargo_toml = SYMJIT_DIR / "Cargo.toml"
    text = cargo_toml.read_text(encoding="utf-8")
    patched, replacements = re.subn(
        r'(?m)^crate-type\s*=\s*\["cdylib"\]\s*$',
        'crate-type = ["rlib"]',
        text,
        count=1,
    )
    if replacements != 1 and 'crate-type = ["rlib"]' not in text:
        raise DependencySetupError(f"Could not patch the symjit crate type in {cargo_toml}")
    cargo_toml.write_text(patched, encoding="utf-8")


def patch_sources() -> None:
    patch_symbolica_community_cargo()
    patch_symbolica_community_module()
    patch_symjit_crate_type()
    patch_symbolica_symjit_constraint()
    pin_symbolica_symjit(SYMBOLICA_DIR)
    pin_symbolica_symjit(SYMBOLICA_COMMUNITY_DIR)
    patch_gammaloop_cargo()
    apply_dependency_patches()


def pin_symbolica_symjit(source_dir: Path) -> None:
    """Refresh lockfiles after switching Symbolica to the managed SymJIT checkout."""

    run(
        ["cargo", "update", "-p", "symjit"],
        cwd=source_dir,
    )


BASE_SMOKE_TEST = """
import importlib.metadata
import ufo_model_loader
from ufo_model_loader.commands import load_model

if importlib.metadata.version("ufo-model-loader") != "0.1.7":
    raise SystemExit("ufo-model-loader 0.1.7 is required")

# Load once before importing the other native extensions. Raw UFO loading is
# isolated in pyAmpliCol proper; this ordering also avoids mixing native symbol
# registries during this single-process installer smoke test.
loader_model, loader_card = load_model("scalars", "full", False)
if len(loader_model.particles) != 3 or not loader_card:
    raise SystemExit("ufo-model-loader bundled-model smoke test failed")

import symbolica
import symbolica.community.idenso
import symbolica.community.spenso
import rusticol

runtime = getattr(rusticol, "Runtime", None)
if runtime is None:
    raise SystemExit("rusticol.Runtime is missing")

required_runtime_methods = {
    "load",
    "evaluate",
    "evaluate_with_prec",
    "profile",
    "stage_diagnostics",
    "metadata",
}
missing_runtime_methods = [
    name for name in sorted(required_runtime_methods)
    if not hasattr(runtime, name)
]
if missing_runtime_methods:
    raise SystemExit(
        f"rusticol.Runtime is missing methods: {missing_runtime_methods}"
    )

versions = getattr(symbolica, "LOCAL_VERSIONS", None)
if not isinstance(versions, dict):
    raise SystemExit("symbolica.LOCAL_VERSIONS is missing or is not a dict")

expected = {"symbolica", "spenso", "idenso"}
missing = expected.difference(versions)
if missing:
    raise SystemExit(f"symbolica.LOCAL_VERSIONS is missing keys: {sorted(missing)}")
"""


def expected_local_versions_for_smoke() -> dict[str, str]:
    try:
        return local_versions()
    except DependencySetupError:
        return {}


def smoke_test_code(*, include_gammaloop: bool) -> str:
    code = BASE_SMOKE_TEST
    expected_versions = expected_local_versions_for_smoke()
    if expected_versions:
        code += (
            "\nexpected_versions = "
            + json.dumps(expected_versions, sort_keys=True)
            + "\n"
            + """
mismatched_versions = {
    name: (expected, versions.get(name))
    for name, expected in expected_versions.items()
    if versions.get(name) != expected
}
if mismatched_versions:
    raise SystemExit(
        f"symbolica.LOCAL_VERSIONS does not match managed sources: {mismatched_versions}"
    )
"""
        )
    if include_gammaloop:
        code += "\nimport gammaloop\n"
    return code


def python_package_is_installed(package: str) -> bool:
    if not venv_python().exists():
        return False

    code = (
        "import importlib.metadata as metadata\n"
        f"metadata.version({package!r})\n"
    )
    completed = run(
        [venv_python(), "-c", code],
        env=venv_environment(),
        capture=True,
        check=False,
    )
    return completed.returncode == 0


def ensure_gammaloop_api_absent() -> None:
    if not python_package_is_installed("gammaloop"):
        return

    print("Removing GammaLoop Python API because it was not requested.")
    run([venv_python(), "-m", "pip", "uninstall", "-y", "gammaloop"], env=venv_environment())


def installed_environment_is_ready(
    *,
    include_gammaloop: bool,
    verbose: bool = False,
) -> bool:
    if not venv_python().exists():
        return False
    if not include_gammaloop and python_package_is_installed("gammaloop"):
        return False

    completed = run(
        [venv_python(), "-c", smoke_test_code(include_gammaloop=include_gammaloop)],
        env=venv_environment(),
        capture=True,
        check=False,
    )
    if verbose:
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr)
    return completed.returncode == 0


def bootstrap_tools_are_ready(*, verbose: bool = False) -> bool:
    if not venv_python().exists():
        return False

    code = "\n".join(f"import {package}" for package in BOOTSTRAP_IMPORTS)
    completed = run(
        [venv_python(), "-c", code],
        env=venv_environment(),
        capture=True,
        check=False,
    )
    if verbose:
        if completed.stdout:
            print(completed.stdout, end="")
        if completed.stderr:
            print(completed.stderr, end="", file=sys.stderr)
    if completed.returncode != 0:
        return False
    maturin = run(
        [venv_python(), "-m", "maturin", "--version"],
        env=venv_environment(),
        capture=True,
        check=False,
    )
    if verbose:
        if maturin.stdout:
            print(maturin.stdout, end="")
        if maturin.stderr:
            print(maturin.stderr, end="", file=sys.stderr)
    return maturin.returncode == 0


def maturin_status() -> str:
    if not venv_python().exists():
        return f"will be bootstrapped into {display_path(VENV_DIR)}"
    completed = run(
        [venv_python(), "-m", "maturin", "--version"],
        env=venv_environment(),
        capture=True,
        check=False,
    )
    if completed.returncode == 0:
        output = (completed.stdout or completed.stderr).strip().splitlines()
        return output[0] if output else "installed in managed venv"
    return f"will be installed/refreshed in {display_path(VENV_DIR)}"


def build_maturin_wheel(
    project_dir: Path,
    output_dir: Path,
    *,
    manifest_path: Path | None = None,
    features: str | None = None,
    profile: str | None = None,
    release: bool = False,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    for wheel in output_dir.glob("*.whl"):
        wheel.unlink()

    env = venv_environment()
    command: list[str | Path] = [
        venv_python(),
        "-m",
        "maturin",
        "build",
        "--interpreter",
        venv_python(),
        "-o",
        output_dir,
    ]
    if release:
        command.append("--release")
    if profile:
        command.extend(["--profile", profile])
    if features:
        command.extend(["--features", features])
    if manifest_path:
        command.extend(["--manifest-path", manifest_path])

    run(command, cwd=project_dir, env=env)

    wheels = sorted(output_dir.glob("*.whl"), key=lambda path: path.stat().st_mtime)
    if not wheels:
        raise DependencySetupError(f"No wheel was produced in {output_dir}")

    return wheels[-1]


def install_wheel(
    wheel: Path,
    *,
    force_reinstall: bool = False,
    no_deps: bool = False,
) -> None:
    env = venv_environment()
    command: list[str | Path] = [venv_python(), "-m", "pip", "install"]
    if force_reinstall:
        command.append("--force-reinstall")
    if no_deps:
        command.append("--no-deps")
    command.append(wheel)
    run(command, env=env)


def build_python_wheel(project_dir: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    for wheel in output_dir.glob("*.whl"):
        wheel.unlink()
    run(
        [
            venv_python(),
            "-m",
            "pip",
            "wheel",
            "--no-build-isolation",
            "--no-deps",
            "--wheel-dir",
            output_dir,
            project_dir,
        ],
        env=venv_environment(),
    )
    wheels = sorted(output_dir.glob("*.whl"), key=lambda path: path.stat().st_mtime)
    if not wheels:
        raise DependencySetupError(f"No wheel was produced in {output_dir}")
    return wheels[-1]


def build_wheels_and_install(*, include_gammaloop: bool) -> None:
    symbolica_wheel = build_maturin_wheel(
        SYMBOLICA_COMMUNITY_DIR,
        SYMBOLICA_WHEEL_DIR,
        release=True,
    )
    install_wheel(symbolica_wheel, force_reinstall=True, no_deps=True)

    loader_wheel = build_python_wheel(
        UFO_MODEL_LOADER_DIR,
        UFO_MODEL_LOADER_WHEEL_DIR,
    )
    install_wheel(loader_wheel, force_reinstall=True, no_deps=True)

    rusticol_wheel = build_maturin_wheel(
        RUSTICOL_DIR,
        RUSTICOL_WHEEL_DIR,
        release=True,
    )
    install_wheel(rusticol_wheel, force_reinstall=True, no_deps=True)

    if not include_gammaloop:
        print("Skipping GammaLoop API build because it was not requested.")
        ensure_gammaloop_api_absent()
        return

    gammaloop_wheel = build_maturin_wheel(
        GAMMALOOP_DIR,
        GAMMALOOP_WHEEL_DIR,
        manifest_path=Path("crates/gammaloop-api/Cargo.toml"),
        features="ufo_support,python_abi",
        profile="release",
    )
    install_wheel(gammaloop_wheel)
    install_wheel(gammaloop_wheel, force_reinstall=True, no_deps=True)


def print_activation_hint() -> None:
    try:
        relative_venv = VENV_DIR.relative_to(Path.cwd())
    except ValueError:
        relative_venv = Path(display_path(VENV_DIR))
    print()
    print("To use this environment in your shell, run:")
    print(f"  source {relative_venv}/bin/activate")


def run_pyamplicol_smoke() -> None:
    package_dir = REPO_ROOT / "src" / "pyamplicol"
    if not package_dir.exists():
        return

    env = venv_environment()
    env["PYTHONPATH"] = str(REPO_ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    if (package_dir / "__main__.py").exists():
        run([venv_python(), "-m", "pyamplicol"], env=env)
    elif (package_dir / "__init__.py").exists():
        run([venv_python(), "-c", "import pyamplicol"], env=env)


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Set up pyamplicol's local Symbolica community dependency environment."
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--recompile",
        action="store_true",
        help="Rebuild and reinstall the Symbolica community wheel without deleting checkouts or target caches.",
    )
    mode.add_argument(
        "--reset",
        action="store_true",
        help="Delete managed dependency checkouts, wheel output, and the venv before rebuilding from scratch.",
    )
    parser.add_argument(
        "--no-gammaloop",
        action="store_true",
        help="Compatibility flag; the GammaLoop Python API is skipped by default.",
    )
    parser.add_argument(
        "--with-gammaloop",
        action="store_true",
        help="Build and install the optional GammaLoop Python API.",
    )
    args = parser.parse_args(argv)
    if args.no_gammaloop and args.with_gammaloop:
        parser.error("--no-gammaloop and --with-gammaloop cannot be used together.")
    return args


def main(argv: list[str] | None = None) -> int:
    raw_args = sys.argv[1:] if argv is None else argv
    reexec_outside_virtualenv_if_needed(raw_args)
    args = parse_args(raw_args)
    include_gammaloop = args.with_gammaloop

    if args.reset:
        reset_managed_state()

    if not include_gammaloop and not args.reset:
        ensure_gammaloop_api_absent()

    ensure_system_prerequisites()

    if (
        not args.recompile
        and not args.reset
        and installed_environment_is_ready(include_gammaloop=include_gammaloop)
    ):
        if not xsimd_headers_are_ready():
            require_tool("git")
            ensure_branch_checkout(
                XSIMD_URL,
                XSIMD_REF,
                XSIMD_DIR,
                update_existing=False,
            )
        if not UFO_MODEL_LOADER_DIR.exists():
            require_tool("git")
            ensure_ref_checkout(
                UFO_MODEL_LOADER_URL,
                UFO_MODEL_LOADER_REV,
                UFO_MODEL_LOADER_DIR,
                update_existing=False,
            )

        if not bootstrap_tools_are_ready():
            print("pyamplicol native dependencies are installed; refreshing Python development tools.")
            ensure_maturin()

        if not include_gammaloop:
            ensure_gammaloop_api_absent()

        print("pyamplicol dependencies are already installed; nothing else to do.")
        write_dependency_manifest(
            gammaloop_requested=include_gammaloop,
            symbolica_installed=True,
            gammaloop_installed=include_gammaloop,
        )
        print_activation_hint()
        return 0

    if is_inside_virtualenv():
        print(
            "Detected an active Python virtual environment. "
            "This installer will ignore it "
            f"and create/use {VENV_DIR} with {base_python_executable()}."
        )

    ensure_maturin()
    if not include_gammaloop:
        ensure_gammaloop_api_absent()
    ensure_sources(update_existing=args.reset)
    patch_sources()
    write_dependency_manifest(
        gammaloop_requested=include_gammaloop,
        symbolica_installed=False,
        gammaloop_installed=False,
    )
    build_wheels_and_install(include_gammaloop=include_gammaloop)

    if not installed_environment_is_ready(include_gammaloop=include_gammaloop, verbose=True):
        raise DependencySetupError("The installed environment failed the import smoke test.")

    write_dependency_manifest(
        gammaloop_requested=include_gammaloop,
        symbolica_installed=True,
        gammaloop_installed=include_gammaloop,
    )
    run_pyamplicol_smoke()
    print_activation_hint()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except DependencySetupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
