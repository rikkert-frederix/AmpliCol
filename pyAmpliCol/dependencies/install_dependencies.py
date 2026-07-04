#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEPS_DIR = Path(__file__).resolve().parent
VENV_DIR = DEPS_DIR / ".venv"
WHEEL_DIR = DEPS_DIR / "wheels"
DEPENDENCY_MANIFEST = DEPS_DIR / "install_manifest.json"
PATCHES_DIR = DEPS_DIR / "patches"
SYMBOLICA_WHEEL_DIR = WHEEL_DIR / "symbolica"
GAMMALOOP_WHEEL_DIR = WHEEL_DIR / "gammaloop"

SYMBOLICA_COMMUNITY_DIR = DEPS_DIR / "symbolica-community"
SYMBOLICA_DIR = DEPS_DIR / "symbolica"
GAMMALOOP_DIR = DEPS_DIR / "gammaloop"
XSIMD_DIR = DEPS_DIR / "xsimd"

SYMBOLICA_COMMUNITY_URL = "https://github.com/symbolica-dev/symbolica-community.git"
SYMBOLICA_URL = "https://github.com/symbolica-dev/symbolica.git"
GAMMALOOP_URL = "https://github.com/alphal00p/gammaloop.git"
XSIMD_URL = "https://github.com/xtensor-stack/xsimd.git"
DEFAULT_GAMMALOOP_REF = "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08"

MIN_RUST_VERSION = (1, 89, 0)
REEXEC_SENTINEL = "PYAMPLICOL_REEXECED_OUTSIDE_VENV"
BOOTSTRAP_REQUIREMENTS = ("pip", "maturin", "pytest", "mypy", "numpy<2.5", "ipykernel")
BOOTSTRAP_IMPORTS = ("pip", "maturin", "pytest", "mypy", "numpy", "ipykernel")

MANAGED_PATHS = (
    VENV_DIR,
    WHEEL_DIR,
    DEPENDENCY_MANIFEST,
    SYMBOLICA_COMMUNITY_DIR,
    SYMBOLICA_DIR,
    GAMMALOOP_DIR,
    XSIMD_DIR,
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


def ensure_system_prerequisites() -> None:
    require_tool("git")
    require_tool("cargo")
    require_tool("rustc")

    rustc = run(["rustc", "--version"], capture=True)
    rust_version = parse_rust_version(rustc.stdout)
    if rust_version and rust_version < MIN_RUST_VERSION:
        required = ".".join(str(part) for part in MIN_RUST_VERSION)
        current = ".".join(str(part) for part in rust_version)
        raise DependencySetupError(
            f"Rust {required} or newer is required, but rustc reports {current}. "
            "Update Rust with `rustup update stable`, then rerun this script."
        )


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
        "main",
        SYMBOLICA_COMMUNITY_DIR,
        update_existing=update_existing,
    )
    ensure_branch_checkout(
        SYMBOLICA_URL,
        "dev",
        SYMBOLICA_DIR,
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
        "master",
        XSIMD_DIR,
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
vakint = { path = "../gammaloop/crates/vakint", features = [
    "symbolica_community_module",
] }
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
vakint = { path = "../gammaloop/crates/vakint" }
"""

    text = replace_toml_section(text, "dependencies", dependencies)
    text = replace_toml_section(text, "patch.crates-io", patches)
    cargo_toml.write_text(text.rstrip() + "\n", encoding="utf-8")


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
                }
            )
    return entries


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
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
    if not (path / ".git").exists():
        return None
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
            "source_path": str(SYMBOLICA_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(SYMBOLICA_DIR),
        },
        "symbolica_community": {
            "requested": True,
            "installed": symbolica_installed,
            "source_path": str(SYMBOLICA_COMMUNITY_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(SYMBOLICA_COMMUNITY_DIR),
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
            "source_path": str(XSIMD_DIR.relative_to(REPO_ROOT)),
            "source_rev": optional_git_head(XSIMD_DIR),
            "usage": "header-only include path for Symbolica complex_4x compiled evaluators",
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
        "vakint": gammaloop_head,
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

    lib_rs.write_text(text, encoding="utf-8")


def patch_sources() -> None:
    patch_symbolica_community_cargo()
    patch_symbolica_community_module()
    patch_gammaloop_cargo()
    apply_dependency_patches()


BASE_SMOKE_TEST = """
import symbolica
import symbolica.community.idenso
import symbolica.community.spenso
import symbolica.community.vakint

versions = getattr(symbolica, "LOCAL_VERSIONS", None)
if not isinstance(versions, dict):
    raise SystemExit("symbolica.LOCAL_VERSIONS is missing or is not a dict")

expected = {"symbolica", "spenso", "idenso", "vakint"}
missing = expected.difference(versions)
if missing:
    raise SystemExit(f"symbolica.LOCAL_VERSIONS is missing keys: {sorted(missing)}")
"""


def smoke_test_code(*, include_gammaloop: bool) -> str:
    code = BASE_SMOKE_TEST
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
    return completed.returncode == 0


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


def build_wheels_and_install(*, include_gammaloop: bool) -> None:
    symbolica_wheel = build_maturin_wheel(
        SYMBOLICA_COMMUNITY_DIR,
        SYMBOLICA_WHEEL_DIR,
        release=True,
    )
    install_wheel(symbolica_wheel, force_reinstall=True, no_deps=True)

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

    if (
        not args.recompile
        and not args.reset
        and installed_environment_is_ready(include_gammaloop=include_gammaloop)
    ):
        if not xsimd_headers_are_ready():
            require_tool("git")
            ensure_branch_checkout(
                XSIMD_URL,
                "master",
                XSIMD_DIR,
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

    ensure_system_prerequisites()
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
