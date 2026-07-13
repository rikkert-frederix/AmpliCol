from __future__ import annotations

import importlib.util
import hashlib
import json
import subprocess
from pathlib import Path


def _load_installer_module():
    path = (
        Path(__file__).resolve().parents[2]
        / "dependencies"
        / "install_dependencies.py"
    )
    spec = importlib.util.spec_from_file_location(
        "pyamplicol_install_dependencies",
        path,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_dependency_smoke_test_checks_rusticol_runtime_api() -> None:
    installer = _load_installer_module()
    smoke = installer.smoke_test_code(include_gammaloop=False)

    assert "import rusticol" in smoke
    assert "import ufo_model_loader" in smoke
    assert 'version("ufo-model-loader") != "0.1.7"' in smoke
    assert 'getattr(rusticol, "Runtime", None)' in smoke
    for method in (
        "load",
        "evaluate",
        "evaluate_with_prec",
        "profile",
        "stage_diagnostics",
        "metadata",
    ):
        assert f'"{method}"' in smoke
    assert "mismatched_versions" in smoke


def test_dependency_manifest_records_rusticol_source() -> None:
    installer = _load_installer_module()

    assert installer.RUSTICOL_DIR.name == "rusticol"
    assert installer.RUSTICOL_WHEEL_DIR.name == "rusticol"
    assert installer.UFO_MODEL_LOADER_DIR.name == "ufo-model-loader"
    assert installer.UFO_MODEL_LOADER_VERSION == "0.1.7"
    assert installer.UFO_MODEL_LOADER_REV == (
        "9cb4deeae40ddd64184049af07ac1d03ce5f6162"
    )


def test_dependency_installer_uses_upstream_symbolica_dev_and_skips_gammaloop_by_default() -> None:
    installer = _load_installer_module()

    assert installer.SYMBOLICA_URL == "https://github.com/symbolica-dev/symbolica.git"
    assert installer.SYMBOLICA_REF == "dev"
    assert installer.SYMJIT_VERSION == "2.19.3"
    assert installer.SYMJIT_REF == "7fb09d1cb2a943c25a6fd71a208af44fcc6d813d"
    assert "import gammaloop" not in installer.smoke_test_code(include_gammaloop=False)


def test_dependency_patch_specs_include_symjit() -> None:
    installer = _load_installer_module()

    specs = {
        dependency: (target_dir, patch_dir)
        for dependency, target_dir, patch_dir in installer.dependency_patch_specs()
    }

    assert specs["symjit"] == (
        installer.SYMJIT_DIR,
        installer.PATCHES_DIR / "symjit",
    )


def test_recorded_dependency_patch_requires_matching_digest(
    monkeypatch,
    tmp_path: Path,
) -> None:
    installer = _load_installer_module()
    patch = tmp_path / "patches" / "fix.patch"
    patch.parent.mkdir()
    patch.write_text("test patch\n", encoding="utf-8")
    manifest = tmp_path / "install_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "dependency_patches": [
                    {
                        "path": "patches/fix.patch",
                        "sha256": hashlib.sha256(patch.read_bytes()).hexdigest(),
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(installer, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(installer, "DEPENDENCY_MANIFEST", manifest)

    assert installer.dependency_patch_is_recorded(patch) is True

    patch.write_text("changed patch\n", encoding="utf-8")
    assert installer.dependency_patch_is_recorded(patch) is False


def test_dependency_prerequisite_report_collects_missing_toolchain(
    monkeypatch,
) -> None:
    installer = _load_installer_module()

    monkeypatch.setattr(installer.shutil, "which", lambda _name: None)

    issues = installer.system_prerequisite_issues()
    text = "\n".join(issues)

    for tool in ("git", "cargo", "rustc", "make", "gfortran", "g++"):
        assert tool in text


def test_dependency_prerequisite_report_validates_rust_version(
    monkeypatch,
) -> None:
    installer = _load_installer_module()

    monkeypatch.setattr(installer.shutil, "which", lambda _name: "/usr/bin/tool")

    def fake_run(command, **_kwargs):
        if command[:2] == ["rustc", "--version"]:
            return subprocess.CompletedProcess(
                command,
                0,
                stdout="rustc 1.70.0 (old)\n",
                stderr="",
            )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(installer, "run", fake_run)

    issues = installer.system_prerequisite_issues()

    assert any("Rust 1.89.0 or newer is required" in issue for issue in issues)


def test_dependency_build_wheels_installs_rusticol_by_default(monkeypatch) -> None:
    installer = _load_installer_module()
    calls: list[tuple[str, str, str, bool | None]] = []

    def fake_build(project_dir, output_dir, **kwargs):
        calls.append(("build", project_dir.name, output_dir.name, kwargs.get("release")))
        return output_dir / f"{project_dir.name}.whl"

    def fake_install(wheel, **_kwargs):
        calls.append(("install", wheel.parent.name, wheel.name, None))

    monkeypatch.setattr(installer, "build_maturin_wheel", fake_build)
    monkeypatch.setattr(installer, "build_python_wheel", fake_build)
    monkeypatch.setattr(installer, "install_wheel", fake_install)
    monkeypatch.setattr(installer, "ensure_gammaloop_api_absent", lambda: None)

    installer.build_wheels_and_install(include_gammaloop=False)

    assert ("build", "symbolica-community", "symbolica", True) in calls
    assert ("build", "ufo-model-loader", "ufo-model-loader", None) in calls
    assert ("build", "rusticol", "rusticol", True) in calls
    assert ("install", "ufo-model-loader", "ufo-model-loader.whl", None) in calls
    assert ("install", "rusticol", "rusticol.whl", None) in calls


def test_bootstrap_tools_ready_requires_maturin_cli(monkeypatch, tmp_path: Path) -> None:
    installer = _load_installer_module()
    python = tmp_path / "python"
    python.touch()
    calls = []

    monkeypatch.setattr(installer, "venv_python", lambda: python)
    monkeypatch.setattr(installer, "venv_environment", lambda: {})

    def fake_run(command, **_kwargs):
        calls.append(command)
        if "maturin" in command:
            return subprocess.CompletedProcess(command, 1, stdout="", stderr="missing")
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(installer, "run", fake_run)

    assert installer.bootstrap_tools_are_ready() is False
    assert any("maturin" in command for command in calls)
