from __future__ import annotations

import importlib.util
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


def test_dependency_installer_uses_upstream_symbolica_dev_and_skips_gammaloop_by_default() -> None:
    installer = _load_installer_module()

    assert installer.SYMBOLICA_URL == "https://github.com/symbolica-dev/symbolica.git"
    assert installer.SYMBOLICA_REF == "dev"
    assert installer.SYMJIT_VERSION == "2.19.3"
    assert installer.SYMJIT_REF == "dee77e14c8bf9c8d5304dbda2e393c52ee3e4cd4"
    assert "import gammaloop" not in installer.smoke_test_code(include_gammaloop=False)


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
    monkeypatch.setattr(installer, "install_wheel", fake_install)
    monkeypatch.setattr(installer, "ensure_gammaloop_api_absent", lambda: None)

    installer.build_wheels_and_install(include_gammaloop=False)

    assert ("build", "symbolica-community", "symbolica", True) in calls
    assert ("build", "rusticol", "rusticol", True) in calls
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
