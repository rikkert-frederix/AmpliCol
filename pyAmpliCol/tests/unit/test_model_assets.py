from __future__ import annotations

from pathlib import Path

import pyamplicol
from pyamplicol.model_assets import (
    BUNDLED_MODEL_NAMES,
    bundled_model_path,
    packaged_models_root,
    verify_model_asset_manifest,
)


def test_bundled_model_assets_are_complete_and_in_sync() -> None:
    package_root = packaged_models_root()
    repository_root = Path(pyamplicol.__file__).resolve().parents[2]
    source_root = repository_root / "assets" / "models"

    assert verify_model_asset_manifest(package_root) == ()
    assert verify_model_asset_manifest(source_root) == ()
    assert (package_root / "MANIFEST.sha256").read_bytes() == (
        source_root / "MANIFEST.sha256"
    ).read_bytes()

    for model_name in BUNDLED_MODEL_NAMES:
        assert bundled_model_path(model_name, "ufo").is_dir()
        assert bundled_model_path(model_name, "json").is_file()


def test_bundled_model_path_rejects_unknown_values() -> None:
    for name, model_format in (("unknown", "ufo"), ("sm", "yaml")):
        try:
            bundled_model_path(name, model_format)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid bundled model lookup should fail")
