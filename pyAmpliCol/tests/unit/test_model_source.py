from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from pyamplicol.model_assets import bundled_model_path
from pyamplicol.model_source import (
    COMPILED_MODEL_KIND,
    CompiledModel,
    compile_model_source,
    detect_model_source,
    load_compiled_model,
    preflight_model,
)
from pyamplicol.ufo_model import CompiledUFOModel


_BUNDLED_MODEL_STRUCTURE_SHA256 = {
    "scalars": "04f4611426abc396ec28feed5fdb0fa90cfa58cd5bc083ad014b13085fad4827",
    "scalar_gravity": "ffab5001bf75cf33f470a06f7b7a332a517d0d16e0fb7edb472a8d19d158841a",
}


def _compiled_ir_structure_sha256(compiled: CompiledModel) -> str:
    ir = compiled.ir.to_dict()
    for kernel in ir["oriented_kernels"]:
        # Symbolica may print equivalent large sums in a different term order.
        expressions = kernel.pop("component_expressions")
        kernel["component_expression_count"] = len(expressions)
        # Evaluation equivalence is derived optimization metadata and does not
        # change the compiled model's physics/topology structure.
        kernel.pop("evaluation_class")
        kernel.pop("evaluation_factor")
        kernel.pop("evaluation_input_order")
        kernel.pop("evaluation_equivalence_verified")
    payload = json.dumps(
        ir,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


@pytest.mark.parametrize("alias", ["BUILTIN_SM", "builtin_sm", "built-in-sm"])
def test_builtin_sm_aliases_compile_to_supported_model(alias: str) -> None:
    compiled = compile_model_source(alias, use_cache=False)

    assert compiled.name == "built-in-sm"
    assert compiled.supported
    assert compiled.capabilities["max_vertex_valence"] == 3
    assert compiled.capabilities["parameter_count"] == len(compiled.parameter_defaults)
    particles = {str(particle["name"]): particle for particle in compiled.model["particles"]}
    assert particles["g"]["spin"] == 3
    assert particles["g"]["antiname"] == "g"
    assert particles["d"]["spin"] == 2
    assert particles["d"]["antiname"] == "d~"
    assert particles["h"]["spin"] == 1


@pytest.mark.parametrize("name", ["sm", "scalars", "scalar_gravity"])
@pytest.mark.parametrize("model_format", ["ufo", "json"])
def test_bundled_models_load_and_pass_preflight(
    name: str,
    model_format: str,
) -> None:
    path = bundled_model_path(name, model_format)

    source_kind, detected_path = detect_model_source(path)
    compiled = compile_model_source(path, use_cache=False)

    assert source_kind == model_format
    assert detected_path == path.resolve()
    assert compiled.name == name
    assert compiled.supported
    assert compiled.capabilities["particle_count"] > 0
    assert compiled.capabilities["vertex_count"] > 0
    assert compiled.parameter_defaults
    if name in _BUNDLED_MODEL_STRUCTURE_SHA256:
        assert (
            _compiled_ir_structure_sha256(compiled)
            == _BUNDLED_MODEL_STRUCTURE_SHA256[name]
        )


def test_compiled_model_round_trip_and_parameter_card(tmp_path: Path) -> None:
    compiled = compile_model_source("BUILTIN_SM", use_cache=False)

    model_path = compiled.write(tmp_path / "standard-model")
    parameter_path = compiled.write_parameter_card(tmp_path / "parameters.json")

    assert model_path.name == "standard-model.pyAmplicol-model.json"
    assert load_compiled_model(model_path) == compiled
    payload = json.loads(parameter_path.read_text(encoding="utf-8"))
    assert payload["alpha_s"] == [0.118, 0.0]


def test_compiled_builtin_uses_symbolica_safe_runtime_parameters() -> None:
    compiled = compile_model_source("BUILTIN_SM", use_cache=False)
    model = CompiledUFOModel(compiled)

    assert all("-" not in name for name in compiled.parameter_defaults)
    assert model.runtime_parameter_names_for_particle(-21) == ()
    assert model.runtime_derived_parameter_domains() == {}


def test_compile_model_source_parses_compiled_json_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    compiled = compile_model_source("BUILTIN_SM", use_cache=False)
    model_path = compiled.write(tmp_path / "standard-model")
    original_read_text = Path.read_text
    read_count = 0

    def counted_read_text(path: Path, *args, **kwargs):
        nonlocal read_count
        if path.resolve() == model_path.resolve():
            read_count += 1
        return original_read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", counted_read_text)

    assert compile_model_source(model_path) == compiled
    assert read_count == 1


def test_loaded_compiled_model_write_preserves_exact_serialized_bytes(
    tmp_path: Path,
) -> None:
    compiled = compile_model_source("BUILTIN_SM", use_cache=False)
    source = compiled.write(tmp_path / "source")
    source.write_bytes(source.read_bytes() + b"\n")

    loaded = compile_model_source(source)
    copied = loaded.write(tmp_path / "copied")

    assert copied.read_bytes() == source.read_bytes()


def test_compiled_model_rejects_dependency_fingerprint_mismatch(
    tmp_path: Path,
) -> None:
    compiled = compile_model_source("BUILTIN_SM", use_cache=False)
    path = compiled.write(tmp_path / "standard-model")
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["producer"]["symbolica"] = "different"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="fingerprint mismatch"):
        load_compiled_model(path)


def test_compiled_model_rejects_restriction_or_simplification_options(
    tmp_path: Path,
) -> None:
    path = compile_model_source("BUILTIN_SM", use_cache=False).write(
        tmp_path / "standard-model"
    )

    with pytest.raises(ValueError, match="already compiled"):
        compile_model_source(path, restriction="none")
    with pytest.raises(ValueError, match="already compiled"):
        compile_model_source(path, simplify=False)


def test_preflight_reports_all_incompatible_features() -> None:
    issues, capabilities = preflight_model(
        {
            "particles": [
                {
                    "name": "chi",
                    "antiname": "chi",
                    "spin": 2,
                    "color": 1,
                    "mass": "MCHI",
                },
                {
                    "name": "x",
                    "antiname": "x~",
                    "spin": 4,
                    "color": 10,
                    "mass": "ZERO",
                },
            ],
            "parameters": [],
            "couplings": [{"name": "GC", "expression": "mystery(x)"}],
            "propagators": [],
            "lorentz_structures": [],
            "vertex_rules": [],
            "functions": [{"name": "callback", "arguments": ["x"]}],
            "form_factors": [{"name": "FF"}],
        }
    )

    codes = {issue.code for issue in issues}
    assert {
        "unsupported-spin",
        "unsupported-color-representation",
        "majorana-fermion",
        "form-factors",
        "unknown-functions",
    } <= codes
    assert capabilities["supported"] is False


def test_preflight_rejects_wrong_standard_function_arity() -> None:
    issues, _capabilities = preflight_model(
        {
            "particles": [],
            "parameters": [],
            "couplings": [],
            "propagators": [],
            "lorentz_structures": [],
            "vertex_rules": [],
            "functions": [{"name": "sec", "arguments": ["x", "y"]}],
            "form_factors": [],
        }
    )

    assert any(
        issue.code == "function-arity" and issue.context == "sec/2"
        for issue in issues
    )


def test_detect_model_source_recognizes_compiled_json(tmp_path: Path) -> None:
    path = tmp_path / "model.json"
    path.write_text(json.dumps({"kind": COMPILED_MODEL_KIND}), encoding="utf-8")

    assert detect_model_source(path) == ("pyamplicol", path.resolve())


def test_cache_is_content_addressed(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    first = compile_model_source("BUILTIN_SM", cache_dir=cache)
    cache_files = tuple(cache.glob("*.pyAmplicol-model.json"))

    second = compile_model_source("BUILTIN_SM", cache_dir=cache)

    assert second == first
    assert len(cache_files) == 1
    assert tuple(cache.glob("*.pyAmplicol-model.json")) == cache_files
