from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Mapping

from .ufo_ir import CompiledModelIR, compile_builtin_model_ir, compile_ufo_model_ir


COMPILED_MODEL_KIND = "pyamplicol-compiled-model"
COMPILED_MODEL_SCHEMA_VERSION = 3
MODEL_COMPILER_VERSION = 2
BUILTIN_SM_ALIASES = frozenset(("builtin_sm", "built-in-sm"))
DEFAULT_MODEL_RESTRICTION = "default"
NO_MODEL_RESTRICTION = "none"

SUPPORTED_FUNCTION_ARITIES = {
    name: frozenset({1})
    for name in (
        "Theta",
        "abs",
        "acos",
        "acosh",
        "acsc",
        "asec",
        "asin",
        "asinh",
        "atan",
        "atanh",
        "complexconjugate",
        "conj",
        "cos",
        "cosh",
        "csc",
        "exp",
        "im",
        "log",
        "log10",
        "re",
        "reglog",
        "reglogm",
        "reglogp",
        "sec",
        "sin",
        "sinh",
        "sqrt",
        "tan",
        "tanh",
    )
}
SUPPORTED_FUNCTION_ARITIES.update(
    {
        "complex": frozenset({2}),
        "cond": frozenset({3}),
        "if": frozenset({3}),
        "pow": frozenset({2}),
    }
)
SUPPORTED_FUNCTIONS = frozenset(SUPPORTED_FUNCTION_ARITIES)

UFO_TENSOR_HEADS = frozenset(
    {
        "C",
        "Epsilon",
        "EpsilonBar",
        "Gamma",
        "Gamma5",
        "Identity",
        "IdentityL",
        "K6",
        "K6Bar",
        "Metric",
        "P",
        "PSlash",
        "ProjM",
        "ProjP",
        "Sigma",
        "T",
        "T6",
        "d",
        "dummy",
        "f",
        "idx",
    }
)


@dataclass(frozen=True)
class ModelCompatibilityIssue:
    severity: str
    code: str
    message: str
    context: str = ""

    def to_dict(self) -> dict[str, str]:
        return {
            "severity": self.severity,
            "code": self.code,
            "message": self.message,
            "context": self.context,
        }

    @staticmethod
    def from_dict(payload: Mapping[str, object]) -> "ModelCompatibilityIssue":
        return ModelCompatibilityIssue(
            severity=str(payload["severity"]),
            code=str(payload["code"]),
            message=str(payload["message"]),
            context=str(payload.get("context", "")),
        )


@dataclass(frozen=True)
class CompiledModel:
    source: Mapping[str, object]
    producer: Mapping[str, object]
    model: Mapping[str, object]
    ir: CompiledModelIR
    parameter_defaults: Mapping[str, tuple[float, float]]
    capabilities: Mapping[str, object]
    issues: tuple[ModelCompatibilityIssue, ...]
    phase_timings: Mapping[str, float]
    conversion_seconds: float
    _serialized_path: Path | None = field(default=None, compare=False, repr=False)

    @property
    def name(self) -> str:
        return str(self.model.get("name", "unnamed-model"))

    @property
    def supported(self) -> bool:
        return not any(issue.severity == "error" for issue in self.issues)

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": COMPILED_MODEL_KIND,
            "schema_version": COMPILED_MODEL_SCHEMA_VERSION,
            "model_compiler_version": MODEL_COMPILER_VERSION,
            "source": dict(self.source),
            "producer": dict(self.producer),
            "model": dict(self.model),
            "ir": self.ir.to_dict(),
            "parameter_defaults": {
                name: [value[0], value[1]]
                for name, value in sorted(self.parameter_defaults.items())
            },
            "capabilities": dict(self.capabilities),
            "issues": [issue.to_dict() for issue in self.issues],
            "phase_timings": dict(self.phase_timings),
            "conversion_seconds": self.conversion_seconds,
        }

    def write(self, path: Path) -> Path:
        output = _compiled_model_output_path(path)
        if self._serialized_path is None or not self._serialized_path.is_file():
            _atomic_write_json(output, self.to_dict(), compact=True)
        else:
            _atomic_copy(self._serialized_path, output)
        return output

    def write_parameter_card(self, path: Path) -> Path:
        _atomic_write_json(
            path,
            {
                name: [value[0], value[1]]
                for name, value in sorted(self.parameter_defaults.items())
            },
        )
        return path

    @staticmethod
    def from_dict(
        payload: Mapping[str, object],
        *,
        validate_fingerprint: bool = True,
        serialized_path: Path | None = None,
    ) -> "CompiledModel":
        if payload.get("kind") != COMPILED_MODEL_KIND:
            raise ValueError("file is not a pyAmpliCol compiled model")
        if int(payload.get("schema_version", -1)) != COMPILED_MODEL_SCHEMA_VERSION:
            raise ValueError("compiled model schema mismatch; regenerate the model")
        if int(payload.get("model_compiler_version", -1)) != MODEL_COMPILER_VERSION:
            raise ValueError("compiled model compiler mismatch; regenerate the model")
        producer = _mapping(payload.get("producer"), "producer")
        if validate_fingerprint and producer != compiler_fingerprint():
            raise ValueError(
                "compiled model dependency fingerprint mismatch; regenerate the model"
            )
        parameter_payload = _mapping(
            payload.get("parameter_defaults"),
            "parameter_defaults",
        )
        parameters = {
            str(name): _complex_pair(value, context=f"parameter {name}")
            for name, value in parameter_payload.items()
        }
        raw_issues = payload.get("issues", ())
        if not isinstance(raw_issues, list):
            raise ValueError("compiled model issues must be a list")
        return CompiledModel(
            source=_mapping(payload.get("source"), "source"),
            producer=producer,
            model=_mapping(payload.get("model"), "model"),
            ir=CompiledModelIR.from_dict(_mapping(payload.get("ir"), "ir")),
            parameter_defaults=parameters,
            capabilities=_mapping(payload.get("capabilities"), "capabilities"),
            issues=tuple(
                ModelCompatibilityIssue.from_dict(_mapping(issue, "issue"))
                for issue in raw_issues
            ),
            phase_timings={
                str(name): float(value)
                for name, value in _mapping(
                    payload.get("phase_timings"),
                    "phase_timings",
                ).items()
            },
            conversion_seconds=float(payload.get("conversion_seconds", 0.0)),
            _serialized_path=(
                None if serialized_path is None else Path(serialized_path).resolve()
            ),
        )


def compiler_fingerprint() -> dict[str, object]:
    return {
        "pyamplicol": _distribution_version("pyamplicol", "0.0.1"),
        "ufo_model_loader": _distribution_version("ufo-model-loader", "missing"),
        "symbolica": _distribution_version("symbolica", "missing"),
        "model_compiler_version": MODEL_COMPILER_VERSION,
        "model_compiler_sha256": _model_compiler_digest(),
        "managed_revisions": _managed_dependency_revisions(),
    }


def detect_model_source(source: str | Path) -> tuple[str, str | Path]:
    source_kind, resolved, _payload = _detect_model_source_payload(source)
    return source_kind, resolved


def _detect_model_source_payload(
    source: str | Path,
) -> tuple[str, str | Path, dict[str, object] | None]:
    text = str(source)
    if text.lower() in BUILTIN_SM_ALIASES:
        return "built-in-sm", "built-in-sm", None
    path = Path(text).expanduser().resolve()
    if path.is_dir():
        if not (path / "__init__.py").is_file():
            raise ValueError(f"UFO model directory has no __init__.py: {path}")
        return "ufo", path, None
    if not path.is_file():
        raise ValueError(f"model source does not exist: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"model JSON could not be read: {path}: {exc}") from exc
    if isinstance(payload, dict) and payload.get("kind") == COMPILED_MODEL_KIND:
        return "pyamplicol", path, payload
    return "json", path, payload if isinstance(payload, dict) else None


def compile_model_source(
    source: str | Path = "BUILTIN_SM",
    *,
    restriction: str = "default",
    simplify: bool = True,
    cache_dir: Path | None = None,
    use_cache: bool = True,
    require_supported: bool = True,
) -> CompiledModel:
    if not restriction or any(character.isspace() for character in restriction):
        raise ValueError("model restriction must be default, none, or a restriction name")
    source_kind, resolved, detected_payload = _detect_model_source_payload(source)
    if source_kind == "pyamplicol":
        if restriction != DEFAULT_MODEL_RESTRICTION or not simplify:
            raise ValueError(
                "restriction and simplification options cannot be applied to an "
                "already compiled pyAmpliCol model"
            )
        if detected_payload is None:
            raise RuntimeError("compiled model detection did not retain its JSON payload")
        return CompiledModel.from_dict(
            detected_payload,
            serialized_path=Path(resolved),
        )
    source_digest = _source_digest(source_kind, resolved)
    fingerprint = compiler_fingerprint()
    cache_key = hashlib.sha256(
        json.dumps(
            {
                "source_digest": source_digest,
                "restriction": restriction,
                "simplify": simplify,
                "producer": fingerprint,
            },
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()
    active_cache = cache_dir or _default_model_cache_dir()
    cache_path = active_cache / f"{cache_key}.pyAmplicol-model.json"
    if use_cache and cache_path.is_file():
        compiled = load_compiled_model(cache_path)
        _raise_for_unsupported(compiled, require_supported=require_supported)
        return compiled

    started = time.perf_counter()
    phase_timings: dict[str, float] = {}
    phase_started = time.perf_counter()
    if source_kind == "built-in-sm":
        model_payload, parameter_defaults = _built_in_model_payload()
    else:
        worker_payload = _load_with_worker(
            Path(resolved),
            restriction=restriction,
            simplify=simplify,
        )
        model_payload = _mapping(worker_payload.get("model"), "worker model")
        card = _mapping(worker_payload.get("parameter_card"), "parameter card")
        parameter_defaults = {
            str(name): _complex_pair(value, context=f"parameter {name}")
            for name, value in card.items()
        }
    phase_timings["model_loading"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    issues, capabilities = preflight_model(model_payload)
    phase_timings["preflight"] = time.perf_counter() - phase_started
    phase_started = time.perf_counter()
    model_ir = (
        compile_builtin_model_ir(model_payload)
        if source_kind == "built-in-sm"
        else compile_ufo_model_ir(model_payload)
    )
    phase_timings["tensor_lowering"] = time.perf_counter() - phase_started
    contact_term_ids = {
        term.id for term in model_ir.vertex_terms if term.valence > 3
    }
    lowered_contact_term_ids = {
        term_id
        for kernel in model_ir.oriented_kernels
        if "::contact-" in kernel.vertex and "final" in kernel.vertex
        for term_id in (kernel.term_ids or (kernel.term_id,))
    }
    unlowered_contact_term_ids = contact_term_ids - lowered_contact_term_ids
    if unlowered_contact_term_ids:
        issues = (
            *issues,
            ModelCompatibilityIssue(
                "error",
                "unsupported-contact-color-lowering",
                "one or more higher-point color tensors could not be lowered",
                ", ".join(str(term_id) for term_id in sorted(unlowered_contact_term_ids)),
            ),
        )
    capabilities = {
        **capabilities,
        "compiled_vertex_term_count": len(model_ir.vertex_terms),
        "compiled_contact_term_count": len(lowered_contact_term_ids),
        "unlowered_contact_term_count": len(unlowered_contact_term_ids),
        "compiled_propagator_count": len(model_ir.propagators),
    }
    conversion_seconds = time.perf_counter() - started
    phase_timings["total"] = conversion_seconds
    compiled = CompiledModel(
        source={
            "kind": source_kind,
            "path": None if source_kind == "built-in-sm" else str(resolved),
            "digest": source_digest,
            "restriction": restriction,
            "simplify": simplify,
        },
        producer=fingerprint,
        model=model_payload,
        ir=model_ir,
        parameter_defaults=parameter_defaults,
        capabilities=capabilities,
        issues=issues,
        phase_timings=phase_timings,
        conversion_seconds=conversion_seconds,
    )
    if use_cache:
        compiled.write(cache_path)
        compiled = replace(compiled, _serialized_path=cache_path.resolve())
    _raise_for_unsupported(compiled, require_supported=require_supported)
    return compiled


def load_compiled_model(path: Path) -> CompiledModel:
    try:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not load compiled model {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"compiled model root must be an object: {path}")
    return CompiledModel.from_dict(payload, serialized_path=Path(path))


def preflight_model(
    model: Mapping[str, object],
) -> tuple[tuple[ModelCompatibilityIssue, ...], dict[str, object]]:
    particles = _list_of_mappings(model.get("particles"), "particles")
    vertices = _list_of_mappings(model.get("vertex_rules"), "vertex_rules")
    functions = _list_of_mappings(model.get("functions", []), "functions")
    form_factors = _list_of_mappings(model.get("form_factors", []), "form_factors")
    issues: list[ModelCompatibilityIssue] = []
    spins = sorted({int(particle.get("spin", 0)) for particle in particles})
    colors = sorted({int(particle.get("color", 0)) for particle in particles})
    unsupported_spins = sorted(set(spins) - {-1, 1, 2, 3, 5})
    if unsupported_spins:
        issues.append(
            ModelCompatibilityIssue(
                "error",
                "unsupported-spin",
                f"unsupported UFO spin codes: {unsupported_spins}",
            )
        )
    unsupported_colors = sorted(set(colors) - {-6, -3, 1, 3, 6, 8})
    if unsupported_colors:
        issues.append(
            ModelCompatibilityIssue(
                "error",
                "unsupported-color-representation",
                f"unsupported UFO color representations: {unsupported_colors}",
            )
        )
    if any(abs(color) == 6 for color in colors):
        issues.append(
            ModelCompatibilityIssue(
                "warning",
                "sextet-process-generation-disabled",
                "sextet metadata is preserved, but processes involving sextets are disabled",
            )
        )
    majorana = sorted(
        str(particle.get("name"))
        for particle in particles
        if int(particle.get("spin", 0)) == 2
        and particle.get("name") == particle.get("antiname")
    )
    if majorana:
        issues.append(
            ModelCompatibilityIssue(
                "error",
                "majorana-fermion",
                "Majorana/FNV fermion flow is not implemented",
                ", ".join(majorana),
            )
        )
    if form_factors:
        issues.append(
            ModelCompatibilityIssue(
                "error",
                "form-factors",
                "UFO form factors are not supported in this milestone",
                ", ".join(sorted(str(item.get("name")) for item in form_factors)),
            )
        )
    declared_functions = {str(function.get("name")) for function in functions}
    invalid_function_arities = []
    for function in functions:
        name = str(function.get("name"))
        if name not in SUPPORTED_FUNCTION_ARITIES:
            continue
        arity = len(_sequence(function.get("arguments")))
        if arity not in SUPPORTED_FUNCTION_ARITIES[name]:
            invalid_function_arities.append(f"{name}/{arity}")
    if invalid_function_arities:
        issues.append(
            ModelCompatibilityIssue(
                "error",
                "function-arity",
                "model function declarations do not match pyAmpliCol's registry",
                ", ".join(sorted(invalid_function_arities)),
            )
        )
    expression_functions = _model_expression_functions(model)
    unknown_functions = sorted(
        (declared_functions | expression_functions)
        - SUPPORTED_FUNCTIONS
        - UFO_TENSOR_HEADS
    )
    if unknown_functions:
        issues.append(
            ModelCompatibilityIssue(
                "error",
                "unknown-functions",
                "model uses functions that are not in pyAmpliCol's hard-coded registry",
                ", ".join(unknown_functions),
            )
        )
    massive_tensors = [
        particle
        for particle in particles
        if int(particle.get("spin", 0)) == 5
        and str(particle.get("mass", "ZERO")).upper() != "ZERO"
    ]
    if massive_tensors:
        issues.append(
            ModelCompatibilityIssue(
                "warning",
                "experimental-massive-spin-2",
                "massive spin-2 support is experimental",
                ", ".join(str(item.get("name")) for item in massive_tensors),
            )
        )
    max_valence = max((len(_sequence(vertex.get("particles"))) for vertex in vertices), default=0)
    capabilities = {
        "supported": not any(issue.severity == "error" for issue in issues),
        "particle_count": len(particles),
        "parameter_count": len(_sequence(model.get("parameters"))),
        "vertex_count": len(vertices),
        "max_vertex_valence": max_valence,
        "spins": spins,
        "color_representations": colors,
        "declared_functions": sorted(declared_functions),
        "expression_functions": sorted(expression_functions),
        "form_factor_count": len(form_factors),
        "has_custom_propagators": any(
            particle.get("propagator")
            and not str(particle.get("propagator")).endswith("propFeynman")
            for particle in particles
        ),
        "color_accuracy_modes": ["lc", "nlc", "full"],
    }
    return tuple(issues), capabilities


def _load_with_worker(
    source: Path,
    *,
    restriction: str,
    simplify: bool,
) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="pyamplicol-ufo-") as temporary:
        output = Path(temporary) / "model.json"
        command = [
            sys.executable,
            "-m",
            "pyamplicol.ufo_worker",
            "--source",
            str(source),
            "--restriction",
            restriction,
            "--output",
            str(output),
        ]
        command.append("--simplify" if simplify else "--no-simplify")
        env = dict(os.environ)
        package_root = str(Path(__file__).resolve().parents[1])
        env["PYTHONPATH"] = package_root + os.pathsep + env.get("PYTHONPATH", "")
        completed = subprocess.run(
            command,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip() or completed.stdout.strip()
            raise ValueError(f"UFO model conversion failed: {detail}")
        payload = json.loads(output.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("UFO worker produced a non-object payload")
        return payload


def _built_in_model_payload() -> tuple[dict[str, object], dict[str, tuple[float, float]]]:
    from .model import AmplicolSMLeadingColorModel

    model = AmplicolSMLeadingColorModel()
    particles = []
    parameters: dict[str, tuple[float, float]] = {
        "alpha_s": (0.118, 0.0),
        "alpha_ew": (1.0 / 132.507, 0.0),
    }
    for _key, particle in sorted(model.particles.items()):
        mass_name = f"mass_{particle.pdg}"
        width_name = f"width_{particle.pdg}"
        parameters[mass_name] = (float(particle.mass), 0.0)
        parameters[width_name] = (float(particle.width), 0.0)
        for pdg in dict.fromkeys((particle.pdg, particle.anti_pdg)):
            is_antiparticle = pdg != particle.pdg
            anti_pdg = particle.pdg if is_antiparticle else particle.anti_pdg
            color = particle.color_rep
            if is_antiparticle and abs(color) in {3, 6}:
                color = -color
            particles.append(
                {
                    "pdg_code": pdg,
                    "name": _built_in_particle_name(pdg),
                    "antiname": _built_in_particle_name(anti_pdg),
                    "spin": _built_in_ufo_spin(pdg, particle.spin),
                    "color": color,
                    "mass": mass_name,
                    "width": width_name,
                    "texname": _built_in_particle_name(pdg),
                    "antitexname": _built_in_particle_name(anti_pdg),
                    "charge": -particle.charge if is_antiparticle else particle.charge,
                    "ghost_number": 0,
                    "lepton_number": 0,
                    "y_charge": 0,
                    "propagating": particle.spin >= 0,
                    "goldstoneboson": False,
                    "propagator": f"builtin_prop_{pdg}",
                }
            )
    return (
        {
            "name": "built-in-sm",
            "restriction": None,
            "orders": [
                {"name": "QCD", "expansion_order": 99, "hierarchy": 1},
                {"name": "QED", "expansion_order": 99, "hierarchy": 2},
            ],
            "parameters": [
                {
                    "name": name,
                    "nature": "external",
                    "parameter_type": "real",
                    "value": [value[0], value[1]],
                    "expression": None,
                    "lhablock": "PYAMPLICOL",
                    "lhacode": [index],
                }
                for index, (name, value) in enumerate(
                    sorted(parameters.items()),
                    start=1,
                )
            ],
            "particles": particles,
            "propagators": [],
            "lorentz_structures": [],
            "couplings": [],
            "vertex_rules": [
                {
                    "name": f"builtin_vertex_{index}",
                    "particles": [
                        _built_in_particle_name(pdg) for pdg in vertex.particles
                    ],
                    "color_structures": ["1"],
                    "lorentz_structures": [f"builtin_kind_{vertex.kind}"],
                    "couplings": [[f"builtin_coupling_{index}"]],
                    "builtin_kind": vertex.kind,
                    "builtin_coupling": list(vertex.coupling),
                }
                for index, vertex in enumerate(model.vertices)
            ],
            "functions": [],
            "form_factors": [],
            "builtin_model": True,
        },
        parameters,
    )


def _built_in_particle_name(pdg: int) -> str:
    from .processes import PDGS

    candidates = [name for name, value in PDGS.items() if int(value) == pdg]
    if candidates:
        return min(candidates, key=lambda name: (len(name), name))
    return f"pdg_{pdg}"


def _built_in_ufo_spin(pdg: int, internal_spin: int) -> int:
    """Translate the legacy kernel spin tag to the UFO 2S+1 convention."""
    if internal_spin < 0:
        return -1
    absolute_pdg = abs(pdg)
    if 1 <= absolute_pdg <= 6 or 11 <= absolute_pdg <= 16:
        return 2
    if absolute_pdg in {21, 22, 23, 24}:
        return 3
    if absolute_pdg == 25:
        return 1
    raise ValueError(
        f"cannot map built-in particle {pdg} with spin tag {internal_spin} to UFO spin"
    )


def _model_expression_functions(model: Mapping[str, object]) -> set[str]:
    expressions: list[str] = []
    for key in ("parameters", "couplings", "propagators", "lorentz_structures"):
        for item in _list_of_mappings(model.get(key, []), key):
            for value in item.values():
                if isinstance(value, str):
                    expressions.append(value)
    for vertex in _list_of_mappings(model.get("vertex_rules", []), "vertex_rules"):
        expressions.extend(
            str(value) for value in _sequence(vertex.get("color_structures"))
        )
    heads = set()
    for expression in expressions:
        for match in re.finditer(
            r"(?:(?:UFO|spenso)::(?:\{\}::)?)?([A-Za-z][A-Za-z0-9_]*)\(",
            expression,
        ):
            heads.add(match.group(1))
    return heads


def _source_digest(kind: str, source: str | Path) -> str:
    digest = hashlib.sha256()
    digest.update(kind.encode("utf-8"))
    if kind == "built-in-sm":
        for name in ("model.py", "processes.py"):
            digest.update(Path(__file__).with_name(name).read_bytes())
        return digest.hexdigest()
    path = Path(source)
    files = [path] if path.is_file() else sorted(
        item
        for item in path.rglob("*")
        if item.is_file() and item.suffix != ".pyc" and "__pycache__" not in item.parts
    )
    for item in files:
        digest.update(str(item.relative_to(path.parent)).encode("utf-8"))
        digest.update(item.read_bytes())
    return digest.hexdigest()


def _managed_dependency_revisions() -> dict[str, object]:
    root = Path(__file__).resolve().parents[2]
    manifest = root / "dependencies" / "install_manifest.json"
    if not manifest.is_file():
        return {}
    try:
        payload = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    revisions: dict[str, object] = {}
    if isinstance(payload, dict):
        for name in ("symbolica", "symbolica_community", "ufo_model_loader"):
            entry = payload.get(name)
            if isinstance(entry, dict):
                revisions[name] = entry.get("source_rev")
    return revisions


def _model_compiler_digest() -> str:
    digest = hashlib.sha256()
    for name in ("model_source.py", "ufo_worker.py", "ufo_ir.py", "ufo_tensors.py"):
        path = Path(__file__).with_name(name)
        digest.update(name.encode("utf-8"))
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _default_model_cache_dir() -> Path:
    configured = os.environ.get("PYAMPLICOL_CACHE_DIR")
    root = Path(configured).expanduser() if configured else Path.home() / ".cache" / "pyamplicol"
    return root / "models"


def _compiled_model_output_path(path: Path) -> Path:
    text = str(path)
    suffix = ".pyAmplicol-model.json"
    return Path(text if text.endswith(suffix) else text + suffix)


def _atomic_write_json(
    path: Path,
    payload: Mapping[str, object],
    *,
    compact: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    serialized = (
        json.dumps(payload, separators=(",", ":"), sort_keys=True)
        if compact
        else json.dumps(payload, indent=2, sort_keys=True)
    )
    temporary.write_text(serialized + "\n", encoding="utf-8")
    temporary.replace(path)


def _atomic_copy(source: Path, path: Path) -> None:
    source = source.resolve()
    path = path.resolve()
    if source == path:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    shutil.copyfile(source, temporary)
    temporary.replace(path)


def _distribution_version(name: str, default: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return default


def _raise_for_unsupported(
    model: CompiledModel,
    *,
    require_supported: bool,
) -> None:
    if require_supported and not model.supported:
        details = "; ".join(
            f"{issue.code}: {issue.message}"
            + (f" ({issue.context})" if issue.context else "")
            for issue in model.issues
            if issue.severity == "error"
        )
        raise ValueError(f"model {model.name!r} is not supported: {details}")


def _mapping(value: object, context: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object")
    return {str(key): item for key, item in value.items()}


def _sequence(value: object) -> list[object]:
    if isinstance(value, (list, tuple)):
        return list(value)
    return []


def _list_of_mappings(value: object, context: str) -> list[dict[str, object]]:
    return [_mapping(item, context) for item in _sequence(value)]


def _complex_pair(value: object, *, context: str) -> tuple[float, float]:
    pair = _sequence(value)
    if len(pair) != 2:
        raise ValueError(f"{context} must be [real, imaginary]")
    return float(pair[0]), float(pair[1])
