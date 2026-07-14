from __future__ import annotations

import json
import math
from dataclasses import replace

import pytest
from symbolica import E, S

from pyamplicol.generic_artifact import (
    GenericProcessManifest,
    _evaluate_current_warmup,
    _generic_warmup_phase_space_point,
    _generic_runtime_schema_payload,
)
from pyamplicol.generic_dag import (
    compile_generic_dag,
    infer_minimal_coupling_order_limits,
    prune_dag_to_amplitude_roots,
)
import pyamplicol.generic_stage_compiler as generic_stage_compiler
import pyamplicol.ufo_model as ufo_model_module

from pyamplicol.generic_stage_compiler import (
    build_generic_stage_compiler_blueprint,
    write_model_parameter_evaluator_artifact,
)
from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.model_assets import bundled_model_path
from pyamplicol.model_processes import build_model_process_ir
from pyamplicol.model_source import CompiledModel, compile_model_source
from pyamplicol.symbolica_evaluator import SymbolicaEvaluatorSettings
from pyamplicol.ufo_ir import (
    _canonicalize_oriented_kernel_component,
    _four_point_contact_color_split,
    _minkowski_dot,
    _transverse_yang_mills_three_vector_components,
    compile_ufo_model_ir,
    eager_color_singlet_vertex_term_components,
)
from pyamplicol.ufo_model import CompiledUFOModel


def _scalar_stage_fixture():
    compiled = compile_model_source(
        bundled_model_path("scalars", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir(
        "scalar_0 scalar_0 > scalar_0",
        compiled.ir,
    )
    dag = compile_generic_dag(process, model=model)
    manifest = GenericProcessManifest(
        dag=dag,
        model=model,
        color_plan=dag.color_plan,
    )
    schema = _generic_runtime_schema_payload(dag, model)
    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        runtime_schema=schema,
        stage_local_parameter_layout=True,
    )
    return model, schema, blueprint


def test_runtime_parameter_clone_reuses_structure_and_overlays_values() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled, {"WZ": 1.25})
    derived_name = next(iter(model._runtime_derived_definitions))
    model.runtime_derived_parameter_defaults_for((derived_name,))

    clone = model.with_runtime_parameters({"WT": 2.5})

    assert clone is not model
    assert clone.compiled is model.compiled
    assert clone._particle_records_by_name is model._particle_records_by_name
    assert clone._parameter_records is model._parameter_records
    assert clone._kernels is model._kernels
    assert clone.particles is model.particles
    assert clone.vertices is model.vertices
    assert clone._runtime_derived_definitions is model._runtime_derived_definitions
    assert (
        clone._runtime_derived_expression_cache
        is model._runtime_derived_expression_cache
    )
    assert clone._runtime_derived_domain_cache is model._runtime_derived_domain_cache
    assert clone._runtime_parameter_domain_cache is model._runtime_parameter_domain_cache
    assert (
        clone._kernel_component_expression_cache
        is model._kernel_component_expression_cache
    )
    assert (
        clone._kernel_coupling_expression_cache
        is model._kernel_coupling_expression_cache
    )
    assert clone._expression_symbol_cache is model._expression_symbol_cache
    assert clone._custom_propagator_expressions is model._custom_propagator_expressions
    assert clone._custom_propagator_templates is model._custom_propagator_templates
    assert clone._color_projection_cache is model._color_projection_cache

    assert clone.width(23) == pytest.approx(1.25)
    assert clone.width(6) == pytest.approx(2.5)
    assert model.width(6) != pytest.approx(2.5)
    assert model._runtime_derived_default_cache
    assert clone._runtime_derived_default_cache == {}


def test_ufo_symbol_substitution_is_simultaneous() -> None:
    x, y, absent = S("x", "y", "absent")

    result = ufo_model_module._replace_symbols(
        x + y,
        {
            x: y + 1,
            y: x + 1,
            absent: 99,
        },
    )

    assert result == x + y + 2


def test_real_external_coupling_feeds_derived_runtime_stage_inputs() -> None:
    _model, schema, blueprint = _scalar_stage_fixture()
    records = {
        str(record["name"]): record
        for record in schema["model_parameters"]
    }

    assert records["lam"]["kind"] == "external_parameter"
    assert records["lam"]["parameter_type"] == "real"
    assert not any(name.startswith("coupling.") for name in records)

    stage = blueprint.stages[0]
    parameter_indices = {
        int(record["parameter_index"]): name
        for name, record in records.items()
    }
    stage_parameters = {
        parameter_indices[component.source_id]: stage.parameter_symbols[
            component.parameter_index
        ]
        for component in stage.input_components
        if component.kind == "model_parameter"
    }
    output_symbols = set(stage.output_expressions[0].get_all_symbols(False))

    assert "lam" not in stage_parameters
    derived_symbols = {
        symbol
        for name, symbol in stage_parameters.items()
        if name.startswith("derived_coupling_")
    }
    assert derived_symbols & output_symbols
    assert stage.model_parameter_count == 1 + len(derived_symbols)

    compiled = compile_model_source(
        bundled_model_path("scalars", "json"),
        use_cache=False,
    )
    definitions = CompiledUFOModel(compiled).runtime_derived_parameter_definitions()
    assert any("UFO::{}::lam" in expression for expression in definitions.values())


def test_ufo_coupling_phase_omits_structurally_zero_stage_components() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir("d d~ > Z g", compiled.ir)
    dag = compile_generic_dag(process, model=model)
    manifest = GenericProcessManifest(
        dag=dag,
        model=model,
        color_plan=dag.color_plan,
    )
    schema = _generic_runtime_schema_payload(dag, model)
    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        runtime_schema=schema,
        stage_local_parameter_layout=True,
    )
    records = {
        int(record["parameter_index"]): record
        for record in schema["model_parameters"]
    }
    domains = model.runtime_derived_parameter_domains()

    assert "imaginary" in domains.values()
    used_derived_records = [
        records[component.source_id]
        for stage in blueprint.stages
        for component in stage.input_components
        if component.kind == "model_parameter"
        and records[component.source_id].get("kind")
        == "derived_parameter_component"
    ]
    imaginary_records = [
        record
        for record in used_derived_records
        if record.get("complex_domain") == "imaginary"
    ]

    assert imaginary_records
    assert all(
        record.get("complex_component") == "imag" for record in imaginary_records
    )


def test_resolved_ufo_w_coupling_uses_only_its_imaginary_stage_slot() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir("u d~ > W+", compiled.ir)
    dag = compile_generic_dag(process, model=model)
    manifest = GenericProcessManifest(
        dag=dag,
        model=model,
        color_plan=dag.color_plan,
    )
    schema = _generic_runtime_schema_payload(dag, model)
    blueprint = build_generic_stage_compiler_blueprint(
        manifest,
        runtime_schema=schema,
        stage_local_parameter_layout=True,
    )
    records = {
        int(record["parameter_index"]): record
        for record in schema["model_parameters"]
    }
    used_derived_records = [
        records[component.source_id]
        for stage in blueprint.stages
        for component in stage.input_components
        if component.kind == "model_parameter"
        and records[component.source_id].get("kind")
        == "derived_parameter_component"
    ]

    assert used_derived_records
    assert all(
        record.get("complex_domain") == "imaginary"
        and record.get("complex_component") == "imag"
        for record in used_derived_records
    )


def test_external_particle_mass_parameters_use_model_owned_names() -> None:
    _model, schema, _blueprint = _scalar_stage_fixture()
    names = {str(record["name"]) for record in schema["model_parameters"]}

    assert "mass_scalar_1" in names
    assert not any(name.startswith("particle.") for name in names)


def test_ufo_warmup_resolves_derived_coupling_defaults(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir("d d~ > Z g", compiled.ir)
    dag = compile_generic_dag(process, model=model)
    points = tuple(
        _generic_warmup_phase_space_point(dag, model, seed=12345 + offset)
        for offset in range(2)
    )
    replace_symbols = ufo_model_module._replace_symbols

    def reject_symbolic_kernel_substitution(expression, substitutions):
        if "kernel_" in expression.to_canonical_string():
            pytest.fail("numeric UFO warmup should use cached Expression.evaluate()")
        return replace_symbols(expression, substitutions)

    monkeypatch.setattr(
        ufo_model_module,
        "_replace_symbols",
        reject_symbolic_kernel_substitution,
    )
    monkeypatch.setattr(
        model,
        "runtime_derived_parameter_defaults",
        lambda: pytest.fail("warmup should request only used derived parameters"),
    )

    maxima, signatures = _evaluate_current_warmup(dag, model, points)

    assert max(maxima.values()) > 0.0
    assert any(signatures[current_id] for current_id in signatures)
    assert model._kernel_component_expression_cache
    assert model._kernel_coupling_expression_cache


def test_ufo_minimal_order_policy_uses_declared_order_hierarchies() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir(
        "d d~ > t t~ z h",
        compiled.ir,
        color_accuracy="lc",
    )
    fixed_helicities = {1: -1, 2: 1, 3: -1, 4: 1, 5: 0, 6: 0}

    limits = infer_minimal_coupling_order_limits(
        process,
        model=model,
        color_accuracy="lc",
    )
    external = prune_dag_to_amplitude_roots(
        compile_generic_dag(
            process,
            model=model,
            max_coupling_orders=limits,
            selected_source_helicities=fixed_helicities,
        )
    )
    built_in = prune_dag_to_amplitude_roots(
        compile_generic_dag(
            "d d~ > t t~ z h",
            model=AmplicolSMLeadingColorModel(),
            color_accuracy="lc",
            selected_source_helicities=fixed_helicities,
        )
    )

    assert limits == {"QCD": 2, "QED": 2}
    assert (
        len(external.currents),
        len(external.interactions),
        len(external.amplitude_roots),
    ) == (
        len(built_in.currents),
        len(built_in.interactions),
        len(built_in.amplitude_roots),
    )


def test_ufo_color_tensor_normalization_is_accuracy_independent() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    generator_vertex = next(
        vertex
        for vertex in model.vertices
        if model.vertex_color_structure(vertex) == "fundamental-generator"
    )
    structure_constant_vertex = next(
        vertex
        for vertex in model.vertices
        if model.vertex_color_structure(vertex) == "adjoint-structure-constant"
    )

    for accuracy in ("lc", "nlc", "full"):
        assert model.vertex_color_weight(
            generator_vertex,
            color_accuracy=accuracy,
        ) == pytest.approx((1.0 / math.sqrt(2.0), 0.0))
        assert model.vertex_color_weight(
            structure_constant_vertex,
            color_accuracy=accuracy,
        ) == pytest.approx((0.0, -1.0 / math.sqrt(2.0)))


def test_compiled_model_persists_oriented_kernel_color_projections(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    round_tripped = CompiledModel.from_dict(compiled.to_dict())

    assert round_tripped.ir.oriented_kernels
    assert all(
        kernel.color_projection_structure is not None
        and kernel.color_projection_coefficient is not None
        for kernel in round_tripped.ir.oriented_kernels
    )

    def reject_reprojection(*_args, **_kwargs):
        raise AssertionError("persisted color projections must not be recomputed")

    monkeypatch.setattr(
        ufo_model_module,
        "classify_trilinear_color_expression",
        reject_reprojection,
    )
    model = CompiledUFOModel(round_tripped)
    vertex = next(
        vertex
        for vertex in model.vertices
        if round_tripped.ir.oriented_kernels[vertex.kind].color_projection_structure
        == "fundamental-generator"
    )

    assert model.vertex_color_structure(vertex) == "fundamental-generator"


def test_custom_ufo_scalar_propagator_is_lowered_and_runtime_parameterized(
    tmp_path,
) -> None:
    payload = json.loads(
        bundled_model_path("scalars", "json").read_text(encoding="utf-8")
    )
    particle = next(item for item in payload["particles"] if item["name"] == "scalar_1")
    propagator = next(
        item for item in payload["propagators"] if item["particle"] == "scalar_1"
    )
    propagator["name"] = "custom_scalar_1_propagator"
    propagator["numerator"] = "2𝑖"
    particle["propagator"] = propagator["name"]
    model_path = tmp_path / "custom-scalars.json"
    model_path.write_text(json.dumps(payload), encoding="utf-8")

    compiled = compile_model_source(
        model_path,
        restriction="none",
        simplify=False,
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    particle_id = int(particle["pdg_code"])
    rule = model.propagator_lowering_rule(particle_id)

    assert rule.backend == "spenso-ufo-custom"
    assert rule.full_tensor_network_ready is True
    assert model.runtime_parameter_names_for_particle(particle_id) == (
        "mass_scalar_1",
    )
    value = model.propagator_component_expression(
        particle_id,
        (3.0,),
        (5.0, 0.0, 0.0, 0.0),
    )[0]
    mass = float(model.mass(particle_id))
    assert complex(value) == pytest.approx(6j / (25.0 - mass * mass))


def test_custom_ufo_vector_and_fermion_propagators_match_default_tensors(
    tmp_path,
) -> None:
    payload = json.loads(
        bundled_model_path("sm", "json").read_text(encoding="utf-8")
    )
    selected = {"Z", "g", "d", "d~", "t", "t~"}
    for particle in payload["particles"]:
        if particle["name"] not in selected:
            continue
        propagator = next(
            item
            for item in payload["propagators"]
            if item["particle"] == particle["name"]
        )
        propagator["name"] += "_custom"
        particle["propagator"] = propagator["name"]
    model_path = tmp_path / "custom-sm.json"
    model_path.write_text(json.dumps(payload), encoding="utf-8")

    custom = CompiledUFOModel(
        compile_model_source(
            model_path,
            restriction="none",
            simplify=False,
            use_cache=False,
        )
    ).with_runtime_parameters({"WZ": 0.0, "WT": 0.0})
    reference = CompiledUFOModel(
        compile_model_source(
            bundled_model_path("sm", "json"),
            use_cache=False,
        )
    ).with_runtime_parameters({"WZ": 0.0, "WT": 0.0})
    momentum = (500.0, 10.0, 20.0, 30.0)

    for name in sorted(selected):
        custom_particle = next(
            particle for particle in custom.compiled.ir.particles if particle.name == name
        )
        reference_particle = next(
            particle
            for particle in reference.compiled.ir.particles
            if particle.name == name
        )
        dimension = custom.dimension(custom_particle.pdg_code)
        current = tuple(complex(index + 1, 0.25 * index) for index in range(dimension))
        actual = custom.propagator_component_expression(
            custom_particle.pdg_code,
            current,
            momentum,
        )
        expected = reference.propagator_component_expression(
            reference_particle.pdg_code,
            current,
            momentum,
        )

        assert custom.propagator_lowering_rule(
            custom_particle.pdg_code
        ).full_tensor_network_ready is True
        assert tuple(complex(value) for value in actual) == pytest.approx(
            tuple(complex(value) for value in expected),
            rel=1.0e-13,
            abs=1.0e-15,
        )
        if custom.is_fermion(custom_particle.pdg_code):
            assert custom.is_chiral_eligible(custom_particle.pdg_code) is False


def test_custom_ufo_spin2_propagator_matches_de_donder_projector() -> None:
    source = json.loads(
        bundled_model_path("scalar_gravity", "json").read_text(encoding="utf-8")
    )
    source_particle = next(
        particle for particle in source["particles"] if particle["name"] == "graviton"
    )
    source_propagator = next(
        propagator
        for propagator in source["propagators"]
        if propagator["particle"] == "graviton"
    )
    propagator_name = "custom_graviton_propagator"
    payload = {
        "name": "custom-spin2-test",
        "orders": [],
        "parameters": [
            {
                "name": "dim",
                "nature": "external",
                "parameter_type": "real",
                "value": [4.0, 0.0],
                "expression": None,
                "lhablock": "DIM",
                "lhacode": [1],
            }
        ],
        "particles": [
            {
                **source_particle,
                "propagator": propagator_name,
            }
        ],
        "couplings": [],
        "propagators": [
            {
                **source_propagator,
                "name": propagator_name,
            }
        ],
        "lorentz_structures": [],
        "vertex_rules": [],
    }
    ir = compile_ufo_model_ir(payload)
    compiled = CompiledModel(
        source={"kind": "test"},
        producer={},
        model=payload,
        ir=ir,
        parameter_defaults={"dim": (4.0, 0.0)},
        capabilities={},
        issues=(),
        phase_timings={},
        conversion_seconds=0.0,
    )
    custom = CompiledUFOModel(compiled)
    reference = CompiledUFOModel(
        replace(
            compiled,
            ir=replace(
                ir,
                propagators=tuple(
                    replace(propagator, custom=False)
                    for propagator in ir.propagators
                ),
            ),
        )
    )
    particle_id = int(source_particle["pdg_code"])
    current = tuple(complex(index + 1, 0.125 * index) for index in range(16))
    momentum = (25.0, 2.0, 3.0, 4.0)

    actual = custom.propagator_component_expression(
        particle_id,
        current,
        momentum,
    )
    expected = reference.propagator_component_expression(
        particle_id,
        current,
        momentum,
    )

    assert custom.runtime_parameter_names_for_particle(particle_id) == ("dim",)
    assert tuple(complex(value) for value in actual) == pytest.approx(
        tuple(complex(value) for value in expected),
        rel=1.0e-13,
        abs=1.0e-14,
    )


def test_momentum_dependent_three_gluon_kernel_builds_pure_gluon_dag() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir("g g > g g", compiled.ir)

    dag = compile_generic_dag(process, model=model)

    assert len(dag.interactions) > 0
    assert len(dag.amplitude_roots) > 0
    assert any(
        model._kernel(interaction.vertex_kind).particles == ("g", "g", "g")
        for interaction in dag.interactions
    )


def test_four_gluon_contact_uses_nonpropagating_model_owned_auxiliaries() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir("g g > g g", compiled.ir)
    dag = compile_generic_dag(
        process,
        model=model,
        selected_color_sector_ids={0},
        max_coupling_orders={"QED": 0},
    )

    auxiliaries = [
        particle
        for particle in compiled.ir.particles
        if particle.auxiliary_kind
        and particle.auxiliary_kind.startswith("ufo-contact:")
    ]
    assert auxiliaries
    assert all(not particle.propagating for particle in auxiliaries)
    contact_interactions = [
        interaction
        for interaction in dag.interactions
        if "contact" in model._kernel(interaction.vertex_kind).vertex
    ]
    assert contact_interactions
    contact_current_ids = {
        interaction.result_id
        for interaction in contact_interactions
        if dag.currents[interaction.result_id].index.auxiliary_kind is not None
    }
    assert contact_current_ids
    for current_id in contact_current_ids:
        current = dag.currents[current_id]
        rule = model.propagator_lowering_rule(current.index.particle_id)
        assert rule.applies_propagator is False
        assert rule.kernel == "ufo_contact_auxiliary_no_propagator"


def test_four_gluon_contact_split_uses_normalized_color_tensor() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    term = next(
        term
        for term in compiled.ir.vertex_terms
        if term.particles == ("g", "g", "g", "g")
    )

    expected = _four_point_contact_color_split(term, 0)
    rewritten_source = replace(term, color_source="writer_specific_color_layout")

    assert _four_point_contact_color_split(rewritten_source, 0) == expected
    assert expected[4:6] == (1, 1)


def test_scalar_contacts_use_complete_balanced_provenance_trees() -> None:
    compiled = compile_model_source(
        bundled_model_path("scalars", "json"),
        use_cache=False,
    )
    contact_term_ids = {
        term.id for term in compiled.ir.vertex_terms if term.valence > 3
    }
    lowered_term_ids = {
        term_id
        for kernel in compiled.ir.oriented_kernels
        if "::contact-" in kernel.vertex and "final" in kernel.vertex
        for term_id in (kernel.term_ids or (kernel.term_id,))
    }
    assert lowered_term_ids == contact_term_ids
    assert compiled.capabilities["compiled_contact_term_count"] == len(contact_term_ids)
    assert compiled.capabilities["unlowered_contact_term_count"] == 0

    ten_point = next(
        term
        for term in compiled.ir.vertex_terms
        if term.valence == 10 and len(set(term.particles)) == 1
    )
    particles = {particle.name: particle for particle in compiled.ir.particles}
    finals = [
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.term_id == ten_point.id
        and kernel.vertex.endswith("::contact-tree-final")
    ]
    assert finals
    root_auxiliaries = [particles[name] for name in finals[0].particles[:2]]
    root_leg_counts = sorted(
        len(str(particle.auxiliary_kind).rsplit(":", 1)[1].split(","))
        for particle in root_auxiliaries
    )
    assert root_leg_counts == [4, 5]

    fragments = [
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.term_id == ten_point.id and "::contact-tree-" in kernel.vertex
    ]
    assert all(
        not kernel.runtime_parameters
        for kernel in fragments
        if kernel.vertex.endswith("::contact-tree-partial")
    )
    assert all(
        kernel.runtime_parameters == (f"derived_coupling_{ten_point.id}",)
        for kernel in fragments
        if kernel.vertex.endswith("::contact-tree-final")
    )
    assert all(
        not particle.propagating and particle.component_dimension == 1
        for particle in particles.values()
        if particle.auxiliary_kind
        and particle.auxiliary_kind.startswith(
            f"ufo-contact-tree:{ten_point.id}:"
        )
    )

    model = CompiledUFOModel(compiled)
    scalar_pdg = next(
        particle.pdg_code
        for particle in compiled.ir.particles
        if particle.name == "scalar_0"
    )
    assert model.vertices_for_inputs(scalar_pdg, scalar_pdg)
    assert model.vertices_for_inputs(
        scalar_pdg,
        scalar_pdg,
        color_accuracy="full",
    ) == model.vertices_for_inputs(scalar_pdg, scalar_pdg)


def test_ten_scalar_contact_has_independent_eager_nary_oracle() -> None:
    compiled = compile_model_source(
        bundled_model_path("scalars", "json"),
        use_cache=False,
    )
    term = next(
        candidate
        for candidate in compiled.ir.vertex_terms
        if candidate.valence == 10
        and candidate.particles == ("scalar_0",) * 10
    )
    input_legs = set(range(term.valence)) - {0}

    components = eager_color_singlet_vertex_term_components(
        term,
        compiled.ir.particles,
        result_leg=0,
        input_components={leg: (E("1"),) for leg in input_legs},
        input_momenta={
            leg: (E("1"), E("0"), E("0"), E("0"))
            for leg in input_legs
        },
        coupling=E("1𝑖"),
    )

    assert components == (E("1𝑖"),)


def test_model_parameter_evaluator_maps_external_inputs_to_derived_slots(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model, schema, _blueprint = _scalar_stage_fixture()
    compiled_stages = []

    def fake_compile(stage, artifact_dir, **kwargs):
        compiled_stages.append(stage)
        assert artifact_dir == tmp_path
        settings = kwargs["symbolica_settings"]
        assert settings.compiled_output_chunk_size is None
        assert settings.output_chunk_strategy == "uniform"
        return {
            "kind": "jit-symbolica-evaluator",
            "input_len": stage.parameter_count,
            "output_len": stage.output_length,
            "evaluator_state_path": "generic_model_parameter_derivation.evaluator.bin",
        }

    monkeypatch.setattr(
        generic_stage_compiler,
        "_compile_stage_evaluator_artifact",
        fake_compile,
    )

    payload = write_model_parameter_evaluator_artifact(
        model,
        schema,
        tmp_path,
        symbolica_settings=SymbolicaEvaluatorSettings(
            compiled_output_chunk_size=128,
            output_chunk_strategy="auto",
        ),
    )

    assert payload is not None
    records = sorted(
        schema["model_parameters"],
        key=lambda record: int(record["parameter_index"]),
    )
    external_indices = [
        int(record["parameter_index"])
        for record in records
        if record["kind"] in {"external_parameter", "external_parameter_component"}
    ]
    assert payload["input_parameter_indices"] == external_indices

    derived_records = {
        (str(record["runtime_name"]), str(record["complex_component"])): int(
            record["parameter_index"]
        )
        for record in records
        if record["kind"] == "derived_parameter_component"
    }
    outputs = payload["outputs"]
    assert outputs
    for output_index, output in enumerate(outputs):
        runtime_name = str(output["runtime_name"])
        assert output["output_index"] == output_index
        assert output["real_parameter_index"] == derived_records[
            (runtime_name, "real")
        ]
        assert output["imag_parameter_index"] == derived_records[
            (runtime_name, "imag")
        ]

    assert len(compiled_stages) == 1
    stage = compiled_stages[0]
    assert stage.stage_kind == "model-parameter-derivation"
    assert stage.parameter_count == len(external_indices)
    assert stage.output_length == len(outputs)
    assert len(stage.output_expressions) == len(outputs)


def test_sm_color_tensors_record_lc_basis_normalization() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    terms = {term.id: term for term in compiled.ir.vertex_terms}

    quark_gluon = next(
        term for term in terms.values() if term.vertex == "V_74"
    )
    three_gluon = next(
        term
        for term in terms.values()
        if term.particles == ("g", "g", "g")
    )
    four_gluon = next(
        term
        for term in terms.values()
        if term.particles == ("g", "g", "g", "g")
    )
    neutral_current = next(
        term for term in terms.values() if term.vertex == "V_79"
    )

    assert quark_gluon.lc_color_normalization_power == 1
    assert three_gluon.lc_color_normalization_power == 1
    assert four_gluon.lc_color_normalization_power == 2
    assert neutral_current.lc_color_normalization_power == 0

    model = CompiledUFOModel(compiled)
    kernel = next(
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.term_id == quark_gluon.id
    )
    vertex = next(vertex for vertex in model.vertices if vertex.kind == kernel.kind)
    assert kernel.color_source == "UFO::{}::T(3,2,1)"
    assert model.vertex_color_structure(vertex) == "fundamental-generator"
    assert model.vertex_color_weight(vertex, color_accuracy="lc") == pytest.approx(
        (2.0**-0.5, 0.0)
    )
    assert model.vertex_color_weight(vertex, color_accuracy="full") == pytest.approx(
        (2.0**-0.5, 0.0)
    )

    three_gluon_kernel = next(
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.term_id == three_gluon.id
    )
    three_gluon_vertex = next(
        vertex
        for vertex in model.vertices
        if vertex.kind == three_gluon_kernel.kind
    )
    assert model.vertex_color_structure(three_gluon_vertex) == (
        "adjoint-structure-constant"
    )
    assert model.vertex_color_weight(
        three_gluon_vertex,
        color_accuracy="lc",
    ) == pytest.approx((0.0, -(2.0**-0.5)))


def test_lc_fundamental_generator_uses_fierz_singlet_branch_only_on_one_line(
) -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir("d d~ > t t~", compiled.ir)

    same_line = compile_generic_dag(
        process,
        model=model,
        selected_color_sector_ids={0},
        max_coupling_orders={"QED": 0},
    )
    connected = compile_generic_dag(
        process,
        model=model,
        selected_color_sector_ids={1},
        max_coupling_orders={"QED": 0},
    )

    same_line_gluons = [
        current
        for current in same_line.currents
        if current.index.particle_id == 21
        and current.index.external_labels == (1, 2)
    ]
    assert same_line_gluons
    assert {
        current.index.color_state.basis_key for current in same_line_gluons
    } == {("lc-fierz-singlet",)}
    same_line_ids = {current.id for current in same_line_gluons}
    same_line_weights = {
        interaction.color_weight
        for interaction in same_line.interactions
        if interaction.result_id in same_line_ids
    }
    assert len(same_line_weights) == 1
    assert next(iter(same_line_weights)) == pytest.approx(
        (1.0 / (3.0 * 2.0**0.5), 0.0)
    )
    assert all(
        "lc-fierz-singlet"
        not in same_line.currents[interaction.result_id].index.color_state.basis_key
        for interaction in same_line.interactions
        if interaction.left_id in same_line_ids or interaction.right_id in same_line_ids
    )

    connected_gluons = [
        current
        for current in connected.currents
        if current.index.particle_id == 21
        and current.index.external_labels == (3, 4)
    ]
    assert connected_gluons
    assert all(not current.index.color_state.basis_key for current in connected_gluons)


def test_shared_lc_color_identity_closure_stays_on_compatible_quark_pairing(
) -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    model = CompiledUFOModel(compiled)
    process = build_model_process_ir(
        "d d~ > t t~ z h",
        compiled.ir,
        color_accuracy="lc",
    )
    same_line = compile_generic_dag(
        process,
        model=model,
        selected_color_sector_ids={0},
    )
    crossed = compile_generic_dag(
        process,
        model=model,
        selected_color_sector_ids={0},
        reference_color_order=(3, 1, 2, 4, 5, 6),
    )

    def root_orders(dag, root) -> tuple[tuple[str, int], ...]:
        totals: dict[str, int] = {}
        for current_id in (root.left_id, root.right_id):
            for name, value in dag.currents[current_id].index.coupling_orders:
                totals[name] = totals.get(name, 0) + value
        return tuple(sorted(totals.items()))

    assert len(same_line.amplitude_roots) == 48
    assert len(crossed.amplitude_roots) == 24
    assert any(
        root_orders(same_line, root) == (("QED", 4),)
        for root in same_line.amplitude_roots
    )
    assert {
        root_orders(crossed, root) for root in crossed.amplitude_roots
    } == {(('QCD', 2), ('QED', 2))}
    assert (len(crossed.currents), len(crossed.interactions)) == (99, 199)


def test_fused_ufo_kernel_factors_repeated_lorentz_blocks() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    kernel = next(
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.particles == ("d", "Z", "d") and len(kernel.term_ids) == 2
    )

    shared_components = kernel.component_expressions[2:]
    assert shared_components
    for component in shared_components:
        for runtime_parameter in kernel.runtime_parameters:
            assert component.count(f"UFO::{{}}::{runtime_parameter}") <= 1


def test_weyl_projection_maps_ufo_dirac_blocks_to_physical_chiralities() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    left_coupling = S("test::left_coupling")
    right_coupling = S("test::right_coupling")
    model = CompiledUFOModel(
        compiled,
        {
            "derived_coupling_80": left_coupling,
            "derived_coupling_81": right_coupling,
        },
    )
    kernel = next(
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.particles == ("d", "Z", "d") and len(kernel.term_ids) == 2
    )
    fermion = (S("test::f0"), S("test::f1"))
    vector = tuple(S(f"test::v{index}") for index in range(4))

    positive = model.vertex_component_expression(
        kernel.kind,
        fermion,
        vector,
        result_particle_id=1,
        result_chirality=1,
        left_chirality=1,
    )
    negative = model.vertex_component_expression(
        kernel.kind,
        fermion,
        vector,
        result_particle_id=1,
        result_chirality=-1,
        left_chirality=-1,
    )
    positive_symbols = {
        symbol
        for component in positive
        for symbol in component.get_all_symbols(False)
    }
    negative_symbols = {
        symbol
        for component in negative
        for symbol in component.get_all_symbols(False)
    }

    assert right_coupling in positive_symbols
    assert left_coupling not in positive_symbols
    assert left_coupling in negative_symbols
    assert right_coupling in negative_symbols


def test_kernel_canonicalization_cancels_momenta_without_expanding_couplings() -> None:
    expression = E("g*(p+q-(p+q)+x+y)")

    canonical = _canonicalize_oriented_kernel_component(expression)

    assert canonical == E("g*(x+y)")
    assert str(canonical).count("g") == 1


def test_yang_mills_kernel_retains_shared_minkowski_contraction() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    kernel = next(
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.vertex == "V_36" and kernel.particles == ("g", "g", "g")
    )
    left = tuple(
        E(f"pyamplicol::kernel_{kernel.kind}_left_{index}")
        for index in range(4)
    )
    right = tuple(
        E(f"pyamplicol::kernel_{kernel.kind}_right_{index}")
        for index in range(4)
    )
    shared_dot = _minkowski_dot(left, right).expand_num().to_canonical_string()

    assert all(shared_dot in component for component in kernel.component_expressions)


def test_massless_yang_mills_kernel_uses_transverse_berends_giele_current() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    kernel = next(
        kernel
        for kernel in compiled.ir.oriented_kernels
        if kernel.vertex == "V_36" and kernel.particles == ("g", "g", "g")
    )
    left = tuple(
        S(f"pyamplicol::kernel_{kernel.kind}_left_{index}")
        for index in range(4)
    )
    right = tuple(
        S(f"pyamplicol::kernel_{kernel.kind}_right_{index}")
        for index in range(4)
    )
    left_momentum = tuple(
        S(f"pyamplicol::kernel_{kernel.kind}_left_momentum_{index}")
        for index in range(4)
    )
    right_momentum = tuple(
        S(f"pyamplicol::kernel_{kernel.kind}_right_momentum_{index}")
        for index in range(4)
    )
    expected = _transverse_yang_mills_three_vector_components(
        left_components=left,
        right_components=right,
        left_momentum=left_momentum,
        right_momentum=right_momentum,
    )
    coupling = S(f"UFO::derived_coupling_{kernel.term_id}")

    assert all(
        (E(actual) - target * coupling).expand() == E("0")
        for actual, target in zip(
            kernel.component_expressions,
            expected,
            strict=True,
        )
    )


def test_sm_goldstone_metadata_and_unitary_gauge_internal_policy() -> None:
    compiled_json = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    compiled_ufo = compile_model_source(
        bundled_model_path("sm", "ufo"),
        use_cache=False,
    )
    expected = {"G0", "G+", "G-"}

    for compiled in (compiled_json, compiled_ufo):
        goldstones = {
            particle.name
            for particle in compiled.ir.particles
            if particle.goldstoneboson
        }
        assert goldstones == expected
        model = CompiledUFOModel(compiled)
        assert model.inactive_goldstone_names == expected
        assert all(
            not (set(kernel.particles) & expected)
            for kernel in compiled.ir.oriented_kernels
            if any(vertex.kind == kernel.kind for vertex in model.vertices)
        )


def test_custom_vector_gauge_keeps_its_matching_goldstone() -> None:
    compiled = compile_model_source(
        bundled_model_path("sm", "json"),
        use_cache=False,
    )
    propagators = tuple(
        replace(propagator, custom=True)
        if propagator.particle == "Z"
        else propagator
        for propagator in compiled.ir.propagators
    )
    custom = replace(compiled, ir=replace(compiled.ir, propagators=propagators))

    model = CompiledUFOModel(custom)

    assert "G0" not in model.inactive_goldstone_names
    assert {"G+", "G-"}.issubset(model.inactive_goldstone_names)
    assert any(
        "G0" in kernel.particles
        for kernel in custom.ir.oriented_kernels
        if any(vertex.kind == kernel.kind for vertex in model.vertices)
    )
