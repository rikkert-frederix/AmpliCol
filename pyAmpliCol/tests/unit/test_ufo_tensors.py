from __future__ import annotations

import pytest
from symbolica import E, S
from symbolica.community.spenso import (
    LibraryTensor,
    Representation,
    TensorLibrary,
    TensorName,
    TensorNetwork,
)

from pyamplicol.model_assets import bundled_model_path
from pyamplicol.model_source import compile_model_source
from pyamplicol.ufo_tensors import (
    classify_trilinear_color_expression,
    normalize_color_expression,
    normalize_lorentz_expression,
    project_trilinear_color_expression,
)


def test_ufo_metric_and_gamma_materialize_with_spenso() -> None:
    library = TensorLibrary.hep_lib_atom()
    metric = normalize_lorentz_expression(
        "UFO::Metric(UFO::idx(1,1),UFO::idx(1,2))",
        [3, 3],
    )
    gamma = normalize_lorentz_expression(
        "UFO::Gamma(UFO::idx(1,3),UFO::idx(1,1),UFO::idx(1,2))",
        [2, 2, 3],
    )

    metric_tensor = TensorNetwork(E(metric.expression), library).result_tensor(library)
    gamma_tensor = TensorNetwork(E(gamma.expression), library).result_tensor(library)

    assert len(metric_tensor) == 16
    assert len(gamma_tensor) == 64
    assert "gamma" in gamma.tensor_heads
    assert "mink" in metric.tensor_heads


def test_ufo_projector_chain_normalizes_shared_dummy_spin_index() -> None:
    normalized = normalize_lorentz_expression(
        "UFO::Gamma(UFO::idx(1,3),UFO::idx(1,1),UFO::dummy(1))"
        "*UFO::ProjM(UFO::dummy(1),UFO::idx(1,2))",
        [2, 2, 3],
    )

    assert "UFO::" not in normalized.expression
    # UFO and pyAmpliCol use opposite chirality signs for the two Weyl blocks.
    # The importer swaps P_L/P_R before spenso moves the projector through
    # gamma, so this UFO P_L chain is represented by pyAmpliCol's P_L label.
    assert "projm" in normalized.expression
    assert "chain" in normalized.expression
    assert "ufo_s_dummy_1" not in normalized.expression


def test_ufo_propagator_momentum_square_is_a_scalar_tensor_contraction() -> None:
    normalized = normalize_lorentz_expression(
        "UFO::P(UFO::idx(1,1))^2-UFO::M^2",
        [3, 3],
    )

    assert "UFO::P" not in normalized.expression
    library = TensorLibrary.hep_lib_atom()
    library.register(
        LibraryTensor.dense(
            TensorName("pyamplicol::ufo_momentum_1")(Representation.mink(4)),
            [2.0, 1.0, 0.0, 0.0],
        )
    )
    network = TensorNetwork(E(normalized.expression), library)
    network.execute(library=library)
    scalar = network.result_tensor(library).scalar().replace(S("UFO::M"), 1)
    assert float(scalar) == pytest.approx(2.0)


def test_ufo_gamma_reverses_output_input_indices_for_spenso() -> None:
    normalized = normalize_lorentz_expression(
        "UFO::Gamma(UFO::idx(1,3),UFO::idx(1,1),UFO::idx(1,2))",
        [2, 2, 3],
    )

    input_position = normalized.expression.index("ufo_s_1_2")
    output_position = normalized.expression.index("ufo_s_1_1")
    assert input_position < output_position


def test_ufo_su3_generator_and_structure_constant_are_typed() -> None:
    generator = normalize_color_expression("UFO::T(3,2,1)", [-3, 3, 8])
    structure_constant = normalize_color_expression("UFO::f(1,2,3)", [8, 8, 8])

    assert "cof(3" in generator.expression
    assert "coad(8" in generator.expression
    assert "spenso::{real,spenso::tensor}::t" in generator.expression
    assert "spenso::{antisymmetric,real,spenso::tensor}::f" in (
        structure_constant.expression
    )


def test_trilinear_color_projection_is_independent_of_source_layout() -> None:
    rewritten_f = normalize_color_expression("-UFO::f(2,1,3)", [8, 8, 8])
    scaled_generator = normalize_color_expression("2*UFO::T(3,2,1)", [-3, 3, 8])
    scaled_identity = normalize_color_expression(
        "2*UFO::Identity(1,2)",
        [-3, 3, 1],
    )
    symmetric = normalize_color_expression("UFO::d(1,2,3)", [8, 8, 8])

    assert project_trilinear_color_expression(
        rewritten_f.expression,
        [8, 8, 8],
    ) == pytest.approx(("adjoint-structure-constant", 1.0 + 0.0j))
    assert project_trilinear_color_expression(
        scaled_generator.expression,
        [-3, 3, 8],
    ) == pytest.approx(("fundamental-generator", 2.0 + 0.0j))
    assert project_trilinear_color_expression(
        scaled_identity.expression,
        [-3, 3, 1],
    ) == pytest.approx(("color-identity", 2.0 + 0.0j))
    assert project_trilinear_color_expression(
        symmetric.expression,
        [8, 8, 8],
    ) == pytest.approx(("adjoint-symmetric-invariant", 1.0 + 0.0j))


def test_trilinear_color_classifier_retains_ufo_source_fallback() -> None:
    assert classify_trilinear_color_expression(
        "writer_specific_normalized_tensor",
        "UFO::f(1, 2, 3)",
        [8, 8, 8],
    ) == pytest.approx(("adjoint-structure-constant", 1.0 + 0.0j))


def test_ufo_tensor_normalization_rejects_mistyped_indices() -> None:
    with pytest.raises(ValueError, match="not a Dirac fermion"):
        normalize_lorentz_expression(
            "UFO::Gamma(UFO::idx(1,3),UFO::idx(1,1),UFO::idx(1,2))",
            [3, 2, 3],
        )
    with pytest.raises(ValueError, match="expected adjoint"):
        normalize_color_expression("UFO::f(1,2,3)", [3, 8, 8])


@pytest.mark.parametrize("name", ["sm", "scalars", "scalar_gravity"])
def test_every_bundled_json_vertex_term_is_normalized(name: str) -> None:
    compiled = compile_model_source(
        bundled_model_path(name, "json"),
        use_cache=False,
    )

    assert compiled.ir.vertex_terms
    assert compiled.ir.max_vertex_valence == compiled.capabilities["max_vertex_valence"]
    assert all("UFO::" not in term.lorentz_expression for term in compiled.ir.vertex_terms)
    assert all("UFO::" not in term.color_expression for term in compiled.ir.vertex_terms)
    assert compiled.capabilities["compiled_vertex_term_count"] == len(
        compiled.ir.vertex_terms
    )
