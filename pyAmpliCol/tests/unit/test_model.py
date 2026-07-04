from __future__ import annotations

import math

from symbolica.community.spenso import (
    LibraryTensor,
    Representation,
    TensorName,
    TensorNetwork,
)

from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.native import (
    _gluon_tensor_to_gluon,
    _gluon_propagator,
    _quark_gluon_to_quark_weyl,
    _quark_propagator_weyl,
    _tensor_gluon_to_gluon,
    _three_gluon,
    _two_gluon_to_tensor,
)
from pyamplicol.symbols import symbols


def test_model_tables_mirror_legacy_particle_and_vertex_counts() -> None:
    model = AmplicolSMLeadingColorModel()

    assert len(model.particles) == 24
    assert len(model.vertices) == 211
    assert model.sin_weak == 0.47143025548407230
    assert model.mass(23) == 91.188
    assert model.width(23) == 2.441404
    assert model.mass(24) == 80.419002445756163


def test_model_couplings_and_leading_color_factor_match_legacy_conventions() -> None:
    model = AmplicolSMLeadingColorModel()
    left, right = model.z_fermion_coupling(1)
    prefactor = 1.0 / (model.sin_weak * model.cos_weak)

    assert math.isclose(left, prefactor * (-0.5 + model.sin_weak**2 / 3.0))
    assert math.isclose(right, prefactor * (model.sin_weak**2 / 3.0))
    assert model.photon_fermion_coupling(1) == (-1.0 / 3.0, -1.0 / 3.0)
    assert model.leading_color_factor([1, -1, 23, 21, 21]) == 27


def test_model_owns_vertex_lowering_rules_for_native_milestone() -> None:
    model = AmplicolSMLeadingColorModel()

    two_gluon_tensor = model.vertex_lowering_rule(1)
    assert two_gluon_tensor.backend == "spenso"
    assert two_gluon_tensor.full_tensor_network_ready is True
    assert two_gluon_tensor.tensor_names == ("two_gluon_to_tensor",)

    quark_gluon = model.vertex_lowering_rule(6)
    assert quark_gluon.backend == "spenso"
    assert quark_gluon.full_tensor_network_ready is True
    assert quark_gluon.expression_head == "quark_gluon_weyl_current"

    fermion_gauge = model.vertex_lowering_rule(10)
    assert fermion_gauge.backend == "spenso"
    assert fermion_gauge.full_tensor_network_ready is True
    assert fermion_gauge.expression_head == "fermion_gauge_weyl_current"

    unknown = model.vertex_lowering_rule(999)
    assert unknown.backend == "unimplemented"


def test_three_gluon_model_expression_matches_native_current_formula() -> None:
    model = AmplicolSMLeadingColorModel()
    left = (0.8 + 0.2j, -1.1 + 0.7j, 0.3 - 0.5j, 2.0 + 0.1j)
    right = (-0.4 + 1.5j, 0.9 - 0.2j, -1.2 + 0.3j, 0.1 - 0.8j)
    left_momentum = (5.0, 1.0, -2.0, 3.0)
    right_momentum = (7.0, -0.5, 0.25, -1.5)

    got = _evaluate_three_gluon_model_expression(
        model,
        left,
        right,
        left_momentum,
        right_momentum,
    )
    expected = _three_gluon(left, left_momentum, right, right_momentum)
    assert max(abs(got[i] - expected[i]) for i in range(4)) < 1.0e-14


def test_gluon_propagator_tensor_matches_native_current_formula() -> None:
    model = AmplicolSMLeadingColorModel()
    gluon = (0.8 + 0.2j, -1.1 + 0.7j, 0.3 - 0.5j, 2.0 + 0.1j)
    momentum = (5.0, 1.0, -2.0, 3.0)

    got = _evaluate_gluon_propagator_tensor(model, gluon, momentum)
    expected = _gluon_propagator(gluon, momentum)
    assert max(abs(got[i] - expected[i]) for i in range(4)) < 1.0e-14


def test_quark_weyl_propagator_tensor_matches_native_current_formula() -> None:
    model = AmplicolSMLeadingColorModel()
    quark = (1.0 + 2.0j, -0.5 + 0.25j)
    momentum = (5.0, 1.0, -2.0, 3.0)

    for chirality in (1, -1):
        got = _evaluate_quark_weyl_propagator_tensor(
            model,
            quark,
            momentum,
            chirality=chirality,
        )
        expected = _quark_propagator_weyl(quark, momentum, chirality)
        assert max(abs(got[i] - expected[i]) for i in range(2)) < 1.0e-14


def test_quark_vector_weyl_tensors_match_native_current_formula() -> None:
    model = AmplicolSMLeadingColorModel()
    quark = (1.0 + 2.0j, -0.5 + 0.25j)
    vector = (0.3 - 0.2j, 1.2 + 0.4j, -0.7 + 0.9j, 2.0 - 1.1j)

    for chirality, tensor_symbol in (
        (1, symbols.quark_vector_weyl_plus),
        (-1, symbols.quark_vector_weyl_minus),
    ):
        got = _evaluate_quark_vector_tensor(model, str(tensor_symbol), quark, vector)
        expected = _quark_gluon_to_quark_weyl(quark, vector, chirality)
        assert max(abs(got[i] - expected[i]) for i in range(2)) < 1.0e-14


def test_auxiliary_gluon_tensors_match_native_current_formulas() -> None:
    model = AmplicolSMLeadingColorModel()
    left = (0.8 + 0.2j, -1.1 + 0.7j, 0.3 - 0.5j, 2.0 + 0.1j)
    right = (-0.4 + 1.5j, 0.9 - 0.2j, -1.2 + 0.3j, 0.1 - 0.8j)
    tensor = (
        0.2 - 0.1j,
        -0.3 + 0.4j,
        0.8 + 0.7j,
        -1.1 + 0.2j,
        0.4 - 0.9j,
        1.3 + 0.5j,
    )

    got_tensor = _evaluate_two_gluon_to_tensor(model, left, right)
    expected_tensor = _two_gluon_to_tensor(left, right)
    assert max(abs(got_tensor[i] - expected_tensor[i]) for i in range(6)) < 1.0e-14

    got_left_aux = _evaluate_tensor_gluon_to_gluon(model, tensor, right)
    expected_left_aux = _tensor_gluon_to_gluon(tensor, right)
    assert (
        max(abs(got_left_aux[i] - expected_left_aux[i]) for i in range(4))
        < 1.0e-14
    )

    got_right_aux = _evaluate_gluon_tensor_to_gluon(model, left, tensor)
    expected_right_aux = _gluon_tensor_to_gluon(left, tensor)
    assert (
        max(abs(got_right_aux[i] - expected_right_aux[i]) for i in range(4))
        < 1.0e-14
    )


def _evaluate_quark_vector_tensor(
    model: AmplicolSMLeadingColorModel,
    tensor_name: str,
    quark: tuple[complex, complex],
    vector: tuple[complex, complex, complex, complex],
) -> tuple[complex, complex]:
    library = model.build_tensor_library()
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    mink = Representation.mink(4)
    library.register(LibraryTensor.dense(TensorName("test_quark")(weyl), quark))
    library.register(LibraryTensor.dense(TensorName("test_vector")(mink), vector))
    expression = (
        TensorName(tensor_name)(
            weyl("alpha"),
            mink("mu"),
            weyl("beta"),
        ).to_expression()
        * TensorName("test_quark")(weyl("alpha")).to_expression()
        * TensorName("test_vector")(mink("mu")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return complex(result[0]), complex(result[1])


def _evaluate_two_gluon_to_tensor(
    model: AmplicolSMLeadingColorModel,
    left: tuple[complex, complex, complex, complex],
    right: tuple[complex, complex, complex, complex],
) -> tuple[complex, complex, complex, complex, complex, complex]:
    library = model.build_tensor_library()
    mink = Representation.mink(4)
    aux = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    library.register(LibraryTensor.dense(TensorName("test_left")(mink), left))
    library.register(LibraryTensor.dense(TensorName("test_right")(mink), right))
    expression = (
        TensorName(str(symbols.two_gluon_to_tensor))(
            mink("left"),
            mink("right"),
            aux("aux"),
        ).to_expression()
        * TensorName("test_left")(mink("left")).to_expression()
        * TensorName("test_right")(mink("right")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return (
        complex(result[0]),
        complex(result[1]),
        complex(result[2]),
        complex(result[3]),
        complex(result[4]),
        complex(result[5]),
    )


def _evaluate_tensor_gluon_to_gluon(
    model: AmplicolSMLeadingColorModel,
    tensor: tuple[complex, complex, complex, complex, complex, complex],
    gluon: tuple[complex, complex, complex, complex],
) -> tuple[complex, complex, complex, complex]:
    library = model.build_tensor_library()
    mink = Representation.mink(4)
    aux = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    library.register(LibraryTensor.dense(TensorName("test_tensor")(aux), tensor))
    library.register(LibraryTensor.dense(TensorName("test_gluon")(mink), gluon))
    expression = (
        TensorName(str(symbols.tensor_gluon_to_gluon))(
            aux("aux"),
            mink("gluon"),
            mink("out"),
        ).to_expression()
        * TensorName("test_tensor")(aux("aux")).to_expression()
        * TensorName("test_gluon")(mink("gluon")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return (
        complex(result[0]),
        complex(result[1]),
        complex(result[2]),
        complex(result[3]),
    )


def _evaluate_gluon_tensor_to_gluon(
    model: AmplicolSMLeadingColorModel,
    gluon: tuple[complex, complex, complex, complex],
    tensor: tuple[complex, complex, complex, complex, complex, complex],
) -> tuple[complex, complex, complex, complex]:
    library = model.build_tensor_library()
    mink = Representation.mink(4)
    aux = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    library.register(LibraryTensor.dense(TensorName("test_gluon")(mink), gluon))
    library.register(LibraryTensor.dense(TensorName("test_tensor")(aux), tensor))
    expression = (
        TensorName(str(symbols.gluon_tensor_to_gluon))(
            mink("gluon"),
            aux("aux"),
            mink("out"),
        ).to_expression()
        * TensorName("test_gluon")(mink("gluon")).to_expression()
        * TensorName("test_tensor")(aux("aux")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return (
        complex(result[0]),
        complex(result[1]),
        complex(result[2]),
        complex(result[3]),
    )


def _evaluate_three_gluon_model_expression(
    model: AmplicolSMLeadingColorModel,
    left: tuple[complex, complex, complex, complex],
    right: tuple[complex, complex, complex, complex],
    left_momentum: tuple[float, float, float, float],
    right_momentum: tuple[float, float, float, float],
) -> tuple[complex, complex, complex, complex]:
    library = model.build_tensor_library()
    mink = Representation.mink(4)
    library.register(LibraryTensor.dense(TensorName("test_left")(mink), left))
    library.register(LibraryTensor.dense(TensorName("test_right")(mink), right))
    library.register(
        LibraryTensor.dense(TensorName("test_left_momentum")(mink), left_momentum)
    )
    library.register(
        LibraryTensor.dense(TensorName("test_right_momentum")(mink), right_momentum)
    )
    expression = (
        model.three_gluon_current_expression(
            left_slot=mink("left"),
            right_slot=mink("right"),
            output_slot=mink("out"),
            left_momentum_tensor_name="test_left_momentum",
            right_momentum_tensor_name="test_right_momentum",
            dummy_prefix="unit_test",
        )
        * TensorName("test_left")(mink("left")).to_expression()
        * TensorName("test_right")(mink("right")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return (
        complex(result[0]),
        complex(result[1]),
        complex(result[2]),
        complex(result[3]),
    )


def _evaluate_gluon_propagator_tensor(
    model: AmplicolSMLeadingColorModel,
    gluon: tuple[complex, complex, complex, complex],
    momentum: tuple[float, float, float, float],
) -> tuple[complex, complex, complex, complex]:
    library = model.build_tensor_library()
    mink = Representation.mink(4)
    library.register(LibraryTensor.dense(TensorName("test_gluon")(mink), gluon))
    library.register(
        LibraryTensor.dense(
            TensorName("test_gluon_propagator")(mink, mink),
            model.gluon_propagator_tensor_data(momentum),
        )
    )
    expression = (
        TensorName("test_gluon_propagator")(
            mink("input"),
            mink("output"),
        ).to_expression()
        * TensorName("test_gluon")(mink("input")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return (
        complex(result[0]),
        complex(result[1]),
        complex(result[2]),
        complex(result[3]),
    )


def _evaluate_quark_weyl_propagator_tensor(
    model: AmplicolSMLeadingColorModel,
    quark: tuple[complex, complex],
    momentum: tuple[float, float, float, float],
    *,
    chirality: int,
) -> tuple[complex, complex]:
    library = model.build_tensor_library()
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    library.register(LibraryTensor.dense(TensorName("test_quark")(weyl), quark))
    library.register(
        LibraryTensor.dense(
            TensorName("test_quark_propagator")(weyl, weyl),
            model.quark_weyl_propagator_tensor_data(
                momentum,
                chirality=chirality,
            ),
        )
    )
    expression = (
        TensorName("test_quark_propagator")(
            weyl("input"),
            weyl("output"),
        ).to_expression()
        * TensorName("test_quark")(weyl("input")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    return complex(result[0]), complex(result[1])
