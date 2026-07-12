from __future__ import annotations

import math

from symbolica.community.spenso import Representation, TensorName, TensorLibrary, TensorNetwork

from pyamplicol.params import ParamBuilder, SymbolicaEvaluatorBundle


def test_param_builder_reuses_symbols_created_when_ranges_are_added() -> None:
    builder = ParamBuilder()
    first = builder.add_parameter_list(("cached", "first"), 2)
    second = builder.add_parameter_list(("cached", "second"), 1)

    symbols = builder.parameter_symbols()

    assert len(symbols) == 3
    assert symbols[0] is first[0]
    assert symbols[1] is first[1]
    assert symbols[2] is second[0]


def test_param_builder_registers_rank1_tensor_for_symbolica_evaluator() -> None:
    library = TensorLibrary.hep_lib_atom()
    mink = Representation.mink(4)
    tensor_name = "pyamplicol::external_momentum_1"
    builder = ParamBuilder()

    builder.register_rank1_tensor(
        library,
        tensor_name=tensor_name,
        representation=mink,
        head=("external_momentum", "1"),
        length=4,
        role="external_momentum",
        real_valued=True,
    )
    assert builder.real_valued_inputs == [0, 1, 2, 3]
    momentum = TensorName(tensor_name)
    expression = (
        momentum(mink("mu")).to_expression()
        * momentum(mink("mu")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    scalar = network.result_scalar()
    evaluator = scalar.evaluator(builder.parameter_symbols())
    bundle = SymbolicaEvaluatorBundle(
        name="external_momentum_square",
        evaluator=evaluator,
        param_builder=builder,
        metadata={"source": "unit-test"},
    )

    bundle.param_builder.set_parameter_values(
        ("external_momentum", "1"),
        [5.0, 1.0, 2.0, 3.0],
    )
    assert math.isclose(bundle.evaluate()[0].real, 25.0 - 1.0 - 4.0 - 9.0)

    loaded = SymbolicaEvaluatorBundle.from_artifact_payload(
        bundle.to_artifact_payload()
    )
    assert loaded.param_builder.real_valued_inputs == [0, 1, 2, 3]
    loaded.param_builder.set_parameter_values(
        ("external_momentum", "1"),
        [7.0, 2.0, 3.0, 6.0],
    )
    assert math.isclose(loaded.evaluate()[0].real, 49.0 - 4.0 - 9.0 - 36.0)


def test_symbolica_evaluator_bundle_supports_complex_rank1_tensor_inputs() -> None:
    library = TensorLibrary.hep_lib_atom()
    mink = Representation.mink(4)
    tensor_name = "pyamplicol::external_polarization_1"
    builder = ParamBuilder()

    builder.register_rank1_tensor(
        library,
        tensor_name=tensor_name,
        representation=mink,
        head=("external_polarization", "1"),
        length=4,
        role="external_polarization",
    )
    polarization = TensorName(tensor_name)
    expression = (
        polarization(mink("mu")).to_expression()
        * polarization(mink("mu")).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    scalar = network.result_scalar()
    evaluator = scalar.evaluator(builder.parameter_symbols())
    bundle = SymbolicaEvaluatorBundle(
        name="external_polarization_square",
        evaluator=evaluator,
        param_builder=builder,
        complex_inputs=True,
        metadata={"source": "unit-test"},
    )

    values = [1.0 + 2.0j, 3.0 - 1.0j, -2.0 + 0.5j, 0.25 - 4.0j]
    bundle.param_builder.set_parameter_values(("external_polarization", "1"), values)
    expected = values[0] ** 2 - values[1] ** 2 - values[2] ** 2 - values[3] ** 2
    assert abs(bundle.evaluate()[0] - expected) < 1.0e-12

    loaded = SymbolicaEvaluatorBundle.from_artifact_payload(
        bundle.to_artifact_payload()
    )
    loaded.param_builder.set_parameter_values(("external_polarization", "1"), values)
    assert loaded.complex_inputs
    assert abs(loaded.evaluate()[0] - expected) < 1.0e-12
