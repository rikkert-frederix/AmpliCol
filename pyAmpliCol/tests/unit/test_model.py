from __future__ import annotations

import math
from types import SimpleNamespace

from symbolica.community.spenso import (
    LibraryTensor,
    Representation,
    TensorName,
    TensorNetwork,
)

from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.native import (
    _antilepton_lepton_to_vector_weyl,
    _gluon_tensor_to_gluon,
    _gluon_propagator,
    _lepton_antilepton_to_vector_weyl,
    _massive_vector_propagator,
    _quark_gluon_to_quark_weyl,
    _quark_gluon_to_quark_coupl_weyl,
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


def test_model_contains_legacy_particle_vocabulary_and_charged_current_vertices() -> None:
    model = AmplicolSMLeadingColorModel()

    expected_pdgs = {
        -24,
        -16,
        -15,
        -14,
        -13,
        -12,
        -11,
        -6,
        -5,
        -4,
        -3,
        -2,
        -1,
        1,
        2,
        3,
        4,
        5,
        6,
        11,
        12,
        13,
        14,
        15,
        16,
        21,
        22,
        23,
        24,
        25,
    }

    for pdg in expected_pdgs:
        assert model.particle(pdg).pdg == pdg or model.particle(pdg).anti_pdg == pdg

    charged_current_vertices = {
        vertex.particles: vertex.coupling
        for vertex in model.vertices
        if vertex.kind in {10, 11} and 24 in {abs(pdg) for pdg in vertex.particles}
    }
    assert charged_current_vertices[(1, 24, 2)] == (
        model.charged_current_coupling(),
        0.0,
    )
    assert charged_current_vertices[(2, -24, 1)] == (
        model.charged_current_coupling(),
        0.0,
    )
    assert charged_current_vertices[(-2, 24, -1)] == (
        model.charged_current_coupling(),
        0.0,
    )
    assert charged_current_vertices[(-1, -24, -2)] == (
        model.charged_current_coupling(),
        0.0,
    )
    neutral_current_vertices = {
        vertex.particles: vertex.coupling
        for vertex in model.vertices
        if vertex.kind in {21, 22} and vertex.particles[2] == 23
    }
    assert (12, -12, 23) not in neutral_current_vertices
    assert (-12, 12, 23) not in neutral_current_vertices
    assert not model.vertices_accepting(12, 23)
    assert not model.vertices_accepting(23, 12)
    assert not model.vertices_accepting(-12, 23)
    assert not model.vertices_accepting(23, -12)


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


def test_model_reports_legacy_vertex_lowering_coverage() -> None:
    model = AmplicolSMLeadingColorModel()

    report = model.vertex_lowering_coverage()

    assert report.ready_kinds == (
        0,
        1,
        2,
        3,
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        14,
        15,
        16,
        17,
        18,
        19,
        20,
        21,
        22,
        23,
        24,
    )
    assert report.pending_kinds == ()
    assert report.unimplemented_kinds == ()
    by_kind = {entry.kind: entry for entry in report.entries}
    assert by_kind[4].backend == "symbolica"
    assert by_kind[4].full_tensor_network_ready is True
    assert by_kind[8].kernel == "fermion_pair_to_vector"
    assert by_kind[16].kernel == "fermion_scalar_to_fermion"
    assert by_kind[16].full_tensor_network_ready is True
    assert by_kind[8].full_tensor_network_ready is True
    assert by_kind[10].vertex_count == 24
    assert by_kind[12].kernel == "three_vector_current"
    assert by_kind[12].backend == "symbolica"
    assert by_kind[12].full_tensor_network_ready is True
    assert by_kind[12].input_roles == ("vector", "vector")
    assert by_kind[12].output_role == "vector"
    assert by_kind[12].coupling_mode == "vertex"
    assert model.vertex_lowering_rule(12).to_json_dict()["kernel"] == "three_vector_current"
    assert by_kind[17].kernel == "two_vector_to_scalar"
    assert by_kind[21].kernel == "fermion_pair_to_vector"
    assert by_kind[21].vertex_count == 9
    assert by_kind[21].full_tensor_network_ready is True
    assert report.to_json_dict()["model"] == "amplicol-sm-leading-color"


def test_model_keeps_auxiliary_u1_vertices_available_for_lc_construction() -> None:
    model = AmplicolSMLeadingColorModel()

    lc_vertices = model.iter_vertices(color_accuracy="lc")
    full_vertices = model.iter_vertices(color_accuracy="full")

    assert len(lc_vertices) == len(full_vertices)
    assert any(99 in vertex.particles for vertex in lc_vertices)
    assert any(99 in vertex.particles for vertex in full_vertices)
    assert [vertex.kind for vertex in model.vertices_for_inputs(1, -1, color_accuracy="lc")] == [8]
    assert [vertex.kind for vertex in model.vertices_for_inputs(1, -1, color_accuracy="full")] == [8]


def test_model_owns_source_states_dimensions_and_quantum_flow_filtering() -> None:
    model = AmplicolSMLeadingColorModel()

    assert [state.chirality for state in model.source_spin_states(1)] == [-1, 1]
    assert [state.helicity for state in model.source_spin_states(23)] == [-1, 0, 1]
    assert model.current_dimension(1, -1) == 2
    assert model.current_dimension(21) == 4
    assert model.auxiliary_kind(-21) == "antisymmetric-tensor"
    assert model.charge_units(1) == -1

    charged_current = model.vertices_accepting(1, 24)[0]
    vector_index = SimpleNamespace(
        particle_id=24,
        chirality=0,
        flavour_flow=(24,),
    )
    left_handed_down = SimpleNamespace(
        particle_id=1,
        chirality=-1,
        flavour_flow=(1,),
    )
    right_handed_down = SimpleNamespace(
        particle_id=1,
        chirality=1,
        flavour_flow=(1,),
    )

    allowed = model.allowed_quantum_flows(
        charged_current,
        left_handed_down,
        vector_index,
    )
    assert len(allowed) == 1
    assert allowed[0].chirality == -1
    assert allowed[0].flavour_flow == (1, 2)
    assert allowed[0].charge_flow == 2
    assert allowed[0].coupling == charged_current.coupling

    assert model.allowed_quantum_flows(
        charged_current,
        right_handed_down,
        vector_index,
    ) == ()

    neutral_vector = model.vertices_accepting(-11, 11)[0]
    positron = SimpleNamespace(
        particle_id=-11,
        chirality=1,
        flavour_flow=(-11,),
    )
    electron = SimpleNamespace(
        particle_id=11,
        chirality=-1,
        flavour_flow=(11,),
    )
    neutral_flows = model.allowed_quantum_flows(
        neutral_vector,
        positron,
        electron,
    )
    assert neutral_flows
    assert neutral_flows[0].flavour_flow == (-11, 11, 22)


def test_model_owns_propagator_lowering_rules() -> None:
    model = AmplicolSMLeadingColorModel()

    gluon = model.propagator_lowering_rule(21)
    assert gluon.applies_propagator is True
    assert gluon.full_tensor_network_ready is True
    assert gluon.kernel == "massless_vector_feynman_gauge"

    z_boson = model.propagator_lowering_rule(23)
    assert z_boson.applies_propagator is True
    assert z_boson.full_tensor_network_ready is True
    assert z_boson.kernel == "massive_vector_unitary_gauge"

    quark = model.propagator_lowering_rule(1, -1)
    assert quark.applies_propagator is True
    assert quark.full_tensor_network_ready is True
    assert quark.kernel == "weyl_fermion"

    auxiliary = model.propagator_lowering_rule(-21)
    assert auxiliary.applies_propagator is False
    assert auxiliary.full_tensor_network_ready is True
    assert auxiliary.kernel == "auxiliary_tensor_embedded_propagator"

    higgsor = model.propagator_lowering_rule(126)
    assert higgsor.applies_propagator is False
    assert higgsor.full_tensor_network_ready is True
    assert higgsor.kernel == "auxiliary_scalar_no_propagator"


def test_model_propagator_component_expressions_match_native_formulas() -> None:
    model = AmplicolSMLeadingColorModel()
    vector = (0.8 + 0.2j, -1.1 + 0.7j, 0.3 - 0.5j, 2.0 + 0.1j)
    momentum = (5.0, 1.0, -2.0, 3.0)
    quark = (1.0 + 2.0j, -0.5 + 0.25j)

    gluon_components = model.propagator_component_expression(
        21,
        vector,
        momentum,
    )
    assert all(isinstance(component, complex) for component in gluon_components)
    _assert_components_close(
        gluon_components,
        _gluon_propagator(vector, momentum),
    )
    _assert_components_close(
        model.propagator_component_expression(23, vector, momentum),
        _massive_vector_propagator(
            vector,
            momentum,
            model.mass(23),
            model.width(23),
        ),
    )
    for chirality in (1, -1):
        _assert_components_close(
            model.propagator_component_expression(
                1,
                quark,
                momentum,
                chirality=chirality,
            ),
            _quark_propagator_weyl(quark, momentum, chirality),
        )

    scalar = (1.7 - 0.2j,)
    denominator = (
        _minkowski_dot(momentum, momentum)
        - model.mass(25) ** 2
        + 1j * model.mass(25) * model.width(25)
    )
    _assert_components_close(
        model.propagator_component_expression(25, scalar, momentum),
        (1j * scalar[0] / denominator,),
    )
    _assert_components_close(
        model.propagator_component_expression(-21, vector[:2], momentum),
        vector[:2],
    )
    _assert_components_close(
        model.propagator_component_expression(126, scalar, momentum),
        scalar,
    )


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


def test_model_vertex_component_expression_covers_ready_fermion_vector_kernels() -> None:
    model = AmplicolSMLeadingColorModel()
    fermion = (1.0 + 2.0j, -0.5 + 0.25j)
    antifermion = (0.75 - 0.3j, -1.1 + 0.6j)
    vector = (0.3 - 0.2j, 1.2 + 0.4j, -0.7 + 0.9j, 2.0 - 1.1j)
    coupling = (2.0, -3.0)

    for chirality in (1, -1):
        _assert_components_close(
            model.vertex_component_expression(
                4,
                vector,
                fermion,
                result_particle_id=1,
                result_chirality=chirality,
            ),
            _quark_gluon_to_quark_weyl(fermion, vector, chirality),
        )
        _assert_components_close(
            model.vertex_component_expression(
                6,
                fermion,
                vector,
                result_particle_id=1,
                result_chirality=chirality,
            ),
            _quark_gluon_to_quark_weyl(fermion, vector, chirality),
        )
        _assert_components_close(
            model.vertex_component_expression(
                10,
                fermion,
                vector,
                result_particle_id=1,
                result_chirality=chirality,
                coupling=coupling,
            ),
            _quark_gluon_to_quark_coupl_weyl(
                fermion,
                vector,
                coupling,
                chirality,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                23,
                vector,
                fermion,
                result_particle_id=1,
                result_chirality=chirality,
                coupling=coupling,
            ),
            _quark_gluon_to_quark_coupl_weyl(
                fermion,
                vector,
                coupling,
                chirality,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                5,
                vector,
                antifermion,
                result_particle_id=-1,
                result_chirality=chirality,
            ),
            _expected_antifermion_vector_weyl(
                antifermion,
                vector,
                chirality,
                coupling=None,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                7,
                antifermion,
                vector,
                result_particle_id=-1,
                result_chirality=chirality,
            ),
            _expected_antifermion_vector_weyl(
                antifermion,
                vector,
                chirality,
                coupling=None,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                11,
                antifermion,
                vector,
                result_particle_id=-1,
                result_chirality=chirality,
                coupling=coupling,
            ),
            _expected_antifermion_vector_weyl(
                antifermion,
                vector,
                chirality,
                coupling=coupling,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                24,
                vector,
                antifermion,
                result_particle_id=-1,
                result_chirality=chirality,
                coupling=coupling,
            ),
            _expected_antifermion_vector_weyl(
                antifermion,
                vector,
                chirality,
                coupling=coupling,
            ),
        )


def test_model_vertex_component_expression_covers_ready_fermion_pair_kernels() -> None:
    model = AmplicolSMLeadingColorModel()
    fermion = (1.0 + 2.0j, -0.5 + 0.25j)
    antifermion = (0.75 - 0.3j, -1.1 + 0.6j)
    coupling = (2.0, -3.0)

    for fermion_chirality, antifermion_chirality in ((-1, 1), (1, -1), (1, 1)):
        _assert_components_close(
            model.vertex_component_expression(
                21,
                fermion,
                antifermion,
                result_particle_id=23,
                result_chirality=0,
                left_chirality=fermion_chirality,
                right_chirality=antifermion_chirality,
                coupling=coupling,
            ),
            _lepton_antilepton_to_vector_weyl(
                fermion,
                antifermion,
                coupling,
                fermion_chirality,
                antifermion_chirality,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                22,
                antifermion,
                fermion,
                result_particle_id=23,
                result_chirality=0,
                left_chirality=antifermion_chirality,
                right_chirality=fermion_chirality,
                coupling=coupling,
            ),
            _antilepton_lepton_to_vector_weyl(
                antifermion,
                fermion,
                coupling,
                antifermion_chirality,
                fermion_chirality,
            ),
        )
        _assert_components_close(
            model.vertex_component_expression(
                9,
                antifermion,
                fermion,
                result_particle_id=21,
                result_chirality=0,
                left_chirality=antifermion_chirality,
                right_chirality=fermion_chirality,
            ),
            _antilepton_lepton_to_vector_weyl(
                antifermion,
                fermion,
                (1.0, 1.0),
                antifermion_chirality,
                fermion_chirality,
            ),
        )


def test_model_vertex_component_expression_covers_coupled_vector_scalar_kernels() -> None:
    model = AmplicolSMLeadingColorModel()
    left = (0.8 + 0.2j, -1.1 + 0.7j, 0.3 - 0.5j, 2.0 + 0.1j)
    right = (-0.4 + 1.5j, 0.9 - 0.2j, -1.2 + 0.3j, 0.1 - 0.8j)
    left_momentum = (5.0, 1.0, -2.0, 3.0)
    right_momentum = (7.0, -0.5, 0.25, -1.5)
    tensor = (
        0.2 - 0.1j,
        -0.3 + 0.4j,
        0.8 + 0.7j,
        -1.1 + 0.2j,
        0.4 - 0.9j,
        1.3 + 0.5j,
    )
    coupling = (2.5, 0.0)

    _assert_components_close(
        model.vertex_component_expression(
            12,
            left,
            right,
            result_particle_id=23,
            result_chirality=0,
            coupling=coupling,
            left_momentum=left_momentum,
            right_momentum=right_momentum,
        ),
        _expected_coupled_three_vector(
            left,
            left_momentum,
            right,
            right_momentum,
            coupling,
        ),
    )
    _assert_components_close(
        model.vertex_component_expression(
            13,
            left,
            right,
            result_particle_id=-23,
            result_chirality=0,
            coupling=coupling,
        ),
        tuple(coupling[0] * component for component in _two_gluon_to_tensor(left, right)),
    )
    _assert_components_close(
        model.vertex_component_expression(
            14,
            tensor,
            right,
            result_particle_id=23,
            result_chirality=0,
            coupling=coupling,
        ),
        tuple(coupling[0] * component for component in _tensor_gluon_to_gluon(tensor, right)),
    )
    _assert_components_close(
        model.vertex_component_expression(
            15,
            left,
            tensor,
            result_particle_id=23,
            result_chirality=0,
            coupling=coupling,
        ),
        tuple(coupling[0] * component for component in _gluon_tensor_to_gluon(left, tensor)),
    )

    scalar = (1.7 - 0.2j,)
    other_scalar = (-0.6 + 0.4j,)
    prefactor = 1j / math.sqrt(2.0)
    _assert_components_close(
        model.vertex_component_expression(
            17,
            left,
            right,
            result_particle_id=25,
            result_chirality=0,
            coupling=coupling,
        ),
        (prefactor * coupling[0] * _minkowski_dot(left, right),),
    )
    _assert_components_close(
        model.vertex_component_expression(
            18,
            scalar,
            right,
            result_particle_id=23,
            result_chirality=0,
            coupling=coupling,
        ),
        tuple(prefactor * coupling[0] * scalar[0] * component for component in right),
    )
    _assert_components_close(
        model.vertex_component_expression(
            19,
            left,
            scalar,
            result_particle_id=23,
            result_chirality=0,
            coupling=coupling,
        ),
        tuple(prefactor * coupling[0] * scalar[0] * component for component in left),
    )
    _assert_components_close(
        model.vertex_component_expression(
            20,
            scalar,
            other_scalar,
            result_particle_id=25,
            result_chirality=0,
            coupling=(3.0, -10.0),
        ),
        (prefactor * 1j * 3.0 * scalar[0] * other_scalar[0],),
    )


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


def _assert_components_close(
    got: tuple[complex, ...],
    expected: tuple[complex, ...],
    *,
    tolerance: float = 1.0e-14,
) -> None:
    assert len(got) == len(expected)
    assert max(abs(complex(got[i]) - complex(expected[i])) for i in range(len(got))) < tolerance


def _expected_antifermion_vector_weyl(
    antifermion: tuple[complex, complex],
    vector: tuple[complex, complex, complex, complex],
    chirality: int,
    *,
    coupling: tuple[float, float] | None,
) -> tuple[complex, complex]:
    tmp1, tmp2, tmp3, tmp4 = _vector_slash_terms(vector)
    prefactor = 1j / math.sqrt(2.0)
    a1, a2 = antifermion
    if chirality == 1:
        factor = prefactor if coupling is None else prefactor * coupling[0]
        return (
            factor * (tmp1 * a1 + tmp4 * a2),
            factor * (tmp2 * a2 + tmp3 * a1),
        )
    if chirality == -1:
        factor = prefactor if coupling is None else prefactor * coupling[1]
        return (
            factor * (tmp2 * a1 - tmp4 * a2),
            factor * (tmp1 * a2 - tmp3 * a1),
        )
    raise AssertionError(f"unexpected chirality: {chirality}")


def _vector_slash_terms(
    vector: tuple[complex, complex, complex, complex],
) -> tuple[complex, complex, complex, complex]:
    v0, v1, v2, v3 = vector
    return v0 + v3, v0 - v3, v1 + 1j * v2, v1 - 1j * v2


def _expected_coupled_three_vector(
    left: tuple[complex, complex, complex, complex],
    left_momentum: tuple[float, float, float, float],
    right: tuple[complex, complex, complex, complex],
    right_momentum: tuple[float, float, float, float],
    coupling: tuple[float, float],
) -> tuple[complex, complex, complex, complex]:
    tmp1 = _minkowski_dot(left, right)
    tmp2_momentum = tuple(
        2.0 * right_momentum[index] + left_momentum[index]
        for index in range(4)
    )
    tmp3_momentum = tuple(
        -2.0 * left_momentum[index] - right_momentum[index]
        for index in range(4)
    )
    tmp2 = _minkowski_dot(left, tmp2_momentum)
    tmp3 = _minkowski_dot(right, tmp3_momentum)
    prefactor = (1j / math.sqrt(2.0)) * coupling[0]
    return tuple(
        prefactor
        * (
            tmp1 * (left_momentum[index] - right_momentum[index])
            + tmp2 * right[index]
            + tmp3 * left[index]
        )
        for index in range(4)
    )


def _minkowski_dot(
    left: tuple[complex, ...],
    right: tuple[complex | float, ...],
) -> complex:
    return left[0] * right[0] - left[1] * right[1] - left[2] * right[2] - left[3] * right[3]
