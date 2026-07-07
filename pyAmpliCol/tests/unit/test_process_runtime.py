from __future__ import annotations

import json
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest
import numpy as np

import pyamplicol.dag_runtime as dag_runtime_module
from pyamplicol.dag_runtime import (
    _RUSTICOL_STANDALONE_CHECK_SCRIPT,
    _build_shared_global_parameter_layout,
    _build_shared_helicity_current_table,
    _charged_leptonic_w_vector_source_specs,
    _fill_shared_source_currents,
    _neutral_dilepton_vector_source_specs,
    _rusticol_layout_manifest,
    _rusticol_model_manifest,
    _rusticol_table_manifest,
    _write_rusticol_process_artifacts,
    _shared_current_table_from_rusticol_manifest,
    SharedCurrentNode,
    SharedCurrentTable,
    SharedSourceCurrent,
    _weighted_abs2_sums,
    _write_zero_gluon_rusticol_process_artifacts,
    _z_gluon_graph,
    ZGluonDAGEvaluator,
)
from pyamplicol.generic_artifact import (
    write_generic_dag_process_artifact,
    write_generic_dag_process_set_artifact,
)
from pyamplicol.legacy_matrix import CurrentKey
from pyamplicol.model import AmplicolSMLeadingColorModel
from pyamplicol.native import (
    LeadingColorZJetsNativeEvaluator,
    NativeEvaluationError,
    _antilepton_lepton_to_vector_weyl,
    _ext_antiquark_weyl,
    _ext_quark_weyl,
    _neutral_vector_propagator,
)
from pyamplicol.process_ir import build_process_set_ir
from pyamplicol.process_runtime import load_process, load_process_manifest


_LEGACY_SCHEMA_V1_TEST = pytest.mark.skip(
    reason=(
        "schema-v1 native/Z-family runtime is deprecated; production coverage "
        "uses schema-v2 generic DAG artifacts"
    )
)


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_process_manifest_helpers_are_process_shaped() -> None:
    model = AmplicolSMLeadingColorModel()
    graph = _z_gluon_graph("d d~ > z g", model)
    table = _build_shared_helicity_current_table(graph, model, gluon_count=1)
    layout = _build_shared_global_parameter_layout(table)

    layout_payload = _rusticol_layout_manifest(table, layout)
    table_payload = _rusticol_table_manifest(table)

    json.dumps(layout_payload)
    json.dumps(table_payload)
    model_payload = _rusticol_model_manifest(model)
    json.dumps(model_payload)

    assert layout_payload["parameter_count"] == layout.parameter_count
    assert len(layout_payload["current_offsets"]) == len(table.currents)
    assert table_payload["currents"][0]["id"] == 0
    assert table_payload["sources"]
    assert table_payload["amplitudes"]
    particles_by_pdg = {
        particle["pdg"]: particle
        for particle in model_payload["particles"]
    }
    assert {1, 21, 23}.issubset(particles_by_pdg)
    assert particles_by_pdg[1]["anti_pdg"] == -1
    assert model_payload["vertices"]
    assert model_payload["mass_w"] == model.mass(24)
    assert "import rusticol" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "rusticol.Runtime.load_legacy(str(root))" in (
        _RUSTICOL_STANDALONE_CHECK_SCRIPT
    )
    assert "--precision" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "--profile" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "--rusticol-folder" in _RUSTICOL_STANDALONE_CHECK_SCRIPT
    assert "runtime.evaluate_with_prec(points, precision)" in (
        _RUSTICOL_STANDALONE_CHECK_SCRIPT
    )
    assert "runtime.profile(batch, precision=precision)" in (
        _RUSTICOL_STANDALONE_CHECK_SCRIPT
    )


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_source_manifest_preserves_composite_source_metadata() -> None:
    model = AmplicolSMLeadingColorModel()
    graph = _z_gluon_graph("d d~ > z g", model)
    table = _build_shared_helicity_current_table(graph, model, gluon_count=1)
    table_payload = _rusticol_table_manifest(table)
    source_payload = table_payload["sources"][0]
    source_payload.update(
        {
            "source_kind": "lepton_pair_vector",
            "partner_leg_label": 5,
            "partner_helicity": -1,
            "partner_chirality": -1,
            "vector_pdg": 24,
            "coupling": [1.23, 0.45],
        }
    )
    table_payload["amplitudes"][0]["coherent_group_id"] = 13

    reloaded = _shared_current_table_from_rusticol_manifest({"table": table_payload})
    source = reloaded.sources[0]

    assert source.source_kind == "lepton_pair_vector"
    assert source.partner_leg_label == 5
    assert source.partner_helicity == -1
    assert source.partner_chirality == -1
    assert source.vector_pdg == 24
    assert source.coupling == pytest.approx((1.23, 0.45))
    assert reloaded.sources[1].source_kind == "external"
    assert reloaded.amplitudes[0].coherent_group_id == 13
    assert reloaded.amplitudes[1].coherent_group_id is None


def test_weighted_abs2_sums_supports_coherent_groups() -> None:
    evaluated = np.asarray(
        [
            [1.0 + 2.0j, 3.0 - 1.0j, 2.0 + 0.0j],
            [0.5 + 0.5j, -0.5 + 0.5j, 1.0j],
        ],
        dtype=np.complex128,
    )
    weights = np.asarray([2.0, 2.0, 3.0], dtype=np.float64)

    raw_sums = _weighted_abs2_sums(
        evaluated,
        weights,
        output_length=3,
        group_ids=(7, 7, None),
    )

    assert raw_sums == pytest.approx(
        (
            2.0 * abs((1.0 + 2.0j) + (3.0 - 1.0j)) ** 2
            + 3.0 * abs(2.0 + 0.0j) ** 2,
            2.0 * abs((0.5 + 0.5j) + (-0.5 + 0.5j)) ** 2
            + 3.0 * abs(1.0j) ** 2,
        )
    )


@_LEGACY_SCHEMA_V1_TEST
def test_fill_shared_source_currents_supports_charged_lepton_pair_vector() -> None:
    model = AmplicolSMLeadingColorModel()
    reference = LeadingColorZJetsNativeEvaluator(model=model)
    particles = reference.canonical_charged_leptonic_w_gluon_point(
        "u d~ > e+ ve",
        vector_pdg=24,
        first_lepton_pdg=-11,
        second_lepton_pdg=12,
        gluon_count=0,
        sqrt_s=1000.0,
    )
    current = SharedCurrentNode(
        id=0,
        key=CurrentKey(24, (3, 4), 0),
        ext_source_bits=1,
        source_ids=(0,),
        is_source=True,
        needs_propagator=False,
        dimension=4,
    )
    source = SharedSourceCurrent(
        current_id=0,
        leg_label=3,
        helicity=1,
        physical_helicity=1,
        chirality=1,
        source_bit=1,
        source_kind="lepton_pair_vector",
        partner_leg_label=4,
        partner_helicity=-1,
        partner_chirality=-1,
        vector_pdg=24,
    )
    table = SharedCurrentTable(
        currents=(current,),
        sources=(source,),
        interactions=(),
        interactions_by_result=((),),
        amplitudes=(),
    )
    values = np.full((1, 4), 99.0 + 0.0j, dtype=np.complex128)

    _fill_shared_source_currents(values, table, particles, model, gluon_count=0)

    first = particles[2]
    second = particles[3]
    current_unpropagated = _antilepton_lepton_to_vector_weyl(
        _ext_antiquark_weyl(first.momentum, 1, 1),
        _ext_quark_weyl(second.momentum, -1, -1),
        (model.charged_current_coupling(), 0.0),
        1,
        -1,
    )
    expected = _neutral_vector_propagator(
        current_unpropagated,
        tuple(
            first.momentum[index] + second.momentum[index]
            for index in range(4)
        ),
        vector_pdg=24,
        model=model,
    )
    np.testing.assert_allclose(values[0], expected, rtol=0.0, atol=1.0e-13)


@_LEGACY_SCHEMA_V1_TEST
def test_shared_current_table_can_use_charged_leptonic_w_source_specs() -> None:
    model = AmplicolSMLeadingColorModel()
    graph = _z_gluon_graph("u d~ > w+ g", model)
    table = _build_shared_helicity_current_table(
        graph,
        model,
        gluon_count=1,
        vector_source_specs=_charged_leptonic_w_vector_source_specs(
            vector_pdg=24,
            first_lepton_label=4,
            second_lepton_label=5,
            model=model,
        ),
    )

    composite_sources = [
        source for source in table.sources if source.source_kind == "lepton_pair_vector"
    ]
    assert len(composite_sources) == 1
    source = composite_sources[0]
    assert source.leg_label == 4
    assert source.partner_leg_label == 5
    assert source.helicity == 1
    assert source.partner_helicity == -1
    assert source.vector_pdg == 24
    assert source.coupling == pytest.approx((model.charged_current_coupling(), 0.0))
    assert table.currents[source.current_id].external_labels == (4, 5)
    assert any(
        current.external_labels[-2:] == (4, 5)
        and len(current.external_labels) > 2
        for current in table.currents
    )
    assert {len(amplitude.helicities) for amplitude in table.amplitudes} == {5}
    assert len(table.amplitudes) == 4

    reference = LeadingColorZJetsNativeEvaluator(model=model)
    particles = reference.canonical_charged_leptonic_w_gluon_point(
        "u d~ > e+ ve g",
        vector_pdg=24,
        first_lepton_pdg=-11,
        second_lepton_pdg=12,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    values = np.zeros((len(table.currents), 6), dtype=np.complex128)
    _fill_shared_source_currents(values, table, particles, model, gluon_count=1)
    assert np.any(np.abs(values[source.current_id, :4]) > 0.0)


@_LEGACY_SCHEMA_V1_TEST
def test_shared_current_table_uses_coherent_neutral_dilepton_sources() -> None:
    model = AmplicolSMLeadingColorModel()
    graph = _z_gluon_graph("d d~ > z g", model)
    table = _build_shared_helicity_current_table(
        graph,
        model,
        gluon_count=1,
        vector_source_specs=_neutral_dilepton_vector_source_specs(
            lepton_pdg=11,
            antilepton_pdg=-11,
            lepton_label=4,
            antilepton_label=5,
            model=model,
        ),
    )

    composite_sources = [
        source for source in table.sources if source.source_kind == "lepton_pair_vector"
    ]
    assert len(composite_sources) == 8
    assert {source.vector_pdg for source in composite_sources} == {22, 23}
    assert {
        table.currents[source.current_id].external_labels
        for source in composite_sources
    } == {(4, 5)}
    assert all(source.coupling is not None for source in composite_sources)
    assert len(table.amplitudes) == 32
    group_counts: dict[int, int] = {}
    for amplitude in table.amplitudes:
        assert amplitude.coherent_group_id is not None
        group_counts[amplitude.coherent_group_id] = (
            group_counts.get(amplitude.coherent_group_id, 0) + 1
        )
    assert len(group_counts) == 16
    assert set(group_counts.values()) == {2}

    reference = LeadingColorZJetsNativeEvaluator(model=model)
    particles = reference.canonical_neutral_dilepton_gluon_point(
        "d d~ > e+ e- g",
        lepton_pdg=11,
        antilepton_pdg=-11,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    values = np.zeros((len(table.currents), 6), dtype=np.complex128)
    _fill_shared_source_currents(values, table, particles, model, gluon_count=1)
    nonzero_sources = [
        source
        for source in composite_sources
        if np.any(np.abs(values[source.current_id, :4]) > 0.0)
    ]
    assert len(nonzero_sources) == 4
    assert {source.vector_pdg for source in nonzero_sources} == {22, 23}


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_manifest_preserves_explicit_lepton_family(
    tmp_path: Path,
) -> None:
    model = AmplicolSMLeadingColorModel()
    graph = _z_gluon_graph("d d~ > z g", model)
    table = _build_shared_helicity_current_table(
        graph,
        model,
        gluon_count=1,
        vector_source_specs=_neutral_dilepton_vector_source_specs(
            lepton_pdg=11,
            antilepton_pdg=-11,
            lepton_label=4,
            antilepton_label=5,
            model=model,
        ),
    )
    layout = _build_shared_global_parameter_layout(table)
    reference = LeadingColorZJetsNativeEvaluator(model=model)
    validation_points = (
        reference.canonical_neutral_dilepton_gluon_point(
            "d d~ > e+ e- g",
            lepton_pdg=11,
            antilepton_pdg=-11,
            gluon_count=1,
            sqrt_s=1000.0,
        ),
    )

    _write_rusticol_process_artifacts(
        tmp_path,
        evaluator_manifest={
            "compiled": {
                "kind": "shared-compiled-sweep",
                "stages": [],
                "amplitude_stage": {
                    "kind": "amplitude-stage",
                    "output_length": len(table.amplitudes),
                    "raw_sum_weights": [
                        amplitude.multiplicity for amplitude in table.amplitudes
                    ],
                },
            },
            "metadata": {},
        },
        evaluator_manifest_name="manifest.json",
        process="d d~ > e+ e- g",
        gluon_count=1,
        table=table,
        layout=layout,
        model=model,
        validation_points=validation_points,
        vector_pdg=23,
        electroweak_coupling_power=2,
        artifact_family="q-qbar-neutral-dilepton-gluons-leading-color",
    )

    manifest_payload = json.loads((tmp_path / "process_manifest.json").read_text())
    assert manifest_payload["family"] == (
        "q-qbar-neutral-dilepton-gluons-leading-color"
    )
    assert manifest_payload["process"] == "d d~ > e+ e- g"
    assert manifest_payload["external_pdg_order"] == [1, -1, 21, 11, -11]
    assert {
        source.get("source_kind", "external")
        for source in manifest_payload["table"]["sources"]
    } == {"external", "lepton_pair_vector"}


def test_load_process_manifest_rejects_wrong_kind(tmp_path: Path) -> None:
    process_dir = tmp_path / "process"
    process_dir.mkdir()
    (process_dir / "process_manifest.json").write_text(
        json.dumps({"kind": "not-a-pyamplicol-process"}),
        encoding="utf-8",
    )

    with pytest.raises(NativeEvaluationError, match="unsupported process artifact kind"):
        load_process_manifest(process_dir)


def test_load_process_manifest_resolves_process_set_default_and_key(tmp_path: Path) -> None:
    root = tmp_path / "process_set"
    first = root / "subprocesses" / "d_dbar_to_z_g"
    second = root / "subprocesses" / "u_ubar_to_z_g"
    first.mkdir(parents=True)
    second.mkdir(parents=True)
    first_payload = _minimal_process_manifest_payload()
    second_payload = _minimal_process_manifest_payload()
    second_payload["process"] = "u u~ > z g"
    (first / "process_manifest.json").write_text(json.dumps(first_payload))
    (second / "process_manifest.json").write_text(json.dumps(second_payload))
    (root / "process_set_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pyamplicol-rusticol-process-set",
                "default_process_key": "d_dbar_to_z_g",
                "processes": [
                    {
                        "key": "d_dbar_to_z_g",
                        "process": "d d~ > z g",
                        "path": "subprocesses/d_dbar_to_z_g",
                    },
                    {
                        "key": "u_ubar_to_z_g",
                        "process": "u u~ > z g",
                        "path": "subprocesses/u_ubar_to_z_g",
                    },
                ],
            }
        )
    )

    assert load_process_manifest(root).process == "d d~ > z g"
    runtime = load_process(root, runtime="python", process_key="u_ubar_to_z_g")
    assert runtime.process == "u u~ > z g"


def test_load_process_manifest_loads_generic_dag_plan_only_artifact(
    tmp_path: Path,
) -> None:
    write_generic_dag_process_artifact("d d~ > z g", tmp_path)

    manifest = load_process_manifest(tmp_path)
    assert manifest.schema_version == 2
    assert manifest.kind == "pyamplicol-generic-dag-process"
    assert manifest.artifact_class == "generic-dag-schema-v2"
    assert manifest.family is None
    assert manifest.gluon_count is None
    assert manifest.process == "d d~ > z g"
    assert manifest.external_pdg_order == (1, -1, 23, 21)

    runtime = load_process(tmp_path, runtime="python")
    metadata = runtime.metadata
    assert metadata["runtime"] == "python-generic-dag-schema-v2"
    assert metadata["artifact_class"] == "generic-dag-schema-v2"
    assert "family" not in metadata
    assert "gluon_count" not in metadata
    assert metadata["schema_version"] == 2
    assert metadata["runtime_available"] is False
    assert metadata["current_count"]
    momenta = [
        [
            [500.0, 0.0, 0.0, 500.0],
            [500.0, 0.0, 0.0, -500.0],
            [500.0, 0.0, 0.0, 0.0],
            [500.0, 0.0, 0.0, 0.0],
        ]
    ]
    diagnostics = runtime.stage_diagnostics(momenta)
    assert diagnostics["runtime"] == "python-generic-dag-schema-v2"
    assert diagnostics["runtime_available"] is False
    assert diagnostics["stages"]
    with pytest.raises(NativeEvaluationError, match="generic DAG evaluator stages"):
        runtime.evaluate(momenta)


def test_load_process_loads_generic_dag_process_set_subprocess(
    tmp_path: Path,
) -> None:
    write_generic_dag_process_set_artifact(
        build_process_set_ir("d d~ > z g | u u~ > z g"),
        tmp_path,
    )

    assert load_process_manifest(tmp_path).process == "d d~ > z g"
    runtime = load_process(tmp_path, runtime="python", process_key="u_ubar_to_z_g")
    assert runtime.process == "u u~ > z g"
    assert runtime.metadata["artifact_class"] == "generic-dag-schema-v2"
    assert "family" not in runtime.metadata
    assert "gluon_count" not in runtime.metadata
    assert runtime.metadata["key"] == "u_ubar_to_z_g"
    with pytest.raises(NativeEvaluationError, match="process 'missing' not found"):
        load_process(tmp_path, runtime="python", process_key="missing")


@pytest.mark.parametrize(
    "process",
    [
        "d d~ > e+ e- g",
        "u d~ > e+ ve g",
        "g g > u u~",
        "d d~ > u u~ s s~",
    ],
)
def test_python_runtime_loads_generic_schema_v2_without_family_whitelist(
    tmp_path: Path,
    process: str,
) -> None:
    process_dir = tmp_path / process.replace(" ", "_").replace(">", "to")
    write_generic_dag_process_artifact(process, process_dir)

    manifest = load_process_manifest(process_dir)
    runtime = load_process(process_dir, runtime="python")
    metadata = runtime.metadata

    assert manifest.process == process
    assert manifest.schema_version == 2
    assert manifest.artifact_class == "generic-dag-schema-v2"
    assert manifest.family is None
    assert manifest.gluon_count is None
    assert metadata["process"] == process
    assert metadata["artifact_class"] == "generic-dag-schema-v2"
    assert "family" not in metadata
    assert "gluon_count" not in metadata
    assert metadata["runtime"] == "python-generic-dag-schema-v2"
    assert metadata["current_count"]
    assert metadata["interaction_count"]
    assert metadata["amplitude_root_count"]


def test_rusticol_load_rejects_incomplete_manifest(tmp_path: Path) -> None:
    _skip_if_rusticol_unavailable()
    process_dir = tmp_path / "process"
    process_dir.mkdir()
    (process_dir / "process_manifest.json").write_text("{}", encoding="utf-8")

    result = _run_rusticol_load_subprocess(process_dir)
    assert result.returncode == 0, result.stderr
    assert "Runtime.load only supports schema-v2 generic DAG artifacts" in result.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_default_load_rejects_valid_schema_v1_reference_artifact(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()
    process_dir = tmp_path / "zero"
    _write_zero_gluon_rusticol_process_artifacts(
        process_dir,
        process="d d~ > z",
        model=AmplicolSMLeadingColorModel(),
    )

    code = f"""
import rusticol

try:
    rusticol.Runtime.load({str(process_dir)!r})
except ValueError as exc:
    assert "schema-v2 generic DAG artifacts" in str(exc)
else:
    raise AssertionError("default Runtime.load accepted a legacy schema-v1 artifact")

legacy = rusticol.Runtime.load_legacy({str(process_dir)!r})
assert legacy.process == "d d~ > z"
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_load_rejects_wrong_kind_and_schema(tmp_path: Path) -> None:
    _skip_if_rusticol_unavailable()
    for name, patch, expected in (
        (
            "wrong-kind",
            {"kind": "not-a-pyamplicol-process", "schema_version": 1},
            "unsupported process artifact kind",
        ),
        (
            "wrong-schema",
            {"kind": "pyamplicol-rusticol-process", "schema_version": 999},
            "unsupported process manifest schema_version",
        ),
    ):
        process_dir = tmp_path / name
        process_dir.mkdir()
        payload = _minimal_process_manifest_payload()
        payload.update(patch)
        (process_dir / "process_manifest.json").write_text(
            json.dumps(payload),
            encoding="utf-8",
        )

        result = _run_rusticol_load_subprocess(process_dir, legacy=True)
        assert result.returncode == 0, result.stderr
        assert expected in result.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_load_rejects_unsupported_family(tmp_path: Path) -> None:
    _skip_if_rusticol_unavailable()
    process_dir = tmp_path / "process"
    process_dir.mkdir()
    payload = _minimal_process_manifest_payload()
    payload.update(
        {
            "process": "u d~ > w+ g g",
            "family": "unsupported-leading-color-family",
            "gluon_count": 2,
            "external_pdg_order": [2, -1, 24, 21, 21],
        }
    )
    (process_dir / "process_manifest.json").write_text(
        json.dumps(payload),
        encoding="utf-8",
    )

    result = _run_rusticol_load_subprocess(process_dir, legacy=True)
    assert result.returncode == 0, result.stderr
    assert "currently supports q-qbar vector/leptonic-vector+gluon" in result.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_w_gluon_process_artifact_loads_in_rusticol(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()

    process = "u d~ > w+ g"
    process_dir = tmp_path / "w_process"
    evaluator = ZGluonDAGEvaluator(
        process,
        symbolica_evaluator_backend="jit",
        symbolica_n_cores=1,
    )
    manifest_path = evaluator.save_evaluator_artifact(process_dir)
    evaluator_manifest = json.loads(manifest_path.read_text())
    manifest_payload = json.loads(
        (process_dir / "process_manifest.json").read_text()
    )

    assert evaluator_manifest["kind"] == "pyamplicol-vector-gluon-shared-dag-compiled"
    assert manifest_payload["family"] == "q-qbar-vector-gluons-leading-color"
    assert manifest_payload["external_pdg_order"] == [2, -1, 21, 24]
    assert manifest_payload["model"]["mass_w"] == evaluator.model.mass(24)

    reference = LeadingColorZJetsNativeEvaluator(evaluator.model).evaluate(
        process,
        particles=LeadingColorZJetsNativeEvaluator(
            evaluator.model
        ).canonical_neutral_vector_gluon_point(
            process,
            vector_pdg=24,
            gluon_count=1,
            sqrt_s=1000.0,
        ),
    ).matrix_element

    code = f"""
import json
from pathlib import Path
import numpy as np
import rusticol

root = Path({str(process_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [
        [float(component) for component in particle["momentum"]]
        for particle in point
    ]
    for point in payload["points"]
], dtype=np.float64)
runtime = rusticol.Runtime.load_legacy(str(root))
print(json.dumps({{"values": list(runtime.evaluate(momenta))}}))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    values = json.loads(result.stdout)["values"]
    assert abs(float(values[0]) - reference) / abs(reference) < 1.0e-12


@_LEGACY_SCHEMA_V1_TEST
def test_reversed_beam_vector_process_artifact_loads_in_rusticol(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()

    process = "d~ d > z g"
    process_dir = tmp_path / "reversed_z_process"
    evaluator = ZGluonDAGEvaluator(
        process,
        symbolica_evaluator_backend="jit",
        symbolica_n_cores=1,
    )
    evaluator.save_evaluator_artifact(process_dir)
    manifest_payload = json.loads(
        (process_dir / "process_manifest.json").read_text()
    )

    assert manifest_payload["external_pdg_order"] == [-1, 1, 21, 23]
    currents = manifest_payload["table"]["currents"]
    assert [
        (source["leg_label"], currents[source["current_id"]]["pdg"])
        for source in manifest_payload["table"]["sources"][:4]
    ] == [(2, -1), (2, -1), (1, 1), (1, 1)]

    native = LeadingColorZJetsNativeEvaluator(evaluator.model)
    reference = native.evaluate(
        process,
        particles=native.canonical_neutral_vector_gluon_point(
            process,
            vector_pdg=23,
            gluon_count=1,
            sqrt_s=1000.0,
        ),
    ).matrix_element

    code = f"""
import json
from pathlib import Path
import numpy as np
import rusticol

root = Path({str(process_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [
        [float(component) for component in particle["momentum"]]
        for particle in point
    ]
    for point in payload["points"]
], dtype=np.float64)
runtime = rusticol.Runtime.load_legacy(str(root))
print(json.dumps({{"values": list(runtime.evaluate(momenta))}}))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    values = json.loads(result.stdout)["values"]
    assert abs(float(values[0]) - reference) / abs(reference) < 1.0e-12


@_LEGACY_SCHEMA_V1_TEST
def test_zero_gluon_process_artifact_loads_in_python_and_rusticol(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()
    model = AmplicolSMLeadingColorModel()
    process_dir = tmp_path / "zero"
    manifest_path = _write_zero_gluon_rusticol_process_artifacts(
        process_dir,
        process="d d~ > z",
        model=model,
    )
    manifest = load_process_manifest(process_dir)
    assert manifest_path.name == "process_manifest.json"
    assert manifest.gluon_count == 0
    assert manifest.compiled_kind == "zero-gluon-symbolic-scalar"
    manifest_payload = json.loads(manifest_path.read_text())
    assert manifest_payload["backend_settings"]["runtime"] == "rusticol"
    assert (
        manifest_payload["backend_settings"]["compiled_kind"]
        == "zero-gluon-symbolic-scalar"
    )
    particles_by_pdg = {
        particle["pdg"]: particle
        for particle in manifest_payload["model"]["particles"]
    }
    assert {1, 23}.issubset(particles_by_pdg)
    assert particles_by_pdg[1]["anti_pdg"] == -1
    assert manifest_payload["model"]["vertices"]
    assert "dependency_fingerprint" in manifest_payload
    assert "symbolica" in manifest_payload["dependency_fingerprint"]
    assert manifest_payload["helicity_filter"]["retained_amplitude_count"] == 1
    assert manifest_payload["helicity_filter"]["raw_sum_weights"] == [1.0]
    assert [particle["pdg"] for particle in manifest_payload["external_particles"]] == [
        1,
        -1,
        23,
    ]
    assert [particle["role"] for particle in manifest_payload["external_particles"]] == [
        "initial",
        "initial",
        "final",
    ]
    assert manifest_payload["momentum_conventions"]["input_shape"] == [
        "batch",
        3,
        4,
    ]
    assert manifest_payload["momentum_conventions"]["incoming_labels"] == [1, 2]
    assert manifest_payload["momentum_conventions"]["final_state_labels"] == [3]

    payload = json.loads((process_dir / "validation_momenta.json").read_text())
    momenta = [
        [
            [float(component) for component in particle["momentum"]]
            for particle in point
        ]
        for point in payload["points"]
    ]
    native = LeadingColorZJetsNativeEvaluator(model)
    reference = native.evaluate(
        "d d~ > z",
        particles=native.canonical_zero_gluon_point("d d~ > z", sqrt_s=1000.0),
    ).matrix_element
    python_runtime = load_process(process_dir, runtime="python")
    python_value = python_runtime.evaluate(momenta)[0]
    assert abs(python_value - reference) / abs(reference) < 1.0e-15
    assert python_runtime.metadata["runtime"] == "python-zero-gluon-symbolic-scalar"
    python_profile = python_runtime.profile(momenta)
    assert python_profile["points"] == 1
    assert python_profile["values"] == [python_value]
    assert python_profile["timing"] is None
    python_diagnostics = python_runtime.stage_diagnostics(momenta)
    assert python_diagnostics["stages"][0]["stage"] == "zero_gluon_symbolic_scalar"
    with pytest.raises(NativeEvaluationError, match=r"\(batch, 3, 4\)"):
        python_runtime.evaluate([momenta[0][:2]])

    code = f"""
import json
from decimal import Decimal
from pathlib import Path
import numpy as np
import pyamplicol
import rusticol

root = Path({str(process_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)
decimal_momenta = [
    [[Decimal(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
]
runtime = rusticol.Runtime.load_legacy(str(root))
value = float(runtime.evaluate(momenta)[0])
value16 = float(runtime.evaluate_with_prec(momenta, 16)[0])
value32 = runtime.evaluate_with_prec(decimal_momenta, 32)[0]
value50 = runtime.evaluate_with_prec(decimal_momenta, 50)[0]
profile = dict(runtime.profile(momenta, precision=16))
batched_momenta = np.repeat(momenta, 3, axis=0)
batched_values = runtime.evaluate(batched_momenta)
batched_profile = dict(runtime.profile(batched_momenta, precision=16))
assert abs(value - {reference!r}) / abs({reference!r}) < 1.0e-14
assert value == value16
assert value32.__class__.__name__ == "Decimal"
assert value50.__class__.__name__ == "Decimal"
assert profile["points"] == 1
assert profile["core_evaluator_time_s"] >= 0.0
assert len(batched_values) == 3
assert all(abs(float(item) - value) / abs(value) < 1.0e-14 for item in batched_values)
assert batched_profile["points"] == 3
assert batched_profile["batch_size"] == 3

bad_momenta = momenta[:, :2, :]
bad_decimal_momenta = [point[:2] for point in decimal_momenta]
for call in (
    lambda: runtime.evaluate(bad_momenta),
    lambda: runtime.evaluate_with_prec(bad_momenta, 16),
    lambda: runtime.profile(bad_momenta, precision=16),
    lambda: runtime.evaluate_with_prec(bad_decimal_momenta, 32),
    lambda: runtime.evaluate_with_prec(bad_decimal_momenta, 50),
):
    try:
        call()
    except ValueError as exc:
        assert "expected 3" in str(exc) or "(batch, 3, 4)" in str(exc)
    else:
        raise AssertionError("rusticol accepted a momentum batch with the wrong leg count")
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout

    rusticol_spec = importlib.util.find_spec("rusticol")
    assert rusticol_spec is not None and rusticol_spec.origin is not None
    rusticol_folder = Path(rusticol_spec.origin).parent
    for args in (
        ("--precision", "16"),
        ("--precision", "32"),
        ("--precision", "50"),
        ("--precision", "16", "--profile", "--target-runtime", "0.001"),
        ("--precision", "16", "--rusticol-folder", str(rusticol_folder)),
    ):
        standalone = subprocess.run(
            [sys.executable, str(process_dir / "check_standalone.py"), *args],
            check=False,
            env=env,
            text=True,
            capture_output=True,
        )
        assert standalone.returncode == 0, standalone.stderr + standalone.stdout
        assert "d d~ > z" in standalone.stdout
        assert "family" in standalone.stdout
        assert "gluons" in standalone.stdout
        assert "backend" in standalone.stdout
        assert "retained amplitudes" in standalone.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_rusticol_process_set_root_loads_default_and_selected_subprocess(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()
    model = AmplicolSMLeadingColorModel()
    root = tmp_path / "process_set"
    first = root / "subprocesses" / "d_dbar_to_z"
    second = root / "subprocesses" / "u_ubar_to_z"
    _write_zero_gluon_rusticol_process_artifacts(
        first,
        process="d d~ > z",
        model=model,
    )
    _write_zero_gluon_rusticol_process_artifacts(
        second,
        process="u u~ > z",
        model=model,
    )
    (root / "process_set_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 1,
                "kind": "pyamplicol-rusticol-process-set",
                "default_process_key": "d_dbar_to_z",
                "processes": [
                    {
                        "key": "d_dbar_to_z",
                        "process": "d d~ > z",
                        "path": "subprocesses/d_dbar_to_z",
                    },
                    {
                        "key": "u_ubar_to_z",
                        "process": "u u~ > z",
                        "path": "subprocesses/u_ubar_to_z",
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    code = f"""
import json
from pathlib import Path
import numpy as np
import rusticol

root = Path({str(root)!r})
default_runtime = rusticol.Runtime.load_legacy(str(root))
selected_by_key = rusticol.Runtime.load_legacy(str(root), process_key="u_ubar_to_z")
selected_by_process = rusticol.Runtime.load_legacy(str(root), "u u~ > z")
assert default_runtime.process == "d d~ > z"
assert selected_by_key.process == "u u~ > z"
assert selected_by_process.process == "u u~ > z"

def validation_points(process_dir):
    payload = json.loads((process_dir / "validation_momenta.json").read_text())
    return np.asarray([
        [[float(component) for component in particle["momentum"]] for particle in point]
        for point in payload["points"]
    ], dtype=np.float64)

default_values = default_runtime.evaluate(
    validation_points(root / "subprocesses" / "d_dbar_to_z")
)
selected_values = selected_by_key.evaluate(
    validation_points(root / "subprocesses" / "u_ubar_to_z")
)
assert len(default_values) == 1
assert len(selected_values) == 1
try:
    rusticol.Runtime.load_legacy(str(root), process_key="missing")
except ValueError as exc:
    assert "not found" in str(exc)
    assert "d_dbar_to_z" in str(exc)
    assert "u_ubar_to_z" in str(exc)
else:
    raise AssertionError("rusticol accepted an unknown process-set key")

try:
    rusticol.Runtime.load_legacy(str(root / "subprocesses" / "d_dbar_to_z"), "u_ubar_to_z")
except ValueError as exc:
    assert "process_key can only be used with a process-set artifact" in str(exc)
else:
    raise AssertionError("rusticol accepted process_key for a single-process artifact")

print(json.dumps({{
    "default": float(default_values[0]),
    "selected": float(selected_values[0]),
}}))
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_python_process_runtime_loads_nonzero_manifest_without_regenerating_graph(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _skip_if_rusticol_unavailable()
    model = AmplicolSMLeadingColorModel()
    process = "d d~ > z g"
    process_dir = tmp_path / "n1"
    stale_root_evaluator = process_dir / "shared_stage_stale.evaluator.bin"
    process_dir.mkdir()
    stale_root_evaluator.write_bytes(b"stale")
    generator = ZGluonDAGEvaluator(
        process,
        model=model,
        batch_size=16,
        symbolica_evaluator_backend="jit",
    )
    generator.save_evaluator_artifact(process_dir)
    assert list(process_dir.glob("*.evaluator.bin")) == []
    evaluator_files = sorted((process_dir / "evaluators").glob("*.evaluator.bin"))
    assert evaluator_files
    manifest_payload = json.loads((process_dir / "process_manifest.json").read_text())
    assert manifest_payload["backend_settings"]["runtime"] == "rusticol"
    assert manifest_payload["backend_settings"]["compiled_kind"] == "shared-compiled-sweep"
    assert (
        manifest_payload["backend_settings"]["evaluator_manifest"]
        == "manifest.json"
    )
    assert "symbolica_evaluator_settings" in manifest_payload["backend_settings"]
    assert "dependency_fingerprint" in manifest_payload
    assert "install_manifest" in manifest_payload["dependency_fingerprint"]
    particles_by_pdg = {
        particle["pdg"]: particle
        for particle in manifest_payload["model"]["particles"]
    }
    assert {1, 21, 23}.issubset(particles_by_pdg)
    assert particles_by_pdg[1]["anti_pdg"] == -1
    assert manifest_payload["model"]["vertices"]
    assert [particle["label"] for particle in manifest_payload["external_particles"]] == [
        1,
        2,
        3,
        4,
    ]
    assert [particle["role"] for particle in manifest_payload["external_particles"]] == [
        "initial",
        "initial",
        "final",
        "final",
    ]
    assert manifest_payload["momentum_conventions"]["input_shape"] == [
        "batch",
        4,
        4,
    ]
    assert manifest_payload["momentum_conventions"]["source_current_crossing"][
        "crossed_incoming_labels"
    ] == [1, 2]
    helicity_filter = manifest_payload["helicity_filter"]
    assert helicity_filter["strategy"] == "shared-current-retained-amplitudes"
    assert helicity_filter["retained_amplitude_count"] == len(
        generator.table.amplitudes
    )
    assert len(helicity_filter["amplitudes"]) == len(generator.table.amplitudes)
    point = LeadingColorZJetsNativeEvaluator(model).canonical_z_gluon_point(
        process,
        gluon_count=1,
        sqrt_s=1000.0,
    )
    reference = generator.evaluate(point).matrix_element
    momenta = [
        [
            [float(component) for component in particle.momentum]
            for particle in point
        ]
    ]

    def fail_graph_rebuild(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        raise AssertionError("process runtime rebuilt the graph from the process string")

    monkeypatch.setattr(dag_runtime_module, "_z_gluon_graph", fail_graph_rebuild)
    runtime = load_process(process_dir, runtime="python")
    metadata = runtime.metadata
    assert metadata["loaded_from_process_manifest"] is True
    assert metadata["current_count"] == generator.metadata.shared_current_count
    value = runtime.evaluate(momenta)[0]
    assert abs(value - reference) / abs(reference) < 1.0e-12
    profile = runtime.profile(momenta)
    assert profile["points"] == 1
    assert profile["values"] == [value]
    assert isinstance(profile["timing"], dict)
    diagnostics = runtime.stage_diagnostics(momenta)
    assert diagnostics["points"] == 1
    assert diagnostics["stages"][0]["stage"] == "sources"

    code = f"""
import json
import numpy as np
import pyamplicol
import rusticol
from decimal import Decimal

runtime = rusticol.Runtime.load_legacy({str(process_dir)!r})
momenta = np.asarray({momenta!r}, dtype=np.float64)
decimal_momenta = [
    [[Decimal(str(component)) for component in particle] for particle in point]
    for point in {momenta!r}
]
single = float(runtime.evaluate(momenta)[0])
value32 = runtime.evaluate_with_prec(decimal_momenta, 32)[0]
value50 = runtime.evaluate_with_prec(decimal_momenta, 50)[0]
batch = np.repeat(momenta, 4, axis=0)
values = runtime.evaluate(batch)
profile = dict(runtime.profile(batch, precision=16))
assert abs(single - {reference!r}) / abs({reference!r}) < 1.0e-12
assert value32.__class__.__name__ == "Decimal"
assert value50.__class__.__name__ == "Decimal"
assert abs(float(value32) - {reference!r}) / abs({reference!r}) < 1.0e-12
assert abs(float(value50) - {reference!r}) / abs({reference!r}) < 1.0e-12
assert values.shape == (4,)
assert all(abs(float(value) - single) / abs(single) < 1.0e-12 for value in values)
assert profile["points"] == 4
assert profile["batch_size"] == 4
assert "current_rss_bytes" in profile
assert "peak_rss_bytes" in profile
diagnostics = dict(runtime.stage_diagnostics(momenta))
diagnostics["stages"] = [dict(stage) for stage in diagnostics["stages"]]
print(json.dumps(diagnostics, sort_keys=True))
"""
    batch_env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    batch_result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=batch_env,
        text=True,
        capture_output=True,
    )
    assert batch_result.returncode == 0, batch_result.stderr + batch_result.stdout
    rust_diagnostics = json.loads(batch_result.stdout)
    assert rust_diagnostics["points"] == diagnostics["points"]
    assert len(rust_diagnostics["stages"]) == len(diagnostics["stages"])
    for rust_stage, python_stage in zip(
        rust_diagnostics["stages"],
        diagnostics["stages"],
        strict=True,
    ):
        assert rust_stage["stage"] == python_stage["stage"]
        assert rust_stage["output_len"] == python_stage["output_len"]
        for key in ("sum_re", "sum_im", "sum_abs2", "max_abs"):
            assert abs(float(rust_stage[key]) - float(python_stage[key])) < 1.0e-9

    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    for args in (
        ("--precision", "16"),
        ("--precision", "32"),
        ("--precision", "16", "--profile", "--target-runtime", "0.001"),
    ):
        standalone = subprocess.run(
            [sys.executable, str(process_dir / "check_standalone.py"), *args],
            check=False,
            env=env,
            text=True,
            capture_output=True,
        )
        assert standalone.returncode == 0, standalone.stderr + standalone.stdout
        assert process in standalone.stdout
        assert "family" in standalone.stdout
        assert "retained amplitudes" in standalone.stdout


@_LEGACY_SCHEMA_V1_TEST
def test_zero_gluon_neutral_dilepton_artifact_loads_in_python_and_rusticol(
    tmp_path: Path,
) -> None:
    _skip_if_rusticol_unavailable()
    model = AmplicolSMLeadingColorModel()
    native = LeadingColorZJetsNativeEvaluator(model)
    process = "d d~ > e+ e-"
    point = native.canonical_neutral_dilepton_gluon_point(
        process,
        lepton_pdg=11,
        antilepton_pdg=-11,
        gluon_count=0,
        sqrt_s=1000.0,
    )
    process_dir = tmp_path / "epem_zero"
    evaluator = ZGluonDAGEvaluator(
        process,
        model=model,
        graph=_z_gluon_graph("d d~ > z", model),
        vector_source_specs=_neutral_dilepton_vector_source_specs(
            lepton_pdg=11,
            antilepton_pdg=-11,
            lepton_label=3,
            antilepton_label=4,
            model=model,
        ),
        rusticol_validation_points=(point,),
        electroweak_coupling_power=2,
        rusticol_artifact_family=(
            "q-qbar-neutral-dilepton-zero-gluon-leading-color"
        ),
        symbolica_evaluator_backend="jit",
        symbolica_n_cores=1,
    )
    manifest_path = evaluator.save_evaluator_artifact(process_dir)
    manifest_payload = json.loads((process_dir / "process_manifest.json").read_text())

    assert manifest_path.name == "manifest.json"
    assert manifest_payload["family"] == (
        "q-qbar-neutral-dilepton-zero-gluon-leading-color"
    )
    assert manifest_payload["compiled"]["kind"] == "shared-compiled-sweep"
    assert manifest_payload["gluon_count"] == 0
    assert manifest_payload["external_pdg_order"] == [1, -1, 11, -11]

    payload = json.loads((process_dir / "validation_momenta.json").read_text())
    momenta = [
        [
            [float(component) for component in particle["momentum"]]
            for particle in validation_point
        ]
        for validation_point in payload["points"]
    ]
    reference = native.evaluate(process, particles=point).matrix_element
    python_value = load_process(process_dir, runtime="python").evaluate(momenta)[0]
    assert abs(python_value - reference) / abs(reference) < 1.0e-12

    code = f"""
import json
from pathlib import Path
import numpy as np
import rusticol

root = Path({str(process_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)
runtime = rusticol.Runtime.load_legacy(str(root))
print(json.dumps({{"values": list(runtime.evaluate(momenta))}}))
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    rust_value = float(json.loads(result.stdout)["values"][0])
    assert abs(rust_value - reference) / abs(reference) < 1.0e-12


def test_python_and_rusticol_stage_diagnostics_match_when_artifact_exists() -> None:
    _skip_if_rusticol_unavailable()
    artifact_dir = _local_rusticol_smoke_artifact()
    if not (artifact_dir / "process_manifest.json").exists():
        pytest.skip("local rusticol smoke process artifact is not available")

    code = f"""
import json
from pathlib import Path
import numpy as np
import rusticol
import pyamplicol

root = Path({str(artifact_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)

rust_runtime = rusticol.Runtime.load_legacy(str(root))
python_runtime = pyamplicol.load_process(root, runtime="python")
rust_diag = dict(rust_runtime.stage_diagnostics(momenta))
python_diag = python_runtime.stage_diagnostics(momenta)

assert rust_diag["points"] == python_diag["points"]
assert len(rust_diag["stages"]) == len(python_diag["stages"])
for rust_stage, python_stage in zip(rust_diag["stages"], python_diag["stages"], strict=True):
    rust_stage = dict(rust_stage)
    assert rust_stage["stage"] == python_stage["stage"]
    assert rust_stage["output_len"] == python_stage["output_len"]
    for key in ("sum_re", "sum_im", "sum_abs2", "max_abs"):
        assert abs(float(rust_stage[key]) - float(python_stage[key])) < 1.0e-9
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


def test_rusticol_api_smoke_when_artifact_exists() -> None:
    _skip_if_rusticol_unavailable()
    artifact_dir = _local_rusticol_smoke_artifact()
    if not (artifact_dir / "process_manifest.json").exists():
        pytest.skip("local rusticol smoke process artifact is not available")

    code = f"""
import json
from decimal import Decimal
from pathlib import Path
import numpy as np
import rusticol

root = Path({str(artifact_dir)!r})
payload = json.loads((root / "validation_momenta.json").read_text())
decimal_momenta = [
    [[Decimal(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
]
momenta = np.asarray([
    [[float(component) for component in particle["momentum"]] for particle in point]
    for point in payload["points"]
], dtype=np.float64)

runtime = rusticol.Runtime.load_legacy(str(root))
values_a = runtime.evaluate(momenta)
values_b = runtime.evaluate(momenta)
values_16 = runtime.evaluate_with_prec(momenta, 16)
values_32 = runtime.evaluate_with_prec(decimal_momenta, 32)
profile = dict(runtime.profile(momenta, precision=16))

assert isinstance(values_a, np.ndarray)
assert values_a.dtype == np.float64
assert np.array_equal(values_a, values_b)
assert np.array_equal(values_a, values_16)
assert values_32 and values_32[0].__class__.__name__ == "Decimal"
for key in (
    "points",
    "batch_size",
    "values",
    "source_fill_time_s",
    "momentum_setup_time_s",
    "parameter_pack_time_s",
    "stage_evaluator_time_s",
    "amplitude_evaluator_time_s",
    "result_reduction_time_s",
    "output_transfer_time_s",
    "total_time_s",
    "core_evaluator_time_s",
    "peak_rss_bytes",
):
    assert key in profile
"""
    env = dict(
        os.environ,
        PYTHONPATH=str(Path(__file__).resolve().parents[2] / "src"),
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


def _local_rusticol_smoke_artifact() -> Path:
    return (
        Path(__file__).resolve().parents[3]
        / "pyAmpliCol"
        / "outputs"
        / "rusticol_high_precision_smoke"
        / "n1"
    )


def _skip_if_rusticol_unavailable() -> None:
    if importlib.util.find_spec("rusticol") is None:
        pytest.skip("rusticol extension is not installed")


def _run_rusticol_load_subprocess(
    process_dir: Path,
    *,
    legacy: bool = False,
) -> subprocess.CompletedProcess[str]:
    load_method = "load_legacy" if legacy else "load"
    code = (
        "import rusticol\n"
        "try:\n"
        f"    rusticol.Runtime.{load_method}({str(process_dir)!r})\n"
        "except Exception as exc:\n"
        "    print(type(exc).__name__ + ': ' + str(exc))\n"
        "else:\n"
        "    raise SystemExit('rusticol load unexpectedly succeeded')\n"
    )
    env = dict(os.environ)
    return subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        env=env,
        text=True,
        capture_output=True,
    )


def _minimal_process_manifest_payload() -> dict[str, object]:
    return {
        "schema_version": 1,
        "kind": "pyamplicol-rusticol-process",
        "process": "d d~ > z g",
        "family": "q-qbar-z-gluons-leading-color",
        "gluon_count": 1,
        "external_pdg_order": [1, -1, 23, 21],
        "model": {
            "alpha_s_me_check": 0.13,
            "alpha_ew": 0.0075,
            "mass_z": 91.188,
        },
        "normalization": {
            "color_factor": 1.0,
            "average_factor": 1.0,
            "identical_factor": 1.0,
            "coupling_factor": 1.0,
        },
        "layout": {
            "parameter_count": 0,
            "current_offsets": [],
            "momentum_offsets_and_labels": [],
        },
        "table": {
            "currents": [],
            "sources": [],
            "amplitudes": [],
        },
        "compiled": {
            "kind": "shared-compiled-sweep",
            "stages": [],
            "amplitude_stage": {
                "output_length": 0,
                "raw_sum_weights": [],
                "amplitude_evaluator": None,
                "raw_sum_evaluator": None,
            },
        },
    }
