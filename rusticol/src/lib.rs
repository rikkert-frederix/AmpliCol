#![recursion_limit = "256"]

use pyo3::IntoPyObjectExt;
use pyo3::buffer::PyBuffer;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyList};
use serde::Deserialize;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;
use symbolica::evaluate::JITCompiledEvaluator;
use symbolica::prelude::{
    BatchEvaluator, CompiledComplexEvaluator, Complex, DoubleFloat, EvaluationDomain,
    ExpressionEvaluator, Float, JITCompilationSettings, Rational, Real, RealLike,
};

const MAX_LC_TOPOLOGY_REPLAY_EXPANDED_POINTS: usize = 8192;

#[derive(Clone, Debug, Deserialize)]
struct ProcessManifest {
    schema_version: u32,
    kind: String,
    process: String,
    family: String,
    gluon_count: usize,
    external_pdg_order: Vec<i32>,
    model: ModelManifest,
    normalization: NormalizationManifest,
    layout: LayoutManifest,
    table: TableManifest,
    compiled: CompiledSweepManifest,
}

#[derive(Clone, Debug, Deserialize)]
struct ProcessSetManifest {
    schema_version: u32,
    kind: String,
    default_process_key: String,
    processes: Vec<ProcessSetEntryManifest>,
}

#[derive(Clone, Debug, Deserialize)]
struct ProcessSetEntryManifest {
    key: String,
    process: String,
    path: String,
    #[serde(default)]
    crossing_alias_of: Option<String>,
    #[serde(default)]
    input_crossing_map: Option<Vec<InputCrossingMapEntry>>,
}

#[derive(Clone, Debug, Deserialize)]
struct InputCrossingMapEntry {
    target_index: usize,
    source_index: usize,
    sign: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericProcessManifestV2 {
    schema_version: u32,
    kind: String,
    process: String,
    key: String,
    color_accuracy: String,
    external_pdg_order: Vec<i32>,
    #[serde(default)]
    lc_topology_reuse: Option<Value>,
    compiled: GenericCompiledManifestV2,
    dag_summary: GenericDagSummaryManifestV2,
    runtime_schema: GenericRuntimeSchemaManifestV2,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericCompiledManifestV2 {
    kind: String,
    runtime_available: bool,
    runtime_unavailable_message: Option<String>,
    #[serde(default)]
    lc_topology_replay: Option<LcTopologyReplayManifestV2>,
    stage_evaluators: Option<GenericStageEvaluatorArtifactsManifestV2>,
}

#[derive(Clone, Debug, Default, Deserialize)]
struct LcTopologyReplayManifestV2 {
    #[serde(default)]
    enabled: bool,
    #[serde(default)]
    mode: String,
    #[serde(default)]
    replayed_sector_count: usize,
    #[serde(default)]
    groups: Vec<LcTopologyReplayGroupManifestV2>,
}

#[derive(Clone, Debug, Default, Deserialize)]
struct LcTopologyReplayGroupManifestV2 {
    representative_sector_id: i64,
    materialized_sector_id: i64,
    #[serde(default)]
    active_sector_ids: Vec<i64>,
    #[serde(default)]
    sector_permutations: Vec<LcTopologyReplaySectorPermutationManifestV2>,
}

#[derive(Clone, Debug, Default, Deserialize)]
struct LcTopologyReplaySectorPermutationManifestV2 {
    sector_id: i64,
    #[serde(default)]
    label_permutation: Vec<LcTopologyReplayLabelPermutationManifestV2>,
}

#[derive(Clone, Debug, Default, Deserialize)]
struct LcTopologyReplayLabelPermutationManifestV2 {
    representative_label: usize,
    sector_label: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericStageEvaluatorArtifactsManifestV2 {
    kind: String,
    runtime_available: bool,
    runtime_unavailable_message: Option<String>,
    parameter_count: usize,
    value_parameter_count: usize,
    momentum_parameter_count: usize,
    real_valued_inputs: Vec<usize>,
    parameter_layout: String,
    stage_count: usize,
    stages: Vec<GenericSerializedStageEvaluatorManifestV2>,
    amplitude_stage: GenericSerializedStageEvaluatorManifestV2,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericSerializedStageEvaluatorManifestV2 {
    stage_index: usize,
    stage_kind: String,
    subset_size: Option<usize>,
    evaluator_label: String,
    parameter_layout: String,
    output_length: usize,
    output_slots: Vec<GenericStageOutputSlotManifestV2>,
    input_value_slot_ids: Vec<usize>,
    output_value_slot_ids: Vec<usize>,
    interaction_ids: Vec<usize>,
    #[serde(default)]
    input_components: Vec<GenericStageInputComponentManifestV2>,
    #[serde(default)]
    parameter_count: usize,
    #[serde(default)]
    value_parameter_count: usize,
    #[serde(default)]
    momentum_parameter_count: usize,
    #[serde(default)]
    real_valued_inputs: Vec<usize>,
    expression_ready: bool,
    blockers: Vec<String>,
    evaluator: EvaluatorManifest,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericStageInputComponentManifestV2 {
    kind: String,
    source_id: usize,
    component: usize,
    global_component: usize,
    parameter_index: usize,
    #[serde(default)]
    real_valued: bool,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericStageOutputSlotManifestV2 {
    value_slot_id: isize,
    current_id: isize,
    variant: String,
    component_start: usize,
    component_stop: usize,
    output_start: usize,
    output_stop: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericDagSummaryManifestV2 {
    current_count: usize,
    source_count: usize,
    interaction_count: usize,
    amplitude_root_count: usize,
    truncated: bool,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericRuntimeSchemaManifestV2 {
    schema_version: u32,
    kind: String,
    process_key: String,
    process: String,
    external_particles: Vec<GenericExternalParticleManifestV2>,
    #[serde(default)]
    model: Option<GenericRuntimeModelManifestV2>,
    #[serde(default)]
    normalization: Option<GenericRuntimeNormalizationManifestV2>,
    parameter_layout: GenericParameterLayoutManifestV2,
    current_storage: GenericCurrentStorageManifestV2,
    value_storage: GenericValueStorageManifestV2,
    source_fill: GenericSourceFillManifestV2,
    momentum_slots: Vec<GenericMomentumSlotManifestV2>,
    stages: Vec<GenericStageManifestV2>,
    amplitude_stage: GenericAmplitudeStageManifestV2,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericExternalParticleManifestV2 {
    label: usize,
    index: usize,
    pdg: i32,
    outgoing_pdg: i32,
    role: String,
    momentum_slot: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericRuntimeModelManifestV2 {
    #[serde(default)]
    particles: Vec<GenericRuntimeParticleManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericRuntimeParticleManifestV2 {
    pdg: i32,
    #[serde(default)]
    mass: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericRuntimeNormalizationManifestV2 {
    #[serde(default = "default_one_f64")]
    color_factor: f64,
    #[serde(default = "default_one_f64")]
    global_coupling_factor: f64,
    #[serde(default = "default_one_f64")]
    average_factor: f64,
    #[serde(default = "default_one_f64")]
    identical_factor: f64,
}

fn default_one_f64() -> f64 {
    1.0
}

#[derive(Clone, Debug, Deserialize)]
struct GenericParameterLayoutManifestV2 {
    source_component_parameter_count: usize,
    momentum_parameter_count: usize,
    parameter_count_if_flattened: usize,
    value_component_count: usize,
    source_components_complex: bool,
    momentum_components_real: bool,
    real_valued_inputs: Vec<usize>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericCurrentStorageManifestV2 {
    component_count: usize,
    number_type: String,
    current_slots: Vec<GenericCurrentSlotManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericCurrentSlotManifestV2 {
    current_id: usize,
    component_start: usize,
    component_stop: usize,
    dimension: usize,
    is_source: bool,
    particle_id: i32,
    external_mask: u64,
    external_labels: Vec<usize>,
    helicity_ancestry: Value,
    chirality: i32,
    spin_state: Value,
    flavour_flow: Vec<i32>,
    charge_flow: i32,
    color_state: Value,
    momentum_mask: u64,
    auxiliary_kind: Option<String>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericValueStorageManifestV2 {
    component_count: usize,
    number_type: String,
    value_slots: Vec<GenericValueSlotManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericValueSlotManifestV2 {
    value_slot_id: usize,
    current_id: usize,
    variant: String,
    component_start: usize,
    component_stop: usize,
    dimension: usize,
    current_component_start: usize,
    current_component_stop: usize,
    is_source: bool,
    applies_propagator: bool,
    particle_id: i32,
    external_mask: u64,
    external_labels: Vec<usize>,
    momentum_mask: u64,
    chirality: i32,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericSourceFillManifestV2 {
    source_count: usize,
    sources: Vec<GenericSourceRecordManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericSourceRecordManifestV2 {
    source_id: usize,
    current_id: usize,
    current_component_start: usize,
    current_component_stop: usize,
    value_slot: GenericValueSlotRefManifestV2,
    source_parameter_start: usize,
    source_parameter_stop: usize,
    leg_label: usize,
    input_momentum_slot: usize,
    side: String,
    crossing: String,
    physical_pdg: i32,
    outgoing_pdg: i32,
    particle_id: i32,
    source_kind: String,
    source_helicity: i32,
    chirality: i32,
    spin_state: Value,
    dimension: usize,
    helicity_ancestry: Value,
    color_state: Value,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericMomentumSlotManifestV2 {
    momentum_slot_id: usize,
    momentum_mask: u64,
    external_labels: Vec<usize>,
    component_start: usize,
    component_stop: usize,
    real_valued: bool,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericStageManifestV2 {
    stage_index: usize,
    stage_kind: String,
    subset_size: usize,
    input_current_ids: Vec<usize>,
    output_current_ids: Vec<usize>,
    input_value_slot_ids: Vec<usize>,
    output_value_slot_ids: Vec<usize>,
    interaction_count: usize,
    interactions: Vec<GenericInteractionManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericInteractionManifestV2 {
    interaction_id: usize,
    vertex_kind: i32,
    vertex_particles: Vec<i32>,
    left_current_id: usize,
    right_current_id: usize,
    result_current_id: usize,
    left_slot: GenericSlotRefManifestV2,
    right_slot: GenericSlotRefManifestV2,
    result_slot: GenericSlotRefManifestV2,
    left_value_slot: GenericValueSlotRefManifestV2,
    right_value_slot: GenericValueSlotRefManifestV2,
    result_value_slots: Vec<GenericValueSlotRefManifestV2>,
    result_requires_propagated_value: bool,
    result_requires_unpropagated_value: bool,
    momentum_slots: GenericInteractionMomentumSlotsManifestV2,
    coupling: Vec<f64>,
    color_weight: Vec<f64>,
    accumulation: String,
    lowering: GenericLoweringManifestV2,
    full_tensor_network_ready: bool,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericLoweringManifestV2 {
    kind: i32,
    backend: String,
    tensor_names: Vec<String>,
    expression_head: String,
    full_tensor_network_ready: bool,
    description: String,
    kernel: String,
    input_roles: Vec<String>,
    output_role: String,
    coupling_mode: String,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericInteractionMomentumSlotsManifestV2 {
    left: usize,
    right: usize,
    result: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericSlotRefManifestV2 {
    current_id: usize,
    component_start: usize,
    component_stop: usize,
    dimension: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericAmplitudeStageManifestV2 {
    stage_kind: String,
    output_count: usize,
    #[serde(default)]
    color_contraction: Option<GenericColorContractionManifestV2>,
    roots: Vec<GenericAmplitudeRootManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericColorContractionManifestV2 {
    supported: bool,
    #[serde(default)]
    reason: Option<String>,
    group_count: usize,
    #[serde(default)]
    includes_color_factor: bool,
    entries: Vec<GenericColorContractionEntryManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericColorContractionEntryManifestV2 {
    left_group_id: i64,
    right_group_id: i64,
    weight: Vec<f64>,
    #[serde(default = "default_symmetry_factor")]
    symmetry_factor: f64,
}

fn default_symmetry_factor() -> f64 {
    1.0
}

#[derive(Clone, Debug, Deserialize)]
struct GenericAmplitudeRootManifestV2 {
    output_index: usize,
    root_id: usize,
    kind: String,
    left_current_id: usize,
    right_current_id: usize,
    left_slot: GenericSlotRefManifestV2,
    right_slot: GenericSlotRefManifestV2,
    left_value_slot: GenericValueSlotRefManifestV2,
    right_value_slot: GenericValueSlotRefManifestV2,
    vertex_kind: Option<i32>,
    vertex_particles: Option<Vec<i32>>,
    coupling: Vec<f64>,
    color_weight: Vec<f64>,
    contraction: String,
    coherent_group_id: Option<Value>,
    helicity_weight: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericValueSlotRefManifestV2 {
    value_slot_id: usize,
    current_id: usize,
    variant: String,
    component_start: usize,
    component_stop: usize,
    dimension: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct ModelManifest {
    alpha_s_me_check: f64,
    alpha_ew: f64,
    mass_z: f64,
    mass_w: Option<f64>,
    width_z: Option<f64>,
    width_w: Option<f64>,
}

#[derive(Clone, Debug, Deserialize)]
struct NormalizationManifest {
    color_factor: f64,
    average_factor: f64,
    identical_factor: f64,
    coupling_factor: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct LayoutManifest {
    parameter_count: usize,
    current_offsets: Vec<usize>,
    momentum_offsets_and_labels: Vec<MomentumOffsetManifest>,
}

#[derive(Clone, Debug, Deserialize)]
struct MomentumOffsetManifest {
    offset: usize,
    labels: Vec<usize>,
}

#[derive(Clone, Debug, Deserialize)]
struct TableManifest {
    currents: Vec<CurrentManifest>,
    sources: Vec<SourceManifest>,
    amplitudes: Vec<AmplitudeManifest>,
}

#[derive(Clone, Debug, Deserialize)]
struct CurrentManifest {
    id: usize,
    pdg: i32,
    dimension: usize,
}

#[derive(Clone, Debug, Deserialize)]
struct SourceManifest {
    current_id: usize,
    leg_label: usize,
    helicity: i32,
    physical_helicity: i32,
    chirality: i32,
    #[serde(default = "default_source_kind")]
    source_kind: String,
    #[serde(default)]
    partner_leg_label: Option<usize>,
    #[serde(default)]
    partner_helicity: Option<i32>,
    #[serde(default)]
    partner_chirality: Option<i32>,
    #[serde(default)]
    vector_pdg: Option<i32>,
    #[serde(default)]
    coupling: Option<[f64; 2]>,
}

#[derive(Clone, Debug, Deserialize)]
struct AmplitudeManifest {
    multiplicity: f64,
    #[serde(default)]
    coherent_group_id: Option<i64>,
}

fn default_source_kind() -> String {
    "external".to_string()
}

#[derive(Clone, Debug, Deserialize)]
struct CompiledSweepManifest {
    kind: String,
    stages: Vec<CurrentStageManifest>,
    amplitude_stage: AmplitudeStageManifest,
    zero_gluon: Option<ZeroGluonManifest>,
}

#[derive(Clone, Debug, Deserialize)]
struct ZeroGluonManifest {
    parameter_names: Vec<String>,
    evaluator_state_path: String,
    z_left: f64,
    z_right: f64,
}

#[derive(Clone, Debug, Deserialize)]
struct CurrentStageManifest {
    output_slots: Vec<OutputSlotManifest>,
    evaluator: EvaluatorManifest,
}

#[derive(Clone, Debug, Deserialize)]
struct AmplitudeStageManifest {
    output_length: usize,
    raw_sum_weights: Vec<f64>,
    #[serde(default)]
    raw_sum_group_ids: Option<Vec<Option<i64>>>,
    amplitude_evaluator: Option<EvaluatorManifest>,
    raw_sum_evaluator: Option<EvaluatorManifest>,
}

fn resolve_process_root(root: &Path, process_key: Option<&str>) -> PyResult<ProcessSetSelection> {
    if root.join("process_manifest.json").exists() {
        if process_key.is_some() {
            return Err(PyValueError::new_err(
                "process_key can only be used with a process-set artifact",
            ));
        }
        return Ok(ProcessSetSelection {
            root: root.to_path_buf(),
            selected_key: None,
            selected_process: None,
            input_crossing_map: None,
            crossing_alias_of: None,
        });
    }
    let manifest_path = root.join("process_set_manifest.json");
    if !manifest_path.exists() {
        return Err(PyValueError::new_err(format!(
            "no process artifact manifest found in {}",
            root.display()
        )));
    }
    let manifest_text = fs::read_to_string(&manifest_path).map_err(|err| {
        PyValueError::new_err(format!(
            "could not read process-set manifest {}: {err}",
            manifest_path.display()
        ))
    })?;
    let manifest: ProcessSetManifest = serde_json::from_str(&manifest_text).map_err(|err| {
        PyValueError::new_err(format!(
            "could not parse process-set manifest {}: {err}",
            manifest_path.display()
        ))
    })?;
    let supported_process_set = (manifest.schema_version == 1
        && manifest.kind == "pyamplicol-rusticol-process-set")
        || (manifest.schema_version == 2 && manifest.kind == "pyamplicol-generic-dag-process-set");
    if !supported_process_set {
        return Err(PyValueError::new_err(format!(
            "unsupported process-set artifact kind {} schema {}",
            manifest.kind, manifest.schema_version
        )));
    }
    let selected = process_key.unwrap_or(&manifest.default_process_key);
    for entry in &manifest.processes {
        if selected == entry.key || selected == entry.process {
            let path = PathBuf::from(&entry.path);
            let selected_root = if path.is_absolute() {
                path
            } else {
                root.join(path)
            };
            return Ok(ProcessSetSelection {
                root: selected_root,
                selected_key: Some(entry.key.clone()),
                selected_process: Some(entry.process.clone()),
                input_crossing_map: entry.input_crossing_map.clone(),
                crossing_alias_of: entry.crossing_alias_of.clone(),
            });
        }
    }
    let available = manifest
        .processes
        .iter()
        .map(|entry| entry.key.as_str())
        .collect::<Vec<_>>()
        .join(", ");
    Err(PyValueError::new_err(format!(
        "process {selected:?} not found in {}; available: {available}",
        root.display()
    )))
}

fn validate_source_metadata(sources: &[SourceManifest]) -> PyResult<()> {
    for source in sources {
        match source.source_kind.as_str() {
            "external" => {}
            "lepton_pair_vector" => {
                if source.partner_leg_label.is_none()
                    || source.partner_helicity.is_none()
                    || source.partner_chirality.is_none()
                    || source.vector_pdg.is_none()
                    || source.coupling.is_none()
                {
                    return Err(PyValueError::new_err(
                        "lepton-pair vector source metadata is incomplete",
                    ));
                }
            }
            other => {
                return Err(PyValueError::new_err(format!(
                    "unsupported source kind {other:?}"
                )));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
fn parse_generic_schema_v2_manifest(manifest: &Value) -> PyResult<Option<GenericRuntimeV2>> {
    let schema_version = manifest
        .get("schema_version")
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let kind = manifest.get("kind").and_then(Value::as_str).unwrap_or("");
    if schema_version == 2 && kind == "pyamplicol-generic-dag-process" {
        let generic: GenericProcessManifestV2 =
            serde_json::from_value(manifest.clone()).map_err(|err| {
                PyValueError::new_err(format!(
                    "could not parse generic DAG process manifest schema v2: {err}"
                ))
            })?;
        validate_generic_schema_v2_manifest(&generic)?;
        return Ok(Some(GenericRuntimeV2::from_manifest(generic)?));
    }
    Ok(None)
}

fn load_generic_schema_v2_manifest(
    manifest: &Value,
    root: &Path,
) -> PyResult<Option<GenericRuntimeV2>> {
    let schema_version = manifest
        .get("schema_version")
        .and_then(Value::as_u64)
        .unwrap_or(0);
    let kind = manifest.get("kind").and_then(Value::as_str).unwrap_or("");
    if schema_version == 2 && kind == "pyamplicol-generic-dag-process" {
        let generic: GenericProcessManifestV2 =
            serde_json::from_value(manifest.clone()).map_err(|err| {
                PyValueError::new_err(format!(
                    "could not parse generic DAG process manifest schema v2: {err}"
                ))
            })?;
        validate_generic_schema_v2_manifest(&generic)?;
        return Ok(Some(GenericRuntimeV2::load_from_manifest(generic, root)?));
    }
    Ok(None)
}

fn validate_generic_schema_v2_manifest(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    if manifest.schema_version != 2 || manifest.kind != "pyamplicol-generic-dag-process" {
        return Err(PyValueError::new_err(format!(
            "unsupported generic DAG process manifest kind {} schema {}",
            manifest.kind, manifest.schema_version
        )));
    }
    if manifest.runtime_schema.schema_version != 2
        || manifest.runtime_schema.kind != "pyamplicol-generic-dag-runtime-schema"
    {
        return Err(PyValueError::new_err(format!(
            "unsupported generic runtime schema kind {} schema {}",
            manifest.runtime_schema.kind, manifest.runtime_schema.schema_version
        )));
    }
    if manifest.runtime_schema.process != manifest.process
        || manifest.runtime_schema.process_key != manifest.key
    {
        return Err(PyValueError::new_err(
            "generic runtime schema process identity does not match manifest",
        ));
    }
    if !matches!(
        manifest.compiled.kind.as_str(),
        "generic-dag-plan-only" | "generic-dag-stage-blueprint"
    ) {
        return Err(PyValueError::new_err(
            "unsupported generic DAG schema-v2 compiled artifact kind",
        ));
    }
    if manifest.compiled.runtime_available && manifest.compiled.stage_evaluators.is_none() {
        return Err(PyValueError::new_err(
            "generic DAG schema-v2 artifact marks runtime available without serialized stage evaluators",
        ));
    }
    if manifest.dag_summary.truncated {
        return Err(PyValueError::new_err(
            "generic DAG schema-v2 artifact was truncated during current construction",
        ));
    }
    if manifest.external_pdg_order.len() != manifest.runtime_schema.external_particles.len() {
        return Err(PyValueError::new_err(
            "generic runtime schema external particle count does not match external_pdg_order",
        ));
    }
    for (index, (pdg, particle)) in manifest
        .external_pdg_order
        .iter()
        .zip(&manifest.runtime_schema.external_particles)
        .enumerate()
    {
        if particle.index != index || particle.momentum_slot != index || particle.pdg != *pdg {
            return Err(PyValueError::new_err(format!(
                "generic external particle metadata mismatch at index {index}"
            )));
        }
        if particle.label != index + 1 {
            return Err(PyValueError::new_err(format!(
                "generic external particle labels must be contiguous; got {} at index {index}",
                particle.label
            )));
        }
        if particle.outgoing_pdg == 0 || (particle.role != "initial" && particle.role != "final") {
            return Err(PyValueError::new_err(format!(
                "generic external particle {index} has invalid role or outgoing PDG"
            )));
        }
    }
    validate_generic_parameter_layout(&manifest.runtime_schema)?;
    validate_generic_current_storage(manifest)?;
    validate_generic_value_storage(manifest)?;
    validate_generic_sources(manifest)?;
    validate_generic_momentum_slots(manifest)?;
    validate_generic_stages(manifest)?;
    validate_generic_amplitudes(manifest)?;
    validate_generic_stage_evaluators(manifest)?;
    validate_lc_topology_replay(manifest)?;
    Ok(())
}

fn validate_lc_topology_replay(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let Some(replay) = manifest.compiled.lc_topology_replay.as_ref() else {
        return Ok(());
    };
    if !replay.enabled {
        return Ok(());
    }
    let mappings = build_lc_topology_replay_mappings(Some(replay))?;
    if mappings.is_empty() {
        return Err(PyValueError::new_err(
            "enabled LC topology replay contains no sector mappings",
        ));
    }
    let expected_legs = manifest.external_pdg_order.len();
    for mapping in &mappings {
        let mut seen = vec![false; expected_legs];
        for (representative_index, sector_index) in mapping {
            if *representative_index >= expected_legs || *sector_index >= expected_legs {
                return Err(PyValueError::new_err(
                    "LC topology replay mapping references an out-of-range external leg",
                ));
            }
            if seen[*representative_index] {
                return Err(PyValueError::new_err(
                    "LC topology replay mapping contains a duplicate representative label",
                ));
            }
            seen[*representative_index] = true;
        }
    }
    Ok(())
}

fn validate_generic_parameter_layout(schema: &GenericRuntimeSchemaManifestV2) -> PyResult<()> {
    let layout = &schema.parameter_layout;
    if !layout.source_components_complex || !layout.momentum_components_real {
        return Err(PyValueError::new_err(
            "generic runtime schema requires complex source components and real momenta",
        ));
    }
    if layout.parameter_count_if_flattened
        != layout.source_component_parameter_count + layout.momentum_parameter_count
    {
        return Err(PyValueError::new_err(
            "generic runtime schema flattened parameter count is inconsistent",
        ));
    }
    if layout.momentum_parameter_count != 4 * schema.momentum_slots.len() {
        return Err(PyValueError::new_err(
            "generic runtime schema momentum parameter count does not match momentum slots",
        ));
    }
    let expected_real_inputs = (layout.source_component_parameter_count
        ..layout.parameter_count_if_flattened)
        .collect::<Vec<_>>();
    if layout.real_valued_inputs != expected_real_inputs {
        return Err(PyValueError::new_err(
            "generic runtime schema real-valued input indices are inconsistent",
        ));
    }
    Ok(())
}

fn positive_json_integer(value: &Value) -> bool {
    if let Some(unsigned) = value.as_u64() {
        return unsigned > 0;
    }
    if let Some(signed) = value.as_i64() {
        return signed > 0;
    }
    if let Some(float_value) = value.as_f64() {
        return float_value.is_finite() && float_value > 0.0;
    }
    if let Some(text) = value.as_str() {
        if let Some(hex) = text.strip_prefix("0x").or_else(|| text.strip_prefix("0X")) {
            return !hex.is_empty()
                && hex.bytes().all(|byte| byte.is_ascii_hexdigit())
                && hex.bytes().any(|byte| byte != b'0');
        }
        return !text.is_empty()
            && text.bytes().all(|byte| byte.is_ascii_digit())
            && text.bytes().any(|byte| byte != b'0');
    }
    false
}

fn validate_generic_current_storage(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let storage = &manifest.runtime_schema.current_storage;
    if storage.number_type != "complex" {
        return Err(PyValueError::new_err(
            "generic current storage must use complex number type",
        ));
    }
    if storage.current_slots.len() != manifest.dag_summary.current_count {
        return Err(PyValueError::new_err(
            "generic current slot count does not match DAG summary",
        ));
    }
    let mut offset = 0usize;
    for (index, slot) in storage.current_slots.iter().enumerate() {
        if slot.current_id != index {
            return Err(PyValueError::new_err(format!(
                "generic current slot id mismatch at index {index}"
            )));
        }
        if slot.component_start != offset
            || slot.component_stop != slot.component_start + slot.dimension
            || slot.dimension == 0
        {
            return Err(PyValueError::new_err(format!(
                "generic current slot {index} has inconsistent component range"
            )));
        }
        if slot.particle_id == 0 || slot.external_mask == 0 || slot.momentum_mask == 0 {
            return Err(PyValueError::new_err(format!(
                "generic current slot {index} has invalid physics identity"
            )));
        }
        if slot.external_labels.is_empty()
            || !positive_json_integer(&slot.helicity_ancestry)
            || slot.color_state.is_null()
        {
            return Err(PyValueError::new_err(format!(
                "generic current slot {index} is missing current-index metadata"
            )));
        }
        if slot.chirality.abs() > 1 || slot.flavour_flow.is_empty() {
            return Err(PyValueError::new_err(format!(
                "generic current slot {index} has invalid quantum-flow metadata"
            )));
        }
        if slot.spin_state.is_null()
            || slot.charge_flow.abs() > 1000
            || slot
                .auxiliary_kind
                .as_ref()
                .map_or(false, |kind| kind.is_empty())
        {
            return Err(PyValueError::new_err(format!(
                "generic current slot {index} has invalid spin/charge/auxiliary metadata"
            )));
        }
        offset = slot.component_stop;
    }
    if storage.component_count != offset {
        return Err(PyValueError::new_err(
            "generic current storage component_count is inconsistent",
        ));
    }
    Ok(())
}

fn validate_generic_value_storage(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let schema = &manifest.runtime_schema;
    let storage = &schema.value_storage;
    if storage.number_type != "complex" {
        return Err(PyValueError::new_err(
            "generic value storage must use complex number type",
        ));
    }
    if schema.parameter_layout.value_component_count != storage.component_count {
        return Err(PyValueError::new_err(
            "generic value component count does not match parameter layout",
        ));
    }
    if storage.value_slots.is_empty() && manifest.dag_summary.current_count != 0 {
        return Err(PyValueError::new_err(
            "generic value storage has no value slots",
        ));
    }
    let current_slots = &schema.current_storage.current_slots;
    let mut offset = 0usize;
    for (index, slot) in storage.value_slots.iter().enumerate() {
        if slot.value_slot_id != index
            || slot.component_start != offset
            || slot.component_stop != slot.component_start + slot.dimension
            || slot.dimension == 0
        {
            return Err(PyValueError::new_err(format!(
                "generic value slot {index} has inconsistent component range"
            )));
        }
        let current = current_slots.get(slot.current_id).ok_or_else(|| {
            PyValueError::new_err(format!(
                "generic value slot {index} references missing current {}",
                slot.current_id
            ))
        })?;
        if current.dimension != slot.dimension || current.is_source != slot.is_source {
            return Err(PyValueError::new_err(format!(
                "generic value slot {index} does not match its current slot"
            )));
        }
        if slot.current_component_start != current.component_start
            || slot.current_component_stop != current.component_stop
            || slot.particle_id != current.particle_id
            || slot.external_mask != current.external_mask
            || slot.external_labels != current.external_labels
            || slot.momentum_mask != current.momentum_mask
            || slot.chirality != current.chirality
        {
            return Err(PyValueError::new_err(format!(
                "generic value slot {index} does not preserve its current identity"
            )));
        }
        match slot.variant.as_str() {
            "source" => {
                if !slot.is_source || slot.applies_propagator {
                    return Err(PyValueError::new_err(format!(
                        "generic source value slot {index} is inconsistent"
                    )));
                }
            }
            "propagated" => {
                if slot.is_source || !slot.applies_propagator {
                    return Err(PyValueError::new_err(format!(
                        "generic propagated value slot {index} is inconsistent"
                    )));
                }
            }
            "unpropagated" => {
                if slot.is_source || slot.applies_propagator {
                    return Err(PyValueError::new_err(format!(
                        "generic unpropagated value slot {index} is inconsistent"
                    )));
                }
            }
            other => {
                return Err(PyValueError::new_err(format!(
                    "generic value slot {index} has unsupported variant {other:?}"
                )));
            }
        }
        offset = slot.component_stop;
    }
    if storage.component_count != offset {
        return Err(PyValueError::new_err(
            "generic value storage component_count is inconsistent",
        ));
    }
    Ok(())
}

fn validate_generic_sources(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let schema = &manifest.runtime_schema;
    if schema.source_fill.source_count != manifest.dag_summary.source_count
        || schema.source_fill.sources.len() != manifest.dag_summary.source_count
    {
        return Err(PyValueError::new_err(
            "generic source count does not match DAG summary",
        ));
    }
    let current_slots = &schema.current_storage.current_slots;
    let mut source_offset = 0usize;
    for (index, source) in schema.source_fill.sources.iter().enumerate() {
        if source.source_id != index {
            return Err(PyValueError::new_err(format!(
                "generic source id mismatch at index {index}"
            )));
        }
        let slot = current_slots.get(source.current_id).ok_or_else(|| {
            PyValueError::new_err(format!(
                "generic source {index} references missing current {}",
                source.current_id
            ))
        })?;
        if !slot.is_source
            || slot.dimension != source.dimension
            || slot.particle_id != source.particle_id
        {
            return Err(PyValueError::new_err(format!(
                "generic source {index} does not match its current slot"
            )));
        }
        if source.current_component_start != slot.component_start
            || source.current_component_stop != slot.component_stop
            || source.source_parameter_start != source_offset
            || source.source_parameter_stop != source_offset + source.dimension
        {
            return Err(PyValueError::new_err(format!(
                "generic source {index} has inconsistent component offsets"
            )));
        }
        validate_value_slot_ref(
            &source.value_slot,
            source.current_id,
            Some("source"),
            &schema.value_storage,
        )?;
        if source.source_kind != "external-wavefunction" {
            return Err(PyValueError::new_err(format!(
                "unsupported generic source kind {:?}",
                source.source_kind
            )));
        }
        if source.side != "initial" && source.side != "final" {
            return Err(PyValueError::new_err(format!(
                "generic source {index} has invalid side {:?}",
                source.side
            )));
        }
        if source.physical_pdg == 0
            || source.outgoing_pdg != source.particle_id
            || source.chirality.abs() > 1
            || source.source_helicity.abs() > 1
            || source.spin_state.is_null()
            || !positive_json_integer(&source.helicity_ancestry)
            || source.color_state.is_null()
        {
            return Err(PyValueError::new_err(format!(
                "generic source {index} has invalid physics metadata"
            )));
        }
        if source.crossing != "identity" && source.crossing != "negate-incoming-momentum" {
            return Err(PyValueError::new_err(format!(
                "generic source {index} has unsupported crossing {:?}",
                source.crossing
            )));
        }
        let particle = schema
            .external_particles
            .get(source.input_momentum_slot)
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "generic source {index} references missing momentum slot {}",
                    source.input_momentum_slot
                ))
            })?;
        if particle.label != source.leg_label {
            return Err(PyValueError::new_err(format!(
                "generic source {index} leg label does not match momentum slot"
            )));
        }
        source_offset = source.source_parameter_stop;
    }
    if source_offset != schema.parameter_layout.source_component_parameter_count {
        return Err(PyValueError::new_err(
            "generic source parameter count does not match source records",
        ));
    }
    Ok(())
}

fn validate_generic_momentum_slots(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let schema = &manifest.runtime_schema;
    for (index, slot) in schema.momentum_slots.iter().enumerate() {
        if slot.momentum_slot_id != index
            || slot.component_start != 4 * index
            || slot.component_stop != slot.component_start + 4
            || !slot.real_valued
            || slot.momentum_mask == 0
            || slot.external_labels.is_empty()
        {
            return Err(PyValueError::new_err(format!(
                "generic momentum slot {index} is inconsistent"
            )));
        }
        for label in &slot.external_labels {
            if *label == 0 || *label > schema.external_particles.len() {
                return Err(PyValueError::new_err(format!(
                    "generic momentum slot {index} has invalid external label {label}"
                )));
            }
        }
    }
    Ok(())
}

fn validate_generic_stages(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let schema = &manifest.runtime_schema;
    let current_count = schema.current_storage.current_slots.len();
    let value_count = schema.value_storage.value_slots.len();
    let momentum_count = schema.momentum_slots.len();
    let mut seen_interactions = 0usize;
    for (stage_offset, stage) in schema.stages.iter().enumerate() {
        if stage.stage_index != stage_offset + 1
            || stage.stage_kind != "current-combine"
            || stage.interaction_count != stage.interactions.len()
            || stage.subset_size < 2
        {
            return Err(PyValueError::new_err(format!(
                "generic stage {} is inconsistent",
                stage.stage_index
            )));
        }
        for id in stage
            .input_current_ids
            .iter()
            .chain(&stage.output_current_ids)
        {
            if *id >= current_count {
                return Err(PyValueError::new_err(format!(
                    "generic stage {} references invalid current {id}",
                    stage.stage_index
                )));
            }
        }
        for value_id in stage
            .input_value_slot_ids
            .iter()
            .chain(&stage.output_value_slot_ids)
        {
            if *value_id >= value_count {
                return Err(PyValueError::new_err(format!(
                    "generic stage {} references invalid value slot {value_id}",
                    stage.stage_index
                )));
            }
        }
        for interaction in &stage.interactions {
            validate_generic_interaction(
                interaction,
                &schema.current_storage,
                &schema.value_storage,
                momentum_count,
            )?;
            if !stage
                .input_current_ids
                .contains(&interaction.left_current_id)
                || !stage
                    .input_current_ids
                    .contains(&interaction.right_current_id)
                || !stage
                    .output_current_ids
                    .contains(&interaction.result_current_id)
            {
                return Err(PyValueError::new_err(format!(
                    "generic interaction {} is not listed in its stage inputs/outputs",
                    interaction.interaction_id
                )));
            }
            if !stage
                .input_value_slot_ids
                .contains(&interaction.left_value_slot.value_slot_id)
                || !stage
                    .input_value_slot_ids
                    .contains(&interaction.right_value_slot.value_slot_id)
                || interaction
                    .result_value_slots
                    .iter()
                    .any(|slot| !stage.output_value_slot_ids.contains(&slot.value_slot_id))
            {
                return Err(PyValueError::new_err(format!(
                    "generic interaction {} value slots are not listed in its stage inputs/outputs",
                    interaction.interaction_id
                )));
            }
        }
        seen_interactions += stage.interactions.len();
    }
    if seen_interactions != manifest.dag_summary.interaction_count {
        return Err(PyValueError::new_err(
            "generic stage interaction count does not match DAG summary",
        ));
    }
    Ok(())
}

fn validate_generic_interaction(
    interaction: &GenericInteractionManifestV2,
    current_storage: &GenericCurrentStorageManifestV2,
    value_storage: &GenericValueStorageManifestV2,
    momentum_count: usize,
) -> PyResult<()> {
    let current_count = current_storage.current_slots.len();
    for current_id in [
        interaction.left_current_id,
        interaction.right_current_id,
        interaction.result_current_id,
    ] {
        if current_id >= current_count {
            return Err(PyValueError::new_err(format!(
                "generic interaction {} references invalid current {current_id}",
                interaction.interaction_id
            )));
        }
    }
    for slot_id in [
        interaction.momentum_slots.left,
        interaction.momentum_slots.right,
        interaction.momentum_slots.result,
    ] {
        if slot_id >= momentum_count {
            return Err(PyValueError::new_err(format!(
                "generic interaction {} references invalid momentum slot {slot_id}",
                interaction.interaction_id
            )));
        }
    }
    validate_slot_ref(&interaction.left_slot, interaction.left_current_id)?;
    validate_slot_ref(&interaction.right_slot, interaction.right_current_id)?;
    validate_slot_ref(&interaction.result_slot, interaction.result_current_id)?;
    validate_value_slot_ref(
        &interaction.left_value_slot,
        interaction.left_current_id,
        Some(input_value_variant(
            current_storage,
            value_storage,
            interaction.left_current_id,
        )?),
        value_storage,
    )?;
    validate_value_slot_ref(
        &interaction.right_value_slot,
        interaction.right_current_id,
        Some(input_value_variant(
            current_storage,
            value_storage,
            interaction.right_current_id,
        )?),
        value_storage,
    )?;
    if interaction.result_value_slots.is_empty() {
        return Err(PyValueError::new_err(format!(
            "generic interaction {} has no result value slots",
            interaction.interaction_id
        )));
    }
    for result_slot in &interaction.result_value_slots {
        validate_value_slot_ref(
            result_slot,
            interaction.result_current_id,
            None,
            value_storage,
        )?;
        if result_slot.variant != "propagated" && result_slot.variant != "unpropagated" {
            return Err(PyValueError::new_err(format!(
                "generic interaction {} has invalid result value variant {:?}",
                interaction.interaction_id, result_slot.variant
            )));
        }
    }
    if interaction.vertex_kind < 0 {
        return Err(PyValueError::new_err(format!(
            "generic interaction {} has invalid vertex kind {}",
            interaction.interaction_id, interaction.vertex_kind
        )));
    }
    if interaction.vertex_particles.len() != 3
        || interaction.coupling.len() != 2
        || interaction.color_weight.len() != 2
        || interaction.accumulation != "sum-into-result-current"
        || interaction.lowering.kind != interaction.vertex_kind
        || interaction.lowering.full_tensor_network_ready != interaction.full_tensor_network_ready
        || interaction.lowering.backend.is_empty()
        || interaction.lowering.expression_head.is_empty()
        || interaction.lowering.kernel.is_empty()
        || interaction.lowering.input_roles.len() != 2
        || interaction.lowering.output_role.is_empty()
        || interaction.lowering.coupling_mode.is_empty()
        || interaction
            .lowering
            .tensor_names
            .iter()
            .any(|name| name.is_empty())
        || interaction.lowering.description.is_empty()
    {
        return Err(PyValueError::new_err(format!(
            "generic interaction {} has inconsistent model/lowering metadata",
            interaction.interaction_id
        )));
    }
    if interaction.result_requires_propagated_value
        && !interaction
            .result_value_slots
            .iter()
            .any(|slot| slot.variant == "propagated")
    {
        return Err(PyValueError::new_err(format!(
            "generic interaction {} declares a propagated result without a propagated slot",
            interaction.interaction_id
        )));
    }
    if interaction.result_requires_unpropagated_value
        && !interaction
            .result_value_slots
            .iter()
            .any(|slot| slot.variant == "unpropagated")
    {
        return Err(PyValueError::new_err(format!(
            "generic interaction {} declares an unpropagated result without an unpropagated slot",
            interaction.interaction_id
        )));
    }
    if !interaction.full_tensor_network_ready {
        return Err(PyValueError::new_err(format!(
            "generic interaction {} is not ready for tensor-network lowering",
            interaction.interaction_id
        )));
    }
    Ok(())
}

fn validate_generic_amplitudes(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let schema = &manifest.runtime_schema;
    let stage = &schema.amplitude_stage;
    if stage.stage_kind != "amplitude-roots"
        || stage.output_count != stage.roots.len()
        || stage.output_count != manifest.dag_summary.amplitude_root_count
    {
        return Err(PyValueError::new_err(
            "generic amplitude stage output count is inconsistent",
        ));
    }
    let current_count = schema.current_storage.current_slots.len();
    for (index, root) in stage.roots.iter().enumerate() {
        if root.output_index != index || root.root_id != index {
            return Err(PyValueError::new_err(format!(
                "generic amplitude root index mismatch at output {index}"
            )));
        }
        if root.left_current_id >= current_count || root.right_current_id >= current_count {
            return Err(PyValueError::new_err(format!(
                "generic amplitude root {index} references invalid current"
            )));
        }
        if root.kind != "direct-contraction" && root.kind != "vertex-closure" {
            return Err(PyValueError::new_err(format!(
                "generic amplitude root {index} has unsupported kind {:?}",
                root.kind
            )));
        }
        if root.coupling.len() != 2
            || root.color_weight.len() != 2
            || root.contraction.is_empty()
            || !root.helicity_weight.is_finite()
            || root
                .coherent_group_id
                .as_ref()
                .map_or(false, Value::is_null)
        {
            return Err(PyValueError::new_err(format!(
                "generic amplitude root {index} has inconsistent physics metadata"
            )));
        }
        if root.kind == "vertex-closure"
            && (root.vertex_kind.is_none()
                || root
                    .vertex_particles
                    .as_ref()
                    .map_or(true, |particles| particles.len() != 3))
        {
            return Err(PyValueError::new_err(format!(
                "generic vertex-closure root {index} is missing vertex metadata"
            )));
        }
        validate_slot_ref(&root.left_slot, root.left_current_id)?;
        validate_slot_ref(&root.right_slot, root.right_current_id)?;
        validate_value_slot_ref(
            &root.left_value_slot,
            root.left_current_id,
            Some(amplitude_value_variant(
                &schema.current_storage,
                root.left_current_id,
            )?),
            &schema.value_storage,
        )?;
        validate_value_slot_ref(
            &root.right_value_slot,
            root.right_current_id,
            Some(amplitude_value_variant(
                &schema.current_storage,
                root.right_current_id,
            )?),
            &schema.value_storage,
        )?;
    }
    Ok(())
}

fn validate_generic_stage_evaluators(manifest: &GenericProcessManifestV2) -> PyResult<()> {
    let Some(stage_evaluators) = &manifest.compiled.stage_evaluators else {
        return Ok(());
    };
    let schema = &manifest.runtime_schema;
    let expected_parameter_count = schema.parameter_layout.value_component_count
        + schema.parameter_layout.momentum_parameter_count;
    let expected_real_inputs = (schema.parameter_layout.value_component_count
        ..expected_parameter_count)
        .collect::<Vec<_>>();
    let header_is_global = stage_evaluators.parameter_layout == "global-value-momentum"
        && stage_evaluators.parameter_count == expected_parameter_count
        && stage_evaluators.value_parameter_count == schema.parameter_layout.value_component_count
        && stage_evaluators.momentum_parameter_count
            == schema.parameter_layout.momentum_parameter_count
        && stage_evaluators.real_valued_inputs == expected_real_inputs;
    let header_is_stage_local = stage_evaluators.parameter_layout == "stage-local-value-momentum"
        && stage_evaluators.parameter_count == 0
        && stage_evaluators.value_parameter_count == 0
        && stage_evaluators.momentum_parameter_count == 0
        && stage_evaluators.real_valued_inputs.is_empty();
    if stage_evaluators.kind != "generic-dag-stage-evaluator-artifacts"
        || (!header_is_global && !header_is_stage_local)
        || stage_evaluators.stage_count != schema.stages.len() + 1
        || stage_evaluators.stages.len() != schema.stages.len()
        || (!stage_evaluators.runtime_available
            && stage_evaluators
                .runtime_unavailable_message
                .as_deref()
                .unwrap_or("")
                .is_empty())
    {
        return Err(PyValueError::new_err(
            "generic stage evaluator artifact header is inconsistent with runtime schema",
        ));
    }
    for (stage, runtime_stage) in stage_evaluators.stages.iter().zip(&schema.stages) {
        validate_generic_serialized_stage_evaluator(
            stage,
            runtime_stage.stage_index,
            "current-combine",
            Some(runtime_stage.subset_size),
            expected_parameter_count,
            Some(runtime_stage),
            None,
        )?;
    }
    validate_generic_serialized_stage_evaluator(
        &stage_evaluators.amplitude_stage,
        0,
        "amplitude-roots",
        None,
        expected_parameter_count,
        None,
        Some(&schema.amplitude_stage),
    )?;
    Ok(())
}

fn validate_generic_serialized_stage_evaluator(
    stage: &GenericSerializedStageEvaluatorManifestV2,
    expected_stage_index: usize,
    expected_stage_kind: &str,
    expected_subset_size: Option<usize>,
    expected_parameter_count: usize,
    runtime_stage: Option<&GenericStageManifestV2>,
    amplitude_stage: Option<&GenericAmplitudeStageManifestV2>,
) -> PyResult<()> {
    let stage_is_global = stage.parameter_layout == "global-value-momentum";
    let stage_is_local = stage.parameter_layout == "stage-local-value-momentum";
    if stage.stage_index != expected_stage_index
        || stage.stage_kind != expected_stage_kind
        || stage.subset_size != expected_subset_size
        || (!stage_is_global && !stage_is_local)
        || !stage.expression_ready
        || !stage.blockers.is_empty()
        || stage.evaluator_label.is_empty()
        || stage.output_length == 0
    {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} is inconsistent",
            stage.evaluator_label
        )));
    }
    validate_generic_stage_output_slots(stage)?;
    let (input_len, output_len) = evaluator_manifest_io_len(&stage.evaluator)?;
    if stage_is_global && input_len != expected_parameter_count {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} has inconsistent global evaluator input length",
            stage.evaluator_label
        )));
    }
    if stage_is_local {
        validate_generic_stage_input_components(stage, expected_parameter_count)?;
        if input_len != stage.parameter_count
            || stage.parameter_count != stage.input_components.len()
            || stage.value_parameter_count + stage.momentum_parameter_count != stage.parameter_count
            || stage
                .real_valued_inputs
                .iter()
                .any(|index| *index >= stage.parameter_count)
        {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} has inconsistent local evaluator input metadata",
                stage.evaluator_label
            )));
        }
    }
    if output_len != stage.output_length {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} has inconsistent evaluator IO",
            stage.evaluator_label
        )));
    }
    if let Some(runtime_stage) = runtime_stage {
        let expected_interactions = runtime_stage
            .interactions
            .iter()
            .map(|interaction| interaction.interaction_id)
            .collect::<Vec<_>>();
        if stage.input_value_slot_ids != runtime_stage.input_value_slot_ids
            || stage.output_value_slot_ids != runtime_stage.output_value_slot_ids
            || stage.interaction_ids != expected_interactions
        {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} does not match runtime stage slots",
                stage.evaluator_label
            )));
        }
    }
    if let Some(amplitude_stage) = amplitude_stage {
        if !stage.output_value_slot_ids.is_empty()
            || !stage.interaction_ids.is_empty()
            || stage.output_length != amplitude_stage.output_count
        {
            return Err(PyValueError::new_err(
                "generic serialized amplitude evaluator does not match amplitude stage",
            ));
        }
    }
    Ok(())
}

fn validate_generic_stage_output_slots(
    stage: &GenericSerializedStageEvaluatorManifestV2,
) -> PyResult<()> {
    if stage.output_slots.is_empty() {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} has no output slots",
            stage.evaluator_label
        )));
    }
    let mut max_output_stop = 0usize;
    for slot in &stage.output_slots {
        if slot.variant.is_empty()
            || slot.component_stop < slot.component_start
            || slot.output_stop <= slot.output_start
            || slot.output_stop > stage.output_length
        {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} has invalid output slot",
                stage.evaluator_label
            )));
        }
        if stage.stage_kind == "amplitude-roots" {
            if slot.value_slot_id != -1 || slot.current_id != -1 || slot.variant != "amplitude-root"
            {
                return Err(PyValueError::new_err(
                    "generic serialized amplitude evaluator has invalid output slot metadata",
                ));
            }
        } else if slot.value_slot_id < 0 || slot.current_id < 0 || slot.variant == "amplitude-root"
        {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} has invalid current output slot metadata",
                stage.evaluator_label
            )));
        }
        max_output_stop = max_output_stop.max(slot.output_stop);
    }
    if max_output_stop != stage.output_length {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} output slots do not cover evaluator outputs",
            stage.evaluator_label
        )));
    }
    Ok(())
}

fn validate_generic_stage_input_components(
    stage: &GenericSerializedStageEvaluatorManifestV2,
    global_parameter_count: usize,
) -> PyResult<()> {
    if stage.input_components.is_empty() {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} has no local input components",
            stage.evaluator_label
        )));
    }
    let mut seen_parameters = BTreeSet::new();
    for component in &stage.input_components {
        if component.kind != "value" && component.kind != "momentum" {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} has invalid local input kind {:?}",
                stage.evaluator_label, component.kind
            )));
        }
        if component.global_component >= global_parameter_count
            || component.parameter_index >= stage.parameter_count
            || component.source_id >= global_parameter_count
            || component.component >= global_parameter_count
            || !seen_parameters.insert(component.parameter_index)
        {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} has invalid local input component map",
                stage.evaluator_label
            )));
        }
        if component.real_valued && component.kind != "momentum" {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} marks a non-momentum input as real",
                stage.evaluator_label
            )));
        }
    }
    if seen_parameters.len() != stage.parameter_count {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} local input map is not dense",
            stage.evaluator_label
        )));
    }
    let real_inputs = stage
        .input_components
        .iter()
        .filter_map(|component| {
            if component.real_valued {
                Some(component.parameter_index)
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if real_inputs != stage.real_valued_inputs {
        return Err(PyValueError::new_err(format!(
            "generic serialized stage evaluator {} local real-input metadata is inconsistent",
            stage.evaluator_label
        )));
    }
    Ok(())
}

fn evaluator_manifest_io_len(manifest: &EvaluatorManifest) -> PyResult<(usize, usize)> {
    match manifest {
        EvaluatorManifest::Jit {
            input_len,
            output_len,
            ..
        }
        | EvaluatorManifest::CompiledComplex {
            input_len,
            output_len,
            ..
        } => Ok((*input_len, *output_len)),
        EvaluatorManifest::Chunked { chunks } => {
            let mut iter = chunks.iter();
            let first = iter.next().ok_or_else(|| {
                PyValueError::new_err("generic serialized evaluator chunk list is empty")
            })?;
            let (input_len, mut output_len) = evaluator_manifest_io_len(first)?;
            for chunk in iter {
                let (chunk_input_len, chunk_output_len) = evaluator_manifest_io_len(chunk)?;
                if chunk_input_len != input_len {
                    return Err(PyValueError::new_err(
                        "generic serialized evaluator chunks have inconsistent input lengths",
                    ));
                }
                output_len += chunk_output_len;
            }
            Ok((input_len, output_len))
        }
    }
}

fn validate_slot_ref(slot: &GenericSlotRefManifestV2, current_id: usize) -> PyResult<()> {
    if slot.current_id != current_id
        || slot.component_stop != slot.component_start + slot.dimension
        || slot.dimension == 0
    {
        return Err(PyValueError::new_err(format!(
            "generic slot reference for current {current_id} is inconsistent"
        )));
    }
    Ok(())
}

fn validate_value_slot_ref(
    slot: &GenericValueSlotRefManifestV2,
    current_id: usize,
    expected_variant: Option<&str>,
    value_storage: &GenericValueStorageManifestV2,
) -> PyResult<()> {
    if slot.current_id != current_id
        || slot.component_stop != slot.component_start + slot.dimension
        || slot.dimension == 0
    {
        return Err(PyValueError::new_err(format!(
            "generic value slot reference for current {current_id} is inconsistent"
        )));
    }
    if let Some(expected) = expected_variant {
        if slot.variant != expected {
            return Err(PyValueError::new_err(format!(
                "generic value slot reference for current {current_id} uses variant {:?}, expected {expected:?}",
                slot.variant
            )));
        }
    }
    let storage_slot = value_storage
        .value_slots
        .get(slot.value_slot_id)
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "generic value slot reference for current {current_id} points outside value storage"
            ))
        })?;
    if storage_slot.current_id != slot.current_id
        || storage_slot.variant != slot.variant
        || storage_slot.component_start != slot.component_start
        || storage_slot.component_stop != slot.component_stop
        || storage_slot.dimension != slot.dimension
    {
        return Err(PyValueError::new_err(format!(
            "generic value slot reference for current {current_id} does not match value storage"
        )));
    }
    Ok(())
}

fn input_value_variant(
    current_storage: &GenericCurrentStorageManifestV2,
    value_storage: &GenericValueStorageManifestV2,
    current_id: usize,
) -> PyResult<&'static str> {
    let current = current_storage
        .current_slots
        .get(current_id)
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "generic input value references missing current {current_id}"
            ))
        })?;
    if current.is_source {
        return Ok("source");
    }
    let mut has_unpropagated = false;
    for slot in &value_storage.value_slots {
        if slot.current_id != current_id {
            continue;
        }
        if slot.variant == "propagated" {
            return Ok("propagated");
        }
        if slot.variant == "unpropagated" {
            has_unpropagated = true;
        }
    }
    if has_unpropagated {
        return Ok("unpropagated");
    }
    Err(PyValueError::new_err(format!(
        "generic input value references current {current_id} without an input value slot"
    )))
}

fn amplitude_value_variant(
    current_storage: &GenericCurrentStorageManifestV2,
    current_id: usize,
) -> PyResult<&'static str> {
    let current = current_storage
        .current_slots
        .get(current_id)
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "generic amplitude value references missing current {current_id}"
            ))
        })?;
    Ok(if current.is_source {
        "source"
    } else {
        "unpropagated"
    })
}

fn validate_amplitude_group_metadata(
    amplitudes: &[AmplitudeManifest],
    stage_group_ids: Option<&[Option<i64>]>,
) -> PyResult<()> {
    let Some(stage_group_ids) = stage_group_ids else {
        return Ok(());
    };
    if amplitudes.len() != stage_group_ids.len() {
        return Err(PyValueError::new_err(
            "amplitude metadata and raw-sum group ids have different lengths",
        ));
    }
    for (index, (amplitude, stage_group_id)) in amplitudes.iter().zip(stage_group_ids).enumerate() {
        if amplitude.coherent_group_id != *stage_group_id {
            return Err(PyValueError::new_err(format!(
                "amplitude coherent group id mismatch at index {index}"
            )));
        }
    }
    Ok(())
}

#[derive(Clone, Debug, Deserialize)]
struct OutputSlotManifest {
    id: usize,
    start: usize,
    stop: usize,
}

#[derive(Clone, Debug, Deserialize)]
#[serde(tag = "kind")]
enum EvaluatorManifest {
    #[serde(rename = "jit-symbolica-evaluator")]
    Jit {
        input_len: usize,
        output_len: usize,
        evaluator_state_path: String,
    },
    #[serde(rename = "compiled-complex-evaluator")]
    CompiledComplex {
        function_name: String,
        input_len: usize,
        output_len: usize,
        library_path: String,
        evaluator_state_path: Option<String>,
        number_type: String,
    },
    #[serde(rename = "chunked-symbolica-evaluator")]
    Chunked { chunks: Vec<EvaluatorManifest> },
}

struct EvaluatorGroup {
    evaluators: Vec<LoadedEvaluator>,
    output_len: usize,
    chunk_scratch_f64: Vec<Complex<f64>>,
}

enum F64Evaluator {
    Compiled(CompiledComplexEvaluator),
    Jit(JITCompiledEvaluator<Complex<f64>>),
    Interpreted(ExpressionEvaluator<Complex<f64>>),
}

struct LoadedEvaluator {
    eval: F64Evaluator,
    exact_eval: Option<ExpressionEvaluator<Complex<Rational>>>,
    double_eval: Option<ExpressionEvaluator<Complex<DoubleFloat>>>,
    arb_eval: Option<(u32, ExpressionEvaluator<Complex<Float>>)>,
    input_len: usize,
    output_len: usize,
}

trait RusticolHighPrecisionNumber:
    Real + RealLike + From<f64> + PartialOrd + Clone + EvaluationDomain
where
    Complex<Self>: Real + EvaluationDomain,
{
    fn evaluate_loaded(
        evaluator: &mut LoadedEvaluator,
        params: &[Complex<Self>],
        out: &mut [Complex<Self>],
        binary_precision: Option<u32>,
    ) -> PyResult<()>;
}

impl RusticolHighPrecisionNumber for DoubleFloat {
    fn evaluate_loaded(
        evaluator: &mut LoadedEvaluator,
        params: &[Complex<Self>],
        out: &mut [Complex<Self>],
        _binary_precision: Option<u32>,
    ) -> PyResult<()> {
        if evaluator.double_eval.is_none() {
            let exact = evaluator.exact_eval.as_ref().ok_or_else(|| {
                PyValueError::new_err(
                    "precision 32 requires evaluator-state artifacts, but this process \
                     artifact has no evaluator_state_path for one or more chunks",
                )
            })?;
            evaluator.double_eval =
                Some(exact.clone().map_coeff(&|c| {
                    Complex::new(DoubleFloat::from(&c.re), DoubleFloat::from(&c.im))
                }));
        }
        evaluator
            .double_eval
            .as_mut()
            .expect("double evaluator initialized")
            .evaluate(params, out);
        Ok(())
    }
}

impl RusticolHighPrecisionNumber for Float {
    fn evaluate_loaded(
        evaluator: &mut LoadedEvaluator,
        params: &[Complex<Self>],
        out: &mut [Complex<Self>],
        binary_precision: Option<u32>,
    ) -> PyResult<()> {
        let binary_precision = binary_precision.ok_or_else(|| {
            PyValueError::new_err("arbitrary-precision evaluation needs a binary precision")
        })?;
        let rebuild = evaluator
            .arb_eval
            .as_ref()
            .map(|(precision, _)| *precision != binary_precision)
            .unwrap_or(true);
        if rebuild {
            let exact = evaluator.exact_eval.as_ref().ok_or_else(|| {
                PyValueError::new_err(
                    "arbitrary-precision evaluation requires evaluator-state artifacts, but \
                     this process artifact has no evaluator_state_path for one or more chunks",
                )
            })?;
            evaluator.arb_eval = Some((
                binary_precision,
                exact.clone().map_coeff_with_prec(
                    &|c| {
                        Complex::new(
                            c.re.to_multi_prec_float(binary_precision),
                            c.im.to_multi_prec_float(binary_precision),
                        )
                    },
                    binary_precision,
                ),
            ));
        }
        evaluator
            .arb_eval
            .as_mut()
            .expect("arbitrary-precision evaluator initialized")
            .1
            .evaluate(params, out);
        Ok(())
    }
}

struct CurrentStage {
    outputs: Vec<(usize, usize, usize, usize)>,
    chunk_outputs: Vec<Vec<(usize, usize)>>,
    evaluator: EvaluatorGroup,
}

struct AmplitudeStage {
    output_length: usize,
    raw_sum_weights: Vec<f64>,
    raw_sum_groups: Vec<RawSumGroup>,
    has_coherent_groups: bool,
    amplitude_evaluator: Option<EvaluatorGroup>,
    raw_sum_evaluator: Option<EvaluatorGroup>,
}

struct RawSumGroup {
    id: i64,
    indices: Vec<usize>,
    weight: f64,
}

struct ColorContractionRuntime {
    group_count: usize,
    entries: Vec<ColorContractionEntry>,
    group_scratch_f64: Vec<Complex<f64>>,
}

struct ColorContractionEntry {
    left_group_index: usize,
    right_group_index: usize,
    weight_re: f64,
    weight_im: f64,
    symmetry_factor: f64,
}

struct ZeroGluonStage {
    parameter_names: Vec<String>,
    evaluator: LoadedEvaluator,
    z_left: f64,
    z_right: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct RuntimeProfile {
    source_fill_s: f64,
    momentum_setup_s: f64,
    stage_evaluator_s: f64,
    output_assign_s: f64,
    amplitude_evaluator_s: f64,
    reduction_s: f64,
    total_s: f64,
}

#[derive(Clone, Copy, Debug, Default)]
struct MemorySnapshot {
    current_rss_bytes: Option<u64>,
    peak_rss_bytes: Option<u64>,
}

#[derive(Clone, Copy, Debug, Default)]
struct ComplexChecksum {
    sum_re: f64,
    sum_im: f64,
    sum_abs2: f64,
    max_abs: f64,
    output_len: usize,
}

struct GenericRuntimeV2 {
    process: String,
    key: String,
    color_accuracy: String,
    external_pdg_order: Vec<i32>,
    external_count: usize,
    parameter_count: usize,
    value_parameter_count: usize,
    momentum_parameter_count: usize,
    current_count: usize,
    source_count: usize,
    interaction_count: usize,
    stage_count: usize,
    amplitude_output_count: usize,
    stage_evaluator_count: usize,
    stage_evaluator_labels: Vec<String>,
    amplitude_evaluator_label: Option<String>,
    lc_topology_reuse_available: bool,
    lc_topology_group_count: usize,
    lc_topology_representative_sector_ids: Vec<i64>,
    lc_topology_replay_enabled: bool,
    lc_topology_replay_sector_count: usize,
    lc_topology_replay_mappings: Vec<Vec<(usize, usize)>>,
    runtime_unavailable_message: Option<String>,
    sources: Vec<GenericSourceRecordManifestV2>,
    momentum_slots: Vec<GenericMomentumSlotManifestV2>,
    external_is_initial: Vec<bool>,
    particle_masses: BTreeMap<i32, f64>,
    normalization_factor: f64,
    stages: Option<Vec<GenericStageRuntimeV2>>,
    amplitude_stage: Option<GenericAmplitudeRuntimeV2>,
    state_scratch_f64: Vec<Complex<f64>>,
    values_scratch_f64: Vec<f64>,
}

struct GenericStageRuntimeV2 {
    outputs: Vec<(usize, usize)>,
    output_spans: Vec<(usize, usize, usize)>,
    input_components: Option<Vec<usize>>,
    parameter_scratch_f64: Vec<Complex<f64>>,
    output_scratch_f64: Vec<Complex<f64>>,
    evaluator: EvaluatorGroup,
}

struct GenericAmplitudeRuntimeV2 {
    output_length: usize,
    raw_sum_weights: Vec<f64>,
    raw_sum_groups: Vec<RawSumGroup>,
    has_coherent_groups: bool,
    color_contraction: Option<ColorContractionRuntime>,
    input_components: Option<Vec<usize>>,
    parameter_scratch_f64: Vec<Complex<f64>>,
    output_scratch_f64: Vec<Complex<f64>>,
    evaluator: EvaluatorGroup,
}

fn build_lc_topology_replay_mappings(
    replay: Option<&LcTopologyReplayManifestV2>,
) -> PyResult<Vec<Vec<(usize, usize)>>> {
    let Some(replay) = replay else {
        return Ok(Vec::new());
    };
    if !replay.enabled {
        return Ok(Vec::new());
    }
    if replay.mode != "external-label-permutation" {
        return Err(PyValueError::new_err(format!(
            "unsupported LC topology replay mode {:?}",
            replay.mode
        )));
    }
    let mut mappings = Vec::new();
    for group in &replay.groups {
        if group.materialized_sector_id != group.representative_sector_id {
            return Err(PyValueError::new_err(
                "LC topology replay currently requires the materialized sector to be the representative sector",
            ));
        }
        if group.sector_permutations.is_empty() {
            return Err(PyValueError::new_err(
                "enabled LC topology replay group contains no sector permutations",
            ));
        }
        if !group.active_sector_ids.is_empty() {
            if !group
                .active_sector_ids
                .contains(&group.representative_sector_id)
            {
                return Err(PyValueError::new_err(
                    "LC topology replay active sector ids do not include the representative sector",
                ));
            }
            let mut seen_active = BTreeSet::new();
            for sector_id in &group.active_sector_ids {
                if !seen_active.insert(*sector_id) {
                    return Err(PyValueError::new_err(
                        "LC topology replay active sector ids contain duplicates",
                    ));
                }
            }
        }
        let mut seen_permutations = BTreeSet::new();
        for permutation in &group.sector_permutations {
            if !group.active_sector_ids.is_empty()
                && !group.active_sector_ids.contains(&permutation.sector_id)
            {
                return Err(PyValueError::new_err(
                    "LC topology replay sector permutation is not listed in active sector ids",
                ));
            }
            if !seen_permutations.insert(permutation.sector_id) {
                return Err(PyValueError::new_err(
                    "LC topology replay sector permutations contain duplicate sector ids",
                ));
            }
            let mut mapping = Vec::new();
            for item in &permutation.label_permutation {
                if item.representative_label == 0 || item.sector_label == 0 {
                    return Err(PyValueError::new_err(
                        "LC topology replay label permutations must use one-based labels",
                    ));
                }
                mapping.push((item.representative_label - 1, item.sector_label - 1));
            }
            mappings.push(mapping);
        }
    }
    if mappings.len() != replay.replayed_sector_count {
        return Err(PyValueError::new_err(format!(
            "LC topology replay declares {} sectors but contains {} permutations",
            replay.replayed_sector_count,
            mappings.len()
        )));
    }
    Ok(mappings)
}

impl GenericRuntimeV2 {
    fn from_manifest(manifest: GenericProcessManifestV2) -> PyResult<Self> {
        let stage_evaluators = manifest.compiled.stage_evaluators.as_ref();
        let topology_reuse = manifest.lc_topology_reuse.as_ref();
        let topology_replay = manifest.compiled.lc_topology_replay.as_ref();
        let topology_replay_mappings = build_lc_topology_replay_mappings(topology_replay)?;
        let external_is_initial = manifest
            .runtime_schema
            .external_particles
            .iter()
            .map(|particle| particle.role == "initial")
            .collect::<Vec<_>>();
        let particle_masses = manifest
            .runtime_schema
            .model
            .as_ref()
            .map(|model| {
                model
                    .particles
                    .iter()
                    .map(|particle| (particle.pdg, particle.mass))
                    .collect::<BTreeMap<_, _>>()
            })
            .unwrap_or_default();
        let color_factor_in_contraction = manifest
            .runtime_schema
            .amplitude_stage
            .color_contraction
            .as_ref()
            .map(|contraction| contraction.supported && contraction.includes_color_factor)
            .unwrap_or(false);
        let normalization_factor = manifest
            .runtime_schema
            .normalization
            .as_ref()
            .map(|normalization| {
                let color_factor = if color_factor_in_contraction {
                    1.0
                } else {
                    normalization.color_factor
                };
                color_factor * normalization.global_coupling_factor
                    / (normalization.average_factor * normalization.identical_factor)
            })
            .unwrap_or(1.0);
        Ok(Self {
            process: manifest.process,
            key: manifest.key,
            color_accuracy: manifest.color_accuracy,
            external_count: manifest.external_pdg_order.len(),
            external_pdg_order: manifest.external_pdg_order,
            parameter_count: manifest
                .runtime_schema
                .parameter_layout
                .value_component_count
                + manifest
                    .runtime_schema
                    .parameter_layout
                    .momentum_parameter_count,
            value_parameter_count: manifest
                .runtime_schema
                .parameter_layout
                .value_component_count,
            momentum_parameter_count: manifest
                .runtime_schema
                .parameter_layout
                .momentum_parameter_count,
            current_count: manifest.dag_summary.current_count,
            source_count: manifest.dag_summary.source_count,
            interaction_count: manifest.dag_summary.interaction_count,
            stage_count: manifest.runtime_schema.stages.len(),
            amplitude_output_count: manifest.runtime_schema.amplitude_stage.output_count,
            stage_evaluator_count: stage_evaluators
                .map(|evaluators| evaluators.stages.len() + 1)
                .unwrap_or(0),
            stage_evaluator_labels: stage_evaluators
                .map(|evaluators| {
                    evaluators
                        .stages
                        .iter()
                        .map(|stage| stage.evaluator_label.clone())
                        .collect()
                })
                .unwrap_or_default(),
            amplitude_evaluator_label: stage_evaluators
                .map(|evaluators| evaluators.amplitude_stage.evaluator_label.clone()),
            lc_topology_reuse_available: topology_reuse
                .and_then(|value| value.get("available"))
                .and_then(Value::as_bool)
                .unwrap_or(false),
            lc_topology_group_count: topology_reuse
                .and_then(|value| value.get("active_topology_group_count"))
                .and_then(Value::as_u64)
                .unwrap_or(0) as usize,
            lc_topology_representative_sector_ids: topology_reuse
                .and_then(|value| value.get("representative_sector_ids"))
                .and_then(Value::as_array)
                .map(|items| items.iter().filter_map(Value::as_i64).collect())
                .unwrap_or_default(),
            lc_topology_replay_enabled: !topology_replay_mappings.is_empty(),
            lc_topology_replay_sector_count: topology_replay_mappings.len(),
            lc_topology_replay_mappings: topology_replay_mappings,
            runtime_unavailable_message: manifest.compiled.runtime_unavailable_message,
            sources: manifest.runtime_schema.source_fill.sources,
            momentum_slots: manifest.runtime_schema.momentum_slots,
            external_is_initial,
            particle_masses,
            normalization_factor,
            stages: None,
            amplitude_stage: None,
            state_scratch_f64: Vec::new(),
            values_scratch_f64: Vec::new(),
        })
    }

    fn load_from_manifest(manifest: GenericProcessManifestV2, root: &Path) -> PyResult<Self> {
        let stage_evaluators = manifest.compiled.stage_evaluators.clone();
        let amplitude_stage_manifest = manifest.runtime_schema.amplitude_stage.clone();
        let mut runtime = Self::from_manifest(manifest)?;
        if let Some(stage_evaluators) = stage_evaluators {
            let stages = stage_evaluators
                .stages
                .iter()
                .map(|stage| GenericStageRuntimeV2::load(stage, root))
                .collect::<PyResult<Vec<_>>>()?;
            runtime.stages = Some(stages);
            runtime.amplitude_stage = Some(GenericAmplitudeRuntimeV2::load(
                &amplitude_stage_manifest,
                &stage_evaluators.amplitude_stage,
                root,
            )?);
        }
        Ok(runtime)
    }

    fn execution_unavailable_error(&self) -> PyErr {
        if self.stages.is_none() || self.amplitude_stage.is_none() {
            return PyValueError::new_err(format!(
                "generic DAG schema-v2 artifact for {} loaded successfully, but it has no \
                 serialized generic stage evaluators",
                self.process
            ));
        }
        let detail = self.runtime_unavailable_message.as_deref().unwrap_or(
            "generic source filling and staged evaluator execution from schema-v2 metadata remain pending",
        );
        PyValueError::new_err(format!(
            "generic DAG schema-v2 artifact for {} loaded successfully, but generic Rusticol \
             execution is not wired yet; {detail}",
            self.process,
        ))
    }

    fn run_f64(&mut self, batch: &[Vec<[f64; 4]>]) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        if self.lc_topology_replay_enabled {
            return self.run_f64_with_lc_topology_replay(batch);
        }
        self.run_f64_materialized(batch)
    }

    fn run_double(
        &mut self,
        batch: &[Vec<[DoubleFloat; 4]>],
    ) -> PyResult<(Vec<DoubleFloat>, RuntimeProfile)> {
        if self.lc_topology_replay_enabled {
            return self.run_generic_with_lc_topology_replay(batch, None);
        }
        self.run_generic_materialized(batch, None)
    }

    fn run_float(
        &mut self,
        batch: &[Vec<[Float; 4]>],
        binary_precision: u32,
    ) -> PyResult<(Vec<Float>, RuntimeProfile)> {
        if self.lc_topology_replay_enabled {
            return self.run_generic_with_lc_topology_replay(batch, Some(binary_precision));
        }
        self.run_generic_materialized(batch, Some(binary_precision))
    }

    fn run_f64_with_lc_topology_replay(
        &mut self,
        batch: &[Vec<[f64; 4]>],
    ) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        let total_start = Instant::now();
        let n_points = batch.len();
        let mut values = vec![0.0; n_points];
        let mut profile = RuntimeProfile::default();
        let mappings = self.lc_topology_replay_mappings.clone();
        let mappings_per_chunk = replay_mappings_per_expanded_batch(n_points);
        for mapping_chunk in mappings.chunks(mappings_per_chunk) {
            let expanded_batch =
                apply_lc_topology_label_permutations(batch, self.external_count, mapping_chunk)?;
            let (expanded_values, sector_profile) = self.run_f64_materialized(&expanded_batch)?;
            for mapping_index in 0..mapping_chunk.len() {
                let offset = mapping_index * n_points;
                for point_index in 0..n_points {
                    values[point_index] += expanded_values[offset + point_index];
                }
            }
            profile.source_fill_s += sector_profile.source_fill_s;
            profile.momentum_setup_s += sector_profile.momentum_setup_s;
            profile.stage_evaluator_s += sector_profile.stage_evaluator_s;
            profile.output_assign_s += sector_profile.output_assign_s;
            profile.amplitude_evaluator_s += sector_profile.amplitude_evaluator_s;
            profile.reduction_s += sector_profile.reduction_s;
        }
        profile.total_s = total_start.elapsed().as_secs_f64();
        Ok((values, profile))
    }

    fn run_f64_materialized(
        &mut self,
        batch: &[Vec<[f64; 4]>],
    ) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        if self.stages.is_none() || self.amplitude_stage.is_none() {
            return Err(self.execution_unavailable_error());
        }
        let total_start = Instant::now();
        let n_points = batch.len();
        let state_len = n_points * self.parameter_count;
        if self.state_scratch_f64.len() != state_len {
            self.state_scratch_f64.resize(state_len, c64(0.0, 0.0));
        } else {
            self.state_scratch_f64.fill(c64(0.0, 0.0));
        }
        let state = &mut self.state_scratch_f64;
        let sources = &self.sources;
        let momentum_slots = &self.momentum_slots;
        let external_count = self.external_count;
        let external_is_initial = &self.external_is_initial;
        let particle_masses = &self.particle_masses;
        let value_parameter_count = self.value_parameter_count;

        let source_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_sources_row(sources, external_count, particle_masses, row_state, point)?;
        }
        let source_fill_s = source_start.elapsed().as_secs_f64();

        let momentum_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_momenta_row(
                momentum_slots,
                value_parameter_count,
                external_count,
                external_is_initial,
                row_state,
                point,
            )?;
        }
        let momentum_setup_s = momentum_start.elapsed().as_secs_f64();

        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        for stage in self.stages.as_mut().expect("generic stages checked") {
            let (eval_s, assign_s) = stage.evaluate_f64_into_state(
                n_points,
                self.parameter_count,
                state.as_mut_slice(),
            )?;
            stage_evaluator_s += eval_s;
            output_assign_s += assign_s;
        }

        let amplitude_start = Instant::now();
        self.amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .evaluate_f64_into_scratch(n_points, state.as_slice())?;
        let amplitude_evaluator_s = amplitude_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        self.amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .reduce_scratch_f64_into(n_points, &mut self.values_scratch_f64)?;
        for value in &mut self.values_scratch_f64 {
            *value *= self.normalization_factor;
        }
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        Ok((
            self.values_scratch_f64.clone(),
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
            },
        ))
    }

    fn run_generic_with_lc_topology_replay<T>(
        &mut self,
        batch: &[Vec<[T; 4]>],
        binary_precision: Option<u32>,
    ) -> PyResult<(Vec<T>, RuntimeProfile)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let total_start = Instant::now();
        let n_points = batch.len();
        let mut values = vec![T::new_zero(); n_points];
        let mut profile = RuntimeProfile::default();
        let mappings = self.lc_topology_replay_mappings.clone();
        let mappings_per_chunk = replay_mappings_per_expanded_batch(n_points);
        for mapping_chunk in mappings.chunks(mappings_per_chunk) {
            let expanded_batch = apply_lc_topology_label_permutations_generic(
                batch,
                self.external_count,
                mapping_chunk,
            )?;
            let (expanded_values, sector_profile) =
                self.run_generic_materialized(&expanded_batch, binary_precision)?;
            for mapping_index in 0..mapping_chunk.len() {
                let offset = mapping_index * n_points;
                for point_index in 0..n_points {
                    values[point_index] += expanded_values[offset + point_index].clone();
                }
            }
            profile.source_fill_s += sector_profile.source_fill_s;
            profile.momentum_setup_s += sector_profile.momentum_setup_s;
            profile.stage_evaluator_s += sector_profile.stage_evaluator_s;
            profile.output_assign_s += sector_profile.output_assign_s;
            profile.amplitude_evaluator_s += sector_profile.amplitude_evaluator_s;
            profile.reduction_s += sector_profile.reduction_s;
        }
        profile.total_s = total_start.elapsed().as_secs_f64();
        Ok((values, profile))
    }

    fn run_generic_materialized<T>(
        &mut self,
        batch: &[Vec<[T; 4]>],
        binary_precision: Option<u32>,
    ) -> PyResult<(Vec<T>, RuntimeProfile)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        if self.stages.is_none() || self.amplitude_stage.is_none() {
            return Err(self.execution_unavailable_error());
        }
        let total_start = Instant::now();
        let n_points = batch.len();
        let mut state = vec![complex_zero::<T>(); n_points * self.parameter_count];
        let sources = &self.sources;
        let momentum_slots = &self.momentum_slots;
        let external_count = self.external_count;
        let external_is_initial = &self.external_is_initial;
        let particle_masses = &self.particle_masses;
        let value_parameter_count = self.value_parameter_count;

        let source_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_sources_row_generic(
                sources,
                external_count,
                particle_masses,
                row_state,
                point,
            )?;
        }
        let source_fill_s = source_start.elapsed().as_secs_f64();

        let momentum_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_momenta_row_generic(
                momentum_slots,
                value_parameter_count,
                external_count,
                external_is_initial,
                row_state,
                point,
            )?;
        }
        let momentum_setup_s = momentum_start.elapsed().as_secs_f64();

        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        for stage in self.stages.as_mut().expect("generic stages checked") {
            let (eval_s, assign_s) = stage.evaluate_generic_into_state(
                n_points,
                self.parameter_count,
                state.as_mut_slice(),
                binary_precision,
            )?;
            stage_evaluator_s += eval_s;
            output_assign_s += assign_s;
        }

        let amplitude_start = Instant::now();
        let raw_sums = self
            .amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .evaluate_raw_sums_generic(n_points, state.as_slice(), binary_precision)?;
        let amplitude_evaluator_s = amplitude_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        let factor = T::from(self.normalization_factor);
        let values = raw_sums
            .into_iter()
            .map(|value| value * factor.clone())
            .collect::<Vec<_>>();
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        Ok((
            values,
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
            },
        ))
    }

    fn evaluate_amplitudes_f64(&mut self, batch: &[Vec<[f64; 4]>]) -> PyResult<Vec<Complex<f64>>> {
        if self.lc_topology_replay_enabled {
            let mappings = self.lc_topology_replay_mappings.clone();
            let mut amplitudes = Vec::new();
            for mapping in &mappings {
                let mapped_batch =
                    apply_lc_topology_label_permutation(batch, self.external_count, mapping)?;
                amplitudes.extend(self.evaluate_amplitudes_f64_materialized(&mapped_batch)?);
            }
            return Ok(amplitudes);
        }
        self.evaluate_amplitudes_f64_materialized(batch)
    }

    fn evaluate_amplitudes_f64_materialized(
        &mut self,
        batch: &[Vec<[f64; 4]>],
    ) -> PyResult<Vec<Complex<f64>>> {
        if self.stages.is_none() || self.amplitude_stage.is_none() {
            return Err(self.execution_unavailable_error());
        }
        let n_points = batch.len();
        let state_len = n_points * self.parameter_count;
        if self.state_scratch_f64.len() != state_len {
            self.state_scratch_f64.resize(state_len, c64(0.0, 0.0));
        } else {
            self.state_scratch_f64.fill(c64(0.0, 0.0));
        }
        let state = &mut self.state_scratch_f64;
        let sources = &self.sources;
        let momentum_slots = &self.momentum_slots;
        let external_count = self.external_count;
        let external_is_initial = &self.external_is_initial;
        let particle_masses = &self.particle_masses;
        let value_parameter_count = self.value_parameter_count;

        for (row, point) in batch.iter().enumerate() {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_sources_row(sources, external_count, particle_masses, row_state, point)?;
        }
        for (row, point) in batch.iter().enumerate() {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_momenta_row(
                momentum_slots,
                value_parameter_count,
                external_count,
                external_is_initial,
                row_state,
                point,
            )?;
        }
        for stage in self.stages.as_mut().expect("generic stages checked") {
            let _ = stage.evaluate_f64_into_state(
                n_points,
                self.parameter_count,
                state.as_mut_slice(),
            )?;
        }
        let amplitude_stage = self
            .amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked");
        amplitude_stage.evaluate_f64_into_scratch(n_points, state.as_slice())?;
        Ok(amplitude_stage.output_scratch_f64.clone())
    }

    fn fill_sources_row(
        sources: &[GenericSourceRecordManifestV2],
        external_count: usize,
        particle_masses: &BTreeMap<i32, f64>,
        row: &mut [Complex<f64>],
        point: &[[f64; 4]],
    ) -> PyResult<()> {
        for source in sources {
            let start = source.value_slot.component_start;
            let stop = source.value_slot.component_stop;
            if stop > row.len() || stop < start {
                return Err(PyValueError::new_err(format!(
                    "generic source {} has an invalid value-slot range",
                    source.source_id
                )));
            }
            Self::write_source_wavefunction(
                source,
                external_count,
                particle_masses,
                point,
                &mut row[start..stop],
            )?;
        }
        Ok(())
    }

    fn fill_sources_row_generic<T>(
        sources: &[GenericSourceRecordManifestV2],
        external_count: usize,
        particle_masses: &BTreeMap<i32, f64>,
        row: &mut [Complex<T>],
        point: &[[T; 4]],
    ) -> PyResult<()>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        for source in sources {
            let start = source.value_slot.component_start;
            let stop = source.value_slot.component_stop;
            if stop > row.len() || stop < start {
                return Err(PyValueError::new_err(format!(
                    "generic source {} has an invalid value-slot range",
                    source.source_id
                )));
            }
            Self::write_source_wavefunction_generic(
                source,
                external_count,
                particle_masses,
                point,
                &mut row[start..stop],
            )?;
        }
        Ok(())
    }

    fn fill_momenta_row(
        momentum_slots: &[GenericMomentumSlotManifestV2],
        value_parameter_count: usize,
        external_count: usize,
        external_is_initial: &[bool],
        row: &mut [Complex<f64>],
        point: &[[f64; 4]],
    ) -> PyResult<()> {
        for slot in momentum_slots {
            let start = value_parameter_count + slot.component_start;
            let stop = value_parameter_count + slot.component_stop;
            if stop > row.len() || stop < start || stop - start != 4 {
                return Err(PyValueError::new_err(format!(
                    "generic momentum slot {} has an invalid component range",
                    slot.momentum_slot_id
                )));
            }
            let mut momentum = [0.0; 4];
            for label in &slot.external_labels {
                let index = label.checked_sub(1).ok_or_else(|| {
                    PyValueError::new_err("generic momentum labels are one-based")
                })?;
                if index >= external_count || index >= external_is_initial.len() {
                    return Err(PyValueError::new_err(format!(
                        "generic momentum slot {} refers to unknown external label {}",
                        slot.momentum_slot_id, label
                    )));
                }
                let sign = if external_is_initial[index] {
                    -1.0
                } else {
                    1.0
                };
                for component in 0..4 {
                    momentum[component] += sign * point[index][component];
                }
            }
            for component in 0..4 {
                row[start + component] = c64(momentum[component], 0.0);
            }
        }
        Ok(())
    }

    fn fill_momenta_row_generic<T>(
        momentum_slots: &[GenericMomentumSlotManifestV2],
        value_parameter_count: usize,
        external_count: usize,
        external_is_initial: &[bool],
        row: &mut [Complex<T>],
        point: &[[T; 4]],
    ) -> PyResult<()>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        for slot in momentum_slots {
            let start = value_parameter_count + slot.component_start;
            let stop = value_parameter_count + slot.component_stop;
            if stop > row.len() || stop < start || stop - start != 4 {
                return Err(PyValueError::new_err(format!(
                    "generic momentum slot {} has an invalid component range",
                    slot.momentum_slot_id
                )));
            }
            let mut momentum: [T; 4] = std::array::from_fn(|_| T::new_zero());
            for label in &slot.external_labels {
                let index = label.checked_sub(1).ok_or_else(|| {
                    PyValueError::new_err("generic momentum labels are one-based")
                })?;
                if index >= external_count || index >= external_is_initial.len() {
                    return Err(PyValueError::new_err(format!(
                        "generic momentum slot {} refers to unknown external label {}",
                        slot.momentum_slot_id, label
                    )));
                }
                for component in 0..4 {
                    if external_is_initial[index] {
                        momentum[component] -= point[index][component].clone();
                    } else {
                        momentum[component] += point[index][component].clone();
                    }
                }
            }
            for component in 0..4 {
                row[start + component] =
                    c_generic(momentum[component].clone(), T::new_zero());
            }
        }
        Ok(())
    }

    fn write_source_wavefunction(
        source: &GenericSourceRecordManifestV2,
        external_count: usize,
        particle_masses: &BTreeMap<i32, f64>,
        point: &[[f64; 4]],
        out: &mut [Complex<f64>],
    ) -> PyResult<()> {
        if source.source_kind != "external-wavefunction" {
            return Err(PyValueError::new_err(format!(
                "generic source kind {:?} is not implemented",
                source.source_kind
            )));
        }
        let index = source
            .leg_label
            .checked_sub(1)
            .ok_or_else(|| PyValueError::new_err("generic source leg labels are one-based"))?;
        if index >= external_count {
            return Err(PyValueError::new_err(format!(
                "generic source {} refers to unknown external label {}",
                source.source_id, source.leg_label
            )));
        }
        let momentum = if source.crossing == "negate-incoming-momentum" {
            negate(point[index])
        } else {
            point[index]
        };
        if source.dimension == 1 {
            if out.len() != 1 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 1 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            out[0] = c64(1.0, 0.0);
            return Ok(());
        }
        if source.dimension == 2 {
            if out.len() != 2 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 2 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            let chirality = source.chirality;
            let wave = if source.particle_id < 0 {
                ext_antiquark_weyl_array(momentum, source.source_helicity, chirality)
            } else {
                ext_quark_weyl_array(momentum, source.source_helicity, chirality)
            };
            out.copy_from_slice(&wave);
            return Ok(());
        }
        if source.dimension == 4 {
            if out.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 4 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            let wave = if is_fermion_pdg(source.particle_id) {
                let mass = particle_mass_from_map(particle_masses, source.particle_id);
                if source.particle_id < 0 {
                    ext_antiquark_dirac_massive(momentum, source.source_helicity, mass)
                } else {
                    ext_quark_dirac_massive(momentum, source.source_helicity, mass)
                }
            } else if source.particle_id.abs() == 21 || source.particle_id == 22 {
                ext_gluon(momentum, source.source_helicity)
            } else {
                ext_massive_vector(
                    momentum,
                    source.source_helicity,
                    particle_mass_from_map(particle_masses, source.particle_id),
                )
            };
            out.copy_from_slice(&wave);
            return Ok(());
        }
        Err(PyValueError::new_err(format!(
            "generic source dimension {} is not implemented",
            source.dimension
        )))
    }

    fn write_source_wavefunction_generic<T>(
        source: &GenericSourceRecordManifestV2,
        external_count: usize,
        particle_masses: &BTreeMap<i32, f64>,
        point: &[[T; 4]],
        out: &mut [Complex<T>],
    ) -> PyResult<()>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        if source.source_kind != "external-wavefunction" {
            return Err(PyValueError::new_err(format!(
                "generic source kind {:?} is not implemented",
                source.source_kind
            )));
        }
        let index = source
            .leg_label
            .checked_sub(1)
            .ok_or_else(|| PyValueError::new_err("generic source leg labels are one-based"))?;
        if index >= external_count {
            return Err(PyValueError::new_err(format!(
                "generic source {} refers to unknown external label {}",
                source.source_id, source.leg_label
            )));
        }
        let momentum = if source.crossing == "negate-incoming-momentum" {
            negate_generic(&point[index])
        } else {
            point[index].clone()
        };
        if source.dimension == 1 {
            if out.len() != 1 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 1 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            out[0] = c_generic(T::from(1.0), T::new_zero());
            return Ok(());
        }
        if source.dimension == 2 {
            if out.len() != 2 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 2 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            let chirality = source.chirality;
            let wave = if source.particle_id < 0 {
                ext_antiquark_weyl_generic(&momentum, source.source_helicity, chirality)
            } else {
                ext_quark_weyl_generic(&momentum, source.source_helicity, chirality)
            };
            out.clone_from_slice(&wave);
            return Ok(());
        }
        if source.dimension == 4 {
            if out.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 4 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            let wave = if is_fermion_pdg(source.particle_id) {
                let mass = particle_mass_from_map(particle_masses, source.particle_id);
                if mass != 0.0 {
                    return Err(PyValueError::new_err(
                        "high-precision generic massive fermion sources are not implemented",
                    ));
                }
                if source.particle_id < 0 {
                    ext_antiquark_dirac_generic(&momentum, source.source_helicity)
                } else {
                    ext_quark_dirac_generic(&momentum, source.source_helicity)
                }
            } else if source.particle_id.abs() == 21 || source.particle_id == 22 {
                ext_gluon_generic(&momentum, source.source_helicity)
            } else {
                ext_massive_vector_generic(
                    &momentum,
                    source.source_helicity,
                    T::from(particle_mass_from_map(particle_masses, source.particle_id)),
                )
            };
            out.clone_from_slice(&wave);
            return Ok(());
        }
        Err(PyValueError::new_err(format!(
            "generic source dimension {} is not implemented",
            source.dimension
        )))
    }
}

fn particle_mass_from_map(particle_masses: &BTreeMap<i32, f64>, particle_id: i32) -> f64 {
    particle_masses
        .get(&particle_id)
        .copied()
        .or_else(|| particle_masses.get(&(-particle_id)).copied())
        .unwrap_or(0.0)
}

fn is_fermion_pdg(particle_id: i32) -> bool {
    let abs_id = particle_id.abs();
    (1..=6).contains(&abs_id) || (11..=16).contains(&abs_id)
}

#[pyclass(module = "rusticol")]
struct Runtime {
    root: PathBuf,
    manifest: Option<ProcessManifest>,
    generic_runtime: Option<GenericRuntimeV2>,
    selected_process_key: Option<String>,
    selected_process: Option<String>,
    input_crossing_map: Option<Vec<InputCrossingMapEntry>>,
    crossing_alias_of: Option<String>,
    stages: Option<Vec<CurrentStage>>,
    amplitude_stage: Option<AmplitudeStage>,
    zero_gluon_stage: Option<ZeroGluonStage>,
    last_profile: RuntimeProfile,
}

struct ProcessSetSelection {
    root: PathBuf,
    selected_key: Option<String>,
    selected_process: Option<String>,
    input_crossing_map: Option<Vec<InputCrossingMapEntry>>,
    crossing_alias_of: Option<String>,
}

fn load_rusticol_runtime(
    process_dir: &str,
    process_key: Option<&str>,
    allow_legacy_schema_v1: bool,
) -> PyResult<Runtime> {
    let process_dir = PathBuf::from(process_dir);
    let requested_root = process_dir.canonicalize().map_err(|err| {
        PyValueError::new_err(format!(
            "could not resolve process directory {}: {err}",
            process_dir.display()
        ))
    })?;
    let selection = resolve_process_root(&requested_root, process_key)?;
    let root = selection.root.clone();
    let manifest_path = root.join("process_manifest.json");
    let manifest_text = fs::read_to_string(&manifest_path).map_err(|err| {
        PyValueError::new_err(format!(
            "could not read process manifest {}: {err}",
            manifest_path.display()
        ))
    })?;
    let manifest_value: Value = serde_json::from_str(&manifest_text).map_err(|err| {
        PyValueError::new_err(format!(
            "could not parse process manifest {}: {err}",
            manifest_path.display()
        ))
    })?;
    if let Some(generic_runtime) = load_generic_schema_v2_manifest(&manifest_value, &root)? {
        return Ok(Runtime {
            root,
            manifest: None,
            generic_runtime: Some(generic_runtime),
            selected_process_key: selection.selected_key,
            selected_process: selection.selected_process,
            input_crossing_map: selection.input_crossing_map,
            crossing_alias_of: selection.crossing_alias_of,
            stages: None,
            amplitude_stage: None,
            zero_gluon_stage: None,
            last_profile: RuntimeProfile::default(),
        });
    }

    if !allow_legacy_schema_v1 {
        let schema_version = manifest_value
            .get("schema_version")
            .and_then(Value::as_u64)
            .unwrap_or(0);
        let kind = manifest_value
            .get("kind")
            .and_then(Value::as_str)
            .unwrap_or("<missing>");
        return Err(PyValueError::new_err(format!(
            "Runtime.load only supports schema-v2 generic DAG artifacts; got kind {kind:?} \
             schema_version {schema_version}. Schema-v1 Rusticol artifacts are legacy; use \
             Runtime.load_legacy(...) for reference artifacts or regenerate with pyAmpliCol \
             generate-process."
        )));
    }

    let manifest: ProcessManifest = serde_json::from_value(manifest_value).map_err(|err| {
        PyValueError::new_err(format!(
            "could not parse process manifest {}: {err}",
            manifest_path.display()
        ))
    })?;
    if manifest.schema_version != 1 {
        return Err(PyValueError::new_err(format!(
            "unsupported process manifest schema_version {}",
            manifest.schema_version
        )));
    }
    if manifest.kind != "pyamplicol-rusticol-process" {
        return Err(PyValueError::new_err(format!(
            "unsupported process artifact kind {}",
            manifest.kind
        )));
    }
    if manifest.family != "q-qbar-z-gluons-leading-color"
        && manifest.family != "q-qbar-vector-gluons-leading-color"
        && manifest.family != "q-qbar-neutral-dilepton-gluons-leading-color"
        && manifest.family != "q-qbar-neutral-dilepton-zero-gluon-leading-color"
        && manifest.family != "q-qbar-charged-leptonic-w-gluons-leading-color"
        && manifest.family != "q-qbar-charged-leptonic-w-zero-gluon-leading-color"
    {
        return Err(PyValueError::new_err(format!(
            "rusticol currently supports q-qbar vector/leptonic-vector+gluon leading-colour artifacts, got {}",
            manifest.family
        )));
    }
    if manifest.compiled.kind == "zero-gluon-symbolic-scalar" {
        if manifest.gluon_count != 0 || manifest.external_pdg_order.len() != 3 {
            return Err(PyValueError::new_err(
                "zero-gluon process artifacts must have exactly three external legs",
            ));
        }
        if manifest.external_pdg_order[0] + manifest.external_pdg_order[1] != 0
            || manifest.external_pdg_order[2] != 23
        {
            return Err(PyValueError::new_err(
                "zero-gluon process artifacts must be q qbar -> Z",
            ));
        }
        if !manifest.normalization.coupling_factor.is_finite()
            || !manifest.model.alpha_ew.is_finite()
        {
            return Err(PyValueError::new_err(
                "model normalization constants are not finite",
            ));
        }
        let zero_manifest = manifest.compiled.zero_gluon.as_ref().ok_or_else(|| {
            PyValueError::new_err("zero-gluon process artifact is missing zero_gluon")
        })?;
        let zero_gluon_stage = ZeroGluonStage::load(zero_manifest, &root)?;
        return Ok(Runtime {
            root,
            manifest: Some(manifest),
            generic_runtime: None,
            selected_process_key: selection.selected_key,
            selected_process: selection.selected_process,
            input_crossing_map: selection.input_crossing_map,
            crossing_alias_of: selection.crossing_alias_of,
            stages: Some(Vec::new()),
            amplitude_stage: None,
            zero_gluon_stage: Some(zero_gluon_stage),
            last_profile: RuntimeProfile::default(),
        });
    }
    if manifest.compiled.kind != "shared-compiled-sweep" {
        return Err(PyValueError::new_err(format!(
            "unsupported compiled sweep kind {}",
            manifest.compiled.kind
        )));
    }
    let has_composite_sources = manifest
        .table
        .sources
        .iter()
        .any(|source| source.source_kind != "external");
    let expected_external_count = if has_composite_sources {
        manifest.gluon_count + 4
    } else {
        manifest.gluon_count + 3
    };
    if manifest.external_pdg_order.len() != expected_external_count {
        return Err(PyValueError::new_err(
            "external PDG order length does not match the shared-current artifact layout",
        ));
    }
    if manifest
        .table
        .currents
        .iter()
        .enumerate()
        .any(|(index, current)| current.id != index || current.dimension > 6 || current.pdg == 0)
    {
        return Err(PyValueError::new_err(
            "current table ids/dimensions are not compatible with rusticol",
        ));
    }
    validate_source_metadata(&manifest.table.sources)?;
    if manifest.table.amplitudes.len() != manifest.compiled.amplitude_stage.raw_sum_weights.len() {
        return Err(PyValueError::new_err(
            "amplitude metadata and raw-sum weights have different lengths",
        ));
    }
    if manifest
        .table
        .amplitudes
        .iter()
        .zip(&manifest.compiled.amplitude_stage.raw_sum_weights)
        .any(|(amplitude, weight)| (amplitude.multiplicity - *weight).abs() > 0.0)
    {
        return Err(PyValueError::new_err(
            "amplitude multiplicities and raw-sum weights differ",
        ));
    }
    validate_amplitude_group_metadata(
        &manifest.table.amplitudes,
        manifest
            .compiled
            .amplitude_stage
            .raw_sum_group_ids
            .as_deref(),
    )?;
    if !manifest.normalization.coupling_factor.is_finite()
        || !manifest.model.alpha_s_me_check.is_finite()
        || !manifest.model.alpha_ew.is_finite()
    {
        return Err(PyValueError::new_err(
            "model normalization constants are not finite",
        ));
    }
    let stages = manifest
        .compiled
        .stages
        .iter()
        .map(|stage| CurrentStage::load(stage, &root, &manifest.layout.current_offsets))
        .collect::<PyResult<Vec<_>>>()?;
    let amplitude_stage = AmplitudeStage::load(&manifest.compiled.amplitude_stage, &root)?;
    Ok(Runtime {
        root,
        manifest: Some(manifest),
        generic_runtime: None,
        selected_process_key: selection.selected_key,
        selected_process: selection.selected_process,
        input_crossing_map: selection.input_crossing_map,
        crossing_alias_of: selection.crossing_alias_of,
        stages: Some(stages),
        amplitude_stage: Some(amplitude_stage),
        zero_gluon_stage: None,
        last_profile: RuntimeProfile::default(),
    })
}

#[pymethods]
impl Runtime {
    #[classmethod]
    #[pyo3(signature = (process_dir, process_key=None))]
    fn load(
        _cls: &Bound<'_, pyo3::types::PyType>,
        process_dir: &str,
        process_key: Option<&str>,
    ) -> PyResult<Self> {
        load_rusticol_runtime(process_dir, process_key, false)
    }

    #[classmethod]
    #[pyo3(signature = (process_dir, process_key=None))]
    fn load_legacy(
        _cls: &Bound<'_, pyo3::types::PyType>,
        process_dir: &str,
        process_key: Option<&str>,
    ) -> PyResult<Self> {
        load_rusticol_runtime(process_dir, process_key, true)
    }

    #[getter]
    fn process(&self) -> PyResult<String> {
        if let Some(generic) = &self.generic_runtime {
            if let Some(selected) = &self.selected_process {
                return Ok(selected.clone());
            }
            return Ok(generic.process.clone());
        }
        Ok(self.legacy_manifest().process.clone())
    }

    #[getter]
    fn gluon_count(&self) -> PyResult<usize> {
        if self.generic_runtime.is_some() {
            return Err(PyValueError::new_err(
                "generic schema-v2 artifacts do not expose a family-level gluon_count",
            ));
        }
        Ok(self.legacy_manifest().gluon_count)
    }

    fn metadata<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        if let Some(generic) = &self.generic_runtime {
            dict.set_item(
                "process",
                self.selected_process.as_ref().unwrap_or(&generic.process),
            )?;
            dict.set_item(
                "key",
                self.selected_process_key.as_ref().unwrap_or(&generic.key),
            )?;
            dict.set_item("representative_process", &generic.process)?;
            dict.set_item("representative_key", &generic.key)?;
            dict.set_item("crossing_alias_of", self.crossing_alias_of.clone())?;
            dict.set_item(
                "input_crossing_map",
                self.input_crossing_map
                    .as_ref()
                    .map(input_crossing_map_to_json)
                    .unwrap_or_default(),
            )?;
            dict.set_item("artifact_class", "generic-dag-schema-v2")?;
            dict.set_item("schema_version", 2)?;
            dict.set_item("color_accuracy", &generic.color_accuracy)?;
            dict.set_item("external_pdg_order", generic.external_pdg_order.clone())?;
            dict.set_item("external_count", generic.external_count)?;
            dict.set_item("parameter_count", generic.parameter_count)?;
            dict.set_item("value_parameter_count", generic.value_parameter_count)?;
            dict.set_item("momentum_parameter_count", generic.momentum_parameter_count)?;
            dict.set_item("current_count", generic.current_count)?;
            dict.set_item("source_count", generic.source_count)?;
            dict.set_item("interaction_count", generic.interaction_count)?;
            dict.set_item("stage_count", generic.stage_count)?;
            dict.set_item("amplitude_output_count", generic.amplitude_output_count)?;
            dict.set_item("stage_evaluator_count", generic.stage_evaluator_count)?;
            dict.set_item(
                "stage_evaluator_labels",
                generic.stage_evaluator_labels.clone(),
            )?;
            dict.set_item(
                "amplitude_evaluator_label",
                generic.amplitude_evaluator_label.clone(),
            )?;
            dict.set_item(
                "lc_topology_reuse_available",
                generic.lc_topology_reuse_available,
            )?;
            dict.set_item("lc_topology_group_count", generic.lc_topology_group_count)?;
            dict.set_item(
                "lc_topology_representative_sector_ids",
                generic.lc_topology_representative_sector_ids.clone(),
            )?;
            dict.set_item(
                "lc_topology_replay_enabled",
                generic.lc_topology_replay_enabled,
            )?;
            dict.set_item(
                "lc_topology_replay_sector_count",
                generic.lc_topology_replay_sector_count,
            )?;
            dict.set_item(
                "runtime_available",
                generic.stages.is_some() && generic.amplitude_stage.is_some(),
            )?;
            dict.set_item("root", self.root.to_string_lossy().to_string())?;
            return Ok(dict);
        }
        dict.set_item("process", &self.legacy_manifest().process)?;
        dict.set_item("family", &self.legacy_manifest().family)?;
        dict.set_item("schema_version", 1)?;
        dict.set_item("gluon_count", self.legacy_manifest().gluon_count)?;
        dict.set_item(
            "parameter_count",
            self.legacy_manifest().layout.parameter_count,
        )?;
        dict.set_item("stage_count", self.stages_ref().len())?;
        dict.set_item("root", self.root.to_string_lossy().to_string())?;
        Ok(dict)
    }

    fn evaluate<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let values = self.evaluate_f64_values(py, momenta)?;
        f64_values_to_numpy_or_list(py, values)
    }

    fn raw_amplitudes<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let external_count = self
            .generic_runtime
            .as_ref()
            .ok_or_else(|| {
                PyValueError::new_err(
                    "raw_amplitudes is only available for schema-v2 generic DAG artifacts",
                )
            })?
            .external_count;
        let batch = self.generic_batch_from_python(py, momenta, external_count)?;
        let generic = self.generic_runtime.as_mut().ok_or_else(|| {
            PyValueError::new_err(
                "raw_amplitudes is only available for schema-v2 generic DAG artifacts",
            )
        })?;
        if generic.lc_topology_replay_enabled {
            return Err(PyValueError::new_err(
                "raw_amplitudes is not available for LC topology replay artifacts yet; \
                 use a fully materialized sector artifact or evaluate()",
            ));
        }
        let output_len = generic.amplitude_output_count;
        let amplitudes = generic.evaluate_amplitudes_f64(&batch)?;
        let rows = amplitudes
            .chunks(output_len)
            .map(|row| {
                row.iter()
                    .map(|value| (value.re, value.im))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        rows.into_py_any(py)
    }

    fn evaluate_with_prec<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
        decimal_digit_precision: u32,
    ) -> PyResult<Py<PyAny>> {
        match PrecisionMode::from_decimal_digits(decimal_digit_precision)? {
            PrecisionMode::F64 => self.evaluate(py, momenta),
            PrecisionMode::DoubleDouble => {
                let (values, profile) = if self.generic_runtime.is_some() {
                    let external_count = self
                        .generic_runtime
                        .as_ref()
                        .expect("checked")
                        .external_count;
                    let batch = self.generic_batch_double_from_python(momenta, external_count)?;
                    self.generic_runtime
                        .as_mut()
                        .expect("checked")
                        .run_double(&batch)?
                } else {
                    let batch = batch_momenta_double(
                        momenta,
                        self.legacy_manifest().external_pdg_order.len(),
                    )?;
                    self.run_double(&batch)?
                };
                self.last_profile = profile;
                decimals_to_python(py, values, decimal_digit_precision)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let (values, profile) = if self.generic_runtime.is_some() {
                    let external_count = self
                        .generic_runtime
                        .as_ref()
                        .expect("checked")
                        .external_count;
                    let batch = self.generic_batch_float_from_python(
                        momenta,
                        binary_precision,
                        external_count,
                    )?;
                    self.generic_runtime
                        .as_mut()
                        .expect("checked")
                        .run_float(&batch, binary_precision)?
                } else {
                    let batch = batch_momenta_float(
                        momenta,
                        binary_precision,
                        self.legacy_manifest().external_pdg_order.len(),
                    )?;
                    self.run_float(&batch, binary_precision)?
                };
                self.last_profile = profile;
                decimals_to_python(py, values, decimal_precision)
            }
        }
    }

    #[pyo3(signature = (momenta, precision = 16, include_values = true))]
    fn profile<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
        precision: u32,
        include_values: bool,
    ) -> PyResult<Bound<'py, PyDict>> {
        let (points, values, profile) = match PrecisionMode::from_decimal_digits(precision)? {
            PrecisionMode::F64 => {
                let (points, values, profile) = if self.generic_runtime.is_some() {
                    let external_count = self
                        .generic_runtime
                        .as_ref()
                        .expect("checked")
                        .external_count;
                    let batch = self.generic_batch_from_python(py, momenta, external_count)?;
                    let points = batch.len();
                    let (values, profile) = self
                        .generic_runtime
                        .as_mut()
                        .expect("checked")
                        .run_f64(&batch)?;
                    (points, values, profile)
                } else {
                    let batch = batch_momenta(
                        py,
                        momenta,
                        self.legacy_manifest().external_pdg_order.len(),
                    )?;
                    let points = batch.len();
                    let (values, profile) = self.run_f64(&batch)?;
                    (points, values, profile)
                };
                let values = if include_values {
                    values.into_py_any(py)?
                } else {
                    py.None()
                };
                (points, values, profile)
            }
            PrecisionMode::DoubleDouble => {
                let (points, values, profile) = if self.generic_runtime.is_some() {
                    let external_count = self
                        .generic_runtime
                        .as_ref()
                        .expect("checked")
                        .external_count;
                    let batch = self.generic_batch_double_from_python(momenta, external_count)?;
                    let points = batch.len();
                    let (values, profile) = self
                        .generic_runtime
                        .as_mut()
                        .expect("checked")
                        .run_double(&batch)?;
                    (points, values, profile)
                } else {
                    let batch = batch_momenta_double(
                        momenta,
                        self.legacy_manifest().external_pdg_order.len(),
                    )?;
                    let points = batch.len();
                    let (values, profile) = self.run_double(&batch)?;
                    (points, values, profile)
                };
                let values = if include_values {
                    decimals_to_python(py, values, precision)?
                } else {
                    py.None()
                };
                (points, values, profile)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let (points, values, profile) = if self.generic_runtime.is_some() {
                    let external_count = self
                        .generic_runtime
                        .as_ref()
                        .expect("checked")
                        .external_count;
                    let batch = self.generic_batch_float_from_python(
                        momenta,
                        binary_precision,
                        external_count,
                    )?;
                    let points = batch.len();
                    let (values, profile) = self
                        .generic_runtime
                        .as_mut()
                        .expect("checked")
                        .run_float(&batch, binary_precision)?;
                    (points, values, profile)
                } else {
                    let batch = batch_momenta_float(
                        momenta,
                        binary_precision,
                        self.legacy_manifest().external_pdg_order.len(),
                    )?;
                    let points = batch.len();
                    let (values, profile) = self.run_float(&batch, binary_precision)?;
                    (points, values, profile)
                };
                (
                    points,
                    if include_values {
                        decimals_to_python(py, values, decimal_precision)?
                    } else {
                        py.None()
                    },
                    profile,
                )
            }
        };
        let dict = PyDict::new(py);
        let memory = memory_snapshot();
        dict.set_item("precision", precision)?;
        dict.set_item("points", points)?;
        dict.set_item("batch_size", points)?;
        dict.set_item("values", values)?;
        dict.set_item("source_fill_time_s", profile.source_fill_s)?;
        dict.set_item("momentum_setup_time_s", profile.momentum_setup_s)?;
        dict.set_item("parameter_pack_time_s", 0.0)?;
        dict.set_item("stage_evaluator_time_s", profile.stage_evaluator_s)?;
        dict.set_item("output_transfer_time_s", profile.output_assign_s)?;
        dict.set_item("output_assign_time_s", profile.output_assign_s)?;
        dict.set_item("amplitude_evaluator_time_s", profile.amplitude_evaluator_s)?;
        dict.set_item("result_reduction_time_s", profile.reduction_s)?;
        dict.set_item("total_time_s", profile.total_s)?;
        dict.set_item("current_rss_bytes", memory.current_rss_bytes)?;
        dict.set_item("peak_rss_bytes", memory.peak_rss_bytes)?;
        dict.set_item(
            "core_evaluator_time_s",
            profile.stage_evaluator_s + profile.amplitude_evaluator_s,
        )?;
        Ok(dict)
    }

    fn stage_diagnostics<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyDict>> {
        self.ensure_legacy_runtime()?;
        let batch = batch_momenta(py, momenta, self.legacy_manifest().external_pdg_order.len())?;
        let n_points = batch.len();
        let parameter_count = self.legacy_manifest().layout.parameter_count;
        let mut state = vec![c64(0.0, 0.0); n_points * parameter_count];

        for (row, point) in batch.iter().enumerate() {
            self.fill_sources(
                &mut state[row * parameter_count..(row + 1) * parameter_count],
                point,
            )?;
        }
        for (row, point) in batch.iter().enumerate() {
            self.fill_momenta(
                &mut state[row * parameter_count..(row + 1) * parameter_count],
                point,
            );
        }

        let stages = PyList::empty(py);
        stages.append(checksum_dict_str(
            py,
            "sources",
            checksum_state_prefix(
                &state,
                n_points,
                parameter_count,
                self.legacy_manifest().table.currents.len() * 6,
            ),
        )?)?;

        let current_offsets = self.legacy_manifest().layout.current_offsets.clone();
        let mut evaluated = Vec::new();
        let stage_count = self.stages_ref().len();
        for stage_index in 0..stage_count {
            let checksum = {
                let stage = &mut self.stages_mut()[stage_index];
                stage
                    .evaluator
                    .evaluate_batch_into(n_points, &state, &mut evaluated)?;
                for row in 0..n_points {
                    let row_state = &mut state[row * parameter_count..(row + 1) * parameter_count];
                    for (column, _, _, state_offset) in &stage.outputs {
                        row_state[*state_offset] =
                            evaluated[row * stage.evaluator.output_len + *column];
                    }
                }
                checksum_stage_outputs(
                    &state,
                    n_points,
                    parameter_count,
                    &current_offsets,
                    &stage.outputs,
                )
            };
            stages.append(checksum_dict_usize(py, stage_index + 1, checksum)?)?;
        }

        let dict = PyDict::new(py);
        dict.set_item("points", n_points)?;
        dict.set_item("stages", stages)?;
        Ok(dict)
    }
}

enum PrecisionMode {
    F64,
    DoubleDouble,
    Arbitrary(u32),
}

impl PrecisionMode {
    fn from_decimal_digits(precision: u32) -> PyResult<Self> {
        if precision == 16 {
            Ok(Self::F64)
        } else if precision == 32 {
            Ok(Self::DoubleDouble)
        } else if precision > 0 {
            Ok(Self::Arbitrary(precision))
        } else {
            Err(PyValueError::new_err("precision must be positive"))
        }
    }
}

impl Runtime {
    fn generic_batch_from_python(
        &self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
        expected_legs: usize,
    ) -> PyResult<Vec<Vec<[f64; 4]>>> {
        let batch = batch_momenta_dynamic(py, momenta, expected_legs)?;
        apply_input_crossing_map(&batch, expected_legs, self.input_crossing_map.as_deref())
    }

    fn generic_batch_double_from_python(
        &self,
        momenta: &Bound<'_, PyAny>,
        expected_legs: usize,
    ) -> PyResult<Vec<Vec<[DoubleFloat; 4]>>> {
        let batch = batch_momenta_double(momenta, expected_legs)?;
        apply_input_crossing_map_generic(
            &batch,
            expected_legs,
            self.input_crossing_map.as_deref(),
        )
    }

    fn generic_batch_float_from_python(
        &self,
        momenta: &Bound<'_, PyAny>,
        binary_precision: u32,
        expected_legs: usize,
    ) -> PyResult<Vec<Vec<[Float; 4]>>> {
        let batch = batch_momenta_float(momenta, binary_precision, expected_legs)?;
        apply_input_crossing_map_generic(
            &batch,
            expected_legs,
            self.input_crossing_map.as_deref(),
        )
    }

    fn legacy_manifest(&self) -> &ProcessManifest {
        self.manifest
            .as_ref()
            .expect("legacy schema-v1 manifest is unavailable for generic runtime")
    }

    fn ensure_legacy_runtime(&self) -> PyResult<()> {
        if let Some(generic) = &self.generic_runtime {
            return Err(generic.execution_unavailable_error());
        }
        Ok(())
    }

    fn stages_ref(&self) -> &Vec<CurrentStage> {
        self.stages
            .as_ref()
            .expect("rusticol runtime stages are unavailable during drop")
    }

    fn stages_mut(&mut self) -> &mut Vec<CurrentStage> {
        self.stages
            .as_mut()
            .expect("rusticol runtime stages are unavailable during drop")
    }

    fn amplitude_stage_mut(&mut self) -> &mut AmplitudeStage {
        self.amplitude_stage
            .as_mut()
            .expect("rusticol amplitude stage is unavailable during drop")
    }

    fn zero_gluon_stage_ref(&self) -> &ZeroGluonStage {
        self.zero_gluon_stage
            .as_ref()
            .expect("rusticol zero-gluon stage is unavailable during drop")
    }

    fn zero_gluon_stage_mut(&mut self) -> &mut ZeroGluonStage {
        self.zero_gluon_stage
            .as_mut()
            .expect("rusticol zero-gluon stage is unavailable during drop")
    }

    fn evaluate_f64_values(
        &mut self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<f64>> {
        if self.generic_runtime.is_some() {
            let external_count = self
                .generic_runtime
                .as_ref()
                .expect("checked")
                .external_count;
            let batch = self.generic_batch_from_python(py, momenta, external_count)?;
            let (values, profile) = self
                .generic_runtime
                .as_mut()
                .expect("checked")
                .run_f64(&batch)?;
            self.last_profile = profile;
            return Ok(values);
        }
        let batch = batch_momenta(py, momenta, self.legacy_manifest().external_pdg_order.len())?;
        let (values, profile) = self.run_f64(&batch)?;
        self.last_profile = profile;
        Ok(values)
    }

    fn run_f64(&mut self, batch: &[[[f64; 4]; 16]]) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        if self.zero_gluon_stage.is_some() {
            return self.run_zero_gluon_f64(batch);
        }
        let total_start = Instant::now();
        let n_points = batch.len();
        let parameter_count = self.legacy_manifest().layout.parameter_count;
        let mut state = vec![c64(0.0, 0.0); n_points * parameter_count];

        let source_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            self.fill_sources(
                &mut state[row * parameter_count..(row + 1) * parameter_count],
                point,
            )?;
        }
        let source_fill_s = source_start.elapsed().as_secs_f64();

        let momentum_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            self.fill_momenta(
                &mut state[row * parameter_count..(row + 1) * parameter_count],
                point,
            );
        }
        let momentum_setup_s = momentum_start.elapsed().as_secs_f64();

        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        for stage in self.stages_mut() {
            let (eval_s, assign_s) =
                stage.evaluate_f64_into_state(n_points, parameter_count, &mut state)?;
            stage_evaluator_s += eval_s;
            output_assign_s += assign_s;
        }

        let amp_start = Instant::now();
        let raw_sums = self
            .amplitude_stage_mut()
            .evaluate_raw_sums(n_points, &state)?;
        let amplitude_evaluator_s = amp_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        let factor = self.legacy_manifest().normalization.color_factor
            * self.legacy_manifest().normalization.coupling_factor
            / (self.legacy_manifest().normalization.average_factor
                * self.legacy_manifest().normalization.identical_factor);
        let values = raw_sums
            .into_iter()
            .map(|raw| raw * factor)
            .collect::<Vec<_>>();
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        let total_s = total_start.elapsed().as_secs_f64();

        Ok((
            values,
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s,
            },
        ))
    }

    fn run_double(
        &mut self,
        batch: &[Vec<[DoubleFloat; 4]>],
    ) -> PyResult<(Vec<DoubleFloat>, RuntimeProfile)> {
        self.run_generic(batch, None)
    }

    fn run_float(
        &mut self,
        batch: &[Vec<[Float; 4]>],
        binary_precision: u32,
    ) -> PyResult<(Vec<Float>, RuntimeProfile)> {
        self.run_generic(batch, Some(binary_precision))
    }

    fn run_generic<T>(
        &mut self,
        batch: &[Vec<[T; 4]>],
        binary_precision: Option<u32>,
    ) -> PyResult<(Vec<T>, RuntimeProfile)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        if self.zero_gluon_stage.is_some() {
            return self.run_zero_gluon_generic(batch, binary_precision);
        }
        let total_start = Instant::now();
        let n_points = batch.len();
        let parameter_count = self.legacy_manifest().layout.parameter_count;
        let mut state = vec![complex_zero::<T>(); n_points * parameter_count];

        let source_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            self.fill_sources_generic(
                &mut state[row * parameter_count..(row + 1) * parameter_count],
                point,
            )?;
        }
        let source_fill_s = source_start.elapsed().as_secs_f64();

        let momentum_start = Instant::now();
        for (row, point) in batch.iter().enumerate() {
            self.fill_momenta_generic(
                &mut state[row * parameter_count..(row + 1) * parameter_count],
                point,
            );
        }
        let momentum_setup_s = momentum_start.elapsed().as_secs_f64();

        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        for stage in self.stages_mut() {
            let eval_start = Instant::now();
            let evaluated =
                stage
                    .evaluator
                    .evaluate_batch_generic(n_points, &state, binary_precision)?;
            stage_evaluator_s += eval_start.elapsed().as_secs_f64();
            let assign_start = Instant::now();
            for row in 0..n_points {
                let row_state = &mut state[row * parameter_count..(row + 1) * parameter_count];
                for (column, _, _, state_offset) in &stage.outputs {
                    row_state[*state_offset] =
                        evaluated[row * stage.evaluator.output_len + *column].clone();
                }
            }
            output_assign_s += assign_start.elapsed().as_secs_f64();
        }

        let amp_start = Instant::now();
        let raw_sums = self.amplitude_stage_mut().evaluate_raw_sums_generic(
            n_points,
            &state,
            binary_precision,
        )?;
        let amplitude_evaluator_s = amp_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        let factor = T::from(
            self.legacy_manifest().normalization.color_factor
                * self.legacy_manifest().normalization.coupling_factor
                / (self.legacy_manifest().normalization.average_factor
                    * self.legacy_manifest().normalization.identical_factor),
        );
        let values = raw_sums
            .into_iter()
            .map(|raw| raw * factor.clone())
            .collect::<Vec<_>>();
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        let total_s = total_start.elapsed().as_secs_f64();

        Ok((
            values,
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s,
            },
        ))
    }

    fn run_zero_gluon_f64(
        &mut self,
        batch: &[[[f64; 4]; 16]],
    ) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        let total_start = Instant::now();
        let source_start = Instant::now();
        let params = self
            .zero_gluon_stage_ref()
            .parameter_rows_f64(self.legacy_manifest(), batch)?;
        let source_fill_s = source_start.elapsed().as_secs_f64();

        let eval_start = Instant::now();
        let output_len = self.zero_gluon_stage_ref().evaluator.output_len;
        let mut evaluated = vec![c64(0.0, 0.0); batch.len() * output_len];
        {
            let stage = self.zero_gluon_stage_mut();
            stage
                .evaluator
                .evaluate_f64_batch(batch.len(), &params, &mut evaluated)?;
        }
        let amplitude_evaluator_s = eval_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        let factor = self.legacy_manifest().normalization.color_factor
            * self.legacy_manifest().normalization.coupling_factor
            / (self.legacy_manifest().normalization.average_factor
                * self.legacy_manifest().normalization.identical_factor);
        let values = (0..batch.len())
            .map(|row| evaluated[row * output_len].re * factor)
            .collect::<Vec<_>>();
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        Ok((
            values,
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s: 0.0,
                stage_evaluator_s: 0.0,
                output_assign_s: 0.0,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
            },
        ))
    }

    fn run_zero_gluon_generic<T>(
        &mut self,
        batch: &[Vec<[T; 4]>],
        binary_precision: Option<u32>,
    ) -> PyResult<(Vec<T>, RuntimeProfile)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let total_start = Instant::now();
        let source_start = Instant::now();
        let params = self
            .zero_gluon_stage_ref()
            .parameter_rows_generic(self.legacy_manifest(), batch)?;
        let source_fill_s = source_start.elapsed().as_secs_f64();

        let eval_start = Instant::now();
        let output_len = self.zero_gluon_stage_ref().evaluator.output_len;
        let mut evaluated = vec![complex_zero::<T>(); batch.len() * output_len];
        {
            let stage = self.zero_gluon_stage_mut();
            let input_len = stage.evaluator.input_len;
            for row in 0..batch.len() {
                let in_start = row * input_len;
                let out_start = row * output_len;
                T::evaluate_loaded(
                    &mut stage.evaluator,
                    &params[in_start..in_start + input_len],
                    &mut evaluated[out_start..out_start + output_len],
                    binary_precision,
                )?;
            }
        }
        let amplitude_evaluator_s = eval_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        let factor = T::from(
            self.legacy_manifest().normalization.color_factor
                * self.legacy_manifest().normalization.coupling_factor
                / (self.legacy_manifest().normalization.average_factor
                    * self.legacy_manifest().normalization.identical_factor),
        );
        let values = (0..batch.len())
            .map(|row| evaluated[row * output_len].re.clone() * factor.clone())
            .collect::<Vec<_>>();
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        Ok((
            values,
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s: 0.0,
                stage_evaluator_s: 0.0,
                output_assign_s: 0.0,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
            },
        ))
    }

    fn vector_mass(&self, pdg: i32) -> PyResult<f64> {
        match pdg {
            22 => Ok(0.0),
            23 => Ok(self.legacy_manifest().model.mass_z),
            24 | -24 => {
                self.legacy_manifest().model.mass_w.ok_or_else(|| {
                    PyValueError::new_err("process manifest model is missing mass_w")
                })
            }
            _ => Err(PyValueError::new_err(format!(
                "unsupported electroweak vector PDG {pdg}"
            ))),
        }
    }

    fn vector_width(&self, pdg: i32) -> PyResult<f64> {
        match pdg {
            22 => Ok(0.0),
            23 => {
                self.legacy_manifest().model.width_z.ok_or_else(|| {
                    PyValueError::new_err("process manifest model is missing width_z")
                })
            }
            24 | -24 => {
                self.legacy_manifest().model.width_w.ok_or_else(|| {
                    PyValueError::new_err("process manifest model is missing width_w")
                })
            }
            _ => Err(PyValueError::new_err(format!(
                "unsupported electroweak vector PDG {pdg}"
            ))),
        }
    }

    fn fill_sources(&self, row: &mut [Complex<f64>], point: &[[f64; 4]; 16]) -> PyResult<()> {
        for source in &self.legacy_manifest().table.sources {
            let offset = self.legacy_manifest().layout.current_offsets[source.current_id];
            for component in 0..6 {
                row[offset + component] = c64(0.0, 0.0);
            }
            let wave = if source.source_kind == "lepton_pair_vector" {
                self.lepton_pair_vector_source(&source, point)?
            } else if source.source_kind == "external" {
                match source.leg_label {
                    1 | 2 => {
                        let current_pdg =
                            self.legacy_manifest().table.currents[source.current_id].pdg;
                        let p = negate(point[source.leg_label - 1]);
                        if current_pdg < 0 {
                            if source.physical_helicity == 1 {
                                ext_antiquark_weyl(p, 1, -1)
                            } else {
                                ext_antiquark_weyl(p, -1, 1)
                            }
                        } else if current_pdg > 0 {
                            if source.chirality == 1 {
                                ext_quark_weyl(p, -1, 1)
                            } else {
                                ext_quark_weyl(p, 1, -1)
                            }
                        } else {
                            return Err(PyValueError::new_err(format!(
                                "unexpected fermion source current PDG {current_pdg}"
                            )));
                        }
                    }
                    label if label >= 3 && label < self.legacy_manifest().gluon_count + 3 => {
                        ext_gluon(point[label - 1], source.helicity).to_vec()
                    }
                    label if label == self.legacy_manifest().gluon_count + 3 => {
                        let vector_pdg = self.legacy_manifest().external_pdg_order[label - 1];
                        if vector_pdg == 22 {
                            ext_gluon(point[label - 1], source.helicity).to_vec()
                        } else {
                            ext_massive_vector(
                                point[label - 1],
                                source.helicity,
                                self.vector_mass(vector_pdg)?,
                            )
                            .to_vec()
                        }
                    }
                    _ => {
                        return Err(PyValueError::new_err(format!(
                            "unexpected source leg label {}",
                            source.leg_label
                        )));
                    }
                }
            } else {
                return Err(PyValueError::new_err(format!(
                    "unsupported source kind {:?}",
                    source.source_kind
                )));
            };
            for (component, value) in wave.into_iter().enumerate() {
                row[offset + component] = value;
            }
        }
        Ok(())
    }

    fn lepton_pair_vector_source(
        &self,
        source: &SourceManifest,
        point: &[[f64; 4]; 16],
    ) -> PyResult<Vec<Complex<f64>>> {
        let partner_label = source
            .partner_leg_label
            .ok_or_else(|| PyValueError::new_err("lepton-pair source is missing partner leg"))?;
        let partner_helicity = source.partner_helicity.ok_or_else(|| {
            PyValueError::new_err("lepton-pair source is missing partner helicity")
        })?;
        let partner_chirality = source.partner_chirality.ok_or_else(|| {
            PyValueError::new_err("lepton-pair source is missing partner chirality")
        })?;
        let vector_pdg = source
            .vector_pdg
            .ok_or_else(|| PyValueError::new_err("lepton-pair source is missing vector PDG"))?;
        let coupling = source
            .coupling
            .ok_or_else(|| PyValueError::new_err("lepton-pair source is missing coupling"))?;
        let first_index = source
            .leg_label
            .checked_sub(1)
            .ok_or_else(|| PyValueError::new_err("source leg labels are one-based"))?;
        let second_index = partner_label
            .checked_sub(1)
            .ok_or_else(|| PyValueError::new_err("source leg labels are one-based"))?;
        if first_index >= self.legacy_manifest().external_pdg_order.len()
            || second_index >= self.legacy_manifest().external_pdg_order.len()
        {
            return Err(PyValueError::new_err(
                "lepton-pair source refers to a leg outside the process",
            ));
        }
        let first_pdg = self.legacy_manifest().external_pdg_order[first_index];
        let second_pdg = self.legacy_manifest().external_pdg_order[second_index];
        let current = if vector_pdg == 22 || vector_pdg == 23 {
            neutral_lepton_pair_current(
                first_pdg,
                point[first_index],
                source.helicity,
                source.chirality,
                second_pdg,
                point[second_index],
                partner_helicity,
                partner_chirality,
                coupling,
            )?
        } else if vector_pdg == 24 || vector_pdg == -24 {
            charged_lepton_pair_current(
                first_pdg,
                point[first_index],
                source.helicity,
                source.chirality,
                second_pdg,
                point[second_index],
                partner_helicity,
                partner_chirality,
                vector_pdg,
                coupling,
            )?
        } else {
            return Err(PyValueError::new_err(format!(
                "unsupported lepton-pair vector source PDG {vector_pdg}"
            )));
        };
        neutral_vector_propagator(
            current,
            add_momenta(point[first_index], point[second_index]),
            vector_pdg,
            self.vector_mass(vector_pdg)?,
            self.vector_width(vector_pdg)?,
        )
    }

    fn fill_sources_generic<T>(&self, row: &mut [Complex<T>], point: &[[T; 4]]) -> PyResult<()>
    where
        T: Real + RealLike + From<f64> + PartialOrd + Clone,
    {
        for source in &self.legacy_manifest().table.sources {
            let offset = self.legacy_manifest().layout.current_offsets[source.current_id];
            for component in 0..6 {
                row[offset + component] = complex_zero::<T>();
            }
            if source.source_kind != "external" {
                return Err(PyValueError::new_err(format!(
                    "high-precision source filling for {:?} is not implemented in rusticol yet",
                    source.source_kind
                )));
            }
            let wave = match source.leg_label {
                1 | 2 => {
                    let current_pdg = self.legacy_manifest().table.currents[source.current_id].pdg;
                    let p = negate_generic(&point[source.leg_label - 1]);
                    if current_pdg < 0 {
                        if source.physical_helicity == 1 {
                            ext_antiquark_weyl_generic(&p, 1, -1)
                        } else {
                            ext_antiquark_weyl_generic(&p, -1, 1)
                        }
                    } else if current_pdg > 0 {
                        if source.chirality == 1 {
                            ext_quark_weyl_generic(&p, -1, 1)
                        } else {
                            ext_quark_weyl_generic(&p, 1, -1)
                        }
                    } else {
                        return Err(PyValueError::new_err(format!(
                            "unexpected fermion source current PDG {current_pdg}"
                        )));
                    }
                }
                label if label >= 3 && label < self.legacy_manifest().gluon_count + 3 => {
                    ext_gluon_generic(&point[label - 1], source.helicity).to_vec()
                }
                label if label == self.legacy_manifest().gluon_count + 3 => {
                    let vector_pdg = self.legacy_manifest().external_pdg_order[label - 1];
                    if vector_pdg == 22 {
                        ext_gluon_generic(&point[label - 1], source.helicity).to_vec()
                    } else {
                        ext_massive_vector_generic(
                            &point[label - 1],
                            source.helicity,
                            T::from(self.vector_mass(vector_pdg)?),
                        )
                        .to_vec()
                    }
                }
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "unexpected source leg label {}",
                        source.leg_label
                    )));
                }
            };
            for (component, value) in wave.into_iter().enumerate() {
                row[offset + component] = value;
            }
        }
        Ok(())
    }

    fn fill_momenta(&self, row: &mut [Complex<f64>], point: &[[f64; 4]; 16]) {
        for entry in &self.legacy_manifest().layout.momentum_offsets_and_labels {
            let mut momentum = [0.0; 4];
            for label in &entry.labels {
                let sign = if *label <= 2 { -1.0 } else { 1.0 };
                for component in 0..4 {
                    momentum[component] += sign * point[*label - 1][component];
                }
            }
            for component in 0..4 {
                row[entry.offset + component] = c64(momentum[component], 0.0);
            }
        }
    }

    fn fill_momenta_generic<T>(&self, row: &mut [Complex<T>], point: &[[T; 4]])
    where
        T: Real + RealLike + From<f64> + Clone,
    {
        for entry in &self.legacy_manifest().layout.momentum_offsets_and_labels {
            let mut momentum = [T::new_zero(), T::new_zero(), T::new_zero(), T::new_zero()];
            for label in &entry.labels {
                let sign = if *label <= 2 {
                    T::from(-1.0)
                } else {
                    T::from(1.0)
                };
                for component in 0..4 {
                    momentum[component] += point[*label - 1][component].clone() * sign.clone();
                }
            }
            for component in 0..4 {
                row[entry.offset + component] =
                    Complex::new(momentum[component].clone(), T::new_zero());
            }
        }
    }
}

impl Drop for Runtime {
    fn drop(&mut self) {
        if let Some(stages) = self.stages.take() {
            std::mem::forget(stages);
        }
        if let Some(amplitude_stage) = self.amplitude_stage.take() {
            std::mem::forget(amplitude_stage);
        }
        if let Some(zero_gluon_stage) = self.zero_gluon_stage.take() {
            std::mem::forget(zero_gluon_stage);
        }
    }
}

impl GenericStageRuntimeV2 {
    fn load(stage: &GenericSerializedStageEvaluatorManifestV2, root: &Path) -> PyResult<Self> {
        let evaluator = EvaluatorGroup::load(&stage.evaluator, root)?;
        let mut outputs = Vec::new();
        for slot in &stage.output_slots {
            let output_len = slot
                .output_stop
                .checked_sub(slot.output_start)
                .ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "generic stage {} has an invalid output slot range",
                        stage.evaluator_label
                    ))
                })?;
            let component_len = slot
                .component_stop
                .checked_sub(slot.component_start)
                .ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "generic stage {} has an invalid component slot range",
                        stage.evaluator_label
                    ))
                })?;
            if output_len != component_len {
                return Err(PyValueError::new_err(format!(
                    "generic stage {} output slot length does not match component length",
                    stage.evaluator_label
                )));
            }
            for component in 0..component_len {
                outputs.push((
                    slot.output_start + component,
                    slot.component_start + component,
                ));
            }
        }
        let input_components = if stage.parameter_layout == "stage-local-value-momentum" {
            let mut map = vec![0usize; stage.parameter_count];
            for component in &stage.input_components {
                map[component.parameter_index] = component.global_component;
            }
            Some(map)
        } else {
            None
        };
        let output_spans = contiguous_output_spans(&outputs);
        Ok(Self {
            outputs,
            output_spans,
            input_components,
            parameter_scratch_f64: Vec::new(),
            output_scratch_f64: Vec::new(),
            evaluator,
        })
    }

    fn evaluate_f64_into_state(
        &mut self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<f64>],
    ) -> PyResult<(f64, f64)> {
        let eval_start = Instant::now();
        if let Some(input_components) = self.input_components.as_ref() {
            let local_parameter_count = input_components.len();
            self.parameter_scratch_f64
                .resize(batch_size * local_parameter_count, c64(0.0, 0.0));
            for row in 0..batch_size {
                let row_state = row * parameter_count;
                let row_params = row * local_parameter_count;
                for (local_index, global_index) in input_components.iter().enumerate() {
                    self.parameter_scratch_f64[row_params + local_index] =
                        state[row_state + *global_index];
                }
            }
            self.evaluator.evaluate_batch_into(
                batch_size,
                &self.parameter_scratch_f64,
                &mut self.output_scratch_f64,
            )?;
        } else {
            self.evaluator
                .evaluate_batch_into(batch_size, state, &mut self.output_scratch_f64)?;
        }
        let evaluator_s = eval_start.elapsed().as_secs_f64();

        let assign_start = Instant::now();
        for row in 0..batch_size {
            let row_state = row * parameter_count;
            let row_eval = row * self.evaluator.output_len;
            if self.output_spans.is_empty() {
                for (column, state_offset) in &self.outputs {
                    state[row_state + *state_offset] = self.output_scratch_f64[row_eval + *column];
                }
            } else {
                for (column_start, state_offset_start, len) in &self.output_spans {
                    let source_start = row_eval + *column_start;
                    let target_start = row_state + *state_offset_start;
                    state[target_start..target_start + *len].copy_from_slice(
                        &self.output_scratch_f64[source_start..source_start + *len],
                    );
                }
            }
        }
        Ok((evaluator_s, assign_start.elapsed().as_secs_f64()))
    }

    fn evaluate_generic_into_state<T>(
        &mut self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<(f64, f64)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let eval_start = Instant::now();
        let evaluated = if let Some(input_components) = self.input_components.as_ref() {
            let local_parameter_count = input_components.len();
            let mut parameter_scratch =
                vec![complex_zero::<T>(); batch_size * local_parameter_count];
            for row in 0..batch_size {
                let row_state = row * parameter_count;
                let row_params = row * local_parameter_count;
                for (local_index, global_index) in input_components.iter().enumerate() {
                    parameter_scratch[row_params + local_index] =
                        state[row_state + *global_index].clone();
                }
            }
            self.evaluator
                .evaluate_batch_generic(batch_size, &parameter_scratch, binary_precision)?
        } else {
            self.evaluator
                .evaluate_batch_generic(batch_size, state, binary_precision)?
        };
        let evaluator_s = eval_start.elapsed().as_secs_f64();

        let assign_start = Instant::now();
        for row in 0..batch_size {
            let row_state = row * parameter_count;
            let row_eval = row * self.evaluator.output_len;
            if self.output_spans.is_empty() {
                for (column, state_offset) in &self.outputs {
                    state[row_state + *state_offset] = evaluated[row_eval + *column].clone();
                }
            } else {
                for (column_start, state_offset_start, len) in &self.output_spans {
                    let source_start = row_eval + *column_start;
                    let target_start = row_state + *state_offset_start;
                    state[target_start..target_start + *len]
                        .clone_from_slice(&evaluated[source_start..source_start + *len]);
                }
            }
        }
        Ok((evaluator_s, assign_start.elapsed().as_secs_f64()))
    }
}

fn contiguous_output_spans(outputs: &[(usize, usize)]) -> Vec<(usize, usize, usize)> {
    if outputs.is_empty() {
        return Vec::new();
    }
    let mut spans = Vec::new();
    let (mut output_start, mut state_start) = outputs[0];
    let mut previous_output = output_start;
    let mut previous_state = state_start;
    let mut len = 1usize;
    for (output, state) in outputs.iter().copied().skip(1) {
        if output == previous_output + 1 && state == previous_state + 1 {
            previous_output = output;
            previous_state = state;
            len += 1;
            continue;
        }
        spans.push((output_start, state_start, len));
        output_start = output;
        state_start = state;
        previous_output = output;
        previous_state = state;
        len = 1;
    }
    spans.push((output_start, state_start, len));
    if spans.len() >= outputs.len() {
        Vec::new()
    } else {
        spans
    }
}

impl GenericAmplitudeRuntimeV2 {
    fn load(
        amplitude_stage: &GenericAmplitudeStageManifestV2,
        stage: &GenericSerializedStageEvaluatorManifestV2,
        root: &Path,
    ) -> PyResult<Self> {
        if stage.stage_kind != "amplitude-roots" {
            return Err(PyValueError::new_err(
                "generic amplitude runtime expected an amplitude-roots stage",
            ));
        }
        let raw_sum_weights = amplitude_stage
            .roots
            .iter()
            .map(|root| root.helicity_weight)
            .collect::<Vec<_>>();
        let raw_sum_group_ids = amplitude_stage
            .roots
            .iter()
            .map(generic_root_group_id)
            .collect::<PyResult<Vec<_>>>()?;
        let has_coherent_groups = raw_sum_group_ids.iter().any(Option::is_some);
        let raw_sum_groups = if has_coherent_groups {
            build_raw_sum_groups(
                amplitude_stage.output_count,
                &raw_sum_weights,
                &raw_sum_group_ids,
            )?
        } else {
            Vec::new()
        };
        let color_contraction = build_color_contraction_runtime(
            amplitude_stage.color_contraction.as_ref(),
            &raw_sum_groups,
        )?;
        Ok(Self {
            output_length: amplitude_stage.output_count,
            raw_sum_weights,
            raw_sum_groups,
            has_coherent_groups,
            color_contraction,
            input_components: if stage.parameter_layout == "stage-local-value-momentum" {
                let mut map = vec![0usize; stage.parameter_count];
                for component in &stage.input_components {
                    map[component.parameter_index] = component.global_component;
                }
                Some(map)
            } else {
                None
            },
            parameter_scratch_f64: Vec::new(),
            output_scratch_f64: Vec::new(),
            evaluator: EvaluatorGroup::load(&stage.evaluator, root)?,
        })
    }

    fn evaluate_f64_into_scratch(
        &mut self,
        batch_size: usize,
        state: &[Complex<f64>],
    ) -> PyResult<()> {
        if let Some(input_components) = self.input_components.as_ref() {
            let local_parameter_count = input_components.len();
            let global_parameter_count = state
                .len()
                .checked_div(batch_size)
                .ok_or_else(|| PyValueError::new_err("generic amplitude batch size is zero"))?;
            self.parameter_scratch_f64
                .resize(batch_size * local_parameter_count, c64(0.0, 0.0));
            for row in 0..batch_size {
                let row_state = row * global_parameter_count;
                let row_params = row * local_parameter_count;
                for (local_index, global_index) in input_components.iter().enumerate() {
                    self.parameter_scratch_f64[row_params + local_index] =
                        state[row_state + *global_index];
                }
            }
            return self.evaluator.evaluate_batch_into(
                batch_size,
                &self.parameter_scratch_f64,
                &mut self.output_scratch_f64,
            );
        }
        self.evaluator
            .evaluate_batch_into(batch_size, state, &mut self.output_scratch_f64)
    }

    fn reduce_scratch_f64_into(
        &mut self,
        batch_size: usize,
        raw_sums: &mut Vec<f64>,
    ) -> PyResult<()> {
        let amplitudes = &self.output_scratch_f64;
        if amplitudes.len() != batch_size * self.output_length {
            return Err(PyValueError::new_err(format!(
                "generic amplitude output buffer has length {}, expected {}",
                amplitudes.len(),
                batch_size * self.output_length
            )));
        }
        raw_sums.clear();
        raw_sums.resize(batch_size, 0.0);
        if let Some(contraction) = self.color_contraction.as_mut() {
            if self.raw_sum_groups.len() != contraction.group_count {
                return Err(PyValueError::new_err(
                    "colour contraction group count does not match coherent groups",
                ));
            }
            contraction
                .group_scratch_f64
                .resize(batch_size * contraction.group_count, c64(0.0, 0.0));
            for row in 0..batch_size {
                let row_offset = row * self.output_length;
                let group_row = row * contraction.group_count;
                for (group_index, group) in self.raw_sum_groups.iter().enumerate() {
                    let mut sum = c64(0.0, 0.0);
                    for index in &group.indices {
                        sum += amplitudes[row_offset + *index];
                    }
                    contraction.group_scratch_f64[group_row + group_index] = sum;
                }
                for entry in &contraction.entries {
                    let left = contraction.group_scratch_f64[group_row + entry.left_group_index];
                    let right = contraction.group_scratch_f64[group_row + entry.right_group_index];
                    let product = left * right.conj();
                    raw_sums[row] += entry.symmetry_factor
                        * (entry.weight_re * product.re - entry.weight_im * product.im);
                }
            }
            return Ok(());
        }
        for row in 0..batch_size {
            let row_offset = row * self.output_length;
            if self.has_coherent_groups {
                for group in &self.raw_sum_groups {
                    let mut sum = c64(0.0, 0.0);
                    for index in &group.indices {
                        sum += amplitudes[row_offset + *index];
                    }
                    raw_sums[row] += group.weight * (sum.re * sum.re + sum.im * sum.im);
                }
                continue;
            }
            for index in 0..self.output_length {
                let value = amplitudes[row_offset + index];
                raw_sums[row] +=
                    self.raw_sum_weights[index] * (value.re * value.re + value.im * value.im);
            }
        }
        Ok(())
    }

    fn evaluate_raw_sums_generic<T>(
        &mut self,
        batch_size: usize,
        state: &[Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<Vec<T>>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let evaluated = if let Some(input_components) = self.input_components.as_ref() {
            let local_parameter_count = input_components.len();
            let global_parameter_count = state
                .len()
                .checked_div(batch_size)
                .ok_or_else(|| PyValueError::new_err("generic amplitude batch size is zero"))?;
            let mut parameter_scratch =
                vec![complex_zero::<T>(); batch_size * local_parameter_count];
            for row in 0..batch_size {
                let row_state = row * global_parameter_count;
                let row_params = row * local_parameter_count;
                for (local_index, global_index) in input_components.iter().enumerate() {
                    parameter_scratch[row_params + local_index] =
                        state[row_state + *global_index].clone();
                }
            }
            self.evaluator
                .evaluate_batch_generic(batch_size, &parameter_scratch, binary_precision)?
        } else {
            self.evaluator
                .evaluate_batch_generic(batch_size, state, binary_precision)?
        };
        if evaluated.len() != batch_size * self.output_length {
            return Err(PyValueError::new_err(format!(
                "generic amplitude output buffer has length {}, expected {}",
                evaluated.len(),
                batch_size * self.output_length
            )));
        }
        let mut raw_sums = vec![T::new_zero(); batch_size];
        for row in 0..batch_size {
            let row_offset = row * self.output_length;
            if self.has_coherent_groups {
                for group in &self.raw_sum_groups {
                    let mut sum_re = T::new_zero();
                    let mut sum_im = T::new_zero();
                    for index in &group.indices {
                        let value = &evaluated[row_offset + *index];
                        sum_re += value.re.clone();
                        sum_im += value.im.clone();
                    }
                    raw_sums[row] +=
                        T::from(group.weight) * (sum_re.clone() * sum_re + sum_im.clone() * sum_im);
                }
                continue;
            }
            for index in 0..self.output_length {
                let value = &evaluated[row_offset + index];
                raw_sums[row] += T::from(self.raw_sum_weights[index])
                    * (value.re.clone() * value.re.clone()
                        + value.im.clone() * value.im.clone());
            }
        }
        Ok(raw_sums)
    }
}

impl CurrentStage {
    fn load(
        manifest: &CurrentStageManifest,
        root: &Path,
        current_offsets: &[usize],
    ) -> PyResult<Self> {
        let mut outputs = Vec::new();
        for slot in &manifest.output_slots {
            for (component, column) in (slot.start..slot.stop).enumerate() {
                let offset = *current_offsets.get(slot.id).ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "stage output references unknown current id {}",
                        slot.id
                    ))
                })? + component;
                outputs.push((column, slot.id, component, offset));
            }
        }
        let evaluator = EvaluatorGroup::load(&manifest.evaluator, root)?;
        let mut chunk_outputs = Vec::with_capacity(evaluator.evaluators.len());
        let mut chunk_start = 0;
        for chunk in &evaluator.evaluators {
            let chunk_stop = chunk_start + chunk.output_len;
            let outputs_for_chunk = outputs
                .iter()
                .filter_map(|(column, _, _, state_offset)| {
                    if *column >= chunk_start && *column < chunk_stop {
                        Some((*column - chunk_start, *state_offset))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            chunk_outputs.push(outputs_for_chunk);
            chunk_start = chunk_stop;
        }
        Ok(Self {
            outputs,
            chunk_outputs,
            evaluator,
        })
    }

    fn evaluate_f64_into_state(
        &mut self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<f64>],
    ) -> PyResult<(f64, f64)> {
        let mut evaluator_time_s = 0.0;
        let mut assign_time_s = 0.0;
        for chunk_index in 0..self.evaluator.evaluators.len() {
            let evaluator = &mut self.evaluator.evaluators[chunk_index];
            if state.len() != batch_size * evaluator.input_len {
                return Err(PyValueError::new_err(format!(
                    "parameter buffer has length {}, expected {}",
                    state.len(),
                    batch_size * evaluator.input_len
                )));
            }
            let chunk_output_len = evaluator.output_len;
            self.evaluator
                .chunk_scratch_f64
                .resize(batch_size * chunk_output_len, c64(0.0, 0.0));

            let eval_start = Instant::now();
            evaluator.evaluate_f64_batch(
                batch_size,
                &*state,
                &mut self.evaluator.chunk_scratch_f64,
            )?;
            evaluator_time_s += eval_start.elapsed().as_secs_f64();

            let assign_start = Instant::now();
            for row in 0..batch_size {
                let state_row = row * parameter_count;
                let chunk_row = row * chunk_output_len;
                for (local_column, state_offset) in &self.chunk_outputs[chunk_index] {
                    state[state_row + *state_offset] =
                        self.evaluator.chunk_scratch_f64[chunk_row + *local_column];
                }
            }
            assign_time_s += assign_start.elapsed().as_secs_f64();
        }
        Ok((evaluator_time_s, assign_time_s))
    }
}

impl AmplitudeStage {
    fn load(manifest: &AmplitudeStageManifest, root: &Path) -> PyResult<Self> {
        let raw_sum_group_ids = manifest
            .raw_sum_group_ids
            .clone()
            .unwrap_or_else(|| vec![None; manifest.output_length]);
        let has_coherent_groups = raw_sum_group_ids.iter().any(Option::is_some);
        let raw_sum_groups = if has_coherent_groups {
            build_raw_sum_groups(
                manifest.output_length,
                &manifest.raw_sum_weights,
                &raw_sum_group_ids,
            )?
        } else {
            Vec::new()
        };
        Ok(Self {
            output_length: manifest.output_length,
            raw_sum_weights: manifest.raw_sum_weights.clone(),
            raw_sum_groups,
            has_coherent_groups,
            amplitude_evaluator: manifest
                .amplitude_evaluator
                .as_ref()
                .map(|m| EvaluatorGroup::load(m, root))
                .transpose()?,
            raw_sum_evaluator: manifest
                .raw_sum_evaluator
                .as_ref()
                .map(|m| EvaluatorGroup::load(m, root))
                .transpose()?,
        })
    }

    fn evaluate_raw_sums(
        &mut self,
        batch_size: usize,
        state: &[Complex<f64>],
    ) -> PyResult<Vec<f64>> {
        if let Some(raw_sum_evaluator) = &mut self.raw_sum_evaluator {
            let evaluated = raw_sum_evaluator.evaluate_batch(batch_size, state)?;
            return Ok((0..batch_size)
                .map(|row| evaluated[row * raw_sum_evaluator.output_len].re)
                .collect());
        }
        let amplitude_evaluator = self
            .amplitude_evaluator
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("amplitude stage has no amplitude evaluator"))?;
        if !self.has_coherent_groups {
            let mut raw_sums = vec![0.0; batch_size];
            let mut output_offset = 0;
            for evaluator in &mut amplitude_evaluator.evaluators {
                if state.len() != batch_size * evaluator.input_len {
                    return Err(PyValueError::new_err(format!(
                        "parameter buffer has length {}, expected {}",
                        state.len(),
                        batch_size * evaluator.input_len
                    )));
                }
                amplitude_evaluator
                    .chunk_scratch_f64
                    .resize(batch_size * evaluator.output_len, c64(0.0, 0.0));
                evaluator.evaluate_f64_batch(
                    batch_size,
                    state,
                    &mut amplitude_evaluator.chunk_scratch_f64,
                )?;
                for row in 0..batch_size {
                    let row_offset = row * evaluator.output_len;
                    for local_index in 0..evaluator.output_len {
                        let value = amplitude_evaluator.chunk_scratch_f64[row_offset + local_index];
                        let weight = self.raw_sum_weights[output_offset + local_index];
                        raw_sums[row] += weight * (value.re * value.re + value.im * value.im);
                    }
                }
                output_offset += evaluator.output_len;
            }
            return Ok(raw_sums);
        }
        let evaluated = amplitude_evaluator.evaluate_batch(batch_size, state)?;
        let mut raw_sums = vec![0.0; batch_size];
        for row in 0..batch_size {
            let row_offset = row * amplitude_evaluator.output_len;
            if self.has_coherent_groups {
                for group in &self.raw_sum_groups {
                    let mut sum = c64(0.0, 0.0);
                    for index in &group.indices {
                        sum += evaluated[row_offset + *index];
                    }
                    raw_sums[row] += group.weight * (sum.re * sum.re + sum.im * sum.im);
                }
                continue;
            }
            for index in 0..self.output_length {
                let value = evaluated[row_offset + index];
                let weight = self.raw_sum_weights[index];
                raw_sums[row] += weight * (value.re * value.re + value.im * value.im);
            }
        }
        Ok(raw_sums)
    }

    fn evaluate_raw_sums_generic<T>(
        &mut self,
        batch_size: usize,
        state: &[Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<Vec<T>>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        if let Some(raw_sum_evaluator) = &mut self.raw_sum_evaluator {
            let evaluated =
                raw_sum_evaluator.evaluate_batch_generic(batch_size, state, binary_precision)?;
            return Ok((0..batch_size)
                .map(|row| evaluated[row * raw_sum_evaluator.output_len].re.clone())
                .collect());
        }
        let amplitude_evaluator = self
            .amplitude_evaluator
            .as_mut()
            .ok_or_else(|| PyRuntimeError::new_err("amplitude stage has no amplitude evaluator"))?;
        let evaluated =
            amplitude_evaluator.evaluate_batch_generic(batch_size, state, binary_precision)?;
        let mut raw_sums = vec![T::new_zero(); batch_size];
        for row in 0..batch_size {
            let row_offset = row * amplitude_evaluator.output_len;
            if self.has_coherent_groups {
                for group in &self.raw_sum_groups {
                    let mut sum_re = T::new_zero();
                    let mut sum_im = T::new_zero();
                    for index in &group.indices {
                        let value = &evaluated[row_offset + *index];
                        sum_re += value.re.clone();
                        sum_im += value.im.clone();
                    }
                    raw_sums[row] +=
                        T::from(group.weight) * (sum_re.clone() * sum_re + sum_im.clone() * sum_im);
                }
                continue;
            }
            for index in 0..self.output_length {
                let value = &evaluated[row_offset + index];
                let weight = T::from(self.raw_sum_weights[index]);
                raw_sums[row] += weight
                    * (value.re.clone() * value.re.clone() + value.im.clone() * value.im.clone());
            }
        }
        Ok(raw_sums)
    }
}

fn build_raw_sum_groups(
    output_length: usize,
    weights: &[f64],
    group_ids: &[Option<i64>],
) -> PyResult<Vec<RawSumGroup>> {
    if weights.len() != output_length || group_ids.len() != output_length {
        return Err(PyValueError::new_err(
            "raw-sum group metadata length does not match amplitude outputs",
        ));
    }
    let mut grouped: BTreeMap<i64, Vec<usize>> = BTreeMap::new();
    let mut groups = Vec::new();
    for index in 0..output_length {
        if let Some(group_id) = group_ids[index] {
            grouped.entry(group_id).or_default().push(index);
        } else {
            groups.push(RawSumGroup {
                id: index as i64,
                indices: vec![index],
                weight: weights[index],
            });
        }
    }
    for (group_id, indices) in grouped {
        let weight = weights[indices[0]];
        if indices
            .iter()
            .any(|index| (weights[*index] - weight).abs() > 0.0)
        {
            return Err(PyValueError::new_err(format!(
                "coherent amplitude group {group_id} has inconsistent raw-sum weights"
            )));
        }
        groups.push(RawSumGroup {
            id: group_id,
            indices,
            weight,
        });
    }
    Ok(groups)
}

fn build_color_contraction_runtime(
    manifest: Option<&GenericColorContractionManifestV2>,
    groups: &[RawSumGroup],
) -> PyResult<Option<ColorContractionRuntime>> {
    let Some(manifest) = manifest else {
        return Ok(None);
    };
    if !manifest.supported {
        return Err(PyValueError::new_err(format!(
            "generic colour contraction is unsupported: {}",
            manifest
                .reason
                .as_deref()
                .unwrap_or("no diagnostic was provided")
        )));
    }
    let group_index_by_id = groups
        .iter()
        .enumerate()
        .map(|(index, group)| (group.id, index))
        .collect::<BTreeMap<_, _>>();
    if group_index_by_id.len() != manifest.group_count {
        return Err(PyValueError::new_err(format!(
            "colour contraction declares {} groups but runtime has {} coherent groups",
            manifest.group_count,
            group_index_by_id.len()
        )));
    }
    let mut entries = Vec::with_capacity(manifest.entries.len());
    for entry in &manifest.entries {
        let left_group_index = *group_index_by_id.get(&entry.left_group_id).ok_or_else(|| {
            PyValueError::new_err(format!(
                "colour contraction references unknown left group {}",
                entry.left_group_id
            ))
        })?;
        let right_group_index = *group_index_by_id
            .get(&entry.right_group_id)
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "colour contraction references unknown right group {}",
                    entry.right_group_id
                ))
            })?;
        let weight_re = entry.weight.first().copied().unwrap_or(0.0);
        let weight_im = entry.weight.get(1).copied().unwrap_or(0.0);
        entries.push(ColorContractionEntry {
            left_group_index,
            right_group_index,
            weight_re,
            weight_im,
            symmetry_factor: entry.symmetry_factor,
        });
    }
    Ok(Some(ColorContractionRuntime {
        group_count: manifest.group_count,
        entries,
        group_scratch_f64: Vec::new(),
    }))
}

fn generic_root_group_id(root: &GenericAmplitudeRootManifestV2) -> PyResult<Option<i64>> {
    let Some(value) = root.coherent_group_id.as_ref() else {
        return Ok(None);
    };
    if value.is_null() {
        return Ok(None);
    }
    if let Some(group_id) = value.as_i64() {
        return Ok(Some(group_id));
    }
    if let Some(group_id) = value.as_u64() {
        return i64::try_from(group_id)
            .map(Some)
            .map_err(|_| PyValueError::new_err("generic coherent group id exceeds i64"));
    }
    if let Some(text) = value.as_str() {
        return text.parse::<i64>().map(Some).map_err(|err| {
            PyValueError::new_err(format!(
                "could not parse generic coherent group id {text:?}: {err}"
            ))
        });
    }
    Err(PyValueError::new_err(format!(
        "generic coherent group id for root {} must be an integer or string",
        root.root_id
    )))
}

impl EvaluatorGroup {
    fn load(manifest: &EvaluatorManifest, root: &Path) -> PyResult<Self> {
        let mut evaluators = Vec::new();
        flatten_evaluators(manifest, root, &mut evaluators)?;
        let output_len = evaluators.iter().map(|e| e.output_len).sum();
        Ok(Self {
            evaluators,
            output_len,
            chunk_scratch_f64: Vec::new(),
        })
    }

    fn evaluate_batch(
        &mut self,
        batch_size: usize,
        params: &[Complex<f64>],
    ) -> PyResult<Vec<Complex<f64>>> {
        let mut out = Vec::new();
        self.evaluate_batch_into(batch_size, params, &mut out)?;
        Ok(out)
    }

    fn evaluate_batch_into(
        &mut self,
        batch_size: usize,
        params: &[Complex<f64>],
        out: &mut Vec<Complex<f64>>,
    ) -> PyResult<()> {
        let expected_output_len = batch_size * self.output_len;
        if out.len() != expected_output_len {
            out.resize(expected_output_len, c64(0.0, 0.0));
        }
        if self.evaluators.len() == 1 {
            let evaluator = &mut self.evaluators[0];
            if params.len() != batch_size * evaluator.input_len {
                return Err(PyValueError::new_err(format!(
                    "parameter buffer has length {}, expected {}",
                    params.len(),
                    batch_size * evaluator.input_len
                )));
            }
            evaluator.evaluate_f64_batch(batch_size, params, out)?;
            return Ok(());
        }
        let mut output_offset = 0;
        for evaluator in &mut self.evaluators {
            if params.len() != batch_size * evaluator.input_len {
                return Err(PyValueError::new_err(format!(
                    "parameter buffer has length {}, expected {}",
                    params.len(),
                    batch_size * evaluator.input_len
                )));
            }
            self.chunk_scratch_f64
                .resize(batch_size * evaluator.output_len, c64(0.0, 0.0));
            evaluator.evaluate_f64_batch(batch_size, params, &mut self.chunk_scratch_f64)?;
            for row in 0..batch_size {
                let src = row * evaluator.output_len;
                let dst = row * self.output_len + output_offset;
                out[dst..dst + evaluator.output_len]
                    .copy_from_slice(&self.chunk_scratch_f64[src..src + evaluator.output_len]);
            }
            output_offset += evaluator.output_len;
        }
        Ok(())
    }

    fn evaluate_batch_generic<T>(
        &mut self,
        batch_size: usize,
        params: &[Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<Vec<Complex<T>>>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let mut out = vec![complex_zero::<T>(); batch_size * self.output_len];
        let mut output_offset = 0;
        for evaluator in &mut self.evaluators {
            if params.len() != batch_size * evaluator.input_len {
                return Err(PyValueError::new_err(format!(
                    "parameter buffer has length {}, expected {}",
                    params.len(),
                    batch_size * evaluator.input_len
                )));
            }
            let mut chunk_out = vec![complex_zero::<T>(); batch_size * evaluator.output_len];
            for row in 0..batch_size {
                let in_start = row * evaluator.input_len;
                let out_start = row * evaluator.output_len;
                T::evaluate_loaded(
                    evaluator,
                    &params[in_start..in_start + evaluator.input_len],
                    &mut chunk_out[out_start..out_start + evaluator.output_len],
                    binary_precision,
                )?;
            }
            for row in 0..batch_size {
                let src = row * evaluator.output_len;
                let dst = row * self.output_len + output_offset;
                out[dst..dst + evaluator.output_len]
                    .clone_from_slice(&chunk_out[src..src + evaluator.output_len]);
            }
            output_offset += evaluator.output_len;
        }
        Ok(out)
    }
}

impl LoadedEvaluator {
    fn evaluate_f64_batch(
        &mut self,
        batch_size: usize,
        params: &[Complex<f64>],
        out: &mut [Complex<f64>],
    ) -> PyResult<()> {
        match &mut self.eval {
            F64Evaluator::Compiled(eval) => eval
                .evaluate_batch(batch_size, params, out)
                .map_err(PyRuntimeError::new_err),
            F64Evaluator::Jit(eval) => {
                if params.len() != batch_size * self.input_len {
                    return Err(PyValueError::new_err(format!(
                        "parameter buffer has length {}, expected {}",
                        params.len(),
                        batch_size * self.input_len
                    )));
                }
                if out.len() != batch_size * self.output_len {
                    return Err(PyValueError::new_err(format!(
                        "output buffer has length {}, expected {}",
                        out.len(),
                        batch_size * self.output_len
                    )));
                }
                // SymJIT's matrix-evaluation entry point is still brittle for a
                // few large complex AArch64 kernels. The scalar entry point is
                // stable on the same serialized evaluators and JIT rows are not
                // the production fast path, so keep C++ evaluators batched and
                // run JIT evaluators row-by-row.
                for row in 0..batch_size {
                    let in_start = row * self.input_len;
                    let out_start = row * self.output_len;
                    eval.evaluate(
                        &params[in_start..in_start + self.input_len],
                        &mut out[out_start..out_start + self.output_len],
                    );
                }
                Ok(())
            }
            F64Evaluator::Interpreted(eval) => {
                if params.len() != batch_size * self.input_len {
                    return Err(PyValueError::new_err(format!(
                        "parameter buffer has length {}, expected {}",
                        params.len(),
                        batch_size * self.input_len
                    )));
                }
                if out.len() != batch_size * self.output_len {
                    return Err(PyValueError::new_err(format!(
                        "output buffer has length {}, expected {}",
                        out.len(),
                        batch_size * self.output_len
                    )));
                }
                for row in 0..batch_size {
                    let in_start = row * self.input_len;
                    let out_start = row * self.output_len;
                    eval.evaluate(
                        &params[in_start..in_start + self.input_len],
                        &mut out[out_start..out_start + self.output_len],
                    );
                }
                Ok(())
            }
        }
    }
}

impl ZeroGluonStage {
    fn load(manifest: &ZeroGluonManifest, root: &Path) -> PyResult<Self> {
        let state_path = artifact_path(root, &manifest.evaluator_state_path);
        let (exact_eval, _jit_complex) = load_evaluator_state(&state_path)?;
        let eval_complex = exact_eval
            .clone()
            .map_coeff(&|c| Complex::new(c.re.to_f64(), c.im.to_f64()));
        let input_len = manifest.parameter_names.len();
        Ok(Self {
            parameter_names: manifest.parameter_names.clone(),
            evaluator: LoadedEvaluator {
                eval: F64Evaluator::Interpreted(eval_complex),
                exact_eval: Some(exact_eval),
                double_eval: None,
                arb_eval: None,
                input_len,
                output_len: 1,
            },
            z_left: manifest.z_left,
            z_right: manifest.z_right,
        })
    }

    fn parameter_rows_f64(
        &self,
        manifest: &ProcessManifest,
        batch: &[[[f64; 4]; 16]],
    ) -> PyResult<Vec<Complex<f64>>> {
        let mut rows = Vec::with_capacity(batch.len() * self.parameter_names.len());
        for point in batch {
            let context = ZeroGluonParameterContext::new_f64(manifest, point, self)?;
            for name in &self.parameter_names {
                rows.push(context.value_f64(name)?);
            }
        }
        Ok(rows)
    }

    fn parameter_rows_generic<T>(
        &self,
        manifest: &ProcessManifest,
        batch: &[Vec<[T; 4]>],
    ) -> PyResult<Vec<Complex<T>>>
    where
        T: Real + RealLike + From<f64> + PartialOrd + Clone,
    {
        let mut rows = Vec::with_capacity(batch.len() * self.parameter_names.len());
        for point in batch {
            let context = ZeroGluonParameterContext::new_generic(manifest, point, self)?;
            for name in &self.parameter_names {
                rows.push(context.value_generic(name)?);
            }
        }
        Ok(rows)
    }
}

struct ZeroGluonParameterContext<T> {
    quark_minus: Vec<Complex<T>>,
    quark_plus: Vec<Complex<T>>,
    antiquark_minus: Vec<Complex<T>>,
    antiquark_plus: Vec<Complex<T>>,
    z_minus: Vec<Complex<T>>,
    z_zero: Vec<Complex<T>>,
    z_plus: Vec<Complex<T>>,
    z_left: T,
    z_right: T,
}

impl ZeroGluonParameterContext<f64> {
    fn new_f64(
        manifest: &ProcessManifest,
        point: &[[f64; 4]; 16],
        stage: &ZeroGluonStage,
    ) -> PyResult<Self> {
        let (quark_momentum, antiquark_momentum) =
            zero_gluon_quark_antiquark_momenta_f64(manifest, point)?;
        let z_momentum = point[2];
        Ok(Self {
            quark_minus: ext_quark_dirac(quark_momentum, -1).to_vec(),
            quark_plus: ext_quark_dirac(quark_momentum, 1).to_vec(),
            antiquark_minus: ext_antiquark_dirac(antiquark_momentum, -1).to_vec(),
            antiquark_plus: ext_antiquark_dirac(antiquark_momentum, 1).to_vec(),
            z_minus: ext_massive_vector(z_momentum, -1, manifest.model.mass_z).to_vec(),
            z_zero: ext_massive_vector(z_momentum, 0, manifest.model.mass_z).to_vec(),
            z_plus: ext_massive_vector(z_momentum, 1, manifest.model.mass_z).to_vec(),
            z_left: stage.z_left,
            z_right: stage.z_right,
        })
    }

    fn value_f64(&self, name: &str) -> PyResult<Complex<f64>> {
        if name == "z_left" {
            return Ok(c64(self.z_left, 0.0));
        }
        if name == "z_right" {
            return Ok(c64(self.z_right, 0.0));
        }
        let (prefix, helicity, component, part) = parse_zero_gluon_parameter_name(name)?;
        let wave = self.wave_f64(prefix, helicity)?;
        let value = wave.get(component).ok_or_else(|| {
            PyValueError::new_err(format!("zero-gluon component index out of range in {name}"))
        })?;
        match part {
            "re" => Ok(c64(value.re, 0.0)),
            "im" => Ok(c64(value.im, 0.0)),
            _ => Err(PyValueError::new_err(format!(
                "unsupported zero-gluon parameter component {part:?} in {name}"
            ))),
        }
    }

    fn wave_f64(&self, prefix: &str, helicity: i32) -> PyResult<&Vec<Complex<f64>>> {
        match (prefix, helicity) {
            ("q", -1) => Ok(&self.quark_minus),
            ("q", 1) => Ok(&self.quark_plus),
            ("aq", -1) => Ok(&self.antiquark_minus),
            ("aq", 1) => Ok(&self.antiquark_plus),
            ("z", -1) => Ok(&self.z_minus),
            ("z", 0) => Ok(&self.z_zero),
            ("z", 1) => Ok(&self.z_plus),
            _ => Err(PyValueError::new_err(format!(
                "unsupported zero-gluon wavefunction label {prefix}_{helicity}"
            ))),
        }
    }
}

impl<T> ZeroGluonParameterContext<T>
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    fn new_generic(
        manifest: &ProcessManifest,
        point: &[[T; 4]],
        stage: &ZeroGluonStage,
    ) -> PyResult<Self> {
        let (quark_momentum, antiquark_momentum) =
            zero_gluon_quark_antiquark_momenta_generic(manifest, point)?;
        let z_momentum = point
            .get(2)
            .ok_or_else(|| PyValueError::new_err("zero-gluon point is missing Z momentum"))?;
        Ok(Self {
            quark_minus: ext_quark_dirac_generic(&quark_momentum, -1).to_vec(),
            quark_plus: ext_quark_dirac_generic(&quark_momentum, 1).to_vec(),
            antiquark_minus: ext_antiquark_dirac_generic(&antiquark_momentum, -1).to_vec(),
            antiquark_plus: ext_antiquark_dirac_generic(&antiquark_momentum, 1).to_vec(),
            z_minus: ext_massive_vector_generic(z_momentum, -1, T::from(manifest.model.mass_z))
                .to_vec(),
            z_zero: ext_massive_vector_generic(z_momentum, 0, T::from(manifest.model.mass_z))
                .to_vec(),
            z_plus: ext_massive_vector_generic(z_momentum, 1, T::from(manifest.model.mass_z))
                .to_vec(),
            z_left: T::from(stage.z_left),
            z_right: T::from(stage.z_right),
        })
    }

    fn value_generic(&self, name: &str) -> PyResult<Complex<T>> {
        if name == "z_left" {
            return Ok(Complex::new(self.z_left.clone(), T::new_zero()));
        }
        if name == "z_right" {
            return Ok(Complex::new(self.z_right.clone(), T::new_zero()));
        }
        let (prefix, helicity, component, part) = parse_zero_gluon_parameter_name(name)?;
        let wave = self.wave_generic(prefix, helicity)?;
        let value = wave.get(component).ok_or_else(|| {
            PyValueError::new_err(format!("zero-gluon component index out of range in {name}"))
        })?;
        match part {
            "re" => Ok(Complex::new(value.re.clone(), T::new_zero())),
            "im" => Ok(Complex::new(value.im.clone(), T::new_zero())),
            _ => Err(PyValueError::new_err(format!(
                "unsupported zero-gluon parameter component {part:?} in {name}"
            ))),
        }
    }

    fn wave_generic(&self, prefix: &str, helicity: i32) -> PyResult<&Vec<Complex<T>>> {
        match (prefix, helicity) {
            ("q", -1) => Ok(&self.quark_minus),
            ("q", 1) => Ok(&self.quark_plus),
            ("aq", -1) => Ok(&self.antiquark_minus),
            ("aq", 1) => Ok(&self.antiquark_plus),
            ("z", -1) => Ok(&self.z_minus),
            ("z", 0) => Ok(&self.z_zero),
            ("z", 1) => Ok(&self.z_plus),
            _ => Err(PyValueError::new_err(format!(
                "unsupported zero-gluon wavefunction label {prefix}_{helicity}"
            ))),
        }
    }
}

fn parse_zero_gluon_parameter_name(name: &str) -> PyResult<(&str, i32, usize, &str)> {
    let parts = name.split('_').collect::<Vec<_>>();
    if parts.len() != 4 {
        return Err(PyValueError::new_err(format!(
            "invalid zero-gluon parameter name {name:?}"
        )));
    }
    let helicity = parse_helicity_token(parts[1])?;
    let component = parts[2].parse::<usize>().map_err(|err| {
        PyValueError::new_err(format!(
            "invalid zero-gluon component index in {name:?}: {err}"
        ))
    })?;
    Ok((parts[0], helicity, component, parts[3]))
}

fn parse_helicity_token(token: &str) -> PyResult<i32> {
    let (sign, value) = token.split_at(1);
    let magnitude = value
        .parse::<i32>()
        .map_err(|err| PyValueError::new_err(format!("invalid helicity token {token:?}: {err}")))?;
    match sign {
        "m" => Ok(-magnitude),
        "p" => Ok(magnitude),
        _ => Err(PyValueError::new_err(format!(
            "invalid helicity token sign in {token:?}"
        ))),
    }
}

fn zero_gluon_quark_antiquark_momenta_f64(
    manifest: &ProcessManifest,
    point: &[[f64; 4]; 16],
) -> PyResult<([f64; 4], [f64; 4])> {
    if manifest.external_pdg_order.len() != 3 {
        return Err(PyValueError::new_err(
            "zero-gluon process requires three external momenta",
        ));
    }
    let current0 = -manifest.external_pdg_order[0];
    let current1 = -manifest.external_pdg_order[1];
    let p0 = negate(point[0]);
    let p1 = negate(point[1]);
    if current0 > 0 && current1 < 0 {
        Ok((p0, p1))
    } else if current1 > 0 && current0 < 0 {
        Ok((p1, p0))
    } else {
        Err(PyValueError::new_err(
            "zero-gluon process requires one incoming quark and one antiquark",
        ))
    }
}

fn zero_gluon_quark_antiquark_momenta_generic<T>(
    manifest: &ProcessManifest,
    point: &[[T; 4]],
) -> PyResult<([T; 4], [T; 4])>
where
    T: Real + Clone,
{
    if manifest.external_pdg_order.len() != 3 || point.len() < 3 {
        return Err(PyValueError::new_err(
            "zero-gluon process requires three external momenta",
        ));
    }
    let current0 = -manifest.external_pdg_order[0];
    let current1 = -manifest.external_pdg_order[1];
    let p0 = negate_generic(&point[0]);
    let p1 = negate_generic(&point[1]);
    if current0 > 0 && current1 < 0 {
        Ok((p0, p1))
    } else if current1 > 0 && current0 < 0 {
        Ok((p1, p0))
    } else {
        Err(PyValueError::new_err(
            "zero-gluon process requires one incoming quark and one antiquark",
        ))
    }
}

fn flatten_evaluators(
    manifest: &EvaluatorManifest,
    root: &Path,
    output: &mut Vec<LoadedEvaluator>,
) -> PyResult<()> {
    match manifest {
        EvaluatorManifest::Jit {
            input_len,
            output_len,
            evaluator_state_path,
        } => {
            let (exact_eval, jit_eval) =
                load_evaluator_state(&artifact_path(root, evaluator_state_path))?;
            let eval = jit_eval.ok_or_else(|| {
                PyValueError::new_err(
                    "jit-symbolica-evaluator artifact has no saved complex JIT payload; \
                     evaluate the Symbolica evaluator once before saving it",
                )
            })?;
            output.push(LoadedEvaluator {
                eval: F64Evaluator::Jit(eval),
                exact_eval: Some(exact_eval),
                double_eval: None,
                arb_eval: None,
                input_len: *input_len,
                output_len: *output_len,
            });
            Ok(())
        }
        EvaluatorManifest::CompiledComplex {
            function_name,
            input_len,
            output_len,
            library_path,
            evaluator_state_path,
            number_type,
        } => {
            if number_type != "complex" {
                return Err(PyValueError::new_err(format!(
                    "rusticol currently supports compiled complex evaluators, got {number_type}"
                )));
            }
            let library = artifact_path(root, library_path);
            let eval = CompiledComplexEvaluator::load(&library, function_name).map_err(|err| {
                PyRuntimeError::new_err(format!(
                    "could not load compiled evaluator {} from {}: {err}",
                    function_name,
                    library.display()
                ))
            })?;
            let exact_eval = evaluator_state_path
                .as_deref()
                .map(|state_path| {
                    load_evaluator_state(&artifact_path(root, state_path)).map(|state| state.0)
                })
                .transpose()?;
            output.push(LoadedEvaluator {
                eval: F64Evaluator::Compiled(eval),
                exact_eval,
                double_eval: None,
                arb_eval: None,
                input_len: *input_len,
                output_len: *output_len,
            });
            Ok(())
        }
        EvaluatorManifest::Chunked { chunks } => {
            for chunk in chunks {
                flatten_evaluators(chunk, root, output)?;
            }
            Ok(())
        }
    }
}

fn load_evaluator_state(
    path: &Path,
) -> PyResult<(
    ExpressionEvaluator<Complex<Rational>>,
    Option<JITCompiledEvaluator<Complex<f64>>>,
)> {
    type SavedEvaluator = (
        bool,
        JITCompilationSettings,
        ExpressionEvaluator<Complex<Rational>>,
        Option<JITCompiledEvaluator<f64>>,
        Option<JITCompiledEvaluator<Complex<f64>>>,
    );
    type LegacySavedEvaluator = (
        bool,
        ExpressionEvaluator<Complex<Rational>>,
        Option<JITCompiledEvaluator<f64>>,
        Option<JITCompiledEvaluator<Complex<f64>>>,
    );

    let bytes = fs::read(path).map_err(|err| {
        PyValueError::new_err(format!(
            "could not read evaluator state {}: {err}",
            path.display()
        ))
    })?;
    match bincode::decode_from_slice::<SavedEvaluator, _>(&bytes, bincode::config::standard()) {
        Ok(((_, _, evaluator, _, jit_complex), _)) => Ok((evaluator, jit_complex)),
        Err(new_err) => {
            let decoded = bincode::decode_from_slice::<LegacySavedEvaluator, _>(
                &bytes,
                bincode::config::standard(),
            )
            .map_err(|_| {
                PyValueError::new_err(format!(
                    "could not decode evaluator state {}: {new_err}",
                    path.display()
                ))
            })?;
            let (_, evaluator, _, jit_complex) = decoded.0;
            Ok((evaluator, jit_complex))
        }
    }
}

fn artifact_path(root: &Path, value: &str) -> PathBuf {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        path
    } else {
        root.join(path)
    }
}

fn batch_momenta(
    py: Python<'_>,
    momenta: &Bound<'_, PyAny>,
    expected_legs: usize,
) -> PyResult<Vec<[[f64; 4]; 16]>> {
    match PyBuffer::<f64>::get(momenta) {
        Ok(buffer) => batch_momenta_from_buffer(py, &buffer, expected_legs),
        Err(buffer_error) => match momenta.extract::<Vec<Vec<Vec<f64>>>>() {
            Ok(points) => batch_momenta_from_nested(points, expected_legs),
            Err(sequence_error) => Err(PyValueError::new_err(format!(
                "momenta must be a C-contiguous f64 buffer or nested sequence with shape \
                 (batch, {expected_legs}, 4); buffer error: {buffer_error}; sequence error: \
                 {sequence_error}",
            ))),
        },
    }
}

fn input_crossing_map_to_json(map: &Vec<InputCrossingMapEntry>) -> Vec<(usize, usize, f64)> {
    map.iter()
        .map(|entry| (entry.target_index, entry.source_index, entry.sign))
        .collect()
}

fn apply_input_crossing_map(
    batch: &[Vec<[f64; 4]>],
    expected_legs: usize,
    input_crossing_map: Option<&[InputCrossingMapEntry]>,
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    let Some(map) = input_crossing_map else {
        return Ok(batch.to_vec());
    };
    if map.len() != expected_legs {
        return Err(PyValueError::new_err(format!(
            "input crossing map has {} entries, expected {expected_legs}",
            map.len()
        )));
    }
    let mut seen = vec![false; expected_legs];
    for entry in map {
        if entry.target_index >= expected_legs || entry.source_index >= expected_legs {
            return Err(PyValueError::new_err(
                "input crossing map references an out-of-range external leg",
            ));
        }
        if seen[entry.target_index] {
            return Err(PyValueError::new_err(
                "input crossing map contains a duplicate target index",
            ));
        }
        seen[entry.target_index] = true;
    }
    if seen.iter().any(|value| !*value) {
        return Err(PyValueError::new_err(
            "input crossing map does not cover every target index",
        ));
    }
    let mut mapped_batch = Vec::with_capacity(batch.len());
    for point in batch {
        let mut mapped = vec![[0.0; 4]; expected_legs];
        for entry in map {
            let source = point[entry.source_index];
            mapped[entry.target_index] = [
                entry.sign * source[0],
                entry.sign * source[1],
                entry.sign * source[2],
                entry.sign * source[3],
            ];
        }
        mapped_batch.push(mapped);
    }
    Ok(mapped_batch)
}

fn validate_input_crossing_map(
    expected_legs: usize,
    input_crossing_map: Option<&[InputCrossingMapEntry]>,
) -> PyResult<Option<&[InputCrossingMapEntry]>> {
    let Some(map) = input_crossing_map else {
        return Ok(None);
    };
    if map.len() != expected_legs {
        return Err(PyValueError::new_err(format!(
            "input crossing map has {} entries, expected {expected_legs}",
            map.len()
        )));
    }
    let mut seen = vec![false; expected_legs];
    for entry in map {
        if entry.target_index >= expected_legs || entry.source_index >= expected_legs {
            return Err(PyValueError::new_err(
                "input crossing map references an out-of-range external leg",
            ));
        }
        if seen[entry.target_index] {
            return Err(PyValueError::new_err(
                "input crossing map contains a duplicate target index",
            ));
        }
        seen[entry.target_index] = true;
    }
    if seen.iter().any(|value| !*value) {
        return Err(PyValueError::new_err(
            "input crossing map does not cover every target index",
        ));
    }
    Ok(Some(map))
}

fn apply_input_crossing_map_generic<T>(
    batch: &[Vec<[T; 4]>],
    expected_legs: usize,
    input_crossing_map: Option<&[InputCrossingMapEntry]>,
) -> PyResult<Vec<Vec<[T; 4]>>>
where
    T: RusticolHighPrecisionNumber,
    Complex<T>: Real + EvaluationDomain,
{
    let Some(map) = validate_input_crossing_map(expected_legs, input_crossing_map)? else {
        return Ok(batch.to_vec());
    };
    let mut mapped_batch = Vec::with_capacity(batch.len());
    for point in batch {
        let mut mapped = vec![std::array::from_fn(|_| T::new_zero()); expected_legs];
        for entry in map {
            let source = &point[entry.source_index];
            for component in 0..4 {
                mapped[entry.target_index][component] =
                    T::from(entry.sign) * source[component].clone();
            }
        }
        mapped_batch.push(mapped);
    }
    Ok(mapped_batch)
}

fn apply_lc_topology_label_permutation(
    batch: &[Vec<[f64; 4]>],
    expected_legs: usize,
    mapping: &[(usize, usize)],
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    let mut seen = vec![false; expected_legs];
    for (representative_index, sector_index) in mapping {
        if *representative_index >= expected_legs || *sector_index >= expected_legs {
            return Err(PyValueError::new_err(
                "LC topology replay label permutation references an out-of-range external leg",
            ));
        }
        if seen[*representative_index] {
            return Err(PyValueError::new_err(
                "LC topology replay label permutation contains a duplicate representative label",
            ));
        }
        seen[*representative_index] = true;
    }
    let mut mapped_batch = Vec::with_capacity(batch.len());
    for point in batch {
        if point.len() != expected_legs {
            return Err(PyValueError::new_err(format!(
                "LC topology replay point has {} external legs, expected {expected_legs}",
                point.len(),
            )));
        }
        let mut mapped = point.clone();
        for (representative_index, sector_index) in mapping {
            mapped[*representative_index] = point[*sector_index];
        }
        mapped_batch.push(mapped);
    }
    Ok(mapped_batch)
}

fn apply_lc_topology_label_permutations(
    batch: &[Vec<[f64; 4]>],
    expected_legs: usize,
    mappings: &[Vec<(usize, usize)>],
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    let mut expanded_batch = Vec::with_capacity(batch.len() * mappings.len());
    for mapping in mappings {
        expanded_batch.extend(apply_lc_topology_label_permutation(
            batch,
            expected_legs,
            mapping,
        )?);
    }
    Ok(expanded_batch)
}

fn apply_lc_topology_label_permutation_generic<T>(
    batch: &[Vec<[T; 4]>],
    expected_legs: usize,
    mapping: &[(usize, usize)],
) -> PyResult<Vec<Vec<[T; 4]>>>
where
    T: RusticolHighPrecisionNumber,
    Complex<T>: Real + EvaluationDomain,
{
    let mut seen = vec![false; expected_legs];
    for (representative_index, sector_index) in mapping {
        if *representative_index >= expected_legs || *sector_index >= expected_legs {
            return Err(PyValueError::new_err(
                "LC topology replay label permutation references an out-of-range external leg",
            ));
        }
        if seen[*representative_index] {
            return Err(PyValueError::new_err(
                "LC topology replay label permutation contains a duplicate representative label",
            ));
        }
        seen[*representative_index] = true;
    }
    let mut mapped_batch = Vec::with_capacity(batch.len());
    for point in batch {
        if point.len() != expected_legs {
            return Err(PyValueError::new_err(format!(
                "LC topology replay point has {} external legs, expected {expected_legs}",
                point.len(),
            )));
        }
        let mut mapped = point.clone();
        for (representative_index, sector_index) in mapping {
            mapped[*representative_index] = point[*sector_index].clone();
        }
        mapped_batch.push(mapped);
    }
    Ok(mapped_batch)
}

fn apply_lc_topology_label_permutations_generic<T>(
    batch: &[Vec<[T; 4]>],
    expected_legs: usize,
    mappings: &[Vec<(usize, usize)>],
) -> PyResult<Vec<Vec<[T; 4]>>>
where
    T: RusticolHighPrecisionNumber,
    Complex<T>: Real + EvaluationDomain,
{
    let mut expanded_batch = Vec::with_capacity(batch.len() * mappings.len());
    for mapping in mappings {
        expanded_batch.extend(apply_lc_topology_label_permutation_generic(
            batch,
            expected_legs,
            mapping,
        )?);
    }
    Ok(expanded_batch)
}

fn replay_mappings_per_expanded_batch(n_points: usize) -> usize {
    if n_points == 0 {
        return 1;
    }
    (MAX_LC_TOPOLOGY_REPLAY_EXPANDED_POINTS / n_points).max(1)
}

fn batch_momenta_dynamic(
    py: Python<'_>,
    momenta: &Bound<'_, PyAny>,
    expected_legs: usize,
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    match PyBuffer::<f64>::get(momenta) {
        Ok(buffer) => batch_momenta_dynamic_from_buffer(py, &buffer, expected_legs),
        Err(buffer_error) => match momenta.extract::<Vec<Vec<Vec<f64>>>>() {
            Ok(points) => batch_momenta_dynamic_from_nested(points, expected_legs),
            Err(sequence_error) => Err(PyValueError::new_err(format!(
                "momenta must be a C-contiguous f64 buffer or nested sequence with shape \
                 (batch, {expected_legs}, 4); buffer error: {buffer_error}; sequence error: \
                 {sequence_error}",
            ))),
        },
    }
}

fn batch_momenta_dynamic_from_buffer(
    py: Python<'_>,
    buffer: &PyBuffer<f64>,
    expected_legs: usize,
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    let shape = buffer.shape();
    if shape.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "momenta must have shape (batch, n_external, 4), got {shape:?}"
        )));
    }
    if shape[2] != 4 {
        return Err(PyValueError::new_err(format!(
            "last momenta axis must have length 4, got {}",
            shape[2]
        )));
    }
    if shape[1] != expected_legs {
        return Err(PyValueError::new_err(format!(
            "momenta must have shape (batch, {expected_legs}, 4), got {shape:?}",
        )));
    }
    let values = buffer
        .as_slice(py)
        .ok_or_else(|| PyValueError::new_err("momenta buffer must be C-contiguous f64 data"))?;
    let mut out = Vec::with_capacity(shape[0]);
    for row in 0..shape[0] {
        let mut point = Vec::with_capacity(shape[1]);
        for leg in 0..shape[1] {
            let mut momentum = [0.0; 4];
            for component in 0..4 {
                let index = (row * shape[1] + leg) * 4 + component;
                momentum[component] = values[index].get();
            }
            point.push(momentum);
        }
        out.push(point);
    }
    Ok(out)
}

fn batch_momenta_dynamic_from_nested(
    points: Vec<Vec<Vec<f64>>>,
    expected_legs: usize,
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    let mut out = Vec::with_capacity(points.len());
    for (row_index, row) in points.into_iter().enumerate() {
        if row.len() != expected_legs {
            return Err(PyValueError::new_err(format!(
                "momenta row {row_index} has {} external legs, expected {expected_legs}",
                row.len(),
            )));
        }
        let mut point = Vec::with_capacity(row.len());
        for (leg_index, leg) in row.into_iter().enumerate() {
            if leg.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "momenta row {row_index} leg {leg_index} has length {}, expected 4",
                    leg.len()
                )));
            }
            point.push([leg[0], leg[1], leg[2], leg[3]]);
        }
        out.push(point);
    }
    Ok(out)
}

fn batch_momenta_from_buffer(
    py: Python<'_>,
    buffer: &PyBuffer<f64>,
    expected_legs: usize,
) -> PyResult<Vec<[[f64; 4]; 16]>> {
    let shape = buffer.shape();
    if shape.len() != 3 {
        return Err(PyValueError::new_err(format!(
            "momenta must have shape (batch, n_external, 4), got {shape:?}"
        )));
    }
    if shape[2] != 4 {
        return Err(PyValueError::new_err(format!(
            "last momenta axis must have length 4, got {}",
            shape[2]
        )));
    }
    if shape[1] != expected_legs {
        return Err(PyValueError::new_err(format!(
            "momenta must have shape (batch, {expected_legs}, 4), got {shape:?}",
        )));
    }
    if shape[1] > 16 {
        return Err(PyValueError::new_err(
            "rusticol currently supports at most 16 external legs in the fixed f64 path",
        ));
    }
    let values = buffer
        .as_slice(py)
        .ok_or_else(|| PyValueError::new_err("momenta buffer must be C-contiguous f64 data"))?;
    let mut out = Vec::with_capacity(shape[0]);
    for row in 0..shape[0] {
        let mut point = [[0.0; 4]; 16];
        for leg in 0..shape[1] {
            for component in 0..4 {
                let index = (row * shape[1] + leg) * 4 + component;
                point[leg][component] = values[index].get();
            }
        }
        out.push(point);
    }
    Ok(out)
}

fn batch_momenta_from_nested(
    points: Vec<Vec<Vec<f64>>>,
    expected_legs: usize,
) -> PyResult<Vec<[[f64; 4]; 16]>> {
    let mut out = Vec::with_capacity(points.len());
    for (row_index, row) in points.into_iter().enumerate() {
        if row.len() != expected_legs {
            return Err(PyValueError::new_err(format!(
                "momenta row {row_index} has {} external legs, expected {expected_legs}",
                row.len(),
            )));
        }
        if row.len() > 16 {
            return Err(PyValueError::new_err(
                "rusticol currently supports at most 16 external legs in the fixed f64 path",
            ));
        }
        let mut point = [[0.0; 4]; 16];
        for (leg_index, leg) in row.into_iter().enumerate() {
            if leg.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "momenta row {row_index} leg {leg_index} has length {}, expected 4",
                    leg.len()
                )));
            }
            point[leg_index] = [leg[0], leg[1], leg[2], leg[3]];
        }
        out.push(point);
    }
    Ok(out)
}

fn batch_momenta_double(
    momenta: &Bound<'_, PyAny>,
    expected_legs: usize,
) -> PyResult<Vec<Vec<[DoubleFloat; 4]>>> {
    let float_points = batch_momenta_float(momenta, 106, expected_legs)?;
    Ok(float_points
        .into_iter()
        .map(|point| {
            point
                .into_iter()
                .map(|leg| {
                    [
                        leg[0].to_double_float(),
                        leg[1].to_double_float(),
                        leg[2].to_double_float(),
                        leg[3].to_double_float(),
                    ]
                })
                .collect()
        })
        .collect())
}

fn batch_momenta_float(
    momenta: &Bound<'_, PyAny>,
    binary_precision: u32,
    expected_legs: usize,
) -> PyResult<Vec<Vec<[Float; 4]>>> {
    let mut out = Vec::new();
    let rows = momenta.try_iter().map_err(|err| {
        PyValueError::new_err(format!(
            "high-precision momenta must be a nested sequence with shape \
             (batch, n_external, 4): {err}"
        ))
    })?;
    for (row_index, row_obj) in rows.enumerate() {
        let row = row_obj?;
        let mut point = Vec::new();
        let legs = row.try_iter().map_err(|err| {
            PyValueError::new_err(format!(
                "momenta row {row_index} is not iterable over legs: {err}"
            ))
        })?;
        for (leg_index, leg_obj) in legs.enumerate() {
            if leg_index >= 16 {
                return Err(PyValueError::new_err(
                    "rusticol currently supports at most 16 external legs",
                ));
            }
            let leg = leg_obj?;
            let components_iter = leg.try_iter().map_err(|err| {
                PyValueError::new_err(format!(
                    "momenta row {row_index} leg {leg_index} is not iterable: {err}"
                ))
            })?;
            let mut components = Vec::with_capacity(4);
            for component_obj in components_iter {
                let component = component_obj?;
                components.push(parse_float_component(&component, binary_precision)?);
            }
            if components.len() != 4 {
                return Err(PyValueError::new_err(format!(
                    "momenta row {row_index} leg {leg_index} has length {}, expected 4",
                    components.len()
                )));
            }
            let leg_array: [Float; 4] =
                components.try_into().map_err(|components: Vec<Float>| {
                    PyValueError::new_err(format!(
                        "momenta row {row_index} leg {leg_index} has length {}, expected 4",
                        components.len()
                    ))
                })?;
            point.push(leg_array);
        }
        if point.len() != expected_legs {
            return Err(PyValueError::new_err(format!(
                "momenta row {row_index} has {} external legs, expected {expected_legs}",
                point.len(),
            )));
        }
        out.push(point);
    }
    Ok(out)
}

fn parse_float_component(value: &Bound<'_, PyAny>, binary_precision: u32) -> PyResult<Float> {
    let text = value.str()?.to_string_lossy().into_owned();
    Float::parse(&text, Some(binary_precision)).map_err(|err| {
        PyValueError::new_err(format!(
            "could not parse high-precision momentum component {text:?}: {err}"
        ))
    })
}

fn decimal_digits_to_bits(decimal_digits: u32) -> u32 {
    (decimal_digits as f64 * std::f64::consts::LOG2_10).ceil() as u32
}

fn decimals_to_python<T: std::fmt::Display + std::fmt::LowerExp>(
    py: Python<'_>,
    values: Vec<T>,
    decimal_digits: u32,
) -> PyResult<Py<PyAny>> {
    let decimal = py.import("decimal")?.getattr("Decimal")?;
    let digits = decimal_digits as usize;
    let mut items = Vec::with_capacity(values.len());
    for value in values {
        items.push(decimal.call1((format!("{value:.digits$e}"),))?.unbind());
    }
    Ok(PyList::new(py, items)?.into_any().unbind())
}

fn f64_values_to_numpy_or_list(py: Python<'_>, values: Vec<f64>) -> PyResult<Py<PyAny>> {
    if let Ok(numpy) = py.import("numpy") {
        let dtype = numpy.getattr("float64")?;
        return numpy
            .getattr("array")?
            .call1((values,))?
            .call_method1("astype", (dtype,))?
            .into_py_any(py);
    }
    values.into_py_any(py)
}

fn checksum_dict_str<'py>(
    py: Python<'py>,
    stage: &str,
    checksum: ComplexChecksum,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = checksum_dict_base(py, checksum)?;
    dict.set_item("stage", stage)?;
    Ok(dict)
}

fn checksum_dict_usize<'py>(
    py: Python<'py>,
    stage: usize,
    checksum: ComplexChecksum,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = checksum_dict_base(py, checksum)?;
    dict.set_item("stage", stage)?;
    Ok(dict)
}

fn checksum_dict_base(py: Python<'_>, checksum: ComplexChecksum) -> PyResult<Bound<'_, PyDict>> {
    let dict = PyDict::new(py);
    dict.set_item("output_len", checksum.output_len)?;
    dict.set_item("sum_re", checksum.sum_re)?;
    dict.set_item("sum_im", checksum.sum_im)?;
    dict.set_item("sum_abs2", checksum.sum_abs2)?;
    dict.set_item("max_abs", checksum.max_abs)?;
    Ok(dict)
}

fn checksum_state_prefix(
    state: &[Complex<f64>],
    batch_size: usize,
    parameter_count: usize,
    prefix_len: usize,
) -> ComplexChecksum {
    let mut checksum = ComplexChecksum {
        output_len: batch_size * prefix_len,
        ..ComplexChecksum::default()
    };
    for row in 0..batch_size {
        let row_start = row * parameter_count;
        for value in &state[row_start..row_start + prefix_len] {
            checksum.add(*value);
        }
    }
    checksum
}

fn checksum_stage_outputs(
    state: &[Complex<f64>],
    batch_size: usize,
    parameter_count: usize,
    current_offsets: &[usize],
    outputs: &[(usize, usize, usize, usize)],
) -> ComplexChecksum {
    let mut checksum = ComplexChecksum {
        output_len: batch_size * outputs.len(),
        ..ComplexChecksum::default()
    };
    for row in 0..batch_size {
        let row_start = row * parameter_count;
        for (_, current_id, component, _) in outputs {
            let value = state[row_start + current_offsets[*current_id] + *component];
            checksum.add(value);
        }
    }
    checksum
}

impl ComplexChecksum {
    fn add(&mut self, value: Complex<f64>) {
        self.sum_re += value.re;
        self.sum_im += value.im;
        let abs2 = value.re * value.re + value.im * value.im;
        self.sum_abs2 += abs2;
        self.max_abs = self.max_abs.max(abs2.sqrt());
    }
}

fn memory_snapshot() -> MemorySnapshot {
    MemorySnapshot {
        current_rss_bytes: current_rss_bytes(),
        peak_rss_bytes: peak_rss_bytes(),
    }
}

#[cfg(target_os = "linux")]
fn current_rss_bytes() -> Option<u64> {
    let statm = fs::read_to_string("/proc/self/statm").ok()?;
    let resident_pages = statm.split_whitespace().nth(1)?.parse::<u64>().ok()?;
    let page_size = unsafe { libc::sysconf(libc::_SC_PAGESIZE) };
    if page_size <= 0 {
        return None;
    }
    resident_pages.checked_mul(page_size as u64)
}

#[cfg(not(target_os = "linux"))]
fn current_rss_bytes() -> Option<u64> {
    None
}

#[cfg(any(target_os = "linux", target_os = "macos"))]
fn peak_rss_bytes() -> Option<u64> {
    let mut usage = std::mem::MaybeUninit::<libc::rusage>::uninit();
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()) };
    if rc != 0 {
        return None;
    }
    let usage = unsafe { usage.assume_init() };
    let maxrss = usage.ru_maxrss;
    if maxrss <= 0 {
        return None;
    }
    #[cfg(target_os = "linux")]
    {
        Some((maxrss as u64).saturating_mul(1024))
    }
    #[cfg(target_os = "macos")]
    {
        Some(maxrss as u64)
    }
}

#[cfg(not(any(target_os = "linux", target_os = "macos")))]
fn peak_rss_bytes() -> Option<u64> {
    None
}

fn c64(re: f64, im: f64) -> Complex<f64> {
    Complex::new(re, im)
}

fn negate(momentum: [f64; 4]) -> [f64; 4] {
    [-momentum[0], -momentum[1], -momentum[2], -momentum[3]]
}

fn add_momenta(left: [f64; 4], right: [f64; 4]) -> [f64; 4] {
    [
        left[0] + right[0],
        left[1] + right[1],
        left[2] + right[2],
        left[3] + right[3],
    ]
}

fn minkowski_square(momentum: [f64; 4]) -> f64 {
    momentum[0] * momentum[0]
        - momentum[1] * momentum[1]
        - momentum[2] * momentum[2]
        - momentum[3] * momentum[3]
}

fn minkowski_dot_momentum(vector: &[Complex<f64>; 4], momentum: [f64; 4]) -> Complex<f64> {
    vector[0] * c64(momentum[0], 0.0)
        - vector[1] * c64(momentum[1], 0.0)
        - vector[2] * c64(momentum[2], 0.0)
        - vector[3] * c64(momentum[3], 0.0)
}

fn neutral_vector_propagator(
    vector: [Complex<f64>; 4],
    momentum: [f64; 4],
    vector_pdg: i32,
    mass: f64,
    width: f64,
) -> PyResult<Vec<Complex<f64>>> {
    if vector_pdg == 22 {
        let denominator = minkowski_square(momentum);
        if denominator == 0.0 {
            return Err(PyValueError::new_err("singular massless vector propagator"));
        }
        let prefactor = c64(0.0, -1.0) / c64(denominator, 0.0);
        return Ok(vector.into_iter().map(|value| value * prefactor).collect());
    }
    let denominator = c64(minkowski_square(momentum) - mass * mass, mass * width);
    if denominator.re == 0.0 && denominator.im == 0.0 {
        return Err(PyValueError::new_err("singular massive vector propagator"));
    }
    let prefactor = c64(0.0, -1.0) / denominator;
    let longitudinal = minkowski_dot_momentum(&vector, momentum) / c64(mass * mass, 0.0);
    Ok((0..4)
        .map(|component| {
            (vector[component] - c64(momentum[component], 0.0) * longitudinal) * prefactor
        })
        .collect())
}

fn lepton_antilepton_to_vector_weyl(
    lepton: Vec<Complex<f64>>,
    antilepton: Vec<Complex<f64>>,
    coupling: [f64; 2],
    lepton_chirality: i32,
    antilepton_chirality: i32,
) -> [Complex<f64>; 4] {
    let prefactor = c64(0.0, 1.0 / 2.0_f64.sqrt());
    let left = c64(coupling[0], 0.0);
    let right = c64(coupling[1], 0.0);
    let l1 = lepton[0];
    let l2 = lepton[1];
    let a1 = antilepton[0];
    let a2 = antilepton[1];
    if lepton_chirality == -1 && antilepton_chirality == 1 {
        let factor = prefactor * left;
        return [
            factor * (l1 * a1 + l2 * a2),
            -factor * (l2 * a1 + l1 * a2),
            c64(0.0, 1.0) * factor * (-l2 * a1 + l1 * a2),
            factor * (-l1 * a1 + l2 * a2),
        ];
    }
    if lepton_chirality == 1 && antilepton_chirality == -1 {
        let factor = prefactor * right;
        return [
            factor * (l1 * a1 + l2 * a2),
            factor * (l1 * a2 + l2 * a1),
            c64(0.0, 1.0) * factor * (-l1 * a2 + l2 * a1),
            factor * (l1 * a1 - l2 * a2),
        ];
    }
    [c64(0.0, 0.0); 4]
}

fn antilepton_lepton_to_vector_weyl(
    antilepton: Vec<Complex<f64>>,
    lepton: Vec<Complex<f64>>,
    coupling: [f64; 2],
    antilepton_chirality: i32,
    lepton_chirality: i32,
) -> [Complex<f64>; 4] {
    let prefactor = c64(0.0, 1.0 / 2.0_f64.sqrt());
    let left = c64(coupling[0], 0.0);
    let right = c64(coupling[1], 0.0);
    let a1 = antilepton[0];
    let a2 = antilepton[1];
    let l1 = lepton[0];
    let l2 = lepton[1];
    if antilepton_chirality == 1 && lepton_chirality == -1 {
        let factor = prefactor * left;
        return [
            factor * (l1 * a1 + l2 * a2),
            -factor * (l2 * a1 + l1 * a2),
            c64(0.0, 1.0) * factor * (-l2 * a1 + l1 * a2),
            factor * (-l1 * a1 + l2 * a2),
        ];
    }
    if antilepton_chirality == -1 && lepton_chirality == 1 {
        let factor = prefactor * right;
        return [
            factor * (l1 * a1 + l2 * a2),
            factor * (l1 * a2 + l2 * a1),
            c64(0.0, 1.0) * factor * (-l1 * a2 + l2 * a1),
            factor * (l1 * a1 - l2 * a2),
        ];
    }
    [c64(0.0, 0.0); 4]
}

fn neutral_lepton_pair_current(
    first_pdg: i32,
    first_momentum: [f64; 4],
    first_helicity: i32,
    first_chirality: i32,
    second_pdg: i32,
    second_momentum: [f64; 4],
    second_helicity: i32,
    second_chirality: i32,
    coupling: [f64; 2],
) -> PyResult<[Complex<f64>; 4]> {
    if first_pdg > 0 && second_pdg < 0 && first_pdg == -second_pdg {
        return Ok(lepton_antilepton_to_vector_weyl(
            ext_quark_weyl(first_momentum, first_helicity, first_chirality),
            ext_antiquark_weyl(second_momentum, second_helicity, second_chirality),
            coupling,
            first_chirality,
            second_chirality,
        ));
    }
    if first_pdg < 0 && second_pdg > 0 && second_pdg == -first_pdg {
        return Ok(antilepton_lepton_to_vector_weyl(
            ext_antiquark_weyl(first_momentum, first_helicity, first_chirality),
            ext_quark_weyl(second_momentum, second_helicity, second_chirality),
            coupling,
            first_chirality,
            second_chirality,
        ));
    }
    Err(PyValueError::new_err(format!(
        "neutral lepton-pair vector source expects l- l+, got {first_pdg} {second_pdg}"
    )))
}

#[allow(clippy::too_many_arguments)]
fn charged_lepton_pair_current(
    first_pdg: i32,
    first_momentum: [f64; 4],
    first_helicity: i32,
    first_chirality: i32,
    second_pdg: i32,
    second_momentum: [f64; 4],
    second_helicity: i32,
    second_chirality: i32,
    vector_pdg: i32,
    coupling: [f64; 2],
) -> PyResult<[Complex<f64>; 4]> {
    if vector_pdg == 24 {
        if matches!(first_pdg, -11 | -13 | -15) && matches!(second_pdg, 12 | 14 | 16) {
            return Ok(antilepton_lepton_to_vector_weyl(
                ext_antiquark_weyl(first_momentum, first_helicity, first_chirality),
                ext_quark_weyl(second_momentum, second_helicity, second_chirality),
                coupling,
                first_chirality,
                second_chirality,
            ));
        }
        if matches!(first_pdg, 12 | 14 | 16) && matches!(second_pdg, -11 | -13 | -15) {
            return Ok(lepton_antilepton_to_vector_weyl(
                ext_quark_weyl(first_momentum, first_helicity, first_chirality),
                ext_antiquark_weyl(second_momentum, second_helicity, second_chirality),
                coupling,
                first_chirality,
                second_chirality,
            ));
        }
        return Err(PyValueError::new_err(format!(
            "W+ lepton-pair vector source expects l+ nu, got {first_pdg} {second_pdg}"
        )));
    }
    if vector_pdg == -24 {
        if matches!(first_pdg, 11 | 13 | 15) && matches!(second_pdg, -12 | -14 | -16) {
            return Ok(lepton_antilepton_to_vector_weyl(
                ext_quark_weyl(first_momentum, first_helicity, first_chirality),
                ext_antiquark_weyl(second_momentum, second_helicity, second_chirality),
                coupling,
                first_chirality,
                second_chirality,
            ));
        }
        if matches!(first_pdg, -12 | -14 | -16) && matches!(second_pdg, 11 | 13 | 15) {
            return Ok(antilepton_lepton_to_vector_weyl(
                ext_antiquark_weyl(first_momentum, first_helicity, first_chirality),
                ext_quark_weyl(second_momentum, second_helicity, second_chirality),
                coupling,
                first_chirality,
                second_chirality,
            ));
        }
        return Err(PyValueError::new_err(format!(
            "W- lepton-pair vector source expects l- nu~, got {first_pdg} {second_pdg}"
        )));
    }
    Err(PyValueError::new_err(format!(
        "charged lepton-pair vector source expects W+/W-, got {vector_pdg}"
    )))
}

fn fortran_sign(value: f64, sign_source: f64) -> f64 {
    value.abs().copysign(sign_source)
}

fn complex_zero<T>() -> Complex<T>
where
    T: Real + Clone,
{
    Complex::new(T::new_zero(), T::new_zero())
}

fn c_generic<T>(re: T, im: T) -> Complex<T> {
    Complex::new(re, im)
}

fn negate_generic<T>(momentum: &[T; 4]) -> [T; 4]
where
    T: Real + Clone,
{
    [
        -momentum[0].clone(),
        -momentum[1].clone(),
        -momentum[2].clone(),
        -momentum[3].clone(),
    ]
}

fn is_zero<T: RealLike>(value: &T) -> bool {
    value.to_f64() == 0.0
}

fn t_min<T>(left: &T, right: &T) -> T
where
    T: RealLike + PartialOrd + Clone,
{
    if left <= right {
        left.clone()
    } else {
        right.clone()
    }
}

fn t_max<T>(left: &T, right: &T) -> T
where
    T: RealLike + PartialOrd + Clone,
{
    if left >= right {
        left.clone()
    } else {
        right.clone()
    }
}

fn fortran_sign_generic<T>(value: &T, sign_source: &T) -> T
where
    T: Real + RealLike + Clone,
{
    let magnitude = value.norm();
    if sign_source.to_f64().is_sign_negative() {
        -magnitude
    } else {
        magnitude
    }
}

fn ext_quark_dirac(momentum: [f64; 4], helicity: i32) -> [Complex<f64>; 4] {
    let [energy, px, py, pz] = momentum;
    if energy > 0.0 {
        let sqp0p3 = if px == 0.0 && py == 0.0 && pz < 0.0 {
            0.0
        } else {
            (energy + pz).max(0.0).sqrt()
        };
        let chi1 = c64(sqp0p3, 0.0);
        let chi2 = if sqp0p3 == 0.0 {
            c64(-(helicity as f64) * (2.0 * energy).sqrt(), 0.0)
        } else {
            c64(helicity as f64 * px / sqp0p3, -py / sqp0p3)
        };
        if helicity == 1 {
            return [chi1, chi2, c64(0.0, 0.0), c64(0.0, 0.0)];
        }
        return [c64(0.0, 0.0), c64(0.0, 0.0), chi2, chi1];
    }
    let sqp0p3 = if px == 0.0 && py == 0.0 && pz > 0.0 {
        0.0
    } else {
        -(-(energy + pz)).max(0.0).sqrt()
    };
    let chi1 = c64(sqp0p3, 0.0);
    let chi2 = if sqp0p3 == 0.0 {
        c64(-(helicity as f64) * (2.0 * energy.abs()).sqrt(), 0.0)
    } else {
        c64(-(helicity as f64) * (-px) / sqp0p3, py / sqp0p3)
    };
    if -helicity == 1 {
        [chi1, chi2, c64(0.0, 0.0), c64(0.0, 0.0)]
    } else {
        [c64(0.0, 0.0), c64(0.0, 0.0), chi2, chi1]
    }
}

fn ext_antiquark_dirac(momentum: [f64; 4], helicity: i32) -> [Complex<f64>; 4] {
    let [energy, px, py, pz] = momentum;
    if energy > 0.0 {
        let sqp0p3 = if px == 0.0 && py == 0.0 && pz < 0.0 {
            0.0
        } else {
            -(energy + pz).max(0.0).sqrt()
        };
        let chi1 = c64(sqp0p3, 0.0);
        let chi2 = if sqp0p3 == 0.0 {
            c64(-(helicity as f64) * (2.0 * energy).sqrt(), 0.0)
        } else {
            c64(-(helicity as f64) * px / sqp0p3, py / sqp0p3)
        };
        if -helicity == 1 {
            return [c64(0.0, 0.0), c64(0.0, 0.0), chi1, chi2];
        }
        return [chi2, chi1, c64(0.0, 0.0), c64(0.0, 0.0)];
    }
    let sqp0p3 = if px == 0.0 && py == 0.0 && pz > 0.0 {
        0.0
    } else {
        (-(energy + pz)).max(0.0).sqrt()
    };
    let chi1 = c64(sqp0p3, 0.0);
    let chi2 = if sqp0p3 == 0.0 {
        c64(-(helicity as f64) * (2.0 * energy.abs()).sqrt(), 0.0)
    } else {
        c64(helicity as f64 * (-px) / sqp0p3, -py / sqp0p3)
    };
    if helicity == 1 {
        [c64(0.0, 0.0), c64(0.0, 0.0), chi1, chi2]
    } else {
        [chi2, chi1, c64(0.0, 0.0), c64(0.0, 0.0)]
    }
}

fn ext_quark_dirac_massive(momentum: [f64; 4], helicity: i32, mass: f64) -> [Complex<f64>; 4] {
    if mass.abs() < 1.0e-8 {
        return ext_quark_dirac(momentum, helicity);
    }
    let [energy, px, py, pz] = momentum;
    let nsf = if energy > 0.0 { 1 } else { -1 };
    let nh = nsf * helicity;
    let pp = (px * px + py * py + pz * pz).sqrt().abs();
    let omega1 = (energy.abs() + pp).sqrt();
    let omega2 = mass / omega1;
    let omega = [omega1, omega2];
    let sf1 = (1 + nsf + (1 - nsf) * nh) as f64 * 0.5;
    let sf2 = (1 + nsf - (1 - nsf) * nh) as f64 * 0.5;
    let ip = ((3 + nh) / 2 - 1) as usize;
    let im = ((3 - nh) / 2 - 1) as usize;
    let sfomeg = [sf1 * omega[ip], sf2 * omega[im]];
    let (signed_px, signed_py, signed_pz) = if energy > 0.0 {
        (px, py, pz)
    } else {
        (-px, -py, -pz)
    };
    let pp3 = (pp + signed_pz).max(0.0);
    let chi1 = if pp == 0.0 {
        c64(1.0, 0.0)
    } else {
        c64((pp3 * 0.5 / pp).sqrt(), 0.0)
    };
    let chi2 = if pp3 == 0.0 || pp == 0.0 {
        c64(-(nh as f64), 0.0)
    } else {
        let denom = (2.0 * pp * pp3).sqrt();
        c64((nh as f64) * signed_px / denom, -signed_py / denom)
    };
    let chi = [chi1, chi2];
    [
        chi[im] * sfomeg[1],
        chi[ip] * sfomeg[1],
        chi[im] * sfomeg[0],
        chi[ip] * sfomeg[0],
    ]
}

fn ext_antiquark_dirac_massive(momentum: [f64; 4], helicity: i32, mass: f64) -> [Complex<f64>; 4] {
    if mass.abs() < 1.0e-8 {
        return ext_antiquark_dirac(momentum, helicity);
    }
    let [energy, px, py, pz] = momentum;
    let nsf = if energy > 0.0 { -1 } else { 1 };
    let nh = nsf * helicity;
    let pp = (px * px + py * py + pz * pz).sqrt().abs();
    let omega1 = (energy.abs() + pp).sqrt();
    let omega2 = mass / omega1;
    let omega = [omega1, omega2];
    let sf1 = (1 + nsf + (1 - nsf) * nh) as f64 * 0.5;
    let sf2 = (1 + nsf - (1 - nsf) * nh) as f64 * 0.5;
    let ip = ((3 + nh) / 2 - 1) as usize;
    let im = ((3 - nh) / 2 - 1) as usize;
    let sfomeg = [sf1 * omega[ip], sf2 * omega[im]];
    let (signed_px, signed_py, signed_pz) = if energy > 0.0 {
        (px, py, pz)
    } else {
        (-px, -py, -pz)
    };
    let pp3 = (pp + signed_pz).max(0.0);
    let chi1 = if pp == 0.0 {
        c64(1.0, 0.0)
    } else {
        c64((pp3 * 0.5 / pp).sqrt(), 0.0)
    };
    let chi2 = if pp3 == 0.0 || pp == 0.0 {
        c64(-(nh as f64), 0.0)
    } else {
        let denom = (2.0 * pp * pp3).sqrt();
        c64((nh as f64) * signed_px / denom, signed_py / denom)
    };
    let chi = [chi1, chi2];
    [
        chi[im] * sfomeg[0],
        chi[ip] * sfomeg[0],
        chi[im] * sfomeg[1],
        chi[ip] * sfomeg[1],
    ]
}

fn ext_quark_dirac_generic<T>(momentum: &[T; 4], helicity: i32) -> [Complex<T>; 4]
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    let energy = momentum[0].clone();
    let px = momentum[1].clone();
    let py = momentum[2].clone();
    let pz = momentum[3].clone();
    if energy.to_f64() > 0.0 {
        let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() < 0.0 {
            T::new_zero()
        } else {
            t_max(&(energy.clone() + pz.clone()), &T::new_zero()).sqrt()
        };
        let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
        let chi2 = if is_zero(&sqp0p3) {
            c_generic(
                -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.clone()).sqrt(),
                T::new_zero(),
            )
        } else {
            c_generic(
                energy.from_i64(helicity as i64) * px.clone() / sqp0p3.clone(),
                -py.clone() / sqp0p3.clone(),
            )
        };
        if helicity == 1 {
            return [chi1, chi2, complex_zero(), complex_zero()];
        }
        return [complex_zero(), complex_zero(), chi2, chi1];
    }
    let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() > 0.0 {
        T::new_zero()
    } else {
        -t_max(&(-(energy.clone() + pz.clone())), &T::new_zero()).sqrt()
    };
    let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
    let chi2 = if is_zero(&sqp0p3) {
        c_generic(
            -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.norm()).sqrt(),
            T::new_zero(),
        )
    } else {
        c_generic(
            -energy.from_i64(helicity as i64) * (-px.clone()) / sqp0p3.clone(),
            py.clone() / sqp0p3.clone(),
        )
    };
    if -helicity == 1 {
        [chi1, chi2, complex_zero(), complex_zero()]
    } else {
        [complex_zero(), complex_zero(), chi2, chi1]
    }
}

fn ext_antiquark_dirac_generic<T>(momentum: &[T; 4], helicity: i32) -> [Complex<T>; 4]
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    let energy = momentum[0].clone();
    let px = momentum[1].clone();
    let py = momentum[2].clone();
    let pz = momentum[3].clone();
    if energy.to_f64() > 0.0 {
        let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() < 0.0 {
            T::new_zero()
        } else {
            -t_max(&(energy.clone() + pz.clone()), &T::new_zero()).sqrt()
        };
        let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
        let chi2 = if is_zero(&sqp0p3) {
            c_generic(
                -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.clone()).sqrt(),
                T::new_zero(),
            )
        } else {
            c_generic(
                -energy.from_i64(helicity as i64) * px.clone() / sqp0p3.clone(),
                py.clone() / sqp0p3.clone(),
            )
        };
        if -helicity == 1 {
            return [complex_zero(), complex_zero(), chi1, chi2];
        }
        return [chi2, chi1, complex_zero(), complex_zero()];
    }
    let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() > 0.0 {
        T::new_zero()
    } else {
        t_max(&(-(energy.clone() + pz.clone())), &T::new_zero()).sqrt()
    };
    let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
    let chi2 = if is_zero(&sqp0p3) {
        c_generic(
            -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.norm()).sqrt(),
            T::new_zero(),
        )
    } else {
        c_generic(
            energy.from_i64(helicity as i64) * (-px.clone()) / sqp0p3.clone(),
            -py.clone() / sqp0p3.clone(),
        )
    };
    if helicity == 1 {
        [complex_zero(), complex_zero(), chi1, chi2]
    } else {
        [chi2, chi1, complex_zero(), complex_zero()]
    }
}

fn ext_massive_vector(momentum: [f64; 4], helicity: i32, mass: f64) -> [Complex<f64>; 4] {
    let [energy, px, py, pz] = momentum;
    let sqh = 0.5f64.sqrt();
    let hel = helicity as f64;
    let nsvahl = helicity.abs() as f64;
    let pt2 = px * px + py * py;
    let pp = energy.min((pt2 + pz * pz).sqrt());
    let pt = pp.min(pt2.sqrt());
    let hel0 = 1.0 - hel.abs();
    if pp == 0.0 {
        return [
            c64(0.0, 0.0),
            c64(-hel * sqh, 0.0),
            c64(0.0, nsvahl * sqh),
            c64(hel0, 0.0),
        ];
    }
    let emp = energy / (mass * pp);
    let wf0 = c64(hel0 * pp / mass, 0.0);
    let wf3 = c64(hel0 * pz * emp + hel * pt / pp * sqh, 0.0);
    let (wf1, wf2) = if pt != 0.0 {
        let pzpt = pz / (pp * pt) * sqh * hel;
        (
            c64(hel0 * px * emp - px * pzpt, -nsvahl * py / pt * sqh),
            c64(hel0 * py * emp - py * pzpt, nsvahl * px / pt * sqh),
        )
    } else {
        (
            c64(-hel * sqh, 0.0),
            c64(0.0, nsvahl * fortran_sign(sqh, pz)),
        )
    };
    [wf0, wf1, wf2, wf3]
}

fn ext_massive_vector_generic<T>(momentum: &[T; 4], helicity: i32, mass: T) -> [Complex<T>; 4]
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    let energy = momentum[0].clone();
    let px = momentum[1].clone();
    let py = momentum[2].clone();
    let pz = momentum[3].clone();
    let sqh = (energy.one() / energy.from_i64(2)).sqrt();
    let hel = energy.from_i64(helicity as i64);
    let nsvahl = energy.from_i64(helicity.abs() as i64);
    let pt2 = px.clone() * px.clone() + py.clone() * py.clone();
    let pp = t_min(&energy, &(pt2.clone() + pz.clone() * pz.clone()).sqrt());
    let pt = t_min(&pp, &pt2.sqrt());
    let hel0 = if helicity == 0 {
        energy.one()
    } else {
        T::new_zero()
    };
    if is_zero(&pp) {
        return [
            complex_zero(),
            c_generic(-hel.clone() * sqh.clone(), T::new_zero()),
            c_generic(T::new_zero(), nsvahl.clone() * sqh.clone()),
            c_generic(hel0, T::new_zero()),
        ];
    }
    let emp = energy.clone() / (mass.clone() * pp.clone());
    let wf0 = c_generic(hel0.clone() * pp.clone() / mass, T::new_zero());
    let wf3 = c_generic(
        hel0.clone() * pz.clone() * emp.clone()
            + hel.clone() * pt.clone() / pp.clone() * sqh.clone(),
        T::new_zero(),
    );
    let (wf1, wf2) = if !is_zero(&pt) {
        let pzpt = pz.clone() / (pp.clone() * pt.clone()) * sqh.clone() * hel.clone();
        (
            c_generic(
                hel0.clone() * px.clone() * emp.clone() - px.clone() * pzpt.clone(),
                -nsvahl.clone() * py.clone() / pt.clone() * sqh.clone(),
            ),
            c_generic(
                hel0 * py.clone() * emp - py.clone() * pzpt,
                nsvahl * px.clone() / pt * sqh,
            ),
        )
    } else {
        (
            c_generic(-hel * sqh.clone(), T::new_zero()),
            c_generic(T::new_zero(), nsvahl * fortran_sign_generic(&sqh, &pz)),
        )
    };
    [wf0, wf1, wf2, wf3]
}

fn ext_gluon(momentum: [f64; 4], helicity: i32) -> [Complex<f64>; 4] {
    let [energy, px, py, pz] = momentum;
    let sqh = 0.5f64.sqrt();
    if energy > 0.0 {
        let hel = helicity as f64;
        let pp = energy;
        let pt = (px * px + py * py).sqrt();
        let wf3 = c64(hel * pt / pp * sqh, 0.0);
        let (wf1, wf2) = if pt != 0.0 {
            let pzpt = pz / (pp * pt) * sqh * hel;
            (
                c64(-px * pzpt, -py / pt * sqh),
                c64(-py * pzpt, px / pt * sqh),
            )
        } else {
            (c64(-hel * sqh, 0.0), c64(0.0, fortran_sign(sqh, pz)))
        };
        return [c64(0.0, 0.0), wf1, wf2, wf3];
    }
    let hel = -helicity as f64;
    let pp = -energy;
    let pt = (px * px + py * py).sqrt();
    let wf3 = c64(hel * pt / pp * sqh, 0.0);
    let (wf1, wf2) = if pt != 0.0 {
        let pzpt = -pz / (pp * pt) * sqh * hel;
        (
            c64(px * pzpt, py / pt * sqh),
            c64(py * pzpt, -px / pt * sqh),
        )
    } else {
        (c64(-hel * sqh, 0.0), c64(0.0, -fortran_sign(sqh, pz)))
    };
    [c64(0.0, 0.0), wf1, wf2, wf3]
}

fn ext_gluon_generic<T>(momentum: &[T; 4], helicity: i32) -> [Complex<T>; 4]
where
    T: Real + RealLike + From<f64> + Clone,
{
    let energy = momentum[0].clone();
    let px = momentum[1].clone();
    let py = momentum[2].clone();
    let pz = momentum[3].clone();
    let sqh = (energy.one() / energy.from_i64(2)).sqrt();
    if energy.to_f64() > 0.0 {
        let hel = energy.from_i64(helicity as i64);
        let pp = energy;
        let pt = (px.clone() * px.clone() + py.clone() * py.clone()).sqrt();
        let wf3 = c_generic(
            hel.clone() * pt.clone() / pp.clone() * sqh.clone(),
            T::new_zero(),
        );
        let (wf1, wf2) = if !is_zero(&pt) {
            let pzpt = pz.clone() / (pp.clone() * pt.clone()) * sqh.clone() * hel.clone();
            (
                c_generic(
                    -px.clone() * pzpt.clone(),
                    -py.clone() / pt.clone() * sqh.clone(),
                ),
                c_generic(-py.clone() * pzpt, px.clone() / pt * sqh),
            )
        } else {
            (
                c_generic(-hel * sqh.clone(), T::new_zero()),
                c_generic(T::new_zero(), fortran_sign_generic(&sqh, &pz)),
            )
        };
        return [complex_zero(), wf1, wf2, wf3];
    }
    let hel = energy.from_i64(-(helicity as i64));
    let pp = -energy;
    let pt = (px.clone() * px.clone() + py.clone() * py.clone()).sqrt();
    let wf3 = c_generic(
        hel.clone() * pt.clone() / pp.clone() * sqh.clone(),
        T::new_zero(),
    );
    let (wf1, wf2) = if !is_zero(&pt) {
        let pzpt = -pz.clone() / (pp.clone() * pt.clone()) * sqh.clone() * hel.clone();
        (
            c_generic(
                px.clone() * pzpt.clone(),
                py.clone() / pt.clone() * sqh.clone(),
            ),
            c_generic(py.clone() * pzpt, -px.clone() / pt * sqh),
        )
    } else {
        (
            c_generic(-hel * sqh.clone(), T::new_zero()),
            c_generic(T::new_zero(), -fortran_sign_generic(&sqh, &pz)),
        )
    };
    [complex_zero(), wf1, wf2, wf3]
}

fn ext_quark_weyl_array(momentum: [f64; 4], helicity: i32, chirality: i32) -> [Complex<f64>; 2] {
    let [energy, px, py, pz] = momentum;
    if energy > 0.0 {
        let sqp0p3 = if px == 0.0 && py == 0.0 && pz < 0.0 {
            0.0
        } else {
            (energy + pz).max(0.0).sqrt()
        };
        let chi1 = c64(sqp0p3, 0.0);
        let chi2 = if sqp0p3 == 0.0 {
            c64(-(helicity as f64) * (2.0 * energy).sqrt(), 0.0)
        } else {
            c64(helicity as f64 * px / sqp0p3, -py / sqp0p3)
        };
        if helicity == 1 && chirality == 1 {
            return [chi1, chi2];
        }
        if helicity == -1 && chirality == -1 {
            return [chi2, chi1];
        }
        return [c64(0.0, 0.0), c64(0.0, 0.0)];
    }
    let sqp0p3 = if px == 0.0 && py == 0.0 && pz > 0.0 {
        0.0
    } else {
        -(-(energy + pz)).max(0.0).sqrt()
    };
    let chi1 = c64(sqp0p3, 0.0);
    let chi2 = if sqp0p3 == 0.0 {
        c64(-(helicity as f64) * (2.0 * energy.abs()).sqrt(), 0.0)
    } else {
        c64(-(helicity as f64) * (-px) / sqp0p3, py / sqp0p3)
    };
    if helicity == -1 && chirality == 1 {
        [chi1, chi2]
    } else if helicity == 1 && chirality == -1 {
        [chi2, chi1]
    } else {
        [c64(0.0, 0.0), c64(0.0, 0.0)]
    }
}

fn ext_quark_weyl(momentum: [f64; 4], helicity: i32, chirality: i32) -> Vec<Complex<f64>> {
    ext_quark_weyl_array(momentum, helicity, chirality).to_vec()
}

fn ext_quark_weyl_generic<T>(momentum: &[T; 4], helicity: i32, chirality: i32) -> Vec<Complex<T>>
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    let energy = momentum[0].clone();
    let px = momentum[1].clone();
    let py = momentum[2].clone();
    let pz = momentum[3].clone();
    if energy.to_f64() > 0.0 {
        let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() < 0.0 {
            T::new_zero()
        } else {
            t_max(&(energy.clone() + pz.clone()), &T::new_zero()).sqrt()
        };
        let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
        let chi2 = if is_zero(&sqp0p3) {
            c_generic(
                -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.clone()).sqrt(),
                T::new_zero(),
            )
        } else {
            c_generic(
                energy.from_i64(helicity as i64) * px.clone() / sqp0p3.clone(),
                -py.clone() / sqp0p3.clone(),
            )
        };
        if helicity == 1 && chirality == 1 {
            return vec![chi1, chi2];
        }
        if helicity == -1 && chirality == -1 {
            return vec![chi2, chi1];
        }
        return vec![complex_zero(), complex_zero()];
    }
    let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() > 0.0 {
        T::new_zero()
    } else {
        -t_max(&(-(energy.clone() + pz.clone())), &T::new_zero()).sqrt()
    };
    let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
    let chi2 = if is_zero(&sqp0p3) {
        c_generic(
            -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.norm()).sqrt(),
            T::new_zero(),
        )
    } else {
        c_generic(
            -energy.from_i64(helicity as i64) * (-px.clone()) / sqp0p3.clone(),
            py.clone() / sqp0p3.clone(),
        )
    };
    if helicity == -1 && chirality == 1 {
        vec![chi1, chi2]
    } else if helicity == 1 && chirality == -1 {
        vec![chi2, chi1]
    } else {
        vec![complex_zero(), complex_zero()]
    }
}

fn ext_antiquark_weyl_array(
    momentum: [f64; 4],
    helicity: i32,
    chirality: i32,
) -> [Complex<f64>; 2] {
    let [energy, px, py, pz] = momentum;
    if energy > 0.0 {
        let sqp0p3 = if px == 0.0 && py == 0.0 && pz < 0.0 {
            0.0
        } else {
            -(energy + pz).max(0.0).sqrt()
        };
        let chi1 = c64(sqp0p3, 0.0);
        let chi2 = if sqp0p3 == 0.0 {
            c64(-(helicity as f64) * (2.0 * energy).sqrt(), 0.0)
        } else {
            c64(-(helicity as f64) * px / sqp0p3, py / sqp0p3)
        };
        if helicity == 1 && chirality == 1 {
            return [chi2, chi1];
        }
        if helicity == -1 && chirality == -1 {
            return [chi1, chi2];
        }
        return [c64(0.0, 0.0), c64(0.0, 0.0)];
    }
    let sqp0p3 = if px == 0.0 && py == 0.0 && pz > 0.0 {
        0.0
    } else {
        (-(energy + pz)).max(0.0).sqrt()
    };
    let chi1 = c64(sqp0p3, 0.0);
    let chi2 = if sqp0p3 == 0.0 {
        c64(-(helicity as f64) * (2.0 * energy.abs()).sqrt(), 0.0)
    } else {
        c64(helicity as f64 * (-px) / sqp0p3, -py / sqp0p3)
    };
    if helicity == -1 && chirality == 1 {
        [chi2, chi1]
    } else if helicity == 1 && chirality == -1 {
        [chi1, chi2]
    } else {
        [c64(0.0, 0.0), c64(0.0, 0.0)]
    }
}

fn ext_antiquark_weyl(momentum: [f64; 4], helicity: i32, chirality: i32) -> Vec<Complex<f64>> {
    ext_antiquark_weyl_array(momentum, helicity, chirality).to_vec()
}

fn ext_antiquark_weyl_generic<T>(
    momentum: &[T; 4],
    helicity: i32,
    chirality: i32,
) -> Vec<Complex<T>>
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    let energy = momentum[0].clone();
    let px = momentum[1].clone();
    let py = momentum[2].clone();
    let pz = momentum[3].clone();
    if energy.to_f64() > 0.0 {
        let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() < 0.0 {
            T::new_zero()
        } else {
            -t_max(&(energy.clone() + pz.clone()), &T::new_zero()).sqrt()
        };
        let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
        let chi2 = if is_zero(&sqp0p3) {
            c_generic(
                -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.clone()).sqrt(),
                T::new_zero(),
            )
        } else {
            c_generic(
                -energy.from_i64(helicity as i64) * px.clone() / sqp0p3.clone(),
                py.clone() / sqp0p3.clone(),
            )
        };
        if helicity == 1 && chirality == 1 {
            return vec![chi2, chi1];
        }
        if helicity == -1 && chirality == -1 {
            return vec![chi1, chi2];
        }
        return vec![complex_zero(), complex_zero()];
    }
    let sqp0p3 = if is_zero(&px) && is_zero(&py) && pz.to_f64() > 0.0 {
        T::new_zero()
    } else {
        t_max(&(-(energy.clone() + pz.clone())), &T::new_zero()).sqrt()
    };
    let chi1 = c_generic(sqp0p3.clone(), T::new_zero());
    let chi2 = if is_zero(&sqp0p3) {
        c_generic(
            -energy.from_i64(helicity as i64) * (energy.from_i64(2) * energy.norm()).sqrt(),
            T::new_zero(),
        )
    } else {
        c_generic(
            energy.from_i64(helicity as i64) * (-px.clone()) / sqp0p3.clone(),
            -py.clone() / sqp0p3.clone(),
        )
    };
    if helicity == -1 && chirality == 1 {
        vec![chi2, chi1]
    } else if helicity == 1 && chirality == -1 {
        vec![chi1, chi2]
    } else {
        vec![complex_zero(), complex_zero()]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn generic_dynamic_momenta_parser_accepts_more_than_sixteen_legs() {
        let point = vec![vec![0.0, 1.0, 2.0, 3.0]; 17];
        let parsed = batch_momenta_dynamic_from_nested(vec![point], 17).unwrap();

        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].len(), 17);
        assert_eq!(parsed[0][16], [0.0, 1.0, 2.0, 3.0]);
    }

    #[test]
    fn legacy_fixed_momenta_parser_keeps_sixteen_leg_cap() {
        pyo3::Python::initialize();
        let point = vec![vec![0.0, 1.0, 2.0, 3.0]; 17];
        let error = batch_momenta_from_nested(vec![point], 17).unwrap_err();

        assert!(
            error
                .to_string()
                .contains("at most 16 external legs in the fixed f64 path")
        );
    }

    #[test]
    fn input_crossing_map_permutates_and_sign_flips_external_momenta() {
        let batch = vec![vec![
            [10.0, 1.0, 2.0, 3.0],
            [20.0, 4.0, 5.0, 6.0],
            [30.0, 7.0, 8.0, 9.0],
        ]];
        let map = vec![
            InputCrossingMapEntry {
                target_index: 0,
                source_index: 2,
                sign: -1.0,
            },
            InputCrossingMapEntry {
                target_index: 1,
                source_index: 0,
                sign: 1.0,
            },
            InputCrossingMapEntry {
                target_index: 2,
                source_index: 1,
                sign: -1.0,
            },
        ];

        let mapped = apply_input_crossing_map(&batch, 3, Some(&map)).unwrap();

        assert_eq!(mapped[0][0], [-30.0, -7.0, -8.0, -9.0]);
        assert_eq!(mapped[0][1], batch[0][0]);
        assert_eq!(mapped[0][2], [-20.0, -4.0, -5.0, -6.0]);
    }

    #[test]
    fn input_crossing_map_rejects_duplicate_or_incomplete_targets() {
        let batch = vec![vec![[0.0, 0.0, 0.0, 0.0]; 2]];
        let duplicate = vec![
            InputCrossingMapEntry {
                target_index: 0,
                source_index: 0,
                sign: 1.0,
            },
            InputCrossingMapEntry {
                target_index: 0,
                source_index: 1,
                sign: 1.0,
            },
        ];

        let error = apply_input_crossing_map(&batch, 2, Some(&duplicate)).unwrap_err();

        assert!(error.to_string().contains("duplicate target index"));
    }

    #[test]
    fn lc_topology_label_permutation_maps_sector_momenta_to_representative_slots() {
        let batch = vec![vec![
            [10.0, 1.0, 0.0, 0.0],
            [20.0, 2.0, 0.0, 0.0],
            [30.0, 3.0, 0.0, 0.0],
            [40.0, 4.0, 0.0, 0.0],
        ]];
        let mapping = vec![(0, 1), (1, 0), (2, 3), (3, 2)];

        let mapped = apply_lc_topology_label_permutation(&batch, 4, &mapping).unwrap();

        assert_eq!(mapped[0][0], batch[0][1]);
        assert_eq!(mapped[0][1], batch[0][0]);
        assert_eq!(mapped[0][2], batch[0][3]);
        assert_eq!(mapped[0][3], batch[0][2]);
    }

    #[test]
    fn lc_topology_label_permutation_rejects_duplicate_representative_slots() {
        let batch = vec![vec![[0.0, 0.0, 0.0, 0.0]; 2]];
        let duplicate = vec![(0, 0), (0, 1)];

        let error = apply_lc_topology_label_permutation(&batch, 2, &duplicate).unwrap_err();

        assert!(error.to_string().contains("duplicate representative label"));
    }

    fn minimal_generic_manifest() -> Value {
        json!({
            "schema_version": 2,
            "kind": "pyamplicol-generic-dag-process",
            "process": "d d~ > z",
            "key": "d_dbar_to_z",
            "color_accuracy": "lc",
            "external_pdg_order": [1, -1, 23],
            "lc_topology_reuse": {
                "available": true,
                "active_topology_group_count": 1,
                "representative_sector_ids": [0],
                "groups": []
            },
            "compiled": {
                "kind": "generic-dag-plan-only",
                "runtime_available": false,
                "runtime_unavailable_message": "pending generic runtime"
            },
            "dag_summary": {
                "current_count": 3,
                "source_count": 2,
                "interaction_count": 1,
                "amplitude_root_count": 1,
                "truncated": false
            },
            "runtime_schema": {
                "schema_version": 2,
                "kind": "pyamplicol-generic-dag-runtime-schema",
                "process_key": "d_dbar_to_z",
                "process": "d d~ > z",
                "external_particles": [
                    {
                        "label": 1,
                        "index": 0,
                        "pdg": 1,
                        "outgoing_pdg": -1,
                        "role": "initial",
                        "momentum_slot": 0
                    },
                    {
                        "label": 2,
                        "index": 1,
                        "pdg": -1,
                        "outgoing_pdg": 1,
                        "role": "initial",
                        "momentum_slot": 1
                    },
                    {
                        "label": 3,
                        "index": 2,
                        "pdg": 23,
                        "outgoing_pdg": 23,
                        "role": "final",
                        "momentum_slot": 2
                    }
                ],
                "parameter_layout": {
                    "source_component_parameter_count": 4,
                    "momentum_parameter_count": 12,
                    "parameter_count_if_flattened": 16,
                    "value_component_count": 8,
                    "source_components_complex": true,
                    "momentum_components_real": true,
                    "real_valued_inputs": [4,5,6,7,8,9,10,11,12,13,14,15]
                },
                "current_storage": {
                    "component_count": 8,
                    "number_type": "complex",
                    "current_slots": [
                        {
                            "current_id": 0,
                            "component_start": 0,
                            "component_stop": 2,
                            "dimension": 2,
                            "is_source": true,
                            "particle_id": -1,
                            "external_mask": 1,
                            "external_labels": [1],
                            "helicity_ancestry": 1,
                            "chirality": -1,
                            "spin_state": -1,
                            "flavour_flow": [-1],
                            "charge_flow": 1,
                            "color_state": {"accuracy": "lc", "sector_id": 0, "line_groups": [0], "basis_key": []},
                            "momentum_mask": 1,
                            "auxiliary_kind": null
                        },
                        {
                            "current_id": 1,
                            "component_start": 2,
                            "component_stop": 4,
                            "dimension": 2,
                            "is_source": true,
                            "particle_id": 1,
                            "external_mask": 2,
                            "external_labels": [2],
                            "helicity_ancestry": 2,
                            "chirality": 1,
                            "spin_state": 1,
                            "flavour_flow": [1],
                            "charge_flow": -1,
                            "color_state": {"accuracy": "lc", "sector_id": 0, "line_groups": [0], "basis_key": []},
                            "momentum_mask": 2,
                            "auxiliary_kind": null
                        },
                        {
                            "current_id": 2,
                            "component_start": 4,
                            "component_stop": 8,
                            "dimension": 4,
                            "is_source": false,
                            "particle_id": 23,
                            "external_mask": 3,
                            "external_labels": [1, 2],
                            "helicity_ancestry": 3,
                            "chirality": 0,
                            "spin_state": 0,
                            "flavour_flow": [-1, 1, 23],
                            "charge_flow": 0,
                            "color_state": {"accuracy": "lc", "sector_id": 0, "line_groups": [0], "basis_key": []},
                            "momentum_mask": 3,
                            "auxiliary_kind": null
                        }
                    ]
                },
                "value_storage": {
                    "component_count": 8,
                    "number_type": "complex",
                    "value_slots": [
                        {
                            "value_slot_id": 0,
                            "current_id": 0,
                            "variant": "source",
                            "component_start": 0,
                            "component_stop": 2,
                            "dimension": 2,
                            "current_component_start": 0,
                            "current_component_stop": 2,
                            "is_source": true,
                            "applies_propagator": false,
                            "particle_id": -1,
                            "external_mask": 1,
                            "external_labels": [1],
                            "momentum_mask": 1,
                            "chirality": -1
                        },
                        {
                            "value_slot_id": 1,
                            "current_id": 1,
                            "variant": "source",
                            "component_start": 2,
                            "component_stop": 4,
                            "dimension": 2,
                            "current_component_start": 2,
                            "current_component_stop": 4,
                            "is_source": true,
                            "applies_propagator": false,
                            "particle_id": 1,
                            "external_mask": 2,
                            "external_labels": [2],
                            "momentum_mask": 2,
                            "chirality": 1
                        },
                        {
                            "value_slot_id": 2,
                            "current_id": 2,
                            "variant": "unpropagated",
                            "component_start": 4,
                            "component_stop": 8,
                            "dimension": 4,
                            "current_component_start": 4,
                            "current_component_stop": 8,
                            "is_source": false,
                            "applies_propagator": false,
                            "particle_id": 23,
                            "external_mask": 3,
                            "external_labels": [1, 2],
                            "momentum_mask": 3,
                            "chirality": 0
                        }
                    ]
                },
                "source_fill": {
                    "source_count": 2,
                    "sources": [
                        {
                            "source_id": 0,
                            "current_id": 0,
                            "current_component_start": 0,
                            "current_component_stop": 2,
                            "value_slot": {
                                "value_slot_id": 0,
                                "current_id": 0,
                                "variant": "source",
                                "component_start": 0,
                                "component_stop": 2,
                                "dimension": 2
                            },
                            "source_parameter_start": 0,
                            "source_parameter_stop": 2,
                            "leg_label": 1,
                            "input_momentum_slot": 0,
                            "side": "initial",
                            "crossing": "negate-incoming-momentum",
                            "physical_pdg": 1,
                            "outgoing_pdg": -1,
                            "particle_id": -1,
                            "source_kind": "external-wavefunction",
                            "source_helicity": -1,
                            "chirality": -1,
                            "spin_state": -1,
                            "dimension": 2,
                            "helicity_ancestry": 1,
                            "color_state": {"accuracy": "lc", "sector_id": 0, "line_groups": [0], "basis_key": []}
                        },
                        {
                            "source_id": 1,
                            "current_id": 1,
                            "current_component_start": 2,
                            "current_component_stop": 4,
                            "value_slot": {
                                "value_slot_id": 1,
                                "current_id": 1,
                                "variant": "source",
                                "component_start": 2,
                                "component_stop": 4,
                                "dimension": 2
                            },
                            "source_parameter_start": 2,
                            "source_parameter_stop": 4,
                            "leg_label": 2,
                            "input_momentum_slot": 1,
                            "side": "initial",
                            "crossing": "negate-incoming-momentum",
                            "physical_pdg": -1,
                            "outgoing_pdg": 1,
                            "particle_id": 1,
                            "source_kind": "external-wavefunction",
                            "source_helicity": 1,
                            "chirality": 1,
                            "spin_state": 1,
                            "dimension": 2,
                            "helicity_ancestry": 2,
                            "color_state": {"accuracy": "lc", "sector_id": 0, "line_groups": [0], "basis_key": []}
                        }
                    ]
                },
                "momentum_slots": [
                    {
                        "momentum_slot_id": 0,
                        "momentum_mask": 1,
                        "external_labels": [1],
                        "component_start": 0,
                        "component_stop": 4,
                        "real_valued": true
                    },
                    {
                        "momentum_slot_id": 1,
                        "momentum_mask": 2,
                        "external_labels": [2],
                        "component_start": 4,
                        "component_stop": 8,
                        "real_valued": true
                    },
                    {
                        "momentum_slot_id": 2,
                        "momentum_mask": 3,
                        "external_labels": [1, 2],
                        "component_start": 8,
                        "component_stop": 12,
                        "real_valued": true
                    }
                ],
                "stages": [
                    {
                        "stage_index": 1,
                        "stage_kind": "current-combine",
                        "subset_size": 2,
                        "input_current_ids": [0, 1],
                        "output_current_ids": [2],
                        "input_value_slot_ids": [0, 1],
                        "output_value_slot_ids": [2],
                        "interaction_count": 1,
                        "interactions": [
                            {
                                "interaction_id": 0,
                                "vertex_kind": 9,
                                "vertex_particles": [-1, 1, 23],
                                "left_current_id": 0,
                                "right_current_id": 1,
                                "result_current_id": 2,
                                "left_slot": {
                                    "current_id": 0,
                                    "component_start": 0,
                                    "component_stop": 2,
                                    "dimension": 2
                                },
                                "right_slot": {
                                    "current_id": 1,
                                    "component_start": 2,
                                    "component_stop": 4,
                                    "dimension": 2
                                },
                                "result_slot": {
                                    "current_id": 2,
                                    "component_start": 4,
                                    "component_stop": 8,
                                    "dimension": 4
                                },
                                "left_value_slot": {
                                    "value_slot_id": 0,
                                    "current_id": 0,
                                    "variant": "source",
                                    "component_start": 0,
                                    "component_stop": 2,
                                    "dimension": 2
                                },
                                "right_value_slot": {
                                    "value_slot_id": 1,
                                    "current_id": 1,
                                    "variant": "source",
                                    "component_start": 2,
                                    "component_stop": 4,
                                    "dimension": 2
                                },
                                "result_value_slots": [
                                    {
                                        "value_slot_id": 2,
                                        "current_id": 2,
                                        "variant": "unpropagated",
                                        "component_start": 4,
                                        "component_stop": 8,
                                        "dimension": 4
                                    }
                                ],
                                "result_requires_propagated_value": false,
                                "result_requires_unpropagated_value": true,
                                "momentum_slots": {
                                    "left": 0,
                                    "right": 1,
                                    "result": 2
                                },
                                "coupling": [1.0, 1.0],
                                "color_weight": [1.0, 0.0],
                                "accumulation": "sum-into-result-current",
                                "lowering": {
                                    "kind": 9,
                                    "backend": "symbolica",
                                    "tensor_names": [],
                                    "expression_head": "quark_gluon_weyl_current",
                                    "full_tensor_network_ready": true,
                                    "description": "Weyl QCD antifermion-fermion vector current",
                                    "kernel": "fermion_pair_to_vector",
                                    "input_roles": ["antifermion", "fermion"],
                                    "output_role": "vector",
                                    "coupling_mode": "fixed"
                                },
                                "full_tensor_network_ready": true
                            }
                        ]
                    }
                ],
                "amplitude_stage": {
                    "stage_kind": "amplitude-roots",
                    "output_count": 1,
                    "roots": [
                        {
                            "output_index": 0,
                            "root_id": 0,
                            "kind": "direct-contraction",
                            "left_current_id": 0,
                            "right_current_id": 1,
                            "left_slot": {
                                "current_id": 0,
                                "component_start": 0,
                                "component_stop": 2,
                                "dimension": 2
                            },
                            "right_slot": {
                                "current_id": 1,
                                "component_start": 2,
                                "component_stop": 4,
                                "dimension": 2
                            },
                            "left_value_slot": {
                                "value_slot_id": 0,
                                "current_id": 0,
                                "variant": "source",
                                "component_start": 0,
                                "component_stop": 2,
                                "dimension": 2
                            },
                            "right_value_slot": {
                                "value_slot_id": 1,
                                "current_id": 1,
                                "variant": "source",
                                "component_start": 2,
                                "component_stop": 4,
                                "dimension": 2
                            },
                            "vertex_kind": null,
                            "vertex_particles": null,
                            "coupling": [1.0, 0.0],
                            "color_weight": [1.0, 0.0],
                            "contraction": "weyl",
                            "coherent_group_id": null,
                            "helicity_weight": 1.0
                        }
                    ]
                }
            }
        })
    }

    fn add_minimal_stage_evaluators(payload: &mut Value) {
        payload["compiled"]["kind"] = json!("generic-dag-stage-blueprint");
        payload["compiled"]["runtime_available"] = json!(true);
        payload["compiled"]["runtime_unavailable_message"] = json!(null);
        payload["compiled"]["stage_evaluators"] = json!({
            "kind": "generic-dag-stage-evaluator-artifacts",
            "runtime_available": true,
            "runtime_unavailable_message": null,
            "parameter_count": 20,
            "value_parameter_count": 8,
            "momentum_parameter_count": 12,
            "real_valued_inputs": [8,9,10,11,12,13,14,15,16,17,18,19],
            "parameter_layout": "global-value-momentum",
            "stage_count": 2,
            "stages": [
                {
                    "stage_index": 1,
                    "stage_kind": "current-combine",
                    "subset_size": 2,
                    "evaluator_label": "generic_stage_1_subset_2",
                    "parameter_layout": "global-value-momentum",
                    "output_length": 4,
                    "output_slots": [
                        {
                            "value_slot_id": 2,
                            "current_id": 2,
                            "variant": "unpropagated",
                            "component_start": 4,
                            "component_stop": 8,
                            "output_start": 0,
                            "output_stop": 4
                        }
                    ],
                    "input_value_slot_ids": [0, 1],
                    "output_value_slot_ids": [2],
                    "interaction_ids": [0],
                    "expression_ready": true,
                    "blockers": [],
                    "first_output_previews": ["p0"],
                    "evaluator": {
                        "kind": "jit-symbolica-evaluator",
                        "input_len": 20,
                        "output_len": 4,
                        "evaluator_state_path": "evaluators/stage_1.evaluator.bin"
                    }
                }
            ],
            "amplitude_stage": {
                "stage_index": 0,
                "stage_kind": "amplitude-roots",
                "subset_size": null,
                "evaluator_label": "generic_amplitude_stage",
                "parameter_layout": "global-value-momentum",
                "output_length": 1,
                "output_slots": [
                    {
                        "value_slot_id": -1,
                        "current_id": -1,
                        "variant": "amplitude-root",
                        "component_start": 0,
                        "component_stop": 1,
                        "output_start": 0,
                        "output_stop": 1
                    }
                ],
                "input_value_slot_ids": [0, 1],
                "output_value_slot_ids": [],
                "interaction_ids": [],
                "expression_ready": true,
                "blockers": [],
                "first_output_previews": ["amp0"],
                "evaluator": {
                    "kind": "jit-symbolica-evaluator",
                    "input_len": 20,
                    "output_len": 1,
                    "evaluator_state_path": "evaluators/amplitude.evaluator.bin"
                }
            }
        });
    }

    fn add_minimal_lc_topology_replay(payload: &mut Value) {
        payload["compiled"]["lc_topology_replay"] = json!({
            "enabled": true,
            "mode": "external-label-permutation",
            "replayed_sector_count": 2,
            "groups": [
                {
                    "representative_sector_id": 0,
                    "materialized_sector_id": 0,
                    "active_sector_ids": [0, 1],
                    "sector_permutations": [
                        {
                            "sector_id": 0,
                            "label_permutation": [
                                {"representative_label": 1, "sector_label": 1},
                                {"representative_label": 2, "sector_label": 2},
                                {"representative_label": 3, "sector_label": 3}
                            ]
                        },
                        {
                            "sector_id": 1,
                            "label_permutation": [
                                {"representative_label": 1, "sector_label": 2},
                                {"representative_label": 2, "sector_label": 1},
                                {"representative_label": 3, "sector_label": 3}
                            ]
                        }
                    ]
                }
            ]
        });
    }

    #[test]
    fn generic_schema_v2_validator_accepts_consistent_manifest() {
        let manifest: GenericProcessManifestV2 =
            serde_json::from_value(minimal_generic_manifest()).unwrap();

        validate_generic_schema_v2_manifest(&manifest).unwrap();
    }

    #[test]
    fn generic_schema_v2_validator_accepts_serialized_stage_evaluators() {
        let mut payload = minimal_generic_manifest();
        add_minimal_stage_evaluators(&mut payload);
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        validate_generic_schema_v2_manifest(&manifest).unwrap();
    }

    #[test]
    fn generic_schema_v2_manifest_parses_to_metadata_runtime() {
        let runtime = parse_generic_schema_v2_manifest(&minimal_generic_manifest())
            .unwrap()
            .expect("generic schema-v2 manifest should be recognized");

        assert_eq!(runtime.process, "d d~ > z");
        assert_eq!(runtime.key, "d_dbar_to_z");
        assert_eq!(runtime.color_accuracy, "lc");
        assert_eq!(runtime.external_pdg_order, [1, -1, 23]);
        assert_eq!(runtime.external_count, 3);
        assert_eq!(runtime.current_count, 3);
        assert_eq!(runtime.source_count, 2);
        assert_eq!(runtime.interaction_count, 1);
        assert_eq!(runtime.stage_count, 1);
        assert_eq!(runtime.stage_evaluator_count, 0);
        assert_eq!(runtime.amplitude_output_count, 1);
        assert!(runtime.lc_topology_reuse_available);
        assert_eq!(runtime.lc_topology_group_count, 1);
        assert_eq!(runtime.lc_topology_representative_sector_ids, [0]);
        assert!(!runtime.lc_topology_replay_enabled);
        assert_eq!(runtime.lc_topology_replay_sector_count, 0);
    }

    #[test]
    fn generic_schema_v2_manifest_parses_lc_topology_replay_metadata() {
        let mut payload = minimal_generic_manifest();
        add_minimal_lc_topology_replay(&mut payload);

        let runtime = parse_generic_schema_v2_manifest(&payload)
            .unwrap()
            .expect("generic schema-v2 manifest should be recognized");

        assert!(runtime.lc_topology_replay_enabled);
        assert_eq!(runtime.lc_topology_replay_sector_count, 2);
        assert_eq!(
            runtime.lc_topology_replay_mappings,
            vec![vec![(0, 0), (1, 1), (2, 2)], vec![(0, 1), (1, 0), (2, 2)],]
        );
    }

    #[test]
    fn generic_schema_v2_validator_rejects_bad_lc_topology_replay_mode() {
        let mut payload = minimal_generic_manifest();
        add_minimal_lc_topology_replay(&mut payload);
        payload["compiled"]["lc_topology_replay"]["mode"] = json!("family-shortcut");
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        let error = validate_generic_schema_v2_manifest(&manifest)
            .expect_err("bad topology replay mode should be rejected");

        assert!(
            error
                .to_string()
                .contains("unsupported LC topology replay mode")
        );
    }

    #[test]
    fn generic_schema_v2_manifest_preserves_stage_evaluator_metadata() {
        let mut payload = minimal_generic_manifest();
        add_minimal_stage_evaluators(&mut payload);

        let runtime = parse_generic_schema_v2_manifest(&payload)
            .unwrap()
            .expect("generic schema-v2 manifest should be recognized");

        assert_eq!(runtime.parameter_count, 20);
        assert_eq!(runtime.value_parameter_count, 8);
        assert_eq!(runtime.momentum_parameter_count, 12);
        assert_eq!(runtime.stage_evaluator_count, 2);
        assert_eq!(runtime.stage_evaluator_labels, ["generic_stage_1_subset_2"]);
        assert_eq!(
            runtime.amplitude_evaluator_label.as_deref(),
            Some("generic_amplitude_stage")
        );
    }

    #[test]
    fn generic_schema_v2_validator_rejects_bad_serialized_stage_evaluator_io() {
        let mut payload = minimal_generic_manifest();
        add_minimal_stage_evaluators(&mut payload);
        payload["compiled"]["stage_evaluators"]["stages"][0]["evaluator"]["output_len"] = json!(3);
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        let error = validate_generic_schema_v2_manifest(&manifest)
            .expect_err("bad serialized stage evaluator IO should be rejected");
        assert!(error.to_string().contains("inconsistent evaluator IO"));
    }

    #[test]
    fn generic_schema_v2_validator_rejects_bad_stage_count() {
        let mut payload = minimal_generic_manifest();
        payload["runtime_schema"]["stages"][0]["interaction_count"] = json!(2);
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        let error = validate_generic_schema_v2_manifest(&manifest)
            .expect_err("bad stage count should be rejected");
        assert!(
            error
                .to_string()
                .contains("generic stage 1 is inconsistent")
        );
    }

    #[test]
    fn generic_schema_v2_validator_accepts_unpropagated_auxiliary_inputs() {
        let mut payload = minimal_generic_manifest();
        payload["dag_summary"]["current_count"] = json!(4);
        payload["dag_summary"]["interaction_count"] = json!(2);
        payload["runtime_schema"]["parameter_layout"]["value_component_count"] = json!(10);
        payload["runtime_schema"]["current_storage"]["component_count"] = json!(10);
        payload["runtime_schema"]["value_storage"]["component_count"] = json!(10);
        payload["runtime_schema"]["current_storage"]["current_slots"]
            .as_array_mut()
            .unwrap()
            .push(json!({
                "current_id": 3,
                "component_start": 8,
                "component_stop": 10,
                "dimension": 2,
                "is_source": false,
                "particle_id": 21,
                "external_mask": 3,
                "external_labels": [1, 2],
                "helicity_ancestry": 3,
                "chirality": 0,
                "spin_state": 0,
                "flavour_flow": [-1, 1, 21],
                "charge_flow": 0,
                "color_state": {"accuracy": "lc", "sector_id": 0, "line_groups": [0], "basis_key": []},
                "momentum_mask": 3,
                "auxiliary_kind": null
            }));
        payload["runtime_schema"]["value_storage"]["value_slots"]
            .as_array_mut()
            .unwrap()
            .push(json!({
                "value_slot_id": 3,
                "current_id": 3,
                "variant": "unpropagated",
                "component_start": 8,
                "component_stop": 10,
                "dimension": 2,
                "current_component_start": 8,
                "current_component_stop": 10,
                "is_source": false,
                "applies_propagator": false,
                "particle_id": 21,
                "external_mask": 3,
                "external_labels": [1, 2],
                "momentum_mask": 3,
                "chirality": 0
            }));
        payload["runtime_schema"]["stages"]
            .as_array_mut()
            .unwrap()
            .push(json!({
                "stage_index": 2,
                "stage_kind": "current-combine",
                "subset_size": 3,
                "input_current_ids": [2, 0],
                "output_current_ids": [3],
                "input_value_slot_ids": [2, 0],
                "output_value_slot_ids": [3],
                "interaction_count": 1,
                "interactions": [
                    {
                        "interaction_id": 1,
                        "vertex_kind": 2,
                        "vertex_particles": [-21, 21, 21],
                        "left_current_id": 2,
                        "right_current_id": 0,
                        "result_current_id": 3,
                        "left_slot": {
                            "current_id": 2,
                            "component_start": 4,
                            "component_stop": 8,
                            "dimension": 4
                        },
                        "right_slot": {
                            "current_id": 0,
                            "component_start": 0,
                            "component_stop": 2,
                            "dimension": 2
                        },
                        "result_slot": {
                            "current_id": 3,
                            "component_start": 8,
                            "component_stop": 10,
                            "dimension": 2
                        },
                        "left_value_slot": {
                            "value_slot_id": 2,
                            "current_id": 2,
                            "variant": "unpropagated",
                            "component_start": 4,
                            "component_stop": 8,
                            "dimension": 4
                        },
                        "right_value_slot": {
                            "value_slot_id": 0,
                            "current_id": 0,
                            "variant": "source",
                            "component_start": 0,
                            "component_stop": 2,
                            "dimension": 2
                        },
                        "result_value_slots": [
                            {
                                "value_slot_id": 3,
                                "current_id": 3,
                                "variant": "unpropagated",
                                "component_start": 8,
                                "component_stop": 10,
                                "dimension": 2
                            }
                        ],
                        "result_requires_propagated_value": false,
                        "result_requires_unpropagated_value": true,
                        "momentum_slots": {
                            "left": 2,
                            "right": 0,
                            "result": 2
                        },
                        "coupling": [1.0, 0.0],
                        "color_weight": [1.0, 0.0],
                        "accumulation": "sum-into-result-current",
                        "lowering": {
                            "kind": 2,
                            "backend": "spenso",
                            "tensor_names": ["tensor_gluon_to_gluon"],
                            "expression_head": "tensor_gluon_to_gluon",
                            "full_tensor_network_ready": true,
                            "description": "auxiliary tensor and gluon to gluon current",
                            "kernel": "tensor_vector_to_vector",
                            "input_roles": ["antisymmetric_tensor", "vector"],
                            "output_role": "vector",
                            "coupling_mode": "fixed"
                        },
                        "full_tensor_network_ready": true
                    }
                ]
            }));
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        validate_generic_schema_v2_manifest(&manifest).unwrap();
    }
}

#[pymodule]
fn rusticol(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Runtime>()?;
    m.add_function(wrap_pyfunction!(build_profile, m)?)?;
    m.add_function(wrap_pyfunction!(build_target, m)?)?;
    Ok(())
}

#[pyfunction]
fn build_profile() -> &'static str {
    env!("RUSTICOL_BUILD_PROFILE")
}

#[pyfunction]
fn build_target() -> &'static str {
    env!("RUSTICOL_BUILD_TARGET")
}
