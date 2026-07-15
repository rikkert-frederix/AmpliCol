#![recursion_limit = "256"]

#[cfg(feature = "python")]
use numpy::IntoPyArray;
#[cfg(feature = "python")]
use pyo3::IntoPyObjectExt;
#[cfg(feature = "python")]
use pyo3::buffer::PyBuffer;
#[cfg(feature = "python")]
use pyo3::exceptions::{PyRuntimeError, PyValueError};
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::{PyAny, PyDict, PyList, PyTuple};
use serde::de::DeserializeOwned;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::Instant;
use symbolica::evaluate::JITCompiledEvaluator;
use symbolica::prelude::{
    BatchEvaluator, CompiledComplexEvaluator, Complex, DoubleFloat, EvaluationDomain,
    ExpressionEvaluator, Float, JITCompilationSettings, Rational, Real, RealLike,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RusticolError {
    message: String,
}

impl RusticolError {
    pub fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl std::fmt::Display for RusticolError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.message)
    }
}

impl std::error::Error for RusticolError {}

#[cfg(not(feature = "python"))]
type PyResult<T> = Result<T, RusticolError>;
#[cfg(not(feature = "python"))]
type PyErr = RusticolError;

#[cfg(not(feature = "python"))]
struct PyValueError;
#[cfg(not(feature = "python"))]
impl PyValueError {
    fn new_err(message: impl std::fmt::Display) -> RusticolError {
        RusticolError::new(message.to_string())
    }
}

#[cfg(not(feature = "python"))]
struct PyRuntimeError;
#[cfg(not(feature = "python"))]
impl PyRuntimeError {
    fn new_err(message: impl std::fmt::Display) -> RusticolError {
        RusticolError::new(message.to_string())
    }
}

const MAX_LC_TOPOLOGY_REPLAY_EXPANDED_POINTS: usize = 8192;
const LC_SECTOR_SELECTOR_PARAMETER: &str = "runtime.lc_sector_id";
const PROCESS_MANIFEST_READ_BUFFER_BYTES: usize = 4 * 1024 * 1024;

#[derive(Clone, Debug, Deserialize)]
struct ArtifactManifestHeader {
    #[serde(default)]
    schema_version: u32,
    #[serde(default)]
    kind: String,
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
    #[serde(default)]
    physics: Option<PhysicsManifestV2>,
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
    #[serde(default)]
    model_parameter_evaluator: Option<GenericModelParameterEvaluatorManifestV2>,
    stage_evaluators: Option<GenericStageEvaluatorArtifactsManifestV2>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericModelParameterEvaluatorManifestV2 {
    kind: String,
    input_parameter_indices: Vec<usize>,
    outputs: Vec<GenericDerivedParameterOutputManifestV2>,
    evaluator: EvaluatorManifest,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericDerivedParameterOutputManifestV2 {
    runtime_name: String,
    output_index: usize,
    real_parameter_index: usize,
    imag_parameter_index: usize,
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
    #[serde(default = "default_lc_topology_replay_weight")]
    weight: f64,
    #[serde(default)]
    label_permutation: Vec<LcTopologyReplayLabelPermutationManifestV2>,
}

#[derive(Clone, Debug, Default, Deserialize)]
struct LcTopologyReplayLabelPermutationManifestV2 {
    representative_label: usize,
    sector_label: usize,
}

fn default_lc_topology_replay_weight() -> f64 {
    1.0
}

#[derive(Clone, Debug, Deserialize)]
struct GenericStageEvaluatorArtifactsManifestV2 {
    kind: String,
    runtime_available: bool,
    runtime_unavailable_message: Option<String>,
    parameter_count: usize,
    value_parameter_count: usize,
    momentum_parameter_count: usize,
    #[serde(default)]
    model_parameter_count: usize,
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
    model_parameter_count: usize,
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
    model_parameters: Vec<GenericRuntimeModelParameterManifestV2>,
    #[serde(default)]
    normalization: Option<GenericRuntimeNormalizationManifestV2>,
    #[serde(default)]
    physics: Option<PhysicsManifestV2>,
    parameter_layout: GenericParameterLayoutManifestV2,
    current_storage: GenericCurrentStorageManifestV2,
    value_storage: GenericValueStorageManifestV2,
    source_fill: GenericSourceFillManifestV2,
    momentum_slots: Vec<GenericMomentumSlotManifestV2>,
    stages: Vec<GenericStageManifestV2>,
    amplitude_stage: GenericAmplitudeStageManifestV2,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsManifestV2 {
    schema_version: u32,
    kind: String,
    process: String,
    process_key: String,
    color_accuracy: String,
    external_particles: Vec<PhysicsExternalParticleManifestV2>,
    helicities: Vec<PhysicsHelicityManifestV2>,
    color_components: Vec<PhysicsColorComponentManifestV2>,
    #[serde(default)]
    model_parameters: Vec<PhysicsModelParameterManifestV2>,
    coverage: PhysicsCoverageManifestV2,
    selectors: PhysicsSelectorsManifestV2,
    reduction: PhysicsReductionManifestV2,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsExternalParticleManifestV2 {
    label: usize,
    index: usize,
    side: String,
    role: String,
    particle: String,
    outgoing_particle: String,
    pdg: i32,
    outgoing_pdg: i32,
    particle_class: String,
    momentum_slot: usize,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsHelicityManifestV2 {
    id: String,
    index: usize,
    helicities: Vec<i32>,
    representative_id: String,
    computed: bool,
    structural_zero: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsColorComponentManifestV2 {
    id: String,
    index: usize,
    kind: String,
    #[serde(default)]
    word: Vec<usize>,
    representative_id: String,
    computed: bool,
    #[serde(default)]
    internal_sector_id: Option<i64>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsModelParameterManifestV2 {
    name: String,
    kind: String,
    parameter_index: usize,
    #[serde(default)]
    default: f64,
    #[serde(flatten)]
    extra: BTreeMap<String, Value>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsCoverageManifestV2 {
    helicities: String,
    color: String,
    color_kind: String,
    #[serde(default)]
    structural_zero_helicity_count: usize,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsSelectorsManifestV2 {
    helicity: bool,
    color_flow: bool,
    contracted_color: bool,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsReductionManifestV2 {
    kind: String,
    groups: Vec<PhysicsReductionGroupManifestV2>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
struct PhysicsReductionGroupManifestV2 {
    group_id: i64,
    representative_helicity_id: String,
    physical_helicity_ids: Vec<String>,
    representative_color_id: String,
    physical_color_ids: Vec<String>,
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
struct GenericRuntimeModelParameterManifestV2 {
    name: String,
    kind: String,
    parameter_index: usize,
    #[serde(default)]
    default: f64,
    #[serde(default)]
    pdg: Option<i32>,
    #[serde(default)]
    runtime_name: Option<String>,
    #[serde(default)]
    complex_component: Option<String>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericRuntimeParticleManifestV2 {
    pdg: i32,
    #[serde(default)]
    mass: f64,
    #[serde(default)]
    mass_parameter: Option<String>,
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
    #[serde(default)]
    qcd_coupling_power: usize,
    #[serde(default)]
    electroweak_coupling_power: usize,
}

fn default_one_f64() -> f64 {
    1.0
}

#[derive(Clone, Debug, Deserialize)]
struct GenericParameterLayoutManifestV2 {
    source_component_parameter_count: usize,
    momentum_parameter_count: usize,
    #[serde(default)]
    model_parameter_count: usize,
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
    #[serde(default)]
    metadata_compacted: bool,
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
    #[serde(default)]
    external_labels: Vec<usize>,
    #[serde(default)]
    helicity_ancestry: Value,
    chirality: i32,
    #[serde(default)]
    spin_state: Value,
    #[serde(default)]
    flavour_flow: Vec<i32>,
    #[serde(default)]
    charge_flow: i32,
    #[serde(default)]
    color_state: Value,
    momentum_mask: u64,
    auxiliary_kind: Option<String>,
}

#[derive(Clone, Debug, Deserialize)]
struct GenericValueStorageManifestV2 {
    component_count: usize,
    number_type: String,
    #[serde(default)]
    metadata_compacted: bool,
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
    #[serde(default)]
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
    #[serde(default)]
    wavefunction_kind: String,
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
    #[serde(default)]
    interactions_compacted: bool,
    #[serde(default)]
    interaction_ids: Vec<usize>,
    #[serde(default)]
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
    #[serde(default)]
    color_sector_id: Option<i64>,
    contraction: String,
    coherent_group_id: Option<Value>,
    helicity_weight: f64,
    #[serde(default)]
    all_sector_weight: Option<f64>,
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
            physics: None,
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
    let supported_process_set =
        manifest.schema_version == 2 && manifest.kind == "pyamplicol-generic-dag-process-set";
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
                physics: entry.physics.clone(),
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

fn read_json_manifest<T: DeserializeOwned>(manifest_path: &Path) -> Result<T, String> {
    let file = fs::File::open(manifest_path).map_err(|err| {
        format!(
            "could not read process manifest {}: {err}",
            manifest_path.display()
        )
    })?;
    serde_json::from_reader(BufReader::with_capacity(
        PROCESS_MANIFEST_READ_BUFFER_BYTES,
        file,
    ))
    .map_err(|err| {
        format!(
            "could not parse process manifest {}: {err}",
            manifest_path.display()
        )
    })
}

fn load_generic_schema_v2_manifest(
    manifest: GenericProcessManifestV2,
    root: &Path,
) -> PyResult<GenericRuntimeV2> {
    validate_generic_schema_v2_manifest(&manifest)?;
    GenericRuntimeV2::load_from_manifest(manifest, root)
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
    let (mappings, _weights) = build_lc_topology_replay_mappings(Some(replay))?;
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
        != layout.source_component_parameter_count
            + layout.momentum_parameter_count
            + layout.model_parameter_count
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
    if schema.model_parameters.len() != layout.model_parameter_count {
        return Err(PyValueError::new_err(
            "generic runtime schema model-parameter count is inconsistent",
        ));
    }
    let mut seen_model_parameters = BTreeSet::new();
    let mut seen_model_parameter_names = BTreeSet::new();
    for parameter in &schema.model_parameters {
        if parameter.name.is_empty()
            || parameter.parameter_index >= layout.model_parameter_count
            || !parameter.default.is_finite()
            || !seen_model_parameters.insert(parameter.parameter_index)
            || !seen_model_parameter_names.insert(parameter.name.clone())
            || parameter
                .runtime_name
                .as_ref()
                .is_some_and(|name| name.is_empty())
            || match (
                &parameter.runtime_name,
                parameter.complex_component.as_deref(),
            ) {
                (Some(_), Some("real" | "imag")) | (None, None) => false,
                _ => true,
            }
        {
            return Err(PyValueError::new_err(
                "generic runtime schema contains invalid model-parameter metadata",
            ));
        }
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
        if slot.chirality.abs() > 1 {
            return Err(PyValueError::new_err(format!(
                "generic current slot {index} has invalid quantum-flow metadata"
            )));
        }
        if !storage.metadata_compacted {
            if slot.external_labels.is_empty()
                || !positive_json_integer(&slot.helicity_ancestry)
                || slot.color_state.is_null()
            {
                return Err(PyValueError::new_err(format!(
                    "generic current slot {index} is missing current-index metadata"
                )));
            }
            if slot.flavour_flow.is_empty() {
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
            || (!storage.metadata_compacted && slot.external_labels != current.external_labels)
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
        let max_helicity = if source.wavefunction_kind == "spin2" {
            2
        } else {
            1
        };
        if source.physical_pdg == 0
            || source.outgoing_pdg != source.particle_id
            || source.chirality.abs() > 1
            || source.source_helicity.abs() > max_helicity
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
    let input_value_variants =
        build_input_value_variants(&schema.current_storage, &schema.value_storage)?;
    let mut seen_interactions = 0usize;
    let mut seen_interaction_ids = BTreeSet::new();
    for (stage_offset, stage) in schema.stages.iter().enumerate() {
        if stage.stage_index != stage_offset + 1
            || stage.stage_kind != "current-combine"
            || stage.subset_size < 2
        {
            return Err(PyValueError::new_err(format!(
                "generic stage {} is inconsistent",
                stage.stage_index
            )));
        }
        if stage.interactions_compacted {
            if !stage.interactions.is_empty()
                || stage.interaction_ids.len() != stage.interaction_count
            {
                return Err(PyValueError::new_err(format!(
                    "generic compact stage {} has inconsistent interaction metadata",
                    stage.stage_index
                )));
            }
            for interaction_id in &stage.interaction_ids {
                if *interaction_id >= manifest.dag_summary.interaction_count {
                    return Err(PyValueError::new_err(format!(
                        "generic compact stage {} references invalid interaction {interaction_id}",
                        stage.stage_index
                    )));
                }
                if !seen_interaction_ids.insert(*interaction_id) {
                    return Err(PyValueError::new_err(format!(
                        "generic compact stage {} repeats interaction {interaction_id}",
                        stage.stage_index
                    )));
                }
            }
        } else if stage.interaction_count != stage.interactions.len()
            || !stage.interaction_ids.is_empty()
        {
            return Err(PyValueError::new_err(format!(
                "generic stage {} is inconsistent",
                stage.stage_index
            )));
        }
        let mut input_current_membership = vec![false; current_count];
        let mut output_current_membership = vec![false; current_count];
        for (ids, membership) in [
            (&stage.input_current_ids, &mut input_current_membership),
            (&stage.output_current_ids, &mut output_current_membership),
        ] {
            for id in ids {
                if *id >= current_count {
                    return Err(PyValueError::new_err(format!(
                        "generic stage {} references invalid current {id}",
                        stage.stage_index
                    )));
                }
                membership[*id] = true;
            }
        }
        let mut input_value_membership = vec![false; value_count];
        let mut output_value_membership = vec![false; value_count];
        for (ids, membership) in [
            (&stage.input_value_slot_ids, &mut input_value_membership),
            (&stage.output_value_slot_ids, &mut output_value_membership),
        ] {
            for value_id in ids {
                if *value_id >= value_count {
                    return Err(PyValueError::new_err(format!(
                        "generic stage {} references invalid value slot {value_id}",
                        stage.stage_index
                    )));
                }
                membership[*value_id] = true;
            }
        }
        for interaction in &stage.interactions {
            if !seen_interaction_ids.insert(interaction.interaction_id) {
                return Err(PyValueError::new_err(format!(
                    "generic stage {} repeats interaction {}",
                    stage.stage_index, interaction.interaction_id
                )));
            }
            validate_generic_interaction(
                interaction,
                &schema.current_storage,
                &schema.value_storage,
                &input_value_variants,
                momentum_count,
            )?;
            if !input_current_membership[interaction.left_current_id]
                || !input_current_membership[interaction.right_current_id]
                || !output_current_membership[interaction.result_current_id]
            {
                return Err(PyValueError::new_err(format!(
                    "generic interaction {} is not listed in its stage inputs/outputs",
                    interaction.interaction_id
                )));
            }
            if !input_value_membership[interaction.left_value_slot.value_slot_id]
                || !input_value_membership[interaction.right_value_slot.value_slot_id]
                || interaction
                    .result_value_slots
                    .iter()
                    .any(|slot| !output_value_membership[slot.value_slot_id])
            {
                return Err(PyValueError::new_err(format!(
                    "generic interaction {} value slots are not listed in its stage inputs/outputs",
                    interaction.interaction_id
                )));
            }
        }
        seen_interactions += stage.interaction_count;
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
    input_value_variants: &[Option<GenericInputValueVariant>],
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
            input_value_variants,
            interaction.left_current_id,
        )?),
        value_storage,
    )?;
    validate_value_slot_ref(
        &interaction.right_value_slot,
        interaction.right_current_id,
        Some(input_value_variant(
            input_value_variants,
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
        + schema.parameter_layout.momentum_parameter_count
        + schema.parameter_layout.model_parameter_count;
    let expected_real_inputs = (schema.parameter_layout.value_component_count
        ..expected_parameter_count)
        .collect::<Vec<_>>();
    let header_is_global = stage_evaluators.parameter_layout == "global-value-momentum"
        && stage_evaluators.parameter_count == expected_parameter_count
        && stage_evaluators.value_parameter_count == schema.parameter_layout.value_component_count
        && stage_evaluators.momentum_parameter_count
            == schema.parameter_layout.momentum_parameter_count
        && stage_evaluators.model_parameter_count == schema.parameter_layout.model_parameter_count
        && stage_evaluators.real_valued_inputs == expected_real_inputs;
    let header_is_stage_local = stage_evaluators.parameter_layout == "stage-local-value-momentum"
        && stage_evaluators.parameter_count == 0
        && stage_evaluators.value_parameter_count == 0
        && stage_evaluators.momentum_parameter_count == 0
        && stage_evaluators.model_parameter_count == 0
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
            || stage.value_parameter_count
                + stage.momentum_parameter_count
                + stage.model_parameter_count
                != stage.parameter_count
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
        let expected_interactions = if runtime_stage.interactions_compacted {
            runtime_stage.interaction_ids.clone()
        } else {
            runtime_stage
                .interactions
                .iter()
                .map(|interaction| interaction.interaction_id)
                .collect::<Vec<_>>()
        };
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
        if component.kind != "value"
            && component.kind != "momentum"
            && component.kind != "model_parameter"
        {
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
        if component.real_valued
            && component.kind != "momentum"
            && component.kind != "model_parameter"
        {
            return Err(PyValueError::new_err(format!(
                "generic serialized stage evaluator {} marks a complex local input as real",
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GenericInputValueVariant {
    Source,
    Propagated,
    Unpropagated,
}

impl GenericInputValueVariant {
    fn as_str(self) -> &'static str {
        match self {
            Self::Source => "source",
            Self::Propagated => "propagated",
            Self::Unpropagated => "unpropagated",
        }
    }
}

fn build_input_value_variants(
    current_storage: &GenericCurrentStorageManifestV2,
    value_storage: &GenericValueStorageManifestV2,
) -> PyResult<Vec<Option<GenericInputValueVariant>>> {
    let mut variants = current_storage
        .current_slots
        .iter()
        .map(|current| {
            current
                .is_source
                .then_some(GenericInputValueVariant::Source)
        })
        .collect::<Vec<_>>();
    for slot in &value_storage.value_slots {
        let current = current_storage
            .current_slots
            .get(slot.current_id)
            .ok_or_else(|| {
                PyValueError::new_err(format!(
                    "generic value slot references missing current {}",
                    slot.current_id
                ))
            })?;
        if current.is_source {
            continue;
        }
        match slot.variant.as_str() {
            "propagated" => {
                variants[slot.current_id] = Some(GenericInputValueVariant::Propagated);
            }
            "unpropagated"
                if variants[slot.current_id] != Some(GenericInputValueVariant::Propagated) =>
            {
                variants[slot.current_id] = Some(GenericInputValueVariant::Unpropagated);
            }
            _ => {}
        }
    }
    Ok(variants)
}

fn input_value_variant(
    input_value_variants: &[Option<GenericInputValueVariant>],
    current_id: usize,
) -> PyResult<&'static str> {
    let variant = input_value_variants.get(current_id).ok_or_else(|| {
        PyValueError::new_err(format!(
            "generic input value references missing current {current_id}"
        ))
    })?;
    variant
        .map(GenericInputValueVariant::as_str)
        .ok_or_else(|| {
            PyValueError::new_err(format!(
                "generic input value references current {current_id} without an input value slot"
            ))
        })
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
    chunk_scratch_native2: Vec<Complex<wide::f64x2>>,
}

enum F64Evaluator {
    Compiled(CompiledComplexEvaluator),
    Jit(JITCompiledEvaluator<Complex<f64>>),
    JitNative2(JITCompiledEvaluator<Complex<wide::f64x2>>),
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

struct RawSumGroup {
    id: i64,
    indices: Vec<usize>,
    weight: f64,
    all_sector_weight: f64,
    sector_ids: Vec<i64>,
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

#[derive(Clone, Debug, Default)]
struct RuntimeProfile {
    source_fill_s: f64,
    momentum_setup_s: f64,
    stage_input_pack_s: f64,
    stage_evaluator_call_s: f64,
    stage_evaluator_s: f64,
    output_assign_s: f64,
    amplitude_input_pack_s: f64,
    amplitude_evaluator_call_s: f64,
    amplitude_evaluator_s: f64,
    reduction_s: f64,
    total_s: f64,
    stage_input_pack_by_stage_s: Vec<f64>,
    stage_evaluator_call_by_stage_s: Vec<f64>,
    stage_output_assign_by_stage_s: Vec<f64>,
}

impl RuntimeProfile {
    fn add_sector(&mut self, sector: &RuntimeProfile) {
        self.source_fill_s += sector.source_fill_s;
        self.momentum_setup_s += sector.momentum_setup_s;
        self.stage_input_pack_s += sector.stage_input_pack_s;
        self.stage_evaluator_call_s += sector.stage_evaluator_call_s;
        self.stage_evaluator_s += sector.stage_evaluator_s;
        self.output_assign_s += sector.output_assign_s;
        self.amplitude_input_pack_s += sector.amplitude_input_pack_s;
        self.amplitude_evaluator_call_s += sector.amplitude_evaluator_call_s;
        self.amplitude_evaluator_s += sector.amplitude_evaluator_s;
        self.reduction_s += sector.reduction_s;
        add_profile_vector(
            &mut self.stage_input_pack_by_stage_s,
            &sector.stage_input_pack_by_stage_s,
        );
        add_profile_vector(
            &mut self.stage_evaluator_call_by_stage_s,
            &sector.stage_evaluator_call_by_stage_s,
        );
        add_profile_vector(
            &mut self.stage_output_assign_by_stage_s,
            &sector.stage_output_assign_by_stage_s,
        );
    }
}

fn add_profile_vector(target: &mut Vec<f64>, source: &[f64]) {
    if target.len() < source.len() {
        target.resize(source.len(), 0.0);
    }
    for (index, value) in source.iter().enumerate() {
        target[index] += value;
    }
}

#[cfg(feature = "python")]
fn parse_optional_color_sector_ids(
    value: Option<&Bound<'_, PyAny>>,
) -> PyResult<Option<BTreeSet<i64>>> {
    let Some(value) = value else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    parse_required_color_sector_ids(value).map(Some)
}

#[cfg(feature = "python")]
fn parse_required_color_sector_ids(value: &Bound<'_, PyAny>) -> PyResult<BTreeSet<i64>> {
    if let Ok(sector_id) = value.extract::<i64>() {
        return Ok(BTreeSet::from([sector_id]));
    }
    if let Ok(sector_ids) = value.extract::<Vec<i64>>() {
        if sector_ids.is_empty() {
            return Err(PyValueError::new_err(
                "color_sector_ids must contain at least one sector id",
            ));
        }
        return Ok(sector_ids.into_iter().collect());
    }
    Err(PyValueError::new_err(
        "color_sector_ids must be an integer or a sequence of integers",
    ))
}

#[derive(Clone, Copy, Debug, Default)]
struct MemorySnapshot {
    current_rss_bytes: Option<u64>,
    peak_rss_bytes: Option<u64>,
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
    model_parameter_count: usize,
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
    lc_topology_replay_weights: Vec<f64>,
    runtime_unavailable_message: Option<String>,
    sources: Vec<GenericSourceRecordManifestV2>,
    momentum_slots: Vec<GenericMomentumSlotManifestV2>,
    external_is_initial: Vec<bool>,
    particle_masses: BTreeMap<i32, f64>,
    particle_mass_parameter_names: BTreeMap<i32, String>,
    normalization_factor: f64,
    normalization_color_factor: f64,
    normalization_average_factor: f64,
    normalization_identical_factor: f64,
    normalization_qcd_coupling_power: usize,
    normalization_electroweak_coupling_power: usize,
    model_parameters: Vec<GenericRuntimeModelParameterManifestV2>,
    model_parameter_name_to_index: BTreeMap<String, usize>,
    model_parameter_runtime_slots: BTreeMap<String, GenericRuntimeParameterSlots>,
    model_parameter_values_f64: Vec<f64>,
    model_parameter_evaluator: Option<GenericModelParameterEvaluatorRuntimeV2>,
    physics: Option<PhysicsRuntimeV2>,
    stages: Option<Vec<GenericStageRuntimeV2>>,
    amplitude_stage: Option<GenericAmplitudeRuntimeV2>,
    state_scratch_f64: Vec<Complex<f64>>,
    values_scratch_f64: Vec<f64>,
}

#[derive(Clone)]
struct PhysicsRuntimeV2 {
    manifest: PhysicsManifestV2,
    helicity_index_by_id: BTreeMap<String, usize>,
    color_index_by_id: BTreeMap<String, usize>,
    reduction_by_group_id: BTreeMap<i64, PhysicsReductionGroupManifestV2>,
}

#[derive(Clone, Debug)]
struct ResolvedValues<T> {
    values: Vec<T>,
    point_count: usize,
    helicity_indices: Vec<usize>,
    color_indices: Vec<usize>,
}

#[derive(Clone, Debug, Serialize)]
pub struct NativeRuntimeMetadata {
    pub abi_version: u32,
    pub schema_version: u32,
    pub process: String,
    pub process_key: String,
    pub representative_process: String,
    pub representative_process_key: String,
    pub crossing_alias_of: Option<String>,
    pub color_accuracy: String,
    pub external_pdg_order: Vec<i32>,
    pub external_count: usize,
    pub current_count: usize,
    pub source_count: usize,
    pub interaction_count: usize,
    pub stage_count: usize,
    pub amplitude_output_count: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct NativeResolvedEvaluation {
    /// Row-major storage with layout `[point][helicity][color]`.
    pub values: Vec<f64>,
    pub point_count: usize,
    pub helicity_ids: Vec<String>,
    pub color_ids: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct NativeExternalParticle {
    pub label: usize,
    pub index: usize,
    pub side: String,
    pub role: String,
    pub particle: String,
    pub outgoing_particle: String,
    pub pdg: i32,
    pub outgoing_pdg: i32,
    pub particle_class: String,
    pub momentum_slot: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct NativeHelicityConfiguration {
    pub id: String,
    pub index: usize,
    pub helicities: Vec<i32>,
    pub representative_id: String,
    pub computed: bool,
    pub structural_zero: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct NativeColorComponent {
    pub id: String,
    pub index: usize,
    pub kind: String,
    pub word: Vec<usize>,
    pub representative_id: String,
    pub computed: bool,
}

#[derive(Clone, Debug, Serialize)]
pub struct NativeModelParameter {
    pub name: String,
    pub kind: String,
    pub parameter_index: usize,
    pub default: f64,
}

impl NativeResolvedEvaluation {
    pub fn shape(&self) -> (usize, usize, usize) {
        (
            self.point_count,
            self.helicity_ids.len(),
            self.color_ids.len(),
        )
    }

    pub fn totals(&self) -> Vec<f64> {
        let component_count = self.helicity_ids.len() * self.color_ids.len();
        self.values
            .chunks(component_count)
            .map(|point| point.iter().sum())
            .collect()
    }
}

/// Python-independent schema-v2 process runtime.
///
/// The input momentum layout is `[point][external particle][E, px, py, pz]`.
/// Instances are mutable and must not be called concurrently; independent
/// instances can be used from separate threads.
pub struct NativeRuntime {
    root: PathBuf,
    runtime: GenericRuntimeV2,
    process: String,
    process_key: String,
    input_crossing_map: Option<Vec<InputCrossingMapEntry>>,
    crossing_alias_of: Option<String>,
    physics: Option<PhysicsManifestV2>,
    warnings_muted: bool,
    warned_kinds: BTreeSet<String>,
    pending_warnings: Vec<String>,
}

impl NativeRuntime {
    pub const ABI_VERSION: u32 = 1;

    pub fn load(
        process_dir: impl AsRef<Path>,
        process_key: Option<&str>,
        model_parameters_path: Option<&Path>,
    ) -> Result<Self, RusticolError> {
        let requested = process_dir.as_ref();
        let requested_root = requested.canonicalize().map_err(|error| {
            RusticolError::new(format!(
                "could not resolve process directory {}: {error}",
                requested.display()
            ))
        })?;
        let selection = resolve_process_root(&requested_root, process_key)
            .map_err(|error| RusticolError::new(error.to_string()))?;
        let manifest_path = selection.root.join("process_manifest.json");
        let manifest = read_json_manifest::<GenericProcessManifestV2>(&manifest_path)
            .map_err(RusticolError::new)?;
        if manifest.schema_version != 2 || manifest.kind != "pyamplicol-generic-dag-process" {
            return Err(RusticolError::new(format!(
                "Rusticol native APIs support only schema-v2 generic DAG artifacts; got kind {:?} schema {}",
                manifest.kind, manifest.schema_version
            )));
        }
        let representative_process = manifest.process.clone();
        let representative_key = manifest.key.clone();
        let mut runtime = load_generic_schema_v2_manifest(manifest, &selection.root)
            .map_err(|error| RusticolError::new(error.to_string()))?;
        if let Some(metadata) = selection.physics.clone() {
            runtime.physics = Some(
                PhysicsRuntimeV2::new(metadata)
                    .map_err(|error| RusticolError::new(error.to_string()))?,
            );
        }
        if let Some(path) = model_parameters_path {
            runtime
                .apply_model_parameter_json_path(path)
                .map_err(|error| RusticolError::new(error.to_string()))?;
        }
        let process = selection
            .selected_process
            .clone()
            .unwrap_or_else(|| representative_process.clone());
        let process_key = selection
            .selected_key
            .clone()
            .unwrap_or_else(|| representative_key.clone());
        let mut physics = runtime.physics.as_ref().map(|value| value.manifest.clone());
        if let Some(metadata) = physics.as_mut() {
            metadata.process = process.clone();
            metadata.process_key = process_key.clone();
        }
        Ok(Self {
            root: selection.root,
            runtime,
            process,
            process_key,
            input_crossing_map: selection.input_crossing_map,
            crossing_alias_of: selection.crossing_alias_of,
            physics,
            warnings_muted: false,
            warned_kinds: BTreeSet::new(),
            pending_warnings: Vec::new(),
        })
    }

    pub fn metadata(&self) -> NativeRuntimeMetadata {
        NativeRuntimeMetadata {
            abi_version: Self::ABI_VERSION,
            schema_version: 2,
            process: self.process.clone(),
            process_key: self.process_key.clone(),
            representative_process: self.runtime.process.clone(),
            representative_process_key: self.runtime.key.clone(),
            crossing_alias_of: self.crossing_alias_of.clone(),
            color_accuracy: self.runtime.color_accuracy.clone(),
            external_pdg_order: self.runtime.external_pdg_order.clone(),
            external_count: self.runtime.external_count,
            current_count: self.runtime.current_count,
            source_count: self.runtime.source_count,
            interaction_count: self.runtime.interaction_count,
            stage_count: self.runtime.stage_count,
            amplitude_output_count: self.runtime.amplitude_output_count,
        }
    }

    pub fn metadata_json(&self) -> Result<String, RusticolError> {
        serde_json::to_string(&self.metadata()).map_err(|error| {
            RusticolError::new(format!("could not serialize runtime metadata: {error}"))
        })
    }

    pub fn physics_json(&self) -> Result<String, RusticolError> {
        let physics = self.physics.as_ref().ok_or_else(|| {
            RusticolError::new(
                "this older schema-v2 artifact has no resolved physics metadata; regenerate it with pyAmpliCol",
            )
        })?;
        serde_json::to_string(physics).map_err(|error| {
            RusticolError::new(format!("could not serialize physics metadata: {error}"))
        })
    }

    pub fn external_count(&self) -> usize {
        self.runtime.external_count
    }

    pub fn external_particles(&self) -> Result<Vec<NativeExternalParticle>, RusticolError> {
        Ok(self
            .physics
            .as_ref()
            .ok_or_else(|| {
                RusticolError::new("external-particle metadata requires artifact regeneration")
            })?
            .external_particles
            .iter()
            .map(|item| NativeExternalParticle {
                label: item.label,
                index: item.index,
                side: item.side.clone(),
                role: item.role.clone(),
                particle: item.particle.clone(),
                outgoing_particle: item.outgoing_particle.clone(),
                pdg: item.pdg,
                outgoing_pdg: item.outgoing_pdg,
                particle_class: item.particle_class.clone(),
                momentum_slot: item.momentum_slot,
            })
            .collect())
    }

    pub fn helicities(&self) -> Result<Vec<NativeHelicityConfiguration>, RusticolError> {
        Ok(self
            .physics
            .as_ref()
            .ok_or_else(|| {
                RusticolError::new("resolved helicities require regenerated physics metadata")
            })?
            .helicities
            .iter()
            .map(|item| NativeHelicityConfiguration {
                id: item.id.clone(),
                index: item.index,
                helicities: item.helicities.clone(),
                representative_id: item.representative_id.clone(),
                computed: item.computed,
                structural_zero: item.structural_zero,
            })
            .collect())
    }

    pub fn color_components(&self) -> Result<Vec<NativeColorComponent>, RusticolError> {
        Ok(self
            .physics
            .as_ref()
            .ok_or_else(|| {
                RusticolError::new("resolved color components require regenerated physics metadata")
            })?
            .color_components
            .iter()
            .map(|item| NativeColorComponent {
                id: item.id.clone(),
                index: item.index,
                kind: item.kind.clone(),
                word: item.word.clone(),
                representative_id: item.representative_id.clone(),
                computed: item.computed,
            })
            .collect())
    }

    pub fn model_parameters(&self) -> Result<Vec<NativeModelParameter>, RusticolError> {
        Ok(self
            .physics
            .as_ref()
            .ok_or_else(|| {
                RusticolError::new("model-parameter metadata requires artifact regeneration")
            })?
            .model_parameters
            .iter()
            .map(|item| NativeModelParameter {
                name: item.name.clone(),
                kind: item.kind.clone(),
                parameter_index: item.parameter_index,
                default: item.default,
            })
            .collect())
    }

    pub fn helicity_ids(&self) -> Result<Vec<String>, RusticolError> {
        Ok(self
            .physics
            .as_ref()
            .ok_or_else(|| {
                RusticolError::new("resolved helicities require regenerated physics metadata")
            })?
            .helicities
            .iter()
            .map(|item| item.id.clone())
            .collect())
    }

    pub fn color_ids(&self) -> Result<Vec<String>, RusticolError> {
        Ok(self
            .physics
            .as_ref()
            .ok_or_else(|| {
                RusticolError::new("resolved color components require regenerated physics metadata")
            })?
            .color_components
            .iter()
            .map(|item| item.id.clone())
            .collect())
    }

    pub fn resolved_shape(
        &self,
        helicity_ids: Option<&[String]>,
        color_ids: Option<&[String]>,
    ) -> Result<(usize, usize), RusticolError> {
        if color_ids.is_some() && self.runtime.color_accuracy != "lc" {
            return Err(RusticolError::new(
                "LC color-flow selection is unavailable for NLC/full artifacts; their resolved color axis is contracted",
            ));
        }
        let selected_helicities = selector_set(helicity_ids, "helicity")?;
        let selected_colors = selector_set(color_ids, "color component")?;
        let physics = self.runtime.physics.as_ref().ok_or_else(|| {
            RusticolError::new(
                "resolved evaluation is unavailable for this older schema-v2 artifact; regenerate it with pyAmpliCol",
            )
        })?;
        let helicity_count = physics
            .selected_helicity_indices(selected_helicities.as_ref())
            .map_err(|error| RusticolError::new(error.to_string()))?
            .len();
        let color_count = physics
            .selected_color_indices(selected_colors.as_ref())
            .map_err(|error| RusticolError::new(error.to_string()))?
            .len();
        Ok((helicity_count, color_count))
    }

    pub fn evaluate_f64(
        &mut self,
        momenta: &[f64],
        point_count: usize,
    ) -> Result<Vec<f64>, RusticolError> {
        let batch = self.prepare_f64_batch(momenta, point_count)?;
        self.runtime
            .run_f64(&batch)
            .map(|(values, _profile)| values)
            .map_err(|error| RusticolError::new(error.to_string()))
    }

    pub fn evaluate_resolved_f64(
        &mut self,
        momenta: &[f64],
        point_count: usize,
        helicity_ids: Option<&[String]>,
        color_ids: Option<&[String]>,
    ) -> Result<NativeResolvedEvaluation, RusticolError> {
        if color_ids.is_some() && self.runtime.color_accuracy != "lc" {
            return Err(RusticolError::new(
                "LC color-flow selection is unavailable for NLC/full artifacts; their resolved color axis is contracted",
            ));
        }
        self.record_resolved_warnings(helicity_ids, color_ids)?;
        let selected_helicities = selector_set(helicity_ids, "helicity")?;
        let selected_colors = selector_set(color_ids, "color component")?;
        let batch = self.prepare_f64_batch(momenta, point_count)?;
        let physics = self.physics.clone().ok_or_else(|| {
            RusticolError::new(
                "resolved evaluation is unavailable for this older schema-v2 artifact; regenerate it with pyAmpliCol",
            )
        })?;
        let (resolved, _profile) = self
            .runtime
            .run_resolved_f64(
                &batch,
                selected_helicities.as_ref(),
                selected_colors.as_ref(),
            )
            .map_err(|error| RusticolError::new(error.to_string()))?;
        let helicity_ids = resolved
            .helicity_indices
            .iter()
            .map(|index| physics.helicities[*index].id.clone())
            .collect();
        let color_ids = resolved
            .color_indices
            .iter()
            .map(|index| physics.color_components[*index].id.clone())
            .collect();
        Ok(NativeResolvedEvaluation {
            values: resolved.values,
            point_count: resolved.point_count,
            helicity_ids,
            color_ids,
        })
    }

    pub fn set_model_parameters(
        &mut self,
        values: &BTreeMap<String, (f64, f64)>,
    ) -> Result<(), RusticolError> {
        self.runtime
            .apply_model_parameter_overrides(values)
            .map_err(|error| RusticolError::new(error.to_string()))
    }

    pub fn set_model_parameter(
        &mut self,
        name: &str,
        real: f64,
        imaginary: f64,
    ) -> Result<(), RusticolError> {
        self.set_model_parameters(&BTreeMap::from([(name.to_string(), (real, imaginary))]))
    }

    pub fn set_model_parameters_json(&mut self, path: &Path) -> Result<(), RusticolError> {
        self.runtime
            .apply_model_parameter_json_path(path)
            .map_err(|error| RusticolError::new(error.to_string()))
    }

    pub fn mute_warnings(&mut self) {
        self.warnings_muted = true;
    }

    pub fn unmute_warnings(&mut self) {
        self.warnings_muted = false;
    }

    pub fn take_warnings(&mut self) -> Vec<String> {
        std::mem::take(&mut self.pending_warnings)
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    fn prepare_f64_batch(
        &self,
        momenta: &[f64],
        point_count: usize,
    ) -> Result<Vec<Vec<[f64; 4]>>, RusticolError> {
        if point_count == 0 {
            return Err(RusticolError::new("point_count must be positive"));
        }
        let values_per_point = self
            .runtime
            .external_count
            .checked_mul(4)
            .ok_or_else(|| RusticolError::new("momentum shape overflow"))?;
        let expected = point_count
            .checked_mul(values_per_point)
            .ok_or_else(|| RusticolError::new("momentum shape overflow"))?;
        if momenta.len() != expected {
            return Err(RusticolError::new(format!(
                "momenta contain {} values, expected {expected} for shape ({point_count}, {}, 4)",
                momenta.len(),
                self.runtime.external_count
            )));
        }
        let mut batch = Vec::with_capacity(point_count);
        for point_values in momenta.chunks_exact(values_per_point) {
            let point = point_values
                .chunks_exact(4)
                .map(|components| [components[0], components[1], components[2], components[3]])
                .collect();
            batch.push(point);
        }
        apply_input_crossing_map(
            batch,
            self.runtime.external_count,
            self.input_crossing_map.as_deref(),
        )
        .map_err(|error| RusticolError::new(error.to_string()))
    }

    fn record_resolved_warnings(
        &mut self,
        helicity_ids: Option<&[String]>,
        color_ids: Option<&[String]>,
    ) -> Result<(), RusticolError> {
        if self.warnings_muted {
            return Ok(());
        }
        let physics = self.runtime.physics.as_ref().ok_or_else(|| {
            RusticolError::new("resolved evaluation requires regenerated physics metadata")
        })?;
        let mut warnings = Vec::new();
        if physics.manifest.coverage.helicities != "complete" {
            warnings.push((
                "incomplete-helicity-coverage",
                "resolved evaluation contains only the helicities represented by this artifact",
            ));
        }
        if physics.manifest.coverage.color != "complete" {
            warnings.push((
                "incomplete-color-coverage",
                "resolved evaluation contains only the color components represented by this artifact",
            ));
        }
        let reduction_only_helicity = helicity_ids.is_some_and(|ids| {
            ids.iter().any(|id| {
                physics
                    .helicity_index_by_id
                    .get(id)
                    .and_then(|index| physics.manifest.helicities.get(*index))
                    .is_some_and(|item| !item.computed)
            })
        });
        let reduction_only_color = color_ids.is_some_and(|ids| {
            ids.iter().any(|id| {
                physics
                    .color_index_by_id
                    .get(id)
                    .and_then(|index| physics.manifest.color_components.get(*index))
                    .is_some_and(|item| !item.computed)
            })
        });
        if reduction_only_helicity || reduction_only_color {
            warnings.push((
                "reduction-only-selection",
                "the selected resolved component reuses an exact symmetry representative",
            ));
        }
        for (kind, message) in warnings {
            if self.warned_kinds.insert(kind.to_string()) {
                self.pending_warnings.push(message.to_string());
            }
        }
        Ok(())
    }
}

fn selector_set(
    ids: Option<&[String]>,
    kind: &str,
) -> Result<Option<BTreeSet<String>>, RusticolError> {
    let Some(ids) = ids else {
        return Ok(None);
    };
    if ids.is_empty() {
        return Err(RusticolError::new(format!(
            "resolved {kind} selection must not be empty"
        )));
    }
    let selected = ids.iter().cloned().collect::<BTreeSet<_>>();
    if selected.len() != ids.len() {
        return Err(RusticolError::new(format!(
            "resolved {kind} selection contains duplicate ids"
        )));
    }
    Ok(Some(selected))
}

impl PhysicsRuntimeV2 {
    fn new(manifest: PhysicsManifestV2) -> PyResult<Self> {
        if manifest.schema_version != 1 || manifest.kind != "pyamplicol-resolved-physics" {
            return Err(PyValueError::new_err(format!(
                "unsupported resolved physics metadata kind {:?} schema {}",
                manifest.kind, manifest.schema_version
            )));
        }
        let mut helicity_index_by_id = BTreeMap::new();
        for (index, helicity) in manifest.helicities.iter().enumerate() {
            if helicity.index != index
                || helicity.helicities.len() != manifest.external_particles.len()
            {
                return Err(PyValueError::new_err(format!(
                    "invalid resolved helicity metadata at index {index}"
                )));
            }
            if helicity_index_by_id
                .insert(helicity.id.clone(), index)
                .is_some()
            {
                return Err(PyValueError::new_err(format!(
                    "duplicate public helicity id {:?}",
                    helicity.id
                )));
            }
        }
        let mut color_index_by_id = BTreeMap::new();
        for (index, color) in manifest.color_components.iter().enumerate() {
            if color.index != index {
                return Err(PyValueError::new_err(format!(
                    "invalid resolved color metadata at index {index}"
                )));
            }
            if color_index_by_id.insert(color.id.clone(), index).is_some() {
                return Err(PyValueError::new_err(format!(
                    "duplicate public color id {:?}",
                    color.id
                )));
            }
        }
        let mut reduction_by_group_id = BTreeMap::new();
        for group in &manifest.reduction.groups {
            if group.physical_helicity_ids.is_empty() || group.physical_color_ids.is_empty() {
                return Err(PyValueError::new_err(format!(
                    "resolved reduction group {} has an empty physical expansion",
                    group.group_id
                )));
            }
            for id in &group.physical_helicity_ids {
                if !helicity_index_by_id.contains_key(id) {
                    return Err(PyValueError::new_err(format!(
                        "resolved reduction group {} references unknown helicity {id:?}",
                        group.group_id
                    )));
                }
            }
            for id in &group.physical_color_ids {
                if !color_index_by_id.contains_key(id) {
                    return Err(PyValueError::new_err(format!(
                        "resolved reduction group {} references unknown color component {id:?}",
                        group.group_id
                    )));
                }
            }
            if reduction_by_group_id
                .insert(group.group_id, group.clone())
                .is_some()
            {
                return Err(PyValueError::new_err(format!(
                    "duplicate resolved reduction group {}",
                    group.group_id
                )));
            }
        }
        Ok(Self {
            manifest,
            helicity_index_by_id,
            color_index_by_id,
            reduction_by_group_id,
        })
    }

    fn selected_helicity_indices(&self, ids: Option<&BTreeSet<String>>) -> PyResult<Vec<usize>> {
        self.select_indices(ids, &self.helicity_index_by_id, "helicity")
    }

    fn selected_color_indices(&self, ids: Option<&BTreeSet<String>>) -> PyResult<Vec<usize>> {
        self.select_indices(ids, &self.color_index_by_id, "color component")
    }

    fn select_indices(
        &self,
        ids: Option<&BTreeSet<String>>,
        available: &BTreeMap<String, usize>,
        kind: &str,
    ) -> PyResult<Vec<usize>> {
        let Some(ids) = ids else {
            return Ok((0..available.len()).collect());
        };
        if ids.is_empty() {
            return Err(PyValueError::new_err(format!(
                "resolved {kind} selection must not be empty"
            )));
        }
        let mut indices = Vec::with_capacity(ids.len());
        for id in ids {
            let index = available.get(id).ok_or_else(|| {
                PyValueError::new_err(format!("unknown resolved {kind} id {id:?}"))
            })?;
            indices.push(*index);
        }
        indices.sort_unstable();
        Ok(indices)
    }
}

#[derive(Clone, Copy, Debug)]
struct GenericRuntimeParameterSlots {
    real: usize,
    imaginary: Option<usize>,
}

struct GenericModelParameterEvaluatorRuntimeV2 {
    input_parameter_indices: Vec<usize>,
    outputs: Vec<GenericDerivedParameterOutputManifestV2>,
    evaluator: EvaluatorGroup,
}

struct GenericStageRuntimeV2 {
    outputs: Vec<(usize, usize)>,
    output_spans: Vec<(usize, usize, usize)>,
    input_components: Option<Vec<usize>>,
    input_spans: Vec<(usize, usize, usize)>,
    parameter_scratch_f64: Vec<Complex<f64>>,
    output_scratch_f64: Vec<Complex<f64>>,
    parameter_scratch_native2: Vec<Complex<wide::f64x2>>,
    output_scratch_native2: Vec<Complex<wide::f64x2>>,
    evaluator: EvaluatorGroup,
}

struct GenericAmplitudeRuntimeV2 {
    output_length: usize,
    raw_sum_weights: Vec<f64>,
    raw_sum_all_sector_weights: Vec<f64>,
    raw_sum_color_sector_ids: Vec<Option<i64>>,
    raw_sum_groups: Vec<RawSumGroup>,
    has_coherent_groups: bool,
    color_contraction: Option<ColorContractionRuntime>,
    input_components: Option<Vec<usize>>,
    input_spans: Vec<(usize, usize, usize)>,
    parameter_scratch_f64: Vec<Complex<f64>>,
    output_scratch_f64: Vec<Complex<f64>>,
    parameter_scratch_native2: Vec<Complex<wide::f64x2>>,
    output_scratch_native2: Vec<Complex<wide::f64x2>>,
    evaluator: EvaluatorGroup,
}

fn build_lc_topology_replay_mappings(
    replay: Option<&LcTopologyReplayManifestV2>,
) -> PyResult<(Vec<Vec<(usize, usize)>>, Vec<f64>)> {
    let Some(replay) = replay else {
        return Ok((Vec::new(), Vec::new()));
    };
    if !replay.enabled {
        return Ok((Vec::new(), Vec::new()));
    }
    if replay.mode != "external-label-permutation" {
        return Err(PyValueError::new_err(format!(
            "unsupported LC topology replay mode {:?}",
            replay.mode
        )));
    }
    let mut mappings = Vec::new();
    let mut weights = Vec::new();
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
            if !permutation.weight.is_finite() || permutation.weight <= 0.0 {
                return Err(PyValueError::new_err(
                    "LC topology replay sector permutation weights must be positive finite numbers",
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
            weights.push(permutation.weight);
        }
    }
    if mappings.len() != replay.replayed_sector_count {
        return Err(PyValueError::new_err(format!(
            "LC topology replay declares {} sectors but contains {} permutations",
            replay.replayed_sector_count,
            mappings.len()
        )));
    }
    Ok((mappings, weights))
}

fn parse_complex_parameter_overrides(
    text: &str,
    path: &Path,
) -> PyResult<BTreeMap<String, (f64, f64)>> {
    let value = serde_json::from_str::<Value>(text).map_err(|err| {
        PyValueError::new_err(format!(
            "could not parse model-parameter JSON {}: {err}",
            path.display()
        ))
    })?;
    let object = value.as_object().ok_or_else(|| {
        PyValueError::new_err(format!(
            "model-parameter JSON {} must contain an object",
            path.display()
        ))
    })?;
    let mut overrides = BTreeMap::new();
    for (name, value) in object {
        let components = value.as_array().ok_or_else(|| {
            PyValueError::new_err(format!(
                "model-parameter JSON entry {name:?} must be [real, imaginary]"
            ))
        })?;
        if components.len() != 2 {
            return Err(PyValueError::new_err(format!(
                "model-parameter JSON entry {name:?} must have exactly two components"
            )));
        }
        let real = components[0].as_f64().ok_or_else(|| {
            PyValueError::new_err(format!(
                "model-parameter JSON entry {name:?} has a non-numeric real component"
            ))
        })?;
        let imaginary = components[1].as_f64().ok_or_else(|| {
            PyValueError::new_err(format!(
                "model-parameter JSON entry {name:?} has a non-numeric imaginary component"
            ))
        })?;
        if !real.is_finite() || !imaginary.is_finite() {
            return Err(PyValueError::new_err(format!(
                "model-parameter JSON entry {name:?} must contain finite values"
            )));
        }
        overrides.insert(name.clone(), (real, imaginary));
    }
    Ok(overrides)
}

fn build_runtime_parameter_slots(
    parameters: &[GenericRuntimeModelParameterManifestV2],
) -> PyResult<BTreeMap<String, GenericRuntimeParameterSlots>> {
    let mut direct = BTreeMap::new();
    let mut complex_components: BTreeMap<String, (Option<usize>, Option<usize>)> = BTreeMap::new();
    for parameter in parameters {
        if parameter.kind == "derived_parameter_component" {
            continue;
        }
        if let Some(runtime_name) = &parameter.runtime_name {
            let slots = complex_components
                .entry(runtime_name.clone())
                .or_insert((None, None));
            match parameter.complex_component.as_deref() {
                Some("real") if slots.0.replace(parameter.parameter_index).is_none() => {}
                Some("imag") if slots.1.replace(parameter.parameter_index).is_none() => {}
                Some(component) => {
                    return Err(PyValueError::new_err(format!(
                        "runtime model parameter {runtime_name:?} has duplicate or invalid component {component:?}"
                    )));
                }
                None => {
                    return Err(PyValueError::new_err(format!(
                        "runtime model parameter {runtime_name:?} is missing component metadata"
                    )));
                }
            }
        } else if direct
            .insert(
                parameter.name.clone(),
                GenericRuntimeParameterSlots {
                    real: parameter.parameter_index,
                    imaginary: None,
                },
            )
            .is_some()
        {
            return Err(PyValueError::new_err(format!(
                "duplicate runtime model parameter name {:?}",
                parameter.name
            )));
        }
    }
    for (name, (real, imaginary)) in complex_components {
        let (Some(real), Some(imaginary)) = (real, imaginary) else {
            return Err(PyValueError::new_err(format!(
                "runtime model parameter {name:?} requires real and imaginary slots"
            )));
        };
        if direct
            .insert(
                name.clone(),
                GenericRuntimeParameterSlots {
                    real,
                    imaginary: Some(imaginary),
                },
            )
            .is_some()
        {
            return Err(PyValueError::new_err(format!(
                "duplicate runtime model parameter name {name:?}"
            )));
        }
    }
    Ok(direct)
}

impl GenericRuntimeV2 {
    fn from_manifest(manifest: GenericProcessManifestV2) -> PyResult<Self> {
        let stage_evaluators = manifest.compiled.stage_evaluators.as_ref();
        let topology_reuse = manifest.lc_topology_reuse.as_ref();
        let topology_replay = manifest.compiled.lc_topology_replay.as_ref();
        let (topology_replay_mappings, topology_replay_weights) =
            build_lc_topology_replay_mappings(topology_replay)?;
        let physics = manifest
            .runtime_schema
            .physics
            .clone()
            .map(PhysicsRuntimeV2::new)
            .transpose()?;
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
        let particle_mass_parameter_names = manifest
            .runtime_schema
            .model
            .as_ref()
            .map(|model| {
                model
                    .particles
                    .iter()
                    .filter_map(|particle| {
                        particle
                            .mass_parameter
                            .as_ref()
                            .map(|name| (particle.pdg, name.clone()))
                    })
                    .collect::<BTreeMap<_, _>>()
            })
            .unwrap_or_default();
        let mut model_parameters = manifest.runtime_schema.model_parameters.clone();
        model_parameters.sort_by_key(|parameter| parameter.parameter_index);
        let model_parameter_values_f64 = model_parameters
            .iter()
            .map(|parameter| parameter.default)
            .collect::<Vec<_>>();
        let model_parameter_name_to_index = model_parameters
            .iter()
            .map(|parameter| (parameter.name.clone(), parameter.parameter_index))
            .collect::<BTreeMap<_, _>>();
        let model_parameter_runtime_slots = build_runtime_parameter_slots(&model_parameters)?;
        let color_factor_in_contraction = manifest
            .runtime_schema
            .amplitude_stage
            .color_contraction
            .as_ref()
            .map(|contraction| contraction.supported && contraction.includes_color_factor)
            .unwrap_or(false);
        let (
            normalization_factor,
            normalization_color_factor,
            normalization_average_factor,
            normalization_identical_factor,
            normalization_qcd_coupling_power,
            normalization_electroweak_coupling_power,
        ) = manifest
            .runtime_schema
            .normalization
            .as_ref()
            .map(|normalization| {
                let color_factor = if color_factor_in_contraction {
                    1.0
                } else {
                    normalization.color_factor
                };
                (
                    color_factor * normalization.global_coupling_factor
                        / (normalization.average_factor * normalization.identical_factor),
                    color_factor,
                    normalization.average_factor,
                    normalization.identical_factor,
                    normalization.qcd_coupling_power,
                    normalization.electroweak_coupling_power,
                )
            })
            .unwrap_or((1.0, 1.0, 1.0, 1.0, 0, 0));
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
                    .momentum_parameter_count
                + manifest
                    .runtime_schema
                    .parameter_layout
                    .model_parameter_count,
            value_parameter_count: manifest
                .runtime_schema
                .parameter_layout
                .value_component_count,
            momentum_parameter_count: manifest
                .runtime_schema
                .parameter_layout
                .momentum_parameter_count,
            model_parameter_count: manifest
                .runtime_schema
                .parameter_layout
                .model_parameter_count,
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
            lc_topology_replay_weights: topology_replay_weights,
            runtime_unavailable_message: manifest.compiled.runtime_unavailable_message,
            sources: manifest.runtime_schema.source_fill.sources,
            momentum_slots: manifest.runtime_schema.momentum_slots,
            external_is_initial,
            particle_masses,
            particle_mass_parameter_names,
            normalization_factor,
            normalization_color_factor,
            normalization_average_factor,
            normalization_identical_factor,
            normalization_qcd_coupling_power,
            normalization_electroweak_coupling_power,
            model_parameters,
            model_parameter_name_to_index,
            model_parameter_runtime_slots,
            model_parameter_values_f64,
            model_parameter_evaluator: None,
            physics,
            stages: None,
            amplitude_stage: None,
            state_scratch_f64: Vec::new(),
            values_scratch_f64: Vec::new(),
        })
    }

    fn load_from_manifest(manifest: GenericProcessManifestV2, root: &Path) -> PyResult<Self> {
        let stage_evaluators = manifest.compiled.stage_evaluators.clone();
        let model_parameter_evaluator = manifest.compiled.model_parameter_evaluator.clone();
        let amplitude_stage_manifest = manifest.runtime_schema.amplitude_stage.clone();
        let mut runtime = Self::from_manifest(manifest)?;
        if let Some(manifest) = model_parameter_evaluator {
            if manifest.kind != "generic-model-parameter-evaluator" {
                return Err(PyValueError::new_err(format!(
                    "unsupported model-parameter evaluator kind {:?}",
                    manifest.kind
                )));
            }
            runtime.model_parameter_evaluator = Some(GenericModelParameterEvaluatorRuntimeV2 {
                input_parameter_indices: manifest.input_parameter_indices,
                outputs: manifest.outputs,
                evaluator: EvaluatorGroup::load(&manifest.evaluator, root)?,
            });
            runtime.refresh_derived_model_parameters()?;
        }
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

    fn apply_model_parameter_json_path(&mut self, path: &Path) -> PyResult<()> {
        let text = fs::read_to_string(path).map_err(|err| {
            PyValueError::new_err(format!(
                "could not read model-parameter JSON {}: {err}",
                path.display()
            ))
        })?;
        let overrides = parse_complex_parameter_overrides(&text, path)?;
        self.apply_model_parameter_overrides(&overrides)
    }

    fn apply_model_parameter_overrides(
        &mut self,
        overrides: &BTreeMap<String, (f64, f64)>,
    ) -> PyResult<()> {
        let mut proposed = self.model_parameter_values_f64.clone();
        for (name, (real, imaginary)) in overrides {
            let Some(slots) = self.model_parameter_runtime_slots.get(name).copied() else {
                return Err(PyValueError::new_err(format!(
                    "model-parameter override {name:?} is not used by process {}",
                    self.process
                )));
            };
            if slots.real >= self.model_parameter_values_f64.len()
                || !real.is_finite()
                || !imaginary.is_finite()
            {
                return Err(PyValueError::new_err(format!(
                    "model-parameter override {name:?} has invalid value [{real}, {imaginary}]",
                )));
            }
            proposed[slots.real] = *real;
            if let Some(index) = slots.imaginary {
                if index >= self.model_parameter_values_f64.len() {
                    return Err(PyValueError::new_err(format!(
                        "model-parameter override {name:?} has an invalid imaginary slot"
                    )));
                }
                proposed[index] = *imaginary;
            } else if *imaginary != 0.0 {
                return Err(PyValueError::new_err(format!(
                    "real model parameter {name:?} cannot receive a nonzero imaginary component"
                )));
            }
        }
        let previous_values = std::mem::replace(&mut self.model_parameter_values_f64, proposed);
        let previous_masses = self.particle_masses.clone();
        let previous_normalization = self.normalization_factor;
        if let Err(error) = self.refresh_derived_model_parameters() {
            self.model_parameter_values_f64 = previous_values;
            self.particle_masses = previous_masses;
            self.normalization_factor = previous_normalization;
            return Err(error);
        }
        self.refresh_particle_mass_parameters();
        self.refresh_normalization_factor();
        Ok(())
    }

    fn refresh_derived_model_parameters(&mut self) -> PyResult<()> {
        let Some(runtime) = self.model_parameter_evaluator.as_mut() else {
            return Ok(());
        };
        let parameters = runtime
            .input_parameter_indices
            .iter()
            .map(|index| {
                self.model_parameter_values_f64
                    .get(*index)
                    .copied()
                    .map(|value| c64(value, 0.0))
                    .ok_or_else(|| {
                        PyValueError::new_err(format!(
                            "model-parameter evaluator input index {index} is out of range"
                        ))
                    })
            })
            .collect::<PyResult<Vec<_>>>()?;
        let evaluated = runtime.evaluator.evaluate_batch(1, &parameters)?;
        for output in &runtime.outputs {
            let value = evaluated.get(output.output_index).ok_or_else(|| {
                PyValueError::new_err(format!(
                    "model-parameter evaluator output {} for {:?} is absent",
                    output.output_index, output.runtime_name
                ))
            })?;
            let Some(real) = self
                .model_parameter_values_f64
                .get_mut(output.real_parameter_index)
            else {
                return Err(PyValueError::new_err(format!(
                    "derived model-parameter real slot {} is out of range",
                    output.real_parameter_index
                )));
            };
            *real = value.re;
            let Some(imaginary) = self
                .model_parameter_values_f64
                .get_mut(output.imag_parameter_index)
            else {
                return Err(PyValueError::new_err(format!(
                    "derived model-parameter imaginary slot {} is out of range",
                    output.imag_parameter_index
                )));
            };
            *imaginary = value.im;
        }
        Ok(())
    }

    fn refresh_particle_mass_parameters(&mut self) {
        for parameter in &self.model_parameters {
            if parameter.kind == "particle_mass" {
                if let Some(pdg) = parameter.pdg {
                    if let Some(value) = self
                        .model_parameter_values_f64
                        .get(parameter.parameter_index)
                    {
                        self.particle_masses.insert(pdg, *value);
                    }
                }
            }
        }
        for (pdg, name) in &self.particle_mass_parameter_names {
            let Some(slots) = self.model_parameter_runtime_slots.get(name) else {
                continue;
            };
            if let Some(value) = self.model_parameter_values_f64.get(slots.real) {
                self.particle_masses.insert(*pdg, *value);
            }
        }
    }

    fn refresh_normalization_factor(&mut self) {
        let alpha_s = self
            .model_parameter_name_to_index
            .get("normalization.alpha_s_me_check")
            .and_then(|index| self.model_parameter_values_f64.get(*index))
            .copied();
        let alpha_ew = self
            .model_parameter_name_to_index
            .get("normalization.alpha_ew")
            .and_then(|index| self.model_parameter_values_f64.get(*index))
            .copied();
        let mut global_coupling_factor = 1.0;
        if self.normalization_qcd_coupling_power > 0 {
            let Some(alpha_s) = alpha_s else {
                return;
            };
            global_coupling_factor *= (4.0 * std::f64::consts::PI * alpha_s)
                .powi(self.normalization_qcd_coupling_power as i32);
        }
        if self.normalization_electroweak_coupling_power > 0 {
            let Some(alpha_ew) = alpha_ew else {
                return;
            };
            global_coupling_factor *= (2.0 * 4.0 * std::f64::consts::PI * alpha_ew)
                .powi(self.normalization_electroweak_coupling_power as i32);
        }
        self.normalization_factor = self.normalization_color_factor * global_coupling_factor
            / (self.normalization_average_factor * self.normalization_identical_factor);
    }

    fn set_lc_sector_selector(&mut self, sector_id: Option<i64>) -> Option<f64> {
        let index = *self
            .model_parameter_name_to_index
            .get(LC_SECTOR_SELECTOR_PARAMETER)?;
        let previous = *self.model_parameter_values_f64.get(index)?;
        if let Some(value) = self.model_parameter_values_f64.get_mut(index) {
            *value = sector_id.map(|id| id as f64).unwrap_or(-1.0);
        }
        Some(previous)
    }

    fn restore_lc_sector_selector(&mut self, previous: Option<f64>) {
        let Some(previous) = previous else {
            return;
        };
        let Some(index) = self
            .model_parameter_name_to_index
            .get(LC_SECTOR_SELECTOR_PARAMETER)
            .copied()
        else {
            return;
        };
        if let Some(value) = self.model_parameter_values_f64.get_mut(index) {
            *value = previous;
        }
    }

    fn run_f64(&mut self, batch: &[Vec<[f64; 4]>]) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        self.run_f64_selected(batch, None)
    }

    fn run_resolved_f64(
        &mut self,
        batch: &[Vec<[f64; 4]>],
        selected_helicity_ids: Option<&BTreeSet<String>>,
        selected_color_ids: Option<&BTreeSet<String>>,
    ) -> PyResult<(ResolvedValues<f64>, RuntimeProfile)> {
        if self.lc_topology_replay_enabled {
            return Err(PyValueError::new_err(
                "resolved evaluation for LC topology-replay artifacts requires regeneration with resolved replay metadata",
            ));
        }
        let physics = self.physics.clone().ok_or_else(|| {
            PyValueError::new_err(
                "resolved evaluation is unavailable for this older schema-v2 artifact; regenerate it with pyAmpliCol",
            )
        })?;
        let (_summed, profile) = self.run_f64_materialized(batch)?;
        let resolved = self
            .amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .reduce_scratch_f64_resolved(
                batch.len(),
                &physics,
                self.normalization_factor,
                selected_helicity_ids,
                selected_color_ids,
            )?;
        Ok((resolved, profile))
    }

    fn run_f64_selected(
        &mut self,
        batch: &[Vec<[f64; 4]>],
        selected_color_sector_ids: Option<&BTreeSet<i64>>,
    ) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        if self.lc_topology_replay_enabled {
            if selected_color_sector_ids.is_some() {
                return Err(PyValueError::new_err(
                    "LC color-sector runtime selection is not available for topology-replay artifacts",
                ));
            }
            return self.run_f64_with_lc_topology_replay(batch);
        }
        let Some(selected) = selected_color_sector_ids else {
            let previous = self.set_lc_sector_selector(None);
            let result = self.run_f64_materialized_selected(batch, None);
            self.restore_lc_sector_selector(previous);
            return result;
        };
        if selected.is_empty() {
            return Err(PyValueError::new_err(
                "LC color-sector runtime selection requires at least one sector id",
            ));
        }
        if !self
            .model_parameter_name_to_index
            .contains_key(LC_SECTOR_SELECTOR_PARAMETER)
        {
            return self.run_f64_materialized_selected(batch, Some(selected));
        }
        let total_start = Instant::now();
        let n_points = batch.len();
        let mut values = vec![0.0; n_points];
        let mut profile = RuntimeProfile::default();
        let previous = self.set_lc_sector_selector(None);
        for sector_id in selected {
            self.set_lc_sector_selector(Some(*sector_id));
            let mut singleton = BTreeSet::new();
            singleton.insert(*sector_id);
            let (sector_values, sector_profile) =
                self.run_f64_materialized_selected(batch, Some(&singleton))?;
            for (value, sector_value) in values.iter_mut().zip(sector_values) {
                *value += sector_value;
            }
            profile.add_sector(&sector_profile);
        }
        self.restore_lc_sector_selector(previous);
        profile.total_s = total_start.elapsed().as_secs_f64();
        return Ok((values, profile));
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
        let weights = self.lc_topology_replay_weights.clone();
        let mappings_per_chunk = replay_mappings_per_expanded_batch(n_points);
        for chunk_start in (0..mappings.len()).step_by(mappings_per_chunk) {
            let chunk_end = usize::min(chunk_start + mappings_per_chunk, mappings.len());
            let mapping_chunk = &mappings[chunk_start..chunk_end];
            let weight_chunk = &weights[chunk_start..chunk_end];
            let expanded_batch =
                apply_lc_topology_label_permutations(batch, self.external_count, mapping_chunk)?;
            let (expanded_values, sector_profile) = self.run_f64_materialized(&expanded_batch)?;
            for mapping_index in 0..mapping_chunk.len() {
                let weight = weight_chunk[mapping_index];
                let offset = mapping_index * n_points;
                for point_index in 0..n_points {
                    values[point_index] += weight * expanded_values[offset + point_index];
                }
            }
            profile.add_sector(&sector_profile);
        }
        profile.total_s = total_start.elapsed().as_secs_f64();
        Ok((values, profile))
    }

    fn run_f64_materialized(
        &mut self,
        batch: &[Vec<[f64; 4]>],
    ) -> PyResult<(Vec<f64>, RuntimeProfile)> {
        self.run_f64_materialized_selected(batch, None)
    }

    fn run_f64_materialized_selected(
        &mut self,
        batch: &[Vec<[f64; 4]>],
        selected_color_sector_ids: Option<&BTreeSet<i64>>,
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
        let model_parameter_start = self.value_parameter_count + self.momentum_parameter_count;
        let model_parameter_values = &self.model_parameter_values_f64;

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

        let model_parameter_start_time = Instant::now();
        for row in 0..n_points {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_model_parameters_row(
                model_parameter_start,
                model_parameter_values,
                row_state,
            )?;
        }
        let model_parameter_setup_s = model_parameter_start_time.elapsed().as_secs_f64();

        let mut stage_input_pack_s = 0.0;
        let mut stage_evaluator_call_s = 0.0;
        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        let mut stage_input_pack_by_stage_s = Vec::new();
        let mut stage_evaluator_call_by_stage_s = Vec::new();
        let mut stage_output_assign_by_stage_s = Vec::new();
        for stage in self.stages.as_mut().expect("generic stages checked") {
            let (pack_s, eval_s, assign_s) = stage.evaluate_f64_into_state(
                n_points,
                self.parameter_count,
                state.as_mut_slice(),
            )?;
            stage_input_pack_s += pack_s;
            stage_evaluator_call_s += eval_s;
            stage_evaluator_s += pack_s + eval_s;
            output_assign_s += assign_s;
            stage_input_pack_by_stage_s.push(pack_s);
            stage_evaluator_call_by_stage_s.push(eval_s);
            stage_output_assign_by_stage_s.push(assign_s);
        }

        let (amplitude_input_pack_s, amplitude_evaluator_call_s) = self
            .amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .evaluate_f64_into_scratch(n_points, state.as_slice())?;
        let amplitude_evaluator_s = amplitude_input_pack_s + amplitude_evaluator_call_s;

        let reduction_start = Instant::now();
        self.amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .reduce_scratch_f64_into_selected(
                n_points,
                &mut self.values_scratch_f64,
                selected_color_sector_ids,
            )?;
        for value in &mut self.values_scratch_f64 {
            *value *= self.normalization_factor;
        }
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        Ok((
            self.values_scratch_f64.clone(),
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s: momentum_setup_s + model_parameter_setup_s,
                stage_input_pack_s,
                stage_evaluator_call_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_input_pack_s,
                amplitude_evaluator_call_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
                stage_input_pack_by_stage_s,
                stage_evaluator_call_by_stage_s,
                stage_output_assign_by_stage_s,
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
        let weights = self.lc_topology_replay_weights.clone();
        let mappings_per_chunk = replay_mappings_per_expanded_batch(n_points);
        for chunk_start in (0..mappings.len()).step_by(mappings_per_chunk) {
            let chunk_end = usize::min(chunk_start + mappings_per_chunk, mappings.len());
            let mapping_chunk = &mappings[chunk_start..chunk_end];
            let weight_chunk = &weights[chunk_start..chunk_end];
            let expanded_batch = apply_lc_topology_label_permutations_generic(
                batch,
                self.external_count,
                mapping_chunk,
            )?;
            let (expanded_values, sector_profile) =
                self.run_generic_materialized(&expanded_batch, binary_precision)?;
            for mapping_index in 0..mapping_chunk.len() {
                let weight = T::from(weight_chunk[mapping_index]);
                let offset = mapping_index * n_points;
                for point_index in 0..n_points {
                    values[point_index] +=
                        weight.clone() * expanded_values[offset + point_index].clone();
                }
            }
            profile.add_sector(&sector_profile);
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
        let model_parameter_start = self.value_parameter_count + self.momentum_parameter_count;
        let model_parameter_values = &self.model_parameter_values_f64;

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

        let model_parameter_start_time = Instant::now();
        for row in 0..n_points {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_model_parameters_row_generic(
                model_parameter_start,
                model_parameter_values,
                row_state,
            )?;
        }
        let model_parameter_setup_s = model_parameter_start_time.elapsed().as_secs_f64();

        let mut stage_input_pack_s = 0.0;
        let mut stage_evaluator_call_s = 0.0;
        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        let mut stage_input_pack_by_stage_s = Vec::new();
        let mut stage_evaluator_call_by_stage_s = Vec::new();
        let mut stage_output_assign_by_stage_s = Vec::new();
        for stage in self.stages.as_mut().expect("generic stages checked") {
            let (pack_s, eval_s, assign_s) = stage.evaluate_generic_into_state(
                n_points,
                self.parameter_count,
                state.as_mut_slice(),
                binary_precision,
            )?;
            stage_input_pack_s += pack_s;
            stage_evaluator_call_s += eval_s;
            stage_evaluator_s += pack_s + eval_s;
            output_assign_s += assign_s;
            stage_input_pack_by_stage_s.push(pack_s);
            stage_evaluator_call_by_stage_s.push(eval_s);
            stage_output_assign_by_stage_s.push(assign_s);
        }

        let (raw_sums, amplitude_input_pack_s, amplitude_evaluator_call_s) = self
            .amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .evaluate_raw_sums_generic(n_points, state.as_slice(), binary_precision)?;
        let amplitude_evaluator_s = amplitude_input_pack_s + amplitude_evaluator_call_s;

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
                momentum_setup_s: momentum_setup_s + model_parameter_setup_s,
                stage_input_pack_s,
                stage_evaluator_call_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_input_pack_s,
                amplitude_evaluator_call_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
                stage_input_pack_by_stage_s,
                stage_evaluator_call_by_stage_s,
                stage_output_assign_by_stage_s,
            },
        ))
    }

    fn run_resolved_generic<T>(
        &mut self,
        batch: &[Vec<[T; 4]>],
        binary_precision: Option<u32>,
        selected_helicity_ids: Option<&BTreeSet<String>>,
        selected_color_ids: Option<&BTreeSet<String>>,
    ) -> PyResult<(ResolvedValues<T>, RuntimeProfile)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        if self.lc_topology_replay_enabled {
            return Err(PyValueError::new_err(
                "resolved evaluation for LC topology-replay artifacts requires regeneration with resolved replay metadata",
            ));
        }
        if self.stages.is_none() || self.amplitude_stage.is_none() {
            return Err(self.execution_unavailable_error());
        }
        let physics = self.physics.clone().ok_or_else(|| {
            PyValueError::new_err(
                "resolved evaluation is unavailable for this older schema-v2 artifact; regenerate it with pyAmpliCol",
            )
        })?;
        let total_start = Instant::now();
        let n_points = batch.len();
        let mut state = vec![complex_zero::<T>(); n_points * self.parameter_count];
        let sources = &self.sources;
        let momentum_slots = &self.momentum_slots;
        let external_count = self.external_count;
        let external_is_initial = &self.external_is_initial;
        let particle_masses = &self.particle_masses;
        let value_parameter_count = self.value_parameter_count;
        let model_parameter_start = self.value_parameter_count + self.momentum_parameter_count;
        let model_parameter_values = &self.model_parameter_values_f64;

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

        let model_parameter_start_time = Instant::now();
        for row in 0..n_points {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_model_parameters_row_generic(
                model_parameter_start,
                model_parameter_values,
                row_state,
            )?;
        }
        let model_parameter_setup_s = model_parameter_start_time.elapsed().as_secs_f64();

        let mut stage_input_pack_s = 0.0;
        let mut stage_evaluator_call_s = 0.0;
        let mut stage_evaluator_s = 0.0;
        let mut output_assign_s = 0.0;
        let mut stage_input_pack_by_stage_s = Vec::new();
        let mut stage_evaluator_call_by_stage_s = Vec::new();
        let mut stage_output_assign_by_stage_s = Vec::new();
        for stage in self.stages.as_mut().expect("generic stages checked") {
            let (pack_s, eval_s, assign_s) = stage.evaluate_generic_into_state(
                n_points,
                self.parameter_count,
                state.as_mut_slice(),
                binary_precision,
            )?;
            stage_input_pack_s += pack_s;
            stage_evaluator_call_s += eval_s;
            stage_evaluator_s += pack_s + eval_s;
            output_assign_s += assign_s;
            stage_input_pack_by_stage_s.push(pack_s);
            stage_evaluator_call_by_stage_s.push(eval_s);
            stage_output_assign_by_stage_s.push(assign_s);
        }

        let reduction_start = Instant::now();
        let (resolved, amplitude_input_pack_s, amplitude_evaluator_call_s) = self
            .amplitude_stage
            .as_mut()
            .expect("generic amplitude stage checked")
            .evaluate_resolved_generic(
                n_points,
                state.as_slice(),
                binary_precision,
                &physics,
                self.normalization_factor,
                selected_helicity_ids,
                selected_color_ids,
            )?;
        let reduction_s = reduction_start.elapsed().as_secs_f64();
        let amplitude_evaluator_s = amplitude_input_pack_s + amplitude_evaluator_call_s;
        Ok((
            resolved,
            RuntimeProfile {
                source_fill_s,
                momentum_setup_s: momentum_setup_s + model_parameter_setup_s,
                stage_input_pack_s,
                stage_evaluator_call_s,
                stage_evaluator_s,
                output_assign_s,
                amplitude_input_pack_s,
                amplitude_evaluator_call_s,
                amplitude_evaluator_s,
                reduction_s,
                total_s: total_start.elapsed().as_secs_f64(),
                stage_input_pack_by_stage_s,
                stage_evaluator_call_by_stage_s,
                stage_output_assign_by_stage_s,
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
        let model_parameter_start = self.value_parameter_count + self.momentum_parameter_count;
        let model_parameter_values = &self.model_parameter_values_f64;

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
        for row in 0..n_points {
            let row_state =
                &mut state[row * self.parameter_count..(row + 1) * self.parameter_count];
            Self::fill_model_parameters_row(
                model_parameter_start,
                model_parameter_values,
                row_state,
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
                row[start + component] = c_generic(momentum[component].clone(), T::new_zero());
            }
        }
        Ok(())
    }

    fn fill_model_parameters_row(
        model_parameter_start: usize,
        model_parameter_values: &[f64],
        row: &mut [Complex<f64>],
    ) -> PyResult<()> {
        let stop = model_parameter_start + model_parameter_values.len();
        if stop > row.len() {
            return Err(PyValueError::new_err(
                "generic model-parameter block exceeds runtime row length",
            ));
        }
        for (index, value) in model_parameter_values.iter().enumerate() {
            row[model_parameter_start + index] = c64(*value, 0.0);
        }
        Ok(())
    }

    fn fill_model_parameters_row_generic<T>(
        model_parameter_start: usize,
        model_parameter_values: &[f64],
        row: &mut [Complex<T>],
    ) -> PyResult<()>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let stop = model_parameter_start + model_parameter_values.len();
        if stop > row.len() {
            return Err(PyValueError::new_err(
                "generic model-parameter block exceeds runtime row length",
            ));
        }
        for (index, value) in model_parameter_values.iter().enumerate() {
            row[model_parameter_start + index] = c_generic(T::from(*value), T::new_zero());
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
            let wave = if source.wavefunction_kind == "fermion"
                || (source.wavefunction_kind.is_empty() && is_fermion_pdg(source.particle_id))
            {
                let mass = particle_mass_from_map(particle_masses, source.particle_id);
                if source.particle_id < 0 {
                    ext_antiquark_dirac_massive(momentum, source.source_helicity, mass)
                } else {
                    ext_quark_dirac_massive(momentum, source.source_helicity, mass)
                }
            } else if (source.wavefunction_kind == "vector"
                && particle_mass_from_map(particle_masses, source.particle_id) == 0.0)
                || (source.wavefunction_kind.is_empty()
                    && (source.particle_id.abs() == 21 || source.particle_id == 22))
            {
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
        if source.dimension == 16 && source.wavefunction_kind == "spin2" {
            if out.len() != 16 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 16 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            let wave = ext_spin2(
                momentum,
                source.source_helicity,
                particle_mass_from_map(particle_masses, source.particle_id),
            )?;
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
            let wave = if source.wavefunction_kind == "fermion"
                || (source.wavefunction_kind.is_empty() && is_fermion_pdg(source.particle_id))
            {
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
            } else if (source.wavefunction_kind == "vector"
                && particle_mass_from_map(particle_masses, source.particle_id) == 0.0)
                || (source.wavefunction_kind.is_empty()
                    && (source.particle_id.abs() == 21 || source.particle_id == 22))
            {
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
        if source.dimension == 16 && source.wavefunction_kind == "spin2" {
            if out.len() != 16 {
                return Err(PyValueError::new_err(format!(
                    "generic source {} expected dimension 16 but slot has length {}",
                    source.source_id,
                    out.len()
                )));
            }
            let wave = ext_spin2_generic(
                &momentum,
                source.source_helicity,
                T::from(particle_mass_from_map(particle_masses, source.particle_id)),
            )?;
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

#[cfg(feature = "python")]
#[pyclass(module = "rusticol", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct ExternalParticle {
    label: usize,
    index: usize,
    side: String,
    role: String,
    particle: String,
    outgoing_particle: String,
    pdg: i32,
    outgoing_pdg: i32,
    particle_class: String,
    momentum_slot: usize,
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol", frozen, get_all, from_py_object)]
#[derive(Clone)]
struct HelicityConfiguration {
    id: String,
    index: usize,
    helicities: Vec<i32>,
    representative_id: String,
    computed: bool,
    structural_zero: bool,
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol", frozen, get_all, from_py_object)]
#[derive(Clone)]
struct ColorFlow {
    id: String,
    index: usize,
    word: Vec<usize>,
    representative_id: String,
    computed: bool,
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct ContractedColorComponent {
    id: String,
    index: usize,
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol", frozen, get_all, skip_from_py_object)]
#[derive(Clone)]
struct ModelParameter {
    name: String,
    kind: String,
    parameter_index: usize,
    default: f64,
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol", frozen, skip_from_py_object)]
#[derive(Clone)]
struct ProcessPhysics {
    process: String,
    process_key: String,
    color_accuracy: String,
    external_particles: Vec<ExternalParticle>,
    helicities: Vec<HelicityConfiguration>,
    color_flows: Vec<ColorFlow>,
    contracted_color_components: Vec<ContractedColorComponent>,
    model_parameters: Vec<ModelParameter>,
    helicity_coverage: String,
    color_coverage: String,
}

#[cfg(feature = "python")]
#[pymethods]
#[cfg(feature = "python")]
#[cfg(feature = "python")]
impl ProcessPhysics {
    #[getter]
    fn process(&self) -> &str {
        &self.process
    }

    #[getter]
    fn process_key(&self) -> &str {
        &self.process_key
    }

    #[getter]
    fn color_accuracy(&self) -> &str {
        &self.color_accuracy
    }

    #[getter]
    fn external_particles(&self) -> Vec<ExternalParticle> {
        self.external_particles.clone()
    }

    #[getter]
    fn helicities(&self) -> Vec<HelicityConfiguration> {
        self.helicities.clone()
    }

    #[getter]
    fn color_flows(&self) -> Vec<ColorFlow> {
        self.color_flows.clone()
    }

    #[getter]
    fn contracted_color_components(&self) -> Vec<ContractedColorComponent> {
        self.contracted_color_components.clone()
    }

    #[getter]
    fn model_parameters(&self) -> Vec<ModelParameter> {
        self.model_parameters.clone()
    }

    #[getter]
    fn helicity_coverage(&self) -> &str {
        &self.helicity_coverage
    }

    #[getter]
    fn color_coverage(&self) -> &str {
        &self.color_coverage
    }
}

#[cfg(feature = "python")]
impl ProcessPhysics {
    fn from_manifest(manifest: &PhysicsManifestV2) -> Self {
        Self {
            process: manifest.process.clone(),
            process_key: manifest.process_key.clone(),
            color_accuracy: manifest.color_accuracy.clone(),
            external_particles: manifest
                .external_particles
                .iter()
                .map(|particle| ExternalParticle {
                    label: particle.label,
                    index: particle.index,
                    side: particle.side.clone(),
                    role: particle.role.clone(),
                    particle: particle.particle.clone(),
                    outgoing_particle: particle.outgoing_particle.clone(),
                    pdg: particle.pdg,
                    outgoing_pdg: particle.outgoing_pdg,
                    particle_class: particle.particle_class.clone(),
                    momentum_slot: particle.momentum_slot,
                })
                .collect(),
            helicities: manifest
                .helicities
                .iter()
                .map(|helicity| HelicityConfiguration {
                    id: helicity.id.clone(),
                    index: helicity.index,
                    helicities: helicity.helicities.clone(),
                    representative_id: helicity.representative_id.clone(),
                    computed: helicity.computed,
                    structural_zero: helicity.structural_zero,
                })
                .collect(),
            color_flows: manifest
                .color_components
                .iter()
                .filter(|color| color.kind == "lc-flow")
                .map(|color| ColorFlow {
                    id: color.id.clone(),
                    index: color.index,
                    word: color.word.clone(),
                    representative_id: color.representative_id.clone(),
                    computed: color.computed,
                })
                .collect(),
            contracted_color_components: manifest
                .color_components
                .iter()
                .filter(|color| color.kind == "contracted")
                .map(|color| ContractedColorComponent {
                    id: color.id.clone(),
                    index: color.index,
                })
                .collect(),
            model_parameters: manifest
                .model_parameters
                .iter()
                .map(|parameter| ModelParameter {
                    name: parameter.name.clone(),
                    kind: parameter.kind.clone(),
                    parameter_index: parameter.parameter_index,
                    default: parameter.default,
                })
                .collect(),
            helicity_coverage: manifest.coverage.helicities.clone(),
            color_coverage: manifest.coverage.color.clone(),
        }
    }
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol")]
struct ResolvedEvaluation {
    values: Py<PyAny>,
    totals: Py<PyAny>,
    shape: (usize, usize, usize),
    helicities: Vec<HelicityConfiguration>,
    color_flows: Vec<ColorFlow>,
    contracted_color_components: Vec<ContractedColorComponent>,
}

#[cfg(feature = "python")]
#[pymethods]
impl ResolvedEvaluation {
    #[getter]
    fn values(&self, py: Python<'_>) -> Py<PyAny> {
        self.values.clone_ref(py)
    }

    #[getter]
    fn shape(&self) -> (usize, usize, usize) {
        self.shape
    }

    #[getter]
    fn helicities(&self) -> Vec<HelicityConfiguration> {
        self.helicities.clone()
    }

    #[getter]
    fn color_flows(&self) -> Vec<ColorFlow> {
        self.color_flows.clone()
    }

    #[getter]
    fn contracted_color_components(&self) -> Vec<ContractedColorComponent> {
        self.contracted_color_components.clone()
    }

    fn total(&self, py: Python<'_>) -> Py<PyAny> {
        self.totals.clone_ref(py)
    }
}

#[cfg(feature = "python")]
#[derive(Clone, Copy)]
enum PublicSelectorKind {
    Helicity,
    ColorFlow,
}

#[cfg(feature = "python")]
fn parse_public_selector_ids(
    value: Option<&Bound<'_, PyAny>>,
    kind: PublicSelectorKind,
) -> PyResult<Option<BTreeSet<String>>> {
    let Some(value) = value else {
        return Ok(None);
    };
    if value.is_none() {
        return Ok(None);
    }
    if let Some(id) = public_selector_item_id(value, kind)? {
        return Ok(Some(BTreeSet::from([id])));
    }
    let mut ids = BTreeSet::new();
    for item in value.try_iter().map_err(|_| {
        PyValueError::new_err(
            "resolved selectors must be an id, a typed metadata object, or an iterable of them",
        )
    })? {
        let item = item?;
        let id = public_selector_item_id(&item, kind)?.ok_or_else(|| {
            PyValueError::new_err("resolved selector iterable contains an unsupported item")
        })?;
        ids.insert(id);
    }
    if ids.is_empty() {
        return Err(PyValueError::new_err(
            "resolved selector iterable must not be empty",
        ));
    }
    Ok(Some(ids))
}

#[cfg(feature = "python")]
fn public_selector_item_id(
    value: &Bound<'_, PyAny>,
    kind: PublicSelectorKind,
) -> PyResult<Option<String>> {
    if let Ok(id) = value.extract::<String>() {
        return Ok(Some(id));
    }
    match kind {
        PublicSelectorKind::Helicity => {
            if let Ok(item) = value.extract::<PyRef<'_, HelicityConfiguration>>() {
                return Ok(Some(item.id.clone()));
            }
        }
        PublicSelectorKind::ColorFlow => {
            if let Ok(item) = value.extract::<PyRef<'_, ColorFlow>>() {
                return Ok(Some(item.id.clone()));
            }
        }
    }
    Ok(None)
}

#[cfg(feature = "python")]
fn parse_python_model_parameter_mapping(
    mapping: &Bound<'_, PyAny>,
) -> PyResult<BTreeMap<String, (f64, f64)>> {
    let items = mapping.call_method0("items").map_err(|_| {
        PyValueError::new_err("model parameter updates must be supplied as a mapping")
    })?;
    let mut result = BTreeMap::new();
    for item in items.try_iter()? {
        let item = item?;
        let pair = item.cast::<PyTuple>().map_err(|_| {
            PyValueError::new_err("model parameter mapping items must be name/value pairs")
        })?;
        if pair.len() != 2 {
            return Err(PyValueError::new_err(
                "model parameter mapping items must contain exactly two values",
            ));
        }
        let name = pair.get_item(0)?.extract::<String>()?;
        let value = pair.get_item(1)?;
        let components = if let Ok(real) = value.extract::<f64>() {
            (real, 0.0)
        } else if let Ok((real, imaginary)) = value.extract::<(f64, f64)>() {
            (real, imaginary)
        } else if value.hasattr("real")? && value.hasattr("imag")? {
            (
                value.getattr("real")?.extract::<f64>()?,
                value.getattr("imag")?.extract::<f64>()?,
            )
        } else {
            return Err(PyValueError::new_err(format!(
                "model parameter {name:?} must be a real, complex, or (real, imag) pair",
            )));
        };
        if result.insert(name.clone(), components).is_some() {
            return Err(PyValueError::new_err(format!(
                "duplicate model parameter update {name:?}",
            )));
        }
    }
    Ok(result)
}

#[cfg(feature = "python")]
fn resolved_f64_to_python(
    py: Python<'_>,
    resolved: ResolvedValues<f64>,
    physics: &PhysicsManifestV2,
) -> PyResult<ResolvedEvaluation> {
    let helicity_count = resolved.helicity_indices.len();
    let color_count = resolved.color_indices.len();
    let mut rows = Vec::with_capacity(resolved.point_count);
    let mut totals = Vec::with_capacity(resolved.point_count);
    for point_index in 0..resolved.point_count {
        let mut helicity_rows = Vec::with_capacity(helicity_count);
        let mut total = 0.0;
        for helicity_offset in 0..helicity_count {
            let start = (point_index * helicity_count + helicity_offset) * color_count;
            let color_values = resolved.values[start..start + color_count].to_vec();
            total += color_values.iter().sum::<f64>();
            helicity_rows.push(color_values);
        }
        rows.push(helicity_rows);
        totals.push(total);
    }
    let typed = ProcessPhysics::from_manifest(physics);
    Ok(ResolvedEvaluation {
        values: rows.into_py_any(py)?,
        totals: totals.into_pyarray(py).into_any().unbind(),
        shape: (resolved.point_count, helicity_count, color_count),
        helicities: resolved
            .helicity_indices
            .iter()
            .map(|index| typed.helicities[*index].clone())
            .collect(),
        color_flows: resolved
            .color_indices
            .iter()
            .filter_map(|index| {
                let color = &physics.color_components[*index];
                (color.kind == "lc-flow").then(|| ColorFlow {
                    id: color.id.clone(),
                    index: color.index,
                    word: color.word.clone(),
                    representative_id: color.representative_id.clone(),
                    computed: color.computed,
                })
            })
            .collect(),
        contracted_color_components: resolved
            .color_indices
            .iter()
            .filter_map(|index| {
                let color = &physics.color_components[*index];
                (color.kind == "contracted").then(|| ContractedColorComponent {
                    id: color.id.clone(),
                    index: color.index,
                })
            })
            .collect(),
    })
}

#[cfg(feature = "python")]
fn resolved_decimals_to_python<T>(
    py: Python<'_>,
    resolved: ResolvedValues<T>,
    physics: &PhysicsManifestV2,
    decimal_digits: u32,
) -> PyResult<ResolvedEvaluation>
where
    T: Real + std::fmt::Display + std::fmt::LowerExp,
{
    let decimal = py.import("decimal")?.getattr("Decimal")?;
    let digits = decimal_digits as usize;
    let helicity_count = resolved.helicity_indices.len();
    let color_count = resolved.color_indices.len();
    let mut point_rows: Vec<Py<PyAny>> = Vec::with_capacity(resolved.point_count);
    let mut totals: Vec<Py<PyAny>> = Vec::with_capacity(resolved.point_count);
    for point_index in 0..resolved.point_count {
        let mut helicity_rows: Vec<Py<PyAny>> = Vec::with_capacity(helicity_count);
        let mut total = T::new_zero();
        for helicity_offset in 0..helicity_count {
            let start = (point_index * helicity_count + helicity_offset) * color_count;
            let mut color_values: Vec<Py<PyAny>> = Vec::with_capacity(color_count);
            for value in &resolved.values[start..start + color_count] {
                total += value.clone();
                color_values.push(
                    decimal
                        .call1((format!("{value:.digits$e}"),))?
                        .into_any()
                        .unbind(),
                );
            }
            helicity_rows.push(PyList::new(py, color_values)?.into_any().unbind());
        }
        point_rows.push(PyList::new(py, helicity_rows)?.into_any().unbind());
        totals.push(
            decimal
                .call1((format!("{total:.digits$e}"),))?
                .into_any()
                .unbind(),
        );
    }
    let typed = ProcessPhysics::from_manifest(physics);
    Ok(ResolvedEvaluation {
        values: PyList::new(py, point_rows)?.into_any().unbind(),
        totals: PyList::new(py, totals)?.into_any().unbind(),
        shape: (resolved.point_count, helicity_count, color_count),
        helicities: resolved
            .helicity_indices
            .iter()
            .map(|index| typed.helicities[*index].clone())
            .collect(),
        color_flows: resolved
            .color_indices
            .iter()
            .filter_map(|index| {
                let color = &physics.color_components[*index];
                (color.kind == "lc-flow").then(|| ColorFlow {
                    id: color.id.clone(),
                    index: color.index,
                    word: color.word.clone(),
                    representative_id: color.representative_id.clone(),
                    computed: color.computed,
                })
            })
            .collect(),
        contracted_color_components: resolved
            .color_indices
            .iter()
            .filter_map(|index| {
                let color = &physics.color_components[*index];
                (color.kind == "contracted").then(|| ContractedColorComponent {
                    id: color.id.clone(),
                    index: color.index,
                })
            })
            .collect(),
    })
}

#[cfg(feature = "python")]
#[pyclass(module = "rusticol")]
struct Runtime {
    root: PathBuf,
    generic_runtime: GenericRuntimeV2,
    selected_process_key: Option<String>,
    selected_process: Option<String>,
    input_crossing_map: Option<Vec<InputCrossingMapEntry>>,
    crossing_alias_of: Option<String>,
    last_profile: RuntimeProfile,
    warnings_muted: bool,
    warned_kinds: BTreeSet<String>,
}

struct ProcessSetSelection {
    root: PathBuf,
    selected_key: Option<String>,
    selected_process: Option<String>,
    input_crossing_map: Option<Vec<InputCrossingMapEntry>>,
    crossing_alias_of: Option<String>,
    physics: Option<PhysicsManifestV2>,
}

#[cfg(feature = "python")]
fn load_rusticol_runtime(
    process_dir: &str,
    process_key: Option<&str>,
    model_parameters_path: Option<&str>,
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
    let header: ArtifactManifestHeader =
        read_json_manifest(&manifest_path).map_err(PyValueError::new_err)?;
    if header.schema_version != 2 || header.kind != "pyamplicol-generic-dag-process" {
        return Err(PyValueError::new_err(format!(
            "Rusticol supports only schema-v2 generic DAG artifacts; got kind {:?} \
             schema_version {}. Schema-v1 execution has been removed; regenerate this process \
             with pyAmpliCol generate-process.",
            header.kind, header.schema_version
        )));
    }
    let generic: GenericProcessManifestV2 =
        read_json_manifest(&manifest_path).map_err(|error| {
            PyValueError::new_err(format!(
                "could not parse generic DAG process manifest schema v2: {error}"
            ))
        })?;
    let mut generic_runtime = load_generic_schema_v2_manifest(generic, &root)?;
    if let Some(metadata) = selection.physics.clone() {
        generic_runtime.physics = Some(PhysicsRuntimeV2::new(metadata)?);
    }
    if let Some(path) = model_parameters_path {
        generic_runtime.apply_model_parameter_json_path(&PathBuf::from(path))?;
    }
    Ok(Runtime {
        root,
        generic_runtime,
        selected_process_key: selection.selected_key,
        selected_process: selection.selected_process,
        input_crossing_map: selection.input_crossing_map,
        crossing_alias_of: selection.crossing_alias_of,
        last_profile: RuntimeProfile::default(),
        warnings_muted: false,
        warned_kinds: BTreeSet::new(),
    })
}

#[cfg(feature = "python")]
#[pymethods]
impl Runtime {
    #[classmethod]
    #[pyo3(signature = (process_dir, process_key=None, model_parameters=None))]
    fn load(
        _cls: &Bound<'_, pyo3::types::PyType>,
        process_dir: &str,
        process_key: Option<&str>,
        model_parameters: Option<&str>,
    ) -> PyResult<Self> {
        load_rusticol_runtime(process_dir, process_key, model_parameters)
    }

    #[getter]
    fn process(&self) -> PyResult<String> {
        Ok(self
            .selected_process
            .clone()
            .unwrap_or_else(|| self.generic_runtime.process.clone()))
    }

    fn metadata<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        let generic = &self.generic_runtime;
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
        dict.set_item("model_parameter_count", generic.model_parameter_count)?;
        dict.set_item(
            "model_parameter_names",
            generic
                .model_parameter_runtime_slots
                .keys()
                .cloned()
                .collect::<Vec<_>>(),
        )?;
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
        Ok(dict)
    }

    #[getter]
    fn physics(&self) -> PyResult<ProcessPhysics> {
        let physics = self.generic_runtime.physics.as_ref().ok_or_else(|| {
            PyValueError::new_err(
                "this older schema-v2 artifact has no resolved physics metadata; regenerate it with pyAmpliCol",
            )
        })?;
        Ok(ProcessPhysics::from_manifest(&physics.manifest))
    }

    fn set_model_parameters(&mut self, mapping: &Bound<'_, PyAny>) -> PyResult<()> {
        let overrides = parse_python_model_parameter_mapping(mapping)?;
        self.generic_runtime
            .apply_model_parameter_overrides(&overrides)
    }

    #[pyo3(signature = (name, real, imaginary=0.0))]
    fn set_model_parameter(&mut self, name: &str, real: f64, imaginary: f64) -> PyResult<()> {
        self.generic_runtime
            .apply_model_parameter_overrides(&BTreeMap::from([(
                name.to_string(),
                (real, imaginary),
            )]))
    }

    fn mute_warnings(&mut self) {
        self.warnings_muted = true;
    }

    fn unmute_warnings(&mut self) {
        self.warnings_muted = false;
    }

    fn evaluate<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let values = self.evaluate_f64_values(py, momenta)?;
        f64_values_to_numpy_or_list(py, values)
    }

    #[pyo3(signature = (momenta, helicities=None, color_flows=None))]
    fn evaluate_resolved(
        &mut self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
        helicities: Option<&Bound<'_, PyAny>>,
        color_flows: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<ResolvedEvaluation> {
        let selected_helicity_ids =
            parse_public_selector_ids(helicities, PublicSelectorKind::Helicity)?;
        let selected_color_ids =
            parse_public_selector_ids(color_flows, PublicSelectorKind::ColorFlow)?;
        self.emit_resolved_warnings(
            py,
            selected_helicity_ids.as_ref(),
            selected_color_ids.as_ref(),
        )?;
        let generic = &self.generic_runtime;
        if selected_color_ids.is_some() && generic.color_accuracy != "lc" {
            return Err(PyValueError::new_err(
                "LC color-flow selection is unavailable for NLC/full artifacts; their resolved color axis is contracted",
            ));
        }
        let external_count = generic.external_count;
        let physics = generic
            .physics
            .as_ref()
            .ok_or_else(|| {
                PyValueError::new_err(
                    "resolved evaluation is unavailable for this older schema-v2 artifact; regenerate it with pyAmpliCol",
                )
            })?
            .manifest
            .clone();
        let batch = self.generic_batch_from_python(py, momenta, external_count)?;
        let (resolved, profile) = self.generic_runtime.run_resolved_f64(
            &batch,
            selected_helicity_ids.as_ref(),
            selected_color_ids.as_ref(),
        )?;
        self.last_profile = profile;
        resolved_f64_to_python(py, resolved, &physics)
    }

    fn evaluate_color_sectors<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
        color_sector_ids: &Bound<'py, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let selected_color_sector_ids = parse_required_color_sector_ids(color_sector_ids)?;
        let values =
            self.evaluate_f64_values_selected(py, momenta, Some(&selected_color_sector_ids))?;
        f64_values_to_numpy_or_list(py, values)
    }

    fn raw_amplitudes<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
    ) -> PyResult<Py<PyAny>> {
        let external_count = self.generic_runtime.external_count;
        let batch = self.generic_batch_from_python(py, momenta, external_count)?;
        let generic = &mut self.generic_runtime;
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
                let external_count = self.generic_runtime.external_count;
                let batch = self.generic_batch_double_from_python(momenta, external_count)?;
                let (values, profile) = self.generic_runtime.run_double(&batch)?;
                self.last_profile = profile;
                decimals_to_python(py, values, decimal_digit_precision)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let external_count = self.generic_runtime.external_count;
                let batch = self.generic_batch_float_from_python(
                    momenta,
                    binary_precision,
                    external_count,
                )?;
                let (values, profile) = self.generic_runtime.run_float(&batch, binary_precision)?;
                self.last_profile = profile;
                decimals_to_python(py, values, decimal_precision)
            }
        }
    }

    #[pyo3(signature = (momenta, decimal_digit_precision, helicities=None, color_flows=None))]
    fn evaluate_resolved_with_prec(
        &mut self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
        decimal_digit_precision: u32,
        helicities: Option<&Bound<'_, PyAny>>,
        color_flows: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<ResolvedEvaluation> {
        if decimal_digit_precision == 16 {
            return self.evaluate_resolved(py, momenta, helicities, color_flows);
        }
        let selected_helicity_ids =
            parse_public_selector_ids(helicities, PublicSelectorKind::Helicity)?;
        let selected_color_ids =
            parse_public_selector_ids(color_flows, PublicSelectorKind::ColorFlow)?;
        self.emit_resolved_warnings(
            py,
            selected_helicity_ids.as_ref(),
            selected_color_ids.as_ref(),
        )?;
        let generic = &self.generic_runtime;
        if selected_color_ids.is_some() && generic.color_accuracy != "lc" {
            return Err(PyValueError::new_err(
                "LC color-flow selection is unavailable for NLC/full artifacts; their resolved color axis is contracted",
            ));
        }
        let external_count = generic.external_count;
        let physics = generic
            .physics
            .as_ref()
            .ok_or_else(|| {
                PyValueError::new_err(
                    "resolved evaluation is unavailable for this older schema-v2 artifact; regenerate it with pyAmpliCol",
                )
            })?
            .manifest
            .clone();
        match PrecisionMode::from_decimal_digits(decimal_digit_precision)? {
            PrecisionMode::F64 => unreachable!("precision 16 returned above"),
            PrecisionMode::DoubleDouble => {
                let batch = self.generic_batch_double_from_python(momenta, external_count)?;
                let (resolved, profile) = self.generic_runtime.run_resolved_generic(
                    &batch,
                    None,
                    selected_helicity_ids.as_ref(),
                    selected_color_ids.as_ref(),
                )?;
                self.last_profile = profile;
                resolved_decimals_to_python(py, resolved, &physics, decimal_digit_precision)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let batch = self.generic_batch_float_from_python(
                    momenta,
                    binary_precision,
                    external_count,
                )?;
                let (resolved, profile) = self.generic_runtime.run_resolved_generic(
                    &batch,
                    Some(binary_precision),
                    selected_helicity_ids.as_ref(),
                    selected_color_ids.as_ref(),
                )?;
                self.last_profile = profile;
                resolved_decimals_to_python(py, resolved, &physics, decimal_precision)
            }
        }
    }

    #[pyo3(signature = (momenta, precision = 16, include_values = true, color_sector_ids = None))]
    fn profile<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
        precision: u32,
        include_values: bool,
        color_sector_ids: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let selected_color_sector_ids = parse_optional_color_sector_ids(color_sector_ids)?;
        if selected_color_sector_ids.is_some() && precision != 16 {
            return Err(PyValueError::new_err(
                "LC color-sector runtime selection is currently supported only at precision 16",
            ));
        }
        let (points, values, profile) = match PrecisionMode::from_decimal_digits(precision)? {
            PrecisionMode::F64 => {
                let external_count = self.generic_runtime.external_count;
                let batch = self.generic_batch_from_python(py, momenta, external_count)?;
                let points = batch.len();
                let (values, profile) = self
                    .generic_runtime
                    .run_f64_selected(&batch, selected_color_sector_ids.as_ref())?;
                let values = if include_values {
                    values.into_py_any(py)?
                } else {
                    py.None()
                };
                (points, values, profile)
            }
            PrecisionMode::DoubleDouble => {
                let external_count = self.generic_runtime.external_count;
                let batch = self.generic_batch_double_from_python(momenta, external_count)?;
                let points = batch.len();
                let (values, profile) = self.generic_runtime.run_double(&batch)?;
                let values = if include_values {
                    decimals_to_python(py, values, precision)?
                } else {
                    py.None()
                };
                (points, values, profile)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let external_count = self.generic_runtime.external_count;
                let batch = self.generic_batch_float_from_python(
                    momenta,
                    binary_precision,
                    external_count,
                )?;
                let points = batch.len();
                let (values, profile) = self.generic_runtime.run_float(&batch, binary_precision)?;
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
        dict.set_item(
            "color_sector_ids",
            selected_color_sector_ids
                .as_ref()
                .map(|ids| ids.iter().copied().collect::<Vec<_>>()),
        )?;
        dict.set_item("source_fill_time_s", profile.source_fill_s)?;
        dict.set_item("momentum_setup_time_s", profile.momentum_setup_s)?;
        dict.set_item("parameter_pack_time_s", 0.0)?;
        dict.set_item("stage_input_pack_time_s", profile.stage_input_pack_s)?;
        dict.set_item(
            "stage_evaluator_call_time_s",
            profile.stage_evaluator_call_s,
        )?;
        dict.set_item("stage_evaluator_time_s", profile.stage_evaluator_s)?;
        dict.set_item("output_transfer_time_s", profile.output_assign_s)?;
        dict.set_item("output_assign_time_s", profile.output_assign_s)?;
        dict.set_item(
            "amplitude_input_pack_time_s",
            profile.amplitude_input_pack_s,
        )?;
        dict.set_item(
            "amplitude_evaluator_call_time_s",
            profile.amplitude_evaluator_call_s,
        )?;
        dict.set_item("amplitude_evaluator_time_s", profile.amplitude_evaluator_s)?;
        dict.set_item("result_reduction_time_s", profile.reduction_s)?;
        dict.set_item("total_time_s", profile.total_s)?;
        dict.set_item(
            "stage_input_pack_by_stage_time_s",
            profile.stage_input_pack_by_stage_s.clone(),
        )?;
        dict.set_item(
            "stage_evaluator_call_by_stage_time_s",
            profile.stage_evaluator_call_by_stage_s.clone(),
        )?;
        dict.set_item(
            "stage_output_assign_by_stage_time_s",
            profile.stage_output_assign_by_stage_s.clone(),
        )?;
        dict.set_item("current_rss_bytes", memory.current_rss_bytes)?;
        dict.set_item("peak_rss_bytes", memory.peak_rss_bytes)?;
        dict.set_item(
            "core_evaluator_time_s",
            profile.stage_evaluator_s + profile.amplitude_evaluator_s,
        )?;
        Ok(dict)
    }
}

#[cfg(feature = "python")]
enum PrecisionMode {
    F64,
    DoubleDouble,
    Arbitrary(u32),
}

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
impl Runtime {
    fn emit_resolved_warnings(
        &mut self,
        py: Python<'_>,
        selected_helicity_ids: Option<&BTreeSet<String>>,
        selected_color_ids: Option<&BTreeSet<String>>,
    ) -> PyResult<()> {
        if self.warnings_muted {
            return Ok(());
        }
        let Some(physics) = self.generic_runtime.physics.as_ref() else {
            return Ok(());
        };
        let mut warnings = Vec::new();
        if physics.manifest.coverage.helicities != "complete" {
            warnings.push((
                "incomplete-helicity-coverage",
                "resolved evaluation contains only the helicities represented by this artifact",
            ));
        }
        if physics.manifest.coverage.color != "complete" {
            warnings.push((
                "incomplete-color-coverage",
                "resolved evaluation contains only the color components represented by this artifact",
            ));
        }
        let reduction_only_helicity = selected_helicity_ids.is_some_and(|ids| {
            ids.iter().any(|id| {
                physics
                    .helicity_index_by_id
                    .get(id)
                    .and_then(|index| physics.manifest.helicities.get(*index))
                    .is_some_and(|item| !item.computed)
            })
        });
        let reduction_only_color = selected_color_ids.is_some_and(|ids| {
            ids.iter().any(|id| {
                physics
                    .color_index_by_id
                    .get(id)
                    .and_then(|index| physics.manifest.color_components.get(*index))
                    .is_some_and(|item| !item.computed)
            })
        });
        if reduction_only_helicity || reduction_only_color {
            warnings.push((
                "reduction-only-selection",
                "the selected resolved component reuses an exact symmetry representative",
            ));
        }
        if warnings.is_empty() {
            return Ok(());
        }
        let warnings_module = py.import("warnings")?;
        for (kind, message) in warnings {
            if self.warned_kinds.insert(kind.to_string()) {
                warnings_module.call_method1("warn", (message,))?;
            }
        }
        Ok(())
    }

    fn generic_batch_from_python(
        &self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
        expected_legs: usize,
    ) -> PyResult<Vec<Vec<[f64; 4]>>> {
        let batch = batch_momenta_dynamic(py, momenta, expected_legs)?;
        apply_input_crossing_map(batch, expected_legs, self.input_crossing_map.as_deref())
    }

    fn generic_batch_double_from_python(
        &self,
        momenta: &Bound<'_, PyAny>,
        expected_legs: usize,
    ) -> PyResult<Vec<Vec<[DoubleFloat; 4]>>> {
        let batch = batch_momenta_double(momenta, expected_legs)?;
        apply_input_crossing_map_generic(&batch, expected_legs, self.input_crossing_map.as_deref())
    }

    fn generic_batch_float_from_python(
        &self,
        momenta: &Bound<'_, PyAny>,
        binary_precision: u32,
        expected_legs: usize,
    ) -> PyResult<Vec<Vec<[Float; 4]>>> {
        let batch = batch_momenta_float(momenta, binary_precision, expected_legs)?;
        apply_input_crossing_map_generic(&batch, expected_legs, self.input_crossing_map.as_deref())
    }

    fn evaluate_f64_values(
        &mut self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<f64>> {
        self.evaluate_f64_values_selected(py, momenta, None)
    }

    fn evaluate_f64_values_selected(
        &mut self,
        py: Python<'_>,
        momenta: &Bound<'_, PyAny>,
        selected_color_sector_ids: Option<&BTreeSet<i64>>,
    ) -> PyResult<Vec<f64>> {
        let external_count = self.generic_runtime.external_count;
        let batch = self.generic_batch_from_python(py, momenta, external_count)?;
        let (values, profile) = self
            .generic_runtime
            .run_f64_selected(&batch, selected_color_sector_ids)?;
        self.last_profile = profile;
        Ok(values)
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
        let (input_components, input_spans) =
            if stage.parameter_layout == "stage-local-value-momentum" {
                let mut map = vec![0usize; stage.parameter_count];
                for component in &stage.input_components {
                    map[component.parameter_index] = component.global_component;
                }
                let spans = contiguous_input_spans(&map);
                (Some(map), spans)
            } else {
                (None, Vec::new())
            };
        let output_spans = contiguous_output_spans(&outputs);
        Ok(Self {
            outputs,
            output_spans,
            input_components,
            input_spans,
            parameter_scratch_f64: Vec::new(),
            output_scratch_f64: Vec::new(),
            parameter_scratch_native2: Vec::new(),
            output_scratch_native2: Vec::new(),
            evaluator,
        })
    }

    fn evaluate_f64_into_state(
        &mut self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<f64>],
    ) -> PyResult<(f64, f64, f64)> {
        if self.evaluator.supports_native2() {
            return self.evaluate_f64_native2_into_state(batch_size, parameter_count, state);
        }
        let mut input_pack_s = 0.0;
        let eval_start;
        if let Some(input_components) = self.input_components.as_ref() {
            let local_parameter_count = input_components.len();
            let pack_start = Instant::now();
            self.parameter_scratch_f64
                .resize(batch_size * local_parameter_count, c64(0.0, 0.0));
            for row in 0..batch_size {
                let row_state = row * parameter_count;
                let row_params = row * local_parameter_count;
                if self.input_spans.is_empty() {
                    for (local_index, global_index) in input_components.iter().enumerate() {
                        self.parameter_scratch_f64[row_params + local_index] =
                            state[row_state + *global_index];
                    }
                } else {
                    for (local_start, global_start, len) in &self.input_spans {
                        let target_start = row_params + *local_start;
                        let source_start = row_state + *global_start;
                        self.parameter_scratch_f64[target_start..target_start + *len]
                            .copy_from_slice(&state[source_start..source_start + *len]);
                    }
                }
            }
            input_pack_s = pack_start.elapsed().as_secs_f64();
            eval_start = Instant::now();
            self.evaluator.evaluate_batch_into(
                batch_size,
                &self.parameter_scratch_f64,
                &mut self.output_scratch_f64,
            )?;
        } else {
            eval_start = Instant::now();
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
        Ok((
            input_pack_s,
            evaluator_s,
            assign_start.elapsed().as_secs_f64(),
        ))
    }

    fn evaluate_f64_native2_into_state(
        &mut self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<f64>],
    ) -> PyResult<(f64, f64, f64)> {
        let pack_start = Instant::now();
        pack_native2_parameters(
            batch_size,
            parameter_count,
            self.input_components.as_deref(),
            state,
            &mut self.parameter_scratch_native2,
        )?;
        let input_pack_s = pack_start.elapsed().as_secs_f64();

        let evaluator_start = Instant::now();
        self.evaluator.evaluate_native2_into(
            batch_size.div_ceil(2),
            &self.parameter_scratch_native2,
            &mut self.output_scratch_native2,
        )?;
        let evaluator_s = evaluator_start.elapsed().as_secs_f64();

        let assign_start = Instant::now();
        let output_len = self.evaluator.output_len;
        for row in 0..batch_size {
            let state_row = row * parameter_count;
            let native_row = row / 2 * output_len;
            let lane = row % 2;
            for (output_column, state_offset) in &self.outputs {
                let value = self.output_scratch_native2[native_row + *output_column];
                state[state_row + *state_offset] =
                    c64(value.re.as_array()[lane], value.im.as_array()[lane]);
            }
        }
        Ok((
            input_pack_s,
            evaluator_s,
            assign_start.elapsed().as_secs_f64(),
        ))
    }

    fn evaluate_generic_into_state<T>(
        &mut self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<(f64, f64, f64)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let mut input_pack_s = 0.0;
        let (evaluated, evaluator_s) =
            if let Some(input_components) = self.input_components.as_ref() {
                let local_parameter_count = input_components.len();
                let pack_start = Instant::now();
                let mut parameter_scratch =
                    vec![complex_zero::<T>(); batch_size * local_parameter_count];
                for row in 0..batch_size {
                    let row_state = row * parameter_count;
                    let row_params = row * local_parameter_count;
                    if self.input_spans.is_empty() {
                        for (local_index, global_index) in input_components.iter().enumerate() {
                            parameter_scratch[row_params + local_index] =
                                state[row_state + *global_index].clone();
                        }
                    } else {
                        for (local_start, global_start, len) in &self.input_spans {
                            let target_start = row_params + *local_start;
                            let source_start = row_state + *global_start;
                            parameter_scratch[target_start..target_start + *len]
                                .clone_from_slice(&state[source_start..source_start + *len]);
                        }
                    }
                }
                input_pack_s = pack_start.elapsed().as_secs_f64();
                let eval_start = Instant::now();
                let evaluated = self.evaluator.evaluate_batch_generic(
                    batch_size,
                    &parameter_scratch,
                    binary_precision,
                )?;
                (evaluated, eval_start.elapsed().as_secs_f64())
            } else {
                let eval_start = Instant::now();
                let evaluated =
                    self.evaluator
                        .evaluate_batch_generic(batch_size, state, binary_precision)?;
                (evaluated, eval_start.elapsed().as_secs_f64())
            };

        self.assign_generic_outputs(
            batch_size,
            parameter_count,
            state,
            evaluated,
            input_pack_s,
            evaluator_s,
        )
    }

    fn assign_generic_outputs<T>(
        &self,
        batch_size: usize,
        parameter_count: usize,
        state: &mut [Complex<T>],
        evaluated: Vec<Complex<T>>,
        input_pack_s: f64,
        evaluator_s: f64,
    ) -> PyResult<(f64, f64, f64)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
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
        Ok((
            input_pack_s,
            evaluator_s,
            assign_start.elapsed().as_secs_f64(),
        ))
    }
}

fn pack_native2_parameters(
    batch_size: usize,
    global_parameter_count: usize,
    input_components: Option<&[usize]>,
    state: &[Complex<f64>],
    target: &mut Vec<Complex<wide::f64x2>>,
) -> PyResult<()> {
    if batch_size == 0 {
        return Err(PyValueError::new_err(
            "native two-lane evaluation requires a non-empty batch",
        ));
    }
    if state.len() != batch_size * global_parameter_count {
        return Err(PyValueError::new_err(format!(
            "state buffer has length {}, expected {}",
            state.len(),
            batch_size * global_parameter_count
        )));
    }
    let local_parameter_count = input_components
        .map(|components| components.len())
        .unwrap_or(global_parameter_count);
    target.resize(
        batch_size.div_ceil(2) * local_parameter_count,
        Complex::new(wide::f64x2::ZERO, wide::f64x2::ZERO),
    );
    for native_row in 0..batch_size.div_ceil(2) {
        let first_row = native_row * 2;
        let second_row = usize::min(first_row + 1, batch_size - 1);
        let target_row = native_row * local_parameter_count;
        for local_index in 0..local_parameter_count {
            let global_index = input_components
                .map(|components| components[local_index])
                .unwrap_or(local_index);
            let first = state[first_row * global_parameter_count + global_index];
            let second = state[second_row * global_parameter_count + global_index];
            target[target_row + local_index] = Complex::new(
                wide::f64x2::new([first.re, second.re]),
                wide::f64x2::new([first.im, second.im]),
            );
        }
    }
    Ok(())
}

fn contiguous_input_spans(input_components: &[usize]) -> Vec<(usize, usize, usize)> {
    if input_components.is_empty() {
        return Vec::new();
    }
    let mut spans = Vec::new();
    let mut local_start = 0usize;
    let mut global_start = input_components[0];
    let mut previous_local = local_start;
    let mut previous_global = global_start;
    let mut len = 1usize;
    for (local, global) in input_components.iter().copied().enumerate().skip(1) {
        if local == previous_local + 1 && global == previous_global + 1 {
            previous_local = local;
            previous_global = global;
            len += 1;
            continue;
        }
        spans.push((local_start, global_start, len));
        local_start = local;
        global_start = global;
        previous_local = local;
        previous_global = global;
        len = 1;
    }
    spans.push((local_start, global_start, len));
    if spans.len() >= input_components.len() {
        Vec::new()
    } else {
        spans
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
        let raw_sum_all_sector_weights = amplitude_stage
            .roots
            .iter()
            .map(|root| root.all_sector_weight.unwrap_or(root.helicity_weight))
            .collect::<Vec<_>>();
        let raw_sum_color_sector_ids = amplitude_stage
            .roots
            .iter()
            .map(|root| root.color_sector_id)
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
                &raw_sum_all_sector_weights,
                &raw_sum_group_ids,
                &raw_sum_color_sector_ids,
            )?
        } else {
            Vec::new()
        };
        let color_contraction = build_color_contraction_runtime(
            amplitude_stage.color_contraction.as_ref(),
            &raw_sum_groups,
        )?;
        let (input_components, input_spans) =
            if stage.parameter_layout == "stage-local-value-momentum" {
                let mut map = vec![0usize; stage.parameter_count];
                for component in &stage.input_components {
                    map[component.parameter_index] = component.global_component;
                }
                let spans = contiguous_input_spans(&map);
                (Some(map), spans)
            } else {
                (None, Vec::new())
            };
        Ok(Self {
            output_length: amplitude_stage.output_count,
            raw_sum_weights,
            raw_sum_all_sector_weights,
            raw_sum_color_sector_ids,
            raw_sum_groups,
            has_coherent_groups,
            color_contraction,
            input_components,
            input_spans,
            parameter_scratch_f64: Vec::new(),
            output_scratch_f64: Vec::new(),
            parameter_scratch_native2: Vec::new(),
            output_scratch_native2: Vec::new(),
            evaluator: EvaluatorGroup::load(&stage.evaluator, root)?,
        })
    }

    fn evaluate_f64_into_scratch(
        &mut self,
        batch_size: usize,
        state: &[Complex<f64>],
    ) -> PyResult<(f64, f64)> {
        if self.evaluator.supports_native2() {
            let global_parameter_count = state
                .len()
                .checked_div(batch_size)
                .ok_or_else(|| PyValueError::new_err("generic amplitude batch size is zero"))?;
            let pack_start = Instant::now();
            pack_native2_parameters(
                batch_size,
                global_parameter_count,
                self.input_components.as_deref(),
                state,
                &mut self.parameter_scratch_native2,
            )?;
            let input_pack_s = pack_start.elapsed().as_secs_f64();

            let eval_start = Instant::now();
            self.evaluator.evaluate_native2_into(
                batch_size.div_ceil(2),
                &self.parameter_scratch_native2,
                &mut self.output_scratch_native2,
            )?;
            self.output_scratch_f64
                .resize(batch_size * self.evaluator.output_len, c64(0.0, 0.0));
            for row in 0..batch_size {
                let native_row = row / 2 * self.evaluator.output_len;
                let output_row = row * self.evaluator.output_len;
                let lane = row % 2;
                for column in 0..self.evaluator.output_len {
                    let value = self.output_scratch_native2[native_row + column];
                    self.output_scratch_f64[output_row + column] =
                        c64(value.re.as_array()[lane], value.im.as_array()[lane]);
                }
            }
            return Ok((input_pack_s, eval_start.elapsed().as_secs_f64()));
        }
        if let Some(input_components) = self.input_components.as_ref() {
            let local_parameter_count = input_components.len();
            let global_parameter_count = state
                .len()
                .checked_div(batch_size)
                .ok_or_else(|| PyValueError::new_err("generic amplitude batch size is zero"))?;
            let pack_start = Instant::now();
            self.parameter_scratch_f64
                .resize(batch_size * local_parameter_count, c64(0.0, 0.0));
            for row in 0..batch_size {
                let row_state = row * global_parameter_count;
                let row_params = row * local_parameter_count;
                if self.input_spans.is_empty() {
                    for (local_index, global_index) in input_components.iter().enumerate() {
                        self.parameter_scratch_f64[row_params + local_index] =
                            state[row_state + *global_index];
                    }
                } else {
                    for (local_start, global_start, len) in &self.input_spans {
                        let target_start = row_params + *local_start;
                        let source_start = row_state + *global_start;
                        self.parameter_scratch_f64[target_start..target_start + *len]
                            .copy_from_slice(&state[source_start..source_start + *len]);
                    }
                }
            }
            let input_pack_s = pack_start.elapsed().as_secs_f64();
            let eval_start = Instant::now();
            self.evaluator.evaluate_batch_into(
                batch_size,
                &self.parameter_scratch_f64,
                &mut self.output_scratch_f64,
            )?;
            return Ok((input_pack_s, eval_start.elapsed().as_secs_f64()));
        }
        let eval_start = Instant::now();
        self.evaluator
            .evaluate_batch_into(batch_size, state, &mut self.output_scratch_f64)?;
        Ok((0.0, eval_start.elapsed().as_secs_f64()))
    }

    fn reduce_scratch_f64_into_selected(
        &mut self,
        batch_size: usize,
        raw_sums: &mut Vec<f64>,
        selected_color_sector_ids: Option<&BTreeSet<i64>>,
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
            if selected_color_sector_ids.is_some() {
                return Err(PyValueError::new_err(
                    "LC color-sector runtime selection is only supported for leading-colour diagonal artifacts",
                ));
            }
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
                    if !raw_sum_group_is_selected(group, selected_color_sector_ids) {
                        continue;
                    }
                    let mut sum = c64(0.0, 0.0);
                    for index in &group.indices {
                        sum += amplitudes[row_offset + *index];
                    }
                    let weight = if selected_color_sector_ids.is_none() {
                        group.all_sector_weight
                    } else {
                        group.weight
                    };
                    raw_sums[row] += weight * (sum.re * sum.re + sum.im * sum.im);
                }
                continue;
            }
            for index in 0..self.output_length {
                if !raw_sum_index_is_selected(
                    self.raw_sum_color_sector_ids.get(index).copied().flatten(),
                    selected_color_sector_ids,
                ) {
                    continue;
                }
                let value = amplitudes[row_offset + index];
                let weight = if selected_color_sector_ids.is_none() {
                    self.raw_sum_all_sector_weights[index]
                } else {
                    self.raw_sum_weights[index]
                };
                raw_sums[row] += weight * (value.re * value.re + value.im * value.im);
            }
        }
        Ok(())
    }

    fn reduce_scratch_f64_resolved(
        &mut self,
        batch_size: usize,
        physics: &PhysicsRuntimeV2,
        normalization_factor: f64,
        selected_helicity_ids: Option<&BTreeSet<String>>,
        selected_color_ids: Option<&BTreeSet<String>>,
    ) -> PyResult<ResolvedValues<f64>> {
        let amplitudes = &self.output_scratch_f64;
        if amplitudes.len() != batch_size * self.output_length {
            return Err(PyValueError::new_err(format!(
                "generic amplitude output buffer has length {}, expected {}",
                amplitudes.len(),
                batch_size * self.output_length
            )));
        }
        let helicity_count = physics.manifest.helicities.len();
        let color_count = physics.manifest.color_components.len();
        let mut full_values = vec![0.0; batch_size * helicity_count * color_count];

        if let Some(contraction) = self.color_contraction.as_mut() {
            if color_count != 1 || physics.manifest.color_components[0].kind != "contracted" {
                return Err(PyValueError::new_err(
                    "resolved NLC/full evaluation requires one contracted color component",
                ));
            }
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
                    let left_group = &self.raw_sum_groups[entry.left_group_index];
                    let right_group = &self.raw_sum_groups[entry.right_group_index];
                    let left_reduction = physics
                        .reduction_by_group_id
                        .get(&left_group.id)
                        .ok_or_else(|| {
                            PyValueError::new_err(format!(
                                "resolved metadata is missing coherent group {}",
                                left_group.id
                            ))
                        })?;
                    let right_reduction = physics
                        .reduction_by_group_id
                        .get(&right_group.id)
                        .ok_or_else(|| {
                            PyValueError::new_err(format!(
                                "resolved metadata is missing coherent group {}",
                                right_group.id
                            ))
                        })?;
                    if left_reduction.physical_helicity_ids != right_reduction.physical_helicity_ids
                    {
                        return Err(PyValueError::new_err(
                            "colour contraction mixed distinct physical helicities",
                        ));
                    }
                    let left = contraction.group_scratch_f64[group_row + entry.left_group_index];
                    let right = contraction.group_scratch_f64[group_row + entry.right_group_index];
                    let product = left * right.conj();
                    let contribution = normalization_factor
                        * entry.symmetry_factor
                        * (entry.weight_re * product.re - entry.weight_im * product.im)
                        / left_reduction.physical_helicity_ids.len() as f64;
                    for helicity_id in &left_reduction.physical_helicity_ids {
                        let helicity_index = physics.helicity_index_by_id[helicity_id];
                        full_values[(row * helicity_count + helicity_index) * color_count] +=
                            contribution;
                    }
                }
            }
        } else {
            if !self.has_coherent_groups {
                return Err(PyValueError::new_err(
                    "resolved evaluation requires coherent amplitude-group metadata",
                ));
            }
            for row in 0..batch_size {
                let row_offset = row * self.output_length;
                for group in &self.raw_sum_groups {
                    let reduction =
                        physics
                            .reduction_by_group_id
                            .get(&group.id)
                            .ok_or_else(|| {
                                PyValueError::new_err(format!(
                                    "resolved metadata is missing coherent group {}",
                                    group.id
                                ))
                            })?;
                    let mut sum = c64(0.0, 0.0);
                    for index in &group.indices {
                        sum += amplitudes[row_offset + *index];
                    }
                    let member_count =
                        reduction.physical_helicity_ids.len() * reduction.physical_color_ids.len();
                    let contribution = normalization_factor
                        * group.all_sector_weight
                        * (sum.re * sum.re + sum.im * sum.im)
                        / member_count as f64;
                    for helicity_id in &reduction.physical_helicity_ids {
                        let helicity_index = physics.helicity_index_by_id[helicity_id];
                        for color_id in &reduction.physical_color_ids {
                            let color_index = physics.color_index_by_id[color_id];
                            full_values[(row * helicity_count + helicity_index) * color_count
                                + color_index] += contribution;
                        }
                    }
                }
            }
        }

        let helicity_indices = physics.selected_helicity_indices(selected_helicity_ids)?;
        let color_indices = physics.selected_color_indices(selected_color_ids)?;
        let mut values =
            Vec::with_capacity(batch_size * helicity_indices.len() * color_indices.len());
        for row in 0..batch_size {
            for helicity_index in &helicity_indices {
                for color_index in &color_indices {
                    values.push(
                        full_values
                            [(row * helicity_count + *helicity_index) * color_count + *color_index],
                    );
                }
            }
        }
        Ok(ResolvedValues {
            values,
            point_count: batch_size,
            helicity_indices,
            color_indices,
        })
    }

    fn evaluate_outputs_generic<T>(
        &mut self,
        batch_size: usize,
        state: &[Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<(Vec<Complex<T>>, f64, f64)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let mut input_pack_s = 0.0;
        let (evaluated, evaluator_call_s) =
            if let Some(input_components) = self.input_components.as_ref() {
                let local_parameter_count = input_components.len();
                let global_parameter_count = state
                    .len()
                    .checked_div(batch_size)
                    .ok_or_else(|| PyValueError::new_err("generic amplitude batch size is zero"))?;
                let pack_start = Instant::now();
                let mut parameter_scratch =
                    vec![complex_zero::<T>(); batch_size * local_parameter_count];
                for row in 0..batch_size {
                    let row_state = row * global_parameter_count;
                    let row_params = row * local_parameter_count;
                    if self.input_spans.is_empty() {
                        for (local_index, global_index) in input_components.iter().enumerate() {
                            parameter_scratch[row_params + local_index] =
                                state[row_state + *global_index].clone();
                        }
                    } else {
                        for (local_start, global_start, len) in &self.input_spans {
                            let target_start = row_params + *local_start;
                            let source_start = row_state + *global_start;
                            parameter_scratch[target_start..target_start + *len]
                                .clone_from_slice(&state[source_start..source_start + *len]);
                        }
                    }
                }
                input_pack_s = pack_start.elapsed().as_secs_f64();
                let eval_start = Instant::now();
                let evaluated = self.evaluator.evaluate_batch_generic(
                    batch_size,
                    &parameter_scratch,
                    binary_precision,
                )?;
                (evaluated, eval_start.elapsed().as_secs_f64())
            } else {
                let eval_start = Instant::now();
                let evaluated =
                    self.evaluator
                        .evaluate_batch_generic(batch_size, state, binary_precision)?;
                (evaluated, eval_start.elapsed().as_secs_f64())
            };
        if evaluated.len() != batch_size * self.output_length {
            return Err(PyValueError::new_err(format!(
                "generic amplitude output buffer has length {}, expected {}",
                evaluated.len(),
                batch_size * self.output_length
            )));
        }
        Ok((evaluated, input_pack_s, evaluator_call_s))
    }

    fn evaluate_raw_sums_generic<T>(
        &mut self,
        batch_size: usize,
        state: &[Complex<T>],
        binary_precision: Option<u32>,
    ) -> PyResult<(Vec<T>, f64, f64)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let (evaluated, input_pack_s, evaluator_call_s) =
            self.evaluate_outputs_generic(batch_size, state, binary_precision)?;
        let mut raw_sums = vec![T::new_zero(); batch_size];
        if let Some(contraction) = self.color_contraction.as_ref() {
            if self.raw_sum_groups.len() != contraction.group_count {
                return Err(PyValueError::new_err(
                    "colour contraction group count does not match coherent groups",
                ));
            }
            let mut group_values = vec![complex_zero::<T>(); batch_size * contraction.group_count];
            for row in 0..batch_size {
                let row_offset = row * self.output_length;
                let group_row = row * contraction.group_count;
                for (group_index, group) in self.raw_sum_groups.iter().enumerate() {
                    let mut sum = complex_zero::<T>();
                    for index in &group.indices {
                        sum.re += evaluated[row_offset + *index].re.clone();
                        sum.im += evaluated[row_offset + *index].im.clone();
                    }
                    group_values[group_row + group_index] = sum;
                }
                for entry in &contraction.entries {
                    let left = &group_values[group_row + entry.left_group_index];
                    let right = &group_values[group_row + entry.right_group_index];
                    let product_re =
                        left.re.clone() * right.re.clone() + left.im.clone() * right.im.clone();
                    let product_im =
                        left.im.clone() * right.re.clone() - left.re.clone() * right.im.clone();
                    raw_sums[row] += T::from(entry.symmetry_factor)
                        * (T::from(entry.weight_re) * product_re
                            - T::from(entry.weight_im) * product_im);
                }
            }
            return Ok((raw_sums, input_pack_s, evaluator_call_s));
        }
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
                    raw_sums[row] += T::from(group.all_sector_weight)
                        * (sum_re.clone() * sum_re + sum_im.clone() * sum_im);
                }
                continue;
            }
            for index in 0..self.output_length {
                let value = &evaluated[row_offset + index];
                raw_sums[row] += T::from(self.raw_sum_all_sector_weights[index])
                    * (value.re.clone() * value.re.clone() + value.im.clone() * value.im.clone());
            }
        }
        Ok((raw_sums, input_pack_s, evaluator_call_s))
    }

    fn evaluate_resolved_generic<T>(
        &mut self,
        batch_size: usize,
        state: &[Complex<T>],
        binary_precision: Option<u32>,
        physics: &PhysicsRuntimeV2,
        normalization_factor: f64,
        selected_helicity_ids: Option<&BTreeSet<String>>,
        selected_color_ids: Option<&BTreeSet<String>>,
    ) -> PyResult<(ResolvedValues<T>, f64, f64)>
    where
        T: RusticolHighPrecisionNumber,
        Complex<T>: Real + EvaluationDomain,
    {
        let (evaluated, input_pack_s, evaluator_call_s) =
            self.evaluate_outputs_generic(batch_size, state, binary_precision)?;
        let helicity_count = physics.manifest.helicities.len();
        let color_count = physics.manifest.color_components.len();
        let mut full_values = vec![T::new_zero(); batch_size * helicity_count * color_count];
        if let Some(contraction) = self.color_contraction.as_ref() {
            if color_count != 1 || physics.manifest.color_components[0].kind != "contracted" {
                return Err(PyValueError::new_err(
                    "resolved NLC/full evaluation requires one contracted color component",
                ));
            }
            let mut group_values = vec![complex_zero::<T>(); batch_size * contraction.group_count];
            for row in 0..batch_size {
                let row_offset = row * self.output_length;
                let group_row = row * contraction.group_count;
                for (group_index, group) in self.raw_sum_groups.iter().enumerate() {
                    let mut sum = complex_zero::<T>();
                    for index in &group.indices {
                        sum.re += evaluated[row_offset + *index].re.clone();
                        sum.im += evaluated[row_offset + *index].im.clone();
                    }
                    group_values[group_row + group_index] = sum;
                }
                for entry in &contraction.entries {
                    let left_group = &self.raw_sum_groups[entry.left_group_index];
                    let right_group = &self.raw_sum_groups[entry.right_group_index];
                    let left_reduction = physics
                        .reduction_by_group_id
                        .get(&left_group.id)
                        .ok_or_else(|| {
                            PyValueError::new_err(format!(
                                "resolved metadata is missing coherent group {}",
                                left_group.id
                            ))
                        })?;
                    let right_reduction = physics
                        .reduction_by_group_id
                        .get(&right_group.id)
                        .ok_or_else(|| {
                            PyValueError::new_err(format!(
                                "resolved metadata is missing coherent group {}",
                                right_group.id
                            ))
                        })?;
                    if left_reduction.physical_helicity_ids != right_reduction.physical_helicity_ids
                    {
                        return Err(PyValueError::new_err(
                            "colour contraction mixed distinct physical helicities",
                        ));
                    }
                    let left = &group_values[group_row + entry.left_group_index];
                    let right = &group_values[group_row + entry.right_group_index];
                    let product_re =
                        left.re.clone() * right.re.clone() + left.im.clone() * right.im.clone();
                    let product_im =
                        left.im.clone() * right.re.clone() - left.re.clone() * right.im.clone();
                    let coefficient = normalization_factor * entry.symmetry_factor
                        / left_reduction.physical_helicity_ids.len() as f64;
                    let contribution = T::from(coefficient)
                        * (T::from(entry.weight_re) * product_re
                            - T::from(entry.weight_im) * product_im);
                    for helicity_id in &left_reduction.physical_helicity_ids {
                        let helicity_index = physics.helicity_index_by_id[helicity_id];
                        full_values[(row * helicity_count + helicity_index) * color_count] +=
                            contribution.clone();
                    }
                }
            }
        } else {
            for row in 0..batch_size {
                let row_offset = row * self.output_length;
                for group in &self.raw_sum_groups {
                    let reduction =
                        physics
                            .reduction_by_group_id
                            .get(&group.id)
                            .ok_or_else(|| {
                                PyValueError::new_err(format!(
                                    "resolved metadata is missing coherent group {}",
                                    group.id
                                ))
                            })?;
                    let mut sum_re = T::new_zero();
                    let mut sum_im = T::new_zero();
                    for index in &group.indices {
                        sum_re += evaluated[row_offset + *index].re.clone();
                        sum_im += evaluated[row_offset + *index].im.clone();
                    }
                    let member_count =
                        reduction.physical_helicity_ids.len() * reduction.physical_color_ids.len();
                    let contribution = T::from(
                        normalization_factor * group.all_sector_weight / member_count as f64,
                    ) * (sum_re.clone() * sum_re + sum_im.clone() * sum_im);
                    for helicity_id in &reduction.physical_helicity_ids {
                        let helicity_index = physics.helicity_index_by_id[helicity_id];
                        for color_id in &reduction.physical_color_ids {
                            let color_index = physics.color_index_by_id[color_id];
                            full_values[(row * helicity_count + helicity_index) * color_count
                                + color_index] += contribution.clone();
                        }
                    }
                }
            }
        }
        let helicity_indices = physics.selected_helicity_indices(selected_helicity_ids)?;
        let color_indices = physics.selected_color_indices(selected_color_ids)?;
        let mut values =
            Vec::with_capacity(batch_size * helicity_indices.len() * color_indices.len());
        for row in 0..batch_size {
            for helicity_index in &helicity_indices {
                for color_index in &color_indices {
                    values.push(
                        full_values
                            [(row * helicity_count + *helicity_index) * color_count + *color_index]
                            .clone(),
                    );
                }
            }
        }
        Ok((
            ResolvedValues {
                values,
                point_count: batch_size,
                helicity_indices,
                color_indices,
            },
            input_pack_s,
            evaluator_call_s,
        ))
    }
}

fn build_raw_sum_groups(
    output_length: usize,
    weights: &[f64],
    all_sector_weights: &[f64],
    group_ids: &[Option<i64>],
    color_sector_ids: &[Option<i64>],
) -> PyResult<Vec<RawSumGroup>> {
    if weights.len() != output_length
        || all_sector_weights.len() != output_length
        || group_ids.len() != output_length
        || color_sector_ids.len() != output_length
    {
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
                all_sector_weight: all_sector_weights[index],
                sector_ids: color_sector_ids[index].into_iter().collect(),
            });
        }
    }
    for (group_id, indices) in grouped {
        let weight = weights[indices[0]];
        let all_sector_weight = all_sector_weights[indices[0]];
        if indices
            .iter()
            .any(|index| (weights[*index] - weight).abs() > 0.0)
        {
            return Err(PyValueError::new_err(format!(
                "coherent amplitude group {group_id} has inconsistent raw-sum weights"
            )));
        }
        if indices
            .iter()
            .any(|index| (all_sector_weights[*index] - all_sector_weight).abs() > 0.0)
        {
            return Err(PyValueError::new_err(format!(
                "coherent amplitude group {group_id} has inconsistent all-sector raw-sum weights"
            )));
        }
        groups.push(RawSumGroup {
            id: group_id,
            sector_ids: unique_color_sector_ids(&indices, color_sector_ids),
            indices,
            weight,
            all_sector_weight,
        });
    }
    Ok(groups)
}

fn unique_color_sector_ids(indices: &[usize], color_sector_ids: &[Option<i64>]) -> Vec<i64> {
    indices
        .iter()
        .filter_map(|index| color_sector_ids.get(*index).copied().flatten())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn raw_sum_group_is_selected(
    group: &RawSumGroup,
    selected_color_sector_ids: Option<&BTreeSet<i64>>,
) -> bool {
    let Some(selected) = selected_color_sector_ids else {
        return true;
    };
    group
        .sector_ids
        .iter()
        .any(|sector_id| selected.contains(sector_id))
}

fn raw_sum_index_is_selected(
    sector_id: Option<i64>,
    selected_color_sector_ids: Option<&BTreeSet<i64>>,
) -> bool {
    let Some(selected) = selected_color_sector_ids else {
        return true;
    };
    sector_id
        .map(|value| selected.contains(&value))
        .unwrap_or(false)
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
            chunk_scratch_native2: Vec::new(),
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

    fn supports_native2(&self) -> bool {
        !self.evaluators.is_empty()
            && self
                .evaluators
                .iter()
                .all(|evaluator| matches!(&evaluator.eval, F64Evaluator::JitNative2(_)))
    }

    fn evaluate_native2_into(
        &mut self,
        native_rows: usize,
        params: &[Complex<wide::f64x2>],
        out: &mut Vec<Complex<wide::f64x2>>,
    ) -> PyResult<()> {
        let expected_output_len = native_rows * self.output_len;
        if out.len() != expected_output_len {
            out.resize(
                expected_output_len,
                Complex::new(wide::f64x2::ZERO, wide::f64x2::ZERO),
            );
        }
        if self.evaluators.len() == 1 {
            return self.evaluators[0].evaluate_native2_batch(native_rows, params, out);
        }

        let mut output_offset = 0;
        for evaluator in &mut self.evaluators {
            self.chunk_scratch_native2.resize(
                native_rows * evaluator.output_len,
                Complex::new(wide::f64x2::ZERO, wide::f64x2::ZERO),
            );
            evaluator.evaluate_native2_batch(
                native_rows,
                params,
                &mut self.chunk_scratch_native2,
            )?;
            for row in 0..native_rows {
                let source_start = row * evaluator.output_len;
                let target_start = row * self.output_len + output_offset;
                out[target_start..target_start + evaluator.output_len].copy_from_slice(
                    &self.chunk_scratch_native2[source_start..source_start + evaluator.output_len],
                );
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
                eval.evaluate_batch(batch_size, params, out)
                    .map_err(PyRuntimeError::new_err)
            }
            F64Evaluator::JitNative2(eval) => {
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
                eval.evaluate_batch(batch_size, params, out)
                    .map_err(PyRuntimeError::new_err)
            }
        }
    }

    fn evaluate_native2_batch(
        &mut self,
        native_rows: usize,
        params: &[Complex<wide::f64x2>],
        out: &mut [Complex<wide::f64x2>],
    ) -> PyResult<()> {
        if params.len() != native_rows * self.input_len {
            return Err(PyValueError::new_err(format!(
                "native parameter buffer has length {}, expected {}",
                params.len(),
                native_rows * self.input_len
            )));
        }
        if out.len() != native_rows * self.output_len {
            return Err(PyValueError::new_err(format!(
                "native output buffer has length {}, expected {}",
                out.len(),
                native_rows * self.output_len
            )));
        }
        match &mut self.eval {
            F64Evaluator::JitNative2(eval) => {
                eval.batch_evaluate(params, out, native_rows);
                Ok(())
            }
            _ => Err(PyRuntimeError::new_err(
                "native two-lane evaluation requested for a non-native evaluator",
            )),
        }
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
            let (jit_settings, exact_eval, jit_eval, jit_native2) =
                load_evaluator_state(&artifact_path(root, evaluator_state_path))?;
            let eval = if native_simd_jit_enabled() {
                F64Evaluator::JitNative2(if native_simd_jit_recompile_enabled() {
                    exact_eval
                        .jit_compile::<Complex<wide::f64x2>>(jit_settings.clone())
                        .map_err(|err| {
                            PyRuntimeError::new_err(format!(
                                "could not recompile native two-lane JIT evaluator from {}: {err}",
                                evaluator_state_path
                            ))
                        })?
                } else {
                    match jit_native2 {
                            Some(eval) => eval,
                            None => exact_eval
                                .jit_compile::<Complex<wide::f64x2>>(jit_settings.clone())
                                .map_err(|err| {
                                    PyRuntimeError::new_err(format!(
                                        "could not compile native two-lane JIT evaluator from {}: {err}",
                                        evaluator_state_path
                                    ))
                                })?,
                        }
                })
            } else {
                F64Evaluator::Jit(match jit_eval {
                    Some(eval) => eval,
                    None => exact_eval
                        .jit_compile::<Complex<f64>>(jit_settings)
                        .map_err(|err| {
                            PyRuntimeError::new_err(format!(
                                "could not compile scalar JIT evaluator from {}: {err}",
                                evaluator_state_path
                            ))
                        })?,
                })
            };
            output.push(LoadedEvaluator {
                eval,
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
                    load_evaluator_state(&artifact_path(root, state_path)).map(|state| state.1)
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
    JITCompilationSettings,
    ExpressionEvaluator<Complex<Rational>>,
    Option<JITCompiledEvaluator<Complex<f64>>>,
    Option<JITCompiledEvaluator<Complex<wide::f64x2>>>,
)> {
    type SavedEvaluatorNative2 = (
        bool,
        JITCompilationSettings,
        ExpressionEvaluator<Complex<Rational>>,
        Option<JITCompiledEvaluator<f64>>,
        Option<JITCompiledEvaluator<Complex<f64>>>,
        Option<JITCompiledEvaluator<Complex<wide::f64x2>>>,
    );
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
    match bincode::decode_from_slice::<SavedEvaluatorNative2, _>(
        &bytes,
        bincode::config::standard(),
    ) {
        Ok(((_, settings, evaluator, _, jit_complex, jit_native2), _)) => {
            Ok((settings, evaluator, jit_complex, jit_native2))
        }
        Err(native_err) => match bincode::decode_from_slice::<SavedEvaluator, _>(
            &bytes,
            bincode::config::standard(),
        ) {
            Ok(((_, settings, evaluator, _, jit_complex), _)) => {
                Ok((settings, evaluator, jit_complex, None))
            }
            Err(new_err) => {
                let decoded = bincode::decode_from_slice::<LegacySavedEvaluator, _>(
                    &bytes,
                    bincode::config::standard(),
                )
                .map_err(|_| {
                    PyValueError::new_err(format!(
                        "could not decode evaluator state {}: {native_err}; {new_err}",
                        path.display()
                    ))
                })?;
                let (_, evaluator, _, jit_complex) = decoded.0;
                Ok((
                    JITCompilationSettings::default(),
                    evaluator,
                    jit_complex,
                    None,
                ))
            }
        },
    }
}

fn native_simd_jit_enabled() -> bool {
    if !cfg!(target_arch = "aarch64") {
        return false;
    }
    std::env::var("RUSTICOL_NATIVE_SIMD_JIT")
        .ok()
        .map(|value| {
            !matches!(
                value.to_ascii_lowercase().as_str(),
                "0" | "false" | "no" | "off"
            )
        })
        .unwrap_or(true)
}

fn native_simd_jit_recompile_enabled() -> bool {
    std::env::var("RUSTICOL_RECOMPILE_NATIVE_SIMD_JIT")
        .ok()
        .map(|value| {
            matches!(
                value.to_ascii_lowercase().as_str(),
                "1" | "true" | "yes" | "on"
            )
        })
        .unwrap_or(false)
}

fn artifact_path(root: &Path, value: &str) -> PathBuf {
    let path = PathBuf::from(value);
    if path.is_absolute() {
        path
    } else {
        root.join(path)
    }
}

fn input_crossing_map_to_json(map: &Vec<InputCrossingMapEntry>) -> Vec<(usize, usize, f64)> {
    map.iter()
        .map(|entry| (entry.target_index, entry.source_index, entry.sign))
        .collect()
}

fn apply_input_crossing_map(
    batch: Vec<Vec<[f64; 4]>>,
    expected_legs: usize,
    input_crossing_map: Option<&[InputCrossingMapEntry]>,
) -> PyResult<Vec<Vec<[f64; 4]>>> {
    let Some(map) = input_crossing_map else {
        return Ok(batch);
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
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

#[cfg(feature = "python")]
fn f64_values_to_numpy_or_list(py: Python<'_>, values: Vec<f64>) -> PyResult<Py<PyAny>> {
    Ok(values.into_pyarray(py).into_any().unbind())
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

fn spin2_outer(left: &[Complex<f64>; 4], right: &[Complex<f64>; 4]) -> [Complex<f64>; 16] {
    std::array::from_fn(|index| left[index / 4] * right[index % 4])
}

fn ext_spin2(momentum: [f64; 4], helicity: i32, mass: f64) -> PyResult<[Complex<f64>; 16]> {
    if mass == 0.0 {
        if ![-2, 2].contains(&helicity) {
            return Err(PyValueError::new_err(
                "massless spin-2 sources only support helicities -2 and 2",
            ));
        }
        let vector = ext_gluon(momentum, helicity / 2);
        return Ok(spin2_outer(&vector, &vector));
    }
    if ![-2, -1, 0, 1, 2].contains(&helicity) {
        return Err(PyValueError::new_err(format!(
            "unsupported massive spin-2 helicity {helicity}"
        )));
    }
    let plus = ext_massive_vector(momentum, 1, mass);
    let minus = ext_massive_vector(momentum, -1, mass);
    let longitudinal = ext_massive_vector(momentum, 0, mass);
    let plus_plus = spin2_outer(&plus, &plus);
    let minus_minus = spin2_outer(&minus, &minus);
    if helicity == 2 {
        return Ok(plus_plus);
    }
    if helicity == -2 {
        return Ok(minus_minus);
    }
    let inverse_sqrt_two = 1.0 / 2.0f64.sqrt();
    if helicity == 1 {
        let first = spin2_outer(&plus, &longitudinal);
        let second = spin2_outer(&longitudinal, &plus);
        return Ok(std::array::from_fn(|index| {
            (first[index] + second[index]) * inverse_sqrt_two
        }));
    }
    if helicity == -1 {
        let first = spin2_outer(&minus, &longitudinal);
        let second = spin2_outer(&longitudinal, &minus);
        return Ok(std::array::from_fn(|index| {
            (first[index] + second[index]) * inverse_sqrt_two
        }));
    }
    let plus_minus = spin2_outer(&plus, &minus);
    let minus_plus = spin2_outer(&minus, &plus);
    let zero_zero = spin2_outer(&longitudinal, &longitudinal);
    let inverse_sqrt_six = 1.0 / 6.0f64.sqrt();
    Ok(std::array::from_fn(|index| {
        (plus_minus[index] + minus_plus[index] + c64(2.0, 0.0) * zero_zero[index])
            * inverse_sqrt_six
    }))
}

fn spin2_outer_generic<T>(left: &[Complex<T>; 4], right: &[Complex<T>; 4]) -> [Complex<T>; 16]
where
    T: Real + RealLike + From<f64> + Clone,
{
    std::array::from_fn(|index| left[index / 4].clone() * right[index % 4].clone())
}

fn ext_spin2_generic<T>(momentum: &[T; 4], helicity: i32, mass: T) -> PyResult<[Complex<T>; 16]>
where
    T: Real + RealLike + From<f64> + PartialOrd + Clone,
{
    if is_zero(&mass) {
        if ![-2, 2].contains(&helicity) {
            return Err(PyValueError::new_err(
                "massless spin-2 sources only support helicities -2 and 2",
            ));
        }
        let vector = ext_gluon_generic(momentum, helicity / 2);
        return Ok(spin2_outer_generic(&vector, &vector));
    }
    if ![-2, -1, 0, 1, 2].contains(&helicity) {
        return Err(PyValueError::new_err(format!(
            "unsupported massive spin-2 helicity {helicity}"
        )));
    }
    let plus = ext_massive_vector_generic(momentum, 1, mass.clone());
    let minus = ext_massive_vector_generic(momentum, -1, mass.clone());
    let longitudinal = ext_massive_vector_generic(momentum, 0, mass.clone());
    if helicity == 2 {
        return Ok(spin2_outer_generic(&plus, &plus));
    }
    if helicity == -2 {
        return Ok(spin2_outer_generic(&minus, &minus));
    }
    let inverse_sqrt_two = mass.one() / mass.from_i64(2).sqrt();
    let weight_two = c_generic(inverse_sqrt_two, T::new_zero());
    if helicity == 1 {
        let first = spin2_outer_generic(&plus, &longitudinal);
        let second = spin2_outer_generic(&longitudinal, &plus);
        return Ok(std::array::from_fn(|index| {
            (first[index].clone() + second[index].clone()) * weight_two.clone()
        }));
    }
    if helicity == -1 {
        let first = spin2_outer_generic(&minus, &longitudinal);
        let second = spin2_outer_generic(&longitudinal, &minus);
        return Ok(std::array::from_fn(|index| {
            (first[index].clone() + second[index].clone()) * weight_two.clone()
        }));
    }
    let plus_minus = spin2_outer_generic(&plus, &minus);
    let minus_plus = spin2_outer_generic(&minus, &plus);
    let zero_zero = spin2_outer_generic(&longitudinal, &longitudinal);
    let two = c_generic(mass.from_i64(2), T::new_zero());
    let inverse_sqrt_six = c_generic(mass.one() / mass.from_i64(6).sqrt(), T::new_zero());
    Ok(std::array::from_fn(|index| {
        (plus_minus[index].clone()
            + minus_plus[index].clone()
            + two.clone() * zero_zero[index].clone())
            * inverse_sqrt_six.clone()
    }))
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

    fn test_physics_runtime(color_accuracy: &str) -> PhysicsRuntimeV2 {
        let contracted = color_accuracy != "lc";
        let color_components = if contracted {
            vec![PhysicsColorComponentManifestV2 {
                id: "contracted".to_string(),
                index: 0,
                kind: "contracted".to_string(),
                word: Vec::new(),
                representative_id: "contracted".to_string(),
                computed: true,
                internal_sector_id: None,
            }]
        } else {
            vec![
                PhysicsColorComponentManifestV2 {
                    id: "flow:0".to_string(),
                    index: 0,
                    kind: "lc-flow".to_string(),
                    word: vec![1, 2],
                    representative_id: "flow:0".to_string(),
                    computed: true,
                    internal_sector_id: Some(0),
                },
                PhysicsColorComponentManifestV2 {
                    id: "flow:1".to_string(),
                    index: 1,
                    kind: "lc-flow".to_string(),
                    word: vec![2, 1],
                    representative_id: "flow:0".to_string(),
                    computed: false,
                    internal_sector_id: Some(1),
                },
            ]
        };
        let physical_color_ids = color_components
            .iter()
            .map(|item| item.id.clone())
            .collect();
        let mut helicities = vec![
            PhysicsHelicityManifestV2 {
                id: "hel:+-".to_string(),
                index: 0,
                helicities: vec![1, -1],
                representative_id: "hel:+-".to_string(),
                computed: true,
                structural_zero: false,
            },
            PhysicsHelicityManifestV2 {
                id: "hel:-+".to_string(),
                index: 1,
                helicities: vec![-1, 1],
                representative_id: "hel:+-".to_string(),
                computed: false,
                structural_zero: false,
            },
        ];
        if !contracted {
            helicities.push(PhysicsHelicityManifestV2 {
                id: "hel:zero".to_string(),
                index: 2,
                helicities: vec![1, 1],
                representative_id: "hel:zero".to_string(),
                computed: false,
                structural_zero: true,
            });
        }
        PhysicsRuntimeV2::new(PhysicsManifestV2 {
            schema_version: 1,
            kind: "pyamplicol-resolved-physics".to_string(),
            process: "x x > y y".to_string(),
            process_key: "x_x_to_y_y".to_string(),
            color_accuracy: color_accuracy.to_string(),
            external_particles: vec![
                PhysicsExternalParticleManifestV2 {
                    label: 1,
                    index: 0,
                    side: "incoming".to_string(),
                    role: "initial".to_string(),
                    particle: "x".to_string(),
                    outgoing_particle: "x~".to_string(),
                    pdg: 1,
                    outgoing_pdg: -1,
                    particle_class: "fermion".to_string(),
                    momentum_slot: 0,
                },
                PhysicsExternalParticleManifestV2 {
                    label: 2,
                    index: 1,
                    side: "incoming".to_string(),
                    role: "initial".to_string(),
                    particle: "x~".to_string(),
                    outgoing_particle: "x".to_string(),
                    pdg: -1,
                    outgoing_pdg: 1,
                    particle_class: "fermion".to_string(),
                    momentum_slot: 1,
                },
            ],
            helicities,
            color_components,
            model_parameters: Vec::new(),
            coverage: PhysicsCoverageManifestV2 {
                helicities: "complete".to_string(),
                color: "complete".to_string(),
                color_kind: if contracted {
                    "contracted".to_string()
                } else {
                    "physical-lc-flows".to_string()
                },
                structural_zero_helicity_count: usize::from(!contracted),
            },
            selectors: PhysicsSelectorsManifestV2 {
                helicity: true,
                color_flow: !contracted,
                contracted_color: contracted,
            },
            reduction: PhysicsReductionManifestV2 {
                kind: "coherent-group-physical-expansion".to_string(),
                groups: vec![PhysicsReductionGroupManifestV2 {
                    group_id: 7,
                    representative_helicity_id: "hel:+-".to_string(),
                    physical_helicity_ids: vec!["hel:+-".to_string(), "hel:-+".to_string()],
                    representative_color_id: if contracted {
                        "contracted".to_string()
                    } else {
                        "flow:0".to_string()
                    },
                    physical_color_ids,
                }],
            },
        })
        .unwrap()
    }

    fn empty_evaluator_group() -> EvaluatorGroup {
        EvaluatorGroup {
            evaluators: Vec::new(),
            output_len: 0,
            chunk_scratch_f64: Vec::new(),
            chunk_scratch_native2: Vec::new(),
        }
    }

    fn test_amplitude_runtime(
        outputs: Vec<Complex<f64>>,
        color_contraction: Option<ColorContractionRuntime>,
    ) -> GenericAmplitudeRuntimeV2 {
        let output_length = outputs.len();
        GenericAmplitudeRuntimeV2 {
            output_length,
            raw_sum_weights: vec![1.0; output_length],
            raw_sum_all_sector_weights: vec![1.0; output_length],
            raw_sum_color_sector_ids: vec![None; output_length],
            raw_sum_groups: vec![RawSumGroup {
                id: 7,
                indices: (0..output_length).collect(),
                weight: 1.0,
                all_sector_weight: 1.0,
                sector_ids: vec![0],
            }],
            has_coherent_groups: true,
            color_contraction,
            input_components: None,
            input_spans: Vec::new(),
            parameter_scratch_f64: Vec::new(),
            output_scratch_f64: outputs,
            parameter_scratch_native2: Vec::new(),
            output_scratch_native2: Vec::new(),
            evaluator: empty_evaluator_group(),
        }
    }

    #[test]
    fn resolved_lc_reduction_expands_symmetries_and_structural_zeros() {
        let physics = test_physics_runtime("lc");
        let mut amplitude = test_amplitude_runtime(vec![c64(2.0, 0.0)], None);

        let resolved = amplitude
            .reduce_scratch_f64_resolved(1, &physics, 4.0, None, None)
            .unwrap();

        assert_eq!(resolved.point_count, 1);
        assert_eq!(resolved.helicity_indices, vec![0, 1, 2]);
        assert_eq!(resolved.color_indices, vec![0, 1]);
        assert_eq!(resolved.values, vec![4.0, 4.0, 4.0, 4.0, 0.0, 0.0]);
        assert_eq!(resolved.values.iter().sum::<f64>(), 16.0);

        let helicities = BTreeSet::from(["hel:-+".to_string()]);
        let colors = BTreeSet::from(["flow:1".to_string()]);
        let selected = amplitude
            .reduce_scratch_f64_resolved(1, &physics, 4.0, Some(&helicities), Some(&colors))
            .unwrap();
        assert_eq!(selected.values, vec![4.0]);
    }

    #[test]
    fn resolved_nlc_and_full_reductions_have_one_contracted_color_component() {
        for color_accuracy in ["nlc", "full"] {
            let physics = test_physics_runtime(color_accuracy);
            let contraction = ColorContractionRuntime {
                group_count: 1,
                entries: vec![ColorContractionEntry {
                    left_group_index: 0,
                    right_group_index: 0,
                    weight_re: 2.0,
                    weight_im: 0.0,
                    symmetry_factor: 1.0,
                }],
                group_scratch_f64: Vec::new(),
            };
            let mut amplitude = test_amplitude_runtime(vec![c64(3.0, 0.0)], Some(contraction));

            let resolved = amplitude
                .reduce_scratch_f64_resolved(1, &physics, 2.0, None, None)
                .unwrap();

            assert_eq!(resolved.helicity_indices, vec![0, 1]);
            assert_eq!(resolved.color_indices, vec![0]);
            assert_eq!(resolved.values, vec![18.0, 18.0]);
            assert_eq!(resolved.values.iter().sum::<f64>(), 36.0);
        }
    }

    #[test]
    fn model_parameter_override_batch_is_atomic() {
        let manifest: GenericProcessManifestV2 =
            serde_json::from_value(minimal_generic_manifest()).unwrap();
        let mut runtime = GenericRuntimeV2::from_manifest(manifest).unwrap();
        runtime.model_parameter_values_f64 = vec![0.118];
        runtime.model_parameter_runtime_slots.insert(
            "normalization.alpha_s_me_check".to_string(),
            GenericRuntimeParameterSlots {
                real: 0,
                imaginary: None,
            },
        );
        let invalid_batch = BTreeMap::from([
            ("normalization.alpha_s_me_check".to_string(), (0.101, 0.0)),
            ("unknown.parameter".to_string(), (1.0, 0.0)),
        ]);

        let error = runtime
            .apply_model_parameter_overrides(&invalid_batch)
            .unwrap_err();

        assert!(error.to_string().contains("unknown.parameter"));
        assert_eq!(runtime.model_parameter_values_f64, vec![0.118]);
        runtime
            .apply_model_parameter_overrides(&BTreeMap::from([(
                "normalization.alpha_s_me_check".to_string(),
                (0.101, 0.0),
            )]))
            .unwrap();
        assert_eq!(runtime.model_parameter_values_f64, vec![0.101]);
    }

    #[test]
    fn resolved_warning_state_is_mutable_and_once_per_native_handle() {
        let mut physics = test_physics_runtime("lc");
        physics.manifest.coverage.helicities = "fixed".to_string();
        physics.manifest.coverage.color = "fixed".to_string();
        let manifest: GenericProcessManifestV2 =
            serde_json::from_value(minimal_generic_manifest()).unwrap();
        let mut generic_runtime = GenericRuntimeV2::from_manifest(manifest).unwrap();
        generic_runtime.physics = Some(physics.clone());
        let mut runtime = NativeRuntime {
            root: PathBuf::new(),
            runtime: generic_runtime,
            process: "x x > y y".to_string(),
            process_key: "x_x_to_y_y".to_string(),
            input_crossing_map: None,
            crossing_alias_of: None,
            physics: Some(physics.manifest.clone()),
            warnings_muted: false,
            warned_kinds: BTreeSet::new(),
            pending_warnings: Vec::new(),
        };
        let helicities = ["hel:-+".to_string()];
        let colors = ["flow:1".to_string()];

        runtime.mute_warnings();
        runtime
            .record_resolved_warnings(Some(&helicities), Some(&colors))
            .unwrap();
        assert!(runtime.take_warnings().is_empty());

        runtime.unmute_warnings();
        runtime
            .record_resolved_warnings(Some(&helicities), Some(&colors))
            .unwrap();
        runtime
            .record_resolved_warnings(Some(&helicities), Some(&colors))
            .unwrap();
        let warnings = runtime.take_warnings();
        assert_eq!(warnings.len(), 3);
        assert!(
            warnings
                .iter()
                .any(|warning| warning.contains("helicities"))
        );
        assert!(warnings.iter().any(|warning| warning.contains("color")));
        assert!(warnings.iter().any(|warning| warning.contains("symmetry")));
        runtime
            .record_resolved_warnings(Some(&helicities), Some(&colors))
            .unwrap();
        assert!(runtime.take_warnings().is_empty());
    }

    #[test]
    fn native_two_lane_parameter_pack_preserves_component_order_and_odd_tail() {
        let state = (0..12)
            .map(|index| c64(index as f64, -(index as f64)))
            .collect::<Vec<_>>();
        let mut packed = Vec::new();

        pack_native2_parameters(3, 4, Some(&[3, 1]), &state, &mut packed).unwrap();

        assert_eq!(packed.len(), 4);
        assert_eq!(packed[0].re.to_array(), [3.0, 7.0]);
        assert_eq!(packed[0].im.to_array(), [-3.0, -7.0]);
        assert_eq!(packed[1].re.to_array(), [1.0, 5.0]);
        assert_eq!(packed[2].re.to_array(), [11.0, 11.0]);
        assert_eq!(packed[3].re.to_array(), [9.0, 9.0]);
    }

    #[test]
    fn derived_model_parameter_components_are_not_user_overrides() {
        let parameters =
            serde_json::from_value::<Vec<GenericRuntimeModelParameterManifestV2>>(json!([
                {
                    "name": "aS.real",
                    "kind": "external_parameter_component",
                    "parameter_index": 0,
                    "runtime_name": "aS",
                    "complex_component": "real"
                },
                {
                    "name": "aS.imag",
                    "kind": "external_parameter_component",
                    "parameter_index": 1,
                    "runtime_name": "aS",
                    "complex_component": "imag"
                },
                {
                    "name": "derived_coupling_75.real",
                    "kind": "derived_parameter_component",
                    "parameter_index": 2,
                    "runtime_name": "derived_coupling_75",
                    "complex_component": "real"
                },
                {
                    "name": "derived_coupling_75.imag",
                    "kind": "derived_parameter_component",
                    "parameter_index": 3,
                    "runtime_name": "derived_coupling_75",
                    "complex_component": "imag"
                }
            ]))
            .unwrap();

        let slots = build_runtime_parameter_slots(&parameters).unwrap();

        assert_eq!(slots.len(), 1);
        assert_eq!(slots["aS"].real, 0);
        assert_eq!(slots["aS"].imaginary, Some(1));
        assert!(!slots.contains_key("derived_coupling_75"));
    }

    #[test]
    fn complex_model_parameter_json_uses_real_imaginary_pairs() {
        let path = Path::new("model-parameters.json");
        let parsed = parse_complex_parameter_overrides(
            r#"{"aS":[0.118,0.0],"complex_mass":[173.0,-1.5]}"#,
            path,
        )
        .unwrap();

        assert_eq!(parsed["aS"], (0.118, 0.0));
        assert_eq!(parsed["complex_mass"], (173.0, -1.5));
    }

    #[test]
    #[cfg(feature = "python")]
    fn generic_dynamic_momenta_parser_accepts_more_than_sixteen_legs() {
        let point = vec![vec![0.0, 1.0, 2.0, 3.0]; 17];
        let parsed = batch_momenta_dynamic_from_nested(vec![point], 17).unwrap();

        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].len(), 17);
        assert_eq!(parsed[0][16], [0.0, 1.0, 2.0, 3.0]);
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

        let mapped = apply_input_crossing_map(batch.clone(), 3, Some(&map)).unwrap();

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

        let error = apply_input_crossing_map(batch, 2, Some(&duplicate)).unwrap_err();

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
    fn generic_schema_v2_validator_accepts_compact_storage_metadata() {
        let mut payload = minimal_generic_manifest();
        payload["runtime_schema"]["current_storage"]["metadata_compacted"] = json!(true);
        payload["runtime_schema"]["value_storage"]["metadata_compacted"] = json!(true);
        for slot in payload["runtime_schema"]["current_storage"]["current_slots"]
            .as_array_mut()
            .unwrap()
        {
            let object = slot.as_object_mut().unwrap();
            object.remove("external_labels");
            object.remove("helicity_ancestry");
            object.remove("spin_state");
            object.remove("flavour_flow");
            object.remove("charge_flow");
            object.remove("color_state");
            object.remove("auxiliary_kind");
        }
        for slot in payload["runtime_schema"]["value_storage"]["value_slots"]
            .as_array_mut()
            .unwrap()
        {
            slot.as_object_mut().unwrap().remove("external_labels");
        }
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        validate_generic_schema_v2_manifest(&manifest).unwrap();
    }

    #[test]
    fn generic_input_value_variant_index_preserves_precedence() {
        let manifest: GenericProcessManifestV2 =
            serde_json::from_value(minimal_generic_manifest()).unwrap();
        let current_storage = &manifest.runtime_schema.current_storage;
        let mut value_storage = manifest.runtime_schema.value_storage.clone();

        let variants = build_input_value_variants(current_storage, &value_storage).unwrap();
        assert_eq!(input_value_variant(&variants, 0).unwrap(), "source");
        assert_eq!(input_value_variant(&variants, 2).unwrap(), "unpropagated");

        let mut propagated = value_storage.value_slots[2].clone();
        propagated.variant = "propagated".to_string();
        value_storage.value_slots.push(propagated);
        let variants = build_input_value_variants(current_storage, &value_storage).unwrap();
        assert_eq!(input_value_variant(&variants, 2).unwrap(), "propagated");
    }

    #[test]
    fn generic_input_value_variant_index_reports_missing_slots() {
        let manifest: GenericProcessManifestV2 =
            serde_json::from_value(minimal_generic_manifest()).unwrap();
        let current_storage = &manifest.runtime_schema.current_storage;
        let mut value_storage = manifest.runtime_schema.value_storage.clone();
        value_storage
            .value_slots
            .retain(|slot| slot.current_id != 2);

        let variants = build_input_value_variants(current_storage, &value_storage).unwrap();
        let error = input_value_variant(&variants, 2).unwrap_err();
        assert!(error.to_string().contains("without an input value slot"));
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
        assert_eq!(runtime.lc_topology_replay_weights, vec![1.0, 1.0]);
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
    fn generic_schema_v2_validator_accepts_compact_stage_interactions() {
        let mut payload = minimal_generic_manifest();
        add_minimal_stage_evaluators(&mut payload);
        payload["runtime_schema"]["stages"][0]["interactions_compacted"] = json!(true);
        payload["runtime_schema"]["stages"][0]["interaction_ids"] = json!([0]);
        payload["runtime_schema"]["stages"][0]["interactions"] = json!([]);
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        validate_generic_schema_v2_manifest(&manifest).unwrap();
    }

    #[test]
    fn generic_schema_v2_validator_rejects_bad_compact_stage_interactions() {
        let mut payload = minimal_generic_manifest();
        payload["runtime_schema"]["stages"][0]["interactions_compacted"] = json!(true);
        payload["runtime_schema"]["stages"][0]["interaction_ids"] = json!([]);
        payload["runtime_schema"]["stages"][0]["interactions"] = json!([]);
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        let error = validate_generic_schema_v2_manifest(&manifest)
            .expect_err("bad compact interaction metadata should be rejected");
        assert!(error.to_string().contains("compact stage 1"));
    }

    #[test]
    fn generic_schema_v2_validator_rejects_out_of_range_compact_interaction() {
        let mut payload = minimal_generic_manifest();
        payload["runtime_schema"]["stages"][0]["interactions_compacted"] = json!(true);
        payload["runtime_schema"]["stages"][0]["interaction_ids"] = json!([1]);
        payload["runtime_schema"]["stages"][0]["interactions"] = json!([]);
        let manifest: GenericProcessManifestV2 = serde_json::from_value(payload).unwrap();

        let error = validate_generic_schema_v2_manifest(&manifest)
            .expect_err("out-of-range compact interaction should be rejected");
        assert!(error.to_string().contains("invalid interaction 1"));
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

#[cfg(feature = "python")]
#[pymodule]
fn rusticol(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Runtime>()?;
    m.add_class::<ProcessPhysics>()?;
    m.add_class::<ExternalParticle>()?;
    m.add_class::<HelicityConfiguration>()?;
    m.add_class::<ColorFlow>()?;
    m.add_class::<ContractedColorComponent>()?;
    m.add_class::<ModelParameter>()?;
    m.add_class::<ResolvedEvaluation>()?;
    m.add_function(wrap_pyfunction!(build_profile, m)?)?;
    m.add_function(wrap_pyfunction!(build_target, m)?)?;
    Ok(())
}

#[cfg(feature = "python")]
#[pyfunction]
fn build_profile() -> &'static str {
    env!("RUSTICOL_BUILD_PROFILE")
}

#[cfg(feature = "python")]
#[pyfunction]
fn build_target() -> &'static str {
    env!("RUSTICOL_BUILD_TARGET")
}
