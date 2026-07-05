use pyo3::IntoPyObjectExt;
use pyo3::buffer::PyBuffer;
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyList};
use serde::Deserialize;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;
use symbolica::evaluate::JITCompiledEvaluator;
use symbolica::prelude::{
    BatchEvaluator, CompiledComplexEvaluator, Complex, DoubleFloat, EvaluationDomain,
    ExpressionEvaluator, Float, JITCompilationSettings, Rational, Real, RealLike,
};

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
struct ModelManifest {
    alpha_s_me_check: f64,
    alpha_ew: f64,
    mass_z: f64,
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
}

#[derive(Clone, Debug, Deserialize)]
struct AmplitudeManifest {
    multiplicity: f64,
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
    amplitude_evaluator: Option<EvaluatorManifest>,
    raw_sum_evaluator: Option<EvaluatorManifest>,
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
    evaluator: EvaluatorGroup,
}

struct AmplitudeStage {
    output_length: usize,
    raw_sum_weights: Vec<f64>,
    amplitude_evaluator: Option<EvaluatorGroup>,
    raw_sum_evaluator: Option<EvaluatorGroup>,
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

#[pyclass(module = "rusticol")]
struct Runtime {
    root: PathBuf,
    manifest: ProcessManifest,
    stages: Option<Vec<CurrentStage>>,
    amplitude_stage: Option<AmplitudeStage>,
    zero_gluon_stage: Option<ZeroGluonStage>,
    last_profile: RuntimeProfile,
}

#[pymethods]
impl Runtime {
    #[classmethod]
    fn load(_cls: &Bound<'_, pyo3::types::PyType>, process_dir: &str) -> PyResult<Self> {
        let process_dir = PathBuf::from(process_dir);
        let root = process_dir.canonicalize().map_err(|err| {
            PyValueError::new_err(format!(
                "could not resolve process directory {}: {err}",
                process_dir.display()
            ))
        })?;
        let manifest_path = root.join("process_manifest.json");
        let manifest_text = fs::read_to_string(&manifest_path).map_err(|err| {
            PyValueError::new_err(format!(
                "could not read process manifest {}: {err}",
                manifest_path.display()
            ))
        })?;
        let manifest: ProcessManifest = serde_json::from_str(&manifest_text).map_err(|err| {
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
        if manifest.family != "q-qbar-z-gluons-leading-color" {
            return Err(PyValueError::new_err(format!(
                "rusticol currently supports q-qbar-z-gluons-leading-color artifacts, got {}",
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
                manifest,
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
        if manifest.external_pdg_order.len() != manifest.gluon_count + 3 {
            return Err(PyValueError::new_err(
                "external PDG order length does not match q qbar -> Z + n gluons",
            ));
        }
        if manifest
            .table
            .currents
            .iter()
            .enumerate()
            .any(|(index, current)| {
                current.id != index || current.dimension > 6 || current.pdg == 0
            })
        {
            return Err(PyValueError::new_err(
                "current table ids/dimensions are not compatible with rusticol",
            ));
        }
        if manifest.table.amplitudes.len()
            != manifest.compiled.amplitude_stage.raw_sum_weights.len()
        {
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
            manifest,
            stages: Some(stages),
            amplitude_stage: Some(amplitude_stage),
            zero_gluon_stage: None,
            last_profile: RuntimeProfile::default(),
        })
    }

    #[getter]
    fn process(&self) -> &str {
        &self.manifest.process
    }

    #[getter]
    fn gluon_count(&self) -> usize {
        self.manifest.gluon_count
    }

    fn metadata<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        dict.set_item("process", &self.manifest.process)?;
        dict.set_item("family", &self.manifest.family)?;
        dict.set_item("gluon_count", self.manifest.gluon_count)?;
        dict.set_item("parameter_count", self.manifest.layout.parameter_count)?;
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

    fn evaluate_with_prec<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
        decimal_digit_precision: u32,
    ) -> PyResult<Py<PyAny>> {
        match PrecisionMode::from_decimal_digits(decimal_digit_precision)? {
            PrecisionMode::F64 => self.evaluate(py, momenta),
            PrecisionMode::DoubleDouble => {
                let batch = batch_momenta_double(momenta, self.manifest.external_pdg_order.len())?;
                let (values, profile) = self.run_double(&batch)?;
                self.last_profile = profile;
                decimals_to_python(py, values, decimal_digit_precision)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let batch = batch_momenta_float(
                    momenta,
                    binary_precision,
                    self.manifest.external_pdg_order.len(),
                )?;
                let (values, profile) = self.run_float(&batch, binary_precision)?;
                self.last_profile = profile;
                decimals_to_python(py, values, decimal_precision)
            }
        }
    }

    #[pyo3(signature = (momenta, precision = 16))]
    fn profile<'py>(
        &mut self,
        py: Python<'py>,
        momenta: &Bound<'py, PyAny>,
        precision: u32,
    ) -> PyResult<Bound<'py, PyDict>> {
        let (points, values, profile) = match PrecisionMode::from_decimal_digits(precision)? {
            PrecisionMode::F64 => {
                let batch = batch_momenta(py, momenta, self.manifest.external_pdg_order.len())?;
                let points = batch.len();
                let (values, profile) = self.run_f64(&batch)?;
                (points, values.into_py_any(py)?, profile)
            }
            PrecisionMode::DoubleDouble => {
                let batch = batch_momenta_double(momenta, self.manifest.external_pdg_order.len())?;
                let points = batch.len();
                let (values, profile) = self.run_double(&batch)?;
                (points, decimals_to_python(py, values, precision)?, profile)
            }
            PrecisionMode::Arbitrary(decimal_precision) => {
                let binary_precision = decimal_digits_to_bits(decimal_precision);
                let batch = batch_momenta_float(
                    momenta,
                    binary_precision,
                    self.manifest.external_pdg_order.len(),
                )?;
                let points = batch.len();
                let (values, profile) = self.run_float(&batch, binary_precision)?;
                (
                    points,
                    decimals_to_python(py, values, decimal_precision)?,
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
        let batch = batch_momenta(py, momenta, self.manifest.external_pdg_order.len())?;
        let n_points = batch.len();
        let parameter_count = self.manifest.layout.parameter_count;
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
                self.manifest.table.currents.len() * 6,
            ),
        )?)?;

        let current_offsets = self.manifest.layout.current_offsets.clone();
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
        let batch = batch_momenta(py, momenta, self.manifest.external_pdg_order.len())?;
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
        let parameter_count = self.manifest.layout.parameter_count;
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
        let mut evaluated = Vec::new();
        for stage in self.stages_mut() {
            let eval_start = Instant::now();
            stage
                .evaluator
                .evaluate_batch_into(n_points, &state, &mut evaluated)?;
            stage_evaluator_s += eval_start.elapsed().as_secs_f64();
            let assign_start = Instant::now();
            for row in 0..n_points {
                let row_state = &mut state[row * parameter_count..(row + 1) * parameter_count];
                for (column, _, _, state_offset) in &stage.outputs {
                    row_state[*state_offset] =
                        evaluated[row * stage.evaluator.output_len + *column];
                }
            }
            output_assign_s += assign_start.elapsed().as_secs_f64();
        }

        let amp_start = Instant::now();
        let raw_sums = self
            .amplitude_stage_mut()
            .evaluate_raw_sums(n_points, &state)?;
        let amplitude_evaluator_s = amp_start.elapsed().as_secs_f64();

        let reduction_start = Instant::now();
        let factor = self.manifest.normalization.color_factor
            * self.manifest.normalization.coupling_factor
            / (self.manifest.normalization.average_factor
                * self.manifest.normalization.identical_factor);
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
        let parameter_count = self.manifest.layout.parameter_count;
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
            self.manifest.normalization.color_factor * self.manifest.normalization.coupling_factor
                / (self.manifest.normalization.average_factor
                    * self.manifest.normalization.identical_factor),
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
            .parameter_rows_f64(&self.manifest, batch)?;
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
        let factor = self.manifest.normalization.color_factor
            * self.manifest.normalization.coupling_factor
            / (self.manifest.normalization.average_factor
                * self.manifest.normalization.identical_factor);
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
            .parameter_rows_generic(&self.manifest, batch)?;
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
            self.manifest.normalization.color_factor * self.manifest.normalization.coupling_factor
                / (self.manifest.normalization.average_factor
                    * self.manifest.normalization.identical_factor),
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

    fn fill_sources(&self, row: &mut [Complex<f64>], point: &[[f64; 4]; 16]) -> PyResult<()> {
        for source in &self.manifest.table.sources {
            let offset = self.manifest.layout.current_offsets[source.current_id];
            for component in 0..6 {
                row[offset + component] = c64(0.0, 0.0);
            }
            let wave = match source.leg_label {
                1 => {
                    let p = negate(point[0]);
                    if source.physical_helicity == 1 {
                        ext_antiquark_weyl(p, 1, -1)
                    } else {
                        ext_antiquark_weyl(p, -1, 1)
                    }
                }
                2 => {
                    let p = negate(point[1]);
                    if source.chirality == 1 {
                        ext_quark_weyl(p, -1, 1)
                    } else {
                        ext_quark_weyl(p, 1, -1)
                    }
                }
                label if label >= 3 && label < self.manifest.gluon_count + 3 => {
                    ext_gluon(point[label - 1], source.helicity).to_vec()
                }
                label if label == self.manifest.gluon_count + 3 => ext_massive_vector(
                    point[label - 1],
                    source.helicity,
                    self.manifest.model.mass_z,
                )
                .to_vec(),
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

    fn fill_sources_generic<T>(&self, row: &mut [Complex<T>], point: &[[T; 4]]) -> PyResult<()>
    where
        T: Real + RealLike + From<f64> + PartialOrd + Clone,
    {
        for source in &self.manifest.table.sources {
            let offset = self.manifest.layout.current_offsets[source.current_id];
            for component in 0..6 {
                row[offset + component] = complex_zero::<T>();
            }
            let wave = match source.leg_label {
                1 => {
                    let p = negate_generic(&point[0]);
                    if source.physical_helicity == 1 {
                        ext_antiquark_weyl_generic(&p, 1, -1)
                    } else {
                        ext_antiquark_weyl_generic(&p, -1, 1)
                    }
                }
                2 => {
                    let p = negate_generic(&point[1]);
                    if source.chirality == 1 {
                        ext_quark_weyl_generic(&p, -1, 1)
                    } else {
                        ext_quark_weyl_generic(&p, 1, -1)
                    }
                }
                label if label >= 3 && label < self.manifest.gluon_count + 3 => {
                    ext_gluon_generic(&point[label - 1], source.helicity).to_vec()
                }
                label if label == self.manifest.gluon_count + 3 => ext_massive_vector_generic(
                    &point[label - 1],
                    source.helicity,
                    T::from(self.manifest.model.mass_z),
                )
                .to_vec(),
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
        for entry in &self.manifest.layout.momentum_offsets_and_labels {
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
        for entry in &self.manifest.layout.momentum_offsets_and_labels {
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
        Ok(Self {
            outputs,
            evaluator: EvaluatorGroup::load(&manifest.evaluator, root)?,
        })
    }
}

impl AmplitudeStage {
    fn load(manifest: &AmplitudeStageManifest, root: &Path) -> PyResult<Self> {
        Ok(Self {
            output_length: manifest.output_length,
            raw_sum_weights: manifest.raw_sum_weights.clone(),
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
        let evaluated = amplitude_evaluator.evaluate_batch(batch_size, state)?;
        let mut raw_sums = vec![0.0; batch_size];
        for row in 0..batch_size {
            let row_offset = row * amplitude_evaluator.output_len;
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
        out.clear();
        out.resize(batch_size * self.output_len, c64(0.0, 0.0));
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
            F64Evaluator::Jit(eval) => eval
                .evaluate_batch(batch_size, params, out)
                .map_err(PyRuntimeError::new_err),
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

fn ext_quark_weyl(momentum: [f64; 4], helicity: i32, chirality: i32) -> Vec<Complex<f64>> {
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
            return vec![chi1, chi2];
        }
        if helicity == -1 && chirality == -1 {
            return vec![chi2, chi1];
        }
        return vec![c64(0.0, 0.0), c64(0.0, 0.0)];
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
        vec![chi1, chi2]
    } else if helicity == 1 && chirality == -1 {
        vec![chi2, chi1]
    } else {
        vec![c64(0.0, 0.0), c64(0.0, 0.0)]
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

fn ext_antiquark_weyl(momentum: [f64; 4], helicity: i32, chirality: i32) -> Vec<Complex<f64>> {
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
            return vec![chi2, chi1];
        }
        if helicity == -1 && chirality == -1 {
            return vec![chi1, chi2];
        }
        return vec![c64(0.0, 0.0), c64(0.0, 0.0)];
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
        vec![chi2, chi1]
    } else if helicity == 1 && chirality == -1 {
        vec![chi1, chi2]
    } else {
        vec![c64(0.0, 0.0), c64(0.0, 0.0)]
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

#[pymodule]
fn rusticol(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<Runtime>()?;
    Ok(())
}
