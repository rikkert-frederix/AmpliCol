from __future__ import annotations

from importlib import metadata

from .logging import (
    DEFAULT_LOG_FORMAT,
    LOGGER_NAME,
    configure_logging,
    disable_logging,
    get_logger,
    progress,
)
from .evaluation import (
    NativeRuntimeEvaluator,
    NativeRuntimeMetadata,
    RuntimeBackend,
)
from .lowering import (
    ColorAlgebraProbe,
    RecursionLoweringPlan,
    SymbolicLoweringReport,
    TensorNetworkProbe,
    TensorNetworkBlueprint,
    VertexLoweringReport,
    VertexLoweringStep,
    build_symbolic_lowering_report,
)
from .matrix import (
    CurrentKey,
    EvaluatorArtifact,
    GenerationResult,
    InteractionNode,
    NativeMatrixElementGenerator,
    RecursionGraph,
    evaluator_artifact_path,
    load_evaluator_artifact,
)
from .model import (
    AmplicolSMLeadingColorModel,
    Model,
    Particle,
    Vertex,
    VertexLoweringRule,
)
from .native import (
    ExternalMomentum,
    HelicityContribution,
    LeadingColorZJetsNativeEvaluator,
    MatrixElementEvaluation,
    NativeEvaluationError,
)
from .params import (
    ParamBuilder,
    ParameterRange,
    SymbolicaEvaluatorBundle,
)
from .process_runtime import (
    ProcessArtifactManifest,
    ProcessRuntimeBackend,
    PythonProcessRuntime,
    load_process,
    load_process_manifest,
)
from .processes import (
    ParsedProcess,
    PhaseSpaceGroup,
    ProcessEnumeration,
    ProcessEnumerator,
    ProcessOptions,
    SubprocessRecord,
    enumerate_processes,
    write_legacy_process_file,
)
from .reference import (
    AmplicolAdapter,
    AmplicolFirstPoint,
    AmplicolWorkflowResult,
    CommandResult,
    TimingRow,
    parse_first_phase_space_point,
    parse_first_matrix_element,
    parse_timing_rows,
    popen_runner,
)
from .symbolic import (
    SymbolicEvaluatorMetadata,
    ZeroGluonSymbolicEvaluator,
    build_zero_gluon_symbolic_evaluator_payload,
)
from .tensor_runtime import (
    NumericTensorNetworkRuntimeMetadata,
    TensorNetworkStrategy,
    ZGluonNumericTensorNetworkEvaluator,
    ZGluonTensorNetworkEvaluator,
)


def _package_version() -> str:
    try:
        return metadata.version("pyamplicol")
    except metadata.PackageNotFoundError:
        return "0.0.1"


__version__ = _package_version()

__all__ = [
    "AmplicolAdapter",
    "AmplicolFirstPoint",
    "AmplicolSMLeadingColorModel",
    "AmplicolWorkflowResult",
    "CommandResult",
    "ColorAlgebraProbe",
    "CurrentKey",
    "DEFAULT_LOG_FORMAT",
    "ExternalMomentum",
    "EvaluatorArtifact",
    "GenerationResult",
    "HelicityContribution",
    "InteractionNode",
    "LOGGER_NAME",
    "LeadingColorZJetsNativeEvaluator",
    "MatrixElementEvaluation",
    "Model",
    "NativeEvaluationError",
    "NativeMatrixElementGenerator",
    "NativeRuntimeEvaluator",
    "NativeRuntimeMetadata",
    "NumericTensorNetworkRuntimeMetadata",
    "ParamBuilder",
    "ParsedProcess",
    "ParameterRange",
    "Particle",
    "PhaseSpaceGroup",
    "ProcessArtifactManifest",
    "ProcessEnumeration",
    "ProcessEnumerator",
    "ProcessOptions",
    "ProcessRuntimeBackend",
    "PythonProcessRuntime",
    "RecursionGraph",
    "RecursionLoweringPlan",
    "RuntimeBackend",
    "SubprocessRecord",
    "SymbolicaEvaluatorBundle",
    "SymbolicEvaluatorMetadata",
    "SymbolicLoweringReport",
    "TensorNetworkBlueprint",
    "TensorNetworkProbe",
    "TensorNetworkStrategy",
    "TimingRow",
    "Vertex",
    "VertexLoweringReport",
    "VertexLoweringRule",
    "VertexLoweringStep",
    "__version__",
    "build_zero_gluon_symbolic_evaluator_payload",
    "configure_logging",
    "build_symbolic_lowering_report",
    "disable_logging",
    "enumerate_processes",
    "evaluator_artifact_path",
    "get_logger",
    "load_evaluator_artifact",
    "load_process",
    "load_process_manifest",
    "parse_first_phase_space_point",
    "parse_first_matrix_element",
    "parse_timing_rows",
    "popen_runner",
    "progress",
    "write_legacy_process_file",
    "ZeroGluonSymbolicEvaluator",
    "ZGluonNumericTensorNetworkEvaluator",
    "ZGluonTensorNetworkEvaluator",
]
