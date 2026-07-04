from __future__ import annotations

import re
import time
from collections import Counter
from dataclasses import dataclass
from typing import Any

from .model import AmplicolSMLeadingColorModel
from .params import ParamBuilder, SymbolicaEvaluatorBundle
from .symbols import symbols

_MAX_RECURSION_EXPRESSION_PREVIEW = 4096


@dataclass(frozen=True)
class TensorNetworkProbe:
    engine: str
    tensor_names: tuple[str, ...]
    expression: str
    output_structure: str
    output_rank: int
    output_size: int
    nonzero_entries: int
    max_abs_entry: float
    weighted_checksum: tuple[float, float]
    first_nonzero_entries: tuple[tuple[int, str], ...]


@dataclass(frozen=True)
class ColorAlgebraProbe:
    engine: str
    input_expression: str
    simplified_expression: str


@dataclass(frozen=True)
class RecursionLoweringPlan:
    engine: str
    current_count: int
    interaction_count: int
    amplitude_count: int
    source_current_count: int
    vertex_kind_counts: tuple[tuple[int, int], ...]
    tensor_route_vertex_kinds: tuple[int, ...]
    color_order: tuple[int, ...]
    expression: str
    expression_length: int
    expression_truncated: bool
    first_assignments: tuple[str, ...]


@dataclass(frozen=True)
class VertexLoweringStep:
    index: int
    vertex_kind: int
    backend: str
    tensor_names: tuple[str, ...]
    expression_head: str
    full_tensor_network_ready: bool
    result_current: str
    left_current: str
    right_current: str


@dataclass(frozen=True)
class VertexLoweringReport:
    total_interactions: int
    full_tensor_network_ready: bool
    backend_counts: tuple[tuple[str, int], ...]
    ready_vertex_kind_counts: tuple[tuple[int, int], ...]
    pending_vertex_kind_counts: tuple[tuple[int, int], ...]
    tensor_names: tuple[str, ...]
    first_steps: tuple[VertexLoweringStep, ...]


@dataclass(frozen=True)
class TensorNetworkBlueprint:
    engine: str
    status: str
    current_count: int
    interaction_count: int
    amplitude_count: int
    expression_built: bool
    expression_executed: bool
    full_me_tensor_network_ready: bool
    propagator_lowering_ready: bool
    ready_interactions: int
    pending_interactions: int
    placeholder_vertex_kinds: tuple[int, ...]
    registered_tensor_names: tuple[str, ...]
    current_leaf_count: int
    parametric_external_current_count: int
    parametric_source_current_parameter_count: int
    parametric_current_momentum_count: int
    parametric_momentum_parameter_count: int
    parametric_parameter_count: int
    expression: str | None
    expression_length: int | None
    expression_truncated: bool | None
    executed_expression: str | None
    executed_expression_length: int | None
    executed_expression_truncated: bool | None
    execution_time_s: float | None


@dataclass(frozen=True)
class SymbolicLoweringReport:
    tensor_library: str
    tensor_network_probe: TensorNetworkProbe
    color_algebra_probe: ColorAlgebraProbe
    recursion_plan: RecursionLoweringPlan | None = None
    vertex_lowering: VertexLoweringReport | None = None
    tensor_network_blueprint: TensorNetworkBlueprint | None = None
    full_me_tensor_network_ready: bool = False


def build_symbolic_lowering_report(
    model: AmplicolSMLeadingColorModel,
    graph: Any | None = None,
) -> SymbolicLoweringReport:
    """Exercise the real Symbolica/spenso/idenso hooks used by ME lowering.

    This is intentionally a small, deterministic probe: it validates that the
    model-owned auxiliary four-gluon tensors are registered in spenso and that
    idenso color simplification is available. It is not the full ME evaluator.
    """

    tensor_probe = _build_auxiliary_tensor_probe(model)
    color_probe = _build_color_probe()
    recursion_plan = None if graph is None else _build_recursion_plan(graph)
    vertex_lowering = None if graph is None else _build_vertex_lowering_report(model, graph)
    tensor_network_blueprint = (
        None if graph is None else _build_tensor_network_blueprint(model, graph)
    )
    full_me_tensor_network_ready = (
        tensor_network_blueprint.full_me_tensor_network_ready
        if tensor_network_blueprint is not None
        else False
    )
    return SymbolicLoweringReport(
        tensor_library="TensorLibrary.hep_lib_atom",
        tensor_network_probe=tensor_probe,
        color_algebra_probe=color_probe,
        recursion_plan=recursion_plan,
        vertex_lowering=vertex_lowering,
        tensor_network_blueprint=tensor_network_blueprint,
        full_me_tensor_network_ready=full_me_tensor_network_ready,
    )


def _build_recursion_plan(graph: Any) -> RecursionLoweringPlan:
    assignments = tuple(
        _interaction_assignment(interaction) for interaction in graph.interactions
    )
    amplitudes = tuple(
        symbols.amplitude(
            _current_atom(left),
            _current_atom(right),
        )
        for left, right in graph.amplitudes
    )
    assignment_sum = _sum_expressions(assignments)
    amplitude_sum = _sum_expressions(amplitudes)
    plan_expression = symbols.matrix_element_plan(assignment_sum, amplitude_sum)
    vertex_kind_counts = Counter(
        int(interaction.vertex_kind) for interaction in graph.interactions
    )
    expression = _clean_symbolica_string(str(plan_expression))
    tensor_route_vertex_kinds = tuple(
        kind for kind in (1, 2, 3) if vertex_kind_counts.get(kind, 0) > 0
    )
    first_assignments = tuple(
        _clean_symbolica_string(str(assignment)) for assignment in assignments[:5]
    )
    return RecursionLoweringPlan(
        engine="symbolica",
        current_count=len(graph.currents),
        interaction_count=len(graph.interactions),
        amplitude_count=len(graph.amplitudes),
        source_current_count=_source_current_count(graph),
        vertex_kind_counts=tuple(sorted(vertex_kind_counts.items())),
        tensor_route_vertex_kinds=tensor_route_vertex_kinds,
        color_order=tuple(int(index) for index in graph.color_order),
        expression=_preview_expression(expression),
        expression_length=len(expression),
        expression_truncated=len(expression) > _MAX_RECURSION_EXPRESSION_PREVIEW,
        first_assignments=first_assignments,
    )


def _interaction_assignment(interaction: Any) -> Any:
    return symbols.assignment(
        _current_atom(interaction.result),
        symbols.vertex(
            int(interaction.vertex_kind),
            _current_atom(interaction.left),
            _current_atom(interaction.right),
            _number(interaction.coupling[0]),
            _number(interaction.coupling[1]),
        ),
    )


def _current_atom(current: Any) -> Any:
    return symbols.current(
        int(current.pdg),
        _label_atom(tuple(int(label) for label in current.external_labels)),
        int(current.chirality),
    )


def _label_atom(labels: tuple[int, ...]) -> Any:
    from symbolica import S

    return S("L" + "_".join(str(label) for label in labels))


def _number(value: float) -> Any:
    from symbolica import Expression

    return Expression.num(value)


def _sum_expressions(expressions: tuple[Any, ...]) -> Any:
    from symbolica import Expression

    total = Expression.num(0)
    for expression in expressions:
        total = total + expression
    return total


def _current_key_tuple(current: Any) -> tuple[int, tuple[int, ...], int]:
    return (
        int(current.pdg),
        tuple(int(label) for label in current.external_labels),
        int(current.chirality),
    )


def _build_vertex_lowering_report(
    model: AmplicolSMLeadingColorModel,
    graph: Any,
) -> VertexLoweringReport:
    steps = tuple(
        _vertex_lowering_step(model, index, interaction)
        for index, interaction in enumerate(graph.interactions, start=1)
    )
    backend_counts = Counter(step.backend for step in steps)
    ready_kind_counts = Counter(
        step.vertex_kind for step in steps if step.full_tensor_network_ready
    )
    pending_kind_counts = Counter(
        step.vertex_kind for step in steps if not step.full_tensor_network_ready
    )
    tensor_names = tuple(
        sorted({tensor_name for step in steps for tensor_name in step.tensor_names})
    )
    return VertexLoweringReport(
        total_interactions=len(steps),
        full_tensor_network_ready=bool(steps)
        and all(step.full_tensor_network_ready for step in steps),
        backend_counts=tuple(sorted(backend_counts.items())),
        ready_vertex_kind_counts=tuple(sorted(ready_kind_counts.items())),
        pending_vertex_kind_counts=tuple(sorted(pending_kind_counts.items())),
        tensor_names=tensor_names,
        first_steps=steps[:8],
    )


def _vertex_lowering_step(
    model: AmplicolSMLeadingColorModel,
    index: int,
    interaction: Any,
) -> VertexLoweringStep:
    rule = model.vertex_lowering_rule(int(interaction.vertex_kind))
    return VertexLoweringStep(
        index=index,
        vertex_kind=int(interaction.vertex_kind),
        backend=rule.backend,
        tensor_names=rule.tensor_names,
        expression_head=rule.expression_head,
        full_tensor_network_ready=rule.full_tensor_network_ready,
        result_current=_clean_symbolica_string(str(_current_atom(interaction.result))),
        left_current=_clean_symbolica_string(str(_current_atom(interaction.left))),
        right_current=_clean_symbolica_string(str(_current_atom(interaction.right))),
    )


def _build_tensor_network_blueprint(
    model: AmplicolSMLeadingColorModel,
    graph: Any,
) -> TensorNetworkBlueprint:
    from symbolica.community.idenso import list_dangling, simplify_color
    from symbolica.community.spenso import TensorNetwork

    max_interactions_to_build = 96
    max_interactions_to_execute = 40
    vertex_lowering = _build_vertex_lowering_report(model, graph)
    ready_interactions = sum(
        count for _, count in vertex_lowering.ready_vertex_kind_counts
    )
    pending_interactions = sum(
        count for _, count in vertex_lowering.pending_vertex_kind_counts
    )
    placeholder_vertex_kinds = tuple(
        kind for kind, _ in vertex_lowering.pending_vertex_kind_counts
    )
    registered_tensor_names = vertex_lowering.tensor_names
    source_currents = _source_currents(graph)
    momentum_currents = _current_momentum_currents(graph)
    current_leaf_count = len(source_currents)
    source_parameter_count = sum(
        _current_dimension(current) for current in source_currents
    )
    momentum_parameter_count = 4 * len(momentum_currents)
    parametric_parameter_count = source_parameter_count + momentum_parameter_count
    propagator_ready = _propagator_lowering_ready(graph)
    full_ready = vertex_lowering.full_tensor_network_ready and propagator_ready

    if len(graph.interactions) > max_interactions_to_build:
        return TensorNetworkBlueprint(
            engine="spenso",
            status="size-guarded",
            current_count=len(graph.currents),
            interaction_count=len(graph.interactions),
            amplitude_count=len(graph.amplitudes),
            expression_built=False,
            expression_executed=False,
            full_me_tensor_network_ready=False,
            propagator_lowering_ready=propagator_ready,
            ready_interactions=ready_interactions,
            pending_interactions=pending_interactions,
            placeholder_vertex_kinds=placeholder_vertex_kinds,
            registered_tensor_names=registered_tensor_names,
            current_leaf_count=current_leaf_count,
            parametric_external_current_count=current_leaf_count,
            parametric_source_current_parameter_count=source_parameter_count,
            parametric_current_momentum_count=len(momentum_currents),
            parametric_momentum_parameter_count=momentum_parameter_count,
            parametric_parameter_count=parametric_parameter_count,
            expression=None,
            expression_length=None,
            expression_truncated=None,
            executed_expression=None,
            executed_expression_length=None,
            executed_expression_truncated=None,
            execution_time_s=None,
        )

    try:
        builder = _GraphTensorExpressionBuilder(model, graph)
        raw_expression = builder.matrix_element_skeleton()
        expression = simplify_color(raw_expression)
        expression_text = _clean_symbolica_string(str(expression))
        if len(graph.interactions) > max_interactions_to_execute:
            return TensorNetworkBlueprint(
                engine="spenso",
                status="execution-size-guarded",
                current_count=len(graph.currents),
                interaction_count=len(graph.interactions),
                amplitude_count=len(graph.amplitudes),
                expression_built=True,
                expression_executed=False,
                full_me_tensor_network_ready=False,
                propagator_lowering_ready=propagator_ready,
                ready_interactions=ready_interactions,
                pending_interactions=pending_interactions,
                placeholder_vertex_kinds=placeholder_vertex_kinds,
                registered_tensor_names=registered_tensor_names,
                current_leaf_count=current_leaf_count,
                parametric_external_current_count=current_leaf_count,
                parametric_source_current_parameter_count=source_parameter_count,
                parametric_current_momentum_count=len(momentum_currents),
                parametric_momentum_parameter_count=momentum_parameter_count,
                parametric_parameter_count=parametric_parameter_count,
                expression=_preview_expression(expression_text),
                expression_length=len(expression_text),
                expression_truncated=(
                    len(expression_text) > _MAX_RECURSION_EXPRESSION_PREVIEW
                ),
                executed_expression=None,
                executed_expression_length=None,
                executed_expression_truncated=None,
                execution_time_s=None,
            )

        start = time.perf_counter()
        library = model.build_tensor_library()
        _register_parametric_source_currents(library, graph)
        _register_parametric_current_momenta(model, library, graph)
        network = TensorNetwork(expression, library)
        network.execute(library=library)
        scalar = network.result_scalar()
        execution_time_s = time.perf_counter() - start
        executed_text = _clean_symbolica_string(str(scalar))
        dangling = list_dangling(scalar)
        status = "scalar-skeleton" if not dangling else "dangling-indices"
    except (RuntimeError, TypeError, ValueError) as exc:
        return TensorNetworkBlueprint(
            engine="spenso",
            status=f"failed: {exc}",
            current_count=len(graph.currents),
            interaction_count=len(graph.interactions),
            amplitude_count=len(graph.amplitudes),
            expression_built=False,
            expression_executed=False,
            full_me_tensor_network_ready=False,
            propagator_lowering_ready=propagator_ready,
            ready_interactions=ready_interactions,
            pending_interactions=pending_interactions,
            placeholder_vertex_kinds=placeholder_vertex_kinds,
            registered_tensor_names=registered_tensor_names,
            current_leaf_count=current_leaf_count,
            parametric_external_current_count=current_leaf_count,
            parametric_source_current_parameter_count=source_parameter_count,
            parametric_current_momentum_count=len(momentum_currents),
            parametric_momentum_parameter_count=momentum_parameter_count,
            parametric_parameter_count=parametric_parameter_count,
            expression=None,
            expression_length=None,
            expression_truncated=None,
            executed_expression=None,
            executed_expression_length=None,
            executed_expression_truncated=None,
            execution_time_s=None,
        )

    return TensorNetworkBlueprint(
        engine="spenso",
        status=status,
        current_count=len(graph.currents),
        interaction_count=len(graph.interactions),
        amplitude_count=len(graph.amplitudes),
        expression_built=True,
        expression_executed=True,
        full_me_tensor_network_ready=full_ready,
        propagator_lowering_ready=propagator_ready,
        ready_interactions=ready_interactions,
        pending_interactions=pending_interactions,
        placeholder_vertex_kinds=placeholder_vertex_kinds,
        registered_tensor_names=registered_tensor_names,
        current_leaf_count=current_leaf_count,
        parametric_external_current_count=current_leaf_count,
        parametric_source_current_parameter_count=source_parameter_count,
        parametric_current_momentum_count=len(momentum_currents),
        parametric_momentum_parameter_count=momentum_parameter_count,
        parametric_parameter_count=parametric_parameter_count,
        expression=_preview_expression(expression_text),
        expression_length=len(expression_text),
        expression_truncated=len(expression_text) > _MAX_RECURSION_EXPRESSION_PREVIEW,
        executed_expression=_preview_expression(executed_text),
        executed_expression_length=len(executed_text),
        executed_expression_truncated=(
            len(executed_text) > _MAX_RECURSION_EXPRESSION_PREVIEW
        ),
        execution_time_s=execution_time_s,
    )


def build_tensor_network_scalar_bundle(
    model: AmplicolSMLeadingColorModel,
    graph: Any,
    *,
    name: str,
    collect_expression_metadata: bool = False,
) -> SymbolicaEvaluatorBundle:
    """Build a Symbolica evaluator bundle from the propagated tensor network."""

    from symbolica.community.idenso import simplify_color
    from symbolica.community.spenso import TensorNetwork

    total_start = time.perf_counter()
    library = model.build_tensor_library()
    param_builder = ParamBuilder()
    _register_parametric_source_currents(library, graph, param_builder)
    _register_parametric_current_momenta(model, library, graph, param_builder)
    raw_expression = _GraphTensorExpressionBuilder(model, graph).matrix_element_skeleton()
    expression = simplify_color(raw_expression)
    network = TensorNetwork(expression, library)
    reduction_start = time.perf_counter()
    network.execute(library=library)
    scalar = network.result_scalar()
    reduction_time_s = time.perf_counter() - reduction_start
    evaluator_start = time.perf_counter()
    evaluator = scalar.evaluator(param_builder.parameter_symbols())
    if param_builder.real_valued_inputs:
        evaluator.set_real_params(param_builder.real_valued_inputs)
    evaluator_build_time_s = time.perf_counter() - evaluator_start
    total_build_time_s = time.perf_counter() - total_start
    expression_length = None
    executed_expression_length = None
    if collect_expression_metadata:
        expression_length = len(_clean_symbolica_string(str(expression)))
        executed_expression_length = len(_clean_symbolica_string(str(scalar)))
    return SymbolicaEvaluatorBundle(
        name=name,
        evaluator=evaluator,
        param_builder=param_builder,
        complex_inputs=True,
        metadata={
            "engine": "symbolica",
            "kernel": "tensor-network-scalar",
            "current_count": len(graph.currents),
            "interaction_count": len(graph.interactions),
            "amplitude_count": len(graph.amplitudes),
            "expression_length": expression_length,
            "executed_expression_length": executed_expression_length,
            "expression_metadata_collected": collect_expression_metadata,
            "tensor_network_reduction_s": reduction_time_s,
            "symbolica_evaluator_build_s": evaluator_build_time_s,
            "total_build_s": total_build_time_s,
        },
    )


def build_interleaved_tensor_network_scalar_bundle(
    model: AmplicolSMLeadingColorModel,
    graph: Any,
    *,
    name: str,
    collect_expression_metadata: bool = False,
) -> SymbolicaEvaluatorBundle:
    """Build a Symbolica evaluator by executing one aggregate TensorNetwork in steps."""

    total_start = time.perf_counter()
    library = model.build_tensor_library()
    param_builder = ParamBuilder()
    _register_parametric_source_currents(library, graph, param_builder)
    _register_parametric_current_momenta(model, library, graph, param_builder)

    build_start = time.perf_counter()
    builder = _GraphTensorExpressionBuilder(model, graph)
    network, interleaved_metadata = builder.matrix_element_interleaved_network(
        library,
        execute_between=True,
    )
    network_build_s = time.perf_counter() - build_start

    final_start = time.perf_counter()
    network.execute(library=library)
    scalar = network.result_scalar()
    final_reduction_s = time.perf_counter() - final_start

    evaluator_start = time.perf_counter()
    evaluator = scalar.evaluator(param_builder.parameter_symbols())
    if param_builder.real_valued_inputs:
        evaluator.set_real_params(param_builder.real_valued_inputs)
    evaluator_build_time_s = time.perf_counter() - evaluator_start
    total_build_time_s = time.perf_counter() - total_start
    return SymbolicaEvaluatorBundle(
        name=name,
        evaluator=evaluator,
        param_builder=param_builder,
        complex_inputs=True,
        metadata={
            "engine": "symbolica",
            "kernel": "interleaved-tensor-network-scalar",
            "strategy": "interleaved",
            "current_count": len(graph.currents),
            "interaction_count": len(graph.interactions),
            "amplitude_count": len(graph.amplitudes),
            "network_build_s": network_build_s,
            "interleaved_execution_s": interleaved_metadata["execution_s"],
            "interleaved_execution_count": interleaved_metadata["execution_count"],
            "interleaved_multiplication_count": interleaved_metadata["multiplication_count"],
            "interleaved_addition_count": interleaved_metadata["addition_count"],
            "tensor_network_reduction_s": final_reduction_s,
            "symbolica_evaluator_build_s": evaluator_build_time_s,
            "total_build_s": total_build_time_s,
            "executed_expression_length": (
                len(_clean_symbolica_string(str(scalar)))
                if collect_expression_metadata
                else None
            ),
            "expression_metadata_collected": collect_expression_metadata,
        },
    )


class _GraphTensorExpressionBuilder:
    def __init__(self, model: AmplicolSMLeadingColorModel, graph: Any) -> None:
        from symbolica.community.spenso import Representation, TensorName

        self.model = model
        self.graph = graph
        self._mink = Representation.mink(4)
        self._aux6 = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
        self._weyl = Representation("pyamplicol::weyl_spinor", 2)
        self._tensor_name = TensorName
        self._interactions_by_result: dict[tuple[int, tuple[int, ...], int], list[Any]] = {}
        for interaction in graph.interactions:
            self._interactions_by_result.setdefault(
                _current_key_tuple(interaction.result),
                [],
            ).append(interaction)

    def matrix_element_skeleton(self) -> Any:
        from symbolica import Expression

        total = Expression.num(0)
        for amplitude_index, (left, right) in enumerate(self.graph.amplitudes, start=1):
            slots = self._slots_for_current(left, f"amp_{amplitude_index}")
            total = total + (
                self._current_expression(left, slots)
                * self._current_expression(right, slots)
            )
        return total

    def matrix_element_interleaved_network(
        self,
        library: Any,
        *,
        execute_between: bool = True,
    ) -> tuple[Any, dict[str, float | int]]:
        from symbolica.community.spenso import TensorNetwork

        metadata: dict[str, float | int] = {
            "execution_s": 0.0,
            "execution_count": 0,
            "multiplication_count": 0,
            "addition_count": 0,
        }
        total = None
        for amplitude_index, (left, right) in enumerate(self.graph.amplitudes, start=1):
            slots = self._slots_for_current(left, f"amp_{amplitude_index}")
            left_network = self._current_network(
                left,
                slots,
                library,
                execute_between=execute_between,
                metadata=metadata,
            )
            right_network = self._current_network(
                right,
                slots,
                library,
                execute_between=execute_between,
                metadata=metadata,
            )
            amplitude_network = self._multiply_networks(
                left_network,
                right_network,
                library,
                execute_between=execute_between,
                metadata=metadata,
            )
            total = (
                amplitude_network
                if total is None
                else self._add_networks(
                    total,
                    amplitude_network,
                    library,
                    execute_between=execute_between,
                    metadata=metadata,
                )
            )
        if total is None:
            total = TensorNetwork.zero()
        return total, metadata

    def _current_network(
        self,
        current: Any,
        output_slots: tuple[Any, ...],
        library: Any,
        *,
        execute_between: bool,
        metadata: dict[str, float | int],
    ) -> Any:
        from symbolica.community.spenso import TensorNetwork

        interactions = self._interactions_by_result.get(_current_key_tuple(current))
        if not interactions:
            return TensorNetwork(self._current_leaf(current, output_slots), library)

        needs_propagator = _current_needs_propagator(self.graph, current)
        result_slots = (
            self._slots_for_current(
                current,
                self._propagator_dummy_prefix(current),
            )
            if needs_propagator
            else output_slots
        )
        total = None
        for interaction_index, interaction in enumerate(interactions, start=1):
            left_slots = self._slots_for_current(
                interaction.left,
                self._dummy_prefix(interaction, interaction_index, "left"),
            )
            right_slots = self._slots_for_current(
                interaction.right,
                self._dummy_prefix(interaction, interaction_index, "right"),
            )
            term = TensorNetwork(
                self._vertex_tensor(interaction, left_slots, right_slots, result_slots),
                library,
            )
            term = self._multiply_networks(
                term,
                self._current_network(
                    interaction.left,
                    left_slots,
                    library,
                    execute_between=execute_between,
                    metadata=metadata,
                ),
                library,
                execute_between=execute_between,
                metadata=metadata,
            )
            term = self._multiply_networks(
                term,
                self._current_network(
                    interaction.right,
                    right_slots,
                    library,
                    execute_between=execute_between,
                    metadata=metadata,
                ),
                library,
                execute_between=execute_between,
                metadata=metadata,
            )
            if needs_propagator:
                term = self._multiply_networks(
                    TensorNetwork(
                        self._propagator_tensor(current, result_slots, output_slots),
                        library,
                    ),
                    term,
                    library,
                    execute_between=execute_between,
                    metadata=metadata,
                )
            total = (
                term
                if total is None
                else self._add_networks(
                    total,
                    term,
                    library,
                    execute_between=execute_between,
                    metadata=metadata,
                )
            )
        if total is None:
            return TensorNetwork.zero()
        return total

    def _multiply_networks(
        self,
        left: Any,
        right: Any,
        library: Any,
        *,
        execute_between: bool,
        metadata: dict[str, float | int],
    ) -> Any:
        result = left * right
        metadata["multiplication_count"] = int(metadata["multiplication_count"]) + 1
        if execute_between:
            result = self._execute_network(result, library, metadata)
        return result

    def _add_networks(
        self,
        left: Any,
        right: Any,
        library: Any,
        *,
        execute_between: bool,
        metadata: dict[str, float | int],
    ) -> Any:
        result = left + right
        metadata["addition_count"] = int(metadata["addition_count"]) + 1
        if execute_between:
            result = self._execute_network(result, library, metadata)
        return result

    def _execute_network(
        self,
        network: Any,
        library: Any,
        metadata: dict[str, float | int],
    ) -> Any:
        from symbolica.community.spenso import TensorNetwork

        start = time.perf_counter()
        network.execute(library=library)
        metadata["execution_s"] = float(metadata["execution_s"]) + (
            time.perf_counter() - start
        )
        metadata["execution_count"] = int(metadata["execution_count"]) + 1
        return TensorNetwork.one() * network.result_tensor(library)

    def _current_expression(self, current: Any, output_slots: tuple[Any, ...]) -> Any:
        from symbolica import Expression

        key = _current_key_tuple(current)
        interactions = self._interactions_by_result.get(key)
        if not interactions:
            return self._current_leaf(current, output_slots)

        needs_propagator = _current_needs_propagator(self.graph, current)
        result_slots = (
            self._slots_for_current(
                current,
                self._propagator_dummy_prefix(current),
            )
            if needs_propagator
            else output_slots
        )
        total = Expression.num(0)
        for interaction_index, interaction in enumerate(interactions, start=1):
            left_slots = self._slots_for_current(
                interaction.left,
                self._dummy_prefix(interaction, interaction_index, "left"),
            )
            right_slots = self._slots_for_current(
                interaction.right,
                self._dummy_prefix(interaction, interaction_index, "right"),
            )
            total = total + (
                self._vertex_tensor(interaction, left_slots, right_slots, result_slots)
                * self._current_expression(interaction.left, left_slots)
                * self._current_expression(interaction.right, right_slots)
            )
        if needs_propagator:
            total = (
                self._propagator_tensor(current, result_slots, output_slots)
                * total
            )
        return total

    def _current_leaf(self, current: Any, output_slots: tuple[Any, ...]) -> Any:
        return self._tensor_name(_current_tensor_name(current))(*output_slots).to_expression()

    def _vertex_tensor(
        self,
        interaction: Any,
        left_slots: tuple[Any, ...],
        right_slots: tuple[Any, ...],
        output_slots: tuple[Any, ...],
    ) -> Any:
        kind = int(interaction.vertex_kind)
        if kind == 0:
            return self.model.three_gluon_current_expression(
                left_slot=left_slots[0],
                right_slot=right_slots[0],
                output_slot=output_slots[0],
                left_momentum_tensor_name=_current_momentum_tensor_name(
                    interaction.left
                ),
                right_momentum_tensor_name=_current_momentum_tensor_name(
                    interaction.right
                ),
                dummy_prefix=self._dummy_prefix(interaction, 0, "three_gluon"),
            )
        if kind == 1:
            return self._tensor_name(str(symbols.two_gluon_to_tensor))(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 2:
            return self._tensor_name(str(symbols.tensor_gluon_to_gluon))(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 3:
            return self._tensor_name(str(symbols.gluon_tensor_to_gluon))(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 6:
            return self._tensor_name(
                _quark_vector_weyl_tensor_name(int(interaction.result.chirality))
            )(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        if kind == 10:
            return _number(
                _weyl_coupling_for_chirality(
                    int(interaction.result.chirality),
                    _coupling_pair(interaction.coupling),
                )
            ) * self._tensor_name(
                _quark_vector_weyl_tensor_name(int(interaction.result.chirality))
            )(
                left_slots[0],
                right_slots[0],
                output_slots[0],
            ).to_expression()
        return self._tensor_name(f"pyamplicol::vertex_kind_{kind}")(
            *left_slots,
            *right_slots,
            *output_slots,
            _number(float(interaction.coupling[0])),
            _number(float(interaction.coupling[1])),
        ).to_expression()

    def _propagator_tensor(
        self,
        current: Any,
        input_slots: tuple[Any, ...],
        output_slots: tuple[Any, ...],
    ) -> Any:
        return self._tensor_name(_propagator_tensor_name(current))(
            input_slots[0],
            output_slots[0],
        ).to_expression()

    def _slots_for_current(self, current: Any, prefix: str) -> tuple[Any, ...]:
        pdg = int(current.pdg)
        if pdg == -21:
            return (self._aux6(f"{prefix}_A"),)
        if pdg == 21 or pdg == 22 or pdg == 23:
            return (self._mink(f"{prefix}_mu"),)
        if 1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16:
            return (self._weyl(f"{prefix}_alpha"),)
        return ()

    def _dummy_prefix(self, interaction: Any, index: int, side: str) -> str:
        labels = "_".join(str(label) for label in interaction.result.external_labels)
        chirality = int(interaction.result.chirality)
        return (
            f"v{int(interaction.vertex_kind)}_{labels}_"
            f"{_signed_label(chirality)}_{index}_{side}"
        )

    def _propagator_dummy_prefix(self, current: Any) -> str:
        labels = "_".join(str(label) for label in current.external_labels)
        return (
            f"prop_{_signed_label(int(current.pdg))}_{labels}_"
            f"{_signed_label(int(current.chirality))}"
        )

    def _current_output_prefix(self, current: Any) -> str:
        labels = "_".join(str(label) for label in current.external_labels)
        return (
            f"cur_{_signed_label(int(current.pdg))}_{labels}_"
            f"{_signed_label(int(current.chirality))}"
        )


def _current_tensor_name(current: Any) -> str:
    labels = "_".join(str(label) for label in current.external_labels)
    return (
        "pyamplicol::"
        f"current_{_signed_label(int(current.pdg))}_{labels}_"
        f"{_signed_label(int(current.chirality))}"
    )


def _signed_label(value: int) -> str:
    if value < 0:
        return f"m{abs(value)}"
    return f"p{value}"


def _quark_vector_weyl_tensor_name(chirality: int) -> str:
    if chirality == 1:
        return str(symbols.quark_vector_weyl_plus)
    if chirality == -1:
        return str(symbols.quark_vector_weyl_minus)
    raise ValueError(f"Weyl vector tensor requires nonzero chirality, got {chirality}")


def _weyl_coupling_for_chirality(
    chirality: int,
    coupling: tuple[float, float],
) -> float:
    if chirality == 1:
        return coupling[1]
    if chirality == -1:
        return coupling[0]
    raise ValueError(f"Weyl coupling requires nonzero chirality, got {chirality}")


def _coupling_pair(coupling: Any) -> tuple[float, float]:
    if len(coupling) != 2:
        raise ValueError(f"expected a two-component coupling, got {coupling}")
    return float(coupling[0]), float(coupling[1])


def _source_current_count(graph: Any) -> int:
    return len(_source_currents(graph))


def _source_currents(graph: Any) -> tuple[Any, ...]:
    interactions_by_result: dict[tuple[int, tuple[int, ...], int], list[Any]] = {}
    for interaction in graph.interactions:
        interactions_by_result.setdefault(
            _current_key_tuple(interaction.result),
            [],
        ).append(interaction)
    source_keys: set[tuple[int, tuple[int, ...], int]] = set()
    seen: set[tuple[int, tuple[int, ...], int]] = set()

    def visit(current: Any) -> None:
        key = _current_key_tuple(current)
        if key in seen:
            return
        seen.add(key)
        interactions = interactions_by_result.get(key)
        if not interactions:
            source_keys.add(key)
            return
        for interaction in interactions:
            visit(interaction.left)
            visit(interaction.right)

    for left, right in graph.amplitudes:
        visit(left)
        visit(right)
    return tuple(
        current
        for current in graph.currents
        if _current_key_tuple(current) in source_keys
    )


def _register_parametric_source_currents(
    library: Any,
    graph: Any,
    builder: ParamBuilder | None = None,
) -> ParamBuilder:
    from symbolica.community.spenso import Representation

    if builder is None:
        builder = ParamBuilder()
    mink = Representation.mink(4)
    aux6 = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    for current in _source_currents(graph):
        pdg = int(current.pdg)
        if pdg == -21:
            representation = aux6
        elif pdg == 21 or pdg == 22 or pdg == 23:
            representation = mink
        elif 1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16:
            representation = weyl
        else:
            continue
        builder.register_rank1_tensor(
            library,
            tensor_name=_current_tensor_name(current),
            representation=representation,
            head=_current_parameter_head(current),
            length=_current_dimension(current),
            role="source_current",
        )
    return builder


def _register_parametric_current_momenta(
    model: AmplicolSMLeadingColorModel,
    library: Any,
    graph: Any,
    builder: ParamBuilder | None = None,
) -> ParamBuilder:
    from symbolica.community.spenso import LibraryTensor, Representation, TensorName

    if builder is None:
        builder = ParamBuilder()
    mink = Representation.mink(4)
    weyl = Representation("pyamplicol::weyl_spinor", 2)
    momentum_symbols_by_current: dict[tuple[int, tuple[int, ...], int], tuple[Any, ...]] = {}
    for current in _current_momentum_currents(graph):
        momentum_symbols = builder.register_rank1_tensor(
            library,
            tensor_name=_current_momentum_tensor_name(current),
            representation=mink,
            head=_current_momentum_parameter_head(current),
            length=4,
            role="current_momentum",
            real_valued=True,
        )
        momentum_symbols_by_current[_current_key_tuple(current)] = momentum_symbols

    for current in _propagating_currents(graph):
        momentum_symbols = momentum_symbols_by_current[_current_key_tuple(current)]
        pdg = int(current.pdg)
        if pdg == 21:
            library.register(
                LibraryTensor.dense(
                    TensorName(_propagator_tensor_name(current))(mink, mink),
                    model.gluon_propagator_tensor_data(momentum_symbols),
                )
            )
        elif _is_weyl_fermion_current(current):
            library.register(
                LibraryTensor.dense(
                    TensorName(_propagator_tensor_name(current))(weyl, weyl),
                    model.quark_weyl_propagator_tensor_data(
                        momentum_symbols,
                        chirality=int(current.chirality),
                    ),
                )
            )
    return builder


def _current_momentum_currents(graph: Any) -> tuple[Any, ...]:
    currents: dict[tuple[int, tuple[int, ...], int], Any] = {}
    for interaction in graph.interactions:
        if int(interaction.vertex_kind) == 0:
            currents.setdefault(_current_key_tuple(interaction.left), interaction.left)
            currents.setdefault(_current_key_tuple(interaction.right), interaction.right)
    for current in _propagating_currents(graph):
        currents.setdefault(_current_key_tuple(current), current)
    return tuple(currents.values())


def _propagating_currents(graph: Any) -> tuple[Any, ...]:
    result_keys = {
        _current_key_tuple(interaction.result) for interaction in graph.interactions
    }
    amplitude_keys = {
        _current_key_tuple(current)
        for amplitude in graph.amplitudes
        for current in amplitude
    }
    return tuple(
        current
        for current in graph.currents
        if _current_key_tuple(current) in result_keys
        and _current_key_tuple(current) not in amplitude_keys
        and _has_supported_propagator(current)
    )


def _current_needs_propagator(graph: Any, current: Any) -> bool:
    propagating_keys = {
        _current_key_tuple(propagating_current)
        for propagating_current in _propagating_currents(graph)
    }
    return _current_key_tuple(current) in propagating_keys


def _has_supported_propagator(current: Any) -> bool:
    pdg = int(current.pdg)
    return pdg == 21 or _is_weyl_fermion_current(current)


def _is_weyl_fermion_current(current: Any) -> bool:
    pdg = int(current.pdg)
    return (1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16) and int(current.chirality) != 0


def _current_dimension(current: Any) -> int:
    pdg = int(current.pdg)
    if pdg == -21:
        return 6
    if pdg == 21 or pdg == 22 or pdg == 23:
        return 4
    if 1 <= abs(pdg) <= 6 or 11 <= abs(pdg) <= 16:
        return 2
    return 0


def _current_parameter_head(current: Any) -> tuple[str, ...]:
    labels = tuple(str(label) for label in current.external_labels)
    return (
        "current",
        _signed_label(int(current.pdg)),
        "_".join(labels) if labels else "empty",
        _signed_label(int(current.chirality)),
    )


def _current_momentum_tensor_name(current: Any) -> str:
    labels = "_".join(str(label) for label in current.external_labels)
    return (
        "pyamplicol::"
        f"current_momentum_{_signed_label(int(current.pdg))}_{labels}_"
        f"{_signed_label(int(current.chirality))}"
    )


def _propagator_tensor_name(current: Any) -> str:
    labels = "_".join(str(label) for label in current.external_labels)
    if int(current.pdg) == 21:
        head = "gluon_propagator"
    elif _is_weyl_fermion_current(current):
        head = "quark_weyl_propagator"
    else:
        head = "propagator"
    return (
        "pyamplicol::"
        f"{head}_{_signed_label(int(current.pdg))}_{labels}_"
        f"{_signed_label(int(current.chirality))}"
    )


def _current_momentum_parameter_head(current: Any) -> tuple[str, ...]:
    labels = tuple(str(label) for label in current.external_labels)
    return (
        "current_momentum",
        _signed_label(int(current.pdg)),
        "_".join(labels) if labels else "empty",
        _signed_label(int(current.chirality)),
    )


def _propagator_lowering_ready(graph: Any) -> bool:
    result_currents = {
        _current_key_tuple(interaction.result): interaction.result
        for interaction in graph.interactions
    }
    amplitude_keys = {
        _current_key_tuple(current)
        for amplitude in graph.amplitudes
        for current in amplitude
    }
    for key, current in result_currents.items():
        if key in amplitude_keys or int(current.pdg) == -21:
            continue
        if not _has_supported_propagator(current):
            return False
    return bool(graph.interactions)


def _build_auxiliary_tensor_probe(
    model: AmplicolSMLeadingColorModel,
) -> TensorNetworkProbe:
    from symbolica.community.spenso import Representation, TensorName, TensorNetwork

    library = model.build_tensor_library()
    mink = Representation.mink(4)
    antisym = Representation("pyamplicol::antisymmetric_lorentz_pair", 6)
    two_gluon_to_tensor = TensorName(str(symbols.two_gluon_to_tensor))
    tensor_gluon_to_gluon = TensorName(str(symbols.tensor_gluon_to_gluon))

    expression = (
        two_gluon_to_tensor(
            mink("mu"),
            mink("nu"),
            antisym("A"),
        ).to_expression()
        * tensor_gluon_to_gluon(
            antisym("A"),
            mink("nu"),
            mink("rho"),
        ).to_expression()
    )
    network = TensorNetwork(expression, library)
    network.execute(library=library)
    result = network.result_tensor(library)
    output_size = len(result)
    structure = result.structure()
    entries = tuple(complex(result[i]) for i in range(output_size))
    nonzero = tuple(
        (index, value) for index, value in enumerate(entries) if abs(value) > 1.0e-15
    )
    weighted_checksum = sum((index + 1) * value for index, value in enumerate(entries))
    return TensorNetworkProbe(
        engine="spenso",
        tensor_names=(
            str(symbols.two_gluon_to_tensor),
            str(symbols.tensor_gluon_to_gluon),
        ),
        expression=_clean_symbolica_string(str(expression)),
        output_structure=_clean_symbolica_string(str(structure)),
        output_rank=2,
        output_size=output_size,
        nonzero_entries=len(nonzero),
        max_abs_entry=max((abs(value) for value in entries), default=0.0),
        weighted_checksum=(weighted_checksum.real, weighted_checksum.imag),
        first_nonzero_entries=tuple(
            (index, _format_complex(value)) for index, value in nonzero[:4]
        ),
    )


def _build_color_probe() -> ColorAlgebraProbe:
    from symbolica import S
    from symbolica.community.idenso import simplify_color
    from symbolica.community.spenso import Representation

    structure_constant = S("spenso::f")
    adjoint = Representation("coad", 8)

    def color_f(i: int, j: int, k: int) -> Any:
        return structure_constant(
            adjoint(i).to_expression(),
            adjoint(j).to_expression(),
            adjoint(k).to_expression(),
        )

    expression = color_f(1, 2, 3) * color_f(3, 2, 1)
    simplified = simplify_color(expression)
    return ColorAlgebraProbe(
        engine="idenso",
        input_expression=_clean_symbolica_string(str(expression)),
        simplified_expression=_clean_symbolica_string(str(simplified)),
    )


def _format_complex(value: complex) -> str:
    real = 0.0 if abs(value.real) < 1.0e-15 else value.real
    imag = 0.0 if abs(value.imag) < 1.0e-15 else value.imag
    return f"{real:.16g}{imag:+.16g}j"


def _clean_symbolica_string(value: str) -> str:
    return re.sub(r"\x1b\[[0-9;]*m", "", value)


def _preview_expression(value: str) -> str:
    if len(value) <= _MAX_RECURSION_EXPRESSION_PREVIEW:
        return value
    suffix = "...<truncated>"
    preview_length = _MAX_RECURSION_EXPRESSION_PREVIEW - len(suffix)
    return f"{value[:preview_length]}{suffix}"


__all__ = [
    "ColorAlgebraProbe",
    "RecursionLoweringPlan",
    "SymbolicLoweringReport",
    "TensorNetworkBlueprint",
    "TensorNetworkProbe",
    "VertexLoweringReport",
    "VertexLoweringStep",
    "build_tensor_network_scalar_bundle",
    "build_symbolic_lowering_report",
]
