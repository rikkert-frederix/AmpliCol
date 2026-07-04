from __future__ import annotations

import base64
import math
import time
from dataclasses import dataclass
from functools import cached_property
from typing import Any, Sequence

from .model import AmplicolSMLeadingColorModel
from .native import (
    ExternalMomentum,
    MatrixElementEvaluation,
    NativeEvaluationError,
    _ext_antiquark,
    _ext_massive_vector,
    _ext_quark,
    _external_current_input,
    _initial_state_average_factor,
)
from .symbols import symbols

ComplexExpression = tuple[Any, Any]


@dataclass(frozen=True)
class SymbolicEvaluatorMetadata:
    engine: str
    kernel: str
    parameter_names: tuple[str, ...]
    expression: str
    expression_length: int
    instruction_count: int
    temporary_count: int
    constant_count: int
    evaluator_state_b64: str
    build_time_s: float

    def to_json_dict(self) -> dict[str, object]:
        return {
            "engine": self.engine,
            "kernel": self.kernel,
            "parameter_names": list(self.parameter_names),
            "parameter_count": len(self.parameter_names),
            "expression": self.expression,
            "expression_length": self.expression_length,
            "instruction_count": self.instruction_count,
            "temporary_count": self.temporary_count,
            "constant_count": self.constant_count,
            "evaluator_state_b64": self.evaluator_state_b64,
            "build_time_s": self.build_time_s,
        }


class ZeroGluonSymbolicEvaluator:
    """Symbolica scalar evaluator for q q~ -> Z raw helicity sums."""

    def __init__(
        self,
        *,
        model: AmplicolSMLeadingColorModel | None = None,
        evaluator_state: bytes | None = None,
        parameter_names: Sequence[str] | None = None,
    ) -> None:
        self.model = model or AmplicolSMLeadingColorModel()
        if evaluator_state is None:
            spec = self._spec
            self._evaluator = spec["evaluator"]
            self._parameter_names = spec["parameter_names"]
        else:
            from symbolica import Evaluator

            self._evaluator = Evaluator.load(evaluator_state)
            self._parameter_names = (
                tuple(parameter_names)
                if parameter_names is not None
                else tuple(_zero_gluon_parameter_names())
            )

    def evaluate(
        self,
        process: str,
        particles: Sequence[ExternalMomentum],
    ) -> MatrixElementEvaluation:
        point = _validate_zero_gluon_point(particles)
        parameter_values = self._parameter_values(process, point)
        raw_sum = float(self._evaluator.evaluate([parameter_values])[0][0])
        pdgs = tuple(particle.pdg for particle in point)
        color_factor = self.model.leading_color_factor(pdgs)
        average_factor = _initial_state_average_factor(pdgs[:2])
        coupling_factor = 2.0 * 4.0 * math.pi * self.model.alpha_ew
        matrix_element = raw_sum * color_factor * coupling_factor / average_factor
        return MatrixElementEvaluation(
            process=process,
            particles=point,
            matrix_element=matrix_element,
            raw_helicity_sum=raw_sum,
            color_factor=color_factor,
            average_factor=average_factor,
            coupling_factor=coupling_factor,
            helicity_contributions=(),
        )

    @cached_property
    def metadata(self) -> SymbolicEvaluatorMetadata:
        spec = self._spec
        evaluator = self._evaluator
        instructions, temporary_count, constants = evaluator.get_instructions()
        return SymbolicEvaluatorMetadata(
            engine="symbolica",
            kernel="symbolica-zero-gluon",
            parameter_names=spec["parameter_names"],
            expression=spec["expression_preview"],
            expression_length=spec["expression_length"],
            instruction_count=len(instructions),
            temporary_count=temporary_count,
            constant_count=len(constants),
            evaluator_state_b64=base64.b64encode(evaluator.save()).decode("ascii"),
            build_time_s=spec["build_time_s"],
        )

    @classmethod
    def from_artifact_payload(
        cls,
        payload: dict[str, object],
        *,
        model: AmplicolSMLeadingColorModel | None = None,
    ) -> ZeroGluonSymbolicEvaluator:
        encoded = payload.get("evaluator_state_b64")
        if not isinstance(encoded, str):
            raise ValueError("symbolic evaluator artifact is missing evaluator_state_b64")
        parameter_names = payload.get("parameter_names")
        if not isinstance(parameter_names, list) or not all(
            isinstance(name, str) for name in parameter_names
        ):
            raise ValueError("symbolic evaluator artifact is missing parameter_names")
        return cls(
            model=model,
            evaluator_state=base64.b64decode(encoded),
            parameter_names=parameter_names,
        )

    @cached_property
    def _spec(self) -> dict[str, Any]:
        return _build_zero_gluon_symbolic_spec()

    def _parameter_values(
        self,
        process: str,
        particles: tuple[ExternalMomentum, ExternalMomentum, ExternalMomentum],
    ) -> list[float]:
        current_inputs = [
            _external_current_input(self.model, particle, is_initial=i < 2)
            for i, particle in enumerate(particles)
        ]
        quark = next(item for item in current_inputs[:2] if item[0] > 0)
        antiquark = next(item for item in current_inputs[:2] if item[0] < 0)
        z_momentum = current_inputs[2][1]
        coupling = self.model.z_fermion_coupling(abs(particles[0].pdg))

        values: dict[str, float] = {
            "z_left": coupling[0],
            "z_right": coupling[1],
        }
        for helicity in (-1, 1):
            q_wf = _ext_quark(quark[1], helicity, self.model.mass(quark[0]))
            aq_wf = _ext_antiquark(
                antiquark[1],
                helicity,
                self.model.mass(antiquark[0]),
            )
            _add_complex_wavefunction_values(values, "q", helicity, q_wf)
            _add_complex_wavefunction_values(values, "aq", helicity, aq_wf)
        for helicity in (-1, 0, 1):
            z_wf = _ext_massive_vector(z_momentum, helicity, self.model.mass(23))
            _add_complex_wavefunction_values(values, "z", helicity, z_wf)

        return [values[name] for name in self._parameter_names]


def build_zero_gluon_symbolic_evaluator_payload() -> dict[str, object]:
    return ZeroGluonSymbolicEvaluator().metadata.to_json_dict()


def _build_zero_gluon_symbolic_spec() -> dict[str, Any]:
    from symbolica import Expression

    start = time.perf_counter()
    parameters: dict[str, Any] = {}
    parameter_names = tuple(_zero_gluon_parameter_names())
    for name in parameter_names:
        parameters[name] = symbols.parameter(f"zero_gluon::{name}")

    total = Expression.num(0)
    left = parameters["z_left"]
    right = parameters["z_right"]
    inv_sqrt2 = Expression.num(1.0 / math.sqrt(2.0))
    for h_quark in (-1, 1):
        for h_antiquark in (-1, 1):
            q = tuple(_complex_parameter(parameters, "q", h_quark, i) for i in range(4))
            aq = tuple(
                _complex_parameter(parameters, "aq", h_antiquark, i) for i in range(4)
            )
            current = _symbolic_fermion_antifermion_to_vector(
                q,
                aq,
                left,
                right,
                inv_sqrt2,
            )
            for h_z in (-1, 0, 1):
                z = tuple(_complex_parameter(parameters, "z", h_z, i) for i in range(4))
                amplitude = _csum(_cmul(current[i], z[i]) for i in range(4))
                total = total + _cabs2(amplitude)

    expression = str(total)
    evaluator = total.evaluator([parameters[name] for name in parameter_names])
    return {
        "evaluator": evaluator,
        "parameter_names": parameter_names,
        "expression_preview": _bounded_expression_preview(expression),
        "expression_length": len(expression),
        "build_time_s": time.perf_counter() - start,
    }


def _bounded_expression_preview(expression: str, limit: int = 4096) -> str:
    suffix = "...<truncated>"
    if len(expression) <= limit:
        return expression
    return expression[: limit - len(suffix)] + suffix


def _symbolic_fermion_antifermion_to_vector(
    quark: tuple[ComplexExpression, ...],
    antiquark: tuple[ComplexExpression, ...],
    left: Any,
    right: Any,
    inv_sqrt2: Any,
) -> tuple[ComplexExpression, ComplexExpression, ComplexExpression, ComplexExpression]:
    l1, l2, l3, l4 = quark
    a1, a2, a3, a4 = antiquark
    return (
        _imul(
            _cscale(
                _cadd(
                    _cscale(_cadd(_cmul(l3, a1), _cmul(l4, a2)), left),
                    _cscale(_cadd(_cmul(l1, a3), _cmul(l2, a4)), right),
                ),
                inv_sqrt2,
            )
        ),
        _imul(
            _cscale(
                _cadd(
                    _cscale(_cadd(_cneg(_cmul(l4, a1)), _cneg(_cmul(l3, a2))), left),
                    _cscale(_cadd(_cmul(l1, a4), _cmul(l2, a3)), right),
                ),
                inv_sqrt2,
            )
        ),
        _cscale(
            _cadd(
                _cscale(_cadd(_cneg(_cmul(l4, a1)), _cmul(l3, a2)), left),
                _cscale(_cadd(_cneg(_cmul(l1, a4)), _cmul(l2, a3)), right),
            ),
            -inv_sqrt2,
        ),
        _imul(
            _cscale(
                _cadd(
                    _cscale(_cadd(_cneg(_cmul(l3, a1)), _cmul(l4, a2)), left),
                    _cscale(_cadd(_cmul(l1, a3), _cneg(_cmul(l2, a4))), right),
                ),
                inv_sqrt2,
            )
        ),
    )


def _zero_gluon_parameter_names() -> list[str]:
    names: list[str] = []
    for prefix, helicities in (("q", (-1, 1)), ("aq", (-1, 1)), ("z", (-1, 0, 1))):
        for helicity in helicities:
            for component in range(4):
                names.extend(
                    [
                        _component_name(prefix, helicity, component, "re"),
                        _component_name(prefix, helicity, component, "im"),
                    ]
                )
    names.extend(["z_left", "z_right"])
    return names


def _component_name(prefix: str, helicity: int, component: int, part: str) -> str:
    sign = "m" if helicity < 0 else "p"
    return f"{prefix}_{sign}{abs(helicity)}_{component}_{part}"


def _complex_parameter(
    parameters: dict[str, Any],
    prefix: str,
    helicity: int,
    component: int,
) -> ComplexExpression:
    return (
        parameters[_component_name(prefix, helicity, component, "re")],
        parameters[_component_name(prefix, helicity, component, "im")],
    )


def _add_complex_wavefunction_values(
    values: dict[str, float],
    prefix: str,
    helicity: int,
    wavefunction: Sequence[complex],
) -> None:
    for index, value in enumerate(wavefunction):
        values[_component_name(prefix, helicity, index, "re")] = value.real
        values[_component_name(prefix, helicity, index, "im")] = value.imag


def _validate_zero_gluon_point(
    particles: Sequence[ExternalMomentum],
) -> tuple[ExternalMomentum, ExternalMomentum, ExternalMomentum]:
    if len(particles) != 3:
        raise NativeEvaluationError("q q~ -> Z requires exactly three external momenta")
    point = tuple(particles)
    if point[2].pdg != 23:
        raise NativeEvaluationError("third external particle must be a Z boson")
    if point[0].pdg + point[1].pdg != 0:
        raise NativeEvaluationError("incoming particles must be a quark/antiquark pair")
    if not 1 <= abs(point[0].pdg) <= 6:
        raise NativeEvaluationError("incoming pair must be quarks")
    return point[0], point[1], point[2]


def _cadd(a: ComplexExpression, b: ComplexExpression) -> ComplexExpression:
    return a[0] + b[0], a[1] + b[1]


def _cneg(a: ComplexExpression) -> ComplexExpression:
    return -a[0], -a[1]


def _cmul(a: ComplexExpression, b: ComplexExpression) -> ComplexExpression:
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def _cscale(a: ComplexExpression, factor: Any) -> ComplexExpression:
    return factor * a[0], factor * a[1]


def _imul(a: ComplexExpression) -> ComplexExpression:
    return -a[1], a[0]


def _cabs2(a: ComplexExpression) -> Any:
    return a[0] * a[0] + a[1] * a[1]


def _csum(values: Any) -> ComplexExpression:
    iterator = iter(values)
    total = next(iterator)
    for value in iterator:
        total = _cadd(total, value)
    return total


__all__ = [
    "SymbolicEvaluatorMetadata",
    "ZeroGluonSymbolicEvaluator",
    "build_zero_gluon_symbolic_evaluator_payload",
]
