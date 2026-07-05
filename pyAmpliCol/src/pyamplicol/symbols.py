from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Any


@dataclass(frozen=True)
class SymbolRegistry:
    """Central registry for reusable Symbolica heads used by pyamplicol."""

    namespace: str = "pyamplicol"

    def symbol(self, name: str) -> Any:
        from symbolica import S

        return S(f"{self.namespace}::{name}")

    def real_symbol(self, name: str) -> Any:
        from symbolica import S

        return S(f"{self.namespace}::{name}", is_real=True)

    def parameter(self, name: str) -> Any:
        return self.symbol(f"param::{name}")

    @cached_property
    def two_gluon_to_tensor(self) -> Any:
        return self.symbol("two_gluon_to_tensor")

    @cached_property
    def tensor_gluon_to_gluon(self) -> Any:
        return self.symbol("tensor_gluon_to_gluon")

    @cached_property
    def gluon_tensor_to_gluon(self) -> Any:
        return self.symbol("gluon_tensor_to_gluon")

    @cached_property
    def quark_vector_weyl_plus(self) -> Any:
        return self.symbol("quark_vector_weyl_plus")

    @cached_property
    def quark_vector_weyl_minus(self) -> Any:
        return self.symbol("quark_vector_weyl_minus")

    @cached_property
    def current(self) -> Any:
        return self.symbol("current")

    @cached_property
    def vertex(self) -> Any:
        return self.symbol("vertex")

    @cached_property
    def assignment(self) -> Any:
        return self.symbol("assign")

    @cached_property
    def amplitude(self) -> Any:
        return self.symbol("amplitude")

    @cached_property
    def matrix_element_plan(self) -> Any:
        return self.symbol("matrix_element_plan")

    @cached_property
    def momentum(self) -> Any:
        return self.symbol("momentum")

    @cached_property
    def current_momentum(self) -> Any:
        return self.symbol("current_momentum")

    @cached_property
    def polarization(self) -> Any:
        return self.symbol("polarization")


symbols = SymbolRegistry()


__all__ = ["SymbolRegistry", "symbols"]
