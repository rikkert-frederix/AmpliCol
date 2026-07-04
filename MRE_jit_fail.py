#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from symbolica import E, Expression, S


# Bug reproduced by this MRE:
# On AArch64, Symbolica/SymJIT successfully builds a JIT evaluator for the
# expressions in MRE_jit_fail_input.txt, but the first evaluate_complex() call
# aborts in symjit-2.18.6/rust/arm/vector.rs:795 with
# `assertion failed: x.abs() < 1048576`. The input file contains one dumped
# pyamplicol shared-current block for d d~ -> Z + 7g, alpha-renamed to p0..pN
# only; it is intentionally not simplified or minimized.
INPUT_FILE = Path(__file__).with_name("MRE_jit_fail_input.txt")


def main() -> int:
    payload = json.loads(INPUT_FILE.read_text(encoding="utf-8"))
    parameter_count = int(payload["parameter_count"])
    output_texts = list(payload["outputs"])
    real_params = [int(index) for index in payload.get("real_params", [])]
    evaluator_kwargs = dict(payload["evaluator_kwargs"])

    params = _symbols(parameter_count)
    outputs = tuple(E(text) for text in output_texts)

    print(
        "Building Symbolica evaluator for "
        f"{len(outputs)} outputs and {len(params)} complex parameters...",
        flush=True,
    )
    evaluator = Expression.evaluator_multiple(outputs, params, **evaluator_kwargs)
    if real_params:
        evaluator.set_real_params(real_params, verbose=False)

    print(
        "Evaluator built. Calling evaluate_complex; on affected AArch64 SymJIT "
        "builds this aborts in symjit/rust/arm/vector.rs:795.",
        flush=True,
    )
    parameter_rows = np.zeros((1, parameter_count), dtype=np.complex128)
    result = evaluator.evaluate_complex(parameter_rows)
    print(f"Unexpectedly succeeded, result shape: {result.shape}", flush=True)
    return 0


def _symbols(count: int) -> list[object]:
    symbols: list[object] = []
    for start in range(0, count, 1000):
        names = [f"p{index}" for index in range(start, min(start + 1000, count))]
        created = S(*names)
        if len(names) == 1:
            symbols.append(created)
        else:
            symbols.extend(created)
    return symbols


if __name__ == "__main__":
    raise SystemExit(main())
