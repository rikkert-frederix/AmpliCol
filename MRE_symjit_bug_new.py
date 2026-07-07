#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from symbolica import Expression, S


# Bug reproduced by this MRE:
# This is the first output expression passed to Symbolica for JIT compilation
# by the pyAmpliCol generic-DAG staged evaluator for:
#
#   d d~ > u u~ s s~ c c~
#
# using leading-colour sector 0 and stage `generic_stage_1_subset_2`.
# With the historical pyAmpliCol managed Symbolica/SymJIT v220 build, Symbolica
# successfully returns the JIT evaluator object, but the first
# evaluate_complex() call segfaults on AArch64:
#
#   symbolica 2.1.0 from the pyAmpliCol managed venv
#   symjit source ref v220, rev 773294f2a9e1b3359cba381e8d5106a854dd3219
#
# The same expression evaluates normally when `jit_compile` is set to False,
# and also succeeds with SymJIT 2.19.2 rev
# 54d6c6171f05b39505d18d1932f0972bfac9e4da. The input file keeps the
# original 428-parameter layout and real-parameter metadata from the
# pyAmpliCol stage compiler.
#
# Reproduction command:
#
#   PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python \
#     MRE_symjit_bug_new.py
#
# Observed behavior on macOS/AArch64:
#
#   Loaded pyAmpliCol generic-DAG SymJIT MRE
#     expression: 2.35702260395516e-1*(c2*c4+c3*c5)/(c32^2-c33^2-c34^2-c35^2)
#     parameters: 428
#     real params: 92
#   Building Symbolica evaluator with jit_compile=True...
#   Evaluator object built. Calling evaluate_complex()...
#   zsh: segmentation fault  ... MRE_symjit_bug_new.py
#
# lldb backtrace captured with:
#
#   lldb --batch \
#     -o 'settings set target.process.stop-on-exec false' \
#     -o 'run' \
#     -k 'thread backtrace all' \
#     -k 'register read x10 x21 x22' \
#     -- pyAmpliCol/dependencies/.venv/bin/python MRE_symjit_bug_new.py
#
# Backtrace:
#
#   Process stopped with EXC_BAD_ACCESS (code=1, address=0x0)
#   * thread #2, queue = 'com.apple.main-thread'
#     * frame #0: 0x0000039a42c88180
#       -> 0x39a42c88180: str d0, [x10, x22, lsl #3]
#          0x39a42c88184: ldr x10, [x21, #0x10]
#          0x39a42c88188: ldr d0, [x19, #0x8]
#          0x39a42c8818c: str d0, [x10, x22, lsl #3]
#       frame #1: core.abi3.so`pyo3::impl_::trampoline::trampoline + 84
#       frame #2: core.abi3.so`pyo3::impl_::trampoline::cfunction_with_keywords + 64
#       frame #3: Python`method_vectorcall_VARARGS_KEYWORDS + 148
#       frame #4: Python`_PyEval_EvalFrameDefault + 43820
#       frame #5: Python`PyEval_EvalCode + 184
#       frame #6: Python`run_eval_code_obj + 88
#       frame #7: Python`run_mod + 132
#       frame #8: Python`pyrun_file + 156
#       frame #9: Python`_PyRun_SimpleFileObject + 288
#       frame #10: Python`_PyRun_AnyFileObject + 80
#       frame #11: Python`pymain_run_file_obj + 164
#       frame #12: Python`pymain_run_file + 72
#       frame #13: Python`Py_RunMain + 760
#       frame #14: Python`pymain_main + 304
#       frame #15: Python`Py_BytesMain + 40
#       frame #16: dyld`start + 2840
#
#   register read x10 x21 x22:
#     x10 = 0x0000000000000000
#     x21 = 0x0000039a42c72800
#     x22 = 0x0000000000000000
#
# No Python exception is raised because the process receives SIGSEGV in the
# native JIT-generated code reached from Evaluator.evaluate_complex().
INPUT_FILE = Path(__file__).with_name("MRE_symjit_bug_new_input.txt")


def main() -> int:
    payload = json.loads(INPUT_FILE.read_text(encoding="utf-8"))
    parameter_count = int(payload["parameter_count"])
    expression_text = str(payload["expression"])
    evaluator_kwargs = dict(payload["evaluator_kwargs"])
    real_params = [int(index) for index in payload.get("real_params", [])]

    params = _symbols(parameter_count)
    expression = Expression.parse(expression_text)

    print("Loaded pyAmpliCol generic-DAG SymJIT MRE", flush=True)
    print(f"  expression: {expression_text}", flush=True)
    print(f"  parameters: {len(params)}", flush=True)
    print(f"  real params: {len(real_params)}", flush=True)
    print("Building Symbolica evaluator with jit_compile=True...", flush=True)

    evaluator = Expression.evaluator_multiple(
        (expression,),
        params,
        **evaluator_kwargs,
    )
    if real_params:
        evaluator.set_real_params(real_params, verbose=False)

    print("Evaluator object built. Calling evaluate_complex()...", flush=True)
    parameter_rows = np.ones((1, parameter_count), dtype=np.complex128)
    result = evaluator.evaluate_complex(parameter_rows)
    print(f"Unexpectedly succeeded: {result}", flush=True)
    return 0


def _symbols(count: int) -> list[object]:
    result: list[object] = []
    for start in range(0, count, 1000):
        names = [f"c{index}" for index in range(start, min(start + 1000, count))]
        created = S(*names)
        if len(names) == 1:
            result.append(created)
        else:
            result.extend(created)
    return result


if __name__ == "__main__":
    raise SystemExit(main())
