# pyAmplicCol Historical Performance Summary

> **Historical note:** this file records the earlier `d d~ -> Z + n g`
> benchmark campaign.  The live refreshed results now live in
> `pyAmpliCol.pdf`, generated from `z_performance_data.json`,
> `z_performance_table.tex`, and the LC/NLC/full-colour result-matrix caches.
> The current production defaults are Rusticol schema-v2 process artifacts,
> SymJIT `opt_level=3`, stage-local evaluator inputs, ten Horner iterations,
> Symbolica's backend-default CPE choice, and enlarged common-pair/Horner
> limits. The live result matrix deliberately uses O1 throughout; repeated
> ten-second checks have not found a consistent, material O3 gain for a current
> matrix cell. Selected-flow
> JIT artifacts use runtime batch 128 and measure stage chunks around base 128;
> documented all-flow artifacts use uniform chunk 8192 and batch 64. Timings
> use `time-process --target-runtime 10`. The rows below are kept for provenance
> and should not be used as the current regenerated table.

Process family: `d d~ -> Z + n g`, for `n = 1,...,9`.

Successful pyAmplicCol rows validated against Fortran AmpliCol with max relative differences typically of order `1e-13`. The rusticol rows shown for `n = 1,...,9` were validated against Fortran AmpliCol probes; the refreshed low-multiplicity Fortran timing references use a direct generated amplitude-library benchmark, with max relative differences no larger than `2.1e-12` in the validation probes. The table below omits the validation column and focuses on generation time, full pyAmplicCol wall runtime, and evaluator-only runtime.

For pyAmplicCol rows, the multiplier in the `Gen` column is relative to AmpliCol generation time for the same `n`. Multipliers in `Wall` and `Eval` are relative to AmpliCol's per-point evaluator timing for the same `n`. Color convention: <span style="color:#1a7f37">green</span> is faster than AmpliCol, <span style="color:#bf8700">orange</span> is slower but below `x2.0`, and <span style="color:#cf222e">red</span> is `x2.0` or worse.

For C++ and ASM rows with multiple measured chunking options, this table reports the existing result with the best evaluator-only runtime, except where an additional O3 C++ row was explicitly added while retaining the previous lower-optimization row. The refreshed C++ O3 save/load rows use batch size 128 after a loaded-artifact batch-size scan; retained older rows use the batch size shown in their notes. The rusticol rows use C++ O3 compiled-complex artifacts; rusticol can also load saved JIT artifacts, but direct rusticol probes at low multiplicity confirmed C++ O3 remains faster. All benchmark measurements summarized here were obtained under a watchdog limiting process-tree RAM usage to less than 30 GB.

The compiled pyAmplicCol DAG marks all known-real momentum-sum parameters as real with Symbolica `set_real_params(...)` before evaluator generation/compilation. In this hot path the couplings are embedded constants, so the relevant runtime input realness is the momentum sector. The refreshed low-multiplicity AmpliCol reference rows use the direct generated-library benchmark: `./amplicol_generate --library=create --amplicol_momenta_probe=10`, `make amplicol_generate_library`, `make amplicol_library_benchmark`, then `./amplicol_library_benchmark N 1 1`, which calls `amp_lib:evaluate_amp` directly and bypasses integration/probe bookkeeping.

The dependency installer now uses upstream Symbolica `dev` plus the local pyAmplicCol patches and pins SymJIT to `2.19.3` at commit `7fb09d1cb2a943c25a6fd71a208af44fcc6d813d`. Managed patches cover the historical AArch64 complex-JIT register-allocation failure, external-call spill offsets at the 4096-byte boundary, stack adjustments above the shifted 12-bit immediate range, and conditional branches beyond the signed 19-bit displacement. The historical rows below use the optimization levels recorded in their notes; the current matrix policy is O1 throughout.

| n | Setup | Gen [s] | Wall [us/pt] | Eval [us/pt] | Notes |
|---:|---|---:|---:|---:|---|
| **1** | **AmpliCol** | **2.33** | **N/A** | **0.0225** | **Fortran AmpliCol; direct generated-library benchmark with momenta warmup** |
| 1 | pyAmplicCol - JIT | 0.016 <span style="color:#1a7f37">(x0.01)</span> | 13.34(29) <span style="color:#cf222e">(x592.52)</span> | 0.607(17) <span style="color:#cf222e">(x26.96)</span> | JIT; batch=16 |
| 1 | pyAmplicCol - C++ | 0.774 <span style="color:#1a7f37">(x0.33)</span> | 13.01(12) <span style="color:#cf222e">(x577.86)</span> | 0.521(10) <span style="color:#cf222e">(x23.14)</span> | C++; O3; chunk=None; batch=16 |
| 1 | pyAmplicCol - ASM | 0.560 <span style="color:#1a7f37">(x0.24)</span> | 14.03(30) <span style="color:#cf222e">(x623.17)</span> | 0.912(18) <span style="color:#cf222e">(x40.51)</span> | ASM; O3; chunk=None; batch=16 |
| 1 | pyAmplicCol - Rusticol | 0.984 <span style="color:#1a7f37">(x0.42)</span> | 0.412(1) <span style="color:#cf222e">(x18.30)</span> | 0.1938(2) <span style="color:#cf222e">(x8.61)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; single-chunk direct output; batch=1024 |
| **2** | **AmpliCol** | **2.35** | **N/A** | **0.101** | **Fortran AmpliCol; direct generated-library benchmark with momenta warmup** |
| 2 | pyAmplicCol - JIT | 0.066 <span style="color:#1a7f37">(x0.03)</span> | 16.15(15) <span style="color:#cf222e">(x160.20)</span> | 1.61(2) <span style="color:#cf222e">(x15.97)</span> | JIT; batch=16 |
| 2 | pyAmplicCol - C++ | 1.86 <span style="color:#1a7f37">(x0.79)</span> | 16.05(14) <span style="color:#cf222e">(x159.21)</span> | 1.16(2) <span style="color:#cf222e">(x11.51)</span> | C++; O3; chunk=None; batch=16 |
| 2 | pyAmplicCol - ASM | 0.940 <span style="color:#1a7f37">(x0.40)</span> | 16.48(20) <span style="color:#cf222e">(x163.47)</span> | 2.16(2) <span style="color:#cf222e">(x21.43)</span> | ASM; O3; chunk=None; batch=16 |
| 2 | pyAmplicCol - Rusticol | 2.00 <span style="color:#1a7f37">(x0.85)</span> | 0.984(1) <span style="color:#cf222e">(x9.76)</span> | 0.653(1) <span style="color:#cf222e">(x6.48)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; single-chunk direct output; batch=1024 |
| **3** | **AmpliCol** | **2.30** | **N/A** | **1.52** | **Fortran AmpliCol; direct generated-library benchmark** |
| 3 | pyAmplicCol - JIT | 0.281 <span style="color:#1a7f37">(x0.12)</span> | 22.67(42) <span style="color:#cf222e">(x14.92)</span> | 4.72(5) <span style="color:#cf222e">(x3.11)</span> | JIT; batch=16 |
| 3 | pyAmplicCol - C++ | 5.69 <span style="color:#cf222e">(x2.47)</span> | 20.90(40) <span style="color:#cf222e">(x13.76)</span> | 3.06(4) <span style="color:#cf222e">(x2.01)</span> | C++; O3; chunk=None; batch=16 |
| 3 | pyAmplicCol - ASM | 1.52 <span style="color:#1a7f37">(x0.66)</span> | 23.87(45) <span style="color:#cf222e">(x15.71)</span> | 6.17(8) <span style="color:#cf222e">(x4.06)</span> | ASM; O3; chunk=None; batch=16 |
| 3 | pyAmplicCol - Rusticol | 5.48 <span style="color:#cf222e">(x2.38)</span> | 2.783(13) <span style="color:#bf8700">(x1.83)</span> | 2.177(10) <span style="color:#bf8700">(x1.43)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; batch=128 |
| **4** | **AmpliCol** | **2.80** | **N/A** | **4.70** | **Fortran AmpliCol; direct generated-library benchmark** |
| 4 | pyAmplicCol - JIT | 0.431 <span style="color:#1a7f37">(x0.15)</span> | 14.23(1) <span style="color:#cf222e">(x3.03)</span> | 11.20(1) <span style="color:#cf222e">(x2.38)</span> | JIT; batch=16; last successful JIT refresh |
| 4 | pyAmplicCol - C++ | 20.7 <span style="color:#cf222e">(x7.40)</span> | 27.78(13) <span style="color:#cf222e">(x5.91)</span> | 7.75(6) <span style="color:#bf8700">(x1.65)</span> | C++; O3; chunk=None; batch=16 |
| 4 | pyAmplicCol - ASM | 2.04 <span style="color:#1a7f37">(x0.73)</span> | 36.30(54) <span style="color:#cf222e">(x7.72)</span> | 15.80(15) <span style="color:#cf222e">(x3.36)</span> | ASM; O3; chunk=None; batch=16 |
| 4 | pyAmplicCol - Rusticol | 19.9 <span style="color:#cf222e">(x7.11)</span> | 7.23(7) <span style="color:#bf8700">(x1.54)</span> | 6.10(6) <span style="color:#bf8700">(x1.30)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; batch=128 |
| **5** | **AmpliCol** | **3.52** | **N/A** | **13.8** | **Fortran AmpliCol; direct generated-library benchmark** |
| 5 | pyAmplicCol - JIT | 0.991 <span style="color:#1a7f37">(x0.28)</span> | 44.33(5) <span style="color:#cf222e">(x3.21)</span> | 40.56(4) <span style="color:#cf222e">(x2.94)</span> | JIT; batch=16; last successful JIT refresh |
| 5 | pyAmplicCol - C++ | 54.3 <span style="color:#cf222e">(x15.41)</span> | 45.06(48) <span style="color:#cf222e">(x3.27)</span> | 20.83(26) <span style="color:#bf8700">(x1.51)</span> | C++; O3; chunk=64; batch=128; chunk compile workers=10; two-stage save/load |
| 5 | pyAmplicCol - C++ (O2) | 50.7 <span style="color:#cf222e">(x14.39)</span> | 51.04(87) <span style="color:#cf222e">(x3.70)</span> | 23.82(42) <span style="color:#bf8700">(x1.73)</span> | C++; O2; chunk=64; batch=16 |
| 5 | pyAmplicCol - ASM | 10.2 <span style="color:#cf222e">(x2.90)</span> | 75.69(88) <span style="color:#cf222e">(x5.49)</span> | 48.74(72) <span style="color:#cf222e">(x3.53)</span> | ASM; O3; chunk=64; batch=16 |
| 5 | pyAmplicCol - Rusticol | 50.2 <span style="color:#cf222e">(x14.25)</span> | 19.63(28) <span style="color:#bf8700">(x1.42)</span> | 17.50(25) <span style="color:#bf8700">(x1.27)</span> | Rusticol PyO3; C++ O3; chunk=64; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **6** | **AmpliCol** | **4.66** | **N/A** | **37.3** | **Fortran AmpliCol; with --library=use** |
| 6 | pyAmplicCol - JIT | 2.96 <span style="color:#1a7f37">(x0.64)</span> | 163.0(5) <span style="color:#cf222e">(x4.37)</span> | 157.0(5) <span style="color:#cf222e">(x4.21)</span> | JIT; batch=16; last successful JIT refresh |
| 6 | pyAmplicCol - C++ | 151 <span style="color:#cf222e">(x32.44)</span> | 94.79(73) <span style="color:#cf222e">(x2.54)</span> | 56.36(49) <span style="color:#bf8700">(x1.51)</span> | C++; O3; chunk=64; batch=128; chunk compile workers=10; two-stage save/load |
| 6 | pyAmplicCol - C++ (O2) | 140 <span style="color:#cf222e">(x30.07)</span> | 108.4(1.8) <span style="color:#cf222e">(x2.91)</span> | 68.5(1.6) <span style="color:#bf8700">(x1.84)</span> | C++; O2; chunk=64; batch=16 |
| 6 | pyAmplicCol - ASM | 21.4 <span style="color:#cf222e">(x4.60)</span> | 183.7(3.3) <span style="color:#cf222e">(x4.93)</span> | 141.2(2.7) <span style="color:#cf222e">(x3.79)</span> | ASM; O3; chunk=64; batch=16 |
| 6 | pyAmplicCol - Rusticol | 141 <span style="color:#cf222e">(x30.29)</span> | 55.63(26) <span style="color:#bf8700">(x1.49)</span> | 50.12(23) <span style="color:#bf8700">(x1.34)</span> | Rusticol PyO3; C++ O3; chunk=64; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **7** | **AmpliCol** | **9.28** | **N/A** | **95.2** | **Fortran AmpliCol; with --library=use** |
| 7 | pyAmplicCol - JIT | 8.85 <span style="color:#1a7f37">(x0.95)</span> | 536(2) <span style="color:#cf222e">(x5.63)</span> | 525(2) <span style="color:#cf222e">(x5.52)</span> | JIT; batch=16; last successful JIT refresh |
| 7 | pyAmplicCol - C++ | 417 <span style="color:#cf222e">(x44.92)</span> | 259.7(1.4) <span style="color:#cf222e">(x2.73)</span> | 194.7(1.2) <span style="color:#cf222e">(x2.05)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 7 | pyAmplicCol - C++ (O1) | 219 <span style="color:#cf222e">(x23.59)</span> | 421.6(2.4) <span style="color:#cf222e">(x4.43)</span> | 348.9(1.9) <span style="color:#cf222e">(x3.67)</span> | C++; O1; chunk=96; batch=16 |
| 7 | pyAmplicCol - ASM | 36.2 <span style="color:#cf222e">(x3.90)</span> | 384.4(2.8) <span style="color:#cf222e">(x4.04)</span> | 320.2(2.4) <span style="color:#cf222e">(x3.36)</span> | ASM; O3; chunk=96; batch=16 |
| 7 | pyAmplicCol - Rusticol | 391 <span style="color:#cf222e">(x42.12)</span> | 197.8(4) <span style="color:#cf222e">(x2.08)</span> | 186.1(4) <span style="color:#bf8700">(x1.95)</span> | Rusticol PyO3; C++ O3; chunk=96; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **8** | **AmpliCol** | **26.1** | **N/A** | **234** | **Fortran AmpliCol; with --library=use** |
| 8 | pyAmplicCol - JIT | 42.3 <span style="color:#bf8700">(x1.62)</span> | 1783(4) <span style="color:#cf222e">(x7.61)</span> | 1763(4) <span style="color:#cf222e">(x7.52)</span> | JIT; batch=16; last successful JIT refresh |
| 8 | pyAmplicCol - C++ | 1.06e3 <span style="color:#cf222e">(x40.69)</span> | 642.6(2.1) <span style="color:#cf222e">(x2.74)</span> | 491.8(1.4) <span style="color:#cf222e">(x2.10)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 8 | pyAmplicCol - C++ (O1) | 577 <span style="color:#cf222e">(x22.15)</span> | 997.2(3.3) <span style="color:#cf222e">(x4.26)</span> | 885.3(2.9) <span style="color:#cf222e">(x3.78)</span> | C++; O1; chunk=96; batch=16 |
| 8 | pyAmplicCol - ASM | 85.6 <span style="color:#cf222e">(x3.29)</span> | 1117.0(4.6) <span style="color:#cf222e">(x4.77)</span> | 1004.9(4.3) <span style="color:#cf222e">(x4.29)</span> | ASM; O3; chunk=96; batch=16 |
| 8 | pyAmplicCol - Rusticol | 1.02e3 <span style="color:#cf222e">(x39.15)</span> | 487.7(9) <span style="color:#cf222e">(x2.08)</span> | 464.4(6) <span style="color:#bf8700">(x1.98)</span> | Rusticol PyO3; C++ O3; chunk=96; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **9** | **AmpliCol** | **39.4** | **N/A** | **567** | **Fortran AmpliCol; with --library=use** |
| 9 | pyAmplicCol - JIT | 80.8 <span style="color:#cf222e">(x2.05)</span> | 6144(17) <span style="color:#cf222e">(x10.84)</span> | 6102(15) <span style="color:#cf222e">(x10.77)</span> | JIT; batch=16; last successful JIT refresh |
| 9 | pyAmplicCol - C++ | 2.83e3 <span style="color:#cf222e">(x71.80)</span> | 1504.9(11.7) <span style="color:#cf222e">(x2.66)</span> | 1252.7(9.6) <span style="color:#cf222e">(x2.21)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 9 | pyAmplicCol - ASM | 349 <span style="color:#cf222e">(x8.85)</span> | 2688.2(5.5) <span style="color:#cf222e">(x4.74)</span> | 2459.1(4.5) <span style="color:#cf222e">(x4.34)</span> | ASM; O3; chunk=64; batch=16 |
| 9 | pyAmplicCol - Rusticol | 2.83e3 <span style="color:#cf222e">(x71.80)</span> | 1164.8(9) <span style="color:#cf222e">(x2.06)</span> | 1120.9(1.0) <span style="color:#bf8700">(x1.98)</span> | Rusticol PyO3; C++ O3; chunk=96; cached NumPy input; reusable stage scratch; batch=128; repackaged from existing C++ O3 artifact in 103 s; two-stage save/load |

## Saved Evaluators

The CLI now supports reusable evaluator artifacts for the shared-current compiled DAG path and for the scalar zero-gluon process artifact:

- `--save-evaluator-dir PATH` materializes and saves an evaluator artifact, including generated C++ sources/shared libraries for compiled-complex backends, a manifest, and serialized Symbolica evaluator byte states where available.
- `--load-evaluator-dir PATH` loads the saved evaluator and skips Symbolica expression lowering and backend compilation.

This is intended for the C++ and ASM performance rows above. It lets runtime-only choices such as `--batch-size` be retimed without regenerating the evaluator. The n=5,6,7,8,9 C++ O3 rows in this table were generated with that two-stage save/load workflow and retimed from saved evaluators at batch size 128.

The rusticol rows use the same saved process directory format. The directory contains the shared-current DAG or scalar zero-gluon manifest, Symbolica stage/scalar evaluators, validation momenta, and the generated `API/python/check_standalone.py` runner. The n=9 rusticol row reuses the previously generated C++ O3 compiled evaluators and records the full C++ O3 generation cost in the table; writing the rusticol process wrapper around that artifact took 103 s. The rusticol runtime can load both saved JIT and compiled-complex stage evaluators, but the rows above use the compiled C++ evaluator libraries for `precision=16` because low-multiplicity probes found C++ O3 faster than saved JIT under rusticol. `evaluate_with_prec(..., 32)` and higher precisions load the serialized Symbolica evaluator states from the same process directory and return Python `Decimal` values. The wall timings above use cached NumPy input arrays in the benchmark wrapper and reusable Rust stage scratch buffers, matching the intended `rusticol.Runtime.evaluate/profile(momenta)` API where momenta are already supplied as arrays.

Process-generic Rusticol probes were also run for `u d~ > w+ g g`,
`d d~ > z z g`, pure-QCD channels, and multi-quark-line channels. These
validate through the generic schema-v2 artifact path rather than through a
`Z+gluons` fallback. Explicit dilepton examples such as `d d~ > e+ e- g g`
generate schema-v2 artifacts as well; in the current one-point comparison the
Fortran reference phase-space setup failed before pyAmpliCol evaluation, so
that process remains a generation/runtime-support check rather than a Fortran
ME-comparison row.

## Rusticol Artifact Load Times

These are direct `rusticol.Runtime.load(process_dir)` timings from the saved process directories used for the rusticol rows. `First load` records the first load in the timing process and includes dynamic-library loading/cache effects. `Warm load` is the mean of the remaining repeated loads in the same process.

| n | Artifact kind | Stages | First load [s] | Warm load [s] | Process directory |
|---:|---|---:|---:|---:|---|
| 1 | shared sweep | 2 | 0.00747 | 0.000204 | `pyAmpliCol/outputs/rusticol_backend_probe/cxx-o3-n1` |
| 2 | shared sweep | 3 | 0.00304 | 0.000345 | `pyAmpliCol/outputs/rusticol_backend_probe/cxx-o3-n2` |
| 3 | shared sweep | 4 | 0.00459 | 0.000657 | `pyAmpliCol/outputs/rusticol_goal_nle5/n3` |
| 4 | shared sweep | 5 | 0.00643 | 0.00141 | `pyAmpliCol/outputs/rusticol_goal_nle5/n4` |
| 5 | shared sweep | 6 | 0.0365 | 0.00879 | `pyAmpliCol/outputs/rusticol_goal_nle5/n5` |
| 6 | shared sweep | 7 | 0.100 | 0.0282 | `pyAmpliCol/outputs/rusticol_goal_n6_n9/n6` |
| 7 | shared sweep | 8 | 0.198 | 0.0758 | `pyAmpliCol/outputs/rusticol_goal_n6_n9/n7` |
| 8 | shared sweep | 9 | 0.669 | 0.293 | `pyAmpliCol/outputs/rusticol_goal_n6_n9/n8` |
| 9 | shared sweep | 10 | 2.52 | 1.12 | `pyAmpliCol/outputs/rusticol_goal_n6_n9/n9` |

## Reproduction Commands

The following commands reproduce the six setups for the representative case `n=3`, i.e. `d d~ > z g g g`. They intentionally do not use the watchdog wrapper or any benchmark collection helper.

AmpliCol:

The direct-library benchmark creates the generated Fortran amplitude library,
compiles it, and then calls `amp_lib:evaluate_amp` through the standalone
`amplicol_library_benchmark` driver. This bypasses integration/probe
bookkeeping and is the convention used for the refreshed low-multiplicity
Fortran reference rows.

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
make clean
make -j8 amplicol_generate
printf 'd d~ > z g g g\n' > processes.txt
./amplicol_generate --library=create --process=processes.txt
make -j8 amplicol_generate_library
make -j8 amplicol_library_benchmark
./amplicol_library_benchmark 1000000 1 1
```

pyAmplicCol - JIT generation and runtime:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol generate-process 'd d~ > z g g g' pyAmpliCol/outputs/reproduce/n3-jit --batch-size 16 --symbolica-evaluator-backend jit --symbolica-n-cores 10 --json
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol time-process pyAmpliCol/outputs/reproduce/n3-jit --batch-size 16 --target-runtime 10 --json
```

pyAmplicCol - C++ generation:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol generate-process 'd d~ > z g g g' pyAmpliCol/outputs/reproduce/n3-cxx-o3 --batch-size 128 --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset runtime-o3 --symbolica-n-cores 10 --symbolica-compiled-chunk-compile-workers 10 --json
```

pyAmplicCol - C++ runtime from saved process artifact:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol time-process pyAmpliCol/outputs/reproduce/n3-cxx-o3 --batch-size 128 --target-runtime 10 --json
```

pyAmplicCol - ASM generation:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol generate-process 'd d~ > z g g g' pyAmpliCol/outputs/reproduce/n3-asm --batch-size 16 --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset generation --symbolica-n-cores 10 --symbolica-compiled-chunk-compile-workers 10 --json
```

pyAmplicCol - ASM runtime from saved process artifact:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol time-process pyAmpliCol/outputs/reproduce/n3-asm --batch-size 16 --target-runtime 10 --json
```

pyAmplicCol - Rusticol default process-artifact generation:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol generate-process 'd d~ > z g g g' pyAmpliCol/outputs/reproduce/n3-rusticol --batch-size 128 --symbolica-compiled-preset runtime-o3 --symbolica-n-cores 10 --symbolica-compiled-chunk-compile-workers 10 --json
```

pyAmplicCol - Rusticol default runtime from saved process artifact:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol time-process pyAmpliCol/outputs/reproduce/n3-rusticol --batch-size 128 --target-runtime 10 --json
```

The generated process directory can also be checked without importing pyAmplicCol:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol/pyAmpliCol/outputs/reproduce/n3-rusticol
python API/python/check_standalone.py --precision 16 --profile
python API/python/check_standalone.py --precision 32 --profile
```

## Dependency Versions

```json
{
  "install_manifest": {
    "dependency_patches": [
      {
        "dependency": "symbolica",
        "path": "dependencies/patches/symbolica/0001-fix-complex-export-aarch64.patch"
      }
    ],
    "rusticol": {
      "installed": true,
      "requested": true,
      "source_path": "../rusticol",
      "source_rev": "92ac9e60932a2a0cf15a72b8f7d1e08e7fd5dccf",
      "usage": "Python, C++, and Fortran runtime for pyAmpliCol schema-v2 process artifacts"
    },
    "symbolica": {
      "installed": true,
      "requested": true,
      "source_path": "dependencies/symbolica",
      "source_ref": "dev",
      "source_rev": "e4167e767147ab8f3b4f039057c396c8fa961f6a",
      "source_url": "https://github.com/symbolica-dev/symbolica.git",
      "symjit_source_ref": "7fb09d1cb2a943c25a6fd71a208af44fcc6d813d",
      "symjit_source_rev": "7fb09d1cb2a943c25a6fd71a208af44fcc6d813d",
      "symjit_source_url": "https://github.com/siravan/symjit.git",
      "symjit_version": "2.19.3"
    },
    "symbolica_community": {
      "installed": true,
      "requested": true,
      "source_path": "dependencies/symbolica-community",
      "source_ref": "main",
      "source_rev": "4daa45e793a3775f13bee33c4ddd2ba7f22ba714",
      "source_url": "https://github.com/symbolica-dev/symbolica-community.git"
    }
  },
  "symbolica_local_versions": {
    "idenso": "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08",
    "spenso": "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08",
    "symbolica": "e4167e767147ab8f3b4f039057c396c8fa961f6a"
  },
  "symbolica_version": "2.1.0"
}
```
