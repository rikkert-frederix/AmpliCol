# pyAmplicCol Performance Summary

Process family: `d d~ -> Z + n g`, for `n = 1,...,9`.

Successful pyAmplicCol rows validated against Fortran AmpliCol with max relative differences typically of order `1e-13`. The rusticol rows shown for `n = 1,...,9` were validated against Fortran AmpliCol/native fixed probes, with max relative differences no larger than `2.1e-12`. The table below omits the validation column and focuses on generation time, full pyAmplicCol wall runtime, and evaluator-only runtime.

For pyAmplicCol rows, the multiplier in the `Gen` column is relative to AmpliCol generation time for the same `n`. Multipliers in `Wall` and `Eval` are relative to AmpliCol's per-point evaluator timing for the same `n`. Color convention: <span style="color:#1a7f37">green</span> is faster than AmpliCol, <span style="color:#bf8700">orange</span> is slower but below `x2.0`, and <span style="color:#cf222e">red</span> is `x2.0` or worse.

For C++ and ASM rows with multiple measured chunking options, this table reports the existing result with the best evaluator-only runtime, except where an additional O3 C++ row was explicitly added while retaining the previous lower-optimization row. The refreshed C++ O3 save/load rows use batch size 128 after a loaded-artifact batch-size scan; retained older rows use the batch size shown in their notes. The rusticol rows use C++ O3 compiled-complex artifacts; rusticol can also load saved JIT artifacts, but direct rusticol probes at low multiplicity confirmed C++ O3 remains faster. All benchmark measurements summarized here were obtained under a watchdog limiting process-tree RAM usage to less than 30 GB.

The compiled pyAmplicCol DAG marks all known-real momentum-sum parameters as real with Symbolica `set_real_params(...)` before evaluator generation/compilation. In this hot path the couplings are embedded constants, so the relevant runtime input realness is the momentum sector.

| n | Setup | Gen [s] | Wall [us/pt] | Eval [us/pt] | Notes |
|---:|---|---:|---:|---:|---|
| **1** | **AmpliCol** | **2.22** | **N/A** | **0.734** | **Fortran AmpliCol** |
| 1 | pyAmplicCol - JIT | 0.016 <span style="color:#1a7f37">(x0.01)</span> | 13.34(29) <span style="color:#cf222e">(x18.17)</span> | 0.607(17) <span style="color:#1a7f37">(x0.83)</span> | JIT; batch=16 |
| 1 | pyAmplicCol - C++ | 0.774 <span style="color:#1a7f37">(x0.35)</span> | 13.01(12) <span style="color:#cf222e">(x17.72)</span> | 0.521(10) <span style="color:#1a7f37">(x0.71)</span> | C++; O3; chunk=None; batch=16 |
| 1 | pyAmplicCol - ASM | 0.560 <span style="color:#1a7f37">(x0.25)</span> | 14.03(30) <span style="color:#cf222e">(x19.12)</span> | 0.912(18) <span style="color:#bf8700">(x1.24)</span> | ASM; O3; chunk=None; batch=16 |
| 1 | pyAmplicCol - Rusticol | 0.984 <span style="color:#1a7f37">(x0.44)</span> | 0.412(1) <span style="color:#1a7f37">(x0.56)</span> | 0.1938(2) <span style="color:#1a7f37">(x0.26)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; single-chunk direct output; batch=1024 |
| **2** | **AmpliCol** | **2.20** | **N/A** | **1.64** | **Fortran AmpliCol** |
| 2 | pyAmplicCol - JIT | 0.066 <span style="color:#1a7f37">(x0.03)</span> | 16.15(15) <span style="color:#cf222e">(x9.85)</span> | 1.61(2) <span style="color:#1a7f37">(x0.98)</span> | JIT; batch=16 |
| 2 | pyAmplicCol - C++ | 1.86 <span style="color:#1a7f37">(x0.85)</span> | 16.05(14) <span style="color:#cf222e">(x9.79)</span> | 1.16(2) <span style="color:#1a7f37">(x0.70)</span> | C++; O3; chunk=None; batch=16 |
| 2 | pyAmplicCol - ASM | 0.940 <span style="color:#1a7f37">(x0.43)</span> | 16.48(20) <span style="color:#cf222e">(x10.05)</span> | 2.16(2) <span style="color:#bf8700">(x1.32)</span> | ASM; O3; chunk=None; batch=16 |
| 2 | pyAmplicCol - Rusticol | 2.00 <span style="color:#1a7f37">(x0.91)</span> | 0.984(1) <span style="color:#1a7f37">(x0.60)</span> | 0.653(1) <span style="color:#1a7f37">(x0.40)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; single-chunk direct output; batch=1024 |
| **3** | **AmpliCol** | **2.31** | **N/A** | **4.55** | **Fortran AmpliCol** |
| 3 | pyAmplicCol - JIT | 0.281 <span style="color:#1a7f37">(x0.12)</span> | 22.67(42) <span style="color:#cf222e">(x4.98)</span> | 4.72(5) <span style="color:#bf8700">(x1.04)</span> | JIT; batch=16 |
| 3 | pyAmplicCol - C++ | 5.69 <span style="color:#cf222e">(x2.46)</span> | 20.90(40) <span style="color:#cf222e">(x4.59)</span> | 3.06(4) <span style="color:#1a7f37">(x0.67)</span> | C++; O3; chunk=None; batch=16 |
| 3 | pyAmplicCol - ASM | 1.52 <span style="color:#1a7f37">(x0.66)</span> | 23.87(45) <span style="color:#cf222e">(x5.25)</span> | 6.17(8) <span style="color:#bf8700">(x1.36)</span> | ASM; O3; chunk=None; batch=16 |
| 3 | pyAmplicCol - Rusticol | 5.48 <span style="color:#cf222e">(x2.37)</span> | 2.783(13) <span style="color:#1a7f37">(x0.61)</span> | 2.177(10) <span style="color:#1a7f37">(x0.48)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; batch=128 |
| **4** | **AmpliCol** | **2.55** | **N/A** | **13.2** | **Fortran AmpliCol** |
| 4 | pyAmplicCol - JIT | 0.252 <span style="color:#1a7f37">(x0.10)</span> | 49.53(42) <span style="color:#cf222e">(x3.75)</span> | 28.37(17) <span style="color:#cf222e">(x2.15)</span> | JIT; batch=16 |
| 4 | pyAmplicCol - C++ | 20.7 <span style="color:#cf222e">(x8.14)</span> | 27.78(13) <span style="color:#cf222e">(x2.10)</span> | 7.75(6) <span style="color:#1a7f37">(x0.59)</span> | C++; O3; chunk=None; batch=16 |
| 4 | pyAmplicCol - ASM | 2.04 <span style="color:#1a7f37">(x0.80)</span> | 36.30(54) <span style="color:#cf222e">(x2.75)</span> | 15.80(15) <span style="color:#bf8700">(x1.20)</span> | ASM; O3; chunk=None; batch=16 |
| 4 | pyAmplicCol - Rusticol | 19.9 <span style="color:#cf222e">(x7.81)</span> | 7.23(7) <span style="color:#1a7f37">(x0.55)</span> | 6.10(6) <span style="color:#1a7f37">(x0.46)</span> | Rusticol PyO3; C++ O3 process artifact; cached NumPy input; reusable stage scratch; batch=128 |
| **5** | **AmpliCol** | **3.22** | **N/A** | **34.3** | **Fortran AmpliCol** |
| 5 | pyAmplicCol - JIT | 0.327 <span style="color:#1a7f37">(x0.10)</span> | 373.6(2.8) <span style="color:#cf222e">(x10.89)</span> | 342.7(2.6) <span style="color:#cf222e">(x9.99)</span> | JIT; batch=16 |
| 5 | pyAmplicCol - C++ | 54.3 <span style="color:#cf222e">(x16.87)</span> | 45.06(48) <span style="color:#bf8700">(x1.31)</span> | 20.83(26) <span style="color:#1a7f37">(x0.61)</span> | C++; O3; chunk=64; batch=128; chunk compile workers=10; two-stage save/load |
| 5 | pyAmplicCol - C++ (O2) | 50.7 <span style="color:#cf222e">(x15.76)</span> | 51.04(87) <span style="color:#bf8700">(x1.49)</span> | 23.82(42) <span style="color:#1a7f37">(x0.69)</span> | C++; O2; chunk=64; batch=16 |
| 5 | pyAmplicCol - ASM | 10.2 <span style="color:#cf222e">(x3.16)</span> | 75.69(88) <span style="color:#cf222e">(x2.21)</span> | 48.74(72) <span style="color:#bf8700">(x1.42)</span> | ASM; O3; chunk=64; batch=16 |
| 5 | pyAmplicCol - Rusticol | 50.2 <span style="color:#cf222e">(x15.60)</span> | 19.63(28) <span style="color:#1a7f37">(x0.57)</span> | 17.50(25) <span style="color:#1a7f37">(x0.51)</span> | Rusticol PyO3; C++ O3; chunk=64; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **6** | **AmpliCol** | **4.55** | **N/A** | **84.1** | **Fortran AmpliCol** |
| 6 | pyAmplicCol - JIT | 1.02 <span style="color:#1a7f37">(x0.22)</span> | 1435(10) <span style="color:#cf222e">(x17.06)</span> | 1380(8) <span style="color:#cf222e">(x16.41)</span> | JIT; batch=16 |
| 6 | pyAmplicCol - C++ | 151 <span style="color:#cf222e">(x33.23)</span> | 94.79(73) <span style="color:#bf8700">(x1.13)</span> | 56.36(49) <span style="color:#1a7f37">(x0.67)</span> | C++; O3; chunk=64; batch=128; chunk compile workers=10; two-stage save/load |
| 6 | pyAmplicCol - C++ (O2) | 140 <span style="color:#cf222e">(x30.74)</span> | 108.4(1.8) <span style="color:#bf8700">(x1.29)</span> | 68.5(1.6) <span style="color:#1a7f37">(x0.81)</span> | C++; O2; chunk=64; batch=16 |
| 6 | pyAmplicCol - ASM | 21.4 <span style="color:#cf222e">(x4.70)</span> | 183.7(3.3) <span style="color:#cf222e">(x2.18)</span> | 141.2(2.7) <span style="color:#bf8700">(x1.68)</span> | ASM; O3; chunk=64; batch=16 |
| 6 | pyAmplicCol - Rusticol | 141 <span style="color:#cf222e">(x31.08)</span> | 55.63(26) <span style="color:#1a7f37">(x0.66)</span> | 50.12(23) <span style="color:#1a7f37">(x0.60)</span> | Rusticol PyO3; C++ O3; chunk=64; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **7** | **AmpliCol** | **9.38** | **N/A** | **207** | **Fortran AmpliCol** |
| 7 | pyAmplicCol - JIT | N/A | N/A | N/A | JIT backend bugged. |
| 7 | pyAmplicCol - C++ | 417 <span style="color:#cf222e">(x44.44)</span> | 259.7(1.4) <span style="color:#bf8700">(x1.25)</span> | 194.7(1.2) <span style="color:#1a7f37">(x0.94)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 7 | pyAmplicCol - C++ (O1) | 219 <span style="color:#cf222e">(x23.35)</span> | 421.6(2.4) <span style="color:#cf222e">(x2.04)</span> | 348.9(1.9) <span style="color:#bf8700">(x1.69)</span> | C++; O1; chunk=96; batch=16 |
| 7 | pyAmplicCol - ASM | 36.2 <span style="color:#cf222e">(x3.86)</span> | 384.4(2.8) <span style="color:#bf8700">(x1.86)</span> | 320.2(2.4) <span style="color:#bf8700">(x1.55)</span> | ASM; O3; chunk=96; batch=16 |
| 7 | pyAmplicCol - Rusticol | 391 <span style="color:#cf222e">(x41.71)</span> | 197.8(4) <span style="color:#1a7f37">(x0.96)</span> | 186.1(4) <span style="color:#1a7f37">(x0.90)</span> | Rusticol PyO3; C++ O3; chunk=96; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **8** | **AmpliCol** | **25.0** | **N/A** | **567** | **Fortran AmpliCol** |
| 8 | pyAmplicCol - JIT | N/A | N/A | N/A | JIT backend bugged. |
| 8 | pyAmplicCol - C++ | 1.06e3 <span style="color:#cf222e">(x42.24)</span> | 642.6(2.1) <span style="color:#bf8700">(x1.13)</span> | 491.8(1.4) <span style="color:#1a7f37">(x0.87)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 8 | pyAmplicCol - C++ (O1) | 577 <span style="color:#cf222e">(x23.09)</span> | 997.2(3.3) <span style="color:#bf8700">(x1.76)</span> | 885.3(2.9) <span style="color:#bf8700">(x1.56)</span> | C++; O1; chunk=96; batch=16 |
| 8 | pyAmplicCol - ASM | 85.6 <span style="color:#cf222e">(x3.42)</span> | 1117.0(4.6) <span style="color:#bf8700">(x1.97)</span> | 1004.9(4.3) <span style="color:#bf8700">(x1.77)</span> | ASM; O3; chunk=96; batch=16 |
| 8 | pyAmplicCol - Rusticol | 1.02e3 <span style="color:#cf222e">(x40.90)</span> | 487.7(9) <span style="color:#1a7f37">(x0.86)</span> | 464.4(6) <span style="color:#1a7f37">(x0.82)</span> | Rusticol PyO3; C++ O3; chunk=96; cached NumPy input; reusable stage scratch; batch=128; two-stage save/load |
| **9** | **AmpliCol** | **48.0** | **N/A** | **2520** | **Fortran AmpliCol** |
| 9 | pyAmplicCol - JIT | N/A | N/A | N/A | JIT backend bugged. |
| 9 | pyAmplicCol - C++ | 2.83e3 <span style="color:#cf222e">(x58.83)</span> | 1504.9(11.7) <span style="color:#1a7f37">(x0.60)</span> | 1252.7(9.6) <span style="color:#1a7f37">(x0.50)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 9 | pyAmplicCol - ASM | 349 <span style="color:#cf222e">(x7.27)</span> | 2688.2(5.5) <span style="color:#bf8700">(x1.07)</span> | 2459.1(4.5) <span style="color:#1a7f37">(x0.98)</span> | ASM; O3; chunk=64; batch=16 |
| 9 | pyAmplicCol - Rusticol | 2.83e3 <span style="color:#cf222e">(x58.83)</span> | 1164.8(9) <span style="color:#1a7f37">(x0.46)</span> | 1120.9(1.0) <span style="color:#1a7f37">(x0.44)</span> | Rusticol PyO3; C++ O3; chunk=96; cached NumPy input; reusable stage scratch; batch=128; repackaged from existing C++ O3 artifact in 103 s; two-stage save/load |

## Saved Evaluators

The CLI now supports reusable evaluator artifacts for the shared-current compiled DAG path and for the scalar zero-gluon process artifact:

- `--save-evaluator-dir PATH` materializes and saves an evaluator artifact, including generated C++ sources/shared libraries for compiled-complex backends, a manifest, and serialized Symbolica evaluator byte states where available.
- `--load-evaluator-dir PATH` loads the saved evaluator and skips Symbolica expression lowering and backend compilation.

This is intended for the C++ and ASM performance rows above. It lets runtime-only choices such as `--batch-size` be retimed without regenerating the evaluator. The n=5,6,7,8,9 C++ O3 rows in this table were generated with that two-stage save/load workflow and retimed from saved evaluators at batch size 128.

The rusticol rows use the same saved process directory format. The directory contains the shared-current DAG or scalar zero-gluon manifest, Symbolica stage/scalar evaluators, validation momenta, and a standalone `check_standalone.py` script. The n=9 rusticol row reuses the previously generated C++ O3 compiled evaluators and records the full C++ O3 generation cost in the table; writing the rusticol process wrapper around that artifact took 103 s. The rusticol runtime can load both saved JIT and compiled-complex stage evaluators, but the rows above use the compiled C++ evaluator libraries for `precision=16` because low-multiplicity probes found C++ O3 faster than saved JIT under rusticol. `evaluate_with_prec(..., 32)` and higher precisions load the serialized Symbolica evaluator states from the same process directory and return Python `Decimal` values. The wall timings above use cached NumPy input arrays in the benchmark wrapper and reusable Rust stage scratch buffers, matching the intended `rusticol.Runtime.evaluate/profile(momenta)` API where momenta are already supplied as arrays.

Process-generic rusticol probes were also run for `u d~ > w+ g g`, `d d~ > e+ e- g g`, and `d d~ > z z g`. The current native graph/lowering route reports these as unsupported with explicit JSON diagnostics (`no native graph available ...`) rather than silently falling back to a Z+gluons-specific path.

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

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol compare-amplicol 'd d~ > z g g g' --amplicol-probe --points 10000 --timing 1 --json
```

pyAmplicCol - JIT:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --points 16 --repetitions 20 --batch-size 16 --symbolica-evaluator-backend jit --symbolica-n-cores 10 --json
```

pyAmplicCol - C++ generation:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --points 128 --repetitions 1 --batch-size 128 --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset runtime-o3 --symbolica-n-cores 10 --symbolica-compiled-chunk-compile-workers 10 --save-evaluator-dir pyAmpliCol/outputs/reproduce/n3-cxx-o3 --json
```

pyAmplicCol - C++ runtime from saved evaluator:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --points 128 --repetitions 10 --batch-size 128 --load-evaluator-dir pyAmpliCol/outputs/reproduce/n3-cxx-o3 --json
```

pyAmplicCol - ASM generation:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --points 16 --repetitions 1 --batch-size 16 --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset generation --symbolica-n-cores 10 --symbolica-compiled-chunk-compile-workers 10 --save-evaluator-dir pyAmpliCol/outputs/reproduce/n3-asm --json
```

pyAmplicCol - ASM runtime from saved evaluator:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --points 16 --repetitions 20 --batch-size 16 --load-evaluator-dir pyAmpliCol/outputs/reproduce/n3-asm --json
```

pyAmplicCol - Rusticol process-artifact generation:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --runtime-backend rusticol --generate-only --batch-size 128 --symbolica-compiled-preset runtime-o3 --symbolica-n-cores 10 --symbolica-compiled-chunk-compile-workers 10 --save-evaluator-dir pyAmpliCol/outputs/reproduce/n3-rusticol --json
```

pyAmplicCol - Rusticol runtime from saved process artifact:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator 'd d~ > z g g g' --runtime-backend rusticol --points 128 --repetitions 10 --batch-size 128 --load-evaluator-dir pyAmpliCol/outputs/reproduce/n3-rusticol --json
```

The generated process directory can also be checked without importing pyAmplicCol:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol/pyAmpliCol/outputs/reproduce/n3-rusticol
python check_standalone.py --precision 16 --profile
python check_standalone.py --precision 32 --profile
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
      "source_rev": "142a48ec621ba1c54b5ed1eacec7a8260067bf0c",
      "usage": "PyO3 runtime for pyAmpliCol eager-DAG process artifacts"
    },
    "symbolica": {
      "installed": true,
      "requested": true,
      "source_path": "dependencies/symbolica",
      "source_ref": "pyamplicol-dev-base",
      "source_rev": "fe3084804c8677b1341e84428d65e0e713b5901c",
      "source_url": "https://github.com/ValentinHirschi/symbolica_mod.git"
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
    "symbolica": "fe3084804c8677b1341e84428d65e0e713b5901c"
  },
  "symbolica_version": "2.1.0"
}
```
