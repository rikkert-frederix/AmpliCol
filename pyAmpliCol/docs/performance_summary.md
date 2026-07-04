# pyAmplicCol Performance Summary

Process family: `d d~ -> Z + n g`, for `n = 1,...,9`.

Successful pyAmplicCol rows validated against Fortran AmpliCol with max relative differences typically of order `1e-13`. The table below omits the validation column and focuses on generation time, full pyAmplicCol wall runtime, and evaluator-only runtime.

For pyAmplicCol rows, the multiplier in the `Gen` column is relative to AmpliCol generation time for the same `n`. Multipliers in `Wall` and `Eval` are relative to AmpliCol's per-point evaluator timing for the same `n`. Color convention: <span style="color:#1a7f37">green</span> is faster than AmpliCol, <span style="color:#bf8700">orange</span> is slower but below `x2.0`, and <span style="color:#cf222e">red</span> is `x2.0` or worse.

For C++ and ASM rows with multiple measured chunking options, this table reports the existing result with the best evaluator-only runtime, except where an additional O3 C++ row was explicitly added while retaining the previous lower-optimization row. The refreshed C++ O3 save/load rows use batch size 128 after a loaded-artifact batch-size scan; retained older rows use the batch size shown in their notes. All benchmark measurements summarized here were obtained under a watchdog limiting process-tree RAM usage to less than 30 GB.

The compiled pyAmplicCol DAG marks all known-real momentum-sum parameters as real with Symbolica `set_real_params(...)` before evaluator generation/compilation. In this hot path the couplings are embedded constants, so the relevant runtime input realness is the momentum sector.

| n | Setup | Gen [s] | Wall [us/pt] | Eval [us/pt] | Notes |
|---:|---|---:|---:|---:|---|
| **1** | **AmpliCol** | **2.22** | **N/A** | **0.734** | **Fortran AmpliCol** |
| 1 | pyAmplicCol - JIT | 0.016 <span style="color:#1a7f37">(x0.01)</span> | 13.34(29) <span style="color:#cf222e">(x18.17)</span> | 0.607(17) <span style="color:#1a7f37">(x0.83)</span> | JIT; batch=16 |
| 1 | pyAmplicCol - C++ | 0.774 <span style="color:#1a7f37">(x0.35)</span> | 13.01(12) <span style="color:#cf222e">(x17.72)</span> | 0.521(10) <span style="color:#1a7f37">(x0.71)</span> | C++; O3; chunk=None; batch=16 |
| 1 | pyAmplicCol - ASM | 0.560 <span style="color:#1a7f37">(x0.25)</span> | 14.03(30) <span style="color:#cf222e">(x19.12)</span> | 0.912(18) <span style="color:#bf8700">(x1.24)</span> | ASM; O3; chunk=None; batch=16 |
| **2** | **AmpliCol** | **2.20** | **N/A** | **1.64** | **Fortran AmpliCol** |
| 2 | pyAmplicCol - JIT | 0.066 <span style="color:#1a7f37">(x0.03)</span> | 16.15(15) <span style="color:#cf222e">(x9.85)</span> | 1.61(2) <span style="color:#1a7f37">(x0.98)</span> | JIT; batch=16 |
| 2 | pyAmplicCol - C++ | 1.86 <span style="color:#1a7f37">(x0.85)</span> | 16.05(14) <span style="color:#cf222e">(x9.79)</span> | 1.16(2) <span style="color:#1a7f37">(x0.70)</span> | C++; O3; chunk=None; batch=16 |
| 2 | pyAmplicCol - ASM | 0.940 <span style="color:#1a7f37">(x0.43)</span> | 16.48(20) <span style="color:#cf222e">(x10.05)</span> | 2.16(2) <span style="color:#bf8700">(x1.32)</span> | ASM; O3; chunk=None; batch=16 |
| **3** | **AmpliCol** | **2.31** | **N/A** | **4.55** | **Fortran AmpliCol** |
| 3 | pyAmplicCol - JIT | 0.281 <span style="color:#1a7f37">(x0.12)</span> | 22.67(42) <span style="color:#cf222e">(x4.98)</span> | 4.72(5) <span style="color:#bf8700">(x1.04)</span> | JIT; batch=16 |
| 3 | pyAmplicCol - C++ | 5.69 <span style="color:#cf222e">(x2.46)</span> | 20.90(40) <span style="color:#cf222e">(x4.59)</span> | 3.06(4) <span style="color:#1a7f37">(x0.67)</span> | C++; O3; chunk=None; batch=16 |
| 3 | pyAmplicCol - ASM | 1.52 <span style="color:#1a7f37">(x0.66)</span> | 23.87(45) <span style="color:#cf222e">(x5.25)</span> | 6.17(8) <span style="color:#bf8700">(x1.36)</span> | ASM; O3; chunk=None; batch=16 |
| **4** | **AmpliCol** | **2.55** | **N/A** | **13.2** | **Fortran AmpliCol** |
| 4 | pyAmplicCol - JIT | 0.252 <span style="color:#1a7f37">(x0.10)</span> | 49.53(42) <span style="color:#cf222e">(x3.75)</span> | 28.37(17) <span style="color:#cf222e">(x2.15)</span> | JIT; batch=16 |
| 4 | pyAmplicCol - C++ | 20.7 <span style="color:#cf222e">(x8.14)</span> | 27.78(13) <span style="color:#cf222e">(x2.10)</span> | 7.75(6) <span style="color:#1a7f37">(x0.59)</span> | C++; O3; chunk=None; batch=16 |
| 4 | pyAmplicCol - ASM | 2.04 <span style="color:#1a7f37">(x0.80)</span> | 36.30(54) <span style="color:#cf222e">(x2.75)</span> | 15.80(15) <span style="color:#bf8700">(x1.20)</span> | ASM; O3; chunk=None; batch=16 |
| **5** | **AmpliCol** | **3.22** | **N/A** | **34.3** | **Fortran AmpliCol** |
| 5 | pyAmplicCol - JIT | 0.327 <span style="color:#1a7f37">(x0.10)</span> | 373.6(2.8) <span style="color:#cf222e">(x10.89)</span> | 342.7(2.6) <span style="color:#cf222e">(x9.99)</span> | JIT; batch=16 |
| 5 | pyAmplicCol - C++ | 54.3 <span style="color:#cf222e">(x16.87)</span> | 45.06(48) <span style="color:#bf8700">(x1.31)</span> | 20.83(26) <span style="color:#1a7f37">(x0.61)</span> | C++; O3; chunk=64; batch=128; chunk compile workers=10; two-stage save/load |
| 5 | pyAmplicCol - C++ (O2) | 50.7 <span style="color:#cf222e">(x15.76)</span> | 51.04(87) <span style="color:#bf8700">(x1.49)</span> | 23.82(42) <span style="color:#1a7f37">(x0.69)</span> | C++; O2; chunk=64; batch=16 |
| 5 | pyAmplicCol - ASM | 10.2 <span style="color:#cf222e">(x3.16)</span> | 75.69(88) <span style="color:#cf222e">(x2.21)</span> | 48.74(72) <span style="color:#bf8700">(x1.42)</span> | ASM; O3; chunk=64; batch=16 |
| **6** | **AmpliCol** | **4.55** | **N/A** | **84.1** | **Fortran AmpliCol** |
| 6 | pyAmplicCol - JIT | 1.02 <span style="color:#1a7f37">(x0.22)</span> | 1435(10) <span style="color:#cf222e">(x17.06)</span> | 1380(8) <span style="color:#cf222e">(x16.41)</span> | JIT; batch=16 |
| 6 | pyAmplicCol - C++ | 151 <span style="color:#cf222e">(x33.23)</span> | 94.79(73) <span style="color:#bf8700">(x1.13)</span> | 56.36(49) <span style="color:#1a7f37">(x0.67)</span> | C++; O3; chunk=64; batch=128; chunk compile workers=10; two-stage save/load |
| 6 | pyAmplicCol - C++ (O2) | 140 <span style="color:#cf222e">(x30.74)</span> | 108.4(1.8) <span style="color:#bf8700">(x1.29)</span> | 68.5(1.6) <span style="color:#1a7f37">(x0.81)</span> | C++; O2; chunk=64; batch=16 |
| 6 | pyAmplicCol - ASM | 21.4 <span style="color:#cf222e">(x4.70)</span> | 183.7(3.3) <span style="color:#cf222e">(x2.18)</span> | 141.2(2.7) <span style="color:#bf8700">(x1.68)</span> | ASM; O3; chunk=64; batch=16 |
| **7** | **AmpliCol** | **9.38** | **N/A** | **207** | **Fortran AmpliCol** |
| 7 | pyAmplicCol - JIT | N/A | N/A | N/A | JIT backend bugged. |
| 7 | pyAmplicCol - C++ | 417 <span style="color:#cf222e">(x44.44)</span> | 259.7(1.4) <span style="color:#bf8700">(x1.25)</span> | 194.7(1.2) <span style="color:#1a7f37">(x0.94)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 7 | pyAmplicCol - C++ (O1) | 219 <span style="color:#cf222e">(x23.35)</span> | 421.6(2.4) <span style="color:#cf222e">(x2.04)</span> | 348.9(1.9) <span style="color:#bf8700">(x1.69)</span> | C++; O1; chunk=96; batch=16 |
| 7 | pyAmplicCol - ASM | 36.2 <span style="color:#cf222e">(x3.86)</span> | 384.4(2.8) <span style="color:#bf8700">(x1.86)</span> | 320.2(2.4) <span style="color:#bf8700">(x1.55)</span> | ASM; O3; chunk=96; batch=16 |
| **8** | **AmpliCol** | **25.0** | **N/A** | **567** | **Fortran AmpliCol** |
| 8 | pyAmplicCol - JIT | N/A | N/A | N/A | JIT backend bugged. |
| 8 | pyAmplicCol - C++ | 1.06e3 <span style="color:#cf222e">(x42.24)</span> | 642.6(2.1) <span style="color:#bf8700">(x1.13)</span> | 491.8(1.4) <span style="color:#1a7f37">(x0.87)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 8 | pyAmplicCol - C++ (O1) | 577 <span style="color:#cf222e">(x23.09)</span> | 997.2(3.3) <span style="color:#bf8700">(x1.76)</span> | 885.3(2.9) <span style="color:#bf8700">(x1.56)</span> | C++; O1; chunk=96; batch=16 |
| 8 | pyAmplicCol - ASM | 85.6 <span style="color:#cf222e">(x3.42)</span> | 1117.0(4.6) <span style="color:#bf8700">(x1.97)</span> | 1004.9(4.3) <span style="color:#bf8700">(x1.77)</span> | ASM; O3; chunk=96; batch=16 |
| **9** | **AmpliCol** | **48.0** | **N/A** | **2520** | **Fortran AmpliCol** |
| 9 | pyAmplicCol - JIT | N/A | N/A | N/A | JIT backend bugged. |
| 9 | pyAmplicCol - C++ | 2.83e3 <span style="color:#cf222e">(x58.83)</span> | 1504.9(11.7) <span style="color:#1a7f37">(x0.60)</span> | 1252.7(9.6) <span style="color:#1a7f37">(x0.50)</span> | C++; O3; chunk=96; batch=128; chunk compile workers=10; two-stage save/load |
| 9 | pyAmplicCol - ASM | 349 <span style="color:#cf222e">(x7.27)</span> | 2688.2(5.5) <span style="color:#bf8700">(x1.07)</span> | 2459.1(4.5) <span style="color:#1a7f37">(x0.98)</span> | ASM; O3; chunk=64; batch=16 |

## Saved Evaluators

The CLI now supports reusable compiled evaluator artifacts for the shared-current compiled DAG path:

- `--save-evaluator-dir PATH` materializes and saves a compiled evaluator artifact, including generated C++ sources, shared libraries, a manifest, and serialized Symbolica evaluator byte states where available.
- `--load-evaluator-dir PATH` loads the saved compiled evaluator and skips Symbolica expression lowering and C++ compilation.

This is intended for the C++ and ASM performance rows above. It lets runtime-only choices such as `--batch-size` be retimed without regenerating the evaluator. The n=5,6,7,8,9 C++ O3 rows in this table were generated with that two-stage save/load workflow and retimed from saved evaluators at batch size 128.

## Reproduction Commands

The following commands reproduce the four setups for the representative case `n=3`, i.e. `d d~ > z g g g`. They intentionally do not use the watchdog wrapper or any benchmark collection helper.

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
    "symbolica": {
      "installed": true,
      "requested": true,
      "source_path": "dependencies/symbolica",
      "source_rev": "e4167e767147ab8f3b4f039057c396c8fa961f6a"
    },
    "symbolica_community": {
      "installed": true,
      "requested": true,
      "source_path": "dependencies/symbolica-community",
      "source_rev": "4daa45e793a3775f13bee33c4ddd2ba7f22ba714"
    }
  },
  "symbolica_local_versions": {
    "idenso": "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08",
    "spenso": "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08",
    "symbolica": "e4167e767147ab8f3b4f039057c396c8fa961f6a",
    "vakint": "db79edc84f6a1580decbcc4ede7ea0b1c79d9a08"
  },
  "symbolica_version": "2.1.0"
}
```
