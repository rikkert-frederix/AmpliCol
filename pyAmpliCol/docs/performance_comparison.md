# pyamplicol Performance Comparison

> **Historical note:** this file records the pre-Rusticol staged-DAG
> performance exploration for the original `d d~ -> Z + n g` milestone. It is
> kept for provenance only and is not the current reproduction guide. Current
> production generation uses schema-v2 generic DAG process artifacts with the
> Rusticol runtime; use `performance_summary.md`, `process_coverage.md`, and
> the result-matrix documentation for the supported workflow.

Process family: `d d~ -> Z + n g`.

All heavy runs for this table were launched under `pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30`.
Validation points per multiplicity: `16`.
pyamplicol timing repetitions: `20`.
Fortran AmpliCol timing sample: `10000`.
pyamplicol batch size: `16` for the original scan; refreshed two-stage C++ O3 rows use batch size `128`.

## Modes

- A: Fortran AmpliCol generated Fortran library, timed with the quiet `--amplicol_probe` timing path.
- X: pyamplicol shared-helicity-current D-mode with Symbolica `compiled-complex` generic C++ output and the `runtime` preset.
- Y: pyamplicol shared-helicity-current D-mode with Symbolica `compiled-complex` emitted assembly output and the `generation` preset.
- Z: pyamplicol shared-helicity-current D-mode with Symbolica JIT compilation. This is the generation-time-oriented mode.

For pyamplicol rows, `wall` is full `evaluate_matrix_elements_many` wall time per phase-space point and `eval` is time spent inside Symbolica evaluator calls only. For Fortran AmpliCol, `eval` is the printed `amplitude evaluation` timing per point.

## Table

| n | mode | status | gen [s] | wall [us/pt] | eval [us/pt] | max rel diff | notes |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | A | ok | 2.218 | n/a | 0.734 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 1 | X | ok | 0.774 | 13.009 +/- 0.116 | 0.521 +/- 0.010 | 2.669e-14 | compiled-complex, preset=runtime, asm=none, O3, chunk=None |
| 1 | Y | ok | 0.560 | 14.032 +/- 0.301 | 0.912 +/- 0.018 | 1.439e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 1 | Z | ok | 0.016 | 13.340 +/- 0.286 | 0.607 +/- 0.017 | 2.614e-14 | jit |
| 2 | A | ok | 2.196 | n/a | 1.638 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 2 | X | ok | 1.864 | 16.048 +/- 0.142 | 1.155 +/- 0.017 | 2.040e-13 | compiled-complex, preset=runtime, asm=none, O3, chunk=None |
| 2 | Y | ok | 0.940 | 16.481 +/- 0.195 | 2.163 +/- 0.024 | 2.036e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 2 | Z | ok | 0.066 | 16.151 +/- 0.148 | 1.612 +/- 0.015 | 1.952e-13 | jit |
| 3 | A | ok | 2.312 | n/a | 4.547 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 3 | X | ok | 5.687 | 20.902 +/- 0.396 | 3.058 +/- 0.035 | 1.363e-13 | compiled-complex, preset=runtime, asm=none, O3, chunk=None |
| 3 | Y | ok | 1.516 | 23.869 +/- 0.446 | 6.170 +/- 0.081 | 1.361e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 3 | Z | ok | 0.281 | 22.668 +/- 0.421 | 4.715 +/- 0.047 | 7.662e-13 | jit |
| 4 | A | ok | 2.550 | n/a | 13.232 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 4 | X | ok | 20.748 | 27.782 +/- 0.131 | 7.752 +/- 0.059 | 2.028e-12 | compiled-complex, preset=runtime, asm=none, O3, chunk=None |
| 4 | Y | ok | 2.039 | 36.298 +/- 0.541 | 15.795 +/- 0.153 | 2.027e-12 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 4 | Z | ok | 0.252 | 49.525 +/- 0.424 | 28.374 +/- 0.166 | 1.712e-12 | jit |
| 5 | A | ok | 3.215 | n/a | 34.284 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 5 | X-O3 | ok | 54.313 | 45.058 +/- 0.480 | 20.830 +/- 0.261 | 6.535e-13 | compiled-complex, preset=runtime-o3, asm=none, O3, chunk=64, batch=128, two-stage save/load, chunk compile workers=10 |
| 5 | X | ok | 50.744 | 51.039 +/- 0.867 | 23.823 +/- 0.421 | 6.535e-13 | compiled-complex, preset=runtime, asm=none, O2, chunk=64 |
| 5 | Y | ok | 3.185 | 91.913 +/- 0.488 | 67.219 +/- 0.414 | 6.545e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 5 | Z | ok | 0.327 | 373.635 +/- 2.778 | 342.747 +/- 2.644 | 4.848e-13 | jit |
| 6 | A | ok | 4.545 | n/a | 84.122 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 6 | X-O3 | ok | 151.200 | 94.794 +/- 0.732 | 56.359 +/- 0.487 | 2.481e-13 | compiled-complex, preset=runtime-o3, asm=none, O3, chunk=64, batch=128, two-stage save/load, chunk compile workers=10 |
| 6 | X | ok | 139.851 | 108.391 +/- 1.819 | 68.481 +/- 1.558 | 2.481e-13 | compiled-complex, preset=runtime, asm=none, O2, chunk=64 |
| 6 | Y | ok | 6.203 | 228.555 +/- 2.586 | 190.724 +/- 2.238 | 2.497e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 6 | Z | ok | 1.015 | 1435.117 +/- 9.545 | 1379.992 +/- 8.463 | 2.492e-13 | jit |
| 7 | A | ok | 9.378 | n/a | 206.639 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 7 | X-O3 | ok | 416.885 | 259.681 +/- 1.396 | 194.669 +/- 1.169 | 7.590e-13 | compiled-complex, preset=runtime-o3, asm=none, O3, chunk=96, batch=128, two-stage save/load, chunk compile workers=10 |
| 7 | X | ok | 218.982 | 421.601 +/- 2.361 | 348.863 +/- 1.872 | 7.590e-13 | compiled-complex, preset=runtime, asm=none, O1, chunk=96 |
| 7 | Y | ok | 14.197 | 629.734 +/- 5.687 | 568.637 +/- 4.896 | 7.581e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 7 | Z | error | n/a | n/a | n/a | n/a | thread '<unnamed>' (23023180) panicked at /Users/vjhirsch/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/symjit-2.18.6/rust/arm/vector.rs:795:43: |
| 8 | A | ok | 25.034 | n/a | 567.249 +/- 0.050 | 0.000e+00 | Fortran AmpliCol |
| 8 | X-O3 | ok | 1055.509 | 642.578 +/- 2.078 | 491.811 +/- 1.444 | 2.151e-13 | compiled-complex, preset=runtime-o3, asm=none, O3, chunk=96, batch=128, two-stage save/load, chunk compile workers=10 |
| 8 | X | ok | 577.208 | 997.182 +/- 3.331 | 885.275 +/- 2.860 | 2.151e-13 | compiled-complex, preset=runtime, asm=none, O1, chunk=96 |
| 8 | Y | ok | 45.191 | 3404.697 +/- 85.266 | 3300.954 +/- 84.972 | 2.136e-13 | compiled-complex, preset=generation, asm=default, O3, chunk=None |
| 8 | Z | error | n/a | n/a | n/a | n/a | thread '<unnamed>' (23044695) panicked at /Users/vjhirsch/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/symjit-2.18.6/rust/arm/vector.rs:795:43: |
| 9 | A | ok | 48.044 | n/a | 2520 | 0.000e+00 | Fortran AmpliCol |
| 9 | X-O3 | ok | 2826.345 | 1504.867 +/- 11.715 | 1252.733 +/- 9.575 | 1.319e-13 | compiled-complex, preset=runtime-o3, asm=none, O3, chunk=96, batch=128, two-stage save/load, chunk compile workers=10 |
| 9 | Z | error | n/a | n/a | n/a | n/a | thread '<unnamed>' (23054874) panicked at /Users/vjhirsch/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/symjit-2.18.6/rust/arm/vector.rs:795:43: |

## Chunking Study

This table compares Symbolica `compiled-complex` generated-code evaluators with and without pyamplicol output chunking. `asm` means `--symbolica-compiled-inline-asm default`; `c++` means `--symbolica-compiled-inline-asm none`. The C++ no-chunk case is bypassed for `n >= 8` because it is not useful operationally at these multiplicities.

| n | backend | chunk | opt | gen [s] | wall [us/pt] | eval [us/pt] | notes |
|---:|---|---:|---:|---:|---:|---:|---|
| 5 | asm | none | O3 | 3.068 | 96.290 +/- 2.218 | 68.862 +/- 1.369 | fresh paired chunking scan |
| 5 | asm | 64 | O3 | 10.184 | 75.688 +/- 0.876 | 48.741 +/- 0.721 | fresh paired chunking scan |
| 5 | c++ | none | O2 | 69.496 | 73.512 +/- 0.334 | 50.473 +/- 0.291 | fresh paired chunking scan |
| 5 | c++ | 64 | O3 | 54.313 | 45.058 +/- 0.480 | 20.830 +/- 0.261 | runtime-o3 preset, batch=128, two-stage save/load, chunk compile workers=10 |
| 5 | c++ | 64 | O2 | 50.744 | 51.039 +/- 0.867 | 23.823 +/- 0.421 | main benchmark X row |
| 6 | asm | none | O3 | 6.089 | 219.508 +/- 1.937 | 183.203 +/- 1.665 | fresh paired chunking scan |
| 6 | asm | 64 | O3 | 21.385 | 183.688 +/- 3.257 | 141.223 +/- 2.733 | fresh paired chunking scan |
| 6 | c++ | none | O2 | 295.517 | 204.467 +/- 3.656 | 165.126 +/- 2.802 | fresh paired chunking scan |
| 6 | c++ | 64 | O3 | 151.200 | 94.794 +/- 0.732 | 56.359 +/- 0.487 | runtime-o3 preset, batch=128, two-stage save/load, chunk compile workers=10 |
| 6 | c++ | 64 | O2 | 139.851 | 108.391 +/- 1.819 | 68.481 +/- 1.558 | main benchmark X row |
| 7 | asm | none | O3 | 14.315 | 654.814 +/- 5.166 | 590.260 +/- 4.070 | fresh paired chunking scan |
| 7 | asm | 96 | O3 | 36.174 | 384.422 +/- 2.829 | 320.193 +/- 2.427 | fresh paired chunking scan |
| 7 | c++ | none | O1 | 298.615 | 623.370 +/- 4.263 | 564.065 +/- 3.724 | fresh paired chunking scan |
| 7 | c++ | 96 | O3 | 416.885 | 259.681 +/- 1.396 | 194.669 +/- 1.169 | runtime-o3 preset, batch=128, two-stage save/load, chunk compile workers=10 |
| 7 | c++ | 96 | O1 | 218.982 | 421.601 +/- 2.361 | 348.863 +/- 1.872 | main benchmark X row |
| 8 | asm | none | O3 | 41.433 | 1443.478 +/- 5.654 | 1375.885 +/- 5.250 | fresh paired chunking scan |
| 8 | asm | 96 | O3 | 85.615 | 1117.026 +/- 4.573 | 1004.861 +/- 4.331 | fresh paired chunking scan |
| 8 | c++ | none | O1 | n/a | n/a | n/a | bypassed per request |
| 8 | c++ | 96 | O3 | 1055.509 | 642.578 +/- 2.078 | 491.811 +/- 1.444 | runtime-o3 preset, batch=128, two-stage save/load, chunk compile workers=10 |
| 8 | c++ | 96 | O1 | 577.208 | 997.182 +/- 3.331 | 885.275 +/- 2.860 | main benchmark X row |
| 9 | asm | 64 | O3 | 348.792 | 2688.228 +/- 5.530 | 2459.077 +/- 4.532 | fresh assembly chunk-size scan |
| 9 | asm | 96 | O3 | 297.279 | 2881.639 +/- 5.260 | 2666.879 +/- 4.691 | fresh assembly chunk-size scan |
| 9 | asm | 128 | O3 | 270.360 | 3126.726 +/- 7.558 | 2908.961 +/- 6.510 | fresh assembly chunk-size scan |
| 9 | c++ | 96 | O3 | 2826.345 | 1504.867 +/- 11.715 | 1252.733 +/- 9.575 | runtime-o3 preset, batch=128, two-stage save/load, chunk compile workers=10 |

Chunking improves runtime for the measured assembly rows as well as the generic C++ rows, although it usually increases generation time for assembly. The refreshed O3 generic C++ backend is faster than assembly for `n=5,6,7,8,9`. For the measured `n=9` assembly scan, chunk 64 gives the best runtime while chunk 128 gives the fastest generation.

## Reproduction Commands

Set up the shell once:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
export PYTHONPATH=pyAmpliCol/src
PY=pyAmpliCol/dependencies/.venv/bin/python
WATCH='pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30 --'
```

The full table can be regenerated with:

```bash
$WATCH $PY pyAmpliCol/scripts/collect_performance_comparison.py \
  --min-gluons 1 --max-gluons 8 --extra-z-gluons 9 \
  --points 16 --repetitions 20 \
  --legacy-timing 10000 --batch-size 16
```

For an individual process, define for example:

```bash
PROCESS='d d~ > z g g g'
```

Fortran AmpliCol validation/timing:

```bash
$WATCH $PY -m pyamplicol compare-amplicol "$PROCESS" \
  --amplicol-probe --points 16 --timing 10000 --json
```

Mode X profile:

```bash
$WATCH $PY -m pyamplicol profile-dag-evaluator "$PROCESS" \
  --points 16 --repetitions 20 --batch-size 16 \
  --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset runtime --json
```

Mode Y profile:

```bash
$WATCH $PY -m pyamplicol profile-dag-evaluator "$PROCESS" \
  --points 16 --repetitions 20 --batch-size 16 \
  --symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset generation --json
```

Mode Z profile:

```bash
$WATCH $PY -m pyamplicol profile-dag-evaluator "$PROCESS" \
  --points 16 --repetitions 20 --batch-size 16 \
  --symbolica-evaluator-backend jit --json
```

The `n=9` assembly chunking scan was run with:

```bash
PROCESS='d d~ > z g g g g g g g g g'
for CHUNK in 64 96 128; do
  $WATCH $PY -m pyamplicol profile-dag-evaluator "$PROCESS" \
    --points 16 --repetitions 20 --batch-size 16 \
    --symbolica-evaluator-backend compiled-complex \
    --symbolica-compiled-preset generation \
    --symbolica-compiled-output-chunk-size "$CHUNK" --json
done
```

The mode-specific validation command is the same as the legacy command, with `--runtime-backend dag` and the mode-specific Symbolica options added.

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
