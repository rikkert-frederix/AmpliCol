# pyAmpliCol Compiled DAG Performance Summary

This document tracks the `--compiled-dag-evaluator` mode. The intended final
table mirrors `performance_summary.md`, but for the alias-backed single
multi-output evaluator route.

Compiled-DAG generated-code profiles should use no final-output chunking by
default. Chunking was tested and is counterproductive for this mode because
each final-output chunk duplicates the same internal alias dependency graph.
The CLI still exposes `--compiled-dag-output-chunk-size` as a debug/memory
experiment, but compiled-DAG presets intentionally do not auto-select chunking.

Current dependency status: the managed installer uses
`valentinHirschi/symbolica_mod` branch `pyamplicol-dev-base`, which is
compatible with the current `spenso`/`idenso` stack and exposes a minimal
Python `aliases=[(handle, body), ...]` evaluator option. The compiled-DAG route
uses this hook to pass current-component alias handles to Symbolica evaluator
construction; metadata should report `symbolica_alias_available = true`. The
fork now also exposes `Expression.alias(handle=None, opaque=True)`, returning
the `(handle, body)` pair accepted by the evaluator alias hook.

Current implementation status:

| Process | Lowering | Backend | Status | Notes |
|---|---|---|---|---|
| `d d~ > z g` | `symbolic` | JIT | validated smoke | Matches the staged native recursion on the canonical point with relative difference below `1e-12`; current-component aliases are passed through Symbolica's Python evaluator alias hook. |
| `d d~ > z g` | `symbolic` | JIT | artifact smoke | Generate/save/load works with the default RAMBO warmup helicity filter metadata persisted in the artifact manifest. |
| `d d~ > z g g` | `symbolic` | JIT | evaluated smoke | The compiled-DAG evaluator runs with the canonical warmup filter and returns the same matrix element convention as the existing native runtime path. |
| `d d~ > z g` | `spenso` | JIT | validated smoke | Current bodies are lowered as bounded Spenso tensor networks using temporary parent-current alias tensors; matches symbolic lowering at the canonical point. |
| `d d~ > z g g` | `spenso` | JIT | validated smoke | Spenso current-body lowering matches symbolic lowering at the canonical point. |
| `d d~ > z g g g` | `spenso` | JIT | validated smoke | Watchdog-backed Spenso-vs-symbolic lowering check on the canonical point gives a maximum relative difference of `1.1e-11` across retained amplitude outputs. |
| `d d~ > z g g g g` | `spenso` | JIT | validated smoke | Watchdog-backed Spenso-vs-symbolic lowering check on the canonical point gives a maximum relative difference of `5.9e-12` across retained amplitude outputs. |
| `d d~ > z g g g g g` | `spenso` | JIT | validated smoke | Watchdog-backed Spenso-vs-symbolic lowering check on the canonical point gives a maximum relative difference of `8.3e-12` across retained amplitude outputs. |
| `d d~ > z g g g g g g` | `spenso` | JIT | validated smoke | Watchdog-backed Spenso-vs-symbolic lowering check on the canonical point gives a maximum relative difference of `9.8e-12` across retained amplitude outputs. |
| `d d~ > z g` | `spenso` | JIT | Fortran comparison smoke | Watchdog-backed `compare-amplicol --amplicol-probe --points 2 --timing 2` gives a maximum relative difference of `8.8e-16` against Fortran AmpliCol. |

The current slice validates the compiled-DAG artifact interface, parameter
layout, real momentum parameter marking, single multi-output amplitude
evaluator, deterministic warmup helicity filter, RAMBO-style on-shell
phase-space sampler, save/load round-trip, and chunked multi-output evaluator
artifact round-trip. It also validates the default Spenso current-body lowering
against the direct symbolic lowering for low multiplicity. The helicity filter
is enabled by default and currently supports two warmup point sources:

- `--compiled-dag-helicity-filter-phase-space rambo` uses deterministic
  RAMBO-style points from `--compiled-dag-helicity-filter-seed`.
- `--compiled-dag-helicity-filter-phase-space canonical` reuses the fixed
  pyAmpliCol sanity points.

Generation progress bars for helicity warmup and alias/current-body
construction are enabled by `--verbose-evaluator-build`. The normal JSON path
stays quiet by default.

Compiled-DAG artifact manifests include Symbolica/spenso/idenso commit
provenance, the parameter layout with real-valued inputs, retained root
metadata, shared-current metadata, evaluator/chunk manifests, and a deterministic
SHA-256 JSON fingerprint.

Full performance rows for `n = 1,...,9` are now present. The `n <= 8` rows
have direct fresh Fortran-probe validation with the source-inlined compiled-DAG
path. The refreshed `n=9` timing row is included, but the fresh Fortran probe
currently fails in the Fortran reference executable with a segmentation fault
before producing probe values.

The optimized low-n path currently uses `--compiled-dag-lowering symbolic`.
The Spenso lowering remains validated and useful as an independent tensor
network cross-check, but direct symbolic current-body lowering gives smaller
Symbolica alias bodies and is faster for the single multi-output evaluator.
The current optimized low-n fast path inlines the physical external
wavefunctions into the alias-backed evaluator. Runtime inputs are therefore
only the real external momenta; source currents and current momenta are built as
Symbolica expressions and optimized together with the current DAG. The old
source-current-parameter path remains available with
`--no-compiled-dag-inline-external-wavefunctions` for debugging. This source
inlining removes the previous Python source-fill bottleneck and brings
wall-time performance to parity for `n <= 5`. A follow-up low-n performance
pass also aliases nonzero external source-wavefunction components and
propagator denominators, marks structurally real aliases as real Symbolica
symbols where applicable, uses `max_common_pair_distance=250` by default for
compiled-DAG evaluator construction, and applies adaptive CPE defaults
(`2` through five gluons, unbounded from six gluons onward). Runtime profiling
now performs a full-batch warm-up before collecting timings so shape-dependent
JIT/setup effects are not included in the measured repetitions. Real-pair
evaluators, current-momentum aliases, selective current aliasing, factor
collection, and per-vertex temporary aliases were benchmarked and disabled by
default because they were slower at `n=5,6`. The single no-chunk evaluator
still hits a code-shape/runtime cliff beyond the parity region.

## Functional Validation Smoke

These rows use the compiled-DAG Spenso lowering, JIT backend, canonical warmup
helicity filter, and a two-point Fortran AmpliCol probe. They are functional
checks, not stable benchmark measurements; two-point wall/evaluator timings are
omitted because they include lazy backend effects and are not representative.

| n in `d d~ > z + n g` | Retained outputs | Alias components | Spenso current bodies | Max relative difference vs Fortran |
|---:|---:|---:|---:|---:|
| 1 | 12 | 44 | 22 | `8.8e-16` |
| 2 | 24 | 124 | 58 | `1.4e-14` |
| 3 | 48 | 348 | 142 | `3.5e-14` |
| 4 | 96 | 836 | 318 | `2.7e-14` |
| 5 | 192 | 1852 | 678 | `7.3e-15` |
| 6 | 384 | 3924 | 1406 | `1.0e-14` |

## n<=9 Profiling Versus Fortran AmpliCol

For pyAmpliCol rows, the multiplier in the `Gen` column is relative to AmpliCol
generation time for the same `n`. Multipliers in `Wall` and `Eval` are relative
to AmpliCol's per-point evaluator timing for the same `n`. Color convention:
<span style="color:#1a7f37">green</span> is faster than AmpliCol,
<span style="color:#bf8700">orange</span> is slower but below `x2.0`, and
<span style="color:#cf222e">red</span> is `x2.0` or worse.

All pyAmpliCol rows below use `--compiled-dag-lowering symbolic`, no final-output
chunking, trusted internally generated canonical benchmark points, and the
source-inlined external wavefunction path described above. The refreshed
`n<=7` rows also include source-component and propagator-denominator aliases.
All compiled-DAG rows build the evaluator with the external momentum parameters
registered as real Symbolica inputs before evaluator optimization/compilation.
The `compiled-complex-4x` backend is not used on this aarch64 machine because
it is guarded as unsafe for the current complex-lane layout. The refresh runs
were executed behind the 30 GiB watchdog; observed peak process-tree RSS was
below `13.39 GiB` for completed probes. Low-n rows now show the three backend
setups matching `performance_summary.md`: JIT, C++ O3, and ASM O3. C++/ASM rows
are no-chunk compiled-DAG rows; the n=5 C++ O3 no-chunk compile reached the
clang++ stage, stayed asleep with no CPU use for about 7.5 minutes, and was
interrupted as an unproductive monolithic compile. An older n=6 C++ O3 no-chunk
probe did complete and is retained for comparison, but larger monolithic C++
rows are marked N/A.

| n | Setup | Gen [s] | Wall [us/pt] | Eval [us/pt] | Notes |
|---:|---|---:|---:|---:|---|
| **1** | **AmpliCol** | **2.22** | **N/A** | **0.734** | **Fortran AmpliCol** |
| 1 | pyAmpliCol - compiled-DAG JIT | 0.0483 <span style="color:#1a7f37">(x0.02)</span> | 0.660 <span style="color:#1a7f37">(x0.90)</span> | 0.108 <span style="color:#1a7f37">(x0.15)</span> | JIT; source and propagator-denominator aliases; batch 256; full-batch warm-up |
| 1 | pyAmpliCol - compiled-DAG C++ | 0.560 <span style="color:#1a7f37">(x0.25)</span> | 0.768 <span style="color:#bf8700">(x1.05)</span> | 0.189 <span style="color:#1a7f37">(x0.26)</span> | C++; O3; chunk=None; batch 256; saved evaluator |
| 1 | pyAmpliCol - compiled-DAG ASM | 0.297 <span style="color:#1a7f37">(x0.13)</span> | 0.978 <span style="color:#bf8700">(x1.33)</span> | 0.403 <span style="color:#1a7f37">(x0.55)</span> | ASM; O3; chunk=None; batch 256; saved evaluator |
| **2** | **AmpliCol** | **2.20** | **N/A** | **1.64** | **Fortran AmpliCol** |
| 2 | pyAmpliCol - compiled-DAG JIT | 0.0941 <span style="color:#1a7f37">(x0.04)</span> | 1.10 <span style="color:#1a7f37">(x0.67)</span> | 0.421 <span style="color:#1a7f37">(x0.26)</span> | JIT; source and propagator-denominator aliases; batch 256; full-batch warm-up |
| 2 | pyAmpliCol - compiled-DAG C++ | 0.726 <span style="color:#1a7f37">(x0.33)</span> | 1.22 <span style="color:#1a7f37">(x0.74)</span> | 0.512 <span style="color:#1a7f37">(x0.31)</span> | C++; O3; chunk=None; batch 256; saved evaluator |
| 2 | pyAmpliCol - compiled-DAG ASM | 0.379 <span style="color:#1a7f37">(x0.17)</span> | 2.24 <span style="color:#bf8700">(x1.36)</span> | 1.51 <span style="color:#1a7f37">(x0.92)</span> | ASM; O3; chunk=None; batch 256; saved evaluator |
| **3** | **AmpliCol** | **2.31** | **N/A** | **4.55** | **Fortran AmpliCol** |
| 3 | pyAmpliCol - compiled-DAG JIT | 0.313 <span style="color:#1a7f37">(x0.14)</span> | 3.18 <span style="color:#1a7f37">(x0.70)</span> | 2.35 <span style="color:#1a7f37">(x0.52)</span> | JIT; source and propagator-denominator aliases; batch 256; full-batch warm-up |
| 3 | pyAmpliCol - compiled-DAG C++ | 2.84 <span style="color:#bf8700">(x1.23)</span> | 2.54 <span style="color:#1a7f37">(x0.56)</span> | 1.70 <span style="color:#1a7f37">(x0.37)</span> | C++; O3; chunk=None; batch 256; saved evaluator |
| 3 | pyAmpliCol - compiled-DAG ASM | 0.715 <span style="color:#1a7f37">(x0.31)</span> | 5.51 <span style="color:#bf8700">(x1.21)</span> | 4.65 <span style="color:#bf8700">(x1.02)</span> | ASM; O3; chunk=None; batch 256; saved evaluator |
| **4** | **AmpliCol** | **2.55** | **N/A** | **13.2** | **Fortran AmpliCol** |
| 4 | pyAmpliCol - compiled-DAG JIT | 0.301 <span style="color:#1a7f37">(x0.12)</span> | 11.8 <span style="color:#1a7f37">(x0.89)</span> | 10.9 <span style="color:#1a7f37">(x0.82)</span> | JIT; source and propagator-denominator aliases; batch 256; full-batch warm-up |
| 4 | pyAmpliCol - compiled-DAG C++ | 40.6 <span style="color:#cf222e">(x15.9)</span> | 16.8 <span style="color:#bf8700">(x1.27)</span> | 15.7 <span style="color:#bf8700">(x1.19)</span> | C++; O3; chunk=None; batch 256; saved evaluator; monolithic C++ slower than JIT |
| 4 | pyAmpliCol - compiled-DAG ASM | 1.02 <span style="color:#1a7f37">(x0.40)</span> | 18.5 <span style="color:#bf8700">(x1.40)</span> | 17.4 <span style="color:#bf8700">(x1.32)</span> | ASM; O3; chunk=None; batch 256; saved evaluator |
| **5** | **AmpliCol** | **3.22** | **N/A** | **34.3** | **Fortran AmpliCol** |
| 5 | pyAmpliCol - compiled-DAG JIT | 0.435 <span style="color:#1a7f37">(x0.14)</span> | 32.8 <span style="color:#1a7f37">(x0.96)</span> | 31.0 <span style="color:#1a7f37">(x0.90)</span> | JIT; source and propagator-denominator aliases; batch 32; full-batch warm-up |
| 5 | pyAmpliCol - compiled-DAG C++ | N/A | N/A | N/A | C++; O3; chunk=None; monolithic compile did not complete productively |
| 5 | pyAmpliCol - compiled-DAG ASM | 1.94 <span style="color:#1a7f37">(x0.60)</span> | 49.5 <span style="color:#bf8700">(x1.44)</span> | 47.6 <span style="color:#bf8700">(x1.39)</span> | ASM; O3; chunk=None; batch 32; saved evaluator |
| **6** | **AmpliCol** | **4.55** | **N/A** | **84.1** | **Fortran AmpliCol** |
| 6 | pyAmpliCol - compiled-DAG JIT | 1.26 <span style="color:#1a7f37">(x0.28)</span> | 89.7 <span style="color:#bf8700">(x1.07)</span> | 87.7 <span style="color:#bf8700">(x1.04)</span> | JIT; adaptive CPE unbounded; source and propagator-denominator aliases; batch 32; full-batch warm-up; max rel diff `5.1e-14` over 5 probe points |
| 6 | pyAmpliCol - compiled-DAG C++ | 814 <span style="color:#cf222e">(x179)</span> | 108 <span style="color:#bf8700">(x1.28)</span> | 106 <span style="color:#bf8700">(x1.26)</span> | C++; O3; chunk=None; older completed probe; generation is not practical |
| 6 | pyAmpliCol - compiled-DAG ASM | 3.99 <span style="color:#1a7f37">(x0.88)</span> | 121 <span style="color:#bf8700">(x1.44)</span> | 120 <span style="color:#bf8700">(x1.42)</span> | ASM; O3; chunk=None; batch 1024 |
| **7** | **AmpliCol** | **9.38** | **N/A** | **207** | **Fortran AmpliCol** |
| 7 | pyAmpliCol - compiled-DAG JIT | 3.90 <span style="color:#1a7f37">(x0.42)</span> | 268 <span style="color:#bf8700">(x1.29)</span> | 265 <span style="color:#bf8700">(x1.28)</span> | JIT; adaptive CPE unbounded; source and propagator-denominator aliases; batch 24; full-batch warm-up; max rel diff `3.7e-13` over 3 probe points |
| 7 | pyAmpliCol - compiled-DAG C++ | N/A | N/A | N/A | C++; O3; chunk=None; not attempted after monolithic n=5 C++ compile failed productively |
| 7 | pyAmpliCol - compiled-DAG ASM | 11.1 <span style="color:#bf8700">(x1.19)</span> | 369 <span style="color:#bf8700">(x1.78)</span> | 364 <span style="color:#bf8700">(x1.76)</span> | ASM; O3; chunk=None; batch 512 |
| **8** | **AmpliCol** | **25.0** | **N/A** | **567** | **Fortran AmpliCol** |
| 8 | pyAmpliCol - compiled-DAG JIT | 16.6 <span style="color:#1a7f37">(x0.67)</span> | 3601 <span style="color:#cf222e">(x6.35)</span> | 3596 <span style="color:#cf222e">(x6.34)</span> | JIT; source and propagator-denominator aliases; batch 512 |
| 8 | pyAmpliCol - compiled-DAG C++ | N/A | N/A | N/A | C++; O3; chunk=None; not attempted after monolithic n=5 C++ compile failed productively |
| 8 | pyAmpliCol - compiled-DAG ASM | 34.7 <span style="color:#bf8700">(x1.39)</span> | 1432 <span style="color:#cf222e">(x2.53)</span> | 1430 <span style="color:#cf222e">(x2.52)</span> | ASM O3; inline external wavefunctions; max rel diff `1.4e-13` over 3 probe points |
| **9** | **AmpliCol** | **48.0** | **N/A** | **2520** | **Fortran AmpliCol** |
| 9 | pyAmpliCol - compiled-DAG JIT | 122 <span style="color:#cf222e">(x2.54)</span> | 8056 <span style="color:#cf222e">(x3.20)</span> | 8051 <span style="color:#cf222e">(x3.19)</span> | JIT; inline external wavefunctions; fresh Fortran probe currently fails with reference-side SIGSEGV |
| 9 | pyAmpliCol - compiled-DAG C++ | N/A | N/A | N/A | C++; O3; chunk=None; not attempted after monolithic n=5 C++ compile failed productively |
| 9 | pyAmpliCol - compiled-DAG ASM | 191 <span style="color:#cf222e">(x3.97)</span> | 10705 <span style="color:#cf222e">(x4.25)</span> | 10700 <span style="color:#cf222e">(x4.25)</span> | ASM; O3; chunk=None; batch 512; peak RSS `8.01 GiB` |

Current conclusion: after inlining physical external wavefunctions and aliasing
source components plus propagator denominators, the compiled-DAG route is at
wall-time and evaluator-only parity for `n <= 5`, and close to parity for
`n=6` (`x1.04` evaluator-only). The fresh `n=7` row improves over the older
profile but remains slower than Fortran (`x1.28` evaluator-only). The refreshed
n=8 JIT and n=9 ASM rows confirm that higher multiplicities are dominated by
evaluator code shape rather than Python overhead. The remaining Python overhead
in the refreshed rows is mainly external-momentum object-to-array packing and
output transfer; source fill and current-momentum setup have moved into the
evaluator.

The current compiled-DAG route is mathematically equivalent to the eager DAG,
but it is not yet the compiled-code equivalent of the eager DAG's staged current
buffer schedule. The eager route compiles bounded current/interaction/amplitude
blocks and materializes current arrays between blocks. The compiled-DAG route
passes one large alias-backed multi-output expression to Symbolica. Generated
C++ for the n=4 no-chunk row lowers this to generic `Z[...]` temporaries
(`605` buffer entries, `11,625` source lines) rather than an AmpliCol-style
current-table schedule. This explains why compiled-DAG C++ is excellent at
`n=1,2,3` but loses to the old eager C++ evaluator from `n=4` onward. The next
performance target is therefore a single Python call whose generated code still
preserves staged current storage internally, instead of a fully flattened
multi-output kernel.

JIT compiled-DAG evaluators support the same two-step save/load workflow as the
generated-code backends. The artifact stores the serialized Symbolica evaluator
state (`symbolica-evaluator-state`) and reloads it with `Evaluator.load(...)`.
The generation-time number reported by a loaded run is therefore not zero
because pyAmpliCol still rebuilds process metadata and restores the helicity
filter, but the Symbolica evaluator itself is loaded instead of rebuilt.

Forcing `--symbolica-no-jit-direct-translation` was tested and is strongly
disfavored for the current compiled-DAG expression shape:

| n | JIT mode | Gen [s] | Wall [us/pt] | Eval [us/pt] | Notes |
|---:|---|---:|---:|---:|---|
| 6 | direct translation | 1.29 | 95.4 | 93.1 | batch 32 |
| 6 | no direct translation | 1.29 | 356 | 352 | batch 32 |
| 7 | direct translation | 4.01 | 449 | 445 | batch 24 |
| 7 | no direct translation | 4.01 | 10421 | 10416 | batch 24 |
| 8 | direct translation | 16.6 | 3601 | 3596 | batch 512 |
| 8 | no direct translation | 16.7 | 48836 | 48831 | batch 512 |

These probes were run behind the 30 GiB watchdog; peak process-tree RSS was
`10.17 GiB`. The higher-multiplicity default should therefore keep JIT direct
translation enabled.

### n=6 and n=7 Backend Probes

These probes use the same symbolic compiled-DAG expression and the 30 GiB
watchdog, but are not all selected main table rows. They test whether generated
code can recover runtime performance at the first multiplicities above parity.

| n | Backend probe | Gen [s] | Wall [us/pt] | Eval [us/pt] | Peak RSS | Outcome |
|---:|---|---:|---:|---:|---:|---|
| 6 | JIT adaptive CPE | 1.26 | 89.7 | 87.7 | 0.34 GiB | Selected n=6 row; source and propagator-denominator aliases; full-batch warm-up. |
| 6 | C++ O3 no chunk | 814 | 108 | 106 | 13.39 GiB | Completes, but generation is expensive and runtime is slower than JIT. |
| 6 | ASM O3 no chunk | 3.99 | 121 | 120 | 0.23 GiB | Completes quickly, but is slower than JIT. |
| 7 | JIT adaptive CPE | 3.90 | 268 | 265 | 1.18 GiB | Selected n=7 row; batch 24; full-batch warm-up. |
| 7 | ASM O3 no chunk | 11.1 | 369 | 364 | 1.07 GiB | Completes quickly, but is slower than JIT. |

### n<=7 Fast-Path Optimization Checks

The following variants were tested while hammering down the compiled-DAG
runtime before proceeding to higher multiplicity. They are disabled by default
unless explicitly re-enabled in code for debugging:

| Variant | n | Gen [s] | Wall [us/pt] | Eval [us/pt] | Decision |
|---|---:|---:|---:|---:|---|
| Source aliases only | 5 | 0.458 | 32.4 | 31.2 | Kept as part of fast path. |
| Source + propagator-denominator aliases | 5 | 0.435 | 32.8 | 31.0 | Kept; current best n=5 JIT path with batch 32. |
| Source + propagator-denominator aliases | 6 | 1.26 | 89.7 | 87.7 | Kept; current best n=6 JIT path with adaptive CPE and batch 32. |
| Real-pair evaluator | 5 | 0.647 | 39.9 | 38.6 | Disabled; real evaluator is slower than Symbolica complex JIT here. |
| Real-pair evaluator | 6 | 2.09 | 130 | 128 | Disabled. |
| Current-momentum aliases | 5 | 0.445 | 32.6 | 31.3 | Disabled; no useful n=5 gain and n=6 regresses/gets unstable. |
| Current-momentum aliases | 6 | 1.20 | 98.0 | 95.9 | Disabled; runtime regresses. |
| Selective current aliases | 5 | 0.521 | 35.9 | 34.1 | Disabled; inlining single-use currents is slower. |
| Selective current aliases | 6 | 1.39 | 89.8 | 87.8 | Disabled; no useful runtime gain. |
| Factor collection | 6 | 1.26 | 96.7 | 94.6 | Disabled; expression shape is slower. |
| Vertex temporary aliases | 5 | 0.520 | 36.2 | 34.4 | Disabled; too many aliases for no gain. |
| Vertex temporary aliases | 6 | 1.34 | 102 | 99.7 | Disabled. |
| n=7 batch 24, common-pair distance 250 | 7 | 3.90 | 268 | 265 | Selected n=7 row. |
| n=7 common-pair distance 100 | 7 | 3.92 | 268 | 265 | Essentially tied with 250; keep compiled-DAG default 250. |
| n=7 common-pair distance 500 | 7 | 3.90 | 272 | 270 | Disabled; runtime regresses. |
| n=7 CPE iterations 1 | 7 | 3.58 | 276 | 274 | Disabled; lower generation time but slower runtime. |
| n=7 CPE iterations 2 | 7 | 3.62 | 289 | 286 | Disabled; runtime regresses. |
| n=7 JIT optimization level 1 | 7 | 3.89 | 269 | 267 | Disabled; slightly slower than optimization level 3. |

### n=8 Backend Probe

The n=8 JIT evaluator remains feasible but hits a larger runtime cliff. JIT
iterations did not help in the earlier sweep, and disabling JIT direct
translation is much worse. Inline ASM is currently the best measured generated
code route.

| Backend probe | Gen [s] | Wall [us/pt] | Eval [us/pt] | Peak RSS | Outcome |
|---|---:|---:|---:|---:|---|
| JIT direct translation | 16.6 | 3601 | 3596 | 10.17 GiB sweep | Feasible, but `x6.3` slower than Fortran evaluator time. |
| JIT no direct translation | 16.7 | 48836 | 48831 | 10.17 GiB sweep | Strongly disfavored. |
| JIT iterations=1, older sweep | 16.6 | 2785 | 2782 | 2.75 GiB sweep | Older pre-final-tightening timing retained only as context. |
| JIT iterations=2 | 16.7 | 2782 | 2779 | 2.75 GiB sweep | No material improvement. |
| ASM O3 no chunk | 34.7 | 1432 | 1430 | 2.75 GiB | Best measured n=8 route so far, but still `x2.5` slower than Fortran. |

### n=9 Reference Status

The refreshed n=9 JIT timing row uses inline external wavefunctions and improves
over the pre-inlining JIT row, but a fresh Fortran comparison could not be
completed in this refresh: a two-point `--amplicol_probe` run segfaulted in the
Fortran `amplicol_generate` executable, and a one-point retry was interrupted
after several minutes without producing reference values. The older
pre-inlining n=9 row had shown agreement at the `3.4e-14` level, but the fresh
source-inlined n=9 row should be treated as timing-only until the Fortran probe
is reliable again at this multiplicity.

## Process-Generic Status

The compiled-DAG route is architecturally separate from the old runtime path,
but the current lowering is still limited to the existing `q q~ -> Z + n g`
recursion graph. The following planned process-generic checks currently fail
early with clear unsupported diagnostics:

| Process | Current status |
|---|---|
| `d d~ > e+ e- g g` | Unsupported in compiled-DAG artifacts: pending the lepton-source shared-current artifact route. |
| `u d~ > w+ g g` | Supported by the eager-DAG/Rusticol vector-plus-gluon route; compiled-DAG support remains outside the current exploratory lowering. |
| `d d~ > z z g` | Unsupported in compiled-DAG artifacts: pending the model-driven multi-boson current lowering. |

## Reproduction Commands

The selected pyAmpliCol rows in the main table were produced with
`profile-dag-evaluator`, `--compiled-dag-lowering symbolic`,
`--compiled-dag-helicity-filter-phase-space canonical`, no final-output chunking,
and the following row-specific settings:

| n | Process | Backend flags | Points/repetitions/batch |
|---:|---|---|---|
| 1 | `'d d~ > z g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 8192 --repetitions 5 --batch-size 256` |
| 2 | `'d d~ > z g g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 8192 --repetitions 5 --batch-size 256` |
| 3 | `'d d~ > z g g g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 8192 --repetitions 5 --batch-size 256` |
| 4 | `'d d~ > z g g g g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 8192 --repetitions 5 --batch-size 256` |
| 5 | `'d d~ > z g g g g g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 8192 --repetitions 5 --batch-size 32` |
| 6 | `'d d~ > z g g g g g g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 8192 --repetitions 5 --batch-size 32` |
| 7 | `'d d~ > z g g g g g g g'` | `--symbolica-evaluator-backend jit --symbolica-iterations 1` | `--points 4096 --repetitions 5 --batch-size 24` |
| 8 | `'d d~ > z g g g g g g g g'` | `--symbolica-evaluator-backend compiled-complex --symbolica-compiled-preset manual --symbolica-compiled-inline-asm default --symbolica-compiled-optimization-level 3` | `--points 1024 --repetitions 2 --batch-size 1024` |
| 9 | `'d d~ > z g g g g g g g g g'` | `--symbolica-evaluator-backend jit` | `--points 512 --repetitions 1 --batch-size 512` |

For example, the n=8 selected row is reproduced with:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
source "$HOME/.bashrc"
PYTHONPATH=pyAmpliCol/src pyAmpliCol/scripts/run_with_memory_watch.py \
  --limit-gb 30 \
  --report-json pyAmpliCol/outputs/compiled-dag/profile-n8-asm-o3-final-watchdog.json \
  -- \
  pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol profile-dag-evaluator \
    --compiled-dag-evaluator \
    --compiled-dag-lowering symbolic \
    --compiled-dag-helicity-filter-phase-space canonical \
    --points 1024 \
    --repetitions 2 \
    --batch-size 1024 \
    --symbolica-evaluator-backend compiled-complex \
    --symbolica-compiled-preset manual \
    --symbolica-compiled-inline-asm default \
    --symbolica-compiled-optimization-level 3 \
    --json \
    'd d~ > z g g g g g g g g'
```

Validation rows through `n=8` use the same process strings with
`compare-amplicol --compiled-dag-evaluator --compiled-dag-lowering symbolic
--compiled-dag-helicity-filter-phase-space canonical --amplicol-probe`.
The refreshed `n=9` row is timing-only because the fresh Fortran probe currently
segfaults in `amplicol_generate`.

Generate and save a JIT compiled-DAG evaluator:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol generate \
  --compiled-dag-evaluator \
  --compiled-dag-lowering symbolic \
  --compiled-dag-helicity-filter-phase-space rambo \
  --symbolica-evaluator-backend jit \
  --save-evaluator-dir pyAmpliCol/outputs/compiled-dag/n1-jit \
  'd d~ > z g' \
  --json
```

Evaluate from the saved artifact:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol evaluate \
  --runtime-backend compiled-dag \
  --compiled-dag-lowering symbolic \
  --batch-size 1024 \
  --load-evaluator-dir pyAmpliCol/outputs/compiled-dag/n1-jit \
  'd d~ > z g' \
  --json
```

Generate a debug-only chunked JIT artifact for output chunking experiments:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol generate \
  --compiled-dag-evaluator \
  --compiled-dag-lowering symbolic \
  --no-compiled-dag-helicity-filter \
  --compiled-dag-output-chunk-size 5 \
  --symbolica-evaluator-backend jit \
  --save-evaluator-dir pyAmpliCol/outputs/compiled-dag/n1-jit-chunked \
  'd d~ > z g' \
  --json
```

Run a low-multiplicity Spenso-vs-symbolic lowering cross-check:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol evaluate \
  --compiled-dag-evaluator \
  --compiled-dag-cross-check-lowering \
  --compiled-dag-helicity-filter-phase-space canonical \
  --symbolica-evaluator-backend jit \
  'd d~ > z g' \
  --json
```

Run a small Fortran AmpliCol comparison probe:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
source "$HOME/.bashrc"
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python \
  pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30 -- \
  pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol compare-amplicol \
    --compiled-dag-evaluator \
    --compiled-dag-lowering symbolic \
    --compiled-dag-helicity-filter-phase-space canonical \
    --amplicol-probe \
    --points 2 \
    --timing 2 \
    --jobs 8 \
    'd d~ > z g' \
    --json
```

Run a compiled-DAG profiler smoke:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
source "$HOME/.bashrc"
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol \
  profile-dag-evaluator \
    --compiled-dag-evaluator \
    --compiled-dag-lowering symbolic \
    --compiled-dag-helicity-filter-phase-space canonical \
    --points 4 \
    --repetitions 1 \
    --batch-size 4096 \
    --json \
    'd d~ > z g'
```

Run the focused unit tests:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
source "$HOME/.bashrc"
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python -m pytest \
  pyAmpliCol/tests/unit/test_compiled_dag_runtime.py \
  pyAmpliCol/tests/unit/test_cli.py::test_cli_compiled_dag_shortcut_selects_runtime_backend \
  pyAmpliCol/tests/unit/test_cli.py::test_cli_compiled_dag_helicity_filter_flags \
  pyAmpliCol/tests/unit/test_cli.py::test_cli_profile_dag_evaluator_supports_compiled_dag_backend \
  -q
```

Heavy validation and benchmark runs for this mode must be launched behind the
30 GiB watchdog. Use `--report-json` to persist the observed peak process-tree
RSS next to the saved evaluator or profile output:

```bash
cd /Users/vjhirsch/HEP_programs/AmpliCol
source "$HOME/.bashrc"
PYTHONPATH=pyAmpliCol/src pyAmpliCol/dependencies/.venv/bin/python \
  pyAmpliCol/scripts/run_with_memory_watch.py \
    --limit-gb 30 \
    --report-json pyAmpliCol/outputs/compiled-dag/watchdog-n3.json \
    -- \
  pyAmpliCol/dependencies/.venv/bin/python -m pyamplicol evaluate \
    --compiled-dag-evaluator \
    --compiled-dag-lowering symbolic \
    --batch-size 4096 \
    'd d~ > z g g g' \
    --json
```
