# pyAmpliCol

pyAmpliCol is the Python generation and validation layer for the modern
AmpliCol matrix-element prototype.  The current production path generates a
shared-current eager-DAG process artifact in Python and evaluates it through the
Rusticol PyO3 runtime using serialized Symbolica evaluators.

The production path is now the leading-colour generic DAG artifact route.  The
process parser/enumerator expands concrete subprocesses, process sets,
anonymous multiparticle slots, built-in `p`/`j` labels, repetition syntax, and
arbitrary quark-line counts.  The model supplies particles, vertices,
propagators, source kernels, and vertex lowerings; the generic DAG compiler
discovers valid currents from that local model information rather than from
process-family branches.  Generated schema-v2 artifacts contain serialized
Symbolica evaluator stages.  `pyamplicol.load_process(..., runtime="python")`
can inspect the same generic artifact metadata and stage plan, while
`runtime="rusticol"` is the production evaluator for serialized stage
execution.

Current validation covers representative `V + gluons`, explicit dilepton,
charged-current, multi-boson, pure-gluon, and multi-quark-line leading-colour
processes.  Unsupported requests fail explicitly with diagnostics that name the
missing layer, such as colour expansion, current closure, or a missing
vertex/propagator lowering.

The external-model validation fixture now contains 180 populated
LC/NLC/full-colour Standard Model cases through final-state multiplicity nine
from the bundled loader JSON.  It covers all three colour modes through `n=5`
and every populated LC reference through `n=9`.  At identical stored
phase-space points and matched parameters, all 180 agree with the
built-in/AmpliCol references at `1e-10` tolerance; the fixture's largest
relative difference is `9.91e-13`.  Incremental portable reports for `n=4`
through `n=9` live under `docs/ufo_sm_n*_validation_report.json`; every
generated process artifact is retained below `.ufo_support_outputs/`.

## Installation

From this directory, install the managed Python and native dependencies:

```sh
python3 ./dependencies/install_dependencies.py
```

This creates `dependencies/.venv`, installs pyAmpliCol in editable mode, builds
the Rusticol PyO3 extension, and installs the managed Symbolica/spenso/idenso
environment.  The GammaLoop Python API is optional and is not installed by
default; request it explicitly with:

```sh
python3 ./dependencies/install_dependencies.py --with-gammaloop
```

## Symbolica License

Symbolica requires a license for evaluator generation.  If you do not already
have one, request a free trial from the managed Python environment:

```sh
source dependencies/.venv/bin/activate
python3
```

```python
from symbolica import *
request_trial_license('NAME', 'EMAIL', 'ORGANIZATION')
```

Save the returned license string in the `SYMBOLICA_LICENSE` environment
variable before running pyAmpliCol:

```sh
export SYMBOLICA_LICENSE='PASTE_THE_RETURNED_LICENSE_HERE'
```

## Quick Start

Generate a reusable process artifact and then time it with the default Rusticol
runtime:

```sh
./pyamplicol.sh generate-process 'd d~ > Z g g g g' outputs/dd_z_4g
./pyamplicol.sh time-process outputs/dd_z_4g
```

## Models And Run Cards

The default model is the hard-coded Standard Model, named `built-in-sm`.
`BUILTIN_SM` and `builtin_sm` remain accepted aliases.  The same commands also
accept a trusted UFO directory, a JSON model serialized by
[`ufo-model-loader`](https://github.com/alphal00p/ufo_model_loader), or an
exact-version `*.pyAmplicol-model.json` artifact:

```sh
./pyamplicol.sh inspect-model --model assets/models/json/sm/sm.json
./pyamplicol.sh compile-model --model assets/models/ufo/sm \
  --output outputs/models/sm.pyAmplicol-model.json
./pyamplicol.sh generate-process \
  --model outputs/models/sm.pyAmplicol-model.json \
  'd d~ > z g' outputs/external_sm_dd_z_g
```

UFO directories are executable Python and must be trusted.  pyAmpliCol imports
them in an isolated worker, applies the requested UFO restriction and
simplification, lowers their tensors into its native model IR, and caches the
result by source contents and dependency fingerprint.  Loader JSON starts at
the same post-import representation.  Compiled pyAmpliCol models are portable
data rather than machine code, but deliberately require the exact matching
compiler/dependency fingerprint; regenerate one after upgrading pyAmpliCol or
its managed dependencies.  The implementation follows the
[UFO 2.0 specification](https://arxiv.org/abs/2304.09883).  The dependency
installer pins the tested `ufo-model-loader` 0.1.7 revision rather than relying
on the GammaLoop reference checkout.

Generic kernel components use named Symbolica functions during model lowering.
Before optimization, process generation inlines component bodies up to 1 KiB,
allowing common-subexpression elimination across current boundaries.  Larger
spin-2 and high-rank bodies remain in Symbolica's stable `functions` map so
stage construction does not replicate them throughout the expression tree.

Every public subcommand accepts `--run-card PATH`.  A card contains one
`[arguments]` table whose keys may be argparse destination names or normalized
public long-option names (`--symbolica-output-chunk-size` becomes
`symbolica_output_chunk_size`), including positionals.  Relative paths are
resolved from the card directory.  Explicit CLI options override card values:

```toml
[arguments]
model = "../assets/models/json/sm/sm.json"
process = "d d~ > z g"
output_dir = "../outputs/external_sm_dd_z_g"
color_accuracy = "lc"
symbolica_iterations = 10
```

```sh
./pyamplicol.sh generate-process --run-card examples/ufo_sm_generate.toml \
  --symbolica-iterations 20
```

Unknown keys, arguments belonging to another subcommand, and missing required
inputs are errors.  The cards under `examples/` cover the built-in, UFO, loader
JSON, compiled-model, runtime-parameter, multiparticle, colour-accuracy,
massless-scalar, and scalar-gravity workflows. The contact-only, uncapped
all-scalar-tree, and scalar-gravity validation ladders and their generated TeX
tables are maintained by `docs/model_results.py`;
all process outputs are retained under `.ufo_support_outputs/model_matrices`.

For normal selected-flow JIT artifacts, the CLI uses batch 128 and measures
stage output chunks around the base size 128.  A candidate replaces uniform
chunk 128 only when its stage microbenchmark is at least 5% faster. On AArch64,
the score counts one shared native input pack and sums evaluator/output-unpack
work across chunks, matching Rusticol's fused stage call. Use
`--symbolica-output-chunk-strategy uniform` for fixed-backend comparisons;
the all-flow documentation benchmarks pair that with output chunk 8192 and
runtime batch 64 to control memory use.

On AArch64, Rusticol uses native two-lane complex SymJIT payloads and packs a
stage's local inputs once for all of its output chunks.  Set
`RUSTICOL_NATIVE_SIMD_JIT=0` only when a scalar-JIT diagnostic comparison is
needed.

The managed dependency install pins SymJIT `2.19.3` at commit
`7fb09d1cb2a943c25a6fd71a208af44fcc6d813d` and applies four AArch64 fixes.
Patch `0001` corrects O3 complex register allocation, including aliased
operands and ABI-preserved SIMD registers. Patch `0002` materializes external
call spill addresses beyond the 12-bit immediate range. Patch `0003` splits
large stack-frame adjustments into encodable instructions. Patch `0004`
expands out-of-range conditional branches into a short inverted branch plus an
unconditional long branch. Their regression tests cover the wrong-result,
invalid spill address, `x < 4096`, and `x.abs() < 1048576` failure modes,
respectively.

Large runnable schema-v2 artifacts compact generation-only current/value
metadata and retain stage interaction IDs plus vertex-kind counts instead of a
second copy of every interaction. Evaluators are lowered, compiled, and written
one stage at a time so prior Symbolica expressions can be released. Small
detailed and plan-only artifacts keep the full records; large diagnostic
sidecars explicitly identify their compact layout. Rusticol accepts both forms,
so existing artifacts remain usable.

The generated process directory is self-contained.  It includes a
`process_manifest.json`, serialized evaluator artifacts under `evaluators/`,
validation momenta, and a standalone checker:

```sh
cd outputs/dd_z_4g
python3 check_standalone.py --precision 16 --profile
```

Masses, widths, normalization inputs, and model couplings remain runtime
parameters.  Their exact names and defaults are listed under
`runtime_schema.model_parameters` in the process manifest.  Override any used
subset with a loader-compatible complex JSON file:

```json
{
  "normalization.alpha_s_me_check": [0.118, 0.0],
  "particle.23.mass": [91.188, 0.0],
  "particle.23.width": [2.441404, 0.0],
  "coupling.10.1_23_1.component_0": [-1.0244420275940371, 0.0]
}
```

```sh
./pyamplicol.sh time-process outputs/dd_z_4g \
  --model-parameters model-parameters.json
```

Rusticol rejects unknown parameter names instead of silently baking or
ignoring them.

## Useful Commands

```sh
./pyamplicol.sh processes 'd d~ > Z g g' --json
./pyamplicol.sh processes 'd d~ > e+ e- [d g] [d g] | p p > Z 2j' --json
./pyamplicol.sh process-plan 'd d~ > Z Z g' outputs/plans/dd_zz_g
./pyamplicol.sh generate-process 'd d~ > e+ e-' outputs/dd_epem_0g
./pyamplicol.sh time-process outputs/dd_epem_0g
./pyamplicol.sh generate-process 'd d~ > e+ e- g' outputs/dd_epem_1g
./pyamplicol.sh time-process outputs/dd_epem_1g
./pyamplicol.sh generate-process 'u d~ > e+ ve g' outputs/ud_epve_1g
./pyamplicol.sh time-process outputs/ud_epve_1g
./pyamplicol.sh generate-process 'd d~ > u u~ s s~' outputs/dd_uuss
./pyamplicol.sh time-process outputs/dd_uuss
./pyamplicol.sh compare-amplicol 'd d~ > Z g g g g' --runtime-backend rusticol --points 10
```

`generate-process` also accepts process sets.  A process-set output contains a
root `process_set_manifest.json` and one nested subprocess artifact per
canonical process key.  The root also has a `check_standalone.py` wrapper that
forwards to the selected subprocess checker:

```sh
./pyamplicol.sh generate-process 'd d~ > Z g | u u~ > Z g' outputs/z_1g_set
./pyamplicol.sh time-process outputs/z_1g_set --process u_ubar_to_z_g
python outputs/z_1g_set/check_standalone.py --process 'u u~ > z g' --precision 16
```

Inclusive labels are expanded at the process-set boundary.  For example,
`p p > Z g` becomes concrete subprocess entries such as `d d~ > g z`,
`u u~ > g z`, and their reversed beam-order partners.  The vector-plus-gluon
artifact route supports both beam orders for these concrete entries.  Entries
whose current lowering is not implemented yet fail with explicit diagnostics
instead of being silently folded into a parent inclusive artifact.

The production artifact defaults to leading colour:

```sh
./pyamplicol.sh generate-process --color-accuracy lc 'd d~ > Z 4*g' outputs/dd_z_4g
```

An LC artifact can retain every ordering while also compiling pruned runtime
sidecars for selected flows.  Omitting `--lc-sector-ids` at timing uses the
main artifact reduction; supplying it loads the matching sidecar when present:

```sh
./pyamplicol.sh generate-process --lc-sector-strategy all \
  --runtime-lc-sector-ids 0 'd d~ > Z 4*g' outputs/dd_z_4g_all
./pyamplicol.sh time-process --lc-sector-ids 0 outputs/dd_z_4g_all
./pyamplicol.sh time-process outputs/dd_z_4g_all
```

NLC and full-colour matrix elements are available through the same generic
staged-DAG/Rusticol route.  They reuse the same colour-ordered amplitude vector
and change only the final sparse colour contraction.  Pure-gluon, one-quark-line
and two-quark-line contractions are matched directly against Fortran AmpliCol;
multi-quark-line contractions use the same open-line trace-overlap rule
generically, including gluon strings on the open lines:

```sh
./pyamplicol.sh generate-process --color-accuracy nlc 'g g > g g g' outputs/gg_3g_nlc
./pyamplicol.sh generate-process --color-accuracy full 'd d~ > t t~ g' outputs/dd_tt_g_full
./pyamplicol.sh generate-process --color-accuracy nlc 'd d~ > u u~ s s~ g' outputs/dd_3q_g_nlc
```

Large inclusive or many-quark-line requests can be controlled with generic
physics budgets, not process-family modes.  Useful examples are coupling-order
caps, explicit LC sector ids, line-pairing representative sectors, and a cap on
external quark-pair count:

Generic current and colour-sector caps are unbounded by default.  Prefer the
RAM watchdog for production benchmark refreshes; pass `--max-currents` or
`--max-color-sectors` only when you intentionally want to truncate exploration.

```sh
./pyamplicol.sh process-plan --coupling-order-policy minimal 'd d~ > Z 4*g' outputs/plans/dd_z_4g
./pyamplicol.sh process-plan --lc-sector-strategy line-pairing-representatives 'd d~ > u u~ s s~ c c~' outputs/plans/dd_3pairs
./pyamplicol.sh generate-process --max-quark-pairs 2 'p p > Z 2j' outputs/pp_z_2j_qcap2
```

For more than two quark-pair lines, Fortran AmpliCol's built-in colour matrix
is not the validation reference.  pyAmpliCol still generates the NLC/full sparse
metric; dedicated raw-amplitude probes are used where a direct Fortran
comparison is needed.

Long validation and benchmark jobs should be run behind the bundled RAM
watchdog so the process is stopped before exceeding the repository guideline of
30 GB.  On macOS the watchdog checks compressed physical footprint in addition
to RSS, so memory moved into the compressor cannot evade the cap:

```sh
python3 ./scripts/run_with_memory_watch.py --limit-gb 30 -- \
  ./pyamplicol.sh compare-amplicol 'd d~ > Z g g g g' \
    --runtime-backend rusticol --points 10
```

## Documentation

See `docs/pyAmpliCol.pdf` for the architecture overview and the live generated
LC/NLC/full-colour and dedicated `d d~ > Z + (n-1)g` performance tables.
The generated JSON/TeX caches next to the PDF are the source of truth for the
current refresh.  Matrix tables are refreshed with `docs/result_matrix.py`; the
dedicated Z-family table is populated with `docs/run_z_performance_table.py`,
which preserves generated processes in `docs/.z_performance_outputs`.
Use `--force-pyamplicol-regeneration` for a full current-code Z refresh: the
driver retains the LC-cache AmpliCol reference but rebuilds every requested
pyAmpliCol selected/all-flow artifact before timing it. Use
`--retime-existing` after runtime-only changes to retain generation times and
rerun both timing blocks from the saved artifacts.
Matrix process artifacts are likewise retained in
`docs/.result_matrix_outputs`, so completed cells can be retimed or profiled
without regenerating their evaluators.
`docs/performance_summary.md` and the older comparison notes are retained as
historical benchmark snapshots, while `docs/process_coverage.md` summarizes the
current parser/model/runtime coverage.
