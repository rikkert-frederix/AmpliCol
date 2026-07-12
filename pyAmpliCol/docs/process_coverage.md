# pyAmpliCol Process Coverage

This document separates three layers that are easy to conflate:

- **Process syntax/enumeration**: parsing a user request into partonic
  subprocesses, colour orders, phase-space groups, multichannel partners, and
  symmetry factors.
- **Canonical process IR**: the typed, all-outgoing description shared by
  process-support diagnostics and the planned model-driven DAG builder.  It
  keeps both physical incoming/final legs and crossed all-outgoing particles,
  records PDG orders, labels gluons/leptons/vectors/Higgs legs, and summarizes
  arbitrary quark-line balance.
- **Matrix-element lowering**: constructing the AmpliCol-style shared-current
  DAG and Symbolica evaluator inputs for a concrete process.
- **Runtime execution**: loading a generated artifact with Python or Rusticol
  and evaluating phase-space points.  Python can load generic schema-v2
  metadata for diagnostics and stage-plan inspection; Rusticol remains the
  production evaluator for serialized schema-v2 stage execution.

## Current Status

| Feature | Status | Notes |
|---|---|---|
| Legacy particle vocabulary | Implemented in process/model metadata | QCD partons, `a`, `z`, `w+`, `w-`, charged leptons, neutrinos, and `h`. |
| Process sets with `PROC | PROC` | Implemented in parser/CLI artifacts | Multi-entry generation writes a root `process_set_manifest.json`, a root standalone checker wrapper, and nested subprocess artifacts. |
| Crossing-equivalent subprocess reuse | Implemented for process-set artifacts | Subprocesses with the same all-outgoing PDG multiset share the first generated representative artifact; the process-set manifest records `crossing_alias_of` and an input momentum crossing/permutation map for Rusticol. |
| Built-in `p` and `j` labels | Implemented | Both use the active massless-QCD flavour scheme, matching legacy AmpliCol. |
| Anonymous labels such as `[d g]` | Implemented in parser/enumerator | Anonymous slots expand by cartesian product; invalid charge/flavour combinations are skipped. |
| Repetition syntax such as `4*g` and `3*[d g]` | Implemented | Repeated anonymous slots are independent. |
| Canonical process IR | Implemented for parser/support layers | `process_ir.py` is the stable process description intended to feed generic current construction and Rusticol schema v2. |
| Generic colour-flow planning | Implemented for production LC and initial NLC/full artifacts | `color_plan.py` enumerates leading-colour open-line sectors for arbitrary balanced quark-pair counts and pure-gluon single-trace sectors.  For NLC/full, it keeps the additional trace/open-line orderings needed by the sparse colour-contraction plan for Fortran-supported 0/1/2 quark-line classes. |
| Generic DAG reachability | Implemented for production LC artifacts | `generic_dag.py` builds source currents, model-vertex interactions, external-subset closures, stage buckets, and lowering-readiness summaries from `CurrentIndex` physics state rather than process-family tags. |
| Generic current-plan facade | Unified with production DAG | `current_plan.py` is now a compatibility view over `GenericDAGCompiler`; it no longer owns a second simplified recursion or weak process-family-shaped current key. |
| Shared runtime value types | Split from native reference evaluator | Common dataclasses and errors live in `core_types.py`, so importing the package root, CLI, generic artifacts, phase-space helpers, or process runtime no longer imports `native.py`. |
| Generic schema-v2 process artifact | Implemented as production artifact | `generic_artifact.py` serializes process IR, colour plan, currents, interactions, closures, runtime schema, Symbolica evaluator stages, validation momenta, and standalone checkers into `pyamplicol-generic-dag-process` artifacts. |
| Symbolica evaluator bridge | Split from retired staged-DAG runtime | `symbolica_evaluator.py` owns evaluator settings, chunking, save/load manifests, and compiled/JIT adapter construction for production schema-v2 stages. Generic stage compilation no longer imports the retired `dag_runtime.py` implementation. |
| Matrix generation facade | Generic-only | `matrix.py` is now a process-generic schema-v2 planning facade.  The previous `CurrentKey`/`RecursionGraph` family graph builder is quarantined in `legacy_matrix.py` for reference tests and migration diagnostics. |
| Arbitrary quark-pair counts in enumeration formulas | Implemented | The Python process layer uses generic quark-line combinatorics instead of explicit 1/2/3 branches. |
| Production ME lowering | Implemented and validated for broad generic LC and the documented NLC/full matrix coverage | The schema-v2 generic DAG path reaches model-owned kernels for multi-boson, scalar/Higgs, pure-gluon, and multi-quark-line processes, builds evaluator-stage expressions, serializes evaluator artifacts, and records explicit blockers when a requested process or colour accuracy is not yet implemented. |
| Python artifact loader | Implemented for schema-v2 generic DAG artifacts | `load_process(..., runtime="python")` resolves process sets, subprocess keys, generic external layouts, and stage-plan metadata without depending on process-family labels. Schema-v1 execution is rejected by default and requires the explicit `allow_legacy_schema_v1=True` reference-only opt-in. |
| Rusticol runtime | Implemented for schema-v2 generic DAG artifacts | `rusticol.Runtime.load()` is the production loader and accepts only schema-v2 generic DAG artifacts. On AArch64, the f64 path uses native two-lane complex SymJIT payloads and shares one packed stage-local input across output chunks; `RUSTICOL_NATIVE_SIMD_JIT=0` selects the scalar diagnostic fallback. Schema-v1 family artifacts are reference-only and require the explicit `rusticol.Runtime.load_legacy()` entry point. |
| Compact runnable manifests | Implemented | After evaluator serialization, generation-only interaction detail is reduced to stage IDs and vertex-kind counts; the full records remain in the trusted local `runtime_schema_diagnostics.pickle.gz` sidecar, plan-only manifests stay detailed, and Rusticol loads both layouts. |
| NLC/full-colour | Initial production support, direct Fortran colour-matrix validation through two quark lines | `--color-accuracy nlc` and `--color-accuracy full` build the same colour-ordered DAG and replace the final LC diagonal contraction by a sparse colour metric.  Pure-gluon, one-quark-line, and two-quark-line classes are matched directly to Fortran AmpliCol colour matrices.  Higher quark-line counts use the generic open-line trace-overlap metric, including gluons on the open lines; those cases require raw-amplitude probes rather than Fortran's built-in colour matrix for validation. |
| Production CLI surface | Generic DAG only | `processes`, `process-plan`, `generate-process`, `time-process`, and `compare-amplicol --runtime-backend rusticol` are the visible workflow. `compare-amplicol` defaults to the generated-library supplied-momenta probe (`--library=create`, `make amplicol_generate_library`, `--library=use`). Legacy native/tensor/Z-family commands remain compatibility stubs only and are hidden from help output. |

## Fast-Path Invariant

The phrase "fast path" refers to the currently validated production artifact
route, not to a permanent physics special case.  For one-quark-line
`q q~ > V + n g` processes, pyAmpliCol generates a fixed shared-current sweep,
serializes bounded Symbolica evaluator blocks, and Rusticol executes the
resulting stage plan without dynamic process interpretation.  This is the route
used for the live `d d~ > Z + (n-1)g` timings in `pyAmpliCol.pdf` and
`z_performance_data.json`; `performance_summary.md` is retained only as a
historical benchmark snapshot.

Generic process support must preserve that property.  Parser, colour-flow, and
current-DAG construction can be process-generic at generation time, but the
generated subprocess artifact should still contain compact fixed arrays:
source-current layouts, stage evaluator manifests, current offsets, amplitude
weights, and coherent-group metadata.  For `d d~ > Z + n g`, the generic
compiler is expected to emit the same effective stage layout as the current
validated artifact route, so timing should match within normal system noise.
If a generalization changes the hot runtime into hash-map lookup, dynamic
dispatch, extra live colour sectors, or different evaluator chunk/CSE
boundaries, that is a regression rather than acceptable genericity.

LC generation supports two complementary sharing modes.  The ordinary CLI
default, `--lc-sector-strategy topology-representatives`, compiles one
replay-safe representative for a selected-flow workload.  The all-flow
benchmark path uses `--lc-sector-strategy all`: it retains every LC amplitude
root but keys reusable currents by physics state and ordered coloured word, not
by sector id.  Current and interaction counts therefore follow the shared
recursion rather than multiplying by the number of colour orderings.

`--runtime-lc-sector-ids` can additionally write a pruned selected-flow
sidecar.  Its colour choice is fixed before Symbolica evaluator construction,
so a one-flow timing does not execute an all-flow evaluator or retain dynamic
selector conditionals.  The result matrix deliberately generates both the
shared all-flow artifact and this selected sidecar: generation is compared
against AmpliCol's complete LC capability, while the two runtime columns time
one flow with a helicity sum and all flows at one fixed helicity separately.

The colour planner still records topology groups and exact label
permutations.  `--lc-topology-replay` remains available as a diagnostic
representative-replay route, but it is not required by the materialized shared
all-sector artifact.  The grouping and shared-current keys are process-generic:
they are derived from process IR, particle state, open colour lines, trace
structure, and singlet attachments rather than a `Z+gluons` family tag.

Generic pruning knobs are allowed when they describe physics budgets rather
than process families.  Examples include maximum QCD/QED/EW coupling order,
selected colour-sector ids or topology representatives, active quark flavour
schemes, maximum external quark-pair counts, ignored particles/vertices, and
limits on virtual non-QCD propagators.  Such options must be implemented at
the model/colour/DAG layer; they must never become tags such as "Z+gluons
mode" or "dilepton mode".

The current generic DAG compiler implements the first of these budget layers.
Each `CurrentIndex` carries an accumulated UFO-style coupling-order vector,
currently `QCD` and `QED` for the AmpliCol SM model.  The model classifies a
local vertex by coupling order, the compiler adds that order to the two parent
current indices, and optional caps such as `--max-coupling-order QED=1`,
`--max-qcd-order 4`, or `--max-qed-order 1` prune over-budget branches before a
new current is inserted.  The same parser accepts generic `--ignore-particles`
and `--ignore-vertex-kinds` constraints.  These are deliberately local model
constraints: no code path asks whether the process is "Z+jets",
"dileptons", or any other family.

The compiler also accepts `--max-lc-current-line-groups N`, which caps how many
leading-colour open-line groups one intermediate current may span.  This is a
generic colour-state constraint for controlling multi-quark-line warmups; it is
not tied to a process family.  For example, on `d d~ > u u~`, the uncapped DAG
contains currents spanning two LC line groups and produces amplitude roots,
whereas `--max-lc-current-line-groups 1` keeps every current within a single
line group and intentionally removes the multi-line closure.

The separate `--max-quark-pairs N` option is a subprocess-level cap for
inclusive or exploratory requests.  It rejects over-budget subprocesses before
current insertion without inspecting the process family.  It is useful when a
user wants, for example, `p p > Z + jets` expansion but only up to a chosen
number of open quark lines.

For selected workloads, topology representatives remain the default generic
generation-time sharing knob.  If a requested topology cannot be replayed
safely, the planner falls back to the contributing sectors rather than
returning an approximation.  For complete all-flow workloads, the shared
all-sector DAG handles open-line and pure single-trace gluon processes
directly.  Pure-gluon reflection symmetry and mixed-process gluon-subcurrent
reflection reuse carry their signs on interaction edges; selected-sector
artifacts disable those all-ordering transformations.

Executable replay validation is available through:

```sh
python pyAmpliCol/scripts/validate_lc_topology_replay.py 'd d~ > z g g'
```

This builds a complete shared artifact and a topology-representative artifact
for the same deterministic point, loads both through Rusticol, and compares
their final matrix elements.  The matrix refresh performs the stronger
production check against AmpliCol as well and records selected-flow and
all-flow values, current counts, generation phases, and runtimes in
`result_matrix_data.json`.

## Legacy Process Families

The process parser and model metadata target the full Fortran AmpliCol
vocabulary:

- QCD partons: `g`, `d`, `u`, `s`, `c`, `b`, `t` and antiparticles.
- Electroweak/Higgs singlets: `a`, `z`, `w+`, `w-`, `h`.
- Charged leptons: `e+`, `e-`, `mu+`, `mu-`, `ta+`, `ta-`.
- Neutrinos: `ve`, `ve~`, `vm`, `vm~`, `vt`, `vt~`.

Examples that should enumerate:

```sh
./pyamplicol.sh processes 'd d~ > e+ e- [d g] [d g]' --json
./pyamplicol.sh processes 'u d~ > w+ 2j' --json
./pyamplicol.sh processes 'p p > z z 2j' --json
./pyamplicol.sh processes 'd d~ > z 3*[d g]' --json
```

Process-set artifact generation expands built-in inclusive labels to concrete
subprocess entries before launching child generators.  Thus `p p > z g` is
represented by subprocess keys such as `d_dbar_to_g_z`, `u_ubar_to_g_z`,
`dbar_d_to_g_z`, and `ubar_u_to_g_z`, not by a single inclusive `p_p_to_z_g`
artifact.  This is the artifact shape expected by Rusticol.  The current
vector-plus-gluon lowering supports both `q q~` and reversed `q~ q` beam order
for these concrete entries, including neutral `gamma/Z` examples and
charged-current examples such as `p p > w+ g`.  The native explicit-lepton
recursion uses the same order-agnostic incoming-current checks. Other
unsupported entries still report clear diagnostics rather than falling back to
an inclusive approximation.

The root process-set manifests written by both `process-plan` and
`generate-process` record a `generic_generation` block with the process
options, process-neutral pruning controls, LC sector strategy, topology replay
flag, and worker count where relevant. New representative entries and crossing
aliases in both planning manifests and executable process-set artifacts also
carry a `generation_request` copy of this metadata, while entries kept during
`--append` retain their original metadata. This makes a process set auditable
without reading each nested subprocess manifest first, and prevents generic
options such as
`--max-quark-pairs`, `--max-qed-order`, or `--lc-sector-ids` from becoming
implicit state.

When several subprocesses in the same process-set request are related by
crossing or by a permutation of identical all-outgoing particle content,
`generate-process` builds only the first representative artifact. Later entries
become manifest aliases pointing at that representative path.  The alias stores
an input map with `source_index`, `target_index`, and a crossing sign, so
Rusticol can accept momenta in the selected subprocess's physical external-leg
order and internally evaluate the shared representative artifact.  Appending to
an existing process set uses only real representative entries to seed this
reuse cache; existing aliases are never selected as new representatives.
Append mode preserves the existing process order and default process key, so
adding subprocesses cannot silently change which entry Rusticol evaluates when
no explicit process key is supplied.  Replace mode rebuilds entries in the new
request order.

Concrete leading-colour processes with implemented model lowerings can be
generated as executable Rusticol process artifacts:

```sh
./pyamplicol.sh generate-process 'd d~ > a g' outputs/dd_a_g
./pyamplicol.sh time-process outputs/dd_a_g
./pyamplicol.sh generate-process 'u d~ > w+ g' outputs/ud_w_1g
./pyamplicol.sh time-process outputs/ud_w_1g
./pyamplicol.sh generate-process 'd d~ > e+ e-' outputs/dd_epem_0g
./pyamplicol.sh time-process outputs/dd_epem_0g
./pyamplicol.sh generate-process 'd d~ > e+ e- g' outputs/dd_epem_1g
./pyamplicol.sh time-process outputs/dd_epem_1g
./pyamplicol.sh generate-process 'u d~ > e+ ve g' outputs/ud_epve_1g
./pyamplicol.sh time-process outputs/ud_epve_1g
./pyamplicol.sh generate-process 'd d~ > u u~ s s~' outputs/dd_uuss
./pyamplicol.sh time-process outputs/dd_uuss
```

Schema-v2 `generate-process` requests for generic processes write executable
artifacts when LC current construction, model lowering, evaluator generation,
and Rusticol runtime loading are all available. Unsupported runtime requests
produce an explicit diagnostic naming the missing layer: current closure,
vertex/propagator lowering, colour expansion, evaluator serialization, or
Rusticol integration.

For process-set generation, parent-level generic options are forwarded to
subprocess child generators. This includes `--color-accuracy`, evaluator
backend/preset options, split-stage settings, and process-enumeration options,
so a multi-entry artifact cannot silently fall back to an LC or legacy
family-specific child mode.

Rusticol's production API intentionally rejects schema-v1 artifacts, and the
Python loader now does the same unless `allow_legacy_schema_v1=True` is passed
explicitly for reference-only diagnostics.  Those artifacts still exist in
tests and migration diagnostics because they were used to validate the original
`q q~ > V + n g` eager-DAG route, but they are loaded only through explicit
legacy entry points such as `rusticol.Runtime.load_legacy()` or the Python
opt-in above.  New process-generation work must target schema-v2 manifests and
must not add new production behavior to the schema-v1 family loader.

## Current LC Validation Snapshot

The generic schema-v2 LC validator reports the full default smoke matrix as
supported in dry-run support mode.  The same 58-process default matrix has
also been run numerically against Fortran AmpliCol generated-library probes
with the 30 GB memory watchdog; all entries passed with maximum relative
difference `1.02e-12` and peak RSS `0.31 GB`.  The validation creates the
Fortran library through `./amplicol_generate --library=create`, warms the
helicity filter and generated library emission with deterministic supplied
momenta, compiles it with `make amplicol_generate_library`, and then compares
against `./amplicol_generate --library=use` on the same supplied point.  This
avoids reference-side phase-space-integration failures while preserving the
generated-library chain used for benchmarking.  This validation uses the
original Fortran-side `process_list.py` writer as the reference process-list
backend, with the legacy `-cc` and `-3` switches enabled for charged-current
conjugates and three-quark-pair subprocesses.  The broad default validation
driver enables those legacy switches automatically for the built-in matrix;
explicit custom process lists keep explicit switch control.
Supplied momenta are reordered to the concrete process-file PDG order before
both library creation and library-use probes, so requests such as
`d d~ > z g` that the Fortran process list canonicalizes as `d d~ > g z`
still compare identical on-shell external assignments.

A representative breakdown of the covered low-multiplicity classes is:

| Class | Processes Checked | Max Relative Difference |
|---|---:|---:|
| neutral vector | `d d~ > z`, `d d~ > z g`, `d d~ > z g g` | `5.99e-16` |
| charged current | `u d~ > w+`, `u d~ > w+ g`, `d u~ > w-`, `d u~ > w- g`, `u d~ > e+ ve g`, `d u~ > e- ve~ g`, `u/d charge-conjugate W/Z/a/H channels` | `1.02e-12` |
| dilepton | `d d~ > e+ e-`, `d d~ > e+ e- g`, `d d~ > e+ e- a` | `4.29e-15` |
| multi-boson/Higgs | `d d~ > z z`, `d d~ > a z g`, `d d~ > h z`, `u d~ > h w+` | `1.02e-12` |
| pure/mixed QCD | `g g > g g`, `g g > g g g`, `g g > u u~ g`, `d d~ > u u~ g` | `2.55e-15` |
| massive fermions | `g g > t t~`, `g g > t t~ g`, `g g > t t~ h`, `d d~ > t t~`, `d d~ > t t~ g` | `5.54e-16` |
| multi-quark lines | `d d~ > u u~ s s~`, via the Fortran-selected concrete subprocess `u u~ > s~ s d~ d` | `1.47e-15` |

The generic current planner has also been stress-tested beyond the old
Fortran-script explicit quark-line cases.  For
`d d~ > u u~ s s~ c c~`, the canonical LC colour planner keeps the 24
line-pairing sectors relevant for the leading-colour representative selection.
A selected-sector process-plan probe builds 520 currents, 2120 interactions,
and 64 amplitude roots with no truncation and no missing vertex or propagator
lowerings.  This is a structural LC planning and lowering-readiness check; full
all-sector Rusticol numerical validation for arbitrary quark-line counts
remains separate work.

These comparisons use deterministic supplied-momenta probes through the
generated Fortran library workflow, `./amplicol_generate --library=create`,
`make amplicol_generate_library`, and `./amplicol_generate --library=use`.
The direct AmpliCol probe normalization now uses the phase-space group external
count `pgl(ichan)%next` rather than the module-level `next` variable.  This is
important in `--library=use` mode, where the module-level value can be zero
after loading a generated library.  The fixed library probe has been
rechecked on focused low-multiplicity representatives:

| Process | Relative Difference vs Fortran Library Probe |
|---|---:|
| `d d~ > z g` | `4.09e-16` |
| `d d~ > e+ e- g` | `2.98e-15` |
| `u d~ > w+ g` | `7.25e-15` |
| `d d~ > u u~` | `1.40e-16` |
| `d d~ > u u~ g` | `1.95e-15` |
| `d d~ > u u~ s s~` | `1.47e-15` |
| `g g > t t~` | `5.54e-16` |
| `g g > t t~ h` | `1.34e-16` |
| `d d~ > t t~` | `1.34e-16` |
| `g g > t t~ g` | `1.84e-16` |
| `d d~ > t t~ g` | `2.01e-16` |
| `u d~ > e+ ve z` | `5.24e-16` |
| `u d~ > e+ ve a` | `6.84e-16` |
| `d u~ > w-` | `1.24e-15` |
| `d u~ > e- ve~` | `2.27e-15` |
| `d u~ > e- ve~ g` | `3.17e-16` |
| `d u~ > e- ve~ z` | `2.07e-15` |
| `d u~ > e- ve~ a` | `3.50e-16` |
| `d u~ > h w-` | `1.02e-12` |
| `d u~ > w- z` | `1.09e-13` |
| `d u~ > w- a` | `7.32e-15` |

## Current NLC/Full-Colour Validation Snapshot

The NLC/full-colour implementation has been smoke-tested against the dedicated
Fortran `amplicol_color_probe` generated-library driver on the first
representative colour classes.  The driver builds the same AmpliCol generated
library used for benchmark-quality LC probes, initializes `init_col` with the
requested colour accuracy, and evaluates supplied momenta through the generated
amplitude library.

| Class | Process | Colour Accuracy | Relative Difference |
|---|---|---:|---:|
| one quark line | `d d~ > z g` | NLC | `1e-16` level |
| one quark line | `d d~ > z g` | full | `1e-16` level |
| pure gluon | `g g > g g` | NLC | `1e-16` level |
| pure gluon | `g g > g g` | full | `1e-16` level |
| two quark lines | `d d~ > u u~` | NLC | `1e-16` level |
| two quark lines | `d d~ > u u~` | full | `1e-16` level |

The broader result-matrix campaign for NLC/full colour is still in progress.
For more than two quark lines at NLC/full, pyAmpliCol now builds and runs the
generic sparse open-line colour metric.  Fortran AmpliCol's built-in colour
matrix still stops there, so those result-matrix reference slots are N/A until
the raw-amplitude validation probe is used for that comparison.

The current validation command for the default matrix is:

```sh
pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30 \
  --report-json pyAmpliCol/outputs/generic_lc_default_matrix_watch.json -- \
  pyAmpliCol/dependencies/.venv/bin/python \
  pyAmpliCol/scripts/validate_generic_lc_against_amplicol.py \
  --output-dir pyAmpliCol/outputs/generic_lc_default_matrix \
  --symbolica-evaluator-backend jit --batch-size 16 \
  --n-cores 4 --jobs 1 --json
```

## Model Lowering Coverage

The model now exposes `vertex_lowering_coverage()` so process-support
diagnostics and generic DAG construction can identify blockers from the
legacy vertex table itself.  The current split is:

- tensor-network or direct-symbolic ready: all legacy vertex kinds `0` through
  `24`, including the QCD U(1) subtraction helper `8` and the massive
  Dirac-fermion scalar/Yukawa current `16`;
- pending: none in the current leading-colour model metadata;
- not yet implemented: none in the current leading-colour model metadata.

Future model extensions should continue adding missing physics through
model-owned lowering rules instead of process-family branches.

The generic colour planner and current planner already use the canonical
all-outgoing process IR.  The colour planner enumerates leading-colour open
quark-line sectors, including arbitrary balanced quark-pair counts and ordered
gluon allocations, plus pure-gluon single-trace sectors.  For NLC/full colour,
the same planner retains the additional non-folded pure-gluon trace orderings
and open-line block permutations needed for sparse colour contractions.  The
current sparse metric contracts generic products of open fundamental strings,
so arbitrary balanced quark-line sectors with gluons are represented by the
same colour-flow trace-overlap rule.  The generic DAG
compiler uses those LC sector ids in `CurrentIndex`, duplicates source/current
tables per sector, rejects cross-sector current combinations, and lets coloured
currents span several open line groups inside one LC sector when model vertices
and colour flow allow it.  Colourless singlet labels remain attachable to any
compatible colour line.  As a result, low-multiplicity multi-quark examples
such as `d d~ > u u~` and guarded larger examples such as
`d d~ > u u~ s s~` now find amplitude closures in the generic planning layer.
The normal unit suite also checks a selected leading-colour sector for the
four-quark-line process `d d~ > u u~ s s~ c c~`; the colour planner enumerates
24 canonical LC line-pairing sectors, and the selected-sector DAG finds
amplitude closures without truncation.  This is a planning/scalability
invariant rather than a direct Fortran colour-matrix validation gate.  NLC/full
colour for more than two quark lines is therefore a pyAmpliCol generic-colour
capability, while independent Fortran comparison for those cases requires the
raw colour-order amplitude probe because Fortran AmpliCol's own `init_col`
matrix only covers zero, one and two quark-line sectors.
The planner also filters out the `99` QCD singlet helper in leading-colour mode
and answers whether a concrete process is reachable through model vertices and
which vertex lowerings block a complete generic evaluator.  These layers now
build symbolic stage expressions for generic current and amplitude stages and
serialize bounded Symbolica evaluator artifacts consumed by Rusticol.

The planner also exposes a generic stage plan.  Current interactions are
bucketed by external-subset size, matching the AmpliCol topological sweep, and
closures are collected into a final amplitude stage.  Each stage carries
current/interaction ids, per-colour-sector work summaries, and the ready,
pending, and unimplemented vertex kinds.  The per-sector summaries include
current, interaction, and closure counts plus vertex-kind readiness, so the
generic compiler can schedule LC sectors without reverse-engineering the raw
current ids.  This is the metadata shape consumed by the process-generic
Python and Rusticol execution schema.

The schema-v2 artifact additionally records a `compiled.stage_compiler`
summary and, when requested, serialized stage evaluators.  It constructs
bounded stage expressions from the runtime value storage and real momentum
parameters, asks the model for each local vertex and propagator expression,
accumulates contributions by `CurrentIndex`, and emits output-slot metadata for
Rusticol.  Current examples such as `d d~ > z z g`, `d d~ > w+ w- g`,
`g g > g g`, and `d d~ > u u~ s s~` have expression-ready stage plans with no
lowering blockers.  Low-multiplicity non-`Z+gluons` smoke checks have also
been pushed through serialized schema-v2 evaluator generation and Rusticol
evaluation: `d d~ > e+ e- g`, `u d~ > w+ g`, and selected-sector
`d d~ > u u~ g`.  These are runtime availability checks, not yet replacements
for the required Fortran-AmpliCol matrix-element agreement campaign.

The schema-v2 planning manifest is available through
`build_generic_process_manifest()` and `write_generic_process_manifest()`.  It
records both the physical external PDG order used by runtime momenta and the
crossed all-outgoing PDG order used by current construction, together with the
leading-colour sector plan consumed by the generic stage compiler and Rusticol
loader.  The `process-plan` manifest is not a production
Rusticol artifact because it intentionally contains no serialized Symbolica
evaluators.  The `generate-process` command now writes the closely related
schema-v2 `pyamplicol-generic-dag-process` artifact shape for concrete
requests, including serialized evaluator stages when generation succeeds.

The same payload can be written from the CLI without invoking evaluator
generation:

```sh
./pyamplicol.sh process-plan 'd d~ > z g' outputs/plans/dd_z_g
./pyamplicol.sh process-plan 'd d~ > z g | d d~ > z z g' outputs/plans/set
./pyamplicol.sh process-plan --color-accuracy full 'd d~ > z z g' outputs/plans/full
```

A single concrete process writes `generic_process_manifest.json` directly in
the output directory.  A request that expands to several partonic subprocesses,
including `p`, `j`, anonymous labels, or explicit `PROC | PROC` sets, writes a
root `generic_process_set_manifest.json` plus nested
`subprocesses/<canonical-key>/generic_process_manifest.json` files.  The
planning command is intentionally broader than `generate-process`: it accepts
LC/NLC/full-colour metadata and records lowering blockers instead of rejecting
processes whose production evaluator is not implemented yet.

Unsupported `generate-process --json` failures include a `support_report`
payload.  When the concrete process is small enough for the guarded planner,
that payload contains `color_plan` and `current_plan` summaries.  The colour
summary reports the leading-colour sector count, truncation status, sector-kind
counts, and whether Idenso basis/metric generation is required for the requested
colour accuracy.  The current summary reports current/interaction/closure
counts, the colour-sector ids actually used by current identities, and the
per-sector vertex-kind worklists, plus the exact pending or unimplemented
vertex kinds.  For example, `d d~ > z z g` and `d d~ > w+ w- g` now report
reachable, tensor-ready generic current graphs using the model-owned
multi-vector and scalar kernels. Multi-quark-line QCD examples report their
open-line colour-sector count, current/interactions, and amplitude closures.
Their remaining hard blocker is no longer family-specific reachability; the
next gate is broad numerical agreement against Fortran AmpliCol, including
normalization and LC colour weights, across representative concrete processes.

Long JIT evaluator builds now emit heartbeat progress while Symbolica is inside
the blocking evaluator-construction call.  The heartbeat cannot expose
Symbolica's internal optimization phase boundaries, but it keeps the
progressbar alive with elapsed `jit initialize` updates and refreshed process
RAM while that call is running.

`generate-process --monitor` exposes the same DAG construction, pruning,
stage-blueprint, evaluator-build, diagnostics, and manifest phases in non-TTY
logs. Rapid per-chunk JIT transitions are coalesced into periodic updates, so
captured logs remain useful without losing the selected autotune layout or
phase-completion timings.

## Colour Accuracy

Leading colour remains the default production setting:

```sh
./pyamplicol.sh generate-process --color-accuracy lc 'd d~ > Z 4*g' outputs/dd_z_4g
```

`--color-accuracy nlc` and `--color-accuracy full` now use the same
colour-ordered generic DAG and replace the LC diagonal contraction by a sparse
colour-contraction block in the Rusticol amplitude reducer.  The contraction
plan is built from the colour sectors emitted by `color_plan.py`:

- LC uses the existing diagonal leading-colour weights.
- NLC keeps the absolute AmpliCol NLC-expanded value, not only the correction
  term.
- Full colour keeps all nonzero sparse colour-metric entries.

The initial NLC/full implementation covers the colour classes for which the
Fortran reference is available and manageable: pure gluon, one quark line, and
two quark lines.  More than two quark lines remains LC-only in production and
reports a structural unsupported diagnostic for NLC/full.  A dedicated
`amplicol_color_probe` driver provides generated-library Fortran reference
values for NLC/full comparisons and exposes colour-matrix diagnostics when
needed.

The long-term target remains Idenso-backed arbitrary-quark-line colour bases,
but this is deliberately staged after Fortran-parity validation of the
implemented NLC/full classes.

## Validation Policy

For matrix-element support beyond `Z + gluons`, validation should start with
low-multiplicity processes with six or fewer final-state particles:

- `d d~ > Z + n g`;
- `d d~ > e+ e- + n g`;
- `u d~ > W+ + n g`;
- `d d~ > Z Z + n g`;
- photon/Z/W/H singlet combinations supported by legacy vertices;
- multi-quark-line QCD/EW examples through the Fortran-supported range.

For quark-line counts beyond Fortran AmpliCol's effective validation range,
use Python-vs-Rusticol parity, Spenso/direct-symbolic lowering parity, Idenso
colour-flow checks, current-table reachability, helicity-filter round trips,
and normalization invariants.

The reusable low-multiplicity LC validation driver is:

```sh
dependencies/.venv/bin/python scripts/run_with_memory_watch.py --limit-gb 30 -- \
  dependencies/.venv/bin/python scripts/validate_generic_lc_against_amplicol.py
```

By default this validates the broad low-multiplicity process matrix listed in
the source constant `DEFAULT_LC_VALIDATION_PROCESSES`, automatically enabling
the legacy charged-current and three-quark-line process-list switches needed by
that matrix. It uses deterministic supplied-momenta Fortran probes through the
generated-library path for point-by-point matrix-element validation. If
Fortran runtime timing is needed, pass `--library-timing-events N`; this reuses
the generated library for a warm `./amplicol_generate --library=use` timing
run.

Explicit custom process lists keep explicit control over the same
process-enumeration switches needed for legacy corners:

```sh
dependencies/.venv/bin/python scripts/validate_generic_lc_against_amplicol.py \
  --include-3qqbar 'd d~ > u u~ s s~'
```

Three-quark-line processes without extra gluons are part of the light
validation matrix. Three-quark-line plus gluon checks are useful stress tests
for colour-sector growth, but are not currently routine gates because generic
current planning becomes much more expensive before sector-growth
optimizations are applied.

## Validation Snapshots

These are focused checks that document the current incremental coverage. The
matrix-element agreement entries below use deterministic supplied-momenta
Fortran probes through the generated-library workflow
`./amplicol_generate --library=create`, `make amplicol_generate_library`, and
`./amplicol_generate --library=use`.

| Process | Layer | Reference | Status |
|---|---|---|---|
| `d d~ > Z + n g` | Rusticol eager-DAG artifact | Fortran AmpliCol | Production benchmark fast path; see `performance_summary.md`. |
| Default LC validation matrix, 58 processes | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Passes under the 30 GB watchdog with max relative difference `1.02e-12` and peak RSS `0.31 GB`. |
| `d d~ > e+ e- + n g`, `n = 0,1,2` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences `3.66e-16`, `1.21e-15`, and `1.34e-15`. |
| `d d~ > mu+ mu- g`, `d d~ > mu+ mu- a`, `u d~ > mu+ vm g` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences between `1.21e-15` and `4.29e-15`, confirming charged-lepton flavour genericity. |
| `u d~ > W+ + n g`, `n = 0,1,2` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences between `8.20e-16` and `1.74e-15`. |
| `u d~ > e+ ve + n g`, `n = 0,1`, plus `u d~ > e+ ve z/a`; charge conjugates `d u~ > e- ve~ + n g`, `n = 0,1`, plus `d u~ > e- ve~ z/a` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences from `3.17e-16` to `2.27e-15`. |
| `d u~ > w-`, `d u~ > h w-`, `d u~ > w- z`, `d u~ > w- a` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences from `1.24e-15` to `1.02e-12`, matching the existing `W+` side within the validation tolerance. |
| `d d~ > Z Z`, `d d~ > Z Z g`, `d d~ > Z Z Z` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences `1.15e-15`, `2.18e-15`, and `1.27e-15`. |
| `d d~ > a Z`, `d d~ > a a`, `d d~ > a a g`, `d d~ > a Z g`, `d d~ > a Z Z` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences from `2.25e-16` to `1.44e-15`. |
| `d d~ > H Z`, `d d~ > H Z g`, `d d~ > H H Z`, `u d~ > H W+`, `u d~ > W+ Z`, `u d~ > W+ a` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences from `1.53e-16` to `1.02e-12`; all remain well below the `1e-8` LC validation tolerance. |
| `g g > g g`, `g g > g g g`, `g g > u u~`, `g g > t t~`, `g g > t t~ g`, `g g > t t~ H`, `d d~ > t t~`, `d d~ > t t~ g` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences from zero to `2.55e-15`, including massive top-pair, top plus gluon, and top-Higgs checks. |
| `g g > u u~ g`, `g g > d d~ g`, `u d~ > u d~ g`, `u u~ > d d~ g` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences between `4.06e-16` and `1.85e-15`, including anti-oriented open-line colour orders and same-flavour two-line exchange. |
| `d d~ > u u~`, `d d~ > u u~ g`, `g g > u u~ d d~`, `d d~ > u u~ s s~` | Generic DAG, Rusticol schema-v2 artifact | Fortran AmpliCol generated-library supplied-momenta probe | Agreement at relative differences between `6.52e-16` and `1.85e-15`, including the light three-quark-line example selected by Fortran as `u u~ > s~ s d~ d`. |
| `d d~ > u u~ s s~ c c~` | Generic DAG selected LC sector | Internal planning/scalability invariant | One representative sector builds 520 currents, 2120 interactions, and 64 amplitude roots without truncation; broader all-sector evaluation is intentionally not a light validation gate. |
| `d d~ > ve ve~`, `d d~ > ve ve~ g` | Unsupported for Fortran-parity LC validation | Fortran AmpliCol generated-library supplied-momenta probe | The legacy lepton vertex table does not contain neutral-current neutrino-`Z` vertices, so pyAmpliCol intentionally omits them in the AmpliCol-parity model and these requests are not validation gates. |
