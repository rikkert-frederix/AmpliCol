# pyAmpliCol Agent Notes

## Non-Negotiable Sandbox Dispatch

This is the first shell-operating rule for this repository: never leave the
thread waiting on a sandbox approval prompt. The `exec_command` parameter
`sandbox_permissions: "require_escalated"` is banned for this repository. Do
not set it, do not request escalation, and do not ask the user for permission
through the tool approval path. The user has provided and started `listener.py`
specifically so sandbox-sensitive commands can be dispatched through the
`run.order`/`run.output` file protocol instead of blocking the one-shot matching
work.

Mechanical fail-closed guard before every shell tool call:

- The `sandbox_permissions` key must be absent from every `exec_command` call.
  Do not set it to `require_escalated`, and do not set it to `use_default`.
- Before sending a shell command, inspect the intended tool payload. If it
  contains `sandbox_permissions`, delete that key and re-plan the command.
- If the reason for adding `sandbox_permissions` would have been "this may need
  approval", do not run that command through escalation. Dispatch it through
  `listener.py` with `run.order`/`run.output` instead.
- Treat the current Codex approval policy as fail-closed: approval escalation is
  unavailable, so any payload containing `sandbox_permissions` is an invalid
  command plan for this repository.

Use this rule before every command that may touch sandbox-sensitive state:

- Ordinary read-only inspection commands and focused tests may run directly.
- Any command that is likely to write `.git` metadata or otherwise hit the
  known restricted-filesystem boundary must go through `listener.py` first.
  This includes `git commit`, `git rebase`, `git merge`, lock cleanup,
  repository maintenance, and any retry after `Operation not permitted`.
- Do not request or wait for a user approval prompt for these repository
  commands. Never retry them with `sandbox_permissions: "require_escalated"`.
  If a direct command unexpectedly fails with `Operation not permitted` or a
  read-only `.git` error, immediately queue the same command through
  `run.order` and inspect `run.output`.
- Queue one command at a time:

```sh
printf 'git commit -F /tmp/commit-message\n' > run.order
sleep 1
tail -200 run.output
```

The listener is only a sandbox-dispatch mechanism. Long or memory-sensitive
tests, matching previews, CDE/vakint probes, and validation fixtures must still
run through `scripts/run_with_memory_watch.py --limit-gb 30`.

## Running Tests

Install the managed dependencies first. The GammaLoop Python API is not built by
default:

```sh
python dependencies/install_dependencies.py
```

Request the optional GammaLoop Python API explicitly when needed:

```sh
python dependencies/install_dependencies.py --with-gammaloop
```

Then run the test suite with pytest from the managed virtual environment:

```sh
source dependencies/.venv/bin/activate
python -m pytest tests
```

Equivalently, without activating the environment:

```sh
dependencies/.venv/bin/python -m pytest tests
```

The pytest suite includes a static typing check that runs `python -m mypy`.
To run it directly:

```sh
dependencies/.venv/bin/python -m mypy
```

Prefer grouped targeted tests during implementation slices, and reserve the
full suite for larger green milestones:

When a slice touches the slow validation fixtures, batch related work first and
then run the validation group once. Do not pay for the full suite after every
small local fix; use focused tests while building the slice, then a broader
targeted gate, and only then a full-suite gate when the milestone is large
enough to justify it.

Heavy AmpliCol or pyamplicol validation runs must always run behind the memory
watchdog with a 30 GiB limit. This includes large `make`/generation/evaluation
loops, AmpliCol reference comparisons, multi-point `ME_TEST` runs, tensor
network reduction, Symbolica optimization/JIT, and performance sweeps. Never
launch those workloads directly from the shell when they could grow large.
Use:

```sh
source "$HOME/.bashrc"
dependencies/.venv/bin/python scripts/run_with_memory_watch.py --limit-gb 30 -- \
  <heavy command and arguments>
```

If the heavy command must run from the repository root instead of `pyAmpliCol`,
invoke the same watchdog by absolute or relative path, for example:

```sh
source "$HOME/.bashrc"
pyAmpliCol/dependencies/.venv/bin/python pyAmpliCol/scripts/run_with_memory_watch.py --limit-gb 30 -- \
  make -j8 amplicol_generate_library
```

Always use the managed virtual environment for pyAmpliCol development and tests.
Do not use the ambient system Python when importing Symbolica, idenso, spenso,
or vakint.

Before running Python commands that import Symbolica, source `~/.bashrc` so the
local `SYMBOLICA_LICENSE` export is available:

```sh
source "$HOME/.bashrc"
dependencies/.venv/bin/python -m pytest tests
```

Run tests or exploratory  workloads that could grow large through the
project memory wrapper with a 30 GiB cap:

```sh
source "$HOME/.bashrc"
dependencies/.venv/bin/python scripts/run_with_memory_watch.py --limit-gb 30 -- \
  dependencies/.venv/bin/python -m pytest tests/integration/matching
```

Use this wrapper for broad pytest groups, slow validation fixtures, CDE/vakint
matching smokes, and any workload that might approach machine RAM limits.
The wrapper also polls `stop.order` in the current working directory by
default. For long exploratory commands, prefer this file-based stop mechanism
over process-management commands that may require sandbox approval: remove any
stale `stop.order` before launching the workload, and create/touch
`stop.order` to ask the wrapper to terminate the wrapped process group.

## Sandbox And Listener Workflow

Do not let work stall on sandbox approval prompts. The expected workflow is:

1. Run ordinary read-only shell commands and ordinary tests directly.
2. For commands that may touch `.git` metadata or another known
   restricted-filesystem boundary, use the user-started `listener.py` route
   immediately. If any other necessary direct command fails with
   `Operation not permitted`, a sandbox write restriction, or the known
   read-only `.git` metadata failure, retry it through the listener rather
   than stalling on an approval prompt.
3. Queue exactly one command by writing it to `run.order`, then read
   `run.output` for its exit code and output. The listener clears `run.order`
   itself and appends history to `run.log`.

Typical listener usage:

```sh
printf 'git commit -F /tmp/commit-message\n' > run.order
sleep 1
tail -200 run.output
```

Use the listener for any `.git` metadata writes that fail in the sandbox,
including `git commit`, `git rebase`, `git merge`, lock cleanup, and similarly
blocked repository-maintenance commands. For long or memory-sensitive
Python/test/matching workloads, still use
`scripts/run_with_memory_watch.py --limit-gb 30`; the listener is a
sandbox-dispatch fallback, not a replacement for the memory watchdog.

## Logging And Progress Output

Use pyamplicol's package logging layer for user-facing progress and debugging
output. Do not add ad hoc `print(...)` calls inside library code. Exported
helpers are available as:

```python
import pyamplicol

pyamplicol.configure_logging()
```

Use `pyamplicol.logging.get_logger(...)` and `pyamplicol.logging.progress(...)` from
implementation modules. Log high-level, notebook-friendly progress at `INFO`
for expensive matching, validation, tensor-reduction, and integral-evaluation
steps. Use `DEBUG` for lower-level internals. Never log full large Symbolica
expressions by default; log stage names, backend choices, counts, and timings.

## Dependency Policy

Using `sympy` and `scipy` is strictly forbidden in this project, including
importing them. All symbolic, algebraic, tensor, and integral work must be done
with Symbolica and the locally built community modules as much as possible.

Use:

- Symbolica for all symbolic and algebraic manipulations.
- idenso for gamma-matrix and colour algebra.
- spenso for tensor-network evaluations when needed.
- symbolica jit_compiled evaluators for efficient low-level implementation of the matrix element

## Matrix-Element Architecture

pyamplicol is a modern matrix-element generator for AmpliCol-style
leading-colour recursion. Keep the implementation architecture-general, but
stage numerical validation first on the `q q~ -> Z + n g` family through at
least six final-state gluons.

Strict ownership rules:

- The process layer owns parsing, crossing conventions, subprocess records,
  colour orders, phase-space groups, multichannel partners, symmetry factors,
  and legacy `processes.txt` export. `processes.txt` is a compatibility output,
  not the internal API.
- A `Model` owns all model-dependent physics: particles, vertices, couplings,
  masses, widths, colour representations, tensor registration, auxiliary
  particles, normalization constants, and Feynman-rule lowering. Do not hard-code
  Feynman rules, coupling constants, colour factors, or auxiliary tensor data
  outside model classes.
- The native matrix-element generator owns recursion graphs, current identity
  keys, chirality propagation, colour-order compatibility, propagators, final
  closure, tensor-network construction, evaluator cache metadata, and timing
  breakdowns.
- The AmpliCol adapter owns all legacy process steering through
  `subprocess.Popen`. Do not call `run.sh` from pyamplicol library code.

The first native implementation strategy is one whole symbolic/tensor network:
external wavefunctions and cached polarizations stay at the evaluator boundary,
while recursion, contractions, colour algebra, simplification, expression
generation, and JIT evaluation are delegated to spenso, idenso, and Symbolica.
Do not switch to nested Python numerical current recursion unless a measured
milestone shows the whole-network strategy is not viable.

For the production `q q~ -> Z + n g` D-mode milestone, stick very close to the
ideas behind legacy AmpliCol conceptually. Build one shared helicity-aware
current table, not one scalar evaluator per full helicity assignment. External
helicity states are source currents with unique `ext_cur`-style source bits;
internal currents are identified by particle type, chirality, external subset,
and source-current ancestry; amplitudes are endpoint pairs equivalent to
legacy `curr2amp`. The evaluator must produce the retained amplitude vector
from one shared recursion sweep.

The scalable D-mode implementation is a staged AmpliCol-like sweep: source
wavefunctions fill a `val_c`-like current table, each current-size layer is
lowered to bounded Symbolica/spenso-backed evaluator outputs, those outputs
are written back into the current table, and a final amplitude evaluator
contracts endpoint currents. Do not reintroduce a Python loop over full
helicity configurations in the hot recursion path. Avoid one fully inlined
monolithic expression for high multiplicity when it exceeds the watchdog or
causes expression blow-up; prefer staged evaluator layers that preserve the
legacy current-table structure.

The `--compiled-dag-evaluator` mode is a separate full-symbolic route. It must
reuse the same AmpliCol-style shared current table, but move runtime current
orchestration to generation time: source current components and real momentum
sums are runtime parameters, while internal current components are opaque
Symbolica aliases and the final retained helicity amplitudes are one
multi-output evaluator. Final complex conjugation, helicity/color factors, and
the ME sum remain outside the evaluator until the multi-output amplitude vector
is validated. Keep this route isolated in compiled-DAG modules while factoring
shared current-table and expression utilities only when that removes real
duplication.

For compiled-DAG lowering, support both `spenso` and `symbolic` routes.
`spenso` is the intended default and should build each current body as a small
tensor network using the fixed model-owned `hep_lib` plus explicit temporary
parent-current tensors whose entries are aliases. `symbolic` constructs the
same current components directly with the verified Symbolica component
builders and is the debugging baseline. Current aliases and temporary parent
tensors are graph construction artifacts; do not register them as new model
physics tensors.

Manual Symbolica aliasing is required for the final compiled-DAG path. Prefer a
managed `symbolica_mod` fork over pyamplicol-side evaluator workarounds, and keep
Symbolica fork changes as small as possible by weaving in existing alias/dev code
rather than rewriting evaluator internals. The currently installable branch is
`pyamplicol-dev-base`, chosen because it is compatible with current
`spenso`/`idenso`; it exposes a minimal Python evaluator hook,
`aliases=[(handle, body), ...]`, plus `Expression.alias(handle=None,
opaque=True)` for constructing such pairs. pyamplicol uses this interface for
current-component alias handles. Compiled-DAG metadata must report whether that hook is active
(`symbolica_alias_available`) and how many current components are represented by
aliases (`opaque_alias_component_count`). Tests for this mode should prove that
aliases survive evaluator construction as reusable DAG slots, and the dependency
installer should remain the source of truth for fetching the selected fork
branch.

Compiled-DAG evaluator construction should apply the warmup helicity filter
before final roots are sent to Symbolica. The filter must be deterministic,
stored in evaluator artifacts, and should use the pure Python/orchestrated
recursion path for signatures rather than compiling another D-mode evaluator.
The default warmup point source is deterministic RAMBO-style on-shell
phase-space generation; fixed canonical pyAmpliCol points remain available for
sanity checks and debugging. Filtering trims roots and records multiplicities;
it must not change the leading-color normalization or move the final
helicity/color sum into the evaluator.

Expose both evaluator-build strategies. `--no-merge-evaluators-strategy` is the
current default and may build a bounded multi-output evaluator in one call,
which has measured faster generation for the current D-mode path.
`--merge-evaluators-strategy` remains available as the memory-saving strategy
and builds bounded evaluator pieces before merging them. `--verbose-evaluator-build`
may enable progress bars and backend verbosity. `--batch-size` defaults to 16
and must control how many phase-space samples are sent to evaluator calls
together. `--no-inlined-helicity-sum` exposes retained raw amplitudes without
changing the shared-current generation strategy.

When reporting performance for the shared D-mode, separate full pyamplicol
wall time from evaluator-only runtime. The milestone comparison against
legacy AmpliCol should use evaluator-only pyamplicol timing for the primary
runtime ratio, after warming up JIT-compiled evaluators. Still report full
runtime overhead separately so slow wavefunction or parameter packing remains
visible.

Rusticol is a PyO3 runtime that links the local Symbolica Rust crate while
pyamplicol also uses Symbolica through the Python extension. Until the linking
model is unified, keep rusticol extension-load/error-path tests isolated in
subprocesses when the same pytest process also constructs Python-side
Symbolica graph objects. Direct runtime use through `rusticol.Runtime.load(...)`
is still the target API, but tests that intentionally exercise failed loads or
unsupported manifests should not poison the main pytest process before later
Symbolica construction tests.

For compiled Symbolica backends, distinguish the valid scalar compiled modes
from SIMD complex mode. On Apple Silicon/aarch64, `compiled-complex-4x` is not
ABI-compatible with Symbolica's current C++ export because the Rust side uses
`Complex<wide::f64x4>` while `xsimd::batch<std::complex<double>,
xsimd::best_arch>` has two complex lanes. Guard this mode instead of accepting
NaN results. Valid modes on this platform are `compiled-complex` with
`--symbolica-compiled-inline-asm default` and `compiled-complex` with
`--symbolica-compiled-inline-asm none`; benchmark both because generic C++ may
run faster but can have much larger generation/compile time.

For colour algebra, use idenso primitives and `simplify_color` whenever possible.
For tensor contractions and Feynman-rule building blocks, register model-owned
tensors in `TensorLibrary.hep_lib_atom()` and execute through spenso tensor
networks. For scalar symbolic evaluators, lower through Symbolica and use its
JIT/evaluator facilities instead of Python loops.

For interleaved tensor-network execution, keep the `hep_lib`/tensor library fixed
throughout the evaluation. Do not register generated intermediate currents as new
library tensors. Build an aggregate `TensorNetwork`, multiply in the next current
or kernel, call `execute()` at the configured interleaving points, then carry the
executed residual network forward. When a residual must participate in later
sums/products, materialize it with `result_tensor(library)` and wrap it as a
local `TensorNetwork` value; this avoids scaled-term sum ambiguities without
mutating the library.

Momentum-dependent Feynman rules still belong to the model. Implement them as
model-owned spenso/Symbolica expressions using metric tensors and named
parametric rank-one momentum tensors supplied by lowering; do not move their
component formula into the recursion graph builder.

Momentum-dependent propagators follow the same boundary: the model owns the
symbolic tensor data/formula, while lowering only decides which graph currents
need a propagator, assigns stable graph-specific tensor names, and binds those
tensors to `ParamBuilder` current-momentum parameters.

When registering dense spenso tensors, inspect or test the stored tensor
structure instead of assuming Python row-major order from the expression call.
spenso may canonicalize slot order by representation while preserving the
expression-level call order. Every nontrivial model-owned tensor should have a
focused contraction regression against the intended AmpliCol current formula.
For the auxiliary gluon tensors used in the four-gluon route, spenso stores the
auxiliary representation before Minkowski slots and metric signs appear on
contracted Minkowski inputs. Keep direct tests against `_two_gluon_to_tensor`,
`_tensor_gluon_to_gluon`, and `_gluon_tensor_to_gluon` when editing these dense
tables.

For full matrix-element tensor-network lowering, apply color simplification to
the raw symbolic tensor expression before constructing/executing the spenso
`TensorNetwork`. Do not execute a network containing explicit color tensors and
then call `simplify_color` afterward, because execution may concretize color
structure too early.

Never expand the expression passed into `TensorNetwork`. The parser input must
remain the factorized raw current/tensor expression produced by recursion,
optionally after pre-network `simplify_color` on that raw expression. Do not use
`network.result_scalar()`, a scalar string preview, or any post-execution
component expression as input to another tensor network. Treat post-execution
scalar stringification as diagnostic-only and keep it opt-in/bounded.

Maintain a fully numerical tensor-network validation mode alongside the
parametric evaluator path. In that mode, register external currents, current
momenta, propagators, and other point-dependent tensors as concrete complex or
real dense tensors in `hep_lib`, then execute the same factorized tensor-network
input directly. This mode is a correctness probe for the parsed network and may
be slow; do not treat it as the production performance path unless measurements
show otherwise.

Do not run a separate Symbolica simplification/optimization pass on
`network.result_scalar()` as a default generation step. The executed scalar can
be large and highly factorized; any intentional simplification should happen on
the raw pre-execution expression, and additional evaluator-level optimization is
handled by Symbolica when building the evaluator.
In particular, do not add `scalar.optimize()` after `network.result_scalar()`
in the tensor-network generation path.

External momenta, external spinors, polarizations, and other point-dependent
rank-one inputs must be represented as parametric tensors in the spenso library.
Their concrete numerical values are filled only through pyamplicol's
`ParamBuilder`, which owns input ordering, contiguous ranges, phase metadata,
and evaluator save/load metadata. Evaluator wrappers should bundle the
Symbolica/spenso evaluator with its `ParamBuilder`; do not scatter ad hoc input
array assembly across generation or CLI code.

Symbolica scalar evaluators receive real-valued input arrays by default through
`Evaluator.evaluate`. Complex spinor and polarization components must use an
explicit complex-input evaluator wrapper that calls `Evaluator.evaluate_complex`
with `ParamBuilder`'s ordered complex values. If a future compiled-real backend
requires real-only inputs, use an explicit real/imaginary parameter packing in
`ParamBuilder` and reconstruct the complex tensor entries symbolically; do not
silently pass Python `complex` values to a real evaluator.

Validation must compare the same quantity as the legacy ME check:
leading-colour helicity/colour summed squared matrix element with the same
coupling normalization, identical-particle factors, initial averages, crossing
signs, massive-Z polarization convention, auxiliary tensor normalization, and
LC colour factor conventions.

Log and cache enough metadata to make failures debuggable: subprocess, colour
order, phase-space group, graph sizes, current counts, amplitude counts,
expression sizes, cache hits, generation time, tensor-network reduction time,
Symbolica optimization/JIT time, artifact load time, and per-point runtime.

## Dependency Patch Policy

The local dependency checkouts are part of the engineering surface for this
project. If Symbolica, spenso, or idenso lacks a primitive, representation hook,
simplification behavior, tensor-library feature, or evaluator interface that
would otherwise force brittle pyamplicol-side workarounds, patch the dependency
directly under `dependencies/symbolica`,
`dependencies/gammaloop/crates/spenso`,
`dependencies/gammaloop/crates/idenso`, or the corresponding Python bindings.

Patch dependencies only when it materially reduces complexity or improves the
native Rust-backed implementation path. Before patching, inspect the relevant
Rust/Python stubs to confirm the feature is missing or insufficient. After
patching, rebuild through `dependencies/install_dependencies.py` or the narrow
dependency build command, and add a focused pyamplicol regression test proving
why the dependency change is needed.

Avoid local Python workaround layers for operations that belong in spenso,
idenso, or Symbolica. Prefer one clean dependency patch over accumulating
translation code, ad hoc simplifiers, or manually expanded tensor/color rules in
pyamplicol.

When requested with `--with-gammaloop`, the installer builds the GammaLoop API
against the local Symbolica checkout with Symbolica's `gmp` feature enabled. GMP
is an accepted dependency for this project.

## Native Symbolica First

Performance matters. Symbolica is Rust-backed; handwritten Python symbolic
algorithms, tree walkers, replacement loops, simplifiers, derivative engines,
polynomial routines, tensor contractions, or integral reducers are forbidden
unless the local Symbolica/idenso/spenso/vakint/GammaLoop APIs have first been
checked and found insufficient for the specific operation.

Before adding or modifying symbolic code, explicitly inspect the Python stubs
and source listed below. Prefer native primitives even when a Python loop seems
small. If Python orchestration remains necessary, keep it at the boundary and
push the actual symbolic operation into Symbolica primitives.

When code starts checking atom types with `Expression.get_type()` or unpacking
children manually, stop and first try to express the operation with Symbolica
patterns: `Expression.match`, `Expression.matches`, `Expression.replace`,
`Expression.replace_multiple`, `Expression.replace_wildcards`, and
`Replacement`. Atom-type dispatch is acceptable only for narrow boundary code
such as external parsers or numeric coercions, or after a pattern-based attempt
has been shown not to express the required semantics.

Treat `Expression.match(pattern, restriction)` as the native way to collect
matching subexpressions: it yields every match from the expression, and the
matched subexpression can be reconstructed with `pattern.replace_wildcards`.
Do not add separate Python tree-walk collectors when a pattern match provides
the same data.

## API Discovery

The local source checkouts under `dependencies/` are the primary API reference.
Inspect implementation details directly in:

- `dependencies/symbolica/src/`
- `dependencies/gammaloop/crates/idenso/`
- `dependencies/gammaloop/crates/spenso/`
- `dependencies/gammaloop/crates/spynso3/`
- `dependencies/gammaloop/crates/vakint/`
- `dependencies/gammaloop/crates/gammaloop-api/`

For Python-facing APIs, also inspect the generated/source stub files:

- `dependencies/symbolica-community/python/symbolica/core.pyi`
- `dependencies/symbolica-community/python/symbolica/community/idenso/__init__.pyi`
- `dependencies/symbolica-community/python/symbolica/community/spenso/__init__.pyi`
- `dependencies/symbolica-community/python/symbolica/community/vakint/__init__.pyi`
- `dependencies/symbolica-community/python/symbolica/community/example_extension/__init__.pyi`
- `dependencies/symbolica/symbolica.pyi`

For GammaLoop's Python API, inspect:

- `dependencies/gammaloop/crates/gammaloop-api/python/gammaloop/__init__.py`
- `dependencies/gammaloop/crates/gammaloop-api/src/python.rs`

## Symbolica Expression Policy

Symbolica expressions are pyAmpliCol's canonical physics representation. Encode
fields, couplings, indices, derivatives, conjugation, group tensors, and EFT
bookkeeping directly with stable Symbolica function heads whenever practical.
Use Python objects only for registries, validated metadata, orchestration, and
services that cannot reasonably live in an expression.

All reusable Symbolica symbols, function heads, and pattern wildcards must be
created once in pyAmpliCol's central symbol registry and referenced from there.
Do not scatter calls such as `S("name")` through the codebase, and do not use
`E("...")` to construct reusable internal symbols, function heads, or pattern
wildcards. String parsing with `E("...")` or `Expression.parse("...")` is fine
for numeric coefficients and genuinely one-off literals such as `1/24`; the
centralization rule is about reusable pyAmpliCol symbols and expression heads, not
about every rational number.

Prefer expression construction such as:

```python
s.phi(s.flavor(s.quark), S("b"))
```

over:

```python
E("phi(flavor(quark), b)")
```

Pattern placeholders follow the same rule and belong in the central registry.
