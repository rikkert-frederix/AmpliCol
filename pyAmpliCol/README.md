# pyAmpliCol

Python interface and tooling foundation for AmpliCol.

Install the managed native dependencies from this directory with:

```sh
python dependencies/install_dependencies.py
```

The GammaLoop Python API is optional and is not installed by default. Request it
explicitly with:

```sh
python dependencies/install_dependencies.py --with-gammaloop
```

After dependency setup, activate the managed environment with:

```sh
source dependencies/.venv/bin/activate
```

Useful development commands:

```sh
pyamplicol processes 'd d~ > z g g' --json
pyamplicol processes 'd d~ > z g g' --legacy-output processes.txt
pyamplicol generate 'd d~ > z g g' --json
pyamplicol evaluate 'd d~ > z' --json
pyamplicol evaluate 'd d~ > z g' --sqrt-s 1000 --json
pyamplicol evaluate 'd d~ > z g g' --sqrt-s 1000 --json
pyamplicol evaluate 'd d~ > z g g g g g g' --sqrt-s 1000 --json
pyamplicol evaluate 'd d~ > z g' --sqrt-s 1000 --runtime-backend numeric-tensor-network --json
pyamplicol compare-amplicol 'd d~ > z g g g g g g' --amplicol-probe --points 10
pyamplicol validate-z-gluon-family --max-gluons 6 --points 10 --runtime-backend dag
```

Current native matrix-element status:

- `q q~ -> Z` is evaluated by a native AmpliCol-convention numerical kernel
  using the ported external wavefunctions, Z couplings, leading-colour factor,
  and initial-state averaging. Its generated artifact also carries a serialized
  Symbolica scalar evaluator for the raw helicity sum; `evaluate` reloads that
  artifact when present and reports the Symbolica/native relative difference.
- `q q~ -> Z + n g` for `n = 1..6` is evaluated by the default
  `native-spenso-symbolica-current-dag` runtime. It reuses the generated
  AmpliCol recursion graph, but delegates vertex, propagator, and tensor
  building blocks to spenso/Symbolica kernels. A staged Python recursion remains
  available as `--runtime-backend python` for reference and debugging. The
  gluon-current recursion includes the color-ordered three-gluon subcurrent,
  the auxiliary antisymmetric tensor route for four-gluon contributions, and
  final identical-gluon normalization.
- Direct AmpliCol probe validation passes for `d d~ -> Z + n g`, `n = 0..6`,
  with 10 deterministic probe points per multiplicity through:

  ```sh
  pyamplicol validate-z-gluon-family --max-gluons 6 --points 10 --runtime-backend dag
  ```

  The latest watchdog run validated all seven multiplicities with maximum
  relative difference `2.027596339685497e-12` after the legacy probe output was
  upgraded to double precision.
- `compare-amplicol --amplicol-probe` reports both legacy AmpliCol command and
  timing-table data, pyamplicol generation/artifact-cache timings, and
  pyamplicol native re-evaluation runtime metrics for the same probe points.
- `validate-z-gluon-family` runs the same direct-probe comparison across the
  milestone family and summarizes generation, legacy command time, per-point
  pyamplicol runtime, and maximum relative difference per multiplicity.
- `q q~ -> Z` is covered by the native AmpliCol-convention numerical kernel,
  unit regression, and the fixed on-shell `--amplicol_fixed_probe` legacy
  adapter path. This avoids the old random 2-to-1 phase-space failure mode.
- Supported `generate` calls now include a `symbolic_lowering` report that
  exercises real spenso/idenso hooks: model-owned auxiliary tensor registration
  and contraction through `TensorLibrary.hep_lib_atom()` plus an idenso SU(N)
  color simplification probe. The same report also records a graph-derived
  Symbolica recursion plan with current, interaction, amplitude, color-order,
  auxiliary-tensor-route, expression-size, and bounded expression-preview
  metadata, plus a model-owned per-vertex lowering map identifying which
  recursion vertices are already spenso-backed and which remain pending full
  tensor-network lowering.
- The lowering report also carries a graph-derived spenso tensor-network
  blueprint. It recursively turns currents into indexed tensor expressions,
  uses registered spenso tensors for the auxiliary `two_gluon_to_tensor`,
  `tensor_gluon_to_gluon`, and `gluon_tensor_to_gluon` vertices, and now also
  uses model-owned chirality-specific Weyl-vector tensors for the quark-gluon
  and Z-current recursion vertices plus a model-owned spenso expression for
  the momentum-dependent color-ordered three-gluon current. Blueprints register
  reachable external source currents and the current momenta needed by
  three-gluon vertices and propagators as `ParamBuilder` rank-one tensors.
  Gluon and chirality-dependent Weyl-quark propagators are registered as
  graph-specific symbolic spenso tensors. Small propagated blueprints execute to
  fully propagated parameterized scalar expressions; larger propagated
  blueprints stay build-only behind a size guard to avoid unbounded expression
  expansion.
- `ParamBuilder` provides the evaluator input boundary for full lowering:
  point-dependent external rank-one tensors can be registered parametrically in
  spenso, bundled with serialized Symbolica evaluators, and reloaded with stable
  input ordering. It covers real-valued tensors such as external momenta through
  Symbolica's real evaluator path, and complex tensors such as spinors or
  polarizations through an explicit `evaluate_complex` bundle mode.
- `--runtime-backend numeric-tensor-network` is an explicit validation mode. It
  registers external currents, current momenta, and graph-specific propagators
  as fully numerical dense tensors in `TensorLibrary.hep_lib_atom()` and then
  executes the same factorized tensor-network input directly. This is slower
  than the current-DAG runtime, but it confirms that the parsed tensor network
  itself contracts to the expected matrix element.
- Cached `generate` calls write both metadata and a versioned evaluator
  manifest. For `q q~ -> Z`, the manifest kernel is `symbolica-zero-gluon` and
  includes a serialized Symbolica evaluator. For `q q~ -> Z + n g`, evaluator
  manifests are emitted for small tensor-network scalar artifacts where the size
  guard allows reduction; larger multiplicities still carry stable graph and
  lowering metadata for the current-DAG runtime and future cached JIT evaluators.
- `evaluate` reports whether such an evaluator manifest was loaded, including
  the zero-gluon Symbolica cross-check when available. `profile` reports
  generation, artifact-cache, artifact-load, Symbolica evaluator build/reduction
  timings when present, selected native runtime backend, and native per-point
  runtime timings.
- The remaining optimization frontier is replacing the current-DAG production
  runtime with a cached full-network/JIT strategy at high multiplicity without
  expanding the parsed tensor-network input or exceeding the watchdog memory
  budget.
