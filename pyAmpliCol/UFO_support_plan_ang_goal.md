# pyAmpliCol Generic UFO Model Support

## Goal Statement

Implement model-general pyAmpliCol support on `pyamplicol_ufo_support`: retain the hard-coded Standard Model as `built-in-sm`, accept UFO directories, `ufo-model-loader` JSON models, and exact-version `*.pyAmplicol-model.json` artifacts, compile them through generic idenso/spenso tensor, color, chirality, propagator, and higher-point-vertex lowering, expose TOML run cards and JSON runtime parameters, validate UFO-SM numerics and production-JIT performance against the built-in SM and AmpliCol, then validate the scalar and scalar-gravity models, update all tests/examples/documentation/PDF matrices, and commit and push both the loader and pyAmpliCol changes.

## Execution Autonomy

- Never escalate or ask the user for permission unless completing a required outcome is genuinely impossible with the available permitted tools.
- Anticipate commands likely to be blocked and use non-destructive, policy-compliant alternatives from the start.
- If a command is blocked, try equivalent approaches such as narrower commands, repository-local moves instead of deletion, language/library APIs, or other permitted tools.
- Do not pause the overall workflow merely because one approach is blocked; continue independent work while finding an alternative.
- Ask for user intervention only after reasonable compliant alternatives have been exhausted and the missing operation is essential.
- Never remove files outside the repository or perform destructive Git operations without explicit authorization.

## 1. Baseline And Upstream Loader

- Record the current branch, dependency revisions, tests, built-in matrix caches, PDF, and representative built-in generation/runtime measurements before refactoring.
- Update the separate `ufo-model-loader` repository on `main`, commit there, and push directly to `origin/main`.
- Bump `ufo-model-loader` from `0.1.6` to `0.1.7`.
- Replace the removed Symbolica `Expression.evaluate_complex(...)` usage with the current Symbolica evaluation API, covering real and complex expressions.
- Fix automatic JSON restriction discovery so `restrict_default.json` and named `restrict_<name>.json` files work consistently.
- Extend loader serialization to preserve `propagating`, `goldstoneboson`, custom propagators, declared function names/arities, and enough form-factor metadata for explicit rejection.
- Continue loading old loader JSON files by supplying backward-compatible defaults for newly introduced fields.
- Read model-supplied custom propagators instead of always synthesizing defaults.
- Add loader tests for UFO/JSON parity, restrictions, complex expressions, metadata round-trips, custom propagators, and all three supplied models.
- Push the loader commit, capture its immutable SHA, and pin it in both pyAmpliCol's `pyproject.toml` and `install_dependencies.py`.
- Make the dependency installer record and smoke-test the pinned loader revision.
- Use published `symbolica==2.1.0` and its `symbolica.community.idenso` and `symbolica.community.spenso` APIs normally. GammaLoop remains reference material unless a concrete dependency defect requires a documented local patch.

## 2. Model Sources And Public CLI

- Keep the hard-coded implementation, expose it canonically as `built-in-sm`, and accept `BUILTIN_SM` and case-insensitive `built-in-sm`.
- Preserve the old Python class name as a compatibility alias; existing process artifacts still require regeneration.
- Detect UFO directories, loader JSON files, compiled pyAmpliCol model files, and the built-in model.
- Add `--model`, `--model-restriction default|none|NAME`, and `--model-simplify/--no-model-simplify` to model-consuming commands.
- Apply `restrict_default` and simplification by default.
- Reject restriction options for already compiled pyAmpliCol models.
- Add `inspect-model` with human and JSON reports covering contents, capabilities, warnings, and all incompatibilities.
- Add `compile-model --model SOURCE --output NAME.pyAmplicol-model.json`.
- Add repeatable `--multiparticle NAME=ITEM,ITEM,...`; ship explicit SM definitions for `p` and `j` without guessing aliases for arbitrary models.
- Source coupling-order definitions from the loaded model while preserving existing order-pruning controls.
- Report model loading, restriction, preflight, tensor lowering, color projection, Weyl projection, contact decomposition, serialization, DAG construction, warmup, evaluator construction, and JIT as separate monitored phases.

## 3. Run Cards And Parameter Cards

- Add `--run-card PATH` to every public subcommand.
- Use an `[arguments]` TOML table containing argparse destination names, including positional inputs.
- Keep the subcommand on the command line; the card does not select a command.
- Resolve paths relative to the run-card directory.
- Apply precedence as parser defaults, card values, then explicitly supplied CLI values.
- Validate card values with argparse's type and choice rules.
- Treat unknown keys, cross-command keys, malformed values, and missing required inputs as hard errors.
- Replace runtime model-parameter TOML with loader-compatible complex JSON:

```json
{
  "aS": [0.118, 0.0],
  "MT": [173.0, 0.0]
}
```

- Use this format for built-in and external models.
- Emit a complete `model-parameters.json` with every process output.
- Permit partial runtime override cards but reject unknown or compile-time-eliminated parameters.
- Keep every surviving external parameter runtime-adjustable. Parameters removed through restriction-based simplification require recompiling an unrestricted model.

## 4. Native Compiled Model

- Introduce an immutable `CompiledModel` IR containing particle resolution, parameters, couplings, orders, wavefunctions, propagators, tensors, oriented kernels, synthetic auxiliaries, color projections, normalization, capabilities, and provenance.
- Serialize it as `*.pyAmplicol-model.json`.
- Treat it as an exact-version artifact and reject mismatched pyAmpliCol or dependency fingerprints with a regeneration instruction.
- Store expressions and tensor data rather than machine code so it remains platform-portable within the matching dependency fingerprint.
- Always include the exact compiled model in every generated process output.
- Maintain a content-addressed model cache keyed by source contents, restrictions, simplification, function registry, compiler version, and dependency revisions.
- Use atomic cache writes and make `compile-model` expose the same conversion explicitly.
- Separate model-conversion timing from process-generation timing. Performance parity starts from the compiled pyAmpliCol model.
- Generate the built-in SM through this same IR while preserving its hard-coded source definitions.

## 5. Generic UFO Lowering

- Import raw UFO modules in an isolated worker process to avoid Python module/path contamination. Document that UFO modules are executable code and must be trusted.
- Preflight the complete model before expensive lowering and report all incompatibilities together.
- Support physical spins 0, 1/2, 1, and 2.
- Support Dirac fermion flow; reject Majorana/FNV interactions, spin-3/2, and higher spins.
- Support process generation for color representations `1`, `±3`, and `8`.
- Preserve sextet metadata but reject processes involving `±6` in this milestone.
- Reject counterterm/loop vertices, opaque callbacks, form factors, and unknown UFO functions.
- Implement a strict name-and-arity registry for standard complex, power, exponential, logarithmic, trigonometric, inverse, hyperbolic, and reciprocal-trigonometric functions.
- Never execute model-supplied function bodies.
- Normalize UFO indices into typed spenso indices without relying on model-writer naming conventions.
- Map standard UFO momentum, metric, gamma, projector, sigma, and color tensors into spenso representations.
- Use `simplify_metrics`, `simplify_gamma`, and `simplify_color`, followed by spenso materialization.
- Build every color x Lorentz x coupling term independently and orient it for each output leg.
- Convert four-component Dirac structures into the existing chiral/Weyl current basis with explicit projectors.
- Derive full color exactly and LC/NLC through the current large-`N_c` definitions applied to generic simplified expressions.
- Respect source gauge and propagator metadata without silent gauge conversion.
- Forbid ghosts and Goldstones as ordinary external states; retain propagating Goldstones internally when required.
- Preserve explicit custom propagators unchanged.
- Use fixed-width denominators `p^2-m^2+i m Gamma` for synthesized propagators.
- Implement scalar, Dirac/Weyl, vector, massless spin-2, and massive spin-2 wavefunctions and propagators.
- Build spin-2 helicities from standard vector-helicity tensor products and Clebsch-Gordan combinations.
- Use the canonical Fierz-Pauli massive projector when no custom propagator exists.
- Mark massive spin-2 support experimental until end-to-end validation with an external model is available.

## 6. Higher-Point Vertices And DAG Integration

- Support arbitrary UFO vertex valence, including ten-point scalar and five-point scalar-gravity contacts.
- Lower each color/Lorentz/coupling term into a deterministic balanced sequence of trivalent operations.
- Represent partial contractions as synthetic non-propagating auxiliary currents carrying unresolved indices.
- Insert each original coupling exactly once and introduce no fake propagator.
- Track original-vertex/term provenance so fragments from different contacts cannot mix.
- Canonicalize identical-particle permutations without duplicating diagrams or losing UFO normalization.
- Use spenso contraction-cost estimates to select balanced contraction order.
- Keep an independent eager n-ary contraction path as a correctness oracle.
- Refactor process parsing, particle classes, source states, closures, DAG construction, and artifacts to use `CompiledModel` rather than SM-global maps.
- Preserve selected-flow/all-flow APIs, coupling pruning, shared LC-order recycling, stage-local parameter layouts, chunking, and current-filter safety.
- Update Rusticol to read complex JSON parameter cards and construct real/imaginary runtime slots.
- Remove the runtime TOML parameter parser.
- Make normalization model-driven so UFO couplings do not receive built-in global coupling factors.
- Detect old process artifacts and instruct users to regenerate them.

## 7. Assets, Examples, And Tests

- Bundle UFO and loader-JSON copies of `sm`, `scalars`, and `scalar_gravity` under `assets/models`.
- Include these resources in wheels and verify repository/package copies with a hash manifest.
- Add example cards for built-in SM, UFO SM, loader JSON, compiled models, runtime parameters, multiparticle aliases, each color mode, scalars, and scalar gravity.
- Add unit tests for model detection, run-card precedence, paths, compiled-model fingerprints, cache invalidation, parameters, preflight, functions, propagators, Weyl/color projection, and contact decomposition.
- Add formula-level massive-spin-2 tests for symmetry, transversality, tracelessness, orthonormality, and completeness at arbitrary precision.
- Add fast integration tests for representative cards and low-multiplicity JIT generation.
- Mark complete multiplicity ladders and performance matrices as slow tests.
- Verify an installed wheel can locate bundled models and run `inspect-model`, `compile-model`, and a low-multiplicity process.

## 8. Numerical And Performance Validation

- Validate UFO and JSON SM against built-in SM and AmpliCol at identical phase-space points, parameters, widths, helicities, color orderings, and normalization.
- Align complex roots by one process-global phase determined from the largest stable root, then compare all mapped roots.
- Compare final reduced observables directly without phase adjustment.
- Let `Q=max(sqrt(abs(s)),1 GeV)` and `d=4-N_ext`. Require:

  ```text
  abs(A-B) / max(abs(A), abs(B), 1e-30*Q^d) <= 1e-10
  ```

  Use `2d` for squared matrix elements.
- Preserve separate exact-zero and current-pruning validation.
- Re-evaluate disputed points with identical double-precision kinematics upcast to at least 50 decimal digits.
- Begin with the existing 94-case LC/NLC/full fixture through `n<=4`.
- Advance multiplicity only when every applicable process class at the lower multiplicity is green.
- Reproduce every currently valid populated built-in matrix cell with UFO-SM.
- Use documented matrix SymJIT O1 selected-flow and all-flow builds as the blocking production-JIT baseline.
- Require process generation from compiled models to remain within 20% of matched built-in results, and runtime to remain within 10%.
- Keep O3, ASM, and C++ comparisons informational.
- Use median-of-three generation measurements below 90 seconds and one isolated matched run for slower builds.
- Use `--target-runtime 10` with repeated runtime samples and report timing spread.
- Run heavy work behind the 30 GB process-tree watchdog.
- Preserve generated process outputs for later timing and inspection.
- Treat external-SM validation through `n<=4` as the first delivery milestone: require every applicable same-point numerical comparison to be green, publish the measured generation/runtime values as informational results, refresh the external-SM PDF tables, then commit and push that baseline before advancing. Performance acceptance is deliberately deferred for this first baseline.
- After the numerical-only `n<=4` baseline, stop progression and investigate whenever the numerical gate, 20% generation gate, or 10% runtime gate fails.
- Extend the external-SM ladder through `n<=6` as the second delivery milestone, again advancing multiplicity only after every lower applicable class is green; refresh the PDF, commit, and push once the complete `n<=6` milestone passes.
- Before attempting any external-SM process above `n=6`, validate increasing-multiplicity scalar and scalar-gravity processes, including a case exercising the ten-point contact. Keep each documented generation attempt to five minutes or less, and build concise model-specific result matrices rather than mirroring the full SM campaign.
- Do not rerun broad built-in or external-SM matrix multiplicities above `n=6` until the final external-SM campaign. Correct blocking implementation or validation defects immediately, but prefer targeted probes over benchmark refreshes during the scalar/scalar-gravity milestone.
- Confirm UFO, loader JSON, compiled model, eager n-ary, and decomposed DAG agreement.
- Confirm LC/NLC/full equality for colorless models.
- Record UFO, loader-JSON, and compiled-model conversion/loading costs separately.
- Commit and push the scalar/scalar-gravity implementation, validation, matrices, and refreshed PDF once that milestone is green; only then resume and finish the external-SM validation above `n=6`.

## 9. Result Matrices And Documentation

- Archive the current PDF, TeX tables, and JSON caches before adding UFO-derived results.
- Retain complete built-in-SM matrices and add complete UFO-SM matrices for every currently populated LC/NLC/full and Z-family cell.
- Populate and publish the external-SM matrices in staged `n<=4`, then `n<=6`, then `n>6` milestones; the scalar and scalar-gravity matrices must be completed between the `n<=6` and `n>6` SM stages.
- Keep C++ matrix comparisons only through `n<=3`; invalidate or hide stale C++ cells above `n=3` rather than regenerating them.
- Keep both built-in and external-SM dedicated Z-family tables complete and current through at least `n=6` throughout the remaining UFO work.
- Add separate concise scalar and scalar-gravity process-result matrices to the PDF, using representative increasing-multiplicity processes and a five-minute maximum generation time per documented case.
- Use separate generated JSON/TeX files per model source and color mode.
- Regenerate results in increasing multiplicity order while preserving selected-flow/helicity-summed and all-flow/fixed-helicity semantics.
- Use matching AmpliCol workloads and never compare different color/helicity reductions.
- Refresh `pyAmpliCol.pdf` after every completed cell or small group.
- Preserve every generated process directory.
- Document model architecture, trust boundaries, restrictions, compiled models, parameter JSON, supported features, contact decomposition, tensor/color/Weyl projection, spin-2 conventions, validation, and conversion costs.
- Update README and Markdown files with installation, loader revision, run cards, CLI examples, parameter format, regeneration requirements, and dependency patches.
- Cite the [UFO 2.0 specification](https://arxiv.org/abs/2304.09883), [`ufo-model-loader`](https://github.com/alphal00p/ufo_model_loader), and Symbolica Community APIs.
- Render and visually inspect the final PDF for widths, page breaks, listings, captions, and stale built-in-only or TOML text.

## 10. Final Acceptance And Delivery

- Run the complete loader test suite against current Symbolica before committing and pushing it.
- Verify dependency installation from a clean managed environment.
- Run pyAmpliCol unit, fast integration, selected slow integration, typing, wheel-install, and result-matrix audit tests.
- Run Rusticol checks/tests for JSON parameters and generic model metadata.
- Run `git diff --check`, compile Python modules, rebuild the PDF, and verify no generation processes remain running.
- Commit and push `ufo-model-loader` 0.1.7 to `origin/main`.
- Commit implementation, tests, assets, examples, generated caches/tables/PDF, outputs, and documentation on `pyamplicol_ufo_support`.
- Commit and push `pyamplicol_ufo_support` after each validated milestone: external SM through `n<=4`; external SM through `n<=6`; scalar and scalar-gravity support and matrices; and the final external-SM `n>6` campaign.
- Report loader and pyAmpliCol commit SHAs, validation coverage, performance gates, unsupported features, dependency patches, and experimental massive-spin-2 limitations.
