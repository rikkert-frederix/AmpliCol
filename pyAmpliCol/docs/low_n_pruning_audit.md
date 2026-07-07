# Low-Multiplicity LC Pruning Audit

This audit checks whether the generic DAG compiler is missing obvious AmpliCol-style current recycling or helicity pruning for representative low-multiplicity processes.  All Fortran numbers below come from the generated-library workflow (`--library=create`, `make amplicol_generate_library`) and inspect the first generated library module for the same selected leading-colour sector used in the result matrix.  pyAmpliCol counts are the selected-sector generic DAG counts after LC sector filtering, helicity pruning, dead-tree pruning, and the default 10-point numerical current warmup passes.

| process | py currents | py interactions | py roots | py helicity weights | Fortran val_c | Fortran int_c | Fortran amps | status |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| `d d~ > z g g g g` | 333 | 1202 | 96 | `1.0` | 330 | 1192 | 94 | aligned within small convention differences |
| `g g > g g g g g` | 308 | 1561 | 56 | `2.0` | 310 | 1665 | 55 | zero-current filter prunes two currents; validation clean |
| `d d~ > u u~ s s~ g` | 51 | 79 | 8 | `2.0` | 73 | 142 | 8 | pyAmpliCol is more pruned and validates |
| `g g > t t~ g g g` | 310 | 1232 | 128 | `1.0` | 312 | 1234 | 128 | zero-current filter removes four chirality-zero currents; validation clean |
| `d d~ > z z g g g` | 518 | 1856 | 144 | `1.0` | 517 | 1856 | 144 | interaction/root recycling matches Fortran |

A temporary `n <= 5` result-matrix validation refresh was run after enabling the numerical zero-current filter and same/opposite-sign identical-current detection: 45 applicable cells validate against Fortran AmpliCol, with maximum relative difference about `3.1e-13` and no validation failures.  The refresh was run behind the 30 GB watchdog and peaked at about 0.50 GB RSS.  The committed performance data file was left untouched because this was a correctness check with a short runtime target, not a performance-table refresh.  The pure-QCD massless rows use global helicity-flip root grouping through explicit `helicity_weight=2.0`; massive-top and EW rows keep unit helicity weights because that pruning is not safe there.

The generation-time numerical current passes are intentionally process-generic.  They evaluate the retained DAG on 10 deterministic RAMBO phase-space points.  The zero-current filter declares only internal currents whose components stay below `max(1e-300, 1e-12 * global_current_max)` to be zero, removes those currents, and runs the reachability prune again.  It can be disabled with `--no-numerical-filter-current` for diagnostics.  The identical-current detector then compares sampled current signatures with a mixed absolute/relative tolerance `max(threshold, 1e-12 * max(|a|, |b|))`, supports both same-sign and opposite-sign equality, and can be disabled with `--no-numerical-current-merging`.  The merge identity includes the AmpliCol-style helicity ancestry; a broader ancestry-blind merge was tested and rejected because it changes massive-top matrix elements at percent level.

With the safe helicity-ancestry key, identical-current merging did not find any mergeable currents in the representative `n <= 5` validation survey.  That is still the desired conservative behavior: the implementation now supports opposite-sign current detection, but it only materializes a merge when the numerical evidence and current identity are both safe.

The detailed count comparison is reproducible with:

```bash
python3 pyAmpliCol/docs/low_n_pruning_audit.py \
  --manifest pyAmpliCol/docs/.result_matrix_outputs/gg_tt_jets/n5/jit/generic_process_manifest.json \
  --fortran-module Library/amp1_1_lib.f03
```

For the first generated Fortran colour-order module of `g g > t t~ g g g`, the post-filter stage/kind delta is:

| stage size | vertex kind | pyAmpliCol | Fortran | delta |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0 | 14 | 14 | +0 |
| 2 | 1 | 14 | 16 | -2 |
| 2 | 6 | 4 | 4 | +0 |
| 3 | 0 | 44 | 44 | +0 |
| 3 | 1 | 48 | 44 | +4 |
| 3 | 2 | 20 | 24 | -4 |
| 3 | 3 | 24 | 24 | +0 |
| 3 | 6 | 12 | 16 | -4 |
| 4 | 0 | 88 | 88 | +0 |
| 4 | 1 | 88 | 88 | +0 |
| 4 | 2 | 56 | 56 | +0 |
| 4 | 3 | 56 | 64 | -8 |
| 4 | 6 | 48 | 48 | +0 |
| 5 | 0 | 112 | 112 | +0 |
| 5 | 2 | 80 | 80 | +0 |
| 5 | 3 | 80 | 96 | -16 |
| 5 | 6 | 128 | 128 | +0 |
| 6 | 6 | 320 | 288 | +32 |

The remaining small deltas are interaction-materialization details rather than missing process-family pruning.  The same filtered DAG validates numerically against Fortran AmpliCol in the low-n matrix.
