# Low-Multiplicity LC Pruning Audit

This audit checks whether the generic DAG compiler is missing obvious AmpliCol-style current recycling or helicity pruning for representative low-multiplicity processes.  All Fortran numbers below come from the generated-library workflow (`--library=create`, `make amplicol_generate_library`) and inspect the first generated library module for the same selected leading-colour sector used in the result matrix.  pyAmpliCol counts are the selected-sector generic DAG counts after LC sector filtering, helicity pruning, and dead-tree pruning.

| process | py currents | py interactions | py roots | py helicity weights | Fortran val_c | Fortran int_c | Fortran amps | status |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| `d d~ > z g g g g` | 333 | 1202 | 96 | `1.0` | 330 | 1192 | 94 | aligned within small convention differences |
| `g g > g g g g g` | 310 | 1665 | 56 | `2.0` | 310 | 1665 | 55 | current/interaction recycling matches Fortran |
| `d d~ > u u~ s s~ g` | 51 | 79 | 8 | `2.0` | 73 | 142 | 8 | pyAmpliCol is more pruned and validates |
| `g g > t t~ g g g` | 314 | 1332 | 128 | `1.0` | 312 | 1234 | 128 | mild interaction-level excess in pyAmpliCol |
| `d d~ > z z g g g` | 518 | 1856 | 144 | `1.0` | 517 | 1856 | 144 | interaction/root recycling matches Fortran |

The representative `n <= 5` validation rows in `result_matrix_data.json` are clean, including pure gluons, massive top-pair channels, EW vector channels, and multi-quark-line channels.  The pure-QCD massless rows use global helicity-flip root grouping through explicit `helicity_weight=2.0`; massive-top and EW rows keep unit helicity weights because that pruning is not safe there.

The only visible low-multiplicity pruning gap is `g g > t t~ g g g`, where pyAmpliCol has the same number of amplitude roots as Fortran but about 8% more interaction entries.  A Fortran-library parse shows the difference is interaction-materialization level, not amplitude-root recycling.  Rechecking the artifact pipeline confirms that this excess survives dead-tree pruning and the helicity-pruning pass; it is not caused by comparing raw compiler counts with pruned runtime counts.  This does not indicate a process-family-specific DAG issue; it is a smaller optimization target for the generic massive-QCD kernels if top-pair runtime becomes a priority.

The detailed count comparison is reproducible with:

```bash
python3 pyAmpliCol/docs/low_n_pruning_audit.py \
  --manifest pyAmpliCol/docs/.result_matrix_outputs/gg_tt_jets/n5/jit/generic_process_manifest.json \
  --fortran-module Library/amp1_1_lib.f03
```

For the first generated Fortran colour-order module, the stage/kind delta is:

| stage size | vertex kind | pyAmpliCol | Fortran | delta |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0 | 16 | 14 | +2 |
| 2 | 1 | 16 | 16 | +0 |
| 2 | 6 | 4 | 4 | +0 |
| 3 | 0 | 48 | 44 | +4 |
| 3 | 1 | 48 | 44 | +4 |
| 3 | 2 | 24 | 24 | +0 |
| 3 | 3 | 24 | 24 | +0 |
| 3 | 6 | 16 | 16 | +0 |
| 4 | 0 | 96 | 88 | +8 |
| 4 | 1 | 96 | 88 | +8 |
| 4 | 2 | 64 | 56 | +8 |
| 4 | 3 | 64 | 64 | +0 |
| 4 | 6 | 48 | 48 | +0 |
| 5 | 0 | 128 | 112 | +16 |
| 5 | 2 | 96 | 80 | +16 |
| 5 | 3 | 96 | 96 | +0 |
| 5 | 6 | 128 | 128 | +0 |
| 6 | 6 | 320 | 288 | +32 |

The pattern points to AmpliCol's generated-library numerical/local-current pruning removing a small set of vertex entries after the generic current topology is known.  pyAmpliCol currently applies structural reachability, colour-order, coupling-order, dead-tree, and safe global-helicity pruning, but it does not yet run an AmpliCol-style numerical zero/equivalence pass over retained interactions.  That pass should be implemented generically, using deterministic probe points and current/value signatures, before attempting process-specific top-pair tuning.
