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

The only visible low-multiplicity pruning gap is `g g > t t~ g g g`, where pyAmpliCol has the same number of amplitude roots as Fortran but about 8% more interaction entries.  A quick Fortran-library parse shows the difference is interaction-materialization level, not amplitude-root recycling.  This does not indicate a process-family-specific DAG issue; it is a smaller optimization target for the generic QCD kernels if top-pair runtime becomes a priority.
