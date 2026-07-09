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

The generated process directory is self-contained.  It includes a
`process_manifest.json`, serialized evaluator artifacts under `evaluators/`,
validation momenta, and a standalone checker:

```sh
cd outputs/dd_z_4g
python3 check_standalone.py --precision 16 --profile
```

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
30 GB:

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
`docs/performance_summary.md` and the older comparison notes are retained as
historical benchmark snapshots, while `docs/process_coverage.md` summarizes the
current parser/model/runtime coverage.
