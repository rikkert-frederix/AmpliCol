# pyAmpliCol

pyAmpliCol is the Python generation and validation layer for the modern
AmpliCol matrix-element prototype.  The current production path generates a
shared-current eager-DAG process artifact in Python and evaluates it through the
Rusticol PyO3 runtime using serialized Symbolica evaluators.

The validated fast path currently targets leading-colour
`q q~ > Z + n gluons` processes.  The surrounding process/model architecture is
kept broader, but other electroweak final states still require additional
current lowering before they can be used with `generate-process`.

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
./pyamplicol.sh generate 'd d~ > Z g g' --json
./pyamplicol.sh evaluate 'd d~ > Z g g' --sqrt-s 1000 --json
./pyamplicol.sh compare-amplicol 'd d~ > Z g g g g' --amplicol-probe --points 10
./pyamplicol.sh validate-z-gluon-family --max-gluons 6 --points 10 --runtime-backend dag
```

Long validation and benchmark jobs should be run behind the bundled RAM
watchdog so the process is stopped before exceeding the repository guideline of
30 GB:

```sh
python3 ./scripts/run_with_memory_watch.py --max-rss-gb 30 -- \
  ./pyamplicol.sh validate-z-gluon-family --max-gluons 6 --points 10
```

## Documentation

See `docs/description.pdf` for the architecture overview and
`docs/performance_summary.md` for the current benchmark summary.
