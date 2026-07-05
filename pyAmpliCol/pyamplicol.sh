#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYAMPLICOL_PYTHON:-${SCRIPT_DIR}/dependencies/.venv/bin/python}"

if [[ ! -x "${PYTHON}" ]]; then
  cat >&2 <<EOF
error: pyAmpliCol Python environment not found at:
  ${PYTHON}

Run this first from ${SCRIPT_DIR}:
  python dependencies/install_dependencies.py
EOF
  exit 127
fi

export PYTHONPATH="${SCRIPT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"
exec "${PYTHON}" -m pyamplicol "$@"
