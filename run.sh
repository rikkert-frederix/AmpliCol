#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 'process string' [process_list.py options]" >&2
  echo "Example: $0 'd d~ > z g g'" >&2
  exit 2
fi

rm -f processes.txt
python3 process_list.py "$@"
make cleanlib
make -j8 amplicol_generate
./amplicol_generate --library=create --process=processes.txt
make -j8 amplicol_generate_library
if [ "${ME_TEST:-0}" != "0" ]; then
  MG5_PATH="${MG5_PATH:-/Users/vjhirsch/MG5/MG5_aMC_v3_6_0}" \
    ./amplicol_generate --me_test="${ME_TEST}" --timing="${ME_TEST}" --process=processes.txt
fi
if [ "${NEVENTS:-0}" != "0" ]; then
  ./amplicol_generate --library=use --nevents="${NEVENTS}" --seed="${SEED:-101}"
fi
