#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
path_code="$repo_root/Utilities"
path_mg="${MG5_PATH:-/Users/vjhirsch/MG5/MG5_aMC_v3_6_0}"

if [ ! -x "$path_mg/bin/mg5_aMC" ]; then
    echo "MG5 executable not found: $path_mg/bin/mg5_aMC" >&2
    echo "Set MG5_PATH to a valid MG5_aMC installation." >&2
    exit 1
fi

version=$(awk '/^version/ {print $3}' "$path_mg/VERSION")
matrix_template="$path_mg/madgraph/iolibs/template_files/matrix_standalone_v4.inc"
color_algebra="$path_mg/madgraph/core/color_algebra.py"
check_sa="$path_mg/madgraph/iolibs/template_files/check_sa.f"

tmpdir=$(mktemp -d)
restore_mg_templates() {
    cp "$tmpdir/matrix_standalone_v4.inc" "$matrix_template"
    cp "$tmpdir/color_algebra.py" "$color_algebra"
    cp "$tmpdir/check_sa.f" "$check_sa"
    rm -rf "$tmpdir"
}
trap restore_mg_templates EXIT

cp "$matrix_template" "$tmpdir/matrix_standalone_v4.inc"
cp "$color_algebra" "$tmpdir/color_algebra.py"
cp "$check_sa" "$tmpdir/check_sa.f"

case "$version" in
    3.6.0|3.6.1)
        cp "$path_code/ME_checks/matrix_standalone_v4_v361.inc" "$matrix_template"
        ;;
    3.6.6)
        cp "$path_code/ME_checks/matrix_standalone_v4_v366.inc" "$matrix_template"
        ;;
    *)
        echo "Unsupported MG5 version: $version" >&2
        echo "Supported versions: 3.6.0, 3.6.1, 3.6.6" >&2
        exit 1
        ;;
esac

cp "$path_code/ME_checks/color_algebra.py" "$color_algebra"
cp "$path_code/ME_checks/check_sa.f" "$check_sa"

cd "$path_mg/bin"
"$path_code/ME_checks/run.sh" "$path_code" "$@"
