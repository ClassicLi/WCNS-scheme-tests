#!/bin/sh
# Check that a reference sidecar log matches the requested mesh/scheme/operator.
#
# Exit 0 when every field matches, 1 when the log is missing or differs.
set -eu

log=$1
mesh=$2
scheme=$3
operator=$4

[ -f "$log" ] || exit 1
grep -qx "mesh_num=$mesh" "$log" || exit 1
grep -qx "scheme=$scheme" "$log" || exit 1
grep -qx "operator=$operator" "$log" || exit 1
