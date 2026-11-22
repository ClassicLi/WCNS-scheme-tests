#!/bin/sh
# Run one ADR case in a temporary directory and copy its CSV/log.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
kind=${1:?usage: run_cases.sh <weno|wcns> <scheme> [operator] [dest]}
scheme=$2
operator=${3:-}
dest=${4:-"$root/results"}
mesh=${ADR_MESH:-827}

bin="$root/bin/$precision/adr_$kind"
if [ ! -x "$bin" ]; then
    echo "error: $bin not found; run make first" >&2
    exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

cp "$root/input.nml" "$work/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = $mesh,/" \
    -e "s/scheme = .*/scheme = '$scheme',/" \
    "$work/input.nml"

if [ "$kind" = "wcns" ]; then
    op=${operator:-HCS4}
    sed -i "s/difference_operator = .*/difference_operator = '$op',/" "$work/input.nml"
else
    op=""
fi

( cd "$work" && "$bin" ) >/dev/null

mkdir -p "$dest/csv" "$dest/logs"
cp "$work"/results/csv/*.csv "$dest/csv/"
cp "$work"/results/logs/*.log "$dest/logs/" 2>/dev/null || true
cp "$work"/results/samples/*.dat "$dest/samples/" 2>/dev/null || true

echo "done: $kind $scheme ${op}"
