#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
bin="$root/bin/$precision/interpolation_test"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cp "$root/input.nml" "$tmp/input.nml"
mkdir -p "$tmp/results/convergence"
sed -i "s/ini_mesh_num = .*/ini_mesh_num = 10,/; s/iter_times = .*/iter_times = 2,/" "$tmp/input.nml"

( cd "$tmp" && "$bin" )

test -f "$tmp/LINEAR3_M1P1_log.dat"
grep -q "L1_error" "$tmp/LINEAR3_M1P1_log.dat"
test -f "$tmp/u_exact.dat"
