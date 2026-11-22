#!/bin/sh
# Reject invalid complex-wave control, scheme, and operator inputs.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/complex_wave"

if [ ! -x "$bin" ]; then
    echo "error: build complex_wave first (bin/$precision/complex_wave)" >&2
    exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

expect_fail() {
    name=$1
    if ( cd "$tmp" && "$bin" >"$tmp/$name.log" 2>&1 ); then
        echo "expected failure: $name" >&2
        cat "$tmp/$name.log" >&2
        exit 1
    fi
    echo "rejected: $name"
}

# Unknown interpolation scheme.
cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/scheme = .*/scheme = 'NOT_A_SCHEME',/" "$tmp/input.nml"
expect_fail unknown_scheme

# Unknown difference operator.
cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/difference_operator = .*/difference_operator = 'NOPE',/" "$tmp/input.nml"
expect_fail unknown_operator

# Period count must be a positive integer.
cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/periods = .*/periods = 0,/" "$tmp/input.nml"
expect_fail periods_zero

cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/periods = .*/periods = -1,/" "$tmp/input.nml"
expect_fail periods_negative

# CFL must be strictly positive.
cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/cfl = .*/cfl = 0.0d0,/" "$tmp/input.nml"
expect_fail cfl_zero

cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/cfl = .*/cfl = -0.1d0,/" "$tmp/input.nml"
expect_fail cfl_negative

# The initial condition needs at least 20 points.
cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/mesh_num = .*/mesh_num = 19,/" "$tmp/input.nml"
expect_fail mesh_below_20

cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/mesh_num = .*/mesh_num = 0,/" "$tmp/input.nml"
expect_fail mesh_zero

echo "complex_wave params: PASS"
