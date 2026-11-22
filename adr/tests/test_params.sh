#!/bin/sh
# Reject invalid scheme, operator, mesh, and dt inputs.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
weno_bin="$root/bin/$precision/adr_weno"
wcns_bin="$root/bin/$precision/adr_wcns"

if [ ! -x "$weno_bin" ] || [ ! -x "$wcns_bin" ]; then
    echo "error: build adr_weno and adr_wcns first" >&2
    exit 1
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

expect_fail() {
    name=$1
    shift
    if ( cd "$tmp" && "$@" >"$tmp/$name.log" 2>&1 ); then
        echo "expected failure: $name" >&2
        cat "$tmp/$name.log" >&2
        exit 1
    fi
    echo "rejected: $name"
}

cp "$root/input.nml" "$tmp/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/scheme = .*/scheme = 'NOT_A_SCHEME',/" \
    "$tmp/input.nml"
expect_fail unknown_weno "$weno_bin"

cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/mesh_num = .*/mesh_num = 826,/" "$tmp/input.nml"
expect_fail even_mesh "$weno_bin"

cp "$root/input.nml" "$tmp/input.nml"
sed -i "s/mesh_num = .*/mesh_num = 3,/" "$tmp/input.nml"
expect_fail mesh_below_5 "$weno_bin"

cp "$root/input.nml" "$tmp/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/dt = .*/dt = -1.0d0,/" \
    "$tmp/input.nml"
expect_fail negative_dt "$weno_bin"

cp "$root/input.nml" "$tmp/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/scheme = .*/scheme = 'LINEAR5_M2P2',/" \
    -e "s/difference_operator = .*/difference_operator = 'NOPE',/" \
    "$tmp/input.nml"
expect_fail unknown_operator "$wcns_bin"

cp "$root/input.nml" "$tmp/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/scheme = .*/scheme = 'NOT_A_SCHEME',/" \
    -e "s/difference_operator = .*/difference_operator = 'HCS4',/" \
    "$tmp/input.nml"
expect_fail unknown_wcns "$wcns_bin"

echo "adr params: PASS"
