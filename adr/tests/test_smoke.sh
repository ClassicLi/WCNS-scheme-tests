#!/bin/sh
# End-to-end smoke tests for adr_weno and adr_wcns.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
weno_bin="$root/bin/$precision/adr_weno"
wcns_bin="$root/bin/$precision/adr_wcns"
base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

run_weno() {
    name=$1
    scheme=$2
    work="$base/$name"
    mkdir -p "$work"
    cp "$root/input.nml" "$work/input.nml"
    sed -i \
        -e "s/mesh_num = .*/mesh_num = 33,/" \
        -e "s/scheme = .*/scheme = '$scheme',/" \
        "$work/input.nml"
    ( cd "$work" && "$weno_bin" )
    csv="$work/results/csv/$scheme.csv"
    test -f "$csv"
    lines=$(wc -l < "$csv")
    test "$lines" -eq 18
    if grep -Eqi 'nan|inf' "$csv"; then
        echo "smoke: NaN/Inf in $csv" >&2
        exit 1
    fi
    echo "adr_weno $scheme smoke: PASS"
}

# Regular WENO cases, including the TENO and Z-type weights.
run_weno weno5_linear WENO5_LINEAR
run_weno weno5_z WENO5_Z
run_weno teno TENO_MOD_SEMI_ENO

# WENO5_DLR is experimental and can produce NaN; only require a clean exit.
work="$base/dlr"
mkdir -p "$work"
cp "$root/input.nml" "$work/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/scheme = .*/scheme = 'WENO5_DLR',/" \
    "$work/input.nml"
( cd "$work" && "$weno_bin" )
echo "adr_weno WENO5_DLR bounds-only: PASS"

# Optional Tecplot sample.
work="$base/sample"
mkdir -p "$work"
cp "$root/input.nml" "$work/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/scheme = .*/scheme = 'WENO5_LINEAR',/" \
    -e "s/write_tecplot_sample = .*/write_tecplot_sample = .true./" \
    "$work/input.nml"
( cd "$work" && "$weno_bin" )
test -f "$work/results/samples/WENO5_LINEAR.dat"
grep -q "ZONE" "$work/results/samples/WENO5_LINEAR.dat"
echo "adr_weno Tecplot sample: PASS"

# WCNS smoke case.
work="$base/wcns"
mkdir -p "$work"
cp "$root/input.nml" "$work/input.nml"
sed -i \
    -e "s/mesh_num = .*/mesh_num = 33,/" \
    -e "s/scheme = .*/scheme = 'LINEAR5_M2P2',/" \
    -e "s/difference_operator = .*/difference_operator = 'HCS4',/" \
    "$work/input.nml"
( cd "$work" && "$wcns_bin" )
csv="$work/results/csv/LINEAR5_M2P2_HCS4.csv"
test -f "$csv"
lines=$(wc -l < "$csv")
test "$lines" -eq 18
if grep -Eqi 'nan|inf' "$csv"; then
    echo "smoke: NaN/Inf in $csv" >&2
    exit 1
fi
echo "adr_wcns smoke: PASS"
