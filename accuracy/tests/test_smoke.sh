#!/bin/sh
# End-to-end smoke tests for the accuracy drivers.
#
# Each case writes a fresh input.nml from a template, so the smoke test is
# independent of whatever the user changed in the shipped input.nml.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin1="$root/bin/$precision/accuracy_1d"
bin2="$root/bin/$precision/accuracy_2d"

if [ ! -x "$bin1" ] || [ ! -x "$bin2" ]; then
    echo "error: build accuracy first (bin/$precision/accuracy_1d, accuracy_2d)" >&2
    exit 1
fi

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

write_nml_1d() {
    path=$1
    init=$2
    scheme=$3
    operator=$4
    cat > "$path" <<EOF
&accuracy_nml
  write_sample = .false.
  output_dir = 'results'
/
&accuracy_1d_nml
  initial_1d = '$init'
  nbase = 10
  levels = 2
  cfl = 0.1d0
  tend = 2.0d0
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$operator'
/
EOF
}

write_nml_2d() {
    path=$1
    case_name=$2
    scheme=$3
    operator=$4
    cat > "$path" <<EOF
&accuracy_nml
  write_sample = .false.
  output_dir = 'results'
/
&accuracy_2d_nml
  case_2d = '$case_name'
  nbase = 10
  levels = 1
  cfl = 0.1d0
  tend = 1.0d0
  sine_rho0 = 1.0d0
  sine_amp = 0.2d0
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$operator'
/
EOF
}

run_1d() {
    name=$1
    init=$2
    scheme=$3
    operator=$4
    work="$base/$name"
    mkdir -p "$work"
    write_nml_1d "$work/input.nml" "$init" "$scheme" "$operator"
    ( cd "$work" && "$bin1" > run.log )
}

run_2d() {
    name=$1
    case_name=$2
    scheme=$3
    operator=$4
    work="$base/$name"
    mkdir -p "$work"
    write_nml_2d "$work/input.nml" "$case_name" "$scheme" "$operator"
    ( cd "$work" && "$bin2" > run.log )
}

check_csv() {
    path=$1
    expected=$2
    if [ ! -f "$path" ]; then
        echo "smoke: missing $path" >&2
        exit 1
    fi
    actual=$(wc -l < "$path")
    if [ "$actual" -ne "$expected" ]; then
        echo "smoke: $path has $actual lines, expected $expected" >&2
        exit 1
    fi
    if grep -E 'NaN|NAN|nan|Infinity|infinity' "$path" > /dev/null; then
        echo "smoke: non-finite value in $path" >&2
        exit 1
    fi
}

run_1d smoke_1d sin2 Z5 HCS6
check_csv "$base/smoke_1d/results/csv/accuracy_1d_sin2_Z5_HCS6.csv" 4
echo "accuracy 1D smoke: PASS"

run_2d smoke_vortex vortex Z3 HCS4
check_csv "$base/smoke_vortex/results/csv/accuracy_2d_vortex_Z3_HCS4.csv" 3
echo "accuracy 2D vortex smoke: PASS"

run_2d smoke_sine sine Z3 HCS4
check_csv "$base/smoke_sine/results/csv/accuracy_2d_sine_Z3_HCS4.csv" 3
echo "accuracy 2D sine smoke: PASS"

echo "accuracy smoke: PASS"
