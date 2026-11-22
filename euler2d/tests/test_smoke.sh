#!/bin/sh
# End-to-end smoke tests for the euler2d driver.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/euler2d"
[ -x "$bin" ] || { echo "error: build euler2d first" >&2; exit 1; }

base=$(mktemp -d); trap 'rm -rf "$base"' EXIT

run_case() {
    case_name=$1; scheme=$2; op=$3; mesh=$4; tend=$5; tend2=$6; nx=$7
    dir="$base/$case_name-$scheme-$op"
    mkdir -p "$dir"
    cat > "$dir/input.nml" <<EOF
&euler2d_nml
  case_2d = '$case_name'
  mesh_num = $mesh
  cfl = 0.1d0
  tend = $tend
  tend2 = $tend2
  output_dir = 'out'
  write_sample = .false.
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$op'
/
EOF
    ( cd "$dir" && "$bin" > run.log 2>&1 )
    test -f "$dir/out/${case_name}_${mesh}x${nx}_t1.dat"
    grep -q 'nan_count=0' "$dir/out/$case_name.log"
    grep -q '^cell_corrections=[0-9]' "$dir/out/$case_name.log"
    # The log records the resolved operator: AUTO with QR (order 3) is HCS4.
    expected_op=$op
    if [ "$op" = AUTO ]; then expected_op=HCS4; fi
    grep -q "^operator=$expected_op$" "$dir/out/$case_name.log"
    echo "PASS $case_name $scheme $op"
}

run_case double_mach       QR AUTO 16 0.001d0 -1.0d0 64
run_case kelvin_helmholtz  QR AUTO 16 0.001d0 -1.0d0 16
run_case rayleigh_taylor   QR AUTO 16 0.001d0 -1.0d0 4
run_case riemann1          QR AUTO 16 0.001d0 -1.0d0 16
run_case riemann2          QR AUTO 16 0.001d0 -1.0d0 16
run_case shock_bubble      QR AUTO 16 0.001d0  0.0005d0 32
run_case shock_vortex      QR AUTO 16 0.001d0 -1.0d0 32
run_case rayleigh_taylor   QR HCS4   16 0.001d0 -1.0d0 4
run_case rayleigh_taylor   QR HCS6   16 0.001d0 -1.0d0 4
run_case rayleigh_taylor   QR WCNSE6 16 0.001d0 -1.0d0 4

test -f "$base/shock_bubble-QR-AUTO/out/shock_bubble_16x32_t2.dat"
echo "euler2d smoke: PASS"
