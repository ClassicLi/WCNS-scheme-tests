#!/bin/sh
# End-to-end smoke tests for the euler1d drivers.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/euler1d"
bin_exact="$root/bin/$precision/euler1d_exact"

if [ ! -x "$bin" ] || [ ! -x "$bin_exact" ]; then
    echo "error: build euler1d first (make)" >&2
    exit 1
fi

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

write_nml() {
    path=$1
    case_name=$2
    mesh=$3
    scheme=$4
    operator=$5
    cat > "$path" <<EOF
&euler1d_nml
  case_1d = '$case_name'
  mesh_num = $mesh
  cfl = 0.1d0
  tend = -1.0d0
  output_dir = 'results'
  write_sample = .false.
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$operator'
/
EOF
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

check_log() {
    path=$1
    if [ ! -f "$path" ]; then
        echo "smoke: missing $path" >&2
        exit 1
    fi
    grep -q 'nan_count=0' "$path" || { echo "smoke: nan_count in $path" >&2; exit 1; }
    grep -q 'corrections=' "$path" || { echo "smoke: corrections missing" >&2; exit 1; }
    grep -q 'steps=' "$path" || { echo "smoke: steps missing" >&2; exit 1; }
}

run_case() {
    case_name=$1
    mesh=$2
    scheme=$3
    operator=$4
    work="$base/$case_name"
    mkdir -p "$work"
    write_nml "$work/input.nml" "$case_name" "$mesh" "$scheme" "$operator"
    ( cd "$work" && "$bin" > run.log )
    check_csv "$work/results/csv/${case_name}_${scheme}_${operator}.csv" $((mesh + 2))
    check_log "$work/results/logs/${case_name}_${scheme}_${operator}.log"
}

run_exact() {
    case_name=$1
    mesh=$2
    left_rho=$3
    right_rho=$4
    work="$base/exact_$case_name"
    mkdir -p "$work"
    write_nml "$work/input.nml" "$case_name" "$mesh" QR HCS4
    ( cd "$work" && "$bin_exact" > run.log )
    csv="$work/results/csv/exact_${case_name}.csv"
    check_csv "$csv" $((mesh + 2))
    first=$(awk -F, 'NR==2 {print $2+0}' "$csv")
    last=$(awk -F, "NR==$((mesh + 2)) {print \$2+0}" "$csv")
    if [ "$first" != "$left_rho" ]; then
        echo "smoke: exact first rho $first != $left_rho" >&2
        exit 1
    fi
    if [ "$last" != "$right_rho" ]; then
        echo "smoke: exact last rho $last != $right_rho" >&2
        exit 1
    fi
}

run_case sod 40 QR HCS4
run_case lax 40 QR HCS4
run_case shu_osher 40 QR HCS4

run_exact sod 40 1 0.125
run_exact lax 40 0.445 0.5

echo "euler1d smoke: PASS"
