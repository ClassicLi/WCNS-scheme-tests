#!/bin/sh
# End-to-end smoke tests for the complex-wave driver.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/complex_wave"

if [ ! -x "$bin" ]; then
    echo "error: build complex_wave first (bin/$precision/complex_wave)" >&2
    exit 1
fi

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

run_case() {
    name=$1
    scheme=$2
    operator=$3
    mesh=$4
    periods=$5
    work="$base/$name"
    mkdir -p "$work"
    cp "$root/input.nml" "$work/input.nml"
    sed -i \
        -e "s/mesh_num = .*/mesh_num = $mesh,/" \
        -e "s/periods = .*/periods = $periods,/" \
        -e "s/scheme = .*/scheme = '$scheme',/" \
        -e "s/difference_operator = .*/difference_operator = '$operator',/" \
        "$work/input.nml"
    ( cd "$work" && "$bin" >/dev/null )
}

# Reference case with the default mesh and one period.
run_case z5_hcs6 Z5 HCS6 200 1
csv="$base/z5_hcs6/results/csv/Z5_HCS6_t1.csv"
test -f "$base/z5_hcs6/results/csv/Z5_HCS6_t0.csv"
test -f "$csv"
test "$(wc -l < "$csv")" -eq 202
test "$(grep -c 'Zone' "$base/z5_hcs6/results/dat/Z5_HCS6.dat")" -eq 2
grep -q 'nan_count=0' "$base/z5_hcs6/results/logs/Z5_HCS6.log"
if grep -Eqi 'nan|inf' "$base"/z5_hcs6/results/csv/Z5_HCS6_t*.csv; then
    echo "smoke: NaN/Inf in Z5_HCS6 snapshots" >&2
    exit 1
fi
echo "complex_wave Z5_HCS6 one period: PASS"

# Three periods must produce exactly t0..t3.
run_case z3_hcs4 Z3 HCS4 40 3
for k in 0 1 2 3; do
    test -f "$base/z3_hcs4/results/csv/Z3_HCS4_t$k.csv"
done
test ! -f "$base/z3_hcs4/results/csv/Z3_HCS4_t4.csv"
test "$(grep -c 'Zone' "$base/z3_hcs4/results/dat/Z3_HCS4.dat")" -eq 4
echo "complex_wave Z3_HCS4 three periods: PASS"

echo "complex_wave smoke: PASS"
