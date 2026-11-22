#!/bin/sh
# Run the quick result set and verify snapshots, logs, and figures.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
cd "$root"

make PRECISION="$precision" all >/dev/null
sh scripts/run_results.sh --quick

test -f results/quick/csv/JS5_HCS6_t1.csv
test -f results/quick/csv/JS3_HCS4_t1.csv
test -f results/quick/figures/JS5_HCS6_evolution.png
test -f results/quick/figures/JS3_HCS4_evolution.png
test -f results/quick/figures/order5_four_periods.png
test -f results/quick/figures/order3_four_periods.png
grep -q 'nan_count=0' results/quick/logs/*.log

if grep -Eqi 'nan|inf' results/quick/csv/*.csv; then
    echo "results: NaN/Inf in quick CSV" >&2
    exit 1
fi

echo "complex_wave results: PASS"
