#!/bin/sh
# Run the quick matrix and verify CSVs and figures.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
cd "$root"

make PRECISION="$precision" all >/dev/null
sh scripts/run_results.sh --quick

test -f results/quick/csv/WENO5_LINEAR.csv
test -f results/quick/figures/quick_dispersion.png
test -f results/quick/figures/quick_dissipation.png
grep -q "nan_count" results/quick/logs/*.log

if grep -Eqi 'nan|inf' results/quick/csv/*.csv; then
    echo "results: NaN/Inf in quick CSV" >&2
    exit 1
fi

echo "adr results: PASS"
