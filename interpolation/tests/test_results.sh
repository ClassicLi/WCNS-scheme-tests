#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

precision=${PRECISION:-qp}
make PRECISION="$precision" >/dev/null

PRECISION="$precision" sh scripts/run_results.sh --quick

test -f results/convergence/LINEAR5_M2P2_udf3.csv
ls results/samples/*.dat >/dev/null 2>&1

python3 scripts/plot_convergence.py --quick
test -f results/convergence/L1_convergence.png

python3 scripts/plot_diagnostics.py --quick
test -f results/diagnostics/rst10Z5_beta_omega.png

python3 -m py_compile scripts/plot_convergence.py scripts/plot_diagnostics.py

echo "results scripts: PASS"
