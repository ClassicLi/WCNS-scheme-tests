#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

test -f README.md
test -f docs/method_zh.md
test -f docs/method_en.md
test -f docs/formula_audit.md

grep -q "WCNS" README.md
if grep -q "WENO test" README.md; then
    echo "README describes the project as a WENO test" >&2
    exit 1
fi

grep -q "MOC322" README.md
grep -q "MOC322" docs/method_zh.md
grep -q "MOC322" docs/method_en.md
grep -q "nonlinear interpolation" README.md
grep -q "same_order" README.md
grep -q "same_order" docs/method_zh.md

grep -q "results-table" README.md
grep -q "L1_convergence.png" README.md
grep -q "rst10Z5_beta_omega.png" README.md
grep -q "PRECISION" README.md

echo "docs: PASS"
