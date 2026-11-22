#!/bin/sh
# Check the accuracy documentation for required content.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

test -f README.md
test -f docs/method_zh.md
test -f docs/method_en.md
test -f docs/formula_audit.md

grep -q "accuracy" README.md
grep -q "Z3" README.md
grep -q "Z5" README.md
grep -q "QR" README.md
grep -q "HCS4" README.md
grep -q "HCS6" README.md
grep -q "vortex" README.md
grep -q "sine" README.md
grep -q "halo" README.md
grep -q "dp" README.md
grep -q "qp" README.md

grep -q "Runge" docs/method_zh.md
grep -q "Runge" docs/method_en.md
grep -q "sin2" docs/method_zh.md
grep -q "sin2" docs/method_en.md
grep -q "HCS4" docs/formula_audit.md
grep -q "HCS6" docs/formula_audit.md
grep -q "Eint" docs/formula_audit.md
grep -q "halo" docs/formula_audit.md
grep -q "LINEAR5" docs/formula_audit.md
if ! grep -q "临界" docs/formula_audit.md && ! grep -qi "critical" docs/formula_audit.md; then
    echo "formula_audit: missing critical-point note" >&2
    exit 1
fi

echo "accuracy docs: PASS"
