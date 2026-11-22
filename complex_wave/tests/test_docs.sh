#!/bin/sh
# Check the complex-wave documentation for required content.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

test -f README.md
test -f docs/method_zh.md
test -f docs/method_en.md
test -f docs/formula_audit.md

grep -q "complex_wave" README.md
grep -q "WCNSE2" README.md
grep -q "HCS6" README.md
grep -q "periods" README.md
grep -q "RK4" README.md
grep -q "dp" README.md
grep -q "JS3" README.md
grep -q "Z5" README.md
grep -q "200" README.md

grep -q "Runge" docs/method_zh.md
grep -q "Runge" docs/method_en.md
grep -q "64/45" docs/method_zh.md
grep -q "64/45" docs/method_en.md

grep -q "HCS6" docs/formula_audit.md
grep -q "Difference_hcs6" docs/formula_audit.md
grep -q "halo" docs/formula_audit.md
if ! grep -q "越界" docs/formula_audit.md && ! grep -qi "out-of-bounds" docs/formula_audit.md; then
    echo "formula_audit: missing halo out-of-bounds note" >&2
    exit 1
fi

# The top-level README must describe six subprojects and per-project precision.
test -f ../README.md
grep -q "六个自包含的子项目" ../README.md
grep -q "six self-contained subprojects" ../README.md
grep -q "complex_wave/.*dp" ../README.md
grep -q "euler1d/.*dp" ../README.md

echo "complex_wave docs: PASS"
