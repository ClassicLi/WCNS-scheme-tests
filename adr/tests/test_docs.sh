#!/bin/sh
# Check the ADR documentation for required content.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

test -f README.md
test -f docs/method_zh.md
test -f docs/method_en.md
test -f docs/formula_audit.md

grep -q "ADR" README.md
grep -q "WENO" README.md
grep -q "WCNS" README.md
grep -q "WCNSE2" README.md
grep -q "HCS4" README.md
grep -q "HCS6" README.md
grep -q "奇数" README.md
grep -q "odd" README.md
grep -q "PRECISION" README.md
grep -q "qp" README.md
grep -q "dp" README.md
grep -q "Pirozzoli" README.md
grep -q "Pirozzoli" docs/method_zh.md
grep -q "THINC" README.md
grep -q "TENO_333_RBF" README.md
grep -q "WENO5_DLR" README.md
grep -q "同一插值" README.md
grep -q "same interpolation" README.md
grep -q "HCS6" docs/formula_audit.md
grep -q "WENO5_Z" docs/formula_audit.md
grep -q "MOC53" docs/formula_audit.md
grep -q "TENO_MOD_SEMI_ENO_RBF" docs/formula_audit.md
grep -q "旧版" docs/formula_audit.md
grep -q "1000" docs/formula_audit.md

echo "adr docs: PASS"
