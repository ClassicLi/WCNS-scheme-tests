#!/bin/sh
# Check the euler1d documentation for required content.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

test -f README.md
test -f docs/method_zh.md
test -f docs/method_en.md
test -f docs/formula_audit.md

grep -q "sod" README.md
grep -q "lax" README.md
grep -q "shu_osher" README.md
grep -q "QR" README.md
grep -q "HCS4" README.md
grep -q "HCS6" README.md
grep -q "WCNSE2" README.md
grep -q "halo" README.md
grep -q "dp" README.md
grep -q "qp" README.md

grep -q "Steger" docs/method_zh.md
grep -q "Steger" docs/method_en.md
grep -q "Runge" docs/method_zh.md
grep -q "Runge" docs/method_en.md
grep -q "x - 0.5" docs/method_zh.md
grep -q "x - 0.5" docs/method_en.md

grep -q "MOC53" docs/formula_audit.md
grep -q "WCNSE2" docs/formula_audit.md
grep -q "HCS4" docs/formula_audit.md
grep -q "HCS6" docs/formula_audit.md
grep -q "halo" docs/formula_audit.md

# Public files must not carry internal implementation-plan references.
if grep -R -E "Task [0-9]" README.md docs Makefile; then
    echo "euler1d docs: internal Task reference in public files" >&2
    exit 1
fi

echo "euler1d docs: PASS"
