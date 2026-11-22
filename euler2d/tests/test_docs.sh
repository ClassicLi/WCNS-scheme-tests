#!/bin/sh
# Documentation checks for the euler2d subproject.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

for f in README.md docs/method_zh.md docs/method_en.md \
         docs/formula_audit.md tests/test_docs.sh; do
    test -f "$f" || { echo "missing file: $f" >&2; exit 1; }
done

# README keyword coverage.
for word in HCS4 HCS6 WCNSE2 QR dp qp halo; do
    grep -q "$word" README.md || {
        echo "README.md missing keyword: $word" >&2
        exit 1
    }
done

# Both method documents must name the split and the time integrator.
for f in docs/method_zh.md docs/method_en.md; do
    for word in Steger Runge; do
        grep -q "$word" "$f" || {
            echo "$f missing keyword: $word" >&2
            exit 1
        }
    done
done

# Audit keyword coverage.
for word in HCS6 MOC53 THINC; do
    grep -q "$word" docs/formula_audit.md || {
        echo "docs/formula_audit.md missing keyword: $word" >&2
        exit 1
    }
done
grep -Fq 'mt = 7' docs/formula_audit.md || {
    echo "docs/formula_audit.md missing: mt = 7" >&2
    exit 1
}

# Public documentation and scripts must not leak internal task labels.
task_pattern='Task [0-9]'
for f in README.md docs/method_zh.md docs/method_en.md \
         docs/formula_audit.md Makefile input.nml; do
    if grep -Eq "$task_pattern" "$f"; then
        echo "internal task label found in $f" >&2
        exit 1
    fi
done
for f in scripts/* tests/*.sh; do
    [ -f "$f" ] || continue
    [ "$f" = "tests/test_docs.sh" ] && continue
    if grep -Eq "$task_pattern" "$f"; then
        echo "internal task label found in $f" >&2
        exit 1
    fi
done

# Public docs and scripts must not contain machine-specific absolute paths.
for f in README.md docs/method_zh.md docs/method_en.md \
         docs/formula_audit.md Makefile input.nml scripts/* tests/*.sh; do
    [ -f "$f" ] || continue
    [ "$f" = "tests/test_docs.sh" ] && continue
    if grep -Eq '(/mnt/|/home/|[A-Za-z]:\\)' "$f"; then
        echo "machine-specific path found in $f" >&2
        exit 1
    fi
done

# Thesis figures: nine PDFs, each with a PNG preview.
for base in Double_mach_QR_1920_480 KH_800_QR RT_240_960_QR \
            RM_type1_800_QR RM_Type2_1600_QR \
            Shock_bubble_1600_800_QR_t02 Shock_bubble_1600_800_QR_t04 \
            SV_360_180_all Double_mach_Z3_1920_480; do
    test -f "docs/figures/$base.pdf" || {
        echo "missing figure: docs/figures/$base.pdf" >&2
        exit 1
    }
    test -f "docs/figures/$base.png" || {
        echo "missing preview: docs/figures/$base.png" >&2
        exit 1
    }
done

pdf_count=$(find docs/figures -maxdepth 1 -type f -name '*.pdf' | wc -l | tr -d ' ')
png_count=$(find docs/figures -maxdepth 1 -type f -name '*.png' | wc -l | tr -d ' ')
if [ "$pdf_count" -lt 8 ] || [ "$png_count" -lt 8 ]; then
    echo "docs/figures has too few figures: $pdf_count PDF, $png_count PNG" >&2
    exit 1
fi

echo "docs: PASS ($pdf_count PDF, $png_count PNG)"
