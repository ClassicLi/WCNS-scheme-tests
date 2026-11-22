#!/bin/sh
# Check the reference-sidecar comparison used by run_results.sh.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

printf 'mesh_num=1600\nscheme=QR\noperator=HCS4\n' > "$tmp/ref.log"
sh "$root/scripts/check_reference.sh" "$tmp/ref.log" 1600 QR HCS4

if sh "$root/scripts/check_reference.sh" "$tmp/ref.log" 400 QR HCS4; then
    echo "reference check: stale mesh accepted" >&2
    exit 1
fi

printf 'mesh_num=1600\nscheme=Z5\noperator=HCS6\n' > "$tmp/ref2.log"
if sh "$root/scripts/check_reference.sh" "$tmp/ref2.log" 1600 QR HCS4; then
    echo "reference check: stale scheme/operator accepted" >&2
    exit 1
fi

if sh "$root/scripts/check_reference.sh" "$tmp/missing.log" 1600 QR HCS4; then
    echo "reference check: missing log accepted" >&2
    exit 1
fi

echo "euler1d reference check: PASS"
