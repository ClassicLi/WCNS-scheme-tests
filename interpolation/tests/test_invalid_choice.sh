#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

precision=${PRECISION:-qp}
bin="bin/$precision/invalid_choice_driver"
make "$bin" PRECISION="$precision" >/dev/null
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT

if "./$bin" >"$tmp" 2>&1; then
    echo "invalid moc33_linear was accepted" >&2
    cat "$tmp" >&2
    exit 1
fi

grep -qi "invalid moc33_linear" "$tmp"
echo "invalid choice: PASS"
