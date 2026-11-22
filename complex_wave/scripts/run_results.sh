#!/bin/sh
# Run the complex-wave result set and generate waveform figures.
#
# Usage: run_results.sh [--quick]
#   --quick  mesh 40, one period, two schemes, results/quick/
#   default  mesh 200, four periods, five schemes, results/
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/complex_wave"
mode=${1:-full}

if [ ! -x "$bin" ]; then
    echo "error: $bin not found; run make first" >&2
    exit 1
fi

case "$mode" in
    --quick)
        mesh=40
        periods=1
        dest="$root/results/quick"
        cases="JS5:HCS6 JS3:HCS4"
        ;;
    *)
        mesh=200
        periods=4
        dest="$root/results"
        cases="JS5:HCS6 Z5:HCS6 JS3:HCS4 Z3:HCS4 QR:HCS4"
        ;;
esac

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

mkdir -p "$dest/csv" "$dest/dat" "$dest/logs" "$dest/figures"

for pair in $cases; do
    scheme=${pair%%:*}
    operator=${pair##*:}
    work="$base/$scheme-$operator"
    mkdir -p "$work"
    cp "$root/input.nml" "$work/input.nml"
    sed -i \
        -e "s/mesh_num = .*/mesh_num = $mesh,/" \
        -e "s/periods = .*/periods = $periods,/" \
        -e "s/scheme = .*/scheme = '$scheme',/" \
        -e "s/difference_operator = .*/difference_operator = '$operator',/" \
        "$work/input.nml"
    ( cd "$work" && "$bin" >/dev/null )
    cp "$work"/results/csv/*.csv "$dest/csv/"
    cp "$work"/results/dat/*.dat "$dest/dat/"
    cp "$work"/results/logs/*.log "$dest/logs/"
    echo "done: $scheme + $operator"
done

python3 "$root/scripts/plot_waveforms.py" --results-dir "$dest"
