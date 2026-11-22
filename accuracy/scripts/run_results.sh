#!/bin/sh
# Run the accuracy result matrix.
#
# Usage: run_results.sh [--quick]
#   --quick  small meshes, two 1D profiles, two schemes, results/quick/
#   default  full 1D and 2D matrices, results/
#
# The 1D profile list can be overridden with INITIALS="sin2 gauss ...".
# Each case writes a complete input.nml from a template, so the matrix does
# not depend on the values currently in the shipped input.nml.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin1="$root/bin/$precision/accuracy_1d"
bin2="$root/bin/$precision/accuracy_2d"
mode=${1:-full}

if [ ! -x "$bin1" ] || [ ! -x "$bin2" ]; then
    echo "error: build accuracy first; run make" >&2
    exit 1
fi

case "$mode" in
    --quick)
        initials="sin2 sin"
        nbase1=10
        levels1=4
        nbase2=10
        levels2=2
        pairs="Z3:HCS4 Z5:HCS6"
        dest="$root/results/quick"
        ;;
    *)
        initials=${INITIALS:-"sin2 sin sin3 sin_phase gauss sqrt rational"}
        nbase1=10
        levels1=8
        nbase2=40
        levels2=4
        pairs="Z3:HCS4 Z5:HCS6 QR:HCS4"
        dest="$root/results"
        ;;
esac

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

mkdir -p "$dest/csv" "$dest/logs"

write_nml_1d() {
    path=$1
    init=$2
    scheme=$3
    operator=$4
    nbase=$5
    levels=$6
    cat > "$path" <<EOF
&accuracy_nml
  write_sample = .false.
  output_dir = 'results'
/
&accuracy_1d_nml
  initial_1d = '$init'
  nbase = $nbase
  levels = $levels
  cfl = 0.1d0
  tend = 2.0d0
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$operator'
/
EOF
}

write_nml_2d() {
    path=$1
    case_name=$2
    scheme=$3
    operator=$4
    nbase=$5
    levels=$6
    cat > "$path" <<EOF
&accuracy_nml
  write_sample = .false.
  output_dir = 'results'
/
&accuracy_2d_nml
  case_2d = '$case_name'
  nbase = $nbase
  levels = $levels
  cfl = 0.1d0
  tend = 1.0d0
  sine_rho0 = 1.0d0
  sine_amp = 0.2d0
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$operator'
/
EOF
}

for init in $initials; do
    for pair in $pairs; do
        scheme=${pair%%:*}
        operator=${pair##*:}
        work="$base/1d-$init-$scheme-$operator"
        mkdir -p "$work"
        write_nml_1d "$work/input.nml" "$init" "$scheme" "$operator" "$nbase1" "$levels1"
        ( cd "$work" && "$bin1" > run.log )
        cp "$work"/results/csv/*.csv "$dest/csv/"
        cp "$work"/results/logs/*.log "$dest/logs/"
        echo "done: 1d $init $scheme $operator"
    done
done

for case_name in vortex sine; do
    for pair in $pairs; do
        scheme=${pair%%:*}
        operator=${pair##*:}
        work="$base/2d-$case_name-$scheme-$operator"
        mkdir -p "$work"
        write_nml_2d "$work/input.nml" "$case_name" "$scheme" "$operator" "$nbase2" "$levels2"
        ( cd "$work" && "$bin2" > run.log )
        cp "$work"/results/csv/*.csv "$dest/csv/"
        cp "$work"/results/logs/*.log "$dest/logs/"
        echo "done: 2d $case_name $scheme $operator"
    done
done

echo "accuracy results written to $dest"
