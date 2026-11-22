#!/bin/sh
# Run the euler1d result matrix and generate comparison figures.
#
# Usage: run_results.sh [--quick|--force]
#   --quick  mesh 100, reference mesh 400, two schemes, results/quick/
#   --force  recompute the Shu-Osher reference solution
#   default  mesh 200, reference mesh 1600, three schemes, results/
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/euler1d"
bin_exact="$root/bin/$precision/euler1d_exact"
mode=${1:-full}

if [ ! -x "$bin" ] || [ ! -x "$bin_exact" ]; then
    echo "error: build euler1d first; run make" >&2
    exit 1
fi

force=0
case "$mode" in
    --quick)
        mesh=100
        ref_mesh=400
        dest="$root/results/quick"
        pairs="QR:HCS4 Z5:HCS6"
        ;;
    --force)
        force=1
        mesh=200
        ref_mesh=1600
        dest="$root/results"
        pairs="MUSCL:WCNSE2 QR:HCS4 Z5:HCS6"
        ;;
    *)
        mesh=200
        ref_mesh=1600
        dest="$root/results"
        pairs="MUSCL:WCNSE2 QR:HCS4 Z5:HCS6"
        ;;
esac

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT

mkdir -p "$dest/csv" "$dest/logs" "$dest/figures"

write_nml() {
    path=$1
    case_name=$2
    mesh_num=$3
    scheme=$4
    operator=$5
    cat > "$path" <<EOF
&euler1d_nml
  case_1d = '$case_name'
  mesh_num = $mesh_num
  cfl = 0.1d0
  tend = -1.0d0
  output_dir = 'results'
  write_sample = .false.
/
&scheme_nml
  scheme = '$scheme'
/
&difference_nml
  difference_operator = '$operator'
/
EOF
}

run_numeric() {
    case_name=$1
    mesh_num=$2
    scheme=$3
    operator=$4
    work="$base/$case_name-$scheme-$operator"
    mkdir -p "$work"
    write_nml "$work/input.nml" "$case_name" "$mesh_num" "$scheme" "$operator"
    ( cd "$work" && "$bin" > run.log )
    cp "$work"/results/csv/*.csv "$dest/csv/"
    cp "$work"/results/logs/*.log "$dest/logs/"
    echo "done: $case_name $scheme $operator"
}

run_exact() {
    case_name=$1
    work="$base/exact-$case_name"
    mkdir -p "$work"
    write_nml "$work/input.nml" "$case_name" "$mesh" QR HCS4
    ( cd "$work" && "$bin_exact" > run.log )
    cp "$work"/results/csv/exact_*.csv "$dest/csv/"
    echo "done: exact $case_name"
}

run_reference() {
    ref_csv="$dest/csv/reference_shu_osher.csv"
    if [ "$force" -eq 0 ] && [ -f "$ref_csv" ]; then
        if sh "$root/scripts/check_reference.sh" \
              "$dest/logs/reference_shu_osher.log" "$ref_mesh" QR HCS4; then
            echo "reuse: reference_shu_osher.csv"
        else
            echo "warning: reference_shu_osher.csv does not match" \
                 "mesh=$ref_mesh QR HCS4; rerun with --force" >&2
            echo "reuse: reference_shu_osher.csv (stale)"
        fi
        return
    fi
    work="$base/reference-shu-osher"
    mkdir -p "$work"
    write_nml "$work/input.nml" shu_osher "$ref_mesh" QR HCS4
    ( cd "$work" && "$bin" > run.log )
    cp "$work/results/csv/shu_osher_QR_HCS4.csv" "$ref_csv"
    cp "$work/results/logs/shu_osher_QR_HCS4.log" \
       "$dest/logs/reference_shu_osher.log"
    echo "done: shu_osher reference mesh=$ref_mesh"
}

for case_name in sod lax shu_osher; do
    for pair in $pairs; do
        scheme=${pair%%:*}
        operator=${pair##*:}
        run_numeric "$case_name" "$mesh" "$scheme" "$operator"
    done
done

run_exact sod
run_exact lax
run_reference

python3 "$root/scripts/plot_profiles.py" --results-dir "$dest"
echo "euler1d results written to $dest"
