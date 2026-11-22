#!/bin/sh
# Run the euler2d result matrix and write fields and logs.
#
# Usage: run_results.sh [--quick|--full|--smoke]
#   --quick  nz=32, tend=0.01, writes results/quick/ (default)
#   --full   case-table nz and case-default tend, writes results/full/
#   --smoke  nz=16, tend=0.001, writes results/smoke/ (committed samples)
#
# Each case runs in its own temporary directory with a dedicated input.nml.
# The fields and key=value logs are copied to the selected output directory.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-dp}
bin="$root/bin/$precision/euler2d"
mode=${1:---quick}

if [ ! -x "$bin" ]; then
    echo "error: build euler2d first; run make -C $root" >&2
    exit 1
fi

case "$mode" in
    --quick)
        dest="$root/results/quick"
        nz=32
        tend='0.01d0'
        sb_tend2='0.005d0'
        ;;
    --full)
        dest="$root/results/full"
        nz=''
        tend='-1.0d0'
        sb_tend2='-1.0d0'
        ;;
    --smoke)
        dest="$root/results/smoke"
        nz=16
        tend='0.001d0'
        sb_tend2='0.0005d0'
        ;;
    *)
        echo "error: unknown mode '$mode'; use --quick, --full, or --smoke" >&2
        exit 2
        ;;
esac

base=$(mktemp -d)
trap 'rm -rf "$base"' EXIT
mkdir -p "$dest"

# Write one input.nml. $8 is the scheme-specific parameter block; the
# time strings already carry the Fortran double-precision exponent.
write_nml() {
    path=$1
    case_name=$2
    scheme=$3
    op=$4
    mesh=$5
    tend_val=$6
    tend2_val=$7
    params=$8
    cat > "$path" <<EOF
&euler2d_nml
  case_2d = '$case_name'
  mesh_num = $mesh
  cfl = 0.1d0
  tend = $tend_val
  tend2 = $tend2_val
  output_dir = 'results'
  write_sample = .false.
/
&scheme_nml
  scheme = '$scheme'
  $params
/
&difference_nml
  difference_operator = '$op'
/
EOF
}

# Run one case and copy its fields and log to the destination directory.
run_case() {
    case_name=$1
    scheme=$2
    op=$3
    full_nz=$4
    params=$5

    if [ -n "$nz" ]; then
        mesh=$nz
    else
        mesh=$full_nz
    fi
    if [ "$case_name" = 'shock_bubble' ]; then
        case_tend2=$sb_tend2
    else
        case_tend2='0.0d0'
    fi

    work="$base/$case_name"
    mkdir -p "$work"
    write_nml "$work/input.nml" "$case_name" "$scheme" "$op" "$mesh" \
              "$tend" "$case_tend2" "$params"
    ( cd "$work" && "$bin" > run.log 2>&1 )
    cp "$work"/results/*.dat "$work"/results/*.log "$dest/"
    echo "done: $case_name ($scheme/$op, nz=$mesh)"
}

# Case table: case|scheme|operator|full nz|scheme parameters.
#
# R2/SB use the thesis MOC33 exponent p=4.  The repository default
# moc33_eps=1e-80 makes alpha**p overflow in constant regions, so these
# cases pin the thesis value moc33_eps=1e-40.
while IFS='|' read -r case_name scheme op full_nz params; do
    [ -n "$case_name" ] || continue
    run_case "$case_name" "$scheme" "$op" "$full_nz" "$params"
done <<'CASES'
double_mach|MOZ|HCS6|200|moz_M_q = 1, moz_N_q = 4, moz_eps = 1.0d-40
kelvin_helmholtz|Z5|HCS6|400|z5_q = 4.0d0, z5_eps = 1.0d-40
rayleigh_taylor|MOC53|HCS6|400|moc53_N_q = 2, moc53_p = 2, moc53_eps = 1.0d-60, moc53_linear = 'GA'
riemann1|MOC53|HCS6|400|moc53_N_q = 2, moc53_p = 2, moc53_eps = 1.0d-60, moc53_linear = 'GA'
riemann2|MOC33|HCS4|1600|moc33_p = 4, moc33_eps = 1.0d-40, moc33_linear = 'MQ'
shock_bubble|MOC33|HCS4|800|moc33_p = 4, moc33_eps = 1.0d-40, moc33_linear = 'MQ'
shock_vortex|Z3|HCS4|180|z3_q = 2.0d0, z3_eps = 1.0d-40
CASES

echo "euler2d results written to $dest"
