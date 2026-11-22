#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
bin="$root/bin/$precision/interpolation_test"
mode=${1:-full}

if [ ! -x "$bin" ]; then
    echo "error: $bin not found; run make first" >&2
    exit 1
fi

case "$mode" in
    --quick)
        schemes="LINEAR5_M2P2"
        funcs="udf3"
        iters=3
        ;;
    --quick-matrix)
        schemes="MUSCL JS3 Z3 JS5 Z5 QR MOC322 MOC33 MOC53 MOZ \
LINEAR3_M1P1 LINEAR3_M2P0 LINEAR3_0P2 LINEAR4_M2P1 LINEAR4_M1P2 \
LINEAR5_M2P2 LINEAR6_M2P3 RBF_GA3 RBF_MQ3 RBF_IQ3 RBF_GA5 RBF_MQ5 RBF_IQ5"
        funcs="udf3"
        iters=3
        ;;
    *)
        schemes="MUSCL JS3 Z3 JS5 Z5 QR MOC322 MOC33 MOC53 MOZ LINEAR5_M2P2 RBF_GA5"
        funcs="udf3 cp1 cp2 udf4"
        iters=10
        ;;
esac

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

mkdir -p "$root/results/convergence" "$root/results/samples"

run_case() {
    scheme=$1
    func=$2

    case "$func" in
        udf3) tft=3; udf=3; ncp=1 ;;
        udf4) tft=3; udf=4; ncp=1 ;;
        cp1)  tft=1; udf=3; ncp=1 ;;
        cp2)  tft=1; udf=3; ncp=2 ;;
        *) echo "unknown function $func" >&2; exit 1 ;;
    esac

    work="$tmp/$scheme-$func"
    mkdir -p "$work/results/convergence"
    cp "$root/input.nml" "$work/input.nml"
    sed -i \
        -e "s/scheme = .*/scheme = '$scheme',/" \
        -e "s/test_func_type = .*/test_func_type = $tft,/" \
        -e "s/n_order_critical_point = .*/n_order_critical_point = $ncp,/" \
        -e "s/udf_id = .*/udf_id = $udf,/" \
        -e "s/iter_times = .*/iter_times = $iters,/" \
        "$work/input.nml"

    ( cd "$work" && "$bin" >/dev/null )
    cp "$work"/*_log.dat "$root/results/convergence/" 2>/dev/null || true
    cp "$work"/results/convergence/*.csv "$root/results/convergence/" 2>/dev/null || true
    cp "$work"/rst10*.dat "$root/results/samples/" 2>/dev/null || true
}

for s in $schemes; do
    for f in $funcs; do
        echo "run $s $f"
        run_case "$s" "$f"
    done
done

echo "results written to $root/results"
