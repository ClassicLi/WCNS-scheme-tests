#!/bin/sh
# Run the curated ADR result matrix and plot the figures.
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
precision=${PRECISION:-qp}
mode=${1:-full}

case "$mode" in
    --quick)
        dest="$root/results/quick"
        mesh=33
        groups="quick"
        ;;
    --all-schemes)
        dest="$root/results/all"
        mesh=827
        groups="weno3 weno5 wcns3 wcns5 operators"
        ;;
    *)
        dest="$root/results"
        mesh=827
        groups="weno3 weno5 wcns3 wcns5 operators"
        ;;
esac

export ADR_MESH=$mesh
make -C "$root" PRECISION="$precision" all >/dev/null

run() {
    sh "$root/scripts/run_cases.sh" "$1" "$2" "${3:-}" "$dest"
}

if [ "$mode" = "--quick" ]; then
    run weno WENO5_LINEAR
    run weno WENO5_Z
    run wcns LINEAR5_M2P2 HCS4
    run wcns JS3 HCS4
elif [ "$mode" = "--all-schemes" ]; then
    for s in WENO5_LINEAR WENO3_LINEAR WENO3_JS WENO5_JS WENO5_M WENO5_Z WENO_C \
             BAEZA_WENO3 WENO3_LLR WENO3_TRIG WENO3_GA_RBF_DEFINE \
             WENO3_GA_RBF_UNDEFINE RBF5_LINEAR TENO_333 TENO_MOD_SEMI_ENO \
             TENO_MOD_SEMI_ENO_RBF; do
        run weno "$s"
    done
    for s in MUSCL JS3 Z3 JS5 Z5 QR MOC322 MOC33 MOC53 MOZ LINEAR3_M1P1 \
             LINEAR3_M2P0 LINEAR3_0P2 LINEAR4_M2P1 LINEAR4_M1P2 LINEAR5_M2P2 \
             LINEAR6_M2P3 RBF_GA3 RBF_MQ3 RBF_IQ3 RBF_GA5 RBF_MQ5 RBF_IQ5; do
        run wcns "$s" HCS4
    done
else
    for s in WENO3_LINEAR WENO3_JS WENO3_LLR BAEZA_WENO3 TENO_MOD_SEMI_ENO WENO_C; do
        run weno "$s"
    done
    for s in WENO5_LINEAR WENO5_JS WENO5_M WENO5_Z TENO_333; do
        run weno "$s"
    done
    for s in LINEAR3_M1P1 JS3 Z3 QR MOC322 MOC33; do
        run wcns "$s" HCS4
    done
    for s in LINEAR5_M2P2 JS5 Z5 MOC53 MOZ RBF_GA5; do
        run wcns "$s" HCS4
    done
    for s in LINEAR5_M2P2 Z5 MOC53; do
        for op in WCNSE2 WCNSE4 WCNSE6 HCS4; do
            run wcns "$s" "$op"
        done
    done
fi

for g in $groups; do
    python3 "$root/scripts/plot_adr.py" --dir "$dest" --group "$g"
done

echo "results written to $dest"
