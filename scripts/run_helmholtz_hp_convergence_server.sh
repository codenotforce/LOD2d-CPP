#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT/results/helmholtz_hp"}
JOBS=${JOBS:-8}
PATCH_THREADS=${PATCH_THREADS:-8}
PROGRESS_INTERVAL=${PROGRESS_INTERVAL:-8}
export OMP_DYNAMIC=FALSE
KAPPA=${KAPPA:-4}
MODE=${MODE:-}
DEEP_FINE_LEVELS=${DEEP_FINE_LEVELS:-"12 14"}
DEEP_H_LEVELS=${DEEP_H_LEVELS:-6,8}
DEGREES=${DEGREES:-1,2,3}
if [[ -z "${ELL_BY_P:-}" ]]; then
    case "$DEGREES" in
    *,*,*) ELL_BY_P=3,3,3 ;;
    *,*) ELL_BY_P=3,3 ;;
    *) ELL_BY_P=3 ;;
    esac
fi
BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_hp_convergence"
IFS=',' read -r -a degree_values <<< "$DEGREES"
IFS=',' read -r -a ell_values <<< "$ELL_BY_P"
if [[ ${#degree_values[@]} -ne ${#ell_values[@]} ]]; then
    printf 'DEGREES and ELL_BY_P must have the same length\n' >&2
    exit 2
fi

if [[ -z "$MODE" ]]; then
    if [[ "${FULL:-0}" == "1" ]]; then
        MODE=full
    else
        MODE=smoke
    fi
fi

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target \
    bench_helmholtz_hp_convergence test_helmholtz_hp_model -j "$JOBS"
"$BUILD_DIR/tests/test_helmholtz_hp_model"

run_case() {
    local name=$1
    shift
    local tmp="$RESULT_DIR/$name.csv.tmp"
    if [[ -f "$RESULT_DIR/$name.done" ]]; then
        printf 'Skipping completed case %s\n' "$name"
        return
    fi
    /usr/bin/time -v "$BENCH" "$@" \
        --threads="$PATCH_THREADS" --progress="$PROGRESS_INTERVAL" \
        --stream --check \
        >"$tmp" 2>"$RESULT_DIR/$name.time"
    mv "$tmp" "$RESULT_DIR/$name.csv"
    touch "$RESULT_DIR/$name.done"
}

case "$MODE" in
smoke)
    run_case smoke_H4_h8 \
        --study=H --k="$KAPPA" --h=8 --p=1 \
        --H-levels=4 --ell-by-p=3
    ;;
full)
    run_case fine_hp \
        --study=fem --k="$KAPPA" --p="$DEGREES" \
        --h-levels=4,6,8,10,12
    for index in "${!degree_values[@]}"; do
        degree=${degree_values[$index]}
        ell=${ell_values[$index]}
        run_case "ell_calibration_p${degree}" \
            --study=ell --k="$KAPPA" --H=4 --h=10 --p="$degree" \
            --ell-levels=1,2,3,4,5
        run_case "H_master_h12_p${degree}" \
            --study=H --k="$KAPPA" --h=12 --p="$degree" \
            --H-levels=2,4,6,8 --ell-by-p="$ell"
        run_case "coupled_gap6_p${degree}" \
            --study=coupled --k="$KAPPA" --p="$degree" \
            --H-levels=2,4,6 --gap=6 --ell-by-p="$ell"
    done
    ;;
deep)
    for fine_level in $DEEP_FINE_LEVELS; do
        for index in "${!degree_values[@]}"; do
            degree=${degree_values[$index]}
            ell=${ell_values[$index]}
            run_case "H_deep_h${fine_level}_p${degree}" \
                --study=H --k="$KAPPA" --h="$fine_level" --p="$degree" \
                --H-levels="$DEEP_H_LEVELS" --ell-by-p="$ell"
        done
    done
    ;;
*)
    printf 'Unknown MODE=%s; expected smoke, full, or deep\n' "$MODE" >&2
    exit 2
    ;;
esac

printf 'hp Helmholtz %s runs completed in %s\n' "$MODE" "$RESULT_DIR"