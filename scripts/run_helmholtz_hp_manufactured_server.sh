#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT/results/helmholtz_hp_manufactured"}
JOBS=${JOBS:-16}
PATCH_THREADS=${PATCH_THREADS:-8}
HP_SOLVER=${HP_SOLVER:-schur}
PROGRESS_INTERVAL=${PROGRESS_INTERVAL:-8}
export OMP_DYNAMIC=FALSE
MODE=${MODE:-smoke}
KAPPA=${KAPPA:-4}
DEGREES=${DEGREES:-1,2,3}
MASTER_FINE_LEVELS=${MASTER_FINE_LEVELS:-"10 12 14"}
H_LEVELS=${H_LEVELS:-4,6,8}
COUPLED_H_LEVELS=${COUPLED_H_LEVELS:-2,4,6}
GAP=${GAP:-6}
BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_hp_convergence"

if [[ -z "${ELL_BY_P:-}" ]]; then
    case "$DEGREES" in
    *,*,*) ELL_BY_P=3,3,3 ;;
    *,*) ELL_BY_P=3,3 ;;
    *) ELL_BY_P=3 ;;
    esac
fi

IFS=',' read -r -a degree_values <<< "$DEGREES"
IFS=',' read -r -a ell_values <<< "$ELL_BY_P"
if [[ ${#degree_values[@]} -ne ${#ell_values[@]} ]]; then
    printf 'DEGREES and ELL_BY_P must have the same length\n' >&2
    exit 2
fi
mkdir -p "$RESULT_DIR"

cat >"$RESULT_DIR/manufactured_solution.txt" <<EOF
u(x,y) = phi(x) phi(y) exp(i k x)
phi(t) = 16 t^2 (1-t)^2
k = $KAPPA
boundary: d_n u - i k u = 0 on the complete boundary

CSV errors:
exact_energy, exact_l2       = ||u_exact - u_LOD||
fine_energy, fine_l2         = ||u_exact - u_h^p||
lod_fine_energy, lod_fine_l2 = ||u_h^p - u_LOD||
EOF

cmake -S "$ROOT" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target \
    bench_helmholtz_hp_convergence test_helmholtz_hp_model -j "$JOBS"
"$BUILD_DIR/tests/test_helmholtz_hp_model"

run_case() {
    local name=$1
    shift
    local csv="$RESULT_DIR/$name.csv"
    local done_file="$RESULT_DIR/$name.done"
    if [[ -f "$done_file" ]]; then
        printf 'Skipping completed case %s\n' "$name"
        return
    fi
    /usr/bin/time -v "$BENCH" "$@" \
        --solver="$HP_SOLVER" --threads="$PATCH_THREADS" \
        --progress="$PROGRESS_INTERVAL" \
        --stream --check \
        >"$csv.tmp" 2>"$RESULT_DIR/$name.time"
    mv "$csv.tmp" "$csv"
    touch "$done_file"
}

run_fixed() {
    for fine_level in $MASTER_FINE_LEVELS; do
        for index in "${!degree_values[@]}"; do
            local degree=${degree_values[$index]}
            local ell=${ell_values[$index]}
            run_case "manufactured_fixed_h${fine_level}_p${degree}" \
                --study=H --k="$KAPPA" --h="$fine_level" --p="$degree" \
                --H-levels="$H_LEVELS" --ell-by-p="$ell"
        done
    done
}

run_coupled() {
    for index in "${!degree_values[@]}"; do
        local degree=${degree_values[$index]}
        local ell=${ell_values[$index]}
        run_case "manufactured_coupled_gap${GAP}_p${degree}" \
            --study=coupled --k="$KAPPA" --p="$degree" \
            --H-levels="$COUPLED_H_LEVELS" --gap="$GAP" \
            --ell-by-p="$ell"
    done
}

case "$MODE" in
smoke)
    run_case manufactured_smoke_p1 \
        --study=H --k="$KAPPA" --h=8 --p=1 \
        --H-levels=4 --ell-by-p=3
    ;;
fixed)
    run_fixed
    ;;
coupled)
    run_coupled
    ;;
all)
    run_fixed
    run_coupled
    ;;
*)
    printf 'Unknown MODE=%s; expected smoke, fixed, coupled, or all\n' \
        "$MODE" >&2
    exit 2
    ;;
esac

printf 'Manufactured-solution %s runs completed in %s\n' \
    "$MODE" "$RESULT_DIR"
