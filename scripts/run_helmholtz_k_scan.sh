#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_k_scan"}
K_VALUES=${K_VALUES:-"4 8 16 32 64"}
FINE_GAP=${FINE_GAP:-8}
KH_TARGET=${KH_TARGET:-1.0}
MODE=${MODE:-two-sided}
STABILITY_MAX_DOFS=${STABILITY_MAX_DOFS:-512}
JOBS=${JOBS:-8}
THREADS=${THREADS:-8}
RESUME=${RESUME:-0}
CONTINUE_ON_ERROR=${CONTINUE_ON_ERROR:-0}

export OMP_NUM_THREADS="$THREADS"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target bench_helmholtz_k -j "$JOBS"

BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_k"
CSV="$RESULT_DIR/scan.csv"
METADATA="$RESULT_DIR/metadata.txt"

if [[ "$RESUME" != "1" || ! -s "$CSV" ]]; then
    "$BENCH" --csv-header > "$CSV"
fi

{
    echo "run_started=$(date --iso-8601=seconds)"
    echo "git_head=$(git -C "$ROOT_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "git_status_begin"
    git -C "$ROOT_DIR" status --short 2>/dev/null || true
    echo "git_status_end"
    echo "compiler=$(c++ --version | head -n 1)"
    echo "cmake=$(cmake --version | head -n 1)"
    echo "kernel=$(uname -a)"
    echo "threads=$THREADS"
    echo "jobs=$JOBS"
    echo "k_values=$K_VALUES"
    echo "fine_gap=$FINE_GAP"
    echo "kH_target=$KH_TARGET"
    echo "mode=$MODE"
    echo "stability_max_dofs=$STABILITY_MAX_DOFS"
    echo "omp_proc_bind=$OMP_PROC_BIND"
    echo "omp_places=$OMP_PLACES"
    echo
} >> "$METADATA"

has_result() {
    local target_k=$1
    awk -F, -v target="$target_k" 'NR > 1 && $1 == target { found=1 } END { exit !found }' "$CSV"
}

for k in $K_VALUES; do
    if [[ "$RESUME" == "1" ]] && has_result "$k"; then
        echo "=== skipping completed k=$k ==="
        continue
    fi

    log="$RESULT_DIR/k${k}.log"
    time_log="$RESULT_DIR/k${k}.time"
    row_file="$RESULT_DIR/k${k}.csv.tmp"
    echo "=== k=$k, kH<=$KH_TARGET, fine-gap=$FINE_GAP, mode=$MODE, threads=$THREADS ===" | tee "$log"

    command=(
        "$BENCH"
        --k="$k"
        --H=auto
        --h=auto
        --fine-gap="$FINE_GAP"
        --ell=auto
        --kH-target="$KH_TARGET"
        --mode="$MODE"
        --stability-max-dofs="$STABILITY_MAX_DOFS"
        --format=csv
    )
    printf 'command:' >> "$log"
    printf ' %q' "${command[@]}" >> "$log"
    printf '\n' >> "$log"

    if /usr/bin/time -v -o "$time_log" "${command[@]}" > "$row_file" 2>> "$log"; then
        row=$(<"$row_file")
        printf '%s\n' "$row" | tee -a "$CSV" "$log"
        cat "$time_log" >> "$log"
        rm -f "$row_file"
    else
        status=$?
        echo "FAILED k=$k exit=$status" | tee -a "$log"
        [[ -f "$time_log" ]] && cat "$time_log" >> "$log"
        rm -f "$row_file"
        if [[ "$CONTINUE_ON_ERROR" != "1" ]]; then
            exit "$status"
        fi
    fi
done

echo "CSV: $CSV"
echo "metadata: $METADATA"
