#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build-pollution-server"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_pollution_server"}
MODE=${MODE:-pilot}
BUILD_JOBS=${BUILD_JOBS:-32}
THREADS=${THREADS:-32}
PILOT_THREADS=${PILOT_THREADS:-"8 16 32 64"}
PILOT_K=${PILOT_K:-32}
K_VALUES=${K_VALUES:-"32 64 128"}
STRICT_K_VALUES=${STRICT_K_VALUES:-"32 64"}
KH_TARGET=${KH_TARGET:-1.0}
FINE_GAP=${FINE_GAP:-6}
STRICT_FINE_GAP=${STRICT_FINE_GAP:-8}
SYMBOLIC_CACHE_SLOTS=${SYMBOLIC_CACHE_SLOTS:-1}
FACTORIZATION_REUSE=${FACTORIZATION_REUSE:-none}
NUMA_POLICY=${NUMA_POLICY:-interleave}
MIN_AVAILABLE_GIB=${MIN_AVAILABLE_GIB:-256}
ENABLE_LTO=${ENABLE_LTO:-1}
RUN_FORCE=${RUN_FORCE:-0}

case "$FACTORIZATION_REUSE" in
  none|identical) ;;
  *) echo "FACTORIZATION_REUSE must be none or identical" >&2; exit 2 ;;
esac
case "$NUMA_POLICY" in
  default|interleave|local) ;;
  *) echo "NUMA_POLICY must be default, interleave, or local" >&2; exit 2 ;;
esac

mkdir -p "$RESULT_DIR"
cmake_args=(
  -S "$ROOT_DIR"
  -B "$BUILD_DIR"
  -DCMAKE_BUILD_TYPE=Release
  -DLOD2D_USE_OPENMP=ON
  -DLOD2D_BUILD_TESTS=ON
  -DLOD2D_BUILD_BENCHMARKS=ON
)
if [[ "$ENABLE_LTO" == 1 ]]; then
  cmake_args+=(-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON)
fi
cmake "${cmake_args[@]}"
cmake --build "$BUILD_DIR" --target bench_helmholtz_k -j "$BUILD_JOBS"
ctest --test-dir "$BUILD_DIR" \
  -R helmholtz_pollution_manufactured_smoke --output-on-failure

BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_k"
HEADER=$("$BENCH" --csv-header)

numa_prefix=()
if [[ "$NUMA_POLICY" != default ]]; then
  if ! command -v numactl >/dev/null 2>&1; then
    echo "numactl is required for NUMA_POLICY=$NUMA_POLICY" >&2
    exit 2
  fi
  if [[ "$NUMA_POLICY" == interleave ]]; then
    numa_prefix=(numactl --interleave=all)
  else
    numa_prefix=(numactl --localalloc)
  fi
fi

available_gib() {
  awk '/^MemAvailable:/ {printf "%d\n", $2 / 1024 / 1024}' /proc/meminfo
}

run_case() {
  local tag=$1
  local threads=$2
  local k=$3
  local gap=$4
  local cache_slots=$5
  local factor_reuse=$6
  local output="$RESULT_DIR/summary_${tag}.csv"
  local temporary="$output.tmp"
  local log="$RESULT_DIR/run_${tag}.log"
  local timing="$RESULT_DIR/time_${tag}.txt"
  local done_file="$RESULT_DIR/${tag}.done"

  if [[ "$RUN_FORCE" != 1 && -f "$done_file" ]]; then
    echo "skip completed case: $tag"
    return
  fi
  local free_gib
  free_gib=$(available_gib)
  if (( free_gib < MIN_AVAILABLE_GIB )); then
    echo "refusing $tag: MemAvailable=${free_gib} GiB < ${MIN_AVAILABLE_GIB} GiB" >&2
    exit 3
  fi

  local command=(
    "$BENCH"
    --k="$k"
    --H=auto
    --h=auto
    --fine-gap="$gap"
    --ell=auto
    --kH-target="$KH_TARGET"
    --source=manufactured
    --mode=two-sided
    --solver=saddle
    --threads="$threads"
    --symbolic-cache-slots="$cache_slots"
    --factorization-reuse="$factor_reuse"
    --stability-max-dofs=0
    --format=csv
    --check
  )
  : > "$log"
  if [[ ${#numa_prefix[@]} -gt 0 ]]; then
    printf '%q ' "${numa_prefix[@]}" >> "$log"
  fi
  printf '%q' "${command[0]}" >> "$log"
  printf ' %q' "${command[@]:1}" >> "$log"
  printf '\n' >> "$log"
  echo "start tag=$tag threads=$threads k=$k gap=$gap cache=$cache_slots reuse=$factor_reuse"
  {
    printf '%s\n' "$HEADER"
    OMP_NUM_THREADS="$threads" \
    OMP_DYNAMIC=FALSE \
    OMP_PROC_BIND=spread \
    OMP_PLACES=cores \
    OMP_MAX_ACTIVE_LEVELS=1 \
    OMP_STACKSIZE=32M \
    MALLOC_ARENA_MAX=8 \
      /usr/bin/time -v -o "$timing" \
      "${numa_prefix[@]}" "${command[@]}"
  } > "$temporary" 2>> "$log"
  if [[ "$(wc -l < "$temporary")" -ne 2 ]]; then
    echo "unexpected CSV row count for $tag" >&2
    exit 4
  fi
  mv "$temporary" "$output"
  touch "$done_file"
  echo "finish tag=$tag"
}

run_pilot() {
  for threads in $PILOT_THREADS; do
    run_case "pilot_t${threads}_k${PILOT_K}_c1_none" \
      "$threads" "$PILOT_K" "$FINE_GAP" 1 none
  done
  run_case "pilot_t${THREADS}_k${PILOT_K}_c4_none" \
    "$THREADS" "$PILOT_K" "$FINE_GAP" 4 none
  run_case "pilot_t${THREADS}_k${PILOT_K}_c4_identical" \
    "$THREADS" "$PILOT_K" "$FINE_GAP" 4 identical
  run_case "pilot_t${THREADS}_k${PILOT_K}_c8_identical" \
    "$THREADS" "$PILOT_K" "$FINE_GAP" 8 identical
}

run_main() {
  for k in $K_VALUES; do
    run_case "main_k${k}_gap${FINE_GAP}_t${THREADS}" \
      "$THREADS" "$k" "$FINE_GAP" \
      "$SYMBOLIC_CACHE_SLOTS" "$FACTORIZATION_REUSE"
  done
}

run_strict() {
  for k in $STRICT_K_VALUES; do
    run_case "strict_k${k}_gap${STRICT_FINE_GAP}_t${THREADS}" \
      "$THREADS" "$k" "$STRICT_FINE_GAP" \
      "$SYMBOLIC_CACHE_SLOTS" "$FACTORIZATION_REUSE"
  done
}

case "$MODE" in
  pilot) run_pilot ;;
  main) run_main ;;
  strict) run_strict ;;
  all)
    run_pilot
    run_main
    run_strict
    ;;
  *) echo "MODE must be pilot, main, strict, or all" >&2; exit 2 ;;
esac

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
  echo "mode=$MODE"
  echo "threads=$THREADS"
  echo "pilot_threads=$PILOT_THREADS"
  echo "k_values=$K_VALUES"
  echo "strict_k_values=$STRICT_K_VALUES"
  echo "kH_target=$KH_TARGET"
  echo "fine_gap=$FINE_GAP"
  echo "strict_fine_gap=$STRICT_FINE_GAP"
  echo "symbolic_cache_slots=$SYMBOLIC_CACHE_SLOTS"
  echo "factorization_reuse=$FACTORIZATION_REUSE"
  echo "numa_policy=$NUMA_POLICY"
  echo "minimum_available_GiB=$MIN_AVAILABLE_GIB"
  echo "enable_lto=$ENABLE_LTO"
  echo "lscpu_begin"
  lscpu
  echo "lscpu_end"
  echo "numa_begin"
  numactl --hardware 2>/dev/null || true
  echo "numa_end"
  echo "memory_begin"
  free -h
  echo "memory_end"
  echo "git_status_begin"
  git -C "$ROOT_DIR" status --short
  echo "git_status_end"
} > "$RESULT_DIR/server_metadata.txt"

{
  printf '%s\n' "$HEADER"
  shopt -s nullglob
  for file in "$RESULT_DIR"/summary_*.csv; do
    tail -n 1 "$file"
  done
} > "$RESULT_DIR/all_results.csv"

echo "Results written to $RESULT_DIR"
