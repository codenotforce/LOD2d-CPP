#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build-H-convergence-server"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_H_convergence_server"}
MODE=${MODE:-pilot}
BUILD_JOBS=${BUILD_JOBS:-32}
THREADS=${THREADS:-32}
PILOT_THREADS=${PILOT_THREADS:-"16 32 64"}
PILOT_SOLVERS=${PILOT_SOLVERS:-"saddle schur"}
K=${K:-32}
FINE_LEVEL=${FINE_LEVEL:-21}
H_LEVELS=${H_LEVELS:-"10 11 12 13 14 15"}
ELL=${ELL:-5}
PILOT_FINE_LEVEL=${PILOT_FINE_LEVEL:-18}
PILOT_H_LEVEL=${PILOT_H_LEVEL:-12}
PILOT_ELL=${PILOT_ELL:-5}
SOLVER=${SOLVER:-schur}
SYMBOLIC_CACHE_SLOTS=${SYMBOLIC_CACHE_SLOTS:-1}
FACTORIZATION_REUSE=${FACTORIZATION_REUSE:-none}
EXPORT_FIELDS=${EXPORT_FIELDS:-0}
FINE_REFERENCE_LEVEL=${FINE_REFERENCE_LEVEL:-15}
NUMA_POLICY=${NUMA_POLICY:-auto}
MIN_AVAILABLE_GIB=${MIN_AVAILABLE_GIB:-200}
ENABLE_LTO=${ENABLE_LTO:-1}
RUN_FORCE=${RUN_FORCE:-0}

case "$MODE" in
  pilot|main|all) ;;
  *) echo "MODE must be pilot, main, or all" >&2; exit 2 ;;
esac
case "$SOLVER" in
  saddle|schur) ;;
  *) echo "SOLVER must be saddle or schur" >&2; exit 2 ;;
esac
case "$FACTORIZATION_REUSE" in
  none|identical) ;;
  *) echo "FACTORIZATION_REUSE must be none or identical" >&2; exit 2 ;;
esac
case "$NUMA_POLICY" in
  auto|default|interleave|local) ;;
  *) echo "NUMA_POLICY must be auto, default, interleave, or local" >&2; exit 2 ;;
esac

mkdir -p "$RESULT_DIR" "$RESULT_DIR/snapshots"
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
cmake --build "$BUILD_DIR" --target bench_helmholtz_H_convergence -j "$BUILD_JOBS"
ctest --test-dir "$BUILD_DIR" \
  -R helmholtz_H_convergence_smoke --output-on-failure

BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_H_convergence"
HEADER=$("$BENCH" --csv-header)

numa_prefix=()
resolved_numa_policy=$NUMA_POLICY
if [[ "$NUMA_POLICY" == auto ]]; then
  if command -v numactl >/dev/null 2>&1; then
    resolved_numa_policy=interleave
  else
    resolved_numa_policy=default
    echo "numactl not found: NUMA_POLICY=auto falls back to default"
  fi
fi
if [[ "$resolved_numa_policy" != default ]]; then
  if ! command -v numactl >/dev/null 2>&1; then
    echo "numactl is required for NUMA_POLICY=$resolved_numa_policy" >&2
    exit 2
  fi
  if [[ "$resolved_numa_policy" == interleave ]]; then
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
  local solver=$3
  local k=$4
  local fine_level=$5
  local H_level=$6
  local ell=$7
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
    --h="$fine_level"
    --H-levels="$H_level"
    --ell="$ell"
    --mode=two-sided
    --solver="$solver"
    --threads="$threads"
    --symbolic-cache-slots="$SYMBOLIC_CACHE_SLOTS"
    --factorization-reuse="$FACTORIZATION_REUSE"
    --export-dir="$RESULT_DIR/snapshots"
    --summary-out="$temporary"
    --check
  )
  if [[ "$EXPORT_FIELDS" == 1 ]]; then
    command+=(--export-fields)
  fi
  if [[ "$H_level" == "$FINE_REFERENCE_LEVEL" ]]; then
    command+=(--fine-reference)
  fi

  : > "$log"
  printf '%q ' "${numa_prefix[@]}" >> "$log"
  printf '%q' "${command[0]}" >> "$log"
  printf ' %q' "${command[@]:1}" >> "$log"
  printf '\n' >> "$log"
  echo "start tag=$tag threads=$threads solver=$solver k=$k H=$H_level h=$fine_level ell=$ell"
  OMP_NUM_THREADS="$threads" \
  OMP_DYNAMIC=FALSE \
  OMP_PROC_BIND=spread \
  OMP_PLACES=cores \
  OMP_MAX_ACTIVE_LEVELS=1 \
  OMP_STACKSIZE=32M \
  MALLOC_ARENA_MAX=8 \
    /usr/bin/time -v -o "$timing" \
    "${numa_prefix[@]}" "${command[@]}" >> "$log" 2>&1
  if [[ "$(wc -l < "$temporary")" -ne 2 ]]; then
    echo "unexpected CSV row count for $tag" >&2
    exit 4
  fi
  mv "$temporary" "$output"
  touch "$done_file"
  echo "finish tag=$tag"
}

run_pilot() {
  for solver in $PILOT_SOLVERS; do
    for threads in $PILOT_THREADS; do
      run_case \
        "pilot_${solver}_t${threads}_k${K}_H${PILOT_H_LEVEL}_h${PILOT_FINE_LEVEL}" \
        "$threads" "$solver" "$K" "$PILOT_FINE_LEVEL" \
        "$PILOT_H_LEVEL" "$PILOT_ELL"
    done
  done
}

run_main() {
  for H_level in $H_LEVELS; do
    run_case \
      "main_${SOLVER}_t${THREADS}_k${K}_H${H_level}_h${FINE_LEVEL}" \
      "$THREADS" "$SOLVER" "$K" "$FINE_LEVEL" "$H_level" "$ELL"
  done
}

case "$MODE" in
  pilot) run_pilot ;;
  main) run_main ;;
  all)
    run_pilot
    run_main
    ;;
esac

if [[ "$MODE" == main || "$MODE" == all ]]; then
  {
    printf '%s\n' "$HEADER"
    for H_level in $H_LEVELS; do
      tail -n 1 \
        "$RESULT_DIR/summary_main_${SOLVER}_t${THREADS}_k${K}_H${H_level}_h${FINE_LEVEL}.csv"
    done
  } | awk -F, -v OFS=, '
    BEGIN {
      OFMT = "%.17g"
      CONVFMT = "%.17g"
    }
    NR == 1 {
      for (i = 1; i <= NF; ++i) column[$i] = i
      print
      next
    }
    {
      H = $(column["H_max"])
      p1e = $(column["p1_energy_abs"])
      lode = $(column["lod_energy_abs"])
      p1l2 = $(column["p1_l2_abs"])
      lodl2 = $(column["lod_l2_abs"])
      if (have_previous) {
        denominator = log(previous_H / H)
        $(column["p1_energy_rate"]) = log(previous_p1e / p1e) / denominator
        $(column["lod_energy_rate"]) = log(previous_lode / lode) / denominator
        $(column["p1_l2_rate"]) = log(previous_p1l2 / p1l2) / denominator
        $(column["lod_l2_rate"]) = log(previous_lodl2 / lodl2) / denominator
      }
      print
      previous_H = H
      previous_p1e = p1e
      previous_lode = lode
      previous_p1l2 = p1l2
      previous_lodl2 = lodl2
      have_previous = 1
    }
  ' > "$RESULT_DIR/all_results.csv.tmp"
  mv "$RESULT_DIR/all_results.csv.tmp" "$RESULT_DIR/all_results.csv"
fi

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
  echo "mode=$MODE"
  echo "k=$K"
  echo "fine_level=$FINE_LEVEL"
  echo "H_levels=$H_LEVELS"
  echo "ell=$ELL"
  echo "threads=$THREADS"
  echo "solver=$SOLVER"
  echo "pilot_threads=$PILOT_THREADS"
  echo "pilot_solvers=$PILOT_SOLVERS"
  echo "symbolic_cache_slots=$SYMBOLIC_CACHE_SLOTS"
  echo "factorization_reuse=$FACTORIZATION_REUSE"
  echo "export_fields=$EXPORT_FIELDS"
  echo "fine_reference_level=$FINE_REFERENCE_LEVEL"
  echo "numa_policy_requested=$NUMA_POLICY"
  echo "numa_policy_resolved=$resolved_numa_policy"
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

echo "Results written to $RESULT_DIR"
