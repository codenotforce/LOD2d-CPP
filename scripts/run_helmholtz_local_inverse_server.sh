#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build-server"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_local_inverse_server"}
MODE=${MODE:-all}
THREADS=${THREADS:-32}
BUILD_JOBS=${BUILD_JOBS:-64}
PILOT_THREADS=${PILOT_THREADS:-"8 16 32 64"}
INITIAL_H=${INITIAL_H:-3}
ELL=${ELL:-3}
K=${K:-2}
STEPS=${STEPS:-6}
FINE_LEVEL=${FINE_LEVEL:-16}
HSCAN_LEVELS=${HSCAN_LEVELS:-"14 15 16"}
HSCAN_STEPS=${HSCAN_STEPS:-5}
RUN_FORCE=${RUN_FORCE:-0}

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOD2D_USE_OPENMP=ON \
  -DLOD2D_BUILD_TESTS=ON \
  -DLOD2D_BUILD_BENCHMARKS=ON
cmake --build "$BUILD_DIR" --target bench_helmholtz_local_inverse -j "$BUILD_JOBS"
ctest --test-dir "$BUILD_DIR" -R helmholtz_local_inverse_smoke --output-on-failure

EXE="$BUILD_DIR/benchmarks/bench_helmholtz_local_inverse"

run_case() {
  local tag=$1
  local threads=$2
  shift 2
  local summary="$RESULT_DIR/summary_${tag}.csv"
  local log="$RESULT_DIR/run_${tag}.log"
  local timing="$RESULT_DIR/time_${tag}.txt"
  local done_file="$RESULT_DIR/${tag}.done"
  if [[ "$RUN_FORCE" != 1 && -f "$done_file" ]]; then
    echo "skip completed case: $tag"
    return
  fi
  local command=("$EXE" "$@" --format=csv --check)
  printf '%q' "${command[0]}" > "$log"
  printf ' %q' "${command[@]:1}" >> "$log"
  printf '\n' >> "$log"
  echo "start case=$tag threads=$threads"
  OMP_NUM_THREADS="$threads" \
  OMP_DYNAMIC=FALSE \
  OMP_PROC_BIND=spread \
  OMP_PLACES=cores \
  OMP_MAX_ACTIVE_LEVELS=1 \
  OMP_STACKSIZE=32M \
  MALLOC_ARENA_MAX=8 \
    /usr/bin/time -v -o "$timing" "${command[@]}" > "$summary" 2>> "$log"
  touch "$done_file"
  echo "finish case=$tag"
}

run_calibration() {
  run_case calibration_full_Lh12 "$THREADS" \
    --k="$K" --initial-H="$INITIAL_H" --fine-level=12 --steps=2 --ell="$ELL" \
    --mark=single-chain --only-final --basis=all --denominators=all
}

run_pilot() {
  for threads in $PILOT_THREADS; do
    run_case "pilot_t${threads}" "$threads" \
      --k="$K" --initial-H="$INITIAL_H" --fine-level=12 --steps=2 --ell="$ELL" \
      --mark=single-chain --only-final --basis=trial --denominators=matched
  done
}

run_hscan() {
  for level in $HSCAN_LEVELS; do
    run_case "hscan_Lh${level}" "$THREADS" \
      --k="$K" --initial-H="$INITIAL_H" --fine-level="$level" \
      --steps="$HSCAN_STEPS" --ell="$ELL" --mark=single-chain --only-final \
      --basis=trial --denominators=element-matched
  done
}

run_feedback() {
  local tag="feedback_patch${ELL}_Lh${FINE_LEVEL}_steps${STEPS}"
  run_case "$tag" "$THREADS" \
    --k="$K" --initial-H="$INITIAL_H" --fine-level="$FINE_LEVEL" \
    --steps="$STEPS" --ell="$ELL" --mark=argmax-patch --neighbor-layers=0 \
    --basis=trial --denominators=element-matched \
    --element-out="$RESULT_DIR/elements_${tag}.csv" \
    --mesh-out="$RESULT_DIR/mesh_${tag}.csv"
}

case "$MODE" in
  pilot) run_pilot ;;
  hscan) run_hscan ;;
  feedback) run_feedback ;;
  all)
    run_calibration
    run_pilot
    run_hscan
    run_feedback
    ;;
  *) echo "MODE must be pilot, hscan, feedback, or all" >&2; exit 2 ;;
esac

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
  echo "mode=$MODE"
  echo "threads=$THREADS"
  echo "pilot_threads=$PILOT_THREADS"
  echo "build_jobs=$BUILD_JOBS"
  echo "initial_H=$INITIAL_H"
  echo "ell=$ELL"
  echo "k=$K"
  echo "steps=$STEPS"
  echo "fine_level=$FINE_LEVEL"
  echo "hscan_levels=$HSCAN_LEVELS"
  echo "hscan_steps=$HSCAN_STEPS"
  echo "lscpu_begin"
  lscpu
  echo "lscpu_end"
  echo "memory_begin"
  free -h
  echo "memory_end"
  echo "git_status_begin"
  git -C "$ROOT_DIR" status --short
  echo "git_status_end"
} > "$RESULT_DIR/server_metadata.txt"

echo "Results written to $RESULT_DIR"
