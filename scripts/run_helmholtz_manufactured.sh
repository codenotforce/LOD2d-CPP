#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_manufactured"}
THREADS=${THREADS:-8}
K=${K:-4}
H_LEVELS=${H_LEVELS:-3,4,5,6}
FINE_LEVEL=${FINE_LEVEL:-10}
ELL_LEVELS=${ELL_LEVELS:-2,3}
FEM_LEVELS=${FEM_LEVELS:-6,8,10}

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target bench_helmholtz_manufactured -j "$THREADS"

COMMAND=(
  "$BUILD_DIR/benchmarks/bench_helmholtz_manufactured"
  "--k=$K"
  "--h=$FINE_LEVEL"
  "--H-levels=$H_LEVELS"
  "--ell-levels=$ELL_LEVELS"
  "--fem-levels=$FEM_LEVELS"
  --check
)

OMP_NUM_THREADS=$THREADS "${COMMAND[@]}" | tee "$RESULT_DIR/validation.txt"
OMP_NUM_THREADS=$THREADS "${COMMAND[@]}" --format=csv > "$RESULT_DIR/validation.csv"

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
  echo "git_status_begin"
  git -C "$ROOT_DIR" status --short
  echo "git_status_end"
  echo "threads=$THREADS"
  echo "k=$K"
  echo "H_levels=$H_LEVELS"
  echo "fine_level=$FINE_LEVEL"
  echo "ell_levels=$ELL_LEVELS"
  echo "fem_levels=$FEM_LEVELS"
} > "$RESULT_DIR/metadata.txt"

echo "Results written to $RESULT_DIR"
