#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_local_inverse"}
THREADS=${THREADS:-8}
K=${K:-2}
INITIAL_H=${INITIAL_H:-3}
FINE_LEVEL=${FINE_LEVEL:-12}
STEPS=${STEPS:-6}
ELL=${ELL:-3}
MARKS=${MARKS:-"fraction single-chain boundary-chain"}
MARK_FRACTION=${MARK_FRACTION:-0.25}
NEIGHBOR_LAYERS=${NEIGHBOR_LAYERS:-0}

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target bench_helmholtz_local_inverse -j "$THREADS"

for mark in $MARKS; do
  summary="$RESULT_DIR/summary_${mark}.csv"
  elements="$RESULT_DIR/elements_${mark}.csv"
  mesh="$RESULT_DIR/mesh_${mark}.csv"
  log="$RESULT_DIR/run_${mark}.log"
  command=(
    "$BUILD_DIR/benchmarks/bench_helmholtz_local_inverse"
    "--k=$K"
    "--initial-H=$INITIAL_H"
    "--fine-level=$FINE_LEVEL"
    "--steps=$STEPS"
    "--ell=$ELL"
    "--mark=$mark"
    "--mark-fraction=$MARK_FRACTION"
    "--neighbor-layers=$NEIGHBOR_LAYERS"
    "--element-out=$elements"
    "--mesh-out=$mesh"
    --format=csv
    --check
  )
  printf '%q' "${command[0]}" > "$log"
  printf ' %q' "${command[@]:1}" >> "$log"
  printf '\n' >> "$log"
  OMP_NUM_THREADS=$THREADS /usr/bin/time -v -o "$RESULT_DIR/time_${mark}.txt" \
    "${command[@]}" > "$summary" 2>> "$log"
done

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
  echo "git_status_begin"
  git -C "$ROOT_DIR" status --short
  echo "git_status_end"
  echo "threads=$THREADS"
  echo "k=$K"
  echo "initial_H=$INITIAL_H"
  echo "fine_level=$FINE_LEVEL"
  echo "steps=$STEPS"
  echo "ell=$ELL"
  echo "marks=$MARKS"
  echo "mark_fraction=$MARK_FRACTION"
  echo "neighbor_layers=$NEIGHBOR_LAYERS"
} > "$RESULT_DIR/metadata.txt"

echo "Results written to $RESULT_DIR"
