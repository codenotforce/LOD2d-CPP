#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_pollution"}
K_VALUES=${K_VALUES:-"4 8 16 32"}
KH_TARGETS=${KH_TARGETS:-"1.0"}
FINE_GAP=${FINE_GAP:-8}
MODE=${MODE:-two-sided}
THREADS=${THREADS:-8}
JOBS=${JOBS:-8}
SOLVER=${SOLVER:-saddle}
SYMBOLIC_CACHE_SLOTS=${SYMBOLIC_CACHE_SLOTS:-1}
FACTORIZATION_REUSE=${FACTORIZATION_REUSE:-none}
STABILITY_MAX_DOFS=${STABILITY_MAX_DOFS:-512}
RESUME=${RESUME:-1}

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target bench_helmholtz_k -j "$JOBS"

BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_k"
HEADER=$($BENCH --csv-header)

export OMP_NUM_THREADS="$THREADS"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
export OMP_PLACES=${OMP_PLACES:-cores}
export OMP_MAX_ACTIVE_LEVELS=1

for target in $KH_TARGETS; do
  target_tag=${target//./p}
  csv="$RESULT_DIR/scan_kH${target_tag}.csv"
  if [[ "$RESUME" != 1 || ! -s "$csv" ]]; then
    printf '%s\n' "$HEADER" > "$csv"
  elif [[ "$(head -n 1 "$csv")" != "$HEADER" ]]; then
    echo "CSV header mismatch in $csv; use a new RESULT_DIR or RESUME=0" >&2
    exit 2
  fi

  has_result() {
    local k=$1
    awk -F, -v target_k="$k" 'NR > 1 && $1 == target_k {found=1} END {exit !found}' "$csv"
  }

  for k in $K_VALUES; do
    if [[ "$RESUME" == 1 ]] && has_result "$k"; then
      echo "skip completed k=$k kH_target=$target"
      continue
    fi
    tag="k${k}_kH${target_tag}_gap${FINE_GAP}"
    row_file="$RESULT_DIR/${tag}.csv.tmp"
    log="$RESULT_DIR/${tag}.log"
    timing="$RESULT_DIR/${tag}.time"
    command=(
      "$BENCH"
      --k="$k"
      --H=auto
      --h=auto
      --fine-gap="$FINE_GAP"
      --ell=auto
      --kH-target="$target"
      --source=manufactured
      --mode="$MODE"
      --solver="$SOLVER"
      --threads="$THREADS"
      --symbolic-cache-slots="$SYMBOLIC_CACHE_SLOTS"
      --factorization-reuse="$FACTORIZATION_REUSE"
      --stability-max-dofs="$STABILITY_MAX_DOFS"
      --format=csv
      --check
    )
    printf '%q' "${command[0]}" > "$log"
    printf ' %q' "${command[@]:1}" >> "$log"
    printf '\n' >> "$log"
    echo "run k=$k kH_target=$target fine_gap=$FINE_GAP"
    /usr/bin/time -v -o "$timing" "${command[@]}" > "$row_file" 2>> "$log"
    row=$(<"$row_file")
    printf '%s\n' "$row" >> "$csv"
    rm -f "$row_file"
    touch "$RESULT_DIR/${tag}.done"
  done
done

{
  echo "date_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "git_commit=$(git -C "$ROOT_DIR" rev-parse HEAD)"
  echo "k_values=$K_VALUES"
  echo "kH_targets=$KH_TARGETS"
  echo "fine_gap=$FINE_GAP"
  echo "mode=$MODE"
  echo "solver=$SOLVER"
  echo "symbolic_cache_slots=$SYMBOLIC_CACHE_SLOTS"
  echo "factorization_reuse=$FACTORIZATION_REUSE"
  echo "threads=$THREADS"
  echo "omp_proc_bind=$OMP_PROC_BIND"
  echo "omp_places=$OMP_PLACES"
  echo "git_status_begin"
  git -C "$ROOT_DIR" status --short
  echo "git_status_end"
} > "$RESULT_DIR/metadata.txt"

echo "Results written to $RESULT_DIR"
