#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build-adaptive-paper-server"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/R2a-reference-adequacy-v4"}
REFERENCE_CACHE_DIR=${REFERENCE_CACHE_DIR:-"$RESULT_DIR/reference-cache"}
SOURCE_ITERATIONS=${SOURCE_ITERATIONS:-"$ROOT_DIR/results/adaptive-paper-extended-calibration-v3/runs/R2a-palod-k16-extended-calibration-v3/R2a_PALOD_k16_r0_e66e8517f8ac5dad/iterations.csv"}
TEMPLATE=${TEMPLATE:-"$ROOT_DIR/experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-reference-audit-v4.json"}
RUNTIME_CONFIG="$RESULT_DIR/runtime-config.json"
JOBS=${JOBS:-16}
MIN_AVAILABLE_GIB=${MIN_AVAILABLE_GIB:-16}

if ! git -C "$ROOT_DIR" diff --quiet ||
   ! git -C "$ROOT_DIR" diff --cached --quiet; then
  echo "Refusing to run from a checkout with tracked modifications" >&2
  exit 2
fi
if [[ ! -f "$SOURCE_ITERATIONS" ]]; then
  echo "Missing source trajectory: $SOURCE_ITERATIONS" >&2
  exit 2
fi
if [[ ! -f "$TEMPLATE" ]]; then
  echo "Missing audit config template: $TEMPLATE" >&2
  exit 2
fi

available_gib() {
  awk '/MemAvailable:/ {printf "%d\n", $2 / 1024 / 1024}' /proc/meminfo
}

available=$(available_gib)
if (( available < MIN_AVAILABLE_GIB )); then
  echo "Refusing to start reference audit: MemAvailable=${available} GiB < ${MIN_AVAILABLE_GIB} GiB" >&2
  exit 3
fi

mkdir -p "$BUILD_DIR" "$RESULT_DIR" "$REFERENCE_CACHE_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOD2D_BUILD_TESTS=ON \
  -DLOD2D_BUILD_BENCHMARKS=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native"
cmake --build "$BUILD_DIR" -j "$JOBS" --target \
  bench_helmholtz_reference_adequacy

BINARY="$BUILD_DIR/benchmarks/bench_helmholtz_reference_adequacy"
GIT_COMMIT=$(git -C "$ROOT_DIR" rev-parse HEAD)
BUILD_HASH="sha256:$(sha256sum "$BINARY" | awk '{print $1}')"
BASELINE="$ROOT_DIR/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256"

python3 - "$TEMPLATE" "$RUNTIME_CONFIG" "$GIT_COMMIT" "$BUILD_HASH" <<'PY'
import json, pathlib, sys
source, destination, commit, build_hash = sys.argv[1:]
data = json.loads(pathlib.Path(source).read_text(encoding="utf-8"))
data["git_commit"] = commit
data["build_hash"] = build_hash
pathlib.Path(destination).write_text(
    json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n",
    encoding="utf-8")
PY

/usr/bin/time -v -o "$RESULT_DIR/reference-adequacy.time" \
  "$BINARY" \
    --config="$RUNTIME_CONFIG" \
    --source-iterations="$SOURCE_ITERATIONS" \
    --output-dir="$RESULT_DIR" \
    --reference-cache-dir="$REFERENCE_CACHE_DIR" \
    --manuscript-baseline="$BASELINE" \
    --check \
    > "$RESULT_DIR/reference-adequacy.stdout" 2>&1

python3 - "$RESULT_DIR/reference_adequacy.json" <<'PY'
import json, pathlib, sys
path = pathlib.Path(sys.argv[1])
data = json.loads(path.read_text(encoding="utf-8"))
print(
    f"{data['status']}: relative_difference={data['relative_reference_difference']:.6g}, "
    f"terminal_fraction={data['terminal_error_fraction']:.6g}")
PY

(
  cd "$RESULT_DIR"
  find . -type f ! -path './SHA256SUMS' -print0 \
    | sort -z | xargs -0 sha256sum > SHA256SUMS
)
