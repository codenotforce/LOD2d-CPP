#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build-adaptive-paper-server"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_adaptive_paper_server"}
REFERENCE_CACHE_DIR=${REFERENCE_CACHE_DIR:-"$RESULT_DIR/reference-cache"}
MODE=${MODE:-pilot}
JOBS=${JOBS:-16}
PATCH_THREADS=${PATCH_THREADS:-4}
MIN_AVAILABLE_GIB=${MIN_AVAILABLE_GIB:-16}
VALIDATE=${VALIDATE:-1}

if ! git -C "$ROOT_DIR" diff --quiet ||
   ! git -C "$ROOT_DIR" diff --cached --quiet; then
  echo "Refusing to run from a checkout with tracked modifications" >&2
  exit 2
fi

case "$MODE" in
  smoke)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/R1-palod-smoke-v4.json
      experiments/helmholtz_adaptive_paper/configs/R1-hlod-fixed-smoke-v4.json
      experiments/helmholtz_adaptive_paper/configs/R1-slod-smoke-v4.json
      experiments/helmholtz_adaptive_paper/configs/R1-ufem-smoke-v4.json
      experiments/helmholtz_adaptive_paper/configs/R1-afem-smoke-v4.json
    )
    ;;
  pilot)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-resource-pilot-v2.json
      experiments/helmholtz_adaptive_paper/configs/S-palod-k16-resource-pilot-v2.json
    )
    ;;
  calibration)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-epoch1-calibration-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-palod-k16-step6-calibration-v4.json
    )
    ;;
  custom)
    DEFAULT_CONFIGS=()
    ;;
  *)
    echo "MODE must be smoke, pilot, calibration, or custom" >&2
    exit 2
    ;;
esac

if [[ -n ${CONFIGS:-} ]]; then
  read -r -a CONFIG_LIST <<< "$CONFIGS"
else
  CONFIG_LIST=("${DEFAULT_CONFIGS[@]}")
fi
if [[ ${#CONFIG_LIST[@]} -eq 0 ]]; then
  echo "No configs selected; set CONFIGS='path1.json path2.json'" >&2
  exit 2
fi

available_gib() {
  awk '/MemAvailable:/ {printf "%d\n", $2 / 1024 / 1024}' /proc/meminfo
}

require_memory_gate() {
  local available
  available=$(available_gib)
  if (( available < MIN_AVAILABLE_GIB )); then
    echo "Refusing to start: MemAvailable=${available} GiB < ${MIN_AVAILABLE_GIB} GiB" >&2
    exit 3
  fi
}

mkdir -p "$BUILD_DIR" "$RESULT_DIR" "$REFERENCE_CACHE_DIR" \
  "$RESULT_DIR/runtime-configs" "$RESULT_DIR/logs"

cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE=Release \
  -DLOD2D_BUILD_TESTS=ON \
  -DLOD2D_BUILD_BENCHMARKS=ON \
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -DNDEBUG -march=native"
cmake --build "$BUILD_DIR" -j "$JOBS" --target \
  test_helmholtz_practical_driver \
  test_helmholtz_reference_retraction \
  test_helmholtz_reference_solution_cache \
  bench_helmholtz_adaptive_paper

if [[ "$VALIDATE" == 1 ]]; then
  "$BUILD_DIR/tests/test_helmholtz_practical_driver"
  "$BUILD_DIR/tests/test_helmholtz_reference_retraction"
  "$BUILD_DIR/tests/test_helmholtz_reference_solution_cache"
fi

BINARY="$BUILD_DIR/benchmarks/bench_helmholtz_adaptive_paper"
GIT_COMMIT=$(git -C "$ROOT_DIR" rev-parse HEAD)
BUILD_HASH="sha256:$(sha256sum "$BINARY" | awk '{print $1}')"
BASELINE="$ROOT_DIR/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256"

printf 'git_commit=%s\nbuild_hash=%s\npatch_threads=%s\n' \
  "$GIT_COMMIT" "$BUILD_HASH" "$PATCH_THREADS" \
  > "$RESULT_DIR/server-build-identity.txt"
lscpu > "$RESULT_DIR/lscpu.txt"
cat /proc/meminfo > "$RESULT_DIR/meminfo-before.txt"

export OMP_NUM_THREADS="$PATCH_THREADS"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

for config_relative in "${CONFIG_LIST[@]}"; do
  template="$ROOT_DIR/$config_relative"
  if [[ ! -f "$template" ]]; then
    echo "Missing config template: $template" >&2
    exit 2
  fi
  stem=$(basename "$template" .json)
  runtime_config="$RESULT_DIR/runtime-configs/${stem}.json"
  case_dir="$RESULT_DIR/runs/$stem"
  done_file="$RESULT_DIR/runs/${stem}.done"
  time_file="$RESULT_DIR/logs/${stem}.time"
  stdout_file="$RESULT_DIR/logs/${stem}.stdout"
  if [[ -f "$done_file" ]]; then
    echo "skip completed $stem"
    continue
  fi
  require_memory_gate
  mkdir -p "$case_dir"

  python3 - "$template" "$runtime_config" "$GIT_COMMIT" "$BUILD_HASH" <<'PY'
import json, pathlib, sys
source, destination, commit, build_hash = sys.argv[1:]
data = json.loads(pathlib.Path(source).read_text(encoding="utf-8"))
data["git_commit"] = commit
data["build_hash"] = build_hash
pathlib.Path(destination).write_text(
    json.dumps(data, sort_keys=True, separators=(",", ":")) + "\n",
    encoding="utf-8")
PY

  echo "start $stem; MemAvailable=$(available_gib) GiB"
  set +e
  /usr/bin/time -v -o "$time_file" \
    "$BINARY" \
      --config="$runtime_config" \
      --output-dir="$case_dir" \
      --reference-cache-dir="$REFERENCE_CACHE_DIR" \
      --manuscript-baseline="$BASELINE" \
      > "$stdout_file" 2>&1
  status=$?
  set -e
  if (( status != 0 )); then
    echo "failed $stem with exit status $status; see $stdout_file" >&2
    exit "$status"
  fi

  python3 - "$case_dir" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
manifests = list(root.glob("*/run.json"))
if len(manifests) != 1:
    raise SystemExit(f"expected one run.json below {root}, found {len(manifests)}")
manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
if manifest["status"] in {"linear_algebra_failure", "unavailable"}:
    raise SystemExit(f"unacceptable terminal status: {manifest['status']}")
for name in ("iterations.csv", "summary.csv", "ell_history.csv", "final_mesh.vtu"):
    if not (manifests[0].parent / name).is_file():
        raise SystemExit(f"missing artifact: {name}")
print(
    f"validated {manifest['run_id']}: status={manifest['status']} "
    f"regime={manifest['convergence_diagnostic']['status']} "
    f"reference_cache_hit={manifest['timing']['reference_solution_cache_hit']}")
PY
  touch "$done_file"
  echo "done $stem"
done

cat /proc/meminfo > "$RESULT_DIR/meminfo-after.txt"
find "$RESULT_DIR" -type f ! -path "$RESULT_DIR/SHA256SUMS" -print0 \
  | sort -z | xargs -0 sha256sum \
  > "$RESULT_DIR/SHA256SUMS"
echo "All selected adaptive-paper runs completed. Results: $RESULT_DIR"
