#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build-adaptive-paper-server"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_adaptive_paper_server"}
REFERENCE_CACHE_DIR=${REFERENCE_CACHE_DIR:-"$RESULT_DIR/reference-cache"}
MODE=${MODE:-pilot}
JOBS=${JOBS:-16}
PATCH_THREADS=${PATCH_THREADS:-16}
MIN_AVAILABLE_GIB=${MIN_AVAILABLE_GIB:-}
MIN_FREE_DISK_GIB=${MIN_FREE_DISK_GIB:-}
EXTERNAL_TIMEOUT_SECONDS=${EXTERNAL_TIMEOUT_SECONDS:-}
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
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-epoch2-level12-calibration-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-palod-k16-epoch2-level12-step6-calibration-v4.json
    )
    ;;
  r2a-krobust)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k2-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-afem-k2-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k4-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-afem-k4-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k8-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-afem-k8-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k16-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-afem-k16-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-palod-k32-krobust-server-level18-step8-v4.json
      experiments/helmholtz_adaptive_paper/configs/R2a-afem-k32-krobust-server-level18-step8-v4.json
    )
    ;;
  s-corner-wave-366g)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-palod-k16-366g-step18-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-afem-k16-16t-step28-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-slod-k16-366g-step10-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-ufem-k16-16t-level20-v4.json
    )
    ;;
  s-corner-wave-medium)
    DEFAULT_CONFIGS=(
#      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-ufem-k16-H6-level20-step14-v4.json
#      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-afem-k16-H6-level20-step14-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-slod-k16-H6-ell2-gap4-h20-step10-v4.json
      experiments/helmholtz_adaptive_paper/configs/S-corner-wave-palod-k16-H6-h12-to-h20-gap4-step10-v4.json
    )
    ;;
  e1-revised-factor)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H4-h12-gap9-step5-refresh-factor-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=96
    MODE_MIN_FREE_DISK_GIB=100
    MODE_TIMEOUT_SECONDS=11700
    ;;
  e1-revised-pilot)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H4-h12-gap9-step8-pilot-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=192
    MODE_MIN_FREE_DISK_GIB=100
    MODE_TIMEOUT_SECONDS=22500
    ;;
  e1-revised-main)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H4-h12-gap9-step12-main-v6.json
      experiments/helmholtz_adaptive_paper/configs/E1-R1-afem-k16-H2-level18-step40-v4.json
      experiments/helmholtz_adaptive_paper/configs/E1-R1-ufem-k16-H2-level18-step16-v4.json
      experiments/helmholtz_adaptive_paper/configs/E1-R1-hlod-fixed-k16-H2-h16-ell3-step15-v4.json
    )
    MODE_MIN_AVAILABLE_GIB=320
    MODE_MIN_FREE_DISK_GIB=200
    MODE_TIMEOUT_SECONDS=36900
    ;;
  e1-revised-unified-factor)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H4-h12-theta03-step6-factor-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=96
    MODE_MIN_FREE_DISK_GIB=100
    MODE_TIMEOUT_SECONDS=11700
    ;;
  e1-revised-unified-pilot)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H4-h12-theta03-step12-pilot-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=192
    MODE_MIN_FREE_DISK_GIB=100
    MODE_TIMEOUT_SECONDS=29700
    ;;
  e1-revised-unified-main)
    # Shortest comparisons run first. PALOD is deliberately last so the
    # long high-memory trajectory begins only after all inexpensive methods
    # have produced valid artifacts.
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E1-R1-afem-k16-H4-level18-theta03-step60-v4.json
      experiments/helmholtz_adaptive_paper/configs/E1-R1-ufem-k16-H4-level18-step14-v4.json
      experiments/helmholtz_adaptive_paper/configs/E1-R1-hlod-fixed-k16-H4-h18-ell3-theta03-step18-v4.json
      experiments/helmholtz_adaptive_paper/configs/E1-R1-palod-reference-epoch-k16-H4-h12-theta03-step20-main-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=320
    MODE_MIN_FREE_DISK_GIB=200
    MODE_TIMEOUT_SECONDS=58500
    ;;
  e2-revised-factor)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E2-S-palod-moving-reference-k16-H3-h12-radius0125-theta03-step2-factor-v6.json
      experiments/helmholtz_adaptive_paper/configs/E2-S-palod-standard-reference-epoch-k16-H3-h12-gap9-step6-refresh-factor-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=96
    MODE_MIN_FREE_DISK_GIB=100
    MODE_TIMEOUT_SECONDS=11700
    ;;
  e2-revised-pilot)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E2-S-palod-moving-reference-k16-H3-h12-radius0125-theta03-step8-pilot-v6.json
      experiments/helmholtz_adaptive_paper/configs/E2-S-palod-standard-reference-epoch-k16-H3-h12-gap9-step8-pilot-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=192
    MODE_MIN_FREE_DISK_GIB=100
    MODE_TIMEOUT_SECONDS=22500
    ;;
  e2-revised-main)
    DEFAULT_CONFIGS=(
      experiments/helmholtz_adaptive_paper/configs/E2-S-afem-k16-H3-level20-step40-v4.json
      experiments/helmholtz_adaptive_paper/configs/E2-S-palod-standard-reference-epoch-k16-H3-h12-gap9-step12-main-v6.json
      experiments/helmholtz_adaptive_paper/configs/E2-S-hlod-fixed-k16-H3-h16-ell3-step12-v4.json
      experiments/helmholtz_adaptive_paper/configs/E2-S-palod-moving-reference-k16-H3-h12-radius0125-theta03-step15-main-v6.json
    )
    MODE_MIN_AVAILABLE_GIB=320
    MODE_MIN_FREE_DISK_GIB=200
    # The AFEM and fixed-LOD comparison configs retain a 24-hour internal
    # guard.  Leave 15 minutes for their structured shutdown/artifact write.
    MODE_TIMEOUT_SECONDS=87300
    ;;
  custom)
    DEFAULT_CONFIGS=()
    ;;
  *)
    echo "unknown MODE=$MODE; use a documented *-revised-factor, *-revised-pilot, *-revised-main, legacy mode, or custom" >&2
    exit 2
    ;;
esac

MIN_AVAILABLE_GIB=${MIN_AVAILABLE_GIB:-${MODE_MIN_AVAILABLE_GIB:-16}}
MIN_FREE_DISK_GIB=${MIN_FREE_DISK_GIB:-${MODE_MIN_FREE_DISK_GIB:-20}}
EXTERNAL_TIMEOUT_SECONDS=${EXTERNAL_TIMEOUT_SECONDS:-${MODE_TIMEOUT_SECONDS:-86400}}

if ! [[ "$PATCH_THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "PATCH_THREADS must be a positive integer" >&2
  exit 2
fi

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

require_disk_gate() {
  local free_gib
  free_gib=$(df -Pk "$RESULT_DIR" | awk 'NR==2 {printf "%d\n", $4 / 1024 / 1024}')
  if (( free_gib < MIN_FREE_DISK_GIB )); then
    echo "Refusing to start: free disk=${free_gib} GiB < ${MIN_FREE_DISK_GIB} GiB" >&2
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
  test_helmholtz_paper_config \
  test_helmholtz_paper_cases \
  test_helmholtz_reference_epoch_driver \
  test_helmholtz_reference_epoch_hierarchy \
  test_helmholtz_singularity_hybrid \
  bench_helmholtz_adaptive_paper

if [[ "$VALIDATE" == 1 ]]; then
  "$BUILD_DIR/tests/test_helmholtz_practical_driver"
  "$BUILD_DIR/tests/test_helmholtz_reference_retraction"
  "$BUILD_DIR/tests/test_helmholtz_reference_solution_cache"
  "$BUILD_DIR/tests/test_helmholtz_paper_config"
  "$BUILD_DIR/tests/test_helmholtz_paper_cases"
  "$BUILD_DIR/tests/test_helmholtz_reference_epoch_driver"
  "$BUILD_DIR/tests/test_helmholtz_reference_epoch_hierarchy"
  "$BUILD_DIR/tests/test_helmholtz_singularity_hybrid"
fi

BINARY="$BUILD_DIR/benchmarks/bench_helmholtz_adaptive_paper"
GIT_COMMIT=$(git -C "$ROOT_DIR" rev-parse HEAD)
BUILD_HASH="sha256:$(sha256sum "$BINARY" | awk '{print $1}')"
BASELINE_V6="$ROOT_DIR/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE.sha256"
BASELINE_V6_PRE_MOVING="$ROOT_DIR/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_PRE_MOVING_REFERENCE.sha256"
BASELINE_V4="$ROOT_DIR/experiments/helmholtz_adaptive_paper/MANUSCRIPT_BASELINE_LEGACY_V4.sha256"

IDENTITY_FILE="$RESULT_DIR/server-build-identity.txt"
if [[ -f "$IDENTITY_FILE" ]]; then
  recorded_commit=$(awk -F= '$1 == "git_commit" {print $2}' "$IDENTITY_FILE")
  recorded_build=$(awk -F= '$1 == "build_hash" {print $2}' "$IDENTITY_FILE")
  recorded_threads=$(awk -F= '$1 == "patch_threads" {print $2}' "$IDENTITY_FILE")
  if [[ "$recorded_commit" != "$GIT_COMMIT" \
        || "$recorded_build" != "$BUILD_HASH" \
        || "$recorded_threads" != "$PATCH_THREADS" ]]; then
    echo "Refusing to mix a different commit/build/thread count in RESULT_DIR=$RESULT_DIR" >&2
    echo "recorded: commit=$recorded_commit build=$recorded_build threads=$recorded_threads" >&2
    echo "current:  commit=$GIT_COMMIT build=$BUILD_HASH threads=$PATCH_THREADS" >&2
    exit 2
  fi
fi
printf 'git_commit=%s\nbuild_hash=%s\npatch_threads=%s\n' \
  "$GIT_COMMIT" "$BUILD_HASH" "$PATCH_THREADS" > "$IDENTITY_FILE"
lscpu > "$RESULT_DIR/lscpu.txt"
cat /proc/meminfo > "$RESULT_DIR/meminfo-before.txt"

export OMP_NUM_THREADS="$PATCH_THREADS"
export OMP_DYNAMIC=FALSE
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export LOD2D_PROGRESS=${LOD2D_PROGRESS:-1}
export LOD2D_PROFILE_LOCALIZATION_STAGES=${LOD2D_PROFILE_LOCALIZATION_STAGES:-1}
export LOD2D_PROFILE_MODEL_STAGES=${LOD2D_PROFILE_MODEL_STAGES:-1}

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
  require_disk_gate
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
  config_manuscript=$(python3 -c \
    'import json,sys; print(json.load(open(sys.argv[1]))["manuscript_sha256"].removeprefix("sha256:"))' \
    "$runtime_config")
  revised_manuscript=$(awk 'NF && $1 !~ /^#/ {print $1; exit}' "$BASELINE_V6")
  pre_moving_manuscript=$(awk 'NF && $1 !~ /^#/ {print $1; exit}' "$BASELINE_V6_PRE_MOVING")
  legacy_manuscript=$(awk 'NF && $1 !~ /^#/ {print $1; exit}' "$BASELINE_V4")
  if [[ "$config_manuscript" == "$revised_manuscript" ]]; then
    baseline="$BASELINE_V6"
  elif [[ "$config_manuscript" == "$pre_moving_manuscript" ]]; then
    baseline="$BASELINE_V6_PRE_MOVING"
  elif [[ "$config_manuscript" == "$legacy_manuscript" ]]; then
    baseline="$BASELINE_V4"
  else
    echo "Config manuscript hash matches neither frozen baseline: $config_manuscript" >&2
    exit 2
  fi
  set +e
  timeout --signal=INT --kill-after=600s "${EXTERNAL_TIMEOUT_SECONDS}s" \
    /usr/bin/time -v -o "$time_file" \
    stdbuf -oL -eL "$BINARY" \
      --config="$runtime_config" \
      --output-dir="$case_dir" \
      --reference-cache-dir="$REFERENCE_CACHE_DIR" \
      --manuscript-baseline="$baseline" \
      > "$stdout_file" 2>&1
  status=$?
  set -e
  if (( status != 0 )); then
    echo "failed $stem with exit status $status; see $stdout_file" >&2
    exit "$status"
  fi

  python3 - "$case_dir" "$time_file" "$stdout_file" "$MODE" <<'PY'
import csv, json, math, pathlib, re, sys
root = pathlib.Path(sys.argv[1])
time_path = pathlib.Path(sys.argv[2])
stdout_path = pathlib.Path(sys.argv[3])
mode = sys.argv[4]
manifests = list(root.glob("*/run.json"))
if len(manifests) != 1:
    raise SystemExit(f"expected one run.json below {root}, found {len(manifests)}")
manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
schema = int(manifest["schema_version"])
revised_mode = "-revised-" in mode
if revised_mode:
    timing_text = time_path.read_text(encoding="utf-8", errors="replace")
    swap_match = re.search(r"^\s*Swaps:\s*(\d+)\s*$", timing_text, re.MULTILINE)
    if not swap_match:
        raise SystemExit(f"missing swap count in {time_path}")
    if int(swap_match.group(1)) != 0:
        raise SystemExit(f"scientific gate rejected swap use: {swap_match.group(1)}")
if schema == 6:
    if manifest["status"] == "Failed":
        raise SystemExit(
            f"reference-epoch run failed: reason={manifest.get('stop_reason', '<missing>')}")
    if manifest["status"] == "WorkLimitReached" and not any(
        manifest.get("stop_reason", "").startswith(reason)
        for reason in (
            "maximum_H_steps reached",
            "insufficient remaining H-step budget for a new reference epoch",
        )
    ):
        raise SystemExit(
            "reference-epoch run hit an unaccepted work limit: "
            f"reason={manifest.get('stop_reason', '<missing>')}")
    required = (
        "iterations.csv", "summary.csv", "epoch_history.csv",
        "mesh_manifest.csv", "corrector_work.csv", "hybrid_reserve.csv")
    for name in required:
        if not (manifests[0].parent / name).is_file():
            raise SystemExit(f"missing schema-v6 artifact: {name}")
    if revised_mode:
        if manifest.get("claim") != "implementation-study":
            raise SystemExit("revised run is missing claim=implementation-study")
        config = manifest["config"]
        method = config["method"]
        hybrid = bool(config.get("singularity_hybrid"))
        variant = manifest.get("algorithm_variant", {})
        expected_conformance = "implementation-study-variant"
        if variant.get("manuscript_conformance") != expected_conformance:
            raise SystemExit(
                "unexpected manuscript conformance: "
                f"{variant.get('manuscript_conformance')!r}")
        if variant.get("candidate_hierarchy_closure_order") != "post-estimator":
            raise SystemExit("candidate hierarchy closure is not post-estimator")
        if variant.get("coarse_admissibility") != \
                "not-enforced;pre-asymptotic-points-retained":
            raise SystemExit("coarse-admissibility study policy is undeclared")

        def number(value):
            if value in (None, "", "NA"):
                return None
            result = float(value)
            return result if math.isfinite(result) else None

        def integer(value):
            parsed = number(value)
            return None if parsed is None else int(parsed)

        def truth(value):
            return str(value).lower() == "true"

        with (manifests[0].parent / "iterations.csv").open(
                newline="", encoding="utf-8") as stream:
            iterations = list(csv.DictReader(stream))
        iterations.sort(key=lambda row: int(row["event_sequence"]))
        solved = [row for row in iterations
                  if row["action"] == "SolveAndEstimate"]
        with (manifests[0].parent / "mesh_manifest.csv").open(
                newline="", encoding="utf-8") as stream:
            mesh_rows = list(csv.DictReader(stream))
        maximum_steps = int(config["work_limits"]["maximum_H_steps"])
        if manifest["status"] == "WorkLimitReached" \
                and manifest.get("stop_reason", "").startswith(
                    "maximum_H_steps reached") \
                and len(solved) != maximum_steps:
            raise SystemExit(
                f"solve count {len(solved)} != maximum_H_steps {maximum_steps}")

        if config["case"] == "R1" and method.startswith("PALOD"):
            checkpoints = [
                row for row in mesh_rows
                if int(row["epoch"]) == 0
                and row["stage"] in {"epoch_start", "committed", "pre_switch"}
            ]
            grouped = {}
            for row in checkpoints:
                key = (int(row["H_step"]), int(row["iteration"]))
                grouped.setdefault(key, set()).add(row["mesh_role"])
            if any(roles != {"coarse", "reference", "candidate"}
                   for roles in grouped.values()):
                raise SystemExit("E1 epoch-0 mesh checkpoint is not a complete triplet")
            solved_steps = {
                int(row["H_step"]) for row in solved if int(row["epoch"]) == 0}
            captured_steps = {key[0] for key in grouped}
            if not solved_steps.issubset(captured_steps):
                raise SystemExit(
                    "E1 epoch-0 mesh triplets miss solved H-steps: "
                    f"{sorted(solved_steps - captured_steps)}")
            reference_rows = [
                row for row in checkpoints if row["mesh_role"] == "reference"]
            reference_versions = {
                row.get("reference_mesh_version", "") for row in reference_rows}
            reference_files = {row["filename"] for row in reference_rows}
            if len(reference_versions) != 1 or "" in reference_versions \
                    or len(reference_files) != 1:
                raise SystemExit(
                    "E1 reference mesh changed inside epoch 0 or lacks version audit")

        # Corrector accounting and cross-epoch ell inheritance are exact
        # state-machine invariants, not empirical convergence heuristics.
        for row in iterations:
            active = integer(row.get("active_correctors"))
            rebuilt = integer(row.get("rebuilt_correctors"))
            reused = integer(row.get("reused_correctors"))
            if active is not None and rebuilt is not None and reused is not None \
                    and active != rebuilt + reused:
                raise SystemExit("active_correctors != rebuilt + reused")
        begin_rows = {
            int(row["epoch"]): row for row in iterations
            if row["action"] == "BeginEpoch"}
        for epoch, begin in sorted(begin_rows.items()):
            if epoch == 0:
                continue
            prior = [row for row in iterations if int(row["epoch"]) == epoch - 1]
            if not prior:
                raise SystemExit(f"epoch {epoch} has no preceding journal rows")
            if integer(begin["ell"]) != integer(prior[-1]["ell"]):
                raise SystemExit(f"ell was not inherited into epoch {epoch}")

        if hybrid:
            if variant.get("name") != "moving-reference-singularity-aware-v1" \
                    or variant.get("reserve_trigger_scope") != "not-applicable" \
                    or variant.get("reserve_selection") != \
                        "not-applicable-immediate-promotion":
                raise SystemExit("moving-reference update policy is undeclared")
            radius = float(config["hybrid_minimum_physical_radius"])
            if abs(radius - 0.125) > 1e-14:
                raise SystemExit(f"unexpected E2 physical radius {radius}")
            theta_h = float(config["theta_H"])
            theta_c = float(config["theta_c"])
            for row in iterations:
                covered = number(row.get("hybrid_covered_physical_radius"))
                ell_s = integer(row.get("hybrid_ell_S"))
                if row["action"] == "SolveAndEstimate" \
                        and (covered is None or ell_s is None):
                    raise SystemExit(
                        "hybrid solve row is missing physical-radius diagnostics")
                if covered is not None and covered + 1e-14 < radius:
                    raise SystemExit("hybrid singular core misses R_star")
                regular = number(row.get("hybrid_regular_indicator_mass"))
                marked_h = number(row.get("hybrid_marked_H_indicator_mass"))
                if regular is not None and regular > 0.0:
                    if marked_h is None or marked_h + 1e-12 * regular \
                            < theta_h * regular \
                            or not truth(row.get("hybrid_full_regular_doerfler")):
                        raise SystemExit("hybrid coarse full-regular Doerfler gate failed")
                for suffix in ("F", "R"):
                    mass = number(row.get(f"indicator_mass_c_{suffix}"))
                    marked = number(row.get(f"marked_mass_c_{suffix}"))
                    if mass is not None and mass > 0.0 \
                            and (marked is None
                                 or marked + 1e-12 * mass < theta_c * mass):
                        raise SystemExit(
                            f"candidate {suffix} regional Doerfler gate failed")

            with (manifests[0].parent / "hybrid_reserve.csv").open(
                    newline="", encoding="utf-8") as stream:
                reserve_rows = list(csv.DictReader(stream))
            closures = [row for row in reserve_rows
                        if row["row_type"] == "moving_reference_closure"]
            expected_closures = max(0, len(solved) - 1)
            if len(closures) != expected_closures:
                raise SystemExit(
                    "moving-reference promotion count does not match nonterminal solves")
            for row in closures:
                if row["status"] != "achieved" \
                        or not truth(row["target_satisfied"]) \
                        or int(row["matching_spill"]) != 0 \
                        or int(row["requested_target_gap"]) != 0:
                    raise SystemExit("moving-reference matching closure gate failed")

            if any(row["action"] == "ComputeCandidateDual"
                   for row in iterations):
                raise SystemExit("moving-reference run performed candidate dual work")

            preflight = re.findall(
                r"\[hybrid-preflight\].*?max_patch_fine_elements=(\d+).*?guard=(\d+)",
                stdout_path.read_text(encoding="utf-8", errors="replace"))
            if not preflight:
                raise SystemExit("hybrid patch preflight diagnostics are missing")
            if any(int(cost) > int(limit) for cost, limit in preflight):
                raise SystemExit("hybrid corrector patch cap was exceeded")
            if not (manifests[0].parent /
                    "mesh_E2_final_hybrid_regions.vtu").is_file():
                raise SystemExit("hybrid final-region mesh artifact is missing")

        actions = [row["action"] for row in iterations]
        if mode.endswith("-factor") and method.startswith("PALOD") \
                and "RefreshReference" not in actions:
            raise SystemExit("factor did not exercise a reference refresh")

        if mode.endswith(("-pilot", "-main")) and method.startswith("PALOD"):
            by_epoch = {}
            for row in solved:
                by_epoch.setdefault(int(row["epoch"]), []).append(row)
            refresh_epochs = {
                int(row["epoch"]) + 1 for row in iterations
                if row["action"] == "RefreshReference"}
            missing_begin = refresh_epochs.difference(begin_rows)
            if missing_begin:
                raise SystemExit(
                    "refreshed epoch has no BeginEpoch row: "
                    f"{sorted(missing_begin)}")
            if config["case"] == "R1":
                candidates = [rows for epoch, rows in by_epoch.items()
                              if epoch > 0 and len(rows) >= 3]
            elif hybrid:
                candidates = [solved] if len(solved) >= 3 else []
            else:
                candidates = []
            if (config["case"] == "R1" or hybrid) and not candidates:
                raise SystemExit("pilot/main has no required three-point epoch segment")
            if candidates:
                errors = [number(row["relative_exact_energy"])
                          for row in candidates[0]]
                if any(value is None for value in errors) or errors[-1] >= errors[0]:
                    raise SystemExit("selected three-point epoch has no net exact-error decay")
    print(
        f"validated {manifest['run_id']}: state={manifest['status']} "
        f"reason={manifest.get('stop_reason', '')}")
    raise SystemExit(0)
if manifest["status"] != "success":
    raise SystemExit(
        f"run did not complete successfully: status={manifest['status']} "
        f"driver_state={manifest.get('driver_state', '<missing>')} "
        f"reason={manifest.get('stop_reason', '<missing>')}")
policy = manifest.get("config", {}).get("trajectory_policy")
if policy in {"fixed_work_horizon", "fixed_empirical_trajectory"} \
        and manifest.get("driver_state") != "TrajectoryComplete":
    raise SystemExit(
        "fixed trajectory did not reach TrajectoryComplete: "
        f"driver_state={manifest.get('driver_state', '<missing>')}")
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
(
  cd "$RESULT_DIR"
  find . -type f ! -path './SHA256SUMS' -print0 \
    | sort -z | xargs -0 sha256sum \
    > SHA256SUMS
)
echo "All selected adaptive-paper runs completed. Results: $RESULT_DIR"
