#!/usr/bin/env bash
set -u

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
build_dir="${1:-${repo_root}/build-release}"
output_dir="${2:-${build_dir}/e0-calibration}"
log_file="${output_dir}.log"
status_file="${output_dir}.exit-status"

mkdir -p "${output_dir}"
rm -f "${status_file}"
OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}" \
    "${build_dir}/benchmarks/bench_helmholtz_e0_calibration" \
    --output-dir="${output_dir}" --check >"${log_file}" 2>&1
status=$?
printf '%s\n' "${status}" >"${status_file}"
exit "${status}"
