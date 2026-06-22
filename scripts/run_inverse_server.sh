#!/usr/bin/env bash
set -euo pipefail

# Server runner for inverse-inequality experiments.
# Modes:
#   MODE=h   : fixed H, ell; sweep h from H_MIN to H_MAX
#   MODE=ell : fixed H, h;   sweep ell from ELL_MIN to ELL_MAX
#   MODE=H   : fixed h, ell; sweep H from H_MIN to H_MAX
# Examples:
#   THREADS=8 MODE=h H=3 ELL=2 H_MIN=6 H_MAX=12 COEFF=unit SOLVER=auto bash scripts/run_inverse_server.sh
#   THREADS=8 MODE=ell H=3 H_FIXED=12 ELL_MIN=1 ELL_MAX=5 COEFF=unit SOLVER=auto bash scripts/run_inverse_server.sh
#   THREADS=8 MODE=H H_FIXED=12 ELL=3 H_MIN=2 H_MAX=4 COEFF=unit SOLVER=auto bash scripts/run_inverse_server.sh

THREADS="${THREADS:-8}"
MODE="${MODE:-h}"
H="${H:-3}"
ELL="${ELL:-2}"
H_MIN="${H_MIN:-6}"
H_MAX="${H_MAX:-12}"
H_FIXED="${H_FIXED:-${H_MAX}}"
ELL_MIN="${ELL_MIN:-1}"
ELL_MAX="${ELL_MAX:-4}"
COEFF="${COEFF:-unit}"
SOLVER="${SOLVER:-auto}"
SPACE="${SPACE:-free}"
BASIS="${BASIS:-lod}"
BUILD_DIR="${BUILD_DIR:-build}"
RESULT_DIR="${RESULT_DIR:-results/inverse_inequality}"

mkdir -p "${RESULT_DIR}"

cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --target bench_inverse_inequality -j "${THREADS}"

run_case() {
    local H_case="$1"
    local h="$2"
    local ell="$3"
    local log="$4"
    echo "=== Running H=${H_case} h=${h} ell=${ell} coeff=${COEFF} basis=${BASIS} solver=${SOLVER} threads=${THREADS} ==="
    set +e
    OMP_NUM_THREADS="${THREADS}" /usr/bin/time -v \
        "./${BUILD_DIR}/benchmarks/bench_inverse_inequality" \
        --H="${H_case}" --h="${h}" --ell="${ell}" \
        --basis="${BASIS}" --coeff="${COEFF}" \
        --solver="${SOLVER}" --threads="${THREADS}" --space="${SPACE}" \
        > "${log}" 2>&1
    local status=$?
    set -e
    cat "${log}"
    if [[ ${status} -ne 0 ]]; then
        echo "Run failed for H=${H_case}, h=${h}, ell=${ell}; continuing to preserve previous results." >&2
    fi
    return "${status}"
}

if [[ "${MODE}" == "h" ]]; then
    summary="${RESULT_DIR}/inverse_H${H}_ell${ELL}_${BASIS}_${COEFF}_h${H_MIN}-${H_MAX}.csv"
    echo "H,h,ell,basis,coeff,solver,space,threads,status,log" > "${summary}"
    for h in $(seq "${H_MIN}" "${H_MAX}"); do
        log="${RESULT_DIR}/inverse_H${H}_h${h}_ell${ELL}_${BASIS}_${COEFF}.log"
        set +e
        run_case "${H}" "${h}" "${ELL}" "${log}"
        status=$?
        set -e
        echo "${H},${h},${ELL},${BASIS},${COEFF},${SOLVER},${SPACE},${THREADS},${status},${log}" >> "${summary}"
    done
elif [[ "${MODE}" == "ell" ]]; then
    summary="${RESULT_DIR}/inverse_H${H}_h${H_FIXED}_${BASIS}_${COEFF}_ell${ELL_MIN}-${ELL_MAX}.csv"
    echo "H,h,ell,basis,coeff,solver,space,threads,status,log" > "${summary}"
    for ell in $(seq "${ELL_MIN}" "${ELL_MAX}"); do
        log="${RESULT_DIR}/inverse_H${H}_h${H_FIXED}_ell${ell}_${BASIS}_${COEFF}.log"
        set +e
        run_case "${H}" "${H_FIXED}" "${ell}" "${log}"
        status=$?
        set -e
        echo "${H},${H_FIXED},${ell},${BASIS},${COEFF},${SOLVER},${SPACE},${THREADS},${status},${log}" >> "${summary}"
    done
elif [[ "${MODE}" == "H" ]]; then
    summary="${RESULT_DIR}/inverse_h${H_FIXED}_ell${ELL}_${BASIS}_${COEFF}_H${H_MIN}-${H_MAX}.csv"
    echo "H,h,ell,basis,coeff,solver,space,threads,status,log" > "${summary}"
    for H_case in $(seq "${H_MIN}" "${H_MAX}"); do
        if (( H_case > H_FIXED )); then
            echo "Skipping H=${H_case} because H must be <= h=${H_FIXED}" >&2
            continue
        fi
        log="${RESULT_DIR}/inverse_H${H_case}_h${H_FIXED}_ell${ELL}_${BASIS}_${COEFF}.log"
        set +e
        run_case "${H_case}" "${H_FIXED}" "${ELL}" "${log}"
        status=$?
        set -e
        echo "${H_case},${H_FIXED},${ELL},${BASIS},${COEFF},${SOLVER},${SPACE},${THREADS},${status},${log}" >> "${summary}"
    done
else
    echo "Unknown MODE=${MODE}; use MODE=h, MODE=ell, or MODE=H" >&2
    exit 2
fi

echo "Summary: ${summary}"
