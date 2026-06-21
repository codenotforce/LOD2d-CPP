#!/usr/bin/env bash
set -euo pipefail

# Server runner for the inverse-inequality h-refinement experiment.
# Defaults target the high-memory EPYC server. Override with environment vars:
#   THREADS=64 H=3 ELL=2 H_MIN=6 H_MAX=12 COEFF=unit SOLVER=auto ./scripts/run_inverse_server.sh

THREADS="${THREADS:-32}"
H="${H:-3}"
ELL="${ELL:-2}"
H_MIN="${H_MIN:-6}"
H_MAX="${H_MAX:-12}"
COEFF="${COEFF:-unit}"
SOLVER="${SOLVER:-auto}"
SPACE="${SPACE:-free}"
BASIS="${BASIS:-lod}"
BUILD_DIR="${BUILD_DIR:-build}"
RESULT_DIR="${RESULT_DIR:-results/inverse_inequality}"

mkdir -p "${RESULT_DIR}"

cmake -S . -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" --target bench_inverse_inequality -j "${THREADS}"

summary="${RESULT_DIR}/inverse_H${H}_ell${ELL}_${BASIS}_${COEFF}_h${H_MIN}-${H_MAX}.csv"
echo "H,h,ell,basis,coeff,solver,space,threads,status,log" > "${summary}"

for h in $(seq "${H_MIN}" "${H_MAX}"); do
    log="${RESULT_DIR}/inverse_H${H}_h${h}_ell${ELL}_${BASIS}_${COEFF}.log"
    echo "=== Running H=${H} h=${h} ell=${ELL} coeff=${COEFF} basis=${BASIS} solver=${SOLVER} threads=${THREADS} ==="
    set +e
    OMP_NUM_THREADS="${THREADS}" /usr/bin/time -v \
        "./${BUILD_DIR}/benchmarks/bench_inverse_inequality" \
        --H="${H}" --h="${h}" --ell="${ELL}" \
        --basis="${BASIS}" --coeff="${COEFF}" \
        --solver="${SOLVER}" --threads="${THREADS}" --space="${SPACE}" \
        > "${log}" 2>&1
    status=$?
    set -e
    cat "${log}"
    echo "${H},${h},${ELL},${BASIS},${COEFF},${SOLVER},${SPACE},${THREADS},${status},${log}" >> "${summary}"
    if [[ ${status} -ne 0 ]]; then
        echo "Run failed for h=${h}; continuing to preserve previous results." >&2
    fi
done

echo "Summary: ${summary}"
