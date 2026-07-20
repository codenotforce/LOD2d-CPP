#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD_DIR=${BUILD_DIR:-"$ROOT_DIR/build"}
RESULT_DIR=${RESULT_DIR:-"$ROOT_DIR/results/helmholtz_schwarz_scale"}
CASES=${CASES:-"5:10:4 7:12:8 9:14:16"}
ELL_VALUES=${ELL_VALUES:-"3"}
SOURCE=${SOURCE:-manufactured}
SOLVER=${SOLVER:-hybrid}
BOUNDARY=${BOUNDARY:-dirichlet}
EXTENSION=${EXTENSION:-weighted}
IMPEDANCE_BETA=${IMPEDANCE_BETA:-1}
LOCAL_SOLVER=${LOCAL_SOLVER:-direct}
LOCAL_INVERSE=${LOCAL_INVERSE:-lu}
FACTORIZATION_REUSE=${FACTORIZATION_REUSE:-none}
THREADS=${THREADS:-8}
JOBS=${JOBS:-8}
RESTART=${RESTART:-100}
MAX_ITERS=${MAX_ITERS:-3000}
TOL=${TOL:-1e-10}
RESUME=${RESUME:-0}
CONTINUE_ON_ERROR=${CONTINUE_ON_ERROR:-0}

export OMP_NUM_THREADS="$THREADS"
export OMP_DYNAMIC=false
export OMP_PROC_BIND=${OMP_PROC_BIND:-close}
export OMP_PLACES=${OMP_PLACES:-cores}

mkdir -p "$RESULT_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target bench_helmholtz_two_level_schwarz -j "$JOBS"

BENCH="$BUILD_DIR/benchmarks/bench_helmholtz_two_level_schwarz"
SUMMARY="$RESULT_DIR/summary.csv"
METADATA="$RESULT_DIR/metadata.txt"

SUMMARY_HEADER="source,solver,boundary,extension,impedance_beta,local_solver,local_inverse,factorization_reuse,outer_tol,outer_restart,outer_max_iters,threads,H_level,h_level,ell,k,kH,coarse_dofs,fine_dofs,subdomains,min_local_dofs,max_local_dofs,min_owned_dofs,max_owned_dofs,local_solver_groups,reused_factorizations,max_reuse_group,build_ms,fine_lu_ms,lod_solve_ms,setup_ms,status,outer_iterations,outer_solve_ms,true_residual,petrov_residual,energy_rel,fine_rel,exact_energy,exact_l2,exact_energy_rel,exact_l2_rel,lod_petrov,lod_energy_rel,lod_exact_energy,lod_exact_l2,lod_exact_energy_rel,lod_exact_l2_rel,corrector_residual,constraint_residual,benchmark_total_ms,wall_seconds,max_rss_kb"

if [[ "$RESUME" == "1" && -s "$SUMMARY" ]]; then
    if [[ "$(head -n 1 "$SUMMARY")" != "$SUMMARY_HEADER" ]]; then
        echo "error: incompatible summary schema in $SUMMARY" >&2
        echo "use a new RESULT_DIR or rerun with RESUME=0" >&2
        exit 2
    fi
else
    echo "$SUMMARY_HEADER" > "$SUMMARY"
fi

{
    echo "run_started=$(date --iso-8601=seconds)"
    echo "git_head=$(git -C "$ROOT_DIR" rev-parse HEAD 2>/dev/null || echo unknown)"
    echo "compiler=$(c++ --version | head -n 1)"
    echo "cmake=$(cmake --version | head -n 1)"
    echo "kernel=$(uname -a)"
    echo "cases=$CASES"
    echo "ell_values=$ELL_VALUES"
    echo "source=$SOURCE"
    echo "solver=$SOLVER"
    echo "boundary=$BOUNDARY"
    echo "extension=$EXTENSION"
    echo "impedance_beta=$IMPEDANCE_BETA"
    echo "local_solver=$LOCAL_SOLVER"
    echo "local_inverse=$LOCAL_INVERSE"
    echo "factorization_reuse=$FACTORIZATION_REUSE"
    echo "threads=$THREADS"
    echo "outer_tol=$TOL"
    echo "outer_restart=$RESTART"
    echo "outer_max_iters=$MAX_ITERS"
    echo "omp_proc_bind=$OMP_PROC_BIND"
    echo "omp_places=$OMP_PLACES"
    echo
} >> "$METADATA"

has_result() {
    local H=$1
    local h=$2
    local ell=$3
    local k=$4
    awk -F, -v source="$SOURCE" -v solver="$SOLVER" \
        -v boundary="$BOUNDARY" -v extension="$EXTENSION" \
        -v impedance_beta="$IMPEDANCE_BETA" -v local_solver="$LOCAL_SOLVER" \
        -v local_inverse="$LOCAL_INVERSE" \
        -v factorization_reuse="$FACTORIZATION_REUSE" -v outer_tol="$TOL" \
        -v outer_restart="$RESTART" -v outer_max_iters="$MAX_ITERS" \
        -v threads="$THREADS" -v H="$H" -v h="$h" \
        -v ell="$ell" -v k="$k" \
        'NR > 1 && $1 == source && $2 == solver && $3 == boundary &&
             $4 == extension && $5 == impedance_beta &&
             $6 == local_solver && $7 == local_inverse &&
             $8 == factorization_reuse && $9 == outer_tol &&
             $10 == outer_restart && $11 == outer_max_iters &&
             $12 == threads && $13 == H && $14 == h && $15 == ell && $16 == k { found=1 }
         END { exit !found }' "$SUMMARY"
}

append_summary() {
    local log=$1
    local time_log=$2
    local H=$3
    local h=$4
    local ell=$5
    local k=$6
    local wall rss
    wall=$(awk -F= '$1 == "wall_seconds" { print $2 }' "$time_log")
    rss=$(awk -F= '$1 == "max_rss_kb" { print $2 }' "$time_log")
    awk -v source="$SOURCE" -v solver="$SOLVER" \
        -v boundary="$BOUNDARY" -v extension="$EXTENSION" \
        -v impedance_beta="$IMPEDANCE_BETA" -v local_solver="$LOCAL_SOLVER" \
        -v local_inverse="$LOCAL_INVERSE" \
        -v factorization_reuse="$FACTORIZATION_REUSE" -v outer_tol="$TOL" \
        -v outer_restart="$RESTART" -v outer_max_iters="$MAX_ITERS" \
        -v threads="$THREADS" -v H="$H" -v h="$h" \
        -v ell="$ell" -v k="$k" -v wall="$wall" -v rss="$rss" '
        function value(prefix,    i) {
            for (i=1; i<=NF; ++i)
                if (index($i, prefix) == 1)
                    return substr($i, length(prefix) + 1)
            return ""
        }
        /^coarse_dofs=/ {
            coarse=value("coarse_dofs="); fine=value("fine_dofs=")
            kH=value("kH="); build=value("build_ms=")
            fine_lu=value("fine_sparse_lu_ms=")
            lod_solve=value("lod_solve_ms=")
            setup=value("schwarz_setup_ms=")
        }
        /^LOD baseline/ {
            lod_petrov=value("petrov=")
            lod_energy=value("energy_rel=")
            lod_exact_E=value("lod_exact_E=")
            lod_exact_L2=value("lod_exact_L2=")
            lod_exact_E_rel=value("lod_exact_E_rel=")
            lod_exact_L2_rel=value("lod_exact_L2_rel=")
            corrector=value("corrector_residual=")
            constraint=value("constraint_residual=")
        }
        /^Schwarz subdomains=/ {
            subdomains=value("subdomains=")
            split(value("local_dofs="), local, /\.\./)
            min_local=local[1]; max_local=local[2]
            split(value("owned_dofs="), owned, /\.\./)
            min_owned=owned[1]; max_owned=owned[2]
            solver_groups=value("solver_groups=")
            reused=value("reused_factorizations=")
            max_reuse=value("max_reuse_group=")
        }
        /^Benchmark total_ms=/ {
            benchmark_total=value("total_ms=")
        }
        $1 == "hybrid" || $1 == "additive" || $1 == "local" || $1 == "identity" {
            status=value("status="); iterations=value("iterations=")
            solve=value("solve_ms="); residual=value("residual=")
            petrov=value("petrov="); energy=value("energy_rel=")
            fine_error=value("fine_rel="); exact_E=value("exact_E=")
            exact_L2=value("exact_L2=")
            exact_E_rel=value("exact_E_rel=")
            exact_L2_rel=value("exact_L2_rel=")
        }
        END {
            OFS=","
            print source,solver,boundary,extension,impedance_beta,local_solver,local_inverse,
                factorization_reuse,outer_tol,outer_restart,outer_max_iters,threads,
                H,h,ell,k,kH,coarse,fine,subdomains,min_local,max_local,
                min_owned,max_owned,solver_groups,reused,max_reuse,
                build,fine_lu,lod_solve,setup,status,iterations,solve,residual,
                petrov,energy,fine_error,exact_E,exact_L2,exact_E_rel,
                exact_L2_rel,lod_petrov,lod_energy,lod_exact_E,lod_exact_L2,
                lod_exact_E_rel,lod_exact_L2_rel,corrector,constraint,benchmark_total,wall,rss
        }' "$log" >> "$SUMMARY"
}

for ell in $ELL_VALUES; do
    for case_spec in $CASES; do
        IFS=: read -r H h k <<< "$case_spec"
        if [[ -z "$H" || -z "$h" || -z "$k" ]]; then
            echo "invalid CASES entry: $case_spec (expected H:h:k)" >&2
            exit 2
        fi
        if [[ "$RESUME" == "1" ]] && has_result "$H" "$h" "$ell" "$k"; then
            echo "=== skipping completed H=$H h=$h ell=$ell k=$k ==="
            continue
        fi

        stem="H${H}_h${h}_ell${ell}_k${k}_${SOLVER}_${BOUNDARY}_${EXTENSION}_b${IMPEDANCE_BETA}_${LOCAL_SOLVER}_${LOCAL_INVERSE}_${FACTORIZATION_REUSE}"
        log="$RESULT_DIR/$stem.log"
        time_log="$RESULT_DIR/$stem.time"
        echo "=== H=$H h=$h ell=$ell k=$k solver=$SOLVER threads=$THREADS ===" | tee "$log"

        command=(
            "$BENCH"
            --H="$H"
            --h="$h"
            --ell="$ell"
            --k="$k"
            --threads="$THREADS"
            --source="$SOURCE"
            --solver="$SOLVER"
            --boundary="$BOUNDARY"
            --extension="$EXTENSION"
            --impedance-beta="$IMPEDANCE_BETA"
            --local-solver="$LOCAL_SOLVER"
            --local-inverse="$LOCAL_INVERSE"
            --factorization-reuse="$FACTORIZATION_REUSE"
            --restart="$RESTART"
            --max-iters="$MAX_ITERS"
            --tol="$TOL"
        )
        printf 'command:' >> "$log"
        printf ' %q' "${command[@]}" >> "$log"
        printf '\n' >> "$log"

        if /usr/bin/time -f 'wall_seconds=%e\nmax_rss_kb=%M'             -o "$time_log" "${command[@]}" >> "$log" 2>&1; then
            cat "$time_log" >> "$log"
            append_summary "$log" "$time_log" "$H" "$h" "$ell" "$k"
        else
            status=$?
            echo "FAILED exit=$status" | tee -a "$log"
            [[ -f "$time_log" ]] && cat "$time_log" >> "$log"
            if [[ "$CONTINUE_ON_ERROR" != "1" ]]; then
                exit "$status"
            fi
        fi
    done
done

echo "summary: $SUMMARY"
echo "metadata: $METADATA"
