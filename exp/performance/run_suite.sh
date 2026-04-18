#!/bin/bash
##############################################################################
# run_suite.sh — Benchmark runner
#
# Measures CPU cycles (via perf stat) for each operation.
# Runs (backend × structure) pairs in parallel.
# Writes one CSV per (backend, structure, operation):
#   results/<backend>_<structure>_<op>.csv  →  columns: n,cycles
#
# Usage:
#   ./run_suite.sh                          # all backends, all structures
#   ./run_suite.sh --backend fcombinat
#   ./run_suite.sh --structure cycles
#   ./run_suite.sh --operation count
#   ./run_suite.sh --jobs 8
##############################################################################

set -euo pipefail
export LC_NUMERIC=C

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/results"
STRUCTURES_DIR="$SCRIPT_DIR/structures"
BACKENDS_DIR="$SCRIPT_DIR/backends"

NUM_TRIALS=3
TIMEOUT_SECONDS=1200   # 20 minutes
JOBS=12

FILTER_BACKEND=""
FILTER_STRUCTURE=""
FILTER_OPERATION=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --backend)    FILTER_BACKEND="$2";    shift 2 ;;
        --structure)  FILTER_STRUCTURE="$2";  shift 2 ;;
        --operation)  FILTER_OPERATION="$2";  shift 2 ;;
        --trials)     NUM_TRIALS="$2";        shift 2 ;;
        --timeout)    TIMEOUT_SECONDS="$2";   shift 2 ;;
        --jobs)       JOBS="$2";              shift 2 ;;
        -h|--help)
            echo "Usage: $0 [--backend NAME] [--structure NAME] [--operation NAME]"
            echo "          [--trials N] [--timeout SECS] [--jobs N]"
            exit 0 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p "$RESULTS_DIR"

# ---------------------------------------------------------------------------
# perf_cycles: run a command and capture its CPU cycle count.
# Sets: CMD_EXIT (int), CMD_OUTPUT (string), CYCLES (int).
# ---------------------------------------------------------------------------
perf_cycles() {
    local tmpperf tmpout tmpstderr
    tmpperf=$(mktemp)
    tmpout=$(mktemp)
    tmpstderr=$(mktemp)

    perf stat -e cycles -x , -o "$tmpperf" -- "$@" >"$tmpout" 2>"$tmpstderr" \
        && CMD_EXIT=0 || CMD_EXIT=$?
    CMD_OUTPUT=$(cat "$tmpout")
    # perf -x , format: value,,event_name,...  (first non-comment line)
    CYCLES=$(grep -m1 '^[0-9]' "$tmpperf" | cut -d, -f1)
    CYCLES=${CYCLES:-0}

    rm -f "$tmpperf" "$tmpout" "$tmpstderr"
}

# ---------------------------------------------------------------------------
# median_of: median of newline-separated numbers from stdin.
# ---------------------------------------------------------------------------
median_of() {
    sort -g | awk '{a[NR]=$1} END { print a[int((NR+1)/2)] }'
}

# ---------------------------------------------------------------------------
# random_rank [0, count-1]
# ---------------------------------------------------------------------------
random_rank() {
    python3 -c "
import sys, random
sys.set_int_max_str_digits(0)
c = int('$1')
print(random.randint(0, c - 1) if c > 0 else 0)
"
}

# ---------------------------------------------------------------------------
# process_struct: measure all operations for one (backend, structure) pair.
# Intended to run as a background job — each invocation is isolated.
# ---------------------------------------------------------------------------
process_struct() {
    local backend_dir="$1"
    local struct_file="$2"

    local backend_name
    backend_name=$(basename "$backend_dir")
    local backend_run="$backend_dir/run.sh"

    if [ ! -x "$backend_run" ]; then
        echo "  SKIP $backend_name (run.sh missing or not executable)"
        return
    fi

    # Source backend info (sets BACKEND_DESC, OPERATIONS)
    local BACKEND_DESC="$backend_name" OPERATIONS=""
    if [ -f "$backend_dir/info.conf" ]; then
        source "$backend_dir/info.conf"
    fi

    # Source structure config (sets NAME, SPEC, LABELED, SYMBOL, N_VALUES)
    source "$struct_file"
    local struct_id
    struct_id=$(basename "$struct_file" .conf)

    local tag="[$backend_name/$struct_id]"

    # Decide which operations to run based on backend support and filter
    local has_count=0 has_unrank=0 has_rank=0 has_draw=0 has_boltzmann=0
    _op_enabled() {
        local op="$1"
        { [ -z "$FILTER_OPERATION" ] || [ "$FILTER_OPERATION" = "$op" ]; } \
            && [[ " $OPERATIONS " == *" $op "* ]]
    }
    _op_enabled count      && has_count=1
    _op_enabled unrank     && has_unrank=1
    _op_enabled rank       && has_rank=1
    _op_enabled draw       && has_draw=1
    _op_enabled boltzmann  && has_boltzmann=1

    # Initialise output CSV files
    local base="$RESULTS_DIR/${backend_name}_${struct_id}"
    [ $has_count -eq 1 ]     && printf 'n,cycles\n' > "${base}_count.csv"
    [ $has_unrank -eq 1 ]    && printf 'n,cycles\n' > "${base}_unrank.csv"
    [ $has_rank -eq 1 ]      && printf 'n,cycles\n' > "${base}_rank.csv"
    [ $has_draw -eq 1 ]      && printf 'n,cycles\n' > "${base}_draw.csv"
    [ $has_boltzmann -eq 1 ] && printf 'n,cycles\n' > "${base}_boltzmann.csv"

    # count and unrank/rank use independent fast-fail flags so a failed
    # unrank never discards valid count data for the same or larger n.
    local count_fast_fail=0 unrank_rank_fast_fail=0

    for n in $N_VALUES; do

        # --- Phase 1: count (independent of unrank/rank) ---
        local count_times=() count_val=""
        local count_failed=0

        if [ $has_count -eq 1 ] && [ $count_fast_fail -eq 0 ]; then
            for trial in $(seq 1 "$NUM_TRIALS"); do
                perf_cycles timeout "$TIMEOUT_SECONDS" \
                    "$backend_run" count "$SPEC" "$n" "$LABELED" "$SYMBOL"
                if [ "$CMD_EXIT" -ne 0 ]; then
                    echo "  $tag n=$n FAIL count (trial $trial)"
                    count_failed=1; count_fast_fail=1; break
                fi
                count_val=$(printf '%s' "$CMD_OUTPUT" | tr -d ' \r\n')
                if [ -z "$count_val" ] || [ "$count_val" = "0" ]; then
                    echo "  $tag n=$n SKIP (count=0)"
                    count_failed=1; count_fast_fail=1; break
                fi
                count_times+=("$CYCLES")
            done
        elif [ $count_fast_fail -eq 1 ]; then
            echo "  $tag n=$n SKIP count (fast-fail)"
        fi

        # --- Phase 2: unrank + rank (needs count_val; independent fast-fail) ---
        local unrank_times=() rank_times=()
        local unrank_rank_failed=0

        if [ $has_unrank -eq 1 ] && [ $unrank_rank_fast_fail -eq 0 ] \
                && [ -n "$count_val" ]; then
            for trial in $(seq 1 "$NUM_TRIALS"); do
                local rk
                rk=$(random_rank "$count_val")
                perf_cycles timeout "$TIMEOUT_SECONDS" \
                    "$backend_run" unrank "$SPEC" "$n" "$LABELED" "$SYMBOL" "$rk"
                if [ "$CMD_EXIT" -ne 0 ]; then
                    echo "  $tag n=$n FAIL unrank (trial $trial)"
                    unrank_rank_failed=1; unrank_rank_fast_fail=1; break
                fi
                unrank_times+=("$CYCLES")
                local obj="$CMD_OUTPUT"

                if [ $has_rank -eq 1 ]; then
                    perf_cycles timeout "$TIMEOUT_SECONDS" \
                        "$backend_run" rank "$SPEC" "$n" "$LABELED" "$SYMBOL" "$obj"
                    if [ "$CMD_EXIT" -ne 0 ]; then
                        echo "  $tag n=$n FAIL rank (trial $trial)"
                        unrank_rank_failed=1; unrank_rank_fast_fail=1; break
                    fi
                    rank_times+=("$CYCLES")
                fi
            done
        elif [ $unrank_rank_fast_fail -eq 1 ]; then
            echo "  $tag n=$n SKIP unrank/rank (fast-fail)"
        fi

        # Compute medians and append to per-operation CSVs
        local mc mu mr
        if [ ${#count_times[@]}  -gt 0 ]; then
            mc=$(printf '%s\n' "${count_times[@]}"  | median_of)
            printf '%d,%s\n' "$n" "$mc" >> "${base}_count.csv"
        fi
        if [ ${#unrank_times[@]} -gt 0 ]; then
            mu=$(printf '%s\n' "${unrank_times[@]}" | median_of)
            printf '%d,%s\n' "$n" "$mu" >> "${base}_unrank.csv"
        fi
        if [ ${#rank_times[@]}   -gt 0 ]; then
            mr=$(printf '%s\n' "${rank_times[@]}"   | median_of)
            printf '%d,%s\n' "$n" "$mr" >> "${base}_rank.csv"
        fi

        # Progress line
        printf '  %s n=%d' "$tag" "$n"
        [ ${#count_times[@]}  -gt 0 ] && printf '  count=%dK'  "$((mc / 1000))"
        [ ${#unrank_times[@]} -gt 0 ] && printf '  unrank=%dK' "$((mu / 1000))"
        [ ${#rank_times[@]}   -gt 0 ] && printf '  rank=%dK'   "$((mr / 1000))"
        printf '\n'
    done

    echo "  $tag done => ${base}_*.csv"
}

# ---------------------------------------------------------------------------
# Main: discover pairs, launch parallel jobs (cap at JOBS).
# ---------------------------------------------------------------------------
echo "=============================================="
echo "  Benchmark Suite (CPU cycles, parallel)"
echo "  Trials: $NUM_TRIALS  Timeout: ${TIMEOUT_SECONDS}s  Jobs: $JOBS"
echo "=============================================="

declare -a job_pids=()

for backend_dir in "$BACKENDS_DIR"/*/; do
    backend_name=$(basename "$backend_dir")
    [ -n "$FILTER_BACKEND" ]   && [ "$FILTER_BACKEND"   != "$backend_name" ] && continue
    [ ! -x "$backend_dir/run.sh" ] && continue

    for struct_file in "$STRUCTURES_DIR"/*.conf; do
        struct_id=$(basename "$struct_file" .conf)
        [ -n "$FILTER_STRUCTURE" ] && [ "$FILTER_STRUCTURE" != "$struct_id" ] && continue

        # Throttle: wait for the oldest job if at the limit
        if [ "${#job_pids[@]}" -ge "$JOBS" ]; then
            wait "${job_pids[0]}" || true
            job_pids=("${job_pids[@]:1}")
        fi

        process_struct "$backend_dir" "$struct_file" &
        job_pids+=($!)
    done
done

# Wait for all remaining jobs
for pid in "${job_pids[@]}"; do
    wait "$pid" || true
done

echo ""
echo "=============================================="
echo "  Done. Results in $RESULTS_DIR/"
echo "=============================================="
