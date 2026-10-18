#!/usr/bin/env bash
# run_experiments.sh — builds the code, regenerates the instances and runs the
# complete experimental protocol of the article (1400 runs).
#
#   E1 main        : 5 instances x 4 algorithms x 50 seeds (T_max = 5000)     = 1000 runs
#                    algorithms: GA, GA+LS, GA+LS-W (ablation), GA+AMLS (proposed)
#   E2 equal_time  : V = 250, 500 x 4 algorithms x 50 seeds, every algorithm
#                    gets the same wall-clock budget = mean time GA needed for
#                    5000 iterations in E1 on this machine                    =  400 runs
#
# Usage:
#   ./run_experiments.sh                 full protocol (deletes ./results first)
#   RESUME=1 ./run_experiments.sh        continue an interrupted run (keeps results,
#                                        skips runs whose summary.csv already exists)
#   NUM_SEEDS=2 ./run_experiments.sh     quick smoke test
#   CXX=g++ ./run_experiments.sh         choose the compiler (default: clang++ if found)
#
# All runs are executed sequentially (one core) so that wall-clock times are not
# distorted by concurrent runs. On macOS the script keeps the computer awake
# with `caffeinate` while it is running.
# (Compatible with the bash 3.2 that ships with macOS.)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_PATH="$SCRIPT_DIR/$(basename "$0")"

# Keep macOS awake (no idle/system/disk sleep) for the duration of the script.
if [ "$(uname -s)" = "Darwin" ] && [ -z "${CAFFEINATED:-}" ] && command -v caffeinate >/dev/null 2>&1; then
    export CAFFEINATED=1
    exec caffeinate -dims /bin/bash "$SCRIPT_PATH" ${1+"$@"}
fi

set -euo pipefail
cd "$SCRIPT_DIR"

RESULTS_DIR="results"
LOG_FILE="$RESULTS_DIR/run.log"
DONE_FILE="$RESULTS_DIR/DONE.log"
INFO_FILE="$RESULTS_DIR/system_info.txt"
NUM_SEEDS="${NUM_SEEDS:-50}"
FIRST_SEED=1000
RESUME="${RESUME:-0}"
if [[ -z "${CXX:-}" ]]; then
    if command -v clang++ >/dev/null 2>&1; then CXX=clang++; else CXX=g++; fi
fi
CXXFLAGS="-std=c++17 -O3 -Wall"

# ── 0. Clean previous results ────────────────────────────────────────────────
if [[ "$RESUME" != "1" ]]; then
    rm -rf "$RESULTS_DIR"
fi
mkdir -p "$RESULTS_DIR"
touch "$LOG_FILE"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg" | tee -a "$LOG_FILE"
}
elapsed_since() {
    local diff=$(( $(date +%s) - $1 ))
    printf '%dh %dm %ds' $((diff/3600)) $(((diff%3600)/60)) $((diff%60))
}

GLOBAL_START=$(date +%s)

# ── 1. System information (for the "Implementation and reproducibility" section)
{
    echo "Date:            $(date '+%Y-%m-%d %H:%M:%S %Z')"
    echo "OS kernel:       $(uname -srm)"
    if [[ "$(uname -s)" == "Darwin" ]]; then
        echo "macOS:           $(sw_vers -productVersion 2>/dev/null) (build $(sw_vers -buildVersion 2>/dev/null))"
        echo "Model:           $(sysctl -n hw.model 2>/dev/null)"
        echo "CPU:             $(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
        echo "CPU cores:       $(sysctl -n hw.ncpu 2>/dev/null) logical;" \
             "performance: $(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null)," \
             "efficiency: $(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null)"
        echo "P-core L2 cache: $(( $(sysctl -n hw.perflevel0.l2cachesize 2>/dev/null || echo 0) / 1048576 )) MB"
        echo "Memory:          $(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 )) GB"
        echo "Power source:    $(pmset -g batt 2>/dev/null | head -1 | sed "s/.*'\(.*\)'.*/\1/")"
        echo "Low power mode:  $(pmset -g 2>/dev/null | awk '/lowpowermode/ {print $2}')"
    else
        echo "CPU:             $(grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 | xargs || lscpu 2>/dev/null | awk -F: '/Vendor ID/ {print $2}' | xargs)"
        echo "CPU cores:       $(nproc 2>/dev/null)"
        echo "Memory:          $(awk '/MemTotal/ {printf "%.1f GB", $2/1048576}' /proc/meminfo 2>/dev/null)"
        echo "OS:              $(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME")"
    fi
    echo "Compiler:        $($CXX --version 2>&1 | head -1)"
    echo "Compiler flags:  $CXXFLAGS"
    echo "Git commit:      $(git rev-parse --short HEAD 2>/dev/null || echo n/a)$(git diff --quiet 2>/dev/null || echo ' (with uncommitted changes)')"
    echo "Seeds:           $FIRST_SEED..$((FIRST_SEED + NUM_SEEDS - 1))"
} > "$INFO_FILE.tmp"
mv "$INFO_FILE.tmp" "$INFO_FILE"

log "=========================================="
log "Enhanced GA Experiments — START (RESUME=$RESUME)"
log "=========================================="
while IFS= read -r line; do log "$line"; done < "$INFO_FILE"
log "=========================================="

# ── 2. Build ──────────────────────────────────────────────────────────────────
log "Building executables with $CXX $CXXFLAGS ..."
$CXX $CXXFLAGS -o generate generate.cpp 2>&1 | tee -a "$LOG_FILE"
$CXX $CXXFLAGS -o schedule schedule.cpp 2>&1 | tee -a "$LOG_FILE"
log "Build complete."

# ── 3. Generate input data and verify it ─────────────────────────────────────
log "Generating input data (5 variants)..."
./generate 2>&1 | tee -a "$LOG_FILE"
EXPECTED_MD5="d40998d35cce1a180a5fcb9a7a80a567 a3f1daec82836fed0765eb642a6ab2e0 074ed355abe1fc61ce05491e05dcbfb0 5df672910d586e9b30f97b3e4def0cf9 e0157a237a79134eb49364b29224a8f0"
md5_of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
ACTUAL_MD5=""
for v in 0 1 2 3 4; do ACTUAL_MD5="$ACTUAL_MD5 $(md5_of data/variant$v/requirements.csv)"; done
if [[ "${ACTUAL_MD5# }" == "$EXPECTED_MD5" ]]; then
    log "Input data verified: identical to the instances used in the article."
else
    log "WARNING: generated instances differ from the reference instances!"
fi

# ── 4. Experiment runner ──────────────────────────────────────────────────────
SEEDS=()
for ((s=0; s<NUM_SEEDS; s++)); do SEEDS+=($((FIRST_SEED + s))); done
MODE_NAMES=(no_ls ls ls_w amls)
FAILED=0
RUNS_DONE=0

# run_set <experiment-name> "<variants>" "<modes>" [extra schedule options...]
run_set() {
    local name="$1" variants="$2" modes="$3"; shift 3
    local extra=("$@")
    local out="$RESULTS_DIR/$name"
    for variant in $variants; do
        for mode in $modes; do
            local mstr="${MODE_NAMES[$mode]}"
            local t0; t0=$(date +%s)
            log "--- [$name] variant=$variant mode=$mstr (${#SEEDS[@]} seeds) ${extra[*]:-} ---"
            for seed in "${SEEDS[@]}"; do
                if [[ "$RESUME" == "1" && -f "$out/v${variant}_${mstr}/seed${seed}/summary.csv" ]]; then
                    continue
                fi
                if ./schedule "$variant" "$mode" "$seed" --out "$out" ${extra[@]+"${extra[@]}"} >> "$LOG_FILE" 2>&1; then
                    RUNS_DONE=$((RUNS_DONE + 1))
                else
                    log "ERROR: [$name] variant=$variant mode=$mstr seed=$seed FAILED"
                    FAILED=$((FAILED + 1))
                fi
            done
            log "[$name] variant=$variant mode=$mstr done in $(elapsed_since "$t0")  [total elapsed $(elapsed_since "$GLOBAL_START")]"
        done
    done
    # Concatenate the per-run summaries into one CSV.
    local first=1
    : > "$out/all_summaries.csv"
    for f in "$out"/v*_*/seed*/summary.csv; do
        if [[ $first -eq 1 ]]; then cat "$f" >> "$out/all_summaries.csv"; first=0
        else tail -n +2 "$f" >> "$out/all_summaries.csv"; fi
    done
}

# Mean total_time_sec of GA (mode no_ls) runs in E1 for a variant.
mean_ga_time() {
    awk -F, -v v="$1" 'NR==1 { for (i=1;i<=NF;i++) col[$i]=i; next }
        $col["variant"]==v && $col["mode"]=="no_ls" { s+=$col["total_time_sec"]; n++ }
        END { if (n>0) printf "%.3f", s/n; else print "0" }' "$RESULTS_DIR/main/all_summaries.csv"
}

# ── E1: main experiment ───────────────────────────────────────────────────────
log "========== E1: main experiment (iteration budget T_max = 5000) =========="
run_set main "0 1 2 3 4" "0 1 2 3"

# ── E2: equal wall-clock budget ──────────────────────────────────────────────
log "========== E2: equal wall-clock budget =========="
for v in 3 4; do
    B=$(mean_ga_time "$v")
    echo "variant=$v time_budget_sec=$B" >> "$RESULTS_DIR/equal_time_budgets.txt"
    log "Equal-time budget for variant $v: $B s (mean GA time for 5000 iterations)"
    run_set equal_time "$v" "0 1 2 3" --time-limit "$B" --max-iter 100000000
done

# ── 5. Completion report ──────────────────────────────────────────────────────
TOTAL_ELAPSED=$(elapsed_since "$GLOBAL_START")
log "=========================================="
log "All experiments finished: $RUNS_DONE runs executed, $FAILED failed"
log "Total wall time: $TOTAL_ELAPSED"
log "=========================================="

if [[ "$FAILED" -eq 0 ]]; then
    {
        echo "All experiments completed successfully."
        echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Total wall time: $TOTAL_ELAPSED"
        echo "Runs executed in this invocation: $RUNS_DONE"
        echo "Seeds: ${SEEDS[*]}"
        echo ""
        echo "Output: results/{main,equal_time}/v<k>_<mode>/seed<N>/{iterations,summary,schedule}.csv"
        echo "        results/{main,equal_time}/all_summaries.csv, results/system_info.txt"
    } | tee "$DONE_FILE"
    log "SUCCESS — see $DONE_FILE"
else
    log "WARNING — $FAILED run(s) failed. $DONE_FILE not written. Re-run with RESUME=1 to retry them."
    exit 1
fi
