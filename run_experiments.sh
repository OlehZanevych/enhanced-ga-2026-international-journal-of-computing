#!/usr/bin/env bash
# run_experiments.sh — builds, generates data, runs all GA experiments
# 5 variants × 2 modes × 50 seeds = 500 runs
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

RESULTS_DIR="results"
LOG_FILE="$RESULTS_DIR/run.log"
DONE_FILE="$RESULTS_DIR/DONE.log"
NUM_SEEDS=50

# ── 0. Clean previous results ────────────────────────────────────────────────
rm -rf "$RESULTS_DIR"
mkdir -p "$RESULTS_DIR"
> "$LOG_FILE"

log() {
    local msg="[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "$msg" | tee -a "$LOG_FILE"
}

elapsed_since() {
    local start=$1
    local now
    now=$(date +%s)
    local diff=$((now - start))
    printf '%dh %dm %ds' $((diff/3600)) $(((diff%3600)/60)) $((diff%60))
}

GLOBAL_START=$(date +%s)

log "=========================================="
log "Enhanced GA Experiments — START"
log "=========================================="
log "Working directory: $SCRIPT_DIR"
log "Seeds per config: $NUM_SEEDS"
log "Total runs: $((5 * 2 * NUM_SEEDS))"
log "Date: $(date '+%Y-%m-%d %H:%M:%S %Z')"
log "Host: $(hostname)"
log "OS: $(uname -srm)"
log "CPU: $(sysctl -n machdep.cpu.brand_string 2>/dev/null || grep -m1 'model name' /proc/cpuinfo 2>/dev/null | cut -d: -f2 || echo 'unknown')"
log "Compiler: $(clang++ --version 2>&1 | head -1)"
log "=========================================="

# ── 1. Build ──────────────────────────────────────────────────────────────────
log "Building executables..."
clang++ -std=c++17 -O3 -Wall -o generate generate.cpp 2>&1 | tee -a "$LOG_FILE"
clang++ -std=c++17 -O3 -Wall -o schedule schedule.cpp 2>&1 | tee -a "$LOG_FILE"
log "Build complete."

# ── 2. Generate input data ────────────────────────────────────────────────────
log "Generating input data (5 variants)..."
./generate 2>&1 | tee -a "$LOG_FILE"
log "Input data generation complete."

# ── 3. Run experiments ────────────────────────────────────────────────────────
TOTAL=$((5 * 2 * NUM_SEEDS))
DONE_CNT=0
FAILED=0

SEEDS=()
for ((s=0; s<NUM_SEEDS; s++)); do
    SEEDS+=($((1000 + s)))
done

for variant in 0 1 2 3 4; do
    for mode in 0 1; do
        MODE_STR=$( [ "$mode" -eq 0 ] && echo "no_ls" || echo "ls" )
        CONFIG_START=$(date +%s)
        log "--- variant=$variant mode=$MODE_STR ($NUM_SEEDS seeds) ---"

        for seed in "${SEEDS[@]}"; do
            DONE_CNT=$((DONE_CNT + 1))

            if ./schedule "$variant" "$mode" "$seed" 2>&1 | tee -a "$LOG_FILE"; then
                : # success
            else
                log "ERROR: variant=$variant mode=$MODE_STR seed=$seed FAILED"
                FAILED=$((FAILED + 1))
            fi
        done

        CONFIG_ELAPSED=$(elapsed_since "$CONFIG_START")
        TOTAL_ELAPSED=$(elapsed_since "$GLOBAL_START")
        log "variant=$variant mode=$MODE_STR: $NUM_SEEDS seeds done in $CONFIG_ELAPSED  [$DONE_CNT/$TOTAL total, elapsed $TOTAL_ELAPSED]"
    done
done

# ── 4. Completion report ──────────────────────────────────────────────────────
TOTAL_ELAPSED=$(elapsed_since "$GLOBAL_START")
log "=========================================="
log "All experiments finished: $((DONE_CNT - FAILED)) succeeded, $FAILED failed out of $TOTAL"
log "Total wall time: $TOTAL_ELAPSED"
log "=========================================="

if [ "$FAILED" -eq 0 ]; then
    {
        echo "All $TOTAL experiments completed successfully."
        echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Total wall time: $TOTAL_ELAPSED"
        echo "Seeds: ${SEEDS[*]}"
        echo ""
        echo "Output structure: results/v<k>_<mode>/seed<N>/{iterations,summary,schedule}.csv"
    } | tee "$DONE_FILE"
    log "SUCCESS — see $DONE_FILE"
else
    log "WARNING — $FAILED experiment(s) failed. $DONE_FILE not written."
    exit 1
fi
