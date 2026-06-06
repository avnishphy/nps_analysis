#!/bin/bash
# ============================================================================
# NPS π0 WF Analysis - Parallel Bash Orchestration
# Processes each run with src/nps_analysis_wfpi0.C in parallel and reconstructs
# master CSV from per-run CSV outputs to avoid in-macro append issues.
# ============================================================================

set -u

RUNLIST="${1:-config/runlist_x60_4b.txt}"
SKIM_DIR="${2:-/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b}"
OUTPUT_DIR="${3:-output/plots/x60_4b/}"
BEAM_ENERGY="${4:-10.538}"
TIMEOUT_SEC="${5:-1800}"
NPROC="${6:-$(nproc)}"
SESSION_TAG="${USER:-user}_$$"

if [[ ! -f "$RUNLIST" ]]; then
    echo "[ERROR] Runlist not found: $RUNLIST"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CSV_FILE="$OUTPUT_DIR/summary_all_runs.csv"
CSV_HEADER="run,accumulated_charge(mC),current_mean_uA,CPUT_LT,Beam_Time(s),total_entries,pass_hms,pass_hms_nps,total_coin_entries,estimated_time_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,hms_track_eff,hms_track_eff_err,s1x_peak,s1x_err,s1y_peak,s1y_err,s2x_peak,s2x_err,s2y_peak,s2y_err,run_status"

echo "$CSV_HEADER" > "$CSV_FILE"

echo "============================================================================"
echo "NPS π0 WF Analysis - Parallel Batch Processing"
echo "Start time: $(date)"
echo "============================================================================"
echo "Runlist:        $RUNLIST"
echo "Skim directory: $SKIM_DIR"
echo "Output dir:     $OUTPUT_DIR"
echo "Beam energy:    $BEAM_ENERGY GeV"
echo "Timeout/run:    $TIMEOUT_SEC s"
echo "Parallel jobs:  $NPROC"
echo "============================================================================"

# Function to process a single run
process_run() {
    local RUN="$1"
    local SKIM_DIR="$2"
    local OUTPUT_DIR="$3"
    local BEAM_ENERGY="$4"
    local TIMEOUT_SEC="$5"
    local TAG="$6"
    local TEMP_LIST="/tmp/runlist_wf_${RUN}_${TAG}.txt"
    local TEMP_OUT="/tmp/nps_wf_output_${RUN}_${TAG}.txt"
    local TEMP_RUN_CSV="/tmp/csv_wf_run_${RUN}_${TAG}.txt"

    echo "$RUN" > "$TEMP_LIST"
    export NPS_SKIM_DIR="$SKIM_DIR"
    export NPS_OUTPUT_DIR="$OUTPUT_DIR"
    export NPS_RUNLIST="$TEMP_LIST"
    export NPS_EBEAM="$BEAM_ENERGY"

    timeout "$TIMEOUT_SEC" root -l -b -q 'src/nps_analysis.C()' > "$TEMP_OUT" 2>&1
    local ROOT_EXIT_CODE=$?

    if grep -q "\[CSV_WRITTEN\]" "$TEMP_OUT" 2>/dev/null; then
        local PER_RUN_CSV=$(grep "\[CSV_WRITTEN\]" "$TEMP_OUT" | tail -1 | awk '{print $NF}')
        if [[ -f "$PER_RUN_CSV" ]]; then
            if [[ $ROOT_EXIT_CODE -ne 0 ]]; then
                sed -i 's/,OK$/,SEGFAULT/' "$PER_RUN_CSV"
            fi
            cp "$PER_RUN_CSV" "$TEMP_RUN_CSV"
            echo "$RUN OK" > "/tmp/wf_status_${RUN}_${TAG}.txt"
        else
            echo "$RUN CSV_FAIL" > "/tmp/wf_status_${RUN}_${TAG}.txt"
        fi
    else
        echo "$RUN NO_CSV" > "/tmp/wf_status_${RUN}_${TAG}.txt"
    fi
    rm -f "$TEMP_LIST" "$TEMP_OUT"
}

# Read runs and launch jobs in parallel
RUNS=()
while IFS= read -r RUN; do
    RUN=$(echo "$RUN" | sed 's/#.*//' | xargs)
    [[ -z "$RUN" ]] && continue
    RUNS+=("$RUN")
done < "$RUNLIST"

TOTAL_RUNS=${#RUNS[@]}
CURRENT_RUN=0

PIDS=()
for RUN in "${RUNS[@]}"; do
    ((CURRENT_RUN++))
    # Wait for a slot if too many jobs
    while (( $(jobs -rp | wc -l) >= NPROC )); do
        sleep 1
    done
    echo "[$CURRENT_RUN/$TOTAL_RUNS] Launching run $RUN..."
    process_run "$RUN" "$SKIM_DIR" "$OUTPUT_DIR" "$BEAM_ENERGY" "$TIMEOUT_SEC" "$SESSION_TAG" &
    PIDS+=("$!")
done

# Wait for all jobs to finish
wait

echo ""
echo "Rebuilding master CSV from successful per-run outputs..."
echo "$CSV_HEADER" > "$CSV_FILE"
for RUN in "${RUNS[@]}"; do
    TEMP_RUN_CSV="/tmp/csv_wf_run_${RUN}_${SESSION_TAG}.txt"
    if [[ -f "$TEMP_RUN_CSV" ]]; then
        tail -1 "$TEMP_RUN_CSV" >> "$CSV_FILE"
        rm -f "$TEMP_RUN_CSV"
    fi
    rm -f "/tmp/wf_status_${RUN}_${SESSION_TAG}.txt"
done
rm -f "/tmp/runlist_wf_*_${SESSION_TAG}.txt" "/tmp/nps_wf_output_*_${SESSION_TAG}.txt" "/tmp/csv_wf_run_*_${SESSION_TAG}.txt" "/tmp/wf_status_*_${SESSION_TAG}.txt" 2>/dev/null

echo "============================================================================"
echo "WF Parallel Analysis Complete!"
echo "End time: $(date)"
echo "============================================================================"
echo "Master CSV:  $CSV_FILE"
echo "============================================================================"
echo "CSV preview:"
if [[ -f "$CSV_FILE" ]]; then
    tail -$((TOTAL_RUNS + 1)) "$CSV_FILE"
fi
