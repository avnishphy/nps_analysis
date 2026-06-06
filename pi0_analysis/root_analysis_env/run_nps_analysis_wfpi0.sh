#!/bin/bash
# ============================================================================
# NPS π0 WF Analysis - Bash Orchestration
# Processes each run with src/nps_analysis_wfpi0.C and reconstructs master CSV
# from per-run CSV outputs to avoid in-macro append issues.
# ============================================================================

set -u

RUNLIST="${1:-config/runlist_x60_4b.txt}"
SKIM_DIR="${2:-/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS}"
OUTPUT_DIR="${3:-output/plots/x60_4b/production_wfpi0}"
BEAM_ENERGY="${4:-10.538}"
TIMEOUT_SEC="${5:-1800}"

if [[ ! -f "$RUNLIST" ]]; then
    echo "[ERROR] Runlist not found: $RUNLIST"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

CSV_FILE="$OUTPUT_DIR/summary_all_runs.csv"
CSV_HEADER="run,accumulated_charge(mC),current_mean_uA,CPUT_LT,Beam_Time(s),total_entries,pass_hms,pass_hms_nps,total_coin_entries,estimated_time_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,hms_track_eff,hms_track_eff_err,s1x_peak,s1x_err,s1y_peak,s1y_err,s2x_peak,s2x_err,s2y_peak,s2y_err,run_status"

echo "$CSV_HEADER" > "$CSV_FILE"

echo "============================================================================"
echo "NPS π0 WF Analysis - Batch Processing"
echo "Start time: $(date)"
echo "============================================================================"
echo "Runlist:        $RUNLIST"
echo "Skim directory: $SKIM_DIR"
echo "Output dir:     $OUTPUT_DIR"
echo "Beam energy:    $BEAM_ENERGY GeV"
echo "Timeout/run:    $TIMEOUT_SEC s"
echo "============================================================================"

N_SUCCESS=0
N_FAILED=0
FAILED=()
SUCCESSFUL_RUNS=()

TOTAL_RUNS=$(grep -v '^#' "$RUNLIST" | grep -v '^\s*$' | wc -l)
CURRENT_RUN=0

while IFS= read -r RUN; do
    RUN=$(echo "$RUN" | sed 's/#.*//' | xargs)
    [[ -z "$RUN" ]] && continue

    ((CURRENT_RUN++))

    BAR_WIDTH=30
    FILLED=$((CURRENT_RUN * BAR_WIDTH / TOTAL_RUNS))
    BAR=$(printf '%0.s=' $(seq 1 $FILLED))
    SPACE=$(printf '%0.s ' $(seq 1 $((BAR_WIDTH - FILLED))))
    PCT=$((CURRENT_RUN * 100 / TOTAL_RUNS))

    echo ""
    echo "[$BAR$SPACE] $PCT% ($CURRENT_RUN/$TOTAL_RUNS) Processing run $RUN..."
    echo "============================================================================"

    TEMP_LIST="/tmp/runlist_wf_${RUN}_$$.txt"
    TEMP_OUT="/tmp/nps_wf_output_${RUN}_$$.txt"
    TEMP_RUN_CSV="/tmp/csv_wf_run_${RUN}_$$.txt"

    echo "$RUN" > "$TEMP_LIST"

    export NPS_SKIM_DIR="$SKIM_DIR"
    export NPS_OUTPUT_DIR="$OUTPUT_DIR"
    export NPS_RUNLIST="$TEMP_LIST"
    export NPS_EBEAM="$BEAM_ENERGY"

    timeout "$TIMEOUT_SEC" root -l -b -q 'src/nps_analysis_wfpi0.C()' > "$TEMP_OUT" 2>&1
    ROOT_EXIT_CODE=$?

    cat "$TEMP_OUT"

    if grep -q "\[CSV_WRITTEN\]" "$TEMP_OUT" 2>/dev/null; then
        PER_RUN_CSV=$(grep "\[CSV_WRITTEN\]" "$TEMP_OUT" | tail -1 | awk '{print $NF}')
        if [[ -f "$PER_RUN_CSV" ]]; then
            if [[ $ROOT_EXIT_CODE -ne 0 ]]; then
                sed -i 's/,OK$/,SEGFAULT/' "$PER_RUN_CSV"
                echo ""
                echo "  ⚠ Run $RUN exited non-zero after CSV write (status: SEGFAULT)"
            else
                echo ""
                echo "  ✓ Run $RUN completed successfully"
            fi
            cp "$PER_RUN_CSV" "$TEMP_RUN_CSV"
            SUCCESSFUL_RUNS+=("$RUN")
            ((N_SUCCESS++))
        else
            echo ""
            echo "  ✗ Run $RUN failed (per-run CSV file not found: $PER_RUN_CSV)"
            ((N_FAILED++))
            FAILED+=("$RUN")
        fi
    else
        echo ""
        echo "  ✗ Run $RUN failed (no CSV_WRITTEN signal)"
        ((N_FAILED++))
        FAILED+=("$RUN")
    fi

    echo "============================================================================"
    rm -f "$TEMP_LIST" "$TEMP_OUT"
done < "$RUNLIST"

echo ""
echo "Rebuilding master CSV from successful per-run outputs..."
echo "$CSV_HEADER" > "$CSV_FILE"
for RUN in "${SUCCESSFUL_RUNS[@]}"; do
    TEMP_RUN_CSV="/tmp/csv_wf_run_${RUN}_$$.txt"
    if [[ -f "$TEMP_RUN_CSV" ]]; then
        tail -1 "$TEMP_RUN_CSV" >> "$CSV_FILE"
        rm -f "$TEMP_RUN_CSV"
    fi
done

rm -f /tmp/runlist_wf_* /tmp/nps_wf_output_* /tmp/csv_wf_run_* 2>/dev/null

echo ""
echo "============================================================================"
echo "WF Analysis Complete!"
echo "End time: $(date)"
echo "============================================================================"
echo "Total runs:  $TOTAL_RUNS"
echo "Successful:  $N_SUCCESS"
echo "Failed:      $N_FAILED"
echo "Master CSV:  $CSV_FILE"
if [[ $N_SUCCESS -gt 0 ]]; then
    echo "Successful runs: ${SUCCESSFUL_RUNS[*]}"
fi
if [[ $N_FAILED -gt 0 ]]; then
    echo "Failed runs: ${FAILED[*]}"
fi
echo "============================================================================"

echo ""
echo "CSV preview:"
if [[ -f "$CSV_FILE" ]]; then
    tail -$((N_SUCCESS + 1)) "$CSV_FILE"
fi
