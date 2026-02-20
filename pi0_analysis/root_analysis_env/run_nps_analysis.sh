#!/bin/bash
# ============================================================================
# NPS π0 Analysis - Bash Orchestration
# Processes each run with the original nps_analysis.C
# ============================================================================

# Clean up any gibberish files at startup
find . -maxdepth 1 -type f \( -size -1k -o -size +100k \) ! -name "*.sh" ! -name ".*" 2>/dev/null | while read f; do
    # Remove files with suspicious names (non-printable characters)
    if ! echo "$f" | grep -q '^\.\/[a-zA-Z0-9._-]*$'; then
        rm -f "$f" 2>/dev/null
    fi
done

RUNLIST="${1:-config/runlist_x60_4b.txt}"
SKIM_DIR="${2:-/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b}"
OUTPUT_DIR="${3:-output/plots/x60_4b}"
BEAM_ENERGY="${4:-10.538}"

mkdir -p "$OUTPUT_DIR"

# CSV file for accumulating results
CSV_FILE="$OUTPUT_DIR/summary_all_runs.csv"
CSV_HEADER="run,accumulated_charge(mC),current_mean_uA,CPUT_LT,Beam_Time(s),total_entries,pass_hms,pass_hms_nps,total_coin_entries,estimated_time_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,hms_track_eff,hms_track_eff_err,s1x_peak,s1x_err,s1y_peak,s1y_err,s2x_peak,s2x_err,s2y_peak,s2y_err,run_status"

# Start fresh
> "$CSV_FILE"
echo "$CSV_HEADER" > "$CSV_FILE"

echo "============================================================================"
echo "NPS π0 Analysis - Batch Processing"
echo "Start time: $(date)"
echo "============================================================================"
echo "Configuration:"
echo "  Runlist: $RUNLIST"
echo "  Skim directory: $SKIM_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Beam energy: $BEAM_ENERGY GeV"
echo "============================================================================"

N_SUCCESS=0
N_FAILED=0
FAILED=()
SUCCESSFUL_RUNS=()

# Count total runs first
TOTAL_RUNS=$(grep -v '^#' "$RUNLIST" | grep -v '^\s*$' | wc -l)
CURRENT_RUN=0

# Process each run
while IFS= read -r RUN; do
    # Skip comments and blanks
    RUN=$(echo "$RUN" | sed 's/#.*//' | xargs)
    [[ -z "$RUN" ]] && continue
    
    ((CURRENT_RUN++))
    
    # Progress indicator with bar
    BAR_WIDTH=30
    FILLED=$((CURRENT_RUN * BAR_WIDTH / TOTAL_RUNS))
    BAR=$(printf '%0.s=' $(seq 1 $FILLED))
    SPACE=$(printf '%0.s ' $(seq 1 $((BAR_WIDTH - FILLED))))
    PCT=$((CURRENT_RUN * 100 / TOTAL_RUNS))
    
    echo ""
    echo "[$BAR$SPACE] $PCT% ($CURRENT_RUN/$TOTAL_RUNS) Processing run $RUN..."
    echo "============================================================================"
    
    # Create temp runlist with just this run
    TEMP_LIST="/tmp/runlist_$$.txt"
    echo "$RUN" > "$TEMP_LIST"
    
    # Save temp CSV from nps_analysis.C
    TEMP_RUN_CSV="/tmp/csv_run_${RUN}_$$.txt"
    
    # Temp file to capture output
    TEMP_OUT="/tmp/nps_output_${RUN}_$$.txt"
    
    # Run nps_analysis.C - show all output (diagnostics + progress)
    # Pass arguments via environment variables to avoid shell quoting issues
    export NPS_SKIM_DIR="$SKIM_DIR"
    export NPS_OUTPUT_DIR="$OUTPUT_DIR"
    export NPS_RUNLIST="$TEMP_LIST"
    export NPS_EBEAM="$BEAM_ENERGY"
    timeout 900 root -l -b -q 'src/nps_analysis.C()' > "$TEMP_OUT" 2>&1
    ROOT_EXIT_CODE=$?  # Capture exit code to detect crashes
    
    # Display the output to the user (shows progress bars and diagnostics)
    cat "$TEMP_OUT"
    
    # Check if per-run CSV was written successfully
    if grep -q "\[CSV_WRITTEN\]" "$TEMP_OUT" 2>/dev/null; then
        # Find the per-run CSV file
        PER_RUN_CSV=$(grep "\[CSV_WRITTEN\]" "$TEMP_OUT" | tail -1 | awk '{print $NF}')
        if [[ -f "$PER_RUN_CSV" ]]; then
            # CSV was written, but check if ROOT crashed (nonzero exit code = crash/segfault)
            if [[ $ROOT_EXIT_CODE -ne 0 ]]; then
                # Segfault/crash detected - update the CSV run_status to SEGFAULT
                # Replace last column (OK) with SEGFAULT
                sed -i 's/,OK$/,SEGFAULT/' "$PER_RUN_CSV"
                echo ""
                echo "  ⚠ Run $RUN segfaulted after CSV write (status: SEGFAULT)"
            else
                echo ""
                echo "  ✓ Run $RUN completed successfully"
            fi
            # Save this run's data
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

# Rebuild the master CSV file with all successful runs
echo ""
echo "Rebuilding master CSV with all successful runs..."
echo "$CSV_HEADER" > "$CSV_FILE"
for RUN in "${SUCCESSFUL_RUNS[@]}"; do
    TEMP_RUN_CSV="/tmp/csv_run_${RUN}_$$.txt"
    if [[ -f "$TEMP_RUN_CSV" ]]; then
        LAST_LINE=$(tail -1 "$TEMP_RUN_CSV")
        echo "$LAST_LINE" >> "$CSV_FILE"
        rm -f "$TEMP_RUN_CSV"
    fi
done

# Clean up any remaining temp files
rm -f /tmp/runlist_* /tmp/nps_output_* /tmp/csv_run_* /tmp/csv_backup_* 2>/dev/null

echo ""
echo "============================================================================"
echo "Analysis Complete!"
echo "End time: $(date)"
echo "============================================================================"
echo ""
echo "Results Summary:"
echo "  Total runs: $TOTAL_RUNS"
echo "  Successful: $N_SUCCESS"
echo "  Failed: $N_FAILED"
echo "  CSV file: $CSV_FILE"
echo ""
if [[ $N_SUCCESS -gt 0 ]]; then
    echo "✓ Successful runs: ${SUCCESSFUL_RUNS[@]}"
fi
if [[ $N_FAILED -gt 0 ]]; then
    echo "✗ Failed runs: ${FAILED[@]}"
fi
echo "============================================================================"
echo ""
echo "CSV Summary:"
tail -$((N_SUCCESS+1)) "$CSV_FILE"
echo "============================================================================"

# Clean up any remaining temp files and corrupted files
rm -f /tmp/runlist_* /tmp/nps_output_* /tmp/csv_run_* /tmp/csv_backup_* 2>/dev/null

# Remove gibberish files (non-standard filenames) from current directory
# These are corrupted temp files created by ROOT during crashes
# They have non-printable characters in filenames and contain CSV data
find . -maxdepth 1 -type f -size 100c,200c ! -name "*.sh" ! -name ".*" 2>/dev/null | while read f; do
    # Check if file looks like CSV data (numeric start)
    if head -c 20 "$f" 2>/dev/null | grep -q '^[0-9]'; then
        rm -f "$f" 2>/dev/null
    fi
done

# Additional cleanup: remove files with non-ASCII in filename (pattern-based)
# The corrupted files have octal escape sequences in names - target those specifically
for f in *; do
    # Skip normal files
    [[ "$f" == "run_nps_analysis.sh" ]] && continue
    [[ "$f" == "src" ]] && continue
    [[ "$f" == "config" ]] && continue
    [[ "$f" == "output" ]] && continue
    [[ "$f" == "." ]] && continue
    [[ "$f" == ".." ]] && continue
    
    # If filename matches pattern [a-zA-Z0-9._-/], keep it
    # Otherwise it's likely corrupted, remove it
    if ! echo "$f" | grep -q '^[a-zA-Z0-9._-/]*$' 2>/dev/null; then
        if [[ -f "$f" ]] && [[ ! -h "$f" ]]; then
            rm -f "$f" 2>/dev/null
        fi
    fi
done
