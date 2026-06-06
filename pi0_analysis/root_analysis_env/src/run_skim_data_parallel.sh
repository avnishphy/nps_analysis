#!/bin/bash

# Number of parallel jobs (default: number of CPU cores)
JOBS=$(nproc)

# Path to the main runlist
RUNLIST="../config/runlist_x60_4b.txt"

# Output directory for temporary runlists
TMPDIR="tmp_runlists"
mkdir -p "$TMPDIR"

# Remove any old chunk files
rm -f "$TMPDIR"/runlist_chunk_*

# Clean runlist to remove empty lines
grep -ve '^\s*$' "$RUNLIST" > "$TMPDIR/runlist_cleaned.txt"
# Split cleaned runlist into $JOBS chunks
split -d -n l/$JOBS "$TMPDIR/runlist_cleaned.txt" "$TMPDIR/runlist_chunk_"

# Launch one job per chunk with progress info
PIDS=()
CHUNK_IDX=0
for chunk in "$TMPDIR"/runlist_chunk_*; do
    CHUNK_IDX=$((CHUNK_IDX+1))
    echo "[INFO] Starting job $CHUNK_IDX for $chunk..."
    (
        root -b -q "skim_data.C(\"/cache/hallc/c-nps/analysis/pass2/replays/updated/\", \"$chunk\", \"/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b/\")"
        echo "[INFO] Finished job $CHUNK_IDX for $chunk."
    ) &
    PIDS+=("$!")
done

# Wait for all jobs to finish, showing progress
for idx in "${!PIDS[@]}"; do
    pid=${PIDS[$idx]}
    wait $pid
    echo "[INFO] Job $((idx+1)) (PID $pid) completed."
done

echo "All parallel jobs completed."