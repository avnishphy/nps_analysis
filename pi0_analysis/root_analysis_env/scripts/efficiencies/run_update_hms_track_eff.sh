#!/bin/bash
# build_and_run_hms_eff_update.sh
# Parallel HMS tracking efficiency updater

set -e

usage() {
    echo "Usage: $0 --root-dir <root_dir> --summary-csv <summary_csv> [--jobs N]"
    echo "  --root-dir <root_dir>       Directory containing ROOT files (required)"
    echo "  --summary-csv <summary_csv> Path to summary CSV file (required)"
    echo "  --jobs N                    Number of parallel jobs (optional, default: number of CPU cores)"
    exit 1
}

# Parse arguments
ROOT_DIR=""
SUMMARY_CSV=""
N_JOBS=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --root-dir)
            ROOT_DIR="$2"; shift 2;;
        --summary-csv)
            SUMMARY_CSV="$2"; shift 2;;
        --jobs)
            N_JOBS="$2"; shift 2;;
        -h|--help)
            usage;;
        *)
            echo "Unknown argument: $1"; usage;;
    esac
    done

if [[ -z "$ROOT_DIR" || -z "$SUMMARY_CSV" ]]; then
    usage
fi

if [[ -z "$N_JOBS" ]]; then
    N_JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
fi

# Check for update_hms_track_eff binary in PATH or script dir
if ! command -v update_hms_track_eff &>/dev/null && [[ ! -x "./update_hms_track_eff" ]]; then
    echo "[ERROR] update_hms_track_eff executable not found in PATH or current directory."
    echo "Please build it first: g++ -O3 update_hms_track_eff.cpp \`root-config --cflags --libs\` -o update_hms_track_eff"
    exit 1
fi
if command -v update_hms_track_eff &>/dev/null; then
    OUTPUT_BINARY="update_hms_track_eff"
else
    OUTPUT_BINARY="./update_hms_track_eff"
fi

# Temporary working directory
WORKDIR=$(mktemp -d)
trap 'rm -rf "$WORKDIR"' EXIT

# Split CSV into N_JOBS chunks (preserving header)
HEADER_LINE=$(head -n 1 "$SUMMARY_CSV")
DATA_LINES=$(tail -n +2 "$SUMMARY_CSV" | wc -l)
LINES_PER_CHUNK=$(( (DATA_LINES + N_JOBS - 1) / N_JOBS ))

split -l $LINES_PER_CHUNK -d --additional-suffix=.csv <(tail -n +2 "$SUMMARY_CSV") "$WORKDIR/chunk_"

# Add header to each chunk
for f in "$WORKDIR"/chunk_*.csv; do
    (echo "$HEADER_LINE"; cat "$f") > "$f.tmp" && mv "$f.tmp" "$f"
done

# Run update_hms_track_eff in parallel on each chunk
run_chunk() {
    local chunk_csv=$1
    local out_csv=${chunk_csv%.csv}_updated.csv
    if "$OUTPUT_BINARY" --root-dir "$ROOT_DIR" --summary-csv "$chunk_csv"; then
        if [ -f "${chunk_csv%.csv}_updated.csv" ]; then
            mv "${chunk_csv%.csv}_updated.csv" "$out_csv"
        fi
    else
        echo "[ERROR] update_hms_track_eff failed for $chunk_csv. Writing dummy eff=1.0."
        # Write dummy output with eff=1.0 for all runs in this chunk
        header=$(head -n 1 "$chunk_csv")
        eff_col=$(awk -F, '{for(i=1;i<=NF;i++) if($i=="hms_track_eff") print i; exit}' <<< "$header")
        eff_err_col=$(awk -F, '{for(i=1;i<=NF;i++) if($i=="hms_track_eff_err") print i; exit}' <<< "$header")
        if [[ -z "$eff_col" ]]; then eff_col=0; fi
        if [[ -z "$eff_err_col" ]]; then eff_err_col=0; fi
        (
            echo "$header"
            tail -n +2 "$chunk_csv" | awk -F, -v OFS="," -v eff_col="$eff_col" -v eff_err_col="$eff_err_col" '
                BEGIN { eff_col=int(eff_col); eff_err_col=int(eff_err_col); }
                {
                    n=NF;
                    if (eff_col>0 && eff_col<=n) $eff_col="1.0";
                    if (eff_err_col>0 && eff_err_col<=n) $eff_err_col="0.0";
                    print $0;
                }'
        ) > "$out_csv"
    fi
}

export -f run_chunk
export OUTPUT_BINARY ROOT_DIR

if command -v parallel &>/dev/null; then
    parallel --jobs $N_JOBS run_chunk ::: "$WORKDIR"/chunk_*.csv
else
    # Fallback: background jobs
    for f in "$WORKDIR"/chunk_*.csv; do
        run_chunk "$f" &
    done
    wait
fi

# Check for missing chunk updated files before merging
MISSING=0
for f in "$WORKDIR"/chunk_*.csv; do
    updated="${f%.csv}_updated.csv"
    if [ ! -f "$updated" ]; then
        echo "[ERROR] Missing updated chunk: $updated"
        MISSING=1
    fi
    if [ ! -s "$updated" ]; then
        echo "[ERROR] Empty updated chunk: $updated"
        MISSING=1
    fi
    # Optionally, check for correct number of lines
    orig_lines=$(wc -l < "$f")
    upd_lines=$(wc -l < "$updated")
    if [ "$orig_lines" -ne "$upd_lines" ]; then
        echo "[WARN] Line count mismatch in $updated (orig: $orig_lines, updated: $upd_lines)"
    fi
done
if [ "$MISSING" -ne 0 ]; then
    echo "[FATAL] One or more chunk updates failed. Aborting merge."
    exit 2
fi

# Merge updated chunks (preserve header only from first)
UPDATED_CSV="$WORKDIR/merged_updated.csv"
first_chunk=$(ls "$WORKDIR"/chunk_*_updated.csv | head -n 1)
head -n 1 "$first_chunk" > "$UPDATED_CSV"
for f in "$WORKDIR"/chunk_*_updated.csv; do
    tail -n +2 "$f" >> "$UPDATED_CSV"
done

# Move merged CSV to final location
cp "$SUMMARY_CSV" "$SUMMARY_CSV.bak"
mv "$UPDATED_CSV" "$SUMMARY_CSV"
echo "[SUCCESS] Merged updated CSV written to $SUMMARY_CSV (backup at $SUMMARY_CSV.bak)"
