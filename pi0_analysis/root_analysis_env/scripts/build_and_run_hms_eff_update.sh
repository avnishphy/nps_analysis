#!/bin/bash
# build_and_run_hms_eff_update.sh
# Compile and run the HMS tracking efficiency updater

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROG_NAME="update_hms_track_eff"
SOURCE_FILE="${SCRIPT_DIR}/${PROG_NAME}.cpp"
OUTPUT_BINARY="${SCRIPT_DIR}/${PROG_NAME}"

# ROOT directories (adjust if needed)
ROOT_DIR="/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b"
SUMMARY_CSV="${SCRIPT_DIR}/../output/plots/x60_4b/summary_all_runs.csv"

echo "============================================"
echo "HMS Tracking Efficiency Updater"
echo "============================================"
echo ""

# Check if ROOT is available
if ! command -v root-config &> /dev/null; then
    echo "[ERROR] ROOT not found. Please set up ROOT first:"
    echo "  source /path/to/root/bin/thisroot.sh"
    exit 1
fi

echo "[INFO] ROOT version: $(root-config --version)"
echo ""

# Compile
echo "[INFO] Compiling ${PROG_NAME}..."
echo "Command: g++ -O3 ${SOURCE_FILE} \$(root-config --cflags --libs) -o ${OUTPUT_BINARY}"
echo ""

g++ -O3 "${SOURCE_FILE}" $(root-config --cflags --libs) -o "${OUTPUT_BINARY}"

if [ -f "${OUTPUT_BINARY}" ]; then
    echo "[SUCCESS] Compilation complete: ${OUTPUT_BINARY}"
    echo ""
else
    echo "[ERROR] Compilation failed"
    exit 1
fi

# Check if CSV exists
if [ ! -f "${SUMMARY_CSV}" ]; then
    echo "[ERROR] Summary CSV not found: ${SUMMARY_CSV}"
    exit 1
fi

if [ ! -d "${ROOT_DIR}" ]; then
    echo "[ERROR] ROOT directory not found: ${ROOT_DIR}"
    exit 1
fi

# Run
echo "============================================"
echo "Running HMS Efficiency Updater"
echo "============================================"
echo ""
echo "[INFO] ROOT files directory: ${ROOT_DIR}"
echo "[INFO] Summary CSV: ${SUMMARY_CSV}"
echo ""

"${OUTPUT_BINARY}" --root-dir "${ROOT_DIR}" --summary-csv "${SUMMARY_CSV}"

echo ""
echo "============================================"
echo "[SUCCESS] Done!"
echo "============================================"
