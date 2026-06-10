#!/bin/bash
# ============================================================================
# End-to-end NPS smearing pipeline
# 1) Compile nps_sim_smearing_new_try.C
# 2) Run section-wise smearing fit (CSV + interpolated maps)
# 3) Regenerate smeared SIMC ROOT output via simc_pi0_analysis.C
# ============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

DATA_FILE="${DATA_FILE:-/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root}"
SIM_FILE="${SIM_FILE:-/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output.root}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/output/plots/x60_4b/production_wfpi0}"
OUT_FILE="${OUT_FILE:-$OUT_DIR/out_smear.root}"
FITTER_SRC="${FITTER_SRC:-scripts/nps_sim_smearing_new_try.C}"
FITTER_BIN="${FITTER_BIN:-scripts/nps_sim_smearing_new_try}"

NX="${NX:-11}"
NY="${NY:-16}"
X_MIN="${X_MIN:--24}"
X_MAX="${X_MAX:-28}"
Y_MIN="${Y_MIN:--34}"
Y_MAX="${Y_MAX:-34}"
OVERLAP_FRAC="${OVERLAP_FRAC:-0.0}"
NSMEAR="${NSMEAR:-80}"

SECTION_MAP="$OUT_DIR/section_map.csv"
CHI2_PDF="$OUT_DIR/chi2_scans.pdf"
OPTIMIZER_SUMMARY="$OUT_DIR/smearing_optimizer_summary.csv"
OPTIMIZER_SEEDS="$OUT_DIR/smearing_optimizer_seeds.csv"
OPTIMIZER_PROFILES="$OUT_DIR/smearing_optimizer_profiles.csv"
CLOSURE_SUMMARY="$OUT_DIR/smearing_closure_summary.csv"
INTERP_FILE="${OUT_FILE%.root}_interpolated.root"
SMEARED_SIMC_OUT="${SMEARED_SIMC_OUT:-$OUT_DIR/simc_pi0_analysis_output_smeared.root}"
UNSMEARED_SIMC_OUT="${UNSMEARED_SIMC_OUT:-$OUT_DIR/simc_pi0_analysis_output.root}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
BACKUP_DIR="$OUT_DIR/smearing_pipeline_backup/$RUN_TAG"

log_step() {
  echo ""
  echo "============================================================================"
  echo "$1"
  echo "============================================================================"
}

require_file() {
  local path="$1"
  local label="$2"
  if [[ ! -f "$path" ]]; then
    echo "[ERROR] Missing $label: $path"
    exit 1
  fi
}

require_command() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[ERROR] $cmd not found in PATH"
    exit 1
  fi
}

archive_if_exists() {
  local path="$1"
  if [[ -e "$path" ]]; then
    mkdir -p "$BACKUP_DIR"
    mv "$path" "$BACKUP_DIR/"
    echo "[INFO] Archived previous output: $path -> $BACKUP_DIR/"
  fi
}

verify_output() {
  local path="$1"
  local label="$2"
  require_file "$path" "$label"
  if [[ ! -s "$path" ]]; then
    echo "[ERROR] Empty $label: $path"
    exit 1
  fi
}

require_command root
require_command root-config
require_file "$FITTER_SRC" "fitter source"
require_file "scripts/simc_pi0_analysis.C" "downstream macro"
require_file "$DATA_FILE" "data file"
require_file "$SIM_FILE" "simulation file"

mkdir -p "$OUT_DIR"
if [[ ! -w "$OUT_DIR" ]]; then
  echo "[ERROR] Output directory is not writable: $OUT_DIR"
  exit 1
fi

log_step "Smearing Pipeline Configuration"
echo "PROJECT_ROOT=$PROJECT_ROOT"
echo "DATA_FILE=$DATA_FILE"
echo "SIM_FILE=$SIM_FILE"
echo "OUT_DIR=$OUT_DIR"
echo "OUT_FILE=$OUT_FILE"
echo "INTERP_FILE=$INTERP_FILE"
echo "SECTION_MAP=$SECTION_MAP"
echo "CHI2_PDF=$CHI2_PDF"
echo "OPTIMIZER_SUMMARY=$OPTIMIZER_SUMMARY"
echo "OPTIMIZER_SEEDS=$OPTIMIZER_SEEDS"
echo "OPTIMIZER_PROFILES=$OPTIMIZER_PROFILES"
echo "CLOSURE_SUMMARY=$CLOSURE_SUMMARY"
echo "SMEARED_SIMC_OUT=$SMEARED_SIMC_OUT"
echo "GRID nx=$NX ny=$NY x=[$X_MIN,$X_MAX] y=[$Y_MIN,$Y_MAX] overlap=$OVERLAP_FRAC Nsmear=$NSMEAR"
echo "Energy mean model: a_plus_bE_plus_clnE"
echo "Exclusive selector branch: is_exclusive_ellipse_combined (data) and is_exclusive_ellipse (sim)"

log_step "Archiving previous canonical outputs"
archive_if_exists "$OUT_FILE"
archive_if_exists "$INTERP_FILE"
archive_if_exists "$SECTION_MAP"
archive_if_exists "$CHI2_PDF"
archive_if_exists "$OPTIMIZER_SUMMARY"
archive_if_exists "$OPTIMIZER_SEEDS"
archive_if_exists "$OPTIMIZER_PROFILES"
archive_if_exists "$CLOSURE_SUMMARY"
archive_if_exists "$SMEARED_SIMC_OUT"

log_step "Step 1/3: Compiling $FITTER_SRC"
g++ "$FITTER_SRC" $(root-config --cflags --libs) -lMathMore -O2 -std=c++17 -fopenmp -I../src -o "$FITTER_BIN"
verify_output "$FITTER_BIN" "compiled fitter binary"

log_step "Step 2/3: Running section smearing fit"
"$FITTER_BIN" \
  "$DATA_FILE" physics \
  "$SIM_FILE" simulation \
  "$OUT_FILE" "$NX" "$NY" "$X_MIN" "$X_MAX" "$Y_MIN" "$Y_MAX" "$OVERLAP_FRAC" "$NSMEAR"

verify_output "$OUT_FILE" "fitter ROOT output"
verify_output "$INTERP_FILE" "interpolated smearing ROOT output"
verify_output "$SECTION_MAP" "section map CSV"
verify_output "$CHI2_PDF" "chi2 PDF"
verify_output "$OPTIMIZER_SUMMARY" "optimizer summary CSV"
verify_output "$OPTIMIZER_SEEDS" "optimizer Sobol seed CSV"
verify_output "$OPTIMIZER_PROFILES" "optimizer profile CSV"
verify_output "$CLOSURE_SUMMARY" "closure consistency CSV"

log_step "Step 3/3: Regenerating smeared simulation ROOT output"
NPS_SMEARING_MODE=section \
NPS_SMEAR_FILE="$OUT_FILE" \
NPS_SMEAR_INTERP_FILE="$INTERP_FILE" \
NPS_SECTION_MAP_FILE="$SECTION_MAP" \
NPS_SIMC_OUTPUT_FILE="$UNSMEARED_SIMC_OUT" \
NPS_SIMC_SMEARED_OUTPUT_FILE="$SMEARED_SIMC_OUT" \
root -l -b -q scripts/simc_pi0_analysis.C

verify_output "$SMEARED_SIMC_OUT" "smeared SIMC ROOT output"

log_step "Pipeline complete"
echo "Smearing fit: $OUT_FILE"
echo "Interpolated maps: $INTERP_FILE"
echo "Section map: $SECTION_MAP"
echo "Diagnostics: $CHI2_PDF"
echo "Optimizer summary: $OPTIMIZER_SUMMARY"
echo "Sobol seeds: $OPTIMIZER_SEEDS"
echo "Profiles: $OPTIMIZER_PROFILES"
echo "Closure consistency: $CLOSURE_SUMMARY"
echo "Smeared SIMC: $SMEARED_SIMC_OUT"
