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
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
CANONICAL_OUT_FILE="${OUT_FILE:-$OUT_DIR/out_smear.root}"
if [[ "$CANONICAL_OUT_FILE" == *.root ]]; then
  OUT_FILE="${CANONICAL_OUT_FILE%.root}_${RUN_TAG}.root"
else
  OUT_FILE="${CANONICAL_OUT_FILE}_${RUN_TAG}.root"
fi
FITTER_SRC="${FITTER_SRC:-scripts/nps_sim_smearing_new_try.C}"
FITTER_BIN="${FITTER_BIN:-scripts/nps_sim_smearing_new_try}"

NX="${NX:-7}"
NY="${NY:-7}"
X_MIN="${X_MIN:--24}"
X_MAX="${X_MAX:-28}"
Y_MIN="${Y_MIN:--34}"
Y_MAX="${Y_MAX:-34}"
OVERLAP_FRAC="${OVERLAP_FRAC:-0.1}"
NSMEAR="${NSMEAR:-80}"

SECTION_MAP="$OUT_DIR/section_map.csv"
CANONICAL_CHI2_PDF="$OUT_DIR/chi2_scans.pdf"
CHI2_PDF="$OUT_DIR/chi2_scans_${RUN_TAG}.pdf"
CHI2_PROGRESS_DIR="$OUT_DIR/chi2_scans_progress_${RUN_TAG}"
ARTIFACT_MANIFEST="$OUT_DIR/smearing_artifacts_${RUN_TAG}.json"
OPTIMIZER_SUMMARY="$OUT_DIR/smearing_optimizer_summary.csv"
OPTIMIZER_SEEDS="$OUT_DIR/smearing_optimizer_seeds.csv"
OPTIMIZER_PROFILES="$OUT_DIR/smearing_optimizer_profiles.csv"
CLOSURE_SUMMARY="$OUT_DIR/smearing_closure_summary.csv"
SWEEP_HISTORY="$OUT_DIR/smearing_sweep_history.csv"
OBJECTIVE_BREAKDOWN="$OUT_DIR/smearing_objective_breakdown.csv"
INTERP_FILE="${OUT_FILE%.root}_interpolated.root"
CANONICAL_INTERP_FILE="${CANONICAL_OUT_FILE%.root}_interpolated.root"
SMEARED_SIMC_OUT="${SMEARED_SIMC_OUT:-$OUT_DIR/simc_pi0_analysis_output_smeared.root}"
UNSMEARED_SIMC_OUT="${UNSMEARED_SIMC_OUT:-$OUT_DIR/simc_pi0_analysis_output.root}"
RUN_COMPARISON="${RUN_COMPARISON:-1}"
COMPARISON_OUT_DIR="${COMPARISON_OUT_DIR:-$OUT_DIR/smearing_comparison}"
COMPARISON_MAX_ENTRIES="${COMPARISON_MAX_ENTRIES:-}"
NPS_SMEARING_SWEEP_ACCEPTANCE="${NPS_SMEARING_SWEEP_ACCEPTANCE:-jacobi_global_accept_rollback}"
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
echo "CANONICAL_OUT_FILE=$CANONICAL_OUT_FILE"
echo "INTERP_FILE=$INTERP_FILE"
echo "CANONICAL_INTERP_FILE=$CANONICAL_INTERP_FILE"
echo "SECTION_MAP=$SECTION_MAP"
echo "CHI2_PDF=$CHI2_PDF"
echo "CANONICAL_CHI2_PDF=$CANONICAL_CHI2_PDF"
echo "CHI2_PROGRESS_DIR=$CHI2_PROGRESS_DIR"
echo "ARTIFACT_MANIFEST=$ARTIFACT_MANIFEST"
echo "OPTIMIZER_SUMMARY=$OPTIMIZER_SUMMARY"
echo "OPTIMIZER_SEEDS=$OPTIMIZER_SEEDS"
echo "OPTIMIZER_PROFILES=$OPTIMIZER_PROFILES"
echo "CLOSURE_SUMMARY=$CLOSURE_SUMMARY"
echo "SWEEP_HISTORY=$SWEEP_HISTORY"
echo "OBJECTIVE_BREAKDOWN=$OBJECTIVE_BREAKDOWN"
echo "SMEARED_SIMC_OUT=$SMEARED_SIMC_OUT"
echo "RUN_COMPARISON=$RUN_COMPARISON"
echo "COMPARISON_OUT_DIR=$COMPARISON_OUT_DIR"
echo "COMPARISON_MAX_ENTRIES=${COMPARISON_MAX_ENTRIES:-full}"
echo "NPS_SMEARING_SWEEP_ACCEPTANCE=$NPS_SMEARING_SWEEP_ACCEPTANCE"
echo "GRID nx=$NX ny=$NY x=[$X_MIN,$X_MAX] y=[$Y_MIN,$Y_MAX] overlap=$OVERLAP_FRAC Nsmear=$NSMEAR"
echo "Energy mean model: a_plus_bE_plus_clnE"
echo "Exclusive selector branch: is_exclusive_ellipse_combined (data) and is_exclusive_ellipse (sim)"

log_step "Archiving previous canonical outputs"
archive_if_exists "$CANONICAL_OUT_FILE"
archive_if_exists "$CANONICAL_INTERP_FILE"
archive_if_exists "$SECTION_MAP"
archive_if_exists "$CANONICAL_CHI2_PDF"
archive_if_exists "$OPTIMIZER_SUMMARY"
archive_if_exists "$OPTIMIZER_SEEDS"
archive_if_exists "$OPTIMIZER_PROFILES"
archive_if_exists "$CLOSURE_SUMMARY"
archive_if_exists "$SWEEP_HISTORY"
archive_if_exists "$OBJECTIVE_BREAKDOWN"
archive_if_exists "$SMEARED_SIMC_OUT"

log_step "Step 1/3: Compiling $FITTER_SRC"
g++ "$FITTER_SRC" $(root-config --cflags --libs) -lMathMore -O2 -std=c++17 -fopenmp -I../src -o "$FITTER_BIN"
verify_output "$FITTER_BIN" "compiled fitter binary"

log_step "Step 2/3: Running section smearing fit"
NPS_SMEARING_SWEEP_ACCEPTANCE="$NPS_SMEARING_SWEEP_ACCEPTANCE" \
RUN_TAG="$RUN_TAG" \
"$FITTER_BIN" \
  "$DATA_FILE" physics \
  "$SIM_FILE" simulation \
  "$OUT_FILE" "$NX" "$NY" "$X_MIN" "$X_MAX" "$Y_MIN" "$Y_MAX" "$OVERLAP_FRAC" "$NSMEAR"

verify_output "$OUT_FILE" "fitter ROOT output"
verify_output "$INTERP_FILE" "interpolated smearing ROOT output"
cp -p "$OUT_FILE" "$CANONICAL_OUT_FILE"
cp -p "$INTERP_FILE" "$CANONICAL_INTERP_FILE"
verify_output "$CANONICAL_OUT_FILE" "latest fitter ROOT output"
verify_output "$CANONICAL_INTERP_FILE" "latest interpolated smearing ROOT output"
verify_output "$SECTION_MAP" "section map CSV"
verify_output "$CHI2_PDF" "chi2 PDF"
verify_output "$CANONICAL_CHI2_PDF" "latest chi2 PDF"
verify_output "$ARTIFACT_MANIFEST" "artifact manifest"
verify_output "$OPTIMIZER_SUMMARY" "optimizer summary CSV"
verify_output "$OPTIMIZER_SEEDS" "optimizer Sobol seed CSV"
verify_output "$OPTIMIZER_PROFILES" "optimizer profile CSV"
verify_output "$CLOSURE_SUMMARY" "closure consistency CSV"
verify_output "$SWEEP_HISTORY" "coupled sweep history CSV"
verify_output "$OBJECTIVE_BREAKDOWN" "objective breakdown CSV"

log_step "Step 3/3: Regenerating smeared simulation ROOT output"
NPS_SMEARING_MODE=section \
NPS_SMEAR_FILE="$OUT_FILE" \
NPS_SMEAR_INTERP_FILE="$INTERP_FILE" \
NPS_SECTION_MAP_FILE="$SECTION_MAP" \
NPS_SIMC_OUTPUT_FILE="$UNSMEARED_SIMC_OUT" \
NPS_SIMC_SMEARED_OUTPUT_FILE="$SMEARED_SIMC_OUT" \
root -l -b -q scripts/simc_pi0_analysis.C

verify_output "$SMEARED_SIMC_OUT" "smeared SIMC ROOT output"

if [[ "$RUN_COMPARISON" == "1" ]]; then
  log_step "Step 4/4: Writing reproducible comparison metrics"
  if command -v python3 >/dev/null 2>&1 && \
     (python3 -c 'import ROOT, numpy, matplotlib' >/dev/null 2>&1 || \
      python3 -c 'import uproot, numpy, matplotlib' >/dev/null 2>&1); then
    comparison_args=()
    if [[ -n "$COMPARISON_MAX_ENTRIES" ]]; then
      comparison_args+=(--max-entries "$COMPARISON_MAX_ENTRIES")
    fi
    python3 scripts/compare_smearing_outputs.py \
      --data "$DATA_FILE" \
      --sim "$UNSMEARED_SIMC_OUT" \
      --smeared "$SMEARED_SIMC_OUT" \
      --out-dir "$COMPARISON_OUT_DIR" \
      "${comparison_args[@]}"
    verify_output "$COMPARISON_OUT_DIR/smearing_comparison_metrics.csv" "comparison metrics CSV"
    verify_output "$COMPARISON_OUT_DIR/smearing_comparison.png" "comparison plot"
  else
    echo "[WARN] Python comparison skipped: python3 with ROOT or uproot, plus numpy and matplotlib, is required."
  fi
fi

log_step "Writing compact run summaries"
python3 scripts/summarize_smearing_run.py \
  --out-dir "$OUT_DIR" \
  --strategy "$NPS_SMEARING_SWEEP_ACCEPTANCE"
verify_output "$OUT_DIR/smearing_runtime_summary.csv" "runtime summary CSV"
verify_output "$OUT_DIR/best_parameters.json" "best-parameters JSON"
verify_output "$OUT_DIR/optimization_config.json" "optimization config JSON"

log_step "Pipeline complete"
echo "Smearing fit: $OUT_FILE"
echo "Latest smearing fit copy: $CANONICAL_OUT_FILE"
echo "Interpolated maps: $INTERP_FILE"
echo "Latest interpolated maps copy: $CANONICAL_INTERP_FILE"
echo "Section map: $SECTION_MAP"
echo "Diagnostics: $CHI2_PDF"
echo "Latest diagnostics copy: $CANONICAL_CHI2_PDF"
echo "Progress diagnostics: $CHI2_PROGRESS_DIR"
echo "Artifact manifest: $ARTIFACT_MANIFEST"
echo "Optimizer summary: $OPTIMIZER_SUMMARY"
echo "Sobol seeds: $OPTIMIZER_SEEDS"
echo "Profiles: $OPTIMIZER_PROFILES"
echo "Closure consistency: $CLOSURE_SUMMARY"
echo "Coupled sweep history: $SWEEP_HISTORY"
echo "Objective breakdown: $OBJECTIVE_BREAKDOWN"
echo "Smeared SIMC: $SMEARED_SIMC_OUT"
echo "Runtime summary: $OUT_DIR/smearing_runtime_summary.csv"
echo "Best parameters JSON: $OUT_DIR/best_parameters.json"
echo "Optimization config JSON: $OUT_DIR/optimization_config.json"
if [[ "$RUN_COMPARISON" == "1" ]]; then
  echo "Comparison output: $COMPARISON_OUT_DIR"
fi
