#!/bin/bash
# ============================================================================
# End-to-end NPS smearing pipeline
# 1) Compile nps_sim_smearing_new.C
# 2) Run section-wise smearing fit (CSV + interpolated maps)
# 3) Regenerate smeared SIMC ROOT output via simc_pi0_analysis.C
# ============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

DATA_FILE="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root"
SIM_FILE="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output.root"
OUT_DIR="$PROJECT_ROOT/output/plots/x60_4b/production_wfpi0"
OUT_FILE="$OUT_DIR/out_smear.root"

if ! command -v root >/dev/null 2>&1; then
  echo "[ERROR] ROOT executable not found in PATH"
  exit 1
fi

if ! command -v root-config >/dev/null 2>&1; then
  echo "[ERROR] root-config not found in PATH"
  exit 1
fi

if [[ ! -f "scripts/nps_sim_smearing_new.C" ]]; then
  echo "[ERROR] Missing source: scripts/nps_sim_smearing_new.C"
  exit 1
fi

if [[ ! -f "scripts/simc_pi0_analysis.C" ]]; then
  echo "[ERROR] Missing macro: scripts/simc_pi0_analysis.C"
  exit 1
fi

if [[ ! -f "$DATA_FILE" ]]; then
  echo "[ERROR] Missing data file: $DATA_FILE"
  exit 1
fi

if [[ ! -f "$SIM_FILE" ]]; then
  echo "[ERROR] Missing simulation file: $SIM_FILE"
  exit 1
fi

mkdir -p "$OUT_DIR"
if [[ ! -w "$OUT_DIR" ]]; then
  echo "[ERROR] Output directory is not writable: $OUT_DIR"
  exit 1
fi

echo "============================================================================"
echo "Step 1/3: Compiling scripts/nps_sim_smearing_new.C"
echo "============================================================================"
g++ scripts/nps_sim_smearing_new.C $(root-config --cflags --libs) -O2 -std=c++17 -fopenmp -I../src -o scripts/nps_sim_smearing_new

echo ""
echo "============================================================================"
echo "Step 2/3: Running section smearing fit"
echo "============================================================================"
./scripts/nps_sim_smearing_new \
  "$DATA_FILE" physics \
  "$SIM_FILE" simulation \
  "$OUT_FILE" 11 16 -24 28 -34 34 0.0 80

echo ""
echo "============================================================================"
echo "Step 3/3: Regenerating smeared simulation ROOT output"
echo "============================================================================"
NPS_SMEARING_MODE=section root -l -b -q scripts/simc_pi0_analysis.C

echo ""
echo "============================================================================"
echo "Pipeline complete"
echo "============================================================================"
