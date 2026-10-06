#!/usr/bin/env bash
# KinC_x36_4 validated subtraction/bootstrap workflow; load Hall C ROOT first.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
OUT="${NPS_BOOTSTRAP_OUTPUT:-${REPO}/output/pi0_bootstrap_reproduction}"
MANIFEST="${NPS_BOOTSTRAP_MANIFEST:-${REPO}/validation/pi0_bootstrap_20261004/run_quality_manifest.csv}"
KIN=KinC_x36_4
CAMPAIGN="${OUT}/nominal/${KIN}"
ROOT_DIR="${CAMPAIGN}/root"
LIB="${OUT}/build/libnps_stat.so"
SIM=/lustre24/expphy/volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root
VERTEX="${REPO}/output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x36_4.root"
B="${NPS_BOOTSTRAP_REPLICAS:-500}"
build() {
  command -v root-config >/dev/null || { echo 'Load /group/nps/singhav/setup.csh through csh first.' >&2; exit 1; }
  mkdir -p "${OUT}/build"
  if [[ ! -f "$LIB" || src/analysis/nps_stat_bridge.cpp -nt "$LIB" || src/analysis/nps_comb_bg_pepsi.h -nt "$LIB" || src/xsec_extract/xsec_linear_solver.h -nt "$LIB" ]]; then
    g++ -O2 -std=c++17 -shared -fPIC src/analysis/nps_stat_bridge.cpp $(root-config --cflags --libs) -lMinuit2 -o "${LIB}.tmp"
    mv "${LIB}.tmp" "$LIB"
  fi
}
xsec() {
  bash src/xsec_extract/run_xsec_pipeline.sh --kin "$KIN" --target LH2 --data-file "$1" \
    --sim-file "$SIM" --vertex_simc_file "$VERTEX" --xsec_config xsec_config_x36_4.json \
    --fit-objective gaussian --fit-variance data --out-dir "$2" --no-pdf --no-png --no-diagnostics
}
response() {
  if [[ ! -f "${OUT}/response_reference/migration_design.csv" ]]; then
    xsec "${ROOT_DIR}/combined_branches_LH2.root" "${OUT}/response_reference"
  fi
}
bootstrap() {
  build; response
  python3 scripts/bootstrap_pi0_data.py --campaign "$CAMPAIGN" --output "${OUT}/bootstrap_$1" \
    --library "$LIB" --manifest "$MANIFEST" --response-dir "${OUT}/response_reference" \
    --replicas "$1" --seed 20261004 --estimator selected
}
case "${1:-help}" in
  nominal)
    bash src/analysis/run_parallel_nps_analysis_main.sh --kin "$KIN" --source waveform \
      --gevnum-cut yes --target LH2 --jobs 3 --types production,Production --no-combine \
      --run-quality-manifest "$MANIFEST" --output-base "${OUT}/nominal" ;;
  combine)
    python3 src/analysis/combine_analysis_branches.py --kin "$KIN" --output-base "${OUT}/nominal" \
      --efficiency-csv output/efficiency_stuff/efficiency_KinC_x36_4.csv \
      --run-quality-manifest "$MANIFEST" --no-analysis-plots ;;
  quick) bootstrap 20 ;;
  validate100) bootstrap 100 ;;
  production) bootstrap "$B" ;;
  extract)
    build; response
    python3 scripts/bootstrap_pi0_data.py --campaign "$CAMPAIGN" --output "${OUT}/nominal_export" \
      --library "$LIB" --manifest "$MANIFEST" --response-dir "${OUT}/response_reference" \
      --replicas 2 --estimator selected --combined-data "${ROOT_DIR}/combined_branches_LH2.root" \
      --write-nominal-input "${ROOT_DIR}/combined_branches_LH2_selected.root"
    xsec "${ROOT_DIR}/combined_branches_LH2_selected.root" "${OUT}/selected_nominal"
    python3 scripts/report_pi0_bootstrap.py --bootstrap "${OUT}/bootstrap_${B}" \
      --extraction "${OUT}/selected_nominal" --campaign "$CAMPAIGN" --manifest "$MANIFEST" --output "${OUT}/statistical_scales"
    python3 scripts/bootstrap_pi0_data.py --campaign "$CAMPAIGN" --output "${OUT}/covariance_export" \
      --library "$LIB" --manifest "$MANIFEST" --response-dir "${OUT}/response_reference" \
      --replicas 2 --estimator selected --combined-data "${ROOT_DIR}/combined_branches_LH2.root" \
      --write-nominal-input "${ROOT_DIR}/combined_branches_LH2_bootstrap.root" \
      --covariance-from "${OUT}/bootstrap_${B}" \
      --mc-covariance "${OUT}/statistical_scales/finite_mc_delta_sigma_covariance.csv"
    xsec "${ROOT_DIR}/combined_branches_LH2_bootstrap.root" "${OUT}/preliminary_xsec" ;;
  *) echo 'Usage: bash scripts/run_pi0_bootstrap.sh {nominal|combine|quick|validate100|production|extract}' ;;
esac
