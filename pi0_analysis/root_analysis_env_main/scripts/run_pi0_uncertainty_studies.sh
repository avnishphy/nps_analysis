#!/usr/bin/env bash
# Follow-up evidence for the validated KinC_x36_4 baseline; load Hall C ROOT.
# Other settings are deliberately refused until coverage/configuration validation.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
BASE="${NPS_PI0_VALIDATED_BASE:-${REPO}/output/pi0_bootstrap_20261004}"
OUT="${NPS_PI0_STUDY_OUTPUT:-${REPO}/output/pi0_uncertainty_reproduction}"
CFG="${REPO}/src/xsec_extract/xsec_config/xsec_config_x36_4.json"
MANIFEST="${REPO}/validation/pi0_bootstrap_20261004/run_quality_manifest.csv"
CAMPAIGN="${BASE}/nominal/KinC_x36_4"
RESPONSE="${BASE}/selected_nominal"
LIB="${BASE}/build/libnps_stat.so"
SIM=/lustre24/expphy/volatile/hallc/nps/singhav/nps_smearing/smear_x36_4/smearing_output/KinC_x36_4/root/simc_pi0_analysis_output_smeared.root
VERTEX="${REPO}/output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x36_4.root"
case "${1:-help}" in
  data)
    [[ ! -e "$OUT/data_2000" ]] || { echo 'Refusing existing data ensemble' >&2; exit 1; }
    python3 scripts/bootstrap_pi0_data.py --campaign "$CAMPAIGN" --output "$OUT/data_2000" \
      --library "$LIB" --manifest "$MANIFEST" --response-dir "$BASE/response_reference" \
      --replicas 2000 --seed 20261004 --estimator selected ;;
  data-check)
    python3 scripts/summarize_pi0_data_2000.py --data "$OUT/data_2000" --previous "$BASE/bootstrap_500" \
      --response "$RESPONSE" --campaign "$CAMPAIGN" --manifest "$MANIFEST" --library "$LIB" --output "$OUT/data_validation"
    python3 scripts/diagnose_pi0_data_failures.py --data "$OUT/data_2000" --response "$RESPONSE" \
      --manifest "$MANIFEST" --campaign "$CAMPAIGN" --library "$LIB" --output "$OUT/data_failure_diagnosis" ;;
  mc)
    command -v root-config >/dev/null || { echo 'Load /group/nps/singhav/setup.csh first.' >&2; exit 1; }
    mkdir -p "$OUT/build"
    python3 src/xsec_extract/generate_xsec_config.py "$CFG" "$OUT/build/xsec_config.h"
    g++ -O2 -std=c++17 -shared -fPIC -I"$OUT/build" -Isrc/xsec_extract scripts/pi0_mc_response_bridge.cpp \
      $(root-config --cflags --libs) -lMinuit2 -o "$OUT/build/libpi0_mc.so.tmp"
    mv "$OUT/build/libpi0_mc.so.tmp" "$OUT/build/libpi0_mc.so"
    python3 scripts/bootstrap_pi0_mc.py --sim "$SIM" --vertex "$VERTEX" --response "$RESPONSE" \
      --library "$LIB" --mc-library "$OUT/build/libpi0_mc.so" --config "$CFG" \
      --first-order "$BASE/statistical_scales/finite_mc_delta_sigma_covariance.csv" --output "$OUT/mc_2000" --replicas 2000
    python3 scripts/complete_pi0_mc_rank_study.py --mc "$OUT/mc_2000" --response "$RESPONSE" --config "$CFG" \
      --first-order "$BASE/statistical_scales/finite_mc_delta_sigma_covariance.csv" --output "$OUT/mc_rank_aware"
    python3 scripts/diagnose_pi0_mc_delta.py --mc "$OUT/mc_2000" --rank-aware "$OUT/mc_rank_aware" \
      --response "$RESPONSE" --first-order "$BASE/statistical_scales/finite_mc_delta_sigma_covariance.csv" --output "$OUT/mc_diagnosis" ;;
  background)
    python3 scripts/validate_pi0_background_shape.py --campaign "$CAMPAIGN" --manifest "$MANIFEST" \
      --config "$CFG" --response "$RESPONSE" --data "$OUT/data_2000" --library "$LIB" --output "$OUT/background_shape" --toys 500 ;;
  coverage)
    python3 scripts/validate_pi0_coverage.py --library "$LIB" --output "$OUT/coverage" --experiments 300 --bootstrap 200
    python3 scripts/validate_pi0_coverage.py --library "$LIB" --output "$OUT/coverage_production_scale" \
      --experiments 300 --bootstrap 200 --signal-yield 3637.166666666666 --accidental-per-bin 0.3766666666666666
    python3 scripts/diagnose_pi0_fit_bias.py --output "$OUT/bias_diagnosis_no_peak_leak"
    python3 scripts/project_pi0_coverage.py --coverage "$OUT/coverage" --response "$RESPONSE" --library "$LIB" --parameter 9 --output "$OUT/coverage_xsec" ;;
  collect|extract)
    extra=()
    [[ "$1" != extract ]] || extra+=(--release)
    python3 scripts/finalize_pi0_uncertainty.py --data "$OUT/data_2000" --mc "$OUT/mc_rank_aware" \
      --response "$RESPONSE" --coverage "$OUT/coverage_production_scale" \
      --output "$OUT/uncertainty_diagnostics" "${extra[@]}" ;;
  *) echo 'Usage: bash scripts/run_pi0_uncertainty_studies.sh {data|data-check|mc|background|coverage|collect|extract}'
     echo 'extract is a release check; current controlled undercoverage blocks a new final extraction.' ;;
esac
