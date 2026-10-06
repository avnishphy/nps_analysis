#!/usr/bin/env bash
# Prerequisite: Hall C ROOT environment loaded through setup.csh.
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_dir="${1:-$(mktemp -d /tmp/nps-stat-safety.XXXXXX)}"
mkdir -p "${test_dir}"
test_dir="$(cd "${test_dir}" && pwd)"
cd "${repo_dir}"
root -l -b -q tests/test_nps_timing_geometry.C > "${test_dir}/geometry.log" 2>&1
root -l -b -q tests/test_nps_signed_timing.C > "${test_dir}/signed_timing.log" 2>&1
root -l -b -q tests/test_nps_zero_boundary.C > "${test_dir}/zero_boundary.log" 2>&1
root -l -b -q 'tests/test_nps_fit_status.C("output/KinC_x36_4/plots/run_6418/combbg_run6418_order4_results.root")' > "${test_dir}/fit_status.log" 2>&1
python3 tests/test_nps_run_status.py > "${test_dir}/run_status.log" 2>&1
python3 src/xsec_extract/generate_xsec_config.py \
  src/xsec_extract/xsec_config/xsec_config_x36_4.json "${test_dir}/xsec_config.h"
g++ -std=c++17 -O1 -I"${test_dir}" -Isrc/xsec_extract \
  tests/test_xsec_signed_accumulation.cpp $(root-config --cflags --libs) \
  -lMinuit2 -o "${test_dir}/signed_accumulation"
"${test_dir}/signed_accumulation" > "${test_dir}/signed_accumulation.log" 2>&1
g++ -std=c++17 -O1 -I"${test_dir}" -Isrc/xsec_extract \
  tests/test_xsec_linear_solver.C $(root-config --cflags --libs) \
  -lMinuit2 -o "${test_dir}/linear_solver"
"${test_dir}/linear_solver" > "${test_dir}/linear_solver.log" 2>&1
g++ -std=c++17 -fsyntax-only -I"${test_dir}" -Isrc/xsec_extract \
  $(root-config --cflags) src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C
bash -n src/analysis/run_parallel_nps_analysis_main.sh
python3 -m py_compile src/analysis/combine_analysis_branches.py
echo "PASS safety tests and extractor compilation; evidence: ${test_dir}"
