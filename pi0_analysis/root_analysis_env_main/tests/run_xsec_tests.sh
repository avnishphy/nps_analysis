#!/usr/bin/env bash
# Run after sourcing /group/nps/singhav/setup.csh (and Modules/init/csh when
# needed). Tests write only to the chosen scratch directory. Production inputs
# are opened read-only; no analysis or detector-production jobs are launched.
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_output="${1:-$(mktemp -d /tmp/nps_xsec_validation.XXXXXX)}"
fixture_root="${2:-/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main}"
mkdir -p "${test_output}"
read -r -a root_flags <<< "$(root-config --cflags --libs)"

# Algebra/units and numerical behavior are checked before production events.
g++ -O2 -std=c++17 "${repo_root}/tests/test_pi0_response_conventions.cpp" -o "${test_output}/test_gk"
"${test_output}/test_gk"
g++ -O2 -std=c++17 "${repo_root}/tests/test_xsec_helicity.cpp" -o "${test_output}/test_helicity"
"${test_output}/test_helicity"
g++ -O2 -std=c++17 "${repo_root}/tests/test_xsec_linear_solver.C" "${root_flags[@]}" -o "${test_output}/test_solver"
"${test_output}/test_solver"
g++ -O2 -std=c++17 "${repo_root}/tests/test_xsec_experimental_points.C" "${root_flags[@]}" -o "${test_output}/test_experimental_points"
"${test_output}/test_experimental_points" "${test_output}/failed_points" > "${test_output}/experimental_points.log"
rg 'PASS experimental points' "${test_output}/experimental_points.log"
g++ -O2 -std=c++17 "${repo_root}/src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C" "${root_flags[@]}" -o "${test_output}/extractor"
"${test_output}/extractor" --help > "${test_output}/help.log"
for obsolete in --normalize_mmiss --normalize-mmiss; do
    if "${test_output}/extractor" "${obsolete}" > "${test_output}/obsolete.log" 2>&1; then
        echo "Obsolete normalization option unexpectedly accepted" >&2
        exit 1
    fi
    rg -q 'removed' "${test_output}/obsolete.log"
done
bash -n "${repo_root}/src/xsec_extract/run_xsec_pipeline.sh"

sim="${fixture_root}/output/KinC_x60_4b/root/simc_pi0_analysis_output_smeared.root"
vertex="${fixture_root}/output/simc/nps_simc_20260824_135058/worksim/simc_gfortran_updated/worksim/nps_excl_pi0_x60_4b.root"
data="${fixture_root}/output/KinC_x60_4b/root/combined_branches_LH2.root"
if [[ ! -r "${sim}" || ! -r "${vertex}" || ! -r "${data}" ]]; then
    echo "PASS standalone tests; SKIP event closure: KinC_x60_4b fixtures unavailable"
    exit 0
fi

# Regression: reduced trees cannot silently fall back to reconstructed truth.
if "${test_output}/extractor" --data-file "${data}" --sim-file "${sim}" \
    --out-dir "${test_output}/missing_truth" --no-png --no-pdf --no-diagnostics \
    > "${test_output}/missing_truth.log" 2>&1; then
    echo "Missing vertex input unexpectedly accepted" >&2
    exit 1
fi
rg -q 'Migration requires' "${test_output}/missing_truth.log"

# Actual data are a smoke/reproducibility check, not a claim of physics closure.
"${test_output}/extractor" --data-file "${data}" --sim-file "${sim}" \
    --vertex_simc_file "${vertex}" --out-dir "${test_output}/current" \
    --no-png --no-pdf --no-diagnostics > "${test_output}/extraction.log" 2>&1
python3 "${repo_root}/tests/validate_migration.py" "${test_output}/current" \
    --json "${test_output}/closure.json" --sim-file "${sim}" --vertex-file "${vertex}"
echo "PASS; results: ${test_output}"
