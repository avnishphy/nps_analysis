#!/usr/bin/env bash
# Prerequisite: source the Hall C/NPS ROOT environment. Writes only to scratch.
set -euo pipefail
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
test_output="${1:-$(mktemp -d /tmp/nps_xsec_recovery.XXXXXX)}"
mkdir -p "${test_output}"
read -r -a root_flags <<< "$(root-config --cflags --libs)"
g++ -O2 -std=c++17 "${repo_root}/tests/test_xsec_bin_recovery.C" "${root_flags[@]}" -o "${test_output}/test_recovery"
"${test_output}/test_recovery" "${test_output}"
python3 "${repo_root}/tests/validate_xsec_bin_recovery.py" "${test_output}"
