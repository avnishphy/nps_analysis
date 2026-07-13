#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIN_DRIVER="${SCRIPT_DIR}/run_parallel_nps_analysis_main.sh"

if [[ ! -x "${MAIN_DRIVER}" ]]; then
  chmod +x "${MAIN_DRIVER}" || true
fi

echo "[WARN] run_parallel_nps_analysis_wfpi0.sh is deprecated. Forwarding to run_parallel_nps_analysis_main.sh (waveform mode)."

exec "${MAIN_DRIVER}" --mode waveform --source waveform "$@"
