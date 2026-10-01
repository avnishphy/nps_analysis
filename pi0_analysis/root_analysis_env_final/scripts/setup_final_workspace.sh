#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
workspace=$(cd -- "${script_dir}/.." && pwd -P)

if [[ "${workspace}" == */root_analysis_env_main ]]; then
  printf 'error: refusing to create runtime directories in frozen MAIN: %s\n' "${workspace}" >&2
  exit 2
fi

mkdir -p -- \
  "${workspace}/build" \
  "${workspace}/logs" \
  "${workspace}/output" \
  "${workspace}/publication/generated" \
  "${workspace}/recovery" \
  "${workspace}/scratch" \
  "${workspace}/validation/runtime"

printf 'workspace directories ready: %s\n' "${workspace}"
