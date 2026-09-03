#!/usr/bin/env bash
set -euo pipefail

# Shared SIMC installation used to build each workflow runtime package.
SIMC_DIR="/u/group/nps/singhav/simc_gfortran_updated"
SIMC_NAME="$(basename "${SIMC_DIR}")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SELF="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"

# Worker mode: SWIF invokes this branch once per staged infile.
if [[ "${1:-}" == "--run-one" ]]; then
  job_dir="${PWD}"
  infile="${job_dir}/$2"
  runtime_tar="${job_dir}/$3"
  base="$(basename "${infile}" .inp)"
  # Unpack private SIMC runtime in job scratch and install this job's infile.
  tar -xzf "${runtime_tar}" -C "${job_dir}"
  work_dir="${job_dir}/${SIMC_NAME}"
  mkdir -p "${work_dir}"/{infiles,worksim,runout,outfiles}
  cp "${infile}" "${work_dir}/infiles/${base}.inp"
  cd "${work_dir}"
  ./run_simc_tree "${base}"
  exit 0
fi

# Submit mode: [workflow] [infile directory] [output directory].
workflow="${1:-nps_simc_$(date +%Y%m%d_%H%M%S)}"
[[ "${workflow}" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || {
  echo "Invalid workflow name: ${workflow}" >&2; exit 1;
}
infile_dir="$(readlink -f "${2:-${REPO_ROOT}/config/simc_infiles}")"
output_dir="$(readlink -m "${3:-${REPO_ROOT}/output/simc/${workflow}}")"
log_dir="/farm_out/${USER}/simc/${workflow}"
tar_dir="/group/nps/${USER}/swif_inputs"
runtime_tar="${tar_dir}/${workflow}_simc_runtime.tar.gz"
mkdir -p "${output_dir}"/{worksim,runout,outfiles} "${log_dir}" "${tar_dir}"
shopt -s nullglob
infiles=("${infile_dir}"/*.inp)
[[ ${#infiles[@]} -gt 0 ]] || { echo "No infiles: ${infile_dir}" >&2; exit 1; }

# Package runtime once; omit old/generated products from the archive.
echo "Creating SIMC runtime: ${runtime_tar}"
tar -C "$(dirname "${SIMC_DIR}")" \
  --exclude="${SIMC_NAME}/.git" \
  --exclude="${SIMC_NAME}/worksim/*" \
  --exclude="${SIMC_NAME}/runout/*" \
  --exclude="${SIMC_NAME}/outfiles/*" \
  -czf "${runtime_tar}" "${SIMC_NAME}"

swif2 create "${workflow}"
# One independent SWIF job per infile.
for infile in "${infiles[@]}"; do
  base="$(basename "${infile}" .inp)"
    # Exact mappings flatten SIMC's internal directory into workflow output.
    # Avoid match: here: it preserves the source path and duplicates directories.
  swif2 add-job "${workflow}" \
    -name "${base}" -account hallc -partition production \
    -cores 1 -ram 2500m -disk 3g -time 4h \
    -stdout "${log_dir}/${base}.out" -stderr "${log_dir}/${base}.err" \
    -input simc_runtime.tar.gz "file:${runtime_tar}" \
    -input "${base}.inp" "file:${infile}" \
    -input submit_simc_swif2.sh "file:${SELF}" \
    -output "${SIMC_NAME}/worksim/${base}.root" "file:${output_dir}/worksim/${base}.root" \
    -output "${SIMC_NAME}/runout/${base}.out" "file:${output_dir}/runout/${base}.out" \
    -output "${SIMC_NAME}/outfiles/${base}.gen" "file:${output_dir}/outfiles/${base}.gen" \
    -output "${SIMC_NAME}/outfiles/${base}.geni" "file:${output_dir}/outfiles/${base}.geni" \
    -output "${SIMC_NAME}/outfiles/${base}.hist" "file:${output_dir}/outfiles/${base}.hist" \
    -output "${SIMC_NAME}/outfiles/${base}_start_random_state.dat" "file:${output_dir}/outfiles/${base}_start_random_state.dat" \
    /bin/bash submit_simc_swif2.sh --run-one "${base}.inp" simc_runtime.tar.gz
done
swif2 run "${workflow}"
echo "Submitted ${#infiles[@]} jobs: ${workflow}"
echo "Outputs: ${output_dir}/{worksim,runout,outfiles}"
