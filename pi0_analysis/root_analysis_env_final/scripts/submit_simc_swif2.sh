#!/usr/bin/env bash
set -euo pipefail

# Shared SIMC installation used to build each workflow runtime package.
SIMC_DIR="/u/group/nps/singhav/simc_gfortran_updated"
SIMC_NAME="$(basename "${SIMC_DIR}")"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
SELF="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"

usage() {
  cat <<'EOF'
Usage:
  submit_simc_swif2.sh [workflow [infile_dir [output_dir]]]
  submit_simc_swif2.sh [--workflow NAME] [--infile-dir DIR]
                       [--infile FILE]... [--output-dir DIR] [--dry-run]

With no infile selection, all .inp files in config/simc_infiles are submitted.
Repeat --infile to submit specific .inp files from any directory. If both
--infile-dir and --infile are given, files from both selections are submitted.
--dry-run prints the selected jobs without packaging SIMC or calling SWIF2.

Examples:
  ./scripts/submit_simc_swif2.sh --workflow simc_test \
    --infile /path/to/custom_excl.inp --infile /path/to/custom_sidis.inp --dry-run
  ./scripts/submit_simc_swif2.sh simc_test /path/to/infile_directory
EOF
}

die() { echo "[ERROR] $*" >&2; exit 1; }

# Worker mode: SWIF invokes this branch once per staged infile.
if [[ "${1:-}" == "--run-one" ]]; then
  [[ $# -eq 3 ]] || die "--run-one requires an infile and runtime tar"
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

# Submit mode: select a directory, individual infiles, or both.
workflow_option=""
infile_dir_option=""
output_dir_option=""
dry_run=0
custom_infiles=()
positional=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help) usage; exit 0 ;;
    --workflow|--infile-dir|--infile|--output-dir)
      option="$1"
      [[ $# -ge 2 && -n "$2" && "$2" != --* ]] || die "${option} requires a value"
      case "$option" in
        --workflow) workflow_option="$2" ;;
        --infile-dir) infile_dir_option="$2" ;;
        --infile) custom_infiles+=("$2") ;;
        --output-dir) output_dir_option="$2" ;;
      esac
      shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -*) die "Unknown option: $1" ;;
    *) positional+=("$1"); shift ;;
  esac
done
[[ ${#positional[@]} -le 3 ]] || die "Too many positional arguments"
[[ -z "$workflow_option" || ${#positional[@]} -eq 0 ]] || die "Use either --workflow or a positional workflow"
[[ -z "$infile_dir_option" || ${#positional[@]} -lt 2 ]] || die "Use either --infile-dir or a positional infile directory"
[[ -z "$output_dir_option" || ${#positional[@]} -lt 3 ]] || die "Use either --output-dir or a positional output directory"

workflow="${workflow_option:-${positional[0]:-nps_simc_$(date +%Y%m%d_%H%M%S)}}"
[[ "$workflow" =~ ^[A-Za-z0-9][A-Za-z0-9_-]*$ ]] || die "Invalid workflow name: $workflow"
infile_dir="${infile_dir_option:-${positional[1]:-}}"
output_dir="$(readlink -m "${output_dir_option:-${positional[2]:-${REPO_ROOT}/output/simc/${workflow}}}")"

infiles=()
if [[ -n "$infile_dir" || ${#custom_infiles[@]} -eq 0 ]]; then
  infile_dir="${infile_dir:-${REPO_ROOT}/config/simc_infiles}"
  [[ -d "$infile_dir" ]] || die "Infile directory not found: $infile_dir"
  infile_dir="$(readlink -f "$infile_dir")"
  shopt -s nullglob
  infiles+=("${infile_dir}"/*.inp)
fi
for infile in "${custom_infiles[@]}"; do
  [[ -f "$infile" && -r "$infile" && "$infile" == *.inp ]] || die "Not a readable .inp file: $infile"
  infiles+=("$(readlink -f "$infile")")
done
[[ ${#infiles[@]} -gt 0 ]] || die "No .inp infiles selected"

# SWIF job names and output paths use the infile basename; collisions would
# otherwise overwrite one job's results with another's.
declare -A selected_bases=()
for infile in "${infiles[@]}"; do
  [[ -f "$infile" && -r "$infile" ]] || die "Not a readable .inp file: $infile"
  base="$(basename "$infile" .inp)"
  [[ "$base" =~ ^[A-Za-z0-9_][A-Za-z0-9_.-]*$ ]] || die "Unsafe infile basename: $base"
  [[ ! -v selected_bases[$base] ]] || die "Duplicate infile basename: $base"
  selected_bases["$base"]=1
done

echo "Selected ${#infiles[@]} SIMC infile(s) for workflow ${workflow}:"
printf '  %s\n' "${infiles[@]}"
echo "Outputs: ${output_dir}/{worksim,runout,outfiles}"
if (( dry_run )); then
  echo "Dry run complete; no SWIF state changed."
  exit 0
fi

log_dir="/farm_out/${USER}/simc/${workflow}"
tar_dir="/group/nps/${USER}/swif_inputs"
runtime_tar="${tar_dir}/${workflow}_simc_runtime.tar.gz"
mkdir -p "${output_dir}"/{worksim,runout,outfiles} "${log_dir}" "${tar_dir}"

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
