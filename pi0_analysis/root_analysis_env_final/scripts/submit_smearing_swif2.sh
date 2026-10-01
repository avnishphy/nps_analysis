#!/usr/bin/env bash
# Submit one SWIF2/OpenMP smearing job per selected NPS kinematic setting.
#
# The current smearing pipeline always rebuilds nominal simulation from the
# three raw Geant channel files.  Each farm job therefore stages the combined
# data ROOT file, exclusive/SIDIS/delta Geant ROOT files, and the matching
# three SIMC .hist files that provide normfac and Ngen.
#
# Example interactive preview followed by approved submission (run from root_analysis_env_final):
#   ./scripts/submit_smearing_swif2.sh --kin KinC_x60_4b --target LH2 --combined-file output/KinC_x60_4b/root/combined_branches_LH2.root --nx 7 --ny 8 --x-min -26 --x-max 26 --y-min -32 --y-max 32 --cores 8 --omp-threads 8 --workflow smear_x60_4b_test
# The script first runs the pipeline on the interactive farm only far enough
# to write the preview PDF. It submits the approved fit to SWIF2 only after y.
# If that preview already exists, the script asks before replacing it.
#
# Recommended workflow:
#   1. Run with --dry-run and inspect every resolved input path.
#   2. Run normally on an interactive farm node. The pipeline creates each
#      preview PDF locally and the submitter prints its path.
#   3. Inspect every PDF while the prompt waits, then answer y to create and
#      run the optimizer/final-production workflow in SWIF2.

set -euo pipefail

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
DEFAULT_ROOT="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final"
PIPELINE_REL="src/simulation_smearing/run_smearing_pipeline.sh"
METADATA_REL="config/nps_simulation_kinematics.csv"

die() { echo "[ERROR] $*" >&2; exit 1; }
need_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }
safe_name() { printf '%s' "$1" | sed 's/[^[:alnum:]_-]/_/g'; }

# Worker mode: all files named below are SWIF-staged into private scratch.
if [[ "${1:-}" == "--run-one" ]]; then
  [[ $# -ge 11 ]] || die "Invalid --run-one invocation"
  runtime_tar="$2"; runtime_root_name="$3"; kin="$4"; result_tar="$5"
  approved="$6"; omp_threads="$7"; beam_energy="$8"; z_nps="$9"
  angle_flag="${10}"; angle_value="${11}"
  shift 11
  pipeline_options=("$@")

  [[ "$runtime_tar" != */* && "$result_tar" != */* ]] || die "Worker tar names must be basenames"
  [[ "$runtime_root_name" =~ ^[[:alnum:]_.-]+$ ]] || die "Unsafe runtime root name"
  [[ "$approved" == 0 || "$approved" == 1 ]] || die "Invalid approval mode"
  [[ "$angle_flag" == --nps-angle-deg || "$angle_flag" == --theta-nps-deg ]] || die "Invalid angle flag"
  [[ -s "$runtime_tar" && -s "${runtime_tar}.sha256" ]] || die "Runtime tar/checksum missing or empty"
  [[ -s combined_input.root ]] || die "Combined data input missing or empty"
  for f in sim_exclusive.root sim_sidis.root sim_delta.root \
           sim_exclusive.hist sim_sidis.hist sim_delta.hist; do
    [[ -s "$f" ]] || die "Staged simulation input missing or empty: $f"
  done

  # Reload through Hall C environment once. Re-entry file avoids csh quoting loss.
  if [[ "${NPS_FARM_ENV_LOADED:-0}" != 1 ]]; then
    reentry="${PWD}/.smearing_farm_reentry.sh"
    {
      printf '#!/usr/bin/env bash\nset -euo pipefail\nexport NPS_FARM_ENV_LOADED=1\nexec %q' "$SCRIPT_PATH"
      printf ' %q' --run-one "$runtime_tar" "$runtime_root_name" "$kin" "$result_tar" \
        "$approved" "$omp_threads" "$beam_energy" "$z_nps" "$angle_flag" "$angle_value" "${pipeline_options[@]}"
      printf '\n'
    } > "$reentry"
    chmod 700 "$reentry"
    export NPS_FARM_REENTRY="$reentry"
    exec csh -f -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; if ( $status != 0 ) exit 97; exec /bin/bash "$NPS_FARM_REENTRY"'
  fi

  need_cmd tar
  sha256sum -c "${runtime_tar}.sha256"
  tar -tzf "$runtime_tar" >/dev/null
  tar -tzf "$runtime_tar" | grep -x "${runtime_root_name}/${PIPELINE_REL}" >/dev/null || die "Runtime tar lacks pipeline"
  tar -xzf "$runtime_tar"

  runtime_root="${PWD}/${runtime_root_name}"
  pipeline="${runtime_root}/${PIPELINE_REL}"
  safe="$(safe_name "$kin")"
  output_base="${PWD}/smearing_output"
  root_dir="${output_base}/${safe}/root"
  mkdir -p "$root_dir"

  sim_file="${root_dir}/simc_pi0_analysis_output.root"
  sim_smeared="${root_dir}/simc_pi0_analysis_output_smeared.root"
  sim_args=(
    --sim-excl-input "${PWD}/sim_exclusive.root"
    --sim-sidis-input "${PWD}/sim_sidis.root"
    --sim-delta-input "${PWD}/sim_delta.root"
    --sim-excl-hist "${PWD}/sim_exclusive.hist"
    --sim-sidis-hist "${PWD}/sim_sidis.hist"
    --sim-delta-hist "${PWD}/sim_delta.hist"
  )

  export OMP_NUM_THREADS="$omp_threads"
  export OMP_DYNAMIC=FALSE
  export OMP_PROC_BIND=spread
  export OMP_PLACES=cores
  export RUN_TAG="farm_${safe}_$(date -u +%Y%m%dT%H%M%SZ)"

  command=(bash "$pipeline" --kin "$kin" --output-base "$output_base" --root-dir "$root_dir" \
    --combined-file "${PWD}/combined_input.root" --sim-file "$sim_file" --sim-smeared-file "$sim_smeared" \
    --out-file "${root_dir}/out_smear.root" --section-map-file "${root_dir}/section_map.csv" \
    --interp-file "${root_dir}/out_smear_interpolated.root" --input-preview-pdf "${root_dir}/smearing_input_histograms_prefit.pdf" \
    --acceptance-config "$([[ -f acceptance_config.conf ]] && echo "${PWD}/acceptance_config.conf" || echo "${runtime_root}/config/acceptance_cuts.conf")" \
    --beam-energy-gev "$beam_energy" --z-nps-cm "$z_nps" "$angle_flag" "$angle_value" \
    "${sim_args[@]}" "${pipeline_options[@]}")

  start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '[worker] command:'; printf ' %q' "${command[@]}"; printf '\n'
  if [[ "$approved" == 1 ]]; then printf 'y\n' | "${command[@]}"; else "${command[@]}" </dev/null; fi
  end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

  manifest="${PWD}/smearing_job_manifest.txt"
  {
    printf 'kin=%s\nmode=%s\nomp_threads=%s\nbeam_energy_gev=%s\nz_nps_cm=%s\nangle_flag=%s\nangle_value=%s\n' \
      "$kin" "$([[ "$approved" == 1 ]] && echo approved || echo preview)" "$omp_threads" "$beam_energy" "$z_nps" "$angle_flag" "$angle_value"
    printf 'start_utc=%s\nend_utc=%s\nruntime_sha256=%s\ncommand=' "$start_utc" "$end_utc" "$(sha256sum "$runtime_tar" | awk '{print $1}')"
    printf '%q ' "${command[@]}"; printf '\n'
  } > "$manifest"

  tmp_result="${result_tar}.tmp.$$"
  tar -czf "$tmp_result" smearing_output "$(basename "$manifest")"
  tar -tzf "$tmp_result" >/dev/null
  mv "$tmp_result" "$result_tar"
  echo "[worker] result tar ready: $result_tar"
  exit 0
fi

usage() {
  cat <<EOF
Usage: $(basename "$0") (--kin KIN [...] | --all-kins) [options]

Examples:
  # Resolve inputs, build/check the runtime tar, and print SWIF commands only.
  $(basename "$0") --kin KinC_x60_4b --workflow smear_x60_4b_check --dry-run

  # Generate the preview locally, inspect it, then answer y to submit to SWIF2.
  $(basename "$0") --kin KinC_x60_4b --workflow smear_x60_4b_fit

  # Generate the local preview and exit without creating a SWIF2 workflow.
  $(basename "$0") --kin KinC_x60_4b --workflow smear_x60_4b_preview --preview-only

Selection/input paths:
  --kin <Kin_old>             Select setting; repeatable
  --all-kins                  Select all entries from metadata CSV
  --metadata <csv>            Default: config/nps_simulation_kinematics.csv
  --target <name>             Default: LH2
  --input-base <path>         Default: <analysis-root>/output
  --combined-file <path>      Single-kin combined ROOT override
  --combined-template <path>  Multi-kin path; placeholders below
  --sim-input-base <path>     Raw Geant ROOT production base
  --sim-hist-base <path>      SIMC .hist normalization-file directory
  --sim-excl-input <path>     Single-kin raw exclusive ROOT override
  --sim-sidis-input <path>    Single-kin raw SIDIS ROOT override
  --sim-delta-input <path>    Single-kin raw delta ROOT override
  --sim-excl-template <path>  Multi-kin raw exclusive template
  --sim-sidis-template <path> Multi-kin raw SIDIS template
  --sim-delta-template <path> Multi-kin raw delta template
  --sim-excl-hist <path>      Single-kin exclusive .hist override
  --sim-sidis-hist <path>     Single-kin SIDIS .hist override
  --sim-delta-hist <path>     Single-kin delta .hist override
  --sim-excl-hist-template <p>  Multi-kin exclusive .hist template
  --sim-sidis-hist-template <p> Multi-kin SIDIS .hist template
  --sim-delta-hist-template <p> Multi-kin delta .hist template

Templates support {input_base}, {sim_input_base}, {sim_hist_base}, {kin},
{short}, and {target}. Here {short} removes the leading "KinC_".
Every selected kin must have one combined ROOT file, three Geant ROOT files,
and three SIMC .hist files. Preflight checks all inputs before SWIF is changed.

Smearing:
  --preview-only              Generate local preview PDF; do not submit SWIF2
  --preview-dir <path>        Default: /volatile/.../nps_smearing_previews/<workflow>
                              Existing previews require overwrite confirmation
  --data-tree/--sim-tree <n>  Defaults: physics/simulation
  --nx/--ny <N>               Defaults: 11/16
  --x-min/--x-max <v>         Defaults: -24/28 cm
  --y-min/--y-max <v>         Defaults: -34/34 cm
  --overlap <v>               Default: 0.0
  --nsmear <N>                Default: 80
  --smear-seed <N>            Default: 42
  --mode <auto|hcana|waveform> Default: auto
  --acceptance-config <path>  Custom config staged into each job
  --beam-energy-gev <v>       Global metadata override
  --z-nps-cm <v>              Global metadata override
  --nps-angle-deg <v>         Global physical-angle override
  --theta-nps-deg <v>         Global signed NPS-to-Hall override

SWIF2/parallel:
  --workflow <name>           Default: timestamped
  --account/--partition <v>   Defaults: hallc/production
  --cores <N>                 SWIF cores per kin; default: 8
  --omp-threads <N>           Fitter threads; default: cores; must be <= cores
  --ram/--disk/--time <v>     Defaults: 32g/100g/2d
  --result-dir/--tar-dir <p>  Output/runtime-tar directories
  --analysis-root <path>      Source tree packaged into runtime tar
  --create-only               Add jobs without starting workflow
  --dry-run                   Package/check/print; skip preview and SWIF mutation
  -h, --help                  Show help

Outputs:
  One <result-dir>/<kin>.tar.gz per job. Each tar contains smearing_output/
  plus smearing_job_manifest.txt. Logs go to /farm_out/${USER}/nps_smearing/<workflow>/.
EOF
}

ANALYSIS_ROOT="$DEFAULT_ROOT"
METADATA=""
TARGET=LH2
WORKFLOW="nps_smearing_$(date +%Y%m%d_%H%M%S)"
ACCOUNT=hallc; PARTITION=production
CORES=8; OMP_THREADS=""; RAM=32g; DISK=100g; WALLTIME=2d
RESULT_DIR=""; TAR_DIR="/group/nps/${USER}/swif_inputs"
INPUT_BASE=""
SIM_INPUT_BASE="/volatile/hallc/nps/singhav/geant4_simc/nps_geant4_20260824_221705"
SIM_HIST_BASE=""
PREVIEW_DIR=""
COMBINED_FILE=""; SIM_EXCL=""; SIM_SIDIS=""; SIM_DELTA=""
SIM_EXCL_HIST=""; SIM_SIDIS_HIST=""; SIM_DELTA_HIST=""
COMBINED_TEMPLATE='{input_base}/{kin}/root/combined_branches_{target}.root'
SIM_EXCL_TEMPLATE='{sim_input_base}/excl/nps_excl_pi0_{short}_geant4.root'
SIM_SIDIS_TEMPLATE='{sim_input_base}/sidis/nps_sidis_pi0_{short}_geant4.root'
SIM_DELTA_TEMPLATE='{sim_input_base}/delta/nps_delta_pi0_{short}_geant4.root'
SIM_EXCL_HIST_TEMPLATE='{sim_hist_base}/nps_excl_pi0_{short}.hist'
SIM_SIDIS_HIST_TEMPLATE='{sim_hist_base}/nps_sidis_pi0_{short}.hist'
SIM_DELTA_HIST_TEMPLATE='{sim_hist_base}/nps_delta_pi0_{short}.hist'
DATA_TREE=physics; SIM_TREE=simulation; NX=11; NY=16; X_MIN=-24; X_MAX=28; Y_MIN=-34; Y_MAX=34
OVERLAP=0.0; NSMEAR=80; SMEAR_SEED=42; MODE=auto
ACCEPTANCE_CONFIG=""; BEAM_OVERRIDE=""; Z_OVERRIDE=""; ANGLE_OVERRIDE=""; THETA_OVERRIDE=""
APPROVED=1; PREVIEW_ONLY=0; ALL_KINS=0; CREATE_ONLY=0; DRY_RUN=0; KINS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --kin) need_value "$@"; KINS+=("$2"); shift 2 ;;
    --all-kins) ALL_KINS=1; shift ;;
    --metadata) need_value "$@"; METADATA="$2"; shift 2 ;;
    --target) need_value "$@"; TARGET="$2"; shift 2 ;;
    --input-base) need_value "$@"; INPUT_BASE="$2"; shift 2 ;;
    --combined-file) need_value "$@"; COMBINED_FILE="$2"; shift 2 ;;
    --combined-template) need_value "$@"; COMBINED_TEMPLATE="$2"; shift 2 ;;
    --sim-file|--sim-file-template) die "$1 is unsupported: the current pipeline rebuilds nominal SIM from raw Geant inputs" ;;
    --sim-input-base) need_value "$@"; SIM_INPUT_BASE="$2"; shift 2 ;;
    --sim-hist-base) need_value "$@"; SIM_HIST_BASE="$2"; shift 2 ;;
    --sim-excl-input) need_value "$@"; SIM_EXCL="$2"; shift 2 ;;
    --sim-sidis-input) need_value "$@"; SIM_SIDIS="$2"; shift 2 ;;
    --sim-delta-input) need_value "$@"; SIM_DELTA="$2"; shift 2 ;;
    --sim-excl-template) need_value "$@"; SIM_EXCL_TEMPLATE="$2"; shift 2 ;;
    --sim-sidis-template) need_value "$@"; SIM_SIDIS_TEMPLATE="$2"; shift 2 ;;
    --sim-delta-template) need_value "$@"; SIM_DELTA_TEMPLATE="$2"; shift 2 ;;
    --sim-excl-hist) need_value "$@"; SIM_EXCL_HIST="$2"; shift 2 ;;
    --sim-sidis-hist) need_value "$@"; SIM_SIDIS_HIST="$2"; shift 2 ;;
    --sim-delta-hist) need_value "$@"; SIM_DELTA_HIST="$2"; shift 2 ;;
    --sim-excl-hist-template) need_value "$@"; SIM_EXCL_HIST_TEMPLATE="$2"; shift 2 ;;
    --sim-sidis-hist-template) need_value "$@"; SIM_SIDIS_HIST_TEMPLATE="$2"; shift 2 ;;
    --sim-delta-hist-template) need_value "$@"; SIM_DELTA_HIST_TEMPLATE="$2"; shift 2 ;;
    --inputs-approved) die "--inputs-approved is no longer needed; approval is requested after local preview generation" ;;
    --preview-only) PREVIEW_ONLY=1; shift ;;
    --preview-dir) need_value "$@"; PREVIEW_DIR="$2"; shift 2 ;;
    --data-tree) need_value "$@"; DATA_TREE="$2"; shift 2 ;;
    --sim-tree) need_value "$@"; SIM_TREE="$2"; shift 2 ;;
    --nx) need_value "$@"; NX="$2"; shift 2 ;;
    --ny) need_value "$@"; NY="$2"; shift 2 ;;
    --x-min) need_value "$@"; X_MIN="$2"; shift 2 ;;
    --x-max) need_value "$@"; X_MAX="$2"; shift 2 ;;
    --y-min) need_value "$@"; Y_MIN="$2"; shift 2 ;;
    --y-max) need_value "$@"; Y_MAX="$2"; shift 2 ;;
    --overlap) need_value "$@"; OVERLAP="$2"; shift 2 ;;
    --nsmear) need_value "$@"; NSMEAR="$2"; shift 2 ;;
    --smear-seed) need_value "$@"; SMEAR_SEED="$2"; shift 2 ;;
    --mode) need_value "$@"; MODE="$2"; shift 2 ;;
    --acceptance-config) need_value "$@"; ACCEPTANCE_CONFIG="$2"; shift 2 ;;
    --beam-energy-gev) need_value "$@"; BEAM_OVERRIDE="$2"; shift 2 ;;
    --z-nps-cm) need_value "$@"; Z_OVERRIDE="$2"; shift 2 ;;
    --nps-angle-deg) need_value "$@"; ANGLE_OVERRIDE="$2"; shift 2 ;;
    --theta-nps-deg) need_value "$@"; THETA_OVERRIDE="$2"; shift 2 ;;
    --workflow) need_value "$@"; WORKFLOW="$2"; shift 2 ;;
    --account) need_value "$@"; ACCOUNT="$2"; shift 2 ;;
    --partition) need_value "$@"; PARTITION="$2"; shift 2 ;;
    --cores) need_value "$@"; CORES="$2"; shift 2 ;;
    --omp-threads) need_value "$@"; OMP_THREADS="$2"; shift 2 ;;
    --ram) need_value "$@"; RAM="$2"; shift 2 ;;
    --disk) need_value "$@"; DISK="$2"; shift 2 ;;
    --time) need_value "$@"; WALLTIME="$2"; shift 2 ;;
    --result-dir) need_value "$@"; RESULT_DIR="$2"; shift 2 ;;
    --tar-dir) need_value "$@"; TAR_DIR="$2"; shift 2 ;;
    --analysis-root) need_value "$@"; ANALYSIS_ROOT="$2"; shift 2 ;;
    --create-only) CREATE_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (use --help)" ;;
  esac
done

[[ "$WORKFLOW" =~ ^[[:alnum:]][[:alnum:]_-]*$ ]] || die "Invalid workflow: $WORKFLOW"
[[ "$CORES" =~ ^[1-9][0-9]*$ ]] || die "--cores must be positive"
OMP_THREADS="${OMP_THREADS:-$CORES}"
[[ "$OMP_THREADS" =~ ^[1-9][0-9]*$ ]] || die "--omp-threads must be positive"
(( OMP_THREADS <= CORES )) || die "OMP threads ($OMP_THREADS) exceed SWIF cores ($CORES)"
[[ "$RAM" =~ ^[1-9][0-9]*[kmgt]?$ && "$DISK" =~ ^[1-9][0-9]*[kmgt]?$ ]] || die "Invalid RAM/disk size"
[[ "$WALLTIME" =~ ^[1-9][0-9]*[smhd]$ ]] || die "Invalid walltime"
[[ "$MODE" == auto || "$MODE" == hcana || "$MODE" == waveform ]] || die "Invalid mode: $MODE"
[[ -z "$ANGLE_OVERRIDE" || -z "$THETA_OVERRIDE" ]] || die "Use one angle override, not both"
for v in NX NY NSMEAR SMEAR_SEED; do [[ "${!v}" =~ ^[0-9]+$ ]] || die "$v must be non-negative"; done
(( NX > 0 && NY > 0 && NSMEAR > 0 )) || die "NX, NY, and NSMEAR must be positive"
(( ALL_KINS == 0 || ${#KINS[@]} == 0 )) || die "Use --all-kins or --kin, not both"
(( PREVIEW_ONLY == 0 || DRY_RUN == 0 )) || die "Use --preview-only or --dry-run, not both"

ANALYSIS_ROOT="$(readlink -f "$ANALYSIS_ROOT")"
METADATA="${METADATA:-${ANALYSIS_ROOT}/${METADATA_REL}}"
[[ -f "$METADATA" ]] || die "Metadata missing: $METADATA"
METADATA="$(readlink -f "$METADATA")"
INPUT_BASE="$(readlink -m "${INPUT_BASE:-${ANALYSIS_ROOT}/output}")"
SIM_INPUT_BASE="$(readlink -m "$SIM_INPUT_BASE")"
SIM_HIST_BASE="$(readlink -m "${SIM_HIST_BASE:-${ANALYSIS_ROOT}/output/simc/nps_simc_20260824_135058/outfiles/simc_gfortran_updated/outfiles}")"
PREVIEW_DIR="$(readlink -m "${PREVIEW_DIR:-/volatile/hallc/nps/${USER}/nps_smearing_previews/${WORKFLOW}}")"
if [[ -n "$ACCEPTANCE_CONFIG" ]]; then [[ -s "$ACCEPTANCE_CONFIG" ]] || die "Acceptance config missing or empty: $ACCEPTANCE_CONFIG"; ACCEPTANCE_CONFIG="$(readlink -f "$ACCEPTANCE_CONFIG")"; fi
[[ -s "${ANALYSIS_ROOT}/${PIPELINE_REL}" ]] || die "Pipeline missing or empty"
[[ -n "$ACCEPTANCE_CONFIG" || -s "${ANALYSIS_ROOT}/config/acceptance_cuts.conf" ]] || \
  die "Default acceptance config missing or empty: ${ANALYSIS_ROOT}/config/acceptance_cuts.conf"

declare -A META_BEAM=() META_Z=() META_ANGLE=() AVAILABLE=()
while IFS=$'\t' read -r k b a z; do
  [[ -n "$k" ]] || continue
  AVAILABLE["$k"]=1; META_BEAM["$k"]="$b"; META_ANGLE["$k"]="$a"; META_Z["$k"]="$z"
done < <(awk -F',' -v target="$TARGET" '
  function t(s){gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s}
  NR==1{for(i=1;i<=NF;i++){h=t($i);if(h=="kin_old")k=i;if(h=="target")q=i;if(h=="ebeam_gev")b=i;if(h=="nps_theta_deg")a=i;if(h=="nps_target_distance_cm")z=i}next}
  t($q)==target{print t($k) "\t" t($b) "\t" t($a) "\t" t($z)}
' "$METADATA")
(( ${#AVAILABLE[@]} > 0 )) || die "No metadata rows for target=$TARGET"

if (( ALL_KINS == 1 )); then mapfile -t KINS < <(printf '%s\n' "${!AVAILABLE[@]}" | sort); fi
(( ${#KINS[@]} > 0 )) || die "Select --kin or --all-kins"
if (( ${#KINS[@]} > 1 )) && [[ -n "$COMBINED_FILE$SIM_EXCL$SIM_SIDIS$SIM_DELTA$SIM_EXCL_HIST$SIM_SIDIS_HIST$SIM_DELTA_HIST" ]]; then
  die "Single-file input overrides require exactly one kin; use templates for multiple/all"
fi

expand_template() {
  local s="$1" kin="$2" short="${2#KinC_}"
  s="${s//\{input_base\}/$INPUT_BASE}"; s="${s//\{sim_input_base\}/$SIM_INPUT_BASE}"
  s="${s//\{sim_hist_base\}/$SIM_HIST_BASE}"
  s="${s//\{kin\}/$kin}"; s="${s//\{short\}/$short}"; s="${s//\{target\}/$TARGET}"
  printf '%s' "$s"
}
is_float() { [[ "$1" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; }
require_swif_readable() {
  local path="$1" label="$2" perms
  perms="$(stat -Lc %A "$path")"
  [[ "${perms:7:1}" == r ]] || \
    die "$label is not world-readable for SWIF2 staging: $path (permissions: $perms)"
}

declare -A JOB_COMBINED=() JOB_EXCL=() JOB_SIDIS=() JOB_DELTA=()
declare -A JOB_EXCL_HIST=() JOB_SIDIS_HIST=() JOB_DELTA_HIST=()
declare -A JOB_BEAM=() JOB_Z=() JOB_ANGLE=() JOB_ANGLE_FLAG=() SEEN=()
selected=(); missing=0
for kin in "${KINS[@]}"; do
  [[ -n "${AVAILABLE[$kin]+x}" ]] || die "Unknown kin for target=$TARGET: $kin"
  [[ -n "${SEEN[$kin]+x}" ]] && continue; SEEN["$kin"]=1; selected+=("$kin")
  combined="${COMBINED_FILE:-$(expand_template "$COMBINED_TEMPLATE" "$kin")}"; [[ -s "$combined" ]] && combined="$(readlink -f "$combined")"
  JOB_COMBINED["$kin"]="$combined"
  beam="${BEAM_OVERRIDE:-${META_BEAM[$kin]}}"; z="${Z_OVERRIDE:-${META_Z[$kin]}}"
  if [[ -n "$THETA_OVERRIDE" ]]; then angle="$THETA_OVERRIDE"; JOB_ANGLE_FLAG["$kin"]=--theta-nps-deg; else angle="${ANGLE_OVERRIDE:-${META_ANGLE[$kin]}}"; JOB_ANGLE_FLAG["$kin"]=--nps-angle-deg; fi
  JOB_BEAM["$kin"]="$beam"; JOB_Z["$kin"]="$z"; JOB_ANGLE["$kin"]="$angle"
  is_float "$beam" && is_float "$z" && is_float "$angle" || die "Invalid beam/z/angle metadata for $kin: $beam/$z/$angle"
  if [[ ! -s "$combined" ]]; then echo "[MISSING/EMPTY] combined kin=$kin: $combined" >&2; missing=1; fi

  excl="${SIM_EXCL:-$(expand_template "$SIM_EXCL_TEMPLATE" "$kin")}"; sidis="${SIM_SIDIS:-$(expand_template "$SIM_SIDIS_TEMPLATE" "$kin")}"; delta="${SIM_DELTA:-$(expand_template "$SIM_DELTA_TEMPLATE" "$kin")}";
  excl_hist="${SIM_EXCL_HIST:-$(expand_template "$SIM_EXCL_HIST_TEMPLATE" "$kin")}"; sidis_hist="${SIM_SIDIS_HIST:-$(expand_template "$SIM_SIDIS_HIST_TEMPLATE" "$kin")}"; delta_hist="${SIM_DELTA_HIST:-$(expand_template "$SIM_DELTA_HIST_TEMPLATE" "$kin")}";
  [[ -s "$excl" ]] && excl="$(readlink -f "$excl")"; [[ -s "$sidis" ]] && sidis="$(readlink -f "$sidis")"; [[ -s "$delta" ]] && delta="$(readlink -f "$delta")"
  [[ -s "$excl_hist" ]] && excl_hist="$(readlink -f "$excl_hist")"; [[ -s "$sidis_hist" ]] && sidis_hist="$(readlink -f "$sidis_hist")"; [[ -s "$delta_hist" ]] && delta_hist="$(readlink -f "$delta_hist")"
  JOB_EXCL["$kin"]="$excl"; JOB_SIDIS["$kin"]="$sidis"; JOB_DELTA["$kin"]="$delta"
  JOB_EXCL_HIST["$kin"]="$excl_hist"; JOB_SIDIS_HIST["$kin"]="$sidis_hist"; JOB_DELTA_HIST["$kin"]="$delta_hist"
  for spec in "exclusive Geant ROOT:$excl" "SIDIS Geant ROOT:$sidis" "delta Geant ROOT:$delta" \
              "exclusive SIMC hist:$excl_hist" "SIDIS SIMC hist:$sidis_hist" "delta SIMC hist:$delta_hist"; do
    [[ -s "${spec#*:}" ]] || { echo "[MISSING/EMPTY] ${spec%%:*} kin=$kin: ${spec#*:}" >&2; missing=1; }
  done
done
KINS=("${selected[@]}")
(( missing == 0 )) || die "Input preflight failed; no tar/workflow created"
if (( PREVIEW_ONLY == 0 )); then
  for kin in "${KINS[@]}"; do
    require_swif_readable "${JOB_COMBINED[$kin]}" "Combined input for $kin"
    require_swif_readable "${JOB_EXCL[$kin]}" "Exclusive Geant input for $kin"
    require_swif_readable "${JOB_SIDIS[$kin]}" "SIDIS Geant input for $kin"
    require_swif_readable "${JOB_DELTA[$kin]}" "Delta Geant input for $kin"
    require_swif_readable "${JOB_EXCL_HIST[$kin]}" "Exclusive SIMC hist for $kin"
    require_swif_readable "${JOB_SIDIS_HIST[$kin]}" "SIDIS SIMC hist for $kin"
    require_swif_readable "${JOB_DELTA_HIST[$kin]}" "Delta SIMC hist for $kin"
  done
  [[ -z "$ACCEPTANCE_CONFIG" ]] || require_swif_readable "$ACCEPTANCE_CONFIG" "Acceptance config"
fi

# Keep all generated preview artifacts outside the analysis inputs. The
# pipeline receives EOF at its own confirmation prompt, so it publishes only
# smearing_input_histograms_prefit.pdf and never starts local optimization.
common_pipeline=(--target "$TARGET" --data-tree "$DATA_TREE" --sim-tree "$SIM_TREE" --nx "$NX" --ny "$NY" \
  --x-min "$X_MIN" --x-max "$X_MAX" --y-min "$Y_MIN" --y-max "$Y_MAX" --overlap "$OVERLAP" \
  --nsmear "$NSMEAR" --smear-seed "$SMEAR_SEED" --mode "$MODE")
declare -A JOB_PREVIEW_PDF=()
tmp_tar=""; worker_tmp=""; preview_runner=""; preview_backup=""; preview_replaced_path=""
cleanup() {
  if [[ -n "${preview_runner:-}" && -f "$preview_runner" ]]; then rm -f -- "$preview_runner" || true; fi
  if [[ -n "${tmp_tar:-}" && -f "$tmp_tar" ]]; then rm -f -- "$tmp_tar" || true; fi
  if [[ -n "${worker_tmp:-}" && -f "$worker_tmp" ]]; then rm -f -- "$worker_tmp" || true; fi
  if [[ -n "${preview_backup:-}" && ( -e "$preview_backup" || -L "$preview_backup" ) ]]; then
    if [[ -n "${preview_replaced_path:-}" ]]; then
      rm -f -- "$preview_replaced_path" || true
      mv -f -- "$preview_backup" "$preview_replaced_path" || true
    else
      rm -f -- "$preview_backup" || true
    fi
  fi
  return 0
}
trap cleanup EXIT

generate_local_preview() {
  local kin="$1" safe preview_root preview_pdf acceptance_config preview_status backup_candidate
  local -a preview_command
  safe="$(safe_name "$kin")"
  preview_root="${PREVIEW_DIR}/${safe}/root"
  preview_pdf="${preview_root}/smearing_input_histograms_prefit.pdf"
  if [[ -e "$preview_pdf" || -L "$preview_pdf" ]]; then
    [[ -f "$preview_pdf" || -L "$preview_pdf" ]] || die "Preview path exists but is not a file: $preview_pdf"
    echo "[preview] Existing PDF: $preview_pdf"
    [[ -t 0 ]] || die "Preview overwrite confirmation requires an interactive terminal"
    printf 'Overwrite and regenerate this input preview PDF? [y/N] '
    IFS= read -r overwrite_response || overwrite_response=""
    case "$overwrite_response" in
      y|Y|yes|Yes|YES)
        backup_candidate="$(mktemp "${preview_pdf}.backup.XXXXXX")"
        if ! mv -f -- "$preview_pdf" "$backup_candidate"; then
          rm -f -- "$backup_candidate"
          die "Could not move the existing preview aside: $preview_pdf"
        fi
        preview_backup="$backup_candidate"
        preview_replaced_path="$preview_pdf"
        echo "[preview] Existing PDF moved aside; regenerating it from current inputs."
        ;;
      *)
        echo "Existing preview preserved; submission cancelled before SWIF2 workflow creation."
        exit 0
        ;;
    esac
  fi
  mkdir -p "$preview_root"
  acceptance_config="${ACCEPTANCE_CONFIG:-${ANALYSIS_ROOT}/config/acceptance_cuts.conf}"
  preview_command=(bash "${ANALYSIS_ROOT}/${PIPELINE_REL}" --kin "$kin" \
    --output-base "$PREVIEW_DIR" --root-dir "$preview_root" \
    --combined-file "${JOB_COMBINED[$kin]}" \
    --sim-file "${preview_root}/simc_pi0_analysis_output.root" \
    --sim-smeared-file "${preview_root}/simc_pi0_analysis_output_smeared.root" \
    --out-file "${preview_root}/out_smear.root" \
    --section-map-file "${preview_root}/section_map.csv" \
    --interp-file "${preview_root}/out_smear_interpolated.root" \
    --input-preview-pdf "$preview_pdf" \
    --acceptance-config "$acceptance_config" \
    --beam-energy-gev "${JOB_BEAM[$kin]}" --z-nps-cm "${JOB_Z[$kin]}" \
    "${JOB_ANGLE_FLAG[$kin]}" "${JOB_ANGLE[$kin]}" \
    --sim-excl-input "${JOB_EXCL[$kin]}" --sim-sidis-input "${JOB_SIDIS[$kin]}" \
    --sim-delta-input "${JOB_DELTA[$kin]}" --sim-excl-hist "${JOB_EXCL_HIST[$kin]}" \
    --sim-sidis-hist "${JOB_SIDIS_HIST[$kin]}" --sim-delta-hist "${JOB_DELTA_HIST[$kin]}" \
    "${common_pipeline[@]}")

  echo "[preview] Generating locally for $kin"
  printf '[preview] command:'; printf ' %q' "${preview_command[@]}"; printf '\n'
  preview_runner="$(mktemp "/tmp/nps_smearing_preview.${safe}.XXXXXX.sh")"
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\n'
    printf 'export OMP_NUM_THREADS=%q OMP_DYNAMIC=FALSE OMP_PROC_BIND=spread OMP_PLACES=cores\n' "$OMP_THREADS"
    printf 'exec'; printf ' %q' "${preview_command[@]}"; printf ' </dev/null\n'
  } > "$preview_runner"
  chmod 700 "$preview_runner"
  export NPS_PREVIEW_REENTRY="$preview_runner"
  set +e
  csh -f -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; if ( $status != 0 ) exit 97; exec /bin/bash "$NPS_PREVIEW_REENTRY"'
  preview_status=$?
  set -e
  rm -f -- "$preview_runner"; preview_runner=""; unset NPS_PREVIEW_REENTRY
  (( preview_status == 0 )) || die "Local preview pipeline failed for $kin with status $preview_status"
  [[ -s "$preview_pdf" ]] || die "Local preview PDF was not created: $preview_pdf"
  if [[ -n "$preview_backup" ]]; then
    rm -f -- "$preview_backup"
    preview_backup=""; preview_replaced_path=""
  fi
  JOB_PREVIEW_PDF["$kin"]="$preview_pdf"
  echo "[preview] ready: $preview_pdf"
}

if (( DRY_RUN == 0 )); then
  if (( PREVIEW_ONLY == 0 )); then
    [[ -t 0 ]] || die "Interactive preview and SWIF2 approval require a terminal; use --preview-only for preview generation without submission"
  fi
  for kin in "${KINS[@]}"; do generate_local_preview "$kin"; done
  echo ""
  echo "Local input-preview PDF(s):"
  for kin in "${KINS[@]}"; do echo "  $kin: ${JOB_PREVIEW_PDF[$kin]}"; done
  if (( PREVIEW_ONLY == 1 )); then
    echo "Preview-only complete; no SWIF2 workflow was created."
    exit 0
  fi

  echo ""
  echo "==== Input preview verification required ===="
  echo "Selected kinematics: ${KINS[*]}"
  echo "Review every PDF listed above while this prompt waits."
  echo "A yes response will submit one pipeline job per kinematic setting to SWIF2."
  printf 'Have all input preview PDFs been reviewed and approved? Submit to SWIF2? [y/N] '
  IFS= read -r approval_response || approval_response=""
  case "$approval_response" in
    y|Y|yes|Yes|YES)
      echo "Input previews approved; continuing with SWIF2 submission."
      ;;
    *)
      echo "Submission cancelled before SWIF2 workflow creation."
      exit 0
      ;;
  esac
fi

RESULT_DIR="${RESULT_DIR:-/volatile/hallc/nps/${USER}/nps_smearing/${WORKFLOW}}"
LOG_DIR="/farm_out/${USER}/nps_smearing/${WORKFLOW}"
runtime_root_name="$(basename "$ANALYSIS_ROOT")"
runtime_tar="${TAR_DIR}/${WORKFLOW}_smearing_runtime.tar.gz"; runtime_sha="${runtime_tar}.sha256"
worker_basename="${WORKFLOW}_smearing_worker.sh"
worker_script="${TAR_DIR}/${worker_basename}"
mkdir -p "$TAR_DIR"; (( DRY_RUN == 1 )) || mkdir -p "$RESULT_DIR" "$LOG_DIR"
[[ ! -e "$runtime_tar" && ! -e "$runtime_sha" && ! -e "$worker_script" ]] || \
  die "Runtime package/worker exists for workflow: $WORKFLOW"
for kin in "${KINS[@]}"; do result="${RESULT_DIR}/$(safe_name "$kin").tar.gz"; [[ ! -e "$result" ]] || die "Result exists: $result"; done

need_cmd tar; need_cmd sha256sum
tmp_tar="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.tar.gz")"
tar -C "$(dirname "$ANALYSIS_ROOT")" \
  --exclude="${runtime_root_name}/src/analysis/*.so" \
  --exclude="${runtime_root_name}/src/analysis/*.d" \
  --exclude="${runtime_root_name}/src/analysis/*.pcm" \
  --exclude="${runtime_root_name}/src/analysis/__pycache__" \
  --exclude="${runtime_root_name}/src/simulation_smearing/*.so" \
  --exclude="${runtime_root_name}/src/simulation_smearing/*.d" \
  --exclude="${runtime_root_name}/src/simulation_smearing/*.pcm" \
  --exclude="${runtime_root_name}/src/simulation_smearing/farm_jobs" \
  --exclude="${runtime_root_name}/src/simulation_smearing/submit*swif2.sh" \
  -czf "$tmp_tar" "${runtime_root_name}/config" "${runtime_root_name}/src/analysis" "${runtime_root_name}/src/simulation_smearing"
tar -tzf "$tmp_tar" >/dev/null
tar -tzf "$tmp_tar" | grep -x "${runtime_root_name}/${PIPELINE_REL}" >/dev/null || die "Generated tar lacks pipeline"
mv "$tmp_tar" "$runtime_tar"; tmp_tar=""
runtime_hash="$(sha256sum "$runtime_tar" | awk '{print $1}')"; printf '%s  analysis_runtime.tar.gz\n' "$runtime_hash" > "$runtime_sha"
# SWIF2's input-staging service reads these as a separate account. mktemp
# creates mode 0600, so publish both package files explicitly before SWIF use.
chmod 0644 "$runtime_tar" "$runtime_sha"

# Freeze the exact worker implementation under a workflow-unique path/name.
# This follows submit_nps_analysis_swif2.sh and avoids stale SWIF input caching.
worker_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.worker.sh")"
cp "$SCRIPT_PATH" "$worker_tmp"
chmod 0755 "$worker_tmp"
mv "$worker_tmp" "$worker_script"; worker_tmp=""
worker_hash="$(sha256sum "$worker_script" | awk '{print $1}')"

echo "Workflow: $WORKFLOW"; echo "Kinematics: ${KINS[*]}"; echo "Mode: $([[ "$APPROVED" == 1 ]] && echo approved || echo preview)"
echo "Runtime: $runtime_tar"; echo "SHA256: $runtime_hash"; echo "Results: $RESULT_DIR"
echo "Worker: $worker_script"; echo "Worker SHA256: $worker_hash"
echo "Parallel: ${#KINS[@]} independent SWIF job(s); cores=$CORES, OMP threads/job=$OMP_THREADS"
echo "Inputs/job: combined data + 3 Geant ROOT + 3 SIMC hist files"

swif() { if (( DRY_RUN == 1 )); then printf '[dry-run]'; printf ' %q' swif2 "$@"; printf '\n'; else swif2 "$@"; fi; }
(( DRY_RUN == 1 )) || need_cmd swif2
swif create "$WORKFLOW"
for kin in "${KINS[@]}"; do
  safe="$(safe_name "$kin")"; result="${RESULT_DIR}/${safe}.tar.gz"
  inputs=(-input analysis_runtime.tar.gz "file:${runtime_tar}" -input analysis_runtime.tar.gz.sha256 "file:${runtime_sha}" \
    -input "$worker_basename" "file:${worker_script}" -input combined_input.root "file:${JOB_COMBINED[$kin]}" \
    -input sim_exclusive.root "file:${JOB_EXCL[$kin]}" -input sim_sidis.root "file:${JOB_SIDIS[$kin]}" \
    -input sim_delta.root "file:${JOB_DELTA[$kin]}" -input sim_exclusive.hist "file:${JOB_EXCL_HIST[$kin]}" \
    -input sim_sidis.hist "file:${JOB_SIDIS_HIST[$kin]}" -input sim_delta.hist "file:${JOB_DELTA_HIST[$kin]}")
  [[ -n "$ACCEPTANCE_CONFIG" ]] && inputs+=(-input acceptance_config.conf "file:${ACCEPTANCE_CONFIG}")
  swif add-job "$WORKFLOW" -name "smear_${safe}" -account "$ACCOUNT" -partition "$PARTITION" \
    -cores "$CORES" -ram "$RAM" -disk "$DISK" -time "$WALLTIME" \
    -stdout "${LOG_DIR}/smear_${safe}.out" -stderr "${LOG_DIR}/smear_${safe}.err" \
    "${inputs[@]}" -output results.tar.gz "file:${result}" \
    /bin/bash "$worker_basename" --run-one analysis_runtime.tar.gz "$runtime_root_name" "$kin" results.tar.gz \
      "$APPROVED" "$OMP_THREADS" "${JOB_BEAM[$kin]}" "${JOB_Z[$kin]}" \
      "${JOB_ANGLE_FLAG[$kin]}" "${JOB_ANGLE[$kin]}" "${common_pipeline[@]}"
done
if (( CREATE_ONLY == 0 )); then swif run "$WORKFLOW"; fi
if (( DRY_RUN == 1 )); then echo "Dry run complete; no SWIF state changed."; elif (( CREATE_ONLY == 1 )); then echo "Created only. Start: swif2 run $WORKFLOW"; else echo "Submitted. Status: swif2 status $WORKFLOW"; fi
