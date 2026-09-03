#!/usr/bin/env bash
# Submit one SWIF2 workflow per NPS kinematic setting, one job per run, and a
# dependent finalizer job that merges summaries and packages the result.

set -euo pipefail

SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
DEFAULT_ANALYSIS_ROOT="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main"
RUNNER_REL="src/analysis/run_parallel_nps_analysis_main.sh"
CONFIG_REL="config/nps_dvcs_all_kins_main.csv"

die() { echo "[ERROR] $*" >&2; exit 1; }
need_value() { [[ $# -ge 2 && -n "$2" ]] || die "$1 requires a value"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Command not found: $1"; }
safe_name() { printf '%s' "$1" | sed 's/[^[:alnum:]_-]/_/g'; }
canonical_target() {
  case "${1,,}" in lh2) echo LH2;; ld2) echo LD2;; dummy) echo Dummy;; *) return 1;; esac
}

farm_reexec() {
  [[ "${NPS_FARM_ENV_LOADED:-0}" == 1 ]] && return 0
  local reentry="${PWD}/.nps_farm_reentry.sh"
  {
    printf '#!/usr/bin/env bash\nset -euo pipefail\nexport NPS_FARM_ENV_LOADED=1\nexec %q' "$SCRIPT_PATH"
    printf ' %q' "$@"
    printf '\n'
  } > "$reentry"
  chmod 700 "$reentry"
  export NPS_FARM_REENTRY="$reentry"
  exec csh -f -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; if ( $status != 0 ) exit 97; exec /bin/bash "$NPS_FARM_REENTRY"'
}

prepare_runtime() {
  local runtime_tar="$1" runtime_root_name="$2"
  require_command tar
  [[ -f "$runtime_tar" ]] || die "Staged runtime tar missing: $runtime_tar"
  [[ -f "${runtime_tar}.sha256" ]] || die "Staged runtime checksum missing: ${runtime_tar}.sha256"
  sha256sum -c "${runtime_tar}.sha256"
  tar -tzf "$runtime_tar" >/dev/null
  tar -tzf "$runtime_tar" | grep -x "${runtime_root_name}/${RUNNER_REL}" >/dev/null || \
    die "Runtime tar lacks ${runtime_root_name}/${RUNNER_REL}"
  tar -xzf "$runtime_tar"
}

# Per-run worker. Every job gets private SWIF scratch, so ROOT ACLiC artifacts
# cannot collide. User-facing outputs/logs are written directly to persistent
# storage; per-run jobs deliberately do not merge or tar shared directories.
if [[ "${1:-}" == "--run-one" ]]; then
  worker_args=("$@")
  [[ $# -ge 8 ]] || die "--run-one requires: TAR ROOT_NAME KIN RUN INPUT_ACCESS OUTPUT_BASE LOG_BASE [runner options]"
  runtime_tar="$2"
  runtime_root_name="$3"
  kin="$4"
  selected_run="$5"
  input_access="$6"
  output_base="$7"
  log_base="$8"
  shift 8
  runner_options=("$@")

  [[ "$runtime_tar" == */* ]] && die "Worker TAR must be a staged basename"
  [[ "$runtime_root_name" =~ ^[[:alnum:]_.-]+$ ]] || die "Unsafe runtime root name"
  [[ "$selected_run" =~ ^[0-9]+$ ]] || die "Invalid run: $selected_run"
  [[ "$input_access" == shared || "$input_access" == stage ]] || die "Invalid input access: $input_access"
  [[ "$output_base" == /volatile/* ]] || die "Worker output base must be under /volatile"
  [[ "$log_base" == /farm_out/* ]] || die "Worker log base must be under /farm_out"
  [[ -f input_manifest.tsv ]] || die "Staged input manifest missing"
  [[ -f run_source_report.csv ]] || die "Staged run-source report missing"

  farm_reexec "${worker_args[@]}"
  require_command "${ROOT_CMD:-root}"
  prepare_runtime "$runtime_tar" "$runtime_root_name"

  runtime_root="${PWD}/${runtime_root_name}"
  runner="${runtime_root}/${RUNNER_REL}"
  efficiency_source="${runtime_root}/output/efficiency_stuff"
  safe_kin="$(safe_name "$kin")"
  metadata_dir="${output_base}/${safe_kin}/metadata/run${selected_run}"
  mkdir -p "$metadata_dir" "${log_base}/${safe_kin}/logs"

  # Main macro has one non-overridable absolute charge-DB constant. Rebase only
  # extracted scratch copy; source tree and runtime tar stay unchanged.
  macro="${runtime_root}/src/analysis/nps_analysis_main.C"
  charge_db="${runtime_root}/config/DataBase_production_runs_newBCMOffset.txt"
  sed -i "s|^constexpr const char\* CHARGE_FALLBACK_DB = .*|constexpr const char* CHARGE_FALLBACK_DB = \"${charge_db}\";|" "$macro"
  export NPS_DEAD_BLOCK_CONFIG="${runtime_root}/config/dead_block_per_runs.csv"
  export NPS_SELECTION_REPORT_CSV="${efficiency_source}/selection_report_${kin}.csv"

  input_count=0
  input_bytes=0
  while IFS=$'\t' read -r run source_kind basename source_path size_bytes; do
    [[ "$run" == run ]] && continue
    [[ -n "$run" && -n "$basename" ]] || continue
    if [[ "$input_access" == stage ]]; then
      [[ -f "${PWD}/${basename}" ]] || die "Staged ROOT input missing: $basename"
    else
      [[ -r "$source_path" ]] || die "Shared ROOT input unreadable: $source_path"
    fi
    input_count=$((input_count + 1))
    input_bytes=$((input_bytes + size_bytes))
  done < input_manifest.tsv
  (( input_count > 0 )) || die "Input manifest contains no ROOT files"
  cp input_manifest.tsv "$metadata_dir/input_manifest.tsv"
  cp run_source_report.csv "$metadata_dir/run_source_report.csv"

  start_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  runtime_sha256="$(sha256sum "$runtime_tar" | awk '{print $1}')"
  runner_command=(bash "$runner" --output-base "$output_base" --log-base "$log_base" --kin "$kin" --run "$selected_run" --jobs 1 --run-only "${runner_options[@]}")
  [[ "$input_access" == stage ]] && runner_command+=(--input-dir "$PWD")
  printf '[worker] command:'
  printf ' %q' "${runner_command[@]}"
  printf '\n'
  runner_status=0
  "${runner_command[@]}" || runner_status=$?

  end_utc="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  manifest="${metadata_dir}/nps_job_manifest.txt"
  {
    printf 'kin=%s\n' "$kin"
    printf 'run=%s\n' "$selected_run"
    printf 'input_access=%s\ninput_files=%s\ninput_bytes=%s\n' "$input_access" "$input_count" "$input_bytes"
    printf 'start_utc=%s\n' "$start_utc"
    printf 'end_utc=%s\n' "$end_utc"
    printf 'runner_status=%s\n' "$runner_status"
    printf 'runtime_sha256=%s\n' "$runtime_sha256"
    printf 'command='
    printf '%q ' "${runner_command[@]}"
    printf '\n'
  } > "$manifest"

  echo "[worker] persistent analysis: ${output_base}/${safe_kin}"
  echo "[worker] persistent logs: ${log_base}/${safe_kin}/logs"
  exit "$runner_status"
fi

# Finalizer runs only after every run job in its kinematic workflow succeeds.
if [[ "${1:-}" == "--finalize-one" ]]; then
  worker_args=("$@")
  [[ $# -ge 7 ]] || die "--finalize-one requires: TAR ROOT_NAME KIN RESULT_ARCHIVE OUTPUT_BASE LOG_BASE [runner options]"
  runtime_tar="$2"
  runtime_root_name="$3"
  kin="$4"
  result_archive="$5"
  output_base="$6"
  log_base="$7"
  shift 7
  runner_options=("$@")

  [[ "$runtime_tar" == */* ]] && die "Worker TAR must be a staged basename"
  [[ "$runtime_root_name" =~ ^[[:alnum:]_.-]+$ ]] || die "Unsafe runtime root name"
  [[ "$result_archive" == /volatile/* ]] || die "Worker result archive must be under /volatile"
  [[ "$output_base" == /volatile/* ]] || die "Worker output base must be under /volatile"
  [[ "$log_base" == /farm_out/* ]] || die "Worker log base must be under /farm_out"
  [[ -f run_source_report.csv ]] || die "Staged run-source report missing"

  farm_reexec "${worker_args[@]}"
  prepare_runtime "$runtime_tar" "$runtime_root_name"

  runtime_root="${PWD}/${runtime_root_name}"
  runner="${runtime_root}/${RUNNER_REL}"
  efficiency_source="${runtime_root}/output/efficiency_stuff"
  safe_kin="$(safe_name "$kin")"
  metadata_dir="${output_base}/${safe_kin}/metadata/finalizer"
  mkdir -p "$metadata_dir" "${output_base}/efficiency_stuff" "${log_base}/${safe_kin}/logs"
  cp "${efficiency_source}/selection_report_${kin}.csv" \
     "${efficiency_source}/efficiency_${kin}.csv" \
     "${output_base}/efficiency_stuff/"
  cp run_source_report.csv "$metadata_dir/run_source_report.csv"

  runner_command=(bash "$runner" --output-base "$output_base" --log-base "$log_base" --kin "$kin" --finalize-only "${runner_options[@]}")
  printf '[finalizer] command:'
  printf ' %q' "${runner_command[@]}"
  printf '\n'
  runner_status=0
  "${runner_command[@]}" || runner_status=$?

  {
    printf 'kin=%s\nrunner_status=%s\nruntime_sha256=%s\n' \
      "$kin" "$runner_status" "$(sha256sum "$runtime_tar" | awk '{print $1}')"
    printf 'command='; printf '%q ' "${runner_command[@]}"; printf '\n'
  } > "${metadata_dir}/nps_finalize_manifest.txt"

  result_tmp="${result_archive}.tmp.$$"
  tar -C "$output_base" -czf "$result_tmp" "$safe_kin"
  tar -tzf "$result_tmp" >/dev/null
  mv "$result_tmp" "$result_archive"
  echo "[finalizer] result archive: $result_archive"
  exit "$runner_status"
fi

usage() {
  cat <<EOF
Usage: $(basename "$0") (--kin KIN... | --all-kins | --run RUN...) [options]

Selection/analysis:
  --kin <Kin_old...>       One/more settings; comma/space lists accepted
  --all-kins               Submit every setting matching --types
  --run <N...>             Restrict runs; space/comma list (also infers kins)
  --target <name...>       LH2, LD2, Dummy, or all; comma/space lists accepted
                           (default: all; matching is case-insensitive)
  --source <value>         updated, production, or waveform (default: auto-probe)
  --input-dir <path>       Explicit source ROOT directory
  --input-access <mode>    shared (default) or stage
  --missing-input <mode>   skip (default) or fail after writing source report
  --updated-dir <path>     Updated HCANA ROOT directory
  --production-dir <path>  Production HCANA ROOT directory
  --waveform-dir <path>    Waveform ROOT directory
  --mode <value>           auto, hcana, or waveform (default: auto)
  --types <csv>            Config Type filter (default: production,Production)
  --gevnum-cut <yes|no>    g.evnum cut choice (default: yes)
  --analysis-jobs <N>      Compatibility flag; must be 1 (one run/job)
  --timeout <seconds>      Per-run timeout; 0 disables (default: 0)
  --no-combine             Skip per-kin combine step
  --combine-target <name>  Restrict combine stage to one selected target

SWIF2/resources:
  --workflow <name>        Workflow name (default: timestamped)
  --account <name>         Default: hallc
  --partition <name>       Default: production
  --cores <N>              Cores per run job (default: 1)
  --ram <size>             RAM per run job (default: 4g)
  --disk <size>            Disk per run job (default: 20g)
  --time <duration>        Default: 12h
  --max-concurrent <N>     Maximum simultaneous run jobs per kin (default: 8)
  --output-dir <path>      Persistent analysis/archive directory under /volatile
  --tar-dir <path>         Runtime tar directory
  --analysis-root <path>   Source tree to package
  --create-only            Create workflow/jobs, but do not run workflow
  --dry-run                Build/verify tar and print SWIF2 commands only
  -h, --help               Show help

Each kinematic gets workflow <WORKFLOW>_<KIN>. Each run is an isolated job; a
dependent finalizer rebuilds summaries, runs combine, and writes <KIN>.tar.gz.
Analysis is written directly to /volatile and logs to /farm_out. SWIF scratch
holds staged inputs and the extracted runtime only.
The workflow-level <WORKFLOW>_run_sources.csv records ready/missing runs and
resolved waveform/updated/production paths; it is also copied into each result.
Runtime tar contains config/, src/analysis/, and top-level efficiency CSV inputs;
generated analysis data and efficiency plots are never packaged. Existing
runtime/result tar files are never overwritten.

Input modes:
  shared: jobs read preflighted ROOT files directly from source paths.
  stage:  ROOT files become SWIF -input entries; requested disk must exceed
          staged input size. Each result contains input_manifest.tsv.
EOF
}

ANALYSIS_ROOT="$DEFAULT_ANALYSIS_ROOT"
WORKFLOW="nps_analysis_$(date +%Y%m%d_%H%M%S)"
SUBMIT_USER="${NPS_SWIF_USER:-${LOGNAME:-${USER:-}}}"
[[ "$SUBMIT_USER" =~ ^[[:alpha:]_][[:alnum:]_-]*$ ]] || die "Cannot determine submit username; set NPS_SWIF_USER"
ACCOUNT="hallc"
PARTITION="production"
CORES=1
RAM="4g"
DISK="20g"
WALLTIME="12h"
MAX_CONCURRENT=8
OUTPUT_DIR=""
TAR_DIR="/group/nps/${SUBMIT_USER}/swif_inputs"
SOURCE=""
INPUT_DIR=""
INPUT_ACCESS="shared"
MISSING_INPUT="skip"
UPDATED_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated"
PRODUCTION_DIR="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production"
WAVEFORM_DIR="/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS"
MODE="auto"
TYPES_CSV="production,Production"
GEVNUM_CUT="yes"
ANALYSIS_JOBS=1
TIMEOUT_SEC=0
COMBINE=1
COMBINE_TARGET=""
ALL_KINS=0
CREATE_ONLY=0
DRY_RUN=0
KINS=()
RUNS=()
TARGETS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --kin)
      shift
      (( $# > 0 )) && [[ "$1" != --* ]] || die "--kin requires one or more settings"
      while (( $# > 0 )) && [[ "$1" != --* ]]; do
        IFS=',' read -r -a kin_tokens <<< "$1"; KINS+=("${kin_tokens[@]}"); shift
      done
      ;;
    --all-kins) ALL_KINS=1; shift ;;
    --run)
      shift
      (( $# > 0 )) && [[ "$1" != --* ]] || die "--run requires one or more run numbers"
      while (( $# > 0 )) && [[ "$1" != --* ]]; do
        IFS=',' read -r -a run_tokens <<< "$1"
        RUNS+=("${run_tokens[@]}")
        shift
      done
      ;;
    --target)
      shift
      (( $# > 0 )) && [[ "$1" != --* ]] || die "--target requires LH2, LD2, Dummy, or all"
      while (( $# > 0 )) && [[ "$1" != --* ]]; do
        IFS=',' read -r -a target_tokens <<< "$1"; TARGETS+=("${target_tokens[@]}"); shift
      done
      ;;
    --source) need_value "$@"; SOURCE="$2"; shift 2 ;;
    --input-dir) need_value "$@"; INPUT_DIR="$2"; shift 2 ;;
    --input-access) need_value "$@"; INPUT_ACCESS="$2"; shift 2 ;;
    --missing-input) need_value "$@"; MISSING_INPUT="$2"; shift 2 ;;
    --updated-dir) need_value "$@"; UPDATED_DIR="$2"; shift 2 ;;
    --production-dir) need_value "$@"; PRODUCTION_DIR="$2"; shift 2 ;;
    --waveform-dir) need_value "$@"; WAVEFORM_DIR="$2"; shift 2 ;;
    --mode) need_value "$@"; MODE="$2"; shift 2 ;;
    --types) need_value "$@"; TYPES_CSV="$2"; shift 2 ;;
    --gevnum-cut) need_value "$@"; GEVNUM_CUT="$2"; shift 2 ;;
    --analysis-jobs) need_value "$@"; ANALYSIS_JOBS="$2"; shift 2 ;;
    --timeout) need_value "$@"; TIMEOUT_SEC="$2"; shift 2 ;;
    --no-combine) COMBINE=0; shift ;;
    --combine-target) need_value "$@"; COMBINE_TARGET="$2"; shift 2 ;;
    --workflow) need_value "$@"; WORKFLOW="$2"; shift 2 ;;
    --account) need_value "$@"; ACCOUNT="$2"; shift 2 ;;
    --partition) need_value "$@"; PARTITION="$2"; shift 2 ;;
    --cores) need_value "$@"; CORES="$2"; shift 2 ;;
    --ram) need_value "$@"; RAM="$2"; shift 2 ;;
    --disk) need_value "$@"; DISK="$2"; shift 2 ;;
    --time) need_value "$@"; WALLTIME="$2"; shift 2 ;;
    --max-concurrent) need_value "$@"; MAX_CONCURRENT="$2"; shift 2 ;;
    --output-dir) need_value "$@"; OUTPUT_DIR="$2"; shift 2 ;;
    --tar-dir) need_value "$@"; TAR_DIR="$2"; shift 2 ;;
    --analysis-root) need_value "$@"; ANALYSIS_ROOT="$2"; shift 2 ;;
    --create-only) CREATE_ONLY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown option: $1 (use --help)" ;;
  esac
done

[[ "$WORKFLOW" =~ ^[[:alnum:]][[:alnum:]_-]*$ ]] || die "Invalid workflow name: $WORKFLOW"
[[ "$CORES" =~ ^[1-9][0-9]*$ ]] || die "--cores must be a positive integer"
[[ "$RAM" =~ ^[1-9][0-9]*[kmgt]?$ ]] || die "Invalid --ram size: $RAM"
[[ "$DISK" =~ ^[1-9][0-9]*[kmgt]?$ ]] || die "Invalid --disk size: $DISK"
[[ "$WALLTIME" =~ ^[1-9][0-9]*[smhd]$ ]] || die "Invalid --time duration: $WALLTIME"
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || die "--max-concurrent must be a positive integer"
[[ "$TIMEOUT_SEC" =~ ^[0-9]+$ ]] || die "--timeout must be a non-negative integer"
[[ "$GEVNUM_CUT" == yes || "$GEVNUM_CUT" == no ]] || die "--gevnum-cut must be yes or no"
[[ "$MODE" == auto || "$MODE" == hcana || "$MODE" == waveform ]] || die "Invalid --mode: $MODE"
[[ -z "$SOURCE" || "$SOURCE" == updated || "$SOURCE" == production || "$SOURCE" == waveform ]] || \
  die "Invalid --source: $SOURCE"
[[ "$INPUT_ACCESS" == shared || "$INPUT_ACCESS" == stage ]] || die "--input-access must be shared or stage"
[[ "$MISSING_INPUT" == skip || "$MISSING_INPUT" == fail ]] || die "--missing-input must be skip or fail"
[[ -z "$INPUT_DIR" || -d "$INPUT_DIR" ]] || die "Input directory missing: $INPUT_DIR"
(( ALL_KINS == 0 || ${#KINS[@]} == 0 )) || die "Use --all-kins or --kin, not both"

ANALYSIS_ROOT="$(readlink -f "$ANALYSIS_ROOT")"
runner_source="${ANALYSIS_ROOT}/${RUNNER_REL}"
config_source="${ANALYSIS_ROOT}/${CONFIG_REL}"
[[ -f "$runner_source" ]] || die "Runner missing: $runner_source"
[[ -f "$config_source" ]] || die "Config missing: $config_source"
[[ "$ANALYSIS_JOBS" =~ ^[1-9][0-9]*$ ]] || die "--analysis-jobs must be positive"
(( ANALYSIS_JOBS == 1 )) || die "--analysis-jobs must be 1: runs are now separate SWIF jobs"

# Normalize/deduplicate run numbers.
declare -A seen_runs=()
unique_runs=()
for run in "${RUNS[@]}"; do
  run="${run//[[:space:]]/}"
  [[ "$run" =~ ^[0-9]+$ ]] || die "Invalid run number: $run"
  if [[ -z "${seen_runs[$run]+x}" ]]; then seen_runs["$run"]=1; unique_runs+=("$run"); fi
done
RUNS=("${unique_runs[@]}")

# Canonical target list; "all" preserves the prior all-target analysis scope.
if (( ${#TARGETS[@]} == 0 )); then
  TARGETS=(LH2 LD2 Dummy)
else
  normalized_targets=(); select_all_targets=0
  declare -A seen_targets=()
  for raw_target in "${TARGETS[@]}"; do
    raw_target="${raw_target//[[:space:]]/}"
    [[ -n "$raw_target" ]] || continue
    if [[ "${raw_target,,}" == all ]]; then select_all_targets=1; continue; fi
    target="$(canonical_target "$raw_target")" || die "Invalid target: $raw_target (use LH2, LD2, Dummy, or all)"
    if [[ -z "${seen_targets[$target]+x}" ]]; then seen_targets["$target"]=1; normalized_targets+=("$target"); fi
  done
  if (( select_all_targets == 1 )); then TARGETS=(LH2 LD2 Dummy); else TARGETS=("${normalized_targets[@]}"); fi
fi
(( ${#TARGETS[@]} > 0 )) || die "No targets selected"
if [[ -n "$COMBINE_TARGET" ]]; then
  raw_combine_target="$COMBINE_TARGET"
  COMBINE_TARGET="$(canonical_target "$raw_combine_target")" || die "Invalid --combine-target: $raw_combine_target"
  [[ " ${TARGETS[*]} " == *" ${COMBINE_TARGET} "* ]] || die "--combine-target $COMBINE_TARGET is outside --target selection"
fi

# List config kinematics matching Type and optional run selection.
run_pattern="|"
for run in "${RUNS[@]}"; do run_pattern+="${run}|"; done
target_pattern="|"; for target in "${TARGETS[@]}"; do target_pattern+="${target,,}|"; done
mapfile -t available_kins < <(awk -F',' -v types="$TYPES_CSV" -v runs="$run_pattern" -v targets="$target_pattern" '
  function trim(s) { gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); gsub(/^"|"$/, "", s); return s }
  BEGIN { n=split(types,a,","); for(i=1;i<=n;i++) wanted[tolower(trim(a[i]))]=1 }
  NR==1 { for(i=1;i<=NF;i++){ h=trim($i); if(h=="Kin_old")kc=i; if(h=="run_number")rc=i; if(h=="Type")tc=i; if(tolower(h)=="target")gc=i } next }
  NR>1 && kc && rc && tc && gc {
    k=trim($kc); r=trim($rc); t=tolower(trim($tc)); g=tolower(trim($gc))
    if(k!="" && (t in wanted) && index(targets,"|" g "|") && (runs=="|" || index(runs,"|" r "|")>0)) seen[k]=1
  }
  END { for(k in seen) print k }
' "$config_source" | sort)
(( ${#available_kins[@]} > 0 )) || die "No kinematics match --types/--run selection"

if (( ALL_KINS == 1 )) || (( ${#KINS[@]} == 0 && ${#RUNS[@]} > 0 )); then
  KINS=("${available_kins[@]}")
elif (( ${#KINS[@]} == 0 )); then
  die "Select --kin, --all-kins, or --run"
fi

declare -A available=() seen_kins=() seen_safe=()
for kin in "${available_kins[@]}"; do available["$kin"]=1; done
selected_kins=()
for kin in "${KINS[@]}"; do
  kin="${kin//[[:space:]]/}"
  [[ -n "$kin" ]] || continue
  [[ -n "${available[$kin]+x}" ]] || die "Kinematic not found for selected types/runs: $kin"
  [[ -n "${seen_kins[$kin]+x}" ]] && continue
  seen_kins["$kin"]=1
  safe="$(safe_name "$kin")"
  [[ -z "${seen_safe[$safe]+x}" ]] || die "Kinematic job-name collision: $kin and ${seen_safe[$safe]}"
  seen_safe["$safe"]="$kin"
  selected_kins+=("$kin")
done
KINS=("${selected_kins[@]}")

# Resolve every ROOT file now. This mirrors runner priority and prevents jobs
# from discovering missing/wrong data only after dispatch.
shopt -s nullglob
RESOLVED_FILES=()
RESOLVED_KIND=""
RESOLVED_SOURCE=""
resolve_run_inputs() {
  local run="$1" dir="" preferred="auto" source_label=""
  local -a matches=()
  RESOLVED_FILES=(); RESOLVED_KIND=""; RESOLVED_SOURCE=""

  if [[ -n "$INPUT_DIR" ]]; then
    dir="$INPUT_DIR"; source_label="input-dir"
    if [[ "$MODE" == waveform || "$SOURCE" == waveform ]]; then preferred=waveform
    elif [[ "$MODE" == hcana || "$SOURCE" == updated || "$SOURCE" == production ]]; then preferred=hcana
    fi
  elif [[ -n "$SOURCE" ]]; then
    case "$SOURCE" in
      waveform) dir="$WAVEFORM_DIR"; preferred=waveform; source_label=waveform ;;
      updated) dir="$UPDATED_DIR"; preferred=hcana; source_label=updated ;;
      production) dir="$PRODUCTION_DIR"; preferred=hcana; source_label=production ;;
    esac
  fi

  try_waveform() {
    local d="$1" label="$2"; matches=("${d%/}"/nps_production_"${run}"_*_wf_calib.root)
    if (( ${#matches[@]} > 0 )); then RESOLVED_FILES=("${matches[@]}"); RESOLVED_KIND=waveform; RESOLVED_SOURCE="$label"; return 0; fi
    return 1
  }
  try_hcana() {
    local d="$1" label="$2" skim="${1%/}/skim_run${run}.root"
    if [[ -f "$skim" ]]; then RESOLVED_FILES=("$skim"); RESOLVED_KIND=hcana; RESOLVED_SOURCE="$label"; return 0; fi
    matches=("${d%/}"/nps_hms_coin_"${run}"_*_1_-1.root)
    if (( ${#matches[@]} > 0 )); then RESOLVED_FILES=("${matches[@]}"); RESOLVED_KIND=hcana; RESOLVED_SOURCE="$label"; return 0; fi
    return 1
  }

  if [[ -n "$dir" ]]; then
    case "$preferred" in
      waveform) try_waveform "$dir" "$source_label" && return 0 ;;
      hcana) try_hcana "$dir" "$source_label" && return 0 ;;
      auto) try_waveform "$dir" "$source_label" && return 0; try_hcana "$dir" "$source_label" && return 0 ;;
    esac
  else
    try_waveform "$WAVEFORM_DIR" waveform && return 0
    try_hcana "$UPDATED_DIR" updated && return 0
    try_hcana "$PRODUCTION_DIR" production && return 0
  fi
  return 1
}

kin_pattern="|"; for kin in "${KINS[@]}"; do kin_pattern+="${kin}|"; done
mapfile -t job_pairs < <(awk -F',' -v kins="$kin_pattern" -v types="$TYPES_CSV" -v runs="$run_pattern" -v targets="$target_pattern" '
  function trim(s){gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); gsub(/^"|"$/, "", s); return s}
  BEGIN{n=split(types,a,",");for(i=1;i<=n;i++)wanted[tolower(trim(a[i]))]=1}
  NR==1{for(i=1;i<=NF;i++){h=trim($i);if(h=="Kin_old")kc=i;if(h=="run_number")rc=i;if(h=="Type")tc=i;if(tolower(h)=="target")gc=i}next}
  NR>1&&kc&&rc&&tc&&gc{k=trim($kc);r=trim($rc);t=tolower(trim($tc));g=tolower(trim($gc));key=k "\t" r;
    if(index(kins,"|" k "|")&&t in wanted&&index(targets,"|" g "|")&&r~/^[0-9]+$/&&(runs=="|"||index(runs,"|" r "|"))&&!seen[key]++)print key "\t" g}
' "$config_source" | sort -k1,1 -k2,2n)
(( ${#job_pairs[@]} > 0 )) || die "No runs matched selected kinematics/types"

mkdir -p "$TAR_DIR"
SOURCE_REPORT="${TAR_DIR}/${WORKFLOW}_run_sources.csv"
[[ ! -e "$SOURCE_REPORT" ]] || die "Run-source report already exists: $SOURCE_REPORT"
source_report_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.run_sources.csv")"
printf 'kin,target,run,input_status,source,source_kind,input_file_count,input_bytes,input_paths\n' > "$source_report_tmp"

write_source_row() {
  local kin="$1" target="$2" run="$3" status="$4" source="$5" kind="$6" count="$7" bytes="$8" paths="$9"
  kin="${kin//\"/\"\"}"; target="${target//\"/\"\"}"; status="${status//\"/\"\"}"
  source="${source//\"/\"\"}"; kind="${kind//\"/\"\"}"; paths="${paths//\"/\"\"}"
  printf '"%s","%s",%s,"%s","%s","%s",%s,%s,"%s"\n' \
    "$kin" "$target" "$run" "$status" "$source" "$kind" "$count" "$bytes" "$paths" >> "$source_report_tmp"
}

declare -A JOB_INPUT_RECORDS=() JOB_INPUT_COUNT=() JOB_INPUT_BYTES=() JOB_RUN_COUNT=() JOB_RUNS=() DEST_SOURCE=()
declare -A RUN_INPUT_RECORDS=() RUN_INPUT_COUNT=() RUN_INPUT_BYTES=()
missing_input_count=0
ready_run_count=0
for pair in "${job_pairs[@]}"; do
  IFS=$'\t' read -r kin run target <<< "$pair"
  target="$(canonical_target "$target")"
  if ! resolve_run_inputs "$run"; then
    write_source_row "$kin" "$target" "$run" missing "${SOURCE:-auto}" "" 0 0 ""
    missing_input_count=$((missing_input_count + 1))
    continue
  fi

  run_records=(); run_paths=""; run_bytes=0; run_file_count=0; run_status=ready
  declare -A run_dest=()
  for file in "${RESOLVED_FILES[@]}"; do
    file="$(readlink -f "$file")"
    if [[ ! -r "$file" ]]; then run_status=unreadable; break; fi
    basename="$(basename "$file")"; key="${kin}|${run}|${basename}"
    if [[ -n "${DEST_SOURCE[$key]+x}" ]]; then
      if [[ "${DEST_SOURCE[$key]}" != "$file" ]]; then run_status=basename_collision; break; fi
      continue
    fi
    if [[ -n "${run_dest[$basename]+x}" ]]; then
      if [[ "${run_dest[$basename]}" != "$file" ]]; then run_status=basename_collision; break; fi
      continue
    fi
    if [[ "$INPUT_ACCESS" == stage ]]; then
      perms="$(stat -c %A "$file")"
      if [[ "${perms:7:1}" != r ]]; then run_status=not_world_readable; break; fi
    fi
    size="$(stat -c %s "$file")"
    run_dest["$basename"]="$file"
    run_records+=("${run}"$'\t'"${RESOLVED_KIND}"$'\t'"${basename}"$'\t'"${file}"$'\t'"${size}")
    run_paths+="${run_paths:+;}${file}"
    run_bytes=$((run_bytes + size)); run_file_count=$((run_file_count + 1))
  done

  if [[ "$run_status" != ready ]]; then
    write_source_row "$kin" "$target" "$run" "$run_status" "$RESOLVED_SOURCE" "$RESOLVED_KIND" 0 0 "$run_paths"
    missing_input_count=$((missing_input_count + 1))
    continue
  fi

  for record in "${run_records[@]}"; do
    IFS=$'\t' read -r record_run record_kind record_basename record_path record_size <<< "$record"
    DEST_SOURCE["${kin}|${run}|${record_basename}"]="$record_path"
    JOB_INPUT_RECORDS["$kin"]+="${record}"$'\n'
    RUN_INPUT_RECORDS["${kin}|${run}"]+="${record}"$'\n'
  done
  JOB_INPUT_COUNT["$kin"]=$(( ${JOB_INPUT_COUNT[$kin]:-0} + run_file_count ))
  JOB_INPUT_BYTES["$kin"]=$(( ${JOB_INPUT_BYTES[$kin]:-0} + run_bytes ))
  JOB_RUN_COUNT["$kin"]=$(( ${JOB_RUN_COUNT[$kin]:-0} + 1 ))
  JOB_RUNS["$kin"]+="${run} "
  RUN_INPUT_COUNT["${kin}|${run}"]="$run_file_count"
  RUN_INPUT_BYTES["${kin}|${run}"]="$run_bytes"
  ready_run_count=$((ready_run_count + 1))
  write_source_row "$kin" "$target" "$run" ready "$RESOLVED_SOURCE" "$RESOLVED_KIND" "$run_file_count" "$run_bytes" "$run_paths"
done

mv "$source_report_tmp" "$SOURCE_REPORT"
chmod 0644 "$SOURCE_REPORT"
if (( missing_input_count > 0 )); then
  if [[ "$MISSING_INPUT" == skip ]]; then
    echo "[WARN] Skipping $missing_input_count run(s) without usable ROOT input; report: $SOURCE_REPORT" >&2
  else
    echo "[WARN] Found $missing_input_count run(s) without usable ROOT input; report: $SOURCE_REPORT" >&2
  fi
fi
if [[ "$MISSING_INPUT" == fail ]] && (( missing_input_count > 0 )); then
  die "Missing/unusable ROOT inputs found; see $SOURCE_REPORT"
fi
(( ready_run_count > 0 )) || die "No usable ROOT inputs; see $SOURCE_REPORT"

active_kins=()
for kin in "${KINS[@]}"; do
  if (( ${JOB_RUN_COUNT[$kin]:-0} > 0 )); then active_kins+=("$kin")
  else echo "[WARN] No usable inputs for $kin; no SWIF job will be created" >&2
  fi
done
KINS=("${active_kins[@]}")

size_to_bytes() {
  local value="$1" number="${1%[kmgt]}" suffix="${1:${#1}-1}" factor=1
  [[ "$suffix" =~ [kmgt] ]] || { printf '%s' "$value"; return; }
  case "$suffix" in k) factor=1024;; m) factor=$((1024**2));; g) factor=$((1024**3));; t) factor=$((1024**4));; esac
  printf '%s' "$((number * factor))"
}
if [[ "$INPUT_ACCESS" == stage ]]; then
  disk_bytes="$(size_to_bytes "$DISK")"
  for kin in "${KINS[@]}"; do
    read -r -a kin_runs <<< "${JOB_RUNS[$kin]}"
    for run in "${kin_runs[@]}"; do
      run_key="${kin}|${run}"
      bytes="${RUN_INPUT_BYTES[$run_key]:-0}"
      (( bytes * 5 <= disk_bytes * 4 )) || die "Staged inputs for $kin run $run use $bytes bytes; request more --disk (20% headroom required)"
    done
  done
fi

OUTPUT_DIR="${OUTPUT_DIR:-/volatile/hallc/nps/${SUBMIT_USER}/nps_analysis/${WORKFLOW}}"
LOG_DIR="/farm_out/${SUBMIT_USER}/nps_analysis/${WORKFLOW}"
[[ "$OUTPUT_DIR" == /volatile/* ]] || die "--output-dir must be under /volatile"
[[ "$LOG_DIR" == /farm_out/* ]] || die "Internal log directory must be under /farm_out"
runtime_root_name="$(basename "$ANALYSIS_ROOT")"
runtime_tar="${TAR_DIR}/${WORKFLOW}_analysis_runtime.tar.gz"
runtime_sha="${runtime_tar}.sha256"
worker_basename="${WORKFLOW}_worker.sh"
worker_script="${TAR_DIR}/${worker_basename}"

require_command tar
require_command sha256sum
mkdir -p "$TAR_DIR"
(( DRY_RUN == 1 )) || mkdir -p "$OUTPUT_DIR" "$LOG_DIR"
[[ ! -e "$runtime_tar" && ! -e "$runtime_sha" && ! -e "$worker_script" ]] || \
  die "Runtime package/worker already exists for workflow: $WORKFLOW"
declare -A RUN_MANIFEST=()
for kin in "${KINS[@]}"; do
  safe="$(safe_name "$kin")"
  result="${OUTPUT_DIR}/${safe}.tar.gz"
  [[ ! -e "$result" ]] || die "Result already exists; refusing overwrite: $result"
  [[ ! -e "${OUTPUT_DIR}/${safe}" ]] || die "Analysis directory already exists; refusing overwrite: ${OUTPUT_DIR}/${safe}"
  [[ -f "${ANALYSIS_ROOT}/output/efficiency_stuff/selection_report_${kin}.csv" ]] || \
    die "Selection-report CSV missing for $kin"
  [[ -f "${ANALYSIS_ROOT}/output/efficiency_stuff/efficiency_${kin}.csv" ]] || \
    die "Efficiency CSV missing for $kin"
  read -r -a kin_runs <<< "${JOB_RUNS[$kin]}"
  for run in "${kin_runs[@]}"; do
    run_key="${kin}|${run}"
    manifest="${TAR_DIR}/${WORKFLOW}_${safe}_run${run}_inputs.tsv"
    [[ ! -e "$manifest" ]] || die "Input manifest already exists: $manifest"
    manifest_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}_${safe}_run${run}.XXXXXX.tsv")"
    printf 'run\tsource_kind\tbasename\tsource_path\tsize_bytes\n%s' "${RUN_INPUT_RECORDS[$run_key]}" > "$manifest_tmp"
    mv "$manifest_tmp" "$manifest"
    chmod 0644 "$manifest"
    RUN_MANIFEST["$run_key"]="$manifest"
  done
done

# Atomic package creation. Keep exact live config, analysis code, and only two
# required efficiency CSV inputs per selected kin; exclude generated products.
tar_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.tar.gz")"
cleanup() {
  [[ -n "${tar_tmp:-}" && -f "$tar_tmp" ]] && rm -f "$tar_tmp"
  [[ -n "${worker_tmp:-}" && -f "$worker_tmp" ]] && rm -f "$worker_tmp"
}
trap cleanup EXIT
archive_members=("${runtime_root_name}/config" "${runtime_root_name}/src/analysis")
for kin in "${KINS[@]}"; do
  archive_members+=(
    "${runtime_root_name}/output/efficiency_stuff/selection_report_${kin}.csv"
    "${runtime_root_name}/output/efficiency_stuff/efficiency_${kin}.csv"
  )
done
tar -C "$(dirname "$ANALYSIS_ROOT")" \
  --exclude="${runtime_root_name}/src/analysis/*.so" \
  --exclude="${runtime_root_name}/src/analysis/*.d" \
  --exclude="${runtime_root_name}/src/analysis/*.pcm" \
  --exclude="${runtime_root_name}/src/analysis/__pycache__" \
  -czf "$tar_tmp" "${archive_members[@]}"
tar -tzf "$tar_tmp" >/dev/null
tar -tzf "$tar_tmp" | grep -x "${runtime_root_name}/${RUNNER_REL}" >/dev/null || die "Generated tar lacks runner"
mv "$tar_tmp" "$runtime_tar"
tar_tmp=""
runtime_hash="$(sha256sum "$runtime_tar" | awk '{print $1}')"
printf '%s  analysis_runtime.tar.gz\n' "$runtime_hash" > "$runtime_sha"
chmod 0644 "$runtime_tar" "$runtime_sha"

# Freeze this exact submit/worker implementation under a workflow-unique source
# and staged name. SWIF may cache a repeatedly reused source path, which can run
# stale worker code after the live script changes.
worker_tmp="$(mktemp "${TAR_DIR}/.${WORKFLOW}.XXXXXX.worker.sh")"
cp "$SCRIPT_PATH" "$worker_tmp"
chmod 0755 "$worker_tmp"
mv "$worker_tmp" "$worker_script"
worker_hash="$(sha256sum "$worker_script" | awk '{print $1}')"

runner_options=(--mode "$MODE" --types "$TYPES_CSV" --target "${TARGETS[@]}" --gevnum-cut "$GEVNUM_CUT" --timeout "$TIMEOUT_SEC" \
  --updated-dir "$UPDATED_DIR" --production-dir "$PRODUCTION_DIR" --waveform-dir "$WAVEFORM_DIR")
finalize_options=("${runner_options[@]}")
if [[ -n "$COMBINE_TARGET" ]]; then
  runner_options+=(--combine-target "$COMBINE_TARGET")
  finalize_options+=(--combine-target "$COMBINE_TARGET")
fi
[[ -n "$SOURCE" ]] && runner_options+=(--source "$SOURCE")
[[ -n "$INPUT_DIR" ]] && runner_options+=(--input-dir "$INPUT_DIR")
if (( COMBINE == 0 )); then
  runner_options+=(--no-combine)
  finalize_options+=(--no-combine)
fi

echo "Workflow base: $WORKFLOW"
echo "Kinematics: ${KINS[*]}"
echo "Targets:     ${TARGETS[*]}"
echo "Runtime:     $runtime_tar"
echo "SHA256:      $(awk '{print $1}' "$runtime_sha")"
echo "Worker:      $worker_script"
echo "Worker SHA:  $worker_hash"
echo "Results:     $OUTPUT_DIR"
echo "Logs:        $LOG_DIR"
echo "Resources:   cores/run=$CORES ram/run=$RAM disk/run=$DISK time=$WALLTIME max-concurrent/kin=$MAX_CONCURRENT"
echo "ROOT input:  $INPUT_ACCESS"
echo "Run sources: $SOURCE_REPORT (ready=$ready_run_count skipped=$missing_input_count)"
for kin in "${KINS[@]}"; do echo "  ${WORKFLOW}_$(safe_name "$kin"): ${JOB_RUN_COUNT[$kin]} run job(s) + finalizer"; done

swif() {
  if (( DRY_RUN == 1 )); then printf '[dry-run]'; printf ' %q' swif2 "$@"; printf '\n';
  else swif2 "$@"; fi
}

(( DRY_RUN == 1 )) || require_command swif2
workflows=()
for kin in "${KINS[@]}"; do
  safe="$(safe_name "$kin")"
  kin_workflow="${WORKFLOW}_${safe}"
  workflows+=("$kin_workflow")
  result="${OUTPUT_DIR}/${safe}.tar.gz"
  read -r -a kin_runs <<< "${JOB_RUNS[$kin]}"
  (( DRY_RUN == 1 )) || mkdir -p "${LOG_DIR}/${safe}"
  swif create "$kin_workflow" -max-concurrent "$MAX_CONCURRENT"

  antecedents=()
  for run in "${kin_runs[@]}"; do
    run_key="${kin}|${run}"
    job_name="run_${run}"
    antecedents+=(-antecedent "$job_name")
    inputs=(-input analysis_runtime.tar.gz "file:${runtime_tar}" \
      -input analysis_runtime.tar.gz.sha256 "file:${runtime_sha}" \
      -input input_manifest.tsv "file:${RUN_MANIFEST[$run_key]}" \
      -input run_source_report.csv "file:${SOURCE_REPORT}" \
      -input "$worker_basename" "file:${worker_script}")
    if [[ "$INPUT_ACCESS" == stage ]]; then
      unset staged_dest
      declare -A staged_dest=()
      while IFS=$'\t' read -r manifest_run source_kind basename source_path size_bytes; do
        [[ "$manifest_run" == run ]] && continue
        [[ -n "$basename" ]] || continue
        [[ -n "${staged_dest[$basename]+x}" ]] && continue
        staged_dest["$basename"]="$source_path"
        inputs+=(-input "$basename" "file:${source_path}")
      done < "${RUN_MANIFEST[$run_key]}"
    fi
    swif add-job "$kin_workflow" \
      -name "$job_name" -account "$ACCOUNT" -partition "$PARTITION" \
      -cores "$CORES" -ram "$RAM" -disk "$DISK" -time "$WALLTIME" \
      -stdout "${LOG_DIR}/${safe}/${job_name}.out" -stderr "${LOG_DIR}/${safe}/${job_name}.err" \
      "${inputs[@]}" \
      /bin/bash "$worker_basename" --run-one \
        analysis_runtime.tar.gz "$runtime_root_name" "$kin" "$run" "$INPUT_ACCESS" "$OUTPUT_DIR" "$LOG_DIR" \
        "${runner_options[@]}"
  done

  final_inputs=(-input analysis_runtime.tar.gz "file:${runtime_tar}" \
    -input analysis_runtime.tar.gz.sha256 "file:${runtime_sha}" \
    -input run_source_report.csv "file:${SOURCE_REPORT}" \
    -input "$worker_basename" "file:${worker_script}")
  final_runner_options=("${finalize_options[@]}" --run "${kin_runs[@]}")
  swif add-job "$kin_workflow" \
    -name "finalize_${safe}" "${antecedents[@]}" \
    -account "$ACCOUNT" -partition "$PARTITION" \
    -cores 1 -ram "$RAM" -disk "$DISK" -time "$WALLTIME" \
    -stdout "${LOG_DIR}/${safe}/finalize.out" -stderr "${LOG_DIR}/${safe}/finalize.err" \
    "${final_inputs[@]}" \
    /bin/bash "$worker_basename" --finalize-one \
      analysis_runtime.tar.gz "$runtime_root_name" "$kin" "$result" "$OUTPUT_DIR" "$LOG_DIR" \
      "${final_runner_options[@]}"
done

if (( CREATE_ONLY == 0 )); then
  for kin_workflow in "${workflows[@]}"; do swif run "$kin_workflow"; done
fi
if (( DRY_RUN == 1 )); then
  echo "Dry run complete; no SWIF2 workflow state changed."
elif (( CREATE_ONLY == 1 )); then
  echo "Workflows created, not started:"
  for kin_workflow in "${workflows[@]}"; do echo "  swif2 run $kin_workflow"; done
else
  echo "Submitted workflows:"
  for kin_workflow in "${workflows[@]}"; do echo "  swif2 status $kin_workflow"; done
fi
